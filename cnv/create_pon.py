#!/usr/bin/env python3
"""Panel of Normals (PoN) for CNVscope.

Reads multiple normalized_outf TSV files, builds per-probe statistics
(dark probe flags, MAD-based variance, and panel medians of log2 DP/DP0,
always stored for valid autosomal probes). Optionally computes truncated
SVD noise components with --pon_svd. With --ratio_sex M|F, chrX/chrY probes
get same-sex panel medians at the two-copy level (a per-sex PoN for
cnvscope --pon_mode ratio).

Outputs:
  <prefix>       — binary PoN file (CNVPON v1: per-probe table + optional V_k matrices)
  <prefix>.tsv   — with --tsv only: the per-probe table as text, plus per-probe capture stats at
                   --capture_floor; for people to read, no tool reads it
  <prefix>.json  — provenance sidecar (not with --legacy)
  <prefix>.qc.{report.tsv,members.tsv,unstable.bed,json,txt}
                 — QC report of pon_qc.py (next to this file) on the binary PoN just written and the members
                   (--qc, the default, with --sample_sex; not with --legacy)

Builder version 2. Against version 1 (md5 27d01655, which --legacy reproduces byte for byte):
  - the probe index is the union of all members' probes, and the PoN does not depend on member order;
  - --min_samples auto (default), counting only members that carry probes;
  - warnings for small panels and for a --ratio_sex PoN without referenced chrX probes;
  - provenance #key=value lines in the .tsv header and the .json sidecar;
  - the QC report of pon_qc.py after the build; it never changes the .pon or the .pon.tsv.
The binary layout (CNVPON v1) is unchanged.
"""

import argparse
import hashlib
import json
import logging
import os
import platform
import resource
import shlex
import struct
import sys
import time
import warnings
import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

PON_VERSION = 1
BUILDER_VERSION = 2

_STUDY = ("the PoN-size study (wes2.1 ponsize/REPORT.md) measured, against a PoN of 106-115 normals, -0.3 (N=50), "
          "-0.4 (N=30) and -1.2 (N=20) points of strict autosomal macro F1, and recommends at least 30 normals, "
          "50 where possible")


def _min_samples_arg(v):
    if v == "auto":
        return v
    try:
        return int(v)
    except ValueError:
        raise argparse.ArgumentTypeError(f"takes 'auto' or an integer, not {v!r}")


def _md5_file(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _write_json(path, obj):
    with open(path, "w") as f:
        json.dump(obj, f, indent=1)
        f.write("\n")


def _build_qc(args, min_samples, has_sex):
    """pon_qc.run_build_qc on the PoN just written; returns its record for the .json sidecar. Never fails the build."""
    if not has_sex:
        why = "needs --sample_sex and --par_bed"
        logger.warning("QC skipped: the QC report %s", why)
        return {"status": f"skipped: {why}"}
    here = os.path.dirname(os.path.realpath(__file__))
    prefix = args.output + ".qc"
    rec = {"prefix": os.path.basename(prefix), "prior": args.qc_prior, "boot": 200, "jobs": args.qc_jobs or 1,
           "pon_qc": os.path.join(here, "pon_qc.py")}
    t = time.time()
    try:
        rec["pon_qc_md5"] = _md5_file(rec["pon_qc"])
        sys.dont_write_bytecode = True  # before the import: pon_qc's own guard cannot stop its .pyc
        if here not in sys.path:
            sys.path.insert(0, here)
        import pon_qc
        logger.info("QC: pon_qc.py (md5 %s) on %s, prior %s, %d bootstrap draws, %d jobs", rec["pon_qc_md5"],
                    args.output, args.qc_prior or "none", rec["boot"], rec["jobs"])
        st = pon_qc.run_build_qc(args.input_list, args.sample_sex, args.par_bed, args.output, min_samples,
                                 ratio_sex=args.ratio_sex, capture_floor=args.capture_floor, prior=args.qc_prior,
                                 out_prefix=prefix, boot=rec["boot"], jobs=rec["jobs"], verbose=args.verbose)
    except Exception as e:  # the PoN is written; a QC failure is reported, not fatal
        logger.error("QC failed after %.0f s: %s: %s", time.time() - t, type(e).__name__, e)
        rec.update(status=f"failed: {type(e).__name__}: {e}", seconds=round(time.time() - t, 1))
        return rec
    rec["seconds"] = round(time.time() - t, 1)
    rec["verdicts"] = {k: st.get(k) for k in ("autosomes", "chrXY_male", "chrX_female")}
    rec["files"] = [os.path.basename(prefix + x) for x in (".report.tsv", ".members.tsv", ".unstable.bed", ".json",
                                                            ".txt")]
    rec["status"] = "done"
    with open(prefix + ".txt") as f:
        for line in f.read().splitlines():
            logger.info("QC: %s", line)
    v = rec["verdicts"]
    bad = any(x in ("warn", "fail") for x in v.values())
    (logger.warning if bad else logger.info)(
        "QC verdicts: autosomes %s, chrX/chrY (male PoN) %s, chrX (female PoN) %s; %.0f s; report %s.*",
        *(str(v[k]).upper() for k in ("autosomes", "chrXY_male", "chrX_female")), rec["seconds"], prefix)
    return rec


def _option_string(parser, args):
    """Every option with its resolved value, in parser order."""
    parts = []
    for a in parser._actions:
        if not a.option_strings or a.dest == "help" or a.option_strings == ["--no_qc"]:
            continue
        v = getattr(args, a.dest)
        if a.dest == "qc":
            parts.append("--qc" if v else "--no_qc")
        elif isinstance(a, argparse._StoreTrueAction):
            if v:
                parts.append(a.option_strings[0])
        elif v is not None:
            parts.append(f"{a.option_strings[0]} {shlex.quote(str(v))}")
    return " ".join(parts)


# Pack (contig_code, start, end) into int64 for fast sort/unique/searchsorted.
# Supports contig_code < 9200, start < 10^9, probe_length < 10^6.
_KEY_C = np.int64(10**15)
_KEY_S = np.int64(10**6)


def _pack_keys(c, s, e):
    return c.astype(np.int64) * _KEY_C + s.astype(np.int64) * _KEY_S + (e - s).astype(np.int64)


def _unpack_keys(keys):
    c = (keys // _KEY_C).astype(np.uint16)
    rem = keys % _KEY_C
    s = (rem // _KEY_S).astype(np.uint32)
    d = (rem % _KEY_S).astype(np.uint32)
    return c, s, s + d


def _fast_col_median(mat, n_valid_per_col):
    """Fast column-wise median via transposed sort.

    Transpose makes each column contiguous (n_samples elements in L1 cache),
    then in-place sort puts NaN at end. Per-column n_valid selects correct
    median position — fully vectorized, no Python loop.
    """
    n_rows, n_cols = mat.shape
    if n_rows == 0:  # no samples (a sex absent from the panel): every median is NaN
        return np.full(n_cols, np.nan)
    mt = np.ascontiguousarray(mat.T)
    mt.sort(axis=1)  # NaN sorts to end; each row now has valid values first
    nv = n_valid_per_col
    if np.all(nv == n_rows):
        k = n_rows // 2
        if n_rows % 2 == 0:
            return (mt[:, k - 1] + mt[:, k]) / 2.0
        else:
            return mt[:, k].copy()
    k_lo = np.maximum((nv - 1) // 2, 0).astype(np.intp)
    k_hi = np.maximum(nv // 2, 0).astype(np.intp)
    idx = np.arange(n_cols)
    result = (mt[idx, k_lo] + mt[idx, k_hi]) / 2.0
    result[nv == 0] = np.nan
    return result


def parse_bed(fname):
    regions = {}
    with open(fname) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            contig, start, end = parts[0], int(parts[1]), int(parts[2])
            regions.setdefault(contig, []).append((start, end))
    for c in regions:
        regions[c].sort()
    return regions


def main():
    parser = argparse.ArgumentParser(
        description="Create a Panel of Normals for CNVscope (panel medians always; "
        "SVD components with --pon_svd)"
    )
    parser.add_argument(
        "--input_list",
        required=True,
        help="File with one normalized_outf path per line",
    )
    parser.add_argument(
        "--output", required=True,
        help="Output PoN file (--tsv also writes <output>.tsv; <output>.json is the provenance, not with --legacy)"
    )
    parser.add_argument(
        "--n_components", type=int, default=20, help="Number of SVD components (K)"
    )
    parser.add_argument(
        "--min_samples",
        type=_min_samples_arg,
        default="auto",
        help="Min valid samples per probe (on chrX/chrY with --ratio_sex: of the PoN's sex). 'auto' (default): "
        "min(20, max(3, floor(n_smaller_sex / 2))) with --sample_sex, else min(20, max(3, floor(N / 4))), counting "
        "only members with probes (a sex with none does not count). --legacy needs a number",
    )
    parser.add_argument(
        "--min_median_dp", type=float, default=0.05, help="Min median DP to keep probe"
    )
    parser.add_argument(
        "--max_sigma2",
        type=float,
        default=0.0,
        help="Exclude probes whose across-sample variance (sigma2_dp) exceeds this "
        "(0 = off). min_median_dp only catches probes that are low in MOST samples; "
        "a cohort with per-sample stochastic dropout has erratic-but-not-low probes "
        "that only the variance channel identifies.",
    )
    parser.add_argument(
        "--exclude_bed", default=None, help="BED file of regions to exclude"
    )
    parser.add_argument(
        "--pon_svd", action="store_true", help="Compute and store SVD components"
    )
    parser.add_argument(
        "--tsv", action="store_true",
        help="Also write <output>.tsv, the per-probe table as text with capture stats, for people to read; no tool "
        "reads it",
    )
    parser.add_argument(
        "--capture_floor",
        type=float,
        default=0.1,
        help="DP at or above which a sample counts as captured, for the TSV capture columns cap_n/cap_share/cap_level, "
        "and the panel median DP from which the QC uses a probe (default 0.1 = cnvscope's pon_min_median)",
    )
    parser.add_argument(
        "--verbose", action="store_true", help="Enable debug logging"
    )
    parser.add_argument(
        "--sample_sex",
        default=None,
        help="TSV file: <normalized_outf_path>\\t<M|F>. Enables sex-aware PoN (v4)",
    )
    parser.add_argument(
        "--par_bed",
        default=None,
        help="BED file of PAR regions. Required when --sample_sex is provided",
    )
    parser.add_argument(
        "--ratio_sex",
        choices=("M", "F"),
        default=None,
        help="Per-sex PoN for --pon_mode ratio: chrX/chrY probes referenced to this "
        "sex's panel samples at the two-copy level (male chrX outside the PARs and "
        "chrY: +1 in log2). Requires --sample_sex and --par_bed; not with --pon_svd",
    )
    parser.add_argument(
        "--qc", dest="qc", action="store_true", default=None,
        help="Run the QC report (pon_qc.py next to this file) on the PoN just written: <output>.qc.* (default; "
        "needs --sample_sex)",
    )
    parser.add_argument("--no_qc", dest="qc", action="store_false", help="Do not run the QC report")
    parser.add_argument(
        "--qc_prior", default=None,
        help="Population CNV prior VCF for the QC's common-locus annotation, e.g. prior/cnvprior.hprc178.vcf.gz, or "
        "BUNDLE/MEMBER, a VCF member of a model bundle, e.g. model/CNV_WES.bundle/cnvprior.vcf.gz; checked at start-up "
        "unless --no_qc (default: none)",
    )
    parser.add_argument(
        "--qc_jobs", type=int, default=None, help="Processes the QC uses to read the .norm files (default 1)"
    )
    parser.add_argument(
        "--legacy",
        action="store_true",
        help="Reproduce create_pon.py 27d01655 byte for byte: early-stopping probe index, a numeric "
        "--min_samples; no provenance lines, no .json",
    )
    args = parser.parse_args()
    t0 = time.time()
    legacy = args.legacy

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        stream=sys.stderr,
    )

    if legacy:
        if args.min_samples == "auto":
            logger.error("--legacy needs a numeric --min_samples (27d01655 has no auto)")
            sys.exit(1)
        if args.qc or args.qc_prior or args.qc_jobs:
            logger.error("--qc, --qc_prior and --qc_jobs are not available with --legacy")
            sys.exit(1)
        logger.info("--legacy: the create_pon.py 27d01655 build")
    else:
        if args.qc is None:
            args.qc = True
        builder_md5 = _md5_file(os.path.realpath(__file__))
        option_string = _option_string(parser, args)
        logger.info("create_pon.py version %d (md5 %s): %s", BUILDER_VERSION, builder_md5, option_string)
        if args.qc_prior and args.qc:  # its form, before any panel file is read: archive headers only, no VCF parse
            here = os.path.dirname(os.path.realpath(__file__))
            sys.dont_write_bytecode = True  # before the import: pon_qc's own guard cannot stop its .pyc
            if here not in sys.path:
                sys.path.insert(0, here)
            try:
                import pon_qc
            except ImportError:  # no pon_qc.py: the QC fails after the build and says so
                pass
            else:
                try:
                    pon_qc.prior_member(args.qc_prior)
                except (ValueError, OSError) as e:
                    logger.error("--qc_prior: %s", e)
                    sys.exit(1)

    if args.sample_sex and not args.par_bed:
        logger.error("--sample_sex requires --par_bed")
        sys.exit(1)
    if args.par_bed and not args.sample_sex:
        logger.error("--par_bed requires --sample_sex")
        sys.exit(1)
    if args.ratio_sex and not args.sample_sex:
        logger.error("--ratio_sex requires --sample_sex and --par_bed")
        sys.exit(1)
    if args.ratio_sex and args.pon_svd:
        logger.error("--ratio_sex cannot be combined with --pon_svd")
        sys.exit(1)

    # Parse sample sex annotations
    sample_sex_map = {}
    if args.sample_sex:
        with open(args.sample_sex) as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.strip().split("\t")
                if len(parts) >= 2 and parts[1] in ("M", "F"):
                    sample_sex_map[parts[0]] = parts[1]
        logger.info("Loaded sex annotations: %d M, %d F",
                     sum(1 for v in sample_sex_map.values() if v == "M"),
                     sum(1 for v in sample_sex_map.values() if v == "F"))
        if args.ratio_sex and not sample_sex_map:
            logger.error("--ratio_sex: no M/F annotations in %s", args.sample_sex)
            sys.exit(1)

    # Parse PAR BED
    par_regions = parse_bed(args.par_bed) if args.par_bed else None

    with open(args.input_list) as f:
        sample_files = [
            line.strip() for line in f if line.strip() and not line.startswith("#")
        ]
    n_samples = len(sample_files)
    if n_samples < 3:
        logger.error("Need at least 3 samples, got %d", n_samples)
        sys.exit(1)

    # Build per-sample sex array
    sample_sex_arr = []
    n_males = 0
    n_females = 0
    if sample_sex_map:
        for fname in sample_files:
            sex = sample_sex_map.get(fname)
            if sex is None:
                logger.error("No sex annotation for sample: %s", fname)
                sys.exit(1)
            sample_sex_arr.append(sex)
        n_males = sum(1 for s in sample_sex_arr if s == "M")
        n_females = sum(1 for s in sample_sex_arr if s == "F")
        if not n_males or not n_females:
            sx = "F" if n_males else "M"
            logger.warning("No %s members in the panel: a --ratio_sex %s PoN from it would have no chrX%s reference",
                           sx, sx, "/chrY" if sx == "M" else "")
        elif n_males < 3 or n_females < 3:
            logger.warning("Low sample count per sex: %d males, %d females (recommend >=3 each)",
                           n_males, n_females)

    # Read all files, converting to numpy immediately (one DataFrame in memory at a time).
    # First file discovers contigs; subsequent files only add new ones (rare for WGS).
    logger.info("Reading %d samples...", n_samples)
    sample_data = []
    contig_list = None
    contig_to_code = None
    n_first = None
    probes_identical = True
    has_probes = []

    for si, fname in enumerate(sample_files):
        logger.debug("  Reading sample %d/%d: %s", si + 1, n_samples, fname)
        df = pd.read_csv(
            fname, sep="\t", usecols=[0, 1, 2, 3, 6], header=0, dtype={0: str}
        )
        df.columns = ["contig", "start", "end", "dp_val", "dp0_val"]
        has_probes.append(len(df) > 0)
        if len(df) == 0 and legacy:
            logger.warning("No probes in %s: it adds no depth but still counts in n_samples", fname)
        elif len(df) == 0:
            logger.warning("No probes in %s: it adds no depth, still counts in the header's n_samples, "
                           "and not in --min_samples auto", fname)
        cat = df["contig"].astype("category")
        if si == 0:
            n_first = len(df)
            contig_list = sorted(cat.cat.categories.tolist())
            contig_to_code = {c: i for i, c in enumerate(contig_list)}
        else:
            if len(df) != n_first:
                probes_identical = False
            for fc in cat.cat.categories:
                if fc not in contig_to_code:
                    contig_to_code[fc] = len(contig_list)
                    contig_list.append(fc)
        code_map = np.array(
            [contig_to_code[c] for c in cat.cat.categories], dtype=np.uint16
        )
        sample_data.append((
            code_map[cat.cat.codes.values],
            df["start"].values.astype(np.uint32),
            df["end"].values.astype(np.uint32),
            df["dp_val"].values.copy(),
            df["dp0_val"].values.copy(),
        ))
        del df, cat

    warn_list = []

    def warn(msg, *a):
        warn_list.append(msg % a)
        logger.warning(msg, *a)

    if legacy:
        min_samples = args.min_samples
    else:
        # Contig codes in name order, so that the PoN does not depend on which member comes first
        canon = sorted(contig_list)
        if canon != contig_list:
            remap = np.array([canon.index(c) for c in contig_list], dtype=np.uint16)
            sample_data = [(remap[sd[0]],) + sd[1:] for sd in sample_data]
            contig_list = canon
            contig_to_code = {c: i for i, c in enumerate(contig_list)}

        n_with = int(sum(has_probes))
        n_males_p = sum(1 for s, h in zip(sample_sex_arr, has_probes) if h and s == "M")
        n_females_p = sum(1 for s, h in zip(sample_sex_arr, has_probes) if h and s == "F")
        if n_with < 3:
            logger.error("Need at least 3 members with probes, got %d", n_with)
            sys.exit(1)
        if args.min_samples != "auto":
            min_samples = args.min_samples
            min_samples_rule = "given"
        elif sample_sex_arr:
            n_small = min(n for n in (n_males_p, n_females_p) if n)  # a sex without members with probes does not count
            min_samples = min(20, max(3, n_small // 2))
            who = ("members of the smaller sex with probes" if n_males_p and n_females_p else
                   f"{'M' if n_males_p else 'F'} members with probes, the only sex with any")
            min_samples_rule = (f"auto: min(20, max(3, floor({n_small} / 2))), {n_small} = {who} ({n_males_p} M, "
                                f"{n_females_p} F)")
        else:
            min_samples = min(20, max(3, n_with // 4))
            min_samples_rule = f"auto: min(20, max(3, floor({n_with} / 4))), {n_with} = members with probes"
        logger.info("--min_samples %d (%s)", min_samples, min_samples_rule)
        if n_with < 30:
            warn("Small panel: %d members with probes (< 30); %s", n_with, _STUDY)
        if sample_sex_arr:
            for sx in [args.ratio_sex] if args.ratio_sex else ["M", "F"]:
                n_sx = n_males_p if sx == "M" else n_females_p
                if n_sx < 30:
                    warn("Small %s panel: %d %s members with probes (< 30) reference the chrX/chrY probes of this "
                         "sex, and the minimum of 30 applies to each sex; %s", sx, n_sx, sx, _STUDY)
                if args.ratio_sex and min_samples > n_sx:
                    warn("--min_samples %d exceeds the %d %s members with probes: every chrX/chrY probe of this "
                         "--ratio_sex %s PoN will be dark", min_samples, n_sx, sx, sx)

    # Build probe index using packed int64 keys (fast sort/unique on flat arrays)
    logger.info("Building probe index...")
    if not legacy:
        # Union of all members' probes; the fast path only when every member has the same probes in the same order
        keys0 = _pack_keys(*sample_data[0][:3])
        probes_identical = True
        all_keys = np.unique(keys0)
        for sd in sample_data[1:]:
            keys_i = _pack_keys(*sd[:3])
            if probes_identical and not np.array_equal(keys_i, keys0):
                probes_identical = False
            if len(all_keys):
                pos = np.minimum(np.searchsorted(all_keys, keys_i), len(all_keys) - 1)
                keys_i = keys_i[all_keys[pos] != keys_i]
            if len(keys_i):
                all_keys = np.union1d(all_keys, keys_i)
    if probes_identical:
        keys0 = _pack_keys(*sample_data[0][:3])
        sort_order = np.argsort(keys0)
        all_keys = keys0[sort_order]
        all_c = sample_data[0][0][sort_order]
        all_s = sample_data[0][1][sort_order]
        all_e = sample_data[0][2][sort_order]
    elif not legacy:
        all_c, all_s, all_e = _unpack_keys(all_keys)
    else:
        all_keys = np.unique(_pack_keys(*sample_data[0][:3]))
        for sd in sample_data[1:]:
            n_before = len(all_keys)
            keys_i = _pack_keys(*sd[:3])
            all_keys = np.unique(np.concatenate([all_keys, keys_i]))
            if len(all_keys) == n_before:
                break
        all_c, all_s, all_e = _unpack_keys(all_keys)

    # BED exclusion (vectorized: loop over BED regions, not probes)
    bed_regions = parse_bed(args.exclude_bed) if args.exclude_bed else None
    if bed_regions:
        keep = np.ones(len(all_keys), dtype=bool)
        for contig, regions in bed_regions.items():
            if contig not in contig_to_code:
                continue
            code = contig_to_code[contig]
            cmask = all_c == code
            if not np.any(cmask):
                continue
            starts = all_s[cmask]
            ends = all_e[cmask]
            indices = np.where(cmask)[0]
            for bs, be in regions:
                overlap = (starts < be) & (ends > bs)
                keep[indices[overlap]] = False
        n_excluded_bed = int(np.sum(~keep))
        all_keys = all_keys[keep]
        all_c = all_c[keep]
        all_s = all_s[keep]
        all_e = all_e[keep]
        probes_identical = False
        logger.info("Excluded %d probes by BED regions", n_excluded_bed)

    n_probes = len(all_keys)
    logger.info("Found %d probes", n_probes)
    if not legacy:
        logger.info("Probe index: the union of the %d members' probes (%s)", n_samples,
                    "identical in every member" if probes_identical else "they differ between members")

    # Align each sample to probe index using searchsorted on int64 keys
    logger.info("Aligning samples to probe index...")
    dp_matrix = np.full((n_samples, n_probes), np.nan, dtype=np.float64)
    dp0_matrix = np.full((n_samples, n_probes), np.nan, dtype=np.float64)

    if probes_identical:
        logger.debug("Fast path: identical probes across all samples")
        for si, (c, s, e, dp, dp0) in enumerate(sample_data):
            dp_matrix[si] = dp[sort_order]
            dp0_matrix[si] = dp0[sort_order]
    else:
        for si, (c, s, e, dp, dp0) in enumerate(sample_data):
            keys_i = _pack_keys(c, s, e)
            idx = np.searchsorted(all_keys, keys_i)
            idx_c = np.minimum(idx, n_probes - 1)
            valid = all_keys[idx_c] == keys_i
            dp_matrix[si, idx[valid]] = dp[valid]
            dp0_matrix[si, idx[valid]] = dp0[valid]
    del sample_data

    # Per-probe statistics via transposed partition (much faster than nanmedian).
    # Columns with NaN are fixed up individually.
    logger.info("Computing per-probe statistics...")
    n_valid = np.sum(~np.isnan(dp_matrix), axis=0)
    n_nan_probes = int(np.sum(n_valid < n_samples))
    logger.debug("  %d/%d probes have NaN values", n_nan_probes, n_probes)
    median_dp = _fast_col_median(dp_matrix, n_valid)

    # Capture stats (TSV only): samples with DP >= floor, and their median DP
    captured = dp_matrix >= args.capture_floor
    cap_n = np.sum(captured, axis=0)
    with np.errstate(invalid="ignore"):
        cap_share = cap_n / n_valid
    cap_level = _fast_col_median(np.where(captured, dp_matrix, np.nan), cap_n)
    del captured
    logger.info(
        "Capture (DP >= %g): %d probes captured in no sample, %d in under half",
        args.capture_floor,
        int(np.sum(cap_n == 0)),
        int(np.sum(cap_share < 0.5)),
    )

    # Mark excluded probes
    excluded = np.zeros(n_probes, dtype=np.int32)
    excluded[n_valid < min_samples] = 1
    excluded[median_dp < args.min_median_dp] = 1
    excluded[np.isnan(median_dp)] = 1

    valid_mask = excluded == 0
    n_valid_probes = int(np.sum(valid_mask))
    logger.info(
        "Valid probes: %d/%d (excluded: %d)",
        n_valid_probes,
        n_probes,
        n_probes - n_valid_probes,
    )

    if n_valid_probes == 0:
        logger.error("No valid probes remain after filtering")
        sys.exit(1)

    if n_samples < 30 and legacy:
        logger.warning(
            "Only %d samples — PoN quality may be limited (recommend >=30)", n_samples
        )

    pseudocount = 1e-6

    # Per-probe robust variance (MAD-based) of log2(DP) across samples
    log2_dp = np.empty_like(dp_matrix)
    np.add(dp_matrix, pseudocount, out=log2_dp)
    np.log2(log2_dp, out=log2_dp)
    probe_medians = _fast_col_median(log2_dp, n_valid)
    log2_dp -= probe_medians
    np.abs(log2_dp, out=log2_dp)
    mad = _fast_col_median(log2_dp, n_valid)
    del log2_dp
    sigma2_dp_full = ((mad * 1.4826) ** 2).astype(np.float32)
    sigma2_dp_full[excluded == 1] = 0
    sigma2_dp_full = np.nan_to_num(sigma2_dp_full, nan=0.0)
    n_high_var = int(np.sum(sigma2_dp_full[valid_mask] > 0.5))
    logger.info(
        "Per-probe variance: median=%.4f, high-variance (>0.5): %d",
        np.median(sigma2_dp_full[valid_mask]),
        n_high_var,
    )

    # Optional second pass: drop erratic probes. median_dp only sees probes that are
    # low in MOST samples, so a cohort whose probes fail stochastically per sample
    # keeps them all; their signature is high across-sample variance, not a low median.
    if args.max_sigma2 > 0:
        newly = (sigma2_dp_full > args.max_sigma2) & (excluded == 0)
        n_newly = int(np.sum(newly))
        excluded[newly] = 1
        sigma2_dp_full[newly] = 0
        valid_mask = excluded == 0
        n_valid_probes = int(np.sum(valid_mask))
        logger.info(
            "max_sigma2=%.3f excluded %d further probes; valid now %d/%d "
            "(variance median=%.4f)",
            args.max_sigma2,
            n_newly,
            n_valid_probes,
            n_probes,
            np.median(sigma2_dp_full[valid_mask]) if n_valid_probes else float("nan"),
        )

    # Sex-specific per-probe statistics for sex chromosomes
    median_dp_male_full = np.full(n_probes, np.nan, dtype=np.float32)
    median_dp_female_full = np.full(n_probes, np.nan, dtype=np.float32)
    sigma2_dp_male_full = np.zeros(n_probes, dtype=np.float32)
    sigma2_dp_female_full = np.zeros(n_probes, dtype=np.float32)
    par_flag = np.zeros(n_probes, dtype=np.uint8)

    if sample_sex_arr:
        # integer even when empty: an absent sex selects no rows and keeps NaN medians and 0 sigma2
        male_idx = np.array([i for i, s in enumerate(sample_sex_arr) if s == "M"], dtype=np.intp)
        female_idx = np.array([i for i, s in enumerate(sample_sex_arr) if s == "F"], dtype=np.intp)

        # Identify sex chromosome probes
        sex_chroms = {"chrX", "X", "chrY", "Y"}
        sex_codes = np.array(
            [contig_to_code[c] for c in sex_chroms if c in contig_to_code], dtype=np.uint16
        )
        sex_mask = np.isin(all_c, sex_codes)

        # Mark PAR probes
        if par_regions:
            for contig, regions in par_regions.items():
                if contig not in contig_to_code:
                    continue
                code = contig_to_code[contig]
                cmask = all_c == code
                if not np.any(cmask):
                    continue
                starts = all_s[cmask]
                ends = all_e[cmask]
                indices = np.where(cmask)[0]
                for bs, be in regions:
                    overlap = (starts < be) & (ends > bs)
                    par_flag[indices[overlap]] = 1
            n_par = int(np.sum(par_flag))
            logger.info("PAR probes: %d", n_par)

        if np.any(sex_mask):
            sex_indices = np.where(sex_mask)[0]
            dp_sex_male = dp_matrix[np.ix_(male_idx, sex_indices)]
            dp_sex_female = dp_matrix[np.ix_(female_idx, sex_indices)]

            n_valid_male = np.sum(~np.isnan(dp_sex_male), axis=0)
            n_valid_female = np.sum(~np.isnan(dp_sex_female), axis=0)

            # Male medians and variance
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                log2_male = np.log2(dp_sex_male + pseudocount)
                med_male = _fast_col_median(log2_male, n_valid_male)
                resid_male = np.abs(log2_male - med_male)
                mad_male = _fast_col_median(resid_male, n_valid_male)
                s2_male = ((mad_male * 1.4826) ** 2).astype(np.float32)
                s2_male = np.nan_to_num(s2_male, nan=0.0)

            # Female medians and variance
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                log2_female = np.log2(dp_sex_female + pseudocount)
                med_female = _fast_col_median(log2_female, n_valid_female)
                resid_female = np.abs(log2_female - med_female)
                mad_female = _fast_col_median(resid_female, n_valid_female)
                s2_female = ((mad_female * 1.4826) ** 2).astype(np.float32)
                s2_female = np.nan_to_num(s2_female, nan=0.0)

            median_dp_male_full[sex_indices] = med_male
            median_dp_female_full[sex_indices] = med_female
            sigma2_dp_male_full[sex_indices] = s2_male
            sigma2_dp_female_full[sex_indices] = s2_female

            # Log sex-specific stats
            chrx_codes = [contig_to_code[c] for c in ("chrX", "X") if c in contig_to_code]
            chry_codes = [contig_to_code[c] for c in ("chrY", "Y") if c in contig_to_code]
            for label, codes in [("chrX", chrx_codes), ("chrY", chry_codes)]:
                if not codes:
                    continue
                cmask = np.isin(all_c[sex_indices], codes)
                if not np.any(cmask):
                    continue
                nonpar = cmask & (par_flag[sex_indices] == 0)
                if np.any(nonpar):
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore", RuntimeWarning)  # All-NaN where a sex has no probes
                        sx_parts = [f"no {lab} members" if not len(ix) else "%s median=%.4f (sigma2=%.4f)" % (
                            lab, np.nanmedian(med[nonpar]), np.nanmedian(s2[nonpar]))
                            for lab, ix, med, s2 in (("male", male_idx, med_male, s2_male),
                                                     ("female", female_idx, med_female, s2_female))]
                    logger.info("  %s non-PAR: %s, %s, mixed sigma2=%.4f", label, *sx_parts,
                                np.median(sigma2_dp_full[sex_indices[nonpar]]))
            del dp_sex_male, dp_sex_female, log2_male, log2_female

    # svd_idx: -1 = excluded (dark probe), -2 = non-autosomal, >=0 = valid autosomal
    # (SVD index with --pon_svd, else 0; with --ratio_sex also valid chrX/chrY)
    autosomes = {f"chr{i}" for i in range(1, 23)} | {str(i) for i in range(1, 23)}
    auto_codes = np.array(
        [contig_to_code[c] for c in autosomes if c in contig_to_code], dtype=np.uint16
    )
    auto_mask = np.isin(all_c, auto_codes)
    svd_idx = np.full(n_probes, -1, dtype=np.int32)
    svd_idx[~auto_mask] = -2
    svd_idx[valid_mask & auto_mask] = 0

    # Mask chrY PAR probes as dark. chrY PAR sequence is identical to chrX PAR;
    # standard aligners deposit those reads on chrX (or assign MAPQ=0), leaving
    # chrY PAR with no usable read-depth signal. Treat them like dark probes so
    # they are masked at inference and excluded from ploidy-block construction.
    if sample_sex_arr and par_regions:
        chry_codes_list = [contig_to_code[c] for c in ("chrY", "Y") if c in contig_to_code]
        if chry_codes_list:
            chry_codes_arr = np.array(chry_codes_list, dtype=np.uint16)
            chry_par_mask = par_flag.astype(bool) & np.isin(all_c, chry_codes_arr)
            n_chry_par = int(np.sum(chry_par_mask))
            if n_chry_par > 0:
                svd_idx[chry_par_mask] = -1
                par_flag[chry_par_mask] = 0
                logger.info(
                    "Masked %d chrY PAR probes as dark (reads multimap to chrX PAR)",
                    n_chry_par,
                )

    k = 0
    k_signal = 0
    n_svd_probes = 0
    median_variance = 0.0
    median_dp_full = np.zeros(n_probes, dtype=np.float32)
    median_dp0_full = np.zeros(n_probes, dtype=np.float32)

    # Panel medians of valid autosomal probes, stored with or without --pon_svd
    ref_mask = valid_mask & auto_mask
    n_ref_probes = int(np.sum(ref_mask))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        median_log2dp = np.nanmedian(np.log2(dp_matrix[:, ref_mask] + pseudocount), axis=0)
        median_log2dp0 = np.nanmedian(np.log2(dp0_matrix[:, ref_mask] + pseudocount), axis=0)
    if n_ref_probes:
        iqr = float(np.percentile(median_log2dp, 75) - np.percentile(median_log2dp, 25))
        median_variance = (iqr / 1.349) ** 2
        logger.debug(
            "Median variance (DP, IQR-based): %.6f (IQR=%.4f)", median_variance, iqr
        )
    else:
        logger.warning("No valid autosomal probes: no panel medians to store")
    median_dp_full[ref_mask] = median_log2dp
    median_dp0_full[ref_mask] = median_log2dp0
    logger.info("Panel medians stored for %d valid autosomal probes", n_ref_probes)

    # --ratio_sex: chrX/chrY probes (chrY PAR stays dark) referenced to same-sex panel samples
    if args.ratio_sex:
        xy_idx = np.flatnonzero(sex_mask)
        same_idx = male_idx if args.ratio_sex == "M" else female_idx
        dp_xy = dp_matrix[np.ix_(same_idx, xy_idx)]
        n_xy = np.sum(~np.isnan(dp_xy), axis=0)
        s2_xy = (sigma2_dp_male_full if args.ratio_sex == "M" else sigma2_dp_female_full)[xy_idx]
        y_codes = [contig_to_code[c] for c in ("chrY", "Y") if c in contig_to_code]
        is_y = np.isin(all_c[xy_idx], y_codes)
        ok = (svd_idx[xy_idx] == -2) & (n_xy >= min_samples)
        ok &= _fast_col_median(dp_xy, n_xy) >= args.min_median_dp
        if args.max_sigma2 > 0:
            ok &= s2_xy <= args.max_sigma2
        if args.ratio_sex == "F":
            ok &= ~is_y  # female chrY: dark
        # one-copy male regions (chrX outside the PARs, chrY) shifted to the two-copy level
        shift = np.where(par_flag[xy_idx] == 0, 1.0, 0.0) if args.ratio_sex == "M" else 0.0
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            med_xy = np.nanmedian(np.log2(dp_xy + pseudocount), axis=0) + shift
            med0_xy = np.nanmedian(
                np.log2(dp0_matrix[np.ix_(same_idx, xy_idx)] + pseudocount), axis=0
            ) + shift
        svd_idx[xy_idx] = np.where(ok, 0, -1)
        excluded[xy_idx] = ~ok
        median_dp_full[xy_idx[ok]] = med_xy[ok]
        median_dp0_full[xy_idx[ok]] = med0_xy[ok]
        sigma2_dp_full[xy_idx] = np.where(ok, s2_xy, 0)
        n_valid_probes = int(np.sum(excluded == 0))
        for label, cmask in (("chrX", ~is_y), ("chrY", is_y)):
            logger.info(
                "--ratio_sex %s: %s %d probes referenced, %d dark",
                args.ratio_sex,
                label,
                int(np.sum(ok & cmask)),
                int(np.sum(~ok & cmask)),
            )
        n_xref = {"chrX": int(np.sum(ok & ~is_y)), "chrY": int(np.sum(ok & is_y))}
        if not legacy and np.any(~is_y) and n_xref["chrX"] == 0:
            n_sx = n_males_p if args.ratio_sex == "M" else n_females_p
            bar = "!" * 100
            logger.warning(bar)
            warn("--ratio_sex %s: 0 of %d chrX probes referenced: EVERY chrX PROBE OF THIS PoN IS DARK and its "
                 "samples get no chrX calls. The panel has %d %s members with probes against --min_samples %d; "
                 "add members of this sex or lower --min_samples.", args.ratio_sex, int(np.sum(~is_y)), n_sx,
                 args.ratio_sex, min_samples)
            logger.warning(bar)

    if args.pon_svd:
        svd_mask = ref_mask
        n_svd_probes = n_ref_probes
        n_non_auto = int(np.sum(~auto_mask))
        logger.info(
            "SVD-eligible probes: %d (excluded %d non-autosomal probes)",
            n_svd_probes,
            n_non_auto,
        )
        svd_idx[svd_mask] = np.arange(n_svd_probes)

        k_max = min(args.n_components, n_samples - 1, n_svd_probes - 1)
        if k_max < 1:
            logger.error("Not enough samples/probes for SVD")
            sys.exit(1)

        def compute_svd(matrix, median_log2, label):
            M = matrix[:, svd_mask]
            L = np.log2(M + pseudocount)
            L_centered = L - median_log2
            L_centered = np.nan_to_num(L_centered, nan=0.0)
            logger.debug("  Running SVD (k_max=%d) for %s...", k_max, label)
            from scipy.sparse.linalg import svds

            np.random.seed(42)
            U, S, Vk = svds(L_centered, k=k_max)
            order = np.argsort(-S)
            S = S[order]
            Vk = Vk[order]
            total_var = np.sum(L_centered**2)
            explained = (S**2) / total_var
            cumulative = np.cumsum(explained)
            for i, e in enumerate(explained):
                logger.debug(
                    "    component %d: %.2f%% (cumulative: %.1f%%)",
                    i + 1,
                    e * 100,
                    cumulative[i] * 100,
                )
            return Vk, explained

        Vk_dp, explained_dp = compute_svd(dp_matrix, median_log2dp, "DP")
        Vk_dp0, _ = compute_svd(dp0_matrix, median_log2dp0, "DP0")

        k = k_max
        threshold = 3.0 / n_samples
        k_signal = max(1, int(np.sum(explained_dp > threshold)))
        logger.info(
            "Storing K=%d components (signal: %d, threshold=%.1f%%)",
            k,
            k_signal,
            threshold * 100,
        )

    # Build contig table (only contigs present in probes)
    used_codes = np.unique(all_c)
    contigs = [contig_list[c] for c in used_codes]
    code_to_cidx = np.zeros(len(contig_list), dtype=np.uint16)
    for i, c in enumerate(used_codes):
        code_to_cidx[c] = i

    # Write binary PoN file
    # Single canonical format (v1). Sex-aware fields are always present;
    # they carry zeros/NaN when --sample_sex was not provided at build time.
    pon_version = PON_VERSION
    out_fname = args.output
    logger.info("Writing %s (v%d)...", out_fname, pon_version)
    with open(out_fname, "wb") as f:
        # Main header: 32 bytes (offsets 0..31)
        header = b"CNVPON"                              # 0..5  : magic
        header += struct.pack("<H", pon_version)        # 6..7  : version
        header += struct.pack("<I", n_svd_probes)       # 8..11 : n_svd_probes
        header += struct.pack("<H", k)                  # 12..13: k_dp
        header += struct.pack("<H", k)                  # 14..15: k_dp0
        header += struct.pack("<I", n_probes)           # 16..19: n_total_probes
        header += struct.pack("<H", len(contigs))       # 20..21: n_contigs
        header += struct.pack("<H", n_samples)          # 22..23: n_samples
        header += struct.pack("<f", median_variance)    # 24..27: median_variance
        header += struct.pack("<H", k_signal)           # 28..29: k_signal
        header += b"\x00" * (32 - len(header))          # 30..31: reserved
        f.write(header)

        # Extended header: 32 bytes (offsets 32..63), uint16 sex counts + reserved
        extended = struct.pack("<H", n_males)           # 32..33: n_males uint16
        extended += struct.pack("<H", n_females)        # 34..35: n_females uint16
        extended += b"\x00" * (32 - len(extended))      # 36..63: reserved
        f.write(extended)

        # Contig table: (uint16 len, chars) per contig
        for c in contigs:
            cb = c.encode("utf-8")
            f.write(struct.pack("<H", len(cb)))
            f.write(cb)

        # Per-probe table: 43 bytes/probe
        probe_dt = np.dtype([
            ("cidx", "<u2"), ("start", "<u4"), ("end", "<u4"),
            ("median_dp", "<f4"), ("median_dp0", "<f4"),
            ("svd_idx", "<i4"), ("sigma2_dp", "<f4"),
            ("median_dp_male", "<f4"), ("median_dp_female", "<f4"),
            ("sigma2_dp_male", "<f4"), ("sigma2_dp_female", "<f4"),
            ("flags", "u1"),
        ])
        probe_table = np.zeros(n_probes, dtype=probe_dt)
        probe_table["cidx"] = code_to_cidx[all_c]
        probe_table["start"] = all_s
        probe_table["end"] = all_e
        probe_table["median_dp"] = median_dp_full
        probe_table["median_dp0"] = median_dp0_full
        probe_table["svd_idx"] = svd_idx
        probe_table["sigma2_dp"] = sigma2_dp_full
        probe_table["median_dp_male"] = median_dp_male_full
        probe_table["median_dp_female"] = median_dp_female_full
        probe_table["sigma2_dp_male"] = sigma2_dp_male_full
        probe_table["sigma2_dp_female"] = sigma2_dp_female_full
        probe_table["flags"] = par_flag
        f.write(probe_table.tobytes())

        # V_k matrices (only if --pon_svd)
        if k > 0:
            f.write(Vk_dp.astype(np.float32).tobytes())
            f.write(Vk_dp0.astype(np.float32).tobytes())

    # Provenance for the .tsv header and the .json sidecar; no timestamps, so that builds are reproducible
    if not legacy:
        def md5_text(text):
            return hashlib.md5(text.encode("utf-8")).hexdigest()

        members = sorted((os.path.basename(fn), sample_sex_arr[i] if sample_sex_arr else ".")
                         for i, fn in enumerate(sample_files))
        order = np.lexsort((all_e, all_s, all_c))  # contig codes follow the name order
        names = np.array(contig_list, dtype=object)[all_c[order]]
        probe_text = "".join(f"{c}\t{s}\t{e}\n" for c, s, e in
                             zip(names, all_s[order].tolist(), all_e[order].tolist()))
        empty = [os.path.basename(fn) for fn, h in zip(sample_files, has_probes) if not h]
        prov = {
            "builder_version": BUILDER_VERSION,
            "builder_md5": builder_md5,
            "builder_options": option_string,
            "member_list_md5": md5_text("".join(f"{b}\t{s}\n" for b, s in members)),
            "n_members_with_probes": n_with,
            "n_males_with_probes": n_males_p,
            "n_females_with_probes": n_females_p,
            "empty_members": ",".join(empty),
            "probe_set_md5": md5_text(probe_text),
            "n_probes": n_probes,
            "min_samples": min_samples,
            "min_samples_rule": min_samples_rule,
            "qc_requested": "yes" if args.qc else "no",
            "n_warnings": len(warn_list),
            "pon_md5": _md5_file(out_fname),
            "runtime": f"python {platform.python_version()}, numpy {np.__version__}, pandas {pd.__version__}",
        }
        del probe_text

    if args.tsv:
        tsv_fname = args.output + ".tsv"
        logger.info("Writing %s...", tsv_fname)
        contig_names = np.array(contig_list)[all_c]
        tsv_dict = {
            "contig": contig_names,
            "start": all_s,
            "end": all_e,
            "median_log2dp": median_dp_full,
            "median_log2dp0": median_dp0_full,
            "n_valid": n_valid.astype(int),
            "excluded": excluded,
            "probe_idx": svd_idx,
            "sigma2_dp": sigma2_dp_full,
            "median_dp_male": median_dp_male_full,
            "median_dp_female": median_dp_female_full,
            "sigma2_dp_male": sigma2_dp_male_full,
            "sigma2_dp_female": sigma2_dp_female_full,
            "par_flag": par_flag,
            "panel_median_dp": median_dp,
            "cap_n": cap_n,
            "cap_share": cap_share,
            "cap_level": cap_level,
        }
        tsv_df = pd.DataFrame(tsv_dict)
        with open(tsv_fname, "w") as f:
            f.write(f"#pon_version={pon_version}\n")
            f.write(f"#n_samples={n_samples}\n")
            f.write(f"#n_males={n_males}\n")
            f.write(f"#n_females={n_females}\n")
            f.write(f"#n_components={k}\n")
            f.write(f"#median_variance={median_variance:.6f}\n")
            f.write(f"#capture_floor={args.capture_floor:g}\n")
            if args.ratio_sex:
                f.write(f"#ratio_sex={args.ratio_sex}\n")
            if not legacy:
                for key, v in prov.items():
                    f.write(f"#{key}={v}\n")
            tsv_df.to_csv(f, sep="\t", index=False, float_format="%.6f")

    if not legacy:
        side = dict(prov)
        side["empty_members"] = empty
        side.update({
            "pon": os.path.basename(out_fname),
            "pon_version": pon_version,
            "tsv": os.path.basename(tsv_fname) if args.tsv else None,
            "tsv_md5": _md5_file(tsv_fname) if args.tsv else None,
            "n_samples": n_samples,
            "n_males": n_males,
            "n_females": n_females,
            "n_valid_probes": n_valid_probes,
            "ratio_sex": args.ratio_sex,
            "referenced_xy": n_xref if args.ratio_sex else None,
            "warnings": warn_list,
        })
        side["qc"] = {"status": "started" if args.qc else "not requested"}
        json_fname = args.output + ".json"
        _write_json(json_fname, side)
        logger.info("Wrote %s", json_fname)
        if args.qc:
            del dp_matrix, dp0_matrix
            if args.tsv:
                del tsv_df
            side["qc"] = _build_qc(args, min_samples, bool(sample_sex_arr))
            _write_json(json_fname, side)

    logger.info(
        "Done. %d valid probes, %d excluded, %d samples",
        n_valid_probes,
        n_probes - n_valid_probes,
        n_samples,
    )
    if k > 0:
        logger.info(
            "  SVD: %d probes, %d components (signal: %d)", n_svd_probes, k, k_signal
        )

    t1 = time.time()
    mm, ut, st = 0, 0, 0
    mem_scale = 1 if platform.system() == 'Darwin' else 1024
    for who in (resource.RUSAGE_SELF, resource.RUSAGE_CHILDREN):
        ru = resource.getrusage(who)
        mm += ru.ru_maxrss * mem_scale
        ut += ru.ru_utime
        st += ru.ru_stime
    logger.info('overall: %d mem %.3f user %.3f sys %.3f real', mm, ut, st, t1 - t0)


if __name__ == "__main__":
    main()
