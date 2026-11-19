#!/usr/bin/env python3
"""pon_qc.py -- quality report for a wes2.1 CNVscope panel of normals (PoN).

usage:
  python3 -I pon_qc.py --input_list depth.<c>.list --sample_sex depth.<c>.sex.list --par_bed config/hs38.PAR.bed \\
      --pon cnv_depth.<c>.male.pon.tsv [cnv_depth.<c>.female.pon.tsv] --min_samples N --out PREFIX \\
      [--prior prior/cnvprior.hprc178.vcf.gz] [--test S.norm ...] [--admit S.norm ... [--admit_sex M F ...]] \\
      [--compare OLD.male.pon.tsv [OLD.female.pon.tsv]] [--bootstrap 200] [--seed 20261110] [--jobs 4]
  python3 -I pon_qc.py --pon NEW.male.pon.tsv [NEW.female.pon.tsv] --compare OLD.male.pon.tsv [...] --out PREFIX

Inputs (--input_list members, --test, --admit) are raw .norm files: CNVscope --output_normalized from a run without
--pon_file (as wes1/data/rd_stats). The .norm that scripts/segment.sh writes is divided by the PoN (masked probes = 1);
such files are flagged pon_divided (> 1% of probes exactly 1.0; raw files <= 0.007%, PoN-divided 8-9%).
The panel inputs are create_pon.py's; each --pon is a --tsv table the builder wrote from that panel. The autosomal
checks read the first table (male and female tables share their autosomal rows), the chrX/chrY counts every table.
Transform and estimators are the builder's: log2(DP + 1e-6), column medians (mean of the two middle values),
sigma2 = (1.4826 MAD)^2.
--prior is the prior .vcf.gz, or BUNDLE/MEMBER, a member of a model bundle (an ar archive), e.g.
CNV_WES.bundle/cnvprior.vcf.gz; the member is copied to a temp file ($TMPDIR), read and removed (no .tbi needed).
PREFIX.json records the path, the member and the md5 of the VCF bytes ("prior").

Checks (thresholds in T, calibrated in ponsize/analysis/qc/ by pon_qc_calib.py)
  0 verification   this script's full-panel median and sigma2 against the stored ones (valid autosomal probes)
  1 composition    N, sexes, empty members, chrX/chrY counts, --min_samples against the builder's auto rule
                   min(20, max(3, floor(n_smaller_sex/2))): members with probes; a sex with none does not count
  2 precision      per-probe SE of the panel median: bootstrap over members (B=200) and 1.2533 sigma / sqrt(n);
                   status on the median of SE/sigma (1.2533/sqrt(N) for a clean panel)
  3 unstable loci  probes with bootstrap SE > 0.2, or with >= 10% of members deviating on one side, merged into
                   regions with the prior's DEL/DUP AF (BED)
  4 members        leave-one-out (LOO) noise, CNV burden (5-probe running median), sex, missing probes
  5 batch          PCA of the noise-scaled LOO profiles on 20k probes; PC1 above a permutation null, 2-means silhouette
  6 test samples   a sample's noise against the panel median, as a percentile and a multiple of the members' LOO noise
  7 admission      --admit: is a called sample fit to join the pool? The member checks against the current panel;
                   reject (does not match the panel), review (a member flag would fire) or admit
  8 compare        --compare: the stored medians and sigma2 of this PoN against a previous version
Verdicts: autosomes (N, relative SE, empty members, --min_samples, member mismatch/outliers, batch); chrX/chrY of the
male PoN and chrX of the female PoN (members of that sex with probes, its chrX_ref, --min_samples, sex flags); na for
a sex with no members with probes when no table of that sex is given.
Noise, burden, PCA and unstable flags use the valid autosomal probes CNVscope does not mask (panel median DP >=
#capture_floor, 0.1); at this setting the test noise equals CNVscope's measured noise.
Writes PREFIX.report.tsv, PREFIX.members.tsv, PREFIX.unstable.bed, PREFIX.tests.tsv (--test), PREFIX.admit.tsv and
PREFIX.admit_runs.tsv (--admit), PREFIX.compare.tsv and PREFIX.compare.bed (--compare), PREFIX.json and PREFIX.txt.
"""
import sys

sys.dont_write_bytecode = True                   # no __pycache__ next to the package

import argparse  # noqa: E402
import gzip  # noqa: E402
import hashlib  # noqa: E402
import json  # noqa: E402
import multiprocessing as mp  # noqa: E402
import os  # noqa: E402
import re  # noqa: E402
import tempfile  # noqa: E402
import time  # noqa: E402
import warnings  # noqa: E402
import zlib  # noqa: E402

import numpy as np
import pandas as pd

PSEUDO = 1e-6
MAD_K = 1.4826
MED_K = 1.2533                       # sqrt(pi/2): SE of a Gaussian median is 1.2533 sigma / sqrt(n)
AUTO = {f'chr{i}' for i in range(1, 23)} | {str(i) for i in range(1, 23)}
XN, YN = ('chrX', 'X'), ('chrY', 'Y')
_KC, _KS = np.int64(10 ** 15), np.int64(10 ** 6)

# Thresholds (status: pass / warn / fail). Calibration: ponsize/analysis/qc/ (pon_qc_calib.py).
T = {
    'n_pass': 50, 'n_warn': 30,              # members behind every autosomal probe
    'sex_pass': 50, 'sex_warn': 30,          # members of one sex behind its PoN's chrX/chrY: >= 50 pass, else warn
                                             # (30-49 below the recommended count, < 30 below the floor); fail only
                                             # when chrX_ref = 0 or the sex has fewer members than --min_samples
    'rel_se_pass': 0.19, 'rel_se_warn': 0.25,  # median bootstrap SE / sigma (1.2533/sqrt(N) on a clean panel):
                                             # N=50 0.178-0.180, N=30 0.235-0.238, N=20 0.293-0.299
    'se_unstable': 0.2,                      # probe flag: bootstrap SE of the median, log2
    'bimodal_dev': 0.5, 'bimodal_frac': 0.10,  # probe flag: >= 10% of members (and >= 2) deviating on one side
    'merge_gap_probes': 1, 'merge_gap_bp': 100_000,
    'common_af': 0.05, 'prior_min_an': 20,   # a region is common when a prior segment has DEL or DUP AF >= 0.05
    'z_out': 3.5,                            # member flag: above median + 3.5 x max(1.4826 MAD, scale floor)
    'noise_scale': 0.15,                     # noise scale floor: 0.15 x the members' median noise
    'noise_extreme': 2.5,                    # noise > 2.5x the members' median: another centre, kit or batch
    'burden_scale': 0.002, 'burden_floor': 0.005,  # burden scale floor; flagged burden also >= 0.5% of probes
    'pon_divided_one': 0.01,                 # > 1% of probes exactly 1.0: a PoN-divided .norm, not raw depth
    'missing_scale': 0.005, 'missing_floor': 0.03,  # missing > 3% of the probes CNVscope uses (members, tests,
                                             # admission); a .norm recovered by pon_norm_recover.py lacks the masked
                                             # (unused) probes only
    'dev_log2': 0.5, 'dev_k': 4.0,           # a member deviates where |5-probe running median of its LOO ratio|
                                             # > max(0.5, 4 x its robust sd): CNV burden and two-state counts
    'run_log2': 0.3, 'run_k': 3.0,           # CNV run: consecutive probes with the 5-probe running median beyond
                                             # max(0.3, 3 x its robust sd) on one side, outside common CNV loci
    'large_cnv_probes': 100,                 # member flag / admission review: a run of >= 100 probes (healthy
                                             # members' longest run: median 20-33, p99 72-111, max 90-174)
    'x_one_copy': -0.5, 'y_present': -2.0,   # sex: chrX / chrY log2 depth against the panel's same-sex reference
    'pc1_warn': 0.15, 'pc1_fail': 0.25,      # PC1 variance share above its permutation null (noise-scaled profiles)
    'sil_warn': 0.60, 'sil_fail': 0.70,      # 2-means silhouette on PC1-3
    'test_warn_ratio': 1.75,                 # test sample: warn above 1.75x the members' median noise
    'outlier_warn': 0.20,                    # share of members with a noise, burden or large_cnv flag
    'turnover_a': 0.027, 'turnover_b': 0.745,  # compare: share of calls that change = a + b x shift_rsd (36 study PoN
                                             # pairs, R^2 0.78, residual sd 0.008)
}


VERBOSE = True


def log(*a):
    if VERBOSE:
        print(time.strftime('%H:%M:%S'), *a, file=sys.stderr, flush=True)


# ---- I/O -----------------------------------------------------------------------------------------------------------
def contig_code(names):
    """Deterministic contig codes (chr1-22 -> 1-22, X 23, Y 24, others hashed): the same in every process."""
    cat = pd.Categorical(names)
    codes = []
    for c in cat.categories:
        b = c[3:] if c.startswith('chr') else c
        if b.isdigit() and 1 <= int(b) <= 22:
            codes.append(int(b))
        elif b in ('X', 'Y', 'M', 'MT'):
            codes.append({'X': 23, 'Y': 24, 'M': 25, 'MT': 25}[b])
        else:
            codes.append(100 + zlib.crc32(c.encode()) % 9000)
    return np.asarray(codes, dtype=np.int64)[cat.codes]


def pack(code, start, end):
    s = np.asarray(start, dtype=np.int64)
    return np.asarray(code, dtype=np.int64) * _KC + s * _KS + (np.asarray(end, dtype=np.int64) - s)


def read_list(path):
    with open(path) as f:
        return [ln.strip() for ln in f if ln.strip() and not ln.startswith('#')]


def read_sex(path):
    m = {}
    with open(path) as f:
        for ln in f:
            if ln.startswith('#') or not ln.strip():
                continue
            p = ln.strip().split('\t')
            if len(p) >= 2 and p[1] in ('M', 'F'):
                m[p[0]] = p[1]
    return m


def read_bed(path):
    reg = {}
    with open(path) as f:
        for ln in f:
            p = ln.strip().split('\t')
            if ln.startswith('#') or len(p) < 3:
                continue
            reg.setdefault(p[0], []).append((int(p[1]), int(p[2])))
    return reg


def in_regions(contig, start, end, reg):
    hit = np.zeros(len(start), dtype=bool)
    for c, rr in reg.items():
        m = np.flatnonzero(contig == c)
        for bs, be in rr:
            hit[m[(start[m] < be) & (end[m] > bs)]] = True
    return hit


def sample_name(path):
    b = os.path.basename(path)
    return b[:-5] if b.endswith('.norm') else b


def read_norm(path):
    """(keys, DP) of a raw .norm, parsed as create_pon.py parses it (float64)."""
    df = pd.read_csv(path, sep='\t', usecols=[0, 1, 2, 3], header=0)
    if len(df) == 0:
        return np.zeros(0, np.int64), np.zeros(0)
    k = pack(contig_code(df.iloc[:, 0].astype(str).values), df.iloc[:, 1].values, df.iloc[:, 2].values)
    return k, df.iloc[:, 3].values.astype(np.float64)


def read_pon_tsv(path):
    hdr, nskip = {}, 0
    with open(path) as f:
        for ln in f:
            if not ln.startswith('#'):
                break
            nskip += 1
            k, _, v = ln[1:].strip().partition('=')
            hdr[k] = v
    cols = ['contig', 'start', 'end', 'median_log2dp', 'n_valid', 'excluded', 'probe_idx', 'sigma2_dp',
            'median_dp_male', 'median_dp_female', 'par_flag', 'panel_median_dp']
    df = pd.read_csv(path, sep='\t', skiprows=nskip, usecols=cols, dtype={'contig': str})
    df['key'] = pack(contig_code(df['contig'].values), df['start'].values, df['end'].values)
    return {'path': path, 'hdr': hdr, 'df': df}


AR_MAGIC = b'!<arch>\n'


def ar_members(path):
    """{name: (data offset, size)} of the ar archive path, read as prior/cnvprior.lua reads a model bundle: 60-byte
    headers, the // table of long names, the / and /SYM64/ symbol tables skipped, 2-byte alignment, first name kept."""
    m, longs, off = {}, None, 8
    with open(path, 'rb') as f:
        total = f.seek(0, os.SEEK_END)
        while off < total:
            f.seek(off)
            h = f.read(60)
            if len(h) < 60 or h[58:60] != b'`\n':
                raise ValueError(f'{path}: bad ar member header at byte {off}')
            ident = h[:16].rstrip(b' ').decode('latin-1')
            sz = re.fullmatch(rb'([0-9]+) *', h[48:58])
            size = int(sz.group(1)) if sz else -1
            if size < 0 or off + 60 + size > total:
                raise ValueError(f'{path}: bad ar member size at byte {off}')
            name = None
            if ident == '//':
                longs = f.read(size)
            elif re.fullmatch('/[0-9]+', ident):   # a long name: an offset into the // table, ended by "/\n"
                i = int(ident[1:])
                j = longs.find(b'/\n', i) if longs is not None else -1
                if j < 0:
                    raise ValueError(f'{path}: bad ar long name {ident} at byte {off}')
                name = longs[i:j].decode('latin-1')
            elif re.fullmatch('[^/]+/', ident):
                name = ident[:-1]
            elif ident not in ('/', '/SYM64/'):
                raise ValueError(f"{path}: unknown ar member name '{ident}' at byte {off}")
            if name is not None and name not in m:
                m[name] = (off + 60, size)
            off += 60 + size + size % 2
    return m


def prior_member(path):
    """(bundle, member) of a prior given as BUNDLE/MEMBER, split at the last /, where BUNDLE is a model bundle (an ar
    archive, told by its first 8 bytes) holding MEMBER; (None, None) for a plain file. Reads archive headers only.
    Raises ValueError for a bare bundle (or BUNDLE/), naming its .vcf.gz members, for a missing member, and for a path
    that is neither, as prior/cnvprior.lua does (which also needs MEMBER.tbi)."""
    try:
        with open(path, 'rb') as f:
            if f.read(8) != AR_MAGIC:
                return None, None
        bundle, member = path, ''
    except OSError as e:
        bundle, _, member = path.rpartition('/')
        try:
            with open(bundle, 'rb') as f:
                ok = f.read(8) == AR_MAGIC
        except OSError:
            ok = False
        if not ok:
            raise ValueError(f'cannot open the prior {path}: {e.strerror}') from None
    m = ar_members(bundle)
    if not member:
        v = sorted(n for n in m if n.endswith('.vcf.gz'))
        raise ValueError(f'{bundle} is a model bundle: give the prior VCF in it as {bundle}/<member>; ' +
                         (f'its .vcf.gz members: {", ".join(v)}' if v else 'it holds no .vcf.gz member'))
    if member not in m:
        raise ValueError(f'model bundle {bundle} has no member {member}')
    return bundle, member


def read_prior(path, source=None):
    """Prior segments per contig: start0, end, AN, DEL_AF, DUP_AF (sorted by start). path is the prior .vcf.gz or
    BUNDLE/MEMBER (prior_member); a member is copied to a temp file (tempfile, so $TMPDIR), read and removed, also on
    error. A dict source receives the path, the member read (None for a plain file) and the md5 of the VCF bytes."""
    bundle, member = prior_member(path)
    h = hashlib.md5()
    if bundle is None:
        with open(path, 'rb') as f:
            for b in iter(lambda: f.read(1 << 20), b''):
                h.update(b)
        out = _read_prior_vcf(path)
    else:
        off, left = ar_members(bundle)[member]
        fd, tmp = tempfile.mkstemp(prefix='pon_qc.', suffix='.' + member)
        try:
            with os.fdopen(fd, 'wb') as o, open(bundle, 'rb') as f:
                f.seek(off)
                while left:
                    b = f.read(min(left, 1 << 20))
                    if not b:
                        raise ValueError(f'{bundle}: short read of member {member}')
                    h.update(b)
                    o.write(b)
                    left -= len(b)
            out = _read_prior_vcf(tmp)
        finally:
            os.remove(tmp)
    if source is not None:
        source.update(path=path, member=member, md5=h.hexdigest())
    return out


def _read_prior_vcf(path):
    """read_prior of a plain .vcf.gz, read with gzip from start to end (no index needed)."""
    out = {}
    with gzip.open(path, 'rt') as f:
        for ln in f:
            if ln.startswith('#'):
                continue
            p = ln.split('\t', 8)
            info = dict(kv.split('=', 1) for kv in p[7].strip().split(';') if '=' in kv)
            an = int(info.get('AN', 0))
            out.setdefault(p[0], []).append((int(p[1]) - 1, int(info['END']), an,
                                            int(info.get('DEL_AC', 0)) / an if an else 0.0,
                                            int(info.get('DUP_AC', 0)) / an if an else 0.0))
    return {c: np.array(sorted(v)) for c, v in out.items()}


_UNIV = None


def _load_one(path):
    k, dp = read_norm(path)
    row = np.full(len(_UNIV), np.nan)
    if len(k):
        i = np.minimum(np.searchsorted(_UNIV, k), len(_UNIV) - 1)
        ok = _UNIV[i] == k
        row[i[ok]] = dp[ok]
    return path, row, len(k)


class NormStore:
    """Raw .norm depth of many files, aligned to one sorted probe universe (float64, NaN = probe absent)."""

    def __init__(self, keys):
        self.keys = np.asarray(keys, dtype=np.int64)
        self.row, self.n_file = {}, {}

    def load(self, paths, jobs=1):
        global _UNIV
        todo = [p for p in dict.fromkeys(paths) if p not in self.row]
        if not todo:
            return
        _UNIV = self.keys
        if jobs > 1 and len(todo) > 1:
            with mp.get_context('fork').Pool(min(jobs, len(todo))) as pool:
                res = pool.map(_load_one, todo, chunksize=1)
        else:
            res = [_load_one(p) for p in todo]
        for p, r, n in res:
            self.row[p], self.n_file[p] = r, n

    def cols(self, keys):
        i = np.searchsorted(self.keys, keys)
        if not np.array_equal(self.keys[np.minimum(i, len(self.keys) - 1)], keys):
            raise ValueError('probe keys outside the store universe')
        return i

    def matrix(self, paths, cols):
        return np.stack([self.row[p][cols] for p in paths]) if paths else np.zeros((0, len(cols)))


# ---- estimators ----------------------------------------------------------------------------------------------------
def fast_col_median(mat, nv):
    """create_pon.py's _fast_col_median: column medians of a samples x probes matrix, NaN = missing."""
    n_rows, n_cols = mat.shape
    mt = np.array(mat.T, order='C', copy=True)       # always a copy: the sort is in place
    mt.sort(axis=1)
    if np.all(nv == n_rows):
        k = n_rows // 2
        return (mt[:, k - 1] + mt[:, k]) / 2.0 if n_rows % 2 == 0 else mt[:, k].copy()
    k_lo = np.maximum((nv - 1) // 2, 0).astype(np.intp)
    k_hi = np.maximum(nv // 2, 0).astype(np.intp)
    idx = np.arange(n_cols)
    res = (mt[idx, k_lo] + mt[idx, k_hi]) / 2.0
    res[nv == 0] = np.nan
    return res


def rmad(x):
    """1.4826 MAD about the median, NaN ignored."""
    x = x[~np.isnan(x)]
    if not len(x):
        return np.nan
    return MAD_K * float(np.median(np.abs(x - np.median(x))))


def mad_sigma2(L, nv):
    """Builder's per-probe sigma2 = (1.4826 MAD)^2 (float64 before its float32 cast)."""
    med = fast_col_median(L, nv)
    return (fast_col_median(np.abs(L - med), nv) * MAD_K) ** 2, med


def med_sorted(S, nv):
    P, N = S.shape
    r = np.arange(P)
    lo = np.clip((nv - 1) // 2, 0, N - 1)
    hi = np.clip(nv // 2, 0, N - 1)
    m = (S[r, lo].astype(np.float64) + S[r, hi]) * 0.5
    m[nv == 0] = np.nan
    return m


def sort_panel(LT):
    """LT: probes x members float32. Sorted values, order and rank (NaN last), valid counts."""
    P, N = LT.shape
    order = np.argsort(LT, axis=1, kind='stable')
    S = np.take_along_axis(LT, order, axis=1)
    dt = np.int16 if N < 32000 else np.int32
    order = order.astype(dt)
    rank = np.empty_like(order)
    np.put_along_axis(rank, order.astype(np.intp), np.broadcast_to(np.arange(N, dtype=dt), (P, N)), axis=1)
    nv = np.sum(~np.isnan(LT), axis=1)
    return S, order, rank, nv


def boot_se(S, order, nv, B=200, seed=20261110, chunk=8192):
    """Per-probe bootstrap SE of the median (members resampled with replacement; exact weighted medians)."""
    P, N = S.shape
    rng = np.random.default_rng(seed)
    cnt = np.zeros((B, N + 1), dtype=np.int16)
    for b in range(B):
        cnt[b, :N] = np.bincount(rng.integers(0, N, size=N), minlength=N)
    m0 = med_sorted(S, nv)
    s1, s2, nb = np.zeros(P), np.zeros(P), np.zeros(P)
    col = np.arange(N)
    for a in range(0, P, chunk):
        z = min(P, a + chunk)
        o = order[a:z].astype(np.intp)
        o[col[None, :] >= nv[a:z, None]] = N          # NaN slots weigh 0
        Sc, r, mc = S[a:z], np.arange(z - a), m0[a:z]
        for b in range(B):
            C = np.cumsum(cnt[b][o], axis=1, dtype=np.int16)
            n = C[:, -1].astype(np.intp)
            ok = n > 0
            ilo = np.minimum((C <= ((n - 1) // 2)[:, None]).sum(axis=1), N - 1)
            ihi = np.minimum((C <= (n // 2)[:, None]).sum(axis=1), N - 1)
            d = (Sc[r, ilo].astype(np.float64) + Sc[r, ihi]) * 0.5 - mc
            d[~ok] = 0.0
            s1[a:z] += d
            s2[a:z] += d * d
            nb[a:z] += ok
    with np.errstate(invalid='ignore', divide='ignore'):
        v = (s2 - s1 * s1 / nb) / (nb - 1)
    return np.sqrt(np.maximum(v, 0.0))


def loo_median(S, rank, nv, i):
    """Median of every probe without member i (the full median where i has no value)."""
    P, N = S.shape
    r = rank[:, i].astype(np.intp)
    has = r < nv
    n1 = nv - has
    lo, hi = (n1 - 1) // 2, n1 // 2
    lo = np.clip(lo + (has & (lo >= r)), 0, N - 1)
    hi = np.clip(hi + (has & (hi >= r)), 0, N - 1)
    rr = np.arange(P)
    m = (S[rr, lo].astype(np.float64) + S[rr, hi]) * 0.5
    m[n1 <= 0] = np.nan
    return m


def runmed5(x, code):
    """5-probe running median within each contig (edges take the nearest full window)."""
    out = np.empty_like(x)
    b = np.flatnonzero(np.diff(code)) + 1
    for a, z in zip(np.r_[0, b], np.r_[b, len(x)]):
        v = x[a:z]
        if len(v) < 5:
            out[a:z] = np.median(v)
            continue
        w = np.median(np.lib.stride_tricks.sliding_window_view(v, 5), axis=1)
        out[a + 2:z - 2] = w
        out[a:a + 2], out[z - 2:z] = w[0], w[-1]
    return out


def probe_prior_af(contig, start, end, prior, min_an):
    """Max DEL/DUP AF over the prior segments a probe overlaps (segments with AN >= min_an)."""
    af = np.zeros(len(start))
    if prior is None:
        return af
    for c in pd.unique(contig):
        seg = prior.get(c)
        if seg is None:
            continue
        m = np.flatnonzero(contig == c)
        saf = np.where(seg[:, 2] >= min_an, np.maximum(seg[:, 3], seg[:, 4]), 0.0)
        a = np.searchsorted(seg[:, 1], start[m], 'right')
        z = np.searchsorted(seg[:, 0], end[m], 'left')
        ok = z > a
        v = np.zeros(len(m))
        v[ok] = np.maximum(saf[np.minimum(a[ok], len(seg) - 1)], saf[np.maximum(z[ok] - 1, 0)])
        for j in np.flatnonzero(z - a > 2):
            v[j] = saf[a[j]:z[j]].max()
        af[m] = v
    return af


def run_list(sgn, ch, min_len=1):
    """Maximal runs of one non-zero sign within a contig: (first, last, sign, length), length >= min_len."""
    if not len(sgn):
        return []
    brk = np.r_[True, (sgn[1:] != sgn[:-1]) | (ch[1:] != ch[:-1])]
    st = np.flatnonzero(brk)
    ln = np.diff(np.r_[st, len(sgn)])
    keep = (sgn[st] != 0) & (ln >= min_len)
    return [(int(a), int(a + n - 1), int(sgn[a]), int(n)) for a, n in zip(st[keep], ln[keep])]


def ratio_stats(R, has, used, code, t):
    """Noise, probe-to-probe noise, burden and the 5-probe running median of one ratio vector (used probes)."""
    hu = has & used
    r = {'noise_all': rmad(R[has]) if has.any() else np.nan}
    if hu.sum() < 100:
        return r, None, hu
    Rh, ch = R[hu], code[hu]
    d = np.diff(Rh)[ch[1:] == ch[:-1]]
    sm = runmed5(Rh, ch)
    ssd = rmad(sm)
    dv = max(t['dev_log2'], t['dev_k'] * ssd)
    r.update({'noise': rmad(Rh), 'noise_dlr': rmad(d) / np.sqrt(2), 'dev_thr': dv,
              'burden': float(np.mean(np.abs(sm) > dv)), 'burden_del': float(np.mean(sm < -dv)),
              'burden_dup': float(np.mean(sm > dv)), 'sm_sd': ssd, 'auto_off': float(np.median(Rh))})
    return r, sm, hu


def cnv_runs(sm, hu, ssd, code, common, pinfo, t, min_len=1):
    """Runs of the running median beyond the run threshold, outside common CNV loci, as dicts."""
    rv = max(t['run_log2'], t['run_k'] * ssd)
    sgn = np.where(sm > rv, 1, np.where(sm < -rv, -1, 0)).astype(np.int8)
    sgn[common[hu]] = 0
    hi = np.flatnonzero(hu)
    out = []
    for a, z, sg, n in run_list(sgn, code[hu], min_len):
        out.append({'contig': pinfo['contig'][hi[a]], 'start': int(pinfo['start'][hi[a]]),
                    'end': int(pinfo['end'][hi[z]]), 'n_probes': n, 'dir': 'dup' if sg > 0 else 'del',
                    'mean_sm': float(np.mean(sm[a:z + 1]))})
    return out, rv


def longest(runs):
    if not runs:
        return {'max_run': 0, 'max_run_dir': '', 'max_run_region': '', 'max_run_mean': np.nan}
    x = max(runs, key=lambda u: u['n_probes'])
    return {'max_run': x['n_probes'], 'max_run_dir': x['dir'],
            'max_run_region': f"{x['contig']}:{x['start']}-{x['end']}", 'max_run_mean': x['mean_sm']}


def kmeans2(X, rng, n_init=10, iters=100):
    best, lab_best = np.inf, np.zeros(len(X), dtype=int)
    for _ in range(n_init):
        c0 = X[rng.integers(len(X))]
        d = ((X - c0) ** 2).sum(1)
        if d.sum() <= 0:
            continue
        C = np.stack([c0, X[rng.choice(len(X), p=d / d.sum())]])
        for _ in range(iters):
            D = ((X[:, None, :] - C[None]) ** 2).sum(2)
            lab = D.argmin(1)
            if lab.min() == lab.max():
                break
            Cn = np.stack([X[lab == k].mean(0) for k in (0, 1)])
            if np.allclose(Cn, C):
                break
            C = Cn
        inertia = D[np.arange(len(X)), lab].sum()
        if lab.min() != lab.max() and inertia < best:
            best, lab_best = inertia, lab.copy()
    return lab_best


def silhouette(X, lab):
    g = X @ X.T
    sq = np.diag(g)
    D = np.sqrt(np.maximum(sq[:, None] + sq[None, :] - 2 * g, 0))
    s = np.zeros(len(X))
    for i in range(len(X)):
        same = lab == lab[i]
        if same.sum() <= 1 or (~same).sum() == 0:
            continue
        a = D[i, same].sum() / (same.sum() - 1)
        b = D[i, ~same].mean()
        s[i] = (b - a) / max(a, b) if max(a, b) > 0 else 0.0
    return float(s.mean())


def nanmed(x):
    x = np.asarray(x)
    x = x[~np.isnan(x)]
    return float(np.median(x)) if len(x) else np.nan


def robust_z(v):
    m = np.nanmedian(v)
    s = MAD_K * np.nanmedian(np.abs(v - m))
    return (v - m) / s if s > 0 else np.zeros_like(v), m


def status(v, warn, fail, high_bad=True):
    if v is None or (isinstance(v, float) and np.isnan(v)):
        return 'na'
    if high_bad:
        return 'fail' if v > fail else 'warn' if v > warn else 'pass'
    return 'fail' if v < fail else 'warn' if v < warn else 'pass'


# ---- the QC --------------------------------------------------------------------------------------------------------
def derive_table(base, store, paths, min_samples, min_median_dp=0.02):
    """A builder-rule autosomal table for a panel without its own PoN (probe index and chrX/chrY rows from base)."""
    df = base['df'].copy()
    auto = df['contig'].isin(AUTO).values
    cols = store.cols(df['key'].values[auto])
    M = store.matrix(paths, cols)
    nv = np.sum(~np.isnan(M), axis=0)
    mdp = fast_col_median(M, nv)
    exc = (nv < min_samples) | (mdp < min_median_dp) | np.isnan(mdp)
    L = np.log2(M + PSEUDO)
    with np.errstate(all='ignore'):
        s2, _ = mad_sigma2(L, nv)
        med = np.nanmedian(L, axis=0)
    df.loc[auto, 'n_valid'] = nv
    df.loc[auto, 'panel_median_dp'] = mdp
    df.loc[auto, 'excluded'] = exc.astype(int)
    df.loc[auto, 'probe_idx'] = np.where(exc, -1, 0)
    df.loc[auto, 'median_log2dp'] = np.where(exc, 0.0, np.float32(med))
    df.loc[auto, 'sigma2_dp'] = np.where(exc, 0.0, np.float32(s2))
    hdr = {'n_samples': str(len(paths)), 'virtual': '1'}
    return {'path': '(derived)', 'hdr': hdr, 'df': df}


def run_qc(paths, sexes, tables, store, par, min_samples=None, prior=None, B=200, seed=20261110, n_pca=20000,
           tests=(), admits=(), admit_sex=None, verify=True, opt=None):
    """Checks 0-7 for one panel. Returns a dict: res, probes, members, regions, tests, admit, runs."""
    t = dict(T, **(opt or {}))
    t0 = time.time()
    N = len(paths)
    names = [sample_name(p) for p in paths]
    sexes = np.asarray(sexes)
    nM, nF = int((sexes == 'M').sum()), int((sexes == 'F').sum())
    df = tables[0]['df']
    contig = df['contig'].values
    va = (df['excluded'].values == 0) & np.isin(contig, list(AUTO))
    A = df[va]
    keys = A['key'].values
    code = contig_code(A['contig'].values)
    P = len(A)
    cf = float(tables[0]['hdr'].get('capture_floor', 0.1))
    used = A['panel_median_dp'].values >= cf          # the valid probes CNVscope does not mask
    res = {'n_members': N, 'n_valid_auto': P, 'n_used_auto': int(used.sum()), 'capture_floor': cf,
           'table': tables[0]['path']}

    # matrix of valid autosomal probes (float64, as the builder reads it)
    DP = store.matrix(paths, store.cols(keys))
    n_file = np.array([store.n_file[p] for p in paths])
    allcols = store.cols(df['key'].values)
    n_in = np.array([np.count_nonzero(~np.isnan(store.row[p][allcols])) for p in paths])
    empty = n_file == 0
    with np.errstate(invalid='ignore'):
        one_share = np.sum(DP == 1.0, axis=1) / np.maximum(np.sum(~np.isnan(DP), axis=1), 1)
    nvd = np.sum(~np.isnan(DP), axis=0)
    L = np.log2(DP + PSEUDO)
    del DP

    # 0 verification against the stored median and sigma2
    if verify:
        with np.errstate(all='ignore'):
            med64 = np.nanmedian(L, axis=0)
        s2, _ = mad_sigma2(L, nvd)
        st_med, st_s2 = A['median_log2dp'].values, A['sigma2_dp'].values
        rnd = lambda x: np.array([float('%.6f' % v) for v in np.float32(x)])
        res['verify'] = {
            'median_max_abs_diff': float(np.nanmax(np.abs(med64 - st_med))),
            'median_mismatch_after_rounding': int(np.sum(np.abs(rnd(med64) - st_med) > 1e-9)),
            'sigma2_max_abs_diff': float(np.nanmax(np.abs(s2 - st_s2))),
            'sigma2_mismatch_after_rounding': int(np.sum(np.abs(rnd(s2) - st_s2) > 1e-9)),
            'n_valid_mismatch': int(np.sum(nvd != A['n_valid'].values)),
        }
        del s2
    st_med = A['median_log2dp'].values.astype(np.float64)
    LT = np.ascontiguousarray(L.T, dtype=np.float32)
    del L
    S, order, rank, nv = sort_panel(LT)
    med = med_sorted(S, nv)
    sig = np.sqrt(np.maximum(mad_sigma2(LT.T, nv)[0], 0)) if not verify else np.sqrt(A['sigma2_dp'].values)
    log(f'  matrix {N} x {P}, sorted ({time.time() - t0:.0f}s)')

    # 1 composition; the rule is the builder's --min_samples auto (a sex without members with probes does not count)
    nMp, nFp = int(((sexes == 'M') & ~empty).sum()), int(((sexes == 'F') & ~empty).sum())
    rule = min(20, max(3, min([n for n in (nMp, nFp) if n] or [0]) // 2))
    lo_ms = int(A['n_valid'].min()) if P else 0
    ex_auto = (df['excluded'].values == 1) & np.isin(contig, list(AUTO)) & (df['panel_median_dp'].values >= 0.05)
    hi_ms = int(df['n_valid'].values[ex_auto].max()) + 1 if ex_auto.any() else 1
    comp = {'n': N, 'n_M': nM, 'n_F': nF, 'empty': [names[i] for i in np.flatnonzero(empty)],
            'n_M_with_probes': nMp, 'n_F_with_probes': nFp,
            'hdr_n_samples': tables[0]['hdr'].get('n_samples'), 'min_samples': min_samples,
            'min_samples_rule': rule, 'min_samples_inferred': [hi_ms, lo_ms], 'tables': []}
    for tb in tables:
        d = tb['df']
        rs = tb['hdr'].get('ratio_sex', '')
        row = {'path': tb['path'], 'ratio_sex': rs, 'n_samples': tb['hdr'].get('n_samples'),
               'n_males': tb['hdr'].get('n_males'), 'n_females': tb['hdr'].get('n_females'),
               'valid': int((d['excluded'] == 0).sum())}
        for lab, nm in (('chrX', XN), ('chrY', YN)):
            m = d['contig'].isin(nm).values
            row[f'{lab}_ref'] = int(np.sum(m & (d['probe_idx'].values >= 0)))
            row[f'{lab}_dark'] = int(np.sum(m & (d['probe_idx'].values < 0)))
        comp['tables'].append(row)
    res['composition'] = comp

    # 2 precision
    se = boot_se(S, order, nv, B=B, seed=seed)
    se_an = MED_K * sig / np.sqrt(np.maximum(nv, 1))
    with np.errstate(invalid='ignore', divide='ignore'):
        rel = np.where(sig > 0, se / sig, np.nan)
    prec = {'B': B, 'n_probes': P, 'n_used': int(used.sum())}
    for sfx, m in (('', np.ones(P, bool)), ('_used', used)):
        for lab, x in (('boot', se[m]), ('analytic', se_an[m])):
            prec[f'se_{lab}_median{sfx}'] = float(np.nanmedian(x))
            prec[f'se_{lab}_p90{sfx}'] = float(np.nanpercentile(x, 90))
            prec[f'se_{lab}_p99{sfx}'] = float(np.nanpercentile(x, 99))
            for c in (0.05, 0.1, 0.2):
                prec[f'frac_se_{lab}_gt_{c}{sfx}'] = float(np.nanmean(x > c))
        prec[f'rel_se_median{sfx}'] = float(np.nanmedian(rel[m]))
    prec['rel_se_clean'] = MED_K / np.sqrt(np.median(nv))
    prec['sigma2_median'] = float(np.median(sig ** 2))
    res['precision'] = prec
    log(f'  bootstrap B={B} ({time.time() - t0:.0f}s)')

    # 4 members (LOO)
    xm, ym = sex_probe_sets(df, par)
    DX = store.matrix(paths, store.cols(df['key'].values[xm])) if xm.any() else np.zeros((N, 0))
    DY = store.matrix(paths, store.cols(df['key'].values[ym])) if ym.any() else np.zeros((N, 0))
    LX, LY = np.log2(DX + PSEUDO), np.log2(DY + PSEUDO)
    refX = sex_ref(LX, sexes, 'F', 'M', +1.0)
    refY = sex_ref(LY, sexes, 'M', None, 0.0)
    rng = np.random.default_rng(seed + 1)
    iu = np.flatnonzero(used)
    sub = np.sort(rng.choice(iu, size=min(n_pca, len(iu)), replace=False))
    prof = np.zeros((N, len(sub)), dtype=np.float32)
    cnt_lo, cnt_hi = np.zeros(P, np.int32), np.zeros(P, np.int32)
    SM = np.full((N, P), np.nan, dtype=np.float32)    # members' 5-probe running medians, for the CNV runs
    pinfo = {'contig': A['contig'].values, 'start': A['start'].values, 'end': A['end'].values}
    rows = []
    for i in range(N):
        R = LT[:, i] - loo_median(S, rank, nv, i)
        has = ~np.isnan(R)
        prof[i] = np.clip(np.nan_to_num(R[sub], nan=0.0), -2, 2)
        r = {'sample': names[i], 'path': paths[i], 'sex': sexes[i], 'n_probes_file': int(n_file[i]),
             'extra_frac': float((n_file[i] - n_in[i]) / n_file[i]) if n_file[i] else np.nan,
             'missing_frac': float(1 - has[used].mean()), 'missing_frac_all': float(1 - has.mean()),
             'one_share': float(one_share[i])}
        st, sm, hu = ratio_stats(R, has, used, code, t)
        r.update(st)
        if sm is not None:
            hi_ = np.flatnonzero(hu)
            cnt_lo[hi_[sm < -st['dev_thr']]] += 1
            cnt_hi[hi_[sm > st['dev_thr']]] += 1
            SM[i, hi_] = sm
            r['x_l2'] = nanmed(LX[i] - refX) - st['auto_off']
            r['y_l2'] = nanmed(LY[i] - refY) - st['auto_off']
        rows.append(r)
    # common CNV loci: two-state in this panel, or DEL/DUP AF >= common_af in the prior
    f_lo, f_hi = cnt_lo / np.maximum(nv, 1), cnt_hi / np.maximum(nv, 1)
    fl_bi = ((f_lo >= t['bimodal_frac']) & (cnt_lo >= 2)) | ((f_hi >= t['bimodal_frac']) & (cnt_hi >= 2))
    paf = probe_prior_af(pinfo['contig'], pinfo['start'], pinfo['end'], prior, t['prior_min_an'])
    common = fl_bi | (paf >= t['common_af'])
    for i, r in enumerate(rows):
        if 'sm_sd' in r:
            hu = ~np.isnan(SM[i])
            runs, rv = cnv_runs(SM[i, hu], hu, r['sm_sd'], code, common, pinfo, t)
            r.update(longest(runs))
            r['run_thr'] = rv
    del SM
    mem = pd.DataFrame(rows).drop(columns=['auto_off', 'sm_sd'], errors='ignore')
    for c in ('noise', 'noise_all', 'noise_dlr', 'dev_thr', 'burden', 'burden_del', 'burden_dup', 'run_thr', 'max_run',
              'max_run_dir', 'max_run_region', 'max_run_mean', 'x_l2', 'y_l2'):
        if c not in mem:
            mem[c] = np.nan
    mem = flag_members(mem, empty, t)
    res['common_cnv_probes'] = {'two_state': int((fl_bi & used).sum()), 'prior': int(((paf >= t['common_af']) & used).sum()),
                                'either': int((common & used).sum()), 'used': int(used.sum())}
    res['members'] = {'n_flagged': int((mem['flags'] != '').sum()),
                      'frac_flagged': float((mem['flags'] != '').mean()),
                      'noise_median': float(np.nanmedian(mem['noise'])), 'noise_max': float(np.nanmax(mem['noise'])),
                      'burden_median': float(np.nanmedian(mem['burden'])),
                      'by_rule': {k: int(mem['flags'].str.split(',').apply(lambda f: k in f).sum()) for k in
                                  ('empty', 'pon_divided', 'noise_extreme', 'noise', 'burden', 'large_cnv',
                                   'missing', 'sex_mismatch', 'sex_ambiguous')},
                      'flagged': {s: f for s, f in zip(mem['sample'], mem['flags']) if f},
                      'sex_flagged': {sx: int(((mem['sex'] == sx) & mem['flags'].str.contains('sex_')).sum())
                                      for sx in ('M', 'F')}}
    log(f'  members ({time.time() - t0:.0f}s)')

    # 3 unstable loci (two-state counts come from the members' smoothed LOO ratios)
    dev = t['bimodal_dev']
    raw_lo = np.sum(S < (med - dev)[:, None].astype(np.float32), axis=1)
    raw_hi = np.sum(S > (med + dev)[:, None].astype(np.float32), axis=1)
    raw_bi = (np.maximum(raw_lo, raw_hi) >= t['bimodal_frac'] * nv) & (np.maximum(raw_lo, raw_hi) >= 2) & used
    fl_se = (se > t['se_unstable']) & used
    flag = fl_se | fl_bi
    probes = {'key': keys, 'code': code, 'start': A['start'].values, 'end': A['end'].values, 'contig': A['contig'].values,
              'med': med, 'stored_med': st_med, 'sigma2': sig ** 2, 'se_boot': se, 'se_an': se_an, 'nv': nv,
              'f_lo': f_lo, 'f_hi': f_hi, 'flag': flag, 'used': used, 'common': common}
    regions = merge_regions(probes, flag, fl_se, fl_bi, t, prior)
    res['unstable'] = {'probes_flagged': int(flag.sum()), 'probes_se': int(fl_se.sum()), 'probes_bimodal': int(fl_bi.sum()),
                       'probes_offmedian_unsmoothed': int(raw_bi.sum()),
                       'frac_probes_flagged': float(flag.mean()), 'regions': len(regions),
                       'regions_multi': int((regions['n_flag'] >= 2).sum()) if len(regions) else 0,
                       'regions_common': int(regions['common'].sum()) if len(regions) else 0,
                       'regions_multi_common': int(((regions['n_flag'] >= 2) & regions['common']).sum()) if len(regions) else 0,
                       'bp': int((regions['end'] - regions['start']).sum()) if len(regions) else 0}

    # 5 batch structure
    res['batch'] = batch_structure(prof[~empty], rng)
    log(f'  batch ({time.time() - t0:.0f}s)')

    # 6 test samples and 7 admission candidates, scored against the stored medians
    ctx = (store, df, keys, st_med, used, code, common, pinfo, mem, empty, refX, refY, xm, ym, t)
    tr = []
    for tp in tests:
        r, _ = score_sample(tp, *ctx)
        r['status'] = test_status(r, t)
        tr.append(r)
    ar, runs = [], []
    asx = admit_sex or {}
    for ap_ in admits:
        r, rr = score_sample(ap_, *ctx, min_len=10)
        r['sex_declared'] = asx.get(ap_, '')
        r['verdict'], r['rules'] = admit_verdict(r, r['sex_declared'], t)
        ar.append(r)
        runs += rr
    res['tests'], res['admit'] = tr, ar
    res['seconds'] = round(time.time() - t0, 1)
    res['status'] = panel_status(res, t)
    return {'res': res, 'probes': probes, 'members': mem, 'regions': regions, 'tests': pd.DataFrame(tr),
            'admit': pd.DataFrame(ar), 'runs': pd.DataFrame(runs)}


def sex_probe_sets(df, par):
    contig, st, en = df['contig'].values, df['start'].values, df['end'].values
    inpar = in_regions(contig, st, en, par) if par else np.zeros(len(df), bool)
    xm = np.isin(contig, XN) & ~inpar & (df['panel_median_dp'].values >= 0.1)
    with np.errstate(invalid='ignore'):
        ym = np.isin(contig, YN) & ~inpar & (df['median_dp_male'].values >= np.log2(0.1))
    return xm, ym


def sex_ref(LZ, sexes, ref_sex, alt_sex, alt_shift):
    """Per-probe median log2 depth of the reference sex (two-copy level for chrX, male level for chrY)."""
    if LZ.shape[1] == 0:
        return np.zeros(0)
    for s, shift in ((ref_sex, 0.0), (alt_sex, alt_shift)):
        if s is None:
            continue
        m = sexes == s
        if m.sum() >= 1:
            sub = LZ[m]
            return fast_col_median(sub, np.sum(~np.isnan(sub), axis=0)) + shift
    return np.full(LZ.shape[1], np.nan)


def infer_sex(x, y, t):
    """Sex from chrX dosage and chrY presence; chrX alone when the depth file has no chrY probe."""
    if np.isnan(x):
        return 'NA'
    one_x = x < t['x_one_copy']
    has_y = (y > t['y_present']) if not np.isnan(y) else one_x
    if one_x and has_y:
        return 'M'
    if not one_x and not has_y:
        return 'F'
    return 'XO-like' if one_x else 'XXY-like'


def flag_members(mem, empty, t):
    ok = ~empty
    flags = [[] for _ in range(len(mem))]
    for i in np.flatnonzero(empty):
        flags[i].append('empty')
    for i in np.flatnonzero(mem['one_share'].fillna(0).values > t['pon_divided_one']):
        flags[i].append('pon_divided')
    for c, lab, scale, floor in (('noise', 'noise', None, None), ('burden', 'burden', t['burden_scale'], t['burden_floor']),
                                 ('missing_frac', 'missing', t['missing_scale'], t['missing_floor'])):
        v = mem[c].values.astype(float)
        z = np.full(len(v), np.nan)
        if ok.sum() >= 3:
            m = np.nanmedian(v[ok])
            sc = max(MAD_K * np.nanmedian(np.abs(v[ok] - m)), t['noise_scale'] * m if scale is None else scale)
            z[ok] = (v[ok] - m) / sc
            hit = ok & (z > t['z_out'])
            if floor is not None:
                hit &= v > floor
            for i in np.flatnonzero(hit):
                flags[i].append(lab)
        mem[f'z_{c}'] = z
        if c == 'burden':
            for i in np.flatnonzero(ok & (mem['max_run'].fillna(0).values >= t['large_cnv_probes'])):
                flags[i].append('large_cnv')
        if c == 'noise' and ok.sum() >= 3:
            for i in np.flatnonzero(ok & (v > t['noise_extreme'] * np.nanmedian(v[ok]))):
                flags[i] = ['noise_extreme' if f == 'noise' else f for f in flags[i]]
                if 'noise_extreme' not in flags[i]:
                    flags[i].append('noise_extreme')
    inf = [infer_sex(x, y, t) for x, y in zip(mem['x_l2'].values, mem['y_l2'].values)]
    mem['sex_inferred'] = inf
    mem['sex_basis'] = ['NA' if np.isnan(x) else 'X only' if np.isnan(y) else 'X+Y'
                        for x, y in zip(mem['x_l2'].values, mem['y_l2'].values)]
    for i, (d, s) in enumerate(zip(mem['sex'].values, inf)):
        if empty[i] or s == 'NA':
            continue
        if s not in ('M', 'F'):
            flags[i].append('sex_ambiguous')
        elif s != d:
            flags[i].append('sex_mismatch')
    mem['flags'] = [','.join(f) for f in flags]
    return mem


def batch_structure(prof, rng, n_null=3):
    """PCA of the LOO profiles, each member scaled by its own robust sd; null = probes permuted across members."""
    X = prof.astype(np.float64)
    sd = np.array([rmad(x) for x in X])
    X /= np.where(sd > 0, sd, 1.0)[:, None]
    n = len(X)
    out = {'n': n}
    if n < 6:
        return out

    def pca(Y):
        Y = Y - Y.mean(0)
        U, s, _ = np.linalg.svd(Y, full_matrices=False)
        lam = s ** 2
        return lam / lam.sum(), U[:, :3] * s[:3]

    frac, sc = pca(X)
    lab = kmeans2(sc, np.random.default_rng(1))
    pc1_null = float(np.mean([pca(rng.permuted(X, axis=0))[0][0] for _ in range(n_null)]))
    out.update({'pc1': float(frac[0]), 'pc2': float(frac[1]), 'pc3': float(frac[2]), 'pc1_null': pc1_null,
                'pc1_struct': float(frac[0] - pc1_null), 'sil': silhouette(sc, lab),
                'clusters': sorted([int((lab == 0).sum()), int((lab == 1).sum())]), 'labels': lab.tolist()})
    return out


def score_sample(path, store, df, keys, st_med, used, code, common, pinfo, mem, empty, refX, refY, xm, ym, t,
                 min_len=10):
    """A non-member's ratio to the stored medians, scored with the member statistics and against the members."""
    store.load([path])
    row = store.row[path]
    x = row[store.cols(keys)]
    n_file = store.n_file[path]
    n_in = int(np.count_nonzero(~np.isnan(row[store.cols(df['key'].values)])))
    R = np.log2(x + PSEUDO) - st_med
    has = ~np.isnan(R)
    st, sm, hu = ratio_stats(R, has, used, code, t)
    one = float(np.mean(x[has] == 1.0)) if has.any() else np.nan
    hu_ = has & used
    pcor = float(np.corrcoef(np.log2(x[hu_] + PSEUDO), st_med[hu_])[0, 1]) if hu_.sum() > 100 else np.nan
    r = {'sample': sample_name(path), 'path': path, 'n_probes_file': n_file, 'missing_frac': float(1 - has[used].mean()),
         'missing_frac_all': float(1 - has.mean()),
         'extra_frac': float((n_file - n_in) / n_file) if n_file else np.nan, 'one_share': one, 'pon_corr': pcor}
    r.update({k: v for k, v in st.items() if k not in ('auto_off', 'sm_sd')})
    ok = ~np.asarray(empty)
    ref, bref = mem['noise'].values[ok].astype(float), mem['burden'].values[ok].astype(float)
    ref, bref = ref[~np.isnan(ref)], bref[~np.isnan(bref)]
    noise = r.get('noise', np.nan)
    m = float(np.nanmedian(ref)) if len(ref) else np.nan
    r.update({'noise_pct': float(100 * np.nanmean(ref <= noise)) if len(ref) and not np.isnan(noise) else np.nan,
              'noise_over_median': float(noise / m) if len(ref) and not np.isnan(noise) else np.nan,
              'noise_over_max': float(noise / np.nanmax(ref)) if len(ref) and not np.isnan(noise) else np.nan,
              'members_noise_median': m})
    if len(ref) >= 3 and not np.isnan(noise):
        r['z_noise'] = (noise - m) / max(MAD_K * np.nanmedian(np.abs(ref - m)), t['noise_scale'] * m)
        mb = float(np.nanmedian(bref))
        r['z_burden'] = (r['burden'] - mb) / max(MAD_K * np.nanmedian(np.abs(bref - mb)), t['burden_scale'])
    runs = []
    if sm is not None:
        rl, rv = cnv_runs(sm, hu, st['sm_sd'], code, common, pinfo, t)
        r.update(longest(rl))
        r['run_thr'] = rv
        runs = [dict(sample=r['sample'], **u) for u in rl if u['n_probes'] >= min_len]
    off = st.get('auto_off', np.nan)
    lx = np.log2(row[store.cols(df['key'].values[xm])] + PSEUDO) if xm.any() else np.zeros(0)
    ly = np.log2(row[store.cols(df['key'].values[ym])] + PSEUDO) if ym.any() else np.zeros(0)
    xl = nanmed(lx - refX) - off if len(lx) and not np.isnan(off) else np.nan
    yl = nanmed(ly - refY) - off if len(ly) and not np.isnan(off) else np.nan
    r.update({'x_l2': xl, 'y_l2': yl, 'sex_inferred': infer_sex(xl, yl, t),
              'sex_basis': 'NA' if np.isnan(xl) else 'X only' if np.isnan(yl) else 'X+Y'})
    return r, runs


def test_status(r, t):
    if np.isnan(r.get('noise', np.nan)) or r['missing_frac'] > t['missing_floor'] or \
            r['noise_over_median'] > t['noise_extreme'] or r.get('one_share', 0) > t['pon_divided_one']:
        return 'fail'
    return 'warn' if r['noise_over_median'] > t['test_warn_ratio'] else 'pass'


def admit_verdict(r, declared, t):
    """reject: the sample does not match the panel; review: a member flag would fire; admit otherwise."""
    rej, rev = [], []
    if not r['n_probes_file'] or np.isnan(r.get('noise', np.nan)):
        return 'reject', 'empty'
    if r.get('one_share', 0) > t['pon_divided_one']:
        return 'reject', 'pon_divided'
    if r['missing_frac'] > t['missing_floor']:
        rej.append('missing')
    if r['noise_over_median'] > t['noise_extreme']:
        rej.append('noise_extreme')
    elif r.get('z_noise', np.nan) > t['z_out']:
        rev.append('noise')
    if r.get('z_burden', np.nan) > t['z_out'] and r['burden'] > t['burden_floor']:
        rev.append('burden')
    if r.get('max_run', 0) >= t['large_cnv_probes']:
        rev.append('large_cnv')
    s_ = r['sex_inferred']
    if s_ in ('XO-like', 'XXY-like'):
        rej.append('sex_ambiguous')
    elif declared in ('M', 'F') and s_ in ('M', 'F') and s_ != declared:
        rev.append('sex_mismatch')
    return ('reject' if rej else 'review' if rev else 'admit'), ','.join(rej + rev)


def compare_tables(new, old, t, xy=()):
    """Stored medians and sigma2 of PoN `new` against `old` over the autosomal probes valid in both."""
    def valid(tb):
        d = tb['df']
        return d[(d['excluded'].values == 0) & d['contig'].isin(AUTO).values]

    a, b = valid(new), valid(old)
    _, ia, ib = np.intersect1d(a['key'].values, b['key'].values, return_indices=True)
    sh = a['median_log2dp'].values[ia] - b['median_log2dp'].values[ib]
    used = ((a['panel_median_dp'].values[ia] >= float(new['hdr'].get('capture_floor', 0.1))) &
            (b['panel_median_dp'].values[ib] >= float(old['hdr'].get('capture_floor', 0.1))))
    s2a, s2b = a['sigma2_dp'].values[ia], b['sigma2_dp'].values[ib]
    hn, ho = new['hdr'], old['hdr']
    out = {'new': new['path'], 'old': old['path'], 'n_samples_new': hn.get('n_samples'),
           'n_samples_old': ho.get('n_samples'), 'valid_new': len(a), 'valid_old': len(b), 'common': len(ia),
           'gained': len(a) - len(ia), 'lost': len(b) - len(ia)}
    for sfx, m in (('', np.ones(len(sh), bool)), ('_used', used)):
        x = sh[m]
        out.update({f'n{sfx}': int(m.sum()), f'shift_mean{sfx}': float(x.mean()), f'shift_sd{sfx}': float(x.std()),
                    f'shift_rsd{sfx}': rmad(x), f'frac_abs_gt_0.1{sfx}': float(np.mean(np.abs(x) > 0.1)),
                    f'frac_abs_gt_0.2{sfx}': float(np.mean(np.abs(x) > 0.2)),
                    f'frac_gt_0.2{sfx}': float(np.mean(x > 0.2)), f'frac_lt_-0.2{sfx}': float(np.mean(x < -0.2)),
                    f'p99_abs{sfx}': float(np.percentile(np.abs(x), 99))})
    pos = s2b > 0
    out.update({'sigma2_median_new': float(np.median(s2a)), 'sigma2_median_old': float(np.median(s2b)),
                'sigma2_ratio_of_medians': float(np.median(s2a) / np.median(s2b)),
                'sigma2_median_probe_ratio': float(np.median(s2a[pos] / s2b[pos]))})
    for tn, to in xy:
        rs = tn['hdr'].get('ratio_sex', '-')
        for lab, nm in (('chrX', XN), ('chrY', YN)):
            def ref(tb):
                d = tb['df']
                return set(d['key'].values[d['contig'].isin(nm).values & (d['probe_idx'].values >= 0)])
            kn, ko = ref(tn), ref(to)
            out.update({f'{lab}_ref_new[{rs}]': len(kn), f'{lab}_ref_old[{rs}]': len(ko),
                        f'{lab}_gained[{rs}]': len(kn - ko), f'{lab}_lost[{rs}]': len(ko - kn)})
    big = np.abs(sh) > 0.2
    shifted = []
    if big.any():
        idx = np.flatnonzero(big)
        ca, sa, ea = a['contig'].values[ia], a['start'].values[ia], a['end'].values[ia]
        new_ = np.r_[True, (ca[idx[1:]] != ca[idx[:-1]]) | (idx[1:] - idx[:-1] > 2) |
                     (sa[idx[1:]] - ea[idx[:-1]] > t['merge_gap_bp'])]
        for g in np.split(idx, np.flatnonzero(new_)[1:]):
            shifted.append({'contig': ca[g[0]], 'start': int(sa[g[0]]), 'end': int(ea[g[-1]]), 'n_shifted': len(g),
                            'mean_shift': float(sh[g].mean()), 'max_abs_shift': float(np.abs(sh[g]).max())})
    out['shifted_regions'] = len(shifted)
    out['expected_call_turnover'] = t['turnover_a'] + t['turnover_b'] * out['shift_rsd']
    return out, pd.DataFrame(shifted, columns=['contig', 'start', 'end', 'n_shifted', 'mean_shift', 'max_abs_shift'])


def merge_regions(p, flag, fl_se, fl_bi, t, prior):
    idx = np.flatnonzero(flag)
    cols = ['contig', 'start', 'end', 'n_flag', 'n_probes', 'reason', 'max_se', 'max_frac_lo', 'max_frac_hi',
            'prior_del_af', 'prior_dup_af', 'prior_an', 'common']
    if not len(idx):
        return pd.DataFrame(columns=cols)
    pos = np.arange(len(flag))
    new = np.ones(len(idx), bool)
    new[1:] = ((p['code'][idx[1:]] != p['code'][idx[:-1]]) | (idx[1:] - idx[:-1] - 1 > t['merge_gap_probes']) |
               (p['start'][idx[1:]] - p['end'][idx[:-1]] > t['merge_gap_bp']))
    gid = np.cumsum(new) - 1
    out = []
    for g in range(gid[-1] + 1):
        ii = idx[gid == g]
        a, z = ii[0], ii[-1]
        reason = '+'.join(x for x, m in (('se', fl_se[ii].any()), ('lo', (p['f_lo'][ii] >= t['bimodal_frac']).any()),
                                          ('hi', (p['f_hi'][ii] >= t['bimodal_frac']).any())) if m)
        out.append({'contig': p['contig'][a], 'start': int(p['start'][a]), 'end': int(p['end'][z]), 'n_flag': len(ii),
                    'n_probes': int(z - a + 1), 'reason': reason, 'max_se': float(p['se_boot'][ii].max()),
                    'max_frac_lo': float(p['f_lo'][ii].max()), 'max_frac_hi': float(p['f_hi'][ii].max())})
    reg = pd.DataFrame(out)
    del_af, dup_af, an = np.zeros(len(reg)), np.zeros(len(reg)), np.zeros(len(reg), int)
    if prior is not None:
        for j, (c, s, e) in enumerate(zip(reg['contig'], reg['start'], reg['end'])):
            seg = prior.get(c)
            if seg is None:
                continue
            a = np.searchsorted(seg[:, 1], s, 'right')
            z = np.searchsorted(seg[:, 0], e, 'left')
            hit = seg[a:z]
            hit = hit[hit[:, 2] >= t['prior_min_an']]
            if len(hit):
                del_af[j], dup_af[j], an[j] = hit[:, 3].max(), hit[:, 4].max(), int(hit[:, 2].min())
    reg['prior_del_af'], reg['prior_dup_af'], reg['prior_an'] = del_af, dup_af, an
    reg['common'] = np.maximum(del_af, dup_af) >= t['common_af']
    return reg[cols]


def panel_status(res, t):
    """Separate verdicts: autosomes, chrX/chrY of the male PoN, chrX of the female PoN (worst of their checks)."""
    c, p, m, b = res['composition'], res['precision'], res['members'], res['batch']
    br = m['by_rule']
    st = {}
    ne = c['n'] - len(c['empty'])
    st['n'] = 'pass' if ne >= t['n_pass'] else 'warn' if ne >= t['n_warn'] else 'fail'
    st['empty'] = 'warn' if c['empty'] else 'pass'
    ms = c['min_samples']
    st['min_samples'] = 'na' if ms is None else 'pass' if ms == c['min_samples_rule'] else 'warn'
    if ms is not None and not (c['min_samples_inferred'][0] <= ms <= c['min_samples_inferred'][1]):
        st['min_samples'] = 'warn'
    st['rel_se'] = status(p['rel_se_median'], t['rel_se_pass'], t['rel_se_warn'])
    st['members_mismatch'] = 'fail' if br['noise_extreme'] + br['missing'] + br.get('pon_divided', 0) else 'pass'
    st['members_outliers'] = ('warn' if (br['noise'] + br['burden'] + br['large_cnv']) / max(c['n'], 1) >
                              t['outlier_warn'] else 'pass')
    if b.get('pc1') is not None:
        st['batch_pc1'] = status(b['pc1_struct'], t['pc1_warn'], t['pc1_fail'])
        st['batch_sil'] = status(b['sil'], t['sil_warn'], t['sil_fail'])
    sexflag = m.get('sex_flagged', {})
    absent = set()
    for sx, lab in (('M', 'male'), ('F', 'female')):
        n = c.get(f'n_{sx}_with_probes', c[f'n_{sx}'])
        st[f'n_{lab}'] = 'pass' if n >= t['sex_pass'] else 'warn'
        if ms is not None and n < ms:
            st[f'n_{lab}'] = 'fail'                      # --min_samples keeps no chrX/chrY probe of this sex
        st[f'members_sex_{lab}'] = 'warn' if sexflag.get(sx, 0) else 'pass'
        st[f'chrX_ref_{lab}'] = 'na'
        for tb in c['tables']:
            if tb['ratio_sex'] == sx:
                st[f'chrX_ref_{lab}'] = 'fail' if tb['chrX_ref'] == 0 else 'pass'
        if n == 0 and st[f'chrX_ref_{lab}'] == 'na':    # no members with probes and no table of this sex
            absent.add(lab)
            st[f'n_{lab}'] = st[f'members_sex_{lab}'] = 'na'
    order = {'pass': 0, 'na': 0, 'warn': 1, 'fail': 2}
    worst = lambda ks: max((st[k] for k in ks if k in st), key=lambda v: order[v])
    st['autosomes'] = worst(['n', 'empty', 'min_samples', 'rel_se', 'members_mismatch', 'members_outliers',
                             'batch_pc1', 'batch_sil'])
    st['chrXY_male'] = 'na' if 'male' in absent else worst(['n_male', 'chrX_ref_male', 'min_samples',
                                                            'members_sex_male'])
    st['chrX_female'] = 'na' if 'female' in absent else worst(['n_female', 'chrX_ref_female', 'min_samples',
                                                                'members_sex_female'])
    return st


# ---- outputs -------------------------------------------------------------------------------------------------------
def report_rows(res, t):
    c, p, u, m, b, st = (res['composition'], res['precision'], res['unstable'], res['members'], res['batch'],
                         res['status'])
    R = []
    add = lambda chk, item, val, rule, s='', note='': R.append(
        {'check': chk, 'item': item, 'value': val, 'rule': rule, 'status': s, 'note': note})
    v = res.get('verify')
    if v:
        ok = v['median_mismatch_after_rounding'] == 0 and v['sigma2_mismatch_after_rounding'] == 0
        add('0_verify', 'median_max_abs_diff', f"{v['median_max_abs_diff']:.2e}", 'equal after float32 + %.6f',
            'pass' if ok else 'fail', f"mismatches: median {v['median_mismatch_after_rounding']}, sigma2 "
            f"{v['sigma2_mismatch_after_rounding']}, n_valid {v['n_valid_mismatch']}")
        add('0_verify', 'sigma2_max_abs_diff', f"{v['sigma2_max_abs_diff']:.2e}", '', 'pass' if ok else 'fail')
    ne = c['n'] - len(c['empty'])
    add('1_composition', 'n_members', c['n'], f">= {t['n_pass']} pass, >= {t['n_warn']} warn", st['n'],
        f'{ne} with probes')
    for sx, lab, basis in (('M', 'male', 'autosomal tiers by analogy, not measured'),
                           ('F', 'female', 'measured down to 10 females: no trend, wide intervals')):
        n = c.get(f'n_{sx}_with_probes', c[f'n_{sx}'])
        if st[f'n_{lab}'] == 'na':
            note = (f'no {lab} members with probes: a {lab} PoN would have no '
                    f'chrX{"/chrY" if sx == "M" else ""} reference')
        elif st[f'n_{lab}'] == 'fail':
            note = f'fewer than --min_samples {c["min_samples"]}: the {lab} PoN keeps no chrX/chrY probe'
        elif n >= t['sex_pass']:
            note = ''
        elif n >= t['sex_warn']:
            note = 'below the recommended count'
        else:
            note = ('below the floor; haploid calls on male chrX/chrY not recommended' if sx == 'M' else
                    'below the floor')
        add('1_composition', f'n_{sx} ({lab} PoN chrX{"/chrY" if sx == "M" else ""})', n,
            f">= {t['sex_pass']} pass, else warn ({basis})", st[f'n_{lab}'], note)
    add('1_composition', 'empty_members', len(c['empty']), '0', st['empty'], ','.join(c['empty']))
    add('1_composition', 'min_samples', c['min_samples'],
        f"= min(20, max(3, floor(n_smaller_sex/2))) = {c['min_samples_rule']}", st['min_samples'],
        f"table implies {c['min_samples_inferred'][0]}..{c['min_samples_inferred'][1]}")
    for tb in c['tables']:
        lab = {'M': 'male', 'F': 'female'}.get(tb['ratio_sex'])
        add('1_composition', f"chrX_ref/dark [{tb['ratio_sex'] or '-'}]", f"{tb['chrX_ref']}/{tb['chrX_dark']}",
            'chrX_ref > 0', st.get(f'chrX_ref_{lab}', '') if lab else '', os.path.basename(tb['path']))
        add('1_composition', f"chrY_ref/dark [{tb['ratio_sex'] or '-'}]", f"{tb['chrY_ref']}/{tb['chrY_dark']}", '', '')
    add('2_precision', 'se_boot_median', f"{p['se_boot_median']:.4f}", 'log2', '')
    add('2_precision', 'se_boot_p90/p99', f"{p['se_boot_p90']:.4f}/{p['se_boot_p99']:.4f}", '', '')
    add('2_precision', 'frac_se_boot_gt_0.05/0.1/0.2',
        '/'.join(f"{p[f'frac_se_boot_gt_{x}']:.4f}" for x in (0.05, 0.1, 0.2)), '', '')
    add('2_precision', 'se_analytic_median', f"{p['se_analytic_median']:.4f}", '1.2533 sigma / sqrt(n)', '')
    add('2_precision', 'frac_se_analytic_gt_0.05/0.1/0.2',
        '/'.join(f"{p[f'frac_se_analytic_gt_{x}']:.4f}" for x in (0.05, 0.1, 0.2)), '', '')
    add('2_precision', 'rel_se_median', f"{p['rel_se_median']:.4f}",
        f"<= {t['rel_se_pass']} pass, <= {t['rel_se_warn']} warn", st['rel_se'],
        f"clean panel of this size: {p['rel_se_clean']:.4f}")
    add('2_precision', 'sigma2_median', f"{p['sigma2_median']:.4f}", '', '')
    add('3_unstable', 'probes_flagged', u['probes_flagged'], f"SE > {t['se_unstable']} or >= {t['bimodal_frac']:.0%} "
        f"of members deviating on one side", '',
        f"SE {u['probes_se']}, two-state {u['probes_bimodal']} (unsmoothed rule {u['probes_offmedian_unsmoothed']})")
    add('3_unstable', 'regions (>=2 probes)', f"{u['regions']} ({u['regions_multi']})", 'BED', '',
        f"{u['regions_common']} at a common prior CNV (AF >= {t['common_af']})")
    br = m['by_rule']
    add('4_members', 'flagged', f"{m['n_flagged']} ({m['frac_flagged']:.1%})", 'see the rules below', '',
        '; '.join(f'{k} {v}' for k, v in br.items() if v))
    add('4_members', 'mismatch (noise_extreme, missing, pon_divided)',
        br['noise_extreme'] + br['missing'] + br.get('pon_divided', 0),
        f"0: noise > {t['noise_extreme']}x median, missing > {t['missing_floor']:.0%}, or a PoN-divided .norm",
        st['members_mismatch'], 'drop them, or use the raw .norm of a run without --pon_file')
    sf = m.get('sex_flagged', {})
    add('4_members', 'sex (mismatch, ambiguous)', br['sex_mismatch'] + br['sex_ambiguous'], '0 (status of the declared sex PoN)',
        'warn' if br['sex_mismatch'] + br['sex_ambiguous'] else 'pass',
        f"declared M {sf.get('M', 0)}, declared F {sf.get('F', 0)}; fix the sex list or drop them")
    add('4_members', 'outliers (noise, burden, large_cnv)', br['noise'] + br['burden'] + br['large_cnv'],
        f"<= {t['outlier_warn']:.0%} of members", st['members_outliers'], 'optional to drop')
    add('4_members', 'noise median/max', f"{m['noise_median']:.4f}/{m['noise_max']:.4f}", 'LOO 1.4826 MAD', '')
    if b.get('pc1') is not None:
        add('5_batch', 'pc1/pc2/pc3', f"{b['pc1']:.4f}/{b['pc2']:.4f}/{b['pc3']:.4f}", 'variance share', '')
        add('5_batch', 'pc1_struct', f"{b['pc1_struct']:.4f}", f"<= {t['pc1_warn']} pass, <= {t['pc1_fail']} warn",
            st['batch_pc1'], f"PC1 {b['pc1']:.4f} minus permutation null {b['pc1_null']:.4f}")
        add('5_batch', 'silhouette', f"{b['sil']:.3f}", f"<= {t['sil_warn']} pass, <= {t['sil_fail']} warn",
            st['batch_sil'], f"2-means on PC1-3, clusters {b['clusters']}")
    for r in res['tests']:
        add('6_test', r['sample'], f"{r['noise']:.4f}", f"noise <= {t['test_warn_ratio']}x members' median pass, "
            f"<= {t['noise_extreme']}x warn; missing <= {t['missing_floor']:.0%}", r['status'],
            f"pct {r['noise_pct']:.0f}, x{r['noise_over_median']:.2f} of median, missing {r['missing_frac']:.3f}, "
            f"sex {r['sex_inferred']}" + (", PoN-divided .norm" if r.get('one_share', 0) > t['pon_divided_one'] else ''))
    for r in res.get('admit', []):
        add('7_admit', r['sample'], r['verdict'], 'reject: noise > 2.5x median, missing > 3%, sex ambiguous; review: a '
            'member flag (noise, burden, large_cnv, sex_mismatch)', r['verdict'],
            f"{r['rules'] or '-'}; noise {r.get('noise', np.nan):.3f} (pct {r.get('noise_pct', np.nan):.0f}, "
            f"x{r.get('noise_over_median', np.nan):.2f}), longest run {r.get('max_run', 0)} {r.get('max_run_region', '')}")
    add('verdict', 'autosomes', st['autosomes'], 'N, rel SE, empty, min_samples, member mismatch/outliers, batch',
        st['autosomes'])
    add('verdict', 'chrX/chrY (male PoN)', st['chrXY_male'], 'males with probes, male chrX_ref, min_samples, sex flags',
        st['chrXY_male'])
    add('verdict', 'chrX (female PoN)', st['chrX_female'], 'females with probes, female chrX_ref, min_samples, sex flags',
        st['chrX_female'])
    return pd.DataFrame(R)


def compare_rows(cmp):
    keys = ['n_samples_new', 'n_samples_old', 'valid_new', 'valid_old', 'common', 'gained', 'lost', 'shift_mean',
            'shift_sd', 'shift_rsd', 'frac_abs_gt_0.1', 'frac_abs_gt_0.2', 'frac_gt_0.2', 'frac_lt_-0.2', 'p99_abs',
            'shift_sd_used', 'shift_rsd_used', 'frac_abs_gt_0.2_used', 'sigma2_median_new', 'sigma2_median_old',
            'sigma2_ratio_of_medians', 'sigma2_median_probe_ratio', 'shifted_regions', 'expected_call_turnover']
    keys += [k for k in cmp if k.startswith(('chrX_', 'chrY_'))]
    fmt = lambda v: f'{v:.4g}' if isinstance(v, float) else v
    return pd.DataFrame([{'check': '8_compare', 'item': k, 'value': fmt(cmp[k]), 'rule': '', 'status': '', 'note': ''}
                         for k in keys if k in cmp])


def summary_text(res, rep):
    c = res['composition']
    st = res['status']
    lines = [f"PoN QC  {c['n']} members ({c['n_M']} M, {c['n_F']} F)  table {os.path.basename(res['table'])}  "
             f"autosomes {st['autosomes'].upper()}  chrX/chrY(M) {st['chrXY_male'].upper()}  "
             f"chrX(F) {st['chrX_female'].upper()}"]
    for _, r in rep.iterrows():
        s = r['status'].upper() if r['status'] else ''
        note = f"  {r['note']}" if r['note'] else ''
        lines.append(f"{s:5s} {r['check']:14s} {r['item']:34s} {r['value']!s:22s}{note}"[:160])
    fl = res['members']['flagged']
    if fl:
        lines.append('flagged members: ' + ', '.join(f'{k} ({v})' for k, v in list(fl.items())[:12]) +
                     (' ...' if len(fl) > 12 else ''))
    return '\n'.join(lines)


def jsonable(o):
    if isinstance(o, dict):
        return {k: jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [jsonable(v) for v in o]
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating, float)):
        return None if np.isnan(o) else float(o)
    if isinstance(o, np.bool_):
        return bool(o)
    return o


def write_outputs(prefix, out, cmp=None, t=T):
    res = out['res']
    rep = report_rows(res, t)
    if cmp is not None:
        rep = pd.concat([rep, compare_rows(cmp[0])], ignore_index=True)
    rep.to_csv(prefix + '.report.tsv', sep='\t', index=False)
    out['members'].drop(columns=['path']).to_csv(prefix + '.members.tsv', sep='\t', index=False, float_format='%.5f')
    with open(prefix + '.unstable.bed', 'w') as f:
        f.write('#' + '\t'.join(out['regions'].columns) + '\n')
        out['regions'].to_csv(f, sep='\t', index=False, header=False, float_format='%.4f')
    for k, sfx in (('tests', '.tests.tsv'), ('admit', '.admit.tsv'), ('runs', '.admit_runs.tsv')):
        x = out[k]
        if len(x):
            x.drop(columns=['path'], errors='ignore').to_csv(prefix + sfx, sep='\t', index=False, float_format='%.5f')
    js = dict(res)
    js['batch'] = {k: v for k, v in res['batch'].items() if k != 'labels'}
    if cmp is not None:
        js['compare'] = cmp[0]
        write_compare(prefix, cmp)
    with open(prefix + '.json', 'w') as f:
        json.dump(jsonable(js), f, indent=1)
    txt = summary_text(res, rep)
    with open(prefix + '.txt', 'w') as f:
        f.write(txt + '\n')
    return txt


def write_compare(prefix, cmp):
    c, shifted = cmp
    compare_rows(c)[['item', 'value']].to_csv(prefix + '.compare.tsv', sep='\t', index=False)
    with open(prefix + '.compare.bed', 'w') as f:
        f.write('#' + '\t'.join(shifted.columns) + '\n')
        shifted.to_csv(f, sep='\t', index=False, header=False, float_format='%.4f')


def compare_text(c):
    nm = lambda p: os.path.join(os.path.basename(os.path.dirname(os.path.abspath(p))), os.path.basename(p))
    lines = [f"PoN compare  new {nm(c['new'])} ({c['n_samples_new']} samples)  old {nm(c['old'])} "
             f"({c['n_samples_old']} samples)"]
    for k in ('common', 'gained', 'lost', 'shift_sd', 'shift_rsd', 'frac_abs_gt_0.1', 'frac_abs_gt_0.2', 'frac_gt_0.2',
              'sigma2_ratio_of_medians', 'sigma2_median_probe_ratio', 'shifted_regions', 'expected_call_turnover'):
        v = c[k]
        lines.append(f"      8_compare      {k:34s} {v:.4f}" if isinstance(v, float) else f"      8_compare      {k:34s} {v}")
    return '\n'.join(lines)


def read_admit_sex(vals, paths):
    """--admit_sex: one M/F per --admit sample, or one file of path<TAB>M|F lines."""
    if len(vals) == 1 and os.path.exists(vals[0]):
        return read_sex(vals[0])
    if vals and len(vals) != len(paths):
        sys.exit(f'--admit_sex: {len(vals)} values for {len(paths)} --admit samples')
    return {p: v.upper()[:1] for p, v in zip(paths, vals)}


def run_build_qc(input_list, sample_sex, par_bed, pon_tsv_paths, min_samples, prior=None, out_prefix='pon_qc',
                 boot=200, jobs=1, seed=20261110, pca_probes=20000, tests=(), admit=(), admit_sex=None, compare=None,
                 verbose=False):
    """QC of a freshly built PoN, for create_pon.py to call after it writes the PoN.

    input_list, sample_sex, par_bed and min_samples are the builder's inputs; pon_tsv_paths are the --tsv tables it
    wrote (male and/or female); prior is the population CNV prior VCF or BUNDLE/MEMBER of a model bundle (read_prior),
    or None (no prior AF annotation), and the JSON's "prior" records the path, member and md5 read; boot is the
    number of bootstrap draws; jobs the processes reading the .norm files. Optional: tests and admit (raw .norm paths),
    admit_sex ({path: 'M'|'F'}), compare (previous .pon.tsv paths). Writes out_prefix.report.tsv, .members.tsv,
    .unstable.bed, .json and .txt (and .tests.tsv, .admit.tsv, .admit_runs.tsv, .compare.tsv, .compare.bed when asked)
    and nothing else; logs to stderr only with verbose. Returns the status dict: one status per check and the verdicts
    'autosomes', 'chrXY_male' and 'chrX_female' (pass, warn or fail; a sex verdict is na when the panel has no members
    of that sex with probes and no table of that sex is given). Importers set sys.dont_write_bytecode = True
    before `import pon_qc` (or run python -B): the guard in this module cannot stop its own .pyc.
    """
    global VERBOSE
    old_v, VERBOSE = VERBOSE, verbose
    try:
        if prior:
            prior_member(prior)  # a bad prior fails here, before any table or .norm is read
        with warnings.catch_warnings(), np.errstate(all='ignore'):
            warnings.simplefilter('ignore', RuntimeWarning)
            tables = [read_pon_tsv(p) for p in pon_tsv_paths]
            tables.sort(key=lambda tb: tb['hdr'].get('ratio_sex') != 'M')
            cmp = None
            if compare:
                old = [read_pon_tsv(p) for p in compare]
                bysex = {tb['hdr'].get('ratio_sex', ''): tb for tb in old}
                old.sort(key=lambda tb: tb['hdr'].get('ratio_sex') != 'M')
                xy = [(tb, bysex[tb['hdr'].get('ratio_sex', '')]) for tb in tables if tb['hdr'].get('ratio_sex', '') in bysex]
                cmp = compare_tables(tables[0], old[0], T, xy)
            paths = read_list(input_list)
            sx = read_sex(sample_sex)
            miss = [p for p in paths if p not in sx]
            if miss:
                raise ValueError(f'no sex annotation for {len(miss)} members, e.g. {miss[0]}')
            store = NormStore(np.unique(np.concatenate([tb['df']['key'].values for tb in tables])))
            log(f'reading {len(paths)} members')
            store.load(paths + list(tests) + list(admit), jobs=jobs)
            src = {}
            pr = read_prior(prior, src) if prior else None
            if prior:
                log(f"prior {src['path']}, md5 {src['md5']}")
            out = run_qc(paths, [sx[p] for p in paths], tables, store, read_bed(par_bed), min_samples=min_samples,
                         prior=pr, B=boot, seed=seed, n_pca=pca_probes, tests=tests,
                         admits=admit, admit_sex=admit_sex or {})
            out['res']['prior'] = src or None
            if cmp is not None:
                out['res']['compare'] = cmp[0]
            d = os.path.dirname(out_prefix)
            if d:
                os.makedirs(d, exist_ok=True)
            write_outputs(out_prefix, out, cmp)
            return dict(out['res']['status'])
    finally:
        VERBOSE = old_v


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--input_list', help='panel: one raw .norm path per line (omit for --compare alone)')
    ap.add_argument('--sample_sex', help='panel: <path>\\t<M|F> per member')
    ap.add_argument('--par_bed', help='PAR bed (config/hs38.PAR.bed)')
    ap.add_argument('--pon', required=True, nargs='+', help='.pon.tsv table(s) built from this panel')
    ap.add_argument('--out', required=True, help='output prefix')
    ap.add_argument('--min_samples', type=int, default=None, help='the --min_samples the PoN was built with')
    ap.add_argument('--prior', default=None, help='population CNV prior VCF (INFO END, AN, DEL_AC, DUP_AC), or '
                    'BUNDLE/MEMBER, a VCF member of a model bundle, e.g. CNV_WES.bundle/cnvprior.vcf.gz')
    ap.add_argument('--test', nargs='*', default=[], help='raw .norm files of samples to check against the panel')
    ap.add_argument('--admit', nargs='*', default=[], help='raw .norm files of called samples to check for the pool')
    ap.add_argument('--admit_sex', nargs='*', default=[], help='M/F per --admit sample, or a path<TAB>M|F file')
    ap.add_argument('--compare', nargs='+', default=None, help='previous .pon.tsv table(s) to compare --pon with')
    ap.add_argument('--bootstrap', type=int, default=200)
    ap.add_argument('--seed', type=int, default=20261110)
    ap.add_argument('--pca_probes', type=int, default=20000)
    ap.add_argument('--jobs', type=int, default=1, help='processes for reading the .norm files')
    a = ap.parse_args()
    if a.input_list:
        if not (a.sample_sex and a.par_bed):
            sys.exit('--input_list needs --sample_sex and --par_bed')
        try:
            run_build_qc(a.input_list, a.sample_sex, a.par_bed, a.pon, a.min_samples, prior=a.prior, out_prefix=a.out,
                         boot=a.bootstrap, jobs=a.jobs, seed=a.seed, pca_probes=a.pca_probes, tests=a.test,
                         admit=a.admit, admit_sex=read_admit_sex(a.admit_sex, a.admit), compare=a.compare, verbose=True)
        except ValueError as e:
            sys.exit(str(e))
        with open(a.out + '.txt') as f:
            print(f.read(), end='')
        return
    if not a.compare:
        sys.exit('--input_list (with --sample_sex and --par_bed) or --compare is required')
    d = os.path.dirname(a.out)
    if d:
        os.makedirs(d, exist_ok=True)
    new_t = sorted((read_pon_tsv(p) for p in a.pon), key=lambda tb: tb['hdr'].get('ratio_sex') != 'M')
    old_t = sorted((read_pon_tsv(p) for p in a.compare), key=lambda tb: tb['hdr'].get('ratio_sex') != 'M')
    bysex = {tb['hdr'].get('ratio_sex', ''): tb for tb in old_t}
    xy = [(tb, bysex[tb['hdr'].get('ratio_sex', '')]) for tb in new_t if tb['hdr'].get('ratio_sex', '') in bysex]
    cmp = compare_tables(new_t[0], old_t[0], T, xy)
    write_compare(a.out, cmp)
    with open(a.out + '.json', 'w') as f:
        json.dump(jsonable({'compare': cmp[0]}), f, indent=1)
    txt = compare_text(cmp[0])
    with open(a.out + '.txt', 'w') as f:
        f.write(txt + '\n')
    print(txt)


if __name__ == '__main__':
    main()
