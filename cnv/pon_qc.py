#!/usr/bin/env python3
"""pon_qc.py -- checks for a wes2.1 CNVscope panel of normals (PoN) and for the samples called with it.

usage:
  python3 -I pon_qc.py --screen SAMPLE.norm --pon LAB.<sex>.pon [--out PREFIX]
  python3 -I pon_qc.py --pon NEW.male.pon [NEW.female.pon] --compare OLD.male.pon [OLD.female.pon] --out PREFIX

--screen   one called sample against the binary PoN its call used. SAMPLE.norm is the call's own depth file (CNVscope
           --output_normalized with --pon_file: PoN-divided, masked probes exactly 1, capped ones 8/m) or the depth file
           of a run without --pon_file (PoN-free); more than 1% of probes exactly 1.0 marks a PoN-divided file.
           On the autosomal probes CNVscope uses (referenced by the PoN, PoN median depth >= 0.1) the noise is the
           1.4826 MAD of the log2 ratio to the PoN, as a multiple of the members' median leave-one-out noise that the
           PoN's per-probe sigma2 predicts. The sex comes from the chrX non-PAR dosage against the autosomes and chrY.
           fail: noise > 2.5x, more than 3% of those probes missing, or a sex other than the PoN's or XO-/XXY-like
           (checked at noise <= 2.5x); warn: noise > 1.75x; pass otherwise. fail: the calls are unreliable and the
           sample is not for the pool; pass or warn: a candidate for the pool. One line on stdout, PREFIX.tsv with
           --out; exit 0 for any verdict, 1 for an input error.
--compare  the stored medians and sigma2 of PoN NEW against OLD on the autosomal probes valid in both (_used: those
           CNVscope uses in both), and the chrX/chrY probes each sex references. Writes PREFIX.compare.tsv,
           PREFIX.compare.bed (regions shifted > 0.2 log2), PREFIX.json and PREFIX.txt.

run_build_qc() is the QC report create_pon.py runs on the PoN it has just written, from the binary PoN and the members'
depth files (thresholds in T, calibrated in ponsize/analysis/qc/ by pon_qc_calib.py):
  0 verification   this script's full-panel median and sigma2 against the stored ones (valid autosomal probes)
  1 composition    N, sexes, empty members, chrX/chrY counts, --min_samples against the builder's auto rule
                   min(20, max(3, floor(n_smaller_sex/2))): members with probes; a sex with none does not count
  2 precision      per-probe SE of the panel median: bootstrap over members (B=200) and 1.2533 sigma / sqrt(n)
  3 unstable loci  probes with bootstrap SE > 0.2, or with >= 10% of members deviating on one side, merged into
                   regions with the prior's DEL/DUP AF (BED)
  4 members        leave-one-out (LOO) noise, CNV burden (5-probe running median), sex, missing probes
  5 batch          PCA of the noise-scaled LOO profiles on 20k probes; PC1 above a permutation null, 2-means silhouette
Verdicts: autosomes; chrX/chrY of a male PoN; chrX of a female PoN (na for a sex with no members with probes, unless the
PoN is of that sex). Members are depth files of runs without --pon_file (a PoN-divided one is flagged pon_divided).
Transform and estimators are the builder's: log2(DP + 1e-6), column medians (mean of the two middle values),
sigma2 = (1.4826 MAD)^2; n_valid and the panel median DP are recomputed from the members as the builder makes them.
Noise, burden, PCA and unstable flags use the valid autosomal probes with panel median DP >= --capture_floor (0.1).
The prior is a .vcf.gz, or BUNDLE/MEMBER, a member of a model bundle (copied to a temp file in $TMPDIR, read and
removed). Writes PREFIX.report.tsv, .members.tsv, .unstable.bed, .json (with the prior's path, member and md5), .txt.
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
import struct  # noqa: E402
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
    'missing_scale': 0.005, 'missing_floor': 0.03,  # missing > 3% of the probes CNVscope uses (members,
                                             # --screen)
    'dev_log2': 0.5, 'dev_k': 4.0,           # a member deviates where |5-probe running median of its LOO ratio|
                                             # > max(0.5, 4 x its robust sd): CNV burden and two-state counts
    'run_log2': 0.3, 'run_k': 3.0,           # CNV run: consecutive probes with the 5-probe running median beyond
                                             # max(0.3, 3 x its robust sd) on one side, outside common CNV loci
    'large_cnv_probes': 100,                 # member flag: a run of >= 100 probes (healthy
                                             # members' longest run: median 20-33, p99 72-111, max 90-174)
    'x_one_copy': -0.5, 'y_present': -2.0,   # sex: chrX / chrY log2 depth against the panel's same-sex reference
    'pc1_warn': 0.15, 'pc1_fail': 0.25,      # PC1 variance share above its permutation null (noise-scaled profiles)
    'sil_warn': 0.60, 'sil_fail': 0.70,      # 2-means silhouette on PC1-3
    'test_warn_ratio': 1.75,                 # --screen: warn above 1.75x the members' median noise
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


# create_pon.py's per-probe record (43 bytes, packed); svd_idx is read as probe_idx, flags as par_flag
PON_DT = np.dtype([('cidx', '<u2'), ('start', '<u4'), ('end', '<u4'), ('median_log2dp', '<f4'),
                   ('median_log2dp0', '<f4'), ('probe_idx', '<i4'), ('sigma2_dp', '<f4'), ('median_dp_male', '<f4'),
                   ('median_dp_female', '<f4'), ('sigma2_dp_male', '<f4'), ('sigma2_dp_female', '<f4'),
                   ('par_flag', 'u1')])
PON_MIN_MEDIAN = 0.1                 # CNVscope's pon_min_median: it masks probes whose PoN median depth is below it


def read_pon(path, sex=None):
    """A binary PoN (CNVPON v1, as create_pon.py writes it): {'path', 'hdr', 'df', 'sex'}.

    df: one row per probe in file order: contig, start, end, the stored fields as float64 (probe_idx is svd_idx: -1
    dark, -2 non-autosomal and not referenced, >= 0 referenced), excluded (1 where probe_idx < 0: the PoN does not use
    the probe) and key. hdr: version, n_samples, n_males, n_females, median_variance, n_svd_probes, k. sex: the
    caller's; else 'M' when the PoN references a chrY probe, 'F' when it references chrX and no chrY, None otherwise.
    The V_k matrices are not read."""
    with open(path, 'rb') as f:
        h = f.read(64)
        if len(h) < 64 or h[:6] != b'CNVPON':
            raise ValueError(f'{path}: not a binary PoN (no CNVPON header)')
        ver, n_svd, k_dp, k_dp0, n_probes, n_contigs, n_samples, mvar = struct.unpack_from('<HIHHIHHf', h, 6)
        if ver != 1:
            raise ValueError(f'{path}: PoN version {ver}; this script reads version 1')
        n_males, n_females = struct.unpack_from('<HH', h, 32)
        names, off = [], 64
        for _ in range(n_contigs):
            b = f.read(2)
            n = struct.unpack('<H', b)[0] if len(b) == 2 else 0
            names.append(f.read(n).decode('utf-8', 'replace'))
            off += 2 + n
        size, want = os.path.getsize(path), off + PON_DT.itemsize * n_probes + 4 * n_svd * (k_dp + k_dp0)
        if size != want:
            raise ValueError(f'{path}: {size} bytes, but its header implies {want}')
        tab = np.frombuffer(f.read(PON_DT.itemsize * n_probes), PON_DT)
    if n_probes and int(tab['cidx'].max()) >= n_contigs:
        raise ValueError(f'{path}: a probe names contig {int(tab["cidx"].max())} of {n_contigs}')
    df = pd.DataFrame({'contig': np.asarray(names, dtype=object)[tab['cidx']],
                       'start': tab['start'].astype(np.int64), 'end': tab['end'].astype(np.int64)})
    for c in PON_DT.names[3:]:
        df[c] = tab[c].astype(np.int64 if c in ('probe_idx', 'par_flag') else np.float64)
    df['excluded'] = (df['probe_idx'].values < 0).astype(np.int64)
    df['key'] = pack(contig_code(df['contig'].values), df['start'].values, df['end'].values)
    if sex is None:
        ref = df['probe_idx'].values >= 0
        sex = ('M' if (ref & df['contig'].isin(YN).values).any() else
               'F' if (ref & df['contig'].isin(XN).values).any() else None)
    hdr = {'version': ver, 'n_samples': n_samples, 'n_males': n_males, 'n_females': n_females,
           'median_variance': mvar, 'n_svd_probes': n_svd, 'k': k_dp}
    return {'path': path, 'hdr': hdr, 'df': df, 'sex': sex}


def cnvscope_on(df):
    """The probes CNVscope does not mask: referenced by the PoN, with PoN median depth >= pon_min_median. The others
    are exactly 1 in a PoN-divided .norm."""
    return (df['probe_idx'].values >= 0) & (df['median_log2dp'].values >= np.log2(PON_MIN_MEDIAN))


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


def member_stats(tb, store, paths, chunk=50000):
    """n_valid (members with a value) and panel_median_dp (their median DP) of every probe of PoN tb, added to its df,
    as create_pon.py makes them over all its members."""
    df = tb['df']
    cols = store.cols(df['key'].values)
    nv, mdp = np.zeros(len(df), np.int64), np.full(len(df), np.nan)
    for a in range(0, len(df), chunk):
        M = store.matrix(paths, cols[a:a + chunk])
        nv[a:a + chunk] = np.sum(~np.isnan(M), axis=0)
        mdp[a:a + chunk] = fast_col_median(M, nv[a:a + chunk])
    df['n_valid'], df['panel_median_dp'] = nv, mdp


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


def status(v, warn, fail, high_bad=True):
    if v is None or (isinstance(v, float) and np.isnan(v)):
        return 'na'
    if high_bad:
        return 'fail' if v > fail else 'warn' if v > warn else 'pass'
    return 'fail' if v < fail else 'warn' if v < warn else 'pass'


# ---- the QC --------------------------------------------------------------------------------------------------------
def run_qc(paths, sexes, tables, store, par, min_samples=None, prior=None, B=200, seed=20261110, n_pca=20000,
           cf=PON_MIN_MEDIAN, verify=True, opt=None):
    """Checks 0-5 for one panel. tables: read_pon() of its PoN, with member_stats(); cf: the capture floor. Returns a
    dict: res, probes, members, regions."""
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

    # 0 verification against the stored float32 median and sigma2
    if verify:
        with np.errstate(all='ignore'):
            med64 = np.nanmedian(L, axis=0)
        s2, _ = mad_sigma2(L, nvd)
        st_med, st_s2 = A['median_log2dp'].values, A['sigma2_dp'].values
        res['verify'] = {
            'median_max_abs_diff': float(np.nanmax(np.abs(med64 - st_med))),
            'median_mismatch_after_rounding': int(np.sum(med64.astype(np.float32) != st_med)),
            'sigma2_max_abs_diff': float(np.nanmax(np.abs(s2 - st_s2))),
            'sigma2_mismatch_after_rounding': int(np.sum(s2.astype(np.float32) != st_s2)),
        }
        del s2
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
        row = {'path': tb['path'], 'ratio_sex': tb['sex'] or '', 'n_samples': tb['hdr'].get('n_samples'),
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
              'med': med, 'sigma2': sig ** 2, 'se_boot': se, 'se_an': se_an, 'nv': nv,
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

    res['seconds'] = round(time.time() - t0, 1)
    res['status'] = panel_status(res, t)
    return {'res': res, 'probes': probes, 'members': mem, 'regions': regions}


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


# ---- the screen ----------------------------------------------------------------------------------------------------
REF_K = 1.0076                    # members' median noise / estimate: median over the study's 22 panels (sd 3.6%)


def ref_noise(pon):
    """The members' median LOO noise that the PoN predicts: REF_K x the median sigma of the autosomal probes CNVscope
    uses x sqrt(1 + pi / (2 (N - 1))), the LOO median's variance."""
    df = pon['df']
    u = cnvscope_on(df) & df['contig'].isin(AUTO).values
    sig = np.sqrt(np.maximum(df['sigma2_dp'].values[u], 0))
    n = max(int(pon['hdr']['n_samples']) - 1, 1)
    return REF_K * float(np.median(sig)) * np.sqrt(1 + np.pi / 2 / n) if len(sig) else np.nan


def screen(norm_path, pon, t=T):
    """--screen: a called sample's depth file against the binary PoN its call used (read_pon). Returns the PREFIX.tsv
    row as a dict; raises ValueError when no probe of the file is in the PoN."""
    df = pon['df']
    k, dp = read_norm(norm_path)
    pk = df['key'].values
    o = np.argsort(pk, kind='stable')
    i = np.minimum(np.searchsorted(pk[o], k), len(pk) - 1)
    ok = pk[o][i] == k
    if not ok.any():
        raise ValueError(f'{norm_path}: none of its {len(k)} probes is in the PoN {pon["path"]}')
    x = np.full(len(df), np.nan)
    x[o[i[ok]]] = dp[ok]
    contig, med = df['contig'].values, df['median_log2dp'].values
    auto, isx, isy = np.isin(contig, list(AUTO)), np.isin(contig, XN), np.isin(contig, YN)
    on = cnvscope_on(df)
    used = on & auto
    have = ~np.isnan(x)
    va = auto & (df['probe_idx'].values >= 0) & have
    one_share = float(np.mean(x[va] == 1.0)) if va.any() else np.nan
    divided = one_share > t['pon_divided_one']
    with np.errstate(divide='ignore', invalid='ignore'):
        if divided:                                   # DP = min(DP_raw / 2^median, 8) / m: masked rows exactly 1
            masked = x == 1.0
            v = x[on & have & ~masked]
            top = v.max() if len(v) else np.nan
            capped = x == top if np.sum(v == top) >= 2 else np.zeros(len(x), bool)   # the ratio cap 8/m
            L = np.where(have & ~masked & ~capped, np.log2(x), np.nan)              # log2 ratio - log2(m)
        else:
            masked = capped = np.zeros(len(x), bool)
            L = np.log2(x + PSEUDO) - med
    hu = used & ~np.isnan(L)
    n_used = int(used.sum())
    n_miss = int(np.sum(used & ~(have & ~masked)))
    noise = rmad(L[hu]) if hu.sum() >= 100 else np.nan
    off = float(np.median(L[hu])) if hu.any() else np.nan
    ref = ref_noise(pon)
    # sex: chrX non-PAR and chrY against the PoN's two-copy reference; without one (PoN-free input only), the stored
    # same-sex medians (female chrX, else male + 1; male chrY at the one-copy level)
    hdr, lo, npar = pon['hdr'], np.log2(PON_MIN_MEDIAN), df['par_flag'].values == 0
    xs, ys, xl, yl, yshift = isx & npar & on, isy & on, L, L, 1.0
    if not xs.any() and not divided:
        alt = df['median_dp_female'].values if hdr['n_females'] else df['median_dp_male'].values + 1.0
        with np.errstate(invalid='ignore'):
            xs = isx & npar & (alt >= lo)
            xl = np.log2(x + PSEUDO) - alt
    if not ys.any() and not divided and hdr['n_males']:
        mm = df['median_dp_male'].values
        with np.errstate(invalid='ignore'):
            ys = isy & npar & (mm >= lo)
        yl, yshift = np.log2(x + PSEUDO) - mm, 0.0
    x_l2 = nanmed(xl[xs]) - off if xs.any() and not np.isnan(off) else np.nan
    y_l2 = nanmed(yl[ys]) - off + yshift if ys.any() and not np.isnan(off) else np.nan
    sx = infer_sex(x_l2, y_l2, t)
    ratio = noise / ref if ref > 0 else np.nan
    fails, notes = [], []
    if np.isnan(ratio):
        fails.append('noise')
        notes.append(f'fewer than 100 of the {n_used} probes CNVscope uses have a value')
    elif ratio > t['noise_extreme']:
        fails.append('noise')
        notes.append(f"noise {ratio:.2f}x the PoN's members (> {t['noise_extreme']}x): another centre, kit or batch")
    if n_used and n_miss / n_used > t['missing_floor']:
        fails.append('missing')
        notes.append(f"{n_miss / n_used:.1%} of the probes CNVscope uses are missing (> {t['missing_floor']:.0%})")
    ps = pon['sex']
    if ps and not ratio <= t['noise_extreme']:       # another centre, kit or batch shifts chrX dosage by up to 0.5 log2
        notes.append('sex not checked at this noise')
    elif ps and sx in ('M', 'F') and sx != ps:
        fails.append('sex_mismatch')
        other = {'M': 'male', 'F': 'female'}[sx]
        notes.append(f"the sample looks {other} and the PoN is {'male' if ps == 'M' else 'female'}: re-call it with "
                     f"the {other} kit files and PoN")
    elif ps and sx in ('XO-like', 'XXY-like'):
        fails.append('sex_aneuploid')
        notes.append(f'{sx} sex chromosomes')
    verdict = 'fail' if fails else 'warn' if ratio > t['test_warn_ratio'] else 'pass'
    if verdict == 'warn':
        notes.append(f"noise {ratio:.2f}x the PoN's members (> {t['test_warn_ratio']}x)")
    return {'sample': sample_name(norm_path), 'verdict': verdict, 'reasons': ','.join(fails),
            'noise_ratio': ratio, 'noise': noise, 'ref_noise': ref, 'missing_frac': n_miss / n_used if n_used else np.nan,
            'n_used': n_used, 'n_missing': n_miss, 'n_capped': int(np.sum(capped & used)),
            'sex_inferred': 'not checked' if ps and not ratio <= t['noise_extreme'] else sx,   # unreliable there
            'sex_basis': 'NA' if np.isnan(x_l2) else 'X only' if np.isnan(y_l2) else 'X+Y',
            'x_l2': x_l2, 'y_l2': y_l2, 'pon_sex': ps or 'none',
            'input': 'PoN-divided' if divided else 'PoN-free', 'one_share': one_share, 'n_probes_file': len(k),
            'pon': pon['path'], 'norm': norm_path, 'note': '; '.join(notes)}


def screen_line(r):
    pct = lambda v: 'NA' if np.isnan(v) else f'{100 * v:.2f}%'
    ratio = 'NA' if np.isnan(r['noise_ratio']) else f"{r['noise_ratio']:.2f}x"
    line = (f"{r['sample']}  {r['verdict'].upper()}  noise {ratio}  missing {pct(r['missing_frac'])}  "
            f"sex {r['sex_inferred']} (PoN {r['pon_sex']})  input {r['input']}")
    return line + (f"  -- {r['note']}" if r['note'] else '')


# ---- compare -------------------------------------------------------------------------------------------------------
def compare_tables(new, old, t, xy=()):
    """Stored medians and sigma2 of PoN `new` against `old` (read_pon) over the autosomal probes valid in both."""
    def valid(tb):
        d = tb['df']
        return d[(d['excluded'].values == 0) & d['contig'].isin(AUTO).values]

    a, b = valid(new), valid(old)
    _, ia, ib = np.intersect1d(a['key'].values, b['key'].values, return_indices=True)
    sh = a['median_log2dp'].values[ia] - b['median_log2dp'].values[ib]
    lo = np.log2(PON_MIN_MEDIAN)
    used = (a['median_log2dp'].values[ia] >= lo) & (b['median_log2dp'].values[ib] >= lo)   # CNVscope uses them
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
        rs = tn['sex'] or '-'
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
        add('0_verify', 'median_max_abs_diff', f"{v['median_max_abs_diff']:.2e}", 'equal after float32',
            'pass' if ok else 'fail', f"mismatches: median {v['median_mismatch_after_rounding']}, sigma2 "
            f"{v['sigma2_mismatch_after_rounding']}")
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


def write_outputs(prefix, out, t=T):
    res = out['res']
    rep = report_rows(res, t)
    rep.to_csv(prefix + '.report.tsv', sep='\t', index=False)
    out['members'].drop(columns=['path']).to_csv(prefix + '.members.tsv', sep='\t', index=False, float_format='%.5f')
    with open(prefix + '.unstable.bed', 'w') as f:
        f.write('#' + '\t'.join(out['regions'].columns) + '\n')
        out['regions'].to_csv(f, sep='\t', index=False, header=False, float_format='%.4f')
    js = dict(res)
    js['batch'] = {k: v for k, v in res['batch'].items() if k != 'labels'}
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


def run_build_qc(input_list, sample_sex, par_bed, pon, min_samples, ratio_sex=None, capture_floor=PON_MIN_MEDIAN,
                 prior=None, out_prefix='pon_qc', boot=200, jobs=1, seed=20261110, pca_probes=20000, verbose=False):
    """QC of a freshly built PoN, for create_pon.py to call after it writes the PoN.

    input_list, sample_sex, par_bed, min_samples, ratio_sex ('M', 'F' or None) and capture_floor are the builder's
    inputs; pon is the binary PoN it wrote; prior is the population CNV prior VCF or BUNDLE/MEMBER of a model bundle
    (read_prior), or None (no prior AF annotation), and the JSON's "prior" records the path, member and md5 read; boot
    is the number of bootstrap draws; jobs the processes reading the .norm files. Writes out_prefix.report.tsv,
    .members.tsv, .unstable.bed, .json and .txt and nothing else; logs to stderr only with verbose. Returns the status
    dict: one status per check and the verdicts 'autosomes', 'chrXY_male' and 'chrX_female' (pass, warn or fail; a sex
    verdict is na when the panel has no members of that sex with probes and the PoN is not of that sex). Importers set
    sys.dont_write_bytecode = True before `import pon_qc` (or run python -B): the guard in this module cannot stop its
    own .pyc.
    """
    global VERBOSE
    old_v, VERBOSE = VERBOSE, verbose
    try:
        if prior:
            prior_member(prior)  # a bad prior fails here, before the PoN or any .norm is read
        with warnings.catch_warnings(), np.errstate(all='ignore'):
            warnings.simplefilter('ignore', RuntimeWarning)
            tb = read_pon(pon)
            tb['sex'] = ratio_sex or None               # the builder's --ratio_sex, not inferred
            paths = read_list(input_list)
            sx = read_sex(sample_sex)
            miss = [p for p in paths if p not in sx]
            if miss:
                raise ValueError(f'no sex annotation for {len(miss)} members, e.g. {miss[0]}')
            store = NormStore(np.unique(tb['df']['key'].values))
            log(f'reading {len(paths)} members')
            store.load(paths, jobs=jobs)
            member_stats(tb, store, paths)
            src = {}
            pr = read_prior(prior, src) if prior else None
            if prior:
                log(f"prior {src['path']}, md5 {src['md5']}")
            out = run_qc(paths, [sx[p] for p in paths], [tb], store, read_bed(par_bed), min_samples=min_samples,
                         prior=pr, B=boot, seed=seed, n_pca=pca_probes, cf=capture_floor)
            out['res']['prior'] = src or None
            d = os.path.dirname(out_prefix)
            if d:
                os.makedirs(d, exist_ok=True)
            write_outputs(out_prefix, out)
            return dict(out['res']['status'])
    finally:
        VERBOSE = old_v


SCREEN_COLS = ['sample', 'verdict', 'reasons', 'noise_ratio', 'noise', 'ref_noise', 'missing_frac', 'n_used',
               'n_missing', 'n_capped', 'sex_inferred', 'sex_basis', 'x_l2', 'y_l2', 'pon_sex', 'input', 'one_share',
               'n_probes_file', 'pon', 'norm', 'note']


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--screen', metavar='SAMPLE.norm', help="a called sample's depth file: its call's own .norm "
                    '(PoN-divided) or the depth file of a run without --pon_file (PoN-free)')
    ap.add_argument('--pon', required=True, nargs='+', metavar='PON', help='binary PoN: with --screen the one the '
                    'call used; with --compare the new PoN(s), male and/or female')
    ap.add_argument('--compare', nargs='+', metavar='OLD', help='the previous binary PoN(s) to compare --pon with')
    ap.add_argument('--out', help='output prefix: PREFIX.tsv for --screen (optional); required with --compare')
    a = ap.parse_args()
    if bool(a.screen) == bool(a.compare):
        ap.error('give --screen or --compare')
    if a.screen and len(a.pon) != 1:
        ap.error('--screen takes one --pon: the PoN the call used')
    if a.compare and not a.out:
        ap.error('--compare needs --out')
    d = os.path.dirname(a.out) if a.out else ''
    try:
        if a.screen:
            with warnings.catch_warnings(), np.errstate(all='ignore'):
                warnings.simplefilter('ignore', RuntimeWarning)
                r = screen(a.screen, read_pon(a.pon[0]))
            if d:
                os.makedirs(d, exist_ok=True)
            if a.out:
                pd.DataFrame([r], columns=SCREEN_COLS).to_csv(a.out + '.tsv', sep='\t', index=False, float_format='%.5f')
            print(screen_line(r))
            return
        new_t = sorted((read_pon(p) for p in a.pon), key=lambda tb: tb['sex'] != 'M')
        old_t = sorted((read_pon(p) for p in a.compare), key=lambda tb: tb['sex'] != 'M')
    except (OSError, ValueError) as e:
        sys.exit(f'pon_qc.py: {e}')
    bysex = {tb['sex']: tb for tb in old_t}
    xy = [(tb, bysex[tb['sex']]) for tb in new_t if tb['sex'] in bysex]
    cmp = compare_tables(new_t[0], old_t[0], T, xy)
    if d:
        os.makedirs(d, exist_ok=True)
    write_compare(a.out, cmp)
    with open(a.out + '.json', 'w') as f:
        json.dump(jsonable({'compare': cmp[0]}), f, indent=1)
    txt = compare_text(cmp[0])
    with open(a.out + '.txt', 'w') as f:
        f.write(txt + '\n')
    print(txt)


if __name__ == '__main__':
    main()
