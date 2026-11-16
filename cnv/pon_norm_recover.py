#!/usr/bin/env python3
"""pon_norm_recover.py -- PoN-free depth from the .norm of a CNVscope run with --pon_file (--pon_mode ratio).

usage: python3 -I pon_norm_recover.py --norm RUN.norm --pon USED.pon.tsv [--recenter X | --log RUN.log]
           --out RECOVERED.norm [--report REPORT.json]

CNVscope's ratio mode writes DP = min(DP_raw / 2^median_log2dp, 8) / m and DP0 = min(DP0_raw / 2^median_log2dp, 8) / m.
Both are divided by the PoN's DP median median_log2dp (not median_log2dp0; tested on the study runs), as the PoN
references each probe: for a --ratio_sex PoN, chrX/chrY take the same-sex medians, and male non-PAR chrX and chrY
carry the +1 shift. m is the run's re-center median.

This tool writes the same .norm format and columns with DP = DP_run * 2^median_log2dp * m, and DP0 likewise.
- Values below 1e-4 are written as 0: CNVscope writes a floor value for zero depth, and real raw depths are >= 1.6e-3.
- The other columns pass through unchanged.
- Rows that cannot be recovered are dropped: DP written as exactly 1 (masked: PoN-excluded, or panel median DP below
  CNVscope's pon_min_median); DP or DP0 at the ratio cap 8/m (saturated); rows the PoN does not reference (excluded,
  or chrX/chrY dark for that sex); rows absent from the PoN.

m:
- --recenter X undoes it exactly, to the printed precision. X is the value in the CNVscope log line
  "PoN applied: ... re-center median=X"; --log parses that line.
- Otherwise m is estimated so that the recovered depth's median log2 over the recovered autosomal probes equals the
  median of the PoN's medians there (the PoN's typical member scale). Over the 54 study samples this errs by a median
  of +0.0013, rms 0.017 and at most 0.038 log2; the log value errs by <= 7e-5.
Use the PoN the run used (its --pon_file), as its .pon.tsv. The report says which way m was set.
"""
import sys

sys.dont_write_bytecode = True                   # no __pycache__ next to the package

import argparse  # noqa: E402
import json  # noqa: E402
import os  # noqa: E402
import re  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

ZERO_FLOOR = 1e-4                                # recovered depth below this is CNVscope's zero floor
AUTO = {f'chr{i}' for i in range(1, 23)} | {str(i) for i in range(1, 23)}


def read_recenter(log_path):
    """The "re-center median=X" value of a CNVscope (or segment.sh) log, or None."""
    with open(log_path) as f:
        m = re.findall(r're-center median=([0-9.eE+-]+)', f.read())
    return float(m[-1]) if m else None


def read_pon(path):
    cols = ['contig', 'start', 'end', 'median_log2dp', 'median_log2dp0', 'probe_idx', 'excluded', 'panel_median_dp']
    nskip = 0
    with open(path) as f:
        for ln in f:
            if not ln.startswith('#'):
                break
            nskip += 1
    p = pd.read_csv(path, sep='\t', skiprows=nskip, usecols=cols, dtype={'contig': str})
    return p


def recover(norm_path, pon_path, recenter=None):
    """Recovered rows of a PoN-divided .norm. Returns (header line, DataFrame of kept rows as text, info dict)."""
    with open(norm_path) as f:
        header = f.readline()
    d = pd.read_csv(norm_path, sep='\t', header=None, skiprows=1, dtype=str, keep_default_na=False)
    n_in = len(d)
    p = read_pon(pon_path)
    key = lambda c, s, e: c.astype(str) + ':' + s.astype(str) + '-' + e.astype(str)
    p.index = key(p['contig'], p['start'], p['end'])
    dk = key(d[0], d[1], d[2])
    j = p.reindex(dk.values)
    absent = j['median_log2dp'].isna().values
    ref = ~absent & (j['probe_idx'].values >= 0)
    dp = pd.to_numeric(d[3], errors='coerce').values
    dp0 = pd.to_numeric(d[6], errors='coerce').values
    masked = dp == 1.0
    keep = ref & ~masked & np.isfinite(dp)
    auto = d[0].isin(AUTO).values
    med = j['median_log2dp'].values.astype(float)       # CNVscope divides DP and DP0 by the DP median
    # m: from the log value, or the PoN's typical member scale over the recovered autosomal probes
    a = keep & auto & (dp > 0)
    if recenter is not None:
        m, how = float(recenter), 'recenter'
    else:
        lr = np.log2(dp[a]) + med[a]
        m, how = float(2.0 ** (np.median(med[a]) - np.median(lr))), 'estimated'
    # the ratio cap: 8/m in the run file; when m is estimated, the file's maximum if it lies near 8/m (m errs <= 2.7%;
    # at worst this drops the one largest genuine value)
    top = float(np.nanmax(np.where(keep, dp, np.nan))) if keep.any() else np.nan
    cap = 8.0 / m if how == 'recenter' else (top if top * m > 7.7 else np.inf)
    capped = keep & ((dp >= cap * (1 - 1e-4)) | (dp0 >= cap * (1 - 1e-4)))
    keep &= ~capped
    with np.errstate(over='ignore', invalid='ignore'):
        rec = dp * np.exp2(med) * m
        rec0 = dp0 * np.exp2(med) * m
    rec = np.where(rec < ZERO_FLOOR, 0.0, rec)
    rec0 = np.where(np.isfinite(rec0) & (rec0 >= ZERO_FLOOR), rec0, 0.0)
    out = d[keep].copy()
    out[3] = [f'{v:.6g}' for v in rec[keep]]
    out[6] = [f'{v:.6g}' for v in rec0[keep]]
    pmd = j['panel_median_dp'].values.astype(float)
    info = {'norm': norm_path, 'pon': pon_path, 'rows_in': int(n_in), 'rows_out': int(keep.sum()), 'recenter': how,
            'm': m, 'log2_m': float(np.log2(m)),
            'dropped_masked': int((ref & masked).sum()),
            'dropped_masked_unexplained': int((ref & masked & (pmd >= 0.1)).sum()), 'dropped_capped': int(capped.sum()),
            'dropped_not_referenced': int((~absent & ~ref).sum()), 'dropped_absent': int(absent.sum()),
            'zero_depth_rows': int((keep & (rec == 0)).sum()),
            'recovered_auto': int((keep & auto).sum()), 'recovered_chrX': int((keep & d[0].isin(['chrX', 'X']).values).sum()),
            'recovered_chrY': int((keep & d[0].isin(['chrY', 'Y']).values).sum())}
    return header, out, info


def write_norm(path, header, rows):
    d = os.path.dirname(path)
    if d:
        os.makedirs(d, exist_ok=True)
    with open(path, 'w') as f:
        f.write(header)
        rows.to_csv(f, sep='\t', header=False, index=False)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--norm', required=True, help='the .norm of a CNVscope run with --pon_file (--pon_mode ratio)')
    ap.add_argument('--pon', required=True, help='the .pon.tsv of the PoN that run used')
    g = ap.add_mutually_exclusive_group()
    g.add_argument('--recenter', type=float, default=None, help='the run\'s "re-center median=" value')
    g.add_argument('--log', default=None, help='the run\'s CNVscope or segment.sh log (re-center median parsed)')
    ap.add_argument('--out', required=True, help='recovered .norm')
    ap.add_argument('--report', default=None, help='optional JSON report')
    a = ap.parse_args()
    rc = a.recenter
    if a.log:
        rc = read_recenter(a.log)
        if rc is None:
            sys.exit(f'no "re-center median=" in {a.log}')
    header, rows, info = recover(a.norm, a.pon, rc)
    write_norm(a.out, header, rows)
    if a.report:
        with open(a.report, 'w') as f:
            json.dump(info, f, indent=1)
    print(f"recovered {info['rows_out']}/{info['rows_in']} rows (dropped: {info['dropped_masked']} masked, "
          f"{info['dropped_capped']} at the ratio cap, {info['dropped_not_referenced']} not referenced by the PoN, "
          f"{info['dropped_absent']} absent from it); "
          f"m = {info['m']:.6g} ({'from --recenter/--log' if info['recenter'] == 'recenter' else 'estimated from the PoN scale'})")


if __name__ == '__main__':
    main()
