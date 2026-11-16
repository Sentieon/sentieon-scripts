#!/usr/bin/env python3
"""Build a CNVscope --blocks BED from a WES/panel target BED.

Groups capture targets into clusters (segmentation units) split at inter-target gaps > --gap, and
emits one line `chrom start end cn_neutral` per cluster. Feeding the result to

    cnvscope ... --blocks CLUSTERS.bed

makes each cluster its own block, so segmentation never crosses the capture gaps between clusters
while probes inside a cluster segment normally.

Sex chromosomes: --blocks and --sex are mutually exclusive in the binary, so the blocks BED is the
single ploidy authority. Pass --sex/--par here and the script folds sex-chromosome ploidy in:
  male:   chrX PAR = 2, chrX/chrY non-PAR = 1; chrY PAR dropped (reads multimap to chrX PAR)
  female: chrX = 2; chrY dropped
A ploidy change (e.g. the PAR boundary) breaks a cluster, since you never segment across it.

Dropping a region can't be expressed in a blocks BED (uncovered probes fall back to cn_neutral),
so with --sex the script also writes a cleaned target BED via --out-targets with the dropped
targets removed. Pass THAT to `cnvscope --target` so no probes are made there. Without it, dropped
regions (female chrY, male chrY-PAR) would be mis-called.

Only BED input is supported (0-based half-open, matching probe coordinates), plain or gzipped.

    python3 make_target_blocks.py targets.bed --gap 10000 -o clusters.bed
    python3 make_target_blocks.py targets.bed --gap 10000 --sex M --par par.hg38.bed \\
        -o clusters.bed --out-targets targets.clean.bed
"""
import argparse
import gzip
import sys
from collections import OrderedDict

CHRX = ('chrX', 'X')
CHRY = ('chrY', 'Y')


def open_maybe_gz(path):
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path, 'r')


def load_bed(path):
    """Read a BED into {chrom: [(start, end), ...]} preserving first-seen chrom order."""
    by_chrom = OrderedDict()
    with open_maybe_gz(path) as fh:
        for line in fh:
            if not line.strip() or line[0] == '#':
                continue
            if line.startswith(('track', 'browser')):
                continue
            f = line.split()
            if len(f) < 3:
                continue
            chrom, start, end = f[0], int(f[1]), int(f[2])
            if end <= start:
                continue
            by_chrom.setdefault(chrom, []).append((start, end))
    return by_chrom


def load_par(path):
    """Read a PAR BED into {chrom: [(start, end), ...]}."""
    par = {}
    with open_maybe_gz(path) as fh:
        for line in fh:
            if not line.strip() or line[0] == '#':
                continue
            f = line.split()
            if len(f) >= 3:
                par.setdefault(f[0], []).append((int(f[1]), int(f[2])))
    return par


def in_par(chrom, start, end, par):
    for s, e in par.get(chrom, ()):
        if start < e and end > s:
            return True
    return False


def parse_chrom_cn(spec):
    """Parse 'chr21=3,chr18=3' into {chrom: cn}."""
    out = {}
    for tok in spec.split(','):
        tok = tok.strip()
        if not tok:
            continue
        if '=' not in tok:
            raise ValueError("--chrom-cn entries must be CHROM=CN: %r" % tok)
        chrom, cn = tok.split('=', 1)
        out[chrom.strip()] = int(cn)
    return out


def target_cn(chrom, start, end, sex, par, base_cn, chrom_cn):
    """cn_neutral for a target, or None to drop it. --sex owns chrX/chrY when set."""
    if sex and (chrom in CHRX or chrom in CHRY):
        par_hit = in_par(chrom, start, end, par)
        if sex == 'M':
            if chrom in CHRY and par_hit:
                return None                 # chrY PAR: reads multimap to chrX PAR
            return 2 if par_hit else 1      # PAR diploid, else haploid
        if chrom in CHRY:
            return None                     # female chrY dropped
        return 2                            # female chrX diploid
    return chrom_cn.get(chrom, base_cn)


def cluster(targets, gap):
    """Merge sorted (start, end, cn) targets into blocks, splitting at gaps > `gap` or cn changes."""
    cs, ce, ccn = targets[0]
    for s, e, cn in targets[1:]:
        if s - ce > gap or cn != ccn:
            yield cs, ce, ccn
            cs, ce, ccn = s, e, cn
        else:
            ce = max(ce, e)
    yield cs, ce, ccn


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('target', help='capture target BED (0-based, plain or .gz)')
    ap.add_argument('-o', '--output', default='-', help='output blocks BED (default stdout)')
    ap.add_argument('--gap', type=int, default=10000,
                    help='split clusters at inter-target gaps larger than this (bp, default 10000)')
    ap.add_argument('--cn', type=int, default=2,
                    help='default cn_neutral for autosomes (default 2)')
    ap.add_argument('--padding', type=int, default=250,
                    help='pad cluster spans each side; match --target_padding (default 250)')
    ap.add_argument('--sex', choices=['M', 'F'],
                    help='fold human sex-chromosome ploidy into the blocks (needs --par for M)')
    ap.add_argument('--par', help='PAR BED (required for --sex M): chrX/chrY PAR = diploid')
    ap.add_argument('--out-targets',
                    help='write the target BED with dropped targets removed; pass to cnvscope --target')
    ap.add_argument('--chrom-cn', default='',
                    help='flat per-chrom cn_neutral overrides for non-sex chroms, e.g. chr21=3')
    args = ap.parse_args()

    if args.gap < 0:
        ap.error('--gap must be >= 0')
    if args.padding < 0:
        ap.error('--padding must be >= 0')
    if args.sex == 'M' and not args.par:
        ap.error('--sex M requires --par (to tell PAR from the haploid region)')
    chrom_cn = parse_chrom_cn(args.chrom_cn)
    for chrom, cn in list(chrom_cn.items()) + [('<default>', args.cn)]:
        if cn < 1 or cn > 4:
            ap.error('cn_neutral for %s must be in [1,4] (got %d)' % (chrom, cn))
    par = load_par(args.par) if args.par else {}

    by_chrom = load_bed(args.target)

    out = sys.stdout if args.output == '-' else open(args.output, 'w')
    tgt_out = open(args.out_targets, 'w') if args.out_targets else None
    n_targets = n_dropped = n_clusters = 0
    out.write('# CNVscope block units from %s (gap=%d padding=%d sex=%s)\n'
              % (args.target, args.gap, args.padding, args.sex or '-'))
    out.write('# chrom\tstart\tend\tcn_neutral\n')
    for chrom, intervals in by_chrom.items():
        # assign cn per target; dropped (cn=None) targets are excluded from both outputs
        targets = []
        for s, e in sorted(intervals):
            n_targets += 1
            cn = target_cn(chrom, s, e, args.sex, par, args.cn, chrom_cn)
            if cn is None:
                n_dropped += 1
                continue
            targets.append((s, e, cn))
            if tgt_out:
                tgt_out.write('%s\t%d\t%d\n' % (chrom, s, e))
        if not targets:
            continue
        prev_end = None
        for cs, ce, cn in cluster(targets, args.gap):
            start = max(0, cs - args.padding)
            end = ce + args.padding
            if prev_end is not None and start < prev_end:
                start = prev_end   # keep padded spans non-overlapping (tiny --gap or a ploidy split)
            out.write('%s\t%d\t%d\t%d\n' % (chrom, start, end, cn))
            prev_end = end
            n_clusters += 1
    if out is not sys.stdout:
        out.close()
    if tgt_out:
        tgt_out.close()

    sys.stderr.write('make_target_blocks: %d targets (%d dropped) -> %d clusters across %d chromosomes\n'
                     % (n_targets, n_dropped, n_clusters, len(by_chrom)))
    if n_dropped and not args.out_targets:
        sys.stderr.write('WARNING: %d targets dropped but --out-targets not set; pass a cleaned '
                         'target BED to cnvscope --target or those regions will be mis-called\n' % n_dropped)


if __name__ == '__main__':
    main()
