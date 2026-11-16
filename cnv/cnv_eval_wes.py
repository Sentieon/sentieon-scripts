"""CNV evaluation and labeling for capture assays: the WGS v2 rules (run7j/cnv_eval.py, 2026-09-10)
measured in probed base pairs. Without --target_bed the measure is plain length and this script
reproduces the WGS evaluator exactly. No depth field is read anywhere."""
import argparse
import sys
import bisect
import heapq
import io
import os
import copy
import numpy as np
import pandas as pd
from tabulate import tabulate
import vcflib
from vcflib import bgzf

pass_filter = ['PASS', 'cnvLength']


# ---------------------------------------------------------------------------
# IntervalList — reused from cnv_eval_vcf2.py
# ---------------------------------------------------------------------------
class IntervalList(object):
    def __init__(self, path) -> None:
        self.regions = self.load_bed(path)

    def load_bed(self, path) -> dict:
        regions = {}
        if path.endswith('.gz'):
            fp = bgzf.open(path, 'rb')
        else:
            fp = io.open(path, 'rb')
        for line in fp:
            if line.startswith(b'#') or line.startswith(b'track'):
                continue
            cols = line.rstrip().decode().split('\t')
            if len(cols) < 3:
                continue
            chrom, start, end = cols[:3]
            regions.setdefault(chrom, []).append([int(start), int(end)])
        fp.close()
        for chrom in regions:
            v = []
            s0 = e0 = None
            for (s, e) in sorted(regions[chrom]):
                if e0 is None:
                    s0 = s
                    e0 = e
                elif s > e0:
                    v.extend((s0, e0))
                    s0 = s
                    e0 = e
                else:
                    e0 = max(e0, e)
            if e0 is not None:
                v.extend((s0, e0))
            regions[chrom] = v
        return regions

    def bridge_gaps(self, max_gap):
        if max_gap <= 0:
            return
        for chrom in self.regions:
            r = self.regions[chrom]
            if len(r) < 4:
                continue
            v = [r[0], r[1]]
            for i in range(2, len(r), 2):
                if r[i] - v[-1] < max_gap:
                    v[-1] = r[i+1]
                else:
                    v.extend((r[i], r[i+1]))
            self.regions[chrom] = v

    def get(self, c, s, e) -> list:
        r = self.regions.get(c)
        v = []
        if r is not None:
            i = bisect.bisect_right(r, s)
            if i % 2:
                v.append((r[i-1], r[i]))
                i += 1
            while i < len(r) and r[i] < e:
                v.append((r[i], r[i+1]))
                i += 2
        return v

    def __contains__(self, cs) -> bool:
        if isinstance(cs, str):
            return cs in self.regions
        c, s = cs
        r = self.regions.get(c)
        return r and bisect.bisect_right(r, s) % 2 != 0

    def intersection(self, c, s, e):
        regs = self.get(c, s, e)
        itsc = []
        for rs, re in regs:
            rs = max(rs, s)
            re = min(re, e)
            itsc.append((rs, re))
        return itsc

    @staticmethod
    def parse_region(region):
        if isinstance(region, IntervalList):
            res = {}
            for chr, locs in region.regions.items():
                cur_regions = []
                for i in range(len(locs)//2):
                    cur_regions.append((locs[i*2], locs[i*2+1]))
                res[chr] = cur_regions
            return res
        if isinstance(region, dict):
            return region
        if isinstance(region, list):
            res = {}
            for c, s, e in region:
                res.setdefault(c, []).append([s, e])
            return res
        res = {}
        for rgn in region.split(','):
            chr, parts = rgn.split(':')
            start, end = [int(p) for p in parts.split('-')]
            res.setdefault(chr, []).append((start-1, end))
        return res

    @staticmethod
    def merge_regions(region_list):
        v = []
        s0 = e0 = None
        for (s, e) in sorted(region_list):
            if e0 is None:
                s0 = s
                e0 = e
            elif s > e0:
                v.extend((s0, e0))
                s0 = s
                e0 = e
            else:
                e0 = max(e0, e)
        if e0 is not None:
            v.extend((s0, e0))
        return v

    def intersect(self, region):
        if not region:
            return self
        result_interval = copy.deepcopy(self)
        regions = {}
        input_region = self.parse_region(region)
        for chr, rgns in input_region.items():
            if chr not in self.regions:
                continue
            res = []
            for r in rgns:
                rs = self.get(chr, r[0], r[1])
                if rs:
                    if rs[0][0] < r[0]:
                        rs[0] = (r[0], rs[0][1])
                    if rs[-1][1] > r[1]:
                        rs[-1] = (rs[-1][0], r[1])
                res += rs
            regions[chr] = self.merge_regions(res)
        result_interval.regions = regions
        return result_interval


# ---------------------------------------------------------------------------
# UnionFind
# ---------------------------------------------------------------------------
class UnionFind:
    def __init__(self, n):
        self.parent = list(range(n))
        self.rank = [0] * n

    def find(self, x):
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra == rb:
            return
        if self.rank[ra] < self.rank[rb]:
            ra, rb = rb, ra
        self.parent[rb] = ra
        if self.rank[ra] == self.rank[rb]:
            self.rank[ra] += 1


# ---------------------------------------------------------------------------
# Utility functions
# ---------------------------------------------------------------------------
def normalize_svlen(v):
    """Normalize SVLEN: handle list (take max abs) and negative values."""
    svlen = v.info.get('SVLEN')
    if svlen is None:
        return
    if isinstance(svlen, list):
        v.info['SVLEN'] = max(abs(x) for x in svlen)
    elif isinstance(svlen, int):
        if svlen < 0:
            v.info['SVLEN'] = -svlen


def clipped_overlap(s1, e1, s2, e2):
    return max(0, min(e1, e2) - max(s1, s2))


def round_half_up(x):
    return int(x + 0.5) if x >= 0 else -int(-x + 0.5)


# ---------------------------------------------------------------------------
# The probed measure — the one WES difference from WGS v2
# ---------------------------------------------------------------------------
target_bed = None   # capture targets (--target_bed); None -> the measure is plain length
call_min_size = None  # FP size gate for calls only; None -> use --min_cnv_size
strict_truth_gate = False  # by design a weak or sub-floor truth is never an FN but is credited when found; --strict_truth_gate fixes the denominator for arm comparisons
min_probed = 0      # bp floor on probed bases of an array for FN eligibility


def probed_len(chrom, s, e):
    """Base pairs of [s, e) the assay measures; plain length without --target_bed."""
    if e <= s:
        return 0
    if target_bed is None:
        return e - s
    return sum(re - rs for rs, re in target_bed.intersection(chrom, s, e))


def probed_overlap(chrom, s1, e1, s2, e2):
    """Probed base pairs shared by two spans."""
    return probed_len(chrom, max(s1, s2), min(e1, e2))


def probed_positive(chrom, s, e):
    """Does the assay measure anything in [s, e)? Always true without --target_bed."""
    return target_bed is None or probed_len(chrom, s, e) > 0


def truth_array(t):
    """The tandem array on the reference, where the depth change is expected."""
    return t.info['REFSTOP'] - t.info['REFSTART']


def truth_shift(t):
    """Expected CN shift averaged over the array: |CNDIFF * SVLEN| / array."""
    array = truth_array(t)
    return abs(t.info['CNDIFF'] * t.info['SVLEN']) / array if array > 0 else 0.0


def truth_probed_shift(t, chrom):
    """Largest expected shift a single capture target sees from this truth:
    |CNDIFF*SVLEN|/|A| * |P n A| / |P| maximised over targets P."""
    if target_bed is None:
        return truth_shift(t)
    ts, te = get_truth_region(t)
    best = 0.0
    for ps, pe in target_bed.get(chrom, ts, te):
        if pe > ps:
            best = max(best, clipped_overlap(ps, pe, ts, te) / (pe - ps))
    return truth_shift(t) * best


def truth_confident(t, min_dist):
    return t.info.get('DIST', 1.0) >= min_dist


def same_direction(call_cn, truth_cndiff, neutral_cn):
    return (call_cn - neutral_cn) * truth_cndiff > 0


def high_conf_overlap_frac(start, end, bed, chrom):
    if bed is None:
        return 1.0
    length = end - start
    if length <= 0:
        return 0.0
    regs = bed.get(chrom, start, end)
    total = sum(min(re, end) - max(rs, start) for rs, re in regs)
    return total / length


def event_inside_highconf(start, end, bed, chrom):
    if bed is None:
        return True
    return high_conf_overlap_frac(start, end, bed, chrom) >= 1.0


def event_any_highconf(start, end, bed, chrom):
    if bed is None:
        return True
    return high_conf_overlap_frac(start, end, bed, chrom) > 0


def probed_inside_highconf(start, end, bed, chrom):
    """All of the span's probed bases lie in the high-confidence bed."""
    if bed is None:
        return True
    if target_bed is None:
        return event_inside_highconf(start, end, bed, chrom)
    total = probed_len(chrom, start, end)
    if total <= 0:
        return False
    inside = sum(probed_len(chrom, rs, re) for rs, re in bed.intersection(chrom, start, end))
    return inside >= total


def get_truth_region(t):
    """The array in 0-based half-open coordinates; REFSTART/REFSTOP are 1-based half-open."""
    return t.info['REFSTART'] - 1, t.info['REFSTOP'] - 1


def get_call_region(c):
    return c.pos, c.end


cn_field = 'CN'
cn_is_diff = False

def get_cn(v, neutral_cn):
    val = None
    if cn_field in v.info:
        val = v.info[cn_field]
    elif cn_field in v.samples[0]:
        val = v.samples[0][cn_field]
    if val is None:
        return neutral_cn
    return neutral_cn + val if cn_is_diff else val


# ---------------------------------------------------------------------------
# Connected component grouping via sweep-line + union-find
# ---------------------------------------------------------------------------
def find_components(truths, calls, chrom):
    n_t = len(truths)
    n_c = len(calls)
    n = n_t + n_c
    if n == 0:
        return []

    # Build intervals: (start, end, global_index)
    intervals = []
    for i, t in enumerate(truths):
        s, e = get_truth_region(t)
        intervals.append((s, e, i))
    for j, c in enumerate(calls):
        s, e = get_call_region(c)
        intervals.append((s, e, n_t + j))

    intervals.sort(key=lambda x: x[0])
    uf = UnionFind(n)

    # Sweep line: maintain active intervals sorted by end
    active = []  # (end, global_index)
    for s, e, idx in intervals:
        # Remove expired
        while active and active[0][0] <= s:
            heapq.heappop(active)
        # Union with all remaining active intervals
        for _, aidx in active:
            uf.union(idx, aidx)
        heapq.heappush(active, (e, idx))

    # Group by component
    comps = {}
    for i in range(n):
        r = uf.find(i)
        comps.setdefault(r, ([], []))
        if i < n_t:
            comps[r][0].append(i)
        else:
            comps[r][1].append(i - n_t)
    return list(comps.values())


# ---------------------------------------------------------------------------
# Evaluation per component
# ---------------------------------------------------------------------------
def span_shift(truth_recs, cs, ce, chrom=None):
    """Expected copy-number shift over [cs, ce) from the truth records it overlaps: each spreads its
    CNDIFF * SVLEN base pairs evenly over its array and the span averages what it covers, counting
    only probed bases. The same number labels a segment (rounded) and judges a call (--tp_shift)."""
    length = probed_len(chrom, cs, ce)
    if length <= 0:
        return 0.0
    shift = 0.0
    for t in truth_recs:
        ts, te = get_truth_region(t)
        ovl = probed_overlap(chrom, ts, te, cs, ce)
        array = te - ts
        if ovl > 0 and array > 0:
            shift += t.info['CNDIFF'] * t.info['SVLEN'] * ovl / (array * length)
    return shift


def eval_component(truths, calls, t_indices, c_indices, neutral_cn, threshold,
                   fn_thresh, min_cnv_size, bed, chrom, tp_shift=0.4):
    comp_truths = [truths[i] for i in t_indices]
    comp_calls = [calls[j] for j in c_indices]

    t_evals = {}  # truth index -> 'TP' or 'FN'
    c_evals = {}  # call index -> 'TP' or 'FP'
    c_shift = {}  # call index -> expected shift over the call

    # Truth-side: same-direction calls must cover `threshold` of the array's probed bases
    for ti, t in zip(t_indices, comp_truths):
        ts, te = get_truth_region(t)
        array = te - ts
        if array <= 0 or t.info['SVLEN'] <= 0:
            t_evals[ti] = 'FN'
            continue
        parray = probed_len(chrom, ts, te)
        if parray <= 0:
            t_evals[ti] = 'FN'  # the assay sees nothing in this array
            continue
        agg_overlap = 0
        for c in comp_calls:
            c_cn = get_cn(c, neutral_cn)
            if not same_direction(c_cn, t.info['CNDIFF'], neutral_cn):
                continue
            cs, ce = get_call_region(c)
            agg_overlap += probed_overlap(chrom, ts, te, cs, ce)
        t_evals[ti] = 'TP' if agg_overlap / parray >= threshold else 'FN'

    # Call-side: the truth's expected shift over the call must reach tp_shift in the call's direction.
    # A call is judged on its own extent, whatever the truth's coverage status: a 5 kb call inside a
    # 100 kb deletion is right where it is even though the deletion is missed.
    for cj, c in zip(c_indices, comp_calls):
        c_cn = get_cn(c, neutral_cn)
        cs, ce = get_call_region(c)
        shift = span_shift(comp_truths, cs, ce, chrom)
        c_shift[cj] = shift
        c_evals[cj] = 'TP' if (c_cn - neutral_cn) * shift > 0 and abs(shift) >= tp_shift else 'FP'

    return t_evals, c_evals, c_shift


# ---------------------------------------------------------------------------
# Filtering gates
# ---------------------------------------------------------------------------
def passes_truth_filter(t, fn_thresh, min_cnv_size, min_shift, bed=None, chrom=None):
    if not truth_confident(t, fn_thresh):
        return False
    if truth_array(t) < min_cnv_size:
        return False
    # a change too small next to its array leaves no depth signal to find
    if truth_shift(t) < min_shift:
        return False
    ts, te = get_truth_region(t)
    if min_probed > 0 and probed_len(chrom, ts, te) < min_probed:
        return False
    # the assay must see the change on some single target, not just on average
    if target_bed is not None and truth_probed_shift(t, chrom) < min_shift:
        return False
    if bed:
        if not event_inside_highconf(ts, te, bed, chrom):
            return False
    return True


def truth_gated(t, fn_thresh, min_cnv_size, min_shift, bed=None, chrom=None):
    """--strict_truth_gate: may this truth be credited TP at all? Off, a truth the FN gate
    excludes is never charged an FN but is still credited when a call lands on it, so
    Truth.TOTAL depends on the call set. On, the same gate decides both. The call side is
    not changed either way."""
    return not strict_truth_gate or passes_truth_filter(
        t, fn_thresh, min_cnv_size, min_shift, bed, chrom)


def passes_call_filter(c, min_cnv_size, neutral_cn, bed=None, chrom=None):
    c_cn = get_cn(c, neutral_cn)
    if c_cn == neutral_cn:
        return False
    cs, ce = get_call_region(c)
    if c.info.get('SVLEN', ce - cs) < (min_cnv_size if call_min_size is None else call_min_size):
        return False
    if not probed_positive(chrom, cs, ce):
        return False  # off-target record: the assay never measured it
    if bed and not probed_inside_highconf(cs, ce, bed, chrom):
        return False
    return True


def passes_reward_filter(t, c, min_cnv_size, bed=None, chrom=None):
    cs, ce = get_call_region(c)
    t_size = truth_array(t)
    c_size = c.info.get('SVLEN', ce - cs)
    if t_size < min_cnv_size and c_size < min_cnv_size:
        return False
    if bed:
        ts, te = get_truth_region(t)
        t_any = event_any_highconf(ts, te, bed, chrom)
        c_any = event_any_highconf(cs, ce, bed, chrom)
        if not t_any and not c_any:
            return False
    return True


# ---------------------------------------------------------------------------
# cnv_type helper (Loss=0, Gain=1, Neutral=-1)
# ---------------------------------------------------------------------------
def cnv_type_val(cn, neutral_cn):
    diff = cn - neutral_cn
    if diff == 0:
        return -1
    return 0 if diff < 0 else 1


def cnv_type_truth(t, neutral_cn):
    diff = t.info['CNDIFF']
    if diff == 0:
        return -1
    return 0 if diff < 0 else 1


# ---------------------------------------------------------------------------
# Segment labeling
# ---------------------------------------------------------------------------
def segment_label(comp_truths, call, neutral_cn, args, chrom=None):
    """Expected CN of a call segment: each overlapping truth spreads its base-pair change
    (CNDIFF * SVLEN) evenly over its array, and the segment averages what it covers over its
    probed bases. Returns (status, label, expected); status is 'label', 'abstain' or 'none'."""
    cs, ce = get_call_region(call)
    seglen = probed_len(chrom, cs, ce)
    if seglen <= 0:
        return 'none', None, None  # nothing measured here
    touching = []
    for t in comp_truths:
        ts, te = get_truth_region(t)
        if probed_overlap(chrom, ts, te, cs, ce) <= 0 or te - ts <= 0:
            continue
        if not truth_confident(t, args.fn_thresh):
            return 'abstain', None, None  # one unreliable truth voids the whole segment
        touching.append(t)
    if not touching:
        return 'none', None, None
    expected = neutral_cn + span_shift(touching, cs, ce, chrom)
    label = round_half_up(expected)
    if args.label_abstain_band > 0 and abs(expected - label) > 0.5 - args.label_abstain_band:
        return 'abstain', None, expected  # the expectation sits on a rounding boundary
    return 'label', max(0, min(label, args.max_cn)), expected


# ---------------------------------------------------------------------------
# Main evaluation loop
# ---------------------------------------------------------------------------
def _in_par(pos, par_intervals_for_chrom):
    if not par_intervals_for_chrom:
        return False
    for s, e in par_intervals_for_chrom:
        if s <= pos < e:
            return True
    return False


def _par_filter_match(pos, par_intervals_for_chrom, par_filter):
    if par_filter is None:
        return True
    in_par = _in_par(pos, par_intervals_for_chrom)
    if par_filter == 'par':
        return in_par
    return not in_par  # 'nonpar'


def read_truth_variants(vcf_obj, chrom, do_filter, neutral_cn,
                         par_intervals_for_chrom=None, par_filter=None):
    variants = []
    for v in vcf_obj.range(chrom):
        if not _par_filter_match(v.pos, par_intervals_for_chrom, par_filter):
            continue
        if do_filter and v.info.get('CNDIFF', 0) == 0:
            continue
        if do_filter and pass_filter and v.filter and set(v.filter).difference(pass_filter):
            continue
        normalize_svlen(v)
        variants.append(v)
    return variants


def read_call_variants(vcf_obj, chrom, do_filter, neutral_cn,
                        par_intervals_for_chrom=None, par_filter=None):
    variants = []
    for v in vcf_obj.range(chrom):
        if not _par_filter_match(v.pos, par_intervals_for_chrom, par_filter):
            continue
        if cn_field in v.samples[0] and cn_field not in v.info:
            v.info[cn_field] = v.samples[0][cn_field]
        if do_filter:
            cn = get_cn(v, neutral_cn)
            if cn == neutral_cn:
                continue
            if pass_filter and v.filter and set(v.filter).difference(pass_filter):
                continue
        normalize_svlen(v)
        variants.append(v)
    return variants


def eval_all(truth_vcf, call_vcf, bed, args):
    default_neutral_cn = args.neutral_cn
    out_truth_stat = []
    out_call_stat = []

    out_call_writer = None
    out_truth_writer = None
    comp_hdr = '##INFO=<ID=COMP,Number=1,Type=String,Description="Connected component ID (chrom:start-end)">'
    pbp_hdr = '##INFO=<ID=PBP,Number=1,Type=Integer,Description="Probed base pairs (bases of the record inside --target_bed)">'
    pshift_hdr = '##INFO=<ID=PSHIFT,Number=1,Type=Float,Description="Largest expected CN shift on a single capture target from this truth; FN eligibility needs |PSHIFT| >= --min_shift">'
    if args.out_call:
        update = ['##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
                  '##INFO=<ID=EVAL,Number=1,Type=String,Description="Evaluation result">',
                  '##INFO=<ID=TRUTH,Number=1,Type=String,Description="Matched truth region">',
                  '##INFO=<ID=TRUTHSIZE,Number=1,Type=String,Description="Matched truth SVLEN">',
                  '##INFO=<ID=SHIFT,Number=1,Type=Float,Description="Expected copy-number shift over the call from the truth (CNDIFF*SVLEN spread over each array); TP needs |SHIFT| >= --tp_shift in the call direction">',
                  comp_hdr]
        if target_bed is not None:
            update.append(pbp_hdr)
        out_call_writer = vcflib.VCF(args.out_call, 'w')
        out_call_writer.copy_header(call_vcf, update=tuple(update))
        out_call_writer.emit_header()
    if args.out_truth:
        update = ['##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
                  '##INFO=<ID=EVAL,Number=1,Type=String,Description="Evaluation result">',
                  '##INFO=<ID=CALL,Number=1,Type=String,Description="Matched call region">',
                  '##INFO=<ID=CALLSIZE,Number=1,Type=String,Description="Matched call SVLEN">',
                  comp_hdr]
        if target_bed is not None:
            update.extend((pbp_hdr, pshift_hdr))
        out_truth_writer = vcflib.VCF(args.out_truth, 'w')
        out_truth_writer.copy_header(truth_vcf, update=tuple(update))
        out_truth_writer.emit_header()

    for chrom in truth_vcf.contigs:
        chrom_neutral_cn = get_chrom_neutral_cn(chrom, args.chrom_cn_map, default_neutral_cn)
        if chrom_neutral_cn < 0:
            continue  # skip (e.g. female chrY)
        par_intervals_for_chrom = args.par_intervals.get(chrom)

        chrom_t_queue = []
        chrom_c_queue = []
        qi1 = qi2 = 0

        for neutral_cn, par_filter in iter_eval_regions(
                chrom, chrom_neutral_cn, par_intervals_for_chrom, default_neutral_cn):
            truths = read_truth_variants(
                truth_vcf, chrom, True, neutral_cn,
                par_intervals_for_chrom, par_filter)
            calls = read_call_variants(
                call_vcf, chrom, True, neutral_cn,
                par_intervals_for_chrom, par_filter)
            components = find_components(truths, calls, chrom)

            # Track which indices get annotated (only those counted in stats)
            annotated_calls = {}   # call index -> (eval, best_truth_or_None, comp_id)
            annotated_truths = {}  # truth index -> (eval, matched_calls_or_None, comp_id)

            for t_indices, c_indices in components:
                if not t_indices and not c_indices:
                    continue

                # Compute component ID from genomic span
                starts, ends = [], []
                for ti in t_indices:
                    ts, te = get_truth_region(truths[ti])
                    starts.append(ts); ends.append(te)
                for cj in c_indices:
                    cs, ce = get_call_region(calls[cj])
                    starts.append(cs); ends.append(ce)
                comp_id = f'{chrom}:{min(starts)}-{max(ends)}'

                # Unmatched calls (no truth in component)
                if not t_indices:
                    for cj in c_indices:
                        c = calls[cj]
                        if passes_call_filter(c, args.min_cnv_size, neutral_cn, bed, chrom):
                            c_cn = get_cn(c, neutral_cn)
                            out_call_stat.append({
                                'call': c_cn,
                                'neutral_cn': neutral_cn,
                                'svlen_orig': c.info.get('SVLEN', c.end - c.pos),
                                'svlen': c.info.get('SVLEN', c.end - c.pos),
                                'truth_svlen': 0,
                                'eval': 'FP'
                            })
                            annotated_calls[cj] = ('FP', None, comp_id, 0.0)
                    continue

                # Unmatched truths (no calls in component)
                if not c_indices:
                    for ti in t_indices:
                        t = truths[ti]
                        if t.info['CNDIFF'] != 0 and passes_truth_filter(t, args.fn_thresh, args.min_cnv_size, args.min_shift, bed, chrom):
                            t_cn = t.info['CNDIFF'] + neutral_cn
                            out_truth_stat.append({
                                'call': t_cn,
                                'neutral_cn': neutral_cn,
                                'svlen': truth_array(t),
                                'eval': 'FN'
                            })
                            annotated_truths[ti] = ('FN', None, comp_id)
                    continue

                # Matched component: evaluate
                t_evals, c_evals, c_shift = eval_component(
                    truths, calls, t_indices, c_indices,
                    neutral_cn, args.min_overlap_thresh,
                    args.fn_thresh, args.min_cnv_size, bed, chrom, args.tp_shift)

                # Truth stats
                for ti in t_indices:
                    t = truths[ti]
                    ev = t_evals[ti]
                    t_cn = t.info['CNDIFF'] + neutral_cn
                    ts, te = get_truth_region(t)

                    if ev == 'TP':
                        # --strict_truth_gate: hold TPs to the call-independent FN gate too
                        if not truth_gated(t, args.fn_thresh, args.min_cnv_size, args.min_shift, bed, chrom):
                            continue
                        any_passes = False
                        for cj in c_indices:
                            c = calls[cj]
                            c_cn = get_cn(c, neutral_cn)
                            if same_direction(c_cn, t.info['CNDIFF'], neutral_cn):
                                if passes_reward_filter(t, c, args.min_cnv_size, bed, chrom):
                                    any_passes = True
                                    break
                        if not any_passes:
                            continue
                        out_truth_stat.append({
                            'call': t_cn,
                            'neutral_cn': neutral_cn,
                            'svlen': truth_array(t),
                            'eval': 'TP'
                        })
                        matched_calls = []
                        for cj in c_indices:
                            c = calls[cj]
                            c_cn = get_cn(c, neutral_cn)
                            cs, ce = get_call_region(c)
                            if same_direction(c_cn, t.info['CNDIFF'], neutral_cn) and clipped_overlap(ts, te, cs, ce) > 0:
                                matched_calls.append(c)
                        annotated_truths[ti] = ('TP', matched_calls, comp_id)
                    else:
                        if t.info['CNDIFF'] == 0:
                            continue
                        if not passes_truth_filter(t, args.fn_thresh, args.min_cnv_size, args.min_shift, bed, chrom):
                            continue
                        out_truth_stat.append({
                            'call': t_cn,
                            'neutral_cn': neutral_cn,
                            'svlen': truth_array(t),
                            'eval': 'FN'
                        })
                        annotated_truths[ti] = ('FN', None, comp_id)

                # Call stats — per-call eval, then event-level counting
                comp_call_tp = False
                comp_call_entries = []  # (cj, eval, best_truth_or_None)
                for cj in c_indices:
                    c = calls[cj]
                    ev = c_evals[cj]
                    c_cn = get_cn(c, neutral_cn)
                    if c_cn == neutral_cn:
                        continue
                    cs, ce = get_call_region(c)
                    if not probed_positive(chrom, cs, ce):
                        continue  # off-target record: the assay never measured it
                    if ev == 'TP':
                        # the same-direction truth contributing most of the call's shift
                        best_t = None
                        best_part = 0.0
                        for ti in t_indices:
                            t = truths[ti]
                            if not same_direction(c_cn, t.info['CNDIFF'], neutral_cn):
                                continue
                            part = abs(span_shift([t], cs, ce, chrom))
                            if part > best_part:
                                best_part = part
                                best_t = t
                        if best_t and passes_reward_filter(best_t, c, args.min_cnv_size, bed, chrom):
                            comp_call_tp = True
                            comp_call_entries.append((cj, 'TP', best_t))
                        else:
                            comp_call_entries.append((cj, 'FP', None))
                    else:
                        comp_call_entries.append((cj, 'FP', None))

                if not comp_call_entries:
                    continue
                event_eval = 'TP' if comp_call_tp else 'FP'
                rep_cn = get_cn(calls[comp_call_entries[0][0]], neutral_cn)
                rep_svlen = sum(calls[cj].info.get('SVLEN', calls[cj].end - calls[cj].pos) for cj, _, _ in comp_call_entries)

                # Skip FP if too small or any call is a border event
                if event_eval == 'FP':
                    if rep_svlen < (args.min_cnv_size if call_min_size is None else call_min_size):
                        continue
                    if bed:
                        all_inside = all(
                            probed_inside_highconf(*get_call_region(calls[cj]), bed, chrom)
                            for cj, _, _ in comp_call_entries)
                        if not all_inside:
                            continue

                truth_svlen = 0
                if comp_call_tp:
                    truth_svlen = max(truth_array(best_t) for _, ev, best_t in comp_call_entries if ev == 'TP')

                out_call_stat.append({
                    'call': rep_cn,
                    'neutral_cn': neutral_cn,
                    'svlen_orig': rep_svlen,
                    'svlen': rep_svlen,
                    'truth_svlen': truth_svlen,
                    'eval': event_eval
                })

                # Mark calls for annotation
                for cj, ev, best_t in comp_call_entries:
                    annotated_calls[cj] = (ev, best_t, comp_id, c_shift.get(cj, 0.0))

            # Emit region's records to chrom queues; only annotate those counted in stats
            for ti, t in enumerate(truths):
                if ti in annotated_truths:
                    ev, matched_calls, cid = annotated_truths[ti]
                    t.info['EVAL'] = ev
                    t.info['COMP'] = cid
                    if matched_calls:
                        t.info['CALL'] = ','.join(f'{mc.pos}-{mc.end}' for mc in matched_calls)
                        t.info['CALLSIZE'] = ','.join(str(mc.info.get('SVLEN', mc.end - mc.pos)) for mc in matched_calls)
                if target_bed is not None:
                    ts, te = get_truth_region(t)
                    t.info['PBP'] = probed_len(chrom, ts, te)
                    t.info['PSHIFT'] = f'{truth_probed_shift(t, chrom):.3f}'
                t.line = None
                heapq.heappush(chrom_t_queue, (t.pos, qi1, t))
                qi1 += 1
            for ci, c in enumerate(calls):
                if ci in annotated_calls:
                    ev, best_t, cid, shift = annotated_calls[ci]
                    c.info['EVAL'] = ev
                    c.info['COMP'] = cid
                    c.info['SHIFT'] = f'{shift:+.3f}'
                    if ev == 'TP' and best_t:
                        c.info['TRUTH'] = f'{best_t.pos}-{best_t.end}'
                        c.info['TRUTHSIZE'] = str(best_t.info['SVLEN'])
                if target_bed is not None:
                    c.info['PBP'] = probed_len(chrom, *get_call_region(c))
                c.line = None
                heapq.heappush(chrom_c_queue, (c.pos, qi2, c))
                qi2 += 1

        # Emit output VCFs in position order
        while out_truth_writer is not None and chrom_t_queue:
            _, _, v = heapq.heappop(chrom_t_queue)
            out_truth_writer.emit(v)
        while out_call_writer is not None and chrom_c_queue:
            _, _, v = heapq.heappop(chrom_c_queue)
            out_call_writer.emit(v)

    if out_call_writer:
        out_call_writer.close()
    if out_truth_writer:
        out_truth_writer.close()

    call_df = pd.DataFrame(out_call_stat) if out_call_stat else None
    truth_df = pd.DataFrame(out_truth_stat) if out_truth_stat else None
    if call_df is not None:
        call_df = call_df[~call_df['eval'].isnull()]
    if truth_df is not None:
        truth_df = truth_df[~truth_df['eval'].isnull()]
    return call_df, truth_df


# ---------------------------------------------------------------------------
# Labeling
# ---------------------------------------------------------------------------
def label_all(truth_vcf, call_vcf, out_call_writer, bed, args):
    default_neutral_cn = args.neutral_cn

    for chrom in truth_vcf.contigs:
        chrom_neutral_cn = get_chrom_neutral_cn(chrom, args.chrom_cn_map, default_neutral_cn)
        if chrom_neutral_cn < 0:
            continue
        par_intervals_for_chrom = args.par_intervals.get(chrom)

        q = []
        qi = 0
        for neutral_cn, par_filter in iter_eval_regions(
                chrom, chrom_neutral_cn, par_intervals_for_chrom, default_neutral_cn):
            truths = read_truth_variants(
                truth_vcf, chrom, True, neutral_cn,
                par_intervals_for_chrom, par_filter)
            calls = read_call_variants(
                call_vcf, chrom, False, neutral_cn,
                par_intervals_for_chrom, par_filter)
            components = find_components(truths, calls, chrom)

            # high-conf check is on the call side; a truth array may legitimately extend outside high-conf
            labeled = {}    # call index -> (label, expected CN)
            abstained = set()
            for t_indices, c_indices in components:
                if not c_indices or not t_indices:
                    continue
                comp_truths = [truths[i] for i in t_indices]
                for cj in c_indices:
                    status, label, expected = segment_label(
                        comp_truths, calls[cj], neutral_cn, args, chrom)
                    if status == 'label':
                        labeled[cj] = (label, expected)
                    elif status == 'abstain':
                        abstained.add(cj)

            # Push region's calls into per-chrom queue
            for ci, c in enumerate(calls):
                cs, ce = get_call_region(c)
                if ci in labeled:
                    label, expected = labeled[ci]
                    c.info['LABEL'] = label
                    c.info['LEXP'] = f'{expected:.3f}'
                elif ci in abstained:
                    pass  # a truth we do not trust, or a coin-flip rounding: say nothing
                elif probed_positive(chrom, cs, ce) and (
                        not bed or probed_inside_highconf(cs, ce, bed, chrom)):
                    c.info['LABEL'] = neutral_cn
                if target_bed is not None:
                    c.info['PBP'] = probed_len(chrom, cs, ce)
                c.line = None
                heapq.heappush(q, (c.pos, qi, c))
                qi += 1

        while q:
            _, _, v = heapq.heappop(q)
            out_call_writer.emit(v)


# ---------------------------------------------------------------------------
# Stats
# ---------------------------------------------------------------------------
def _compute_loss_gain(c_df, t_df):
    if t_df is not None and len(t_df) > 0:
        ncn = t_df['neutral_cn'] if 'neutral_cn' in t_df.columns else 2
        t_df_i = [t_df[t_df['call'] < ncn], t_df[t_df['call'] > ncn]]
        tp1 = [t_df_i[i][t_df_i[i]['eval'] == 'TP'].shape[0] for i in range(2)]
        fn1 = [t_df_i[i][t_df_i[i]['eval'] == 'FN'].shape[0] for i in range(2)]
    else:
        tp1 = [0, 0]
        fn1 = [0, 0]
    if c_df is not None and len(c_df) > 0:
        ncn = c_df['neutral_cn'] if 'neutral_cn' in c_df.columns else 2
        c_df_i = [c_df[c_df['call'] < ncn], c_df[c_df['call'] > ncn]]
        tp2 = [c_df_i[i][c_df_i[i]['eval'] == 'TP'].shape[0] for i in range(2)]
        fp1 = [c_df_i[i][c_df_i[i]['eval'] == 'FP'].shape[0] for i in range(2)]
    else:
        tp2 = [0, 0]
        fp1 = [0, 0]
    return tp1, fn1, tp2, fp1


def _stats_rows(tp1, fn1, tp2, fp1, types):
    sensitivity = [tp1[i]/(tp1[i] + fn1[i]) if tp1[i] + fn1[i] > 0 else 0 for i in range(2)]
    precision = [tp2[i]/(tp2[i] + fp1[i]) if tp2[i] + fp1[i] > 0 else 0 for i in range(2)]
    f1s = [2 * sensitivity[i] * precision[i] / (sensitivity[i] + precision[i]) if sensitivity[i] + precision[i] > 0 else 0 for i in range(2)]
    t_total = sum(tp1) + sum(fn1)
    c_total = sum(tp2) + sum(fp1)
    total_sensitivity = sum(tp1)/t_total if t_total else 0
    total_precision = sum(tp2)/c_total if c_total else 0
    f1 = 2 * total_sensitivity * total_precision / (total_sensitivity + total_precision) if total_sensitivity + total_precision > 0 else 0
    rows = [['All', t_total, sum(tp1), sum(fn1), c_total, sum(tp2), sum(fp1), total_sensitivity, total_precision, f1]]
    for i in range(2):
        rows.append([types[i], tp1[i] + fn1[i], tp1[i], fn1[i], tp2[i] + fp1[i], tp2[i], fp1[i], sensitivity[i], precision[i], f1s[i]])
    return rows


def get_stats(c_df, t_df):
    tp1, fn1, tp2, fp1 = _compute_loss_gain(c_df, t_df)
    results = _stats_rows(tp1, fn1, tp2, fp1, ["Loss", "Gain"])

    # Per-ploidy breakdown (only if multiple neutral_cn values present)
    has_ncn = lambda df: df is not None and len(df) > 0 and 'neutral_cn' in df.columns
    ploidy_vals = set()
    if has_ncn(t_df):
        ploidy_vals |= set(t_df['neutral_cn'].unique())
    if has_ncn(c_df):
        ploidy_vals |= set(c_df['neutral_cn'].unique())
    if len(ploidy_vals) > 1:
        for ncn_val in sorted(ploidy_vals):
            t_sub = t_df[t_df['neutral_cn'] == ncn_val] if has_ncn(t_df) else None
            c_sub = c_df[c_df['neutral_cn'] == ncn_val] if has_ncn(c_df) else None
            tp1s, fn1s, tp2s, fp1s = _compute_loss_gain(c_sub, t_sub)
            suffix = f"({int(ncn_val)})"
            rows = _stats_rows(tp1s, fn1s, tp2s, fp1s, [f"Loss{suffix}", f"Gain{suffix}"])
            results.extend(rows[1:])  # skip "All" for sub-groups

    cols = ['Type', 'Truth.TOTAL', 'Truth.TP', 'Truth.FN', 'CALL.TOTAL', 'CALL.TP', 'CALL.FP', 'Sensitivity', 'Precision', 'F1.Score']
    return pd.DataFrame(results, columns=cols).set_index('Type')


def get_stats_by_length(call_df, truth_df, cnv_lengths):
    max_sv = max(truth_df['svlen'].max() if len(truth_df) else 0,
                 call_df['svlen'].max() if len(call_df) else 0)
    if max_sv > cnv_lengths[-1]:
        cnv_lengths.append(max_sv)
    cnv_lengths = [0] + cnv_lengths
    all_results = []
    for i in range(1, len(cnv_lengths)):
        min_s, max_s = cnv_lengths[i-1], cnv_lengths[i]
        if call_df is not None and len(call_df) > 0:
            call_df_copy = call_df.copy()
            call_df_copy['eval_svlen'] = call_df_copy.apply(
                lambda r: r['truth_svlen'] if r.get('truth_svlen', 0) else r['svlen'], axis=1)
            c_df = call_df_copy[(call_df_copy['eval_svlen'] >= min_s) & (call_df_copy['eval_svlen'] <= max_s)]
        else:
            c_df = None
        if truth_df is not None and len(truth_df) > 0:
            t_df = truth_df[(truth_df['svlen'] >= min_s) & (truth_df['svlen'] <= max_s)]
        else:
            t_df = None
        results_df_s = get_stats(c_df, t_df)
        results_df_s['size_range_start'] = min_s
        results_df_s['size_range_stop'] = max_s
        all_results.append(results_df_s)
    return pd.concat(all_results)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def build_chrom_cn_neutral(sex, default_cn):
    chrom_cn = {}
    if not sex:
        return chrom_cn
    sex_chroms = {'chrX': True, 'X': True, 'chrY': True, 'Y': True}
    for name in sex_chroms:
        is_chrY = name in ('chrY', 'Y')
        if sex == 'M':
            chrom_cn[name] = 1
        else:
            chrom_cn[name] = -1 if is_chrY else default_cn
    return chrom_cn


def get_chrom_neutral_cn(chrom, chrom_cn_map, default_cn):
    return chrom_cn_map.get(chrom, default_cn)


def parse_par_bed(path):
    intervals = {}
    if not path:
        return intervals
    with open(path) as f:
        for line in f:
            if not line.strip() or line.startswith('#'):
                continue
            flds = line.rstrip().split('\t')
            if len(flds) < 3:
                continue
            chrom, start, end = flds[0], int(flds[1]), int(flds[2])
            intervals.setdefault(chrom, []).append((start, end))
    return intervals


def iter_eval_regions(chrom, chrom_neutral_cn, par_intervals_for_chrom, default_neutral_cn):
    if not par_intervals_for_chrom or chrom_neutral_cn == default_neutral_cn:
        yield (chrom_neutral_cn, None)
        return
    yield (default_neutral_cn, 'par')
    yield (chrom_neutral_cn, 'nonpar')


def main(args):
    global cn_field, cn_is_diff, pass_filter
    global target_bed, call_min_size, strict_truth_gate, min_probed
    cn_field = args.cn_field
    cn_is_diff = cn_field.endswith('DIFF')
    if args.all_records:
        pass_filter = None
    call_min_size = args.call_min_size
    strict_truth_gate = args.strict_truth_gate
    min_probed = args.min_probed
    if not args.label and cn_field == 'HMMCN' and not strict_truth_gate:
        print('warning: raw records without --strict_truth_gate: Truth.TOTAL depends on the call set; '
              'use --strict_truth_gate for arm comparisons', file=sys.stderr)
    args.chrom_cn_map = build_chrom_cn_neutral(args.sex, args.neutral_cn)
    args.par_intervals = parse_par_bed(args.par)
    call = vcflib.VCF(args.call, 'r')
    truth = vcflib.VCF(args.truth, 'r')
    bed = IntervalList(args.high_conf_bed) if args.high_conf_bed else None
    target_bed = IntervalList(args.target_bed) if args.target_bed else None

    if args.label:
        if not args.out_call:
            print("Error: --out_call required for label mode", file=sys.stderr)
            return 1
        update = ['##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
                  '##INFO=<ID=LABEL,Number=1,Type=Integer,Description="Label">',
                  '##INFO=<ID=LEXP,Number=1,Type=Float,Description="Expected copy number from truth, before rounding">']
        if target_bed is not None:
            update.append('##INFO=<ID=PBP,Number=1,Type=Integer,Description="Probed base pairs (bases of the record inside --target_bed)">')
        out_call = vcflib.VCF(args.out_call, 'w')
        out_call.copy_header(call, update=tuple(update))
        out_call.emit_header()
        label_all(truth, call, out_call, bed, args)
        out_call.close()
    else:
        call_df, truth_df = eval_all(truth, call, bed, args)
        results_df = get_stats(call_df, truth_df)
        if args.out:
            results_df.to_csv(args.out, sep='\t')
        print(tabulate(results_df, headers='keys', tablefmt='psql', floatfmt=('.4g')), flush=True)
        if args.cnvlengths and args.out_by_length:
            cnv_lengths = sorted([int(l) for l in args.cnvlengths.split(',')])
            # the floor starts the first bin; every edge above it is kept
            cnv_lengths = [args.min_cnv_size] + [l for l in cnv_lengths if l > args.min_cnv_size]
            if call_df is not None and truth_df is not None:
                all_results_df = get_stats_by_length(call_df, truth_df, cnv_lengths)
                all_results_df.to_csv(args.out_by_length, sep='\t')

        if args.stratifications:
            all_strats_df = None
            strat_basedir = os.path.dirname(args.stratifications)
            with open(args.stratifications) as f:
                for line in f:
                    if not line.strip():
                        continue
                    flds = line.strip().split('\t')
                    if len(flds) >= 2:
                        strat_name = flds[0]
                        strat_bed = IntervalList(os.path.join(strat_basedir, flds[1]))
                        if bed:
                            strat_bed = strat_bed.intersect(bed)
                        strat_call_df, strat_truth_df = eval_all(truth, call, strat_bed, args)
                        strat_results = get_stats(strat_call_df, strat_truth_df).reset_index()
                        strat_results.insert(0, 'Stratification', strat_name)
                        all_strats_df = strat_results if all_strats_df is None else pd.concat([all_strats_df, strat_results], ignore_index=True)
            if all_strats_df is not None and args.out:
                all_strats_df.to_csv(args.out + '.stratification', sep='\t', index=None)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='CNV evaluation and labeling with many-to-many matching, under a probed-base measure')
    parser.add_argument("truth", help="Truth VCF file")
    parser.add_argument("call", help="Call VCF file")
    parser.add_argument("out", help="Result output TSV file", nargs='?', default='')
    parser.add_argument('--high_conf_bed', help='High confidence region bed file')
    parser.add_argument('--target_bed', default=None, help='Capture-target BED. Every span is measured in the base pairs it shares with these targets: the label, the call shift, the truth coverage and the FN gate. Omit for WGS, where the measure is plain length.')
    parser.add_argument('--min_probed', default=0, type=int, help="Minimum probed base pairs in a truth array for FN eligibility (0 disables)")
    parser.add_argument('--call_min_size', default=None, type=int, help="Minimum SVLEN for a call to be charged as an FP; overrides --min_cnv_size on the call side only (0 charges every on-target call)")
    parser.add_argument('--strict_truth_gate', default=False, action='store_true', help="Credit a truth as TP only if it is also FN-eligible, so Truth.TOTAL does not depend on the call set; use for arm comparisons. Default: a weak, small or diluted truth is never an FN but counts when a call finds it; the call side is the same either way")
    parser.add_argument('--generous_truth_gate', default=False, action='store_true', help="The default behaviour; accepted for compatibility")
    parser.add_argument('--high_conf_max_gap', help="Bridge BED gaps smaller than this", default=200, type=int)
    parser.add_argument('--stratifications', help='Stratification file')
    parser.add_argument('--out_call', help='Output annotated call VCF')
    parser.add_argument('--out_truth', help='Output annotated truth VCF')
    parser.add_argument('--min_overlap_thresh', help="Fraction of a truth array's probed bases that same-direction calls must cover for the truth to count as found", default=0.3, type=float)
    parser.add_argument('--tp_shift', help="A call is TP when the truth's expected CN shift over its probed bases reaches this in the call's direction; 0.5 minus --label_abstain_band keeps eval and labels consistent", default=0.4, type=float)
    parser.add_argument('--label', help="Label mode: annotate calls with truth CN", default=False, action='store_true')
    parser.add_argument('--fn_thresh', '--min_dist', help="Minimum truth DIST to trust a record (FN eligibility and labeling)", default=0.95, type=float)
    parser.add_argument('--min_shift', help="Minimum |CNDIFF*SVLEN|/array for a truth to be FN-eligible; with --target_bed some single target must also see this much", default=0.5, type=float)
    parser.add_argument('--max_cn', help="Maximum copy number a label may take", default=10, type=int)
    parser.add_argument('--label_abstain_band', help="Skip the label when the expected CN is this close to a .5 boundary (0 disables)", default=0.1, type=float)
    parser.add_argument('--min_cnv_size', help="Minimum CNV size for penalty counting", default=5000, type=int)
    parser.add_argument('--cnvlengths', help="Comma-separated CNV length bins", default='10000,50000,100000')
    parser.add_argument('--out_by_length', help='Output stats by CNV length bins')
    parser.add_argument('--neutral_cn', default=2, type=int, help='Neutral copy number (default: 2)')
    parser.add_argument('--cn_field', default='CN', help='INFO/FORMAT field for copy number (default: CN, use HMMCN for raw CNVscope output)')
    parser.add_argument('--all_records', action='store_true', default=False, help='Ignore FILTER field (include all records)')
    parser.add_argument('--sex', default=None, choices=['M', 'F'], help='Sample sex for per-chromosome cn_neutral (M or F)')
    parser.add_argument('--par', default=None, help='PAR BED file: split chrX/chrY into PAR (cn_neutral=2) and non-PAR regions for sex-aware eval')

    args = parser.parse_args()
    sys.exit(main(args))
