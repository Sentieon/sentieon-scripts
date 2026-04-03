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
def clipped_overlap(s1, e1, s2, e2):
    return max(0, min(e1, e2) - max(s1, s2))


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


def get_truth_region(t):
    return t.info['REFSTART'], t.info['REFSTOP']


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
def eval_component(truths, calls, t_indices, c_indices, neutral_cn, threshold,
                   fn_thresh, min_cnv_size, bed, chrom):
    comp_truths = [truths[i] for i in t_indices]
    comp_calls = [calls[j] for j in c_indices]

    t_evals = {}  # truth index -> 'TP' or 'FN'
    c_evals = {}  # call index -> 'TP' or 'FP'

    # Truth-side: aggregate coverage from same-direction calls
    for ti, t in zip(t_indices, comp_truths):
        ts, te = get_truth_region(t)
        t_svlen = t.info['SVLEN']
        if t_svlen <= 0:
            t_evals[ti] = 'FN'
            continue
        agg_overlap = 0
        for c in comp_calls:
            c_cn = get_cn(c, neutral_cn)
            if not same_direction(c_cn, t.info['CNDIFF'], neutral_cn):
                continue
            cs, ce = get_call_region(c)
            agg_overlap += clipped_overlap(ts, te, cs, ce)
        coverage = agg_overlap / t_svlen
        t_evals[ti] = 'TP' if coverage >= threshold else 'FN'

    # Call-side: TP if it overlaps any TP truth in same direction
    for cj, c in zip(c_indices, comp_calls):
        c_cn = get_cn(c, neutral_cn)
        cs, ce = get_call_region(c)
        matched = False
        for ti, t in zip(t_indices, comp_truths):
            if t_evals.get(ti) != 'TP':
                continue
            if not same_direction(c_cn, t.info['CNDIFF'], neutral_cn):
                continue
            ts, te = get_truth_region(t)
            if clipped_overlap(ts, te, cs, ce) > 0:
                matched = True
                break
        c_evals[cj] = 'TP' if matched else 'FP'

    return t_evals, c_evals


# ---------------------------------------------------------------------------
# Filtering gates
# ---------------------------------------------------------------------------
def passes_truth_filter(t, fn_thresh, min_cnv_size, bed=None, chrom=None):
    event_size = min(t.info.get('PERIOD', t.info['SVLEN']), t.info['SVLEN'])
    if event_size < min_cnv_size:
        return False
    if t.info.get('DIST', 1.0) < fn_thresh:
        return False
    if bed:
        ts, te = get_truth_region(t)
        if not event_inside_highconf(ts, te, bed, chrom):
            return False
    return True


def passes_call_filter(c, min_cnv_size, neutral_cn, bed=None, chrom=None):
    c_cn = get_cn(c, neutral_cn)
    if c_cn == neutral_cn:
        return False
    cs, ce = get_call_region(c)
    if c.info.get('SVLEN', ce - cs) < min_cnv_size:
        return False
    if bed and not event_inside_highconf(cs, ce, bed, chrom):
        return False
    return True


def passes_reward_filter(t, c, min_cnv_size, bed=None, chrom=None):
    cs, ce = get_call_region(c)
    t_size = min(t.info.get('PERIOD', t.info['SVLEN']), t.info['SVLEN'])
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
# PERIOD-aware labeling
# ---------------------------------------------------------------------------
def effective_label(truth, call, neutral_cn, max_cn=10):
    ts, te = get_truth_region(truth)
    cs, ce = get_call_region(call)
    ovl = clipped_overlap(ts, te, cs, ce)
    if ovl <= 0:
        return neutral_cn
    period = truth.info.get('PERIOD', truth.info['SVLEN'])
    if period <= 0:
        period = truth.info['SVLEN']
    svlen = truth.info['SVLEN']
    if svlen <= 0:
        return neutral_cn
    n_periods_total = svlen / period
    n_periods_covered = max(1, round(ovl / period))
    cndiff = truth.info['CNDIFF']
    effective_cndiff = round(cndiff * n_periods_total / n_periods_covered)
    return max(0, min(neutral_cn + effective_cndiff, max_cn))


# ---------------------------------------------------------------------------
# Main evaluation loop
# ---------------------------------------------------------------------------
def read_truth_variants(vcf_obj, chrom, do_filter, neutral_cn):
    variants = []
    for v in vcf_obj.range(chrom):
        if do_filter and v.info.get('CNDIFF', 0) == 0:
            continue
        if do_filter and v.filter and set(v.filter).difference(pass_filter):
            continue
        variants.append(v)
    return variants


def read_call_variants(vcf_obj, chrom, do_filter, neutral_cn):
    variants = []
    for v in vcf_obj.range(chrom):
        if cn_field in v.samples[0] and cn_field not in v.info:
            v.info[cn_field] = v.samples[0][cn_field]
        if do_filter:
            cn = get_cn(v, neutral_cn)
            if cn == neutral_cn:
                continue
            if v.filter and set(v.filter).difference(pass_filter):
                continue
        variants.append(v)
    return variants


def eval_all(truth_vcf, call_vcf, bed, args):
    neutral_cn = args.neutral_cn
    out_truth_stat = []
    out_call_stat = []

    out_call_writer = None
    out_truth_writer = None
    comp_hdr = '##INFO=<ID=COMP,Number=1,Type=String,Description="Connected component ID (chrom:start-end)">'
    if args.out_call:
        update = ('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
                  '##INFO=<ID=EVAL,Number=1,Type=String,Description="Evaluation result">',
                  '##INFO=<ID=TRUTH,Number=1,Type=String,Description="Matched truth region">',
                  '##INFO=<ID=TRUTHSIZE,Number=1,Type=String,Description="Matched truth SVLEN">',
                  comp_hdr)
        out_call_writer = vcflib.VCF(args.out_call, 'w')
        out_call_writer.copy_header(call_vcf, update=update)
        out_call_writer.emit_header()
    if args.out_truth:
        update = ('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
                  '##INFO=<ID=EVAL,Number=1,Type=String,Description="Evaluation result">',
                  '##INFO=<ID=CALL,Number=1,Type=String,Description="Matched call region">',
                  '##INFO=<ID=CALLSIZE,Number=1,Type=String,Description="Matched call SVLEN">',
                  comp_hdr)
        out_truth_writer = vcflib.VCF(args.out_truth, 'w')
        out_truth_writer.copy_header(truth_vcf, update=update)
        out_truth_writer.emit_header()

    for chrom in truth_vcf.contigs:
        truths = read_truth_variants(truth_vcf, chrom, True, neutral_cn)
        calls = read_call_variants(call_vcf, chrom, True, neutral_cn)
        components = find_components(truths, calls, chrom)

        chrom_t_queue = []
        chrom_c_queue = []
        qi1 = qi2 = 0

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
                            'svlen_orig': c.info.get('SVLEN', c.end - c.pos),
                            'svlen': c.info.get('SVLEN', c.end - c.pos),
                            'truth_svlen': 0,
                            'eval': 'FP'
                        })
                        annotated_calls[cj] = ('FP', None, comp_id)
                continue

            # Unmatched truths (no calls in component)
            if not c_indices:
                for ti in t_indices:
                    t = truths[ti]
                    if t.info['CNDIFF'] != 0 and passes_truth_filter(t, args.fn_thresh, args.min_cnv_size, bed, chrom):
                        t_cn = t.info['CNDIFF'] + neutral_cn
                        out_truth_stat.append({
                            'call': t_cn,
                            'svlen': t.info['SVLEN'],
                            'eval': 'FN'
                        })
                        annotated_truths[ti] = ('FN', None, comp_id)
                continue

            # Matched component: evaluate
            t_evals, c_evals = eval_component(
                truths, calls, t_indices, c_indices,
                neutral_cn, args.min_overlap_thresh,
                args.fn_thresh, args.min_cnv_size, bed, chrom)

            # Truth stats
            for ti in t_indices:
                t = truths[ti]
                ev = t_evals[ti]
                t_cn = t.info['CNDIFF'] + neutral_cn
                ts, te = get_truth_region(t)

                if ev == 'TP':
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
                        'svlen': t.info['SVLEN'],
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
                    if not passes_truth_filter(t, args.fn_thresh, args.min_cnv_size, bed, chrom):
                        continue
                    out_truth_stat.append({
                        'call': t_cn,
                        'svlen': t.info['SVLEN'],
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
                if ev == 'TP':
                    best_t = None
                    best_ovl = 0
                    for ti in t_indices:
                        t = truths[ti]
                        if t_evals[ti] != 'TP':
                            continue
                        if not same_direction(c_cn, t.info['CNDIFF'], neutral_cn):
                            continue
                        ts, te = get_truth_region(t)
                        ovl = clipped_overlap(ts, te, cs, ce)
                        if ovl > best_ovl:
                            best_ovl = ovl
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
                if rep_svlen < args.min_cnv_size:
                    continue
                if bed:
                    all_inside = all(
                        event_inside_highconf(*get_call_region(calls[cj]), bed, chrom)
                        for cj, _, _ in comp_call_entries)
                    if not all_inside:
                        continue

            truth_svlen = 0
            if comp_call_tp:
                for ti in t_indices:
                    t = truths[ti]
                    if t_evals[ti] == 'TP':
                        truth_svlen = max(truth_svlen, t.info['SVLEN'])

            out_call_stat.append({
                'call': rep_cn,
                'svlen_orig': rep_svlen,
                'svlen': rep_svlen,
                'truth_svlen': truth_svlen,
                'eval': event_eval
            })

            # Mark calls for annotation
            for cj, ev, best_t in comp_call_entries:
                annotated_calls[cj] = (ev, best_t, comp_id)

        # Emit all records to VCFs; only annotate those counted in stats
        for ti, t in enumerate(truths):
            if ti in annotated_truths:
                ev, matched_calls, cid = annotated_truths[ti]
                t.info['EVAL'] = ev
                t.info['COMP'] = cid
                if matched_calls:
                    t.info['CALL'] = ','.join(f'{mc.pos}-{mc.end}' for mc in matched_calls)
                    t.info['CALLSIZE'] = ','.join(str(mc.info.get('SVLEN', mc.end - mc.pos)) for mc in matched_calls)
            t.line = None
            heapq.heappush(chrom_t_queue, (t.pos, qi1, t))
            qi1 += 1
        for ci, c in enumerate(calls):
            if ci in annotated_calls:
                ev, best_t, cid = annotated_calls[ci]
                c.info['EVAL'] = ev
                c.info['COMP'] = cid
                if ev == 'TP' and best_t:
                    c.info['TRUTH'] = f'{best_t.pos}-{best_t.end}'
                    c.info['TRUTHSIZE'] = str(best_t.info['SVLEN'])
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
    neutral_cn = args.neutral_cn

    for chrom in truth_vcf.contigs:
        truths = read_truth_variants(truth_vcf, chrom, True, neutral_cn)
        calls = read_call_variants(call_vcf, chrom, False, neutral_cn)
        components = find_components(truths, calls, chrom)

        # Determine labels: only from truths fully inside high-conf
        labeled = {}  # call index -> label
        for t_indices, c_indices in components:
            if not c_indices or not t_indices:
                continue
            # Filter to truths 100% inside high-conf
            comp_truths = []
            for i in t_indices:
                t = truths[i]
                ts, te = get_truth_region(t)
                if bed and not event_inside_highconf(ts, te, bed, chrom):
                    continue
                comp_truths.append(t)
            if not comp_truths:
                continue
            for cj in c_indices:
                c = calls[cj]
                cs, ce = get_call_region(c)
                best_label = None
                best_ovl = 0
                for t in comp_truths:
                    ts, te = get_truth_region(t)
                    ovl = clipped_overlap(ts, te, cs, ce)
                    if ovl <= 0:
                        continue
                    c_svlen = c.info.get('SVLEN', ce - cs)
                    if c_svlen > 0 and ovl / c_svlen < args.min_overlap_thresh:
                        continue
                    if ovl > best_ovl:
                        best_ovl = ovl
                        best_label = effective_label(t, c, neutral_cn)
                if best_label is not None:
                    labeled[cj] = best_label

        # Emit all calls: matched → truth-derived label, unmatched inside high-conf → neutral
        q = []
        for ci, c in enumerate(calls):
            if ci in labeled:
                c.info['LABEL'] = labeled[ci]
            elif not bed or event_inside_highconf(*get_call_region(c), bed, chrom):
                c.info['LABEL'] = neutral_cn
            c.line = None
            heapq.heappush(q, (c.pos, ci, c))
        while q:
            _, _, v = heapq.heappop(q)
            out_call_writer.emit(v)


# ---------------------------------------------------------------------------
# Stats
# ---------------------------------------------------------------------------
def get_stats(c_df, t_df):
    if t_df is not None and len(t_df) > 0:
        t_df_i = [t_df[t_df['call'] < 2], t_df[t_df['call'] > 2]]
        tp1 = [t_df_i[i][t_df_i[i]['eval'] == 'TP'].shape[0] for i in range(2)]
        fn1 = [t_df_i[i][t_df_i[i]['eval'] == 'FN'].shape[0] for i in range(2)]
        sensitivity = [tp1[i]/(tp1[i] + fn1[i]) if tp1[i] + fn1[i] > 0 else 0 for i in range(2)]
    else:
        tp1 = [0, 0]
        fn1 = [0, 0]
        sensitivity = [1, 1]
    if c_df is not None and len(c_df) > 0:
        c_df_i = [c_df[c_df['call'] < 2], c_df[c_df['call'] > 2]]
        tp2 = [c_df_i[i][c_df_i[i]['eval'] == 'TP'].shape[0] for i in range(2)]
        fp1 = [c_df_i[i][c_df_i[i]['eval'] == 'FP'].shape[0] for i in range(2)]
        precision = [tp2[i]/(tp2[i] + fp1[i]) if tp2[i] + fp1[i] > 0 else 0 for i in range(2)]
    else:
        tp2 = [0, 0]
        fp1 = [0, 0]
        precision = [1, 1]
    f1s = [2 * sensitivity[i] * precision[i] / (sensitivity[i] + precision[i]) if sensitivity[i] + precision[i] > 0 else 0 for i in range(2)]
    t_total = sum(tp1) + sum(fn1)
    c_total = sum(tp2) + sum(fp1)
    total_sensitivity = sum(tp1)/t_total if t_total else 0
    total_precision = sum(tp2)/c_total if c_total else 0
    f1 = 2 * total_sensitivity * total_precision / (total_sensitivity + total_precision) if total_sensitivity + total_precision > 0 else 0
    results = [['All', t_total, sum(tp1), sum(fn1), c_total, sum(tp2), sum(fp1), total_sensitivity, total_precision, f1]]
    types = ["Loss", "Gain"]
    for i in range(2):
        results.append([types[i], tp1[i] + fn1[i], tp1[i], fn1[i], tp2[i] + fp1[i], tp2[i], fp1[i], sensitivity[i], precision[i], f1s[i]])
    results_df = pd.DataFrame(results, columns=['Type', 'Truth.TOTAL', 'Truth.TP', 'Truth.FN', 'CALL.TOTAL', 'CALL.TP', 'CALL.FP', 'Sensitivity', 'Precision', 'F1.Score']).set_index('Type')
    return results_df


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
def main(args):
    global cn_field, cn_is_diff
    cn_field = args.cn_field
    cn_is_diff = cn_field.endswith('DIFF')
    call = vcflib.VCF(args.call, 'r')
    truth = vcflib.VCF(args.truth, 'r')
    bed = IntervalList(args.high_conf_bed) if args.high_conf_bed else None

    if args.label:
        if not args.out_call:
            print("Error: --out_call required for label mode", file=sys.stderr)
            return 1
        update = ('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
                  '##INFO=<ID=LABEL,Number=1,Type=Integer,Description="Label">')
        out_call = vcflib.VCF(args.out_call, 'w')
        out_call.copy_header(call, update=update)
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
            if args.min_cnv_size <= cnv_lengths[0]:
                cnv_lengths[0] = args.min_cnv_size
            else:
                cnv_lengths = [args.min_cnv_size] + cnv_lengths
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
    parser = argparse.ArgumentParser(description='CNV evaluation and labeling with many-to-many matching')
    parser.add_argument("truth", help="Truth VCF file")
    parser.add_argument("call", help="Call VCF file")
    parser.add_argument("out", help="Result output TSV file", nargs='?', default='')
    parser.add_argument('--high_conf_bed', help='High confidence region bed file')
    parser.add_argument('--stratifications', help='Stratification file')
    parser.add_argument('--out_call', help='Output annotated call VCF')
    parser.add_argument('--out_truth', help='Output annotated truth VCF')
    parser.add_argument('--min_overlap_thresh', help="Minimum overlap threshold", default=0.3, type=float)
    parser.add_argument('--label', help="Label mode: annotate calls with truth CN", default=False, action='store_true')
    parser.add_argument('--fn_thresh', help="Minimum DIST to count truth as FN", default=0.95, type=float)
    parser.add_argument('--min_cnv_size', help="Minimum CNV size for penalty counting", default=5000, type=int)
    parser.add_argument('--cnvlengths', help="Comma-separated CNV length bins", default='10000,50000,100000')
    parser.add_argument('--out_by_length', help='Output stats by CNV length bins')
    parser.add_argument('--neutral_cn', default=2, type=int, help='Neutral copy number (default: 2)')
    parser.add_argument('--cn_field', default='CN', help='INFO/FORMAT field for copy number (default: CN, use HMMCN for raw CNVscope output)')

    args = parser.parse_args()
    sys.exit(main(args))
