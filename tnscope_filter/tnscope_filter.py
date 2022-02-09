#!/usr/bin/env python
from __future__ import print_function
import argparse
import heapq
import itertools
import sys
import vcflib
import operator
import os
from vcflib.compat import *

class FilterRegistry(object):
    registry = []   # list of all filters

    @classmethod
    def add(cls, id, descr, params, func):
        params = params and params.split(',') or []
        cls.registry.append((id, descr, params, func))

    @classmethod
    def active(cls, args):
        return sorted([(id, descr, params, func)
            for id, descr, params, func in cls.registry if all(
                getattr(args, param, None) is not None for param in params)])

# decorator for each individual filter
def Filter(id, descr, params=''):
    def decorator(func):
        FilterRegistry.add(id, descr, params, func)
        return func
    return decorator

class TNscopeFilter(object):

    params = ( # name, defval or type, descr, metavar
        ('clear',                   'none', 'Existing filters to clear'),
        ('min_qual',                float,  'Minimum quality score'),
        ('min_depth',               int,    'Minimum depth'),
        ('min_qd',                  float,  'Minimum quality by depth'),
        ('min_read_pos_ranksum',    float,  'Minimum ReadPosRankSum'),
        ('min_base_qual_ranksum',   float,  'Minimum BaseQRankSum'),
        ('max_str_length',          int,    'Maximum STR length'),
        ('min_neighbor_base_qual',  float,  'Minimum NBQPS'),
        ('min_map_qual_ranksum',    float,  'Minimum MQRankSumPS'),
        ('min_depth_high_conf',     float,  'Minimum ALTHC'),
        ('max_pv',                  float,  'Maximum PV'),
        ('max_pv2',                 float,  'Maximum PV2'),
        ('max_str_pv',              float,  'Maximum PV in STR regions'),
        ('max_foxog',               float,  'Maximum fraction of alt reads indicating OxoG'),
        ('max_sor',                 float,  'Maximum Symmetric Odds Ratio'),
        ('max_ecnt',                int,    'Maximum number of events in this haplotype'),
        ('min_tumor_af',            float,  'Minimum tumor allele fraction'),
    )

    presets = {
        'amplicon': {
            'clear': 'triallelic_site',
            'max_pv': 0.1,
            'max_pv2': 0.1,
            'max_str_pv': 0.05,
            'max_foxog': 1,
            'max_sor': 3,
            'min_qual': 40,
            'max_str_length': 10,
            'max_ecnt': 10,
        },
        'ctdna': {
            'clear': 'triallelic_site',
            'min_qual': 9.225,
            'min_read_pos_ranksum': -2.78,
            'min_base_qual_ranksum': -5.85,
            'min_depth_high_conf':3.5,
        },
        'ctdna_umi': {
            'clear': 'triallelic_site',
            'min_neighbor_base_qual': 47.7, 
            'min_qual': 11.2,
            'min_read_pos_ranksum': -2.14,
            'min_base_qual_ranksum': -5.02,
            'min_depth_high_conf':4.5,
        },
        'tissue_panel': {
            'clear': 'triallelic_site',
            'max_pv': 0.1,
            'max_pv2': 0.1,
            'max_str_pv': 0.05,
            'max_foxog': 1,
            'max_sor': 3,
            'min_qual': 200,
            'max_str_length': 10,
            'max_ecnt': 10,
            'min_read_pos_ranksum': -8,
        }
    }

    @Filter('low_qual', 'Low quality score', 'min_qual')
    def low_qual(self, v):
        thr = self.args.min_qual
        val = v.qual
        return thr is not None and val is not None and val < thr

    @Filter('low_depth', 'Low coverage depth', 'min_depth')
    def low_depth(self, v):
        thr = self.args.min_depth
        val = v.samples[self.t_smid].get('AD')
        return thr is not None and isinstance(val, list) and sum(val) < thr
    
    @Filter('low_allele_frac', 'Low allele fraction in tumor', 'min_tumor_af')
    def low_depth(self, v):
        thr = self.args.min_tumor_af
        val = v.samples[self.t_smid].get('AF')
        return thr is not None and val is not None and val < thr

    @Filter('low_qual_by_depth', 'Low quality by depth', 'min_qd')
    def low_qual_by_depth(self, v):
        thr = self.args.min_qd
        val = v.samples[self.t_smid].get('AD')
        if not isinstance(val, list) or v.qual is None or thr is None:
            return False
        return v.qual < thr * sum(val)

    @Filter('read_pos_bias', 'Alt vs. Ref read position bias', 'min_read_pos_ranksum')
    def read_pos_bias(self, v):
        thr = self.args.min_read_pos_ranksum
        val = v.samples[self.t_smid].get('ReadPosRankSumPS')
        return thr is not None and val is not None and val < thr

    @Filter('base_qual_bias', 'Alt Vs. Ref base quality bias', 'min_base_qual_ranksum')
    def base_qual_bias(self, v):
        thr = self.args.min_base_qual_ranksum
        val = v.samples[self.t_smid].get('BaseQRankSumPS')
        return thr is not None and val is not None and val < thr

    @Filter('neighbor_base_qual', 'Low mean neighboring base quality', 'min_neighbor_base_qual')
    def neighbor_base_qual(self, v):
        thr = self.args.min_neighbor_base_qual
        val = v.samples[self.t_smid].get('NBQPS')
        return thr is not None and val is not None and val < thr

    @Filter('map_qual_bias', 'Alt Vs. Ref mapping qualities bias', 'min_map_qual_ranksum')
    def z_score_wilcoxon(self, v):
        thr = self.args.min_z_score_wilcoxon
        val = v.samples[self.t_smid].get('MQRankSumPS')
        return thr is not None and val is not None and val < thr

    @Filter('depth_high_conf', 'Depth high conf', 'min_depth_high_conf')
    def depth_high_conf(self, v):
        thr = self.args.min_depth_high_conf
        val = v.samples[self.t_smid].get('ALTHC')
        return thr is not None and val is not None and val < thr

    @Filter('short_tandem_repeat', 'Short tandem repeat', 'max_str_length')
    def short_tandem_repeat(self, v):
        thr = self.args.max_str_length
        val = v.info.get('RPA')
        return thr is not None and isinstance(val, list) and val[0] >= thr

    @Filter('insignificant', 'Insignificant call', 'max_pv,max_pv2')
    def insignificant(self, v):
        thr1 = self.args.max_pv
        thr2 = self.args.max_pv2
        val1 = v.info.get('PV')
        val2 = v.info.get('PV2')
        return thr1 is not None and thr2 is not None and val1 is not None and val1 > thr1 and val2 is not None and val2 > thr2

    @Filter('insignificant_str', 'Insignificant call in STR regions', 'max_str_pv')
    def insignificant_str(self, v):
        thr = self.args.max_str_pv
        val = v.info.get('PV')
        return v.info.get('STR') is not None and thr is not None and val is not None and val > thr

    @Filter('orientation_bias', 'Orientation bias', 'max_foxog')
    def depth_high_conf(self, v):
        thr = self.args.max_foxog
        val = v.samples[self.t_smid].get('FOXOG')
        return thr is not None and val is not None and val >= thr

    @Filter('strand_bias', 'Strand bias', 'max_sor')
    def depth_high_conf(self, v):
        thr = self.args.max_sor
        val = v.info.get('SOR')
        return thr is not None and val is not None and val > thr

    @Filter('noisy_region', 'Variant in noisy region', 'max_ecnt')
    def depth_high_conf(self, v):
        thr = self.args.max_ecnt
        val = v.info.get('ECNT')
        return thr is not None and val is not None and val > thr

    @classmethod
    def add_arguments(cls, parser, preset=None):
        vals = {}
        if preset is not None:
            vals = cls.presets.get(preset, vals)
        if cls.presets:
            v = cls.presets.keys()
            h = 'Parameter preset (choices: '+','.join(v)+')'
            parser.add_argument('-x', '--preset', choices=v, type=str, help=h)
        for k,v,h in cls.params:
            v = vals.get(k, v)
            if isinstance(v, type):
                h = argparse.SUPPRESS if vals else h
                parser.add_argument('--'+k, type=v, help=h)
            elif v is not None:
                h += ' (default: %(default)s)'
                parser.add_argument('--'+k, default=v, type=type(v), help=h)

    @staticmethod
    def grouper(*vcfs):
        q = []
        for k, vcf in enumerate(vcfs):
            i = iter(vcf)
            v = next(i, None)
            if v:
                heapq.heappush(q, (v.pos, v.end, k, v, i))
        while q:
            grp = [[] for _ in vcfs]
            pos, end, k, v, i = heapq.heappop(q)
            grp[k].append(v)
            v = next(i, None)
            if v:
                heapq.heappush(q, (v.pos, v.end, k, v, i))
            while q and q[0][0] == pos:
                pos, end, k, v, i = heapq.heappop(q)
                grp[k].append(v)
                v = next(i, None)
                if v:
                    heapq.heappush(q, (v.pos, v.end, k, v, i))
            yield (pos, grp)

    def __init__(self, args, t_smid, n_smid):
        self.args = args
        self.t_smid = t_smid
        self.n_smid = n_smid
        self.filters = FilterRegistry.active(args)
        if args.clear == 'all':
            self.clear = None
        elif args.clear == 'none' or args.clear == '':
            self.clear = set()
            self.clear.add('PASS')
        else:
            self.clear = set(args.clear.split(','))
            self.clear.add('PASS')

    def apply(self, invcf, outvcf):
        self.copy_header(invcf, outvcf)
        for chrom in invcf.contigs.keys():
            for pos, grp in self.grouper(invcf.range(chrom)):
                if len(grp[0]) > 1:
                    grp[0].sort(key=operator.attrgetter('qual'), reverse=True)
                for v in grp[0]:
                    outvcf.emit(self.apply_filters(v))

    def apply_filters(self, v):
        filters = set(id for id, _, _, func in self.filters if func(self, v))
        if self.clear is not None:
            filters = set(v.filter) - self.clear | filters
        v.filter = sorted(filters)
        flds = v.line.split('\t')
        flds[6] = v.filter and ';'.join(v.filter) or 'PASS'
        v.line = '\t'.join(flds)
        return v

    def copy_header(self, invcf, outvcf):
        fmt = '##FILTER=<ID=%s,Description="%s">'
        idg = itertools.groupby(self.filters, key=operator.itemgetter(0))
        toupdate = [fmt % next(g)[:2] for id, g in idg]
        if self.clear is None:
            clear = set(invcf.filters.keys())
        else:
            clear = set(invcf.filters.keys()) & self.clear
        clear.discard('PASS')
        toremove = [fmt % (id, '') for id in clear]
        outvcf.copy_header(invcf, toupdate, toremove)
        outvcf.emit_header()

def main(args):
    if not os.path.exists(args.vcf):
        print('Error: input file %s does not exist' % args.vcf, file=sys.stderr)
        return -1

    invcf = vcflib.VCF(args.vcf, 'r')
    try:
        t_smid = invcf.samples.index(args.tumor_sample)
    except ValueError:
        print('Error: tumor sample "%s" not in input file %s' %
            (args.tumor_sample, args.vcf), file=sys.stderr)
        return -1
    try:
        if args.normal_sample:
            n_smid = invcf.samples.index(args.normal_sample)
        else:
            n_smid = -1
    except ValueError:
        print('Error: normal sample "%s" not in input file %s' %
            (args.normal_sample, args.vcf), file=sys.stderr)
        return -1
    outvcf = vcflib.VCF(args.output, 'w')

    filter = TNscopeFilter(args, t_smid, n_smid)
    filter.apply(invcf, outvcf)
    outvcf.close()
    invcf.close()
    return 0

class MixedHelpFormatter(argparse.HelpFormatter):
    def _format_usage(self, usage, actions, groups, prefix):
        if prefix is None:
            prefix = 'usage: sentieon pyexec '
        return argparse.HelpFormatter._format_usage(
            self, usage, actions, groups, prefix)

    def _metavar_formatter(self, action, default):
        if action.metavar is None and action.type is not None:
            action.metavar = action.type.__name__.upper()
        return argparse.HelpFormatter._metavar_formatter(
            self, action, default)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument('-x', '--preset', const='?', nargs='?')
    args, argv = parser.parse_known_args()
    parser = argparse.ArgumentParser(formatter_class=MixedHelpFormatter)
    parser.add_argument('output', help='Output vcf file name')
    parser.add_argument('-v', '--vcf', required=True, help='Input vcf file name')
    parser.add_argument('--tumor_sample', required=True, help='Tumor sample name', type=str)
    parser.add_argument('--normal_sample', help='Normal sample name', type=str)
    TNscopeFilter.add_arguments(parser, args.preset)
    if args.preset == '?':
        argv = ['-x'] + argv
    elif args.preset is not None:
        argv = ['-x', args.preset] + argv
    sys.exit(main(parser.parse_args(argv)))

# vim: ts=4 sw=4 expandtab
