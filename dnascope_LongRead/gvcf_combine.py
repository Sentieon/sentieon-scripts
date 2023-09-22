from __future__ import print_function
import heapq
import sys
import argparse
import vcflib
import copy
import uuid
import os
import time
import resource

from vcflib.compat import *

class Combiner(vcflib.Shardable, vcflib.ShardResult):

    def __init__(self, gvcf_input, vcf_input, gvcf_output):
        self.vcf_in = vcflib.VCF(vcf_input, 'r')
        self.gvcf_in = vcflib.VCF(gvcf_input, 'r')
        self.contigs = self.gvcf_in.contigs
        extras = self.extra_headers(self.gvcf_in, self.vcf_in)
        self.vcftmp = None
        self.stdout = False
        if gvcf_output == '-':
            dir = os.getenv('SENTIEON_TMPDIR') 
            if not dir:
                dir = '/tmp'
            gvcf_output = os.path.join(dir, str(uuid.uuid4()) + ".vcf")
            self.stdout = True
        self.gvcf_out = vcflib.VCF(gvcf_output, 'w')
        self.gvcf_out.copy_header(self.gvcf_in, extras)
        self.gvcf_out.emit_header()
        if self.stdout:
            chr_line = ""
            for line in self.gvcf_out.headers:
                if line.startswith('#CHROM'):
                    chr_line = line
                    continue
                sys.stdout.write(line + '\n')
            sys.stdout.write(chr_line + '\n')
    
    @staticmethod
    def extra_headers(vcf1, vcf2):
        filter_diff = ['##FILTER=<ID=%s,' %f for f in vcf2.filters.keys() if f not in vcf1.filters]
        info_diff = ['##INFO=<ID=%s,' %f for f in vcf2.infos.keys() if f not in vcf1.infos]
        fmt_diff = ['##FORMAT=<ID=%s,' %f for f in vcf2.formats.keys() if f not in vcf1.formats]
        diffs = filter_diff + info_diff + fmt_diff
        return [h for h in vcf2.headers if any([h.startswith(s) for s in diffs])]

    def __del__(self):
        self.vcf_in.close()
        self.gvcf_in.close()
        self.gvcf_out.close()
        if self.stdout and not hasattr(self, 'shard'):
            os.unlink(self.gvcf_out.path)
            os.unlink(self.gvcf_out.path + '.idx')

    @staticmethod
    def grouper(gvcfi, vcfi):
        q = []
        iters = [iter(i) for i in (gvcfi, vcfi)]
        for k, i in enumerate(iters):
            v = next(i, None)
            if v:
                heapq.heappush(q, (v.pos, v.end, k, v, i))
        while q:
            grp = [None, None]
            pos, end, k, v, i = heapq.heappop(q)
            grp[k] = v
            vv = next(i, None)
            if vv:
                heapq.heappush(q, (vv.pos, vv.end, k, vv, i))
            if q and q[0][2] != k:
                pos2, end2, k2, v2, i2 = heapq.heappop(q)
                if k < k2:
                    g, v = v, v2
                else:
                    g = v2
                if len(g.alt) == 1:
                    if Combiner.ovl(g, v):
                        for gg in Combiner.split_g(g, v):
                            heapq.heappush(q, (gg.pos, gg.end, 0, gg, iters[0]))
                        heapq.heappush(q, (v.pos, v.end, 1, v, iters[1]))
                        continue
                    heapq.heappush(q, (pos2, end2, k2, v2, i2))
                elif pos2 == pos:
                    grp[k2] = v2
                    nv = next(i2, None)
                    if nv:
                        heapq.heappush(q, (nv.pos, nv.end, k2, nv, i2))
                else:
                    heapq.heappush(q, (pos2, end2, k2, v2, i2))
            yield (pos, grp)

    @staticmethod
    def ovl(g, v):
        if not len(g.alt) == 1:
            return False
        return v.pos >= g.pos and v.pos < g.end

    @staticmethod
    def split_g(g, v):
        gs = []
        if g.pos < v.pos:
            g1 = copy.deepcopy(g)
            g1.end = v.pos
            g1.line = None
            gs.append(g1)
        if g.end > v.end:
            g.pos = v.end
            g.line = None
            gs.append(g)
        return gs
    
    def __shard__(self, cse):
        self.shard = cse
        self.vcftmp = self.gvcf_out.__shard__(cse)
        return self    

    def __getstate__(self):
        odict = self.__dict__.copy()
        return odict
        
    def __getdata__(self):
        return self.vcftmp.__getdata__()
    
    def __accum__(self, data):
        if not self.stdout:
            self.gvcf_out.__accum__(data)
        else:
            with open(data) as f:
                for line in f:
                    sys.stdout.write(line)
        os.unlink(data)

    def combine(self, shard=None):
        shard = shard or self.shard
        chrom, start, end = shard
        vcfs = [vcf.__shard__(shard) for vcf in (self.gvcf_in, self.vcf_in)]
        for pos, grp in self.grouper(*vcfs):
            g, v = grp
            if v is None:
                if len(g.alt) == 1:
                    if pos < start:
                        g.pos = start
                        g.line = None
                    if g.end > end:
                        g.info['END'] = end
                        g.line = None
                    self.vcftmp.emit(g)
                elif pos >= start:
                    g.samples[0]['GT'] = '0/0'
                    g.line = None
                    self.vcftmp.emit(g)
            if v and pos >= start:
                v.alt.append('<NON_REF>')
                v.line = None
                self.vcftmp.emit(v)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("gvcf", help="Input GVCF file")
    parser.add_argument("vcf", help="Input VCF file")
    parser.add_argument("gvcf_out", help="Output GVCF file")
    parser.add_argument("-t", "--threads", help="Number of threads", type=int)
    
    args = parser.parse_args()
    combiner = Combiner(args.gvcf, args.vcf, args.gvcf_out)

    t0 = time.time()
    nthr = args.threads
    step = 10*1000*1000

    sharder = vcflib.Sharder(nthr)
    contigs = ((c,0,int(t['length'])) for c,t in iteritems(combiner.contigs))
    shards = sharder.cut(contigs, step)
    _ = sharder.run(shards, Combiner.combine, [], combiner)

    t1 = time.time()
    mm, ut, st = 0, 0, 0
    for who in (resource.RUSAGE_SELF, resource.RUSAGE_CHILDREN):
        ru = resource.getrusage(who)
        mm += ru.ru_maxrss * 1024
        ut += ru.ru_utime
        st += ru.ru_stime
    print('overall: %d mem %.3f user %.3f sys %.3f real' %
        (mm, ut, st, t1-t0), file=sys.stderr)

    return 0

if __name__ == '__main__':
    sys.exit(main())

# vim: ts=4 sw=4 expandtab
