#!/usr/bin/env python

'''
Copyright (c) Sentieon Inc. All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are met:
  
  * Redistributions of source code must retain the above copyright notice, this
    list of conditions and the following disclaimer.
  
  * Redistributions in binary form must reproduce the above copyright notice,
    this list of conditions and the following disclaimer in the documentation
    and/or other materials provided with the distribution.
  
  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
  OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
'''

import argparse
import heapq
import json
import struct
import sys
import urllib.parse
import zlib

try:
    import botocore.session
    def s3open(path):
        if path.startswith('s3://'):
            return S3File(path)
        return open(path, 'rb')
except ImportError:
    def s3open(path):
        return open(path, 'rb')

class S3File(object):
    def __init__(self, path):
        url = urllib.parse.urlparse(path)
        self.bucket = url.netloc
        self.key = url.path[1:]
        self.iosize = 1024*1024
        self.reset()
        self.offset = 0
        self.client = botocore.session.get_session().create_client('s3')
        self.filesize = self.size()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, tb):
        self.close()

    def close(self):
        self.reset()
        self.client = None

    def reset(self):
        self.buffer = None
        self.buflen = 0
        self.bufptr = 0

    def seek(self, off, whence=0):
        if whence < 0 or whence > 2:
            raise ValueError
        if whence == 1:
            off += self.offset
        elif whence == 2:
            off = self.filesize + off
        if off >= self.offset - self.buflen and off < self.offset:
            self.bufptr = off - self.offset + self.buflen
        else:
            self.reset()
            self.offset = off

    def tell(self):
        return self.offset - self.buflen + self.bufptr

    def read(self, size=None):
        if self.client is None:
            raise RuntimeError('File closed')
        if size is None:
            size = sys.maxsize
        data = b''
        while size > 0:
            if self.bufptr == self.buflen:
                if self.fill() <= 0:
                    break
            n = min(self.buflen - self.bufptr, size)
            data += self.buffer[self.bufptr:self.bufptr+n]
            self.bufptr += n
            size -= n
        return data

    def fill(self):
        self.reset()
        s, e = (self.offset, self.offset+self.iosize-1)
        if s >= self.filesize:
            return 0
        e = min(e, self.filesize)
        range = 'bytes=%d-%d' % (s, e)
        r = self.client.get_object(Bucket=self.bucket, Key=self.key,
            Range=range)
        self.buffer = r['Body'].read()
        self.buflen = len(self.buffer)
        self.bufptr = 0
        self.offset += self.buflen
        return self.buflen

    def size(self):
        r = self.client.head_object(Bucket=self.bucket, Key=self.key)
        return r['ContentLength']

class BAMIndex(object):
    SHIFTS = (14, 17, 20, 23, 26, 29)
    MAXBIN = ((1 << SHIFTS[-1]-SHIFTS[0]+3) - 1) // 7 + 1
    MAGIC = 0x01494142

    def __init__(self, bamf):
        self.chrs = None
        self.minoff = sys.maxsize
        self.maxoff = 0
        base = bamf
        while base:
            try:
                self.load(base + '.bai')
                base = None
                break
            except:
                pass
            if not base.endswith('.bam'):
                break
            base = base[:-4]
        if base is not None:
            raise RuntimeError('Failed to load the index of ' + bamf)

    def load(self, idxf):
        with s3open(idxf) as fp:
            s4 = struct.Struct('<I')
            s8 = struct.Struct('<Q')
            data = fp.read()
            off = 0
            magic, = s4.unpack_from(data, off); off += s4.size
            if magic != self.MAGIC:
                raise RuntimeError('Not a bam index file')
            chrs = []
            n_ref, = s4.unpack_from(data, off); off += s4.size
            for _ in range(n_ref):
                bins = {}
                n_bin, = s4.unpack_from(data, off); off += s4.size
                for _ in range(n_bin):
                    bin, = s4.unpack_from(data, off); off += s4.size
                    chunks = []
                    n_chunk, = s4.unpack_from(data, off); off += s4.size
                    for _ in range(n_chunk):
                        s, = s8.unpack_from(data, off); off += s8.size
                        e, = s8.unpack_from(data, off); off += s8.size
                        chunks.append((s, e))
                    if bin < self.MAXBIN and n_chunk > 0:
                        self.minoff = min(self.minoff, chunks[0][0])
                        self.maxoff = max(self.maxoff, chunks[-1][1])
                    bins[bin] = chunks
                intvs = []
                n_intv, = s4.unpack_from(data, off); off += s4.size
                for _ in range(n_intv):
                    o, = s8.unpack_from(data, off); off += s8.size
                    intvs.append(o)
                if n_intv == 0:
                    intvs.append(0)
                chrs.append((bins, intvs))
            self.chrs = chrs

    def save(self, idxf):
        with open(idxf, 'wb') as fp:
            s4 = struct.Struct('<I')
            s8 = struct.Struct('<Q')
            fp.write(s4.pack(self.MAGIC))
            fp.write(s4.pack(len(self.chrs)))
            for bins, intvs in self.chrs:
                fp.write(s4.pack(len(bins)))
                for bin in sorted(bins.keys()):
                    chunks = bins[bin]
                    fp.write(s4.pack(bin))
                    fp.write(s4.pack(len(chunks)))
                    for s,e in chunks:
                        fp.write(s8.pack(s))
                        fp.write(s8.pack(e))
                fp.write(s4.pack(len(intvs)))
                for o in intvs:
                    fp.write(s8.pack(o))

    def query(self, tid, s, e):
        ranges = []
        if tid < 0:
            ranges.append((self.maxoff, sys.maxsize))
            return ranges
        if tid >= len(self.chrs):
            return ranges
        ci = self.chrs[tid]
        i = s >> self.SHIFTS[0]
        minoff = i >= len(ci[1]) and ci[1][-1] or ci[1][i]
        for shift in reversed(self.SHIFTS):
            bo = ((1 << 29-shift) - 1) // 7
            bs = bo + (s >> shift)
            be = bo + (e-1 >> shift)
            for bi in range(bs, be+1):
                if bi not in ci[0]:
                    continue
                for chunk in ci[0][bi]:
                    if chunk[1] > minoff:
                        ranges.append(chunk)
        return ranges

    def dump(self):
        print('minoff', self.minoff, 'maxoff', self.maxoff)
        for tid, (bins, intvs) in enumerate(self.chrs):
            if len(bins) == 0:
                continue
            print('tid', tid)
            print(json.dumps((bins, intvs), indent=4))

    def density(self):
        density = []
        for bins, intvs in self.chrs:
            dens = [0] * len(intvs)
            minoff, maxoff = (1<<63)-1, 0
            for bin, chunks in bins.items():
                minoff = min(minoff, chunks[0][0]>>16)
                maxoff = max(maxoff, chunks[-1][1]>>16)
                size = sum([(e>>16) - (s>>16) for s, e in chunks])
                bi, bs = 0, 0
                for shift in self.SHIFTS:
                    bo = ((1 << 29-shift) - 1) // 7
                    if bin >= bo and bin <= bo * 8:
                        bs = 1 << shift-14
                        bi = (bin - bo) * bs
                        break
                if bs == 0:
                    continue
                if len(dens) <= bi:
                    continue
                if len(dens) < bi + bs:
                    bs = len(dens) - bi
                for i in range(bs):
                    dens[bi + i] += size / bs
            size = max(maxoff - minoff, 0)
            density.append((size, dens))
        return density

class BGZF(object):
    def __init__(self, fp):
        self.fp = fp
        self.pos = 0
        self.block = None
        self.blkoff = 0
        self.blkend = 0

    def seek(self, off, whence):
        if whence < 0 or whence == 1 and off != 0 or whence >= 2:
            raise ValueError
        off = off * whence + self.pos
        if self.pos >> 16 != off >> 16:
            self.block = None
        self.pos = off

    def tell(self):
        return self.pos

    def read(self, size):
        if size is None:
            raise ValueError
        data = b''
        blk = self.pos >> 16
        off = self.pos & 65535
        while size > 0:
            if self.blkoff != blk or self.block is None:
                self.read_block(blk)
            end = min(off + size, len(self.block))
            data += self.block[off:end]
            size -= end - off
            if end < len(self.block):
                off = end
                break
            blk = self.blkend
            off = 0
        self.pos = blk << 16 | off
        return data

    def close(self):
        self.fp.close()

    def read_block(self, offset):
        self.fp.seek(offset)
        header = self.fp.read(18)
        length = (header[16] | header[17] << 8) + 1
        block = self.fp.read(length-18)
        zs = zlib.decompressobj(-15)
        self.block = zs.decompress(block[:-8]) + zs.flush()
        self.blkoff = offset
        self.blkend = offset + length

class BAMHeader(object):
    def __init__(self, bamf):
        self.chrs = []
        self.tids = {}
        self.texts = None
        if bamf:
            try:
                self.load(bamf)
                return
            except:
                pass
        base = bamf
        while base:
            try:
                self.load_text(base + '.hdr')
                base = None
                break
            except:
                pass
            if not base.endswith('.bam'):
                break
            base = base[:-4]
        if base is not None:
            raise RuntimeError('Failed to load the header of ' + bamf)

    def load(self, bamf):
        bgzf = BGZF(s3open(bamf))
        d = bgzf.read(4)
        if d != b'BAM\001':
            raise RuntimeError('Not a bam file')
        d = bgzf.read(4)
        n = d[0] | d[1] << 8 | d[2] << 16 | d[3] << 24
        self.texts = bgzf.read(n).decode('utf-8')
        d = bgzf.read(4)
        n = d[0] | d[1] << 8 | d[2] << 16 | d[3] << 24
        chrs = []
        tids = {}
        for tid in range(n):
            d = bgzf.read(4)
            n = d[0] | d[1] << 8 | d[2] << 16 | d[3] << 24
            d = bgzf.read(n)
            name = d[:-1].decode('utf-8')
            d = bgzf.read(4)
            n = d[0] | d[1] << 8 | d[2] << 16 | d[3] << 24
            chrs.append((name, n))
            tids[name] = tid
        bgzf.close()
        self.chrs = chrs
        self.tids = tids

    def load_text(self, hdrf):
        with open(hdrf, 'r') as fp:
            chrs = []
            for line in fp:
                flds = line.rstrip('\r\n').split('\t')
                if flds[0] == '@SQ':
                    attr = dict(f.split(':',1) for f in flds[1:])
                    chrs.append((attr['SN'], int(attr['LN'])))
        self.chrs = chrs
        self.tids = dict((c[0], i) for i,c in enumerate(chrs))

def parse_shard(arg):
    shards = []
    for a in arg.split(','):
        c, se = a.rsplit(':', 1)
        s, e = map(int, se.split('-'))
        shards.append((c, s-1, e))
    return (shards, arg)

def main():
    fraglen = 1000000

    parser = argparse.ArgumentParser(description='A memory usage estimator.')
    parser.add_argument('-i', '--input', metavar='BAM', action='append',
       help='Input BAM file name, required', required=True)
    parser.add_argument('-t', '--thread_count', metavar='NUM', dest='nthr',
       help='Number of threads, required', type=int, required=True)
    parser.add_argument('-f', '--frag_len', metavar='LEN', dest='fraglen',
       help='Fragment length, default %d' % fraglen, type=int, default=fraglen)
    parser.add_argument('-q', '--qcal', action='store_true',
       help='Whether qcal correction will be applied')
    parser.add_argument('-d', '--densest', metavar='N',
       help='List the densest N fragments', type=int, default=0)
    parser.add_argument('shard', nargs='*',
       help='A list of shards. Each shard may be a comma seperated list of '
       'genomic locations, default one shard per chromosome')
    args = parser.parse_args()

    chroms = {}
    density = {}
    for bamf in args.input:
        hdr = BAMHeader(bamf)
        idx = BAMIndex(bamf)
        for (tchr, tlen), (size, dens) in zip(hdr.chrs, idx.density()):
            d = density.setdefault(tchr,[])
            if len(d) < len(dens):
                d.extend([0] * (len(dens)-len(d)))
            for i in range(len(dens)):
                d[i] += dens[i]
            chroms[tchr] = max(chroms.get(tchr,0), tlen)
    chr_tid = dict((chr[0],i) for i, chr in enumerate(hdr.chrs))

    multiplier = 4 * 2.0
    if args.qcal:
        multiplier *= 2
    fraglen = args.fraglen
    nthr = args.nthr
    if len(args.shard) > 0:
        shards = map(parse_shard, args.shard)
    else:
        shards = [([(c,0,n)], c) for c,n in chroms.items()]
    try:
        shards.sort(key=lambda shard: chr_tid[shard[1]])
    except KeyError:
        shards.sort(key=lambda shard: shard[1])
    nmax = [(0,'',0)] * args.densest
    for shard, name in shards:
        dmax = [0] * nthr
        for c,s,e in shard:
            if c not in density:
                continue
            for f in range(s, e, fraglen):
                i = f // 16384
                j = (f + fraglen + 16383) // 16384
                for k in range(i,j):
                    if k >= len(density[c]):
                        break
                    d = density[c][k]
                    heapq.heappushpop(dmax, d)
                    heapq.heappushpop(nmax, (d,c,k*16384))
        d = sum(dmax) * multiplier / (1024*1024*1024)
        print('%s\t%f' % (name, d))

    if args.densest > 0:
        print('\nThe densest %d fragments:' % args.densest)
        for d,c,s in heapq.nlargest(args.densest, nmax):
            if c is None: continue
            frag = '%s:%d-%d' % (c, s+1, s+16384)
            print('%-32s%f' % (frag, d/(1024*1024*1024.)))

if __name__ == '__main__':
    sys.exit(main())

# vim: ts=4 sw=4 expandtab
