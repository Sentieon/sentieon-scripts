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

Script defining the merging strategy for merging MNPs using a codon
file.

'''

class merge_by_codon:

    codon_dict = {}
    ignore_INDELs = False

    def __init__(self, codon_file, ignore_INDELs=0):
        self.ignore_INDELs = bool(int(ignore_INDELs))
	#load the codon file and create a dictionary of all loci and the codon they belong to
	#the codon file is a tab separated list of locus\tCodonID
        with open(codon_file) as fp:
            for line in fp:
                line_contents = line.rstrip().split('\t')
                self.codon_dict[line_contents[0]] = line_contents[1]

    def get_codon(self, v):
	#return CodonID of all codons the variant covers
        codon_set = set()
        pos = v.pos + 1
        if len(v.ref) > 1:
            #DEL, need to check all codons it stradles
            DEL_length = len(v.ref)
            seqs = [v.chrom + ':' + str(i) for i in range(pos+1, pos+DEL_length)]
            codon_set = set([self.codon_dict[s] for s in seqs if s in self.codon_dict])
        else:
            locus = v.chrom + ':' + str(pos)
            if locus in self.codon_dict:
                codon_set.add(self.codon_dict[locus])
        return codon_set

    def is_merge(self, v1, v2):
        #determine if one of the 2 variants is INDEL, then ignore if ignore_INDELs=True
        isINDEL = False
        if len(v1.ref) > 1 or len(v2.ref) > 1:
            isINDEL = True
        for alt in v1.alt:
            if len(alt) > 1:
                isINDEL = True
        for alt in v2.alt:
            if len(alt) > 1:
                isINDEL = True
        if self.ignore_INDELs and isINDEL:
            return False
        #determine if the variants share a codon
        return bool(self.get_codon(v1).intersection(self.get_codon(v2))) 
