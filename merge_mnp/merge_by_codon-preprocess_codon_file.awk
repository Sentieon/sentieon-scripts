# 
# Copyright (c) Sentieon Inc. All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
# 
# * Redistributions of source code must retain the above copyright notice, this
#   list of conditions and the following disclaimer.
# 
# * Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#  
# 
# Script to pre-process a BED file of codons and create a list of loci and
# and the codon they belong to.

{
	#read exon line
	transcript=$4;chr=$1;start=$2;end=$3;
	#Need to do this because the NCBI RefSeq Curated table in https://genome.ucsc.edu/cgi-bin/hgTable has changed the transcriptID format
	split(transcript,v,".");transcript=v[1]
	#check if new transcript
	if (transcript != current_transcript)
	{
		#make sure last codon finished in the exon, otherwise error out
		if (codon_end != previous_end)
		{
			print "Error, transcript "current_transcript" in "chr" did not end in a complete codon; the last bases of the transcript will not be included in a codon." > "/dev/stderr"
			#Do not exit as the latest version of the NCBI RefSeq Curated table in https://genome.ucsc.edu/cgi-bin/hgTable has transcripts containing a number of bases not multiple of 3
			#exit
		}
		#set new info
		codon_start=start+1
		current_transcript = transcript
		#print "debug: Calculating codons for transcript "current_transcript
		codon_id=0
	}
	#otherwise this is a new exon from the same transcript
	else
	{
		#print "debug: codon_start "codon_start" codon_mid "codon_mid" codon_end "codon_end
		# ALL the exon starts are +1 for the codon, so the first codon on an exon always is on start+1, not only the first one in the transcript
		if (codon_start == -1) codon_start = start + 1
		else if (codon_mid == -1) codon_mid = start + 1
			else if (codon_end == -1) codon_end = start + 1
	}
	#print "debug: new exon "chr":"start"-"end
	previous_end = end
	beyond_exon=0
	#codons fully contained in the exon
	while (!beyond_exon)
	{
		if (codon_start >= start) codon_mid=codon_start+1
		if (codon_mid > end)
		{
			codon_mid = -1
			codon_end = -1
			beyond_exon = 1
		}
		else
		{
			if (codon_mid >= start) codon_end=codon_mid+1
		}
		if (codon_end > end)
		{
			codon_end = -1
			beyond_exon = 1
		}
		if ((codon_end != -1) && (codon_mid != -1))
		{
			codon_id = codon_id + 1
			codon_name=transcript"-Codon"codon_id
			print chr":"codon_start"\t"codon_name
			print chr":"codon_mid"\t"codon_name
			print chr":"codon_end"\t"codon_name
			codon_start=codon_end+1
		}
		#check if codon ended in the exon
		if (codon_start > end)
		{
			codon_start = -1
			beyond_exon = 1
		}
	}
}

