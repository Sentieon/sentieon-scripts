# MNP
Merge neighboring variants along the same haplotype into multi-nucleotide polymorphisms (MNPs)

Combines neighboring PASS variants on the same haplotype while supporting user-defined merging
strategies. Variants in the same haplotype are defined as those with common PID/PGT FORMAT fields
(as used by Sentieon) or with the same PS FORMAT field (as used by
[whatshap](https://whatshap.readthedocs.io/en/latest/guide.html)). Original variants that are
combined will be marked as "MERGED".

Requirements:
This script requires access to the vcflib library contained in the Sentieon
software package located in $SENTIEON_INSTALL_DIR/lib/python/sentieon.

### Usage ###
```
usage: merge_mnp.py [-h] [--out_file OUT_FILE] [--max_distance MAX_DISTANCE]
                    vcf_file reference [merge [merge ...]]

Merge phased SNPs into MNPs.

positional arguments:
  vcf_file              Input VCF file
  reference             The reference genome
  merge                 Optional merge options. First argument is the merge
                        file/class name. The rest are the arguments to
                        initialize the merge object.

optional arguments:
  -h, --help            show this help message and exit
  --out_file OUT_FILE   Ouptut VCF file. If not specified, it will be output
                        to stdout.
  --max_distance MAX_DISTANCE
                        Maximum distance between two variants to be merged.
```

To merge phased variants only if they belong to the same codon, you need to use the
merge_by_codon.py class definition file. When using the merge_by_codon.py class, the
first argument after the class is the codon file as a tab separated list of
locus\tCodonID, while the second argument is a boolean [0/1] to determine whether to
ignore INDELs in the merging.
To create the codon file, you can use the script merge_by_codon-preprocess_codon_file.awk
to pre-process it.

### Example Usage - Standard Processing ###
```
export PYTHONPATH="$SENTIEON_INSTALL_DIR/lib/python/sentieon":"$PYTHONPATH"
python merge_mnp.py input.vcf[.gz] ucsc.hg19.fasta > output.vcf
$SENTIEON_INSTALL_DIR/bin/sentieon util vcfconvert output.vcf output.vcf.gz
```

### Example Usage - Merging only variants within the same codon ###
#### Creating a codons.txt file ####
The `codons.txt` file is a tab-delimited user-supplied file containing a locus (in `<CHROM>:<POS>` format)
and a codon-ID, that is unique for each codon. Users might use their own resources to create this file,
but we provide the following directions and awk script for convenience:

- Download a BED file of RefSeq transcripts from the UCSC Table Browser.
  - Go to https://genome.ucsc.edu/cgi-bin/hgTables and:
    - Set group at "Genes and Gene Predictions".
    - Set track to "NCBI RefSeq?".
    - Set table to "RefSeq? Curated (ncbiRefSeqCurated)"
    - Set output format to "BED - browser extensible data"
  - Or use the direct link for hg38 - https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=961881833_BPsmSsIoBe8MFwDPEAlJlo9HVZCm&clade=mammal&org=Human&db=hg38&hgta_group=genes&hgta_track=refSeqComposite&hgta_table=ncbiRefSeqCurated&hgta_regionType=genome&hgta_outputType=bed
  - Press Get output
  - In the next panel, select "Coding Exons"
  - Press get BED and download the output file
- Pass the downloaded BED file to the supplied awk script:
  `cat <BED> | awk -f mnp_with_codons-create_codon_file.awk > codons.txt`

#### Usage ####
```
export PYTHONPATH="$SENTIEON_INSTALL_DIR/lib/python/sentieon":"$PYTHONPATH"
python merge_mnp.py input.vcf[.gz] ucsc.hg19.fasta merge_by_codon codons.txt 1> output_noINDEL.vcf
$SENTIEON_INSTALL_DIR/bin/sentieon util vcfconvert output_noINDEL.vcf output_noINDEL.vcf.gz
```

### Example Usage - Using an unphased VCF and a BAM ###
When using a VCF from Sentieon TNscope, the PID/PGT FORMAT fields will be used to determine
whether variants are in the same haplotype. If you have an unphased VCF, it is still possible
to use the script by pre-processing it to add the phasing information to it; one possible tool
to do this is [whatshap](https://whatshap.readthedocs.io/en/latest/guide.html) which given
an unphased VCF and the corresponding BAM file will add the phasing information.
```
whatshap phase -o phased.vcf --reference=reference.fasta input.vcf input.bam
$SENTIEON_INSTALL_DIR/bin/sentieon util vcfconvert phased.vcf phased.vcf.gz
export PYTHONPATH="$SENTIEON_INSTALL_DIR/lib/python/sentieon":"$PYTHONPATH"
python merge_mnp.py phased.vcf.gz reference.fasta > output.vcf
$SENTIEON_INSTALL_DIR/bin/sentieon util vcfconvert output.vcf output.vcf.gz
```
