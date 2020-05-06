# MNP
Merge neighboring variants along the same haplotype into multi-nucleotide polymorphisms (MNPs)

Combines neighboring PASS variants on the same haplotype while supporting user-defined merging
strategies. Original variants that are combined will be marked as "MERGED".

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
```
zcat GRCh38_ucsc_refseq_good_transcripts.noALT.bed.zip |awk -f merge_by_codon-preprocess_codon_file.awk > codons.txt
export PYTHONPATH="$SENTIEON_INSTALL_DIR/lib/python/sentieon":"$PYTHONPATH"
python merge_mnp.py input.vcf[.gz] ucsc.hg19.fasta merge_by_codon codons.txt 1> output_noINDEL.vcf
$SENTIEON_INSTALL_DIR/bin/sentieon util vcfconvert output_noINDEL.vcf output_noINDEL.vcf.gz
```
