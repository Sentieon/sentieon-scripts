# DNAscope LongRead pipeline scripts

Sentieon DNAscope LongRead is a pipeline for calling germline SNVs and indels from long-read sequence data. The DNAscope LongRead pipeline is able to take advantage of longer read lengths to perform quick and accurate variant calling using specially calibrated machine learning models.

This folder provides a shell script implementing the DNAscope LongRead pipeline along with Python scripts for VCF manipulation. The pipeline expects aligned reads as input in BAM or CRAM format and will output variants in the VCF (and gVCF) formats.

DNAscope LongRead is implemented using the Sentieon software package, which requires a valid license for use. Please contact info@sentieon.com for access to the Sentieon software and an evaluation license.

## Prerequisites

In order to run the pipeline scripts in this folder, you need to use the Sentieon software package version 202308 or higher. The pipeline also requires [Python] versions >2.7 or >3.3, [bcftools] version 1.10 or higher, and [bedtools]. The `sentieon`, `python`, `bcftools`, and `bedtools` executables will be accessed through the user's `PATH` environment variable.

## Input data requirements

### Aligned reads

The pipeline will accept Oxford Nanopore (ONT) long reads that have been aligned to the reference genome with `minimap2`. When aligning reads with Sentieon minimap2, using the Sentieon model for ONT is recommended. When aligning reads with the open-source minimap2, `-x map-ont` is recommended.

### The Reference genome

DNAscope LongRead will call variants present in the sample relative to a high quality reference genome sequence. Besides the reference genome file, a samtools fasta index file (.fai) needs to be present. We recommend aligning to a reference genome without alternate contigs.


## Usage

A single command is run to call variants:
```sh
dnascope_ONT.sh [-h] -r REFERENCE -i INPUT_BAM -m MODEL_BUNDLE [-d DBSNP_VCF] [-b DIPLOID_BED] [-t NUMBER_THREADS] [-g]  [--] VARIANT_VCF
```

The Sentieon LongRead pipeline for ONT requires the following arguments:
- `-r REFERENCE`: the location of the reference FASTA file. You should make sure that the reference is the same as the one used in the mapping stage.
- `-i INPUT_BAM`: the location of the input BAM or CRAM file.
- `-m MODEL_BUNDLE`: the location of the model bundle for ONT.

The Sentieon LongRead pipeline for ONT accepts the following optional arguments:
- `-d DBSNP_VCF`: the location of the Single Nucleotide Polymorphism database (dbSNP) used to label known variants. Only one file is supported. Supplying this file will annotate variants with their dbSNP refSNP ID numbers.
- `-b DIPLOID_BED`: interval in the reference to restrict variant calling, in BED file format. Supplying this file will limit variant calling to the intervals inside the BED file.
- `-t NUMBER_THREADS`: number of computing threads that will be used by the software to run parallel processes. The argument is optional; if omitted, the pipeline will use as many threads as the server has.
- `-g`: output variants in the gVCF format, in addition to the VCF output file. The tool will output a bgzip compressed gVCF file with a corresponding index file.
- `-h`: print the command-line help and exit.

The Sentieon LongRead pipeline for ONT requires the following positional arguments:
- `VARIANT_VCF`: the location and filename of the variant calling output. The tool will output a bgzip compressed VCF file with a corresponding index file.

## Pipeline output

The Sentieon LongRead pipeline for ONT will output a bgzip compressed file (.vcf.gz) containing variant calls in the standard VCFv4.2 format along with a tabix index file (.vcf.gz.tbi). If the `-g` option is used, the pipeline will also output a bgzip compressed file (.g.vcf.gz) containing variant calls in the gVCF format along with a tabix index file (.g.vcf.gz.tbi).

## Other considerations

Currently, the pipeline is only recommended for use with samples from diploid organisms. For samples with both diploid and haploid chromosomes, the `-b DIPLOID_BED` option can be used to limit variant calling to diploid chromosomes.

## References
**[Sentieon DNAscope LongRead – A highly Accurate, Fast, and Efficient Pipeline for Germline Variant Calling from PacBio HiFi reads]** - A preprint describing the DNAscope LongRead pipeline for calling variants from PacBio HiFi data.


[Python]: https://www.python.org/
[bcftools]: http://samtools.github.io/bcftools/bcftools.html
[bedtools]: https://bedtools.readthedocs.io/en/latest/

[Sentieon DNAscope LongRead – A highly Accurate, Fast, and Efficient Pipeline for Germline Variant Calling from PacBio HiFi reads]: https://www.biorxiv.org/content/10.1101/2022.06.01.494452v1