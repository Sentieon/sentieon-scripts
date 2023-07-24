# Example Pipelines With the Sentieon tools #

Example Sentieon pipelines implemented in shell script. You may modify these pipelines to work with your existing bioinformatics infrastructure.

See Sentieon's [Argument Correspondence application note][correspondence] for information on how parameters in the Sentieon tools correspond to parameters in the GATK.

## Example pipelines ##
- **Germline variant calling**
  - [whole-genome sequencing][wgs] - Standard germline variant calling suitable for whole-genome samples.
  - [whole-genome sequencing with DNAscope][dnascope] - Higher accuracy germline variant calling for whole-genome samples.
  - [whole-exome sequencing][wes] - Standard germline variant calling pipeline across target intervals.
  - [whole-exome sequencing with DNAscope][dnascope-wes] - Higher accuracy germline variant calling for whole-exome samples.
  - [multiple input fastq][multi] - Standard germline variant calling pipeline, modified to support multiple pirrs of input FASTQ files.
  - [joint calling][joint] - Standard germline variant calling pipeline with joint genotyping of multiple samples.
  - [CCDG/functional equivalence pipeline][ccdg] - A Centers for Common Disease Genomics (CCDG) functional equivalent pipeline implemented in the Sentieon tools. See Sentieon's [Functional Equivalent Pipeline application note][fe-app] for more information.
  - [RNAseq][rna] - Germline variant calling from aligned RNAseq data. See the [Sentieon manual][rna-doc] for more information on this pipeline.
- **Somatic variant calling**
  - *TNscope* - Improved accuracy somatic variant calling.
    - [default][tnscope] - Standard somatic variant calling with TNscope.
    - [50x to 100x][50x] - Somatic variant calling with TNscope. Recommended settings for 50x to 100x coverage.
    - [300x or greater][300x] - Somatic variant calling with TNscope. Recommended settings for coverage over 300x.
    - [unique molecular identifier (UMI)][umi] - Somatic variant calling and filtering of UMI-tagged reads.
  - *TNseq* - Pipelines matching best practice implementations for somatic variant calling.
    - [tumor-normal][tn-paired] - Somatic variant calling with TNsnv and TNhaplotyper from paired tumor-normal samples.
    - [tumor-only][tumor-only] - Somatic variant calling with TNsnv and TNhaplotyper from tumor-only samples.
    - [TNhaplotyper2][tnhap2] - Somatic variant calling with TNhaplotyper2.

[correspondence]: https://support.sentieon.com/appnotes/arguments/
[wgs]: germline/DNAseq/wgs.sh
[dnascope]: germline/DNAscope/wgs.sh
[wes]: germline/DNAseq/wes-interval.sh
[dnascope-wes]: germline/DNAscope/wes-interval.sh
[multi]: germline/DNAseq/multi-FASTQ.sh
[joint]: germline/DNAseq/joint-calling.sh
[ccdg]: germline/DNAseq/ccdg_functional-equivalent.sh
[fe-app]: https://support.sentieon.com/appnotes/functional_equivalent/
[rna]: germline/DNAseq/RNAseq-calling.sh
[rna-doc]: https://support.sentieon.com/manual/RNA_call/rna/
[tnscope]: somatic/TNscope/default.sh
[50x]: somatic/TNscope/50x_to_100x.sh
[300x]: somatic/TNscope/300x_or_greater.sh
[umi]: somatic/TNscope/umi.sh
[tn-paired]: somatic/TNseq/tumor_normal.sh
[tumor-only]: somatic/TNseq/tumor_only.sh
[tnhap2]: somatic/TNseq/TNhaplotyper2.sh
