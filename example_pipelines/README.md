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
    - [ml-model][ml-model] - TNscope somatic variant calling using a machine learning model.
    - [amplicon][tn-amplicon] - Somatic variant calling from amplicon data using a matched normal sample.
    - [ctDNA][tn-ctDNA] - Somatic variant calling from cfDNA tumor-only samples.
    - [ctDNA with UMI][tn-ctDNA-umi] - Somatic variant calling from cfDNA tumor-only samples with UMI-tagged reads.
    - [target capture][tn-targeted] - Somatic variant calling from hybridization catpure data using a matched normal sample.
  - *TNseq* - Pipelines matching best practice implementations for somatic variant calling.
    - [tumor-normal][tn-paired] - Somatic variant calling with TNsnv and TNhaplotyper from paired tumor-normal samples.
    - [tumor-only][tumor-only] - Somatic variant calling with TNsnv and TNhaplotyper from tumor-only samples.

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
[ml-model]: somatic/TNscope/ml-model.sh
[tn-amplicon]: somatic/TNscope/Somatic_amplicon_panel.sh
[tn-ctDNA]: somatic/TNscope/Somatic_ctDNA_without_UMI.sh
[tn-ctDNA-umi]: somatic/TNscope/Somatic_ctDNA_with_UMI.sh
[tn-targeted]: somatic/TNscope/Somatic_hybridization_panel.sh
[tn-paired]: somatic/TNseq/tumor_normal.sh
[tumor-only]: somatic/TNseq/tumor_only.sh
