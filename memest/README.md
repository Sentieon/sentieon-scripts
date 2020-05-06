# Memest #
Memory usage estimation from a BAM header and index (.bai) files

Estimates the memory usage of Sentieon Haplotyper, Sentieon TNhaplotyper and Sentieon TNscope using only a BAM header and a BAM index file. Memest supports direct access to files in AWS S3, and downloading only the BAM header and index. Input to Memest is one or more sorted and indexed BAM files and, optionally, a list of shards to calculate memory across.


### Usage ###
```
usage: memest.py [-h] -i BAM -t NUM [-f LEN] [-q] [shard [shard ...]] [-d N]

A memory usage estimator.

positional arguments:
  shard                 A list of shards. Each shard may be a comma seperated
                        list of genomic locations, default one shard per
                        chromosome.

optional arguments:
  -h, --help            show this help message and exit
  -i BAM, --input BAM   Input BAM file name, required
  -t NUM, --thread_count NUM
                        Number of threads, required
  -f LEN, --frag_len LEN
                        Fragment length, default 1000000
  -q, --qcal            Whether qcal correction will be applied
  -d N, --densest N     List the densest N fragments
```

### Output ###

Memest ouputs the estimated memory usage of each chromosome or shard in GB, one chromosome/shard per line.

```
1	0.256133
2	1.317063
3	0.071602
...
```

### Example Usage - Standard Processing ###

Estimating the memory usage of three input BAM files running on a single thread.
```
python memest.py -i sample1.bam -i sample2.bam -i sample3.bam \
    -t 1 > memest_output.txt
```

Estimating the memory usage of a tumor-normal paired sample in AWS S3 including base-quality score recalibration and running on 16 threads.

```
python memest.py \
    -i s3://my-bucket/ngs/bams/sample1/tumor.bam \
    -i s3://my-bucket/ngs/bams/sample1/normal.bam \
    -q -t 16 > sample1_memest.txt
```

### Example Usage - Distributed Processing ###

Distributed processing divides the genome up into small, sequential, non-overlapping parts called shards and then processing each of these shards independently on different machines. Memest can estimate memeory usage across each shard from many BAM files.

Estimating the memory usage of many BAM files in AWS S3 across 100 Mb shards with 16 threads. Notice that the third shard is a comma-seperated list of two genomic locations.

```
python memest.py \
    -i s3://my-bucket/sample1.bam -i s3://my-bucket/sample2.bam ... \
    -t 16 chr1:1-100000000 chr1:100000001-200000000 \
    chr1:200000001-249250621,chr2:1-50749379
```
