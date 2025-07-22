# 1 安装Nextflow和nf-core

```bash
# python=3.13.5
conda install nextflow nf-core
```

或者也可以创建一个新的conda环境，避免冲突

```bash
# use python 3.12 as example
conda create --name nf-core python=3.12 nf-core nextflow
conda activate nf-core
```

作者建议nf-core需要经常性的更新

```bash
# update nf-core
conda update nf-core
```

安装完成后，试试第一个流程hello，看看相关依赖是否都正常工作

```bash
# hello, world!
nextflow hello

# Later you will see "hello" in multi-languages
```

如果看到下面的输出，就表示nextflow和nf-core安装完成了

>(RNAseq) tl5024@iyun50:~/Lregale/analyse$ nextflow run hello
>
> N E X T F L O W   ~  version 25.04.6
>
>Pulling nextflow-io/hello ...
> downloaded from https://github.com/nextflow-io/hello.git
>Launching `https://github.com/nextflow-io/hello` [spontaneous_euclid] DSL2 - >revision: 2ce0b0e294 [master]
>
>executor >  local (4)
>[8e/33eb36] process > sayHello (1) [100%] 4 of 4 ✔
>Hola world!
>
>Hello world!
>
>Ciao world!
>
>Bonjour world!

# 2 查找可用的pipeline

可以在[nf-core官网](https://nf-co.re/pipelines)上查找

懒得开浏览器也可以直接在终端里查找

```bash
# nf-core version 3.32, a low version will report an error
nf-core pipelines list
```

不过终端里的查找同时也会显示你是否pull过以及是否有最新的release

这里以rnaseq 3.19.0 的[pipeline](https://github.com/nf-core/rnaseq)为例

下面是这个pipeline的示意图，如果看不到说明你需要些某科学的上网技巧

如果看不到下面的动图，那么相比做分析先学会特殊上网技巧才是你的当务之急

![nf-core pipeline of RNA-seq](https://github.com/nf-core/rnaseq/raw/master/docs/images/nf-core-rnaseq_metro_map_grey_animated.svg)

# 3 RNA-seq

首先测试一下rnaseq流程

```bash
# test
nextflow run nf-core/rnaseq -profile test,docker --outdir "/media/desk16/tl5024/Lreg_rnaseq"
```

```bash
# nf-core rnaseq pipeline只接受fq.gz或者fastq.gz文件
# 递归压缩RNAseq目录下所有fastq文件为fastq.gz
find /media/desk16/tl5024/RNAseq/ -type f -name "*.fastq" | parallel -j 4 'pigz -p 10 "{}"'
```

因为参数很多，所以可以通过json设置

下面是[nf-params.json](https://github.com/ruellia-zhu/myCodes/blob/main/nf-core-rnaseq/nf-params.json)的一个示例

```json
{
    "input": "/media/desk16/tl5024/Lreg_rnaseq/samplesheet.csv",
    "outdir": "/media/desk16/tl5024/Lreg_rnaseq",
    "email": "yifeidangnian@gmail.com",
    "fasta": "/media/desk16/tl5024/Lregale/lrv2_mc.fa",
    "gff": "/media/desk16/tl5024/Lregale/lrv2.gff",
    "gtf_extra_attributes": "gene_name",
    "gtf_group_features": "gene_id",
    "featurecounts_group_type": "gene_biotype",
    "featurecounts_feature_type": "exon",
    "igenomes_base": "s3://ngi-igenomes/igenomes/",
    "trimmer": "trimgalore",
    "min_trimmed_reads": 10000,
    "bam_csi_index": true,
    "ribo_database_manifest": "${projectDir}/workflows/rnaseq/assets/rrna-db-defaults.txt",
    "umi_dedup_tool": "umitools",
    "umitools_extract_method": "string",
    "umitools_grouping_method": "directional",
    "aligner": "star_rsem",
    "pseudo_aligner_kmer_size": 31,
    "min_mapped_reads": 5.0,
    "kallisto_quant_fraglen": 200,
    "kallisto_quant_fraglen_sd": 200,
    "stranded_threshold": 0.8,
    "unstranded_threshold": 0.1,
    "extra_fqlint_args": "--disable-validator P001",
    "deseq2_vst": true,
    "rseqc_modules": "bam_stat,inner_distance,infer_experiment,junction_annotation,junction_saturation,read_distribution,read_duplication",
    "bracken_precision": "S",
    "skip_bbsplit": true,
    "skip_preseq": true,
    "custom_config_version": "master",
    "custom_config_base": "https://raw.githubusercontent.com/nf-core/configs/master",
    "publish_dir_mode": "copy",
    "max_multiqc_email_size": "25.MB",
    "validate_params": true,
    "pipelines_testdata_base_path": "https://raw.githubusercontent.com/nf-core/test-datasets/7f1614baeb0ddf66e60be78c3d9fa55440465ac8/"
}
```

其中的samplesheet.csv是一个逗号分割的csv，包括4列，典型如下：

```csv
sample,fastq_1,fastq_2,strandedness
leaf,/media/desk16/tl5024/RNAseq/ERR15075318/ERR15075318_1.fastq.gz,/media/desk16/tl5024/RNAseq/ERR15075318/ERR15075318_2.fastq.gz,auto
stem,/media/desk16/tl5024/RNAseq/ERR15076288/ERR15076288_1.fastq.gz,/media/desk16/tl5024/RNAseq/ERR15076288/ERR15076288_2.fastq.gz,auto
bulb,/media/desk16/tl5024/RNAseq/ERR15076289/ERR15076289_1.fastq.gz,/media/desk16/tl5024/RNAseq/ERR15076289/ERR15076289_2.fastq.gz,auto
flower,/media/desk16/tl5024/RNAseq/ERR15076290/ERR15076290_1.fastq.gz,/media/desk16/tl5024/RNAseq/ERR15076290/ERR15076290_2.fastq.gz,auto
root-6dpi,/media/desk16/tl5024/RNAseq/SRR28536713/SRR28536713_1.fastq.gz,/media/desk16/tl5024/RNAseq/SRR28536713/SRR28536713_2.fastq.gz,auto
root-4dpi,/media/desk16/tl5024/RNAseq/SRR28536714/SRR28536714_1.fastq.gz,/media/desk16/tl5024/RNAseq/SRR28536714/SRR28536714_2.fastq.gz,auto
root-2dpi,/media/desk16/tl5024/RNAseq/SRR28536715/SRR28536715_1.fastq.gz,/media/desk16/tl5024/RNAseq/SRR28536715/SRR28536715_2.fastq.gz,auto
root-CK,/media/desk16/tl5024/RNAseq/SRR28536716/SRR28536716_1.fastq.gz,/media/desk16/tl5024/RNAseq/SRR28536716/SRR28536716_2.fastq.gz,auto
leaf-CK,/media/desk16/tl5024/RNAseq/SRR7031451/SRR7031451_1.fastq.gz,/media/desk16/tl5024/RNAseq/SRR7031451/SRR7031451_2.fastq.gz,auto
leaf-Be4h,/media/desk16/tl5024/RNAseq/SRR7031452/SRR7031452_1.fastq.gz,/media/desk16/tl5024/RNAseq/SRR7031452/SRR7031452_2.fastq.gz,auto
leaf-Be24h,/media/desk16/tl5024/RNAseq/SRR7031454/SRR7031454_1.fastq.gz,/media/desk16/tl5024/RNAseq/SRR7031454/SRR7031454_2.fastq.gz,auto
```

第1列是样品ID，如果一致会被认为是技术重复自动合并

第2、3列是双端fq.gz的路径，我这里用的绝对路径防止bug

第4列一般auto即可

怎么快速获得这个samplesheet.csv呢？我这有个脚本

参见[generate_sheetsample.py](https://github.com/ruellia-zhu/myCodes/blob/main/nf-core-rnaseq/generate_samplesheet.py)

用法很简单，就2参数，-i <样品存放路径> -o </path/to/samplesheet.csv>

```bash
python3 generate_samplesheet.py -i '~/RNAseq' -o '~/analysis/samplesheet.csv'
```

大基因组的chr/contig过长，可能超过Java的整数限制2^31-1，因此对于大基因组需要设置json:

```json
{
  ...其他参数
  "bam_csi_index": true,
  ...其他参数
}
```

设置好samplesheet.csv和nf-params.json后，开始运行

记住-name后面的run name每次都得不一样才行

-profile可以接conda docker singularity，都行，我这用的docker稳一点

```bash
nextflow run nf-core/rnaseq -r 3.19.0 -name "Lregale_rnaseq20250723-5" -work-dir "/media/desk16/tl5024/Lreg_rnaseq" -params-file "/media/desk16/tl5024/Lreg_rnaseq/nf-params.json" -profile docker
```

运行后会需要一段时间pull docker image，所以会表现为假死很久，多等会就行，没报错就算胜利

顺利的话过一段时间就可以看到服务器在跑了，terminal里可以看到各个process

挂机的话，可以限制一下运行的资源

```bash
nextflow run nf-core/rnaseq -r 3.19.0 -name "Lregale_rnaseq20250723-5" -work-dir "/media/desk16/tl5024/Lreg_rnaseq" -params-file "/media/desk16/tl5024/Lreg_rnaseq/nf-params.json" -profile docker --max_cpus 20 --max_memory 256.GB --max_time '72.h'
```
