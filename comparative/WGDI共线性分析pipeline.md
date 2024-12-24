# 比较基因组共线性分析（一）：共线性点阵图到ks计算

环境：(wsl) Ubuntu 24.04 x64, python 3.12, AMD R7-5700X3D 8C16T, 32G RAM

## 0 前言 

### 基因组共线性和基因组同线性

所谓基因组共线性（Genome Collinearity）指的是不同物种或同一物种的不同个体中，基因组的某些区域在基因顺序和结构上保持一致或相似的现象。

共线性主要是描述同一条染色体上基因的位置关系，也就是指由同一祖先型分化而来的不同物种间基因的类型以及相对顺序的保守性，换言之：

> 共线性 = 基因的同源性 + 基因的排列顺序

这里需要做出区分的是同线性（synteny），同线性指的是两个区域有一定数量的同源基因，而对基因排列顺序没有要求，共线性是同线性的一种特殊情况，要求同源基因的排列顺序也相似

### 基因组共线性与物种之间的分化时间

![共线性与物种之间的分化时间](image-1.png)

共线性片段的大小与物种之间的分化时间有很大关系。分化时间较短的物种之间，积累的变异较少，因此它们保留了更多从祖先遗传下来的特征；反之，分化时间较长的物种之间，积累的变异较多，因此导致它们共有的特征变少，即共线性片段较短。（上图）

### 直系同源和旁系同源

基因同源（homology）又可以分为直系同源（orthologs）和旁系同源（paralogs）。直系同源基因指存在于祖先基因组种，随后因为物种的分化，分别遗传给不同的后代，这些基因在结构和功能上往往有很高的相似性（尤其是如果它们对生物体至关重要、不可替代时），例如拟南芥和水稻种都有JAZ1，AtJAZ1和OsJAZ1就属于直系同源基因。旁系同源基因指的是同一基因组种由于基因复制而产生的同源基因，旁系同源基因在最初功能一致，但是由于没有选择压力，有的拷贝很可能发生分化，而其中一个拷贝的分化受到限制以确保其能实现所需的功能。

仅仅靠序列信息区分直系同源基因和旁系同源基因十分困难，因为除了上述功能上的分化，从最初的基因复制事件之后，还可能发生了多次基因复制、丢失、水平转移事件等等。

### 中性演化假说

根据中性演化假说，基因的变化大多数是中性突变，即不会影响生物生存的突变。

Ka指的是非同义替换位点数（aa发生改变），Ks指的是同义替换位点数（aa不发生改变）

如果一对同源基因分开的时间越早，Ks就越多，Ka/Ks就越低

### 共线性分析的应用

测序发展的初期，人们可获得的序列信息十分有限，无法全面分析基因功能。单个物种的基因组序列也无法挖掘真正的进化事件（基因丢失、基因获得等等）。

近十几年测序技术快速发展，大规模的全基因组测序（WGS）成为现实，比较基因组的出现，更进一步推动了近缘物种或品种的WGS，其中比较基因组的很大一部分工作就是在全基因组比对上，共线性分析是比较基因组中经典的分析策略。

共线性分析允许分析物种之间的不同尺度的分子进化事件。大尺度的进化事件包括：基因组重排事件、基因组复制事件（经常说的WGD）；小尺度进化事件包括：针对基因组水平的碱基替换速率以及插入、缺失事件。从共线性片段中可以识别出物种之间的小尺度和大尺度突变事件，这些事件可以揭示物种之间的进化关系，识别物种之间的基因功能保守性，追踪基因组进化历史，辅助基因注释和同源基因鉴定。

对同一个物种也可以做共线性分析，这些分子进化事件可以揭示基因组倍增历史，研究基因的扩增和分化，定位重要基因家族的分布（许多重要基因家族在基因组内分布于共线性区域），理解特定基因家族的起源、扩张和分化，理解基因组重排（倒位、易位，……）与染色体演化，指导基因组注释。

### WGDI：Whole Genome Duplication Intergrated Analysis

WGDI是一篇发表在Molecular Plant上的植物比较基因组工具，详情见文章，作者徐洲更

[WGDI: A user-friendly toolkit for evolutionary  analyses of whole-genome duplications and  ancestral karyotypes](10.1016/j.molp.2022.10.018)

以下是机翻的Abstract

> 在地球上大多数主要的生物谱系中都发现了全基因组重复（WGDs）和随后的核型变化的证据。为了阐明基因组分析中复杂的基因共线性多层模式，需要方便准确的工具包。为了满足这一需求，我们开发了WGDI（全基因组复制集成分析），这是一种基于Python的命令行工具，有助于全面分析递归多倍体事件和跨物种基因组比对。WGDI支持三个主要工作流程（多倍体推断，基因组同源性的层次推断和祖先染色体核型分析），可以基于高质量的染色体水平基因组改进WGD的检测和WGD相关事件的表征。值得注意的是，它可以提取完整的同线性块，并有助于重建详细的核型进化。该工具包可在GitHub免费获得(https://github.com/SunPengChuan/wgdi)。作为其应用的一个例子，WGDI令人信服地阐明了WGDs后Aquilegia coerulea和Vitis vinifera的核型进化，并拒绝了Aquilegia作为核心双子叶植物异源多倍体起源的亲本谱系的假设。

使用方法参见[WGDI的Github主页](https://github.com/SunPengChuan/wgdi)以及文档

### WGDI的三个关键features：

1. Polyploid Inference（多倍体推断）
Identifies and confirms polyploid events with high accuracy.
2. Genomic Homology Inference（基因组同源性推断）
Traces the evolutionary history of duplicated regions across species, with a focus on distinguishing subgenomes.
3. Ancestral Karyotyping（祖先染色体组）
Reconstructs protochromosomes and traces common chromosomal rearrangements to understand chromosome evolution.

## 1 配置WGDI分析环境

WGDI最好是在x86_64 linux环境下部署，Windows下因为mafft无法正常安装导致需要重新指定alignment软件，macOS下因为现在新出的mac都是M系列芯片，都属于ARM64架构，导致paml无法正常安装，综上就在服务器里或者wsl里跑，不要折磨自己

首先创建conda环境，或者直接在vscode里创建，我以python 3.12为例

```shell
conda create -n wgdi python=3.12
conda activate wgdi
```

安装wgdi及其依赖，当然你用pip安装也可以

```shell
conda install -c bioconda wgdi
```

新建一个wgdi工作目录，用于存储所需脚本

```shell
mkdir wgdi
cd wgdi
```

然后准备好自己要分析的基因组，包括以下文件：

- 基因组.fasta
- CDS.fasta
- 蛋白.fasta
- 基因结构注释.gff

需要注意的是CDS和cDNA之间的区别，CDS只包括（含）起始密码子到终止密码子的序列，而cDNA还包括5' UTR和3' UTR

## 2 数据预处理

数据预处理是WGDI进行共线性分析最漫长的步骤

因为上述的基因组文件并不能直接用于分析

而需要进行一些转换才能得到WGDI所需的文件形式

WGDI需要三种输入文件：
>1. BLAST的-outfmt 6输出的文件
>2. 基因的位置信息，以tab分隔，分别为chr,id,start,end,strand,order,old_id，并非是我们熟知的那个gff格式
>3. 染色体长度信息和染色体上的基因个数，格式为chr,length,gene number

同时，对于每个基因只需要一个转录本，对于有多条转录本的基因，通常使用最长的转录本代表此基因

幸好作者已经提供了一个脚本[generate_conf.py](https://github.com/xuzhougeng/myscripts/blob/master/comparative/generate_conf.py)，我们有救了

使用vim编辑器复制粘贴，或者直接wget

```shell
wget -c -t 0 https://github.com/xuzhougeng/myscripts/blob/master/comparative/generate_conf.py
```

如果是在Windows中出现wget二义性错误，参考我[这篇文章](https://yuntaobioinformatics.wordpress.com/2024/12/23/%e8%a7%a3%e5%86%b3powershell%e4%bd%bf%e7%94%a8wget-curl%e5%87%ba%e7%8e%b0%e4%ba%8c%e4%b9%89%e6%80%a7%e9%94%99%e8%af%af/)

下面是该脚本，使用vim的话，直接在目录下打开vim，复制，然后:w generate_conf.py保存

```python
#!/usr/bin/env python3

# GFF must have CDS feature
# GFF must have ID and Parent in column 9

import re
import argparse
from collections import defaultdict
from collections import OrderedDict

def get_parser():
    """Get options"""
    parser = argparse.ArgumentParser()
    parser.add_argument('fasta',
                        help="fasta file name")
    parser.add_argument('gff',
                        help="gff file name")
    parser.add_argument('-p','--prefix', type=str, default="output",
                        help="prefix for ouput ")

    return parser



# get the fasta  len
def get_fasta_len(fasta):
    fasta_dict = OrderedDict()
    handle = open(fasta, "r")
    active_sequence_name = ""
    for line in handle:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"): 
            active_sequence_name = line[1:]
            active_sequence_name = active_sequence_name.split(" ")[0]
        if active_sequence_name not in fasta_dict:
            fasta_dict[active_sequence_name] = 0
            continue
        sequence = line
        fasta_dict[active_sequence_name] += len(sequence)
    handle.close()
    return fasta_dict

# parse the gff 
def parse_gff(gff):

    gene_dict = OrderedDict()
    tx_pos_dict = defaultdict(list)
    CDS_dict = defaultdict(list)

    handle = open(gff, "r")

    for line in handle:
        if line.startswith("#"):
            continue
        content = line.split("\t")
        if len(content) <= 8:
            continue
        #print(content)
        if content[2] == "transcript" or content[2] == "mRNA":

            # use regual expression to extract the gene ID
            # match the pattern ID=xxxxx; or ID=xxxxx

            tx_id = re.search(r'ID=(.*?)[;\n]',content[8]).group(1)
            tx_parent = re.search(r'Parent=(.*?)[;\n]',content[8])
            if tx_parent is None:
                tx_parent = tx_id
            else:
                tx_parent = tx_parent.group(1)
    
            tx_parent = tx_parent.strip() # remove the 'r' and '\n'
            
            # if the parent of transcript is not in the gene_dict, create it rather than append
            if tx_parent in gene_dict:
                gene_dict[tx_parent].append(tx_id)
            else:
                gene_dict[tx_parent] = [tx_id]
            tx_pos_dict[tx_id] = [content[0],content[3], content[4], content[6] ]
        # GFF must have CDS feature
        if content[2] == 'CDS':
            width = int(content[4]) - int(content[3])
            CDS_parent = re.search(r'Parent=(.*?)[;\n]',content[8]).group(1)
            CDS_parent = CDS_parent.strip() # strip the '\r' and '\n'
            CDS_dict[CDS_parent].append(width)
    handle.close()
    return [gene_dict, tx_pos_dict, CDS_dict]

if __name__ == "__main__":
    parser = get_parser()
    args = parser.parse_args()
    fa_dict = get_fasta_len( args.fasta)
    gene_dict, tx_pos_dict, CDS_dict= parse_gff( args.gff )
    gene_count = {}

    # outfile
    len_file = open(args.prefix + ".len", "w")
    gff_file = open(args.prefix + ".gff", "w")

    for gene, txs in gene_dict.items():
        tmp = 0
        for tx in txs:
            tx_len = sum(CDS_dict[tx])
            if tx_len > tmp:
                lst_tx = tx
                tmp = tx_len
        tx_chrom = tx_pos_dict[lst_tx][0]
        if tx_chrom not in gene_count:
            gene_count[tx_chrom] = 1
        else:
            gene_count[tx_chrom] += 1
        tx_start = tx_pos_dict[lst_tx][1]
        tx_end   = tx_pos_dict[lst_tx][2]
        tx_strand = tx_pos_dict[lst_tx][3]

        print("{chrom}\t{gene}\t{start}\t{end}\t{strand}\t{order}\t{tx}".format(\
                chrom=tx_chrom,gene=gene,start=tx_start,end=tx_end,strand=tx_strand,order=gene_count[tx_chrom],tx=lst_tx), file=gff_file )

    for chrom,lens in fa_dict.items():
        print("{chrom}\t{lens}\t{count}".format(\
                chrom=chrom,lens=lens,count=gene_count.get(chrom,0)), file=len_file)
    len_file.close()
    gff_file.close()
```

下载好或者复制好generate_conf.py后，使用方式如下

```shell
python generate_conf.py -p palm /mnt/e/phoenix_dactylifera_genome/palm_genome.fa /mnt/e/phoenix_dactylifera_genome/palm.gff
# 输出文件为palm.gff和palm.len，在当前目录下，耐心等会儿
```

我这里用的是NCBI下载的油棕基因组，先不吐槽为什么只有NCBI特色的gbff格式不给常用的gff了，NCBI上基因组的基因序列ID大都命名的又臭又长，恨不得把所有信息都塞到fasta id里

用head palm.gff，发现格式有点奇怪
```
NC_052392.1     gene-LOC103714474       125152  133752  +       1       rna-XM_008801737.4
NC_052392.1     gene-LOC103714500       133778  137710  -       2       rna-XM_017844588.3
NC_052392.1     gene-LOC103714499       138535  140896  -       3       rna-XM_017844587.3
NC_052392.1     gene-LOC103696426       280954  298392  +       4       rna-XM_008778048.4
NC_052392.1     gene-LOC103696429       300565  302489  +       5       rna-XM_008778050.4
NC_052392.1     gene-LOC120104173       366690  377636  -       6       rna-XM_039114789.1
NC_052392.1     gene-LOC120104175       382497  384833  -       7       rna-XM_039114797.1
NC_052392.1     gene-LOC120104176       412404  415094  -       8       rna-XM_039114803.1
NC_052392.1     gene-LOC120106207       434342  438594  -       9       rna-XM_039119151.1
NC_052392.1     gene-LOC120104182       446167  448209  -       10      rna-XM_039114811.1
```

于是需要用sed删除奇怪的前缀“NC_052392.1     gene-”和“rna-”

```shell
sed -i -e 's/NC_052392.1	gene-//' -e 's/rna-//' palm.gff
```

再次head palm.gff，这下强迫症舒服了

```
LOC103714474    125152  133752  +       1       XM_008801737.4
LOC103714500    133778  137710  -       2       XM_017844588.3
LOC103714499    138535  140896  -       3       XM_017844587.3
LOC103696426    280954  298392  +       4       XM_008778048.4
LOC103696429    300565  302489  +       5       XM_008778050.4
LOC120104173    366690  377636  -       6       XM_039114789.1
LOC120104175    382497  384833  -       7       XM_039114797.1
LOC120104176    412404  415094  -       8       XM_039114803.1
LOC120106207    434342  438594  -       9       XM_039119151.1
LOC120104182    446167  448209  -       10      XM_039114811.1
```

另外，cds.fa里基因的命名特别长，如下所示，也需要处理一下，保证gff len以及后面的blastp结果里，基因id是一致的

处理方法自己写个shell/python/R脚本啥的吧，简单地很，这里地方太小了，我写不下

```
>lcl|NC_052392.1_cds_XP_017700082.1_1 [gene=LOC103714474] [db_xref=GeneID:103714474] [protein=glucose-1-phosphate adenylyltransferase large subunit 4, chloroplastic/amyloplastic-like] [protein_id=XP_017700082.1] [location=join(125829..126080,126400..126528,127107..127280,127383..127475,127784..127867,127962..128017,128534..128627,128792..128904,129354..129426,129616..129696,129797..129883,129969..130073,130649..130755,130835..130895,133206..133307)] [gbkey=CDS]
ATGGCGGTGGCGCTGCCCATGCACACTCTGGCGCTCGGGTCGGGCGGCGGCGATCACTGCTTGCCGAGCCTCCCTCGGCT
TCGGGGCAGAGATCTGGGGTCTGCTGGGTTGTCATCCTTCCTTGGACGGAGAGTTTGCTGCCCGGGAGGGTCGCCCAGGA
ATCATGGAGCTGGTCGACTCAGAAGGGACCTTGGCAAGCCGCCCGCTGCCGCCATCAACTCTGTTCTTGCTGACGTTGCC
AAAAATTTCAAGGCACAGCCCCTTGAGGCACCGGCGTTCGAAAGGCCGTTAGCTGATCCGGGAACCGTTGCCTCAATCAT
ATTAGGTGGAGGAGCCGGAACTCGACTTTTTCCTCTCACTCGAAGAAGGGCCAAACCAGCTGTACCGATTGGAGGTTGTT
ACAGGCTTATTGATGTTCCAATGAGCAATTGCATCAACAGTGGAATAAATAAAATTTATGTTCTAACTCAATTTAACTCT
CAGTCTCTGAATCGCCATCTTGCTCGGACATATAACTTGGGTAATGGCGTAAACTTTGGTGATGGATTTGTAGAGGTACT
AGCAGCAACACAAACACCTGGAGAATTTGGGAAGAGGTGGTTTCAGGGAACAGCAGATGCTGTTAGGCAATTCGGATGGC
TATTTGAGGATGCCAAACTTAGACATATAGAGAACATGCTGATTTTGTCTGGCGATCATCTGTATCGAATGAACTACATG
GACTTCGTGCAGAAACACATTAACTCTGGGGCTGATATATCCGTTTCTTGTGTCCCAGTGGATGACAGTCGTGCTTCGGA
TTTTGGATTGATGAAGGTTGACAAAATGGGGCGCATCCATCAGTTTCTTGAAAAGCCAAAGGGTGAAAATCTGAGGACTA
TGCAAGTGGACACAACAGTCTTAGGATTGCCTCCTGAAGACGCAAAAATGCATAAATACATCGCATCAATGGGAATATAC
GTGTTCAAGACAGATGTTCTTCTGAAGCTATTAAGATGGCGCTATCCTTCTGCTAATGATTTTGGCTCCGAGATCATACC
AATGGCAGCTAAAGACTACAATGTGCAGGCATATCTATTCGATGGATACTGGGAGGATATTGGAACTATCAAATCATTTT
TTGATGCAAATTTGGCTCTCACCGATCAGCCTCCCAAGTTTCATTTTCATGATCCTAAGAAGCCAATCTTTACATCACCA
CGGTTTTTACCACCAACCAAGATAGAAAAATGCAGGGTTGTAGACTCTATAATTTCACATGGATGCTTCTTGACACAATG
CAGTGTCAAACATTCCATTGTTGGTGTTCGCTCAAGATTAGAGTATGGGGCAGAGCTTAAGGATACCATGATGATGGGTG
CAGACTATTATCAGACTGAGGCAGAGAGAGCATCTTTCTTGGCTGAGGGGAAGGTTCCAGTTGGTGTTGGAGAGAACACT
AGGATTAGGAACTGCATCATCGACAAAAATGCCAGGATCGGAAAGGATGTGATCATCGCAAATGCAGATAATGTGGAGGA
AGCTGATAAACCATCAGATGGTTTCTACATTCGTTCTGGAATCACTGTGGTGCTGAAGAATTCTGTTATTCCAGATGGGA
CCATTATCTAG
```

有时候NCBI下载的基因组注释文件是gbff格式，需要转换为gff

如果安装bioperl报错，需要另起一个conda环境，略

```bash
# 创建个新的conda环境，为了防止古董perl出现奇怪的兼容性问题，用老一点版本的python
conda create -n bioperl python=3.9
conda activate bioperl
conda install -c bioconda perl-bioperl
# 流传很久的神秘上古perl脚本
wget https://fastapi.metacpan.org/source/CJFIELDS/BioPerl-1.7.8/bin/bp_genbank2gff3 -O bp_genbank2gff3.pl
perl bp_genbank2gff3.pl your_genbank_file.gbff
# your_genbank_file为待处理的gbff格式文件
```

通过blastp或者DIAMOND进行蛋白之间的相互比对，输出格式为-outfmt 6

```shell
# conda install bioconda::blast
makeblastdb -in example.protein.fasta -dbtype prot
blastp -num_threads 16 -db example.protein.fasta -query example.protein.fasta -outfmt 6 -evalue 1e-5 -num_alignments 20  -out example.blastp.txt
```

## 3 绘制点阵图

WGDI的使用逻辑是：创建配置文件，修改配置文件，运行程序

### Step1 创建配置文件

WGDI提供了创建配置文件的方法

```bash
# 创建配置文件
wgdi -d \? > example_dot.conf
# 编辑配置文件
vim example_dot.conf
```

打开看到配置文件应该如下所示（这并不是真正的参数，只是创建时的默认参数）

```bash
[dotplot]
blast = blast file  # Step2最后获得的blast outfmt 6文件
gff1 =  gff1 file   # 如果是种内共线性，那么gff1和gff2是同一个文件，下面类似
gff2 =  gff2 file
lens1 = lens1 file
lens2 = lens2 file
genome1_name =  Genome1 name
genome2_name =  Genome2 name
multiple  = 1   # 最好的同源基因数, 用输出结果中会用红点表示
score = 100     # blast输出的score 过滤 
evalue = 1e-5   # blast输出的evalue 过滤 
repeat_number = 20  # genome2相对于genome1的最多同源基因数
position = order
blast_reverse = false
ancestor_left = none
ancestor_top = none
markersize = 0.5  # 点的大小
figsize = 10,10   # 图片大小
savefig = savefile(.png,.pdf)
```

### Step2 修改配置文件

下面是我自己用的一个种内共线性分析的配置

反正不行就调参，都已经调包了，再多调调参也无所谓

人人都笑调参侠，人人都是调参侠

```bash
blast  =  example.blastp.txt
gff1  =  example.gff
gff2  =  example.gff
lens1  =  example.len
lens2  =  example.len
genome1_name  =  example
genome2_name  =  example
multiple  =  1
score  =  100
evalue  =  1e-5
repeat_number  =  5
position  =  order
blast_reverse  =  false
ancestor_left  =  none
ancestor_top  =  none
markersize  =  0.5
figsize  =  10,10
savefig  =  example.dot.pdf
```

### Step3 运行程序

确认配置修改完成后，运行

```bash
wgdi -d example_dot.conf
```

目录下就可以找到输出的文件了，文件名为savefig设置的值

非常不建议保存pdf，因为这使我的电脑卡死，看似只有3MB不大，但是打开来要渲染半天，一动就未响应

另外要注意的是，gff、len和blast文件里，基因ID必须格式一致，完全匹配，否则会报错

这一点十分重要，在前期准备文件时敬请注意

### 点阵图的解读？

![alt text](image.png)

这是用拟南芥做的一个共线性点阵图（种内）

很明显看到一堆红红蓝蓝还有灰灰（？）的点

#### 颜色说明

- 红色点：表示基因组中一个基因的最优同源（即相似性最高的同源基因）。
- 蓝色点：表示基因组中一个基因的次优同源（相似性次高的同源基因）。
- 灰色点：表示其他非最佳或次优的同源关系，或未能检测到显著同源关系的区域。

当然还有一条非常显眼的对角线，这不是你的屏幕被划了一刀

聪明的人应该已经猜到了：图中对角线出现的片段是自身比对

这很好，不过还是猜错了，因为WGDI已经过滤掉自身比对自身的结果了

## 4 共线性分析

WGDI开发了-icl模块（Improved version of ColinearScan）进行共线性分析，逻辑和绘制点阵图一致

### Step1 创建配置文件

```bash
wgdi -icl \? >> example_col.conf
# 稍微注意下，conf文件别重名了，否则会覆盖，如果要返工就麻烦大了
```

### Step2 修改配置文件

```
[collinearity]
gff1 = example.gff
gff2 = example.gff
lens1 = example.len
lens2 = example.len
blast = example.blastp.txt
blast_reverse = false
multiple  = 1
process = 8
evalue = 1e-5
score = 100
grading = 50,40,25
mg = 40,40
pvalue = 0.2
repeat_number = 10
positon = order
savefile = example.collinearity.txt
```

evalue和score是在blast文件中筛选同源基因对的

multiple = 1表示确定一个最佳比对基因

repeat_number表示取多少个高同源的比对结果

grading是构建共线性块的一个参数，也就是点阵图上的红、蓝、灰点分别检索下一个同源基因对的范围，详细看论文

### Step3 运行共线性分析

```bash
wgdi -icl example_col.conj
# 运行结束后得到example.collinearity.txt
# 以#起始的行记录共线性区域的Metadata
```

然后对#开头的行进行一下统计分析

```bash
# 提取#开头的结果行
grep '^#' Ldavi.collinearity.txt > Ldavi.collinearity.results.txt
# 得到的结果用shell或者excel处理下就能得到统计表格了，此处略
# 展示统计学的魅力时刻到了
# 不会也别急，接着向下看
```

## 5 根据共线性结果计算ks

一样的逻辑，这里需要用到一开始准备的cds和蛋白文件了

### Step1 创建配置文件

```bash
wgdi -ks \? >> example_ks.conj
```

### Step2 修改配置文件

```bash
[ks]
cds_file =   example.cds.fa
pep_file =   example.pep.fa
align_software = muscle
pairs_file = example.collinearity.txt
ks_file = example.ks
```

### Step3 运行ks计算

```bash
wgdi -ks example_ks.conj
```

### Step4 数据整合

```bash
wgdi -bi ? >> example_bi.conf
```

编辑配置文件

```
[blockinfo]
blast = example.blastp.txt
gff1 =  example.gff
gff2 =  example.gff
lens1 = example.len
lens2 = example.len
collinearity = example.collinearity.txt
score = 100
evalue = 1e-5
repeat_number = 20
position = order
ks = example.ks
ks_col = ks_NG86
savefile = example_block_information.csv
```

运行

```bash
wgdi -bi example_bi.conj
```

输出文件为csv格式，可以直接用excel打开

>1. id：共线性结果的唯一标识ID
>2. chr1, start1, end1：参考基因组（点阵图左侧）的共线性范围
>3. chr2, start2, end2：参考基因组（点阵图上侧）的共线性范围
>4. pvalue：共线性结果评估，一般认为p < 0.01更合理
>5. length：共线性片段长度
>6. ks_median：共线性片段上所有基因对的ks的中位数（主要用于评判ks分布）
>7. ks_average：~~~的平均值
>8. homo1, homo2, homo3, ..., homoX：与multiple参数有关，multiple = N，就有N列
>主要规则是：基因对如果在点阵图上为红色，赋值1；蓝色，赋值0；灰色，赋值-1。该值越接近1，说明共线性的点大部分为红色，这是共线性片段筛选的依据。
>9. block1, block2：分别为共线性片段上基因order的位置
>10. ks：共线性片段的基因对的ks值
>11. density1, density2：共线性片段的基因分布密集程度，值越小说明基因分布越稀疏

这个表格就是后续分析的万恶之源了，下次再写，今天下班了

参考资料：
https://www.jianshu.com/p/a50548e81ac0
https://wgdi.readthedocs.io/en/latest/index.html