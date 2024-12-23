# 用WGDI分析比较基因组
2024.12.23
朱云涛
环境：(wsl) Ubuntu 24.04 x64, python 3.12, AMD R7-5700X3D 8C16T, 32G RAM

## 0 前言 关于WGDI

WGDI是一篇发表在Molecular Plant上的植物比较基因组工具，详情见文章

[WGDI: A user-friendly toolkit for evolutionary  analyses of whole-genome duplications and  ancestral karyotypes](10.1016/j.molp.2022.10.018)

以下是机翻的Abstract

> 在地球上大多数主要的生物谱系中都发现了全基因组重复（WGDs）和随后的核型变化的证据。为了阐明基因组分析中复杂的基因共线性多层模式，需要方便准确的工具包。为了满足这一需求，我们开发了WGDI（全基因组复制集成分析），这是一种基于Python的命令行工具，有助于全面分析递归多倍体事件和跨物种基因组比对。WGDI支持三个主要工作流程（多倍体推断，基因组同源性的层次推断和祖先染色体核型分析），可以基于高质量的染色体水平基因组改进WGD的检测和WGD相关事件的表征。值得注意的是，它可以提取完整的同线性块，并有助于重建详细的核型进化。该工具包可在GitHub免费获得(https://github.com/SunPengChuan/wgdi)。作为其应用的一个例子，WGDI令人信服地阐明了WGDs后Aquilegia coerulea和Vitis vinifera的核型进化，并拒绝了Aquilegia作为核心双子叶植物异源多倍体起源的亲本谱系的假设。

使用方法参见[WGDI的Github主页](https://github.com/SunPengChuan/wgdi)

WGDI的三个关键特性：
1. Polyploid Inference（多倍体推断）
Identifies and confirms polyploid events with high accuracy.
2. Genomic Homology Inference（基因组同源性推断）
Traces the evolutionary history of duplicated regions across species, with a focus on distinguishing subgenomes.
3. Ancestral Karyotyping（祖先染色体组）
Reconstructs protochromosomes and traces common chromosomal rearrangements to understand chromosome evolution.

## 1 配置WGDI分析环境

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

### 2 数据预处理

数据预处理是WGDI进行共线性分析最漫长的步骤

因为上述的基因组文件并不能直接用于分析

而需要进行一些转换才能得到WGDI所需的文件形式

WGDI需要三种输入文件：
>1. BLAST的-outfmt 6输出的文件
>2. 基因的位置信息，以tab分隔，分别为chr,id,start,end,strand,order,old_id，并非是我们熟知的那个gff格式
>3. 染色体长度信息和染色体上的基因个数，格式为chr,length,gene number

同时，对于每个基因只需要一个转录本，对于有多条转录本的基因，通常使用最长的转录本代表此基因

幸好作者已经提供了一个脚本[generate_conf.py](https://github.com/xuzhougeng/myscripts/blob/master/comparative/generate_conf.py)

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

再次head palm.gff，这下舒服了

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

另外，cds.fa里基因的命名特别长，如下所示
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

