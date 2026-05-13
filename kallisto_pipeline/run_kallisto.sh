#!/bin/bash

# 欢迎界面
echo "============================================================"
echo "Kallisto RNA-Seq 批量定量分析软件 v2.5 2026.5.13"
echo "作者：朱云涛（北京林业大学）"
echo "------------------------------------------------------------"
echo "用法：sh $0 <INDEX> <INDIR> [THREADS]"
echo "INDEX 可选预设: Lsor2021（索邦百合华大版本）/Lsor2025（索邦百合2025版本）/Ldavidii（兰州百合）/Lregale（岷江百合）"
echo "    如果 INDEX 不在预设中，则使用该路径的自定义索引"
echo "INDIR: 包含多个样品子文件夹的目录，每个子文件夹内包含成对的测序文件*.fq.gz/*.fastq.gz"
echo "THREADS: 可选参数，指定线程数，默认为20，太高了也并不会显著增加并行（有I/O瓶颈），可以不改"
echo "用例1（使用预设）：sh $0 Ldavidii /RNA-seq/30个兰州百合转录组"
echo "用例2（使用自定义）：sh $0 /path/to/kallisto_index.idx /RNA-seq/60个转录组"
echo "用例3（指定线程数）：sh $0 Lsor /RNA-seq/60个索邦百合转录组 36"
echo "------------------------------------------------------------"
echo "目前已有预设索引说明："
echo "1、Lsor2021：索邦百合2021年华大基因组装版本，之前的Lsor，现在已经被Lsor2025取代"
echo "2、Lsor2025：索邦百合2025年10月组装的参考转录本，目前索邦用这个"
echo "3、Ldavidii：兰州百合基因组mRNA构建的索引"
echo "4、Lregale：岷江百合基因组mRNA构建的索引"
echo "5、LsorRoot2026：索邦百合2026年1月组装的根尖参考转录本"
echo "------------------------------------------------------------"
echo "V2.5新增功能："
echo "新增加了预设索引LilySecretKiss，对应2026年5月组装的亚洲百合‘秘密之吻’的参考转录本"
echo "V2.4新增功能："
echo "新增加了预设索引LsorRoot2026，对应2026年1月组装的索邦百合根尖参考转录本"
echo "V2.3新增功能："
echo "新增加了生成html交互式可视化报告功能，实验性，开发中"
echo "V2.2新增功能："
echo "新增加了预设索引Lsor2025，对应2025年10月组装的索邦百合茎尖参考转录本"
echo "V2.1新增功能："
echo "1、自动识别fq.gz和fastq.gz格式的测序文件"
echo "2、整合所有样品的TPM和Count矩阵"
echo "3、生成详细的分析摘要报告"
echo "4、修正了样品文件名中含有空格时跳过样品的bug"
echo "------------------------------------------------------------"
echo "其他说明："
echo "1、结果将保存到与INDIR同级的kallisto_result目录下，结构如下："
echo "   - details/: 各样品的详细kallisto结果"
echo "   - summary/: TPM矩阵、Count矩阵和分析摘要"
echo "2、使用前请确保已经切换到RNA-seq或其他配置好的conda环境"
echo "3、线程数可通过第三个参数自定义，默认为20"
echo "4、HTML交互式可视化需要Python依赖: numpy pandas scipy scikit-learn"
echo "   如缺少依赖，脚本会生成基础HTML报告（不含交互图表）"

# 判断参数个数
if [ $# -lt 2 ]; then
  echo "参数数量不对！"
  echo "用法: sh $0 <INDEX> <INDIR> [THREADS]"
  echo "INDEX 可选预设: Lsor2021/Lsor2025/Ldavidii/Lregale；或者自定义"
  echo "THREADS: 可选参数，指定线程数，默认为20"
  exit 1
fi

INDEX=$1
INDIR=$2
# 第三个参数为线程数，未提供则默认为20
if [ $# -ge 3 ] && [ -n "$3" ]; then
  THREADS=$3
  echo "使用自定义线程数: $THREADS"
else
  THREADS=20
  echo "使用默认线程数: $THREADS"
fi

# 检查输入目录是否存在
if [ ! -d "$INDIR" ]; then
  echo "错误：输入目录不存在！($INDIR)"
  exit 1
fi

# 根据 INDEX 选择预设路径
if [ "$INDEX" = "Lsor2021" ]; then
  INDEX_PATH=$HOME/Lsor/Lsor.kallisto.idx
  echo "使用预设的kallisto索引：Lsor（$INDEX_PATH）"
elif [ "$INDEX" = "Lsor2025" ]; then
  INDEX_PATH="/media/desk16/tl5024/LsorJingjianRef/LsorJingjian.kallisto.idx"
  echo "使用预设的kallisto索引：Lsor2025（$INDEX_PATH）"
elif [ "$INDEX" = "Ldavidii" ]; then
  INDEX_PATH=$HOME/Ldavidii/Ldavidii.kallisto.idx
  echo "使用预设的kallisto索引：Ldavidii（$INDEX_PATH）"
elif [ "$INDEX" = "Lregale" ]; then
  INDEX_PATH=$HOME/Lregale/Lregale.kallisto.idx
  echo "使用预设的kallisto索引：Lregale（$INDEX_PATH）"
elif [ "$INDEX" = "LsorRoot2026" ]; then
  INDEX_PATH=$HOME/LsorRoot2026/LsorRoot2026.kallisto.idx
  echo "使用预设的kallisto索引：LsorRoot2026（$INDEX_PATH）"
elif [ "$INDEX" = "LilySK" ]; then
  INDEX_PATH=$HOME/LilySecretKiss/LilySK.kallisto.idx
  echo "使用预设的kallisto索引：LilySecretKiss（$INDEX_PATH）"
else
  echo "使用自定义的kallisto索引： '$INDEX'"
  INDEX_PATH="$INDEX"
  echo "提示：预设索引可选: Lsor2021/Lsor2025/Ldavidii/Lregale/LsorRoot2026/LilySK"
fi

# 检查文件是否存在
if [ ! -f "$INDEX_PATH" ]; then
  echo "错误：指定的索引文件不存在！($INDEX_PATH)"
  exit 1
fi

# 获取输入目录的父目录
PARENT_DIR=$(dirname "$INDIR")

# 设置输出目录
OUTDIR="$PARENT_DIR/kallisto_result"
DETAILS_DIR="$OUTDIR/details"
SUMMARY_DIR="$OUTDIR/summary"

mkdir -p "$DETAILS_DIR"
mkdir -p "$SUMMARY_DIR"

# 验证线程数是否有效
if ! echo "$THREADS" | grep -q '^[0-9][0-9]*$' || [ "$THREADS" -le 0 ]; then
  echo "警告: 线程数 '$THREADS' 无效，使用默认值20"
  THREADS=20
fi

# 初始化样品信息记录
SAMPLE_INFO="$SUMMARY_DIR/sample_info.txt"
printf "序号\t样品名称\tRead1路径\tRead2路径\t分析状态\n" > "$SAMPLE_INFO"

# 获取kallisto版本信息
KALLISTO_VERSION=$(kallisto version 2>&1 | head -n 1 || echo "未知版本")

# 记录分析命令
COMMAND_USED="kallisto quant -i $INDEX_PATH -o <output_dir> -t $THREADS -b 100 <read1> <read2>"

echo "已经创建好结果目录，分析即将开始......"
echo "使用索引: $INDEX_PATH"
echo "输入目录: $INDIR"
echo "结果主目录: $OUTDIR"
echo "详细结果目录: $DETAILS_DIR"
echo "汇总结果目录: $SUMMARY_DIR"
echo "使用线程数: $THREADS"
echo "Kallisto版本: $KALLISTO_VERSION"

# 遍历输入目录下的样品文件夹
sample_count=0
success_count=0
for sample in "$INDIR"/*; do
  if [ -d "$sample" ]; then
    name=$(basename "$sample")
    
    # 跳过结果目录（以kallisto开头的文件夹）
    case "$name" in
      kallisto*)
        echo "跳过结果目录: $name"
        continue
        ;;
    esac
    
    sample_count=$((sample_count + 1))
    
    # 尝试找到配对的测序文件（支持fq.gz和fastq.gz）
    fq1=""
    fq2=""
    
    # 优先查找_R1和_R2格式（最常见）
    fq1_candidates=$(ls "$sample"/*_R1*.f*q.gz 2>/dev/null)
    fq2_candidates=$(ls "$sample"/*_R2*.f*q.gz 2>/dev/null)
    
    # 如果没有找到_R1/_R2，则查找_1/_2格式
    if [ -z "$fq1_candidates" ] || [ -z "$fq2_candidates" ]; then
      fq1_candidates=$(ls "$sample"/*_1.fq.gz 2>/dev/null)
      fq2_candidates=$(ls "$sample"/*_2.fq.gz 2>/dev/null)
    fi
    
    # 如果还没找到fq.gz，则查找fastq.gz格式
    if [ -z "$fq1_candidates" ] || [ -z "$fq2_candidates" ]; then
      fq1_candidates=$(ls "$sample"/*_1.fastq.gz 2>/dev/null)
      fq2_candidates=$(ls "$sample"/*_2.fastq.gz 2>/dev/null)
    fi
    
    # 取第一个匹配的文件（安全处理含空格的路径）
    if [ -n "$fq1_candidates" ]; then
      fq1=$(echo "$fq1_candidates" | head -n1)
    else
      fq1=""
    fi
    
    if [ -n "$fq2_candidates" ]; then
      fq2=$(echo "$fq2_candidates" | head -n1)
    else
      fq2=""
    fi
    
    if [ -n "$fq1" ] && [ -n "$fq2" ] && [ -f "$fq1" ] && [ -f "$fq2" ]; then
      echo "正在处理样品 $name ..."
      echo "  Read1: $fq1"
      echo "  Read2: $fq2"
      
      # 运行kallisto
      if kallisto quant -i "$INDEX_PATH" -o "$DETAILS_DIR/$name" -t "$THREADS" -b 100 "$fq1" "$fq2"; then
        printf "%s\t%s\t%s\t%s\t成功\n" "$sample_count" "$name" "$fq1" "$fq2" >> "$SAMPLE_INFO"
        success_count=$((success_count + 1))
        echo "样品 $name 分析完成"
      else
        printf "%s\t%s\t%s\t%s\t失败\n" "$sample_count" "$name" "$fq1" "$fq2" >> "$SAMPLE_INFO"
        echo "警告: 样品 $name 分析失败"
      fi
    else
      echo "警告: 样品 $name 缺少配对的测序文件，跳过"
      printf "%s\t%s\t未找到\t未找到\t跳过\n" "$sample_count" "$name" >> "$SAMPLE_INFO"
    fi
  fi
done

echo "============================================================"
echo "Kallisto定量分析完成！"
echo "总样品数: $sample_count"
echo "成功分析: $success_count"
echo "正在生成汇总报告..."

# 生成TPM和Count矩阵
echo "正在整合TPM和Count矩阵..."
export DETAILS_DIR
export SUMMARY_DIR
python3 - << 'EOF'
import os
import pandas as pd
import glob

# 设置路径
details_dir = os.environ['DETAILS_DIR']
summary_dir = os.environ['SUMMARY_DIR']

# 获取所有成功的样品目录
sample_dirs = glob.glob(os.path.join(details_dir, '*'))
# 过滤掉结果目录和非目录文件
sample_dirs = [d for d in sample_dirs if os.path.isdir(d) and 
               not os.path.basename(d).startswith('multiQC') and
               not os.path.basename(d).startswith('kallisto')]

if not sample_dirs:
    print("未找到任何样品结果目录")
    exit(1)

# 读取第一个样品的abundance.tsv获取基因列表
first_sample = sample_dirs[0]
abundance_file = os.path.join(first_sample, 'abundance.tsv')

if not os.path.exists(abundance_file):
    print(f"未找到abundance.tsv文件: {abundance_file}")
    exit(1)

# 读取基因信息
df_first = pd.read_csv(abundance_file, sep='\t')

# 初始化TPM和Count矩阵，只保留target_id作为第一列
tpm_matrix = df_first[['target_id']].copy()
count_matrix = df_first[['target_id']].copy()

# 重命名第一列为gene_id
tpm_matrix.rename(columns={'target_id': 'gene_id'}, inplace=True)
count_matrix.rename(columns={'target_id': 'gene_id'}, inplace=True)

# 按样品名称排序（保证相同组的样品在一起）
sample_names = [os.path.basename(d) for d in sample_dirs]
sample_names.sort()

# 读取每个样品的数据
for sample_name in sample_names:
    sample_dir = os.path.join(details_dir, sample_name)
    abundance_file = os.path.join(sample_dir, 'abundance.tsv')
    
    if os.path.exists(abundance_file):
        df_sample = pd.read_csv(abundance_file, sep='\t')
        tpm_matrix[sample_name] = df_sample['tpm']
        count_matrix[sample_name] = df_sample['est_counts']
    else:
        print(f"警告: 未找到样品 {sample_name} 的abundance.tsv文件")

# 保存矩阵
tpm_output = os.path.join(summary_dir, 'TPM.csv')
count_output = os.path.join(summary_dir, 'estCounts.csv')

tpm_matrix.to_csv(tpm_output, index=False)
count_matrix.to_csv(count_output, index=False)

print(f"TPM矩阵已保存到: {tpm_output}")
print(f"Count矩阵已保存到: {count_output}")
print(f"共处理 {len(sample_names)} 个样品")
EOF

# 生成分析摘要
echo "正在生成分析摘要..."
cat > "$SUMMARY_DIR/summary.md" << EOF
# Kallisto RNA-Seq 批量定量分析v2.3摘要报告

by 朱云涛 at $(date)
欢迎访问[我的博客](https://yuntaobioinformatics.wordpress.com/)
Welcome to visit [My Blog](https://yuntaobioinformatics.wordpress.com/)

## 1 分析基本信息

- **Kallisto版本**: $KALLISTO_VERSION
- **完成时间**: $(date)
- **分析命令**: \`$COMMAND_USED\`
- **使用的索引**: $INDEX_PATH
- **输入目录**: $INDIR
- **输出目录**: $OUTDIR

## 2 分析样品统计

- **总样品数**: $sample_count
- **成功分析**: $success_count
- **失败样品**: $((sample_count - success_count))

## 3 分析结果详情

EOF

# 添加样品信息表格到摘要
echo "| 序号 | 样品名称 | Read1路径 | Read2路径 | 分析状态 |" >> "$SUMMARY_DIR/summary.md"
echo "|------|----------|-----------|-----------|----------|" >> "$SUMMARY_DIR/summary.md"
# 使用sed处理制表符分隔的文件
tail -n +2 "$SAMPLE_INFO" | sed 's/\t/ | /g' | sed 's/^/| /' | sed 's/$/ |/' >> "$SUMMARY_DIR/summary.md"

# 生成比对统计摘要
echo "" >> "$SUMMARY_DIR/summary.md"
echo "## 4 Kallisto比对统计" >> "$SUMMARY_DIR/summary.md"
echo "" >> "$SUMMARY_DIR/summary.md"
echo "| 样品名称 | 处理的reads数 | 比对上的reads数 | 比对率(%) |" >> "$SUMMARY_DIR/summary.md"
echo "|----------|---------------|-----------------|-----------|" >> "$SUMMARY_DIR/summary.md"

# 从每个样品的run_info.json提取统计信息
for sample in "$INDIR"/*; do
  if [ -d "$sample" ]; then
    name=$(basename "$sample")
    run_info="$DETAILS_DIR/$name/run_info.json"
    if [ -f "$run_info" ]; then
      # 使用python解析JSON
      python3 - << EOF >> "$SUMMARY_DIR/summary.md"
import json
import os

run_info_file = "$run_info"
sample_name = "$name"

try:
    with open(run_info_file, 'r') as f:
        data = json.load(f)
    
    n_processed = data.get('n_processed', 0)
    n_pseudoaligned = data.get('n_pseudoaligned', 0)
    
    if n_processed > 0:
        alignment_rate = (n_pseudoaligned / n_processed) * 100
    else:
        alignment_rate = 0
    
    print(f"| {sample_name} | {n_processed:,} | {n_pseudoaligned:,} | {alignment_rate:.2f} |")
    
except Exception as e:
    print(f"| {sample_name} | - | - | - |")
EOF
    fi
  fi
done

# 生成HTML摘要报告

echo "正在生成HTML摘要报告..."
echo "------------------------------------------------------------"
echo "说明：HTML报告包含交互式可视化功能，需要以下Python依赖："
echo "  - numpy"
echo "  - pandas"
echo "  - scipy"
echo "  - scikit-learn"
echo ""
echo "如果缺少依赖，请使用以下命令安装："
echo "  conda install numpy pandas scipy scikit-learn"
echo "  或"
echo "  pip install numpy pandas scipy scikit-learn"
echo "------------------------------------------------------------"

export SUMMARY_DIR
export sample_count
export success_count
export KALLISTO_VERSION
export INDEX_PATH
export INDIR
export OUTDIR
python3 - << 'HTMLEOF'
import os
import json
from datetime import datetime

# 检查并导入可视化所需的库
try:
    import numpy as np
    import pandas as pd
    from scipy.cluster.hierarchy import linkage, dendrogram, leaves_list
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler
    viz_available = True
    print("✓ 可视化依赖已满足，将生成交互式图表")
except ImportError as e:
    viz_available = False
    print(f"缺少可视化依赖: {e}")
    print("  将生成基础HTML报告（不包含交互式图表）")
    print("  如需完整功能，请安装: conda install numpy pandas scipy scikit-learn")

summary_dir = os.environ['SUMMARY_DIR']
sample_count = os.environ['sample_count']
success_count = os.environ['success_count']
kallisto_version = os.environ['KALLISTO_VERSION']
index_path = os.environ['INDEX_PATH']
indir = os.environ['INDIR']
outdir = os.environ['OUTDIR']

# 读取样品信息
sample_info_file = os.path.join(summary_dir, 'sample_info.txt')
samples = []
if os.path.exists(sample_info_file):
    with open(sample_info_file, 'r') as f:
        lines = f.readlines()[1:]  # 跳过表头
        for line in lines:
            parts = line.strip().split('\t')
            if len(parts) >= 5:
                samples.append({
                    'num': parts[0],
                    'name': parts[1],
                    'read1': parts[2],
                    'read2': parts[3],
                    'status': parts[4]
                })

# 读取比对统计信息
details_dir = os.path.join(os.path.dirname(summary_dir), 'details')
alignment_stats = []
for sample in samples:
    if sample['status'] == '成功':
        run_info_file = os.path.join(details_dir, sample['name'], 'run_info.json')
        if os.path.exists(run_info_file):
            with open(run_info_file, 'r') as f:
                data = json.load(f)
                n_processed = data.get('n_processed', 0)
                n_pseudoaligned = data.get('n_pseudoaligned', 0)
                alignment_rate = (n_pseudoaligned / n_processed * 100) if n_processed > 0 else 0
                alignment_stats.append({
                    'name': sample['name'],
                    'n_processed': f"{n_processed:,}",
                    'n_pseudoaligned': f"{n_pseudoaligned:,}",
                    'rate': f"{alignment_rate:.2f}"
                })

# 读取TPM矩阵用于可视化
tpm_file = os.path.join(summary_dir, 'TPM.csv')
tpm_data = None
correlation_data = None
pca_data = None

if viz_available and os.path.exists(tpm_file):
    try:
        df_tpm = pd.read_csv(tpm_file)
        df_tpm = df_tpm.set_index('gene_id')
        
        # 1. TPM分布数据（每个样品的TPM分布）
        tpm_distributions = {}
        for col in df_tpm.columns:
            tpm_values = df_tpm[col][df_tpm[col] > 0]  # 只取表达的基因
            if len(tpm_values) > 0:
                # 计算TPM的对数分布（避免0）
                log_tpm = np.log10(tpm_values + 0.1)
                tpm_distributions[col] = log_tpm.tolist()
        
        # 2. 样品相关性矩阵（层次聚类排序）
        corr_matrix = df_tpm.corr()
        # 进行层次聚类
        linkage_matrix = linkage(corr_matrix, method='average')
        ordered_indices = leaves_list(linkage_matrix)
        ordered_samples = [corr_matrix.columns[i] for i in ordered_indices]
        # 重新排序相关性矩阵
        corr_matrix_ordered = corr_matrix.loc[ordered_samples, ordered_samples]
        
        correlation_data = {
            'samples': ordered_samples,
            'matrix': corr_matrix_ordered.values.tolist()
        }
        
        # 3. PCA分析（选择变异最大的前2000个基因）
        # 计算每个基因的方差
        gene_vars = df_tpm.var(axis=1)
        top_genes = gene_vars.nlargest(min(2000, len(gene_vars))).index
        df_top = df_tpm.loc[top_genes]
        
        # 标准化并进行PCA
        scaler = StandardScaler()
        df_scaled = scaler.fit_transform(df_top.T)
        
        pca = PCA(n_components=min(3, df_top.shape[1]))
        pca_result = pca.fit_transform(df_scaled)
        
        pca_data = {
            'samples': df_tpm.columns.tolist(),
            'pc1': pca_result[:, 0].tolist(),
            'pc2': pca_result[:, 1].tolist(),
            'explained_variance': [f"{var*100:.2f}" for var in pca.explained_variance_ratio_[:2]]
        }
        
        tpm_data = tpm_distributions
        
    except Exception as e:
        print(f"警告: 生成可视化数据时出错: {str(e)}")
        tpm_data = None
        correlation_data = None
        pca_data = None

# 生成HTML报告
html_content = f'''<!DOCTYPE html>
<html lang="zh-CN">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Kallisto RNA-Seq 批量定量分析报告</title>
    <script src="https://cdn.plot.ly/plotly-2.26.0.min.js"></script>
    <style>
        * {{
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }}
        body {{
            font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif;
            line-height: 1.6;
            color: #333;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            padding: 20px;
        }}
        .container {{
            max-width: 1200px;
            margin: 0 auto;
            background: white;
            border-radius: 10px;
            box-shadow: 0 10px 40px rgba(0,0,0,0.1);
            overflow: hidden;
        }}
        .header {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 40px;
            text-align: center;
        }}
        .header h1 {{
            font-size: 2.5em;
            margin-bottom: 10px;
            text-shadow: 2px 2px 4px rgba(0,0,0,0.2);
        }}
        .header p {{
            font-size: 1.1em;
            opacity: 0.9;
        }}
        .content {{
            padding: 40px;
        }}
        .section {{
            margin-bottom: 40px;
        }}
        .section h2 {{
            color: #667eea;
            font-size: 1.8em;
            margin-bottom: 20px;
            padding-bottom: 10px;
            border-bottom: 3px solid #667eea;
        }}
        .info-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
            gap: 20px;
            margin-bottom: 30px;
        }}
        .info-card {{
            background: #f8f9fa;
            padding: 20px;
            border-radius: 8px;
            border-left: 4px solid #667eea;
        }}
        .info-card h3 {{
            color: #667eea;
            font-size: 0.9em;
            text-transform: uppercase;
            margin-bottom: 10px;
        }}
        .info-card p {{
            font-size: 1.5em;
            font-weight: bold;
            color: #333;
        }}
        .stats-summary {{
            display: flex;
            justify-content: space-around;
            margin: 30px 0;
            flex-wrap: wrap;
        }}
        .stat-box {{
            text-align: center;
            padding: 20px;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            border-radius: 10px;
            min-width: 150px;
            margin: 10px;
            box-shadow: 0 4px 15px rgba(102, 126, 234, 0.3);
        }}
        .stat-box h3 {{
            font-size: 2.5em;
            margin-bottom: 5px;
        }}
        .stat-box p {{
            font-size: 0.9em;
            opacity: 0.9;
        }}
        table {{
            width: 100%;
            border-collapse: collapse;
            margin: 20px 0;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
        }}
        th {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 15px;
            text-align: left;
            font-weight: 600;
        }}
        td {{
            padding: 12px 15px;
            border-bottom: 1px solid #e0e0e0;
        }}
        tr:hover {{
            background-color: #f5f5f5;
        }}
        .status-success {{
            color: #28a745;
            font-weight: bold;
        }}
        .status-failed {{
            color: #dc3545;
            font-weight: bold;
        }}
        .status-skipped {{
            color: #ffc107;
            font-weight: bold;
        }}
        .footer {{
            background: #f8f9fa;
            padding: 30px;
            text-align: center;
            color: #666;
            border-top: 1px solid #e0e0e0;
        }}
        .footer a {{
            color: #667eea;
            text-decoration: none;
        }}
        .footer a:hover {{
            text-decoration: underline;
        }}
        .path-info {{
            background: #f8f9fa;
            padding: 15px;
            border-radius: 5px;
            font-family: monospace;
            font-size: 0.9em;
            margin: 10px 0;
            word-break: break-all;
        }}
        .chart-container {{
            background: white;
            padding: 20px;
            border-radius: 8px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
            margin: 20px 0;
        }}
        .tab-container {{
            margin: 20px 0;
        }}
        .tab-buttons {{
            display: flex;
            border-bottom: 2px solid #e0e0e0;
            margin-bottom: 20px;
        }}
        .tab-button {{
            padding: 12px 24px;
            background: none;
            border: none;
            cursor: pointer;
            font-size: 1em;
            color: #666;
            transition: all 0.3s;
            border-bottom: 3px solid transparent;
        }}
        .tab-button.active {{
            color: #667eea;
            border-bottom-color: #667eea;
            font-weight: bold;
        }}
        .tab-button:hover {{
            color: #667eea;
        }}
        .tab-content {{
            display: none;
        }}
        .tab-content.active {{
            display: block;
        }}
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <h1>🧬 Kallisto RNA-Seq 批量定量分析报告</h1>
            <p>版本 2.3 | 生成时间: {datetime.now().strftime("%Y年%m月%d日 %H:%M:%S")}</p>
        </div>
        
        <div class="content">
            <!-- 统计概览 -->
            <div class="section">
                <h2>📊 分析统计概览</h2>
                <div class="stats-summary">
                    <div class="stat-box">
                        <h3>{sample_count}</h3>
                        <p>总样品数</p>
                    </div>
                    <div class="stat-box">
                        <h3>{success_count}</h3>
                        <p>成功分析</p>
                    </div>
                    <div class="stat-box">
                        <h3>{int(sample_count) - int(success_count)}</h3>
                        <p>失败/跳过</p>
                    </div>
                </div>
            </div>

            <!-- 基本信息 -->
            <div class="section">
                <h2>ℹ️ 分析基本信息</h2>
                <div class="info-grid">
                    <div class="info-card">
                        <h3>Kallisto 版本</h3>
                        <p>{kallisto_version}</p>
                    </div>
                    <div class="info-card">
                        <h3>使用线程数</h3>
                        <p>20</p>
                    </div>
                </div>
                <div class="info-card">
                    <h3>使用的索引</h3>
                    <div class="path-info">{index_path}</div>
                </div>
                <div class="info-card">
                    <h3>输入目录</h3>
                    <div class="path-info">{indir}</div>
                </div>
                <div class="info-card">
                    <h3>输出目录</h3>
                    <div class="path-info">{outdir}</div>
                </div>
            </div>

            <!-- 样品详情 -->
            <div class="section">
                <h2>📋 样品分析详情</h2>
                <table>
                    <thead>
                        <tr>
                            <th>序号</th>
                            <th>样品名称</th>
                            <th>Read1路径</th>
                            <th>Read2路径</th>
                            <th>分析状态</th>
                        </tr>
                    </thead>
                    <tbody>
'''

for sample in samples:
    status_class = 'status-success' if sample['status'] == '成功' else ('status-failed' if sample['status'] == '失败' else 'status-skipped')
    html_content += f'''
                        <tr>
                            <td>{sample['num']}</td>
                            <td><strong>{sample['name']}</strong></td>
                            <td style="font-size: 0.85em; word-break: break-all;">{sample['read1']}</td>
                            <td style="font-size: 0.85em; word-break: break-all;">{sample['read2']}</td>
                            <td class="{status_class}">{sample['status']}</td>
                        </tr>
'''

html_content += '''
                    </tbody>
                </table>
            </div>

            <!-- 比对统计 -->
            <div class="section">
                <h2>🎯 Kallisto 比对统计</h2>
                <table>
                    <thead>
                        <tr>
                            <th>样品名称</th>
                            <th>处理的reads数</th>
                            <th>比对上的reads数</th>
                            <th>比对率(%)</th>
                        </tr>
                    </thead>
                    <tbody>
'''

for stat in alignment_stats:
    html_content += f'''
                        <tr>
                            <td><strong>{stat['name']}</strong></td>
                            <td>{stat['n_processed']}</td>
                            <td>{stat['n_pseudoaligned']}</td>
                            <td><strong>{stat['rate']}%</strong></td>
                        </tr>
'''

html_content += f'''
                    </tbody>
                </table>
            </div>

            <!-- 输出文件 -->
            <div class="section">
                <h2>📁 输出文件</h2>
                <div class="info-grid">
                    <div class="info-card">
                        <h3>TPM 矩阵</h3>
                        <p>TPM.csv</p>
                    </div>
                    <div class="info-card">
                        <h3>Count 矩阵</h3>
                        <p>estCounts.csv</p>
                    </div>
                    <div class="info-card">
                        <h3>Markdown 摘要</h3>
                        <p>summary.md</p>
                    </div>
                    <div class="info-card">
                        <h3>HTML 摘要</h3>
                        <p>summary.html</p>
                    </div>
                </div>
            </div>
'''

# 添加可视化部分
if tpm_data and correlation_data and pca_data:
    html_content += '''
            <!-- 数据可视化 -->
            <div class="section">
                <h2>📊 数据可视化分析</h2>
                <div class="tab-container">
                    <div class="tab-buttons">
                        <button class="tab-button active" onclick="showTab('tpm-dist')">TPM分布</button>
                        <button class="tab-button" onclick="showTab('correlation')">样品相关性</button>
                        <button class="tab-button" onclick="showTab('pca')">PCA分析</button>
                    </div>
                    
                    <div id="tpm-dist" class="tab-content active">
                        <div class="chart-container">
                            <div id="tpm-plot" style="width:100%; height:500px;"></div>
                        </div>
                    </div>
                    
                    <div id="correlation" class="tab-content">
                        <div class="chart-container">
                            <div id="corr-plot" style="width:100%; height:600px;"></div>
                        </div>
                    </div>
                    
                    <div id="pca" class="tab-content">
                        <div class="chart-container">
                            <div id="pca-plot" style="width:100%; height:600px;"></div>
                        </div>
                    </div>
                </div>
            </div>
'''

html_content += '''
        </div>

        <div class="footer">
            <p><strong>作者：朱云涛（北京林业大学）</strong></p>
            <p>欢迎访问 <a href="https://yuntaobioinformatics.wordpress.com/" target="_blank">我的博客</a></p>
            <p style="margin-top: 10px; font-size: 0.9em;">Kallisto RNA-Seq 批量定量分析软件 v2.3</p>
        </div>
    </div>
'''

# 添加JavaScript代码
if tpm_data and correlation_data and pca_data:
    # 准备数据
    tpm_json = json.dumps(tpm_data)
    corr_json = json.dumps(correlation_data)
    pca_json = json.dumps(pca_data)
    
    html_content += f'''
    <script>
        // Tab切换功能
        function showTab(tabId) {{
            // 隐藏所有tab内容
            document.querySelectorAll('.tab-content').forEach(content => {{
                content.classList.remove('active');
            }});
            // 移除所有按钮的active状态
            document.querySelectorAll('.tab-button').forEach(button => {{
                button.classList.remove('active');
            }});
            // 显示选中的tab
            document.getElementById(tabId).classList.add('active');
            // 激活对应按钮
            event.target.classList.add('active');
        }}
        
        // TPM分布数据
        var tpmData = {tpm_json};
        
        // 1. 绘制TPM分布箱线图
        var tpmTraces = [];
        for (var sample in tpmData) {{
            tpmTraces.push({{
                y: tpmData[sample],
                type: 'box',
                name: sample,
                boxmean: 'sd'
            }});
        }}
        
        var tpmLayout = {{
            title: {{
                text: 'TPM分布箱线图（log10标度）',
                font: {{size: 20, color: '#667eea'}}
            }},
            yaxis: {{
                title: 'log10(TPM + 0.1)',
                gridcolor: '#e0e0e0'
            }},
            xaxis: {{
                title: '样品',
                tickangle: -45
            }},
            plot_bgcolor: '#f8f9fa',
            paper_bgcolor: 'white',
            showlegend: false,
            hovermode: 'closest'
        }};
        
        Plotly.newPlot('tpm-plot', tpmTraces, tpmLayout, {{responsive: true}});
        
        // 2. 绘制相关性热图
        var corrData = {corr_json};
        
        var corrTrace = [{{
            z: corrData.matrix,
            x: corrData.samples,
            y: corrData.samples,
            type: 'heatmap',
            colorscale: [
                [0, '#0d47a1'],
                [0.5, '#ffffff'],
                [1, '#b71c1c']
            ],
            zmid: 0.5,
            text: corrData.matrix.map(row => 
                row.map(val => val.toFixed(3))
            ),
            hovertemplate: '%{{y}} vs %{{x}}<br>相关性: %{{text}}<extra></extra>',
            colorbar: {{
                title: '相关系数',
                titleside: 'right'
            }}
        }}];
        
        var corrLayout = {{
            title: {{
                text: '样品间相关性热图（层次聚类排序）',
                font: {{size: 20, color: '#667eea'}}
            }},
            xaxis: {{
                tickangle: -45,
                side: 'bottom'
            }},
            yaxis: {{
                tickangle: 0
            }},
            plot_bgcolor: 'white',
            paper_bgcolor: 'white',
            width: null,
            height: 600
        }};
        
        Plotly.newPlot('corr-plot', corrTrace, corrLayout, {{responsive: true}});
        
        // 3. 绘制PCA图
        var pcaData = {pca_json};
        
        var pcaTrace = [{{
            x: pcaData.pc1,
            y: pcaData.pc2,
            mode: 'markers+text',
            type: 'scatter',
            text: pcaData.samples,
            textposition: 'top center',
            marker: {{
                size: 12,
                color: '#667eea',
                line: {{
                    color: '#764ba2',
                    width: 2
                }}
            }},
            hovertemplate: '<b>%{{text}}</b><br>PC1: %{{x:.2f}}<br>PC2: %{{y:.2f}}<extra></extra>'
        }}];
        
        var pcaLayout = {{
            title: {{
                text: 'PCA分析（基于前2000个高变异基因）',
                font: {{size: 20, color: '#667eea'}}
            }},
            xaxis: {{
                title: `PC1 (${{pcaData.explained_variance[0]}}%)`,
                gridcolor: '#e0e0e0',
                zeroline: true,
                zerolinecolor: '#999',
                zerolinewidth: 2
            }},
            yaxis: {{
                title: `PC2 (${{pcaData.explained_variance[1]}}%)`,
                gridcolor: '#e0e0e0',
                zeroline: true,
                zerolinecolor: '#999',
                zerolinewidth: 2
            }},
            plot_bgcolor: '#f8f9fa',
            paper_bgcolor: 'white',
            hovermode: 'closest',
            showlegend: false
        }};
        
        Plotly.newPlot('pca-plot', pcaTrace, pcaLayout, {{responsive: true}});
    </script>
'''

html_content += '''
</body>
</html>
'''

# 保存HTML文件
html_file = os.path.join(summary_dir, 'summary.html')
with open(html_file, 'w', encoding='utf-8') as f:
    f.write(html_content)

print(f"HTML摘要已生成: {html_file}")
HTMLEOF

echo "============================================================"
echo "所有分析完成！"
echo "完成时间: `date '+%Y-%m-%d %H:%M:%S'`"
echo "结果主目录: $OUTDIR"
echo "详细结果: $DETAILS_DIR"
echo "汇总结果: $SUMMARY_DIR"
echo "- TPM矩阵: $SUMMARY_DIR/TPM.csv"
echo "- Count矩阵: $SUMMARY_DIR/estCounts.csv"
echo "- Markdown摘要: $SUMMARY_DIR/summary.md"
echo "- HTML摘要: $SUMMARY_DIR/summary.html"
echo "============================================================"
