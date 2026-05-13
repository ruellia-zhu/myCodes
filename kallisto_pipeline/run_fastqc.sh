#!/bin/bash

# 欢迎界面
echo "============================================================"
echo "FastQC 批量质控分析软件 v1.0 2025.10.10"
echo "作者：朱云涛（北京林业大学）"
echo "------------------------------------------------------------"
echo "用法：sh $0 <INDIR> [THREADS]"
echo "INDIR: 包含多个样品子文件夹的目录，每个子文件夹内包含测序文件*.fq.gz/*.fastq.gz"
echo "       INDIR一般为存放Rawdata的目录"
echo "THREADS: 可选参数，指定线程数，默认为8，FastQC内部多线程主要用于单个样品内的并行"
echo "用例：sh $0 /RNA-seq/30个兰州百合转录组/rawdata 16"
echo "------------------------------------------------------------"
echo "功能："
echo "1、自动识别fq.gz和fastq.gz格式的测序文件"
echo "2、为每个样品生成FastQC报告(HTML + ZIP)"
echo "3、输出详细的分析摘要报告"
echo "4、自动跳过结果目录，避免重复处理"
echo "5、自动调用MultiQC汇总FastQC结果"
echo "------------------------------------------------------------"
echo "其他说明："
echo "1、结果将保存到与INDIR同级的fastqc_result目录下"
echo "   例如：输入/RNA-seq/rawdata，输出/RNA-seq/fastqc_result"
echo "   每个样品有独立的子文件夹存放FastQC生成的报告"
echo "2、使用前请确保已经切换到RNA-seq或其他配置好的conda环境"
echo "3、线程数可通过第二个参数自定义，默认为8"

# 判断参数个数
if [ $# -lt 1 ]; then
  echo "参数数量不对！"
  echo "用法: sh $0 <INDIR> [THREADS]"
  echo "INDIR: 包含多个样品子文件夹的目录，每个子文件夹内包含测序文件*.fq.gz/*.fastq.gz"
  echo "THREADS: 可选参数，指定线程数，默认为8"
  exit 1
fi

INDIR=$1
THREADS=${2:-8}

# 检查fastqc是否安装
if ! command -v fastqc >/dev/null 2>&1; then
  echo "错误：未检测到fastqc命令，请先安装或激活包含FastQC的环境"
  exit 1
fi

# 检查输入目录是否存在
if [ ! -d "$INDIR" ]; then
  echo "错误：输入目录不存在！($INDIR)"
  exit 1
fi

# 获取输入目录的父目录
PARENT_DIR=$(dirname "$INDIR")

# 输出目录设置为与输入目录同级
OUTDIR="$PARENT_DIR/fastqc_result"
mkdir -p "$OUTDIR"

# 创建MultiQC目录
MULTIQC_DIR="$OUTDIR/multiQC"
mkdir -p "$MULTIQC_DIR"

# 验证线程数是否为有效数字
if ! echo "$THREADS" | grep -q '^[0-9][0-9]*$' || [ "$THREADS" -le 0 ]; then
  echo "警告: 线程数 '$THREADS' 无效，使用默认值8"
  THREADS=8
fi

# 初始化样品信息记录
SAMPLE_INFO="$MULTIQC_DIR/sample_info.txt"
printf "序号\t样品名称\t文件数量\t文件列表\t分析状态\n" > "$SAMPLE_INFO"

# 获取FastQC版本信息
FASTQC_VERSION=$(fastqc --version 2>&1 | head -n 1 || echo "未知版本")

# 构建示例命令显示
audited_command="fastqc -t $THREADS -o fastqc_output_dir sample_1.fq.gz sample_2.fq.gz"

echo "已经创建好结果目录，分析即将开始......"
echo "输入目录: $INDIR"
echo "FastQC结果目录: $OUTDIR"
echo "MultiQC报告目录: $MULTIQC_DIR"
echo "使用线程数: $THREADS"
echo "FastQC版本: $FASTQC_VERSION"
echo "------------------------------------------------------------"

sample_count=0
success_count=0

declare -A sample_status

for sample in "$INDIR"/*; do
  if [ -d "$sample" ]; then
    name=$(basename "$sample")

    # 跳过结果目录
    case "$name" in
      kallisto*|fastp*|fastqc*|cleandata|multiQC)
        echo "跳过结果目录: $name"
        continue
        ;;
    esac

    sample_count=$((sample_count + 1))

    # 收集样品内所有fastq/fq.gz文件
    mapfile -t fastq_files < <(find "$sample" -maxdepth 1 -type f \( -name "*.fq.gz" -o -name "*.fastq.gz" \) | sort)

    if [ ${#fastq_files[@]} -eq 0 ]; then
      echo "警告: 样品 $name 未找到任何fastq文件，跳过"
      printf "%s\t%s\t0\t未找到\t跳过\n" "$sample_count" "$name" >> "$SAMPLE_INFO"
      sample_status[$name]="跳过"
      continue
    fi

    echo "正在处理样品 $name ..."

    SAMPLE_OUT_DIR="$OUTDIR/$name"
    mkdir -p "$SAMPLE_OUT_DIR"

    # 运行FastQC
    if fastqc -t "$THREADS" -o "$SAMPLE_OUT_DIR" "${fastq_files[@]}"; then
      file_list=$(printf "%s" "${fastq_files[*]}")
      printf "%s\t%s\t%d\t%s\t成功\n" "$sample_count" "$name" "${#fastq_files[@]}" "$file_list" >> "$SAMPLE_INFO"
      sample_status[$name]="成功"
      success_count=$((success_count + 1))
      echo "样品 $name FastQC分析完成"
    else
      file_list=$(printf "%s" "${fastq_files[*]}")
      printf "%s\t%s\t%d\t%s\t失败\n" "$sample_count" "$name" "${#fastq_files[@]}" "$file_list" >> "$SAMPLE_INFO"
      sample_status[$name]="失败"
      echo "警告: 样品 $name FastQC分析失败"
    fi
  fi
done

echo "============================================================"
echo "FastQC分析完成！"
echo "总样品数: $sample_count"
echo "成功分析: $success_count"
echo "正在生成汇总报告..."

SUMMARY_MD="$MULTIQC_DIR/summary.md"

cat > "$SUMMARY_MD" << EOF
# FastQC RNA-Seq 批量质控分析 v1.0 摘要报告

by 朱云涛 at $(date)
欢迎访问[我的博客](https://yuntaobioinformatics.wordpress.com/)
Welcome to visit [My Blog](https://yuntaobioinformatics.wordpress.com/)

## 1 分析基本信息

- **FastQC版本**: $FASTQC_VERSION
- **完成时间**: $(date)
- **FastQC命令示例**: \
  \
  \`$audited_command\`
- **输入目录**: $INDIR
- **结果目录**: $OUTDIR

## 2 样品分析统计

- **总样品数**: $sample_count
- **成功分析**: $success_count
- **失败样品**: $((sample_count - success_count))

| 序号 | 样品名称 | 文件数量 | 文件列表 | 分析状态 |
|------|----------|----------|----------|----------|
EOF

# 将样品信息写入Markdown表格
tail -n +2 "$SAMPLE_INFO" | while IFS=$'\t' read -r idx sample_name file_num file_list status; do
  safe_list=${file_list//|/\\|}
  echo "| $idx | $sample_name | $file_num | $safe_list | $status |" >> "$SUMMARY_MD"
done

echo "\n## 3 FastQC指标统计" >> "$SUMMARY_MD"
echo "" >> "$SUMMARY_MD"
echo "| 样品名称 | FASTQ文件 | PASS数量 | WARN数量 | FAIL数量 |" >> "$SUMMARY_MD"
echo "|-----------|-----------|-----------|-----------|-----------|" >> "$SUMMARY_MD"

# 从FastQC zip内的summary.txt提取通过/警告/失败数量
for sample_dir in "$OUTDIR"/*; do
  [ -d "$sample_dir" ] || continue
  sample_name=$(basename "$sample_dir")
  [ "$sample_name" = "multiQC" ] && continue

  for zip_file in "$sample_dir"/*_fastqc.zip; do
    [ -f "$zip_file" ] || continue
    fastq_base=$(basename "$zip_file" "_fastqc.zip")
    read pass warn fail < <(unzip -p "$zip_file" summary.txt 2>/dev/null | awk '
      BEGIN { pass=0; warn=0; fail=0 }
      { if ($1=="PASS") pass++; else if ($1=="WARN") warn++; else if ($1=="FAIL") fail++; }
      END { printf "%d %d %d", pass, warn, fail }
    ')
    echo "| $sample_name | $fastq_base | ${pass:-0} | ${warn:-0} | ${fail:-0} |" >> "$SUMMARY_MD"
  done
done

echo "" >> "$SUMMARY_MD"
echo "> 提示：PASS/WARN/FAIL数量来源于FastQC summary.txt。更多细节请查看对应HTML报告。" >> "$SUMMARY_MD"

# 运行MultiQC生成报告
echo "正在运行MultiQC生成质量控制报告..."
if command -v multiqc >/dev/null 2>&1; then
  multiqc "$OUTDIR" -o "$MULTIQC_DIR" --title "FastQC RNA-Seq Quality Control Report" --comment "Generated by FastQC批量质控分析软件 v1.0" --force
  echo "MultiQC报告已生成到: $MULTIQC_DIR"
else
  echo "警告: MultiQC未安装，跳过质量控制报告生成"
  echo "可以使用以下命令安装: conda install -c bioconda multiqc"
  echo "你是不是忘记切换conda环境了？conda activate RNAseq"
fi

echo "============================================================"
echo "所有分析完成！"
echo "FastQC结果目录: $OUTDIR"
echo "汇总报告: $MULTIQC_DIR"
echo "- 分析摘要: $SUMMARY_MD"
echo "- MultiQC报告: $MULTIQC_DIR/multiqc_report.html"
echo "============================================================"
