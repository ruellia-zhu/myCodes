import os
import subprocess
import argparse

if __name__ == "__main__":
    # SRA2FASTQ tool是将SRA文件批量转换为FASTQ格式的工具，遍历目录下所有子目录的SRA文件并转换
    print("========================== SRA2FASTQ tool is running... ==========================")
    print("Author: Yuntao Zhu")
    print("Version: 2.0 2025-07-29")
    print("This script will search for all SRA files in the specified directory.")
    print("The output FASTQ files will be stored in the same directory as the SRA files.")
    print("Usage: python sra2fastq.py <directory> [--delete-sra/-d] [--cpu/-c <threads>] [--gz/-g]")
    print("Example: python sra2fastq.py /path/to/sra_files -d -c 16 -g")
    print("===================================================================================")

def parse_args():
    parser = argparse.ArgumentParser(description="Convert SRA files to FASTQ format using fasterq-dump.")
    parser.add_argument("directory", type=str, help="The directory to search for SRA files.")
    parser.add_argument("-d", "--delete-sra", action="store_true", help="Delete original SRA files after successful conversion.")
    parser.add_argument("-c", "--cpu", type=int, default=8, help="Number of CPU threads to use (default: 8).")
    parser.add_argument("-g", "--gz", action="store_true", help="Compress output files to fastq.gz format.")
    return parser.parse_args()

# 转换 SRA 文件为 FASTQ 格式
def convert_sra_to_fastq(directory, threads=10, delete_sra=False, compress_gz=False):
    # 遍历指定目录下的所有子目录
    for root, dirs, files in os.walk(directory):
        for dir_name in dirs:
            # 构建子目录路径
            dir_path = os.path.join(root, dir_name)
            
            # 查找 SRA 文件
            for file in os.listdir(dir_path):
                if file.endswith(".sra"):
                    sra_file_path = os.path.join(dir_path, file)
                    
                    # 构建 fasterq-dump 命令（不再使用--gzip参数）
                    cmd = f"fasterq-dump --threads {threads} --outdir {dir_path} --split-3 {sra_file_path}"
                    
                    format_info = "fastq.gz" if compress_gz else "fastq"
                    print(f"Converting {sra_file_path} to {format_info} with {threads} threads...")
                    
                    # 执行fasterq-dump命令
                    result = subprocess.run(cmd, shell=True)
                    
                    # 检查转换是否成功
                    if result.returncode == 0:
                        print(f"Converted {sra_file_path} to fastq in {dir_path}")
                        
                        # 如果需要压缩，使用pigz压缩生成的fastq文件
                        if compress_gz:
                            sra_basename = os.path.splitext(file)[0]  # 去掉.sra扩展名
                            
                            # 查找生成的fastq文件并压缩
                            fastq_files = []
                            for fastq_file in os.listdir(dir_path):
                                if fastq_file.startswith(sra_basename) and fastq_file.endswith(".fastq"):
                                    fastq_files.append(os.path.join(dir_path, fastq_file))
                            
                            # 使用pigz压缩所有fastq文件
                            for fastq_file in fastq_files:
                                print(f"Compressing {fastq_file} with pigz using {threads} threads...")
                                pigz_cmd = f"pigz -p {threads} {fastq_file}"
                                pigz_result = subprocess.run(pigz_cmd, shell=True)
                                
                                if pigz_result.returncode == 0:
                                    print(f"Successfully compressed {fastq_file}")
                                else:
                                    print(f"Failed to compress {fastq_file}")
                        
                        # 如果转换成功且用户选择删除SRA文件
                        if delete_sra:
                            try:
                                os.remove(sra_file_path)
                                print(f"Deleted original SRA file: {sra_file_path}")
                            except OSError as e:
                                print(f"Error deleting SRA file {sra_file_path}: {e}")
                    else:
                        print(f"Failed to convert {sra_file_path}")

if __name__ == "__main__":
    args = parse_args()    
    # 使用用户指定的参数
    convert_sra_to_fastq(args.directory, threads=args.cpu, delete_sra=args.delete_sra, compress_gz=args.gz)
