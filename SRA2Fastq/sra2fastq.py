import os
import subprocess
import argparse

if __name__ == "__main__":
    # SRA2FASTQ tool是将SRA文件批量转换为FASTQ格式的工具，遍历目录下所有子目录的SRA文件并转换
    print("========================== SRA2FASTQ tool is running... ==========================")
    print("Author: Yuntao Zhu")
    print("This script will search for all SRA files in the specified directory.")
    print("The output FASTQ files will be stored in the same directory as the SRA files.")
    print("Usage: python sra2fastq.py <directory> [--delete-sra/-d] [--cpu/-c <threads>]")
    print("Example: python sra2fastq.py "/media/desk16/tl5024/RNAseq" -d -c 16")
    print("===================================================================================")

def parse_args():
    parser = argparse.ArgumentParser(description="Convert SRA files to FASTQ format using fasterq-dump.")
    parser.add_argument("directory", type=str, help="The directory to search for SRA files.")
    parser.add_argument("-d", "--delete-sra", action="store_true", help="Delete original SRA files after successful conversion.")
    parser.add_argument("-c", "--cpu", type=int, default=8, help="Number of CPU threads to use (default: 8).")
    return parser.parse_args()

# 转换 SRA 文件为 FASTQ 格式
def convert_sra_to_fastq(directory, threads=8, delete_sra=False):
    # 遍历指定目录下的所有子目录
    for root, dirs, files in os.walk(directory):
        for dir_name in dirs:
            # 构建子目录路径
            dir_path = os.path.join(root, dir_name)
            
            # 查找 SRA 文件
            for file in os.listdir(dir_path):
                if file.endswith(".sra"):
                    sra_file_path = os.path.join(dir_path, file)
                    
                    # 使用 fasterq-dump 转换 SRA 文件为 FASTQ 格式
                    # 输出文件将存储在相同的子目录下，并设置线程数
                    cmd = f"fasterq-dump --threads {threads} --outdir {dir_path} --split-3 {sra_file_path}"
                    
                    print(f"Converting {sra_file_path} with {threads} threads...")
                    
                    # 执行命令
                    result = subprocess.run(cmd, shell=True)
                    
                    # 检查转换是否成功
                    if result.returncode == 0:
                        print(f"Converted {sra_file_path} to FASTQ in {dir_path}")
                        
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
    # 使用用户指定的线程数
    convert_sra_to_fastq(args.directory, threads=args.cpu, delete_sra=args.delete_sra)
