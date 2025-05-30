#!/usr/bin/env python3

import argparse

def filter_duplicate_proteins(input_file, output_file):
    """
    从输入文件中筛选出唯一的蛋白质，并写入输出文件。

    Args:
        input_file (str): 输入文件路径。
        output_file (str): 输出文件路径。
    """
    unique_proteins = {}
    with open(input_file, 'r') as f:
        lines = f.readlines()
        for line in lines:
            if line.startswith(">"):
                protein_id = line.split()[0]
                wp_id = protein_id.split("_**_")[-1].split("\t")[0]  # 获取蛋白质ID中的特定部分
                if wp_id not in unique_proteins:
                    unique_proteins[wp_id] = True

    with open(output_file, 'w') as f:
        l = 0
        for line in lines:
            if line.startswith(">"):
                protein_id = line.split()[0]
                wp_id = protein_id.split("_**_")[-1].split("\t")[0]  # 获取蛋白质ID中的特定部分
                if wp_id in unique_proteins:
                    l = 1
                    f.write(">" + wp_id + "\n")  # 只保留 WP 部分
                    del unique_proteins[wp_id]
            elif l == 1 :
                f.write(line)
                l = 0



def main():
    parser = argparse.ArgumentParser(description="Tool for filtering duplicate proteins, merging FASTA files, and running Clustal Omega for multiple sequence alignment.")
    parser.add_argument('--input', '-i', type=str, help='input file path (optional)')
    parser.add_argument('--output', '-o', type=str, help='output file path (optional)')
    #parser.add_argument('--log_file', '-l', type=str, help='Log file path (optional)')
    #parser.add_argument('--config', type=str, default='config.yaml', help='Configuration file path (default: config.yaml)')

    args = parser.parse_args()

    filter_duplicate_proteins(args.input, args.output)  # 过滤重复的蛋白质ID并写入临时文件

if __name__ == "__main__":
    main()

