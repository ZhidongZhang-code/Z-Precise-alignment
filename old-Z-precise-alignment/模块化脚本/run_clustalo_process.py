#!/usr/bin/env python3

import argparse
import os
import sys
import yaml 
import subprocess
from FileOperation import remove_file
from log_config import setup_logging

# 加载YAML配置文件的函数
def load_config(config_file):
    with open(config_file, 'r') as file:
        return yaml.safe_load(file)

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

def merge_fasta_files(file1, file2, output_file):
    """
    将两个 FASTA 格式的文件合并成一个输出文件。

    Args:
        file1 (str): 第一个输入文件路径。
        file2 (str): 第二个输入文件路径。
        output_file (str): 输出文件路径。
    """
    with open(output_file, 'w') as output:
        with open(file1, 'r') as f1:
            output.write(f1.read())
        with open(file2, 'r') as f2:
            output.write(f2.read())

def run_clustalo(input_file, output_file, num_threads, logger, config):
    """
    运行 Clustal Omega 进行多序列比对。

    Args:
        input_file (str): 输入文件路径。
        output_file (str): 输出文件路径。
        num_threads (int): 使用的线程数，默认为 8。
    """
    logger.info(f"Starting clustalo for query file: {input_file} with {num_threads} threads")
    command = [
        config['clustalo_path'],
        '--infile', input_file,
        '-o', output_file,
        '--threads', str(num_threads),
        '--force'
    ]    
    
    try:
        subprocess.run(command, check=True)
        logger.info("clustalo command executed successfully.")
    except subprocess.CalledProcessError as e:
        logger.error(f"Error executing clustalo command: {e}")
        sys.exit(1)  # 出错时退出程序

#def remove_file(file):
#    """
#    删除指定文件。
#
#    Args:
#        file (str): 要删除的文件路径。
#    """
#    if os.path.exists(file):
#        os.remove(file)


def find_start_end_positions(alignment, sequence_id):
    """
    查找序列在多序列比对结果中的起始和结束位置。

    Args:
        alignment (dict): 多序列比对结果，格式为序列ID到序列的映射字典。
        sequence_id (str): 要查找的序列ID。

    Returns:
        tuple: 起始位置和结束位置的元组。
    """
    start_position = None
    end_position = None

    # 找到序列的起始位置
    for i in range(len(alignment[sequence_id])):
        if alignment[sequence_id][i] != '-':
            start_position = i
            break

    # 找到序列的结束位置
    for i in range(len(alignment[sequence_id]) - 1, -1, -1):
        if alignment[sequence_id][i] != '-':
            end_position = i
            break

    return start_position, end_position

def find_end_position(alignment, sequence_id):
    """
    查找序列在多序列比对结果中的结束位置。

    Args:
        alignment (dict): 多序列比对结果，格式为序列ID到序列的映射字典。
        sequence_id (str): 要查找的序列ID。

    Returns:
        int: 结束位置。
    """
    sequence = alignment[sequence_id]
    end_position = None
    for i in range(len(sequence) - 1, -1, -1):
        if sequence[i] != '-':
            end_position = i + 3  # 找到最后一个非'-'的位置后，往后数三位作为结束位置
            break
    return end_position

def trim_alignment(alignment, start_position, end_position, target_sequence_id, padding=3):
    """
    裁剪多序列比对结果，只保留指定范围内的序列。

    Args:
        alignment (dict): 多序列比对结果，格式为序列ID到序列的映射字典。
        start_position (int): 起始位置。
        end_position (int): 结束位置。
        target_sequence_id (str): 要裁剪的序列ID。
        padding (int): 裁剪时的填充大小，默认为 3。

    Returns:
        dict: 裁剪后的多序列比对结果。
    """
    trimmed_alignment = {}
    for seq_id, seq in alignment.items():
        trimmed_seq = seq[start_position - padding:end_position + padding + 1]  # 加一是因为切片时右边界不包含
        trimmed_alignment[seq_id] = trimmed_seq
    # 将指定序列移动到字典的第一个位置
    target_sequence = trimmed_alignment.pop(target_sequence_id)
    trimmed_alignment = {target_sequence_id: target_sequence, **trimmed_alignment}
    return trimmed_alignment

def read_alignment_from_file(file):
    """
    从文件中读取多序列比对结果。

    Args:
        file (str): 输入文件路径。

    Returns:
        dict: 多序列比对结果，格式为序列ID到序列的映射字典。
    """
    alignment = {}
    with open(file, 'r') as f:
        current_seq_id = None
        current_seq = ''
        for line in f:
            if line.startswith('>'):
                if current_seq_id:
                    alignment[current_seq_id] = current_seq
                current_seq_id = line.strip()
                current_seq = ''
            else:
                current_seq += line.strip()
        if current_seq_id:
            alignment[current_seq_id] = current_seq
    return alignment

def write_alignment_to_file(alignment, output_file):
    """
    将多序列比对结果写入文件。

    Args:
        alignment (dict): 多序列比对结果，格式为序列ID到序列的映射字典。
        output_file (str): 输出文件路径。
    """
    with open(output_file, 'w') as f:
        for seq_id, seq in alignment.items():
            f.write(seq_id + '\n')
            f.write(seq + '\n')

def main():
    parser = argparse.ArgumentParser(description="Tool for filtering duplicate proteins, merging FASTA files, and running Clustal Omega for multiple sequence alignment.")
    parser.add_argument("--target_sequence_file", help="File containing the target sequence(150nt).")
    parser.add_argument('--extractprotein', type=str, default='alignment_extract_protein.fasta', help='File where the complete protein sequences are saved (default: alignment_extract_protein.fasta)', metavar='PROTEIN_FILE')
    parser.add_argument("--clustalo_out", type=str ,help="clustalo output file .")
    parser.add_argument('--log_file', '-l', type=str, help='Log file path (optional)')
    parser.add_argument('--config', type=str, default='config.yaml', help='Configuration file path (default: config.yaml)')

    args = parser.parse_args(args)

    output_protein_file = args.extractprotein
    target_sequence_file = args.target_sequence_file
    temp_filtered_file = "temp_filtered.fasta"
    merged_file = "merged.fasta"
    clustal_output_file = "clu-out.fasta"
    num_threads = 8

    # 加载配置
    config = load_config(args.config)

    # 使用 setup_logging 配置日志
    logger = setup_logging(log_file_path=args.log_file)

    # 读取目标序列文件，获取序列ID
    with open(target_sequence_file, 'r') as f:
        target_sequence_id = f.readline().strip()

    remove_file(merged_file)  # 删除之前的 merged.fasta 文件
    filter_duplicate_proteins(output_protein_file, temp_filtered_file)  # 过滤重复的蛋白质ID并写入临时文件
    merge_fasta_files(temp_filtered_file, target_sequence_file, merged_file)
    run_clustalo(merged_file, clustal_output_file, num_threads, logger, config)
    remove_file(merged_file)  # 删除临时合并的文件
    remove_file(temp_filtered_file)  # 删除临时过滤的文件

    if args.clustalo_out : 
        clustal_output_file = args.clustalo_out

    # 从比对结果中裁剪指定序列并写入文件
    alignment = read_alignment_from_file(clustal_output_file)
    start_position, end_position = find_start_end_positions(alignment, target_sequence_id)
    end_position = find_end_position(alignment, target_sequence_id)
    trimmed_alignment = trim_alignment(alignment, start_position, end_position, target_sequence_id)
    write_alignment_to_file(trimmed_alignment, clustal_output_file)

if __name__ == "__main__":
    main()

