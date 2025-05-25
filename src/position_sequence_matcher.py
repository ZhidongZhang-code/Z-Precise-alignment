"""
此脚本专为处理FASTA格式的蛋白质序列文件设计，旨在识别并提取包含特定氨基酸位置模式的序列。
它允许用户指定一个或多个氨基酸的位置以及相应的氨基酸类型，然后从一个大型的FASTA文件中
查找并提取所有符合这些条件的序列。提取的序列将被保存到指定的输出文件中。

主要步骤包括：
1. 设置参考位置索引：根据第一条序列（通常作为参考序列）和用户指定的位置与氨基酸对应关系，
   确定在序列中的实际索引位置。
2. 提取并追加序列：对于每个序列，如果它在指定的参考位置含有正确的氨基酸，则提取该序列
   并将其添加到结果列表中。每个提取的序列将包括序列ID和其后续的10个非"-"字符（氨基酸）。
3. 从FASTA文件中查找序列：读取FASTA文件，并基于上述逻辑查找所有符合条件的序列。
4. 将结果写入输出文件：所有提取的序列将被写入到用户指定的输出文件中，每个序列一行，
   序列ID前带有'>'标记。

使用示例：
    

注意：在使用本脚本前，请确保输入的FASTA文件格式正确，并且已经准备好了指定的位置和氨基酸信息。
"""

import os
import argparse

def parse_position_arg(position_str):
    """
    解析命令行参数中的位置和氨基酸信息。
    参数 position_str 是一个字符串列表，每个元素格式为 '位置:氨基酸'。
    返回一个字典，键为位置（整数），值为对应的氨基酸（字符串）。
    """
    positions = {}
    for pos_aa in position_str:
        pos, aa = pos_aa.split(':')
        positions[int(pos)] = aa
    return positions

def set_reference_indices(sequence, positions, reference_indices):
    """
    根据第一条序列（参考序列）设置参考位置索引。
    
    Args:
        sequence (str): 当前处理的序列。
        positions (dict): 指定的位置和对应的氨基酸。
        reference_indices (dict): 存储序列中特定位置的索引和对应的氨基酸。
    """
    count = 0  # 计数非"-"（gap）字符
    for i, char in enumerate(sequence):
        if char != '-':
            count += 1
            if count in positions:
                reference_indices[i] = positions[count]  # 保存实际索引位置和期望的氨基酸

def extract_and_append_sequence(seq_id, sequence, reference_indices, rf_sequences):
    """
    提取满足条件的序列并添加到结果列表中。
    
    Args:
        seq_id (str): 序列的ID。
        sequence (str): 序列数据。
        reference_indices (dict): 参考位置索引和期望的氨基酸。
        rf_sequences (list): 存储满足条件的序列ID和序列数据的列表。
    """
    # 检查序列是否包含所有指定的参考位置和氨基酸
    if all(sequence[idx] == aa for idx, aa in reference_indices.items()):
        first_index = min(reference_indices.keys())
        rf_seq, count = "", 0
        for i in range(first_index, len(sequence)):
            if sequence[i] != '-':
                rf_seq += sequence[i]
                count += 1
                if count == 10:  # 只提取10个非gap字符
                    break
        if count == 10:
            rf_sequences.append((seq_id, rf_seq))  # 添加到结果列表

def find_sequences_with_rf(fasta_file, positions):
    """
    从FASTA文件中找到包含特定参考位置氨基酸的序列。
    
    Args:
        fasta_file (str): 输入的FASTA文件路径。
        positions (dict): 指定的位置和对应的氨基酸。
        
    Returns:
        list: 包含满足条件的序列ID和序列数据的列表。
    """
    with open(fasta_file, 'r') as f:
        lines = f.readlines()

    rf_sequences = []
    reference_indices = {}
    current_id = None
    current_sequence = ""
    first_id = None

    for line in lines:
        line = line.strip()
        if line.startswith(">"):
            if current_id and current_id != first_id:
                extract_and_append_sequence(current_id, current_sequence, reference_indices, rf_sequences)
            current_id = line[1:]
            if first_id is None:
                first_id = current_id
            current_sequence = ""
        else:
            if not reference_indices and current_sequence == "":
                set_reference_indices(line, positions, reference_indices)
            current_sequence += line

    if current_id and current_id != first_id:
        extract_and_append_sequence(current_id, current_sequence, reference_indices, rf_sequences)

    return rf_sequences

def main(fasta_file, positions, output_file='output.txt'):
    """
    主函数：处理输入FASTA文件，并将结果写入输出文件。
    
    Args:
        fasta_file (str): 输入的FASTA文件路径。
        positions (dict): 指定的位置和对应的氨基酸。
        output_file (str): 结果输出文件的路径。
    """
    # 输入fasta文件路径
    fasta_file = args.input
    # 输入序列位置信息
    positions = parse_position_arg(args.position)
    # 获取所有满足条件的序列的ID以及相应的锚点序列信息
    rf_sequences = find_sequences_with_rf(fasta_file, positions)
    #保存到文件
    with open(output_file, 'w') as out_file:
        for sequence_id, anchor_seq in rf_sequences:
            out_file.write(f">{sequence_id} {anchor_seq}\n")
    print("结果已保存到output.txt文件中")

    parser = argparse.ArgumentParser(description="Process protein sequences")
    parser.add_argument("--input", help="input fasta file")
    parser.add_argument("--position", nargs='+', type=str, help="positions and corresponding amino acids")
    args = parser.parse_args()
#    matches_protein = 'matches-protein.fasta'
#    matches_cds = 'matches-cds.fasta'

if __name__ == "__main__":

    main()