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

def find_sequences_with_rf(fasta_file, positions):
    """
    在FASTA文件中查找具有特定参考帧（RF）位置的序列。
    fasta_file 是FASTA格式的文件路径。
    positions 是一个字典，包含要查找的氨基酸位置和类型。
    返回满足条件的序列ID和对应的锚点序列。
    """
    with open(fasta_file, 'r') as f:
        lines = f.readlines()

    rf_sequences = []
    reference_indices = {}  # 存储参考位置索引和期望的氨基酸
    current_id = None
    current_sequence = ""
    first_id = None  # 用于存储第一条序列的ID

    for line in lines:
        line = line.strip()
        if line.startswith(">"):
            # 处理上一个序列
            if current_id and current_id != first_id:  # 排除第一条序列
                extract_and_append_sequence(current_id, current_sequence, reference_indices, rf_sequences)
            # 准备读取新序列
            current_id = line[1:]  # 移除">"获取ID
            if first_id is None:
                first_id = current_id  # 记录第一条序列的ID
            current_sequence = ""
        else:
            if not reference_indices and current_sequence == "":  # 如果是第一条序列，确定参考位置
                set_reference_indices(line, positions, reference_indices)
            current_sequence += line  # 继续构建当前序列

    # 处理最后一个序列，排除第一条序列
    if current_id and current_id != first_id:
        extract_and_append_sequence(current_id, current_sequence, reference_indices, rf_sequences)

    return rf_sequences

def set_reference_indices(sequence, positions, reference_indices):
    """
    根据第一条序列确定参考位置索引。
    sequence 是序列字符串。
    positions 是位置和氨基酸的映射字典。
    reference_indices 是用于存储实际索引位置和期望的氨基酸。
    """
    count = 0
    for i, char in enumerate(sequence):
        if char != '-':
            count += 1
            if count in positions:
                reference_indices[i] = positions[count]  # 保存实际索引位置和期望的氨基酸

def extract_and_append_sequence(seq_id, sequence, reference_indices, rf_sequences):
    """
    检查序列是否匹配参考索引中的氨基酸，并提取锚点序列。
    seq_id 是序列的ID。
    sequence 是序列字符串。
    reference_indices 是参考位置索引和期望的氨基酸。
    rf_sequences 是结果列表，存储匹配的序列ID和锚点序列。
    """
    if all(sequence[idx] == aa for idx, aa in reference_indices.items()):
        # 提取接下来的10个非"-"氨基酸序列
        first_index = min(reference_indices.keys())
        rf_seq, count = "", 0
        for i in range(first_index, len(sequence)):
            if sequence[i] != '-':
                rf_seq += sequence[i]
                count += 1
                if count == 10:
                    break
        if count == 10:  # 确保提取了足够的氨基酸
            rf_sequences.append((seq_id, rf_seq))


def read_names_and_queries(filename):
    """
    从文件中读取名称和查询序列。
    filename 是包含名称和查询序列的文件路径。
    返回两个列表：names 包含名称，queries 是名称到查询序列的映射字典。
    """
    names = []
    queries = {}
    with open(filename, 'r') as file:
        for line in file:
            if line.startswith('>'):
                name, query = line.strip().split()
                name = name[1:]
                names.append(name)
                queries[name] = query
    return names, queries

def find_first_aa_position(target_sequence, query_sequence):
    """
    在目标序列中查找查询序列的第一个氨基酸位置。
    target_sequence 是目标序列字符串。
    query_sequence 是查询序列字符串。
    返回查询序列在目标序列中的起始位置（从1开始计数）。
    """
    return target_sequence.find(query_sequence) + 1  # 返回位置，从1开始计数

def parse_cds_fasta(filename):
    """
    解析CDS FASTA文件，返回一个字典，键是CDS的ID，值是序列。
    filename 是CDS FASTA文件的路径。
    """
    cds_entries = {}
    with open(filename, 'r') as file:
        current_id = None
        current_sequence = ''
        for line in file:
            if line.startswith('>'):
                if current_id is not None:
                    cds_entries[current_id] = current_sequence
                current_id = line.strip().split(' ')[0][1:]
                current_sequence = ''
            else:
                current_sequence += line.strip()
        if current_id is not None:
            cds_entries[current_id] = current_sequence
    return cds_entries

def extract_nucleotide_sequence(cds_entries, protein_id, genome_id, position):
    """
    从CDS条目中提取与给定蛋白质ID、基因组ID和位置相匹配的核苷酸序列。
    cds_entries 是CDS条目的字典，protein_id 是蛋白质ID，genome_id 是基因组ID，position 是位置。
    返回匹配的核苷酸序列，如果没有找到匹配，则返回None。
    """
    for cds_id, sequence in cds_entries.items():
        if protein_id in cds_id and genome_id in cds_id:
            cds_position = int(position) * 3  # 一个氨基酸由三个核苷酸组成
            start = max(0, cds_position - 75)
            end = min(len(sequence), cds_position + 75)
            nucleotide_sequence = sequence[start:end]
            return nucleotide_sequence
    return None


def main(args):
    # 输入fasta文件路径
    fasta_file = args.input
    positions = parse_position_arg(args.position)
    # 获取所有满足条件的序列的ID以及相应的锚点序列信息
    rf_sequences = find_sequences_with_rf(fasta_file, positions)
    # 保存结果到文件
    with open('output.txt', 'w') as out_file:
        for sequence_id, anchor_seq in rf_sequences:
            out_file.write(f">{sequence_id} {anchor_seq}\n")

    # 读取output.txt文件以获取所有name
    names = []
    with open('output.txt', 'r') as file:
        for line in file:
            if line.startswith('>'):
                name = line.split()[0][1:]  # 移除开头的'>'字符并提取name
                names.append(name)

    # 打开output_protein.fasta文件，并查找所有匹配的序列
    matches = []
    with open('output_protein.fasta', 'r') as fasta_file:
        record = None
        for line in fasta_file:
            if line.startswith('>'):
                # 如果这一行是一个新的记录的开始，检查之前的记录是否匹配
                if record:
                    # 检查是否有任何name出现在当前记录的ID中
                    for name in names:
                        if name in record[0]:
                            matches.append(record)
                            break
                # 开始新的记录
                record = (line, [])
            elif record:
                # 如果这是序列的一部分，添加到当前记录
                record[1].append(line)

        # 检查最后一个记录是否匹配
        if record:
            for name in names:
                if name in record[0]:
                    matches.append(record)
                    break

    # 将所有匹配的序列保存到一个文件中
    with open(matches_protein, 'w') as match_file:
        for match in matches:
            match_file.write(match[0])  # 写入序列的头部
            match_file.writelines(match[1])  # 写入序列数据

    # 读取名称和查询序列
    names, queries = read_names_and_queries('output.txt')

    # 打开output_protein.fasta文件，并查找所有匹配的序列
    matches = []
    with open('output_protein.fasta', 'r') as fasta_file:
        record = None
        for line in fasta_file:
            if line.startswith('>'):
                # 如果这一行是一个新的记录的开始，检查之前的记录是否匹配
                if record:
                    for name in names:
                        if name in record[0]:
                            # 在匹配的序列中查找查询序列的位置
                            position = find_first_aa_position(''.join(record[1]), queries[name])
                            matches.append((name, record[0], position))
                            break
                # 开始新的记录
                record = (line, [])
            elif record:
                # 如果这是序列的一部分，添加到当前记录
                record[1].append(line)

        # 检查最后一个记录是否匹配
        if record:
            for name in names:
                if name in record[0]:
                    # 在匹配的序列中查找查询序列的位置
                    position = find_first_aa_position(''.join(record[1]), queries[name])
                    matches.append((name, record[0], position))
                    break

    # 将结果保存到文件
    with open('output_matches.txt', 'w') as output_file:
        for match in matches:
            name = match[1].split('_**_')[0][1:].strip()
            output_file.write(f"{match[0]},{name},{match[2]}\n")

    # 解析output_cds.fasta文件
    cds_entries = parse_cds_fasta('output_cds.fasta')

    # 打开文件以写入匹配结果
    with open(matches_cds, 'w') as result_file:
        # 读取output_matches.txt中的匹配结果
        with open('output_matches.txt', 'r') as matches_file:
            for line in matches_file:
                protein_id, genome_id, position = line.strip().split(',')
                nucleotide_sequence = extract_nucleotide_sequence(cds_entries, protein_id, genome_id, position)
                if nucleotide_sequence is not None:
                    result_file.write(f">{genome_id}_**_{protein_id}  {int(position)-75}--{int(position)+72}\n{nucleotide_sequence}\n")
                else:
                    result_file.write(f"No match found for protein ID {protein_id}, genome ID {genome_id}, position {position}\n")

    # 删除不再需要的文件
    os.remove('output.txt')
    os.remove('output_matches.txt')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process protein sequences")
    parser.add_argument("--input", help="input fasta file")
    parser.add_argument("--position", nargs='+', type=str, help="positions and corresponding amino acids")
    args = parser.parse_args()
    matches_protein = 'matches-protein.fasta'
    matches_cds = 'matches-cds.fasta'
    main(args)
