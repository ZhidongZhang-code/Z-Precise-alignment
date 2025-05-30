#!/usr/bin/env python3

"""
根据输入文件，提取目标位点附近上游24个，下游25个氨基酸位置的核苷酸序列
用来构建bowtie2的索引

###输入文件格式
WP_003018763.1,GCF_029104325.1_ASM2910432v1,383
WP_006234258.1,GCF_000169035.1_ASM16903v1,369
WP_013451841.1,GCF_000183405.1_ASM18340v1,359
WP_013479277.1,GCF_000175215.2_ASM17521v2,374
WP_006042661.1,GCF_000153285.1_ASM15328v1,373


"""

def parse_cds_fasta(filename):
    """
    解析output_cds.fasta文件
    返回一个字典，其中键是CDS的ID，值是序列
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
    从解析的CDS条目中提取特定蛋白质ID和基因组ID对应的核苷酸序列。

    参数:
    cds_entries (dict): 解析的CDS条目，键是CDS的ID，值是序列。
    protein_id (str): 需要提取序列的蛋白质ID。
    genome_id (str): 需要提取序列的基因组ID。
    position (int): 蛋白质序列中感兴趣的氨基酸的位置。

    返回:
    str 或 None: 如果找到匹配的序列，则返回核苷酸序列字符串；否则返回None。
    """
    for cds_id, sequence in cds_entries.items():
        if protein_id in cds_id and genome_id in cds_id:
            cds_position = int(position) * 3  # 一个氨基酸由三个核苷酸组成
            start = max(0, cds_position - 75)
            end = min(len(sequence), cds_position + 75)
            nucleotide_sequence = sequence[start:end]
            return cds_position,nucleotide_sequence
    return None

def process_cds_fasta(cds_filename, matches_filename, output_filename):
    """
    处理CDS FASTA文件和匹配结果文件，将匹配的核苷酸序列保存到指定的输出文件。

    参数:
    cds_filename (str): CDS FASTA文件的路径。
    matches_filename (str): 匹配结果文件的路径。
    output_filename (str): 输出文件的路径，用于保存匹配的核苷酸序列信息。

    返回:
    None: 将匹配的核苷酸序列信息保存到指定的输出文件。
    """
    
    cds_entries = parse_cds_fasta(cds_filename)
    with open(matches_filename, 'r') as matches_file, open(output_filename, 'w') as output_file:
        for line in matches_file:
            protein_id, genome_id, position = line.strip().split(',')
            result = extract_nucleotide_sequence(cds_entries, protein_id, genome_id, position)
            if result is not None:
                cds_position, nucleotide_sequence = result
                position_int = int(position)  # 将字符串转换为整数
                output_file.write(f">{genome_id}_**_{protein_id}||{position_int - 24}:{position_int + 25}\n{nucleotide_sequence}\n")
            else:
                # 如果你想在没有找到匹配项时做一些操作，可以在这里添加代码
                pass
