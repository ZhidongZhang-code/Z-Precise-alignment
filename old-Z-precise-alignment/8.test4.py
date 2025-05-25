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
    从CDS条目中提取与给定蛋白质ID、基因组ID和位置相匹配的核苷酸序列
    """
    for cds_id, sequence in cds_entries.items():
        if protein_id in cds_id and genome_id in cds_id:
            cds_position = int(position) * 3  # 一个氨基酸由三个核苷酸组成
            start = max(0, cds_position - 75)
            end = min(len(sequence), cds_position + 75)
            nucleotide_sequence = sequence[start:end]
            return nucleotide_sequence
    return None

def main():
    # 解析output_cds.fasta文件
    cds_entries = parse_cds_fasta('output_cds.fasta')

    # 读取output_matches.txt中的匹配结果
    with open('output_matches.txt', 'r') as matches_file:
        for line in matches_file:
            protein_id, genome_id, position = line.strip().split(',')
            nucleotide_sequence = extract_nucleotide_sequence(cds_entries, protein_id, genome_id, position)
            if nucleotide_sequence is not None:
                print(f"Match found for protein ID {protein_id}, genome ID {genome_id}, position {position}: {nucleotide_sequence}")
            else:
                print(f"No match found for protein ID {protein_id}, genome ID {genome_id}, position {position}")

if __name__ == "__main__":
    main()
