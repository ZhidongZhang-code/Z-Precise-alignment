def find_sequences_with_rf(fasta_file, positions):
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
    count = 0
    for i, char in enumerate(sequence):
        if char != '-':
            count += 1
            if count in positions:
                reference_indices[i] = positions[count]  # 保存实际索引位置和期望的氨基酸

def extract_and_append_sequence(seq_id, sequence, reference_indices, rf_sequences):
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

# 示例使用
fasta_file = 'clu-out.fasta'
positions = {25: 'R', 28: 'R'}  # 使用字典定义位置和期望的氨基酸

rf_sequences = find_sequences_with_rf(fasta_file, positions)

# 保存结果到文件
with open('output.txt', 'w') as out_file:
    for sequence_id, anchor_seq in rf_sequences:
        out_file.write(f">{sequence_id} {anchor_seq}\n")

print("结果已保存到output.txt文件中")

