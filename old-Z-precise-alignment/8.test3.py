def read_names_and_queries(filename):
    """
    从文件中读取名称和查询序列
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
    在目标序列中查找查询序列的第一个氨基酸位置
    """
    return target_sequence.find(query_sequence) + 1  # 返回位置，从1开始计数

def main():
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

if __name__ == "__main__":
    main()
