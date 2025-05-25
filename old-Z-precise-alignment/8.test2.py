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
with open('matches.fasta', 'w') as match_file:
    for match in matches:
        match_file.write(match[0])  # 写入序列的头部
        match_file.writelines(match[1])  # 写入序列数据
