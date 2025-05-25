"""
此脚本用于处理生物序列数据，特别是从特定的FASTA文件中找到与输入查询序列匹配的序列，并确定查询序列在匹配序列中的起始位置。脚本的主要工作流程如下：

1. 从一个文本文件（例如output.txt）中读取序列名称和对应的查询序列。该文件的每一行都包含一个序列名称和它的查询序列，格式为">序列名称 查询序列"。

2. 打开并逐行读取另一个FASTA格式的文件（例如output_protein.fasta），该文件包含多个生物序列。FASTA格式通常由一个以">"开头的标题行和随后的序列行组成。

3. 对于FASTA文件中的每个序列，脚本检查它是否包含任何输入文件中指定的查询序列。如果找到匹配的序列，脚本会确定查询序列在该匹配序列中的起始位置。

4. 所有找到的匹配项及其详细信息（包括序列名称、匹配的序列名称、以及查询序列的起始位置）将被保存到一个新的输出文件（例如output_matches.txt）中，用于后续分析。

此脚本可用于生物信息学领域，特别是在基因序列分析、蛋白质序列比对等方面，帮助研究者快速定位特定序列的位置，从而进行深入的生物学特性研究。

使用方法：
1. 确保输入文件和FASTA文件按照上述格式准备好。
2. 运行脚本，指定输入文件和FASTA文件的路径。
3. 查看输出文件，获取匹配序列的详细信息。

注意：此脚本需要Python环境支持，且依赖于标准的文件读写操作。
"""

import argparse
from log_config import setup_logging

def read_names_and_queries(filename):
    """
    从文件中读取名称和查询序列。
    """
    names = []
    queries = {}
    with open(filename, 'r') as file:
        for line in file:
            if line.startswith('>'):
                name, query = line.strip().split()
                name = name[1:]  # 移除">"
                names.append(name)
                queries[name] = query
    return names, queries

def find_first_aa_position(target_sequence, query_sequence):
    """
    在目标序列中查找查询序列的第一个氨基酸位置。
    """
    return target_sequence.find(query_sequence) + 1  # 返回位置，从1开始计数

def process_site_sequences(input_filename, fasta_filename, output_filename):
    """
    处理序列，找出匹配的序列及其位置，并将结果写入文件。
    """
    # 读取名称和查询序列
    names, queries = read_names_and_queries(input_filename)

    # 打开FASTA文件，并查找所有匹配的序列
    matches = []
    with open(fasta_filename, 'r') as fasta_file:
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
            else:
                # 如果这是序列的一部分，添加到当前记录
                record[1].append(line.strip())

        # 检查最后一个记录是否匹配
        if record:
            for name in names:
                if name in record[0]:
                    position = find_first_aa_position(''.join(record[1]), queries[name])
                    matches.append((name, record[0], position))
                    break

    # 将结果保存到文件
    with open(output_filename, 'w') as output_file:
        for match in matches:
            output_file.write(f"{match[0]},{match[1].split('_**_')[0][1:].strip()},{match[2]}\n")

# 防止脚本直接运行
if __name__ == "__main__":
    print("This script is not intended to be run directly.")
