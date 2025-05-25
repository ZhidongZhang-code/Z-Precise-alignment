"""
此脚本用于从FASTA格式的蛋白质序列文件中，基于特定的名称列表，查找并提取匹配的序列记录。它首先从一个文本文件中读取序列名称，
然后在一个较大的FASTA文件中查找这些名称对应的序列，并将所有找到的匹配序列保存到另一个文件中。这个过程支持生物信息学分析，
尤其是在需要从大量序列数据中筛选特定序列时。

使用步骤如下：
1. 将需要查找的序列名称列表保存在一个文本文件中，每个名称位于以'>'开头的行。
2. 准备一个FASTA格式的文件，其中包含了大量的蛋白质序列。
3. 指定一个输出文件路径，用于保存匹配的序列。

脚本包含四个主要函数：
- read_names_from_file: 从文件中读取名称列表。
- find_matching_sequences: 在FASTA文件中查找与名称列表匹配的序列。
- save_matches_to_file: 将找到的匹配序列保存到指定的文件中。
- main: 将上述步骤组织起来，从读取名称到保存匹配序列的完整流程。

此脚本的设计目的是为了提高生物信息学数据处理的效率，特别是在处理大规模序列数据时，通过自动化筛选特定序列，来支持研究和分析。

使用示例：
    

注意：在使用本脚本前，请确保你已经有了包含名称的文本文件和FASTA格式的蛋白质序列文件。
"""

import argparse
from log_config import setup_logging

def read_names_from_file(file_path):
    """
    从文件中读取并返回所有以'>'开头的行的第一个单词，去掉'>'。

    Args:
        file_path (str): 文件的路径。

    Returns:
        list: 包含所有提取的名称的列表。
    """
    names = []
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('>'):
                name = line.split()[0][1:]
                names.append(name)
    return names

def find_matching_sequences(input_fasta_file, names):
    """
    在FASTA文件中查找名称匹配的序列。

    Args:
        input_fasta_file (str): FASTA文件的路径。
        names (list): 需要匹配的名称列表。

    Returns:
        list: 包含匹配的序列记录的列表，每个记录是一个元组，包含头部和序列列表。
    """
    matches = []
    record = None
    with open(input_fasta_file, 'r') as fasta_file:
        for line in fasta_file:
            if line.startswith('>'):
                if record:
                    for name in names:
                        if name in record[0]:
                            matches.append(record)
                            break
                record = (line, [])
            elif record:
                record[1].append(line)

        # 检查最后一个记录是否匹配
        if record:
            for name in names:
                if name in record[0]:
                    matches.append(record)
                    break
    return matches

def save_matches_to_file(matches, output_file_path):
    """
    将所有匹配的序列保存到文件中。

    Args:
        matches (list): 匹配的序列记录列表。
        output_file_path (str): 输出文件的路径。
    """
    with open(output_file_path, 'w') as match_file:
        for match in matches:
            match_file.write(match[0])
            match_file.writelines(match[1])
#    logger.info(f'找到的符合条件的蛋白质序列保存在：{output_file_path}')

def main():

    parser = argparse.ArgumentParser(description="Process protein sequences")
    parser.add_argument("--input_file", help="input fasta file")
    parser.add_argument("--protein_db", help="input protein fasta file")
    parser.add_argument("--output_file", help="output protein fasta file")
    args = parser.parse_args()
    
    # 使用 setup_logging 配置日志
    logger = setup_logging(log_file_path=args.log_file)

    # 定义文件路径
    input_names_file = args.input_file
    input_fasta_file = args.protein_db
    output_file_path = args.output_file

    names = read_names_from_file(input_names_file)
    matches = find_matching_sequences(input_fasta_file, names)
    save_matches_to_file(matches, output_file_path,logger)

if __name__ == "__main__":

    main()
