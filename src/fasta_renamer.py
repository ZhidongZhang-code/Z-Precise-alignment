#!/usr/bin/env python3

import os
import sys
import argparse
import logging
from log_config import setup_logging

###将属于同一类型的蛋白质放在同一个文件中，并且加上_**_
###结尾必须为fasta
###使用方式python fatsa_renamer  --input Mura1.fasta Mura2.fasta --output Mura.fasta


def rename(input_files, output_file_path):
    """
    批量处理输入的 FASTA 文件，为每个文件中的序列头部添加文件名作为前缀。

    :param input_files: 输入文件的路径列表，期望为 FASTA 格式。
    :param output_file_path: 处理后的序列将被追加到此输出文件中。
    """
    for input_file_path in input_files:
        if not input_file_path.endswith(".fasta"):
            logging.warning(f"Ignoring non-fasta file: {input_file_path}")
            continue

        if not os.path.isfile(input_file_path):
            logging.error(f"File not found: {input_file_path}")
            continue

        file_name = os.path.splitext(os.path.basename(input_file_path))[0]
        modified_lines = []

        with open(input_file_path, 'r') as input_file:
            for line in input_file:
                line = line.strip()
                if line.startswith('>') and '_**_' not in line:
                    modified_lines.append(f'>{file_name}_**_{line[1:]}')
                else:
                    modified_lines.append(line)
        # 将修改后的内容写入输出文件，并且增加换行符
        modified_content = '\n'.join(modified_lines) + '\n'

        with open(output_file_path, 'a') as output_file:
            output_file.write(modified_content)

        logging.info(f"Processed: {input_file_path}")

    logging.info(f"All files have been processed. Output saved to: {output_file_path}")

def main():
    parser = argparse.ArgumentParser(description='Process input FASTA files and modify headers.')
    parser.add_argument('--fasta_input', nargs='+', type=str, help='Input FASTA file paths')
    parser.add_argument('--fasta_rename', type=str, help='Output FASTA file path')
    parser.add_argument('--log_file', type=str, help='Log file path (optional)')
    args = parser.parse_args()

    output_file_path = args.fasta_rename
    input_files = args.fasta_input

    # 使用 setup_logging 配置日志
    logger = setup_logging(log_file_path=args.log_file)

    rename(input_files, output_file_path)
    logger.info('fasta_renamer analysis completed.')

if __name__ == "__main__":
    main()
