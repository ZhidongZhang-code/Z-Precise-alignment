#!/usr/bin/env python3
import os
import sys
import argparse
import logging
from log_config import setup_logging

# 支持的FASTA文件扩展名
VALID_EXTENSIONS = ('.fa', '.fasta', '.fna', '.faa')

def rename(input_file, output_file_path, prefix):
    """
    处理输入的 FASTA 文件，为序列头部添加指定的前缀。
    :param input_file: 输入文件的路径，期望为 FASTA 格式。
    :param output_file_path: 处理后的序列将被写入此输出文件中。
    :param prefix: 要添加到序列头部的前缀。
    """
    if not input_file.endswith(VALID_EXTENSIONS):
        logging.error(f"Input file must have one of the following extensions: {', '.join(VALID_EXTENSIONS)}")
        return

    if not os.path.isfile(input_file):
        logging.error(f"File not found: {input_file}")
        return

    modified_lines = []

    with open(input_file, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('>') and '_**_' not in line:
                modified_lines.append(f'>{prefix}_**_{line[1:]}')
            else:
                modified_lines.append(line)

    # 将修改后的内容写入输出文件
    modified_content = '\n'.join(modified_lines) + '\n'
    with open(output_file_path, 'w') as output_file:
        output_file.write(modified_content)

    logging.info(f"Processed: {input_file}")
    logging.info(f"Output saved to: {output_file_path}")

def main():
    parser = argparse.ArgumentParser(description='Process a FASTA file (.fa, .fasta, .fna, .faa) and modify headers with a custom prefix.')
    parser.add_argument('--fasta_input', type=str, required=True, help='Input FASTA file path (.fa, .fasta, .fna, .faa)')
    parser.add_argument('--output', type=str, required=True, help='Output FASTA file path')
    parser.add_argument('--prefix', type=str, required=True, help='Prefix to add to sequence headers')
    parser.add_argument('--log_file', type=str, help='Log file path (optional)')

    args = parser.parse_args()

    output_file_path = args.output
    input_file = args.fasta_input
    prefix = args.prefix

    # 使用 setup_logging 配置日志
    logger = setup_logging(log_file_path=args.log_file)

    rename(input_file, output_file_path, prefix)
    logger.info('fasta_renamer analysis completed.')

if __name__ == "__main__":
    main()

