#!/usr/bin/env python3

import sys
import re
import argparse
import logging
import yaml 
from log_config import setup_logging

###序列比对结果输出文件标准化，输入的序列比对结果：-outfmt', '6 qseqid sseqid qlen slen qstart qend sstart send evalue bitscore length pident
###筛选cov大于70%，id大于30的结果

##输入文件格式：
#0                          1                                             2        3        4       5       6      7            8           9        10     11      
#qseqid                     sseqid                                       qlen    slen    qstart  qend    sstart  send       evalue      bitscore length pident
#K00097_**_eco:b0052     GCF_000190995.1_ASM19099v1_**_WP_000241242.1    329     329     1       329     1       329     1.9e-182        644.8   329     325

# 加载YAML配置文件的函数
def load_config(config_file):
    with open(config_file, 'r') as file:
        return yaml.safe_load(file)
    
    
def filenormalization(input_file, output_file, coverage_threshold, identity_threshold, logger):
#    coverage_threshold = config['coverage_threshold']
#    identity_threshold = config['identity_threshold']

    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        new_lines = []
        for line in infile:
            line = line.strip()
            if re.match(r'^\w+', line):
                test = line.split('\t')
                if (abs(float(test[5]) - float(test[4]) + 1) / float(test[2]) >= coverage_threshold and
                    abs(float(test[7]) - float(test[6]) + 1) / float(test[3]) >= coverage_threshold and
                    float(test[11]) >= identity_threshold):
                    new_lines.append(line + '\n')

        for line in new_lines:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                col1, col2 = parts[0], parts[1]
                col1 = col1.split("_**_")[0]    #去掉KO信息
                #col2 = "_".join(col2.split("_", 2)[:2]) #去掉蛋白信息
                new_line = "{}\t{}\t{}\n".format(col2, col1, '\t'.join(parts[2:]))
                outfile.write(new_line)
        logger.info("File normalization processed successfully.")
    # 统计文件行数，用于日志记录
    with open(output_file, 'r') as f:
        line_count = sum(1 for _ in f)
    #logger.info(f"Total lines written to output file: {line_count}")
    logger.info(f"blast processing complete. Results have been written to {output_file}, with total lines: {line_count}.\n")

def main():
    parser = argparse.ArgumentParser(description='File normalization script for alignment result files.')
    parser.add_argument('--alignmentout', type=str, required=True, help='Input file name , alignment result')
    parser.add_argument('--normailzationout', type=str, required=True, help='normailzation output file name')
    parser.add_argument('--coverage_threshold', type=float, help='Coverage threshold (default 0.7)')
    parser.add_argument('--identity_threshold', type=float, help='Identity threshold (default 30)')
    parser.add_argument('--log', type=str, help='Log file path (optional)')
    parser.add_argument('--config', type=str, default='config.yaml', help='Configuration file path (default: config.yaml)')
    args = parser.parse_args()

    # 加载配置
    config = load_config(args.config)

    # 配置日志并获取 logger 实例
    logger = setup_logging(args.log)

    # 使用配置文件中的默认值，如果命令行中提供了值，则覆盖
    coverage_threshold = args.coverage_threshold if args.coverage_threshold is not None else config['coverage_threshold']
    identity_threshold = args.identity_threshold if args.identity_threshold is not None else config['identity_threshold']

    # 调用文件处理函数
    filenormalization(args.alignmentout, args.normailzationout, coverage_threshold, identity_threshold, logger)


if __name__ == "__main__":
    main()
