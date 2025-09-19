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
#0 1 2 3 4 5 6 7 8 9 10 11
#qseqid sseqid qlen slen qstart qend sstart send evalue bitscore length pident
#K00097_**_eco:b0052 GCF_000190995.1_ASM19099v1_**_WP_000241242.1 329 329 1 329 1 329 1.9e-182 644.8 329 325

# 加载YAML配置文件的函数
def load_config(config_file):
    with open(config_file, 'r') as file:
        return yaml.safe_load(file)

def calculate_coverage(start, end, length):
    """计算覆盖率"""
    return abs(float(end) - float(start) + 1) / float(length)

def filenormalization(input_file, output_file, coverage_threshold, identity_threshold, logger):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        # 写入header，增加qcov和scov两列
        header = "qseqid\tsseqid\tqlen\tslen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tlength\tpident\tqcov\tscov\n"
        outfile.write(header)
        
        new_lines = []
        for line in infile:
            line = line.strip()
            if re.match(r'^\w+', line):
                parts = line.split('\t')
                if len(parts) >= 12:
                    qcov = calculate_coverage(parts[4], parts[5], parts[2])
                    scov = calculate_coverage(parts[6], parts[7], parts[3])
                    
                    if (qcov >= coverage_threshold and 
                        scov >= coverage_threshold and 
                        float(parts[11]) >= identity_threshold):
                        # 添加qcov和scov到行末尾
                        new_line = line + f"\t{qcov:.3f}\t{scov:.3f}\n"
                        new_lines.append(new_line)
        
        for line in new_lines:
            outfile.write(line)
        
        logger.info("File normalization processed successfully.")
    
    # 统计文件行数，用于日志记录
    with open(output_file, 'r') as f:
        line_count = sum(1 for _ in f)
    
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
