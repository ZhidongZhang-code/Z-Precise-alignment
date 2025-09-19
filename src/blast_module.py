#!/usr/bin/env python3

import subprocess
import argparse
import logging
import sys
import yaml
from log_config import setup_logging

# 加载YAML配置文件的函数
def load_config(config_file):
    with open(config_file, 'r') as file:
        return yaml.safe_load(file)

# 运行BLAST或DIAMOND的主函数
def run_blast(query_file, output_file, blast_params, logger):
    logger.info(f'{blast_params["blast_type"]} 分析开始.')

    # 根据不同类型生成BLAST或DIAMOND的命令
    if blast_params['use_diamond']:
        #print(blast_params['blast_path'])
        # DIAMOND 命令
        command = [
            blast_params['diamond_path'],
            blast_params['blast_type'],  # 提取blastp或blastx
            '--db', blast_params['diamond_db_path'],
            '--query', query_file,
            '-o', output_file,
            '--swipe',
            '-f', '6',  # 指定输出格式为tabular（数字6代表tabular格式）
            'qseqid', 'sseqid', 'qlen', 'slen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'length', 'pident',
            '--evalue', blast_params['evalue'],
            '--max-target-seqs', blast_params['max_target_seqs'],
            '--threads', str(blast_params['num_threads'])
        ]
    else:
        # 标准BLAST 命令
        command = [
            blast_params['blast_path'],  # 从配置中读取BLAST的路径
            '-db', blast_params['db_path'],
            '-query', query_file,
            '-out', output_file,
            '-outfmt', '6 qseqid sseqid qlen slen qstart qend sstart send evalue bitscore length pident',
            '-evalue', blast_params['evalue'],
            '-max_target_seqs', blast_params['max_target_seqs'],
            '-num_threads', str(blast_params['num_threads'])  # BLAST使用 -num_threads 参数
        ]

        # 如果选择的是blastp，则添加矩阵参数
        if blast_params['blast_type'] == 'blastp':
            command.extend(['-matrix', blast_params['matrix']])

    try:
        logger.info(f"执行命令: {' '.join(command)}")
        subprocess.run(command, check=True)
        logger.info(f"{blast_params['blast_type']} 命令执行成功.")
    except subprocess.CalledProcessError as e:
        logger.error(f"执行 {blast_params['blast_type']} 命令时出错: {e}")
        sys.exit(1)
        
    logger.info(f'{blast_params["blast_type"]} 分析完成.')

# 主运行函数
def run():
    parser = argparse.ArgumentParser(description='运行BLAST或DIAMOND分析程序.')
    parser.add_argument('--blast_type', choices=['blastp', 'blastx'], required=True, help='选择使用BLAST或DIAMOND程序')
    parser.add_argument('--use_diamond', action='store_true', help='Use DIAMOND instead of BLAST for sequence alignment')
    parser.add_argument('--query', '-q', required=True, help='查询文件的路径')
    parser.add_argument('--outblastput', '-o', required=True, help='输出文件的路径')
    parser.add_argument('--evalue', '-e', default='1e-5', help='E-value 阈值')
    parser.add_argument('--max_target_seqs', default='100000000000', help='最大目标序列数')
    parser.add_argument('--matrix', default='BLOSUM62', help='打分矩阵名称（仅用于blastp）')
    parser.add_argument('--num_threads', '-t', default=8, help='使用的线程数')
    parser.add_argument('--log_file', '-l', help='日志文件路径（可选）')
    parser.add_argument('--config', default='config.yaml', help='配置文件的路径')
    args = parser.parse_args()

    # 加载配置文件
    config = load_config(args.config)

    # 设置日志
    logger = setup_logging(log_file_path=args.log_file)

    # 整理BLAST或DIAMOND参数
    blast_params = {
        'blast_path': config[args.blast_type + '_path'],  # 根据用户选择的类型从配置文件中获取路径
        'db_path': config['db_path'],
        'evalue': args.evalue,
        'max_target_seqs': args.max_target_seqs,
        'matrix': args.matrix,
        'num_threads': args.num_threads,
        'blast_type': args.blast_type,
    }

    # 执行BLAST或DIAMOND分析
    run_blast(args.query, args.outblastput, blast_params, logger)

if __name__ == "__main__":
    run()
