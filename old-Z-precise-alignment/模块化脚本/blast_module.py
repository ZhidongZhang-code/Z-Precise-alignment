#!/usr/bin/env python3

import subprocess
import argparse
import logging
import sys
import yaml
from log_config import setup_logging

#blastp -db /public/group_share_data/TonyWuLab/zzd/db/medusa-2023/all-hum-medusa-db-genmoic-protein \
#-query MurC.fasta -out blast-MurC.csv \
#-outfmt "6 qseqid sseqid qlen slen qstart qend sstart send evalue bitscore length pident" \
#-evalue 1e-5 -max_target_seqs 50000 -matrix BLOSUM62 -num_threads 8



# 加载YAML配置文件的函数
def load_config(config_file):
    with open(config_file, 'r') as file:
        return yaml.safe_load(file)

def run_blast(query_file, output_file, blast_params, logger):
    logger.info('BLAST analysis started.')

    # 根据 blast_type 构建命令行参数
    command = [
        blast_params['blast_path'],  # 根据配置文件和用户选择确定BLAST程序的路径
        '-db', blast_params['db_path'],
        '-query', query_file,
        '-out', output_file,
        '-outfmt', '6 qseqid sseqid qlen slen qstart qend sstart send evalue bitscore length pident',
        '-evalue', blast_params['evalue'],
        '-max_target_seqs', blast_params['max_target_seqs'],
        '-num_threads', str(blast_params['num_threads'])
    ]

    # 如果是 blast，添加矩阵参数
    if blast_params['blast_type'] == 'blastp':
        command.extend(['-matrix', blast_params['matrix']])

    try:
        subprocess.run(command, check=True)
        logger.info(f"{blast_params['blast_type']} command executed successfully.")
    except subprocess.CalledProcessError as e:
        logger.error(f"Error executing {blast_params['blast_type']} command: {e}")
        sys.exit(1)
        
    logger.info(f'{blast_params["blast_type"]} analysis completed.')

def run():
    parser = argparse.ArgumentParser(description='Run BLAST analysis.')
    #parser.add_argument('--diamond_type',  default='diamond', help='Type of diamond BLAST program to use')
    parser.add_argument('--blast_type', choices=['blastp', 'blastx'], default='blastp',required=True, help='Type of BLAST program to use')
    parser.add_argument('--query', '-q', required=True, help='Path to the query file')
    parser.add_argument('--outblastput', '-o', required=True, help='Path to the output file')
    parser.add_argument('--evalue', '-e', default='1e-5', help='E-value threshold')
    parser.add_argument('--max_target_seqs', default='5000', help='Maximum number of target sequences')
    parser.add_argument('--matrix', default='BLOSUM62', help='Scoring matrix name (only used for blastp)')
    parser.add_argument('--num_threads', '-t', default=8, help='Number of threads')
    parser.add_argument('--log_file', '-l', help='Log file path (optional)')
    parser.add_argument('--config', default='config.yaml', help='Configuration file path')
    args = parser.parse_args()

    # 加载配置
    config = load_config(args.config)

    # 使用 setup_logging 配置日志
    logger = setup_logging(log_file_path=args.log_file)

    # 将命令行参数整理成字典
    blast_params = {
        'blast_path': config[args.blast_type + '_path'],  # 根据选择的BLAST类型确定路径
        'db_path': config['db_path'],
        'evalue': args.evalue,
        'max_target_seqs': args.max_target_seqs,
        'matrix': args.matrix,
        'num_threads': args.num_threads,
        'blast_type': args.blast_type,
    }

    run_blast(args.query, args.outblastput, blast_params, logger)

if __name__ == "__main__":
    run()
