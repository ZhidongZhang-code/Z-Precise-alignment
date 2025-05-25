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

def run_pfam(input,pfam_params, logger):
    logger.info('pfam analysis started.')

    # 根据 blast_type 构建命令行参数
    command = [
        "python3",
        pfam_params['pfam_path'],  
        input, 
        pfam_params['pfamdb_path'], 
        '-out', pfam_params['out'],
        '-cpu', str(pfam_params['cpu'])
    ]

    try:
        subprocess.run(command, check=True)
        logger.info(f'pfam_scan command executed successfully.')
    except subprocess.CalledProcessError as e:
        logger.error(f'Error executing pfam_scan command: {e}')
        sys.exit(1)
        
    logger.info(f'{pfam_params["pfam_scan"]} analysis completed.')

def run():
    parser = argparse.ArgumentParser(description='Run pfam analysis.')
    parser.add_argument('--extractprotein', type=str, default='alignment_extract_protein.fasta', help='File where the complete protein sequences are saved (default: alignment_extract_protein.fasta)', metavar='PROTEIN_FILE')
    #parser.add_argument('--blast_type', choices=['blastp', 'blastx'], default='blastp',required=True, help='Type of BLAST program to use')
    #parser.add_argument('--query', '-q', required=True, help='Path to the query file')
    #parser.add_argument('--pfamoutput', '-o',default='tmp_pfam_output.csv', help='Path to the output file')
    #parser.add_argument('--evalue', '-e', default='1e-5', help='E-value threshold')
    #parser.add_argument('--max_target_seqs', default='5000', help='Maximum number of target sequences')
    #parser.add_argument('--matrix', default='BLOSUM62', help='Scoring matrix name (only used for blastp)')
    parser.add_argument('--cpu', '-t', default=8, help='Number of threads')
    parser.add_argument('--log_file', '-l', help='Log file path (optional)')
    parser.add_argument('--config', default='config.yaml', help='Configuration file path')
    args = parser.parse_args()

    # 加载配置
    config = load_config(args.config)

    # 使用 setup_logging 配置日志
    logger = setup_logging(log_file_path=args.log_file)

    # 临时文件
    tmp_pfam_output = "tmp_pfam_output.csv"

    # 将命令行参数整理成字典
    pfam_params = {
        'pfam_path': config["pfam_scan"],  
        'pfamdb_path': config['pfam_dirextory'],
        'cpu': args.cpu,
        'out': tmp_pfam_output ,
    }

    run_pfam(args.extractprotein,pfam_params, logger)

if __name__ == "__main__":
    run()
