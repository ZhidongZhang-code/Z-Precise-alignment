#!/usr/bin/env python3

import os
import re
import argparse
import yaml 
from concurrent.futures import ThreadPoolExecutor
from threading import Lock
from log_config import setup_logging


###############
###传入为临时blast文件，只带筛选

###############

# 加载YAML配置文件的函数
def load_config(config_file):
    with open(config_file, 'r') as file:
        return yaml.safe_load(file)



# 初始化锁
write_lock_cds = Lock()
write_lock_protein = Lock()

def delete_output_files(alignment_cds , alignment_protein):
    """初始化输出文件"""
    if os.path.exists(alignment_cds):
        os.remove(alignment_cds)
    if os.path.exists(alignment_protein):
        os.remove(alignment_protein)
""" 
def findextractcdsprotein(input_file, directory, num_threads ,
                          alignment_cds = 'alignment_extract_cds.fasta', 
                          alignment_protein = 'alignment_extract_protein.fasta', 
                          config_file='config.yaml'):
    """ 
def findextractcdsprotein(input_file, num_threads, alignment_cds, alignment_protein, blast_identity,blast_cover,logger, config):
    """
    处理输入文件，提取符合条件的标识符，并使用线程池并行处理序列文件。

    参数：
        input_file：包含数据的输入文件路径
        num_threads：要使用的线程数
    """

    logger.info('find extract cds and protein analysis started.')
    
    # 清理旧的输出文件，避免追加到已存在的文件
    delete_output_files(alignment_cds, alignment_protein)
    
    identifiers = set()  # 用于存储唯一的标识符
    directory = config['GCF_directory']  # 从配置文件中读取路径
    with open(input_file, 'r') as infile:
        for line in infile:
            line = line.strip()
            if re.match(r'^\w+', line):
                columns = line.split('\t')
                # 根据特定条件筛选标识符
                if (abs(float(columns[5]) - float(columns[4]) + 1) / int(columns[2]) >= blast_cover and      #（qstar-qend）/qlen
                        abs(float(columns[7]) - float(columns[6]) + 1) / int(float(columns[3])) >= blast_cover and      #（sstar-send）sqlen
                        float(columns[11]) >= blast_identity):      #identity
                    identifier = columns[1]
                    identifiers.add(identifier)
    
    logger.info(f'find extract cds and protein analysis completed.')
    logger.info(f'complete coding sequences(cds) sequence is saved in {alignment_cds}')
    logger.info(f'complete protein sequence is saved in {alignment_protein}')

    with ThreadPoolExecutor(max_workers=num_threads) as executor:
        futures = []
        for identifier in identifiers:
            gcf_identifier = re.search(r'(GCF_\d+\.\d+)', identifier).group(1)
#            wp_identifier = re.search(r'(WP_\d+\.\d+)', identifier).group(1)
            folder_path = os.path.join(directory, gcf_identifier)
            if os.path.exists(folder_path):
                cds_file_path = os.path.join(folder_path, "cds_from_genomic.fna")
                protein_file_path = os.path.join(folder_path, "protein.faa")

                if os.path.exists(cds_file_path) and os.path.exists(protein_file_path):
                    futures.append(executor.submit(process_sequence, identifier, cds_file_path, alignment_cds, write_lock_cds))
                    futures.append(executor.submit(process_sequence, identifier, protein_file_path, alignment_protein, write_lock_protein))

        # Optional: Wait for all futures to complete if needed
        # for future in futures:
        #     future.result()

def process_sequence(identifier, file_path, output_file, lock):
    """
    处理单个序列文件，提取序列并写入到输出文件中。

    参数：
        identifier：标识符 GCF_000027225.1_ASM2722v1_**_WP_012989815.1
        file_path：序列文件路径
        output_file：输出文件路径
        lock：线程锁
    """
    header, sequence = extract_sequence(identifier, file_path)
    if sequence:
        with lock:
            with open(output_file, 'a') as outfile:
                outfile.write(f"{header}\n{sequence}\n")

def extract_sequence(identifier, file_path):
    """
    从序列文件中提取特定标识符对应的序列。

    参数：
        identifier：标识符
        file_path：序列文件路径

    返回：
        header：包含标识符的序列标题
        sequence：提取的序列内容
    """
    header = ""
    sequence = ""
    gcf_identifier_search = re.search(r'(GCF_\d+\.\d+)', identifier)
    wp_identifier_search = re.search(r'_\*\*_(.+?)$', identifier)  # 除了GCF_000022005.1_ASM2200v1_**_WP_002518002.1还有GCF_000022005.1_ASM2200v1_**_YP_002518002.1，筛选条件要改成这种
    
    # 确保搜索到了GCF和WP标识符
    if gcf_identifier_search and wp_identifier_search:
        gcf_identifier = gcf_identifier_search.group(1)
        wp_identifier = wp_identifier_search.group(1)  # 提取 _**_ 后面的内容
        
        try:
            with open(file_path, 'r') as file:
                is_sequence_started = False
                for line in file:
                    # 如果行中同时包含GCF和提取的标识符，则开始提取序列
                    if gcf_identifier in line and wp_identifier in line:
                        header = line.strip()
                        is_sequence_started = True
                    elif is_sequence_started and line.startswith(">"):
                        break
                    elif is_sequence_started:
                        sequence += line.strip()
        except IOError as e:
            print(f"Error file {file_path}: {e}")
    else:
        print(f"Identifier format mismatch in {identifier}")
    return header, sequence

def main():
    parser = argparse.ArgumentParser(description="python3 find-extract-cds-protein.py input_file , Output two files output_cds.fasta and output_protein.fasta ")
    parser.add_argument("--input_file", help="Input file containing data.")
    parser.add_argument('--extractcds', type=str, default='alignment_extract_cds.fasta', help='complete coding sequences(cds) sequence is saved (default: alignment_extract_cds.fasta)')
    parser.add_argument('--extractprotein', type=str, default='alignment_extract_protein.fasta', help='complete protein sequence is saved (default: alignment_extract_protein.fasta)')
    parser.add_argument("--num_threads", type=int, default=8, help="Number of threads to use.")
    parser.add_argument('--log_file', '-l', type=str, help='Log file path (optional)')
    parser.add_argument('--coverage_threshold', '-cov',type=float, help='Coverage threshold (default 0.7)')
    parser.add_argument('--identity_threshold', '-id',type=float, help='Identity threshold (default 30)')
    parser.add_argument('--config', type=str, default='config.yaml', help='Configuration file path (default: config.yaml)')
    args = parser.parse_args()

    alignment_cds = args.extractcds   #'alignment_extract_cds.fasta'
    alignment_protein = args.extractprotein #'alignment_extract_protein.fasta'

    # 初始化原始文件
    delete_output_files(alignment_cds , alignment_protein)      
    # 加载配置
    config = load_config(args.config)

    # 使用配置文件中的默认值，如果命令行中提供了值，则覆盖
    coverage_threshold = args.coverage_threshold if args.coverage_threshold is not None else config['coverage_threshold']
    identity_threshold = args.identity_threshold if args.identity_threshold is not None else config['identity_threshold']


    ### 加载配置文件、log
    logger = setup_logging(args.log_file)
    config = load_config(args.config)  

    findextractcdsprotein(args.input_file, args.num_threads, alignment_cds ,alignment_protein, identity_threshold,coverage_threshold,logger, config)

if __name__ == "__main__":
    main()

