#!/usr/bin/env python3

#对pfam结果进行筛选，覆盖原始结果文件

import pandas as pd
import re
import argparse
import shutil
import yaml
from log_config import setup_logging
import sys

# 加载YAML配置文件的函数
def load_config(config_file):
    with open(config_file, 'r') as file:
        return yaml.safe_load(file)
    
def pfam_compare_files_and_keep_df1_only(pfamfile, blastoutfile, keys , output_file, logger):
#
    # 加载PFAM文件
    pfam_df = pd.read_csv(pfamfile, sep=',')
    # 筛选出hmm_name为Mur_ligase_C的记录，并处理seq_id列
    filtered_pfam_df = pfam_df[pfam_df['hmm_name'].isin(keys)].copy()
    #filtered_pfam_df['processed_seq_id'] = filtered_pfam_df['seq_id'].apply(lambda x: '_'.join(x.split('_')[:2]))
    filtered_pfam_df['processed_seq_id'] = filtered_pfam_df['seq_id']
    # 加载BLAST文件
    blast_df = pd.read_csv(blastoutfile, sep='\t')
    # 筛选BLAST文件中存在于处理后PFAM seq_id中的行
    #filtered_blast_df = blast_df[blast_df.iloc[:,1].isin(filtered_pfam_df['processed_seq_id'])]
    filtered_blast_df = blast_df[blast_df.iloc[:,1].apply(lambda x: x.split('_**_')[-1]).isin(filtered_pfam_df['processed_seq_id'])]
    # 找出被删除的行
    #removed_blast_df = blast_df[~blast_df.iloc[:,1].isin(filtered_pfam_df['processed_seq_id'])]
    removed_blast_df = blast_df[~blast_df.iloc[:,1].apply(lambda x: x.split('_**_')[-1]).isin(filtered_pfam_df['processed_seq_id'])]
    # 保存筛选后的结果，包含标题行
    filtered_blast_df.to_csv(output_file, index=False, sep='\t')

    # 统计行数
    rows_in_filtered = len(filtered_blast_df.index)
    rows_in_removed = len(removed_blast_df.index)
    #difference = rows_in_file1 - rows_in_output

    #打印日志
    logger.info(f"After screening, {rows_in_removed} sequences were removed not with {keys}")
    logger.info(f"After screening, {rows_in_filtered} sequences were retained centain {keys}")
    #logger.info(f"Rows in output file: {rows_in_output}")
    #logger.info(f"Difference in rows: {difference}")

def main():
    # 设置命令行参数
    parser = argparse.ArgumentParser(description='Compare two files and keep rows from file1 that match patterns from file2, including the header of file1.')
    parser.add_argument('--pfamfile', type=str, required=True, help='Output file name for hmmsearch results.')
    parser.add_argument('--blastoutfile', type=str, required=True, help='blast Normalization output file name', metavar='NORM_OUTPUT')
    parser.add_argument('--pfamkey', nargs='+', type=str, help='One or more pfam domain keys , Example: "Mur_ligase_C Another_HMM" ', required=True)
    parser.add_argument('--log_file', '-l', type=str, help='Optional log file path', metavar='LOG_FILE')
    parser.add_argument('--output_file', help='The output file path to save the filtered results.')

    # 解析命令行参数
    args = parser.parse_args()

    # 加载配置
    logger = setup_logging(args.log_file)
    #config = load_config(args.config)  

    # 临时文件
    tmp_pfam_excat = 'tmp_pfam_excat.csv'


    # 调用函数执行比较和保存结果
    #excat_hmm_postition(args.target_sequence_file, tmp_hmmsearch_excat)
    pfam_compare_files_and_keep_df1_only(args.pfamfile, args.blastoutfile, args.pfamkey, tmp_pfam_excat,logger)
    # 覆盖原始的normailzationout
    shutil.move(tmp_pfam_excat, args.blastoutfile)

    #删除临时文件


if __name__ == '__main__':
    main()

