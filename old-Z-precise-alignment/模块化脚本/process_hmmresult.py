#!/usr/bin/env python3

#对hmmseach结果进行筛选，覆盖原始结果文件

import pandas as pd
import re
import argparse
import shutil
import yaml
from log_config import setup_logging
import sys

# hmmbuild tmp_inputfasta_clustao.hmm tmp_inputfasta_clustao.fasta
# hmmsearch --cpu 8 -o hmmsearch_output.txt tmp_inputfasta_clustao.hmm /public/group_share_data/TonyWuLab/zzd/db/medusa-2023/all-hum-medusa-db-genmoic-protein

# 加载YAML配置文件的函数
def load_config(config_file):
    with open(config_file, 'r') as file:
        return yaml.safe_load(file)

def excat_hmm_postition(input_file, output_file,hmm_evalue,hmm_score):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        # 初始化字典来存储数据
        data = {}
        
        for line in infile:
            line = line.strip()
            # 匹配细菌名称行
            match_bacteria = re.match(r'^>> (\S+)\s+.*', line)
            if match_bacteria:
                bacteria = match_bacteria.group(1)
                # 初始化或重置细菌数据
                if bacteria not in data:
                    data[bacteria] = {'domains': 0, 'score': 0, 'evalue': 0,
                                      'hmmfrom': 50000, 'hmmto': 0,
                                      'alignfrom': 50000, 'alignto': 0}
            
            # 匹配数据行
            match_data = re.match(r'(\d) \!\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\d+)\s+.*\s+(\d+)\s+(\d+)\s+.*\s+(\d+)\s+(\d+).*\S+', line)
            if match_data:
                domain, score, evalue, hmmfrom, hmmto, alignfrom, alignto = match_data.groups()[0], float(match_data.groups()[1]), float(match_data.groups()[3]), int(match_data.groups()[5]), int(match_data.groups()[6]), int(match_data.groups()[7]), int(match_data.groups()[8])
                # 更新数据
                data[bacteria]['domains'] += 1
                data[bacteria]['score'] += score
                data[bacteria]['evalue'] += evalue
                data[bacteria]['hmmfrom'] = min(data[bacteria]['hmmfrom'], hmmfrom)
                data[bacteria]['hmmto'] = max(data[bacteria]['hmmto'], hmmto)
                data[bacteria]['alignfrom'] = min(data[bacteria]['alignfrom'], alignfrom)
                data[bacteria]['alignto'] = max(data[bacteria]['alignto'], alignto)

        # 写入表头
        outfile.write("Bac_protein\tdomains\tscore\tmean_evalue\thmm_from\thmm_to\talign_from\talign_to\n")
        
        # 写入数据
        for bacteria, info in data.items():
            if info['evalue'] <= float(hmm_evalue) and info['score'] >= hmm_score :
                mean_evalue = info['evalue'] / info['domains']
                hmm_cov = (abs(info['hmmto'] - info['hmmfrom']) + 1) / 63
                outfile.write(f"{bacteria}\t{info['domains']}\t{info['score']}\t{mean_evalue}\t{info['hmmfrom']}\t{info['hmmto']}\t{info['alignfrom']}\t{info['alignto']}\n")


def compare_files_and_keep_df1_only(file1, file2, output_file, logger):
#与blast结果相比较，将只属于blast结果、不属于hmm的物种的GCF删除

    # 读取文件1，假设第一行是标题行
    print(file1)
    df1 = pd.read_csv(file1, sep='\t', header=0)
    
    # 读取文件2，并只处理第一列
    df2 = pd.read_csv(file2, sep='\t', header=0, usecols=['Bac_protein'])
    
    # 使用正则表达式提取每个条目的第二个_之前的内容
    #df2['Bac_protein'] = df2['Bac_protein'].apply(lambda x: re.match(r'(.*?_.*?)_.*', x).group(1))
    df2['Bac_protein'] = df2['Bac_protein']
    
    # 对df2的'Bac_protein'列去重
    #unique_df2_values = df2['Bac_protein'].drop_duplicates()
    unique_df2_values = df2['Bac_protein']

    # 筛选df1中第二列存在于df2处理后并去重的'Bac_protein'列的行
    filtered_df1 = df1[df1[df1.columns[1]].isin(unique_df2_values)]

    # 保存筛选后的结果，包含标题行
    filtered_df1.to_csv(output_file, index=False, sep='\t')
    
    # 统计行数
    rows_in_file1 = len(df1.index)
    rows_in_output = len(filtered_df1.index)
    difference = rows_in_file1 - rows_in_output

    #打印日志
    #print(f"Merged result saved to {output_file}")
    logger.info(f"Merged result saved to {output_file}")
    #print(f"Rows in file1: {rows_in_file1}")
    logger.info(f"blast organison in {file1}: {rows_in_file1}")
    #print(f"Rows in output file: {rows_in_output}")
    logger.info(f"HMM filtered in output {file2}: {rows_in_output}")
    #print(f"Difference in rows: {difference}")
    logger.info(f"Difference in rows: {difference}")

def main():
    # 设置命令行参数
    parser = argparse.ArgumentParser(description='Compare two files and keep rows from file1 that match patterns from file2, including the header of file1.')
    parser.add_argument('--hmmsearch_output', type=str, help='Output file name for hmmsearch results.', required=True)
    parser.add_argument('--blastoutput',  type=str, help='Path to the BLAST output file', metavar='BLAST_OUTPUT')
    parser.add_argument('--hmm_evalue', default='1e-4',type=str, help='hmmsearch evalue  (default:1e-4)', metavar='evalue')
    parser.add_argument('--hmm_score', default='200',type=int, help='hmmsearch score  (default:200)', metavar='score')
    parser.add_argument('--config', type=str, default='config.yaml', help='Configuration file path (default: config.yaml)', metavar='CONFIG_FILE')
    parser.add_argument('--log_file', '-l', type=str, help='Optional log file path', metavar='LOG_FILE')
    parser.add_argument('--hmm_screen_output_file', help='The output file path to save the filtered results.')

    # 解析命令行参数
    args = parser.parse_args()

    # 加载配置
    config = load_config(args.config)

    # 配置日志并获取 logger 实例
    logger = setup_logging(args.log_file)

    # 临时文件
    tmp_hmmsearch_excat = 'tmp_hmmsearch_excat.csv'
    tmp_with_hmm_normailzationout = 'tmp_with_hmm_normailzationout.csv'

    # 使用配置文件中的默认值，如果命令行中提供了值，则覆盖
    hmm_evalue = args.hmm_evalue if args.hmm_evalue is not None else config['hmm_evalue']
    hmm_score = args.hmm_score if args.hmm_score is not None else config['hmm_score']

    # 调用函数执行比较和保存结果
    excat_hmm_postition(args.hmmsearch_output, tmp_hmmsearch_excat, hmm_evalue, hmm_score)
    compare_files_and_keep_df1_only(args.blastoutput, tmp_hmmsearch_excat, tmp_with_hmm_normailzationout, logger)
    # 覆盖原始的normailzationout
    shutil.move(tmp_with_hmm_normailzationout, args.blastoutput)


if __name__ == '__main__':
    main()

