#!/usr/bin/env python3

"""
makeformat.py - 处理并转换生物数据为表格格式

版本：1.1.0 (2025年9月19日)
作者：ZhidongZhang-code

更新日志：
----------
2025年9月19日：
- 添加 --unique_by 参数，用于基于分类水平的数据去重
- 将 GCF 编号设置为第一列
- 改进 GCF ID 匹配逻辑，支持忽略版本号
- 重新组织列顺序，优化数据展示
- 增强基于 sum 值的去重逻辑

使用方式：
    python makeformat.py --inputfile input.txt \
                        --normailzationout output.csv \
                        --config config.yaml \
                        [--unique_by {Species|Genus|Family|Order|Class|Phylum}]
"""

import argparse
import pandas as pd
import logging
import yaml 
import re
from log_config import setup_logging

def clean_gcf_id(gcf_id):
    """Remove version number from GCF ID (everything after the dot)"""
    return re.sub(r'\.\d+$', '', gcf_id)

#1制作表格
#0                   1         2        3        4       5       6      7            8           9        10     11
#qseqid              sseqid    qlen    slen    qstart  qend    sstart  send       evalue      bitscore length pident
#GCF_000190995.1     K00097    329     329     1       329     1       329     1.9e-182        644.8     329     25

# 加载YAML配置文件的函数
def load_config(config_file):
    with open(config_file, 'r') as file:
        return yaml.safe_load(file)

def makeformat(inputfile, outputfile, logger, config, unique_by=None):

    logger.info('BLASTP makeformat started.')
    # 初始化一个字典来存储唯一的行
    unique_lines = {}

    with open(inputfile, 'r', encoding='utf-8') as fr:
        # 跳过header行
        next(fr)
        for line in fr:
            sp = line.strip().split("\t")
            # 确保e-value满足条件
            if float(sp[8]) < 0.00001:
                # 解析新的格式
                query_parts = sp[0].split('_**_')  # 分割查询序列ID（uniport部分）
                target_parts = sp[1].split('_**_')  # 分割目标序列ID（GCF部分）
                
                # 获取正确的ID部分 - 现在GCF是目标序列（sseqid）
                sseqid = target_parts[0]  # 从目标序列中获取GCF ID
                qseqid = query_parts[0]  # 从查询序列中获取uniport ID
                
                # 构造唯一标识符 - 保持与原来的格式一致
                unique_key = f"{qseqid}_**_{sseqid}"  # qseqid_**_sseqid
                # 如果这个唯一标识符还没有被添加到字典中，则添加它
                if unique_key not in unique_lines:
                    unique_lines[unique_key] = line.strip()  # 或者其他需要保留的数据形式

    # 接下来的部分处理DataFrame的创建
    index_list = []
    column_list = []
    index_columns = {}
    for key in unique_lines.keys():
        sp = key.split('_**_')  # K00097_**_GCF_000190995.1
        sseqid, qseqid = sp[0], sp[1]
        index_list.append(sseqid)
        column_list.append(qseqid)
        if sseqid not in index_columns:
            index_columns[sseqid] = []
        index_columns[sseqid].append(qseqid)

    df = pd.DataFrame('', index=list(set(index_list)), columns=list(set(column_list)))

    for index, column in zip(index_list, column_list):
        df.loc[index, column] = 'yes'

    df2 = pd.DataFrame(df.values.T, index=df.columns, columns=df.index)
    # 创建存储所有分类水平信息的字典
    taxonomy_info = {}
    medusa_annotation = config['medusa-annotation']
    # 定义分类水平的列名
    taxonomy_columns = ['strain_Name', 'Species', 'Genus', 'Family', 'Order', 'Class', 'Phylum']
    
    with open(medusa_annotation, 'r', encoding='utf-8') as fr:
        header = next(fr).strip().split(',')  # 读取并保存header
        # 获取每个分类水平在CSV中的索引位置
        col_indices = {col: header.index(col) for col in taxonomy_columns}
        
        for line in fr:
            sp1 = line.strip().split(',')
            if len(sp1) > max(col_indices.values()):
                # 清理GCF编号（移除版本号）
                cleaned_gcf = clean_gcf_id(sp1[1])
                if cleaned_gcf not in taxonomy_info:
                    # 为每个GCF存储所有分类水平的信息
                    taxonomy_info[cleaned_gcf] = {col: sp1[col_indices[col]] for col in taxonomy_columns}
    
    # 为每个分类水平创建新列
    for col in taxonomy_columns:
        df2[col] = ''
        for index, row in df2.iterrows():
            # 清理index中的GCF编号
            cleaned_index = clean_gcf_id(index)
            if cleaned_index in taxonomy_info:
                df2.loc[index, col] = taxonomy_info[cleaned_index][col]
    
    # 计算sum列
    df2["sum"] = (df2 == "yes").sum(axis=1)
    
    # 将index转换为列
    df2.reset_index(inplace=True)
    df2.rename(columns={'index': 'GCF'}, inplace=True)
    
    # 设置列的顺序
    column_order = ['GCF', 'strain_Name', 'Species', 'Genus', 'Family', 'Order', 'Class', 'Phylum', 'sum']
    
    # 获取其他列（非标准列，即'yes'数据列）
    other_columns = [col for col in df2.columns if col not in column_order]
    
    # 合并列顺序
    final_column_order = column_order + other_columns
    
    # 重新排序列
    df2 = df2[final_column_order]
    
    # 如果指定了unique_by参数，进行去重处理
    if unique_by:
        logger.info(f'Deduplicating based on {unique_by} level...')
        # 根据指定的分类水平和sum值进行排序
        df2 = df2.sort_values(['sum'], ascending=False)
        # 保留每个分类水平中sum最大的记录
        df2 = df2.drop_duplicates(subset=[unique_by], keep='first')
        logger.info(f'After deduplication, {len(df2)} records remained.')
    
    # 保存到CSV，不再需要index=True因为GCF已经是普通列了
    df2.to_csv(outputfile, sep='\t', index=False)
    logger.info('BLASTP makeformat end.')
    logging.info(f"File format saved to {outputfile}")

def main():
    parser = argparse.ArgumentParser(description="Process and transform biological data into a tabular format.")
    parser.add_argument("--inputfile", type=str, help="The path to the normailzation  input file.")
    parser.add_argument("--normailzationout", type=str, help="The path to the  output file.")
    parser.add_argument('--log_file', '-l', type=str, help='Log file path (optional)')
    parser.add_argument('--config', type=str, default='config.yaml', help='Configuration file path (default: config.yaml)')
    parser.add_argument('--unique_by', type=str, choices=['Species', 'Genus', 'Family', 'Order', 'Class', 'Phylum'],
                        help='Specify the taxonomic level for deduplication. Will keep the entry with highest sum for each unique value.')
    args = parser.parse_args()

    # 加载配置
    config = load_config(args.config)

    # 配置日志并获取 logger 实例
    logger = setup_logging(args.log_file)

    makeformat(args.inputfile, args.normailzationout, logger, config, args.unique_by)

if __name__ == "__main__":
    main()
