import argparse
import pandas as pd
import logging
import yaml 
from log_config import setup_logging

#1制作表格
#0                   1         2        3        4       5       6      7            8           9        10     11
#qseqid              sseqid    qlen    slen    qstart  qend    sstart  send       evalue      bitscore length pident
#GCF_000190995.1     K00097    329     329     1       329     1       329     1.9e-182        644.8     329     25

# 加载YAML配置文件的函数
def load_config(config_file):
    with open(config_file, 'r') as file:
        return yaml.safe_load(file)

def makeformat(inputfile, outputfile, logger,config):

    logger.info('BLASTP makeformat started.')
    # 初始化一个字典来存储唯一的行
    unique_lines = {}

    with open(inputfile, 'r', encoding='utf-8') as fr:
        for line in fr:
            sp = line.strip().split("\t")
            # 确保e-value满足条件
            if float(sp[8]) < 0.00001:
                # 构造唯一标识符
                unique_key = f"{sp[1]}_**_{sp[0]}"  # sseqid_**_qseqid
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
    tmp = {}
    GCF_annotation = config['GCF_annotation']
    with open(GCF_annotation, 'r', encoding='utf-8') as fr:
        for line in fr:
            sp1 = line.strip().split('\t')
            if sp1[0] not in tmp:
                tmp[sp1[0]] = sp1[1]
    df2['baname'] = ''
    for index, row in df2.iterrows():
        df2.loc[index, "baname"] = tmp.get(index, '')
    c = df2.pop("baname")
    df2.insert(0, 'baname', c)
    df2["sum"] = (df2 == "yes").sum(axis=1)
    d = df2.pop('sum')
    df2.insert(1, 'sum', d)

    # 保存到CSV
    df2.to_csv(outputfile, sep='\t', index=True, header=True)
    logger.info('BLASTP makeformat end.')
    logging.info(f"File format saved to {outputfile}")

def main():
    parser = argparse.ArgumentParser(description="Process and transform biological data into a tabular format.")
    parser.add_argument("--inputfile", type=str, help="The path to the normailzation  input file.")
    parser.add_argument("--normailzationout", type=str, help="The path to the  output file.")
    parser.add_argument('--log_file', '-l', type=str, help='Log file path (optional)')
    parser.add_argument('--config', type=str, default='config.yaml', help='Configuration file path (default: config.yaml)')
    args = parser.parse_args()

    # 加载配置
    config = load_config(args.config)

    # 配置日志并获取 logger 实例
    logger = setup_logging(args.log_file)

    makeformat(args.inputfile, args.normailzationout, logger,config)

if __name__ == "__main__":
    main()