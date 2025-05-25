import argparse
import pandas as pd
import logging
#1制作表格
#0                   1         2        3        4       5       6      7            8           9        10     11
#qseqid              sseqid    qlen    slen    qstart  qend    sstart  send       evalue      bitscore length pident
#GCF_000190995.1     K00097    329     329     1       329     1       329     1.9e-182        644.8     329     25

def setup_logging(log_file_path):
    # 创建一个logger
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)

    if log_file_path:
        # 如果指定了日志文件路径，只创建用于输出到文件的handler
        file_handler = logging.FileHandler(log_file_path)
        file_handler.setLevel(logging.INFO)
        file_formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        file_handler.setFormatter(file_formatter)
        logger.addHandler(file_handler)
    else:
        # 否则，只创建用于输出到控制台的handler
        console_handler = logging.StreamHandler()
        console_handler.setLevel(logging.INFO)
        console_formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        console_handler.setFormatter(console_formatter)
        logger.addHandler(console_handler)

        
def makeformat(inputfile, outputfile):
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

    df = pd.DataFrame('', index=set(index_list), columns=set(column_list))

    for index, column in zip(index_list, column_list):
        df.loc[index, column] = 'yes'

    df2 = pd.DataFrame(df.values.T, index=df.columns, columns=df.index)
    tmp = {}
    with open('/public/group_share_data/TonyWuLab/zzd/db/medusa-2023/all-hum-medusa-db-GCF-name.txt', 'r', encoding='utf-8') as fr:
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
    logging.info(f"File format saved to {outputfile}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process and transform biological data into a tabular format.")
    parser.add_argument("inputfile", type=str, help="The path to the input file.")
    parser.add_argument("outputfile", type=str, help="The path to the output file.")
    parser.add_argument('--log', type=str, help='Log file path (optional)')
    args = parser.parse_args()
    makeformat(args.inputfile, args.outputfile)
