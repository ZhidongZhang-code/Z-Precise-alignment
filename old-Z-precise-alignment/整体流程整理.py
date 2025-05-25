import os
import sys
import re
import logging
import argparse
import subprocess
import pandas as pd
from concurrent.futures import ThreadPoolExecutor
from threading import Lock

###编程日志###
###先不用rename，序列怎么输入还没想好
###还需要一个用来保存log的参数
###脚本模块化
###模块继承，多态  使用多继承来选择使用blast还是diamond
###GUI界面


def setup_logging(log_file_path=None, log_level=logging.INFO):
    # 创建一个具有唯一名称的logger，避免全局配置冲突
    logger = logging.getLogger('MyScriptLogger')
    logger.setLevel(log_level)

    # 清除已存在的handlers，避免重复添加
    if logger.handlers:
        for handler in logger.handlers:
            logger.removeHandler(handler)

    # 根据是否指定了日志文件路径来设置日志处理程序
    if log_file_path:
        file_handler = logging.FileHandler(log_file_path)
    else:
        file_handler = logging.StreamHandler()

    file_handler.setLevel(log_level)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    return logger

def rename(input_files, output_file_path):

    for input_file_path in input_files:
        if not input_file_path.endswith(".fasta"):
            logging.warning(f"Ignoring non-fasta file: {input_file_path}")
            continue

        if not os.path.isfile(input_file_path):
            logging.error(f"File not found: {input_file_path}")
            continue

        file_name = os.path.splitext(os.path.basename(input_file_path))[0]
        modified_lines = []

        with open(input_file_path, 'r') as input_file:
            for line in input_file:
                line = line.strip()
                if line.startswith('>') and '_**_' not in line:
                    modified_lines.append(f'>{file_name}_**_{line[1:]}')
                else:
                    modified_lines.append(line)

        modified_content = '\n'.join(modified_lines)

        with open(output_file_path, 'a') as output_file:
            output_file.write(modified_content)

        logging.info(f"Processed: {input_file_path}")

    logging.info(f"All files have been processed. Output saved to: {output_file_path}")


def run_blastp( query_file, output_file, num_threads,logger):
    ###-outfmt "6 qseqid sseqid qlen slen qstart qend sstart send evalue bitscore length pident" \
    ###-evalue 1e-5 -max_target_seqs 50000 -num_threads 8
    logger.info(f"Starting BLASTp for query file: {query_file} with {num_threads} threads")
    command = [
        '/public/home/TonyWuLab/zhangzhd/anaconda3/envs/blast/bin/blastp',
        '-db', '/public/group_share_data/TonyWuLab/zzd/db/medusa-2023/all-hum-medusa-db-genmoic-protein',
        '-query', query_file,
        '-out', output_file,
        '-outfmt', '6 qseqid sseqid qlen slen qstart qend sstart send evalue bitscore length pident',
        '-evalue', '1e-5',
        '-max_target_seqs', '5000',
        '-num_threads', str(num_threads)
    ]

    try:
        subprocess.run(command, check=True)
        logger.info("Blastp command executed successfully.")
    except subprocess.CalledProcessError as e:
        logger.error(f"Error executing Blastp command: {e}")
        sys.exit(1)  # 出错时退出程序

def filenormalization(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        #####截取合适的cov与identity的值
        new_lines = []
        for line in infile:
            line = line.strip()
            if re.match(r'^\w+', line):
                test = line.split('\t')
                if (abs(float(test[5]) - float(test[4]) + 1) / int(test[2]) >= 0.7 and      #（qstar-qend）/qlen
                        abs(float(test[7]) - float(test[6]) + 1) / int(test[3]) >= 0.7 and      #（sstar-send）sqlen
                        float(test[11]) >= 30):      #identity
                    new_lines.append(line + '\n')
        ######对于_**_内容，取重要的部分
        #K00097_**_eco:b0052     GCF_000190995.1_ASM19099v1_**_WP_000241242.1
    #GCF_000190995.1         K00097
        for line in new_lines:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                col1, col2 = parts[0], parts[1]
                col1 = col1.split("_**_")[0]
                col2 = "_".join(col2.split("_", 2)[:2])
                new_line = "{}\t{}\t{}\n".format(col2, col1, '\t'.join(parts[2:]))
                outfile.write(new_line)
        logging.info("File normalization processed successfully.")

def makeformat(inputfile, outputfile,logger):

    logger.info(f"Starting makeformat for BLASTp result")
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
    logger.info(f"File format saved to {outputfile}")

# 初始化锁
write_lock_cds = Lock()
write_lock_protein = Lock()


def delete_output_files():
    """初始化输出文件"""
    if os.path.exists("output_cds.fasta"):
        os.remove("output_cds.fasta")
    if os.path.exists("output_protein.fasta"):
        os.remove("output_protein.fasta")

def findextractcdsprotein(input_file, directory, num_threads,logger):
    """
    处理输入文件，提取符合条件的标识符，并使用线程池并行处理序列文件。

    参数：
        input_file：包含数据的输入文件路径
        directory：包含 GCF 编号文件夹的目录路径
        num_threads：要使用的线程数
    """

    logger.info("Starting to extract and process find-extract-cds-protein.")
    identifiers = set()  # 用于存储唯一的标识符

    with open(input_file, 'r') as infile:
        for line in infile:
            line = line.strip()
            if re.match(r'^\w+', line):
                columns = line.split('\t')
                # 根据特定条件筛选标识符
                if (abs(float(columns[5]) - float(columns[4]) + 1) / int(columns[2]) >= 0.7 and      #（qstar-qend）/qlen
                        abs(float(columns[7]) - float(columns[6]) + 1) / int(columns[3]) >= 0.7 and      #（sstar-send）sqlen
                        float(columns[11]) >= 30):      #identity
                    identifier = columns[1]
                    identifiers.add(identifier)

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
                    futures.append(executor.submit(process_sequence, identifier, cds_file_path, "output_cds.fasta", write_lock_cds))
                    futures.append(executor.submit(process_sequence, identifier, protein_file_path, "output_protein.fasta", write_lock_protein))
            else:
                logger.warning(f"Folder {folder_path} for identifier {identifier} does not exist.")
    for future in futures:
        future.result()  # This line ensures all futures have completed.
    logger.info("Finished processing all sequences.")
    

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
    wp_identifier_search = re.search(r'_\*\*_(.+?)$', identifier)  # 修改此处的正则表达式
    
    # 确保搜索到了GCF和_**_之后标识符
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
            print(f"Error reading file {file_path}: {e}")
    else:
        print(f"Identifier format mismatch in {identifier}")
    return header, sequence

def filter_duplicate_proteins(input_file, output_file):
    """
    从输入文件中筛选出唯一的蛋白质，并写入输出文件。

    Args:
        input_file (str): 输入文件路径。
        output_file (str): 输出文件路径。
    """
    unique_proteins = {}
    with open(input_file, 'r') as f:
        lines = f.readlines()
        for line in lines:
            if line.startswith(">"):
                protein_id = line.split()[0]
                wp_id = protein_id.split("_**_")[-1].split("\t")[0]  # 获取蛋白质ID中的特定部分
                if wp_id not in unique_proteins:
                    unique_proteins[wp_id] = True

    with open(output_file, 'w') as f:
        l = 0
        for line in lines:
            if line.startswith(">"):
                protein_id = line.split()[0]
                wp_id = protein_id.split("_**_")[-1].split("\t")[0]  # 获取蛋白质ID中的特定部分
                if wp_id in unique_proteins:
                    l = 1
                    f.write(">" + wp_id + "\n")  # 只保留 WP 部分
                    del unique_proteins[wp_id]
            elif l == 1 :
                f.write(line)
                l = 0

def merge_fasta_files(file1, file2, output_file):
    """
    将两个 FASTA 格式的文件合并成一个输出文件。

    Args:
        file1 (str): 第一个输入文件路径。
        file2 (str): 第二个输入文件路径。
        output_file (str): 输出文件路径。
    """
    with open(output_file, 'w') as output:
        with open(file1, 'r') as f1:
            output.write(f1.read())
        with open(file2, 'r') as f2:
            output.write(f2.read())

def run_clustalo(input_file, output_file, num_threads ,logger):
    """
    运行 Clustal Omega 进行多序列比对。

    Args:
        input_file (str): 输入文件路径。
        output_file (str): 输出文件路径。
        num_threads (int): 使用的线程数，默认为 8。
    """
    logger.info(f"Starting clustalo for query file: {input_file} with {num_threads} threads")
    command = [
        '/public/home/TonyWuLab/zhangzhd/anaconda3/envs/clustalo/bin/clustalo',
        '--infile', input_file,
        '-o', output_file,
        '--threads', str(num_threads),
        '--force'
    ]    
    
    try:
        subprocess.run(command, check=True)
        logger.info("clustalo command executed successfully.")
    except subprocess.CalledProcessError as e:
        logger.error(f"Error executing clustalo command: {e}")
        sys.exit(1)  # 出错时退出程序

def remove_file(file):
    """
    删除指定文件。

    Args:
        file (str): 要删除的文件路径。
    """
    if os.path.exists(file):
        os.remove(file)

def find_start_end_positions(alignment, sequence_id):
    """
    查找序列在多序列比对结果中的起始和结束位置。

    Args:
        alignment (dict): 多序列比对结果，格式为序列ID到序列的映射字典。
        sequence_id (str): 要查找的序列ID。

    Returns:
        tuple: 起始位置和结束位置的元组。
    """
    start_position = None
    end_position = None

    # 找到序列的起始位置
    for i in range(len(alignment[sequence_id])):
        if alignment[sequence_id][i] != '-':
            start_position = i
            break

    # 找到序列的结束位置
    for i in range(len(alignment[sequence_id]) - 1, -1, -1):
        if alignment[sequence_id][i] != '-':
            end_position = i
            break

    return start_position, end_position

def find_end_position(alignment, sequence_id):
    """
    查找序列在多序列比对结果中的结束位置。

    Args:
        alignment (dict): 多序列比对结果，格式为序列ID到序列的映射字典。
        sequence_id (str): 要查找的序列ID。

    Returns:
        int: 结束位置。
    """
    sequence = alignment[sequence_id]
    end_position = None
    for i in range(len(sequence) - 1, -1, -1):
        if sequence[i] != '-':
            end_position = i + 3  # 找到最后一个非'-'的位置后，往后数三位作为结束位置
            break
    return end_position

def trim_alignment(alignment, start_position, end_position, target_sequence_id, padding=3):
    """
    裁剪多序列比对结果，只保留指定范围内的序列。

    Args:
        alignment (dict): 多序列比对结果，格式为序列ID到序列的映射字典。
        start_position (int): 起始位置。
        end_position (int): 结束位置。
        target_sequence_id (str): 要裁剪的序列ID。
        padding (int): 裁剪时的填充大小，默认为 3。

    Returns:
        dict: 裁剪后的多序列比对结果。
    """
    trimmed_alignment = {}
    for seq_id, seq in alignment.items():
        trimmed_seq = seq[start_position - padding:end_position + padding + 1]  # 加一是因为切片时右边界不包含
        trimmed_alignment[seq_id] = trimmed_seq
    # 将指定序列移动到字典的第一个位置
    target_sequence = trimmed_alignment.pop(target_sequence_id)
    trimmed_alignment = {target_sequence_id: target_sequence, **trimmed_alignment}
    return trimmed_alignment

def read_alignment_from_file(file):
    """
    从文件中读取多序列比对结果。

    Args:
        file (str): 输入文件路径。

    Returns:
        dict: 多序列比对结果，格式为序列ID到序列的映射字典。
    """
    alignment = {}
    with open(file, 'r') as f:
        current_seq_id = None
        current_seq = ''
        for line in f:
            if line.startswith('>'):
                if current_seq_id:
                    alignment[current_seq_id] = current_seq
                current_seq_id = line.strip()
                current_seq = ''
            else:
                current_seq += line.strip()
        if current_seq_id:
            alignment[current_seq_id] = current_seq
    return alignment

def write_alignment_to_file(alignment, output_file):
    """
    将多序列比对结果写入文件。

    Args:
        alignment (dict): 多序列比对结果，格式为序列ID到序列的映射字典。
        output_file (str): 输出文件路径。
    """
    with open(output_file, 'w') as f:
        for seq_id, seq in alignment.items():
            f.write(seq_id + '\n')
            f.write(seq + '\n')

def main():
    parser = argparse.ArgumentParser(description='处理输入的FASTA文件并修改标题。')
    parser.add_argument('--input', type=str, required=True, help='输入FASTA文件路径')
    parser.add_argument('--blastdatabase', '-db', type=str, default='/public/group_share_data/TonyWuLab/zzd/db/medusa-2023/all-hum-medusa-db-genmoic-protein', help='Path to the BLAST database')
    parser.add_argument('--num_threads', '-t', type=int, default=8, help='Number of threads (default: 8)')
    parser.add_argument('--alignmentout', '-o', type=str, required=True, help='Path to the output file format')
    parser.add_argument("--GCFdirectory", default="/public/group_share_data/TonyWuLab/zzd/db/medusa-taxon-genmoic/hum-medusa-db-all",help="Directory containing GCF numbered folders(/public/group_share_data/TonyWuLab/zzd/db/medusa-taxon-genmoic/hum-medusa-db-all).")
    parser.add_argument("--target_sequence_file", type=str, required=True, help="File containing the target sequence(150nt).")
    parser.add_argument('--log', type=str, help='Log file path (optional)')
    parser.add_argument('--log_level', type=str, default='INFO', choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'], help='Set the log level (default: INFO)')
    args = parser.parse_args()

#    output_file_path_rename = 'tmp-rename.fasta' ###中间结果需要删除
    output_file_path_blast = 'tmp-rename-blastp.csv' ###中间结果需要删除
    output_file_path_normalization = 'tmp-normalization.csv' ###中间结果---有多少的细菌含有目标基因
    output_file_path_alignmentout = args.alignmentout
    temp_filtered_file = "tmp_filtered.fasta" ###中间结果需要删除
    merged_file = "merged.fasta"   ###中间结果需要删除
    output_protein_file = "output_protein.fasta"
    clustal_output_file = "clu-out.fasta" ###clustalo比对结果文件
    input_files = args.input
    target_sequence_file = args.target_sequence_file
    # 读取目标序列文件，获取序列ID
    with open(args.target_sequence_file, 'r') as f:
        target_sequence_id = f.readline().strip()

#    for input_file_path in input_files:
#        if not input_file_path.endswith(".fasta"):
#            print(f"忽略非fasta文件: {input_file_path}")
#            continue

#        if not os.path.isfile(input_file_path):
#            print(f"文件不存在: {input_file_path}")
#            continue

#        rename(input_file_path, output_file_path_rename)
#        print(f"已处理: {input_file_path}")
    
    # 将字符串日志级别转换为logging模块的级别
    numeric_level = getattr(logging, args.log_level.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError(f'Invalid log level: {args.log_level}')

    # 设置日志
    logger = setup_logging(log_file_path=args.log, log_level=numeric_level)

    # 使用logger来记录日志
    logger.info("Starting script")

#    run_blastp(input_files, output_file_path_blast, args.num_threads,logger) ##输出的中间文件是output_file_path_blast
    filenormalization(output_file_path_blast, output_file_path_normalization)
    makeformat(output_file_path_normalization, output_file_path_alignmentout,logger)
    findextractcdsprotein(output_file_path_blast, args.GCFdirectory, args.num_threads,logger)
    remove_file(merged_file)  # 删除之前的 merged.fasta 文件
    filter_duplicate_proteins(output_protein_file, temp_filtered_file)  # 过滤重复的蛋白质ID并写入临时文件
    merge_fasta_files(temp_filtered_file, target_sequence_file, merged_file)
    run_clustalo(merged_file, clustal_output_file, args.num_threads,logger)
    remove_file(merged_file)  # 删除临时合并的文件
    remove_file(temp_filtered_file)  # 删除临时过滤的文件

    # 从比对结果中裁剪指定序列并写入文件
    alignment = read_alignment_from_file(clustal_output_file)
    start_position, end_position = find_start_end_positions(alignment, target_sequence_id)
    end_position = find_end_position(alignment, target_sequence_id)
    trimmed_alignment = trim_alignment(alignment, start_position, end_position, target_sequence_id)
    write_alignment_to_file(trimmed_alignment, clustal_output_file)
    
#    print(f"已处理: {output_file_path}")

if __name__ == "__main__":
    main()