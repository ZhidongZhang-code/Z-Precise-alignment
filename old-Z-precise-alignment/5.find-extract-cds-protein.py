import os
import re
import argparse
from concurrent.futures import ThreadPoolExecutor
from threading import Lock

# 初始化锁
write_lock_cds = Lock()
write_lock_protein = Lock()

def delete_output_files():
    """初始化输出文件"""
    if os.path.exists("output_cds.fasta"):
        os.remove("output_cds.fasta")
    if os.path.exists("output_protein.fasta"):
        os.remove("output_protein.fasta")

def findextractcdsprotein(input_file, directory, num_threads):
    """
    处理输入文件，提取符合条件的标识符，并使用线程池并行处理序列文件。

    参数：
        input_file：包含数据的输入文件路径
        directory：包含 GCF 编号文件夹的目录路径
        num_threads：要使用的线程数
    """
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
            print(f"Error 
                   file {file_path}: {e}")
    else:
        print(f"Identifier format mismatch in {identifier}")
    return header, sequence

def main():
    parser = argparse.ArgumentParser(description="python3 5.find-extract-cds-protein.py input_file , Output two files output_cds.fasta and output_protein.fasta ")
    parser.add_argument("input_file", help="Input file containing data.")
    parser.add_argument("--directory", default="/public/group_share_data/TonyWuLab/zzd/db/medusa-taxon-genmoic/hum-medusa-db-all",
                        help="Directory containing GCF numbered folders.")
    parser.add_argument("--num_threads", type=int, default=8, help="Number of threads to use.")
    args = parser.parse_args()

    delete_output_files()  # 删除输出文件

    findextractcdsprotein(args.input_file, args.directory, args.num_threads)

if __name__ == "__main__":
    main()
