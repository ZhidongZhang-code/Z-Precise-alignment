import time
import pandas as pd
import argparse

def split_ids(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()

    gcf_ids = []
    wp_ids = []

    for line in lines:
        part1 = line.split(' ')[0]
        parts = part1.strip().split('_')
        if len(parts) >= 4:
            gcf_id = '_'.join(parts[:2])[1:]
            wp_id = '_'.join(parts[-2:])
            gcf_ids.append(gcf_id)
            wp_ids.append(wp_id)

    return gcf_ids, wp_ids

def write_to_file(gcf_ids, wp_ids, tmp_table):
    with open(tmp_table, 'w') as file:
        for gcf_id, wp_id in zip(gcf_ids, wp_ids):
            file.write(gcf_id + '\t' + wp_id + '\n')



def find_gcf_id_for_protein(protein_id, fasta_id_path):
    """
    查找给定蛋白ID对应的GCF_ID。

    Parameters:
    - protein_id: 蛋白ID。
    - fasta_id_path: 包含GCF_ID和蛋白ID对应关系的文件路径。

    Returns:
    - 对应的GCF_ID，如果找不到则返回None。
    """
    result=[]
    with open(fasta_id_path, 'r') as file:
        for line in file:
            gcf_id, pid = line.strip().split('\t')
            if pid == protein_id:
                result.append(gcf_id)
    return result



def process_Multifiles(input_csv, output_csv, position_aa_pairs, fasta_id_path):
    """
    根据指定的位置和氨基酸配对处理文件。

    Parameters:
    - input_csv: 输入CSV文件的路径。
    - output_csv: 输出CSV文件的路径。
    - position_aa_pairs: 位置和氨基酸对的列表，格式为POS:AA。
    - fasta_id_path: 包含GCF_ID和蛋白ID对应关系的文件路径。
    """
    pois_list = [pair.split(':')[0] for pair in position_aa_pairs]
    fileA = pd.read_csv(input_csv, sep='\t')
    #print(pois_list)
    for pois in pois_list:
        tmp_file = f'tmp_output_{pois}.txt'
        fileA[pois] = 'no'  # 初始化标记列为'no'

        with open(tmp_file, 'r') as f:
            protein_ids = [line.strip().split('>')[1].split(' ')[0] for line in f]
        for protein_id in protein_ids:
            gcf_ids = find_gcf_id_for_protein(protein_id, fasta_id_path)
            #print("gcf_ids length:",len(gcf_ids))
            fileA.loc[fileA['GCF_ID'].isin(gcf_ids),pois]=f'yes:{protein_id}'

    fileA.to_csv(output_csv, index=False)


if __name__ == '__main__':
    # 设置argparse来解析命令行参数
    parser = argparse.ArgumentParser(description='Process tmp_output files based on specified positions.')
    parser.add_argument('--input_csv', type=str, required=True, help='Path to the input CSV file')
    parser.add_argument('--output_csv', type=str, required=True, help='Path to the output CSV file')
    parser.add_argument("--position", required=True, nargs='+', type=str, help="Positions and corresponding amino acids.", metavar="POS:AA")
    parser.add_argument('--extractprotein', required=True, type=str, help='Output file for complete protein sequences.', metavar='PROTEIN_FILE')

    #临时文件
    tmp_table = "tmp_table.txt"

    args = parser.parse_args()

    gcf_ids, wp_ids = split_ids(args.extractprotein)
    write_to_file(gcf_ids, wp_ids, tmp_table)

    # 调用主函数处理文件
    process_Multifiles(args.input_csv, args.output_csv, args.position, tmp_table)