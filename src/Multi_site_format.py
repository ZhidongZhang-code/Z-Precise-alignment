import time
import pandas as pd
import re
import argparse
import os

def split_ids(filename):
    gcf_ids = []
    wp_ids = []
    with open(filename, 'r') as file:
        for line in file:
            if line.startswith('>'):
                # 提取 GCF ID
                gcf_match = re.search(r'(GCF_\d+\.\d+)', line)
                # 提取 WP ID
                wp_match = re.search(r'(WP_\d+\.\d+)', line)
                if gcf_match and wp_match:
                    gcf_ids.append(gcf_match.group(1))
                    wp_ids.append(wp_match.group(1))
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
    result = []
    try:
        with open(fasta_id_path, 'r') as file:
            for line in file:
                if '\t' in line:
                    gcf_id, pid = line.strip().split('\t')
                    if pid == protein_id:
                        result.append(gcf_id)
    except FileNotFoundError:
        print(f"警告: 找不到文件 {fasta_id_path}")
    return result

def process_Multifiles(input_csv, output_csv, position_aa_pairs, fasta_id_path, logger=None):
    """
    根据指定的位置和氨基酸配对处理文件。
    Parameters:
    - input_csv: 输入CSV文件的路径。
    - output_csv: 输出CSV文件的路径。
    - position_aa_pairs: 位置和氨基酸对的列表，格式为POS:AA。
    - fasta_id_path: 包含GCF_ID和蛋白ID对应关系的文件路径。
    - logger: 日志记录器（可选）
    """
    if logger:
        logger.info("开始处理多位点文件...")
    
    pois_list = [pair.split(':')[0] for pair in position_aa_pairs]
    
    try:
        fileA = pd.read_csv(input_csv, sep='\t')
        if logger:
            logger.info(f"成功读取输入文件: {input_csv}")
    except Exception as e:
        error_msg = f"读取输入文件失败 {input_csv}: {e}"
        if logger:
            logger.error(error_msg)
        raise Exception(error_msg)
    
    for pois in pois_list:
        tmp_file = f'tmp_output_{pois}.txt'
        
        # 检查临时文件是否存在
        if not os.path.exists(tmp_file):
            error_msg = f"找不到位点 {pois} 的临时文件: {tmp_file}"
            if logger:
                logger.error(error_msg)
            raise FileNotFoundError(error_msg)
        
        fileA[pois] = 'no'  # 初始化标记列为'no'
        
        try:
            with open(tmp_file, 'r') as f:
                protein_ids = []
                for line in f:
                    line = line.strip()
                    if line.startswith('>'):
                        # 提取蛋白质ID
                        protein_id = line.split('>')[1].split(' ')[0]
                        protein_ids.append(protein_id)
            
            if logger:
                logger.info(f"位点 {pois}: 找到 {len(protein_ids)} 个蛋白质ID")
            
            # 处理每个蛋白质ID
            for protein_id in protein_ids:
                gcf_ids = find_gcf_id_for_protein(protein_id, fasta_id_path)
                if gcf_ids:
                    # 更新对应的行
                    mask = fileA['GCF_ID'].isin(gcf_ids)
                    fileA.loc[mask, pois] = f'yes:{protein_id}'
                    if logger:
                        logger.debug(f"更新了 {mask.sum()} 行，蛋白质ID: {protein_id}")
                        
        except Exception as e:
            error_msg = f"处理位点 {pois} 时出错: {e}"
            if logger:
                logger.error(error_msg)
            raise Exception(error_msg)
    
    try:
        fileA.to_csv(output_csv, index=False)
        if logger:
            logger.info(f"结果已保存到: {output_csv}")
    except Exception as e:
        error_msg = f"保存输出文件失败 {output_csv}: {e}"
        if logger:
            logger.error(error_msg)
        raise Exception(error_msg)

if __name__ == '__main__':
    # 设置argparse来解析命令行参数
    parser = argparse.ArgumentParser(description='Process tmp_output files based on specified positions.')
    parser.add_argument('--input_csv', type=str, required=True, help='Path to the input CSV file')
    parser.add_argument('--output_csv', type=str, required=True, help='Path to the output CSV file')
    parser.add_argument("--position", required=True, nargs='+', type=str, help="Positions and corresponding amino acids.", metavar="POS:AA")
    parser.add_argument('--extractprotein', required=True, type=str, help='Output file for complete protein sequences.', metavar='PROTEIN_FILE')
    
    # 临时文件
    tmp_table = "tmp_table.txt"
    args = parser.parse_args()
    
    try:
        gcf_ids, wp_ids = split_ids(args.extractprotein)
        write_to_file(gcf_ids, wp_ids, tmp_table)
        # 调用主函数处理文件
        process_Multifiles(args.input_csv, args.output_csv, args.position, tmp_table)
        print("处理完成！")
    except Exception as e:
        print(f"处理过程中发生错误: {e}")
        raise
    finally:
        # 清理临时文件
        if os.path.exists(tmp_table):
            os.remove(tmp_table)

