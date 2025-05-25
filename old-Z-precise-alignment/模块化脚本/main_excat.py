#!/usr/bin/env python3

import argparse
import yaml
import shutil
from Multi_site_format import process_Multifiles,split_ids,write_to_file
from position_sequence_matcher import parse_position_arg, find_sequences_with_rf
from sequence_pattern_extractor import read_names_from_file, find_matching_sequences, save_matches_to_file
from excat_site_protein import process_site_sequences
from excat_site_cds import process_cds_fasta
from FileOperation import remove_file
from species_checker_format import SpeciesChecker
from log_config import setup_logging

###日志
###修改cds输出格式，并且要注意会有没有找到的 ：No match found for protein ID WP_040344983.1, genome ID GCF_000165795.1_ASM16579v1, position 364 
###还没想好怎么处理，在excat_site_cds的模块中  （完成）
###规范化输出文件，该删除的删除   （完成）
###异常处理！不是所有人都看代码，如果出错要让大家知道错误发生在什么上面
###


###新增规则：25：R表示在25位置为R，25！R表示在25位置不为R
###新增规则：输出结果表示为选中位点存在/不存在两种情况
###python src/main_excat.py --input clu-extract-out.fasta --position 25:R 28:R --extractprotein blast_normailzation_excate_protein.fasta --extractcds blast_normailzation_excate_cds.fasta --site_protein_output site-protein.fasta --site_cds_output site-cds.fasta --formatoutput 111 --config ./src/config.yaml

# 加载YAML配置文件的函数
def load_config(config_file):
    with open(config_file, 'r') as file:
        return yaml.safe_load(file)
    
def main():
    parser = argparse.ArgumentParser(description="Process protein sequences")
    
    # 创建参数分组
    input_group = parser.add_argument_group('Input options')
    output_group = parser.add_argument_group('Output options')
    config_group = parser.add_argument_group('Configuration')

    #输入参数
    input_group.add_argument("--input",'-i', required=True, help="Input ClustalO FASTA file.", metavar="FASTA")
    input_group.add_argument("--position",required=True, nargs='+', type=str, help="Positions and corresponding amino acids.", metavar="POS:AA")
    input_group.add_argument('--extractprotein', required=True, type=str, help='Output file for complete protein sequences.', metavar='PROTEIN_FILE')
    input_group.add_argument('--extractcds', required=True, type=str, help='Output file for complete CDS sequences.', metavar='CDS_FILE') 
    #输出参数
    output_group.add_argument("--site_protein_output", required=True,  type=str, help="Output file for matched protein sequences.", metavar="PROTEIN_OUT")
    output_group.add_argument("--site_cds_output", required=True,  type=str, help="Output file for matched CDS sequences.", metavar="CDS_OUT")
    output_group.add_argument('--formatoutput', type=str, required=True, help='Path to output file')
    # 配置选项
    config_group.add_argument('--config', type=str, required=True, default='config.yaml', help='Configuration file path (default: config.yaml)', metavar='CONFIG_FILE')
    config_group.add_argument('--log_file', '-l', help='Log file path (optional)')

    args = parser.parse_args()

    #临时文件
    tmp_position_swquence = "tmp_output_original.txt"
    tmp_matchesq =  "tmp_matches.fasta"
    tmp_multifiles = "tmp_multifiles.csv"
    tmp_table = "tmp_table.txt"

    # 加载配置
    config = load_config(args.config)

    #日志、cfg文件模块加载
    logger = setup_logging(args.log_file)

    # 定义文件路径
    fasta_file = args.input
    #input_names_file = args.input_file
    protein_db_file = args.extractprotein
    #output_file_path = args.output_file 

    #  config传入文件
    format_path = config['medusa-annotation']

    # 输入序列位置信息  根据输入位点（25：R） 有这个位置的序列的下10个氨基酸
    # >WP_011802784.1 RYSRTRDCFE
    positions = parse_position_arg(args.position)
    # 获取所有满足条件的序列的ID以及相应的锚点序列信息
    all_positions = parse_position_arg(args.position)
    # 对每个位置进行处理
    for pos, (aa, negation) in all_positions.items():
        # 每次循环中只使用一个位置的信息
        positions = {pos: (aa, negation)}        
        # 获取所有满足条件的序列的ID以及相应的锚点序列信息
        rf_sequences = find_sequences_with_rf(fasta_file, positions)        
        # 保存到文件，为每个位置创建一个新的输出文件
        output_file = f"tmp_output_{pos}.txt"
        with open(output_file, 'w') as out_file:
            for sequence_id, anchor_seq in rf_sequences:
                out_file.write(f">{sequence_id} {anchor_seq}\n")
    #对于all_positions，也需要
    rf_sequences = find_sequences_with_rf(fasta_file, all_positions) 
    # 保存到文件
    with open(tmp_position_swquence, 'w') as out_file:
        for sequence_id, anchor_seq in rf_sequences:
            out_file.write(f">{sequence_id} {anchor_seq}\n")
    #print(f"结果已保存到{tmp_position_swquence}文件中")
    # 读取输入序列的id
    names = read_names_from_file(tmp_position_swquence)
    # 根据names信息在蛋白质的结果里面找到对应的序列，做成表格
    matches = find_matching_sequences(protein_db_file, names)
    # 将所有匹配的序列保存到文件中
    save_matches_to_file(matches, args.site_protein_output)#
    #根据位点蛋白序列提取核苷酸位点信息  WP_003018763.1,GCF_024638115.1_ASM2463811v1,383 
    process_site_sequences(tmp_position_swquence , args.extractprotein , tmp_matchesq)
    #根据位点核苷酸序列提取
    process_cds_fasta(args.extractcds, tmp_matchesq, args.site_cds_output)
    # 制作表格
    checker = SpeciesChecker(args.extractprotein, args.site_protein_output, format_path, args.formatoutput, logger)
    checker.generate_report()
    #如果是多位点文件，增加表格中位点的信息
    #gcf_ids, wp_ids = split_ids(args.extractprotein)
    #write_to_file(gcf_ids, wp_ids, tmp_table)
    process_Multifiles(args.formatoutput, tmp_multifiles, args.position, tmp_table)
    shutil.move(tmp_multifiles, args.formatoutput)


    #移除中间文件
    #remove_file(tmp_matchesq)
    #remove_file(tmp_position_swquence)
    #制作表格
    


if __name__ == "__main__":

    main()
