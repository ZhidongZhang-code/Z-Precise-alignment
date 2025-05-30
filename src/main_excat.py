#!/usr/bin/env python3

import argparse
import yaml
import shutil
import os
import tempfile
import atexit
from Multi_site_format import process_Multifiles,split_ids,write_to_file
from position_sequence_matcher import parse_position_arg, find_sequences_with_rf
from sequence_pattern_extractor import read_names_from_file, find_matching_sequences, save_matches_to_file
from excat_site_protein import process_site_sequences
from excat_site_cds import process_cds_fasta
from FileOperation import remove_file
from species_checker_format import SpeciesChecker,OneSpeciesChecker
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

# 全局变量存储临时文件列表
temp_files = []

def cleanup_temp_files():
    """清理临时文件（可选择性禁用）"""
    # 如果不想删除临时文件，直接返回
    for temp_file in temp_files:
        if os.path.exists(temp_file):
            try:
                os.remove(temp_file)
                print(f"已删除临时文件: {temp_file}")
            except Exception as e:
                print(f"删除临时文件失败 {temp_file}: {e}")

def create_temp_file_in_current_dir(filename):
    """在当前目录创建指定名称的临时文件"""
    temp_path = os.path.join(os.getcwd(), filename)
    with open(temp_path, 'w') as f:
        pass
    temp_files.append(temp_path)
    return temp_path

# 加载YAML配置文件的函数
def load_config(config_file):
    with open(config_file, 'r') as file:
        return yaml.safe_load(file)

def main():
    # 注册清理函数，确保程序退出时清理临时文件
    atexit.register(cleanup_temp_files)

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

    #输出参数 - 恢复--site_protein_output参数
    output_group.add_argument("--site_protein_output", required=True, type=str, help="Output file for matched protein sequences.", metavar="PROTEIN_OUT")
    output_group.add_argument("--site_cds_output", required=True,  type=str, help="Output file for matched CDS sequences.", metavar="CDS_OUT")
    output_group.add_argument('--formatoutput', type=str, required=False, help='Path to output file')

    # 配置选项
    config_group.add_argument('--config', type=str, required=True, default='config.yaml', help='Configuration file path (default: config.yaml)', metavar='CONFIG_FILE')
    config_group.add_argument('--log_file', '-l', help='Log file path (optional)')
    # 新增参数：是否保留临时文件
    config_group.add_argument('--keep_temp', action='store_true', help='Keep temporary files in current directory')

    args = parser.parse_args()

    try:
        # 创建临时文件的方式
        tmp_position_sequence = create_temp_file_in_current_dir("tmp_output_original.txt")
        tmp_matches = create_temp_file_in_current_dir("tmp_matches.fasta")
        tmp_multifiles = create_temp_file_in_current_dir("tmp_multifiles.csv")
        tmp_table = create_temp_file_in_current_dir("tmp_table.txt")

        # 加载配置
        config = load_config(args.config)

        #日志、cfg文件模块加载
        logger = setup_logging(args.log_file)
        logger.info("开始处理序列...")

        # 定义文件路径
        fasta_file = args.input
        protein_db_file = args.extractprotein

        #  config传入文件
        format_path = config['medusa-annotation']

        # 输入序列位置信息  根据输入位点（25：R） 有这个位置的序列的下10个氨基酸
        positions = parse_position_arg(args.position)
        # 获取所有满足条件的序列的ID以及相应的锚点序列信息
        all_positions = parse_position_arg(args.position)

        # 创建位置文件映射字典，用于Multi_site_format模块
        position_files = {}

        # 对每个位置进行处理
        for pos, (aa, negation) in all_positions.items():
            # 每次循环中只使用一个位置的信息
            positions = {pos: (aa, negation)}
            # 获取所有满足条件的序列的ID以及相应的锚点序列信息
            rf_sequences = find_sequences_with_rf(fasta_file, positions)

            # 使用Multi_site_format期望的文件名格式
            output_file = f"tmp_output_{pos}.txt"
            if args.keep_temp:
                # 在当前目录创建位置特定的临时文件
                output_file = os.path.join(os.getcwd(), output_file)

            position_files[pos] = output_file
            temp_files.append(output_file)  # 加入清理列表

            with open(output_file, 'w') as out_file:
                for sequence_id, anchor_seq in rf_sequences:
                    out_file.write(f">{sequence_id} {anchor_seq}\n")
            logger.info(f"位置 {pos} 处理完成，找到 {len(rf_sequences)} 个匹配序列")

        # 对于all_positions，也需要处理
        rf_sequences = find_sequences_with_rf(fasta_file, all_positions)
        # 保存到文件
        with open(tmp_position_sequence, 'w') as out_file:
            for sequence_id, anchor_seq in rf_sequences:
                out_file.write(f">{sequence_id} {anchor_seq}\n")

        logger.info(f"总共找到 {len(rf_sequences)} 个符合所有位置条件的序列")

        # 根据位点蛋白序列提取核苷酸位点信息  WP_003018763.1,GCF_024638115.1_ASM2463811v1,383
        process_site_sequences(tmp_position_sequence, args.extractprotein, tmp_matches)

        # 根据位点核苷酸序列提取，输出到最终文件
        process_cds_fasta(args.extractcds, tmp_matches, args.site_cds_output)
        logger.info(f"CDS序列已保存到: {args.site_cds_output}")

        # 恢复蛋白质序列提取功能 - 输出到用户指定的文件
        names = read_names_from_file(tmp_position_sequence)
        matches = find_matching_sequences(protein_db_file, names)
        save_matches_to_file(matches, args.site_protein_output)
        logger.info(f"蛋白质序列已保存到: {args.site_protein_output}")

        # 制作表格 - 使用用户指定的蛋白输出文件
        if args.formatoutput:
            checker = SpeciesChecker(args.extractprotein, args.site_protein_output, format_path, args.formatoutput, logger)
            checker.generate_report()

            gcf_ids, wp_ids = split_ids(args.extractprotein)
            write_to_file(gcf_ids, wp_ids, tmp_table)
            process_Multifiles(args.formatoutput, tmp_multifiles, args.position, tmp_table, logger)
            shutil.move(tmp_multifiles, args.formatoutput)

            logger.info(f"最终输出文件: {args.formatoutput}")

        logger.info("处理完成！")
        logger.info(f"最终输出文件: {args.formatoutput}")
        logger.info(f"蛋白质序列文件: {args.site_protein_output}")
        logger.info(f"CDS序列文件: {args.site_cds_output}")

        # 如果保留临时文件，打印文件位置信息
        if args.keep_temp:
            logger.info("临时文件已保留在当前目录:")
            for temp_file in temp_files:
                if os.path.exists(temp_file):
                    logger.info(f"  {temp_file}")

    except Exception as e:
        logger.error(f"处理过程中发生错误: {e}")
        raise

    # 根据参数决定是否清理临时文件
    if not args.keep_temp:
        cleanup_temp_files()

if __name__ == "__main__":
    main()


