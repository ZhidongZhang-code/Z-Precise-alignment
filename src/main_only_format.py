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
    input_group.add_argument("--input",'-i', required=True, help="Input Blast form file.", metavar="Table")
    #输出参数
    output_group.add_argument('--formatoutput', type=str, required=True, help='Path to output file')
    # 配置选项
    config_group.add_argument('--config', type=str, required=True, default='config.yaml', help='Configuration file path (default: config.yaml)', metavar='CONFIG_FILE')
    config_group.add_argument('--log_file', '-l', help='Log file path (optional)')

    args = parser.parse_args()

    #临时文件
    tmp_multifiles = "tmp_multifiles.csv"
    tmp_table = "tmp_table.txt"

    # 加载配置
    config = load_config(args.config)

    #日志、cfg文件模块加载
    logger = setup_logging(args.log_file)

    # 定义文件路径
    #input_names_file = args.input_file
    #output_file_path = args.output_file 

    #  config传入文件
    format_path = config['medusa-annotation']

    # 制作表格
    checker = OneSpeciesChecker(args.input, format_path, args.formatoutput, logger)
    checker.generate_report()
    #如果是多位点文件，增加表格中位点的信息
    #gcf_ids, wp_ids = split_ids(args.extractprotein)
    #write_to_file(gcf_ids, wp_ids, tmp_table)
    #process_Multifiles(args.formatoutput, tmp_multifiles, args.position, tmp_table)
    #shutil.move(tmp_multifiles, args.formatoutput)


    #移除中间文件
    #remove_file(tmp_matchesq)
    #remove_file(tmp_position_swquence)
    #制作表格
    


if __name__ == "__main__":

    main()
