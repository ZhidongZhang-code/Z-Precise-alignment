#!/usr/bin/env python3

import os
import argparse
import csv
import shutil
from blast_module import setup_logging, run_blast, load_config  
from FileOperation import remove_file
from fasta_renamer import rename 
from alignment_normalization import filenormalization
from find_extract_cds_protein import findextractcdsprotein
from makeformat import makeformat
from run_clustalo_process import read_alignment_from_file,find_start_end_positions,find_end_position,trim_alignment,write_alignment_to_file,filter_duplicate_proteins,merge_fasta_files,run_clustalo,remove_file
from hmm_module import run_HMMprcess
from process_hmmresult import excat_hmm_postition, compare_files_and_keep_df1_only
from pfam_module import run_pfam
from process_pfamresult import pfam_compare_files_and_keep_df1_only


###blastout算是一个默认参数，最后需要删除

###日志：
###运行完成后删除中间文件
###软件打包之后在安装的时候需要开启权限，否则无法运行
###每一个模块运行后，要给出运行这个模块的命令行脚本
###软件运行时怎么才能不选择python sh.py而且直接sh.py  完成
###脚本还需要一个例子，在运行-h的时候，可以看到例子  完成
###还是没有解决rename.py的问题，这个脚本是不是应该单独拿出来？
###除了blastp还应该可以支持blastx 由于数据库是蛋白库，所以不能支持blastn    完成
###blastp的参数应该作为一个参数能在脚本中进行选择（cov、id、e值、k）    完成
###如果某一部分出了问题，需要报错在那一部分出的错误 完成
###hmm模块还没有装载好，并且还是考虑这个是否必须使用    完成 
###在config中选择是否使用HMM，默认为使用    完成 
###pfam模块也要加载 (1648个序列8线程跑半小时，有些慢)  完成
###pfam模块不需要每一次都跑一边，因为pfam的输入序列是固定的蛋白数据库中的，所有可以直接先跑完，跑出来pfam的那个表格，然后用blaast的相似序列直接去查找就可以了
###模块的异常处理
###参数列表检查是不是存在必须要有的数据库（参考pfam脚本）
###做GUI，在乌班图这种有图形化界面的linux系统上使用
###2024/10/17 增加了diamond进行快速分析 但是不属于config（因为比属于必选参数），放在alignment.py中作为加快运行速度的选项



__version__ = '1.0'

def main():
    parser = argparse.ArgumentParser(description='Control script to run BLAST analysis')
    parser = argparse.ArgumentParser(description=f'NBrun.py v.{__version__}')

    # 创建参数分组
    input_group = parser.add_argument_group('Input options')
    output_group = parser.add_argument_group('Output options')
    config_group = parser.add_argument_group('Configuration')
    blast_group = parser.add_argument_group('Blast(default)')
    HMMsearch_group = parser.add_argument_group('HMMsearch(default)')
    clustal_group = parser.add_argument_group('clustal')
    pfam_group = parser.add_argument_group('pfam')

    # 输入选项
    input_group.add_argument('--fasta_input', type=str, required=True, help='Input FASTA file path', metavar='FASTA_FILE')
    input_group.add_argument('--target_sequence_file', type=str,  help="File containing the target sequence (150nt)", metavar='TARGET_SEQ_FILE')

    # 输出选项
    #output_group.add_argument('--blastoutput', '-o', type=str, help='Path to the BLAST output file', metavar='BLAST_OUTPUT')
    output_group.add_argument('--normailzationout', type=str, required=True, help='Normalization output file name', metavar='NORM_OUTPUT')
    output_group.add_argument('--extractcds', type=str, default='alignment_extract_cds.fasta', help='File where the complete CDS sequences are saved (default: alignment_extract_cds.fasta)', metavar='CDS_FILE')
    output_group.add_argument('--extractprotein', type=str, default='alignment_extract_protein.fasta', help='File where the complete protein sequences are saved (default: alignment_extract_protein.fasta)', metavar='PROTEIN_FILE')
    output_group.add_argument("--clustalo_extract_out", type=str ,help="clustalo extract output file.")

    # 配置选项
    config_group.add_argument('--num_threads', '-t', type=int, default=8, help='Number of threads to use (default: 8)', metavar='THREADS')
    config_group.add_argument('--cpu', type=str, default=8, help='Number of CPUs to use for hmmsearch.')
    config_group.add_argument('--log_file', '-l', type=str, help='Optional log file path', metavar='LOG_FILE')
    config_group.add_argument('--config', type=str, default='config.yaml', help='Configuration file path (default: config.yaml)', metavar='CONFIG_FILE')
    config_group.add_argument('--version', '-v', action='version',version=f'{__version__}',help="Print version information and exit")

    # pfam参数设置
    pfam_group.add_argument('--pfamkey', nargs='+',type=str, help='One or more pfam domain keys , Example: "Mur_ligase_C Another_HMM" ')
    pfam_group.add_argument('--pfamevalue', type=float, help='E-value threshold')

    # clustalo参数设置
    clustal_group.add_argument("--clustalo_out", type=str ,help="clustalo output file.")

    # blast默认参数设置
    blast_group.add_argument('--blast_type', choices=['blastp', 'blastx','diamond blastp', 'diamond blastx'], default='blastp', help='Type of BLAST program to use   (default:blastp)')
    blast_group.add_argument('--use_diamond', action='store_true', help='Use DIAMOND instead of BLAST for sequence alignment')
    blast_group.add_argument('--evalue', '-e', default='1e-5', help='E-value threshold  (default:1e-5)')
    blast_group.add_argument('--coverage_threshold', '-cov',type=float, help='Coverage threshold (default 0.7)')
    blast_group.add_argument('--identity_threshold', '-id',type=float, help='Identity threshold (default 30)')
    blast_group.add_argument('--max_target_seqs', default='5000', help='Maximum number of target sequences  (default:5000)')
    blast_group.add_argument('--matrix', default='BLOSUM62', help='Scoring matrix name  (default:BLOSUM62)')

    # hmmsearch默认参数设置
    HMMsearch_group.add_argument('--hmm_evalue',type=str, help='hmmsearch evalue  (default:1e-4)', metavar='evalue')
    HMMsearch_group.add_argument('--hmm_score',type=int, help='hmmsearch score  (default:0)', metavar='score')

    args = parser.parse_args()

    ### 加载配置文件、log
    logger = setup_logging(args.log_file)
    config = load_config(args.config)  

    ### 特殊参数处理
    #blast 使用配置文件中的默认值，如果命令行中提供了值，则覆盖
    coverage_threshold = args.coverage_threshold if args.coverage_threshold is not None else config['coverage_threshold']
    identity_threshold = args.identity_threshold if args.identity_threshold is not None else config['identity_threshold']
    
    #HMM 使用配置文件中的默认值，如果命令行中提供了值，则覆盖
    hmm_evalue = args.hmm_evalue if args.hmm_evalue is not None else config['hmm_evalue']
    hmm_score = args.hmm_score if args.hmm_score is not None else config['hmm_score']

    ### 中间临时文件
    tmp_blastout = "tmp_blastout.csv"
    normailzation_out = 'tmp-normailzation-out.csv'
    temp_filtered_file = "temp_filtered.fasta"
    temp_filtered_topfam_file = "temp_filtered_topfam.fasta"
    temp_extractcds = 'tmp_extractcds.fasta'
    temp_extractprotein = 'tmp_extractprotein.fasta'
    tmp_inputfasta_clustao = "tmp_inputfasta_clustao.fasta"
    tmp_inputfasta_clustao_hmmbuild = "tmp_inputfasta_clustao.hmm"
    tmp_hmmsearch_output = 'hmmsearch_output.txt'
    tmp_hmmsearch_excat = 'tmp_hmmsearch_excat.csv'
    tmp_with_hmm_normailzationout = 'tmp_with_hmm_normailzationout.csv'
    tmp_pfam_output = "tmp_pfam_output.csv"
    merged_file = "tmp_merged.fasta"
    clustal_output_file = "tmp_clu-out.fasta"
    #tmp_pfam_result = "tmp_pfam_result.csv"
    tmp_pfam_result = "tmp_pfam_output.csv"
    tmp_pfam_excat = 'tmp_pfam_excat.csv'
    
    #初始化临时文件
    remove_file(merged_file)  # 删除之前的 merged.fasta 文件
    remove_file(tmp_blastout)  
    remove_file(args.extractcds)  
    remove_file(args.extractprotein)  
    remove_file(temp_filtered_file)  
    remove_file(temp_filtered_topfam_file)  
    remove_file(args.normailzationout) 

    # 判断中间文件是否保留
    if args.clustalo_out :
        clustal_output_file = args.clustalo_out
    if config["clustal"] == "yes" : 
        # 读取目标序列文件，获取序列ID
        with open(args.target_sequence_file, 'r') as f:
            target_sequence_id = f.readline().strip()

    ### 将命令行参数整理成字典
    #blast
    blast_params = {
        'blast_path': config[args.blast_type + '_path'],  # 根据选择的BLAST类型确定路径
        'diamond_path':config['diamond_' + args.blast_type + '_path'],
        'db_path': config['db_path'],
        'diamond_db_path': config['diamond_db_path'],
        'evalue': args.evalue,
        'max_target_seqs': args.max_target_seqs,
        'matrix': args.matrix,
        'use_diamond': args.use_diamond,
        'num_threads': args.num_threads,
        'blast_type': args.blast_type,
    }
    #HMM
    hmmsearch_params = {
        'hmmsearch' : config['hmmsearch_path'] ,
        'hmmscan' : config['hmmscan_path'] ,
        'hmmbuild': config['hmmbuild'] ,
        'clustalo' : config['clustalo_path'] ,
        'hmm_db' : config['db_path']
    }
    #pfam
    pfam_params = {
        'pfam_path': config['pfam_scan'],  
        'pfamdb_path': config['pfam_dirextory'],
        'cpu': args.cpu,
        'out': tmp_pfam_output ,
    }

    ### 整体流程
    #rename(args.fasta_input, args.fasta_output, logger)
    ##BLAST序列比对
    #run_blast(args.fasta_rename, tmp_blastout, args.num_threads, logger, config)
    run_blast(args.fasta_input, tmp_blastout, blast_params, logger)
    ##HMM进一步筛选比对结果
    if config["hmmsearch"] == "yes" : 
        #添加一步hmmsearch进行进一步的筛选
        run_HMMprcess(args.fasta_input, args.cpu, tmp_hmmsearch_output,hmmsearch_params, tmp_inputfasta_clustao, tmp_inputfasta_clustao_hmmbuild, logger)
        #调用函数执行比较和保存结果
        excat_hmm_postition(tmp_hmmsearch_output, tmp_hmmsearch_excat, hmm_evalue, hmm_score)
        compare_files_and_keep_df1_only(tmp_blastout, tmp_hmmsearch_excat, tmp_with_hmm_normailzationout, logger)
        #覆盖原始的normailzationout
        #shutil.move(tmp_with_hmm_normailzationout, tmp_blastout)
    ###pfam进行进一步筛选
    if config["pfam"] == "yes" : 
        #将这个蛋白结果提取出来，用于pfam的搜索
        findextractcdsprotein(tmp_blastout, args.num_threads, args.extractcds ,args.extractprotein,identity_threshold,coverage_threshold, logger, config)
        #进行筛选，将序列去重
        filter_duplicate_proteins(args.extractprotein, temp_filtered_topfam_file)  # 过滤重复的蛋白质ID并写入临时文件
        #使用pfam进行进一步筛选
        run_pfam(temp_filtered_topfam_file,pfam_params, logger)   #调用pfam.py，结果保存在tmp_pfam_result
        pfam_compare_files_and_keep_df1_only(tmp_pfam_result, tmp_blastout, args.pfamkey, tmp_pfam_output,logger)
        #shutil.move(tmp_pfam_excat, args.blastoutput)
        #shutil.move(tmp_pfam_output, tmp_blastout)
    filenormalization(tmp_blastout, args.normailzationout, coverage_threshold, identity_threshold, logger)
    # makeformat(normailzation_out, args.normailzationout, logger,config)
    findextractcdsprotein(tmp_blastout, args.num_threads, args.extractcds ,args.extractprotein,identity_threshold,coverage_threshold, logger, config)
    ##clustalo多序列比对
    if config["clustal"] == "yes" : 
        filter_duplicate_proteins(args.extractprotein, temp_filtered_file)  # 过滤重复的蛋白质ID并写入临时文件
        merge_fasta_files(temp_filtered_file, args.target_sequence_file, merged_file)   #融合、制作多序列比对文件
        run_clustalo(merged_file, clustal_output_file, args.num_threads, logger, config)
        remove_file(merged_file)  # 删除临时合并的文件
        remove_file(temp_filtered_file)  # 删除临时过滤的文件
    ##提取达标序列
        alignment = read_alignment_from_file(clustal_output_file)
        start_position, end_position = find_start_end_positions(alignment, target_sequence_id)
        end_position = find_end_position(alignment, target_sequence_id)
        trimmed_alignment = trim_alignment(alignment, start_position, end_position, target_sequence_id)
        write_alignment_to_file(trimmed_alignment, args.clustalo_extract_out)

    ### 临时文件删除
    # 判断是否有需要保存的文件
        if args.clustalo_out is None:
            remove_file(clustal_output_file)



if __name__ == "__main__":
    main()

