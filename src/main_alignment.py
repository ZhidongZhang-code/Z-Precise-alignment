import argparse
import csv
import shutil
from blast_module import setup_logging, run_blast, load_config  
from FileOperation import remove_file
from fasta_renamer import rename 
from alignment_normalization import filenormalization
from find_extract_cds_protein import findextractcdsprotein
from makeformat import makeformat
from run_clustalo_process import *
from hmm_module import run_HMMprcess
from process_hmmresult import excat_hmm_postition, compare_files_and_keep_df1_only
from pfam_module import *
from process_pfamresult import pfam_compare_files_and_keep_df1_only

###日志：
###blastout算是一个默认参数，最后需要删除
###验证结果是否正确：可以与urda的150nt比较
###还是没有解决rename.py的问题，这个脚本是不是应该单独拿出来？
###除了blastp还应该可以支持blastx 由于数据库是蛋白库，所以不能支持blastn    完成
###blastp的参数应该作为一个参数能在脚本中进行选择（cov、id、e值、k）    完成
###如果某一部分出了问题，需要报错在那一部分出的错误 完成
###hmm模块还没有装载好，并且还是考虑这个是否必须使用    完成 
###在config中选择是否使用HMM，默认为使用    完成 
###pfam模块也要加载
###pfam是不是要提供全部的type用来参考？
###每一个模块运行后，要给出运行这个模块的命令行脚本
###参数列表检查是不是存在必须要有的数据库（参考pfam脚本）
###还是得做GUI，在乌班图这种有图形化界面的linux系统上使用
###测试完成后临时文件的删除

__version__ = '1.0'

def main():
    parser = argparse.ArgumentParser(description='Control script to run BLAST analysis')
    parser = argparse.ArgumentParser(description=f'pfam_scan.py v.{__version__}')

    # 创建参数分组
    input_group = parser.add_argument_group('Input options')
    output_group = parser.add_argument_group('Output options')
    config_group = parser.add_argument_group('Configuration')
    blast_group = parser.add_argument_group('Blast(default)')
    HMMsearch_group = parser.add_argument_group('HMMsearch(default)')
    pfam_group = parser.add_argument_group('pfam')

    # 输入选项
    input_group.add_argument('--fasta_input', type=str, help='Input FASTA file path', metavar='FASTA_FILE')
    input_group.add_argument('--target_sequence_file', type=str, required=True, help="File containing the target sequence (150nt)", metavar='TARGET_SEQ_FILE')

    # 输出选项
    output_group.add_argument('--blastoutput', '-o', type=str, help='Path to the BLAST output file', metavar='BLAST_OUTPUT')
    output_group.add_argument('--normailzationout', type=str, required=True, help='Normalization output file name', metavar='NORM_OUTPUT')
    output_group.add_argument('--extractcds', type=str, default='alignment_extract_cds.fasta', help='File where the complete CDS sequences are saved (default: alignment_extract_cds.fasta)', metavar='CDS_FILE')
    output_group.add_argument('--extractprotein', type=str, default='alignment_extract_protein.fasta', help='File where the complete protein sequences are saved (default: alignment_extract_protein.fasta)', metavar='PROTEIN_FILE')

    # 配置选项
    config_group.add_argument('--num_threads', '-t', type=int, default=8, help='Number of threads to use (default: 8)', metavar='THREADS')
    config_group.add_argument('--cpu', type=str, default=8, help='Number of CPUs to use for hmmsearch.')
    config_group.add_argument('--log_file', '-l', type=str, help='Optional log file path', metavar='LOG_FILE')
    config_group.add_argument('--config', type=str, default='config.yaml', help='Configuration file path (default: config.yaml)', metavar='CONFIG_FILE')
    config_group.add_argument('--version', '-v', action='version',version=f'{__version__}',help="Print version information and exit")

    # blast默认参数设置
    blast_group.add_argument('--blast_type', choices=['blastp', 'blastx'], default='blastp', help='Type of BLAST program to use   (default:blastp)')
    blast_group.add_argument('--evalue', '-e', default='1e-5', help='E-value threshold  (default:1e-5)')
    blast_group.add_argument('--coverage_threshold', '-cov',type=float, help='Coverage threshold (default 0.7)')
    blast_group.add_argument('--identity_threshold', '-id',type=float, help='Identity threshold (default 30)')
    blast_group.add_argument('--max_target_seqs', default='5000', help='Maximum number of target sequences  (default:5000)')
    blast_group.add_argument('--matrix', default='BLOSUM62', help='Scoring matrix name  (default:BLOSUM62)')

    # hmmsearch默认参数设置
    HMMsearch_group.add_argument('--hmm_evalue', default='1e-4',type=str, help='hmmsearch evalue  (default:1e-4)', metavar='evalue')
    HMMsearch_group.add_argument('--hmm_score', default='200',type=int, help='hmmsearch score  (default:200)', metavar='score')

    # pfam参数设置
    pfam_group.add_argument('--pfamkey', nargs='+', type=str, help='One or more pfam domain keys , Example: "Mur_ligase_C Another_HMM" ', required=True)
    pfam_group.add_argument('--pfamevalue', type=float, help='E-value threshold')

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
    tmp_inputfasta_clustao = "tmp_inputfasta_clustao.fasta"
    tmp_inputfasta_clustao_hmmbuild = "tmp_inputfasta_clustao.hmm"
    tmp_hmmsearch_output = 'hmmsearch_output.txt'
    tmp_hmmsearch_excat = 'tmp_hmmsearch_excat.csv'
    tmp_with_hmm_normailzationout = 'tmp_with_hmm_normailzationout.csv'
    tmp_pfam_output = "tmp_pfam_output.csv"
    merged_file = "merged.fasta"
    clustal_output_file = "clu-out.fasta"
    tmp_pfam_result = "tmp_pfam_result.csv"
    tmp_pfam_excat = 'tmp_pfam_excat.csv'

    # 读取目标序列文件，获取序列ID
    with open(args.target_sequence_file, 'r') as f:
        target_sequence_id = f.readline().strip()

    ### 将命令行参数整理成字典
    #blast
    blast_params = {
        'blast_path': config[args.blast_type + '_path'],  # 根据选择的BLAST类型确定路径
        'db_path': config['db_path'],
        'evalue': args.evalue,
        'max_target_seqs': args.max_target_seqs,
        'matrix': args.matrix,
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
    #run_blast(args.query, tmp_blastout, blast_params, logger)
    ##HMM进一步筛选比对结果
    # if config["hmmsearch"] == "yes" : 
    #     #添加一步hmmsearch进行进一步的筛选
    #     run_HMMprcess(args.fasta_input, args.cpu, tmp_hmmsearch_output,hmmsearch_params, tmp_inputfasta_clustao, tmp_inputfasta_clustao_hmmbuild)
    #     #调用函数执行比较和保存结果
    #     excat_hmm_postition(tmp_hmmsearch_output, tmp_hmmsearch_excat, hmm_evalue, hmm_score)
    #     compare_files_and_keep_df1_only(args.blastoutput, tmp_hmmsearch_excat, tmp_with_hmm_normailzationout, logger)
    #     #覆盖原始的normailzationout
    #     shutil.move(tmp_with_hmm_normailzationout, args.blastoutput)
    ###将这个蛋白结果提取出来，用于pfam的搜索
    #findextractcdsprotein(args.blastoutput, args.num_threads, args.extractcds ,args.extractprotein, logger, config)
    ###pfam进行进一步筛选
    if config["pfam"] == "yes" : 
        #使用pfam进行进一步筛选
        run_pfam(args.extractprotein,pfam_params, logger)   #调用pfam.py，结果保存在tmp_pfam_result
        pfam_compare_files_and_keep_df1_only(tmp_pfam_result, args.blastoutfile, args.pfamkey, tmp_pfam_excat,logger)
        shutil.move(tmp_pfam_excat, args.blastoutfile)
    filenormalization(args.blastoutfile, args.normailzationout, coverage_threshold, identity_threshold, logger)
    makeformat(normailzation_out, args.normailzationout, logger,config)
    findextractcdsprotein(tmp_blastout, args.num_threads, args.extractcds ,args.extractprotein, logger, config)
    ##clustalo多序列比对
    remove_file(merged_file)  # 删除之前的 merged.fasta 文件
    filter_duplicate_proteins(args.extractprotein, temp_filtered_file)  # 过滤重复的蛋白质ID并写入临时文件
    merge_fasta_files(temp_filtered_file, args.target_sequence_file, merged_file)   #融合、制作多序列比对文件
    run_clustalo(merged_file, clustal_output_file, args.num_threads, logger, config)
    #remove_file(merged_file)  # 删除临时合并的文件
    #remove_file(temp_filtered_file)  # 删除临时过滤的文件
    ##提取达标序列
    alignment = read_alignment_from_file(clustal_output_file)
    start_position, end_position = find_start_end_positions(alignment, target_sequence_id)
    end_position = find_end_position(alignment, target_sequence_id)
    trimmed_alignment = trim_alignment(alignment, start_position, end_position, target_sequence_id)
    write_alignment_to_file(trimmed_alignment, clustal_output_file)

    ### 临时文件删除


if __name__ == "__main__":
    main()
