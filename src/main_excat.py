import argparse
from position_sequence_matcher import parse_position_arg, find_sequences_with_rf
from sequence_pattern_extractor import read_names_from_file, find_matching_sequences, save_matches_to_file
from excat_site_protein import process_site_sequences
from excat_site_cds import process_cds_fasta
from FileOperation import remove_file

###日志
###修改cds输出格式，并且要注意会有没有找到的 ：No match found for protein ID WP_040344983.1, genome ID GCF_000165795.1_ASM16579v1, position 364
###还没想好怎么处理，在excat_site_cds的模块中
###规范化输出文件，该删除的删除   （完成）
###

###python ./src/main_excat.py --input clu-out.fasta --position 25:R 28:R --extractprotein alignment_extract_protein.fasta --extractcds alignment_extract_cds.fasta --site_protein_output site_protein.output --site_cds_output site_cds.output

def main():
    parser = argparse.ArgumentParser(description="Process protein sequences")
    
    # 创建参数分组
    input_group = parser.add_argument_group('Input options')
    output_group = parser.add_argument_group('Output options')

    #输入参数
    input_group.add_argument("--input",'-i', required=True, help="Input ClustalO FASTA file.", metavar="FASTA")
    input_group.add_argument("--position",required=True, nargs='+', type=str, help="Positions and corresponding amino acids.", metavar="POS:AA")
    input_group.add_argument('--extractprotein', required=True, type=str, help='Output file for complete protein sequences.', metavar='PROTEIN_FILE')
    input_group.add_argument('--extractcds', required=True, type=str, help='Output file for complete CDS sequences.', metavar='CDS_FILE') 
    #输出参数
    output_group.add_argument("--site_protein_output", required=True,  type=str, help="Output file for matched protein sequences.", metavar="PROTEIN_OUT")
    output_group.add_argument("--site_cds_output", required=True,  type=str, help="Output file for matched CDS sequences.", metavar="CDS_OUT")

    args = parser.parse_args()

    #临时文件
    tmp_position_swquence = "tmp_position_swquence.fasta"
    tmp_matchesq =  "tmp_matches.fasta"

    #日志、cfg文件模块加载
    #logger = setup_logging(args.log_file)

    # 定义文件路径
    fasta_file = args.input
    #input_names_file = args.input_file
    protein_db_file = args.extractprotein
    #output_file_path = args.output_file 

    # 输入序列位置信息  根据输入位点（25：R） 有这个位置的序列的下10个氨基酸
    # >WP_011802784.1 RYSRTRDCFE
    positions = parse_position_arg(args.position)
    # 获取所有满足条件的序列的ID以及相应的锚点序列信息
    rf_sequences = find_sequences_with_rf(fasta_file, positions)
    #将得到的序列位置信息保存到文件 >WP_011802784.1 RYSRTRDCFE
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
    #移除中间文件
    #remove_file(tmp_matchesq)
    #remove_file(tmp_position_swquence)



if __name__ == "__main__":

    main()