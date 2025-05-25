import subprocess
import argparse
import argparse
import yaml
from log_config import setup_logging
import logging
from FileOperation import remove_file

# 加载YAML配置文件的函数
def load_config(config_file):
    with open(config_file, 'r') as file:
        return yaml.safe_load(file)

def run_command(command):
    try:
        subprocess.run(command, check=True)
        print(f"Command executed successfully: {' '.join(command)}")
    except subprocess.CalledProcessError as e:
        print(f"Error executing command: {' '.join(command)}")
        print(e)

def run_HMMprcess(clustalo_input, cpu, hmmsearch_output,hmmsearch_params,tmp_inputfasta_clustao,tmp_inputfasta_clustao_hmmbuild):
    # print(f"clustalo_input: {clustalo_input}")
    # print(f"cpu: {cpu}")
    # print(f"hmmsearch_output: {hmmsearch_output}")
    # print(f"hmmsearch_params: {hmmsearch_params}")
    # print(f"tmp_inputfasta_clustao: {tmp_inputfasta_clustao}")
    # print(f"tmp_inputfasta_clustao_hmmbuild: {tmp_inputfasta_clustao_hmmbuild}")

    # Step 1: Run clustalo
    remove_file(tmp_inputfasta_clustao)
    clustalo_command = [ hmmsearch_params['clustalo'], "-i", clustalo_input, "-o", tmp_inputfasta_clustao]
    run_command(clustalo_command)

    # Step 2: Run hmmbuild
    hmmbuild_command = [hmmsearch_params['hmmbuild'], tmp_inputfasta_clustao_hmmbuild, tmp_inputfasta_clustao]
    run_command(hmmbuild_command)

    # Step 3: Run hmmsearch
    hmmsearch_command = [hmmsearch_params['hmmsearch'], "--cpu", cpu, "-o", hmmsearch_output, tmp_inputfasta_clustao_hmmbuild, hmmsearch_params['hmm_db']]
    run_command(hmmsearch_command)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Run clustalo, hmmbuild, and hmmsearch in sequence.')
    parser.add_argument('--fasta_input', type=str, help='Input FASTA file path', metavar='FASTA_FILE')
    parser.add_argument('--cpu', type=str, help='Number of CPUs to use for hmmsearch.', required=True)
    parser.add_argument('--hmmsearch_output', type=str, help='Output file name for hmmsearch results.', required=True)
    parser.add_argument('--config', type=str, default='config.yaml', help='Configuration file path (default: config.yaml)', metavar='CONFIG_FILE')

    args = parser.parse_args()

    # 临时文件
    tmp_inputfasta_clustao = "tmp_inputfasta_clustao.fasta"
    tmp_inputfasta_clustao_hmmbuild = "tmp_inputfasta_clustao.hmm"

    # 加载配置
    config = load_config(args.config)

    # 使用 setup_logging 配置日志
    logger = setup_logging(log_file_path=args.log_file)

    hmmsearch_params = {
        'hmmsearch' : config['hmmsearch_path'] ,
        'hmmscan' : config['hmmscan_path'] ,
        'hmmbuild': config['hmmbuild'] ,
        'clustalo' : config['clustalo_path'] ,
        'hmm_db' : config['db_path']
    }
    #print(hmmsearch_params)
    #独立脚本模块文件
    tmp_hmmsearch_output = args.hmmsearch_output

    run_HMMprcess(args.input, args.cpu, tmp_hmmsearch_output,hmmsearch_params,tmp_inputfasta_clustao,tmp_inputfasta_clustao_hmmbuild)

    # 临时文件删除
    #remove_file(tmp_inputfasta_clustao)
    #remove_file(tmp_inputfasta_clustao_hmmbuild)