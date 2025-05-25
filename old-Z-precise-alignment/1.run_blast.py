import subprocess
import argparse
import logging
import sys

def setup_logging(log_file_path):
    # 创建一个logger
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)

    # 检查是否已经存在文件或控制台处理器，避免重复添加
    has_file_handler, has_console_handler = False, False
    for handler in logger.handlers:
        if isinstance(handler, logging.FileHandler):
            has_file_handler = True
        elif isinstance(handler, logging.StreamHandler):
            has_console_handler = True

    if log_file_path and not has_file_handler:
        # 如果指定了日志文件路径，只创建用于输出到文件的handler
        file_handler = logging.FileHandler(log_file_path)
        file_handler.setLevel(logging.INFO)
        file_formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        file_handler.setFormatter(file_formatter)
        logger.addHandler(file_handler)
    elif not log_file_path and not has_console_handler:
        # 否则，只创建用于输出到控制台的handler
        console_handler = logging.StreamHandler()
        console_handler.setLevel(logging.INFO)
        console_formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        console_handler.setFormatter(console_formatter)
        logger.addHandler(console_handler)

def run_blastp(query_file, output_file, num_threads):
    command = [
        '/public/home/TonyWuLab/zhangzhd/anaconda3/envs/blast/bin/blastp',
        '-db', '/public/group_share_data/TonyWuLab/zzd/db/medusa-2023/all-hum-medusa-db-genmoic-protein',
        '-query', query_file,
        '-out', output_file,
        '-outfmt', '6 qseqid sseqid qlen slen qstart qend sstart send evalue bitscore length pident',
        '-evalue', '1e-5',
        '-max_target_seqs', '5000',
        '-num_threads', str(num_threads)
    ]

    try:
        subprocess.run(command, check=True)
        logging.info("Blastp command executed successfully.")
    except subprocess.CalledProcessError as e:
        logging.error(f"Error executing Blastp command: {e}")
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(description='Run blastp with specified parameters')
    parser.add_argument('--query', '-q', type=str, required=True, help='Path to the query file')
    parser.add_argument('--output', '-o', type=str, required=True, help='Path to the output file')
    parser.add_argument('--num_threads', '-t', type=int, default=8, help='Number of threads (default: 8)')
    parser.add_argument('--log_file', '-l', type=str, default='blastp.log', help='Log file path (default: blastp.log)')
    args = parser.parse_args()

    # 配置日志
    setup_logging(args.log_file)

    run_blastp(args.query, args.output, args.num_threads)

if __name__ == "__main__":
    main()
