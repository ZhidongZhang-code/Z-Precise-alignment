import os
import sys
import argparse
import logging

def setup_logging(log_file_path):
    # 创建一个logger
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)

    if log_file_path:
        # 如果指定了日志文件路径，只创建用于输出到文件的handler
        file_handler = logging.FileHandler(log_file_path)
        file_handler.setLevel(logging.INFO)
        file_formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        file_handler.setFormatter(file_formatter)
        logger.addHandler(file_handler)
    else:
        # 否则，只创建用于输出到控制台的handler
        console_handler = logging.StreamHandler()
        console_handler.setLevel(logging.INFO)
        console_formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        console_handler.setFormatter(console_formatter)
        logger.addHandler(console_handler)

def rename(input_files, output_file_path):
    for input_file_path in input_files:
        if not input_file_path.endswith(".fasta"):
            logging.warning(f"Ignoring non-fasta file: {input_file_path}")
            continue

        if not os.path.isfile(input_file_path):
            logging.error(f"File not found: {input_file_path}")
            continue

        file_name = os.path.splitext(os.path.basename(input_file_path))[0]
        modified_lines = []

        with open(input_file_path, 'r') as input_file:
            for line in input_file:
                line = line.strip()
                if line.startswith('>') and '_**_' not in line:
                    modified_lines.append(f'>{file_name}_**_{line[1:]}')
                else:
                    modified_lines.append(line)

        modified_content = '\n'.join(modified_lines)

        with open(output_file_path, 'a') as output_file:
            output_file.write(modified_content)

        logging.info(f"Processed: {input_file_path}")

    logging.info(f"All files have been processed. Output saved to: {output_file_path}")

def main():
    parser = argparse.ArgumentParser(description='Process input FASTA files and modify headers.')
    parser.add_argument('--output', type=str, required=True, help='Output FASTA file path')
    parser.add_argument('--input', nargs='+', type=str, required=True, help='Input FASTA file paths')
    parser.add_argument('--log', type=str, help='Log file path (optional)')
    args = parser.parse_args()

    setup_logging(args.log)

    output_file_path = args.output
    input_files = args.input

    rename(input_files, output_file_path)

if __name__ == "__main__":
    main()
