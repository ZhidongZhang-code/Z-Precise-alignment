import sys
import re
import argparse
import logging

##输入文件格式：
#0                          1                                             2        3        4       5       6      7            8           9        10     11      
#qseqid                     sseqid                                       qlen    slen    qstart  qend    sstart  send       evalue      bitscore length pident
#K00097_**_eco:b0052     GCF_000190995.1_ASM19099v1_**_WP_000241242.1    329     329     1       329     1       329     1.9e-182        644.8   329     325

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

def filenormalization(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        #####截取合适的cov与identity的值
        new_lines = []
        for line in infile:
            line = line.strip()
            if re.match(r'^\w+', line):
                test = line.split('\t')
                if (abs(float(test[5]) - float(test[4]) + 1) / int(test[2]) >= 0.4 and      #（qstar-qend）/qlen
                        abs(float(test[7]) - float(test[6]) + 1) / int(test[3]) >= 0.4 and      #（sstar-send）sqlen
                        float(test[11]) >= 20):      #identity
                    new_lines.append(line + '\n')
        ######对于_**_内容，取重要的部分
        #K00097_**_eco:b0052     GCF_000190995.1_ASM19099v1_**_WP_000241242.1
    #GCF_000190995.1         K00097
        for line in new_lines:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                col1, col2 = parts[0], parts[1]
                col1 = col1.split("_**_")[0]
                col2 = "_".join(col2.split("_", 2)[:2])
                new_line = "{}\t{}\t{}\n".format(col2, col1, '\t'.join(parts[2:]))
                outfile.write(new_line)
        logging.info("File normalization processed successfully.")

def main():
    parser = argparse.ArgumentParser(description='File normalization script')
    parser.add_argument('--input', type=str, required=True, help='Input file name')
    parser.add_argument('--output', type=str, required=True, help='Output file name')
    parser.add_argument('--log', type=str, help='Log file path (optional)')
    args = parser.parse_args()

    input_file = args.input
    output_file = args.output

    filenormalization(input_file, output_file)
    logging.info(f"Processing complete. Results have been written to {output_file}")

if __name__ == "__main__":
    main()