####
#前置步骤，标准化蛋白库的id，
#GCF_923868895.1_FutureBRF_Halomonas_Bluephagenesis_TD01_Genome_v1_**_WP_009722895.1
#GCF_923868895.1_**_WP_009722895.1
####



import re

def normalize_ids(input_file, output_file):
    # 定义正则表达式模式来匹配所需的ID部分
    pattern = re.compile(r"(GCF_\d+\.\d+)_.*_([A-Z]{2}_\d+\.\d+)")

    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            if line.startswith('>'):
                # 提取并连接匹配的ID部分
                match = pattern.search(line)
                if match:
                    normalized_id = f"{match.group(1)}_**_{match.group(2)}"
                    outfile.write(f">{normalized_id}\n")
                else:
                    # 如果没有匹配，保持原样
                    outfile.write(line)
            else:
                # 保持非ID行原样
                outfile.write(line)

# 使用文件路径调用函数
input_file = '../all-hum-medusa-db-genmoic-protein'
output_file = 'all-hum-medusa-db-genmoic-protein-normalizename'
normalize_ids(input_file, output_file)
