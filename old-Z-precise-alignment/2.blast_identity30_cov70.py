import sys
import re

##输入文件格式：
#0                          1                                             2        3        4       5       6      7            8           9        10     11      
#qseqid                     sseqid                                       qlen    slen    qstart  qend    sstart  send       evalue      bitscore length positive
#K00097_**_eco:b0052     GCF_000190995.1_ASM19099v1_**_WP_000241242.1    329     329     1       329     1       329     1.9e-182        644.8   329     325

input_file = sys.argv[1]
output_file = sys.argv[2]

def process_file(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            line = line.strip()
            if re.match(r'^\w+', line):
                test = line.split('\t')
                if (abs(float(test[5]) - float(test[4]) + 1) / int(test[2]) >= 0.4 and   #（qstar-qend）/qlen
                        abs(float(test[7]) - float(test[6]) + 1) / int(test[3]) >= 0.4 and  #（sstar-send）sqlen
                        float(test[9]) >= 20): #identity
                    outfile.write(line + '\n')

# 调用函数，传入输入文件名和输出文件名
input_file = sys.argv[1]
output_file = sys.argv[2]
process_file(input_file, output_file)