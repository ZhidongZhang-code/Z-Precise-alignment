import re
import sys

def process_file(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        # 初始化字典来存储数据
        data = {}
        
        for line in infile:
            line = line.strip()
            # 匹配细菌名称行
            match_bacteria = re.match(r'^>> (\S+)\s+.*', line)
            if match_bacteria:
                bacteria = match_bacteria.group(1)
                # 初始化或重置细菌数据
                if bacteria not in data:
                    data[bacteria] = {'domains': 0, 'score': 0, 'evalue': 0,
                                      'hmmfrom': 50000, 'hmmto': 0,
                                      'alignfrom': 50000, 'alignto': 0}
            
            # 匹配数据行
            match_data = re.match(r'(\d) \!\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\d+)\s+.*\s+(\d+)\s+(\d+)\s+.*\s+(\d+)\s+(\d+).*\S+', line)
            if match_data:
                domain, score, evalue, hmmfrom, hmmto, alignfrom, alignto = match_data.groups()[0], float(match_data.groups()[1]), float(match_data.groups()[3]), int(match_data.groups()[5]), int(match_data.groups()[6]), int(match_data.groups()[7]), int(match_data.groups()[8])
                # 更新数据
                data[bacteria]['domains'] += 1
                data[bacteria]['score'] += score
                data[bacteria]['evalue'] += evalue
                data[bacteria]['hmmfrom'] = min(data[bacteria]['hmmfrom'], hmmfrom)
                data[bacteria]['hmmto'] = max(data[bacteria]['hmmto'], hmmto)
                data[bacteria]['alignfrom'] = min(data[bacteria]['alignfrom'], alignfrom)
                data[bacteria]['alignto'] = max(data[bacteria]['alignto'], alignto)

        # 写入表头
        outfile.write("Bac_protein\tdomains\tscore\tmean_evalue\thmm_from\thmm_to\talign_from\talign_to\thmm_cov\n")
        
        # 写入数据
        for bacteria, info in data.items():
            if info['evalue'] <= 1e-4 and info['score'] >= 200:
                mean_evalue = info['evalue'] / info['domains']
                hmm_cov = (abs(info['hmmto'] - info['hmmfrom']) + 1) / 63
                outfile.write(f"{bacteria}\t{info['domains']}\t{info['score']}\t{mean_evalue}\t{info['hmmfrom']}\t{info['hmmto']}\t{info['alignfrom']}\t{info['alignto']}\t{hmm_cov}\n")

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: python script.py input_file output_file")
        sys.exit(1)
        
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    process_file(input_file, output_file)

