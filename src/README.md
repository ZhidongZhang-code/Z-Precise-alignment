#需要安装#
conda install pyyaml -y
conda install pandas -y

项目名称：[项目名称]

项目描述：本项目融合了多种序列比对技术,与活性位点相结合,鉴定生境中同源不同功基因。

安装

本软件支持以下几种安装方式:
1. 下载包手动安装
从release页面下载最新版本,解压后运行:
```shell
python setup.py install
```

（#注意安装数据库还没有做好）

快速上手：
1：使用main_alignment.py进行序列比对，在config中筛选需要什么软件
python ./src/main_alignment.py  --fasta_input MurC.fasta  --normailzationout blast_normailzation.csv --extractcds blast_normailzation_excate_cds.fasta --extractprotein blast_normailzation_excate_protein.fasta --num_threads 8 --cpu 1 --log_file 2024_3_13.log --config ./src/config.yaml --target_sequence_file EcMurC-383-50bp.fasta --pfamkey Mur_ligase_C  --clustalo_extract_out clu-extract-out.fasta --clustalo_out clu-out.fasta
2：使用main_excat.py挑选感兴趣的活性位点
python src/main_excat.py --input clu-extract-out.fasta --position 25:R 28:R --extractprotein blast_normailzation_excate_protein.fasta --extractcds blast_normailzation_excate_cds.fasta --site_protein_output site-protein.fasta --site_cds_output site-cds.fasta --formatoutput 111 --config ./src/config.yaml --log 2024_3_12.log
3：或者直接使用比对的结果构建表格
python ./src/main_only_format.py --input alignment_extract_protein.fasta --formatoutput 1 --config ./src/config.yaml


使用

首先，运行脚本./src/main_alignment.py，对于目标序列、活性位点信息，在人类肠道微生物数据库中进行比对，提取蛋白序列、CDS序列、标准化的文件表、clustalo比对的文件等。

```shell
#python ./src/main_alignment.py  --fasta_input MurC.fasta  --normailzationout blast_normailzation.csv --extractcds blast_normailzation_excate_cds.fasta --extractprotein blast_normailzation_excate_protein.fasta --num_threads 8 --cpu 1 --log_file 2024_3_13.log --config ./src/config.yaml --target_sequence_file EcMurC-383-50bp.fasta --pfamkey Mur_ligase_C  --clustalo_extract_out clu-extract-out.fasta --clustalo_out clu-out.fasta
```

```shell
python ./src/main_alignment.py --help
```

#输入--help后会显示帮助信息

```shell
#输入选项：
--fasta_input  目标查询的fasta文件
--target_sequence_file  目标活性位点序列文件，一般是目标活性位点上下游30bp
#输出选项：
--normailzationout 输出一个标准化的文件表，包含序列的物种
--extractcds 输出一个提取的cds文件
--extractprotein 输出一个提取的蛋白质文件
--clustalo_extract_out 输出一个clustalo比对的文件
#控制选项
--num_threads 线程数
--cpu 使用的cpu数
--log_file 输出日志文件
--config 配置文件
--version 输出版本信息
#blast选项
--blast_type 选择blast的类型，默认为blastp
--evalue 选择blast的evalue，默认为1e-5
--coverage_threshold 选择blast的coverage，默认为0.7 可以在cfg文件中修改
--identity_threshold 选择blast的identity，默认为30 可以在cfg文件中修改
--max_target_seqs 选择blast的max_target_seqs，默认为5000
--matrix 选择blast的matrix，默认为BLOSUM62
#HMMsearch选项
--hmm_evalue 选择hmmsearch的evalue，默认为1e-4
--hmm_score 选择hmmsearch的score，默认为200 可以在cfg文件中修改
#clustalo选项
--clustalo_out 输出clustalo的完整比对结果文件
#pfam选项
--pfamkey 一个或多个pfam域键，例如:"Mur_ligase_C Another_HMM"
--pfamevalue 选择pfam的evalue，默认为1e-5
```

在运行过程中，会生成一个日志文件，记录运行的时间、参数、输出文件等信息。

对于结果得到的多序列比对结果，使用./src/main_extract.py脚本，选择合适的活性位点，提取出相应的序列

```shell
#python src/main_excat.py --input clu-extract-out.fasta --position 25:R 28:R --extractprotein blast_normailzation_excate_protein.fasta --extractcds blast_normailzation_excate_cds.fasta --site_protein_output site-protein.fasta --site_cds_output site-cds.fasta --formatoutput 111 --config ./src/config.yaml --log 2024_3_12.log
```

```shell
python src/main_excat.py --help
```

#输入--help后会显示帮助信息

```shell
#输入选项：
--input 输入的clustalo比对文件
--position 选择的活性位点，支持多个位点的选择，例如：25:R 28:R，表示选择第25位和第28位的氨基酸。使用:代表该位点为对应的氨基酸，！代表该位点为非对应的氨基酸
--extractprotein 输入的蛋白质文件，一般是main_alignment.py输出的extractprotein
--extractcds 输入的cds文件，一般是main_alignment.py输出的extractcds
#输出选项：
--site_protein_output 输出的蛋白质文件
--site_cds_output 输出的cds文件
--formatoutput 输出一个文件，里面有序列的物种信息，是一个csv文件
#控制选项
--config 配置文件
--log_file 输出日志文件

