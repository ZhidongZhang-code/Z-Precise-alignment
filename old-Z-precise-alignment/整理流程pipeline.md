#0：命名规范化，需要将同一个KO的序列放在同一个txt中
python script.py --output output_file.fasta --input input_file1.fasta input_file2.fasta ...
#1.1:运行diamond
diamond blastp --db /public/group_share_data/TonyWuLab/zzd/db/medusa-2023/all-hum-medusa-db-genmoic-protein-diamond.dmnd \
--query output.fasta  --out diamondout.csv \
--outfmt 6 qseqid sseqid qlen slen qstart qend sstart send evalue bitscore length pident \
--max-target-seqs 50000 --threads 8
#1.1：运行blast
blastp -db /public/group_share_data/TonyWuLab/zzd/db/medusa-2023/all-hum-medusa-db-genmoic-protein \
-query output.fasta -out diamondout.csv \
-outfmt "6 qseqid sseqid qlen slen qstart qend sstart send evalue bitscore length pident" \
-evalue 1e-5 -max_target_seqs 50000 -num_threads 8
#2.1：blast输出结果标准化(id>70,cov>30)  2-3-normalization.py  
#2.2:hmmer输出结果标准化
#3：制作表格 4-form.py
#4: 获取经过序列比对（blast、hmm）筛选得到的基因组的蛋白序列信息、可编码蛋白区域信息 \
    python3 5.find-extract-cds-protein.py --num_threads 8  blast-MurC.csv   需要内存10G
#5：计算pfam：python3 ~/software/pfam_scan-main/pfam_scan.py output_protein.fasta /public/group_share_data/TonyWuLab/zzd/db/pfamdb/ -out 111-pfam-py.csv -cpu 8 #注意，pfam_scan.py效果比pfam_scan.pl的好，出来的csv比pl的好操作，以后就用这个
#6:提取氨基酸活性位点前24后25bp的序列（共150nt，50bp）-----需要手动筛选
#7：通过多序列比对得到目标序列上下游150nt的比对信息   python3 7-clusato.py EcMurC-383-50bp.fasta
#8：按照给定条件进行筛选，将那些满足要求的蛋白序列（完整的）与核苷酸序列（150nt）与基因组序列提取出来