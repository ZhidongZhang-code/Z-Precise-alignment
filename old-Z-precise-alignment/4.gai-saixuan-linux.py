import pandas as pd
#from Bio import Entrez
import re
import os
import sys

#0                   1         2        3        4       5       6      7            8           9        10     11      
#qseqid              sseqid    qlen    slen    qstart  qend    sstart  send       evalue      bitscore length positive
#GCF_000190995.1     K00097    329     329     1       329     1       329     1.9e-182        644.8   329     325

inputfile = sys.argv[1]
outputfile = sys.argv[2]

fw = open("./tmp-1","w",encoding='utf-8')
map = {}
map2 = {}
with open(inputfile,'r',encoding='utf-8') as fr:
    for line in fr:
        sp1 = line.strip("\n").split("\t")
        if float(sp1[8])<0.00001 : #e值
            sp2 = line.split('\t')  #GCF_000190995.1     K00097
            sp4 = sp2[1] + '_**_' + sp2[0]  #K00097_**_GCF_000190995.1
            if sp2[1] not in map :   #筛选出来共有多少KO 
                map[sp4] = ''
                fw.write(sp4+'\n')
#map中包含含有这个KO的所有物种

map2 = {}

fw2 = open("./tmp-2","w",encoding='utf-8')
with open("./tmp-1",'r',encoding='utf-8') as fw:
    for line in fw:
        sp = line.strip("\n").split("_**_")
        sp1 = sp[0]
        sp2 = sp[1]
        sp3 =sp1 + '_**_' + sp2
        if sp3 not in map2: #筛选去掉重复的
            map2[sp3] = ''
            fw2.write(sp3+'\n')

i = 0
index_list = []
column_list = []
index_columns = {}
tmp = {}
for line in map2:
    i += 1
    sp = line.strip('\n').split('_**_')  ##K00097_**_GCF_000190995.1  
    index_list.append(sp[0])
    column_list.append(sp[1])
    if sp[0] not in index_columns:
        index_columns[sp[0]] = []
    index_columns[sp[0]].append(sp[1])

df = pd.DataFrame('', index=list(set(index_list)), columns=list(set(column_list)))   

i = 0
for index, column in zip(index_list, column_list):
    df.loc[index, column] = 'yes'
    i += 1
df2 = pd.DataFrame(df.values.T, index=df.columns, columns=df.index)
tmp = {}
with open('/public/group_share_data/TonyWuLab/zzd/db/medusa-2023/all-hum-medusa-db-GCF-name.txt','r',encoding='utf-8') as fr:
    for line in fr:
        sp1 = line.strip('\n').split('\t')
        if sp1[0] not in tmp:
            tmp[sp1[0]] = sp1[1]
df2['baname'] = ''
for index,row in df2.iterrows():
    df2.loc[index,"baname"] = tmp.get(index)
c = df2.pop("baname")
df2.insert(0, 'baname', c)
df2["sum"] = (df2 == "yes").sum(axis=1)
d = df2.pop('sum')
df2.insert(1, 'sum', d)
df2.to_csv(outputfile,sep='\t',index=True,header=True)
os.remove("./tmp-1")
os.remove("./tmp-2")

