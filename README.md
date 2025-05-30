# EXCAT 生物信息学分析工具集

## 工具概述

EXCAT 包含两个主要工具，分别用于不同的生物信息学分析场景：

1. **main_alignment.py** - 综合序列分析流程
   - 整合 BLAST/HMMER/Pfam/Clustal Omega 等多种分析工具
   - 实现从序列比对到多序列比对的完整流程

2. **main_excat.py** - 特异性位点序列提取工具
   - 根据特定氨基酸位点提取蛋白质和CDS序列
   - 支持正选和反选位点模式

## 安装要求

### 系统要求
- Linux/Unix 系统
- Python 3.6 或更高版本

### 依赖安装
```bash
pip install pyyaml
```

### 生物信息学工具依赖
需要预先安装以下工具并配置PATH：
- BLAST+ 或 DIAMOND
- HMMER (hmmsearch, hmmbuild)
- Pfam scan
- Clustal Omega

## 配置文件说明

需创建 `config.yaml` 配置文件，示例内容：

```yaml
# 工具路径配置
blastp_path: /path/to/blastp
diamond_path: /path/to/diamond
db_path: /path/to/blast_db
hmmsearch_path: /path/to/hmmsearch
clustalo_path: /path/to/clustalo

# 参数默认值
coverage_threshold: 0.7
identity_threshold: 30
hmm_evalue: 1e-4
hmm_score: 0

# 模块开关
hmmsearch: "yes"  # 是否使用HMM
pfam: "yes"       # 是否使用Pfam
clustal: "yes"    # 是否使用Clustal
```

## 工具一：main_alignment.py

### 功能描述
实现从原始序列到多序列比对的完整分析流程，包含：
1. 序列相似性搜索(BLAST/DIAMOND)
2. HMMER筛选(可选)
3. Pfam结构域分析(可选)
4. 多序列比对(Clustal Omega)

### 使用示例
```bash
python main_alignment.py --fasta_input query.faa \
       --target_sequence_file target.fasta \
       --normailzationout results.csv \
       --extractcds cds_out.fasta \
       --extractprotein proteins_out.fasta \
       --config config.yaml \
       --num_threads 16
```

### 常用参数说明
| 参数 | 说明 |
|------|------|
| `--fasta_input` | 输入FASTA文件(必需) |
| `--target_sequence_file` | 目标序列文件(150nt) |
| `--normailzationout` | 标准化结果输出文件 |
| `--blast_type` | 选择blastp/blastx/diamond blastp等 |
| `--pfamkey` | 指定Pfam结构域关键词 |
| `--num_threads` | 线程数(默认8) |

## 工具二：main_excat.py

### 功能描述
从多序列比对结果中提取特定氨基酸位点的序列：
- 支持指定多个位点条件
- 可提取对应CDS序列
- 支持正选(25:R)和反选(25!R)模式

### 使用示例
```bash
python main_excat.py --input alignment.fasta \
       --position 25:R 28:R \
       --extractprotein protein_db.fasta \
       --extractcds cds_db.fasta \
       --site_protein_output matched_proteins.fasta \
       --site_cds_output matched_cds.fasta \
       --config config.yaml
```

### 位点语法说明
- `25:R` - 第25位必须是精氨酸(R)
- `25!R` - 第25位不能是精氨酸(R)
- 可同时指定多个位点条件

### 输出文件
1. 位点匹配的蛋白质序列
2. 对应的CDS序列
3. 格式化报告(可选)

## 注意事项

1. 临时文件默认自动清理，可通过`--keep_temp`保留
2. 运行前请确保所有依赖工具已正确安装
3. 建议首次使用时使用小数据集测试参数
4. 详细日志可通过`--log_file`指定输出路径

## 版本信息
当前版本：1.0

如需帮助或发现问题，请联系开发者。