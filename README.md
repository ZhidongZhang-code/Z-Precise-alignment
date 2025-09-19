# Z-Precise-Alignment

[![Version](https://img.shields.io/badge/version-1.0-blue.svg)](https://github.com/yourusername/Z-Precise-Alignment)
[![Python](https://img.shields.io/badge/python-3.6+-green.svg)](https://www.python.org/)
[![License](https://img.shields.io/badge/license-MIT-yellow.svg)](LICENSE)

## 项目简介

Z-Precise-Alignment 是一个高性能的生物信息学序列比对与功能注释工具套件。该工具专门设计用于精确分析蛋白质序列，通过整合多种比对算法和功能域预测方法，实现对目标序列的全面分析和注释。

### 核心特性

- 🚀 **多引擎序列比对**：支持 BLAST、DIAMOND 等多种比对工具，可根据需求选择最优算法
- 🔍 **精确功能域注释**：集成 HMM、Pfam 等功能域数据库，提供准确的蛋白质功能预测
- 📊 **多序列比对分析**：使用 Clustal Omega 进行多序列比对，识别保守区域和关键位点
- ⚡ **高性能计算**：支持多线程并行处理，优化大规模序列分析性能
- 🔧 **模块化设计**：灵活的模块化架构，支持定制化分析流程
- 📝 **完整的结果追踪**：详细的日志记录和中间文件管理

## 系统要求

### 运行环境
- 操作系统：Linux (推荐 Ubuntu 18.04+) / macOS / Windows (WSL)
- Python：3.6 或更高版本
- 内存：建议 8GB 以上
- 存储：至少 50GB 可用空间（用于数据库存储）

### 必需依赖
```bash
# 基础依赖
conda install -c conda-forge pyyaml pandas numpy biopython

# 序列比对工具
conda install -c bioconda blast diamond clustalo

# 功能域分析工具
conda install -c bioconda hmmer pfam_scan
```

## 快速安装

### 1. 克隆仓库
```bash
git clone https://github.com/yourusername/Z-Precise-Alignment.git
cd Z-Precise-Alignment
```

### 2. 创建 Conda 环境
```bash
conda create -n precise-align python=3.8
conda activate precise-align
```

### 3. 安装依赖
```bash
# 安装 Python 依赖
pip install -r requirements.txt

# 安装生物信息学工具
conda install -c bioconda blast diamond clustalo hmmer
```

### 4. 配置数据库
```bash
# 下载并配置必需的数据库
./scripts/setup_databases.sh
```

## 配置说明

在运行分析前，需要配置 `cfg/config.yaml` 文件：

```yaml
# 工具路径配置
blastp_path: "/path/to/blastp"
diamond_path: "/path/to/diamond"
clustalo_path: "/path/to/clustalo"
hmmsearch_path: "/path/to/hmmsearch"
pfam_scan: "/path/to/pfam_scan.py"

# 数据库路径
db_path: "/path/to/protein_database"
pfam_directory: "/path/to/pfam_database"
GCF_directory: "/path/to/genome_directory"

# 分析参数（可在命令行覆盖）
coverage_threshold: 0.7
identity_threshold: 30
hmm_evalue: 0.0001
```

## 使用指南

### 基础用法

#### 1. 序列比对与注释（main_alignment.py）

最基本的序列比对分析：
```bash
python src/main_alignment.py \
    --fasta_input input.fasta \
    --normailzationout output_normalized.csv \
    --config cfg/config.yaml
```

完整的分析流程（包含所有功能模块）：
```bash
python src/main_alignment.py \
    --fasta_input input.fasta \
    --normailzationout output_normalized.csv \
    --extractcds output_cds.fasta \
    --extractprotein output_protein.fasta \
    --target_sequence_file target.fasta \
    --clustalo_out alignment.fasta \
    --clustalo_extract_out extracted_alignment.fasta \
    --blast_type blastp \
    --num_threads 16 \
    --coverage_threshold 0.8 \
    --identity_threshold 40 \
    --config cfg/config.yaml
```

使用 DIAMOND 加速分析：
```bash
python src/main_alignment.py \
    --fasta_input large_dataset.fasta \
    --normailzationout output.csv \
    --use_diamond \
    --blast_type "diamond blastp" \
    --num_threads 32 \
    --config cfg/config.yaml
```

#### 2. 位点特异性分析（main_excat.py）

提取特定位点的序列：
```bash
python src/main_excat.py \
    --input clustal_alignment.fasta \
    --position 25:R 28:K 100:D \
    --extractprotein proteins.fasta \
    --extractcds cds.fasta \
    --site_protein_output site_proteins.fasta \
    --site_cds_output site_cds.fasta \
    --formatoutput results_table.txt \
    --config cfg/config.yaml
```

### 高级功能

#### 启用 HMM 筛选
在 `config.yaml` 中设置：
```yaml
hmmsearch: "yes"
hmm_evalue: 0.0001
hmm_score: 50
```

#### 启用 Pfam 功能域筛选
```bash
python src/main_alignment.py \
    --fasta_input input.fasta \
    --normailzationout output.csv \
    --pfamkey "Mur_ligase_C" "PBP_dimer" \
    --pfamevalue 1e-10 \
    --config cfg/config.yaml
```

#### 保留中间文件用于调试
```bash
python src/main_alignment.py \
    --fasta_input input.fasta \
    --normailzationout output.csv \
    --keep_temp \
    --config cfg/config.yaml
```

### 参数说明

#### 主要输入参数
- `--fasta_input`: 输入的 FASTA 序列文件
- `--config`: 配置文件路径（默认：config.yaml）

#### 输出参数
- `--normailzationout`: 标准化的比对结果
- `--extractcds`: 提取的完整 CDS 序列
- `--extractprotein`: 提取的完整蛋白质序列
- `--clustalo_out`: Clustal Omega 多序列比对结果
- `--clustalo_extract_out`: 提取的比对片段

#### 比对参数
- `--blast_type`: BLAST 程序类型（blastp/blastx/diamond）
- `--evalue`: E-value 阈值（默认：1e-5）
- `--coverage_threshold`: 覆盖度阈值（默认：0.7）
- `--identity_threshold`: 相似度阈值（默认：30）
- `--max_target_seqs`: 最大目标序列数（默认：5000）

#### 性能参数
- `--num_threads`: 并行线程数（默认：8）
- `--use_diamond`: 使用 DIAMOND 加速
- `--keep_temp`: 保留临时文件

## 输出文件格式

### 1. 标准化比对结果（CSV）
包含以下字段：
- Query ID：查询序列标识
- Subject ID：目标序列标识
- Identity：序列相似度百分比
- Coverage：序列覆盖度
- E-value：期望值
- Bit Score：比特分数

### 2. 提取序列文件（FASTA）
标准 FASTA 格式，包含完整的 CDS 或蛋白质序列

### 3. 多序列比对结果
Clustal 格式的多序列比对文件，可用于下游系统发育分析

## 示例数据

项目提供了示例数据用于测试：
```bash
# 运行示例分析
cd demo
bash make.sh
```

## 模块说明

### 核心模块
- `blast_module.py`：BLAST/DIAMOND 序列比对
- `hmm_module.py`：HMM 模型构建与搜索
- `pfam_module.py`：Pfam 功能域注释
- `run_clustalo_process.py`：多序列比对处理
- `alignment_normalization.py`：结果标准化

### 辅助模块
- `find_extract_cds_protein.py`：序列提取
- `process_hmmresult.py`：HMM 结果处理
- `process_pfamresult.py`：Pfam 结果处理
- `species_checker_format.py`：物种信息注释

## 性能优化建议

1. **大规模数据集**：使用 DIAMOND 替代 BLAST，可提升 100-1000 倍速度
2. **多线程优化**：根据 CPU 核心数设置 `--num_threads`
3. **内存管理**：处理大文件时使用 `--keep_temp` 监控中间文件
4. **数据库索引**：定期更新和索引数据库以保证查询效率

## 故障排除

### 常见问题

1. **数据库路径错误**
   - 检查 `config.yaml` 中的路径配置
   - 确保数据库文件存在且有读取权限

2. **内存不足**
   - 减少 `--max_target_seqs` 参数值
   - 分批处理大型输入文件

3. **依赖缺失**
   - 运行 `python check_dependencies.py` 检查环境
   - 按提示安装缺失的组件

## 开发团队

- 主要开发者：[Your Name]
- 贡献者：[Contributors]

## 引用

如果您在研究中使用了 Z-Precise-Alignment，请引用：

```bibtex
@software{zprecisealign2024,
  title={Z-Precise-Alignment: A High-Performance Sequence Analysis Toolkit},
  author={Your Name},
  year={2024},
  url={https://github.com/yourusername/Z-Precise-Alignment}
}
```

## 许可证

本项目采用 MIT 许可证。详见 [LICENSE](LICENSE) 文件。

## 联系方式

- 项目主页：https://github.com/yourusername/Z-Precise-Alignment
- 问题反馈：[Issues](https://github.com/yourusername/Z-Precise-Alignment/issues)
- 邮件：your.email@example.com

## 更新日志

### v1.0 (2024-10)
- 初始版本发布
- 支持 BLAST/DIAMOND 序列比对
- 集成 HMM 和 Pfam 功能域分析
- 添加多序列比对功能
- 实现模块化架构