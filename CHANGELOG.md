# 更新日志

本文件记录项目的所有重要更新。

## [1.1.0] - 2025年9月19日

### 新增功能
- 在 `makeformat.py` 中添加新的命令行参数 `--unique_by`，用于基于分类水平的数据去重
  - 支持根据 Species（种）、Genus（属）、Family（科）、Order（目）、Class（纲）或 Phylum（门）进行去重
  - 对于每个唯一的分类水平值，保留 sum 值最高的记录
  - 添加了去重过程的日志记录功能

### 功能改进
- 修改了 `makeformat.py` 的输出格式：
  - 将 GCF 编号从索引转换为第一列
  - 重新组织列顺序为：GCF、strain_Name、Species、Genus、Family、Order、Class、Phylum、sum
  - 改进了 GCF ID 的匹配逻辑，现在可以忽略版本号进行匹配（例如：GCF_000710365.2 可以匹配 GCF_000710365）
  - 优化了 DataFrame 的数据组织方式

### 技术细节
- 新增 `clean_gcf_id()` 函数处理 GCF 版本号
- 使用 pandas DataFrame 操作实现高效的数据去重
- 更新了日志记录功能，可追踪去重结果
- 修改了 CSV 输出格式（移除了索引列）

### 使用示例
```bash
# 按照属（Genus）级别进行去重
python src/makeformat.py --inputfile input.txt --normailzationout output.csv --config cfg/config.yaml --unique_by Genus
```