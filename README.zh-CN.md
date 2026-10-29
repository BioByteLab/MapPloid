# MapPloid 0.10.4

<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="assets/brand/mapploid-wordmark-white.svg">
    <source media="(prefers-color-scheme: light)" srcset="assets/brand/mapploid-wordmark-transparent.svg">
    <img src="assets/brand/mapploid-wordmark.svg" width="580" alt="MapPloid">
  </picture>
</p>

[English](README.md) | **简体中文**

通过参考基因组高可比对区域中的等位基因频率推断植物倍性。

[图标文件与使用说明](assets/brand/README.zh-CN.md)

## 安装

```bash
conda env create -f environment.yml
conda activate mapploid
python -m pip install --no-deps .
mapploid --version
```

环境包含 SAMtools、GenMap 和 Python 依赖。
许可证：[MIT](LICENSE)。[软件可用性](docs/availability.zh-CN.md)。

## 快速开始

```bash
conda activate mapploid
./mapploid prepare --reference genome.fa --kmer 125 --out reference_assets
./mapploid call --bam sample.bam --assets reference_assets --out sample_result
./mapploid call --bam-list bams.txt --assets reference_assets --out batch_results
```

在浏览器中打开 `report.html`。报告可离线使用，图形已内嵌。
后续分析使用 `result.tsv`。

需要 Python 3.10+、NumPy、SciPy、pysam、matplotlib，以及 PATH 中的 SAMtools 和 GenMap。
在本版本目录运行 `python -m pip install .` 安装，或在已有依赖的环境中使用本目录的 `mapploid` 启动器。
`--version` 显示 0.10.4。

## 一键示例（无需下载）

激活 MapPloid 环境后，在本版本目录运行：

```bash
python examples/quickstart/run.py
```

示例生成 36 kb 人工参考序列及已排序、已索引的合成三倍体 BAM，利用提供的人工 mappability 分数准备参考资源，
再执行公开的 `call` 命令。使用默认过滤、25,000 个 SNP 的最低门槛和 100 次 bootstrap。
脚本检查结果为 `CALLED 3x` 后正常退出。示例无需运行 GenMap，也无需联网；Python 依赖和 SAMtools 必须事先安装。
提供的分数是测试数据，并非对全 A 人工参考序列测得的 mappability。该示例验证软件流程，不评价生物学准确率，也不验证 GenMap。

结果位于 `examples/quickstart/output/result/report.html` 和 `result.tsv`；
PDF/PNG 位于 `examples/quickstart/output/result/figures/`。
重复运行相同命令会复用已完成结果；添加 `--out /path/to/new-example` 可重新运行。
[快速入门教程](examples/quickstart/README.zh-CN.md) 包含快速测试、参考准备、单样本判定和批量分析。

## 中英文说明

所有说明均提供英文与简体中文版本，页首可切换语言。
命令、代码及代码注释使用英文；命令行帮助、运行日志、报告标签、字段和诊断标识也使用英文，便于检索及程序处理。
对应的中文解释见 [命令行参考](docs/CLI.zh-CN.md) 和 [报告阅读指南](docs/REPORT_GUIDE.zh-CN.md)。

## 结果

```text
sample_result/
  result.tsv
  report.html
  regions.bed
  details/
    run.json
    run.log
```

一个小型隐藏完成标记及目录锁用于校验结果和续跑。

| result.tsv 列名 | 含义 |
|---|---|
| sample | 样本名称 |
| predicted_ploidy | 2x、3x、4x 或 UNCERTAIN；运行失败时为 NA |
| status | CALLED、UNCERTAIN 或 FAILED |
| selected_sites | 主模型实际拟合使用的 SNP 数 |
| delta_bic | 次优模型 BIC 减去最优模型 BIC |
| analyzed_region_mb | 最终分析区域大小，单位 Mb |

缺失或未计算的表格值为 NA；JSON 中非有限数值使用 null。
机器可读的数值保留原精度。输出结构版本为 **3**，与软件版本独立。
旧版结果目录不能直接复用。

报告包含：

- 最终判定，以及三项判定要求：至少 25,000 个拟合 SNP、delta BIC 至少为 10、三个主模型全部收敛。
- Bootstrap 稳定性、有效染色体数量及一致率，以及观测分布与拟合分布的形状差异（JS 散度）。
- 一张三行拟合图：分别比较 2x、3x、4x 模型的观测与拟合次要等位基因频率分布，标注相对 BIC 和最佳候选。
  三行使用实际拟合位点集、相同分箱及统一坐标尺度。
- 可展开的详细警告、各参考序列的拟合结果和实际过滤参数。

诊断警告不会否决倍性判定。Bootstrap 支持度是重采样一致率，并非判定正确的概率。
UNCERTAIN 结果可能有最佳候选，报告只将其标记为候选。
没有可用位点或覆盖时，报告说明原因，不生成空图。
没有有效染色体时，一致率不可评估，而不是 0%。

`details/run.json` 记录输入身份、请求与实际设置、随机种子、各阶段摘要、模型拟合、诊断指标和绘图数据。
`reason`、`warnings` 保留在该文件的 `result` 下及报告的折叠诊断中，不占用主结果表的列。
`details/run.log` 记录执行消息及异常详情。
BAM 身份由路径、大小和修改时间标识；索引和参考资源清单另有 SHA256 校验值。
这不等于对整个 BAM 计算校验值。

### 可选输出

```bash
# Preserve stage tables, bootstrap replicates and intermediate site data.
./mapploid call --bam sample.bam --assets reference_assets --out detailed --keep-details

# Also export the three-row evidence figure in PDF and 300-dpi PNG.
./mapploid call --bam sample.bam --assets reference_assets --out exported --export-figures
```

默认精简模式在运行中保留阶段检查点，只有最终输出写入并通过校验后才清理阶段目录。
中断或失败的任务保留检查点；正常完成但为 UNCERTAIN 的分析也执行精简清理。
`--keep-details` 可保留这些目录。

相同命令可用于续跑：已完成结果通过校验后直接复用，未完成阶段从保留的检查点继续。
输入、参数、输出选项或软件版本改变时，需要新输出目录。同一输出目录只能由一个进程使用。
清理样本结果不会删除参考资源。

## 批量分析

`bams.txt` 每行一个 BAM 路径，空行和以 `#` 开头的行会被忽略。
相对路径以列表文件所在目录为基准。样本必须使用同一参考，且 BAM 文件名不同。
样本名可包含字母、数字、点、下划线和连字符。

```text
batch_results/
  result.tsv       Same six columns, plus relative result_file
  report.html      Search/filter samples and open individual reports
  samples/
    sample_A/
    sample_B/
  details/
    run.json
    run.log
```

总表在单样本六列基础上增加相对路径 `result_file`；总览支持搜索、筛选和打开单样本报告。
每处理完一个样本即更新总览。运行失败记为 FAILED，其余样本继续；存在失败样本时退出码为 1。
UNCERTAIN 表示分析已完成，并非运行失败。中断时退出码为 130。
补齐缺失 BAM 后可重复执行相同批量命令以重试失败样本。
报告中的倍性统计是预测结果数量，不是准确率。

## 输入及过滤

BAM 必须按坐标排序，并带 BAI 或 CSI 索引。参考序列名称和长度必须匹配。
`--assets` 是 `prepare --out` 生成的目录，不是 FASTA 或 BED 文件。

参考准备默认值：`--kmer 150`、GenMap 允许 2 个错配、`--min-region-length 10000`、`--threads 8`。
可重复指定 `--exclude-bed`；`--genmap-bedgraph` 可提供参数匹配的既有 GenMap 分数。
扣除参考缺口及排除区域后保留严格分数为 1 的区域；短区域不会跨缺口拼接。
准备参数改变时使用新资源目录。原始 FASTA 必须保持可访问。

本版本保留 manuscript v20 源代码的读长检查：最多检查 10,000 条主比对记录，
若读长小于参考资源的 k-mer 长度，则报输入错误。

| call 选项 | 默认值 |
|---|---|
| --allele-min / --allele-max | 0.15 / 0.85 |
| --min-depth / --max-depth | 自动，按分析阶段估计 |
| --min-mapq / --min-baseq | 20 / 20 |
| --min-allele-reads | 2 |
| --pair-purity | 0.95 |

等位基因频率以计数最高的两个等位基因计数之和为分母，两者均需满足指定频率区间。
位点深度使用质量过滤后的 A/C/G/T 总计数，不等于整个样本的测序深度。
深度上下界均包含边界值。实际边界见报告及运行记录。

## 判定及诊断规则

分析从约 20 Mb 开始；未满足三项判定要求时扩展至约 50 Mb 和 100 Mb。
较小参考使用可用区域。主拟合要求至少 25,000 个 SNP，最多使用 100,000 个。
超过上限时，按各染色体位点数比例分配配额并进行确定性抽样。

三项判定失败分别对应 LOW_SITE_COUNT、LOW_DELTA_BIC、MODEL_NOT_CONVERGED，结果为 UNCERTAIN。
NO_ELIGIBLE_REGIONS、NO_COVERAGE、NO_VALID_SITES 也产生 UNCERTAIN，且没有可用的模型拟合。
运行错误产生 FAILED。

优先诊断项：

- LOW_BOOTSTRAP_SUPPORT：100 次计划区块 bootstrap 中的一致率低于 0.85；未得到明确结果的重复仍计入分母。
- INSUFFICIENT_CHROMOSOME_EVIDENCE：有效参考序列少于 3 条。
- CHROMOSOME_HETEROGENEITY：有效参考序列对全局候选的一致率低于 70%。
  有效序列要求至少 200 个拟合 SNP、模型收敛且 BIC 有限、delta BIC 至少为 10。
- POOR_THEORETICAL_FIT：JS 散度大于 0.05 或不是有限值。

详细诊断还保留 TOO_FEW_BLOCKS（有位点的 5 Mb 区块少于 5 个）、RHO_AT_BOUNDARY、
A_AT_BOUNDARY、A_AT_UPPER_BOUNDARY、MODEL_ASSUMPTION_SENSITIVE。这些均仅作为警告。

## 兼容性与开发

0.10.1、0.10.2、0.10.3 的参考资源仍可读取，样本输出需要新目录。
0.10.4 更新展示方式，并移除自动生成的二倍体解释警告；似然、位点选择及三项判定门槛保持不变。
详见 [版本说明](RELEASE_NOTES_0.10.4.zh-CN.md) 和 [验证记录](validation/README.zh-CN.md)。

```bash
PYTHONPATH=src python -B -m unittest discover -s tests -v
```
