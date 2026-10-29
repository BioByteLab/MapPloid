# Quickstart：快速测试与标准流程

[English](README.md) | **简体中文**

以下命令均在 `mapploid_v0.10.4/` 目录执行，先激活已安装依赖的环境：

```bash
conda activate mapploid
```

需要 Python 3.10+、NumPy、SciPy、pysam、matplotlib 和 SAMtools；真实参考基因组的准备还需要 GenMap。首次安装可在本版本目录执行 `python -m pip install .`。

## 1. 一键测试（约半分钟，无需下载）

```bash
python examples/quickstart/run.py
```

脚本依次生成小型数据、准备参考区域、判定倍性，并打印每一步对应的命令。
预期看到 `Quickstart passed: CALLED 3x`，拟合位点数为 **35,255**。
本次实测首次运行 12.3 秒（此前为 24 秒）；复用通常约 1 秒，耗时因机器而异。

- [打开样例报告](output/result/report.html)：查看倍性结果和三行拟合图。
- [主结果表](output/result/result.tsv)：六列，可直接用于后续分析。
- [PDF](output/result/figures/evidence.pdf) / [PNG](output/result/figures/evidence.png)：可导出的图，PNG 为 300 dpi。
- [详细记录](output/result/details/run.json)：参数、诊断原因和警告。

重复运行同一命令会复用结果。需要重新体验完整计算时，指定新的目录：

```bash
python examples/quickstart/run.py --out /tmp/mapploid_demo_new
```

样例是 36 kb 的人工参考序列和 30x 合成 BAM；提供了人工 mappability 分数，因此不运行 GenMap。它用于验证软件流程，不能用于评价生物学准确率。判定使用默认过滤、25,000 位点门槛和 100 次 bootstrap；小参考可能出现 `TOO_FEW_BLOCKS` 等诊断警告。

## 2. 学会标准流程：准备一次参考 → 分析样本 → 查看结果

### 第一步：准备输入

需要与 BAM 比对时一致的 `genome.fa`，以及按坐标排序、带 BAI/CSI 索引的 BAM。参考序列名称和长度必须一致。若 BAM 尚未排序或索引：

```bash
samtools sort -o sample.sorted.bam sample.bam
samtools index sample.sorted.bam
```

后续示例使用 `sample.sorted.bam`；已有合格 BAM 时直接替换为其路径。

### 第二步：准备参考区域（同一参考及参数仅做一次）

以下以 150 bp reads 为例；`--kmer` 不能大于实际使用的 reads 长度，修剪后的短 reads 也需考虑。

```bash
python mapploid prepare --reference genome.fa --kmer 150 --threads 8 --out reference_assets
```

此步骤通过 GenMap 计算 mappability 并生成分析区域。耗时取决于参考基因组大小。`reference_assets/` 是下一步的 `--assets`，可供同一参考的多个样本复用；原始 FASTA 需保留。

### 第三步：分析一个样本

```bash
python mapploid call --bam sample.sorted.bam --assets reference_assets --out sample_result --export-figures
```

打开 `sample_result/report.html`，或读取 `sample_result/result.tsv`。
`--export-figures` 额外生成 `figures/evidence.pdf` 和 `figures/evidence.png`。

| 关键字段 | 怎么读 |
|---|---|
| predicted_ploidy | 2x / 3x / 4x；不确定时为 UNCERTAIN，运行失败时为 NA |
| status | CALLED = 得到判定；UNCERTAIN = 证据不足；FAILED = 运行失败 |
| selected_sites | 实际拟合的 SNP 数 |
| delta_bic | 最优与次优模型的 BIC 差，越大表示区分越明显；不等同于准确率 |

原因和警告见 `details/run.json`，运行错误见 `details/run.log`。更换输入或参数请使用新输出目录；同一命令可用于中断后的续跑。默认保留精简结果，需要中间表时添加 `--keep-details`。

### 多个样本：复用参考，改用 BAM 列表

创建 `bams.txt`，每行一个 BAM 路径，例如：

```text
/data/sample_A.bam
/data/sample_B.bam
```

```bash
python mapploid call --bam-list bams.txt --assets reference_assets --out batch_results --export-figures
```

样本应使用同一参考，且 BAM 文件名不能重复。相对路径以列表文件所在目录为基准。
总表和总览在 `batch_results/`，单样本结果在 `batch_results/samples/SAMPLE/`。

### 用刚生成的小样例手动练习

一键测试后，下面两条命令可以单独执行，帮助理解 `prepare --out` 与 `call --assets` 的关系（会复用已有结果）：

```bash
python mapploid prepare --reference examples/quickstart/output/inputs/reference.fa --kmer 50 --threads 1 --genmap-bedgraph examples/quickstart/output/inputs/scores.bg --out examples/quickstart/output/assets
python mapploid call --bam examples/quickstart/output/inputs/synthetic_3x.bam --assets examples/quickstart/output/assets --out examples/quickstart/output/result --export-figures
```

这里的 `--genmap-bedgraph` 只用于提供样例已有的分数。真实数据通常使用上面的标准 `prepare` 命令自行计算。

进一步阅读：[命令行参考](../../docs/CLI.zh-CN.md) · [报告阅读指南](../../docs/REPORT_GUIDE.zh-CN.md)。
