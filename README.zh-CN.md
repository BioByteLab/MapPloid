# MapPloid 0.10.5

<p align="center">
  <img src="assets/brand/mapploid-wordmark.svg" width="580" alt="MapPloid">
</p>

[English](README.md) | **简体中文**

利用参考基因组高可比对区域的等位基因频率推断样本倍性。

## 安装

安装环境： Linux系统、 Python 3.10+、 NumPy 1.26+、
SciPy 1.11+、pysam 0.22+、matplotlib 3.8+、SAMtools 和 GenMap。

```bash
conda env create -f environment.yml
conda activate mapploid
python -m pip install --no-deps .
mapploid --version
```

已有环境可执行 `python -m pip install .`，并确保 SAMtools 和 GenMap 位于 PATH。


## 运行样例

激活环境后，在仓库根目录执行：

```bash
python examples/quickstart/run.py
```

脚本生成合成输入，调用 GenMap 准备参考区域并检查 `CALLED 3x`，无需下载数据。
样例用于检查软件流程，无任何生物学意义。
打开 `examples/quickstart/output/result/report.html`；同目录包含 `result.tsv`，
`figures/` 下提供 PDF/PNG 图。

## MapPloid基础流程

BAM 必须按坐标排序并建立索引，参考序列名称及长度须匹配。
k-mer 长度不能大于实际读长；下例使用 125。

```bash
mapploid prepare --reference genome.fa --kmer 125 --out reference_assets
mapploid call --bam sample.bam --assets reference_assets --out sample_result
mapploid call --bam-list bams.txt --assets reference_assets --out batch_results
```

`bams.txt` 每行一个 BAM 路径，相对路径以列表文件所在目录为基准。
样本须使用同一参考，且 BAM 文件名不同。参考及准备参数相同时可复用参考资源。
重复相同命令可续跑或复用校验通过的结果。
## 结果示例

| 文件 | 用途 |
|---|---|
| `report.html` | 结果报告；批量报告可打开各样本结果 |
| `result.tsv` | 倍性、状态、拟合 SNP 数、delta BIC 和分析区域大小 |
| `regions.bed` | 单样本分析区域 |
| `details/run.json` | 参数、模型拟合、原因及警告 |
| `details/run.log` | 执行消息及错误 |

`CALLED` 表示判为 2x、3x 或 4x；`UNCERTAIN` 表示无法明确判定；`FAILED` 表示运行错误。
添加 `--export-figures` 导出 PDF/PNG，添加 `--keep-details` 保留中间数据。


## 文档

- [快速入门与标准流程](examples/quickstart/README.zh-CN.md)
- [命令、输入要求及过滤参数](docs/CLI.zh-CN.md)
- [结果、判定规则及诊断解释](docs/REPORT_GUIDE.zh-CN.md)
- [更新记录与兼容性](CHANGELOG.zh-CN.md)
