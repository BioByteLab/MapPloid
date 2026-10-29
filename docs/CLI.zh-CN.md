# 命令行参考

[English](CLI.md) | **简体中文**

[主文档](../README.zh-CN.md) · [快速入门](../examples/quickstart/README.zh-CN.md)

以下命令在本版本目录执行；安装后也可直接使用 `mapploid` 命令，无需 `python` 前缀。终端帮助使用英文；中英文说明覆盖相同的公开参数。代码、注释及命令均使用英文。

```bash
python mapploid --help
python mapploid prepare --help
python mapploid call --help
python mapploid --version
python examples/quickstart/run.py --help
```

## 通用选项

| 选项 | 含义 |
|---|---|
| -h, --help | 显示当前命令的帮助并退出 |
| --version | 显示软件版本并退出；放在子命令之前 |

## 准备参考区域

```bash
python mapploid prepare --reference genome.fa --kmer 150 --out reference_assets
```

| 选项 | 必填 / 默认值 | 含义 |
|---|---|---|
| --reference FASTA | 必填 | 未压缩的参考 FASTA，需与 BAM 比对参考一致 |
| --out DIR | 必填 | 参考资源目录，供后续 `call --assets` 使用 |
| --kmer BP | 150 | GenMap k-mer 长度，单位 bp；不能大于实际使用的读长 |
| --min-region-length BP | 10000 | 扣除排除区域后，连续高可比对区域的最小长度 |
| --threads N | 8 | 参考准备使用的线程数 |
| --exclude-bed BED | 无 | 排除 BED 文件中的区域；可重复指定 |
| --genmap-bedgraph FILE | 无 | 复用同一参考、相同 k-mer 长度及 2 个错配设置生成的 GenMap bedGraph 分数 |

未提供分数时运行 GenMap。原始 FASTA 需保持可访问。仅对相同参考及准备参数复用资源，改变这些设置时使用新目录。

## 分析单样本或批量样本

```bash
python mapploid call --bam sample.bam --assets reference_assets --out sample_result
python mapploid call --bam-list bams.txt --assets reference_assets --out batch_results
```

| 选项 | 必填 / 默认值 | 含义 |
|---|---|---|
| --bam FILE | 与 `--bam-list` 二选一 | 一个按坐标排序、带 BAI/CSI 索引的 BAM |
| --bam-list FILE | 与 `--bam` 二选一 | 每行一个 BAM 路径，各 BAM 均需排序并索引 |
| --assets DIR | 必填 | `prepare --out` 创建的目录 |
| --out DIR | 必填 | 单样本输出目录或批量汇总目录 |
| --sample NAME | 去掉 `.bam` 后缀的 BAM 文件名 | 单样本名称，不能与 `--bam-list` 一起使用 |
| --keep-details | 关闭 | 完成后保留阶段表及中间位点数据 |
| --export-figures | 关闭 | 有可用模型证据时，将报告图导出为 PDF 和 300 dpi PNG |

BAM 列表忽略空行和以 `#` 开头的行，相对路径以列表文件所在目录为基准。
BAM 文件名必须不同，且样本使用同一参考。批量单样本结果位于 `samples/SAMPLE/`，批量 TSV 增加 `result_file` 列。

### 可选位点过滤

| 选项 | 默认值 | 含义 |
|---|---|---|
| --allele-min FLOAT | 0.15 | 计数最高的两个等位基因各自的最低频率 |
| --allele-max FLOAT | 0.85 | 计数最高的两个等位基因各自的最高频率 |
| --min-depth N | 自动 | 质量过滤后的最低位点深度，包含边界值 |
| --max-depth N | 自动 | 质量过滤后的最高位点深度，包含边界值 |
| --min-mapq N | 20 | 最低比对质量，范围 0–255 |
| --min-baseq N | 20 | 最低碱基质量，范围 0–255 |
| --min-allele-reads N | 2 | 每个保留等位基因的最少支持 reads 数 |
| --pair-purity FLOAT | 0.95 | 两个主要等位基因在有效碱基总计数中的最低占比 |

频率以两个主要等位基因的计数之和为分母。
位点深度使用质量过滤后的全部 A/C/G/T 计数，不等于整个样本的测序深度。
自动深度边界按分析阶段估计，并记入报告。公开 `call` 命令使用内置模型判定设置，详见 [报告阅读指南](REPORT_GUIDE.zh-CN.md)。

## Quickstart 脚本

```bash
python examples/quickstart/run.py
python examples/quickstart/run.py --out /tmp/mapploid_demo_new
```

`--out DIR` 指定样例工作目录，默认是 `examples/quickstart/output/`。
脚本创建输入、准备参考资源并运行合成三倍体判定。重复执行相同命令可复用结果，使用新目录可完整重算。
`-h` / `--help` 显示脚本帮助。

## 退出码及续跑

| 退出码 | 含义 |
|---|---|
| 0 | 正常完成，包括已完成但判为 UNCERTAIN 的分析 |
| 1 | 运行失败；批量分析中表示至少一个样本失败 |
| 2 | 参数解析器拒绝了无效的命令行参数 |
| 130 | 执行过程中被中断 |

相同命令可继续中断任务或复用校验通过的结果。
输入、参数、输出选项或软件版本改变时，需要新输出目录。
失败或中断的样本保留检查点；缺失 BAM 时，补齐文件后重复执行相同批量命令即可。
执行详情见 `details/run.log`，结构化原因和警告见 `details/run.json`，解释见 [报告阅读指南](REPORT_GUIDE.zh-CN.md)。

## 输入检查

BAM 必须按坐标排序，并带 BAI 或 CSI 索引。参考序列名称和长度必须匹配。
`--assets` 是 `prepare --out` 生成的目录，不是 FASTA 或 BED 文件。

参考准备默认值：`--kmer 150`、GenMap 允许 2 个错配、`--min-region-length 10000`、`--threads 8`。
可重复指定 `--exclude-bed`；`--genmap-bedgraph` 可提供参数匹配的既有 GenMap 分数。
扣除参考缺口及排除区域后保留严格分数为 1 的区域；短区域不会跨缺口拼接。
准备参数改变时使用新资源目录。原始 FASTA 必须保持可访问。

本版本保留 manuscript v20 源代码的读长检查：最多检查 10,000 条主比对记录，
若读长小于参考资源的 k-mer 长度，则报输入错误。

两个主要等位基因均需通过频率区间过滤。样本名可包含字母、数字、点、下划线及连字符。

精简模式仅在最终输出通过校验后清理阶段目录，正常完成的 UNCERTAIN 分析同样执行清理；`--keep-details` 保留阶段数据。同一输出目录只能由一个进程使用，样本清理不会删除参考资源。
