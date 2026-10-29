# 报告与诊断阅读指南

[English](REPORT_GUIDE.md) | **简体中文**

[主文档](../README.zh-CN.md) · [命令行参考](CLI.zh-CN.md)

在浏览器打开 `report.html`，可离线使用。报告、图内标签、日志和机器可读标识使用英文，本指南提供对应的中英文解释。
后续分析读取 `result.tsv`，完整诊断值查看 `details/run.json`。

## 先看结果

| 标签 | 含义 |
|---|---|
| CALLED | 三项判定要求均通过，报告倍性为 2x、3x 或 4x |
| UNCERTAIN | 未得到明确判定；显示的最佳候选仍只是候选 |
| FAILED | 运行失败；查看 `details/run.log` 及 `details/run.json` 中的原因 |
| Best candidate | 可用模型中 BIC 最优的候选，不一定得到明确倍性判定 |
| Analyzed region (Mb) | 最终分析区域大小，不是基因组大小 |
| Available qualifying SNPs | 通过位点过滤的可用 SNP 数 |
| Fitting SNPs | 主拟合实际使用的位点数，已考虑最大位点数抽样 |
| Not calculated / NOT_CALCULATED | 指标未计算，不能当作零 |
| Not assessable | 证据不足，无法评估，例如没有有效染色体 |

主 TSV 包含 `sample`、`predicted_ploidy`、`status`、`selected_sites`、`delta_bic`、`analyzed_region_mb` 六列。
批量 TSV 增加 `result_file`。缺失值为 `NA`，JSON 中非有限数值为 `null`。
原因和警告位于各样本 `details/run.json` 的 `result` 下。

## 判定要求

| 报告标签 | 要求 |
|---|---|
| Fitting SNPs | 实际拟合至少使用 25,000 个 SNP |
| BIC separation | 次优模型 BIC 减去最优模型 BIC，至少为 10 |
| Main-model convergence | 三个主模型（2x、3x、4x）全部收敛 |
| PASS / FAIL | 通过 / 未通过该要求 |

未满足要求时，分析从约 20 Mb 扩展至 50 Mb，再到 100 Mb，以可用区域为限。
运行失败与证据不足是不同情况。以下诊断警告均不增加额外判定门槛。

## 诊断证据

| 报告标签 | 解释 |
|---|---|
| Bootstrap stability | 计划区块 bootstrap 重复与全局候选的一致率；未得到明确结果的重复仍计入分母 |
| Valid replicates | 得到明确结果的 bootstrap 重复数与计划重复数 |
| Chromosome consistency | 有效参考序列数量及其对全局候选的一致率 |
| Fit shape (JS divergence) | 观测与拟合分布的形状差异，越小表示越接近 |
| WARNING | 未达到诊断阈值，不会否决倍性判定 |

Bootstrap 支持度不是判定正确的概率。
有效参考序列要求至少 200 个拟合 SNP、模型收敛且 BIC 有限、delta BIC 至少为 10。
没有有效参考序列时无法评估一致率，不代表一致率为 0%。

## 三行证据图

三行分别对应 2x、3x、4x 候选，使用相同拟合 SNP、分箱和坐标尺度。

| 图内标签 | 含义 |
|---|---|
| Observed and fitted distributions | 观测位点频率与模型预测的对比 |
| Observed | 灰色观测分布 |
| Fitted model | 彩色拟合分布 |
| Minor allele fraction | 次要等位基因计数除以两个主要等位基因的计数之和 |
| Sites (%) | 每个分箱中的拟合位点百分比 |
| ΔBIC | 当前模型 BIC 减去最优模型 BIC；最优模型为零 |
| BEST | 拟合模型中的最佳候选，不独立表示已得到明确判定 |
| not converged | 对应模型未收敛 |

图中各模型的 ΔBIC 与 TSV 中唯一的 `delta_bic` 不同，后者比较最优与次优模型。
曲线使用各模型自己的拟合参数、观测深度分布和似然中的离散截断，不做平滑。
`--export-figures` 导出 PDF 和 300 dpi PNG；不导出独立 SVG，离线 HTML 内部仍嵌入矢量图。
没有可用模型证据时不生成空图。

## 原因及警告

### 不确定判定的原因

| 标识 | 解释及建议检查项 |
|---|---|
| LOW_SITE_COUNT | 拟合 SNP 太少，检查覆盖、保留位点数和可用参考区域 |
| LOW_DELTA_BIC | 最优与次优模型区分不足，检查频率分布和区域一致性 |
| MODEL_NOT_CONVERGED | 至少一个主模型未收敛，解释候选前先检查模型诊断 |
| NO_ELIGIBLE_REGIONS | 没有可用分析区域，检查参考准备及排除区域 |
| NO_COVERAGE | 分析区域没有有效 A/C/G/T 观测，检查 BAM、参考及覆盖 |
| NO_VALID_SITES | 存在覆盖但无位点通过过滤，检查实际深度及等位基因过滤参数 |

### 优先诊断警告

| 标识 | 解释 |
|---|---|
| LOW_BOOTSTRAP_SUPPORT | 100 次计划区块 bootstrap 的一致率低于 0.85 |
| INSUFFICIENT_CHROMOSOME_EVIDENCE | 有效参考序列少于 3 条 |
| CHROMOSOME_HETEROGENEITY | 有效参考序列对全局候选的支持比例低于 70% |
| POOR_THEORETICAL_FIT | JS 散度超过 0.05，或无法计算为有限值 |

### 其他诊断警告

| 标识 | 解释 |
|---|---|
| TOO_FEW_BLOCKS | 可用于重采样的有位点 5 Mb 区块少于 5 个 |
| RHO_AT_BOUNDARY | 拟合离散参数接近允许的边界 |
| A_AT_BOUNDARY | 四倍体外侧峰权重接近下界 |
| A_AT_UPPER_BOUNDARY | 四倍体中央峰权重接近零 |
| MODEL_ASSUMPTION_SENSITIVE | 放宽四倍体外侧峰约束后，最佳候选发生变化 |

## 详细部分与批量总览

| 报告标签 | 内容 |
|---|---|
| Detailed diagnostics | 原因、警告，以及可用时各参考序列的拟合结果 |
| Reference-sequence fits | 序列名称、SNP 数、候选、delta BIC 和状态 |
| Inputs and analysis settings | BAM、参考资源、k-mer 长度、最小区域长度、可用区域大小、种子和请求的过滤参数 |
| Actual site filters | 最终分析阶段实际使用的过滤参数 |
| Input notices | 输入相关提示，仅在存在时显示 |
| Result table / Run record / Run log / Analysis regions | 分别链接到 TSV、详细 JSON、运行日志及 BED 区域 |
| Sample overview | 批量进度、各状态数量及已判定倍性的数量 |
| Search sample / All statuses | 按样本名搜索、按状态筛选 |

批量统计表示预测结果，不表示准确率。点击样本名可打开单样本报告。
个别样本失败不会中止其余样本。续跑行为见 [命令行参考中的退出码及续跑说明](CLI.zh-CN.md)。

## 运行记录与兼容性

主拟合使用 25,000–100,000 个 SNP；超过上限时按各染色体位点数比例分配配额并进行确定性抽样。
分析区域大小以 Mb 表示。输出结构版本为 3，旧版结果目录不能直接复用。

`details/run.json` 记录输入身份、请求与实际设置、随机种子、阶段摘要、模型拟合、诊断和绘图数据。
BAM 身份由路径、大小和修改时间标识，索引和参考资源清单另有 SHA256 校验值；这不等于整个 BAM 的校验值。
`details/run.log` 记录执行消息及异常详情，小型隐藏完成标记和目录锁用于结果校验及续跑。
