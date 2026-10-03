# Availability and requirements for manuscript section 6

| 项目 | 内容 |
| --- | --- |
| 项目名称 | MapPloid |
| 版本 | 0.10.1 |
| 拟发布主页 | https://github.com/BioByteLab/MapPloid （创建并上传成功后才能作为已公开地址） |
| 操作系统 | Linux；Ubuntu 环境验证 |
| 编程语言 | Python |
| 运行依赖 | Python >=3.10、NumPy >=1.26、SciPy >=1.11、pysam >=0.22；GenMap、SAMtools |
| 许可证 | MIT |
| 非学术用途限制 | 无额外非学术用途限制；遵循MIT许可证 |
| 安装、最小示例与预期输出 | README.md、environment.yml、scripts/validate_synthetic.py、validation/README.md |
| 软件存档与持久标识符 | 尚未分配；不能将未创建的DOI填写为已存在 |

本发布包只包含软件源码、测试、文档和合成示例生成器；不包含原始测序数据、真实样本BAM、参考基因组或尚未发表的论文。

软件核心源码来自 v19 配套 v0.10.1，未修改模型或过滤参数。合成验证脚本仅修正了从干净目录运行时创建父目录的问题。历史绝对路径依赖的内部开发脚本未纳入公开包。
