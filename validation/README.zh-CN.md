# MapPloid 0.10.4 验证

[English](README.md) | **简体中文**

2026-10-29 从干净的发布准备目录检查公开源码，记录见[单元测试](github_tests.log)和
[全新 quickstart](github_quickstart.log)。这些是软件检查，不代表生物学准确率或通用测序深度要求。

激活文档中的环境后，在仓库根目录运行：

```bash
PYTHONPATH=src python -m unittest discover -s tests -v
python examples/quickstart/run.py
python scripts/validate_output.py
```

测试覆盖模型行为、参考准备、BAM 计数、流程复用、单样本及批量处理。
quickstart 使用合成三倍体 BAM 和提供的人工 mappability 分数，检查 `CALLED 3x`，
导出 PDF/PNG 图和离线 HTML 报告，不执行 GenMap。

可选输出验证脚本包含空 BAM 和故意缺失的 BAM，用于检查 `UNCERTAIN`、`FAILED`、
批次继续和输出复用。生成的输入和报告保存在本地，并从 Git 中排除。
旧 `scripts/validate_synthetic.py` 入口由当前 quickstart 和输出验证脚本替代；
历史 v0.10.1 检查不计为 v0.10.4 的新验证。

此前本地验证日志和审计记录继续保留用于追溯，其中绝对路径对应原开发运行。
以上 GitHub 发布准备结果用于检查可移植源码。
图标视觉检查见[审阅记录](../assets/brand/review/REVIEW.zh-CN.md)。

58 项测试通过，零跳过，耗时 61.340 秒。全新 quickstart 在 12.4 秒内得到
`CALLED 3x`，拟合 SNP 为 35,255 个。Wheel 打包通过，见[打包日志](github_wheel.log)。
