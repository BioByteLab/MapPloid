# MapPloid 图标

[English](README.md) | **简体中文**

![MapPloid](mapploid-wordmark.svg)

采用已同意的 A 方案：两个圆润的频率峰构成 M，下方配短基线。
图形呼应等位基因频率建模，不表示某一种倍性或实际测得的分布。
颜色为深青 `#087F8C`、暖金 `#D89B45` 和米白 `#F7FAF9`；名称保留 `MapPloid` 大小写。

| 用途 | SVG 源文件 | PNG |
|---|---|---|
| 方形项目图标 | `mapploid-icon.svg` | `mapploid-icon-512.png` |
| 透明图标 | `mapploid-icon-transparent.svg` | `mapploid-icon-transparent-512.png` |
| README 横版标识 | `mapploid-wordmark.svg` | `mapploid-wordmark.png` |
| 透明横版标识 | `mapploid-wordmark-transparent.svg` | `mapploid-wordmark-transparent.png` |
| 浅底黑色单色版 | `mapploid-icon-black.svg`、`mapploid-wordmark-black.svg` | 对应 PNG 文件 |
| 深底白色单色版 | `mapploid-icon-white.svg`、`mapploid-wordmark-white.svg` | 对应 PNG 文件 |

四种图标各提供 32、64、256、512、1024 像素 PNG；四种横版各提供 1160 × 400 像素 PNG。
SVG 由可编辑路径构成，没有嵌入位图；文字由 DejaVu Sans Bold 转为轮廓，显示时无需安装字体。
32 像素下金色点缀较弱，以整体轮廓识别为主。使用时保留比例和留白，深色背景优先使用白色版，头像优先使用米白底方形版。

项目中英文 README 已接入随主题切换的 SVG 横版。文件随 GitHub 仓库提供。

初稿由内置 imagegen 工具生成，完整[提示词](source/imagegen-prompt.txt)已保存。
未修改的原始 PNG 保留在原本地项目归档中；公开仓库提供最终 SVG/PNG 图标。
最终 SVG 使用原生曲线重绘，清理位图边缘杂点、统一配色、缩短金色段并调整居中；它是矢量化设计版本，并非逐像素描摹。

在版本目录执行 `python assets/brand/build.py` 可重新生成，需 matplotlib 和 ImageMagick `convert`。
脚本会覆盖命名的 SVG/PNG 导出文件，保留初稿；仅用于图标构建，不是 MapPloid 运行依赖。

详见[视觉检查记录](review/REVIEW.zh-CN.md)和[图片元数据](review/image_metadata.json)。
