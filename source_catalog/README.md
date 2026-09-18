# source_catalog 来源代码说明

本目录保存 385 份其他人员空间分析代码的逐字节副本，用于查找历史写法，不是一个可整体运行的 workflow。

## 为什么不直接给这些文件插入注释

`SOURCE_CATALOG.tsv` 记录每份文件的原路径、大小和 SHA-256。直接修改代码会破坏溯源，因此：

- 原文件保持不变；
- 每个文件的中文用途、检测平台、子分析标签和使用建议写在 [`../docs/source_code_guide.tsv`](../docs/source_code_guide.tsv)；
- 正式修正版函数放在 `../src/`；
- 用户原始完整脚本放在 `../legacy/original/`。

## 文件夹含义

| 文件夹 | 主要分析点 | 常见混入内容 |
|---|---|---|
| `spatial_preprocessing` | 数据搬运、格式转换、图像和坐标准备 | Space Ranger、Xenium、qptiff、h5ad、8/16 µm 映射 |
| `segmentation_qc` | Xenium 分割与 QC | 表达过滤、h5ad 转换、注释和绘图 |
| `spatial_clustering` | Seurat/Harmony 聚类 | cNMF、图像旋转、空间图 |
| `spatial_annotation` | 注释、比例、打分 | SingleR、marker、距离和空间图 |
| `rctd_deconvolution` | RCTD | BANKSY、聚类、客户图和后处理 |
| `spatial_neighborhood` | 邻域和距离 | BANKSY、CN、kNN、最近邻、proximity、PCF |
| `spatial_cellchat` | 空间 CellChat | 通讯网络和空间参数 |

一个文件可能包含多个分析点。选择代码时不要只看文件夹名，应先查 `source_code_guide.tsv` 的 `detected_topics` 和 `platform_hint`。

## 正式使用规则

- Visium HD 主线优先使用 `../src/`。
- Xenium 分割/细胞边界代码转到同级 Python Xenium 模板参考。
- CellChat 使用独立 `R-cellchat-workflow-template`。
- R Markdown 和 Shell 只作为来源参考；正式运行重新写成简短 `workflow/stepXX.*.R`。
