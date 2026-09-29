# 平台判定：Visium HD 与 Xenium

## 正式入口结论

本仓库的正式 R 入口定位为 **10x Visium HD**，不是 Xenium。

用户提供的原始 `spatial_dat_analysis.R` 有以下明确证据：

- 输入目录为 `spaceranger_output/<sample>/outs`；
- 读取 `outs/binned_outputs/square_008um/`；
- 调用 `Load10X_Spatial(..., bin.size = 8)`；
- 表达 assay 为 `Spatial.008um`；
- 图像对象为 `VisiumV2`；
- 使用 `SketchData()`、`ProjectData()` 和 RCTD 处理大量 8 µm bins。

这些都是 Visium HD/Space Ranger binned output 特征。Xenium 通常读取 `cell_feature_matrix.h5`、`cells.csv.gz`、`transcripts.*` 和细胞边界，不使用 `Spatial.008um` 或 `square_008um`。

## 来源目录组成

`SOURCE_CATALOG.tsv` 共 385 条，是混合来源素材：

- 来源路径含 Xenium：235 条；
- 来源路径含 Visium 或 Space Ranger：67 条；
- 其余通用空间分析路径：83 条；
- 代码内容检索到 Xenium 特征：186 份；
- 代码内容检索到 Visium HD 特征：92 份。

因此，`source_catalog/` 不能整体视为一个平台的正式流程。正式 `src/` 只采用 Visium HD 逻辑；Xenium 代码保留用于溯源，正式 Xenium 分析继续使用同级 `Python_Xenium`。
