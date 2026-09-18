# 方法网页与核心文献

检索日期：2026-08-03。正式代码优先采用 Seurat/方法包当前官方接口，并用方法论文限定结果解释。

## Visium HD、sketch 与完整投射

1. Seurat 5.5.1：*Analysis, visualization, and integration of Visium HD spatial datasets with Seurat*  
   https://satijalab.org/seurat/articles/visiumhd_analysis_vignette

   关键依据：

   - `Load10X_Spatial(bin.size = c(8, 16))` 会建立 `Spatial.008um`/`Spatial.016um`；
   - 官方示例使用标准 `NormalizeData()`；
   - 大规模 8 µm 数据使用 `SketchData(method = "LeverageScore")`；
   - `ProjectData()` 将 sketch cluster、PCA 和 UMAP 投射到完整数据；
   - marker 在投射后的完整 `Spatial.008um` assay 上计算；
   - BANKSY 用于空间定义的组织域，官方示例为 `lambda = 0.8`、`k_geom = 50`。

2. Seurat reference：`Load10X_Spatial`  
   https://satijalab.org/seurat/reference/load10x_spatial

3. Seurat reference：`SketchData`  
   https://satijalab.org/seurat/reference/sketchdata

4. Seurat reference：`ProjectData`  
   https://satijalab.org/seurat/reference/projectdata

## 空间变异基因

5. Seurat reference：`FindSpatiallyVariableFeatures`  
   https://satijalab.org/seurat/reference/findspatiallyvariablefeatures

   当前接口支持 `markvariogram` 和 `moransi`。正式模板使用 Moran's I，并明确它衡量空间自相关而非组间差异。

## RCTD

6. Cable DM et al. *Robust decomposition of cell type mixtures in spatial transcriptomics.* Nature Biotechnology, 2022.  
   PMID: 33603203  
   DOI: https://doi.org/10.1038/s41587-021-00830-w

7. spacexr 官方仓库  
   https://github.com/dmcable/spacexr

   当前官方示例仍使用 `SpatialRNA`、`Reference`、`create.RCTD` 和 `run.RCTD`。`doublet_mode = "doublet"` 的主要分类在 `results_df`，细胞类型权重应同时保留。

## BANKSY 空间域

8. Singhal V et al. *BANKSY unifies cell typing and tissue domain segmentation for scalable spatial omics data analysis.* Nature Genetics, 2024.  
   PMID: 38413725  
   DOI: https://doi.org/10.1038/s41588-024-01664-3

9. BANKSY 官方仓库  
   https://github.com/prabhakarlab/Banksy

10. SeuratWrappers BANKSY 接口  
    https://github.com/satijalab/seurat-wrappers/blob/master/docs/banksy.md

BANKSY 的 `lambda` 控制空间邻域权重，`k_geom` 控制局部邻域尺度。正式模板把结果命名为“空间组织域”，避免误写成细胞类型。

## 空间邻域

11. Palla G et al. *Squidpy: a scalable framework for spatial omics analysis.* Nature Methods, 2022.  
    PMID: 35102346  
    DOI: https://doi.org/10.1038/s41592-021-01358-2

该文与常用空间分析框架把空间图、邻域富集、共现和空间统计作为独立分析模块。模板的 R kNN 入口采用相同的“按切片建空间图、比较局部组成与全局背景”思想，但不声称提供置换检验 P 值。

## CellChat

12. Jin S et al. *Inference and analysis of cell-cell communication using CellChat.* Nature Communications, 2021.  
    PMID: 33597522  
    DOI: https://doi.org/10.1038/s41467-021-21246-9

CellChat 是独立配体-受体分析点。空间距离可用于限制候选通讯，但几何邻近本身不能证明通讯，因此不并入 Visium HD 核心聚类文件。

## 本模板据此做出的主要更新

- 使用官方 `NormalizeData()` 默认量级，而不是把平均 UMI 写死为 scale factor；
- marker 改为完整投射 bins + 完整 Spatial assay；
- Harmony UMAP 与 ProjectData 的 PCA 投射模型分开，避免用没有 feature loadings 的 Harmony reduction 做投射；
- RCTD 同时输出类型分配和权重；
- BANKSY、Moran's I、kNN 邻域各自成为独立分析点；
- Xenium segmentation QC、PCF 和 CellChat 保留为来源/独立分析点，不混入 Visium HD 主线。
