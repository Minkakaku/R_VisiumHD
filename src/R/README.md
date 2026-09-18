# R/Seurat Visium HD 复用代码

| 文件 | 代码负责什么 |
|---|---|
| `00_runtime.R` | 日志、RDS、CSV/TSV、PNG/PDF 成对输出和中文说明 |
| `01_platform_io.R` | 平台识别、项目样本表和 Space Ranger 8 µm bins 输入 |
| `02_qc.R` | 线粒体/核糖体比例、复杂度、逐样本阈值和过滤 |
| `03_clustering.R` | LogNormalize、SketchData、自动 PC、Harmony、聚类和 ProjectData |
| `04_markers.R` | 完整投射 bins 的 marker、Top marker 表和 DotPlot |
| `05_annotation.R` | 人工 cluster→cell type 映射和完整 UMAP |
| `06_spatial_plots.R` | 空间标签图、基因表达图和坐标/metadata 表 |
| `07_rctd.R` | 单细胞参考、逐图像 RCTD、类型分配和权重 |
| `08_banksy.R` | BANKSY 邻域表达/梯度和空间组织域 |
| `09_spatial_statistics.R` | Moran's I 空间变异基因 |
| `10_neighborhood.R` | 每图像独立的 kNN 标签邻域组成 |
| `project_config.R` | 合并公共/项目 YAML，并建立项目输入输出路径 |
| `steps/*.R` | 将函数组合成 9 个可由 workflow 调用的分析点 |

`load_all.R` 只加载这里的正式代码，不加载 `source_catalog/` 或 `legacy/original/`。阅读时先看根目录/项目 YAML，再看 `steps/`，需要修改算法时再进入对应编号函数文件。
