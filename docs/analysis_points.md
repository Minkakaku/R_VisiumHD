# 空间分析点与代码归类

本说明回答两个问题：每个正式结果目录代表什么分析点，以及 `source_catalog/` 的大文件夹内部还混有哪些子分析。

## 一、Visium HD 正式主线

| 分析点 | 输入 | 主要方法 | 核心输出 | 为什么放在这里 |
|---|---|---|---|---|
| 01 平台与 QC | Space Ranger 样本目录/`outs` | 平台文件结构、nFeature/nCount、MT%、ribo%、复杂度 | 平台表、阈值表、过滤前后汇总与 QC 图 | 先排除错误平台和低信息 bins |
| 02 聚类 | 过滤后完整对象 | LogNormalize、HVG、LeverageScore sketch、PCA 拐点、可选 Harmony、ProjectData | sketch/full UMAP、完整投射 cluster | 对 Visium HD 大对象兼顾速度和稀有结构 |
| 03 marker | 完整投射对象 | 每 cluster 最多抽样 1000 bins，`Spatial.008um` Wilcoxon marker | all/top marker、dotplot 和作图表 | 与 Seurat Visium HD 官方教程一致，不把 sketch 当最终表达层 |
| 04 人工注释 | marker + cluster | 显式 cluster→cell type 表 | 映射表、完整 UMAP、每样本 bin 数 | 保留人工判断和版本记录 |
| 05 空间图 | 注释对象 | SpatialDimPlot、SpatialFeaturePlot、q05/q95 色阶 | PNG/PDF、坐标和 metadata | 核对标签与组织结构是否一致 |
| 06 RCTD | 注释对象 + scRNA 参考 | Reference、SpatialRNA、RCTD doublet mode | 主要类型、权重、空间图和 RDS | bin 可能含混合细胞，解卷积比强制单标签更合理 |
| 07 BANKSY | 标准化完整对象 | 邻域平均表达 + 表达梯度、PCA、聚类 | 空间组织域及每样本域数量 | 识别连续组织域，而不是普通细胞类型 |
| 08 空间统计 | 标准化完整对象 | Seurat Moran's I | 全量/Top 空间变异基因和空间图 | 寻找非随机空间聚集的基因 |
| 09 邻域组成 | 注释完整对象 | 每图像独立 kNN、背景校正 log2 enrichment | 邻域 CSV 和热图 | 描述哪些标签在空间上更常相邻 |

## 二、source_catalog 的 7 个主文件夹

| 来源文件夹 | 文件数 | 内部常见子分析 | 对正式 Visium HD 主线的处理 |
|---|---:|---|---|
| `spatial_preprocessing` | 61 | Space Ranger 命令、Xenium 文件复制、图像/qptiff、坐标、8/16 µm 映射、Seurat↔h5ad 转换 | 只复用通用输入思想；文件操作和项目路径保留为来源 |
| `segmentation_qc` | 75 | Xenium 细胞/细胞核边界、转录本数量、面积、分割后 h5ad、QC 图 | 主要属于 Xenium，不进入 Visium HD 正式 src |
| `spatial_clustering` | 31 | Seurat 聚类、Harmony、cNMF/NMF、空间图旋转/翻转、客户图形 | 标准聚类已重写；项目绘图和 cNMF 保留为可选来源 |
| `spatial_annotation` | 50 | SingleR、人工注释、module score、细胞比例、空间基因图、距离 | 人工注释和空间图已正式化；自动注释/打分按项目选择 |
| `rctd_deconvolution` | 86 | RCTD 核心、参考构建、权重后处理、空间图，也有脚本混入 BANKSY/聚类 | 正式主线只保留独立 RCTD；混合脚本拆开理解 |
| `spatial_neighborhood` | 81 | cellular neighborhood、kNN、最近邻、proximity、距离、PCF、BANKSY | BANKSY 与通用 kNN 已正式化；PCF/接触统计按研究问题选用 |
| `spatial_cellchat` | 1 | 空间 CellChat | 属于独立分析点，正式代码使用同级 CellChat 仓库 |

文件夹数量合计为 385。一个脚本可能命中多个子分析，因此子分析统计会重叠。

## 三、容易混淆的代码边界

### segmentation_qc

这里的“分割”主要指 Xenium 的细胞或细胞核边界。Visium HD 8 µm 数据是规则 bins，没有同等含义的单细胞 polygon，所以不能把 Xenium 的面积、边界和 transcript-per-cell QC 直接套用。

### spatial_clustering 与 BANKSY

普通 Seurat 聚类主要按表达相似性分群；BANKSY 把邻域平均和梯度加入特征，主要用于组织域。两者都产生 cluster，但生物学含义不同。

### annotation 与 RCTD

人工 annotation 给 cluster 一个可读标签；RCTD 给每个 bin 估计参考细胞类型权重。前者便于描述，后者更适合处理混合 bins，二者不应互相替代。

### neighborhood、distance、PCF 与 CellChat

- kNN/CN：谁在谁周围；
- nearest/proximity：两类标签之间的最近距离或阈值内邻近；
- PCF：不同距离尺度上的空间共现；
- CellChat：根据配体-受体表达推断通讯。

空间靠近不等于分子通讯，CellChat 结果也不能替代几何邻域统计。

## 四、逐文件说明

[`source_code_guide.tsv`](source_code_guide.tsv) 为 385 个来源文件逐条提供：

- 来源文件及原路径；
- 主要文件夹和 workflow 角色；
- 检测到的 Visium HD/Xenium/混合平台；
- 文件内子分析标签；
- 中文用途说明；
- 建议是正式参考、Xenium 专用参考还是独立分析点；
- 运行前注意事项。

该索引是旁注，不改来源脚本字节，因此 `SOURCE_CATALOG.tsv` 中的 SHA-256 仍可继续校验。
