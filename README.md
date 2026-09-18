# R-visium-hd-workflow-template

正式定位：**R/Seurat Visium HD 多步骤主流程**。它读取 Space Ranger binned output，默认分析 8 µm bins；不是 Xenium 分割细胞流程。

原始 `spatial_dat_analysis.R` 使用 `binned_outputs/square_008um`、`Load10X_Spatial(bin.size = 8)`、`Spatial.008um`、`VisiumV2`、`SketchData()` 和 `ProjectData()`，因此平台明确为 Visium HD。

## 项目结构

```text
.
├── src/R/                         # 可复用函数与 9 个组合分析点
├── config/                        # 公共默认参数和参考文件名
├── projects/project1/
│   ├── config/                    # 样本、QC 阈值和人工注释
│   ├── workflow/                  # 按实际运行顺序编号的短入口
│   ├── data/{raw,processed,reference}/
│   ├── results/                   # 01–09 分析结果
│   ├── logs/
│   └── docs/reports/
├── envs/
├── tests/
├── legacy/original/               # 原始脚本，逐字节保留
└── source_catalog/                # 其他人员来源代码，仅供查阅
```

## 正式分析点

| 顺序 | 结果目录 | 入口 | 作用 |
|---|---|---|---|
| 01 | `01_qc` | `run_step01_load_qc()` | 平台识别、8 µm bins 加载、逐样本 QC 和过滤 |
| 02 | `02_clustering` | `run_step02_clustering()` | LogNormalize、sketch、自动 PC、可选 Harmony、聚类和完整 bins 投射 |
| 03 | `03_markers` | `run_step03_markers()` | 完整投射 bins + `Spatial.008um` marker |
| 04 | `04_annotation` | `run_step04_annotation()` | cluster 到生物学标签的人工映射 |
| 05 | `05_spatial` | `run_step05_spatial_plots()` | 注释和基因表达空间图、坐标与 metadata |
| 06 | `06_rctd` | `run_step06_rctd()` | 可选 RCTD 解卷积及完整权重 |
| 07 | `07_banksy_domains` | `run_step07_banksy()` | 可选 BANKSY 连续组织域 |
| 08 | `08_spatial_statistics` | `run_step08_spatial_statistics()` | Moran's I 空间变异基因 |
| 09 | `09_neighborhood` | `run_step09_neighborhood()` | 每个切片独立的 kNN 邻域组成 |

## 开始项目

1. 复制 `projects/project1/` 并改成项目名。
2. 把 Space Ranger 结果放到项目 `data/raw/<sample>/outs/`。
3. 修改项目 `config/samples.tsv`、`qc_thresholds.tsv` 和 `config.yaml`。
4. 从仓库根目录按编号运行项目 `workflow/`。
5. marker 完成后人工复核 `cluster_annotation.tsv`，再执行注释。

单个入口保持简短：

```r
source("src/R/project_config.R")
load_project("projects/project1")
run_step02_clustering()
```

已经运行过的 workflow 不覆盖；参数变化或补充分析新建下一个编号脚本。模板不提供跨越人工检查点的 `run_all`。

## 环境与测试

```bash
conda env create -f envs/environment.yaml
conda activate r-visium-hd-workflow
Rscript tests/test_project_config.R
```

RCTD 另装 `spacexr`；BANKSY 步骤另装 `Banksy` 和 `SeuratWrappers`。这些依赖只在对应可选步骤执行时检查。

## 方法边界

- marker 使用 ProjectData 后的完整 bins 和完整 `Spatial.008um` assay，不只使用 sketch。
- ProjectData 使用带 feature loadings 的 PCA；Harmony 可用于多样本 sketch 聚类，但不作为投射模型。
- 8 µm bin 不等同于 Xenium 分割细胞。
- BANKSY cluster 表示空间组织域，不直接等同于细胞类型。
- Moran's I 是空间自相关，不是 CON/EXP 差异表达。
- 组间统计必须以样本为重复单位，不能把 bins 当作生物学重复。

详细方法依据见 [docs/method_references.md](docs/method_references.md)，分析点与来源代码边界见 [docs/analysis_points.md](docs/analysis_points.md)。

## 来源与原始代码

`source_catalog/` 的逐文件用途见 [docs/source_code_guide.tsv](docs/source_code_guide.tsv)。正式 workflow 不加载来源代码或 Xenium 脚本。

- `legacy/original/spatial_dat_analysis.R`：SHA-256 `1ca295163196b8369de2f2f0ff90b483fb9b5f5dc8ad8f149f8732ad8eace2a1`
- `legacy/original/muti_group_cellchat_Spatial.R`：SHA-256 `e25e4a906602748efb314665c58df9103ad3a93bfd177b0194bff6127d4e52d7`
