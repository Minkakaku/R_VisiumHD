# R_VisiumHD 示例项目

这是未运行的项目示例。复制整个目录到 `projects/<实际项目名>/`，再修改 config；不复制共享 src。

项目参数覆盖仓库公共默认。请记录实际数据来源、样本分组、参考版本、关键参数、运行时间和结论；数据、图表、日志不提交。样本表以项目 config 中实际加载的表为唯一来源。

步骤：
- `00_setup.R`
- `01_Load_QC.R`
- `02_Clustering.R`
- `03_Markers.R`
- `04_Annotation.R`
- `05_Spatial_plots.R`
- `06_RCTD.R`
- `07_BANKSY.R`
- `08_Spatial_statistics.R`
- `09_Neighborhood.R`

先准备并检查输入；QC 和人工注释后再执行下游。独立分析方法是可选入口，不要求全部运行。详细环境与命令见仓库 README。
