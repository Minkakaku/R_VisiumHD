# 公共配置

- `default.yaml`：Visium HD 8 µm 主流程的公共 QC、聚类、marker、RCTD、BANKSY、空间统计和邻域参数。
- `reference.yaml`：项目参考对象的默认文件名。

项目差异写入 `projects/<项目>/config/`。其中 `config.yaml` 覆盖公共参数，`samples.tsv` 保存样本与 Space Ranger 路径，`qc_thresholds.tsv` 保存逐样本 QC 阈值，`cluster_annotation.tsv` 保存人工注释。
