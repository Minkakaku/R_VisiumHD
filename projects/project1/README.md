# project1

这是 Visium HD 模板中的示例项目层。复制为新项目名后，只修改本目录中的配置和 workflow。

- `config/config.yaml`：项目特异的空间展示基因等参数。
- `config/samples.tsv`：样本、分组和 Space Ranger 路径；相对路径以本项目目录为起点。
- `config/qc_thresholds.tsv`：逐样本 QC 阈值。
- `config/cluster_annotation.tsv`：投射 cluster 的人工生物学注释。
- `data/raw/`：Space Ranger 输入，推荐 `raw/<sample>/outs/`。
- `data/processed/`：按步骤编号的 RDS 中间对象。
- `data/reference/`：RCTD 使用的单细胞参考对象。
- `results/`：01–09 各分析点的图、表和中文说明。

从仓库根目录按编号运行 workflow。03 生成 marker 后先复核 `cluster_annotation.tsv`，再运行 04；06 RCTD 和 07 BANKSY 均为可选步骤。
