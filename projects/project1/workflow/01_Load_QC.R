# 核对平台，读取 8 µm bins，并按样本执行 QC。
source("src/R/project_config.R")
load_project("projects/project1")
run_step01_load_qc()
