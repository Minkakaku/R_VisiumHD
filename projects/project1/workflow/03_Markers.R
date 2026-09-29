# 在完整投射 bins 的 Spatial.008um assay 上计算 marker。
script <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
project <- dirname(dirname(script))
source(file.path(project, "../../src/R/project_config.R"))
load_project(project)
run_step03_markers()
