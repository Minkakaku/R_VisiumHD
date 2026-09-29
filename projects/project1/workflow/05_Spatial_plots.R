# 输出注释与基因表达空间图、坐标和 metadata 表。
script <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
project <- dirname(dirname(script))
source(file.path(project, "../../src/R/project_config.R"))
load_project(project)
run_step05_spatial_plots()
