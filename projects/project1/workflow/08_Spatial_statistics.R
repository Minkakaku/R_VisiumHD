# 使用 Moran's I 识别空间变异基因。
script <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
project <- dirname(dirname(script))
source(file.path(project, "../../src/R/project_config.R"))
load_project(project)
run_step08_spatial_statistics()
