# 可选：使用项目 data/reference/ 下的单细胞参考执行 RCTD。
script <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
project <- dirname(dirname(script))
source(file.path(project, "../../src/R/project_config.R"))
load_project(project)
run_step06_rctd()
