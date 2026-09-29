# 可选：识别连续空间组织域；BANKSY domain 不等同于细胞类型。
script <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
project <- dirname(dirname(script))
source(file.path(project, "../../src/R/project_config.R"))
load_project(project)
run_step07_banksy()
