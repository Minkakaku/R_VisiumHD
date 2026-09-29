# 初始化项目目录并检查 YAML 配置；不运行分析。
script <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
project <- dirname(dirname(script))
source(file.path(project, "../../src/R/project_config.R"))
load_project_config(project)
message("Project ready: ", PROJECT_DIR)
