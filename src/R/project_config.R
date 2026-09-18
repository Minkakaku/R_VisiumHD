# 加载公共 YAML，并用项目 YAML 递归覆盖同名参数。
# 配置为 NULL 或空值时使用回退值。
.config_or <- function(x, fallback) if (is.null(x) || length(x) == 0) fallback else x
# 按 UTF-8 读取 YAML，空文件视为空配置列表。
.read_yaml <- function(path) .config_or(yaml::read_yaml(path, fileEncoding = "UTF-8"), list())

# 功能：加载一个 Visium HD 项目的参数和输入输出路径。
load_project_config <- function(project_dir) {
  project_dir <- normalizePath(project_dir, winslash = "/", mustWork = TRUE)
  workspace_dir <- normalizePath(file.path(project_dir, "..", ".."), winslash = "/", mustWork = TRUE)
  config_files <- c(
    file.path(workspace_dir, "config", "default.yaml"),
    file.path(workspace_dir, "config", "reference.yaml"),
    file.path(project_dir, "config", "config.yaml")
  )
  # 按公共默认、参考、项目配置的顺序递归覆盖；项目层具有最高优先级。
  config <- Reduce(
    function(current, path) utils::modifyList(current, .read_yaml(path), keep.null = TRUE),
    config_files,
    init = list()
  )

  config_dir <- file.path(project_dir, "config")
  data_dir <- file.path(project_dir, "data")
  data_raw_dir <- file.path(data_dir, "raw")
  data_processed_dir <- file.path(data_dir, "processed")
  data_reference_dir <- file.path(data_dir, "reference")
  results_dir <- file.path(project_dir, "results")
  logs_dir <- file.path(project_dir, "logs")
  report_dir <- file.path(project_dir, "docs", "reports")

  result_dirs <- list(
    qc = file.path(results_dir, "01_qc"),
    cluster = file.path(results_dir, "02_clustering"),
    marker = file.path(results_dir, "03_markers"),
    annotation = file.path(results_dir, "04_annotation"),
    spatial = file.path(results_dir, "05_spatial"),
    rctd = file.path(results_dir, "06_rctd"),
    banksy = file.path(results_dir, "07_banksy_domains"),
    spatial_stats = file.path(results_dir, "08_spatial_statistics"),
    neighborhood = file.path(results_dir, "09_neighborhood")
  )

  cluster <- config$cluster
  cluster$dims <- if (is.null(cluster$dims)) NULL else as.integer(unlist(cluster$dims, use.names = FALSE))
  cluster$regress_vars <- unlist(cluster$regress_vars, use.names = FALSE)
  banksy <- config$banksy
  banksy$dims <- as.integer(unlist(banksy$dims, use.names = FALSE))

  project <- c(list(
    workspace_dir = workspace_dir,
    project_dir = project_dir,
    config_dir = config_dir,
    data_dir = data_dir,
    data_raw_dir = data_raw_dir,
    data_processed_dir = data_processed_dir,
    data_reference_dir = data_reference_dir,
    results_dir = results_dir,
    logs_dir = logs_dir,
    report_dir = report_dir
  ), result_dirs)

  values <- list(
    CONFIG = config,
    WORKSPACE_DIR = workspace_dir,
    PROJECT_DIR = project_dir,
    CONFIG_DIR = config_dir,
    DATA_DIR = data_dir,
    DATA_RAW_DIR = data_raw_dir,
    DATA_PROCESSED_DIR = data_processed_dir,
    DATA_REFERENCE_DIR = data_reference_dir,
    RESULTS_DIR = results_dir,
    REPORT_DIR = report_dir,
    LOGS_DIR = logs_dir,
    QC_DIR = result_dirs$qc,
    CLUSTER_DIR = result_dirs$cluster,
    MARKER_DIR = result_dirs$marker,
    ANNOTATION_DIR = result_dirs$annotation,
    SPATIAL_DIR = result_dirs$spatial,
    RCTD_DIR = result_dirs$rctd,
    BANKSY_DIR = result_dirs$banksy,
    SPATIAL_STATS_DIR = result_dirs$spatial_stats,
    NEIGHBORHOOD_DIR = result_dirs$neighborhood,
    SAMPLE_SHEET_PATH = file.path(config_dir, "samples.tsv"),
    QC_THRESHOLD_PATH = file.path(config_dir, "qc_thresholds.tsv"),
    ANNOTATION_PATH = file.path(config_dir, "cluster_annotation.tsv"),
    REFERENCE_RDS = file.path(data_reference_dir, config$reference_rds_name),
    OBJECT_QC = file.path(data_processed_dir, "01.visium_hd_qc.rds"),
    OBJECT_CLUSTERED = file.path(data_processed_dir, "02.visium_hd_clustered.rds"),
    OBJECT_ANNOTATED = file.path(data_processed_dir, "03.visium_hd_annotated.rds"),
    OBJECT_RCTD = file.path(data_processed_dir, "04.visium_hd_rctd.rds"),
    OBJECT_BANKSY = file.path(data_processed_dir, "05.visium_hd_banksy.rds"),
    RCTD_RESULT_RDS = file.path(data_processed_dir, "04.rctd_results.rds"),
    PLATFORM = config$platform,
    SEED = as.integer(config$seed),
    PACKAGES = unlist(config$packages, use.names = FALSE),
    OPTIONAL_PACKAGES = unlist(config$optional_packages, use.names = FALSE),
    VISIUM_HD = config$visium_hd,
    QC_FEATURE_PATTERNS = config$qc_feature_patterns,
    QC_DEFAULT = config$qc_default,
    CLUSTER = cluster,
    MARKERS = config$markers,
    SPATIAL_FEATURES = unlist(config$spatial_features, use.names = FALSE),
    SPATIAL_PLOT = config$spatial_plot,
    RCTD = config$rctd,
    BANKSY = banksy,
    SPATIAL_STATS = config$spatial_stats,
    NEIGHBORHOOD = config$neighborhood,
    PROJECT = project
  )
  # 将解析后的配置注册到全局环境，供现有 R 步骤函数共同使用。
  list2env(values, envir = .GlobalEnv)
  options(future.globals.maxSize = config$future_globals_max_size_gb * 1024^3)
  invisible(lapply(unique(unlist(project, use.names = FALSE)), dir.create,
    recursive = TRUE, showWarnings = FALSE))
  invisible(project)
}

# 功能：加载项目配置，再注册全部正式 Visium HD 分析函数。
load_project <- function(project_dir) {
  load_project_config(project_dir)
  sys.source(file.path(WORKSPACE_DIR, "src", "R", "load_all.R"), envir = .GlobalEnv)
  invisible(PROJECT)
}
