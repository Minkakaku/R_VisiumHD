# 最小测试：检查 YAML 覆盖、项目路径和样本相对路径解析。
source("src/R/project_config.R")
project <- load_project_config("projects/project1")
sys.source("src/R/01_platform_io.R", envir = .GlobalEnv)
sys.source("src/R/00_runtime.R", envir = .GlobalEnv)
samples <- read_spatial_sample_sheet()
tsv_path <- tempfile(fileext = ".tsv")
save_table(data.frame(sample = "sample1", group = "CON"), tsv_path)
tsv_roundtrip <- read.delim(tsv_path, stringsAsFactors = FALSE)
unlink(tsv_path)

stopifnot(
  basename(project$project_dir) == "project1",
  PLATFORM == "visium_hd",
  VISIUM_HD$bin_size == 8,
  VISIUM_HD$assay == "Spatial.008um",
  identical(BANKSY$dims, 1:30),
  all(startsWith(samples$spaceranger_dir, PROJECT_DIR)),
  identical(tsv_roundtrip$group, "CON"),
  dirname(OBJECT_QC) == DATA_PROCESSED_DIR,
  startsWith(QC_DIR, RESULTS_DIR)
)

cat("project config: OK\n")
