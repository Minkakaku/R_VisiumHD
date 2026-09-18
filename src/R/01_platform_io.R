# 平台识别与 Visium HD 输入
# 本文件只负责识别 Space Ranger/Xenium 目录、读取样本表和合并样本。

# 接受样本目录或直接的 outs/；如果子目录 outs/ 存在则自动进入。
resolve_outs_dir <- function(path) {
  path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  candidate <- file.path(path, "outs")
  # basename(path) 会提取路径的最后一层文件夹名
  if (basename(path) == "outs" || !dir.exists(candidate)) path else candidate
}

# 根据官方输出文件结构区分 Visium HD、Xenium 和普通 Visium。
detect_spatial_platform <- function(path) {
  outs <- resolve_outs_dir(path)
  binned <- file.path(outs, "binned_outputs")
  xenium_files <- c(
    "cell_feature_matrix.h5", "cells.csv.gz",
    "transcripts.csv.gz", "transcripts.zarr.zip"
  )

  if (dir.exists(binned) &&
      any(grepl("^square_[0-9]+um$", list.files(binned)))) {
    return("visium_hd")
  }
  if (any(file.exists(file.path(outs, xenium_files)))) return("xenium")
  if (dir.exists(file.path(outs, "spatial")) &&
      file.exists(file.path(outs, "filtered_feature_bc_matrix.h5"))) {
    return("visium")
  }
  "unknown"
}

# 读取样本、分组和 Space Ranger 路径；相对路径以项目目录为起点。
read_spatial_sample_sheet <- function(path = SAMPLE_SHEET_PATH) {
  sheet <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  required <- c("sample", "group", "spaceranger_dir")
  if (!all(required %in% names(sheet))) {
    stop(
      "samples.tsv requires: ", paste(required, collapse = ", "),
      call. = FALSE
    )
  }
  is_absolute <- grepl("^([A-Za-z]:[/\\\\]|/)", sheet$spaceranger_dir)
  sheet$spaceranger_dir[!is_absolute] <- file.path(
    PROJECT_DIR,
    sheet$spaceranger_dir[!is_absolute]
  )
  sheet
}

# 在加载大矩阵前记录每个样本的平台判定，防止误把 Xenium 混入主流程。
audit_sample_platforms <- function(sheet) {
  data.frame(
    sample = sheet$sample,
    spaceranger_dir = sheet$spaceranger_dir,
    detected_platform = vapply(
      sheet$spaceranger_dir, detect_spatial_platform, character(1)
    ),
    stringsAsFactors = FALSE
  )
}

# 读取一个样本的 Space Ranger binned output，并选择 8 µm assay。
load_visium_hd_sample <- function(sample_row) {
  sample_name <- as.character(sample_row[["sample"]])
  group_name <- as.character(sample_row[["group"]])
  outs_dir <- resolve_outs_dir(sample_row[["spaceranger_dir"]])
  detected <- detect_spatial_platform(outs_dir)

  if (detected != "visium_hd") {
    stop(
      sample_name, " detected as ", detected,
      "; this formal entry expects Visium HD.",
      call. = FALSE
    )
  }

  obj <- Seurat::Load10X_Spatial(
    data.dir = outs_dir,
    bin.size = VISIUM_HD$bin_size,
    slice = sample_name
  )
  Seurat::DefaultAssay(obj) <- VISIUM_HD$assay
  obj[[VISIUM_HD$sample_col]] <- sample_name
  obj[[VISIUM_HD$group_col]] <- group_name
  obj$platform <- PLATFORM

  # 加样本前缀保证合并后每个 bin 的条形码唯一。
  Seurat::RenameCells(obj, add.cell.id = sample_name)
}

# 按样本表逐行加载，返回以样本名命名的 Seurat 对象列表。
load_visium_hd_samples <- function(sheet) {
  objects <- lapply(seq_len(nrow(sheet)), function(i) {
    run_step(
      paste("Load Visium HD", sheet$sample[i]),
      load_visium_hd_sample(sheet[i, , drop = FALSE])
    )
  })
  names(objects) <- sheet$sample
  objects
}

# 合并已经完成 QC 过滤的对象；单样本直接返回原对象。
merge_spatial_samples <- function(objects) {
  if (length(objects) == 1) return(objects[[1]])
  merge(objects[[1]], y = objects[-1], merge.data = TRUE)
}
