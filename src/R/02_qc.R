# Visium HD bin 质量控制
# 8 µm bin 不是 Xenium 分割细胞，因此这里检查表达复杂度和组织覆盖，不做细胞边界 QC。

# 返回当前 assay 对应的 Seurat metadata 列名。
spatial_qc_columns <- function(assay = VISIUM_HD$assay) {
  list(
    feature = paste0("nFeature_", assay),
    count = paste0("nCount_", assay),
    mt = "percent.mt",
    ribo = "percent.ribo",
    complexity = "log10GenesPerUMI"
  )
}

# 计算线粒体、核糖体比例和每个 UMI 对应的基因复杂度。
add_spatial_qc <- function(obj, assay = VISIUM_HD$assay) {
  Seurat::DefaultAssay(obj) <- assay
  cols <- spatial_qc_columns(assay)
  obj[[cols$mt]] <- Seurat::PercentageFeatureSet(
    obj, pattern = QC_FEATURE_PATTERNS$mitochondrial
  )
  obj[[cols$ribo]] <- Seurat::PercentageFeatureSet(
    obj, pattern = QC_FEATURE_PATTERNS$ribosomal
  )
  meta <- obj[[]]
  obj[[cols$complexity]] <-
    log10(pmax(meta[[cols$feature]], 1)) /
    log10(pmax(meta[[cols$count]], 2))
  obj
}

# 生成默认阈值，并用 qc_thresholds.tsv 中的逐样本设置覆盖。
read_qc_thresholds <- function(sheet, path = QC_THRESHOLD_PATH) {
  thresholds <- data.frame(
    sample = sheet$sample,
    nFeature_min = QC_DEFAULT$nFeature_min,
    nFeature_max = QC_DEFAULT$nFeature_max,
    nCount_min = QC_DEFAULT$nCount_min,
    nCount_max = QC_DEFAULT$nCount_max,
    percent_mt_max = QC_DEFAULT$percent_mt_max,
    log10GenesPerUMI_min = QC_DEFAULT$log10GenesPerUMI_min,
    stringsAsFactors = FALSE
  )
  if (!file.exists(path)) return(thresholds)

  custom <- read.delim(path, stringsAsFactors = FALSE, na.strings = c("", "NA"))
  value_cols <- setdiff(names(thresholds), "sample")
  for (i in seq_len(nrow(custom))) {
    idx <- match(custom$sample[i], thresholds$sample)
    if (is.na(idx)) next
    for (column in intersect(value_cols, names(custom))) {
      value <- custom[[column]][i]
      if (!is.na(value)) thresholds[[column]][idx] <- value
    }
  }
  thresholds
}

# 使用显式逻辑向量过滤一个样本，避免动态 eval(parse())。
filter_visium_hd_sample <- function(obj, threshold) {
  cols <- spatial_qc_columns()
  meta <- obj[[]]
  keep <-
    meta[[cols$feature]] >= threshold$nFeature_min &
    meta[[cols$feature]] <= threshold$nFeature_max &
    meta[[cols$count]] >= threshold$nCount_min &
    meta[[cols$count]] <= threshold$nCount_max &
    meta[[cols$mt]] <= threshold$percent_mt_max &
    meta[[cols$complexity]] >= threshold$log10GenesPerUMI_min
  subset(obj, cells = rownames(meta)[keep])
}

# 汇总过滤前后 bin 数量和主要 QC 指标。
spatial_qc_summary <- function(obj, sample_name, stage) {
  cols <- spatial_qc_columns()
  meta <- obj[[]]
  data.frame(
    sample = sample_name,
    stage = stage,
    bins = ncol(obj),
    median_nFeature = median(meta[[cols$feature]], na.rm = TRUE),
    median_nCount = median(meta[[cols$count]], na.rm = TRUE),
    median_percent_mt = median(meta[[cols$mt]], na.rm = TRUE),
    median_percent_ribo = median(meta[[cols$ribo]], na.rm = TRUE),
    median_log10GenesPerUMI = median(
      meta[[cols$complexity]], na.rm = TRUE
    )
  )
}

# 绘制单样本过滤前/后 QC 分布；超大对象只抽样用于作图。
save_spatial_qc_plot <- function(obj, sample_name, stage) {
  cols <- spatial_qc_columns()
  plot_obj <- if (ncol(obj) > 100000) {
    subset(obj, downsample = 100000)
  } else {
    obj
  }

  p <- Seurat::VlnPlot(
    plot_obj,
    features = c(
      cols$feature, cols$count, cols$mt,
      cols$ribo, cols$complexity
    ),
    ncol = 5,
    pt.size = 0
  ) + patchwork::plot_annotation(title = paste(sample_name, stage))

  save_plot(
    p,
    file.path(QC_DIR, paste0(sample_name, "_", stage, "_qc")),
    width = 18, height = 5
  )
}

# 对所有样本执行 QC、保存图表，并返回过滤后的对象列表。
run_visium_hd_qc <- function(objects, sheet) {
  thresholds <- read_qc_thresholds(sheet)
  filtered <- vector("list", length(objects))
  names(filtered) <- names(objects)
  summaries <- list()

  for (sample_name in names(objects)) {
    obj <- add_spatial_qc(objects[[sample_name]])
    save_spatial_qc_plot(obj, sample_name, "before_filter")
    summaries[[length(summaries) + 1]] <-
      spatial_qc_summary(obj, sample_name, "before_filter")

    threshold <- thresholds[
      match(sample_name, thresholds$sample), , drop = FALSE
    ]
    obj <- filter_visium_hd_sample(obj, threshold)
    save_spatial_qc_plot(obj, sample_name, "after_filter")
    summaries[[length(summaries) + 1]] <-
      spatial_qc_summary(obj, sample_name, "after_filter")
    filtered[[sample_name]] <- obj
  }

  save_table(do.call(rbind, summaries), file.path(QC_DIR, "qc_summary.csv"))
  save_table(
    thresholds,
    file.path(QC_DIR, "qc_thresholds_used.csv")
  )
  filtered
}
