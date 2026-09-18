# 空间变异基因与 Moran's I
# 该分析寻找在组织坐标上呈现非随机聚集的基因，补充普通 cluster marker。

# 对一个图像运行 Seurat Moran's I，并返回完整统计表。
spatial_variable_feature_table <- function(obj, image) {
  assay <- VISIUM_HD$assay
  features <- head(
    Seurat::VariableFeatures(obj, assay = assay),
    SPATIAL_STATS$nfeatures
  )

  analysed <- Seurat::FindSpatiallyVariableFeatures(
    obj,
    assay = assay,
    layer = "data",
    features = features,
    image = image,
    selection.method = SPATIAL_STATS$method,
    nfeatures = length(features),
    verbose = FALSE
  )
  info <- Seurat::SVFInfo(
    analysed,
    method = SPATIAL_STATS$method,
    status = TRUE,
    assay = assay
  )
  info$gene <- rownames(info)
  info$image <- image
  info[order(info$rank), , drop = FALSE]
}

# 逐图像计算空间变异基因，并保存 Top gene 空间表达图。
run_spatial_statistics <- function(obj) {
  tables <- lapply(Seurat::Images(obj), function(image) {
    log_msg("Spatial variable features: ", image)
    table <- spatial_variable_feature_table(obj, image)
    top_genes <- head(
      table$gene[!is.na(table$rank)],
      SPATIAL_STATS$top_n
    )
    save_spatial_feature_plots(
      obj,
      features = top_genes,
      images = image,
      output_dir = SPATIAL_STATS_DIR,
      stem_prefix = "top_spatial_features"
    )
    table
  })
  do.call(rbind, tables)
}

# 保存 Moran's I 表、Top gene 表和解释说明。
save_spatial_statistics_outputs <- function(table) {
  save_table(
    table,
    file.path(
      SPATIAL_STATS_DIR,
      "spatial_variable_features_all.csv"
    )
  )
  top <- table |>
    dplyr::filter(!is.na(rank)) |>
    dplyr::group_by(image) |>
    dplyr::slice_min(
      order_by = rank,
      n = SPATIAL_STATS$top_n,
      with_ties = FALSE
    ) |>
    dplyr::ungroup()
  save_table(
    top,
    file.path(
      SPATIAL_STATS_DIR,
      "spatial_variable_features_top.csv"
    )
  )

  save_results_note(
    file.path(
      SPATIAL_STATS_DIR,
      "README_results_CN.md"
    ),
    c(
      "# Visium HD 空间变异基因结果说明",
      "",
      paste0(
        "- 方法：Seurat FindSpatiallyVariableFeatures；",
        SPATIAL_STATS$method, "。"
      ),
      paste0(
        "- 每个图像最多检测 ",
        SPATIAL_STATS$nfeatures, " 个高变基因。"
      ),
      "- rank 越小表示空间聚集模式越突出。",
      "- 该结果描述空间自相关，不等同于组间差异表达。",
      "- 正式组间比较应以样本为重复单位，不能把 bins 当作生物学重复。"
    )
  )
}
