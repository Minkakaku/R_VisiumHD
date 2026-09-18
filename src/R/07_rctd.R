# RCTD 参考解卷积
# RCTD 用单细胞参考表达谱估计每个 Visium HD bin 的细胞类型组成。

# 从单细胞 Seurat 对象构建 spacexr::Reference。
build_rctd_reference <- function(reference_obj) {
  assay <- RCTD$reference_assay
  reference_obj <- SeuratObject::JoinLayers(
    reference_obj, assay = assay
  )
  counts <- SeuratObject::LayerData(
    reference_obj, assay = assay, layer = "counts"
  )
  labels <- reference_obj[[]][[RCTD$reference_celltype_col]]
  keep <- !is.na(labels) & labels != "Unknown"

  spacexr::Reference(
    counts[, keep],
    factor(labels[keep]),
    Matrix::colSums(counts[, keep])
  )
}

# 提取 RCTD 权重为统一长表；doublet 两列同时关联各自预测细胞类型。
extract_rctd_weights <- function(rctd, image) {
  full <- as.matrix(rctd@results$weights)
  full_long <- data.frame(
    cell_id = rep(rownames(full), times = ncol(full)),
    cell_type = rep(colnames(full), each = nrow(full)),
    component = "all_types",
    weight = as.numeric(full),
    image = image,
    weight_type = "weights",
    stringsAsFactors = FALSE
  )

  # doublet 权重的列表示第一/第二候选，其细胞类型名称来自 results_df。
  doublet <- as.matrix(rctd@results$weights_doublet)
  assignments <- rctd@results$results_df
  doublet_long <- lapply(colnames(doublet), function(component) {
    data.frame(
      cell_id = rownames(doublet),
      cell_type = as.character(assignments[rownames(doublet), component]),
      component = component,
      weight = as.numeric(doublet[, component]),
      image = image,
      weight_type = "weights_doublet",
      stringsAsFactors = FALSE
    )
  })

  dplyr::bind_rows(full_long, doublet_long)
}

# 对每个 Visium HD 图像分别构建 SpatialRNA 并运行 RCTD。
run_rctd_visium_hd <- function(spatial_obj, reference_obj) {
  require_optional("spacexr")
  spatial_obj <- SeuratObject::JoinLayers(
    spatial_obj, assay = VISIUM_HD$assay
  )
  counts <- SeuratObject::LayerData(
    spatial_obj,
    assay = VISIUM_HD$assay,
    layer = "counts"
  )
  reference <- build_rctd_reference(reference_obj)
  results <- list()
  metadata <- list()
  weights <- list()

  for (image in Seurat::Images(spatial_obj)) {
    cells <- intersect(
      SeuratObject::Cells(spatial_obj[[image]]),
      colnames(counts)
    )
    coords <- spatial_coordinates(spatial_obj, image)
    rownames(coords) <- coords$cell_id
    coords <- coords[cells, 1:2, drop = FALSE]
    query_counts <- counts[, cells, drop = FALSE]
    query <- spacexr::SpatialRNA(
      coords,
      query_counts,
      Matrix::colSums(query_counts)
    )

    rctd <- spacexr::create.RCTD(
      query,
      reference,
      max_cores = RCTD$max_cores,
      UMI_min = RCTD$umi_min
    )
    rctd <- spacexr::run.RCTD(
      rctd,
      doublet_mode = RCTD$doublet_mode
    )
    results[[image]] <- rctd

    sample_meta <- rctd@results$results_df
    sample_meta$cell_id <- rownames(sample_meta)
    sample_meta$image <- image
    metadata[[image]] <- sample_meta
    weights[[image]] <- extract_rctd_weights(rctd, image)
  }

  combined <- do.call(rbind, metadata)
  rownames(combined) <- combined$cell_id
  combined$cell_id <- NULL
  spatial_obj <- Seurat::AddMetaData(
    spatial_obj, metadata = combined
  )

  list(
    object = spatial_obj,
    results = results,
    metadata = combined,
    weights = do.call(rbind, weights)
  )
}

# 保存 RCTD 对象、分配表、权重表、空间图和中文说明。
save_rctd_outputs <- function(rctd_output) {
  saveRDS(rctd_output$results, RCTD_RESULT_RDS)

  assignments <- rctd_output$metadata
  assignments$cell_id <- rownames(assignments)
  save_table(
    assignments,
    file.path(RCTD_DIR, "rctd_assignments.csv")
  )
  save_table(
    rctd_output$weights,
    file.path(RCTD_DIR, "rctd_weights.csv")
  )

  save_spatial_annotation_plots(
    rctd_output$object,
    group_by = "first_type",
    output_dir = RCTD_DIR,
    stem_prefix = "spatial_RCTD"
  )

  save_results_note(
    file.path(RCTD_DIR, "README_results_CN.md"),
    c(
      "# Visium HD RCTD 结果说明",
      "",
      paste0("- 单细胞参考对象：", REFERENCE_RDS, "。"),
      paste0(
        "- 参考细胞类型列：",
        RCTD$reference_celltype_col, "。"
      ),
      paste0(
        "- UMI_min = ", RCTD$umi_min,
        "；doublet_mode = ", RCTD$doublet_mode, "。"
      ),
      "- rctd_assignments.csv 保存每个 8 µm bin 的主要类型判断。",
      paste0(
        "- rctd_weights.csv 为长表：cell_type 是细胞类型，",
        "component 区分全类型权重或 doublet 第一/第二候选。"
      ),
      "- RCTD 结果依赖参考对象的细胞类型完整性和测序平台差异。"
    )
  )
  invisible(rctd_output)
}
