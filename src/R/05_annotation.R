# 投射 cluster 的人工生物学注释
# marker 只提供证据；最终标签通过项目 cluster_annotation.tsv 显式记录，便于审计。

# 根据完整 bins 的 cluster_full 生成待填写的注释模板。
write_annotation_template <- function(obj, path = ANNOTATION_PATH) {
  clusters <- sort(unique(as.character(obj$cluster_full)))
  template <- data.frame(cluster = clusters, cell_type = "")
  save_table(template, path)
  invisible(path)
}

# 应用人工 cluster→cell_type 映射，并保存实际使用的映射表。
apply_spatial_annotation <- function(obj, path = ANNOTATION_PATH) {
  if (!file.exists(path)) {
    write_annotation_template(obj, path)
    stop("Edit ", path, " and run the annotation step again.", call. = FALSE)
  }

  mapping <- read.delim(path, stringsAsFactors = FALSE)
  if (any(is.na(mapping$cell_type) | mapping$cell_type == "")) {
    stop("Fill every cell_type in ", path, " before annotation.", call. = FALSE)
  }

  obj$cell_type <- mapping$cell_type[
    match(as.character(obj$cluster_full), as.character(mapping$cluster))
  ]
  if (any(is.na(obj$cell_type))) {
    stop("Some cluster_full values are missing from ", path, ".", call. = FALSE)
  }

  Seurat::Idents(obj) <- "cell_type"
  save_table(
    mapping,
    file.path(ANNOTATION_DIR, "cluster_celltype_mapping.csv")
  )
  obj
}

# 优先使用完整 bins 的投射 UMAP；非 sketch 对象使用聚类阶段记录的 UMAP。
default_annotation_reduction <- function(obj) {
  settings <- obj@misc$spatial_pipeline
  reduction <- settings$full_umap %||% "umap"
  if (reduction %in% Seurat::Reductions(obj)) reduction else "umap"
}

# 保存以生物学标签显示的完整数据 UMAP。
save_annotation_umap <- function(obj) {
  reduction <- default_annotation_reduction(obj)
  p <- Seurat::DimPlot(
    obj,
    reduction = reduction,
    group.by = "cell_type",
    label = TRUE,
    repel = TRUE,
    raster = TRUE
  ) + ggplot2::coord_fixed()

  save_plot(
    p,
    file.path(ANNOTATION_DIR, "UMAP_cell_type"),
    width = 10,
    height = 8
  )
}
