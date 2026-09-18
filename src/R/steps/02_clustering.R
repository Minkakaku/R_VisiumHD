# 分析点 02：sketch 聚类、可选 Harmony 和完整 bins 投射

run_step02_clustering <- function() {
  obj <- read_object(OBJECT_QC)
  obj <- run_visium_hd_clustering(obj)
  save_object(obj, OBJECT_CLUSTERED)
  settings <- obj@misc$spatial_pipeline

  # sketch UMAP 显示实际用于建图和聚类的子集。
  p_cluster <- Seurat::DimPlot(
    obj,
    reduction = settings$cluster_umap,
    group.by = "seurat_clusters",
    label = TRUE,
    repel = TRUE,
    raster = TRUE
  ) + ggplot2::coord_fixed()
  p_sample <- Seurat::DimPlot(
    obj,
    reduction = settings$cluster_umap,
    group.by = VISIUM_HD$sample_col,
    raster = TRUE
  ) + ggplot2::coord_fixed()
  save_plot(
    p_cluster,
    file.path(CLUSTER_DIR, "UMAP_sketch_clusters"),
    width = 9, height = 8
  )
  save_plot(
    p_sample,
    file.path(CLUSTER_DIR, "UMAP_sketch_samples"),
    width = 9, height = 8
  )

  # sketch 模式额外输出 ProjectData 后的全部 8 µm bins UMAP。
  if (settings$sketch &&
      settings$full_umap %in% Seurat::Reductions(obj)) {
    p_full_cluster <- Seurat::DimPlot(
      obj,
      reduction = settings$full_umap,
      group.by = "cluster_full",
      label = TRUE,
      repel = TRUE,
      raster = TRUE
    ) + ggplot2::coord_fixed()
    p_full_sample <- Seurat::DimPlot(
      obj,
      reduction = settings$full_umap,
      group.by = VISIUM_HD$sample_col,
      raster = TRUE
    ) + ggplot2::coord_fixed()
    save_plot(
      p_full_cluster,
      file.path(CLUSTER_DIR, "UMAP_full_projected_clusters"),
      width = 9, height = 8
    )
    save_plot(
      p_full_sample,
      file.path(CLUSTER_DIR, "UMAP_full_projected_samples"),
      width = 9, height = 8
    )
  }

  save_results_note(
    file.path(CLUSTER_DIR, "README_results_CN.md"),
    c(
      "# Visium HD 聚类结果说明",
      "",
      paste0(
        "- sketch bins = ", CLUSTER$sketch_ncells,
        "；HVG = ", CLUSTER$nfeatures, "。"
      ),
      paste0(
        "- 自动选择 PC：",
        paste(settings$dims, collapse = ", "), "。"
      ),
      paste0(
        "- 聚类降维：", settings$cluster_reduction,
        "；resolution = ", settings$resolution, "。"
      ),
      "- PCAtools 和 nFactors 结果见 PCA_elbow_selection.csv。",
      "- ProjectData 使用带 feature loadings 的 PCA 投射完整 8 µm bins。",
      "- Harmony 如启用，只用于多样本 sketch 聚类；投射模型仍使用 PCA。"
    )
  )
  invisible(obj)
}
