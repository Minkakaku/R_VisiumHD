# 投射后完整 bins 的 marker 分析
# Seurat Visium HD 官方教程在完整 Spatial assay 和投射 cluster 上找 marker，而不是只分析 sketch assay。

# 按 cluster 分层抽样完整 bins，兼顾代表性和运行速度。
prepare_visium_hd_marker_object <- function(obj) {
  assay <- MARKERS$assay
  Seurat::DefaultAssay(obj) <- assay
  Seurat::Idents(obj) <- MARKERS$cluster_col

  marker_obj <- subset(
    obj,
    downsample = MARKERS$downsample_per_cluster
  )
  SeuratObject::JoinLayers(marker_obj, assay = assay)
}

# 在抽样后的完整 Spatial assay 上计算每个投射 cluster 的正向 marker。
find_visium_hd_markers <- function(marker_obj) {
  Seurat::DefaultAssay(marker_obj) <- MARKERS$assay
  Seurat::Idents(marker_obj) <- MARKERS$cluster_col
  Seurat::FindAllMarkers(
    marker_obj,
    assay = MARKERS$assay,
    only.pos = MARKERS$only_pos,
    test.use = MARKERS$test_use,
    min.pct = MARKERS$min_pct,
    logfc.threshold = MARKERS$logfc_threshold
  )
}

# 保存每个 cluster 的 Top marker 表和论文式 dotplot。
save_marker_dotplot <- function(marker_obj, markers) {
  top_markers <- markers |>
    dplyr::filter(p_val_adj < 0.05) |>
    dplyr::group_by(cluster) |>
    dplyr::slice_max(
      order_by = avg_log2FC,
      n = MARKERS$top_n,
      with_ties = FALSE
    ) |>
    dplyr::ungroup()
  features <- unique(top_markers$gene)

  dot <- Seurat::DotPlot(
    marker_obj,
    features = features,
    group.by = MARKERS$cluster_col,
    assay = MARKERS$assay
  )
  plot_data <- dot$data
  save_table(
    plot_data,
    file.path(MARKER_DIR, "top_marker_dotplot_data.csv")
  )
  save_table(
    top_markers,
    file.path(MARKER_DIR, "top_markers.csv")
  )

  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(features.plot, id)
  ) +
    ggplot2::geom_point(
      ggplot2::aes(size = pct.exp, color = avg.exp.scaled)
    ) +
    ggplot2::scale_color_gradientn(
      colours = HONGF_NPG_SEQUENTIAL
    ) +
    ggplot2::labs(
      x = NULL, y = NULL,
      size = "Percent expressed",
      color = "Scaled expression"
    ) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(color = "#E6E6E6", linewidth = 0.32),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(
        color = "#4A4A4A", fill = NA, linewidth = 0.72
      ),
      axis.ticks = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(
        angle = 90, hjust = 1, vjust = 0.5
      ),
      plot.title = ggplot2::element_blank(),
      plot.subtitle = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(3, 3, 3, 3)
    )

  save_plot(
    p,
    file.path(MARKER_DIR, "top_marker_dotplot"),
    width = max(10, length(features) * 0.28),
    height = 7,
    dpi = 600
  )
  invisible(top_markers)
}
