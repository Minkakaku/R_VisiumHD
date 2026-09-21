# 注释后的空间 kNN 邻域组成
# 该步骤回答“某类 8 µm bin 周围更常见哪些标签”，不把 bins 当作独立重复做显著性检验。

# 计算一个图像中 focal label 到 neighbor label 的 kNN 组成和富集倍数。
spatial_knn_enrichment <- function(
  obj,
  image,
  group_by = NEIGHBORHOOD$group_by,
  k = NEIGHBORHOOD$k
) {
  require_optional("RANN")

  coords <- spatial_coordinates(obj, image)
  rownames(coords) <- coords$cell_id
  cells <- intersect(rownames(coords), rownames(obj[[]]))
  labels <- as.character(obj[[]][cells, group_by])
  keep <- !is.na(labels) & labels != ""
  coords <- as.matrix(coords[cells[keep], 1:2, drop = FALSE])
  labels <- labels[keep]

  k_use <- min(k, nrow(coords) - 1L)
  nn_index <- RANN::nn2(
    data = coords,
    k = k_use + 1L
  )$nn.idx[, -1, drop = FALSE]

  label_levels <- sort(unique(labels))
  counts <- matrix(
    0,
    nrow = length(label_levels),
    ncol = length(label_levels),
    dimnames = list(label_levels, label_levels)
  )

  # 按邻居序号累加列联表，避免构建数百万行的边表。
  focal <- factor(labels, levels = label_levels)
  for (j in seq_len(ncol(nn_index))) {
    neighbor <- factor(
      labels[nn_index[, j]],
      levels = label_levels
    )
    counts <- counts + as.matrix(table(focal, neighbor))
  }

  out <- as.data.frame(
    as.table(counts),
    responseName = "neighbor_pairs",
    stringsAsFactors = FALSE
  )
  colnames(out)[1:2] <- c("focal_type", "neighbor_type")
  focal_total <- rowSums(counts)
  background <- table(factor(labels, levels = label_levels))

  out$neighbor_fraction <-
    out$neighbor_pairs / focal_total[as.character(out$focal_type)]
  out$background_fraction <-
    as.numeric(background[as.character(out$neighbor_type)]) / length(labels)
  out$log2_enrichment <- log2(
    (out$neighbor_fraction + NEIGHBORHOOD$pseudocount) /
      (out$background_fraction + NEIGHBORHOOD$pseudocount)
  )
  out$image <- image
  out$sample <- as.character(
    unique(obj[[]][cells[keep], VISIUM_HD$sample_col])[1]
  )
  out$group <- as.character(
    unique(obj[[]][cells[keep], VISIUM_HD$group_col])[1]
  )
  out
}

# 对所有图像分别计算邻域，确保不同切片之间不建立人工邻边。
run_spatial_neighborhood <- function(obj) {
  tables <- lapply(Seurat::Images(obj), function(image) {
    log_msg("Spatial kNN neighborhood: ", image)
    spatial_knn_enrichment(obj, image)
  })
  do.call(rbind, tables)
}

# 保存邻域表、蓝白红富集热图和中文解释。
save_spatial_neighborhood_outputs <- function(table) {
  save_table(
    table,
    file.path(
      NEIGHBORHOOD_DIR,
      "spatial_knn_neighborhood_enrichment.csv"
    )
  )

  p <- ggplot2::ggplot(
    table,
    ggplot2::aes(
      neighbor_type,
      focal_type,
      fill = log2_enrichment
    )
  ) +
    ggplot2::geom_tile() +
    ggplot2::facet_wrap(~image) +
    ggplot2::scale_fill_gradient2(
      low = HONGF_NPG_DIVERGING[1],
      mid = "white",
      high = HONGF_NPG_DIVERGING[3],
      midpoint = 0
    ) +
    ggplot2::coord_fixed() +
    ggplot2::labs(
      x = "Neighbor label",
      y = "Focal label",
      fill = "log2 enrichment"
    ) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 60, hjust = 1
      ),
      strip.background = ggplot2::element_blank()
    )

  save_plot(
    p,
    file.path(
      NEIGHBORHOOD_DIR,
      "spatial_knn_neighborhood_enrichment"
    ),
    width = 10,
    height = 8,
    dpi = 600
  )

  save_results_note(
    file.path(
      NEIGHBORHOOD_DIR,
      "README_results_CN.md"
    ),
    c(
      "# Visium HD kNN 邻域结果说明",
      "",
      paste0(
        "- 邻域标签列：", NEIGHBORHOOD$group_by,
        "；k = ", NEIGHBORHOOD$k, "。"
      ),
      "- 每个图像单独建图，不跨样本连接邻居。",
      "- log2 enrichment > 0 表示比该标签的全图背景比例更常见。",
      "- 热图只描述空间共定位强度，不提供生物学重复层面的显著性。",
      "- 8 µm bin 邻域不等同于 Xenium 单细胞接触或配体-受体通讯。"
    )
  )
}
