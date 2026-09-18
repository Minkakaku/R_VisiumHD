# Visium HD sketch 聚类与完整数据投射
# 核心顺序：归一化 → HVG → sketch → PCA/自动拐点 → 可选 Harmony → 聚类 → ProjectData。

# 提取各 PC 的方差、方差比例和累计方差。
pca_variance_table <- function(obj, reduction = "pca") {
  pc_sd <- Seurat::Stdev(obj, reduction = reduction)
  variance <- pc_sd^2
  variance_percent <- variance / sum(variance) * 100
  data.frame(
    pc = seq_along(variance),
    variance = variance,
    variance_percent = variance_percent,
    cumulative_percent = cumsum(variance_percent)
  )
}

# 用 PCAtools 和 nFactors 给出三个自动拐点结果。
detect_pca_elbow <- function(obj, reduction = "pca") {
  variance <- pca_variance_table(obj, reduction)
  scree <- nFactors::nScree(x = variance$variance)
  selections <- data.frame(
    method = c("PCAtools", "nFactors_noc", "nFactors_naf"),
    pc = as.integer(c(
      PCAtools::findElbowPoint(variance$variance_percent),
      scree$Components$noc,
      scree$Components$naf
    ))
  )
  list(
    recommended_pc = selections$pc[selections$method == "PCAtools"],
    selections = selections,
    variance = variance
  )
}

# 保存拐点表、方差表和带三种自动结果的折线图。
save_pca_elbow <- function(obj, path_no_ext, reduction = "pca") {
  result <- detect_pca_elbow(obj, reduction)
  save_table(result$variance, paste0(path_no_ext, "_variance.csv"))
  save_table(result$selections, paste0(path_no_ext, "_selection.csv"))

  method_colors <- c(
    PCAtools = "#B2182B",
    nFactors_noc = "#2166AC",
    nFactors_naf = "#f6ae00"
  )
  p <- ggplot2::ggplot(
    result$variance,
    ggplot2::aes(pc, variance_percent)
  ) +
    ggplot2::geom_line(linewidth = 0.55, color = "#4D4D4D") +
    ggplot2::geom_point(size = 1.6, color = "#4D4D4D") +
    ggplot2::geom_vline(
      data = result$selections,
      ggplot2::aes(
        xintercept = pc, color = method, linetype = method
      ),
      linewidth = 0.7
    ) +
    ggplot2::scale_color_manual(values = method_colors) +
    ggplot2::labs(
      x = "Principal component",
      y = "Explained variance (%)",
      title = "Visium HD PCA variance and automatic elbow",
      subtitle = paste0(
        "PCAtools recommendation: PC", result$recommended_pc
      ),
      color = NULL,
      linetype = NULL
    ) +
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5),
      plot.subtitle = ggplot2::element_text(hjust = 0.5)
    )
  save_plot(p, path_no_ext)
  result
}

# 用户未手动指定 dims 时，以 PCAtools 拐点为中心并限制上下界。
resolve_cluster_dims <- function(elbow) {
  if (!is.null(CLUSTER$dims)) return(as.integer(CLUSTER$dims))
  n_dims <- max(
    CLUSTER$min_dims,
    min(CLUSTER$max_dims, elbow$recommended_pc)
  )
  seq_len(n_dims)
}

# 执行 Seurat 官方推荐的 Visium HD sketch 聚类流程。
run_visium_hd_clustering <- function(obj) {
  full_assay <- VISIUM_HD$assay
  Seurat::DefaultAssay(obj) <- full_assay

  # 官方教程使用标准 LogNormalize；scale factor 保留在 config 中可调整。
  obj <- run_step(
    "NormalizeData",
    Seurat::NormalizeData(
      obj,
      assay = full_assay,
      normalization.method = CLUSTER$normalization_method,
      scale.factor = CLUSTER$scale_factor,
      verbose = FALSE
    )
  )
  obj <- run_step(
    "FindVariableFeatures full assay",
    Seurat::FindVariableFeatures(
      obj,
      assay = full_assay,
      nfeatures = CLUSTER$nfeatures,
      verbose = FALSE
    )
  )

  # 大对象用 leverage-score sketch；小对象直接在完整 assay 上聚类。
  use_sketch <- ncol(obj) > CLUSTER$sketch_ncells
  if (use_sketch) {
    obj <- run_step(
      "SketchData",
      Seurat::SketchData(
        obj,
        ncells = CLUSTER$sketch_ncells,
        sketched.assay = "sketch",
        method = "LeverageScore",
        var.name = "leverage.score",
        seed = SEED
      )
    )
    cluster_assay <- "sketch"
  } else {
    cluster_assay <- full_assay
  }

  # sketch assay 重新选 HVG、缩放和 PCA，供后续聚类使用。
  Seurat::DefaultAssay(obj) <- cluster_assay
  obj <- Seurat::FindVariableFeatures(
    obj, nfeatures = CLUSTER$nfeatures, verbose = FALSE
  )
  regress_vars <- intersect(CLUSTER$regress_vars, colnames(obj[[]]))
  obj <- run_step(
    "ScaleData",
    Seurat::ScaleData(
      obj,
      vars.to.regress = if (length(regress_vars)) regress_vars else NULL,
      verbose = FALSE
    )
  )
  obj <- run_step(
    "RunPCA",
    Seurat::RunPCA(
      obj,
      npcs = CLUSTER$npcs,
      reduction.name = "pca",
      verbose = FALSE
    )
  )

  elbow <- save_pca_elbow(
    obj, file.path(CLUSTER_DIR, "PCA_elbow")
  )
  dims <- resolve_cluster_dims(elbow)
  cluster_reduction <- "pca"

  # 多样本时可在 sketch PCA 上做 Harmony；单样本跳过。
  if (CLUSTER$batch_method == "harmony" &&
      length(unique(obj[[]][[CLUSTER$batch_col]])) > 1) {
    require_optional("harmony")
    obj <- run_step(
      "RunHarmony",
      harmony::RunHarmony(
        obj,
        group.by.vars = CLUSTER$batch_col,
        reduction = "pca",
        dims.use = dims,
        reduction.save = "harmony"
      )
    )
    cluster_reduction <- "harmony"
  }

  # ProjectData 需要带 feature loadings 的 PCA，所以投射模型始终基于 pca。
  obj <- run_step(
    "RunUMAP projection model",
    Seurat::RunUMAP(
      obj,
      reduction = "pca",
      dims = dims,
      reduction.name = "umap",
      return.model = use_sketch,
      seed.use = SEED
    )
  )

  # 如果聚类使用 Harmony，另建一个 Harmony UMAP 显示 sketch 聚类结构。
  cluster_umap <- "umap"
  if (cluster_reduction != "pca") {
    obj <- run_step(
      "RunUMAP clustering view",
      Seurat::RunUMAP(
        obj,
        reduction = cluster_reduction,
        dims = dims,
        reduction.name = "umap.cluster",
        seed.use = SEED
      )
    )
    cluster_umap <- "umap.cluster"
  }

  obj <- run_step(
    "FindNeighbors",
    Seurat::FindNeighbors(
      obj,
      reduction = cluster_reduction,
      dims = dims,
      verbose = FALSE
    )
  )
  obj <- run_step(
    "FindClusters",
    Seurat::FindClusters(
      obj,
      resolution = CLUSTER$resolution,
      random.seed = SEED,
      verbose = FALSE
    )
  )

  # 将 sketch 学到的 cluster 标签、PCA 和 UMAP 投射到全部 8 µm bins。
  if (use_sketch) {
    obj <- run_step(
      "ProjectData",
      Seurat::ProjectData(
        obj,
        assay = full_assay,
        full.reduction = "full.pca",
        sketched.assay = cluster_assay,
        sketched.reduction = "pca",
        umap.model = "umap",
        dims = dims,
        refdata = list(cluster_full = "seurat_clusters")
      )
    )
    full_umap <- "full.umap"
  } else {
    obj$cluster_full <- as.character(Seurat::Idents(obj))
    full_umap <- cluster_umap
  }

  # 将关键参数写入对象，后续结果说明直接读取，避免手工抄错。
  obj@misc$spatial_pipeline <- list(
    platform = PLATFORM,
    bin_size = VISIUM_HD$bin_size,
    full_assay = full_assay,
    cluster_assay = cluster_assay,
    cluster_reduction = cluster_reduction,
    cluster_umap = cluster_umap,
    projection_reduction = "pca",
    full_umap = full_umap,
    dims = dims,
    resolution = CLUSTER$resolution,
    sketch = use_sketch
  )
  obj
}
