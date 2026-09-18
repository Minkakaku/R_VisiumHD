# BANKSY 空间组织域
# BANKSY 把邻域平均表达和表达梯度加入聚类，适合识别 Visium HD 的连续组织域。

# 从每个 Seurat image 提取坐标，写入 BANKSY 需要的 metadata 列。
add_banksy_coordinates <- function(obj) {
  coords <- lapply(Seurat::Images(obj), function(image) {
    out <- spatial_coordinates(obj, image)
    rownames(out) <- out$cell_id
    out <- out[, 1:2, drop = FALSE]
    colnames(out) <- c("banksy_x", "banksy_y")
    out
  })
  coords <- do.call(rbind, coords)
  obj$banksy_x <- coords[colnames(obj), "banksy_x"]
  obj$banksy_y <- coords[colnames(obj), "banksy_y"]
  obj
}

# 按 Seurat Visium HD 官方示例运行 BANKSY、PCA 和组织域聚类。
run_banksy_domains <- function(obj) {
  require_optional("Banksy")
  require_optional("SeuratWrappers")

  full_assay <- VISIUM_HD$assay
  obj <- add_banksy_coordinates(obj)
  obj <- SeuratObject::JoinLayers(obj, assay = full_assay)
  Seurat::DefaultAssay(obj) <- full_assay

  obj <- run_step(
    "RunBanksy",
    SeuratWrappers::RunBanksy(
      obj,
      lambda = BANKSY$lambda,
      verbose = TRUE,
      assay = full_assay,
      slot = "data",
      features = "variable",
      k_geom = BANKSY$k_geom,
      dimx = "banksy_x",
      dimy = "banksy_y",
      group = VISIUM_HD$sample_col
    )
  )

  Seurat::DefaultAssay(obj) <- "BANKSY"
  obj <- run_step(
    "RunPCA BANKSY",
    Seurat::RunPCA(
      obj,
      assay = "BANKSY",
      reduction.name = "pca.banksy",
      features = rownames(obj[["BANKSY"]]),
      npcs = BANKSY$npcs,
      verbose = FALSE
    )
  )
  obj <- run_step(
    "FindNeighbors BANKSY",
    Seurat::FindNeighbors(
      obj,
      reduction = "pca.banksy",
      dims = BANKSY$dims,
      verbose = FALSE
    )
  )
  obj <- run_step(
    "FindClusters BANKSY",
    Seurat::FindClusters(
      obj,
      cluster.name = BANKSY$cluster_col,
      resolution = BANKSY$resolution,
      random.seed = SEED,
      verbose = FALSE
    )
  )

  Seurat::DefaultAssay(obj) <- full_assay
  obj@misc$banksy <- BANKSY
  obj
}

# 保存 BANKSY 空间域图、每样本域数量表和中文说明。
save_banksy_outputs <- function(obj) {
  save_spatial_annotation_plots(
    obj,
    group_by = BANKSY$cluster_col,
    output_dir = BANKSY_DIR,
    stem_prefix = "spatial_domain"
  )

  meta <- obj[[]]
  summary <- dplyr::count(
    meta,
    dplyr::across(dplyr::all_of(c(
      VISIUM_HD$sample_col,
      VISIUM_HD$group_col,
      BANKSY$cluster_col
    ))),
    name = "bins"
  )
  save_table(
    summary,
    file.path(BANKSY_DIR, "banksy_domain_bin_counts.csv")
  )

  save_results_note(
    file.path(BANKSY_DIR, "README_results_CN.md"),
    c(
      "# BANKSY 空间组织域结果说明",
      "",
      paste0("- 输入 assay：", VISIUM_HD$assay, "。"),
      paste0(
        "- lambda = ", BANKSY$lambda,
        "；k_geom = ", BANKSY$k_geom,
        "；resolution = ", BANKSY$resolution, "。"
      ),
      "- lambda 控制邻域信息权重；k_geom 越大，得到的组织域通常越大。",
      "- banksy_domain_bin_counts.csv 保存每个样本的域组成。",
      "- BANKSY domain 是空间组织域，不应直接当作单细胞类型。"
    )
  )
}
