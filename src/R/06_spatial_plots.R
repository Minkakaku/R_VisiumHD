# 空间坐标、注释和基因表达图
# 图形遵循统一配色，并与精确坐标/metadata 表放在同一个结果目录。

# 统一使用 NPG 浅红到深酒红连续表达梯度。
feature_cols <- HONGF_NPG_SEQUENTIAL

# 统一提取空间坐标和 cell_id，兼容 VisiumV2 行名与 FOV 的 cell 列。
spatial_coordinates <- function(obj, image) {
  coords <- Seurat::GetTissueCoordinates(obj[[image]])
  if ("cell" %in% colnames(coords)) {
    cell_id <- as.character(coords$cell)
    coords$cell <- NULL
  } else {
    cell_id <- rownames(coords)
  }
  coords$cell_id <- cell_id
  coords
}

# 合并一个图像的组织坐标和 Seurat metadata，作为空间图支撑表。
spatial_image_table <- function(obj, image) {
  coords <- spatial_coordinates(obj, image)
  meta <- obj[[]]
  meta$cell_id <- rownames(meta)
  dplyr::left_join(coords, meta, by = "cell_id")
}

# 按图像保存注释、RCTD 或 BANKSY 等离散标签的空间分布。
save_spatial_annotation_plots <- function(
  obj,
  group_by = "cell_type",
  output_dir = SPATIAL_DIR,
  stem_prefix = "spatial"
) {
  for (image in Seurat::Images(obj)) {
    labels <- obj[[group_by, drop = TRUE]]
    labels <- if (is.factor(labels)) levels(droplevels(labels)) else sort(unique(as.character(labels)))
    p <- Seurat::SpatialDimPlot(
      obj,
      images = image,
      group.by = group_by,
      pt.size.factor = SPATIAL_PLOT$point_size_factor,
      cols = unname(npg_palette(labels))
    ) +
      ggplot2::ggtitle(sub("\\.[0-9]{3}um$", "", image))

    save_plot(
      p,
      file.path(
        output_dir,
        paste0(stem_prefix, "_", group_by, "_", image)
      ),
      width = SPATIAL_PLOT$width,
      height = SPATIAL_PLOT$height
    )
    save_table(
      spatial_image_table(obj, image),
      file.path(
        output_dir,
        paste0(stem_prefix, "_metadata_", image, ".csv")
      )
    )
  }
}

# 按指定图像保存连续基因表达图；q05/q95 降低极端值对色阶的影响。
save_spatial_feature_plots <- function(
  obj,
  features = SPATIAL_FEATURES,
  images = Seurat::Images(obj),
  output_dir = SPATIAL_DIR,
  stem_prefix = "spatial_features"
) {
  Seurat::DefaultAssay(obj) <- VISIUM_HD$assay
  features <- intersect(features, rownames(obj))

  for (image in images) {
    p <- Seurat::SpatialFeaturePlot(
      obj,
      images = image,
      features = features,
      pt.size.factor = SPATIAL_PLOT$point_size_factor,
      keep.scale = "all",
      min.cutoff = "q05",
      max.cutoff = "q95",
      ncol = 2
    ) +
      patchwork::plot_annotation(
        title = sub("\\.[0-9]{3}um$", "", image)
      ) &
      ggplot2::scale_fill_gradientn(colours = feature_cols)

    save_plot(
      p,
      file.path(output_dir, paste0(stem_prefix, "_", image)),
      width = 12,
      height = max(7, ceiling(length(features) / 2) * 5)
    )
  }
}

# 记录空间图所用对象、降维、PC 和解释边界。
write_spatial_results_note <- function(obj) {
  settings <- obj@misc$spatial_pipeline
  save_results_note(
    file.path(SPATIAL_DIR, "README_results_CN.md"),
    c(
      "# Visium HD 空间图结果说明",
      "",
      paste0(
        "- 平台：Visium HD；Space Ranger bin size = ",
        settings$bin_size, " µm。"
      ),
      paste0("- 表达 assay：", settings$full_assay, "。"),
      paste0(
        "- 聚类 assay：", settings$cluster_assay,
        "；聚类降维：", settings$cluster_reduction, "。"
      ),
      paste0(
        "- 使用 PC：", paste(settings$dims, collapse = ", "),
        "；resolution = ", settings$resolution, "。"
      ),
      "- 每张正式图均成对保存为 600 dpi PNG 和矢量 PDF。",
      "- spatial_metadata_*.csv 保存坐标、样本、cluster 和细胞类型。",
      "- Visium HD 的观测单位为 8 µm bins，不等同于 Xenium 单细胞边界。"
    )
  )
}
