# 分析点 04：应用人工生物学注释并汇总每样本 bin 数

run_step04_annotation <- function() {
  obj <- read_object(OBJECT_CLUSTERED)
  obj <- apply_spatial_annotation(obj)
  save_object(obj, OBJECT_ANNOTATED)
  save_annotation_umap(obj)

  # 按样本、分组和注释汇总；这里只做描述，不把 bins 当作重复检验。
  meta <- obj[[]]
  summary <- dplyr::count(
    meta,
    dplyr::across(dplyr::all_of(c(
      VISIUM_HD$sample_col,
      VISIUM_HD$group_col,
      "cell_type"
    ))),
    name = "bins"
  )
  save_table(
    summary,
    file.path(ANNOTATION_DIR, "annotation_bin_counts.csv")
  )

  save_results_note(
    file.path(ANNOTATION_DIR, "README_results_CN.md"),
    c(
      "# Visium HD 人工注释结果说明",
      "",
      paste0("- 输入对象：", OBJECT_CLUSTERED, "。"),
      paste0("- 人工映射表：", ANNOTATION_PATH, "。"),
      "- cluster_celltype_mapping.csv 保存 cluster 到细胞类型的映射。",
      "- annotation_bin_counts.csv 保存各样本、分组和标签的 8 µm bin 数。",
      "- UMAP_cell_type.png/.pdf 使用完整投射 UMAP。",
      "- 8 µm bin 不等同于单细胞；标签应结合 marker 和组织位置复核。"
    )
  )
  invisible(obj)
}
