# 分析点 03：投射 cluster 的完整 Spatial assay marker

run_step03_markers <- function() {
  obj <- read_object(OBJECT_CLUSTERED)

  # 与 Seurat Visium HD 官方教程一致：完整 bins 分层抽样后找 marker。
  marker_obj <- prepare_visium_hd_marker_object(obj)
  markers <- run_step(
    "FindAllMarkers full projected bins",
    find_visium_hd_markers(marker_obj)
  )
  save_table(markers, file.path(MARKER_DIR, "all_markers.csv"))
  save_marker_dotplot(marker_obj, markers)

  # 首次运行 marker 后生成待填写的人工注释表。
  if (!file.exists(ANNOTATION_PATH)) {
    write_annotation_template(obj, ANNOTATION_PATH)
  }

  save_results_note(
    file.path(MARKER_DIR, "README_results_CN.md"),
    c(
      "# Visium HD marker 结果说明",
      "",
      paste0("- marker assay：", MARKERS$assay, "。"),
      paste0(
        "- cluster 列：", MARKERS$cluster_col,
        "；每群最多抽样 ",
        MARKERS$downsample_per_cluster, " 个 bins。"
      ),
      paste0(
        "- min.pct = ", MARKERS$min_pct,
        "；logFC threshold = ",
        MARKERS$logfc_threshold, "。"
      ),
      "- all_markers.csv 保存完整 marker 表。",
      "- top_markers.csv 和 top_marker_dotplot_data.csv 保存作图依据。",
      paste0(
        "- 根据 marker 编辑 ",
        ANNOTATION_PATH,
        " 后再执行注释步骤。"
      )
    )
  )
  invisible(markers)
}
