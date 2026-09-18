# 分析点 05：生物学注释空间图、基因表达图和坐标表

run_step05_spatial_plots <- function() {
  obj <- read_object(OBJECT_ANNOTATED)
  save_spatial_annotation_plots(
    obj,
    group_by = "cell_type",
    output_dir = SPATIAL_DIR,
    stem_prefix = "spatial"
  )
  save_spatial_feature_plots(obj)
  write_spatial_results_note(obj)
  invisible(obj)
}
