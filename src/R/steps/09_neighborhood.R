# 分析点 09：注释后 8 µm bin 的 kNN 邻域组成

run_step09_neighborhood <- function() {
  obj <- read_object(OBJECT_ANNOTATED)
  table <- run_spatial_neighborhood(obj)
  save_spatial_neighborhood_outputs(table)
  invisible(table)
}
