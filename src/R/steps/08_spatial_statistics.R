# 分析点 08：Moran's I 空间变异基因

run_step08_spatial_statistics <- function() {
  obj <- read_object(OBJECT_ANNOTATED)
  table <- run_spatial_statistics(obj)
  save_spatial_statistics_outputs(table)
  invisible(table)
}
