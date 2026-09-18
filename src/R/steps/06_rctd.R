# 分析点 06：使用单细胞参考进行 RCTD 解卷积（可选）

run_step06_rctd <- function() {
  obj <- read_object(OBJECT_ANNOTATED)
  reference <- read_object(REFERENCE_RDS)
  output <- run_step(
    "RCTD",
    run_rctd_visium_hd(obj, reference)
  )
  save_rctd_outputs(output)
  save_object(output$object, OBJECT_RCTD)
  invisible(output$object)
}
