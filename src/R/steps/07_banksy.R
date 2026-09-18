# 分析点 07：BANKSY 空间组织域（可选）

run_step07_banksy <- function() {
  obj <- read_object(OBJECT_ANNOTATED)
  obj <- run_banksy_domains(obj)
  save_banksy_outputs(obj)
  save_object(obj, OBJECT_BANKSY)
  invisible(obj)
}
