# 分析点 01：平台核对、加载和 QC
# 新项目的 workflow 脚本只需加载 config/src 后调用本函数。

run_step01_load_qc <- function() {
  sheet <- read_spatial_sample_sheet()

  # 先写平台审计表，再加载大表达矩阵。
  platform_table <- audit_sample_platforms(sheet)
  save_table(
    platform_table,
    file.path(QC_DIR, "platform_detection.csv")
  )
  if (any(platform_table$detected_platform != PLATFORM)) {
    stop(
      "All formal inputs must be detected as visium_hd. See ",
      file.path(QC_DIR, "platform_detection.csv"), ".",
      call. = FALSE
    )
  }

  # 逐样本过滤后再合并，避免把未过滤对象错误并入。
  objects <- load_visium_hd_samples(sheet)
  objects <- run_visium_hd_qc(objects, sheet)
  merged <- merge_spatial_samples(objects)
  save_object(merged, OBJECT_QC)

  save_results_note(
    file.path(QC_DIR, "README_results_CN.md"),
    c(
      "# Visium HD QC 结果说明",
      "",
      "- 输入：Space Ranger outs/ 或其上级样本目录。",
      paste0(
        "- 平台：Visium HD；bin size = ",
        VISIUM_HD$bin_size,
        " µm；assay = ", VISIUM_HD$assay, "。"
      ),
      "- platform_detection.csv 记录每个样本的平台判定。",
      "- qc_summary.csv 记录过滤前后 bin 数和主要 QC 指标。",
      "- qc_thresholds_used.csv 记录本次实际使用的逐样本阈值。",
      "- 每个 QC 图同时保存 300 dpi PNG 和矢量 PDF。"
    )
  )
  invisible(merged)
}
