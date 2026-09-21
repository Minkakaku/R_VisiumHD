# 通用运行与输出工具
# 所有正式分析点都通过这里统一写对象、表格、图和中文说明。

# 当左侧为空时使用默认值，常用于读取对象中的可选参数。
`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

HONGF_NPG <- c(
  "#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F",
  "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85"
)
HONGF_NPG_SEQUENTIAL <- c("#FFF5F0", "#F39B7F", "#E64B35", "#7A0019")
HONGF_NPG_DIVERGING <- c("#3C5488", "#F7F7F7", "#E64B35")

# 按固定 NPG 顺序生成分类配色；超过 10 类时从 NPG 锚点扩展。
npg_palette <- function(labels) {
  labels <- as.character(labels)
  colors <- HONGF_NPG
  if (length(labels) > length(colors)) {
    extra_n <- length(labels) - length(colors)
    extra <- grDevices::colorRampPalette(HONGF_NPG)(2L * extra_n + 1L)[seq.int(2L, 2L * extra_n, by = 2L)]
    colors <- c(colors, extra)
  }
  stats::setNames(colors[seq_along(labels)], labels)
}

options(
  ggplot2.discrete.colour = function(n) unname(npg_palette(seq_len(n))),
  ggplot2.discrete.fill = function(n) unname(npg_palette(seq_len(n)))
)

# 输出带时间戳的运行信息，便于在服务器日志中定位步骤。
log_msg <- function(...) {
  message("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", ...)
}

# 包装一个耗时步骤，并记录开始、结束和分钟数。
run_step <- function(label, expr) {
  start <- Sys.time()
  log_msg("START ", label)
  out <- force(expr)
  log_msg("DONE ", label, " (",
    round(as.numeric(difftime(Sys.time(), start, units = "mins")), 2),
    " min)")
  out
}

# 将 Seurat 或其他 R 对象保存到 data/processed/。
save_object <- function(obj, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(obj, path)
  log_msg("Saved ", path)
  invisible(path)
}

# 读取上一步保存的 RDS 对象。
read_object <- function(path) {
  readRDS(path)
}

# 将结果表保存为 UTF-8 CSV；配置模板使用 TSV。
save_table <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (tolower(tools::file_ext(path)) == "tsv") {
    write.table(x, path, sep = "\t", row.names = FALSE, quote = FALSE, fileEncoding = "UTF-8")
  } else {
    write.csv(x, path, row.names = FALSE, fileEncoding = "UTF-8")
  }
  invisible(path)
}

# 同时保存 600 dpi PNG 和矢量 PDF；正式图只调用这一入口。
save_plot <- function(p, path_no_ext, width = 8, height = 6, dpi = 600) {
  dir.create(dirname(path_no_ext), recursive = TRUE, showWarnings = FALSE)
  if (inherits(p, "patchwork")) {
    p <- p & ggplot2::theme(plot.subtitle = ggplot2::element_blank())
  } else if (inherits(p, "ggplot")) {
    p <- p + ggplot2::labs(subtitle = NULL) + ggplot2::theme(plot.subtitle = ggplot2::element_blank())
  }
  ggplot2::ggsave(
    paste0(path_no_ext, ".png"), p,
    width = width, height = height, dpi = dpi,
    limitsize = FALSE, bg = "white"
  )
  ggplot2::ggsave(
    paste0(path_no_ext, ".pdf"), p,
    width = width, height = height,
    limitsize = FALSE, bg = "white"
  )
  invisible(path_no_ext)
}

# 默认输出非栅格 UMAP，并以标签的稳定顺序应用 NPG 分类色。
npg_dimplot <- function(object, group.by = NULL, ..., raster = FALSE) {
  values <- if (is.null(group.by)) Seurat::Idents(object) else object[[group.by, drop = TRUE]]
  labels <- if (is.factor(values)) levels(droplevels(values)) else sort(unique(as.character(values)))
  Seurat::DimPlot(object, group.by = group.by, ..., raster = raster, cols = unname(npg_palette(labels)))
}

# 写入当前分析点的中文结果说明。
save_results_note <- function(path, lines) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(lines, path, useBytes = TRUE)
  invisible(path)
}

# 仅在执行可选分析点时检查其依赖包。
require_optional <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("Install optional package before this step: ", package, call. = FALSE)
  }
}
