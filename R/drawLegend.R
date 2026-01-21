#' @name drawLegend
#' @author Junjun Lao
#' @title 在图旁边补一个“细胞类型-簇编号”对照图例
#'
#' @description
#' 常见场景：你已经有一个 `DimPlot()`/自定义 ggplot，但是图例里看不到“cluster 编号 ↔ celltype 名称”的对应关系。
#' 本函数会从 `object@meta.data` 中提取两列（细胞类型列、簇编号列），生成一个小图例并与原图拼接。
#'
#' @param object Seurat 对象（用于提供 meta.data）。
#' @param plot ggplot 对象（主图）。
#' @param cellType 字符串；meta.data 中“细胞类型”的列名。
#' @param clusters 字符串；meta.data 中“簇编号/cluster”的列名。
#' @param ncol 整数；图例分成几列（常用 1 或 2），默认 1。
#' @param col 颜色向量；为每个 cellType 指定颜色。默认 NULL 表示使用 `scales::hue_pal()` 自动生成。
#' @param pt.size 数值；图例点大小，默认 8。
#' @param text.size 数值；图例文字大小，默认 4。
#'
#' @return 返回 `cowplot::plot_grid()` 拼接后的图。
#' @export
#'
#' @examples
#' \dontrun{
#' library(Seurat)
#' library(scRNAtoolVis)
#'
#' obj <- readRDS(system.file("extdata", "seuratTest.RDS", package = "scRNAtoolVis"))
#' obj$celltype <- as.character(Seurat::Idents(obj))
#'
#' p <- Seurat::DimPlot(obj, reduction = "umap", group.by = "celltype") + ggplot2::NoLegend()
#' drawLegend(
#'   object = obj,
#'   plot = p,
#'   cellType = "celltype",
#'   clusters = "seurat_clusters",
#'   ncol = 1
#' )
#' }

globalVariables(c("x", "y"))

drawLegend <- function(
    object = NULL,
    plot = NULL,
    cellType = NULL,
    clusters = NULL,
    ncol = 1,
    col = NULL,
    pt.size = 8,
    text.size = 4) {
  # ============================================================================
  # 1) 参数检查
  # ============================================================================
  if (is.null(object) || is.null(plot)) {
    stop("请同时提供 `object`（Seurat 对象）与 `plot`（ggplot 主图）。")
  }
  if (is.null(cellType) || is.null(clusters)) {
    stop("请提供 `cellType` 与 `clusters`（它们必须是 meta.data 中的列名）。")
  }

  # 逐行注释：取出元信息并检查列是否存在
  meta <- object@meta.data
  if (!cellType %in% colnames(meta)) stop(sprintf("meta.data 中不存在列：%s", cellType))
  if (!clusters %in% colnames(meta)) stop(sprintf("meta.data 中不存在列：%s", clusters))

  # ============================================================================
  # 2) 准备图例数据（去重后形成“cellType ↔ cluster”对照表）
  # ============================================================================
  leg.data <- object@meta.data %>%
    dplyr::select(.data[[cellType]], .data[[clusters]]) %>%
    unique()

  colnames(leg.data) <- c("cellType", "clusters")

  # ============================================================================
  # 3) 调整顺序（尽量按 cellType 因子顺序展示）
  # ============================================================================
  # 逐行注释：
  # - 如果 cellType 是 factor，则按 levels 的顺序重排
  # - 如果 cellType 是 character（levels 为 NULL），则保持当前顺序（避免重排后变成 0 行）
  if (!is.null(levels(leg.data$cellType))) {
    leg.data <- leg.data[match(levels(leg.data$cellType), leg.data$cellType), , drop = FALSE]
  }

  # ============================================================================
  # 4) 给图例中的点/文字分配网格坐标（支持多列）
  # ============================================================================
  if (ncol > 1) {
    leg.data$x <- rep(
      1:ncol,
      c(
        ceiling(nrow(leg.data) / ncol),
        nrow(leg.data) - ceiling(nrow(leg.data) / ncol)
      )
    )

    leg.data$y <- c(
      1:ceiling(nrow(leg.data) / ncol),
      1:(nrow(leg.data) - ceiling(nrow(leg.data) / ncol))
    )
  } else {
    leg.data$x <- 1
    leg.data$y <- seq_len(nrow(leg.data))
  }

  # order
  leg.data$cellType <- factor(leg.data$cellType, levels = rev(levels(leg.data$cellType)))

  # ============================================================================
  # 5) 绘制图例（一个单独的 ggplot）
  # ============================================================================
  if (is.null(col)) {
    color <- rev(scales::hue_pal()(nrow(leg.data)))
  } else {
    color <- rev(col)
  }

  pleg <- ggplot2::ggplot(leg.data, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_point(
      ggplot2::aes(color = cellType),
      show.legend = FALSE,
      size = pt.size
    ) +
    ggplot2::geom_text(ggplot2::aes(label = clusters)) +
    ggplot2::geom_text(
      ggplot2::aes(label = cellType),
      hjust = 0,
      nudge_x = 0.2,
      size = text.size
    ) +
    ggplot2::scale_color_manual(values = color) +
    ggplot2::scale_y_reverse() +
    ggplot2::xlim(0, ncol + 1) +
    ggplot2::theme_void()

  # ============================================================================
  # 6) 拼接主图与图例
  # ============================================================================
  cowplot::plot_grid(plot, pleg)
}
