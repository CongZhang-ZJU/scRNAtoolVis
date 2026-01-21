#' @name cellRatioPlot
#' @author Junjun Lao
#' @title 细胞比例堆叠柱形图（可选流带连接）
#'
#' @description
#' 用于展示不同样本（或分组）中各细胞类型的**相对比例**（百分比）。
#' 默认会叠加 `ggalluvial::geom_flow()` 的流带效果，便于观察“同一细胞类型”在不同样本间的变化趋势。
#'
#' @param object Seurat 对象（需要在 `object@meta.data` 中包含分组列）。
#' @param sample.name 字符串；样本/分组列名（位于 `object@meta.data`）。
#' @param celltype.name 字符串；细胞类型列名（位于 `object@meta.data`）。
#' @param sample.order 字符向量；样本显示顺序（用于设置因子 levels），默认 NULL 表示按原始顺序。
#' @param col.width 数值；柱形宽度，默认 0.7。
#' @param flow.alpha 数值；流带透明度，默认 0.25。
#' @param flow.curve 数值；流带弯曲程度（knot.pos），默认 0。
#' @param fill.col 颜色向量；为每个细胞类型指定颜色，默认 NULL 表示自动生成。
#'
#' @return 返回一个 ggplot 对象。
#' @export
#'
#' @examples
#' \dontrun{
#' library(Seurat)
#' library(scRNAtoolVis)
#'
#' obj <- readRDS(system.file("extdata", "seuratTest.RDS", package = "scRNAtoolVis"))
#' # 构造示例分组列：样本（orig.ident）与细胞类型（idents）
#' obj$sample <- obj$orig.ident
#' obj$celltype <- as.character(Seurat::Idents(obj))
#'
#' cellRatioPlot(
#'   object = obj,
#'   sample.name = "sample",
#'   celltype.name = "celltype",
#'   sample.order = unique(obj$sample)
#' )
#' }

globalVariables(c("n", "num"))

cellRatioPlot <- function(
    object = NULL,
    sample.name = NULL,
    celltype.name = NULL,
    sample.order = NULL,
    col.width = 0.7,
    flow.alpha = 0.25,
    flow.curve = 0,
    fill.col = NULL) {
  # ============================================================================
  # 1) 参数检查（尽量在一开始给出明确的中文报错）
  # ============================================================================
  # 逐行注释：必须提供 sample.name 与 celltype.name，否则无法分组统计
  if (is.null(sample.name) || is.null(celltype.name)) {
    stop("请同时提供 `sample.name` 与 `celltype.name`（它们必须是 meta.data 中的列名）。")
  }

  # 逐行注释：从 Seurat 对象中取出元信息表（每行对应一个细胞）
  meta <- object@meta.data

  # 逐行注释：检查列名是否存在
  if (!sample.name %in% colnames(meta)) {
    stop(sprintf("`sample.name` 指定的列不存在：%s", sample.name))
  }
  if (!celltype.name %in% colnames(meta)) {
    stop(sprintf("`celltype.name` 指定的列不存在：%s", celltype.name))
  }

  # ============================================================================
  # 2) 设置样本顺序（可选）
  # ============================================================================
  # 逐行注释：
  # - 原实现写死了 meta$sample.name，容易与真实列名不一致导致无效
  # - 这里改成动态列名：meta[[sample.name]]
  if (!is.null(sample.order)) {
    meta[[sample.name]] <- factor(meta[[sample.name]], levels = sample.order)
  }

  # ============================================================================
  # 3) 统计每个 sample × celltype 的数量与相对比例
  # ============================================================================
  # 逐行注释：
  # - group_by：按样本与细胞类型分组
  # - summarise：统计该组细胞数 num
  # - mutate：将 num 转成该样本内的相对比例 rel_num（0~1）
  ratio.info <- meta %>%
    dplyr::group_by(.data[[sample.name]], .data[[celltype.name]]) %>%
    dplyr::summarise(num = n()) %>%
    dplyr::mutate(rel_num = num / sum(num))

  # ============================================================================
  # 4) 颜色准备
  # ============================================================================
  # 逐行注释：
  # - 若用户未指定 fill.col，则自动生成一组区分度较高的配色
  if (is.null(fill.col)) {
    fill.col <- jjAnno::useMyCol("paired", n = length(unique(meta[, celltype.name])))
  } else {
    fill.col <- fill.col
  }

  # ============================================================================
  # 5) 绘图
  # ============================================================================
  # 逐行注释：
  # - x：样本（sample.name）
  # - y：相对比例（rel_num）
  # - fill：细胞类型（celltype.name）
  p <-
    ggplot2::ggplot(
      ratio.info,
      ggplot2::aes_string(x = sample.name, y = "rel_num")
    ) +
    ggplot2::geom_col(
      ggplot2::aes_string(fill = celltype.name),
      width = col.width
    ) +
    ggalluvial::geom_flow(
      ggplot2::aes_string(
        stratum = celltype.name,
        alluvium = celltype.name,
        fill = celltype.name
      ),
      width = 0.5,
      alpha = flow.alpha,
      knot.pos = flow.curve
    ) +
    ggplot2::theme_bw() +
    ggplot2::coord_cartesian(expand = 0) +
    ggplot2::scale_y_continuous(labels = scales::label_percent()) +
    ggplot2::scale_fill_manual(
      values = fill.col,
      name = "Cell Type"
    ) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(size = ggplot2::rel(1.2), color = "black"),
      axis.title = ggplot2::element_text(size = ggplot2::rel(1.5), color = "black"),
      legend.text = ggplot2::element_text(size = ggplot2::rel(1.2), color = "black"),
      legend.title = ggplot2::element_text(size = ggplot2::rel(1.5), color = "black")
    ) +
    ggplot2::xlab("") +
    ggplot2::ylab("Cell percent ratio")

  return(p)
}
