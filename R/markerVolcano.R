#' @name markerVolcano
#' @author Junjun Lao
#' @title Marker 基因“火山分面图”（按 cluster 分面）
#'
#' @description
#' 输入通常为 Seurat `FindAllMarkers()` 的结果（也可读取包内示例 csv），
#' x 轴为 \(\Delta\%\)（pct.1 - pct.2），y 轴为 avg_log2FC，并按 cluster 分面。
#' 默认会为每个 cluster 标注 top 上调/下调基因。
#'
#' @param markers 数据框；marker 结果。至少需要列：\code{cluster}, \code{gene}, \code{avg_log2FC}, \code{pct.1}, \code{pct.2}。
#' @param ownGene 字符向量；若提供，则只标注这些基因（覆盖 topn 逻辑），默认 NULL。
#' @param topn 整数；每个 cluster 标注的 top 基因数量（正/负各 topn），默认 5。
#' @param log2FC 数值；横线阈值（±log2FC），默认 0.25。
#' @param labelCol 颜色向量；不同 cluster 的标签颜色，默认 NULL（建议传 `ggsci::pal_npg()(9)` 等）。
#' @param hlineSize 数值；横线粗细，默认 1。
#' @param hlineColor 字符串；横线颜色，默认 "grey50"。
#' @param pforce 数值；正向基因标签 repel 力度，默认 5。
#' @param nforce 数值；负向基因标签 repel 力度，默认 2.5。
#' @param nudge_x 数值；标签水平偏移基准，默认 0.8。
#' @param pnudge_y 数值；正向标签 y 偏移，默认 0.25。
#' @param nnudge_y 数值；负向标签 y 偏移，默认 0。
#' @param base_size 主题基准字号，默认 14。
#' @param facetColor 分面边框颜色，默认 NA。
#' @param facetFill 分面背景填充色，默认 "white"。
#' @param ylab y 轴标题，默认 "Log2-Fold Change"。
#' @param nrow 分面行数，默认 1。
#'
#' @return 返回一个 ggplot 对象。
#' @export
#'
#' @examples
#' \dontrun{
#' library(scRNAtoolVis)
#' test <- system.file("extdata", "pbmc.markers.csv", package = "scRNAtoolVis")
#' markers <- read.csv(test)
#'
#' markerVolcano(
#'   markers = markers,
#'   topn = 5,
#'   labelCol = ggsci::pal_npg()(9)
#' )
#' }
#'
#' @examples
#' test <- system.file("extdata", "pbmc.markers.csv", package = "scRNAtoolVis")
#' markers <- read.csv(test)
#'
#' markerVolcano(
#'   markers = markers,
#'   topn = 5,
#'   labelCol = ggsci::pal_npg()(9)
#' )
# define variables
globalVariables(c("avg_log2FC", "cluster", "gene", "pct.1", "pct.2"))

# define function
markerVolcano <- function(
    markers = NULL,
    ownGene = NULL,
    topn = 5,
    log2FC = 0.25,
    labelCol = NULL,
    hlineSize = 1,
    hlineColor = "grey50",
    pforce = 5,
    nforce = 2.5,
    nudge_x = 0.8,
    pnudge_y = 0.25,
    nnudge_y = 0,
    base_size = 14,
    facetColor = NA,
    facetFill = "white",
    ylab = "Log2-Fold Change",
    nrow = 1) {
  # ============================================================================
  # 1) 选取要标注的基因
  # ============================================================================
  # 逐行注释：
  # - 若 ownGene=NULL：每个 cluster 选 top 上调与 top 下调
  # - 若 ownGene 非空：只标注指定基因
  if (is.null(ownGene)) {
    # top genes
    toppos <- markers %>%
      dplyr::group_by(cluster) %>%
      dplyr::top_n(n = topn, wt = avg_log2FC)
    topneg <- markers %>%
      dplyr::group_by(cluster) %>%
      dplyr::top_n(n = -topn, wt = avg_log2FC)

    # merge
    topgene <- rbind(toppos, topneg)
  } else {
    topgene <- markers %>% dplyr::filter(gene %in% ownGene)
    toppos <- topgene %>% dplyr::filter(avg_log2FC > 0)
    topneg <- topgene %>% dplyr::filter(avg_log2FC < 0)
  }

  # ============================================================================
  # 2) 绘图（x = pct.1 - pct.2；y = avg_log2FC；按 cluster 分面）
  # ============================================================================
  ggplot2::ggplot(
    markers,
    ggplot2::aes(x = pct.1 - pct.2, y = avg_log2FC)
  ) +
    ggplot2::geom_point(color = "grey80") +
    ggplot2::geom_hline(
      yintercept = c(-log2FC, log2FC),
      lty = "dashed",
      size = hlineSize,
      color = hlineColor
    ) +
    ggrepel::geom_text_repel(
      data = toppos,
      ggplot2::aes(
        x = pct.1 - pct.2,
        y = avg_log2FC,
        label = gene,
        color = cluster
      ),
      show.legend = FALSE,
      direction = "y",
      hjust = 1,
      nudge_y = pnudge_y,
      force = pforce,
      nudge_x = -nudge_x - (toppos$pct.1 - toppos$pct.2)
    ) +
    ggrepel::geom_text_repel(
      data = topneg,
      ggplot2::aes(
        x = pct.1 - pct.2,
        y = avg_log2FC,
        label = gene,
        color = cluster
      ),
      show.legend = FALSE,
      direction = "y",
      hjust = 0,
      nudge_y = nnudge_y,
      force = nforce,
      nudge_x = nudge_x - (topneg$pct.1 - topneg$pct.2)
    ) +
    ggplot2::geom_point(
      data = topgene,
      show.legend = FALSE,
      ggplot2::aes(
        x = pct.1 - pct.2,
        y = avg_log2FC,
        color = cluster
      )
    ) +
    ggplot2::scale_color_manual(name = "", values = labelCol) +
    # x y breaks label
    # scale_y_continuous(limits = c(-6,10),breaks = seq(-6,10,2)) +
    # scale_x_continuous(limits = c(-1,1),breaks = seq(-1,1,0.5)) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      strip.background = ggplot2::element_rect(color = facetColor, fill = facetFill)
    ) +
    ggplot2::xlab(expression(Delta ~ "Percentage Difference")) +
    ggplot2::ylab(ylab) +
    ggplot2::facet_wrap(~cluster, nrow = nrow, scales = "fixed")
}
