#' @name jjVolcano
#' @author Junjun Lao
#' @title 按 cluster 展示 marker 的“火山棒棒糖图”（正负 avg_log2FC）
#'
#' @description
#' 输入通常来自 Seurat 的差异分析结果（例如 `FindAllMarkers()` 输出），
#' 本函数会按 cluster 分组，在 y 轴展示 `avg_log2FC`，并自动标注每个 cluster 的 top 上调/下调基因。
#'
#' @param diffData 数据框；差异分析结果。至少需要包含列：
#' \code{cluster}, \code{gene}, \code{avg_log2FC}, \code{p_val}, \code{p_val_adj}。
#' @param myMarkers 字符向量；若提供，则只标注这些基因（覆盖 topGeneN 逻辑），默认 NULL。
#' @param order.by 字符串；用于挑选 top 基因的排序列名，默认 "avg_log2FC"。
#' @param log2FC.cutoff 数值；|avg_log2FC| 过滤阈值，默认 0.25。
#' @param pvalue.cutoff 数值；p_val 过滤阈值，默认 0.05。
#' @param adjustP.cutoff 数值；用于“adjustP”配色类型的阈值，默认 0.01。
#' @param topGeneN 整数；每个 cluster 标注的 top 基因数量（正/负各 topGeneN），默认 5。
#' @param col.type 字符串；点的着色策略："updown" 或 "adjustP"，默认 "updown"。
#' @param back.col 字符串；背景条颜色，默认 "grey93"。
#' @param pSize 数值；点大小，默认 0.75。
#' @param aesCol 颜色向量；两种颜色（下调/上调），默认 c("#0099CC","#CC3333")。
#' @param legend.position 图例位置，默认 c(0.7,0.9)。
#' @param base_size 主题基准字号，默认 14。
#' @param tile.col 颜色向量；cluster 底部 tile 的配色，默认 `jjAnno::useMyCol("paired", n=9)`。
#' @param cluster.order 字符向量；指定 cluster 显示顺序，默认 NULL。
#' @param polar 逻辑值；是否转成极坐标展示，默认 FALSE。
#' @param expand 数值向量；极坐标时 y 方向扩展倍数，默认 c(-1,1)。
#' @param flip 逻辑值；是否翻转坐标，默认 FALSE。
#' @param celltypeSize 数值；cluster 名称文字大小，默认 3。
#' @param ... 其他参数会传给 `ggrepel::geom_text_repel()`（用于微调标注）。
#'
#' @return 返回一个 ggplot 对象。
#' @export
#'
#' @examples
#' \dontrun{
#' library(scRNAtoolVis)
#' data("pbmc.markers")
#' jjVolcano(diffData = pbmc.markers, topGeneN = 5)
#' }
globalVariables(c("p_val", "p_val_adj", "type", "type2"))
jjVolcano <- function(
    diffData = NULL,
    myMarkers = NULL,
    order.by = c("avg_log2FC"), # c("avg_log2FC","p_val")
    log2FC.cutoff = 0.25,
    pvalue.cutoff = 0.05,
    adjustP.cutoff = 0.01,
    topGeneN = 5,
    col.type = "updown",
    back.col = "grey93",
    pSize = 0.75,
    aesCol = c("#0099CC", "#CC3333"),
    legend.position = c(0.7, 0.9),
    base_size = 14,
    tile.col = jjAnno::useMyCol("paired", n = 9),
    cluster.order = NULL,
    polar = FALSE,
    expand = c(-1, 1),
    flip = FALSE,
    celltypeSize = 3,
    ...) {
  # ============================================================================
  # 1) 数据过滤
  # ============================================================================
  # 逐行注释：按 log2FC 与 p 值过滤，减少噪音点
  diff.marker <- diffData %>%
    dplyr::filter(abs(avg_log2FC) >= log2FC.cutoff & p_val < pvalue.cutoff)

  # ============================================================================
  # 2) 赋予类别标签（用于上色）
  # ============================================================================
  # 逐行注释：
  # - type：上调/下调（sigUp/sigDown）
  # - type2：adjusted p value 是否显著（用于另一种配色策略）
  diff.marker <- diff.marker %>%
    dplyr::mutate(type = ifelse(avg_log2FC >= log2FC.cutoff, "sigUp", "sigDown")) %>%
    dplyr::mutate(type2 = ifelse(p_val_adj < adjustP.cutoff,
      paste("adjust Pvalue < ", adjustP.cutoff, sep = ""),
      paste("adjust Pvalue >= ", adjustP.cutoff, sep = "")
    ))

  # ============================================================================
  # 3) cluster 顺序（可选）
  # ============================================================================
  if (!is.null(cluster.order)) {
    diff.marker$cluster <- factor(diff.marker$cluster,
      levels = cluster.order
    )
  }

  # ============================================================================
  # 4) 背景条数据：每个 cluster 的 min/max（用于画浅色背景）
  # ============================================================================
  purrr::map_df(unique(diff.marker$cluster), function(x) {
    tmp <- diff.marker %>%
      dplyr::filter(cluster == x)

    new.tmp <- data.frame(
      cluster = x,
      min = min(tmp$avg_log2FC) - 0.2,
      max = max(tmp$avg_log2FC) + 0.2
    )
    return(new.tmp)
  }) -> back.data

  # ============================================================================
  # 5) 选取每个 cluster 的 top 基因（正/负各 topGeneN）
  # ============================================================================
  top.marker.tmp <- diff.marker %>%
    dplyr::group_by(cluster)

  # order
  # if(length(order.by) == 1){
  #   top.marker.max <- top.marker.tmp %>%
  #     dplyr::slice_max(n = topGeneN,order_by = get(order.by))
  #
  #   top.marker.min <- top.marker.tmp %>%
  #     dplyr::group_by(cluster) %>%
  #     dplyr::slice_min(n = topGeneN,order_by = get(order.by))
  #
  # }else{
  #   top.marker.max <- top.marker.tmp %>%
  #     dplyr::arrange(dplyr::desc(get(order.by[1])),get(order.by[2])) %>%
  #     dplyr::slice_head(n = topGeneN)
  #
  #   top.marker.min <- top.marker.tmp %>%
  #     dplyr::arrange(dplyr::desc(get(order.by[1])),get(order.by[2])) %>%
  #     dplyr::slice_tail(n = topGeneN)
  # }

  top.marker.max <- top.marker.tmp %>%
    dplyr::slice_max(n = topGeneN, order_by = get(order.by))

  top.marker.min <- top.marker.tmp %>%
    dplyr::slice_min(n = topGeneN, order_by = get(order.by))

  # combine
  top.marker <- rbind(top.marker.max, top.marker.min)

  # ============================================================================
  # 6) 是否使用用户自定义要标注的基因
  # ============================================================================
  if (!is.null(myMarkers)) {
    top.marker <- diff.marker %>%
      dplyr::filter(gene %in% myMarkers)
  } else {
    top.marker <- top.marker
  }

  # ====================================================================
  # plot
  p1 <- ggplot2::ggplot(
    diff.marker,
    ggplot2::aes(x = cluster, y = avg_log2FC)
  ) +
    # add back cols
    ggplot2::geom_col(
      data = back.data,
      ggplot2::aes(x = cluster, y = min), fill = back.col
    ) +
    ggplot2::geom_col(
      data = back.data,
      ggplot2::aes(x = cluster, y = max), fill = back.col
    )

  # ap1 <- paste("adjust Pvalue >= ", adjustP.cutoff, sep = '')
  # ap2 <- paste("adjust Pvalue < ", adjustP.cutoff, sep = '')

  # color type
  if (col.type == "updown") {
    p2 <- p1 +
      # add point
      ggplot2::geom_jitter(ggplot2::aes(color = type), size = pSize) +
      ggplot2::scale_color_manual(values = c("sigDown" = aesCol[1], "sigUp" = aesCol[2]))
  } else if (col.type == "adjustP") {
    p2 <- p1 +
      # add point
      ggplot2::geom_jitter(ggplot2::aes(color = type2), size = pSize) +
      ggplot2::scale_color_manual(values = c(aesCol[2], aesCol[1]))
  }

  # theme details
  p3 <- p2 +
    ggplot2::scale_y_continuous(n.breaks = 6) +
    ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      legend.position = legend.position,
      legend.title = ggplot2::element_blank(),
      legend.background = ggplot2::element_blank()
    ) +
    ggplot2::xlab("Clusters") + ggplot2::ylab("Average log2FoldChange") +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 5)))

  # add tile
  p4 <- p3 +
    ggplot2::geom_tile(ggplot2::aes(x = cluster, y = 0, fill = cluster),
      color = "black",
      height = log2FC.cutoff * 2,
      alpha = 0.3,
      show.legend = FALSE
    ) +
    ggplot2::scale_fill_manual(values = tile.col) +
    # add gene label
    ggrepel::geom_text_repel(
      data = top.marker,
      ggplot2::aes(x = cluster, y = avg_log2FC, label = gene),
      max.overlaps = 50,
      ...
    )

  # whether coord_plolar
  if (polar == TRUE) {
    p5 <- p4 +
      geomtextpath::geom_textpath(ggplot2::aes(x = cluster, y = 0, label = cluster)) +
      ggplot2::scale_y_continuous(
        n.breaks = 6,
        expand = ggplot2::expansion(mult = expand)
      ) +
      ggplot2::theme_void(base_size = base_size) +
      ggplot2::theme(
        legend.position = legend.position,
        legend.title = ggplot2::element_blank()
      ) +
      ggplot2::coord_polar(clip = "off", theta = "x")
  } else {
    # whether flip plot
    if (flip == TRUE) {
      p5 <- p4 +
        ggplot2::scale_y_continuous(n.breaks = 6) +
        ggplot2::geom_label(ggplot2::aes(x = cluster, y = 0, label = cluster),size = celltypeSize) +
        ggplot2::theme(
          axis.line.y = ggplot2::element_blank(),
          axis.text.y = ggplot2::element_blank(),
          axis.ticks.y = ggplot2::element_blank()
        ) +
        ggplot2::coord_flip()
    } else {
      p5 <- p4 +
        ggplot2::scale_y_continuous(n.breaks = 6) +
        ggplot2::geom_text(ggplot2::aes(x = cluster, y = 0, label = cluster),size = celltypeSize) +
        ggplot2::theme(
          axis.line.x = ggplot2::element_blank(),
          axis.text.x = ggplot2::element_blank(),
          axis.ticks.x = ggplot2::element_blank()
        )
    }
  }
  return(p5)
}

###############################
#' This is a test data for this package
#' test data description
#'
#' @name pbmc.markers
#' @docType data
#' @author Junjun Lao
"pbmc.markers"
