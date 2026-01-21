#' @name clusterCornerAxes
#' @author Junjun Lao
#' @title 在 UMAP/tSNE 聚类图上添加“角标坐标轴”（Corner Axes）
#'
#' @description
#' 用于绘制聚类散点图（按 cluster 上色），并在左下角加一组“角标坐标轴”，
#' 让无坐标轴的简洁图也能表达坐标方向（常用于论文图）。
#'
#' @param object Seurat 对象，必须包含对应降维（例如 `RunUMAP()` 产生的 "umap"）。
#' @param reduction 字符串；降维名称（umap/tsne 等），默认 "umap"。
#' @param groupFacet 字符串；当 `noSplit = FALSE` 时用于分面（facet）的 meta.data 列名。
#' @param clusterCol 字符串；meta.data 中用于上色/分组的列名，默认 "seurat_clusters"。
#' @param pSize 数值；点大小，默认 1。
#' @param aspect.ratio 数值；图形宽高比，默认 NULL。
#' @param noSplit 逻辑值；是否不分面（TRUE 不分面，FALSE 按 groupFacet 分面），默认 TRUE。
#' @param nrow 整数；分面行数（仅 noSplit=FALSE 时生效），默认 1。
#' @param relLength 数值；角标轴长度相对比例（0~1），默认 0.25。
#' @param relDist 数值；角标轴标签与轴线的相对距离，默认 0.1。
#' @param axes 字符串；"mul" 显示两条轴，"one" 只在某个分面上显示一组轴，默认 "mul"。
#' @param show.legend 是否显示图例，默认 TRUE。
#' @param legendPos 图例位置，默认 "right"。
#' @param keySize 图例点大小，默认 5。
#' @param cellLabel 是否在每个 cluster 的中心位置标注文字，默认 FALSE。
#' @param cellLabelSize 标注文字大小，默认 6。
#' @param cellLabelColor 标注文字颜色，默认 "black"。
#' @param lineTextcol 角标轴线与文字颜色，默认 "black"。
#' @param stripCol 分面条带背景色，默认 "white"。
#' @param arrowType 箭头类型（open/closed），默认 "closed"。
#' @param cornerTextSize 角标文字大小，默认 3。
#' @param base_size 主题基准字号，默认 14。
#' @param themebg 主题风格："default" 或 "bwCorner"，默认 "default"。
#' @param addCircle 是否给每个 cluster 加外轮廓（类似圈住一团点），默认 FALSE。
#' @param addCircle.legacy 是否使用 legacy 算法绘制轮廓（影响 nbin/nsm/addsm/sfac/qval 参数），默认 FALSE。
#' @param cicAlpha 轮廓填充透明度，默认 0.1。
#' @param cicDelta 轮廓扩张距离（仅 addCircle.legacy=FALSE 时生效），默认 0.1。
#' @param cicLineSize 轮廓线宽，默认 1。
#' @param cicLineColor 轮廓线颜色，默认 "grey50"。
#' @param cicLineLty 轮廓线型，默认 "dashed"。
#' @param nbin legacy：用于拟合 hull 的点数，默认 100。
#' @param nsm legacy：卷积点数（应小于 nbin），默认 10。
#' @param addsm legacy：额外卷积次数，默认 1。
#' @param qval legacy：扩张系数，默认 1。
#' @param sfac legacy：分位数阈值（<1），默认 1.5。
#'
#' @importFrom ggunchull stat_unchull0
#'
#' @return 返回一个 ggplot 对象。
#' @export
#' @examples
#' test <- system.file("extdata", "seuratTest.RDS", package = "scRNAtoolVis")
#'
#' tmp <- readRDS(test)
#'
#' # umap
#' clusterCornerAxes(
#'   object = tmp, reduction = "umap",
#'   noSplit = TRUE
#' )
#'
#' # arrowType
#' clusterCornerAxes(
#'   object = tmp, reduction = "umap",
#'   noSplit = TRUE, arrowType = "open"
#' )
#'
#' # facet by metadata column "orig.ident"
#' clusterCornerAxes(
#'   object = tmp, reduction = "umap",
#'   noSplit = FALSE, groupFacet = "orig.ident",
#'   relLength = 0.5
#' )
#'
#' # retain only one axes
#' clusterCornerAxes(
#'   object = tmp, reduction = "umap",
#'   noSplit = FALSE, groupFacet = "orig.ident",
#'   relLength = 0.5,
#'   axes = "one"
#' )
#'
#' # line color
#' clusterCornerAxes(
#'   object = tmp, reduction = "umap",
#'   noSplit = FALSE, groupFacet = "orig.ident",
#'   relLength = 0.5,
#'   lineTextcol = "grey50"
#' )
#'
#' # tsne
#' clusterCornerAxes(
#'   object = tmp, reduction = "tsne",
#'   noSplit = FALSE, groupFacet = "orig.ident",
#'   relLength = 0.5
#' )
#'
# define variables
globalVariables(c("x1", "y1", "linegrou", "angle", "lab", ".data"))

# define function
clusterCornerAxes <- function(
    object = NULL,
    reduction = "umap",
    groupFacet = NULL,
    clusterCol = "seurat_clusters",
    pSize = 1,
    aspect.ratio = NULL,
    noSplit = TRUE,
    nrow = 1,
    relLength = 0.25,
    relDist = 0.1,
    axes = "mul",
    show.legend = TRUE,
    legendPos = "right",
    keySize = 5,
    cellLabel = FALSE,
    cellLabelSize = 6,
    cellLabelColor = "black",
    lineTextcol = "black",
    stripCol = "white",
    arrowType = "closed",
    cornerTextSize = 3,
    base_size = 14,
    themebg = "default",
    addCircle = FALSE,
    addCircle.legacy = FALSE,
    cicDelta = 0.1,
    cicAlpha = 0.1,
    cicLineSize = 1,
    cicLineColor = "grey50",
    cicLineLty = "dashed",
    nbin = 100,
    nsm = 10,
    addsm = 1,
    qval = 1,
    sfac = 1.5) {
  # ============================================================================
  # 1) 参数检查（降维与 meta.data 列）
  # ============================================================================
  # 逐行注释：检查降维是否存在
  available_reductions <- Seurat::Reductions(object)
  if (length(available_reductions) == 0) {
    stop(
      "当前 Seurat 对象中没有任何降维结果（reductions）。\n",
      "请先运行例如 `RunUMAP()`/`RunTSNE()`，或把 `reduction` 改为对象中已有的 reduction。"
    )
  }
  if (!reduction %in% available_reductions) {
    stop(sprintf("未在对象中找到 reduction='%s'。可用 reductions: %s",
                 reduction, paste(available_reductions, collapse = ", ")))
  }

  # 逐行注释：当需要分面时，groupFacet 必须存在于 meta.data
  if (isFALSE(noSplit) && (is.null(groupFacet) || groupFacet == "")) {
    stop("当 `noSplit = FALSE` 时，必须提供 `groupFacet`（meta.data 列名）用于分面。")
  }

  # 逐行注释：clusterCol 必须存在
  if (!clusterCol %in% colnames(object@meta.data)) {
    stop(sprintf("meta.data 中不存在 clusterCol 指定的列：%s", clusterCol))
  }

  # ============================================================================
  # 2) 组织数据：降维坐标 + 元信息
  # ============================================================================
  # 逐行注释：提取降维坐标（两列，通常是 UMAP_1/UMAP_2 或 tSNE_1/tSNE_2）
  reduc <- data.frame(Seurat::Embeddings(object, reduction = reduction))

  # 逐行注释：提取元信息（每个细胞的 meta.data）
  meta <- object@meta.data

  # 逐行注释：合并成一个数据框，便于 ggplot2 直接使用
  pc12 <- cbind(reduc, meta)

  #######################################
  # text data
  # 逐行注释：计算每个 cluster 的“中心位置”（用中位数更稳健，避免离群点影响）
  namePos <- pc12 %>%
    dplyr::group_by(.data[[clusterCol]]) %>%
    dplyr::summarise(
      posMedia1 = stats::median(get(colnames(pc12)[1])),
      posMedia2 = stats::median(get(colnames(pc12)[2]))
    )

  #######################################

  # data range
  # 逐行注释：取坐标的最小值（向下取整），用于确定左下角角标轴的起点
  range <- floor(min(min(pc12[, 1]), min(pc12[, 2])))

  # get bottom-left coord
  # 逐行注释：角标轴起点（lower）向左下偏移一段距离 relDist
  lower <- range - relDist * abs(range)

  # label reldist to axes
  # 逐行注释：角标文字相对轴线的偏移距离
  labelRel <- relDist * abs(lower)

  # get relative line length
  # 逐行注释：角标轴终点（线段长度由 relLength 决定）
  linelen <- abs(relLength * lower) + lower

  # mid point
  # 逐行注释：用于放置角标轴文字的中点位置
  mid <- abs(relLength * lower) / 2 + lower

  # give reduction type
  # 逐行注释：根据 reduction 名称决定轴标签写 UMAP 还是 t-SNE
  if (startsWith(reduction, "umap")) {
    axs_label <- paste("UMAP", 2:1, sep = "")
  } else if (startsWith(reduction, "tsne")) {
    axs_label <- paste("t-SNE", 2:1, sep = "")
  } else {
    print("Please give correct type(umap or tsne)!")
  }

  if (axes == "mul") {
    # axies data
    # 逐行注释：两条轴线（竖线 + 横线），用 linegrou 区分分组
    axes <- data.frame(
      "x1" = c(lower, lower, lower, linelen),
      "y1" = c(lower, linelen, lower, lower),
      "linegrou" = c(1, 1, 2, 2)
    )
    # axies label
    # 逐行注释：两条轴的文字标签（一个竖着写，一个横着写）
    label <- data.frame(
      "lab" = c(axs_label),
      "angle" = c(90, 0),
      "x1" = c(lower - labelRel, mid),
      "y1" = c(mid, lower - labelRel)
    )
  } else if (axes == "one") {
    # 逐行注释：只保留一组角标轴（只在某个分面/组上显示）
    firstFacet <- unique(pc12[, groupFacet])[1]
    # axies data
    axes <- data.frame(
      "x1" = c(lower, lower, lower, linelen),
      "y1" = c(lower, linelen, lower, lower),
      "linegrou" = c(1, 1, 2, 2),
      "group" = rep(firstFacet, 2)
    )
    # axies label
    label <- data.frame(
      "lab" = c(axs_label),
      "angle" = c(90, 0),
      "x1" = c(lower - labelRel, mid),
      "y1" = c(mid, lower - labelRel),
      "group" = rep(firstFacet, 2)
    )

    # rename group name
    # 逐行注释：让列名与 facet 用的 groupFacet 一致，便于 ggplot2 在分面时匹配数据
    colnames(axes)[4] <- groupFacet
    colnames(label)[5] <- groupFacet
  } else {
    print("Please give correct args(mul or one)!")
  }

  ######################################################
  # plot
  p <- ggplot2::ggplot(
    pc12,
    ggplot2::aes_string(x = colnames(pc12)[1], y = colnames(pc12)[2])
  ) +
    ggplot2::geom_point(
      ggplot2::aes_string(color = clusterCol),
      size = pSize,
      show.legend = show.legend
    ) +
    ggplot2::theme_classic(base_size = base_size) +
    ggplot2::labs(x = "", y = "") +
    ggplot2::theme(
      strip.background = ggplot2::element_rect(colour = NA, fill = stripCol),
      aspect.ratio = aspect.ratio,
      legend.position = legendPos,
      plot.title = ggplot2::element_text(hjust = 0.5),
      axis.line = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank()
    ) +
    ggplot2::geom_line(
      data = axes,
      ggplot2::aes(x = x1, y = y1, group = linegrou),
      color = lineTextcol,
      arrow = ggplot2::arrow(
        length = ggplot2::unit(0.1, "inches"),
        ends = "last",
        type = arrowType
      )
    ) +
    ggplot2::geom_text(
      data = label,
      ggplot2::aes(x = x1, y = y1, angle = angle, label = lab),
      color = lineTextcol,
      fontface = "italic",
      size = cornerTextSize
    ) +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = keySize)))

  ######################################################
  # add text label
  if (cellLabel == FALSE) {
    plabel <- p
  } else {
    plabel <- p +
      ggrepel::geom_text_repel(
        data = namePos,
        ggplot2::aes_string(x = "posMedia1", y = "posMedia2", label = clusterCol),
        show.legend = FALSE,
        size = cellLabelSize,
        color = cellLabelColor
      )
  }

  ######################################################
  # add circle line
  if (addCircle == FALSE) {
    p0 <- plabel
  } else if (addCircle.legacy) {
    p0 <- plabel +
      ggunchull::stat_unchull0(
        ggplot2::aes_string(fill = clusterCol),
        alpha = cicAlpha,
        size = cicLineSize,
        color = cicLineColor,
        lty = cicLineLty,
        show.legend = FALSE,
        nbin = nbin,
        nsm = nsm,
        addsm = addsm,
        sfac = sfac,
        qval = qval
      )
  } else {
    p0 <- plabel +
      ggunchull::stat_unchull(
        ggplot2::aes_string(fill = clusterCol),
        alpha = cicAlpha,
        size = cicLineSize,
        color = cicLineColor,
        lty = cicLineLty,
        show.legend = FALSE,
        delta = cicDelta
      )
  }

  ######################################################
  # facet plot
  if (noSplit == TRUE) {
    p1 <- p0
  } else {
    p1 <- p0 + ggplot2::facet_wrap(facets = groupFacet, nrow = nrow)
  }

  ######################################################
  # theme style
  if (themebg == "bwCorner") {
    p2 <- p1 +
      ggplot2::theme_bw(base_size = base_size) +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        axis.text = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        aspect.ratio = 1,
        strip.background = ggplot2::element_rect(colour = NA, fill = stripCol)
      )
  } else if (themebg == "default") {
    p2 <- p1
  }

  return(p2)
}
