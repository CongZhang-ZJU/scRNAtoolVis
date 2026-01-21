#' @name featureCornerAxes
#' @author Junjun Lao
#' @title 给 FeaturePlot 加“角标坐标轴”（Corner Axes）
#'
#' @description
#' 基于 Seurat 的降维坐标（UMAP/tSNE）与基因表达，绘制类似 FeaturePlot 的散点图，
#' 并在左下角添加“角标坐标轴”，便于在去掉坐标轴时仍能表达坐标方向。
#'
#' @param object Seurat 对象。
#' @param reduction 字符串；降维名称（umap/tsne），默认 "umap"。
#' @param features 字符向量；要绘制的基因列表。
#' @param slot 读取表达矩阵的层（兼容 Seurat v4/v5），默认 "data"。
#' - Seurat v4：等价于 `FetchData(..., slot = slot)`
#' - Seurat v5：等价于 `FetchData(..., layer = slot)`
#' @param groupFacet 字符串；用于分面的 meta.data 列名。若为 NULL，则只按基因分面。
#' @param relLength "num", the corner axis line relative length to plot axis(0-1).
#' @param relDist "num", the relative distance of corner axis label to axis.
#' @param aspect.ratio "num", plot width and height ratio, default NULL.
#' @param low "string", point color with low expression.
#' @param high "string", point color with high expression.
#' @param axes "string", show multiple corner axis or only one (mul/one), default "mul".
#' @param legendPos "string", legend position same as ggplot theme function, default "right".
#' @param stripCol "string", facet background color, defaults "white".
#' @param pSize "num", point size.
#' @param arrowType "string", arrow type (open/closed), default "closed".
#' @param lineTextcol "string", facet background color, default "white".
#' @param cornerTextSize "num", the corner label text size, default is 5.
#' @param base_size "num", theme base size, default is 14.
#' @param themebg Another theme style, default is "default", or "bwCorner".
#' @param show.legend Whether show legend, default "TRUE".
#' @param cornerVariable Which group  corner axis to be added when "axes" set to "one", default is the first group.
#' @param nLayout = NULL Similar to the ncol/nrow for the layout, default is the gene numbers.
#' @param minExp Minimum expression value defined, default is NULL.
#' @param maxExp Maxmum expression value defined, default is NULL.
#' @return Return a ggplot.
#' @export
#' @examples
#'
#' test <- system.file("extdata", "seuratTest.RDS", package = "scRNAtoolVis")
#'
#' tmp <- readRDS(test)
#'
#' # umap
#' featureCornerAxes(
#'   object = tmp,
#'   reduction = "umap",
#'   groupFacet = "orig.ident",
#'   relLength = 0.5,
#'   relDist = 0.2,
#'   features = c("Actb", "Ythdc1", "Ythdf2"),
#'   slot = "data"
#' )
#'
#' # one axes
#' featureCornerAxes(
#'   object = tmp, reduction = "umap",
#'   groupFacet = "orig.ident",
#'   features = c("Actb", "Ythdc1", "Ythdf2"),
#'   relLength = 0.5, relDist = 0.2,
#'   axes = "one",
#'   lineTextcol = "grey50"
#' )
#'
#' # tsne
#' featureCornerAxes(
#'   object = tmp, reduction = "tsne",
#'   groupFacet = "orig.ident",
#'   relLength = 0.5, relDist = 0.2,
#'   features = c("Actb", "Ythdc1", "Ythdf2")
#' )
#'
#'
# define variables
globalVariables(c("x1", "y1", "linegrou", "angle", "lab", "gene_name", "scaledValue"))

# define function
featureCornerAxes <- function(
    object = NULL,
    reduction = "umap",
    features = NULL,
    slot = "data",
    groupFacet = "orig.ident",
    minExp = NULL,
    maxExp = NULL,
    relLength = 0.25,
    relDist = 0.1,
    aspect.ratio = NULL,
    low = "lightgrey",
    high = "red",
    axes = "mul",
    show.legend = TRUE,
    legendPos = "right",
    stripCol = "white",
    cornerVariable = NULL,
    nLayout = NULL,
    pSize = 1,
    arrowType = "closed",
    lineTextcol = "black",
    cornerTextSize = 3,
    base_size = 14,
    themebg = "default") {
  # 逐行注释：
  # - 与 featurePlot 同理：先检查 reduction 是否存在，给出中文错误信息
  available_reductions <- Seurat::Reductions(object)
  if (length(available_reductions) == 0) {
    stop(
      "当前 Seurat 对象中没有任何降维结果（reductions）。\n",
      "请先运行例如 `RunPCA()` / `RunUMAP()` / `RunTSNE()`，或把 `reduction` 改为对象中已有的 reduction。"
    )
  }
  if (!reduction %in% available_reductions) {
    stop(
      sprintf("未在对象中找到 reduction='%s'。可用 reductions: %s",
              reduction, paste(available_reductions, collapse = ", "))
    )
  }

  # make PC data
  reduc <- data.frame(Seurat::Embeddings(object, reduction = reduction))

  # metadata
  meta <- object@meta.data

  # combine
  pc12 <- cbind(reduc, meta)

  # get gene expression
  # 逐行注释：
  # - 同 featurePlot：统一使用兼容封装，确保 Seurat v5 不触发 slot defunct
  geneExp <- .scRNAtoolVis_fetch_data(
    object = object,
    vars = features,
    slot = slot
  )

  # cbind
  mer <- cbind(pc12, geneExp)

  # merge data
  megredf <- reshape2::melt(
    mer,
    id.vars = colnames(pc12),
    variable.name = "gene_name",
    value.name = "scaledValue"
  )

  # data range
  range <- floor(min(min(pc12[, 1]), min(pc12[, 2])))

  # get bottom-left coord
  lower <- range - relDist * abs(range)

  # label reldist to axes
  labelRel <- relDist * abs(lower)

  # get relative line length
  linelen <- abs(relLength * lower) + lower

  # mid point
  mid <- abs(relLength * lower) / 2 + lower

  # give reduction type
  if (startsWith(reduction, "umap")) {
    axs_label <- paste("UMAP", 2:1, sep = "")
  } else if (startsWith(reduction, "tsne")) {
    axs_label <- paste("t-SNE", 2:1, sep = "")
  } else {
    print("Please give correct type(umap or tsne)!")
  }

  if (axes == "mul") {
    # axises data
    axes <- data.frame(
      "x1" = c(lower, lower, lower, linelen),
      "y1" = c(lower, linelen, lower, lower),
      "linegrou" = c(1, 1, 2, 2)
    )
    # axises label
    label <- data.frame(
      "lab" = c(axs_label),
      "angle" = c(90, 0),
      "x1" = c(lower - labelRel, mid),
      "y1" = c(mid, lower - labelRel)
    )
  } else if (axes == "one") {
    # add specific group corner
    if (is.null(cornerVariable)) {
      lev <- levels(pc12[, groupFacet])
      if (!is.null(lev)) {
        firstFacet <- factor(lev[1], levels = lev)
      } else {
        firstFacet <- unique(pc12[, groupFacet])[1]
      }
    } else {
      lev <- levels(pc12[, groupFacet])
      if (!is.null(lev)) {
        firstFacet <- factor(cornerVariable, levels = lev)
      } else {
        firstFacet <- cornerVariable
      }
    }

    # axises data
    axes <- data.frame(
      "x1" = c(lower, lower, lower, linelen),
      "y1" = c(lower, linelen, lower, lower),
      "linegrou" = c(1, 1, 2, 2),
      "group" = rep(firstFacet, 2)
    )
    # axises label
    label <- data.frame(
      "lab" = c(axs_label),
      angle = c(90, 0),
      "x1" = c(lower - labelRel, mid),
      "y1" = c(mid, lower - labelRel),
      "group" = rep(firstFacet, 2)
    )

    # rename group name
    colnames(axes)[4] <- groupFacet
    colnames(label)[5] <- groupFacet
  } else {
    print("Please give correct args(mul or one)!")
  }

  ####################################
  # set color value range
  if (is.null(minExp) && is.null(maxExp)) {
    minexp <- 0
    maxexp <- round(max(megredf$scaledValue) + 1, digits = 0)
  } else {
    minexp <- minExp
    maxexp <- maxExp
  }

  ####################################################
  # plot
  pmain <- ggplot2::ggplot(
    megredf,
    ggplot2::aes(x = megredf[, 1], y = megredf[, 2])
  ) +
    ggplot2::geom_point(
      ggplot2::aes(color = scaledValue),
      size = pSize,
      show.legend = show.legend
    ) +
    ggplot2::theme_classic(base_size = base_size) +
    ggplot2::scale_color_gradient(
      name = "", low = low, high = high,
      limits = c(minexp, maxexp),
      na.value = high
    ) +
    ggplot2::labs(x = "", y = "") +
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
      fontface = "italic",
      color = lineTextcol,
      size = cornerTextSize
    ) +
    ggplot2::theme(
      strip.background = ggplot2::element_rect(colour = NA, fill = stripCol),
      strip.text = ggplot2::element_text(size = base_size),
      strip.text.y = ggplot2::element_text(angle = 0),
      aspect.ratio = aspect.ratio,
      legend.position = legendPos,
      plot.title = ggplot2::element_text(hjust = 0.5),
      axis.line = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank()
    )

  ######################################
  # plot layout
  if (is.null(nLayout)) {
    nLayout <- length(features)
  } else {
    nLayout <- nLayout
  }

  ######################################
  # facet plot
  if (is.null(groupFacet)) {
    p1 <- pmain +
      ggplot2::facet_wrap(facets = "gene_name", ncol = nLayout)
  } else {
    p1 <- pmain +
      # ggplot2::facet_grid(facets = c("gene_name", groupFacet))
      ggplot2::facet_grid(rows = ggplot2::vars(.data[["gene_name"]]),
                          cols = ggplot2::vars(.data[[groupFacet]]))
  }

  ######################################
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

  # output
  return(p2)
}
