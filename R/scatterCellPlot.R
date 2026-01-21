#' @name scatterCellPlot
#' @author Jun Zhang
#' @title 使用 grid 绘制“降维散点 + 细胞数量条形图 + 图例”的组合图
#'
#' @description
#' 该函数不依赖 ggplot2 的 facet/legend 系统，而是用 `grid` 手动布局：
#' - 左侧：降维散点（UMAP/tSNE 等）
#' - 中间：每个细胞类型的细胞数量条形图（可选）
#' - 右侧：颜色图例（点+文本）
#'
#' @param object Seurat 对象，必须包含指定的降维结果（例如 "umap"）。
#' @param color 颜色向量；每个 celltype 一种颜色。默认 NULL 表示随机生成。
#' @param dim 字符串；降维名称（Seurat reductions），默认 "umap"。
#' @param rm.axis 是否移除坐标轴（TRUE 时会用箭头标注轴方向），默认 FALSE。
#' @param cell.id 字符串；可选，meta.data 中的列名。若提供，则图例左侧会额外显示该列值。
#' @param bar.width 数值；条形图区域宽度，默认 0.2。
#' @param point.size 数值；散点大小（grid 点大小，单位 pt），默认 1。
#' @param rm.barplot 是否移除中间条形图，默认 FALSE。
#' @param legend.psize 数值；图例点大小（单位 char），默认 1.5。
#' @param arrow.len 数值；当 `rm.axis=TRUE` 时箭头长度，默认 0.2。
#'
#' @return 本函数直接在当前图形设备上绘制，不返回对象（返回值为 NULL/不可用）。
#'
#' @examples
#' \dontrun{
#' library(Seurat)
#' library(scRNAtoolVis)
#'
#' obj <- readRDS(system.file("extdata", "seuratTest.RDS", package = "scRNAtoolVis"))
#' scatterCellPlot(object = obj, dim = "umap")
#' }
#'
#' @importFrom grid grid.newpage pushViewport popViewport viewport grid.rect grid.xaxis grid.yaxis grid.points grid.segments grid.text arrow gpar
#'
#' @export

globalVariables(c("idents"))

scatterCellPlot <- function(
    object = NULL,
    color = NULL,
    dim = "umap",
    rm.axis = FALSE,
    cell.id = NULL,
    bar.width = 0.2,
    point.size = 1,
    rm.barplot = FALSE,
    legend.psize = 1.5,
    arrow.len = 0.2) {
  # ============================================================================
  # 1) 提取数据
  # ============================================================================
  # 逐行注释：检查 Seurat 对象是否包含指定的降维结果
  available_reductions <- Seurat::Reductions(object)
  if (length(available_reductions) == 0) {
    stop(
      "当前 Seurat 对象中没有任何降维结果（reductions）。\n",
      "请先运行例如 `RunUMAP()`/`RunTSNE()`，或把 `dim` 改为对象中已有的 reduction。"
    )
  }
  if (!dim %in% available_reductions) {
    stop(sprintf("未在对象中找到 reduction='%s'。可用 reductions: %s",
                 dim, paste(available_reductions, collapse = ", ")))
  }

  # make PC data
  reduc <- data.frame(Seurat::Embeddings(object, reduction = dim))

  # metadata
  meta <- object@meta.data

  # combine
  pc12 <- cbind(reduc, meta)
  pc12$idents <- as.character(Seurat::Idents(object))

  # summary celltype numbers
  if (is.null(cell.id)) {
    cell_num <- pc12 |>
      dplyr::group_by(idents) |>
      dplyr::summarise(n = dplyr::n()) |>
      dplyr::arrange(n)
  } else {
    cell_num <- pc12 |>
      dplyr::group_by(idents, .data[[cell.id]]) |>
      dplyr::summarise(n = dplyr::n()) |>
      dplyr::arrange(n)
  }

  # ============================================================================
  # 2_draw plot
  # ============================================================================
  # 逐行注释：根据是否移除坐标轴，调整轴标签/箭头的留白位置
  if (rm.axis == FALSE) {
    lab.shift <- unit(-2.5, "lines")
  } else {
    lab.shift <- unit(-1, "lines")
  }

  grid.newpage()
  pushViewport(
    viewport(
      x = unit(0.1, "npc"), y = unit(0.5, "npc"),
      width = unit(0.5, "npc"),
      height = unit(0.7, "npc"),
      just = "left",
      xscale = grDevices::extendrange(range(pc12[, 1]), f = 0.05),
      yscale = grDevices::extendrange(range(pc12[, 2]), f = 0.05),
    )
  )
  grid.rect()
  if (rm.axis == FALSE) {
    # grid.xaxis()
    # grid.yaxis()
    jjPlot::grid.xaxis2(label.space = 0.5)
    jjPlot::grid.yaxis2(label.space = 0.25)
  }

  celltype <- cell_num$idents

  if (is.null(color)) {
    # create colors
    cols <- circlize::rand_color(n = length(celltype))
  } else {
    cols <- color
  }

  # draw points
  # i = 1
  for (i in seq_along(celltype)) {
    # tmp <- pc12 |>
    #   dplyr::filter(idents == celltype[i])
    tmp <- pc12[which(pc12$idents %in% celltype[i]), ]

    grid.points(
      x = tmp[, 1], y = tmp[, 2], pch = 19, size = unit(point.size, "pt"),
      gp = gpar(col = cols[i])
    )
  }

  # arrow
  if (rm.axis == TRUE) {
    grid.segments(
      x0 = 0.025, x1 = arrow.len, y0 = 0.05, y1 = 0.05,
      arrow = arrow(length = unit(2, "mm"), type = "closed"),
      gp = gpar(fill = "black")
    )
    grid.text(
      label = paste0(toupper(dim), " 1"),
      x = (arrow.len + 0.025) / 2, y = 0.025,
      gp = gpar(fontsize = 6, fontface = "bold.italic")
    )
    grid.segments(
      x0 = 0.05, x1 = 0.05, y0 = 0.025, y1 = arrow.len,
      arrow = arrow(length = unit(2, "mm"), type = "closed"),
      gp = gpar(fill = "black")
    )
    grid.text(
      label = paste0(toupper(dim), " 2"),
      x = 0.025, y = (arrow.len + 0.025) / 2, rot = 90,
      gp = gpar(fontsize = 6, fontface = "bold.italic")
    )
  } else {
    # labs
    grid.text(label = paste0(toupper(dim), " dimension 1"), x = 0.5, y = lab.shift)
    grid.text(label = paste0(toupper(dim), " dimension 2"), x = lab.shift, y = 0.5, rot = 90)
  }

  popViewport()

  # ============================================================================
  # barplot
  if (isFALSE(rm.barplot)) {
    pushViewport(
      viewport(
        x = unit(0.61, "npc"), y = unit(0.5, "npc"),
        width = unit(bar.width, "npc"),
        height = unit(0.7, "npc"),
        just = "left",
        yscale = c(0, nrow(cell_num) + 0.75),
        xscale = c(0, max(cell_num$n) + 0.1 * max(cell_num$n))
      )
    )

    if (rm.axis == FALSE) {
      # grid.xaxis()
      jjPlot::grid.xaxis2(
        label.space = 0.5,
        at = c(0, max(cell_num$n)),
        labels = as.character(c(0, max(cell_num$n)))
      )
    }
    grid.rect(
      x = rep(0, nrow(cell_num)), y = seq_len(nrow(cell_num)),
      width = cell_num$n, height = unit(0.08, "npc"),
      just = "left",
      gp = gpar(fill = cols, col = NA),
      default.units = "native"
    )
    grid.rect(gp = gpar(fill = "transparent"))
    grid.text(label = "Number of cells", x = 0.5, y = lab.shift)
    popViewport()
  }

  # ============================================================================
  # legend
  if (isTRUE(rm.barplot)) {
    bar.width <- 0
  }

  pushViewport(
    viewport(
      x = unit(0.61 + bar.width, "npc"), y = unit(0.5, "npc"),
      width = unit(0.2, "npc"),
      height = unit(0.7, "npc"),
      just = "left",
      yscale = c(0, nrow(cell_num) + 0.75)
    )
  )

  grid.points(
    x = rep(0.1, nrow(cell_num)), y = seq_len(nrow(cell_num)), pch = 19,
    gp = gpar(col = cols), size = unit(legend.psize, "char")
  )
  if (!is.null(cell.id)) {
    grid.text(
      label = as.character(unlist(cell_num[, cell.id])),
      x = 0.1, y = seq_len(nrow(cell_num)),
      default.units = "native"
    )
  }
  grid.text(
    label = cell_num$idents,
    x = 0.2, y = seq_len(nrow(cell_num)),
    just = "left",
    gp = gpar(fontsize = 10),
    default.units = "native"
  )
  # grid.rect(gp = gpar(fill = "transparent"))
  popViewport()
}
