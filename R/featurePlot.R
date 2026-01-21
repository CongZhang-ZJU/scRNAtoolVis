#' @name featurePlot
#' @author Jun Zhang
#' @title 多基因 FeaturePlot（grid 手工布局版本）
#'
#' @description
#' 从 Seurat 对象中提取降维坐标与基因表达，按给定网格布局绘制多个散点子图。
#' 兼容 Seurat v4/v5：表达读取在 v5 下会自动走 `layer`（本函数保留 `slot` 参数名以兼容旧用法）。
#'
#' @param object Seurat 对象。
#' @param dim 字符串；降维名称（reduction），默认 "umap"。
#' @param genes 字符向量；要绘制的基因/feature 列表。
#' @param slot 读取表达矩阵的层（兼容 Seurat v4/v5），默认 "data"。
#' - Seurat v4：等价于 `FetchData(..., slot = slot)`
#' - Seurat v5：等价于 `FetchData(..., layer = slot)`
#' @param nrow,ncol 整数；子图网格行数/列数。若其中一个为 NULL，会自动根据基因数推断。
#' @param quantile.val 数值；用于截断极端值的分位数（例如 0.99），默认 1 表示不截断。
#' @param color 颜色向量；低-高渐变色，默认使用内置配色。
#' @param rm.axis 逻辑值；是否移除坐标轴，默认 FALSE。
#' @param rm.legend 逻辑值；是否移除每个子图的颜色条，默认 FALSE。
#' @param add.rect 逻辑值；是否为每个子图加边框，默认 FALSE。
#' @param add.corArrow 逻辑值；是否在角落画“方向箭头”（类似 Corner Axes），默认 FALSE。
#' @param add.strip 逻辑值；是否给每个子图加顶部条带，默认 FALSE。
#' @param corLabel.dist 数值；角标文字与箭头的距离。
#' @param arrow.len 数值；箭头长度。
#' @param arrow.label.size 数值；箭头文字字号。
#' @param plot.size 数值；每个子图的相对大小。
#' @param keep.oneCor 逻辑值；是否只保留一套角标箭头，默认 FALSE。
#' @param xlab,ylab 字符串；坐标轴标题，默认 NULL 表示自动。
#' @param respect 逻辑值；是否在布局中保持长宽比，默认 TRUE。
#' @param point.size 数值；点大小（pt），默认 1。
#'
#' @examples
#' \dontrun{
#' library(Seurat)
#' library(scRNAtoolVis)
#'
#' obj <- readRDS(system.file("extdata", "seuratTest.RDS", package = "scRNAtoolVis"))
#' featurePlot(object = obj, dim = "umap", genes = c("Actb", "Ythdc1"), slot = "data")
#' }
#'
#' @import dplyr
#' @import grDevices
#' @import ggplot2
#' @import Seurat
#' @importFrom grid unit viewport pushViewport popViewport grid.rect grid.text
#' grid.points arrow gpar grid.layout
#' @importFrom stats quantile
#'
#' @export

globalVariables(c("col_rg", "tmp_col"))

featurePlot <- function(
    object = NULL,
    dim = "umap",
    genes = NULL,
    slot = "data",
    nrow = NULL,
    ncol = NULL,
    quantile.val = 1,
    color = NULL,
    rm.axis = FALSE,
    rm.legend = FALSE,
    add.rect = FALSE,
    add.corArrow = FALSE,
    add.strip = FALSE,
    corLabel.dist = 0.08,
    arrow.len = 0.2,
    arrow.label.size = 6,
    plot.size = 0.6,
    keep.oneCor = FALSE,
    xlab = NULL,
    ylab = NULL,
    respect = TRUE,
    point.size = 1) {
  # ============================================================================
  # 1_extract data
  # ============================================================================
  # 逐行注释：
  # - Seurat 的降维结果存放在 reductions 里（例如 pca/umap/tsne）
  # - 如果对象里没有对应 reduction，`Seurat::Embeddings()` 会直接报错
  # - 这里提前做一次更友好的检查，并给出中文提示
  available_reductions <- Seurat::Reductions(object)
  if (length(available_reductions) == 0) {
    stop(
      "当前 Seurat 对象中没有任何降维结果（reductions）。\n",
      "请先运行例如 `RunPCA()` / `RunUMAP()` / `RunTSNE()`，或把 `dim` 改为对象中已有的 reduction。"
    )
  }

  # 逐行注释：
  # - 若用户指定的 dim 不存在，给出可用列表
  if (!dim %in% available_reductions) {
    stop(
      sprintf("未在对象中找到 reduction='%s'。可用 reductions: %s",
              dim, paste(available_reductions, collapse = ", "))
    )
  }

  # make PC data
  reduc <- data.frame(Seurat::Embeddings(object, reduction = dim))

  # metadata
  meta <- object@meta.data

  # combine
  pc12 <- cbind(reduc, meta)
  # 逐行注释：显式使用 Seurat::Idents，避免与其他包函数重名
  pc12$idents <- Seurat::Idents(object)

  # get gene expression
  # 逐行注释：
  # - Seurat v5 下 `FetchData(slot=...)` 已经被废弃并会直接报错
  # - 因此这里统一走 `.scRNAtoolVis_fetch_data()`：
  #   - v5：自动使用 `layer = slot`
  #   - v4：自动回退到 `slot = slot`
  geneExp <- .scRNAtoolVis_fetch_data(
    object = object,
    vars = genes,
    slot = slot
  )

  # cbind
  mer <- cbind(pc12, geneExp)

  # get nrow and ncol
  if (is.null(nrow)) {
    nrow <- ifelse(is.null(ncol), 1, ceiling(length(genes) / ncol))
  }

  if (is.null(ncol)) {
    ncol <- ifelse(is.null(nrow), length(genes), ceiling(length(genes) / nrow))
  }

  gene_mtx <- suppressWarnings(matrix(genes, nrow = nrow, ncol = ncol))

  # assign colors
  if (is.null(color)) {
    cols <- c("grey90", "#57C5B6", "#159895", "#1A5F7A", "#002B5B")
  } else {
    cols <- color
  }

  # ============================================================================
  # 2_draw plot
  # ============================================================================
  if (rm.axis == FALSE) {
    lab.shift <- unit(-2.5, "lines")
  } else {
    lab.shift <- unit(-1, "lines")
  }

  # CANVAS FOR PLOT
  grid.newpage()
  pushViewport(
    viewport(
      x = 0.5, y = 0.5,
      width = 0.9, height = 0.9,
      xscale = range(mer[, 1]), yscale = range(mer[, 2]),
      layout = grid.layout(nrow = nrow, ncol = ncol, respect = respect)
    )
  )

  # loop
  for (i in 1:nrow) {
    for (j in 1:ncol) {
      # check genes numbers
      if (i * j > length(genes)) {
        break
      }

      # ===========================================================
      # 1_panel grid
      pushViewport(
        viewport(layout.pos.row = i, layout.pos.col = j)
      )

      if (add.rect == TRUE) {
        grid.rect()
      }

      # process data
      quantile_val <- quantile(mer[, gene_mtx[i, j]], probs = quantile.val)
      mer <- mer |>
        dplyr::mutate(tmp_col = if_else(.data[[gene_mtx[i, j]]] > quantile_val,
          quantile_val,
          .data[[gene_mtx[i, j]]]
        ))

      tmp_data <- mer |>
        dplyr::arrange(tmp_col)
      col_p <- colorRampPalette(cols)(100)
      cut_range <- cut(tmp_data[, "tmp_col"], 100)

      labs <- levels(cut_range)
      names(labs) <- col_p

      tmp_data <- tmp_data |>
        dplyr::mutate(col_rg = as.character(cut_range)) |>
        dplyr::mutate(col_f = ifelse(col_rg %in% labs, names(labs)[match(col_rg, labs)], "black"))

      # ===========================================================
      # 2_scatter plot
      pushViewport(
        viewport(
          x = 0.5, y = 0.5, width = plot.size, height = plot.size,
          xscale = extendrange(range(tmp_data[, 1]), f = 0.05),
          yscale = extendrange(range(tmp_data[, 2]), f = 0.05)
        )
      )

      # whether add corner arrows
      if (keep.oneCor == TRUE) {
        if (j == 1) {
          if (add.corArrow == TRUE) {
            grid.segments(
              x0 = 0, x1 = arrow.len, y0 = 0, y1 = 0,
              arrow = arrow(length = unit(2, "mm"), type = "closed"),
              gp = gpar(fill = "black")
            )
            grid.text(
              label = paste0(toupper(dim), " 1"), x = arrow.len / 2, y = -corLabel.dist,
              gp = gpar(fontsize = arrow.label.size, fontface = "bold.italic")
            )
            grid.segments(
              x0 = 0, x1 = 0, y0 = 0, y1 = arrow.len,
              arrow = arrow(length = unit(2, "mm"), type = "closed"),
              gp = gpar(fill = "black")
            )
            grid.text(
              label = paste0(toupper(dim), " 2"),
              x = -corLabel.dist, y = arrow.len / 2, rot = 90,
              gp = gpar(fontsize = arrow.label.size, fontface = "bold.italic")
            )
          } else {
            grid.rect()
          }
        }
      } else {
        if (add.corArrow == TRUE) {
          grid.segments(
            x0 = 0, x1 = arrow.len, y0 = 0, y1 = 0,
            arrow = arrow(length = unit(2, "mm"), type = "closed"),
            gp = gpar(fill = "black")
          )
          grid.text(
            label = paste0(toupper(dim), " 1"), x = arrow.len / 2, y = -corLabel.dist,
            gp = gpar(fontsize = arrow.label.size, fontface = "bold.italic")
          )
          grid.segments(
            x0 = 0, x1 = 0, y0 = 0, y1 = arrow.len,
            arrow = arrow(length = unit(2, "mm"), type = "closed"),
            gp = gpar(fill = "black")
          )
          grid.text(
            label = paste0(toupper(dim), " 2"),
            x = -corLabel.dist, y = arrow.len / 2, rot = 90,
            gp = gpar(fontsize = arrow.label.size, fontface = "bold.italic")
          )
        } else {
          grid.rect()
        }
      }

      grid.points(
        x = tmp_data[, 1], y = tmp_data[, 2], pch = 19, size = unit(point.size, "pt"),
        gp = gpar(col = tmp_data$col_f)
      )

      # whether draw axis
      if (add.corArrow == FALSE) {
        if (rm.axis == FALSE) {
          # grid.xaxis()
          # grid.yaxis()
          jjPlot::grid.xaxis2(label.space = 0.5)
          jjPlot::grid.yaxis2(label.space = 0.25)
        }
      }

      # add strip
      if (add.strip == TRUE) {
        grid.rect(
          x = 0.5, y = 1, width = 1,
          height = 0.15, gp = gpar(fill = "grey85"),
          just = "bottom"
        )
      }

      grid.text(
        label = gene_mtx[i, j], x = 0.5, y = unit(1 + 0.15 / 2, "npc"),
        gp = gpar(fontface = "bold.italic")
      )
      if (add.corArrow == FALSE) {
        # axis labels
        if (!is.null(xlab) || !is.null(ylab)) {
          axis.label.x <- xlab
          axis.label.y <- ylab
        } else {
          axis.label.x <- paste0(toupper(dim), " dimension 1")
          axis.label.y <- paste0(toupper(dim), " dimension 2")
        }

        grid.text(label = axis.label.x, x = 0.5, y = lab.shift)
        grid.text(label = axis.label.y, x = lab.shift, y = 0.5, rot = 90)
      }

      popViewport()

      # ===========================================================
      # 3_draw legend
      if (rm.legend == FALSE) {
        pushViewport(
          viewport(
            x = 0.5 + plot.size / 2 + 0.01, y = 0.5,
            width = 0.025, height = unit(plot.size, "npc"),
            just = "left",
            yscale = range(tmp_data[, gene_mtx[i, j]])
          )
        )
        # grid.rect(x = 0.5, y = unit(seq(0.25,0.75, length = 100), "npc"),
        #           width = unit(1, "npc"), height = unit(0.5, "npc"),
        #           just = "centre",default.units = "npc",
        #           gp = gpar(col = NA, fill = col_p))
        # grid.rect(gp = gpar(fill = NA))
        # # grid.yaxis(main = FALSE)
        # jjPlot::grid.yaxis2(side = "right",tick.len = 0.25)

        jjPlot::grid.colorkey(
          x = tmp_data[, gene_mtx[i, j]],
          color = cols,
          pos = "v",
          ticks.side = "right"
        )

        popViewport()
      }
      popViewport()
    }
  }
}
