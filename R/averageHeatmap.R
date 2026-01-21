#' @name averageHeatmap
#' @author Junjun Lao
#' @title 平均表达热图（按分组聚合）
#'
#' @description
#' 计算指定基因在不同分组（cluster/ident/celltype 等）中的平均表达量，并绘制 Z-score 热图。
#' 兼容 Seurat v4/v5：Seurat v5 会自动使用 `layer`（本函数仍保留 `slot` 参数名以兼容旧用法）。
#'
#' @param object Seurat 对象。
#' @param markerGene 字符向量；要展示的基因列表（marker genes）。
#' @param group.by 字符串；分组依据（例如 "ident"/"seurat_clusters"/"celltype"），默认 "ident"。
#' @param assays 字符串；使用哪个 assay，默认 "RNA"。
#' @param slot 字符串；表达矩阵层（兼容 v4/v5）。常用 "data"/"counts"/"scale.data"，默认 "data"。
#' @param htCol 颜色向量；热图颜色，默认 c("#0099CC", "white", "#CC0033")。
#' @param colseed 整数；随机颜色种子（用于 cluster 注释色），默认 666。
#' @param htRange 数值向量；热图取值范围（对应 htCol），默认 c(-2, 0, 2)。
#' @param annoCol 逻辑值；是否使用自定义注释颜色，默认 FALSE。
#' @param myanCol 颜色向量；自定义注释颜色（当 annoCol=TRUE 时使用），默认 NULL。
#' @param annoColType 字符串；随机颜色亮度风格，默认 "light"。
#' @param annoColTypeAlpha 数值；随机颜色透明度，默认 0。
#' @param row_title 字符串；行标题，默认 "Cluster top Marker genes"。
#' @param clusterAnnoName 逻辑值；是否显示顶部注释名称，默认 TRUE。
#' @param showRowNames 逻辑值；是否显示行名（基因名），默认 TRUE。
#' @param row_names_side 字符串；行名位置，默认 "left"。
#' @param markGenes 字符向量；需要在右侧标注的基因，默认 NULL。
#' @param border 逻辑值；是否显示边框，默认 FALSE。
#' @param fontsize 数值；行名字号，默认 10。
#' @param column_names_rot 数值；列名旋转角度，默认 45。
#' @param width,height 数值；热图主体宽高（cm），默认 NULL 表示自动。
#' @param cluster.order 字符向量；列（cluster）顺序，默认 NULL。
#' @param cluster_columns,cluster_rows 逻辑值；是否聚类列/行，默认 FALSE。
#' @param gene.order 字符向量；行（gene）顺序，默认 NULL。
#' @param ... 其他参数会传给 `ComplexHeatmap::rowAnnotation()` 与 `ComplexHeatmap::Heatmap()`。
#'
#' @return 返回一个 ComplexHeatmap 对象。
#' @export
#'
#' @examples
#' \dontrun{
#' library(Seurat)
#' library(scRNAtoolVis)
#'
#' obj <- readRDS(system.file("extdata", "htdata.RDS", package = "scRNAtoolVis"))
#' averageHeatmap(object = obj, markerGene = c("LYZ", "MS4A1"), slot = "data")
#' }
#'
#' @examples
#' httest <- system.file("extdata", "htdata.RDS", package = "scRNAtoolVis")
#' pbmc <- readRDS(httest)
#'
#' # load markergene
#' markergene <- system.file("extdata", "top5pbmc.markers.csv", package = "scRNAtoolVis")
#' markers <- read.table(markergene, sep = ",", header = TRUE)
#'
#' # plot
#' averageHeatmap(
#'   object = pbmc,
#'   markerGene = markers$gene
#' )
#'
#' # change color
#' averageHeatmap(
#'   object = pbmc,
#'   markerGene = markers$gene,
#'   htCol = c("#339933", "#FFCC00", "#FF0033")
#' )
#'
# define function
averageHeatmap <- function(
    object = NULL,
    markerGene = NULL,
    group.by = "ident",
    assays = "RNA",
    slot = "data",
    htCol = c("#0099CC", "white", "#CC0033"),
    colseed = 666,
    htRange = c(-2, 0, 2),
    annoCol = FALSE,
    myanCol = NULL,
    annoColType = "light",
    annoColTypeAlpha = 0,
    row_title = "Cluster top Marker genes",
    clusterAnnoName = TRUE,
    showRowNames = TRUE,
    row_names_side = "left",
    markGenes = NULL,
    border = FALSE,
    fontsize = 10,
    column_names_rot = 45,
    width = NULL,
    height = NULL,
    cluster.order = NULL,
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    gene.order = NULL,
    ...) {
  # get cells mean gene expression
  # 逐行注释：
  # - AverageExpression 在 v5 推荐使用 layer；在 v4 使用 slot
  # - 为避免“仅靠 formals 判断不准/不同版本行为差异”，这里采用与兼容层一致的策略：
  #   先尝试 layer，失败且错误信息表明 layer 不被支持时，再回退 slot
  mean_gene_exp <- tryCatch(
    as.matrix(
      data.frame(
        Seurat::AverageExpression(
          object,
          features = markerGene,
          group.by = group.by,
          assays = assays,
          layer = slot
        )
      )
    ),
    error = function(e) {
      msg <- conditionMessage(e)
      if (grepl("unused argument|formal argument", msg, ignore.case = TRUE)) {
        return(
          as.matrix(
            data.frame(
              Seurat::AverageExpression(
                object,
                features = markerGene,
                group.by = group.by,
                assays = assays,
                slot = slot
              )
            )
          )
        )
      }
      stop(e)
    }
  )

  # filter gene
  mean_gene_exp <- mean_gene_exp[rowSums(mean_gene_exp) > 0,]

  # add colnames
  # name1 <- gsub(
  #   pattern = paste0(assays, ".", sep = ""),
  #   replacement = "",
  #   colnames(mean_gene_exp)
  # )
  #
  # colnames(mean_gene_exp) <- gsub(
  #   pattern = "\\.",
  #   replacement = " ", name1
  # )

  colnames(mean_gene_exp) <- levels(Seurat::Idents(object))

  # Z-score
  htdf <- t(scale(t(mean_gene_exp), scale = TRUE, center = TRUE))

  # cluster order
  if (!is.null(cluster.order)) {
    htdf <- htdf[, cluster.order]
  }

  # gene order
  if (!is.null(gene.order)) {
    htdf <- htdf[gene.order, ]
  }

  # color
  col_fun <- circlize::colorRamp2(htRange, htCol)

  # anno color
  if (annoCol == FALSE) {
    set.seed(colseed)
    anno_col <- circlize::rand_color(
      ncol(htdf),
      luminosity = annoColType,
      transparency = annoColTypeAlpha
    )
    print(c("Your cluster annotation color is:", anno_col))
  } else if (annoCol == TRUE) {
    # give your own color vectors
    anno_col <- myanCol
  } else {
    print("Give TRUE or FALSE paramters!")
  }
  names(anno_col) <- colnames(htdf)

  # top annotation
  column_ha <- ComplexHeatmap::HeatmapAnnotation(
    cluster = colnames(htdf),
    show_legend = FALSE,
    show_annotation_name = clusterAnnoName,
    col = list(cluster = anno_col)
  )

  # whether mark your genes on plot
  if (!is.null(markGenes)) {
    # all genes
    rowGene <- rownames(htdf)

    # tartget gene
    annoGene <- markGenes

    # get target gene index
    index <- match(annoGene, rowGene)

    # some genes annotation
    geneMark <- ComplexHeatmap::rowAnnotation(
      gene = ComplexHeatmap::anno_mark(
        at = index,
        labels = annoGene,
        labels_gp = grid::gpar(
          fontface = "italic",
          fontsize = fontsize
        )
      ),
      ...
    )

    right_annotation <- geneMark
  } else {
    right_annotation <- NULL
  }

  # control heatmap width and height
  if (is.null(width) || is.null(height)) {
    # plot
    ComplexHeatmap::Heatmap(
      htdf,
      show_row_dend = TRUE,
      show_column_dend = TRUE,
      name = "Z-score",
      cluster_columns = cluster_columns,
      cluster_rows = cluster_rows,
      row_title = row_title,
      # column_title = "Clusters",
      right_annotation = right_annotation,
      show_row_names = showRowNames,
      row_names_gp = grid::gpar(
        fontface = "italic",
        fontsize = fontsize
      ),
      row_names_side = row_names_side,
      border = border,
      column_names_side = "top",
      column_names_rot = column_names_rot,
      top_annotation = column_ha,
      col = col_fun,
      ...
    )
  } else {
    # plot
    ComplexHeatmap::Heatmap(
      htdf,
      show_row_dend = TRUE,
      show_column_dend = TRUE,
      name = "Z-score",
      cluster_columns = FALSE,
      cluster_rows = FALSE,
      row_title = row_title,
      # column_title = "Clusters",
      right_annotation = right_annotation,
      show_row_names = showRowNames,
      row_names_gp = grid::gpar(
        fontface = "italic",
        fontsize = fontsize
      ),
      row_names_side = row_names_side,
      border = border,
      column_names_side = "top",
      column_names_rot = column_names_rot,
      top_annotation = column_ha,
      col = col_fun,
      width = ggplot2::unit(width, "cm"),
      height = ggplot2::unit(height, "cm"),
      ...
    )
  }
}
