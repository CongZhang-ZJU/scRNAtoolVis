#' Seurat v4/v5 兼容工具函数
#'
#' 说明（重要）：
#' - Seurat v5 引入了 Assay5 / layers（例如 "counts"、"data"、"scale.data"）
#' - Seurat v4 主要使用 slot（例如 slot="data"）
#' - 为了让本包同时兼容 Seurat v4/v5，这里统一提供“读取表达矩阵/FetchData”的兼容封装
#'
#' 注意：
#' - 这些函数是包内部使用的，不导出（不在 NAMESPACE 里暴露）
#' - 所有参数都尽量保留与 Seurat 常用接口一致的命名习惯，便于理解与维护
#'
#' @keywords internal
#' @noRd

# ------------------------------------------------------------------------------
# 统一读取 Seurat 表达矩阵（GetAssayData 的 slot/layer 兼容）
# ------------------------------------------------------------------------------
.scRNAtoolVis_get_assay_data <- function(object,
                                        assay = NULL,
                                        slot = "data",
                                        layer = NULL) {
  # 逐行注释：
  # - assay 为空时，默认使用 Seurat 当前 DefaultAssay
  assay <- assay %||% Seurat::DefaultAssay(object = object)

  # 逐行注释：
  # - 如果用户没有显式给 layer，则保持与旧参数 slot 一致（layer=slot）
  # - 这样用户继续传 slot="data" 时，在 v5 下会自动映射到 layer="data"
  layer <- layer %||% slot

  # 逐行注释：
  # - 这里采用“先尝试 layer，再回退 slot”的策略
  # - 原因：在 Seurat v5 中 slot 可能已经变成 defunct（直接报错）
  # - 仅靠 formals() 判断并不可靠（很多函数是 S3 泛型，形参不一定包含 layer）
  res <- tryCatch(
    Seurat::GetAssayData(object = object, assay = assay, layer = layer),
    error = function(e) e
  )

  # 逐行注释：如果 layer 方式成功，直接返回
  if (!inherits(res, "error")) {
    return(res)
  }

  # 逐行注释：
  # - 只有在“layer 参数不被支持”的情况下，才回退到 slot
  # - 典型报错信息包含 "unused argument" 或 "formal argument"
  msg <- conditionMessage(res)
  if (grepl("unused argument|formal argument", msg, ignore.case = TRUE)) {
    return(Seurat::GetAssayData(object = object, assay = assay, slot = slot))
  }

  # 逐行注释：如果不是“layer 不支持”的错误，原样抛出，方便定位真实问题
  stop(res)
}

# ------------------------------------------------------------------------------
# 统一 FetchData（FetchData 的 slot/layer 兼容）
# ------------------------------------------------------------------------------
.scRNAtoolVis_fetch_data <- function(object,
                                    vars,
                                    slot = "data",
                                    layer = NULL,
                                    ...) {
  # 逐行注释：
  # - layer 为空时，默认与 slot 保持一致
  layer <- layer %||% slot

  # 逐行注释：
  # - 同样采用“先尝试 layer，再回退 slot”的策略
  # - 在 Seurat v5 中，slot 参数可能已经 defunct（会直接中断）
  res <- tryCatch(
    Seurat::FetchData(object = object, vars = vars, layer = layer, ...),
    error = function(e) e
  )

  # 逐行注释：如果 layer 方式成功，直接返回
  if (!inherits(res, "error")) {
    return(res)
  }

  # 逐行注释：
  # - 只有在“layer 参数不被支持”的情况下，才回退到 slot
  # - Seurat v4 常见报错：unused argument (layer = ...)
  msg <- conditionMessage(res)
  if (grepl("unused argument|formal argument", msg, ignore.case = TRUE)) {
    return(Seurat::FetchData(object = object, vars = vars, slot = slot, ...))
  }

  # 逐行注释：其他错误直接抛出，避免掩盖真实问题
  stop(res)
}

# ------------------------------------------------------------------------------
# 将表达矩阵稳定转为 base matrix（兼容 dgCMatrix 等稀疏矩阵）
# ------------------------------------------------------------------------------
.scRNAtoolVis_as_matrix <- function(x) {
  # 逐行注释：
  # - Seurat 常用稀疏矩阵（dgCMatrix）直接转 matrix 可能触发警告
  # - 这里统一用 suppressWarnings 包裹，不影响结果但减少无意义的提示
  suppressWarnings(as.matrix(x))
}

