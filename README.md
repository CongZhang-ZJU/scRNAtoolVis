# scRNAtoolVis 

<img src="man/scRNAtoolVis-logo.png" align="right" height="200" />

 Some useful function to make your scRNA-seq plot more beautiful.



## Installation

```R
install.packages('devtools')
devtools::install_github('junjunlab/scRNAtoolVis')

# if not install ggunchull
devtools::install_github("sajuukLyu/ggunchull", type = "source")

library(scRNAtoolVis)
```

## Seurat v5 兼容说明（重要）

从 Seurat v5 开始，表达矩阵从传统的 `slot`（v4）逐步迁移到 `layer`（v5，Assay5）。  
本仓库已做兼容处理：**包内所有需要读取表达矩阵的位置，都会自动适配 Seurat v4/v5**。

- **你仍然可以像以前一样传 `slot = "data"`**  
  在 Seurat v5 下它会被自动映射为 `layer = "data"`（同理 `counts/scale.data` 也一致）。

### 最小示例（Seurat v5）

```r
library(Seurat)
library(scRNAtoolVis)

# 示例数据（包内自带）
pbmc <- readRDS(system.file("extdata", "htdata.RDS", package = "scRNAtoolVis"))

# 1) DotPlot（本包的 jjDotPlot 兼容 v5 的 layer）
jjDotPlot(object = pbmc, gene = c("LYZ", "MS4A1"), slot = "data")

# 2) 平均表达热图（averageHeatmap 在 v5 下会自动走 layer）
averageHeatmap(object = pbmc, markerGene = c("LYZ", "MS4A1"), slot = "data")

# 3) tracksPlot（已避免直接 @data，兼容 v5 layers）
tracksPlot(object = pbmc, genes = c("LYZ", "MS4A1"), slot = "data")
```

### FeaturePlot/CornerAxes 示例（需要对象包含 umap/tsne 等降维）

```r
library(Seurat)
library(scRNAtoolVis)

obj <- readRDS(system.file("extdata", "seuratTest.RDS", package = "scRNAtoolVis"))

# featurePlot：显式指定 slot，Seurat v5 下会自动按 layer 读取
featurePlot(object = obj, dim = "umap", genes = c("Actb", "Ythdc1"), slot = "data")

# featureCornerAxes：同样支持 slot 参数
featureCornerAxes(object = obj, reduction = "umap", features = c("Actb", "Ythdc1"), slot = "data")
```

## Citation

> Jun Zhang (2022). *scRNAtoolVis: Useful Functions to Make Your scRNA-seq Plot More Cool!.*  <https://github.com/junjunlab/scRNAtoolVis>, <https://junjunlab.github.io/scRNAtoolVis-manual/>.


## More examples see

<https://junjunlab.github.io/scRNAtoolVis-manual/>

---

![image](https://user-images.githubusercontent.com/64965509/198531385-00b0587d-e202-4417-b11d-53cd419594e6.png)

## Related blogs

> - [**scRNAtoolVis 尝试一下?**](https://mp.weixin.qq.com/s?__biz=MzkyMTI1MTYxNA==&mid=2247500252&idx=1&sn=d194022b7b394f1976ef8ad6cb3f8540&chksm=c184fbadf6f372bb37c23610d2beac1eae0421678bf4470876f0b33b59e30fc2034a7cb69a85&token=1253522169&lang=zh_CN#rd)
> - [**scRNAtoolVis 绘制单细胞 Marker 基因均值表达热图**](https://mp.weixin.qq.com/s?__biz=MzkyMTI1MTYxNA==&mid=2247500529&idx=1&sn=daed331779f34479cf301ccdc5d57232&chksm=c184f880f6f37196ad1f0ac4409cc6a6429c8a2654742066aa3df86c22f1e583136d3d5283ee&token=1253522169&lang=zh_CN#rd)
> - [**scRNAtoolVis 0.0.3 版本更新**](https://mp.weixin.qq.com/s?__biz=MzkyMTI1MTYxNA==&mid=2247501432&idx=1&sn=93a7bb3506c911845ef7d10a2fcfdea3&chksm=c184fc09f6f3751f0c2b94f0a6d69692efcc916a45a5409f93b23f90f8db809ca750f87ba234&token=1253522169&lang=zh_CN#rd)
> - [**jjDotPlot 优雅的可视化单细胞基因表达**](https://mp.weixin.qq.com/s?__biz=MzkyMTI1MTYxNA==&mid=2247505159&idx=1&sn=e81b7f890e93419be154b1cf432ae84f&chksm=c184ef76f6f366608d6b2b2d214aecaa4d6c57e976f9fe43a93201a3f6795051395ffc2addf4&token=1253522169&lang=zh_CN#rd)
> - [**jjDotPlot 对基因和亚群排序**](https://mp.weixin.qq.com/s?__biz=MzkyMTI1MTYxNA==&mid=2247505344&idx=1&sn=40b5934ca798e6518172b6c5743c9d95&chksm=c184efb1f6f366a7db3d924f54e9e763c78366021a6d263a502f4d4f4a806250743a74283d5f&token=1253522169&lang=zh_CN#rd)
> - [**averageHeatmap 调整细胞亚群顺序**](https://mp.weixin.qq.com/s?__biz=MzkyMTI1MTYxNA==&mid=2247505603&idx=1&sn=3f85e23cc0eaffa5d92fb4aa03ca0d14&chksm=c184ecb2f6f365a41a428f5600fa9797263bb4432d51be7dc452f213bc1c9be70d50ab669159&token=1253522169&lang=zh_CN#rd)
> - [**jjVolcano 一行代码绘制单细胞火山图**](https://mp.weixin.qq.com/s?__biz=MzkyMTI1MTYxNA==&mid=2247506316&idx=1&sn=307fffe550e25987b148f843e169cbcd&chksm=c184e3fdf6f36aebbe7e948d029f831f03c2906e02a271ee85900a2faac262d2c80bd1b7d2b5&token=1253522169&lang=zh_CN#rd)
> - [**单细胞火山图的旋转和环形**](https://mp.weixin.qq.com/s?__biz=MzkyMTI1MTYxNA==&mid=2247506337&idx=1&sn=5d84b4cbddf053456561d806f3d04737&chksm=c184e3d0f6f36ac65a1cee3d353de3715c30460af0acfbc35446eac61b5b5e9170d7fb1e6315&token=1253522169&lang=zh_CN#rd)
> - [**averageHeatmap 对单细胞 marker 基因聚类**](https://mp.weixin.qq.com/s?__biz=MzkyMTI1MTYxNA==&mid=2247506807&idx=1&sn=f649c782a8d21765f185f03dec0fd9c5&chksm=c184e106f6f3681003142cbb48c5bc406da8458e0153bd00e805ae8fb9c01c6b25c6ff4d6156&token=1253522169&lang=zh_CN#rd)
> - [**听说你想绘制 scanpy 的 tracksPlot?**](https://mp.weixin.qq.com/s?__biz=MzkyMTI1MTYxNA==&mid=2247509478&idx=1&sn=71977a66e74d9dc4fa6aa1b2adbc50be&chksm=c1849f97f6f31681543b475f4df381c9e462dfd617260d5fbeb8956731033c97eb4c20a49f9b&token=143921488&lang=zh_CN#rd)
> - [**用 grid 手搓一个单细胞散点图+细胞数量条形图**](https://mp.weixin.qq.com/s?__biz=MzkyMTI1MTYxNA==&mid=2247509735&idx=1&sn=212e45f9e6c8cfd13bc943b6ac806dc2&chksm=c1849c96f6f31580624e0d4bab65b47d7b3e4cf7fc604881fb9c27a32249bd071dd1eae20d66&token=868587677&lang=zh_CN#rd)
> - [**听说你要绘制 scanpy 版的散点图?**](https://mp.weixin.qq.com/s?__biz=MzkyMTI1MTYxNA==&mid=2247509799&idx=1&sn=b3bc4bbbb8c1e3a98be3f5815cbcb4f0&chksm=c1849d56f6f3144000e22d4b9ebf989ad79349c450928c2a2346749b0187d857fca15d9e7183&token=613129250&lang=zh_CN#rd)
