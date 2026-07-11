library(DEseq)
#转化R环境只需要在设置里面搜索r.rterm.windows，将对应的R.exe填入对应的文件就行
#设置打开方式：左上角 文件（File）→ 首选项（Preferences）→ 设置（Settings）

sce <- merge(scRNAlist[[1]], y= scRNAlist[ -1 ]) 
sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^mt-")
sce <- subset(sce, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & nCount_RNA<15000 & percent.mt < 10)
sce <-  sce %>% 
  NormalizeData() %>%  
  FindVariableFeatures(nfeatures = 2000) %>%  
  ScaleData() %>%  
  RunPCA(features = VariableFeatures(.), npcs = 60)

sce <- IntegrateLayers(
  object = sce, 
  method = HarmonyIntegration, 
  orig.reduction = "pca", 
  new.reduction = "harmony"
)

sce_int <- JoinLayers(sce)
Idents(sce_int) <- "RNA_snn_res.0.4"
sce.markers <- FindAllMarkers(object = sce_int, 
                              only.pos = TRUE, #只保留上调基因
                              logfc.threshold = 0.25,#log2FC阀值
                              min.pct = 0.1  #基因最少在%1的基因表达才纳入计算
)

print("helllo")
R.home("bin")
