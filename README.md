#C:/Users/ION/Documents에 Donor 1~4 폴더들에 각각 h5파일들 넣어 세팅

#mkdir Donor1
#curl.exe -L "https://cf.10xgenomics.com/samples/cell-exp/9.0.0/5k_Human_Donor1_PBMC_3p_gem-x_5k_Human_Donor1_PBMC_3p_gem-x/5k_Human_Donor1_PBMC_3p_gem-x_5k_Human_Donor1_PBMC_3p_gem-x_count_sample_filtered_feature_bc_matrix.h5" -o "Donor1\filtered_feature_bc_matrix.h5"

#작업 경로 설정
setwd("C:/Users/ION/Documents")

#패키지 불러오기
library(Seurat)
library(tidyverse)
library(patchwork)

#Donor1 파일 경로 지정, !!변수있음!!
base_dir <- "C:/Users/ION/Documents"

#변수
sample_id <- "Donor1"

data_file <- file.path(
  base_dir,
  sample_id,
  "filtered_feature_bc_matrix.h5"
)

out_dir <- file.path(
  base_dir,
  sample_id,
  "outputs"
)

dir.create(out_dir, showWarnings = FALSE)

#파일 경로 확인
data_file
file.exists(data_file)

#10x h5 데이터 불러오기
counts <- Read10X_h5(data_file)
class(counts)
dim(counts)

#Seurat object 만들기
donor <- CreateSeuratObject(
  counts = counts,
  project = sample_id,
  min.cells = 3,
  min.features = 200
)
donor
head(donor@meta.data)

#metadata 추가 !!변수있음!!
donor$sample_id <- "Donor1"
donor$sex <- "Male"
donor$age_group <- "18-35"
head(donor@meta.data)

#mitochondrial read 비율 계산
donor[["percent.mt"]] <- PercentageFeatureSet(
  donor,
  pattern = "^MT-"
)
summary(donor$percent.mt)

#QC 전 plot 확인
p_vln_before <- VlnPlot(
  donor,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3,
  pt.size = 0.1
)

print(p_vln_before)

ggsave(
  filename = file.path(out_dir, paste0(sample_id, "_VlnPlot_before_filtering.png")),
  plot = p_vln_before,
  width = 12,
  height = 8,
  dpi = 300
)

#QC scatter plot
p_scatter_before_1 <- FeatureScatter(
  donor,
  feature1 = "nCount_RNA",
  feature2 = "percent.mt"
)

p_scatter_before_2 <- FeatureScatter(
  donor,
  feature1 = "nCount_RNA",
  feature2 = "nFeature_RNA"
)

p_scatter_before <- p_scatter_before_1 + p_scatter_before_2

print(p_scatter_before)

ggsave(
  filename = file.path(out_dir, paste0(sample_id, "_ScatterPlot_before_filtering.png")),
  plot = p_scatter_before,
  width = 10,
  height = 5,
  dpi = 300
)

#QC filtering default 세팅으로 한번
min_features <- 300
max_features <- 6500
min_counts <- 500
max_counts <- 30000
max_percent_mt <- 10

#QC filtering 실행
donor_filtered <- subset(
  donor,
  subset =
    nFeature_RNA > min_features &
    nFeature_RNA < max_features &
    nCount_RNA > min_counts &
    nCount_RNA < max_counts &
    percent.mt < max_percent_mt
)
ncol(donor)
ncol(donor_filtered)

removed_percent <- (1 - ncol(donor_filtered) / ncol(donor)) * 100
removed_percent

#Filtering 후 QC 확인
VlnPlot(
  donor_filtered,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3,
  pt.size = 0.1
)

#QC filtering plot 보고 직접 수정하기 !!중요한 변수!!
min_features <- 2000
max_features <- 5000
min_counts <- 3000
max_counts <- 20000
max_percent_mt <- 7.5

#수정한 QC filtering 실행
donor_filtered <- subset(
  donor,
  subset =
    nFeature_RNA > min_features &
    nFeature_RNA < max_features &
    nCount_RNA > min_counts &
    nCount_RNA < max_counts &
    percent.mt < max_percent_mt
)
ncol(donor)
ncol(donor_filtered)

removed_percent <- (1 - ncol(donor_filtered) / ncol(donor)) * 100
removed_percent

#수정한 Filtering 후 QC 확인
p_vln_after <- VlnPlot(
  donor_filtered,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3,
  pt.size = 0.1
)

print(p_vln_after)

ggsave(
  filename = file.path(out_dir, paste0(sample_id, "_VlnPlot_after_filtering.png")),
  plot = p_vln_after,
  width = 12,
  height = 8,
  dpi = 300
)

#수정한 Filtering 후 scatter plot 확인
p_scatter_after_1 <- FeatureScatter(
  donor_filtered,
  feature1 = "nCount_RNA",
  feature2 = "percent.mt"
)

p_scatter_after_2 <- FeatureScatter(
  donor_filtered,
  feature1 = "nCount_RNA",
  feature2 = "nFeature_RNA"
)

p_scatter_after <- p_scatter_after_1 + p_scatter_after_2

print(p_scatter_after)

ggsave(
  filename = file.path(out_dir, paste0(sample_id, "_ScatterPlot_after_filtering.png")),
  plot = p_scatter_after,
  width = 10,
  height = 5,
  dpi = 300
)

#Normalization
donor_filtered <- NormalizeData(
  donor_filtered,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

#Variable feature 찾기
donor_filtered <- FindVariableFeatures(
  donor_filtered,
  selection.method = "vst",
  nfeatures = 2000
)
top10 <- head(VariableFeatures(donor_filtered), 10)
top10

#Variable feature plot
p_variable <- VariableFeaturePlot(donor_filtered)

p_variable_labeled <- LabelPoints(
  plot = p_variable,
  points = top10,
  repel = TRUE
)

print(p_variable_labeled)

ggsave(
  filename = file.path(out_dir, paste0(sample_id, "_VariableFeaturePlot.png")),
  plot = p_variable_labeled,
  width = 8,
  height = 6,
  dpi = 300
)

#Scaling
all_genes <- rownames(donor_filtered)

donor_filtered <- ScaleData(
  donor_filtered,
  features = all_genes
)

#PCA 실행
donor_filtered <- RunPCA(
  donor_filtered,
  features = VariableFeatures(object = donor_filtered),
  npcs = 50
)
print(donor_filtered[["pca"]], dims = 1:5, nfeatures = 5)

#ElbowPlot으로 PC 수 결정 !!변수있음!!
p_elbow <- ElbowPlot(
  donor_filtered,
  ndims = 50
)

print(p_elbow)

ggsave(
  filename = file.path(out_dir, paste0(sample_id, "_ElbowPlot.png")),
  plot = p_elbow,
  width = 7,
  height = 5,
  dpi = 300
)
#변수
dims_use <- 1:20

#Neighbor graph 생성
donor_filtered <- FindNeighbors(
  donor_filtered,
  dims = dims_use
)
names(donor_filtered@graphs)

#Clustering 기본 0.5로
donor_filtered <- FindClusters(
  donor_filtered,
  resolution = 0.5
)
table(Idents(donor_filtered))

#UMAP 실행
donor_filtered <- RunUMAP(
  donor_filtered,
  dims = dims_use
)

#UMAP 저장
p_umap_cluster <- DimPlot(
  donor_filtered,
  reduction = "umap",
  label = TRUE
)

print(p_umap_cluster)

ggsave(
  filename = file.path(out_dir, paste0(sample_id, "_UMAP_clusters_resolution_0.5.png")),
  plot = p_umap_cluster,
  width = 8,
  height = 6,
  dpi = 300
)

#UMAP 보고 Clustering resolution 수정 !!중요한 변수!!
donor_filtered <- FindClusters(
  donor_filtered,
  resolution = 0.3
)
table(Idents(donor_filtered))

#UMAP 저장
p_umap_cluster <- DimPlot(
  donor_filtered,
  reduction = "umap",
  label = TRUE
)

print(p_umap_cluster)

ggsave(
  filename = file.path(out_dir, paste0(sample_id, "_UMAP_clusters_resolution_change.png")),
  plot = p_umap_cluster,
  width = 8,
  height = 6,
  dpi = 300
)

#Cluster marker gene 찾기
markers <- FindAllMarkers(
  donor_filtered,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)
head(markers)

#Marker 저장
write.csv(
  markers,
  file.path(out_dir, paste0(sample_id, "_cluster_markers.csv")),
  row.names = FALSE
)

#PBMC marker 확인
pbmc_markers <- c(
  "CD3D", "CD3E", "TRAC",
  "IL7R", "CCR7",
  "CD8A", "CD8B",
  "NKG7", "GNLY",
  "MS4A1", "CD79A",
  "LYZ", "S100A8", "S100A9",
  "FCGR3A", "MS4A7",
  "FCER1A", "CST3",
  "PPBP", "PF4"
)

p_dot_cluster <- DotPlot(
  donor_filtered,
  features = pbmc_markers
) + RotatedAxis()

print(p_dot_cluster)

ggsave(
  filename = file.path(out_dir, paste0(sample_id, "_DotPlot_cluster_markers.png")),
  plot = p_dot_cluster,
  width = 14,
  height = 6,
  dpi = 300
)

#FeaturePlot으로 marker 위치 확인
#major cell type 확인용
p_feature_major <- FeaturePlot(
  donor_filtered,
  features = c("CD3D", "MS4A1", "LYZ", "NKG7"),
  reduction = "umap",
  ncol = 2
)

print(p_feature_major)

ggsave(
  filename = file.path(out_dir, paste0(sample_id, "_FeaturePlot_major_markers.png")),
  plot = p_feature_major,
  width = 10,
  height = 8,
  dpi = 300
)

#monocyte와 dendritic cell 확인용
p_feature_mono_dc <- FeaturePlot(
  donor_filtered,
  features = c("LYZ", "S100A8", "S100A9", "FCGR3A", "MS4A7", "FCER1A", "CST3"),
  reduction = "umap",
  ncol = 3
)

print(p_feature_mono_dc)

ggsave(
  filename = file.path(out_dir, paste0(sample_id, "_FeaturePlot_monocyte_DC_markers.png")),
  plot = p_feature_mono_dc,
  width = 10,
  height = 8,
  dpi = 300
)

#T cell 세부 확인용
p_feature_tcell <- FeaturePlot(
  donor_filtered,
  features = c("CD3D", "CD3E", "TRAC", "IL7R", "CCR7", "CD8A", "CD8B", "NKG7", "GNLY"),
  reduction = "umap",
  ncol = 3
)

print(p_feature_tcell)

ggsave(
  filename = file.path(out_dir, paste0(sample_id, "_FeaturePlot_Tcell_markers.png")),
  plot = p_feature_tcell,
  width = 10,
  height = 8,
  dpi = 300
)

#Cell type 이름 직접 붙이기 dotplot과 featureplot 참고 !!중요한 변수!!
levels(donor_filtered)

cluster_to_celltype <- c(
  "0"  = "T cell",
  "1"  = "CD14 Monocyte",
  "2"  = "T cell",
  "3"  = "T cell",
  "4"  = "T cell",
  "5"  = "CD8/Cytotoxic T cell",
  "6"  = "NK cell",
  "7"  = "FCGR3A Monocyte",
  "8"  = "T cell",
  "9"  = "B cell",
  "10" = "T cell",
  "11" = "Dendritic cell"
)

donor_filtered$celltype_manual <- cluster_to_celltype[
  as.character(Idents(donor_filtered))
]

table(donor_filtered$celltype_manual)

#Cell type UMAP보기
p_celltype <- DimPlot(
  donor_filtered,
  reduction = "umap",
  group.by = "celltype_manual",
  label = TRUE,
  repel = TRUE
)

p_celltype

#Cell type UMAP 저장
ggsave(
  filename = file.path(out_dir, paste0(sample_id, "_UMAP_celltype_manual.png")),
  plot = p_celltype,
  width = 8,
  height = 6,
  dpi = 300
)

#Cell type 비율 계산
celltype_count <- donor_filtered@meta.data %>%
  count(celltype_manual) %>%
  mutate(percent = n / sum(n) * 100)

celltype_count

write.csv(
  celltype_count,
  file.path(out_dir, paste0(sample_id, "_celltype_composition.csv")),
  row.names = FALSE
)

#seurat object로 저장
saveRDS(
  donor_filtered,
  file.path(out_dir, paste0(sample_id, "_Seurat_processed.rds"))
)
