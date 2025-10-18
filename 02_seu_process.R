require("Seurat")
require("tidyverse")
require("glue")
require("SeuratObject")
require("optparse")
require("scuttle")
require("scRNAtoolVis")
require("ggsci")
require("paletteer")
require("ComplexHeatmap")
require("circlize")
require("ggpubr")
require("ggrepel")
require("writexl")
require("viridis")

## define my plot theme in ggplot2 ----
my_theme1 <- 
  ggplot2::theme_classic(base_size = 8) + 
  theme(legend.key.size = unit(3, "mm"), axis.text = element_text(size = 7), 
        axis.line = element_blank(), legend.background = element_blank(), 
        panel.border = element_rect(color = "black", linetype = "solid", 
                                    linewidth = .5))
## argments -----
option_list <- list(
  rdsfn = 
    optparse::make_option(c("--rdsPath", "-rds_fn"), type = "character", default = NULL, 
                          help  = "Path of seurat object after QC, `rds` format."), 
  working_directory = 
    optparse::make_option(c("--workingDirectory", "-wkdir"), type = "character", default = NULL, 
                          help  = "Set working directory, if it not exist, 
                                   the wkdir would be created at first."), 
  meta_info_path = 
    optparse::make_option(c("--sampleInfoPath", "-meta_fn"), type = "character", default = NULL, 
                          help  = "Set path of excel file, which contains necessary info for each sample.
                                   \n CAUTION 1: it must exist before running this code.
                                   \n CAUTION 2: `sample_names` must exist in column names, describing data file name."), 
  species_input = 
    optparse::make_option(c("--speciesInput", "-species"), type = "character", default = NULL, 
                          help  = "Which species did the data describe?"), 
  split_input = 
    optparse::make_option(c("--splitInput", "-split_fct"), type = "character", default = NULL, 
                          help  = "Which component shold the data be split by? \nThis term must be one of the column names in meta data."), 
  integration_method = 
    optparse::make_option(c("--integrateMethod", "-intg_method"), type = "character", default = NULL, 
                          help  = "Which integration method shold be selected? \nOnly could be one of 'harmony', 'rpca' or 'cca'."), 
  group_input = 
    optparse::make_option(c("--groupInput", "-grp_fct"), type = "character", default = NULL, 
                          help  = "Which component shold be adopted for `group_by` argument? \nThis term must be one of the column names in meta data."), 
  label_feat_csv = 
    optparse::make_option(c("--labelFeatureCSV", "-lab_csv"), type = "character", default = NULL, 
                          help  = "Which labeling features be adopted for `dotplot` or `featureplot`?"), 
  column4htp = 
    optparse::make_option(c("--columnForHeatmap", "-column4htp"), type = "character", default = "RNA_snn_res.0.2", 
                          help  = "Which component shold be adopted for heatmap column? \nThis term must be one of the column names in meta data.")
  
)

args <- parse_args(OptionParser(option_list=option_list))

rds_fn0 <- args$rdsfn
wkdir <- args$working_directory
datdir <- args$data_directory
meta_fn <- args$meta_info_path
species <- args$species_input
split_fct <- args$split_input
intg_method <- args$integration_method
grp_fct <- args$group_input
label_feat_csv <- args$label_feat_csv
column4htp <- args$column4htp

seu_mrg <- read_rds(rds_fn0)

## step2: Seurat processing -----
seu_process_func <- function(seu_obj = NULL, reduction_method = "pca") {
  seu_obj1 <- seu_obj %>% Seurat::NormalizeData() %>% Seurat::ScaleData() %>% 
    Seurat::FindVariableFeatures() %>% Seurat::RunPCA(npcs = 50) %>% 
    Seurat::RunUMAP(dims = 1:20, reduction = reduction_method, reduction.name = "umap") %>% 
    Seurat::RunTSNE(dims = 1:10, reduction = reduction_method, reduction.name = "tsne") %>% 
    Seurat::FindNeighbors(dims = 1:20) %>% Seurat::FindClusters(resolution = c(.1, .2, .3, .4, .5, .6))
  return(seu_obj1)
}
seu_mrg <- seu_process_func(seu_obj = seu_mrg, reducttion_method = "pca")
dim_plst0 <- lapply(c("orig.ident", split_fct, grp_fct, paste0("RNA_snn_res.", c(.1, .2, .3, .4, .5, .6))), \(grp){
  Seurat::DimPlot(seu_mrg, group.by = grp, pt.size = .1, alpha = .6, 
                  label = TRUE, label.size = 2, repel = TRUE) + 
    my_theme1 + labs(title = grp) + ggplot2::coord_equal()
})
pdf("original_seu_mrg_dimplot0.pdf", width = 5, height = 5)
print(dim_plst0)
dev.off()

mks_lst0 <- list()
for(reso in c(.1, .2, .3, .4, .5, .6)) {
  Idents(seu_mrg) <- seu_mrg@meta.data[[glue("RNA_snn_res.{reso}")]]
  mks_lst[[glue("RNA_snn_res.{reso}")]] <- seu_mrg %>% 
    Seurat::FindAllMarkers() %>% tibble() %>% mutate(resolution = glue("RNA_snn_res.{reso}"))
}
writexl::write_excel(mks_lst0, "original_seu_mrg_mks_lst.xlsx")

if(!is.null(intg_method)) {
  reduct_name_new <- case_when(
    intg_method %in% c("harmony", "hmy") ~ "intg_hmy", 
    intg_method == "rpca" ~ "intg_rpca", 
    intg_method == "cca" ~ "intg_cca", 
    TRUE ~ stop("This method, `{intg_method}` is not supported by Seurat, please check.")
  )
}

if(!is.null(split_fct)) {
  seu_mrg <- Seurat::SplitObject(seu_mrg, split.by = split_fct)
  if(is.null(intg_method)) {
    stop("You have confirm the object should be split, but none of the integration methods was set. Please check out.")
  } else {
    if (intg_method == "harmony") {
      seu_mrg <- Seurat::IntegrateLayers(method = HarmonyIntegration, orig.reduction = "pca", new.reduction = reduct_name_new)
  } else if (intg_method == "rpca") {
      seu_mrg <- Seurat::IntegrateLayers(method = RPCAIntegration, orig.reduction = "pca", new.reduction = reduct_name_new)
  } else if (intg_method == "cca") {
      seu_mrg <- Seurat::IntegrateLayers(method = CCAIntegration, orig.reduction = "pca", new.reduction = reduct_name_new)
  }
    seu_mrg <- seu_mrg %>% 
      Seurat::RunUMAP(reduction = reduct_name_new, dims = 1:30, reduction.name = "umap") %>% 
      Seurat::RunTSNE(reduction = reduct_name_new, dims = 1:10, reduction.name = "tsne") %>% 
      Seurat::FindNeighbors(dims = 1:20, reduction = reduct_name_new) %>% 
      Seurat::FindClusters(resolution = c(.1, .2, .3, .4, .5, .6))
    seu_mrg <- JoinLayers(seu_mrg)
    dim_plst0 <- lapply(c("orig.ident", split_fct, grp_fct, paste0("RNA_snn_res.", c(.1, .2, .3, .4, .5, .6))), \(grp){
      Seurat::DimPlot(seu_mrg, group.by = grp, pt.size = .1, alpha = .6, 
                      label = TRUE, label.size = 2, repel = TRUE, reduction = "umap") + 
        my_theme1 + labs(title = grp) + ggplot2::coord_equal()
    })
    pdf("integrat_seu_mrg_dimplot0.pdf", width = 5, height = 5)
    print(dim_plst0)
    dev.off()
    
    mks_lst0 <- list()
    for(reso in c(.1, .2, .3, .4, .5, .6)) {
      Idents(seu_mrg) <- seu_mrg@meta.data[[glue("RNA_snn_res.{reso}")]]
      mks_lst[[glue("RNA_snn_res.{reso}")]] <- seu_mrg %>% 
        Seurat::FindAllMarkers() %>% tibble() %>% mutate(resolution = glue("RNA_snn_res.{reso}"), integration_method = intg_method)
    }
    writexl::write_excel(mks_lst0, "integrat_seu_mrg_mks_lst.xlsx")
    
  }
} 

## dotplot & featureplot for labeling features ----
label_feat_tbl <- readr::read_csv(label_feat_csv)
lab_feat_lst <- lapply(unique(label_feat_tbl[["type"]]), \(tp) {
  label_feat_tbl %>% dplyr::filter(type == tp) %>% .[["feature"]]
}) %>% setNames(nm = unique(label_feat_tbl[["type"]]))

dot_plst1 <- lapply(c(.1, .2, .3, .4, .5, .6), \(reso) {
  Seurat::DotPlot(seu_mrg, features = lab_feat_lst, group.by = glue("RNA_snn_res.{reso}"), 
                  dot.scale = 3, assay = "RNA") + 
    my_theme1 + Seurat::RotatedAxis() + labs(title = glue("RNA_snn_res.{reso}")) + 
    viridis::scale_color_viridis() + labs(x = "", y = "") + 
    theme(legend.position = "bottom", 
          plot.margin = margin(l = -0.2, t = .2, r = .1, b = .2, unit = "cm"), 
          legend.margin = margin(t = -0.5, unit = "cm"), axis.text.x = element_text(size = 7))
})
pdf("lab_feat_seu_mrg_dotplot.pdf", width = 7, height = 4)
print(dot_plst1)
dev.off()
feat_plst1 <- lapply(unlist(lab_feat_lst), \(feat) {
  Seurat::FeaturePlot(seu_mrg, features = nm, pt.size = .1, alpha = .6, order = T, 
                      reduction = "umap") + 
    my_theme1 + ggplot2::coord_equal() + viridis::scale_color_viridis() + 
    labs(color = "", title = nm)
})
pdf("lab_feat_seu_mrg_featureplot.pdf", width = 5, height = 5)
print(feat_plst1)
dev.off()

## heatmap, manually -----
seu_mrg$seurat_cluster <- seu_mrg@meta.data[[column4htp]] %>% as.character() %>% fct()
Idents(seu_mrg) <- seu_mrg$seurat_cluster
mks4plot <- Serurat::FindAllMarkers(seu_mrg)

htp_top_mks <- mks4plot %>% dplyr::filter(p_val_adj < .2, pct.1 > .1, avg_log2FC > 0) %>% 
  dplyr::group_by(cluster) %>% dplyr::slice_head(n = 100)
meta_sel <- seu_mrg@meta.data %>% rownames_to_column("barcode") %>% dplyr::group_by(seurat_cluster) %>% 
  dplyr::slice_sample(n = 50) %>% ungroup() %>% dplyr::arrange(seurat_cluster)
dat_mtx <- LayerData(seu_mrg[["RNA"]], layer = "data") %>% .[unique(htp_top_mks$gene), meta_sel$barcode] %>% 
  as.matrix() %>% t() %>% scale() %>% t()
col_fun1 <- circlize::colorRamp2(breaks = c(-3, 0, 3), 
                                 colors = c("lightblue", "white", "darkred"))
anot_top1 <- ComplexHeatmap::HeatmapAnnotation(type = meta_sel$seurat_cluster, show_legend = F, 
                                               show_annotation_name = F, simple_anno_size = unit(4, "mm"), 
                                               col = list(type = col_seurat_cluster))
ord1 <- which(rownames(dat_mtx) %in% unlist(lab_lst))
anot_right <- ComplexHeatmap::rowAnnotation(
  foo = anno_mark(at = ord1, labels = rownames(dat_mtx)[ord1],
                  labels_gp = gpar(fontsize = 6)) 
)
htp1 <- ComplexHeatmap::Heatmap(dat_mtx, top_annotation = anot_top1, name = "Z-score", 
                                right_annotation = anot_right, border = T, 
                                column_gap = unit(0, "mm"), column_split = meta_sel$seurat_cluster, 
                                column_title_rot = 45, column_title_gp = gpar(fontsize = 6), 
                                cluster_rows = F, cluster_columns = F, 
                                show_row_dend = F, show_row_names = F, 
                                width = unit(6, "cm"), height = unit(8, "cm"), raster_by_magick = T, 
                                show_column_names = F, col = col_fun1)
pdf(glue("seu_mrg_top100_mks_htp_{column4htp}.pdf"), width = 4, height = 4)
print(htp1)
dev.off()


## scRNAtoolVis for marker list -----
mk_vol1 <- markerVocalno(markers = mks4plot,
              topn = 5, base_size = 7,
              labelCol = ggsci::pal_d3("category20b"))
ggsave(glue("seu_mrg_mk_vol1_{column4htp}.pdf"), mk_vol1, width = 8, height = 3)
mk_vol2 <- jjVolcano(diffData = mks4plot, base_size = 7, log2FC.cutoff = 0.5, 
                     col.type = "adjustP", topGeneN = 3, tile.col = ggsci::pal_d3("category20b"))
ggsave(glue("seu_mrg_mk_jjvol2_{column4htp}.pdf"), mk_vol2, width = 8, height = 3)

## save this Seurat object -----
write_rds(seu_mrg, "processed_seu_mrg.rds", compress = "gz")


