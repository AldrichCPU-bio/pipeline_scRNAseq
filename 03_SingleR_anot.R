require("SingeR")
require("Seurat")
require("optparse")
require("tidyverse")
require("scater")
require("scuttle")
require("celldex")
require("BiocParallel")
require("LoupeR")
require("hdf5r")

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
  reffn = 
    optparse::make_option(c("--refScePath", "-ref_fn"), type = "character", default = NULL, 
                          help  = "Path of reference to applied in `SingleR`, 
                                   \nwhich must be `SingleCellExperiment` class and `rds` file format."), 
  clst_anot = 
      optparse::make_option(c("--clusterNeedAnnotat", "-clst_anot"), type = "character", default = "RNA_snn_res.0.2",
                            help  = "Which component shold be adopted for `SingleR` cluster mode? \nThis term must be one of the column names in meta data."), 
  species_input =
    optparse::make_option(c("--speciesInput", "-species"), type = "character", default = NULL,
                          help  = "Which species did the data describe?"), 
  label_feat_csv =
    optparse::make_option(c("--labelFeatureCSV", "-lab_csv"), type = "character", default = NULL,
                          help  = "Which labeling features be adopted for `dotplot` or `featureplot`?")
)

args <- parse_args(OptionParser(option_list=option_list))

rds_fn0 <- args$rdsfn
wkdir <- args$working_directory
meta_fn <- args$meta_info_path
ref_fn <- args$reffn
clst_anot <- args$clst_anot
species <- args$species_input
label_feat_csv <- args$label_feat_csv

wkdir %>% fs::dir_create() %>% setwd()
seu_mrg <- read_rds(rds_fn0)

if(!is.null(ref_fn)) {
  if(file.exists(ref_fn)) {
    ref_sce <- read_rds(ref_fn)
  } else {
    stop("The reference file path not existed, please check again.")
  }
} else if (species == "human") {
  print("`HumanPrimaryCellAtlasData` from `celldex` package would be the reference")
  ref_sce <- celldex::HumanPrimaryCellAtlasData(ensembl=FALSE)
} else if (species == "mouse") {
  print("`MouseRNAseqData` from `celldex` package would be the reference")
  ref_sce <- celldex::MouseRNAseqData(ensembl=FALSE)
} else {
  stop("You have list your species is: ", species, ", but none of the reference found, please check again.")
} 

## SingleR usage for seurat object annotation -----
### first, convert to `SingleCellExperiment` object
sce <- seu_mrg %>% Seurat::as.SingleCellExperiment() %>% scater::logNormCounts()
### second, perform SingleR annotation with reference SCE object 
if(grepl("ref_sce", ls())) {
  print("The reference of `SingleR` has been loaded, automatic annotation processing start.")
  predictions <- SingleR::SingleR(sce, cluster = seu_mrg@meta.data[[clst_anot]], 
                   ref = ref_sce, labels = ref_sce$label.main, de.method="wilcox", 
                   BPPARAM = BiocParallel::MulticoreParam(8))
}
write_rds(predictions, "SingleR_predictions.rds")
### transfer the predicted label to Seurat object -----
seu_mrg$SingleR_predict <- predictions$labels

seu_mrg$seurat_cluster <- seu_mrg@meta.data[["SingleR_predict"]] %>% as.character() %>% fct()
Idents(seu_mrg) <- seu_mrg$seurat_cluster
mks4plot <- Serurat::FindAllMarkers(seu_mrg)
write_csv(mks4plot, "mks_SingleR_predict.csv")

## checking SingleR annotation with UMAP -----
dim1 <- 
  Seurat::DimPlot(seu_mrg, group.by = "SingleR_predict", pt.size = .1, alpha = .6, 
                  label = TRUE, label.size = 2, repel = TRUE) + 
    my_theme1 + labs(title = "SingleR_predict") + ggplot2::coord_equal()
ggsave("umap_SingleR_predict.pdf", dim1, width = 5, height = 5)

## checking SingleR annotation with provided markers ----
label_feat_tbl <- readr::read_csv(label_feat_csv)
lab_feat_lst <- lapply(unique(label_feat_tbl[["type"]]), \(tp) {
  label_feat_tbl %>% dplyr::filter(type == tp) %>% .[["feature"]]
}) %>% setNames(nm = unique(label_feat_tbl[["type"]]))

dot1 <- Seurat::DotPlot(seu_mrg, features = lab_feat_lst, 
                        group.by = "SingleR_predict", 
                        dot.scale = 3, assay = "RNA") + 
  my_theme1 + Seurat::RotatedAxis() + labs(title = "SingleR_predict") + 
  viridis::scale_color_viridis() + labs(x = "", y = "") + 
  theme(legend.position = "bottom", 
        plot.margin = margin(l = -0.2, t = .2, r = .1, b = .2, unit = "cm"), 
        legend.margin = margin(t = -0.5, unit = "cm"), 
        axis.text.x = element_text(size = 7))
ggsave("lab_feat_SingleR_anot_dotplot.pdf", width = 7, height = 4)

## save the Seurat object & convert to cloupe file ----
readr::write_rds(seu_mrg, "SingleR_anot_seu_mrg.rds", compress = "gz")
loupeR::create_loupe_from_seurat(
  obj = seu_mrg, output_dir = wkdir, force = TRUE, ## forced to overwrite
  output_name = "seu_obj_SingleR_anot"
)

