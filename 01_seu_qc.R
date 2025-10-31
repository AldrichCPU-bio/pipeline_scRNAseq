require("Seurat")
require("tidyverse")
require("glue")
require("SeuratObject")
require("optparse")
require("scDblFinder")
require("scuttle")

option_list <- list(
  working_directory = 
    optparse::make_option(c("--workingDirectory", "-w"), type = "character", default = NULL, 
                          help  = "Set working directory, if it not exist, 
                                   the wkdir would be created at first."), 
  data_directory = 
    optparse::make_option(c("--dataDirectory", "-d"), type = "character", default = NULL, 
                          help  = "Set data directory, feature names, barcodes & gene matrix would be searched in this path. 
                                   \n CAUTION: it must exist before running this code."), 
  meta_info_path = 
    optparse::make_option(c("--sampleInfoPath", "-m"), type = "character", default = NULL, 
                          help  = "Set path of excel file, which contains necessary info for each sample.
                                   \n CAUTION 1: it must exist before running this code.
                                   \n CAUTION 2: `sample_names` must exist in column names, describing data file name."), 
  species_input = 
    optparse::make_option(c("--speciesInput", "-s"), type = "character", default = NULL, 
                          help  = "Which species did the data describe?"), 
  mitochondrial_pattern = 
    optparse::make_option(c("--mitochondrialGenePattern"), type = "character", default = NULL, 
                          help  = "What the Mito gene's pattern?"), 
  ribosome_pattern = 
    optparse::make_option(c("--ribosomeGenePattern"), type = "character", default = NULL, 
                          help  = "What the Ribosome gene's pattern?")
)

args <- parse_args(OptionParser(option_list=option_list))

wkdir <- args$working_directory
datdir <- args$data_directory
meta_fn <- args$meta_info_path
species <- args$species_input
mt_parten <- args$mitochondrial_pattern
ribo_parten <- args$ribosome_pattern

## check point: did data directory input exist?
if (!file.exist(datdir)) {
  stop(
    "Data directory not exist. Please check out."
  )
}
if (!file.exist(meta_fn)) {
  stop(
    "Meta data/sample info not exist. Please check out."
  )
} else if (!grepl(".xlsx$", meta_fn) {
  stop(
    "Meta data/sample info file is not a `.xlsx` format. Please check out."
  )
}

## check point: is there any defined gene pattern?
if (species == "human") {
  mt_parten <- "^MT-"
  ribo_parten <- "^RP(L|S)"
} else if (species == "mouse") {
  mt_parten <- "^mt-"
  ribo_parten <- "^Rp(l|s)"
} else {
  if(any(is.null(mt_parten), is.null(ribo_parten))) {
    stop(
      "species you intend to analyse is, ", species, 
      ", \nbut mitochondrial or ribosome genes' pattern are blank, \nplease check.")
    )
  }
  print("species you intend to analyse is, ", species, 
        ", \nand mitochondrial & ribosome genes' pattern are prepared. \nGood luck.")
}

## basic info checking and importing
wkdir %>% fs::dir_create() %>% setwd()
sample_info <- readxl::read_excel(meta_fn)
print(sample_info)

## define my plot theme in ggplot2 ----
my_theme1 <- 
  ggplot2::theme_classic(base_size = 8) + 
  theme(legend.key.size = unit(3, "mm"), axis.text = element_text(size = 7), 
        axis.line = element_blank(), legend.background = element_blank(), 
        panel.border = element_rect(color = "black", linetype = "solid", 
                                    linewidth = .5))

## files imported confirmation -----
sample_names <- sample_info$sample_names
mtx_fns <- glue("{datdir}/{sample_names}_matrix.mtx.gz")
feat_fns <- c(
  glue("{datdir}/{sample_names}_features.tsv.gz"), 
  glue("{datdir}/{sample_names}_genes.tsv.gz")
)
barcode_fns <- glue("{datdir}/{sample_names}_barcodes.tsv.gz")
dat_fn_tbl <- tibble(
  sample_name = sample_names, 
  mtx_fn = mtx_fns[file.exists(mtx_fns)], 
  barcode_fn = barcode_fns[file.exists(barcode_fns)], 
  feat_fn = feat_fns[file.exists(feat_fns)]
) %>% dplyr::filter(rowSums(is.na(.)) == 0, rowSums(. == "") == 0)

## check point: is there any matched file path?
if(nrow(dat_fn_tbl) == 0) {
  stop("there's no availiable data path, please make sure your data path and sample names exist.")
}

## data importing -----
seu_lst0 <- lapply(dat_fn_tbl$sample_name, \(samp) {
  cnt <- Seurat::ReadMtx(
    mtx = dat_fn_tbl %>% dplyr::filter(sample_name == samp) %>% .$mtx_fn, 
    cells = dat_fn_tbl %>% dplyr::filter(sample_name == samp) %>% .$barcode_fn, 
    features = dat_fn_tbl %>% dplyr::filter(sample_name == samp) %>% .$feat_fn
  )
  seu0 <- SeuratObject::CreateSeuratObject(
    counts = cnt, assay = "RNA", project = samp, 
    min.cells = 10, min.features = 3
  )
  meta_data <- dplyr::filter(sample_info, sample_names == samp)
  seu0@meta.data <- 
    seu0@meta.data %>% tibble::rownames_to_column("barcode") %>% 
    dplyr::left_join(y = meta_data, by = c("orig.ident" = "sample_names")) %>% as.data.frame() %>% 
    column_to_rownames("barcode") %>% 
    mutate(orig.ident = as.character(samp))
  return(seu0)
})

## join the list & basic statistics -----
if(nrow(dat_fn_tbl) == 1) {
  seu_mrg <- seu_lst0[[1]]
} else {
  seu_mrg <- merge(seu_lst0[[1]], seu_lst0[-1])
  seu_mrg <- SeuratObject::JoinLayers(seu_mrg)
}
seu_mrg$pct_mito <- Seurat::PercentageFeatureSet(
  seu_mrg, pattern = mt_parten
)
seu_mrg$pct_ribo <- Seurat::PercentageFeatureSet(
  seu_mrg, pattern = ribo_parten
)
sce0 <- seu_mrg %>% Seurat::as.SingleCellExperiment() %>% 
  scuttle::logNormCounts() %>% 
  scDblFinder::scDblFinder(samples = seu_mrg$orig.ident)

seu_mrg$scDblFinder.class <- sce0$scDblFinder.class
seu_mrg$scDblFinder.score <- sce0$scDblFinder.score

## plot raw QC parameters: 
##  nCounts, nFeatures, pct_mito, pct_ribo & scDblFinder.score
qc_feat <- c("nCount_RNA", "nFeature_RNA", "pct_mito", 
             "pct_ribo", "scDblFinder.score")
qc_vln1 <- Seurat::VlnPlot(seu_mrg, group.by = qc_feat, ncol = 2,  
                          features = feat, pt.size = 0, log = FALSE) & 
  my_theme1 & Seurat::RotatedAxis() & Seurat::NoLegend() & 
  labs(x = "") & theme(plot.title = element_text(size = 8))
ggsave("qc_vln1.pdf", qc_vln1, width = 0.4 * nrow(dat_fn_tbl), 
       height = 7, width = nrow(dat_fn_tbl) * 1.6)
qc_vln11 <- Seurat::VlnPlot(seu_mrg, group.by = qc_feat, ncol = 2,  
                           features = feat, pt.size = 0, log = TRUE) & 
  my_theme1 & Seurat::RotatedAxis() & Seurat::NoLegend() & 
  labs(x = "") & theme(plot.title = element_text(size = 8))
ggsave("qc_log_vln11.pdf", qc_vln1, width = 0.4 * nrow(dat_fn_tbl), 
       height = 7, width = nrow(dat_fn_tbl) * 1.6)

## QC filtering ----
### common threshold, of course not suitable for every situation
cel_idx1 <- (seu_mrg$pct_mito < 30) & (seu_mrg$nFeature_RNA > 100) & 
  (seu_mrg$scDblFinder.class == "singlet")

### according to OSCA book (bioconductor)
### 5MAD threshold would be adopted, for safety
mad_threhold_func <- function(x) {
  mad_val <- abs(x - median(x)) %>% median()
  thrd_hi <- median(x) + 3 * mad_val
  thrd_lo <- median(x) - 3 * mad_val
  return(c(thrd_hi, thrd_lo))
}
mad_call_func <- function(seu_obj = NULL, nm = NULL) {
  input <- seu_obj@meta.data[[nm]] %>% as.numeric()
  thrd_vals <- mad_threhold_func(input)
  print(paste0("QC parameter: ", nm))
  print(paste0("the high threshold is ", thrd_vals[1]))
  print(paste0("low threshold is ", thrd_vals[2]))
  return(thrd_vals)
}

cel_idx2 <- 
  (seu_mrg$pct_mito < mad_call_func(seu_mrg, "pct_mito")[1]) & 
  (seu_mrg$pct_mito > mad_call_func(seu_mrg, "pct_mito")[2]) & 
  (seu_mrg$nFeature_RNA < mad_call_func(seu_mrg, "nFeature_RNA")[1]) & 
  (seu_mrg$nFeature_RNA > mad_call_func(seu_mrg, "nFeature_RNA")[2]) 

seu_mrg$remain_general <- cel_idx1
seu_mrg$remain_3mad <- cel_idx2

seu_mrg_qc <- seu_mrg[, seu_mrg$remain_general & seu_mrg$remain_3mad] ## a tough threshold
## save the objects ------
rds_dir <- glue("{wkdir}/rds") %>% fs::dir_create()
write_rds(seu_mrg, glue("{rds_dir}/seu_mrg.rds"), compress = "gz")
write_rds(seu_mrg_qc, glue("{rds_dir}/seu_mrg_qc.rds"), compress = "gz")


