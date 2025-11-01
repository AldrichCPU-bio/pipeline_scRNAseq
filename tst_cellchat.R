require(CellChat)
require(Seurat)
require(tidyverse)
require(ComplexHeatmap)
require(circlize)
require(paletteer)
require(future)
require(glue)
require(optparse)

options(future.globals.maxSize = 1.024e12)

option_list <- list(
  rdsfn = 
    optparse::make_option(c("--rdsPath", "-f"), type = "character", default = NULL, 
                          help  = "Path of seurat object after QC, `rds` format."), 
  working_directory = 
    optparse::make_option(c("--workingDirectory", "-w"), type = "character", default = NULL, 
                          help  = "Set working directory, if it not exist, 
                                   the wkdir would be created at first."), 
  group_input = 
    optparse::make_option(c("--groupInput", "-g"), type = "character", default = NULL, 
                          help  = "Which component shold be adopted for `group_by` argument? \nThis term must be one of the column names in meta data.")
  
)

args <- parse_args(OptionParser(option_list=option_list))

rds_fn0 <- args$rdsfn
wkdir <- args$working_directory
grp_fct <- args$group_input


seu1_mrg2 <- read_rds(rds_fn0)

wkdir %>% fs::dir_create() %>% setwd()

seu1_mrg2@meta.data <- 
  seu1_mrg2@meta.data %>% rownames_to_column("barcodes") %>% 
  dplyr::rename_with(.fn = function(x) str_replace(group_input, "col_grp")) %>% 
  as.data.frame() %>% column_to_rownames("barcodes")
seu1_mrg3 = seu1_mrg2
cel_idx <- seu1_mrg3@meta.data %>% as.data.frame() %>% 
  rownames_to_column("barcodes") %>% group_by(col_grp) %>% 
  dplyr::slice_sample(., prop = .15) %>% .$barcodes
# norm_data <- SeuratObject::LayerData(seu1_mrg3, layer = "data")
# meta_data <- seu1_mrg3@meta.data %>% tibble::rownames_to_column("samples")
# 'meta' data must have a column named `samples` & force it as a factor
colnames(seu1_mrg3@meta.data)[4] <- "samples"
cellchat <- 
  CellChat::createCellChat(
    object = seu1_mrg3[, cel_idx], group.by = "anot_cisplatin_resist")
levels(cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group

CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)

CellChatDB.use <- subsetDB(
  CellChatDB, search = "Secreted Signaling", key = "annotation") # use Secreted Signaling
cellchat@DB <- CellChatDB.use

# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multisession", workers = 8) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
#> The number of highly variable ligand-receptor pairs used for signaling inference is 692


# project gene expression data onto PPI (Optional: when running it, 
# USER should set `raw.use = FALSE` in the function `computeCommunProb()` 
# in order to use the projected data)
# cellchat <- projectData(cellchat, PPI.human)
future::plan("sequential")
ptm = Sys.time()
cellchat <- computeCommunProb(cellchat, type = "triMean")

cellchat <- filterCommunication(cellchat, min.cells = 10)

cellchat <- computeCommunProbPathway(cellchat)

cellchat <- aggregateNet(cellchat)
execution.time = Sys.time() - ptm
print(as.numeric(execution.time, units = "secs"))
write_rds(cellchat, "cellchat.rds", compress = "gz")

ptm = Sys.time()
groupSize <- as.numeric(table(cellchat@idents))
pdf("net_aggregation.pdf", width = 11, height = 5)
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= F, 
                 title.name = "Number of interactions") 
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize,
                 weight.scale = T, label.edge= F, 
                 title.name = "Interaction weights/strength") 
dev.off()

mat <- cellchat@net$weight
# par(mfrow = c(3,3), xpd=TRUE)
pdf("each_type_weight.pdf", width = 4, height = 4)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), 
                 dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, 
                   weight.scale = T, 
                   edge.weight.max = max(mat), 
                   title.name = rownames(mat)[i]) %>% print()
}
dev.off()

# Hierarchy plot
# Here we define `vertex.receive` so that 
# the left portion of the hierarchy plot shows 
# signaling to fibroblast and the right portion shows 
# signaling to immune cells 
glue("pathways") %>% fs::dir_create()
for(pathways.show in cellchat@netP$pathways) {
  pdf(glue("pathways/pathway_net_aggr_{pathways.show}.pdf"), 
      width = 5, height = 5)
  vertex.receiver = c(1, 5, 6, 8) # a numeric vector. 
  netVisual_aggregate(cellchat, signaling = pathways.show, 
                      pt.title = 8, layout = "hierarchy", 
                      vertex.receiver = vertex.receiver) 
  # Circle plot
  par(mfrow=c(1,1))
  netVisual_aggregate(cellchat, signaling = pathways.show, 
                      layout = "circle") 
  
  # netVisual_chord_cell(cellchat, signaling = pathways.show, 
  #                      group = group.cellType, 
  #                      title.name = 
  #                        paste0(pathways.show, " signaling network"))
  
  netAnalysis_contribution(cellchat, signaling = pathways.show)
  
  dev.off()
  
}
for(pathways.show in cellchat@netP$pathways) {
  pdf(glue("pathways/pathway_heatmap_{pathways.show}.pdf"), 
      width = 8, height = 5)
  # Heatmap
  par(mfrow=c(1,1))
  netVisual_heatmap(cellchat, signaling = pathways.show, 
                    color.heatmap = "Reds")
  #> Do heatmap based on a single object
}

# Access all the signaling pathways showing significant communications
pathways.show.all <- cellchat@netP$pathways
# check the order of cell identity to set suitable vertex.receiver
levels(cellchat@idents)
vertex.receiver = c(1, 5, 6, 8) # a numeric vector. 
for (i in 1:length(pathways.show.all)) {
  # Visualize communication network associated with both signaling pathway and individual L-R pairs
  pdf(glue("{pathways.show.all[i]}_net_vertex.pdf"), width = 8, height = 8)
  netVisual(
    cellchat, signaling = pathways.show.all[i], 
    vertex.receiver = vertex.receiver, layout = "hierarchy")
  dev.off()
  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
  gg <- netAnalysis_contribution(
    cellchat, signaling = pathways.show.all[i])
  ggsave(filename=paste0("pathways/", pathways.show.all[i], 
                         "_L-R_contribution.pdf"), 
         plot=gg, width = 3, height = 2, units = 'in', dpi = 300)
}

## modification needed ----
netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11), remove.isolate = FALSE)


