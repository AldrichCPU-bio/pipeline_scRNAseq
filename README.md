# pipeline_scRNAseq
this is the collection of how to perform scRNAseq analysis, mostly based on `Seurat`
\n
Content:
  01_seu_qc.R: import data from `cellranger` output, QC with percentage of mitochondria, number of feature and doublet finding from `scDblFinder`, using pre-defined threshold and 3 MAD range.
  
  02_seu_process.R: processing with `Seurat` package, including normalization, scaling, variables finding, dimensional reduction, neighbor finding and unsupervised clustering. markers under each resolution also calculated, and heatmap, feature plot and dot plot with user provided markers were allowed. 
  
  03_SingleR_anot.R: automatic annotation with `SingleR` package.
