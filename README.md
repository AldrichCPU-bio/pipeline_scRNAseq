# pipeline_scRNAseq
this is the collection of how to perform scRNAseq analysis, mostly based on `Seurat`


## Content:

  **01_seu_qc.R**: import data from `cellranger` output, QC with percentage of mitochondria, number of feature and doublet finding from `scDblFinder`, using pre-defined threshold and 3 MAD range.
  
  **02_seu_process.R**: processing with `Seurat` package, including normalization, scaling, variables finding, dimensional reduction, neighbor finding and unsupervised clustering. markers under each resolution also calculated, and heatmap, feature plot and dot plot with user provided markers were allowed. 
  
  **03_SingleR_anot.R**: automatic annotation with `SingleR` package. `cloupe` file also produced with `loupeR` package.


## Individual usage: 

```{bash}
Rscript 01_seu_qc.R --help
  Usage: scrna_step1.R [options]
  
  Options:
          -w WORKINGDIRECTORY, --workingDirectory=WORKINGDIRECTORY
                  Set working directory, if it not exist,
                                     the wkdir would be created at first.
  
          -d DATADIRECTORY, --dataDirectory=DATADIRECTORY
                  Set data directory, feature names, barcodes & gene matrix would be searched in this path.
  
   CAUTION: it must exist before running this code.
  
          -m SAMPLEINFOPATH, --sampleInfoPath=SAMPLEINFOPATH
                  Set path of excel file, which contains necessary info for each sample.
  
   CAUTION 1: it must exist before running this code.
  
   CAUTION 2: `sample_names` must exist in column names, describing data file name.
  
          -s SPECIESINPUT, --speciesInput=SPECIESINPUT
                  Which species did the data describe?
  
          --mitochondrialGenePattern=MITOCHONDRIALGENEPATTERN
                  What the Mito gene's pattern?
  
          --ribosomeGenePattern=RIBOSOMEGENEPATTERN
                  What the ribosome gene's pattern?
  
          -h, --help
                  Show this help message and exit
```

