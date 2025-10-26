# pipeline_scRNAseq
this is the collection of how to perform scRNAseq analysis, mostly based on `Seurat`


## Content:

  **01_seu_qc.R**: import data from `cellranger` output, QC with percentage of mitochondria, number of feature and doublet finding from `scDblFinder`, using pre-defined threshold and 3 MAD range.
  
  **02_seu_process.R**: processing with `Seurat` package, including normalization, scaling, variables finding, dimensional reduction, neighbor finding and unsupervised clustering. markers under each resolution also calculated, and heatmap, feature plot and dot plot with user provided markers were allowed. 
  
  **03_SingleR_anot.R**: automatic annotation with `SingleR` package. `cloupe` file also produced with `loupeR` package.


## Individual usage: 

Step1: QC, using `Seurat`

```
Rscript 01_seu_qc.R --help
```

```{bash}

  Usage: 01_seu_qc.R [options]
  
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

Step2: processing, dimension reducing, graph neighbors finding, unsupervised clustering

```
Rscript 02_seu_process.R --help
```

```{bash}

  Usage: 02_seu_process.R [options]


  Options:
          -f RDSPATH, --rdsPath=RDSPATH
                  Path of seurat object after QC, `rds` format.
  
          -w WORKINGDIRECTORY, --workingDirectory=WORKINGDIRECTORY
                  Set working directory, if it not exist,
                                     the wkdir would be created at first.
  
          -m SAMPLEINFOPATH, --sampleInfoPath=SAMPLEINFOPATH
                  Set path of excel file, which contains necessary info for each sample.
  
   CAUTION 1: it must exist before running this code.
  
   CAUTION 2: `sample_names` must exist in column names, describing data file name.
  
          -s SPECIESINPUT, --speciesInput=SPECIESINPUT
                  Which species did the data describe?
  
          -p SPLITINPUT, --splitInput=SPLITINPUT
                  Which component shold the data be split by?
  This term must be one of the column names in meta data.
  
          -i INTEGRATEMETHOD, --integrateMethod=INTEGRATEMETHOD
                  Which integration method shold be selected?
  Only could be one of 'harmony', 'rpca' or 'cca'.
  
          -g GROUPINPUT, --groupInput=GROUPINPUT
                  Which component shold be adopted for `group_by` argument?
  This term must be one of the column names in meta data.
  
          -l LABELFEATURECSV, --labelFeatureCSV=LABELFEATURECSV
                  Which labeling features be adopted for `dotplot` or `featureplot`?
  
          --columnForHeatmap=COLUMNFORHEATMAP
                  Which component shold be adopted for heatmap column?
  This term must be one of the column names in meta data.
  
          -h, --help
                  Show this help message and exit


```

Step3: cell annotation, with automatic method using `SingleR`

```
Rscript 03_SingleR_anot.R --helpRscript 03_SingleR_anot.R --help
```

```{bash}

  Usage: 03_SingleR_anot.R [options]
  
  
  Options:
          -f RDSPATH, --rdsPath=RDSPATH
                  Path of seurat object after QC, `rds` format.
  
          -w WORKINGDIRECTORY, --workingDirectory=WORKINGDIRECTORY
                  Set working directory, if it not exist,
                                     the wkdir would be created at first.
  
          -m SAMPLEINFOPATH, --sampleInfoPath=SAMPLEINFOPATH
                  Set path of excel file, which contains necessary info for each sample.
  
   CAUTION 1: it must exist before running this code.
  
   CAUTION 2: `sample_names` must exist in column names, describing data file name.
  
          --refScePath=REFSCEPATH
                  Path of reference to applied in `SingleR`,
  
  which must be `SingleCellExperiment` class and `rds` file format.
  
          -c CLUSTERNEEDANNOTAT, --clusterNeedAnnotat=CLUSTERNEEDANNOTAT
                  Which component shold be adopted for `SingleR` cluster mode?
  This term must be one of the column names in meta data.
  
          -s SPECIESINPUT, --speciesInput=SPECIESINPUT
                  Which species did the data describe?
  
          -l LABELFEATURECSV, --labelFeatureCSV=LABELFEATURECSV
                  Which labeling features be adopted for `dotplot` or `featureplot`?
  
          -h, --help
                  Show this help message and exit

```
