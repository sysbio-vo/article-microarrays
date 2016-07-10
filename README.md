## Assessment of alternative pipelines for cross-platform microarray data integration using RNA-seq data

To reproduce the results you will need the following structure of directories:

```
article-microarrays
│   README.md
└───allsamples_exprs
└───combined
└───combined_exprs
└───exprs
└───general
└───pdata
└───plots
└───preprocessed
└───raws
|   └───GSE19615
|   └───GSE37614
|   └───GSE58644
|   └───GSE60785
|   └───GSE65194
└───rnaseq
└───RNA-seq
└───scores
└───scripts
```

#### RNA-seq data

Use scripts in `RNA-seq` folder in order to download and aseemble the data (change scripts in case you use some task manager etc), then use `rnaseq_assemble.R` to combine and average across samples. Alternatively, use already processed data in  `rnaseq` folder.

#### Microarray data

Download raw Affymetrix and non-normalized Illumina data into corresponding subfolders in `raw` folder. Use `preprocessing_affymetrix.R` and `preprocessing_illumina.R` scripts to pre-process and normalize the data, instead use zipped datasets in `preprocessed` folder.

If you wish to test tranformation and normalization methods for Illumina raw bead-level data, use `raw_illumina_preprocessing.R` script and `raw/047_20150203_Tchou_CAFs.zip` file, which you should put into `raw/GSE37614` folder. Alternatively, you can create plot using already generated files for this dataset in folder `exprs`. See the script for details.
