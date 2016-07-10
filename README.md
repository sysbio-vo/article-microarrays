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

Not all the scripts are meant to be launched in one go, read each script before executing. If you have troubles installing Brainarray CDFs manually, use `install.brainarray.R` script.

#### RNA-seq data

Use scripts in `RNA-seq/` folder in order to download and assemble the data (change scripts in case you use some task manager etc), then use `rnaseq_assemble.R` to combine and average across samples. Alternatively, use already processed data in  `rnaseq/` folder.

#### Microarray data

Download raw Affymetrix and non-normalized Illumina data into corresponding subfolders in `raw` folder. Use `preprocessing_affymetrix.R` and `preprocessing_illumina.R` scripts to pre-process and normalize the data, instead use zipped datasets in `preprocessed` folder.

Use `generate_common_genes_list.R` to generate common genes list needed in most of the scripts.

Use `combine.R` script to generate combined microarray datasets, which are needed for futher analysis. This script also generates intermediate files, needed for `arrays_comparison.R`.

If you wish to test tranformation and normalization methods for Illumina raw bead-level data, use `raw_illumina_preprocessing.R` script and `raw/047_20150203_Tchou_CAFs.zip` file, which you should put into `raw/GSE37614` folder. Alternatively, you can create plots using already generated files for this dataset in folder `exprs/`. See the script for details.

For array versus array comparison use `arrays_comparison.R` script, for studying many-to-one probesets groups variation use `variation_analysis.R`

#### Pipelines

Use `brainarray.R` and  `max_mean_random_scores.R` to process datasets accroding to each pipeline and average across samples. Alternatively, use files in `exprs/` folder instead.

Use `pipelines_analysis.R` to generate plots for datasets processed with different pipelines comparison, and `combined_data_analysis.R` to compared datasets, combined within each pipeline.

Feel free to create issues if you have troubles executing the scripts or installing libraries used in the scripts.
