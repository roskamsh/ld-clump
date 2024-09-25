# ld-clump

This is a nextflow pipeline which takes as input data from eQTLGen, BGEN files for a cohort-of-interest, and runs LD clumping for a set of variants associated with a specific Gene Symbol. The purpose of this is to generate a final list of eQTLs that are independent. This is necessary due to the fact that the input from eQTLGen are cis-eQTLs which, by design, are close together and likely to be in LD with one another.

## Setup

### Software requirements
To run this pipeline, you must have nextflow and singularity installed. To install these individually, please see https://www.nextflow.io/docs/latest/getstarted.html#installation and https://docs.sylabs.io/guides/3.0/user-guide/installation.html#installation. This can also be done through conda using the `env.yaml` file contained in this repository. Please install conda on your operating system following the instructions here: https://docs.anaconda.com/miniconda/#quick-command-line-install.

Once conda is installed, you can create the environment using the following command:

```
conda env create --file env.yaml
```

Once installed, please activate your environment so that it can be used to run the pipeline.

```
conda activate nextflow
```

### Configuration

The configuration for a given run is specified in a run-specific config file, which you create and is dataset-specific. An example is shown in `genomicc.config`. The following parameters are required:

* `BGEN_FILES [ required ]` : absolute path to the BGEN files for your cohort-of-interest. This should be specified using brace expansion for each chromosome as well as each file type (bgen, bgen.bgi, sample).
* ASSEMBLY [ required ] : either grch37 or grch38 (also will accept hg19 or hg38), which is used to define problem areas of the genome to exclude (See: https://github.com/gabraham/flashpca/blob/master/exclusion_regions_hg19.txt). 
* eQTLGEN_DATA [optional, default: 2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt from eQTLGen]: Downloaded data from eQTLGen, in the form of a tab-separated text file. This contains all SNPs you would like to run LD clumping on, which are associated with a given transcription factor. This should be formatted the same as the output from eQTLGen (See: https://molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz). The SNP IDs must be in the same format as the IDs in your BGEN_FILES. If this is not provided, the file will be automatically pulled from eQTLGen.
* R2_THRESHOLD [optional, default: 0.8] : R-squared threshold to use for idenftifying which SNPs are in LD with one another. Default is 0.8. 
* INFO_THRESHOLD [ optional, default: 0.9 ] : Info score threshold set on imputation quality. Anything below this value will not be considered for hard-calling genotypes from the BGEN files. The default is 0.9.
* MAF_THRESHOLD [ optional, default: 0.01 ] : Minor allele frequency below which SNPs will be filtered out from your final estimands file. 
* OUTDIR [ optional, default: output ] : Output directory for results.
* SNPSTATS_CACHE [ optional, default: OUTDIR/info_scores ] : Cache for saving snp stats output from QCtoolv2. This contains the output from the process bgen::generate_info_score(), and includes the following crucial information for every SNP in your BGEN_FILES: rsid, info, minor_allele_frequency.

## Run

### Run profiles

The following profiles are available in this pipeline:
* eddie
* ultra2

To run this pipeline, just run the following command:

```
nextflow run main.nf -profile your_profile -c path/to/custom.config -with-trace
```

## Docker containers

Docker containers for this workflow can be found either in the `containers` directory in this repository, or at https://github.com/TARGENE/ld-block-removal/containers

