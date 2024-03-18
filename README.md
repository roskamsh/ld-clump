# ld-clump

This is a nextflow pipeline which takes as input data from eQTLGen, BGEN files for a cohort-of-interest, and runs LD clumping for a set of variants associated with a specific Gene Symbol. The purpose of this is to generate a final list of eQTLs that are independent. This is necessary due to the fact that the input from eQTLGen are cis-eQTLs which, by design, are close together and likely to be in LD with one another.

## Setup

### Software requirements
To run this pipeline, you must have nextflow installed. Please see build instructions here: https://www.nextflow.io/docs/latest/getstarted.html#installation

Additionally, you must have singularity installed or loaded as a module, as this pipeline uses singularity to pull from DockerHub for the software required at each step in the pipeline. On Eddie, this can be done with the following:

```
module load roslin/singularity/3.10.1
```

Otherwise, please see: https://docs.sylabs.io/guides/3.0/user-guide/installation.html#installation

### Configuration

The configuration for a given run is specified in a run-specific config file, which you create and is dataset-specific. An example is shown in `genomicc.config`. The following parameters are required:

* TFs : a CSV file which contains two columns, the first being `TF`, which contains the GeneSymbol for a given transcription factor of interest. The second column is `CHR` and contains the chromosome number where that gene is positioned on the genome.
* SNPs : a CSV file which contains all SNPs you would like to run LD clumping on, which are associated with a given transcription factor. This should be formatted the same as the output from eQTLGen (See: https://molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz). This can be a "master list" and doesn't necessarily need to be filtered for only the genes you are interested in.
* BGEN_FILES : absolute path to the BGEN files for your cohort-of-interest. This should be specified using brace expansion for each chromosome as well as each file type (bgen, bgen.bgi, sample). 

## Run

To run this pipeline, just run the following command:

```
nextflow run main.nf -profile eddie -c path/to/custom.config 
```
