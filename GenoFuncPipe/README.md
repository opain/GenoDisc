# GenoFunc: Genomic Analysis of Functional Annotations Pipeline

This is a snakemake pipeline for running a range of genome-wide association study (GWAS) summary statistic-based analyses. Analyses can be split into three parts:

Part 1: Descriptive analyses

  * GWAS summary statistic quality control
  * Confounding and SNP-heritability estimation
  * Identification of indpendent associations

Part 2: Gene finding analyses

  * SNP-based finemapping
  * Positional mapping of associated variants
  * Inference of gene expression/protein levels associated with the GWAS phenotype

Part 3: Drug finding analyses

  * Identify enriched bipartite drug-gene sets
  * Identify drugs that reverse diseases gene expression
  * Identify enriched ATC classifications

The results of all analyses are summarised in an [.html report](https://opain.github.io/GenoFunc/example_report.html).

### Using the pipeline

Users can upload their GWAS summary statistics [here](https://opain.shinyapps.io/NEUROHACK_GenoFunc_demo/) for analysis via our King's College London server.

Alternatively, users can download the pipeline and run analyses [locally](#locally) (see instructions below).

***

<a name="locally"/>

## Running locally

### Getting started

#### Step 1: Clone the repository.

```bash
git clone https://github.com/opain/GenoFunc.git
cd GenoFunc/GenoFuncPipe
```

#### Step 2: Install [Anaconda](https://conda.io/en/latest/miniconda.html).

**Linux:**
```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
sh Miniconda3-latest-Linux-x86_64.sh
```

Install Python 3.8, Snakemake 5.32, and the basic project dependencies. Note. I am installing these packages in an environment called 'base', if you already have an environment called 'base', you may need to create a new environment to avoid conflicts.

```bash
conda activate base
conda install python=3.8
conda install -c conda-forge mamba
mamba install -c bioconda -c conda-forge snakemake-minimal==5.32.2
mamba install pandas
```

> Note. If ghostscript is not already installed on your system, you will need to install it. You can check this by typing 'ghostscript' into the terminal.

#### Step 3: Prepare a [snakemake profile](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles) for parallel computing.
I have provided an example for users using a slurm scheduler called 'example_slurm_profile_config.yaml'. Slurm users should create a folder called 'slurm' in '$HOME/.config/snakemake', and then copy in the example_slurm_profile_config.yaml, renaming it to config.yaml.

### Running pipeline using test data

#### Step 4: Download test data.

```bash
wget -O test_data.tar.gz https://zenodo.org/record/6093584/files/test_data.tar.gz?download=1
tar -xf test_data.tar.gz
rm test_data.tar.gz
```

#### Step 5: Run pipeline for Coronary Artery Disease GWAS (COAD01).

```bash
snakemake --profile slurm --use-conda results/COAD01/reports/COAD01_report.html
```

> Note. If you receive an error saying 'MissingOutputException', you should try adding '--latency-wait 20' to the snakemake command, which tells the pipeline to wait 20 seconds between steps, thereby allowing filesystem latency.

> Note. Please be patient when running the pipeline for the first time. Expect the 'downloading and installing remote packages' to take ~1 hour. It has to create the conda environment in first instance, which involves installing python and R and many packages. Expect this to take ~1 hour.

### Running pipeline using your own data

You must specify a [gwas_list file](#gwas_list_file) listing GWAS summary statistics for the pipeline to use.

In addition, some external datasets cannot be downloaded automatically due to data restrictions. If you would like to infer altered protein levels associated with the GWAS phenotype using ROSMAP or Banner datasets, these must be downloaded in advance from [here](https://www.synapse.org/#!Synapse:syn23627957).

The location of those files must be specified in the [config.yaml](config.yaml) file. 

<a name="gwas_list_file">
 
#### gwas_list file format

- name: Short name for the GWAS
- path: File path to the GWAS summary statistics (uncompressed or gzipped)
- population: The super population that the GWAS was performed in (AFR/AMR/EAS/EUR/SAS)
- sampling: The proportion of the GWAS sample that were cases (if binary, otherwise NA)
- prevalence: The population prevelance of the phenotype (if binary, otherwise NA)
- mean: The phenotype mean in the general population (if continuous, otherwise NA)
- sd: The phenotype sd in the general population (if continuous, otherwise NA)
- label: A human readable name for the GWAS phenotype. Wrap in quotes if multiple words. For example, "Body Mass Index".

> Note. This file should be space delimited.

#### GWAS summary statistic format

The following column names are expected in the GWAS summary statistics files:

- SNP: RSID for variant (either SNP or CHR and BP required)
- CHR: Chromosome number (either SNP or CHR and BP required)
- ORIGBP: Base pair position (either SNP or CHR and BP required)
- A1: Allele 1 (effect allele) (required)
- A2: Allele 2 (required)
- P: P-value of association (required)
- OR: Odds ratio effect size (required if binary)
- BETA: BETA effect size (required if continuous)
- SE: Standard error of log(OR) or BETA (required)
- N: Total sample size (required)
- FREQ: Allele frequency in GWAS sample (optional)
- INFO: Imputation quality (optional)

### Output

All the results for a given GWAS will be stored within the folder 'results/\<GWAS ID>'. An .html file summarising the results can be found in 'results/\<GWAS ID>/reports/\<GWAS ID>_report.html'.

### Running parts of the pipeline

By default, analyses using all methods and external datasets will be run. This behaviour can be changed by modifying the config.yaml file. For example, if you do not want to run LD score regression, you change 'ldsc: T' to 'ldsc: F'. 

***

### Troubleshooting

If using the --profile flag, the log files will be saved in the GenoPredPipe/logs folder. If running interactively (i.e. -j1), the error should be printed on the screen.
 
***

### Acknowledgements

This project was initiated by Team MND during the NEUROHACK 2022 hackathon. Our team won 1st prize and was awarded a pilot grant to continue the project. The funding was used to purchase the server which hosts GenoFunPipe for users to submit GWAS for analysis online.

Oliver Pain is supported by a Sir Henry Wellcome Postdoctoral Fellowship [222811/Z/21/Z]. 

The authors acknowledge use of the research computing facility at King’s College London, Rosalind (https://rosalind.kcl.ac.uk), which is delivered in partnership with the NIHR Maudsley BRC, and part-funded by capital equipment grants from the Maudsley Charity (award 980) and Guy’s & St. Thomas’ Charity (TR130505). The views expressed are those of the authors and not necessarily those of the NHS, the NIHR or the Department of Health and Social Care.

***

Please contact Oliver Pain (oliver.pain@kcl.ac.uk) for any questions or comments.

***
<img src="https://user-images.githubusercontent.com/82537630/161148867-bc2cd3c3-799a-4deb-84bc-3864bfba4f00.png" width="300" height="229" />&nbsp;<img src="https://user-images.githubusercontent.com/82537630/161148282-0eb5c7fd-03d2-4cc0-9724-2391c29a6a53.png" width="424" height="96" />


