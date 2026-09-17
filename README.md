# TF genomic analysis pipeline

This repository contains the scripts used to process and analyse genome-wide transcription factor (TF)–DNA binding data.

The pipeline was developed to reprocess publicly available *in vitro* and *in vivo* TF-binding datasets using a consistent workflow. It covers the main steps from read mapping and peak calling to binding-model construction, model evaluation, and the analysis of TF-binding-site (TFBS) syntax and genomic positioning.

The principal Bash functions are collected in `compil_functions.sh`. Their usage is documented in `compil_usage.sh`, and a step-by-step tutorial is available in the [tutorial directory](tutorial/README.md).

## Main features

The repository provides functions and scripts for:

- downloading sequencing data from the Sequence Read Archive (SRA);
- mapping single-end and paired-end FASTQ files with Bowtie2;
- calling peaks with MACS3;
- generating consensus peak sets and comparing biological replicates;
- comparing binding profiles between experiments or conditions;
- generating matched negative sequence sets;
- discovering and training TF-binding models, including PWM, TFFM, DNAshape, KSM, and SeqConv models;
- evaluating model performance using ROC curves;
- estimating genome-wide binding-score distributions and model thresholds;
- analysing homotypic TFBS spacing and orientation; and
- analysing TFBS positioning and orientation relative to transcription start sites (TSSs).

## Repository structure

The main components of the repository are:

- `compil_functions.sh`: Bash functions implementing the principal analysis steps;
- `compil_usage.sh`: help messages describing the arguments accepted by each function;
- `config.sh`: paths to programs, scripts, reference files, and Conda environments;
- `bin/`: scripts called by the main Bash functions;
- [`tutorial/README.md`](tutorial/README.md): detailed instructions and a complete worked example;
- `conda_TFgenomics.yml`: definition of the main Conda environment; and
- `conda_MACS3.yml`: definition of the environment used for MACS3.

## Requirements and installation

### 1. Download the repository

```bash
git clone https://github.com/Bioinfo-LPCV-RDF/TF_genomic_analysis_v2.git
cd TF_genomic_analysis_v2
```


### 2. Create the Conda environments

The pipeline uses two Conda environments. They can be created from the supplied YAML files using environment names:

```bash
conda env create -f conda_TFgenomics.yml --name tf-genomic-analysis
conda env create -f conda_MACS3.yml --name macs3
```

Alternatively, they can be installed at user-defined locations:

```bash
conda env create -f conda_TFgenomics.yml --prefix <path_to_main_environment>
conda env create -f conda_MACS3.yml --prefix <path_to_MACS3_environment>
```

Activate the main environment with either:

```bash
conda activate tf-genomic-analysis
```

or:

```bash
conda activate <path_to_main_environment>
```

### 3. Configure paths and external software

Some dependencies are not installed through the supplied Conda environments and must be installed separately. All required external programs and reference files are listed in `config.sh`.

Before running the pipeline, edit `config.sh` to define the paths appropriate for your system.

The pipeline was developed using the *A. thaliana* TAIR10 reference genome. Other genomes can be used when the corresponding FASTA sequence, chromosome-size file, annotation, and any required blacklist file are supplied.

## Usage

Load the functions and their help messages into the current Bash session:

```bash
source compil_functions.sh
source compil_usage.sh
```

Help for an individual function can then be displayed with `-h` or `--help`, for example:

```bash
main_mapping_Fastq -h
main_peakcalling --help
compute_motif -h
compute_ROCS -h
```

The general form of a command is:

```bash
function_name [arguments]
```

Several functions accept Bash arrays for lists of samples, files, matrices, thresholds, or colours. The detailed syntax and complete argument descriptions are provided in `compil_usage.sh` and in the tutorial.

## Typical workflow

A complete analysis does not necessarily require every function. A typical workflow consists of:

1. downloading or preparing the sequencing reads;
2. mapping reads to the reference genome;
3. calling peaks independently for each replicate;
4. constructing consensus peak sets and assessing replicate agreement;
5. generating positive and matched negative sequence sets;
6. training one or more TFBS models;
7. evaluating and comparing the models; and
8. analysing TFBS spacing, orientation, genomic positioning, or condition-specific binding.

Detailed commands and examples for each step are provided in the [tutorial](tutorial/README.md).

## Outputs and reproducibility

Output directories and filenames are created by the individual functions and depend on the selected analysis. Wherever applicable, the pipeline records intermediate files, model scores, graphical outputs, and logs containing the arguments used for the analysis.

To facilitate reproducibility, we recommend:

- keeping a copy of the completed `config.sh` file;
- recording the versions of external programs that are not managed through Conda;
- retaining the original Conda YAML files; and
- setting explicit random seeds when the relevant functions provide this option.

The code reflects the workflow used for the TransAt study. Although the scripts are made available to support transparency and reuse, they may require adaptation to local computing environments, scheduler configurations, file-naming conventions, or species other than *A. thaliana*.

## Associated data

The processed peak files, coverage tracks, TFBS models, and accompanying metadata generated for TransAt are available from [Zenodo](https://doi.org/10.5281/zenodo.20556976).


## Citation

If you use this pipeline, please cite:

> Jegou A., Lucas J., Dreuillet M., Parcy F. and Blanc-Mathieu R. *A comprehensive Arabidopsis transcription factor binding atlas reveals pervasive positional and syntactic organization of their DNA binding*.

## Licence

This project is licensed under the [Creative Commons Attribution 4.0 International licence](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0).

## Contact

For questions about the pipeline, please contact:

- Alice Jegou: [alice.jegou@cea.fr](mailto:alice.jegou@cea.fr)
- Jérémy Lucas: [jeremy.lucas@cea.fr](mailto:jeremy.lucas@cea.fr)
- Romain Blanc-Mathieu: [romain.blancmathieu@cea.fr](mailto:romain.blancmathieu@cea.fr)

