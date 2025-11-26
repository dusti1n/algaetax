# algaetax: taxonomic data query

**algaetax is a bioinformatics tool** for automated querying, extraction, and standardization of taxonomic data from multiple biological databases. It streamlines the retrieval of taxonomy information from sources such as NCBI, PR2, and AlgaeBase, and organizes the results into a unified, structured format suitable for downstream analyses. By ensuring consistency and reproducibility across heterogeneous data sources, algaetax facilitates efficient integration of taxonomy data into bioinformatics pipelines and ecological or molecular studies.

## Table of Contents
1. [Dependencies](#dependencies)  
2. [Using algaetax](#using-algaetax)  
3. [Install Conda](#install-conda)  
4. [Install Environment](#install-environment)  
5. [Configuration Setup](#configuration-setup)
6. [AlgaeBase API](#algaebase-api)
7. [Run algaetax](#run-algaetax)  
8. [References](#references)
9. [Acknowledgement](#acknowledgement)
10. [Troubleshooting](#troubleshooting)

---

## Dependencies

- [**Linux (recommended)**](https://ubuntu.com/)  
  The workflow was developed and tested on Ubuntu Linux. Linux provides a stable, high-performance environment for running computationally intensive bioinformatics pipelines and ensures compatibility with most scientific software.

- [**Conda**](https://conda.io/en/latest/index.html)  
  Cross-platform package and environment manager used to install all required software in isolated, reproducible environments. Conda ensures that the correct dependency versions are used and simplifies sharing of the computational setup.

- [**Python**](https://www.python.org/)  
  Interpreted programming language required as the foundation for many bioinformatics tools. Python provides scripting capabilities, modularity, and flexibility for integrating multiple workflow components.

- **Conda Environment Packages**  
  When creating the Conda environment for *algaetax*, specific packages are required. These dependencies are listed in the `algaetax.yaml` file and will be installed automatically when setting up the environment.

## Using algaetax

Install, configure, and run algaetax to retrieve and standardize taxonomic data from multiple databases. Set up the required environment, adjust configuration parameters to your dataset and computational resources, and execute the workflow to obtain reproducible results. Proper setup ensures efficient processing and consistent taxonomic data retrieval with algaetax.

### Example run of algaetax

In this example, algaetax was executed with a small biological test dataset to demonstrate its core functionality. The tool takes a list of organism or taxon names as input, applies a filtering and preprocessing step to clean and standardize the data, and then automatically queries several reference databases available within algaetax — including NCBI, P2, and AlgaeBase. The demonstration illustrates how algaetax retrieves valid and up-to-date taxonomic names as well as hierarchical classifications from these sources. This example run highlights the typical workflow of algaetax, from data preparation to automated taxonomic verification, using a compact test dataset for clarity.

<p><img src="documentation/images/algaetax_run.gif" alt="algaetax demo" width="700"></p>

<p><b>Fig. 1:</b> Example run of algaetax with a small biological test dataset.</p>

### Install Conda

#### Install **Miniconda3** on Linux

Install Miniconda to manage the packages and environments required by algaetax.
It provides a lightweight and reproducible setup that ensures all dependencies are installed correctly and isolated for a consistent workflow.

```bash
# Download the Miniconda3 installer
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
```

```bash
# Make the installer executable
chmod +x Miniconda3-latest-Linux-x86_64.sh
```

```bash
# Start the installation (accept the license during setup)
bash Miniconda3-latest-Linux-x86_64.sh
```

### Install Environment

Setting up a dedicated Conda environment ensures that algaetax runs with all required dependencies and compatible package versions. By installing the environment from the provided algaetax.yaml file, you can easily reproduce the exact setup used for development and testing. This isolated environment prevents version conflicts with other Python projects and provides a stable, controlled workspace for running the workflow reliably across different systems.

#### Create a working directory

Set up a dedicated folder for your **algaetax** project to keep files and results organized.

```bash
# Example: create and enter a working directory
mkdir -p ~/algaetax_project
cd ~/algaetax_project
```

#### Download **algaetax**  
Clone the repository from GitHub or download it as a ZIP archive. After downloading, make sure to navigate into the algaetax directory, as the following setup steps need to be executed from within this folder.

```bash
git clone https://github.com/dusti1n/algaetax.git
```

#### Create Conda environment  

Before creating the environment, navigate to the root directory of algaetax, where the `algaetax.yaml` file is located. This file includes all required dependencies and version specifications to ensure that the workflow runs smoothly and reproducibly. By creating the environment from this file, you set up a clean and consistent workspace for running algaetax without dependency conflicts.

```bash
# Create the Conda ENV
conda env create -f algaetax.yaml
```

```bash
# Activate the environment
conda activate algaetax
```

---

### Configuration Setup

Before running **algaetax**, you should review and adjust the main parameters in the `config.yaml` file located in the project root directory. This configuration file controls how **algaetax** processes your data — including input file paths, filter settings, selected databases, and API credentials. Open `config.yaml` with a text editor (e.g., VS Code, Nano) and modify the values to match your dataset structure, system setup, and analysis requirements. Proper configuration ensures that the workflow runs efficiently and produces consistent, reproducible results.

#### Configuration parameters  

| Parameter | Description | Example |
|------------|--------------|----------|
| `input_data` | Path to input Excel file | `input_data/dustin/input_list_trim_2.xlsx` |
| `output_dir` | Output folder for results | `results_algaetax` |
| `taxa_column_number` | Column number for taxa names (1-based) | `4` |
| `header_row` | Header row number (`false` if none) | `1` |
| `backup_config` | Save config copy to output | `true` |
| `backup_input` | Backup input file before run | `true` |
| `export_not_found_taxa_list` | Export list of unfound taxa | `true` |
| `filter.blacklist_file` | Path to blacklist file | `blacklist.txt` |
| `filter.backup_blacklist` | Backup blacklist to output | `true` |
| `filter.skiplist_file` | Path to skiplist file | `skiplist.txt` |
| `filter.backup_skiplist` | Backup skiplist to output | `true` |
| `database.NCBI` | Use NCBI database | `true` |
| `database.PR2` | Use local PR2 database | `true` |
| `database.ALGB` | Use AlgaeBase (API key required) | `true` |
| `database_path.db_pr2` | Path to local PR2 file | `database/pr2_version_5.1.0_taxonomy.xlsx` |
| `db_ncbi.api_key` | NCBI API key (optional) | `abc123xyz` |
| `db_ncbi.ncbi_email` | Email for NCBI queries | `anonymous@ncbi.com` |
| `db_algaebase.api_key` | AlgaeBase API key | `xyz987abc` |
| `db_algaebase.api_url` | AlgaeBase API endpoint | `https://api.algaebase.org/v1.3/species` |

<p><b>Table 1:</b> Main configuration parameters.</p>

## AlgaeBase API

If you want to use the AlgaeBase database, a valid API key is required. Without it, you will not be able to access or use the AlgaeBase DB. If you do not have an API key, set the `ALGB` option to `false` in the configuration file.

## Run algaetax

After activating the Conda environment and setting the configuration file, run algaetax to retrieve and standardize taxonomic data from multiple sources, ensuring consistent taxonomy and reproducible results.

### Start analysis  

#### Run the main script to start the algaetax workflow 

During this step, the tool reads the configuration file, loads the input data, and queries the selected databases (e.g., NCBI, PR2, AlgaeBase). It retrieves and standardizes the corresponding taxonomic information, then saves the results to the output directory. The process duration depends on dataset size and the number of databases used.

**Important:** Activate the **algaetax** Conda environment before running commands.

Run **algaetax** using the default `config.yaml` located in the root directory:

```bash
# Run algaetax using the default config.yaml in the root directory
python3 algaetax.py
```

You can also run algaetax with a custom configuration file:

```bash
# Run algaetax with a custom configuration file
python3 algaetax.py --configfile <path>/<config.yaml>
```

Example for a custom configuration file:

```bash
# Run algaetax with a specific configuration file
python3 algaetax.py --configfile config_presets/tax_config_ncbi_pr2_algb.yaml
```

### Outputfiles 

After the analysis is completed, **algaetax** automatically saves all generated files, logs, and backups in the directory defined under `output_dir` in your `config.yaml`. Each run creates a structured set of output files that document both the retrieved taxonomy data and the overall workflow execution.

| **Output type** | **Description** |
|------------------|-----------------|
| Taxonomy tables | CSV files with standardized taxonomic data. |
| Reports | Summaries of matched taxa and query stats. |
| Unmatched taxa (`taxa_not_found.csv`) | List of taxa without database matches. |
| Backups | Copies of config, input data, and blacklist. |
| Log files | Records of all steps, timestamps, and warnings. |

<p><b>Table 2:</b> Overview of main output files generated by algaetax.</p>

These outputs ensure that every analysis performed with **algaetax** is transparent, traceable, and fully reproducible, providing a clear record of all processed data, database interactions, and workflow steps for reliable documentation and future reference.

---

## References  

The following references include the main software libraries, environment tools, and taxonomic databases used by algaetax. They represent the core scientific and technical foundations of the workflow, supporting data processing, taxonomy retrieval, and reproducibility across analyses.

#### Environment 

- Van Rossum, G., & Drake, F. L. (2009). *Python 3 Reference Manual.* CreateSpace, Scotts Valley, CA.  
- Anaconda, Inc. (2012). *Conda: Package, dependency and environment management for any language.* [https://docs.conda.io](https://docs.conda.io)

#### Software

- Cock, P.J.A., et al. (2009). *Biopython: Freely available Python tools for computational molecular biology and bioinformatics.* *Bioinformatics*, 25(11), 1422–1423. [https://biopython.org](https://biopython.org)  
- McKinney, W. (2010). *Data Structures for Statistical Computing in Python.* *Proc. 9th Python in Science Conf.*, 51–56. [https://pandas.pydata.org](https://pandas.pydata.org)  
- OpenPyXL Developers. (2023). *OpenPyXL: A Python library to read/write Excel files.* [https://openpyxl.readthedocs.io](https://openpyxl.readthedocs.io)  
- Simonov, K., et al. (2022). *PyYAML: YAML parser and emitter for Python.* [https://pyyaml.org](https://pyyaml.org)  
- Reitz, K., & Contributors. (2011). *Requests: HTTP for Humans.* [https://requests.readthedocs.io](https://requests.readthedocs.io)

#### Databases

- Federhen, S. (2012). *The NCBI Taxonomy Database.* *Nucleic Acids Research*, 40(D1), D136–D143. [https://www.ncbi.nlm.nih.gov/taxonomy](https://www.ncbi.nlm.nih.gov/taxonomy)  
- Guillou, L., et al. (2013). *The Protist Ribosomal Reference Database (PR2): A catalog of unicellular eukaryote SSU rRNA sequences with curated taxonomy.* *Nucleic Acids Research*, 41(D1), D597–D604. [https://pr2-database.org](https://pr2-database.org)  
- Guiry, M.D., & Guiry, G.M. (2024). *AlgaeBase: World-wide electronic publication, National University of Ireland, Galway.* [https://www.algaebase.org](https://www.algaebase.org)

---

## Acknowledgement

The tool **algaetax** was developed in the frame of a cooperation between the [Phycology working group of the University of Duisburg-Essen](https://www.uni-due.de/phycology/) and the [Department U2 Microbial Ecology of the Bundesanstalt für Gewässerkunde](https://www.bafg.de/DE/3_Beraet/4_Exp_oekologie/Planktonzoenosen_U2/planktonzoenosen_node.html).

---

## Troubleshooting

Running **algaetax** may occasionally lead to unexpected behavior or failed executions, especially during setup or database queries. Below is a list of common issues, their likely causes, and hints for troubleshooting. If issues persist, test **algaetax** with a small dataset to verify dependencies and database access before running large analyses.

### Common causes of workflow errors

- **Missing config file.** Ensure `config.yaml` exists and is valid.  
- **Wrong file paths.** Input and output paths must match the configuration.  
- **Empty or invalid input.** Check that the input list contains valid taxa.  
- **Database errors.** Network issues may block access to NCBI, PR2, or AlgaeBase.  
- **Unrecognized taxa.** See `taxa_not_found.csv` for unmatched entries.  
- **Conda issues.** Recreate the environment if dependencies fail.  
- **Permission denied.** Check read/write access to working directories. 
- **Invalid AlgaeBase API key.** AlgaeBase queries require a valid key.

### Installation problems

- **Conda not found.** Verify Conda or Miniconda is installed and in your `PATH`.  
- **Env creation failed.** Remove old environments and recreate from YAML.  
- **Dependency conflicts.** Avoid mixing Conda with `pip` or `mamba` installs. 

### Runtime or output problems

- **Empty results.** Input taxa may be filtered or unmatched — check logs.  
- **Run interrupted.** Restart if the workflow stopped or network dropped.  
- **Wrong database setup.** Confirm databases in `config.yaml` are correct.  
- **Errors in logs.** Inspect `algaetax.log` for details on failed queries.
