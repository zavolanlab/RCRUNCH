# RCRUNCH
<div align="left">
    <img width="10%" align="left" style="margin-right: 20px;" src=images/rcrunch_logo_square.png>
</div> 

RCRUNCH is a workflow that identifies binding sites of RNA Binding Proteins (RBPs) and infers the binding motifs from (e)CLIP data. To accommodate the different types of targets that RBPs have, RCRUNCH offers a few variants and options that can be set by the user. RCRUNCH relies on the [Snakemake](https://github.com/snakemake/snakemake) workflow management system to coordinate the execution of the different components. It can be executed either locally or in a HPC system (e.g. SLURM).

## RCRUNCH components

RCRUNCH consists of the following components:

### <span style="color:purple">Read preprocessing</span>
* 3' or 5' adapter removal
* alignment of reads to the reference genome
* elimination of PCR duplicates (or UMIs)
* optional removal of reads that originate from abundant non-coding RNAs (e.g. tRNAs) 


### <span style="color:green">Splice-Junction-aware (transcriptomic) approach</span>
If the user chooses the Splice-Junction-aware approach (which we call the "TR" (transcriptomic) for simplicity) of RCRUNCH, some additional steps are performed to identify reads that map across splice junctions. That is, after all the preprocessing steps, the remaining alignments for foreground (CLIP) samples are used to select the most expressed transcript isoform for each gene and construct a dataset-specific transcriptome. Then the genome and transcriptome alignment files are jointly analyzed to identify the highest scoring alignment for each read. Peaks are then detected either on the genome or the transcriptome (see RCRUNCH model), treating individual transcripts as chromosomes. This approach allows for the detection and proper quantification of RBP binding sites in the vicinity or even spanning splice junctions.

### <span style="color:red">RCRUNCH model</span>

At the heart of RCRUNCH lies the RCRUNCH model for the detection of RBP-binding regions. Genome/transcriptome-wide identification of peaks corresponding to individual binding sites for an RBP is time consuming. For this reason RCRUNCH implements a two-step process:
1. Identify broader genomic regions that are enriched in reads in the foreground (CLIP) compared to the background sample
2. Identify individual peaks within these selected broader windows

> ???? Please read the "Methods" Section of the manuscript for an extensive description of RCRUNCH.

### <span style="color:blue">Motif analysis</span>
The last part of RCRUNCH is the de-novo prediction of binding motifs and the computation of enrichment scores for known (e.g. from [ATtRACT](https://attract.cnic.es/search)) and de-novo motifs for the RBP of interest.

<div align="left">
    <img width="50%" align="center" src=images/rcrunch_components.png>
</div> 


## Installation

### Requirements:

The following dependencies need to be installed on your system to be able to install and run RCRUNCH:
- [git](https://git-scm.com/)
- [conda](https://docs.conda.io/en/latest/) and/or [mamba](https://github.com/mamba-org/mamba)
- [singularity](https://sylabs.io/singularity/)

> ??? Currently, we only support Linux execution. Tested on CentOS 7.5, with miniconda 4.7.12 and Singularity 3.5.2.

### 1. Clone the repository

Go to the desired directory on your file system, then clone the repository and move into the RCRUNCH directory with:
```bash
git clone https://git.scicore.unibas.ch/zavolan_group/pipelines/RCRUNCH.git
cd RCRUNCH
```

### 2. Install Conda/Mamba

Workflow dependencies can be conveniently installed with the [Conda](https://docs.conda.io/projects/conda/en/latest/index.html)
package manager. If you haven't already, please install [Miniconda](https://docs.conda.io/en/latest/miniconda.html) 
for your system (Linux). Be sure to select Python 3 option. 
The workflow was built and tested with `miniconda 4.7.12`.
Other versions are not guaranteed to work as expected.

In addition to Miniconda, you will need the [Mamba](https://github.com/mamba-org/mamba) package manager, which -if you don't have it yet- needs to be installed in
the `base` conda environment with:

```bash
conda install mamba -n base -c conda-forge
```

### 3. Dependencies installation

For improved reproducibility and reusability of the workflow, each individual step of the workflow runs in its own [Singularity](https://sylabs.io/singularity/) container. As a consequence, running this workflow has very few individual dependencies. The container execution requires Singularity to be installed on the system where the workflow is executed. The functional installation of Singularity requires root privileges, and Conda currently only provides Singularity for Linux architectures. If you do not have root privileges on the machine you want to run the workflow please install Singularity separately and in privileged mode, depending on your system. You may have to ask an authorized person (e.g., a systems administrator) to do that. This will almost certainly be required if you want to run the workflow on a high-performance computing (HPC) cluster (which is highly recommended).

After installing Singularity, install the remaining dependencies with:


```bash
mamba env create -f install/environment.yml
```

This will create an environmnent `rcrunch`

### 4. Activate the environmnent

You can activate the Conda environment with:
```bash
conda activate rcrunch
```

### 5. Execution of RCRUNCH


### 5a. Fill in the necessary config file

In order to run RCRUNCH, please fill in the organism related data and the experiment-dependent parameters for the different samples in the file `config.yaml`.

> ??? For your convenience an empty [config.yaml](config.yaml) file is available to fill in.


### 5b. Dry run and DAG generation (optional)

You can generate a dry run of the pipeline by running:
```bash
snakemake \
-np \
-s workflow/Snakefile \
--profile profiles/local-singularity/ \
--configfile config.yaml
```


A directed acyclic graph (dag) of the run can be generated by running:
```bash
snakemake --dag -np --use-singularity | dot -Tpng > dag.png
```

or a rule graph can be generated by running:
```bash
snakemake --rulegraph --configfile config.yaml | dot -Tpng > rulegraph.png
```

### 5c. RCRUNCH execution

Finally you can trigger the pipeline by running:
```bash
bash run_local_singularity.sh
```

- We also support execution in the **SLURM** workload manager.
If you use SLURM start running RCRUNCH like this:

```bash
bash run_slurm_singularity.sh
```

> ??? Note: It might be useful to run rcrunch in the background, with the aid of, for example, `nohup`, which will create a nohup.out file containing the info of the run:
```bash
nohup bash run_local_singularity.sh &
```

> ??? Note: You can use any one of the: `nohup`, `screen` or `tmux`.

## Testing
To ensure that the version is working properly you can test it by:
1. Activate the Conda environment with:
```bash
conda activate rcrunch
```
2. Run:
```bash
bash test/test_singularity_execution/test_local.sh
```
or for **SLURM** workload manager,
```bash
bash test/test_singularity_execution/test_slurm.sh
```


## Output architecture
