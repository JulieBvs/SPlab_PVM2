# Deciphering the heterogeneity of the Lyve1+ perivascular Macrophages in the mouse brain

## Article information

**Title:** Deciphering the heterogeneity of the Lyve1+ perivascular Macrophages in the mouse brain

**Authors:** C. Siret1, M. van Lessen2, J. Bavais1,3, H. W. Jeong4, S. K. Reddy Samawar5, K. Kapupara5, S. Wang1, M. Simic1, L. de Fabritus1, A. Tchoghandjian6, M. Fallet1, H. Huang1,7, S. Sarrazin1, M.H. Sieweke1,7, R. Stumm8, L. Sorokin5, R. H. Adams4, S. Schulte-Merker2, F. Kiefer4,9, S.A. van de Pavert1%

1 Aix-Marseille Univ, CNRS, INSERM, Centre d’Immunologie de Marseille-Luminy (CIML), Marseille, France
2 Institute for Cardiovascular Organogenesis and Regeneration, Faculty of Medicine, University of Münster, Münster, Germany
3 Turing Centre for Living systems, Marseille, France
4 Max Planck Institute for Molecular Biomedicine, Münster, Germany
5 Institute of Physiological Chemistry and Pathobiochemistry and Cells in Motion Interfaculty Centre (CIMIC), University of Münster, Münster, Germany
6 Aix-Marseille Univ, CNRS, INP, Inst Neurophysiopathol, Marseille, France
7 Center for Regenerative Therapies Dresden (CRTD), Technische Universität Dresden, 01307 Dresden, Germany
8 Institute of Pharmacology and Toxicology, Jena University Hospital, Jena, Germany
9 University of Münster, European Institute for Molecular Imaging (EIMI) and Cells in Motion Interfaculty Centre (CIMIC), Münster, Germany 

% Corresponding author: E-mail: vandepavert@ciml.univ-mrs.fr

**Summary:**
TODO: Add abstract

DOI : TODO: Add DOI
---
---

## Goal of the github
This github project contains the instructions and material to reproduce the analysis reported in the article (and more).
Source code (scripts and dockerfiles) are available in the github repository. Required data and builded Docker image is available on download. Instructions to reproduce the analysis are provided below.

To reproduce the analysis, you have to first, prepare the environments (see "Prepare the Environments" section below), then execute the analysis step by step (see "Run the analysis" section below).

---
---

## Description of the datasets

As described in the article, there are 8 datasets of myelin-depleted brain cortex single-cell mRNA sequencing in this study. Four of these datasets are replicates from juvenile mice (postnatal day 10), two are replicated from adult mice (7-11 weeks) and the last two ones are replicated from aged mice (18 months).

When downloading the code and data, you will obtains 8 files with names as below:

    All_brain_cortex_non_neural_cells
    ├── GSM3904816_Adult-1_gene_counts.tsv
    ├── GSM3904817_Adult-2_gene_counts.tsv
    ├── GSM3904819_Aged-1_gene_counts.tsv
    ├── GSM3904820_Aged-2_gene_counts.tsv
    ├── GSM3904823_Juvenile-1_gene_counts.tsv
    ├── GSM3904824_Juvenile-2_gene_counts.tsv
    ├── GSM3904825_Juvenile-3_gene_counts.tsv
    └── GSM3904826_Juvenile-4_gene_counts.tsv

---
---

## Prepare the environments

In order to prepare the environment for analysis execution, it is required to:

- Clone the github repository and set the WORKING_DIR environment variable
- Download the docker image tar file
- Install Docker
- Load the docker image on your system
- Download the count tables

Below you will find detailed instruction for each of these steps.

### Clone the github repository

Use you favorite method to clone this repository in a chosen folder. This will create a folder "All_brain_cortex_non_neural_cells" with all the source code. 

Then, you must set an environment variable called WORKING_DIR with a value set to the path to this folder.

For instance, if you have chosen to clone the Git repository in "/home/bavaisj/workspace", then the WORKING_DIR variable will be set to "/home/bavaisj/workspace/All_brain_cortex_non_neural_cells"

**On linux:**

    export WORKING_DIR=/home/bavaisj/workspace/All_brain_cortex_non_neural_cells

### Download the raw data

Count tables can be downloaded from Gene Expression Omnibus (GEO) and uncompressed. The GEO accession number is [GSE133283](https://zenodo.org/badge/DOI/10.5281/zenodo.3946154.svg).

To download and uncompress the data, use the following code:

**On linux:**

    cd $WORKING_DIR
    mkdir 00_RawData/
    cd 00_RawData
    
    wget 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM3904816&format=file&file=GSM3904816%5FAdult%2D1%5Fgene%5Fcounts%2Etsv%2Egz' -O GSM3904816_Adult-1_gene_counts.tsv.gz
    gzip -d GSM3904816_Adult-1_gene_counts.tsv.gz
    
    wget 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM3904817&format=file&file=GSM3904817%5FAdult%2D2%5Fgene%5Fcounts%2Etsv%2Egz' -O GSM3904817_Adult-2_gene_counts.tsv.gz
    gzip -d GSM3904817_Adult-2_gene_counts.tsv.gz

    wget 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM3904819&format=file&file=GSM3904819%5FAged%2D1%5Fgene%5Fcounts%2Etsv%2Egz' -O GSM3904819_Aged-1_gene_counts.tsv.gz
    gzip -d GSM3904819_Aged-1_gene_counts.tsv.gz

    wget 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM3904820&format=file&file=GSM3904820%5FAged%2D2%5Fgene%5Fcounts%2Etsv%2Egz' -O GSM3904820_Aged-2_gene_counts.tsv.gz
    gzip -d GSM3904820_Aged-2_gene_counts.tsv.gz

    wget 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM3904823&format=file&file=GSM3904823%5FJuvenile%2D1%5Fgene%5Fcounts%2Etsv%2Egz' -O GSM3904823_Juvenile-1_gene_counts.tsv.gz
    gzip -d GSM3904823_Juvenile-1_gene_counts.tsv.gz

    wget 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM3904824&format=file&file=GSM3904824%5FJuvenile%2D2%5Fgene%5Fcounts%2Etsv%2Egz' -O GSM3904824_Juvenile-2_gene_counts.tsv.gz
    gzip -d GSM3904824_Juvenile-2_gene_counts.tsv.gz

    wget 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM3904826&format=file&file=GSM3904826%5FJuvenile%2D4%5Fgene%5Fcounts%2Etsv%2Egz' -O GSM3904826_Juvenile-4_gene_counts.tsv.gz
    gzip -d GSM3904826_Juvenile-4_gene_counts.tsv.gz

Once done, you may obtain the following structure.

    All_brain_cortex_non_neural_cells
    └── 00_RawData
        ├── GSM3904816_Adult-1_gene_counts.tsv
        ├── GSM3904817_Adult-2_gene_counts.tsv
        ├── GSM3904819_Aged-1_gene_counts.tsv
        ├── GSM3904820_Aged-2_gene_counts.tsv
        ├── GSM3904823_Juvenile-1_gene_counts.tsv
        ├── GSM3904824_Juvenile-2_gene_counts.tsv
        ├── GSM3904825_Juvenile-3_gene_counts.tsv
        └── GSM3904826_Juvenile-4_gene_counts.tsv

### Download the reference files

The study uses references (genome annotations) you have to download. The annotations used during the study are available on Zenodo [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3949849.svg)](https://doi.org/10.5281/zenodo.3949849). Use the following command to download the tarball file and uncompress it.

Note: Since the reference files are used for the 4 single-cell samples analysis, they must be present in all the sample folder in the same 01_Reference subfolder. Instead of copying the files, we will create symbolic links:

**On linux:**

    cd $WORKING_DIR
    wget https://zenodo.org/record/3949849/files/SPlab_BecomingLTi_01_Reference.tar.gz?download=1 -O SPlab_BecomingLTi_01_Reference.tar.gz
    tar zxvf SPlab_BecomingLTi_01_Reference.tar.gz
    ln -s Embryo_Stage13.5_FetalLiver/01_Reference Embryo_Stage13.5_Periphery_CellRangerV3/01_Reference
    ln -s Embryo_Stage13.5_FetalLiver/01_Reference Embryo_Stage14.5_FetalLiver/01_Reference
    ln -s Embryo_Stage13.5_FetalLiver/01_Reference Embryo_Stage14.5_Periphery_CellRangerV3/01_Reference

These commands will create 4 sub-folders named 01_Reference:

    BecomingLTi
    ├── Embryo_Stage13.5_FetalLiver
    │   └── 01_Reference
    ├── Embryo_Stage13.5_Periphery_CellRangerv3
    │   └── 01_Reference
    ├── Embryo_Stage14.5_FetalLiver
    │   └── 01_Reference
    └── Embryo_Stage14.5_Periphery_CellRangerv3
        └── 01_Reference

### Download the Docker and Singularity images

Docker image tar file and Singularity img files are stored on Zenodo [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3949849.svg)](https://doi.org/10.5281/zenodo.3949849). Open a shell command and change dir to the root of the cloned Git repository (WORKING_DIR). Then execute the following commands to download the tarball file and untar  it:

**On linux:**

    cd $WORKING_DIR
    wget https://zenodo.org/record/3949849/files/SPlab_BecomingLTi_02_containers.tar.gz?download=1 -O SPlab_BecomingLTi_02_containers.tar.gz
    tar zxvf SPlab_BecomingLTi_02_containers.tar.gz

These commands will create 2 sub-folders named 02_Container:

    BecomingLTi
    ├── Embryo_Bulk_Stage13.5_2tissues
    │   └── 02_Container
    └── Embryo_Stage13.5_FetalLiver
        └── 02_Container

The first one contains a Docker image tar file used for the bulk RNA-seq analysis. The second one contains the Singularity images for the single-cell RNA-seq analysis. Since the singularity images are used for the 4 single-cell samples analysis, they must be present in all the sample folder in the same 02_Container subfolder. Instead of copying the image files, we will create symbolic links:

**On linux:**

    cd $WORKING_DIR
    ln -s Embryo_Stage13.5_FetalLiver/02_Container Embryo_Stage13.5_Periphery_CellRangerV3/02_Container
    ln -s Embryo_Stage13.5_FetalLiver/02_Container Embryo_Stage14.5_FetalLiver/02_Container
    ln -s Embryo_Stage13.5_FetalLiver/02_Container Embryo_Stage14.5_Periphery_CellRangerV3/02_Container

### Install Docker and Singularity

You need to install Docker and Singularity v2.6 on your system.

- To install Docker, follow the instructions here : https://docs.docker.com/get-docker/

- To install Singularity v2.6, follow the instructions here : https://sylabs.io/guides/2.6/admin-guide/

### Load docker images on the system

In order to execute analysis of the bulk RNA-seq, you must load the provided docker image onto your Docker. Docker must be installed on your system. 
See https://docs.docker.com/install/ for details on Docker installation.
Open a shell command and type:

**On linux:**

    docker load -i $WORKING_DIR/Embryo_Bulk_Stage13.5_2tissues/02_Container/splab_ilcyou_deg_gsea.tar

This command may take some time. If you encounter an issue loading some docker image layer, try again. Sometimes issue would be resolved. 

### Install Snakemake

If you want to take advantage of the workflow management we used for the single-cell RNA-seq analysis, you have to install Snakemake. See the official instruction and use your prefered solution:

https://snakemake.readthedocs.io/en/stable/getting_started/installation.html

---
---

## Run the analysis

There are two types of analysis in this study : bulk RNA-seq and single-cell RNA-seq. The bulk RNA-seq analysis uses the Docker image you loaded. The single-cell RNA-seq analysis uses the Singularity images and optionnaly Snakemake.

### Run the bulk RNA-seq analysis

The RNA-seq analysis are in two steps (step1 and step2). The first step make the QC, study the differentially expressed genes and their functionnal enrichment. The second step study the pattern of evolution of group of genes along the cell types (see article methods).

To run the step1 analysis, use the following command:

**On Linux:**

    docker run -v $WORKING_DIR:$WORKING_DIR -e WORKING_DIR=$WORKING_DIR splab_ilcyou_deg_gsea 'cd $WORKING_DIR/Embryo_Bulk_Stage13.5_2tissues/03_Script/step1;Rscript launch_reports_compilation.R'

To run the step2 analysis, use the following command:

**On Linux:**

     docker run -v $WORKING_DIR:$WORKING_DIR -e WORKING_DIR=$WORKING_DIR splab_ilcyou_deg_gsea 'cd $WORKING_DIR/Embryo_Bulk_Stage13.5_2tissues/03_Script/step2;Rscript launch_reports_compilation.R'

Each analysis will generate a result in $WORKING_DIR/Embryo_Bulk_Stage13.5_2tissues/05_output/step1 or $WORKING_DIR/Embryo_Bulk_Stage13.5_2tissues/05_output/step2.
In the output of the analysis, you will find a HTML file that contains the report of the analysis, with all figures. Some extra file are generated to export data in plain text.


### Run the single-cell RNA-seq analysis

The study contains 4 samples of single-cell RNA-seq data. Each sample have 5 step of analysis you will find the R script files in the subfolder 03_Script. The 5 steps are:

 * 01_QC : General quality control and bad cell removal
 * 02_GlobalHeterogeneity : First study of cell heterogeneity and sample contamination by undesired cell types
 * 03_GlobalHeterogeneity_NoContamination : Study of cell heterogeniety in absence of contamination
 * 04_Dynamics_Monocle : analysis of the cellular process dynamics using pseudotime analysis by Monocle
 * 05_Dynamics_RNAVelocity : analysis of the cellular process dynamics using RNA velocity (Velocyto)

Each step of analysis generates its own HTML report file and several output files. Some output files of some steps are used by other steps, making a complete workflow of analysis.

The simpliest way to run the complete single-cell analysis of a sample is to use the Snakemake workflow dedicated to each sample. The workflow is controled by a snakefile stored in the 04_Workflow subfolder of each sample folder. This workflow uses Singularity images (see above) to control the software environment for each analysis step. So you need both Snakemake and Singularity installed on your system to use this workflow.

In order to use the snakemake workflow, please type first the following commands:

     cd $WORKING_DIR
     ln -s Embryo_Stage13.5_FetalLiver/04_Workflow/snakefile.yml Embryo_Stage13.5_FetalLiver/snakefile.yml
     ln -s Embryo_Stage13.5_Periphery_CellRangerV3/04_Workflow/snakefile.yml Embryo_Stage13.5_Periphery_CellRangerV3/snakefile.yml
     ln -s Embryo_Stage14.5_FetalLiver/04_Workflow/snakefile.yml Embryo_Stage14.5_FetalLiver/snakefile.yml
     ln -s Embryo_Stage14.5_Periphery_CellRangerV3/04_Workflow/snakefile.yml Embryo_Stage14.5_Periphery_CellRangerV3/snakefile.yml

To run the analysis for the Embryo_Stage13.5_FetalLiver (for instance), then run the following commands:

Note: you have to manually change the "$WORKING_DIR" string in the snakemake command below by the value of the environment variable (i.e the path where you clone the project) because snakemake may not interpret the variable name correctly:

     cd $WORKING_DIR/Embryo_Stage13.5_FetalLiver
     snakemake -r --snakefile snakefile.yml --use-singularity --singularity-args "-B $WORKING_DIR:$WORKING_DIR"
     
To execute the analysis of the other sample, simply change folder to the target sample and run again the same snakemake command.

