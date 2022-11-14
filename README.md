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
Perivascular macrophages (pvMs) are associated with cerebral vasculature and mediate brain drainage and immune regulation. Here, using reporter mouse models, whole brain and section immunofluorescence, flow cytometry, and single cell RNA sequencing, besides the Lyve1+F4/80+ CD206+CX3CR1+ pvMs we identified a CX3CR1– pvM population that shares phagocytic functions and location. Furthermore, the brain parenchyma vasculature mostly hosts Lyve1+MHCII–  pvMs with low to intermediate CD45 expression. Using the double Cx3cr1GFP x Cx3cr1-Cre;RosatdT reporter mice for finer mapping of the lineages, we establish that CD45lowCX3CR1– pvMs are derived from CX3CR1+ precursors and require PU.1 during their ontogeny. In parallel, the Cxcr4-CreErt2;Rosa26tdT lineage tracing model supports a bone marrow-independent replenishment of all Lyve1+ pvMs in the adult mouse brain. Lastly, flow cytometry and 3D immunofluorescence analysis uncover increased percentage of pvMs following photothrombotic induced stroke. Our results thus show that the parenchymal pvM population is more heterogenous than previously described, and includes a CD45low and CX3CR1negative pvM population.

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

Count tables can be downloaded from Gene Expression Omnibus (GEO) and uncompressed. The GEO accession number is [GSE133283](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE133283).

To download and uncompress the data, use the following code:

**On linux:**

    cd $WORKING_DIR
    mkdir 01_RawData/
    cd 01_RawData
    
    wget 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM3904817&format=file&file=GSM3904817%5FAdult%2D2%5Fgene%5Fcounts%2Etsv%2Egz' -O GSM3904817_Adult-2_gene_counts.tsv.gz
    gzip -d GSM3904817_Adult-2_gene_counts.tsv.gz

    wget 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM3904818&format=file&file=GSM3904818%5FAdult%2D3%5Fgene%5Fcounts%2Etsv%2Egz' -O GSM3904818_Adult-3_gene_counts.tsv.gz
    gzip -d GSM3904818_Adult-3_gene_counts.tsv.gz

    wget 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM3904819&format=file&file=GSM3904819%5FAged%2D1%5Fgene%5Fcounts%2Etsv%2Egz' -O GSM3904819_Aged-1_gene_counts.tsv.gz
    gzip -d GSM3904819_Aged-1_gene_counts.tsv.gz

    wget 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM3904820&format=file&file=GSM3904820%5FAged%2D2%5Fgene%5Fcounts%2Etsv%2Egz' -O GSM3904820_Aged-2_gene_counts.tsv.gz
    gzip -d GSM3904820_Aged-2_gene_counts.tsv.gz

    wget 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM3904823&format=file&file=GSM3904823%5FJuvenile%2D1%5Fgene%5Fcounts%2Etsv%2Egz' -O GSM3904823_Juvenile-1_gene_counts.tsv.gz
    gzip -d GSM3904823_Juvenile-1_gene_counts.tsv.gz

    wget 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM3904824&format=file&file=GSM3904824%5FJuvenile%2D2%5Fgene%5Fcounts%2Etsv%2Egz' -O GSM3904824_Juvenile-2_gene_counts.tsv.gz
    gzip -d GSM3904824_Juvenile-2_gene_counts.tsv.gz

    wget 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM3904825&format=file&file=GSM3904825%5FJuvenile%2D3%5Fgene%5Fcounts%2Etsv%2Egz' -O GSM3904825_Juvenile-3_gene_counts.tsv.gz
    gzip -d GSM3904825_Juvenile-3_gene_counts.tsv.gz

    wget 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM3904826&format=file&file=GSM3904826%5FJuvenile%2D4%5Fgene%5Fcounts%2Etsv%2Egz' -O GSM3904826_Juvenile-4_gene_counts.tsv.gz
    gzip -d GSM3904826_Juvenile-4_gene_counts.tsv.gz

Once done, you may obtain the following structure.

    All_brain_cortex_non_neural_cells
    └── 01_RawData
        ├── GSM3904816_Adult-2_gene_counts.tsv
        ├── GSM3904817_Adult-3_gene_counts.tsv
        ├── GSM3904819_Aged-1_gene_counts.tsv
        ├── GSM3904820_Aged-2_gene_counts.tsv
        ├── GSM3904823_Juvenile-1_gene_counts.tsv
        ├── GSM3904824_Juvenile-2_gene_counts.tsv
        ├── GSM3904825_Juvenile-3_gene_counts.tsv
        └── GSM3904826_Juvenile-4_gene_counts.tsv


### Download the Docker

Docker image tar file is stored on Zenodo [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7319580.svg)](https://doi.org/10.5281/zenodo.7319580). Open a shell command and change dir to the root of the cloned Git repository (WORKING_DIR). Then execute the following commands to download the tarball file and untar it:

**On linux:**

    cd $WORKING_DIR
    mkdir 02_Container/
    cd 02_Container
    wget https://zenodo.org/record/7319580/files/splab_pvm2_seurat4.1.0_integration.tar?download=1 -O splab_pvm2_seurat4.1.0_integration.tar

These commands will create one sub-folder named 02_Container and containing the docker image tar file:

    All_brain_cortex_non_neural_cells
    └── 02_Container
        └── splab_pvm2_seurat4.1.0_integration.tar


### Install Docker

You need to install Docker on your system following the instructions here : https://docs.docker.com/get-docker/

### Load docker images on the system

In order to execute analysis, you must load the provided docker image onto your Docker. Docker must be installed on your system. 
See https://docs.docker.com/install/ for details on Docker installation.
Open a shell command and type:

**On linux:**

    docker load -i $WORKING_DIR/All_brain_cortex_non_neural_cells/02_Container/splab_pvm2_seurat4.1.0_integration.tar

This command may take some time. If you encounter an issue loading some docker image layer, try again. Sometimes issue would be resolved. 


## Run the analysis

### Quality control
First, you need to run the quality control step for each datasets following these command lines:

**On Linux:**

    docker run -v $WORKING_DIR:$WORKING_DIR -e WORKING_DIR=$WORKING_DIR splab_pvm2_seurat4.1.0_integration 'cd $WORKING_DIR/All_brain_cortex_non_neural_cells/03_Script/01_QC_Adult_rep2;Rscript launch_reports_compilation.R'
    docker run -v $WORKING_DIR:$WORKING_DIR -e WORKING_DIR=$WORKING_DIR splab_pvm2_seurat4.1.0_integration 'cd $WORKING_DIR/All_brain_cortex_non_neural_cells/03_Script/01_QC_Adult_rep3;Rscript launch_reports_compilation.R'
    docker run -v $WORKING_DIR:$WORKING_DIR -e WORKING_DIR=$WORKING_DIR splab_pvm2_seurat4.1.0_integration 'cd $WORKING_DIR/All_brain_cortex_non_neural_cells/03_Script/01_QC_Aged_rep1;Rscript launch_reports_compilation.R'
    docker run -v $WORKING_DIR:$WORKING_DIR -e WORKING_DIR=$WORKING_DIR splab_pvm2_seurat4.1.0_integration 'cd $WORKING_DIR/All_brain_cortex_non_neural_cells/03_Script/01_QC_Aged_rep2;Rscript launch_reports_compilation.R'
    docker run -v $WORKING_DIR:$WORKING_DIR -e WORKING_DIR=$WORKING_DIR splab_pvm2_seurat4.1.0_integration 'cd $WORKING_DIR/All_brain_cortex_non_neural_cells/03_Script/01_QC_Juvenile_rep1;Rscript launch_reports_compilation.R'
    docker run -v $WORKING_DIR:$WORKING_DIR -e WORKING_DIR=$WORKING_DIR splab_pvm2_seurat4.1.0_integration 'cd $WORKING_DIR/All_brain_cortex_non_neural_cells/03_Script/01_QC_Juvenile_rep2;Rscript launch_reports_compilation.R'
    docker run -v $WORKING_DIR:$WORKING_DIR -e WORKING_DIR=$WORKING_DIR splab_pvm2_seurat4.1.0_integration 'cd $WORKING_DIR/All_brain_cortex_non_neural_cells/03_Script/01_QC_Juvenile_rep3;Rscript launch_reports_compilation.R'
    docker run -v $WORKING_DIR:$WORKING_DIR -e WORKING_DIR=$WORKING_DIR splab_pvm2_seurat4.1.0_integration 'cd $WORKING_DIR/All_brain_cortex_non_neural_cells/03_Script/01_QC_Juvenile_rep4;Rscript launch_reports_compilation.R'

### Merge datasets from same time point
Then, datasets replicates from same time point are merged giving three seurat object:

**On Linux:**

    docker run -v $WORKING_DIR:$WORKING_DIR -e WORKING_DIR=$WORKING_DIR splab_pvm2_seurat4.1.0_integration 'cd $WORKING_DIR/All_brain_cortex_non_neural_cells/03_Script/02_GlobalHeterogeneity_Adult_Merge;Rscript launch_reports_compilation.R'
    docker run -v $WORKING_DIR:$WORKING_DIR -e WORKING_DIR=$WORKING_DIR splab_pvm2_seurat4.1.0_integration 'cd $WORKING_DIR/All_brain_cortex_non_neural_cells/03_Script/02_GlobalHeterogeneity_Aged_Merge;Rscript launch_reports_compilation.R'
    docker run -v $WORKING_DIR:$WORKING_DIR -e WORKING_DIR=$WORKING_DIR splab_pvm2_seurat4.1.0_integration 'cd $WORKING_DIR/All_brain_cortex_non_neural_cells/03_Script/02_GlobalHeterogeneity_Juvenile_Merge;Rscript launch_reports_compilation.R'

### Integration
Those three seurat object are integrated and analysed following steps described in the Methods of the paper.

**On Linux:**

    docker run -v $WORKING_DIR:$WORKING_DIR -e WORKING_DIR=$WORKING_DIR splab_pvm2_seurat4.1.0_integration 'cd $WORKING_DIR/All_brain_cortex_non_neural_cells/03_Script/04_Sample_Integration;Rscript launch_reports_compilation.R'

### Lyve1+ subcluster analysis
Finally, Lyve1+ cluster is extracted and re-analysed.

**On Linux:**

    docker run -v $WORKING_DIR:$WORKING_DIR -e WORKING_DIR=$WORKING_DIR splab_pvm2_seurat4.1.0_integration 'cd $WORKING_DIR/All_brain_cortex_non_neural_cells/03_Script/05_Lyve1_Cluster_From_Sample_Integration;Rscript launch_reports_compilation.R'
    docker run -v $WORKING_DIR:$WORKING_DIR -e WORKING_DIR=$WORKING_DIR splab_pvm2_seurat4.1.0_integration 'cd $WORKING_DIR/All_brain_cortex_non_neural_cells/03_Script/06_Lyve1_Cluster_From_Sample_Integration_Without_Apoptosis;Rscript launch_reports_compilation.R'






