

# Modeling the competing effects of the immune system and EMT on tumor development

This repository contains the code used for model development, simulation, and data analysis associated with the paper:

Daniel R Bergman, Matthew K Karikomi, Qing Nie, Adam L MacLean (2019).  Modeling the competing effects of the immune system and EMT on tumor development. *bioRxiv*, [doi:10.1101/615971](https://doi.org/10.1101/615971).


## Description 

We develop an agent-based model to describe the relationships between cancer, the immune system, and EMT. Agents in the model are epithelial tissue cells with oncogenic potential. Cells can either be mutation-free or, resulting from DNA damage during the cell cycle, may harbor any combination of three possible pathway mutations.

We model immune cell populations as continuous variables (using differential equations), these include natural killer cells, cytotoxic T cells, and T regulatory cells. TGF-\(\beta\) is also modeled as a continuous variable within the tumor microenvironment. Tissue cells can transition between epithelial and mesenchymal phenotypes in response to cell-intrinsic and extrinsic factors in a plastic manner. We provide details of the simulation of patient cohorts to study the effects of different EMT and immune/inflammation regimes.


## Contents
- model/ contains the scripts (in MATLAB) required to run the model.
- analyzedata/ contains the scripts (in R) required to query [TCGA](https://www.cancer.gov/about-nci/organization/ccg/research/structural-genomics/tcga) and perform data analysis.


## Usage

Follow these steps to run the model:

1. Start MATLAB
2. In MATLAB cd to the repository tumor-immune-emt-code/

