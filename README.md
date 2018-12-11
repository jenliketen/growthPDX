# growthPDX

## Overview

The PDX Growth Signatures project aims to use machine learning methods to discover novel biomarkers associated with tumor growth rates in patient-derived xenograft (PDX) mouse models. Studies have shown that the rate of engraftment of PDXs is correlated with survival outcomes, such as this recent article about [head-and-neck cancer](https://www.cell.com/cell-reports/pdf/S2211-1247(18)31567-5.pdf). As PDXs have been widely used in drug development and discovery, we propose to study the gene signatures associated with PDX growth rates to develop robust models for cancer risk stratification and drug targeting.

## Getting Started

The project is divided into the following components:

   1. Correlation analysis between gene expression and PDX growth features (doubling time, survival, and slope of growth curve; mainly doubling time)
   2. Differential expression analysis looking at samples with high vs. low doubling times
   3. Gene-set enrichment analysis to identify pathways that are significantly altered between the two classes in doubling time
   4. Machine learning for biomarker discovery `UPCOMING`

Please head over to the [Wiki](https://github.com/jenliketen/growthPDX/wiki) for work progress updates, results, and important assumptions made throughout the course of the project.

## Installation
   * Install R version 3.5.1 from [CRAN](https://cran.r-project.org/)
   * Install the following R packages:
    
    # devtools
    install.packages("devtools")
   
    # Xeva
    devtools::install_github("bhklab/Xeva")
    
    # BiocManager
    if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    
    # Biobase
    BiocManager::install("Biobase", version = "3.8")
    
    # ggplot2
    install.packages("ggplot2")
    
    # BBmisc
    install.packages("BBmisc")
    
    # edgeR
    BiocManager::install("edgeR", version = "3.8")
    
    # limma
    BiocManager::install("limma", version = "3.8")
    
    # piano
    BiocManager::install("piano", version = "3.8")
