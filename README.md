---
title: "Parallel Evolution Focused on Copy Number Abnormalities of Chromosome 1 Drives the Progression of Multiple Myeloma at The Single-Cell Level"
output: html_document
date: "2025-07-17"
version: beta 0.0.1
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

# Parallel Evolution Focused on Copy Number Abnormalities of Chromosome 1 Drives the Progression of Multiple Myeloma at The Single-Cell Level

![image](https://github.com/user-attachments/assets/577824fc-6b27-45f0-80fd-e613bd0a3e84)

## A. Introduction

We performed single cell sequencing of genomic DNA derived from 2,521 cells isolated from patients (n=5) with newly diagnosed multiple myeloma (NDMM) associated with copy number gain of chromosome 1q. The data highlights parallel evolution focused on copy number abnormalities of chromosome 1 as a being a key driver mechanism underlying disease evolution and progression to multiple myeloma (MM).

DLP+7_daVinch is a R package to analyze SCNAs identified through single-cell whole genome sequencing (scWGS) with HMMcopy and generate clustering heatmaps, ploidy plots, and phylogenetic tree. 

## B. Install R package

```{r Installation}
install.packages("devtools")
library(devtools)

devtools::install_github("sanghoonleepitt/iGenSigRx")
library(HMMCopyClusterPloidyv1.2)
```

## Step1. Preprocess HMMCopy CNA state data  

Read the HMMcopyRead.csv.gz file and process it.

  - Remove outlier bins that have NA values > 10 % of cells 
  - Remove artifact cells with Gini coefficient > median(GiniCoefficient)+ mad(GiniCoefficient)

```{Preprecessing}
####### ++++++++++++  Patient #3, 01_206_143839A
HMMcopyReadFile <- "YourDirectory/01_206_143839A/hmmcopy_reads.csv.gz"  
MySampleID="01_206_143839A";  ThresholdLowState=1; ThresholdHighState=10

HMMcopy_AfterRmvCell_WithoutNABinRegion <- CNVStateProcess(HMMcopyReadFile, RemoveArtifactCell=TRUE, BinRegionDelete=TRUE, MySampleID, ThresholdLowState, ThresholdHighState)
dim(HMMcopy_AfterRmvCell_WithoutNABinRegion)  #  656 5228

## This stepp will generate
## - "HMMcopy_StateIdeal_RmvCellByGiniMean_WithoutBinRegionAnnot_01_213_143839A.rds"
## - "HMMcopy_StateIdeal_RmvArtifactCNV_ByMeanGini_01_235_143921A.rds" 
```

## Step2. Run PCA and Kmean.  Make a clustering heatmap for all chromosomes.  - This heatmap still includes artifact or low quality cells of clusters.  
      Give 4 or 5 KmeanClusterSize and see the clusters first. You can choose which cluster to remove and how to order clusters. 
      This step is just to see the clustering with artifict cells.  But, this step is required to make ploidy plot input. 

```{PCA and Kmean cluatering}
Processed_HMMcopyStateData <- readRDS(paste0(VignetteDir , "HMMcopy_StateIdeal_RmvCellByGiniMeanLowHighState_NoBinRegionAnnot_01_206_143839A.rds"))
MySampleID="01_206_143839A";  KmeanClusterSize=5

PCAKmeanClustering(Processed_HMMcopyStateData,  MySampleID=MySampleID, KmeanClusterSize)

### This will generate three files 
# "CompHmapCNV_RmvByGiniMeanLowHighState_AllChr_01_040_143929A_242Cells_5clu.pdf"  Overall heatmap with clusters. Still normal or low quality cell of cluster. 
# "RowSplitHeatmapDraw_WithoutArtifactLowHighState_AllChr_01_040_143929A_5clu.rds" ## It is a list. Names are cluster number and elements are row numbers. 
# ClusterNumber_10LowCell_01_040_143929A_.txt ## show th cluster number that has < 10 cells and excluded 


```
<img width="457" height="457" alt="image" src="https://github.com/user-attachments/assets/8b5def35-b6c5-4421-8ba0-203ed44d09ff" />


## Step3. You got clustering heatmap idea in step2. Remove normal or low quality cell clusters and make a final heatmap. 
            This is the final heatmap for all chromosomes and and chromosome1. Stack the clusters by clone order. KM_ClusterSize is determined after looking at the heatmap in Step2.   

```{Remove normal cells and Reorder Cluster}
## Patient #3, 01_206_143839A
PCAKmeanClusterInfoFile <- paste0(VignetteDir, "RowSplitHeatmapDraw_WithoutArtifactLowHighState_AllChr_01_206_143839A_5clu.rds")
CNVStateData_AfterRmvHighGiniMeanFile <- paste0(VignetteDir, "HMMcopy_StateIdeal_RmvCellByGiniMeanLowHighState_Kmean_01_206_143839A.rds")
MySampleID <- "01_206_143839A"
ClusterNumbToRmv="group2"
ClusterReorder <- c("group4","group5","group1","group3")

RemoveNormalClusterReorderHeatmap(PCAKmeanClusterInfoFile, CNVStateData_AfterRmvHighGiniMeanFile,MySampleID,ClusterNumbToRmv, ClusterReorder)

## It will generate
# "HMMcopy_StateIdeal_AfterRmvArtifactCell_RmvReorderCluster_01_206_143839A.rds"
# CompHmapCNV_RmvByGiniMeanLowHighState_LowSt_AllChr_01_206_143839A_651cells_RordClone.pdf
# CompHmapCNV_RmvByGiniMeanLowHighState_LowSt_Chr1only_01_206_143839A_651cells_RordClone.pdf

```
<img width="456" height="458" alt="image" src="https://github.com/user-attachments/assets/ca818eda-2261-4c7c-8313-6a28e4c27638" />

<img width="367" height="367" alt="image" src="https://github.com/user-attachments/assets/7611cb10-8a20-46ac-9016-dcf43aad1722" />


## Step4. Extract clone number for each cell. This is to split the merged .bam and run 

```{xtract clone number for each cell}
## 3rd patients, 01_206_143839A
PCAKmeanClusterInfoFile <- paste0(VignetteDir, "RowSplitHeatmapDraw_WithoutArtifactLowHighState_AllChr_01_206_143839A_6clu.rds")
CNVStateData_AfterRmvHighGiniMeanFile <- paste0(VignetteDir, "HMMcopy_StateIdeal_RmvCellByGiniMean_01_206_143839A_LowHighState1.rds")
MySampleID <- "01_206_143839A"
ClusterNumbToRmv=6
ClusterReorder <- c(1,3,5,4,2)

RemoveNormalClusterReorder_CloneNumb(PCAKmeanClusterInfoFile, CNVStateData_AfterRmvHighGiniMeanFile, MySampleID,ClusterNumbToRmv, ClusterReorder) 


```

## Step5. Heatmap for chromosome 1 pericentromeric region

```{ }
## Pt. #3 01_206_143839A
PCAKmeanClusterInfoFile <- paste0(VignetteDir, "RowSplitHeatmapDraw_WithoutArtifactLowHighState_AllChr_01_206_143839A_5clu.rds")
HMMcopy_AfterRmvCell_WithoutNABinRegionFile <- paste0(VignetteDir, "HMMcopy_StateIdeal_AfterRmvArtifactCell_RmvReorderCluster_01_206_143839A.rds")
HMMcopy_StateIdeal_RmvArtifactCNV_ByMeanGiniStateLowHighFile <- paste0(VignetteDir, "HMMcopy_StateIdeal_RmvCellByGiniMean_01_206_143839A_LowHighState1.rds") # This file is from Step1 output.It has Bin Region
MySampleID <- "01_206_143839A" 

Chr1Centromere_Heatmap(PCAKmeanClusterInfoFile, HMMcopy_AfterRmvCell_WithoutNABinRegionFile, HMMcopy_StateIdeal_RmvArtifactCNV_ByMeanGiniStateLowHighFile, MySampleID)

## This step will generate "CompHeatmap_Chr1CentromereRegion_01_206_143839A_OOOcells.pdf"

```

## Step6. Ploidy plot for all chromosomes.

```{Ploidyplot}
library(regioneR)
library(karyoploteR)
library(CopyNumberPlots)

### Patient #3 - 01_206_143839A
CNVSegmentFile=paste0(DLPpHMMCopyOutDir,"/01_206_143839A/hmmcopy_segments.csv.gz");
PCAKmeanClusterInfoAfRmvClusterFile=paste0(VignetteDir, "/RowSplitHeatmapDraw_WithoutArtifactLowHighState_AllChr_01_206_143839A_5clu.rds");
CNVStateData_AfterRmvHighGiniMeanFile=paste0(VignetteDir, "/HMMcopy_StateIdeal_AfterRmvArtifactCell_RmvReorderCluster_01_206_143839A.rds" );   #"/HMMcopy_StateIdeal_RmvCellByGiniMean_WithoutBinRegionAnnot_01_206_143839A.rds");
ClusterNumbToRmv=NULL;  ClusterNumbToPick=NULL; MySampleID="01_206_143839A"

Karyoplote_BySegment(CNVSegmentFile=CNVSegmentFile,
                     PCAKmeanClusterInfoAfRmvClusterFile=PCAKmeanClusterInfoAfRmvClusterFile,PCAKmeanClusterInfoAfRmvClusterChr1File=NULL,
                     CNVStateData_AfterRmvHighGiniMeanFile,
                     CNVStateData_AfterRmvHighGiniMeanClusterFile=NULL,
                     ClusterNumbToRmv,  ClusterNumbToPick, ClusterNumbToPick_InChr1=ClusterNumbToPick_InChr1,MySampleID=MySampleID )

## This step will generate "OutKaryoplote_AllCell_BeforeRmvHighGiniMean_BySegMedian_ByCluster_01_206_143839A_PloidyLineOOOc.pdf"
```
<img width="809" height="378" alt="image" src="https://github.com/user-attachments/assets/3629a474-6076-412e-8b40-570083394f12" />


## Step7. Subclone phylogenetic tree

```{r phylogenetic tree}
  library(igraph)
  library(scales)
  library(ggtree) ## bioconductor
  library(ggrepel)
  library(glue)
  library(stringr)
  library(cowplot)
  library(ggtext)
  library(data.table)
  library(reshape2)
  library(ggpubr)

## Pt3, 01_206  clone evolution
TreeData <- data.frame("sample" = c("01_206","01_206","01_206"), #name of sample, useful to read in all samples at once if processing multiple
                           "from" = c("x",   "Clone1", "Clone2_3" ), #parental Clones. 'x' is a virtual normal. 
                           "to"= c("Clone1", "Clone2_3", "Clone4" ))  #child Clone
                           #"score" = c() #a metric of CNV burden 
# Pt3 01_206 cell numbers 
NumbCellData <- data.frame("sample" = c("01_206","01_206", "01_206", "01_206"), #name of sample, useful to read in all samples at once if processing multiple
                                  "name" = c("x", "Clone1", "Clone2_3", "Clone4"), #Clone
                                  "counts"= c(0, 16, 607, 28)) #number of cells in Clone
MySampleID <- "01_206";
## Make Phylogenetic tree
PhylogeneticTree(df=TreeData, num_cells_in=NumbCellData, MySampleID=MySampleID)

## This step will generate "OutPhylogeneticTree_01_206.pdf"

```
<img width="292" height="310" alt="image" src="https://github.com/user-attachments/assets/15577b5f-248b-46d9-a800-7a342eacb4f8" />



contact: Sanghoon Lee, lees130@nyulangone.org

The end
