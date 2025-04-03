---
title: "Parallel Evolution Focused on Copy Number Abnormalities of Chromosome 1 Drives the Progression of Multiple Myeloma at The Single-Cell Level"
output: html_document
date: "2025-02-28"
version: beta 0.0.1
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```
# Parallel Evolution Focused on Copy Number Abnormalities of Chromosome 1 Drives the Progression of Multiple Myeloma at The Single-Cell Level

![image](https://github.com/user-attachments/assets/577824fc-6b27-45f0-80fd-e613bd0a3e84)

## A. Introduction

We performed single cell sequencing of genomic DNA derived from 1,216 cells isolated from patients (n=7) with newly diagnosed multiple myeloma (NDMM) associated with copy number gain of chromosome 1q. The data highlights parallel evolution focused on copy number abnormalities of chromosome 1 as a being a key driver mechanism underlying disease evolution and progression to multiple myeloma (MM).

DLP+7_daVinch is a R package to analyze SCNAs identified through single-cell whole genome sequencing (scWGS) with HMMcopy and generate hierarchical clustering heatmaps, ploidy plots, and phylogenetic tree. 

## B. Install R package

```{Installation}
install.packages("devtools")
library(devtools)

devtools::install_github("sanghoonleepitt/iGenSigRx")
library(HMMCopyClusterPloidyv1.2)
```

## Step1. Preprocess HMMCopy CNA state data  

### Read the HMMcopyRead.csv.gz file and process it.

  - Remove outlier bins that have NA values > 20 % of cells 
  - Remove artifact cells with Gini coefficient > median(GiniCoefficient)+ mad(GiniCoefficient)

```{Preprecessing}
####### ++++++++++++  Patient #3, 01_206_143839A
HMMcopyReadFile <- "YourDirectory/01_206_143839A/hmmcopy_reads.csv.gz"  
MySampleID="01_206_143839A"

HMMcopy_AfterRmvCell_WithoutNABinRegion <- CNVStateProcess(HMMcopyReadFile,RemoveArtifactCell=TRUE, BinRegionDelete=TRUE, MySampleID)
dim(HMMcopy_AfterRmvCell_WithoutNABinRegion) #  656 5228

## This stepp will generate
## - "HMMcopy_StateIdeal_RmvCellByGiniMean_WithoutBinRegionAnnot_01_213_143839A.rds"
## - "HMMcopy_StateIdeal_RmvArtifactCNV_ByMeanGini_01_235_143921A.rds" 
```

## Step2. A draft hierarchical clustering heatmap for all chromosomes. 
This heatmap still includes artifact cells of clusters, such as cluster #1 in the heatmap below. 
This step is just to see the hierarchcial clustering with artifict cells.  But, this step is required to make ploidy plot input. 

```{Hiearchical cluatering}
HMMcopy_AfterRmvCell_WithoutNABinRegion <- readRDS("HMMcopy_StateIdeal_RmvArtifactCNV_ByMeanGini_01_040_143929A.rds") # This is made from Step1. 
HierarchialClusteringHeatmap(Processed_HMMcopyStateData=HMMcopy_AfterRmvCell_WithoutNABinRegion, MySampleID="01_206_143839A", AllCell="WithoutArtifact", ChromosomeToSubset="All",ClusterSize=9)

## This stepp will generate
## "CompHeatmap_CNV_WithoutArtifact_AllChr_01_206_143839A_9Cluster.pdf"
## "RowSplitHeatmapDraw_WithoutArtifact_AllChr_01_206_143839A.rds"
```
![image](https://github.com/user-attachments/assets/24b0a3c8-21e2-4aa7-9f21-311e6a6f3420)

## Step3. Remove artifact cells by CNA state Low or High, and asign cluster size
KM_ClusterSize is determined after looking at the heatmap in Step2.   
Remove cells that have CNA state 0 or 1 in 20% of total bins (n=5259) and have CNA state ≥ 5 in more than 10 bins.
Remove clusters that have cells less than 10.

```{remove artifact cells}
## Patient #3, 01_206_143839A
CNVStateData_AfterRmvHighGiniMeanFile <- "HMMcopy_StateIdeal_RmvCellByGiniMean_WithoutBinRegionAnnot_01_206_143839A.rds"  # This file is generated in Step1.
KM_ClusterSize=6
ThresholdLowState=1; ThresholdHighState=10 # default 10
MySampleID <- "01_206_143839A"
RemoveArtifactCellHeatmap(CNVStateData_AfterRmvHighGiniMeanFile, MySampleID, KM_ClusterSize, ThresholdLowState, ThresholdHighState)

## This stepp will generate
## CompHmapCNV_RmvByGiniMeanLowHighState_AllChr_09_025_143855A_LowSt0.18_527Cells_5clu.pdf
## "RowSplitHeatmapDraw_WithoutArtifactLowHighState_AllChr_",MySampleID,"_",KM_ClusterSize,"clu.rds"
## "ClusterNumber_15LowCell_01_235_143921A_9_8_6_4.txt" This one tells which clusters of cells < 15 were removed. 
```
![image](https://github.com/user-attachments/assets/11bb850c-a1a0-405f-ba1b-ddb86a32b694)
 


## Step4. Final hierarchical clustering heatmap for all chromosome
Check the heatmaps made in Step3 and remove normal cluster cells. 
Make a heatmap for cell clusters in chromosome 1 only following the same cluster order in the heatmap of all chromosome. 

```{Make a refined heatmap}
HierarchialClusterInfoFile <- "RowSplitHeatmapDraw_WithoutArtifactLowHighState_AllChr_01_206_143839A_6clu.rds"
CNVStateData_AfterRmvHighGiniMeanFile <-"HMMcopy_StateIdeal_RmvCellByGiniMean_01_206_143839A_LowHighState1.rds"
MySampleID <- "01_206_143839A"
ClusterNumbToRmv=6
ClusterReorder <- c(1,3,5,4,2)

## This stepp will generate
## HMMcopy_StateIdeal_AfterRmvArtifactCell_RmvReorderCluster_MySampleID.rds
## CompHmapCNV_RmvByGiniMeanLowHighState_LowSt_AllChr_01_213_143839A_389cells_RordClone.pdf
## CompHmapCNV_RmvByGiniMeanLowHighState_LowSt_Chr1only_01_213_143839A_389cells_RordClone.pdf
```

![image](https://github.com/user-attachments/assets/58a54abc-0a12-4286-8daf-29e81d571d89)

![image](https://github.com/user-attachments/assets/a71c086e-d38a-48ee-b328-c0ce6437f9b1)


## Step5. Heatmap for chromosome 1 pericentromeric region

```{ }
HMMcopy_StateIdeal_RmvArtifactCNV_ByMeanGini_File<- "HMMcopy_StateIdeal_RmvArtifactCNV_ByMeanGini_01_206_143839A.rds"
HMMcopy_StateIdeal_RmvArtifactCNV_ByMeanGini_StateLowHighReorderFile <-  "HMMcopy_StateIdeal_AfterRmvArtifactCell_RmvReorderCluster_01_206_143839A.rds"
MySampleID <- "01_206_143839A"
Chr1Centromere_Heatmap(HMMcopy_StateIdeal_RmvArtifactCNV_ByMeanGini_File, HMMcopy_StateIdeal_RmvArtifactCNV_ByMeanGini_StateLowHighReorderFile, MySampleID)

## This step will generate "CompHeatmap_Chr1CentromereRegion_01_235_143921A_OOOcells.pdf"

```

## Step6. Ploidy plot for all chromosomes.

```{Ploidyplot}
library(regioneR)
library(karyoploteR)
library(CopyNumberPlots)

CNVSegmentFile="/YourDirectory/01_206_143839A/hmmcopy_segments.csv.gz";
HierarchialClusterInfoFile="RowSplitHeatmapDraw_WithoutArtifact_AllChr_01_206_143839A.rds";
HierarchialClusterInfoAfRmvClusterFile="RowSplitHeatmapDraw_WithoutArtifactLowHighState_AllChr_01_206_143839A_6clu.rds";
CNVStateData_AfterRmvHighGiniMeanFile="HMMcopy_StateIdeal_AfterRmvArtifactCell_RmvReorderCluster_01_206_143839A.rds";  
ClusterNumbToRmv=NULL;  ClusterNumbToPick=NULL; MySampleID="01_206_143839A"

Karyoplote_BySegment(CNVSegmentFile=CNVSegmentFile,HierarchialClusterInfoFile=HierarchialClusterInfoFile,
                     HierarchialClusterInfoAfRmvClusterFile=HierarchialClusterInfoAfRmvClusterFile,HierarchialClusterInfoAfRmvClusterChr1File=NULL,
                     CNVStateData_AfterRmvHighGiniMeanFile=CNVStateData_AfterRmvHighGiniMeanFile,
                     CNVStateData_AfterRmvHighGiniMeanClusterFile=CNVStateData_AfterRmvHighGiniMeanClusterFile,
                     ClusterNumbToRmv,  ClusterNumbToPick, ClusterNumbToPick_InChr1=ClusterNumbToPick_InChr1,MySampleID=MySampleID )

## This step will generate "OutKaryoplote_AllCell_BeforeRmvHighGiniMean_BySegMedian_ByCluster_01_206_143839A_PloidyLine651c.pdf"
```
![image](https://github.com/user-attachments/assets/91cd7c5e-a5f4-4b63-b834-f376feecce40)

## Step7. Subclone phylogenetic tree

```{Phylogenetic tree}
## Pt3, 01_206  clone evolution
TreeData <- data.frame("sample" = c("01_206","01_206","01_206","01_206","01_206"), #name of sample, useful to read in all samples at once if processing multiple
                           "from" = c("x",   "Clone1", "x",      "Clone3","Clone4"), #parental Clone
                           "to"= c("Clone1", "Clone2", "Clone3", "Clone4","Clone5"))  #child Clone
                           #"score" = c() #a metric of CNV burden 
# Pt3 01_206 cell numbers 
NumbCellData <- data.frame("sample" = c("01_206","01_206", "01_206", "01_206","01_206","01_206"), #name of sample, useful to read in all samples at once if processing multiple
                                  "name" = c("x", "Clone1", "Clone2", "Clone3", "Clone4", "Clone5"), #Clone
                                  "counts"= c(0,19, 133,59,415, 25)) #number of cells in Clone
MySampleID <- "01_206";

## Make Phylogenetic tree
PhylogeneticTree(df=TreeData, num_cells_in=NumbCellData, MySampleID=MySampleID)

## This step will generate "OutPhylogeneticTree_01_206.pdf"

```
![image](https://github.com/user-attachments/assets/c1a676f5-0839-40b8-818b-8702d3bd141f)


contact: Sanghoon Lee, lees130@nyulangone.org

The end
