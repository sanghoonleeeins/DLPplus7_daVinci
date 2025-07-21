
#' Title  CountNumbNA
#'
#' @param x
#'
#' @returns Count of NAs
#' @export
#'
#' @examples CountNumbNA(x)
CountNumbNA <- function(x) {
  # x<- HMMcopy_StateIdeal_Tp[,1]
  NumbNA <- table(is.na(x))[names(table(is.na(x)))=="TRUE"]
  return(NumbNA)
}


#' Title CNVStateProcess
#'
#' @param HMMcopyReadFile original HMMcopy output file
#' @param RemoveArtifactCell  boolean
#' @param BinRegionDelete   boolean
#' @param MySampleID  manually typed in.
#' @param ThresholdLowState  Threshold to remove cells by the percentage of bins with state 0 or 1
#' @param ThresholdHighState  Threshold to remove cells by the number of bins with state > 5
#'
#' @return Processed CNV state data for making hierarchical clustering
#' @export
#'
#' @examples CNVStateProcess(HMMcopyReadFile, RemoveArtifactCell, BinRegionDelete, MySampleID)
CNVStateProcess <- function(HMMcopyReadFile, RemoveArtifactCell=TRUE, BinRegionDelete=TRUE, MySampleID="MySampleID", ThresholdLowState=0.2, ThresholdHighState=10) {
       # MySampleID <- gsub("(.*)04_DLPplus_1qGain_8ptData\\/|\\/hmmcopy_reads.csv.gz","", HMMcopyReadFile)
       # MySampleID <- gsub("-","_", MySampleID)

       print(paste0("HMMcopyReadFile: ", HMMcopyReadFile))
       HMMcopyRead <- data.table::fread(HMMcopyReadFile, header=TRUE, stringsAsFactors=FALSE);

       # I need only 'idea=TRUE' bins and remove chrX and chrY
       HMMcopy_CopyStateIdeal <- dplyr::filter(HMMcopyRead, ideal=="TRUE", !chr %in% c("chrX","chrY")); dim(HMMcopy_CopyStateIdeal) # 6155368      17 # Pt5. [1] 6560407      17

       ## I need 'StartEndChrRead","copy", "state", or "cell_id" columns.
       HMMcopy_CopyStateIdealSelect <- dplyr::mutate(HMMcopy_CopyStateIdeal, StartEndChr=paste0(chr,"_", start,"_",end)) # %>%    # "_",seq(1:nrow(HMMcopy_CopyStateIdeal
       HMMcopy_CopyStateIdealSelect <- dplyr::select(HMMcopy_CopyStateIdealSelect, c(StartEndChr, copy, state, cell_id));  # 5678019       4

       ######## =========== ######## =========== ######## =========== ######## ===========
       ## Step3. state or copy value per cell, and full_join
       ######## =========== ######## =========== ######## =========== ######## ===========
       AllCellID <- names(table(HMMcopy_CopyStateIdealSelect$cell_id)); length(AllCellID); AllCellID[1:3] # 1216  # "01-040-143929A-R03-C07" "01-040-143929A-R03-C08" "01-040-143929A-R03-C09"
       LoopNumb<-0; HMMcopy_CopyStateIdealSelect_AllCellID <- data.frame()
       for(EachCellID in AllCellID) {
           # EachCellID <- AllCellID[500]; print(EachCellID) # "01-040-143929A-R03-C07"
           LoopNumb <- LoopNumb+1; #print(LoopNumb);

           if(LoopNumb==1) print("It takes 1~2 minutes")
           if(LoopNumb %% 100 == 0)  print(paste0("Out of ",length(AllCellID), " cells, currently ", LoopNumb, " are being processed"));

           HMMcopy_CopyStateIdealSelect_ByCellID <-  dplyr::filter(HMMcopy_CopyStateIdealSelect, cell_id==EachCellID) #  %>%
           HMMcopy_CopyStateIdealSelect_ByCellID <-  dplyr::select(HMMcopy_CopyStateIdealSelect_ByCellID, -cell_id); # dim(HMMcopy_CopyStateIdealSelect_ByCellID) # 5404 3
           colnames(HMMcopy_CopyStateIdealSelect_ByCellID)[2:3] <- paste0(EachCellID, "_",colnames(HMMcopy_CopyStateIdealSelect_ByCellID)[2:3])    # <-paste0(EachCellID) #  , "_state")
           # print(paste0("dimension of new Cell ID CNV: ", dim(HMMcopy_CopyStateIdealSelect_AllCellID))) # 5404 3

           if(LoopNumb==1) {
             HMMcopy_CopyStateIdealSelect_AllCellID <- HMMcopy_CopyStateIdealSelect_ByCellID
           } else {
             suppressWarnings(suppressMessages( HMMcopy_CopyStateIdealSelect_AllCellID <- dplyr::full_join(HMMcopy_CopyStateIdealSelect_AllCellID, HMMcopy_CopyStateIdealSelect_ByCellID)))
           }
       }
       print(paste0("The number of cells in the ", MySampleID, ": ", LoopNumb))  #  1216
       # dim(HMMcopy_CopyStateIdealSelect_AllCellID) ## 5482 2433

       ## ================ Put "0" after 'chr' to make 'chr01' ==================== ##
       HMMcopy_CopyStateIdealSelect_AllCellID$StartEndChr <- ifelse(grepl("chr\\d{1}_", HMMcopy_CopyStateIdealSelect_AllCellID$StartEndChr), gsub("chr", "chr0",HMMcopy_CopyStateIdealSelect_AllCellID$StartEndChr),
                                                                    HMMcopy_CopyStateIdealSelect_AllCellID$StartEndChr)
       # dim(HMMcopy_CopyStateIdealSelect_AllCellID); HMMcopy_CopyStateIdealSelect_AllCellID[1:2,1:5]  #   5475 2433
       #               StartEndChr 01-040-143929A-R03-C07_copy 01-040-143929A-R03-C07_state 01-040-143929A-R03-C08_copy 01-040-143929A-R03-C08_state
       #                   <char>                       <num>                        <int>                       <num>                        <int>
       # 1: chr01_2000001_2500000                    5.877135                            5                    5.579707                            5

       ## ================  sort the row by Chr#, Start position  ================== ##
       ## HMMcopy_CopyStateIdealSelect_AllCellID_Sort <- HMMcopy_CopyStateIdealSelect_AllCellID[with(HMMcopy_CopyStateIdealSelect_AllCellID, order(StartEndChr)),]  # This is not good.
       HMMcopy_CopyStateIdealSelect_AllCellID_Proc <-  dplyr::mutate(HMMcopy_CopyStateIdealSelect_AllCellID,
                                                                     Chr=sapply(strsplit(HMMcopy_CopyStateIdealSelect_AllCellID$StartEndChr, split="_"),"[", 1 ),
                                                                     Start=sapply(strsplit(HMMcopy_CopyStateIdealSelect_AllCellID$StartEndChr, split="_"),"[", 2 ))
       HMMcopy_CopyStateIdealSelect_AllCellID_Proc$Start <- as.numeric(as.character(HMMcopy_CopyStateIdealSelect_AllCellID_Proc$Start))
       # min(HMMcopy_CopyStateIdealSelect_AllCellID_Proc$Start[HMMcopy_CopyStateIdealSelect_AllCellID_Proc$Chr=="chr01"]) # 1500001
       ## ++++++++  Before sorting
       # View(HMMcopy_CopyStateIdealSelect_AllCellID_Proc[HMMcopy_CopyStateIdealSelect_AllCellID_Proc$Chr=="chr01",c(1:5, 2430:2435)])

       ## +++++++++  After sorting.
       sort(HMMcopy_CopyStateIdealSelect_AllCellID_Proc$Start[HMMcopy_CopyStateIdealSelect_AllCellID_Proc$Chr=="chr01"])[1:5] # [1] 1500001 2000001 3000001 3500001 4000001 This is goal.
       HMMcopy_CopyStateIdealSelect_AllCellID_Sort <- dplyr::arrange(HMMcopy_CopyStateIdealSelect_AllCellID_Proc, Chr, as.numeric(Start)); # dim(HMMcopy_CopyStateIdealSelect_AllCellID_Sort) # 5475 2434  ##

       # HMMcopy_CopyStateIdealSelect_AllCellID_Sort[1:6,c(1:3,2430:2435)]; HMMcopy_CopyStateIdealSelect_AllCellID_Sort[5470:5475,c(1:3,2430:2435)]
       # dim(HMMcopy_CopyStateIdealSelect_AllCellID_Sort);  # 5475 2435 #  Number of samples 1217
       # HMMcopy_CopyStateIdealSelect_AllCellID_Sort[1:10,1:3]; HMMcopy_CopyStateIdealSelect_AllCellID_Sort[2400:2420,1:3]

       ## 1st column to rownames, and remove "Chr" and "Start" column
       HMMcopy_CopyStateIdealSelect_AllCellID_CNRegion <- tibble::column_to_rownames(HMMcopy_CopyStateIdealSelect_AllCellID_Sort, "StartEndChr")
       HMMcopy_CopyStateIdealSelect_AllCellID_CNRegion <- dplyr::select(HMMcopy_CopyStateIdealSelect_AllCellID_CNRegion, -c(Chr, Start))
       # min(HMMcopy_CopyStateIdealSelect_AllCellID_CNRegion); max(HMMcopy_CopyStateIdealSelect_AllCellID_CNRegion) ## 0  # 165.79

      ## transpose the data: CelLID becomes rownames and Region becomes colnames.
       HMMcopy_CopyStateIdeal_Tp <- data.frame(t(HMMcopy_CopyStateIdealSelect_AllCellID_CNRegion) ); #  dim(HMMcopy_CopyStateIdeal_Tp) # 2432 5475

      ### ++++++++++ Store data here by saveRDS
      # saveRDS(HMMcopy_CopyStateIdeal_Tp, file="HMMcopy_CopyStateIdeal_Tp_ConvertNA2.rds")

       ######## =========== ######## =========== ######## =========== ######## ===========
       ## Step4.  Split the data to copy or state.  "CellID_copy" or "CellID_state" is now rownames.
       ######## =========== ######## =========== ######## =========== ######## ===========
       # ## 1) Copy value data.
       # HMMcopy_CopyIdeal_Tp <- dplyr::filter(data.frame(HMMcopy_CopyStateIdeal_Tp), grepl("copy", rownames(HMMcopy_CopyStateIdeal_Tp)));
       # rownames(HMMcopy_CopyIdeal_Tp) <- gsub("_copy", "", rownames(HMMcopy_CopyIdeal_Tp))
        # dim(HMMcopy_CopyIdeal_Tp); HMMcopy_CopyIdeal_Tp[1:2,1:5] # 1216 5475

      ## 2) CN state data.
      HMMcopy_StateIdeal_Tp <- dplyr::filter(data.frame(HMMcopy_CopyStateIdeal_Tp), grepl("state", rownames(HMMcopy_CopyStateIdeal_Tp)));
       rownames(HMMcopy_StateIdeal_Tp) <- gsub("_state", "", rownames(HMMcopy_StateIdeal_Tp))
       # dim(HMMcopy_StateIdeal_Tp); HMMcopy_StateIdeal_Tp[1:2,1:5] # 1216 5475
       rownames(HMMcopy_StateIdeal_Tp) <- paste0("X", gsub("-","_", rownames(HMMcopy_StateIdeal_Tp) )); # tail(rownames(HMMcopy_StateIdeal_Tp) ) # 1216

       ######## =========== ######## =========== ######## =========== ######## ===========
       ## Step5.  Count the number of "NA" in CNV bins and remove outlier bins. Then, convert NA to 2 for the CNV state
       ######## =========== ######## =========== ######## =========== ######## ===========
       # Number_NA <- unlist(apply(HMMcopy_StateIdeal_Tp, 2, function(x) {
       #   # x<- HMMcopy_StateIdeal_Tp[,1]
       #   table(is.na(x))[names(table(is.na(x)))=="TRUE"]
       # }
       # ))
       Number_NA <- unlist(apply(HMMcopy_StateIdeal_Tp, 2, CountNumbNA))
       Number_NA[1:5];length(Number_NA) # plot(sort(Number_NA)) # boxplot(Number_NA)
       # Number_NA[grepl("chr01",names(Number_NA))]
       ## Get 75% quantile of the Number_NA
       # BoxQuantile_TotalCell <-quantile(1:nrow(HMMcopy_StateIdeal_Tp), probs = c(0,0.25,0.5,0.75,1)); print(BoxQuantile_TotalCell)

       ## Outlier bins and remove those bins   - "chrY" bins all have high Gini Coef, so they are removed
       Bin_ToRemoveByGini <- colnames(HMMcopy_StateIdeal_Tp)[Number_NA > nrow(HMMcopy_StateIdeal_Tp)*0.20  ] # Number_NA > BoxQuantile_TotalCell[names(BoxQuantile_NABin)=="25%"]
       length(Bin_ToRemoveByGini); Bin_ToRemoveByGini[1:3]; #Pt5, 52 #223 #  [1] "chr01_1000001_1500000"     "chr01_154500001_155000000" "chr01_156000001_156500000"
       HMMcopy_StateIdeal_RmvArtifactBin <-  dplyr::select(HMMcopy_StateIdeal_Tp, -all_of(Bin_ToRemoveByGini) );  dim(HMMcopy_StateIdeal_RmvArtifactBin) # 1216 5259

       # Number_NA_AftRmvArtifact <- apply(HMMcopy_StateIdeal_RmvArtifactBin, 2, function(x) {
       #   # x<- HMMcopy_StateIdeal_Tp[,1]
       #   table(is.na(x))[names(table(is.na(x)))=="TRUE"]
       # }  )
       Number_NA_AftRmvArtifact <- apply(HMMcopy_StateIdeal_RmvArtifactBin, 2, CountNumbNA)
       Number_NA_AftRmvArtifact[1:5]; # plot(sort(Number_NA_AftRmvArtifact)); boxplot(Number_NA_AftRmvArtifact)

       ## Convert NA to 0 or 2 for the CNV states, after removing artifact bins.
       table(is.na(HMMcopy_StateIdeal_RmvArtifactBin)) # FALSE: 5486067  TRUE: 908877
       HMMcopy_StateIdeal_Tp_NonNA <- HMMcopy_StateIdeal_RmvArtifactBin
       HMMcopy_StateIdeal_Tp_NonNA[is.na(HMMcopy_StateIdeal_Tp_NonNA)] <- 2
       table(is.na(HMMcopy_StateIdeal_Tp_NonNA)) # FALSE: 6394944

       ######## =========== ######## =========== ######## =========== ######## ===========
       ### Step6. Calculate Gini coefficient for all cells - remove outlier cells    get_gini() is from SCOPE package
       ######## =========== ######## =========== ######## =========== ######## ===========
       library(SCOPE)
       GiniCoefficient <- SCOPE::get_gini(t(HMMcopy_StateIdeal_Tp_NonNA)); length(GiniCoefficient); head(GiniCoefficient) # [1]  0.2158 0.1670 0.2868 0.1464 0.1890 0.1928
       #     0%    25%    50%    75%   100%
       # 0.0000 0.1367 0.1705 0.1856 0.4410
       median(GiniCoefficient); mad(GiniCoefficient);
       Median_MAD_GiniCoef <- median(GiniCoefficient)+mad(GiniCoefficient); print(Median_MAD_GiniCoef)  # 0.1705 # 0.02905896 # 0.199559  # 09_025, Median(Gini)+1MAD: 0.1337 # 0.16
       # plot(sort(GiniCoefficient)); abline(h=Median_MAD_GiniCoef, col="red") ; abline(h= median(GiniCoefficient), col="blue")
       # boxplot(GiniCoefficient);  abline(h=Median_MAD_GiniCoef, col="red") ; abline(h= median(GiniCoefficient), col="blue")
       # table(GiniCoefficient >= Median_MAD_GiniCoef ) # FALSE 1042  TRUE 174

       ## Artifactual Cell IDs
       CellID_ToRemoveByGini <- rownames(HMMcopy_StateIdeal_Tp_NonNA)[GiniCoefficient >= Median_MAD_GiniCoef] ##  BoxQuantile_Gini[names(BoxQuantile_Gini)=="50%"]
       length(CellID_ToRemoveByGini); CellID_ToRemoveByGini[1:3]; # pt1 254 #  "X01_040_143929A_R04_C18" "X01_040_143929A_R07_C36" "X01_040_143929A_R10_C11" # 09_025: 430

       ######## =========== ######## =========== ######## =========== ######## ===========
       ## Step7. Let's remove cells (rows) with very high mean state > Median+MAD
       ######## =========== ######## =========== ######## =========== ######## ===========
       StateMean_PerCell <- apply(HMMcopy_StateIdeal_Tp_NonNA, 1, mean); length(StateMean_PerCell)# 1216
       min(StateMean_PerCell); max(StateMean_PerCell) # 0.9592 # 10.149  # 4th smp: 1.377638 # 8.821259
       median(StateMean_PerCell); mad(StateMean_PerCell); # 2.479 # 0.705.
       Median_MAD_MeanState <- median(StateMean_PerCell)+ mad(StateMean_PerCell); print(Median_MAD_MeanState)  # 2.479 # 0.70.  # 3.184
       # plot(sort(StateMean_PerCell)); abline(h=Median_MAD_MeanState, col="red");  abline(h= median(StateMean_PerCell), col="blue")

       # boxplot(StateMean_PerCell);  abline(h=Median_MAD_MeanState, col="red")
       # table(StateMean_PerCell >= Median_MAD_MeanState ) # FALSE 697  TRUE 519

       ## Artifactual Cell IDs by StateMean
       CellID_ToRemoveByMean <- rownames(HMMcopy_StateIdeal_Tp_NonNA)[StateMean_PerCell>= Median_MAD_MeanState ] # BoxQuantile[names(BoxQuantile)=="75%"]
       length(CellID_ToRemoveByMean); CellID_ToRemoveByMean[1:3]; # 519 # "X01_040_143929A_R03_C07" "X01_040_143929A_R03_C08" "X01_040_143929A_R03_C09"

       ## compare Cell IDs between CellID_ToRemovByMean and CellID_ToRemoveByGini
       table(CellID_ToRemoveByMean %in% CellID_ToRemoveByGini) # FALSE: 221, TRUE:284
       CellID_ToRemoveByMeanGini <- unique(c(CellID_ToRemoveByMean, CellID_ToRemoveByGini)); length(CellID_ToRemoveByMeanGini) # Pt1. 511 #  09_025: 651 by Med(Gini)+MAD;  471 by Med(Gini)

       ## ========  Remove artifactual Cells by both Gini Coef. and StateMean  ========= ##
       HMMcopy_StateIdeal_RmvArtifactCNV_ByMeanGini <- dplyr::filter(HMMcopy_StateIdeal_Tp_NonNA, !rownames(HMMcopy_StateIdeal_Tp_NonNA) %in% CellID_ToRemoveByMeanGini)
       dim(HMMcopy_StateIdeal_RmvArtifactCNV_ByMeanGini); HMMcopy_StateIdeal_RmvArtifactCNV_ByMeanGini[1:2,1:3] # Pt1. 698 5073 Pt2. 886 5082  #Pt3. 705 4922. # 09_025: 565 5105
       #                            chr01_1500001_2000000 chr01_2000001_2500000 chr01_3000001_3500000
       # X01_056_143929A_R03_C39                     1                     1                     2

       # saveRDS(HMMcopy_StateIdeal_RmvArtifactCNV_ByMeanGini, file=paste0("HMMcopy_StateIdeal_RmvArtifactCNV_ByMeanGini_", MySampleID,".rds"))

       ######## =========== ######## =========== ######## =========== ######## ===========
       ## Step8. Remove cells by CNA state Low or High: Cells have CNA state 0 or 1 in 20% of total bins (n=5259) and Cells have CNA state ≥ 5 in more than 10 bins
       ######## =========== ######## =========== ######## =========== ######## ===========
       print(paste0("Number of cells after removing artifact by clusters: ", nrow(HMMcopy_StateIdeal_RmvArtifactCNV_ByMeanGini) ) )  #Pt 09_025 565

       ########## ++++++++++ Remove cells that have state=1 in the overall bins > 15% of all bins. On Nov. 15th 2024
       # aa<- t(HMMcopy_StateIdeal_AfterRmvArtifactCell_ColNameProc[473, ])
       # table(aa)

       ### count the number of bins with state=0 or 1.
       Count_LowState01 <- apply(HMMcopy_StateIdeal_RmvArtifactCNV_ByMeanGini, 1, function(x) {
         # x=HMMcopy_StateIdeal_AfterRmvArtifactCell_ColNameProc[19,]
         EachRowTp <- t(x)
         StateTable <- table(EachRowTp) # don't put EachRowTp[,1]
         ## The number of bins with state 0 or 1
         CountState0 <- StateTable[names(StateTable)==0]; #print(CountState0)
         CountState1 <- StateTable[names(StateTable)==1]; # print(CountState1)

         if ( length(CountState0)==0 & length(CountState1)==0 ) {
           CountState0_1 <- 0
         } else if (length(CountState0)!=0 & length(CountState1)==0 ) {
           CountState0_1 <- CountState0
         } else if (length(CountState0)==0 & length(CountState1)!=0) {
           CountState0_1 <- CountState1
         } else if (length(CountState0)!=0 & length(CountState1)!=0) {
           CountState0_1 <- CountState0 + CountState1
         }
       })
       length(Count_LowState01) # total 565 cells

       ### count the number of bins with state > 5
       Count_HighState5_10 <- apply(HMMcopy_StateIdeal_RmvArtifactCNV_ByMeanGini, 1, function(x) {
         # x=HMMcopy_StateIdeal_AfterRmvArtifactCell_ColNameProc[150,]
         EachRowTp <- t(x)
         StateTable <- table(EachRowTp); # print(StateTable) # don't put EachRowTp[,1]
         ## The number of bins with state 5, 6, 7, ~ 10
         CountStateHigh <- sum(StateTable[names(StateTable)>=5]); #print(CountState0)
       })
       length(Count_HighState5_10) # total 565 cells

       ## 10~20% bins of ncol(HMMcopy_StateIdeal_AfterRmvArtifactCell_ColNameProc)
       NumbCellLowState <- ncol(HMMcopy_StateIdeal_RmvArtifactCNV_ByMeanGini)*ThresholdLowState; # ## Originaly, I used "*0.15" print(NumbCellLowState)
       ## Cell IDs to remove
       CellID_LowStateToRmv <- names(Count_LowState01)[Count_LowState01 >= NumbCellLowState];
       print(paste0("number of cell IDs to remov by Low State: ", length(CellID_LowStateToRmv))) # 13

       ## Find cells of High State > 10 bins. # default is 10
       # ThresholdHighState <- 10 #  ncol(HMMcopy_StateIdeal_RmvArtifactCNV_ByMeanGini)*0.15; #  print(NumbCellLowState)
       ## Cell IDs to remove
       CellID_HighStateToRmv <- names(Count_HighState5_10)[Count_HighState5_10 >= ThresholdHighState];
       print(paste0("number of cell IDs to remov by high State: ", length(CellID_HighStateToRmv))) # 54

       ## Calculate rowSums and remove cell of high rowSums > 13000
       RowSum_CNVState <- rowSums(HMMcopy_StateIdeal_RmvArtifactCNV_ByMeanGini)
       RowSum_CNVStateSort <- sort(RowSum_CNVState, decreasing=TRUE);  RowSum_CNVStateSort[1:20]
       CellID_HighRowSumStateToRmv <- names(RowSum_CNVStateSort)[ RowSum_CNVStateSort > 13000]
       print(paste0("Number of cells rowSums > 13000: ", length(CellID_HighRowSumStateToRmv) ))   # 8

       ## remove cells of low or high CNV state, or HighRowSum State
       HMMcopy_StateIdeal_AfterRmvArtifactCell_RmvLowHighState <- HMMcopy_StateIdeal_RmvArtifactCNV_ByMeanGini[!(rownames(HMMcopy_StateIdeal_RmvArtifactCNV_ByMeanGini) %in%
                                                                            c(CellID_LowStateToRmv, CellID_HighStateToRmv, CellID_HighRowSumStateToRmv)),]
       dim(HMMcopy_StateIdeal_AfterRmvArtifactCell_RmvLowHighState) # Pt 09_ 025, 534 5105
       Numb_Cell <- nrow(HMMcopy_StateIdeal_AfterRmvArtifactCell_RmvLowHighState)
       print(paste0("After remove outlier cells of Low or High CNV state: ", nrow(HMMcopy_StateIdeal_AfterRmvArtifactCell_RmvLowHighState))) # 519
       saveRDS(HMMcopy_StateIdeal_AfterRmvArtifactCell_RmvLowHighState, paste0("HMMcopy_StateIdeal_RmvCellByGiniMean","_", MySampleID,"_LowHighState",
                                ThresholdLowState, ".rds"))  # "HMMcopy_StateIdeal_RmvCellByGiniMean_09_025_143855A_LowHighState0.2.rds"

       #### +++++++++++++++++ #### +++++++++++++++++ #### +++++++++++++++++ #### +++++++++++++++++ #### +++++++++++++++++ #### +++++++++++++++++
      #### ++++++++++++ #### ++++++++++++ #### ++++++++++++ #### ++++++++++++ #### ++++++++++++ #### ++++++++++++ #### ++++++++++++
       ###### Process bin names in column names and delete Bin Regions.
      #### +++++++++++ #### +++++++++++ #### +++++++++++ #### +++++++++++ #### +++++++++++ #### +++++++++++ #### +++++++++++ #### +++++++++++
       if(RemoveArtifactCell==FALSE & BinRegionDelete==TRUE) {
             ####  ======== Save rds file before removing High Gini coef and Mean state cell IDs.  =======
             ColNameProcess_BfRmvCell <- gsub("_(.*)","",colnames(HMMcopy_StateIdeal_Tp_NonNA)); length(ColNameProcess_BfRmvCell) # 5475
             ColNameProcess_BfRmvCell[duplicated(ColNameProcess_BfRmvCell)] <- ""

             HMMcopy_StateIdeal_Tp_NonNA_WithoutBinRegion <- HMMcopy_StateIdeal_Tp_NonNA
             colnames(HMMcopy_StateIdeal_Tp_NonNA_WithoutBinRegion) <- ColNameProcess_BfRmvCell
             saveRDS(HMMcopy_StateIdeal_Tp_NonNA_WithoutBinRegion, file=paste0("HMMcopy_StateIdeal_BfRmvCell_NonNA_WithoutBinRegion", MySampleID, ".rds"))
             return(HMMcopy_StateIdeal_Tp_NonNA_WithoutBinRegion)
       } else if (RemoveArtifactCell==FALSE & BinRegionDelete==FALSE){
             return(HMMcopy_StateIdeal_Tp_NonNA)
       } else if (RemoveArtifactCell==TRUE) {
         if(BinRegionDelete==TRUE) {
             ### Process column names, CNV bin names.  Keep just chr01, chr02, .... chrY  One time.  Delete all other column names.
             ColNameProcess <- gsub("_(.*)","",colnames(HMMcopy_StateIdeal_AfterRmvArtifactCell_RmvLowHighState)); length(ColNameProcess) # 5475
             table(duplicated(ColNameProcess)) ## FALSE 23  TRUE 5236
             ColNameProcess[duplicated(ColNameProcess)] <- ""

             colnames(HMMcopy_StateIdeal_AfterRmvArtifactCell_RmvLowHighState) <- ColNameProcess

             saveRDS(HMMcopy_StateIdeal_AfterRmvArtifactCell_RmvLowHighState, file=paste0("HMMcopy_StateIdeal_RmvCellByGiniMeanLowHighState_NoBinRegionAnnot_", MySampleID, ".rds"))
             return(HMMcopy_StateIdeal_AfterRmvArtifactCell_RmvLowHighState)
         } else {
             return(HMMcopy_StateIdeal_AfterRmvArtifactCell_RmvLowHighState)
         }
       }
}


PCAKmeanClustering <- function(Processed_HMMcopyStateData,  MySampleID="MySampleID", KmeanClusterSize=5 ) {
      ##############################################################################################
      ## Methods2,  using PCA and kmean()
      ##############################################################################################
      library("FactoMineR")
      res.pca <- PCA(t(Processed_HMMcopyStateData), graph = FALSE)
      var <- get_pca_var(res.pca)
      var$coord[1:3,1:5]

      set.seed(123)
      res.km <- kmeans(var$coord, centers=KmeanClusterSize, nstart=25)
      grp <- as.factor(res.km$cluster);
      # Color variables by groups
      # fviz_pca_var(res.pca, col.var = grp, palette = c("#0073C2FF", "#EFC000FF", "#868686FF"),\ legend.title = "Cluster")

      ClusterInfo_DF <- data.frame(grp) %>% tibble::rownames_to_column("CellID"); dim(ClusterInfo_DF) # 565 2
      ClusterInfo_Sort <- ClusterInfo_DF[with(ClusterInfo_DF, order(ClusterInfo_DF$grp, decreasing=FALSE)), ];
      table(ClusterInfo_Sort$grp)
      # 1   2   3   4   5
      # 25  16 137  14 344
      ClusterInfo_SortDF <- as.data.frame(table(ClusterInfo_Sort$grp))

      Processed_HMMcopyStateData_RownNameCol <- Processed_HMMcopyStateData %>% data.frame %>% tibble::rownames_to_column("CellID")
      Processed_HMMcopyStateData_ReorderByKmean <- dplyr::inner_join(ClusterInfo_Sort, Processed_HMMcopyStateData_RownNameCol) %>% dplyr::select(-grp) %>%
                              tibble::column_to_rownames("CellID")
      dim(Processed_HMMcopyStateData_ReorderByKmean); Processed_HMMcopyStateData_ReorderByKmean[1:10,1:5] # 534 5105
      Numb_Cell <- nrow(Processed_HMMcopyStateData_ReorderByKmean)

      ### remove column names "Var.#"
      colnames(Processed_HMMcopyStateData_ReorderByKmean) <- gsub("Var(.*)","",colnames(Processed_HMMcopyStateData_ReorderByKmean))

      saveRDS(Processed_HMMcopyStateData_ReorderByKmean, paste0("HMMcopy_StateIdeal_RmvCellByGiniMeanLowHighState_Kmean_", MySampleID,".rds"))

      ##### How to split heatmap by rows as I want.   Lecture: insidehealth.nyumc.org/my.policy
      # # row_split = rep("group1", nrow(Processed_HMMcopyStateData_ReorderByKmean)) # 534
      # row_split[1:16] = "group1"
      # row_split[17:149] = "group2"
      # row_split[150:462] = "group3"
      # row_split[463:500] = "group4"
      # row_split[501:534] = "group5"

      row_split<-c()
      for(EachClusterNumb in 1:KmeanClusterSize ) {
            # EachClusterNumb <- 1
            # print(paste0("current cluster numb: ", EachClusterNumb))
            if(EachClusterNumb==1) {
                  RowStartNumb<-1; #  print(paste0("RowStartNumb: ", RowStartNumb))
                  RowEndNumb <- ClusterInfo_SortDF$Freq[ClusterInfo_SortDF$Var1==EachClusterNumb];  # print(paste0("RowEndNumb: ", RowEndNumb))
            } else {
                  RowStartNumb <- RowEndNumb + 1; #  print(paste0("RowStartNumb: ", RowStartNumb))
                  RowEndNumb <- RowEndNumb + ClusterInfo_SortDF$Freq[ClusterInfo_SortDF$Var1==EachClusterNumb]; #  print(paste0("RowEndNumb: ", RowEndNumb))
            }
            row_split[RowStartNumb:RowEndNumb] = paste0("group", EachClusterNumb)
      }
      table(row_split)

      ### make hierarchical clustering heatmap
      col_fun_CNV = circlize::colorRamp2(c(0:10),
                                         c( "#2EA3DE","#72C5EF","#D5E3EB","#FFCC99","lightsalmon","#FF9933","darkorange","red1","red2","#CC0066", "#99004C")) # # 12 bins  ,"maroon4"
      RowSplitHeatmap_Cluster <- ComplexHeatmap::Heatmap(as.matrix(Processed_HMMcopyStateData_ReorderByKmean), col=col_fun_CNV, show_row_names=FALSE,show_column_names=TRUE,
                                                         cluster_columns=FALSE, cluster_rows=FALSE, row_split=row_split, use_raster=TRUE ) # row_km=5,

      pdf(paste0("CompHmapCNV_RmvByGiniMeanLowHighState_AllChr_", MySampleID,"_", Numb_Cell, "Cells_",KmeanClusterSize,"clu.pdf"), width=10,height=10)  # height=3
                        RowSplitHeatmapDraw <- ComplexHeatmap::draw(RowSplitHeatmap_Cluster) # This is a list
                        ## To save a heatmap to PDF, I need "ComplexHeatmap::draw" Otherwise, I get an error; Error in grid.Call.graphics(C_downviewport, name$name, strict) :
      dev.off()

      ##### +++++++++   =============       ##### +++++++++   ============= ###      ##### +++++++++   ============= ###
      #### Remove cluster with cells of the number < 10
      RowIndexExtract <- ComplexHeatmap::row_order(RowSplitHeatmapDraw); length(RowIndexExtract) # 5  ## This is a list.
      # $group5
      # [1] 501 502 503 504 505 506 507 508 509 510 511 512 513 514 515 516 517 518 519 520 521 522 523 524 525 526 527 528 529 530 531 532 533 534

      ### Count the number of cells in each cluster and remove the cluster.
      CellNumberEachCluster <- unlist(lapply(RowIndexExtract, length)); print(CellNumberEachCluster)
      # group1 group2 group3 group4 group5
      # 16    133    313     38     34

      ## Find which cluster has the elements < 10
      ClusterNumbLowCell <- names(CellNumberEachCluster[CellNumberEachCluster < 10]); print(ClusterNumbLowCell) # ClusterNumbLowCell <- c(5)
      if (length(ClusterNumbLowCell)==0) {
              RowIndexExtract_NoLowCellCluster <- RowIndexExtract
      } else if(length(ClusterNumbLowCell) >=1 ) {
              RowIndexExtract_NoLowCellCluster <- RowIndexExtract[!names(RowIndexExtract) %in% ClusterNumbLowCell]
      }
      fwrite(data.frame(ClusterNumbLowCell), paste0("ClusterNumber_10LowCell_", MySampleID, "_",paste(ClusterNumbLowCell, collapse="_"), ".txt"),col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)
      saveRDS(RowIndexExtract_NoLowCellCluster, file=paste0("RowSplitHeatmapDraw_WithoutArtifactLowHighState_AllChr_",MySampleID,"_",length(RowIndexExtract_NoLowCellCluster),"clu.rds"))  # problem is here.

}


#' Title
#'
#' @param Processed_HMMcopyStateData HMMcopyStateData processed from Step1
#' @param MySampleID  Manually typed in
#' @param AllCell With or without artifact cells
#' @param ChromosomeToSubset  all chromosome
#' @param ClusterSize arbitrary cluster size.
#'
#' @return save RowSplitHeatmapDraw_WithArtifact rds file
#' @export
#'
#' @examples HierarchialClusteringHeatmap(Processed_HMMcopyStateData=HMMcopy_AfterRmvCell_WithoutNABinRegion, MySampleID="01_206_143839A", AllCell="WithoutArtifact", ChromosomeToSubset="All",ClusterSize=9)
HierarchialClusteringHeatmap <- function(Processed_HMMcopyStateData=HMMcopy_StateIdeal_RmvArtifactCNV_WithoutBinRegion,  MySampleID="MySampleID",
                                         AllCell="WithArtifact", ChromosomeToSubset="All", ClusterSize=9) {  # "All", "chr01", "chr02"
  ## Processed_HMMcopyStateData=HMMcopy_StateIdeal_RmvArtifactCNV_Chr1Only;ChromosomeToSubset=FALSE; SampleID="01_040_143929A"; ClusterSize=9;
  col_fun_CNV = circlize::colorRamp2(c(0:10),
                                     c( "#2EA3DE","#72C5EF","#D5E3EB","#FFCC99","lightsalmon","#FF9933","darkorange","red1","red2","#CC0066", "#99004C")) # # 12 bins  ,"maroon4"

  ### 9 clusters in row or I can extract only chr1 across all cells, and make a hierarchical clustering heatmap
  RowSplitHeatmap_Cluster <- ComplexHeatmap::Heatmap(as.matrix(Processed_HMMcopyStateData), col=col_fun_CNV, show_row_names=FALSE,show_column_names=FALSE,
                                                     cluster_columns=FALSE, cluster_rows=TRUE, row_km=ClusterSize, use_raster=TRUE)
  RowSplitHeatmapDraw <- ComplexHeatmap::draw(RowSplitHeatmap_Cluster) # This is a list

  Numb_Cell <- nrow(Processed_HMMcopyStateData)

  if(AllCell=="WithArtifact") {
      if (ChromosomeToSubset=="All") {
          MyWidth=10; MyHeight=10;
          # pdf(paste0("CompHeatmap_CNV_WithArtifact_AllChr_",SampleID,"_", ClusterSize, "Cluster.pdf"), width=MyWidth,height=MyHeight)
          tidyHeatmap::save_pdf(print(RowSplitHeatmap_Cluster), filename=paste0("./CompHeatmap_CNV_WithArtifact_AllChr_",MySampleID,"_", ClusterSize, "Cluster.pdf"), width=MyWidth,height=MyHeight)
          saveRDS(RowSplitHeatmapDraw, file=paste0("RowSplitHeatmapDraw_WithArtifact_AllChr_", MySampleID,".rds"))
      } else if (ChromosomeToSubset=="chr01")  {
          MyWidth=14; MyHeight=10;
          # pdf(paste0("CompHeatmap_CNV_WithArtifact_", ChromosomeToSubset,"_",SampleID,"_", ClusterSize, "Cluster.pdf"), width=MyWidth,height=MyHeight)
          tidyHeatmap::save_pdf(print(RowSplitHeatmap_Cluster), filename=paste0("./CompHeatmap_CNV_WithArtifact_",ChromosomeToSubset, "_",MySampleID,"_", ClusterSize, "Cluster.pdf"), width=MyWidth,height=MyHeight)
          saveRDS(RowSplitHeatmapDraw, file=paste0("RowSplitHeatmapDraw_WithArtifact_",ChromosomeToSubset,"_", MySampleID,".rds"))
      }
  } else if(AllCell=="WithoutArtifact" ) {
    if (ChromosomeToSubset=="All") {
      MyWidth=10; MyHeight=10;
      # pdf(paste0("CompHeatmap_CNV_WithoutArtifact_AllChr_",MySampleID,"_", ClusterSize, "Cluster.pdf"), width=MyWidth,height=MyHeight)
      #      RowSplitHeatmap_Cluster
      # dev.off()
      tidyHeatmap::save_pdf(RowSplitHeatmap_Cluster, filename=paste0("./CompHeatmap_CNV_WithoutArtifact_AllChr_",Numb_Cell,"cells_",MySampleID,"_", ClusterSize, "Cluster.pdf"), width=MyWidth,height=MyHeight)
      saveRDS(RowSplitHeatmapDraw, file=paste0("RowSplitHeatmapDraw_WithoutArtifact_AllChr_", MySampleID,".rds"))
    } else if (ChromosomeToSubset=="chr01")  {
      MyWidth=14; MyHeight=10;
      # pdf(paste0("CompHeatmap_CNV_WithoutArtifact_",ChromosomeToSubset," _",MySampleID,"_", ClusterSize, "Cluster.pdf"), width=MyWidth,height=MyHeight)
      #      RowSplitHeatmapDraw  # This will work. Test   #. RowSplitHeatmap_Cluster <= this doesn't work to save a pdf file.
      # dev.off()
      tidyHeatmap::save_pdf(print(RowSplitHeatmap_Cluster), filename=paste0("./CompHeatmap_CNV_WithoutArtifact_",ChromosomeToSubset,"_",Numb_Cell,"cells_",MySampleID,"_", ClusterSize, "Cluster.pdf"), width=MyWidth,height=MyHeight)
      saveRDS(RowSplitHeatmapDraw, file=paste0("RowSplitHeatmapDraw_WithoutArtifact_",ChromosomeToSubset,"_", MySampleID,".rds"))
    }
  }
  # RowSplitHeatmapDraw <- ComplexHeatmap::draw(RowSplitHeatmap_Cluster);
  # print(RowSplitHeatmapDraw)
  # dev.off()
}



#' Title
#'
#' @param HMMcopy_StateIdeal_AfterRmvArtifactCell_ColNameProc  Output file of Step3. Input for Step4.
#' @param HierarchialClusterInfoFile   cluster numbers for each cell
#' @param ClusterNumbToRmv  cluster number that I want to remove
#' @param ClusterReorder    Cluster reorder number
#' @param ClusterNumbToPick  If I want to pick a couple of clusters.
#' @param ForPloidyPlot  Is RemoveReorderCell_ByClusterNumb is for making heatmap or ploidyplot?
#'
#' @return HMMcopy_StateIdeal_AfterRmvArtifactCell_ColNameProc
#' @export
#'
#' @examples RemoveReorderCell_ByClusterNumb(HMMcopy_StateIdeal_AfterRmvArtifactCell_ColNameProc, HierarchialClusterInfoFile, ClusterNumbToRmv=NULL,ClusterReorder=NULL, ClusterNumbToPick=NULL)
RemoveReorderCell_ByClusterNumb <- function(HMMcopy_StateIdeal_AfterRmvArtifactCell_ColNameProc, HierarchialClusterInfoFile, ClusterNumbToRmv=NULL,
                                            ClusterReorder=NULL, ClusterNumbToPick=NULL, MySampleID, ForPloidyPlot=FALSE) {
    #### Finding cell IDs by cluster number.
    # RowSplitHeatmapDraw_NoRmvCell <- readRDS(HierarchialClusterInfoFile)
    # RowIndexExtract <-ComplexHeatmap::row_order(RowSplitHeatmapDraw_NoRmvCell); sum(length(unlist(RowIndexExtract)))  # "6" "7" "5" "3" "4" "1" "2"
    RowIndexExtract <- readRDS(HierarchialClusterInfoFile)
    if (!ForPloidyPlot) {
        CellNumberEachCluster <- unlist(lapply(RowIndexExtract, length)); print(CellNumberEachCluster)
        fwrite(data.frame(CellNumberEachCluster), paste0("CellNumberEachCluster_", MySampleID,"_",sum(CellNumberEachCluster) ,  "cells.txt"), col.names=TRUE, row.names=TRUE, sep="\t", quote=FALSE)
    }

    IndexCell_ReorderClusterAll <- c();
    if(is.null(ClusterNumbToPick) & !is.null(ClusterReorder)) { ## This is reordering.
        for (EachClusterNumb in ClusterReorder) {
            # EachClusterNumb <- ClusterReorder[1]; print(EachClusterNumb)
            IndexCell_EachCluster <- unlist(RowIndexExtract[names(RowIndexExtract)==as.character( EachClusterNumb)])
            IndexCell_ReorderClusterAll <- c(IndexCell_ReorderClusterAll, IndexCell_EachCluster)
        }
        length(IndexCell_ReorderClusterAll) # 387

        HMMcopy_State_ClusterToReordered <- HMMcopy_StateIdeal_AfterRmvArtifactCell_ColNameProc[IndexCell_ReorderClusterAll,]; dim(HMMcopy_State_ClusterToReordered)  # Pt 09_025, 497 5105
        # HMMcopy_State_Cluster6789 <- HMMcopy_StateIdeal_NoRmvArtifactCell_ColNameProc[c(RowIndexExtract$"6",
        #                              RowIndexExtract$"7",RowIndexExtract$"8",RowIndexExtract$"9"), ]; dim(HMMcopy_State_Cluster6789) # 448 5259

        return(HMMcopy_State_ClusterToReordered)
    } else if (!is.null(ClusterNumbToPick)  & is.null(ClusterNumbToRmv)) {
        HMMcopy_StateIdeal_ClusterToPick_rbind <- data.frame; LoopNumb <- 0;
        for (EachClusterNumb in ClusterNumbToPick) {
          # EachClusterNumb <- ClusterNumbToPick[1]; print(EachClusterNumb)
          LoopNumb <- LoopNumb+1;
          IndexCell_EachCluster <- unlist(RowIndexExtract[names(RowIndexExtract)==as.character(EachClusterNumb)])

          ### Pick up the cells of the cluster
          HMMcopy_StateIdeal_ClusterToPick <- HMMcopy_StateIdeal_AfterRmvArtifactCell_ColNameProc[IndexCell_EachCluster,]; dim(HMMcopy_StateIdeal_ClusterToPick)  #  6 5259
          rownames(HMMcopy_StateIdeal_ClusterToPick) <- paste0(rownames(HMMcopy_StateIdeal_ClusterToPick), "_Cluster", EachClusterNumb)

          ### rbind the picked-up cluster data.
          if(LoopNumb==1) {
            HMMcopy_StateIdeal_ClusterToPick_rbind <-  HMMcopy_StateIdeal_ClusterToPick
          } else {
            HMMcopy_StateIdeal_ClusterToPick_rbind <- rbind(HMMcopy_StateIdeal_ClusterToPick_rbind, HMMcopy_StateIdeal_ClusterToPick)
          }
      } # end of for loop
      return(HMMcopy_StateIdeal_ClusterToPick_rbind)
    } else if (is.null(ClusterNumbToPick) & is.null(ClusterNumbToRmv) & is.null(ClusterReorder)) {
              HMMcopy_StateIdeal_AfterRmvArtifactCell_ColNameProc <- HMMcopy_StateIdeal_AfterRmvArtifactCell_ColNameProc
              return(HMMcopy_StateIdeal_AfterRmvArtifactCell_ColNameProc)
    }
    # else if (is.null(ClusterNumbToPick) & !is.null(ClusterNumbToRmv) & is.null(ClusterReorder)) {  ## to remove cells of ClusterNumbToRmv.  <- I don't need this.
    #         for (EachClusterNumb in ClusterNumbToRmv) {
    #           # EachClusterNumb <- ClusterNumbToRmv[1]; print(EachClusterNumb)
    #           IndexCell_EachCluster <- unlist(RowIndexExtract[names(RowIndexExtract)==as.character(EachClusterNumb)])
    #           IndexCell_ReorderClusterAll <- c(IndexCell_ReorderClusterAll, IndexCell_EachCluster)
    #         }
    #         length(IndexCell_ReorderClusterAll) # 387
    #
    #         HMMcopy_StateIdeal_AfterRmvArtifactCell_ColNameProc <- HMMcopy_StateIdeal_AfterRmvArtifactCell_ColNameProc
    #         return(HMMcopy_StateIdeal_AfterRmvArtifactCell_ColNameProc)
    # }
}

#' Title
#'
#' @param HierarchialClusterInfoFile    cluster number information of each cell
#' @param ClusterReorder    cluster reorder numbers
#'
#' @return IndexCell_ReorderClusterAll
#' @export
#'
#' @examples CountCell_EachCluster(HierarchialClusterInfoFile="RowSplitHeatmapDraw_WithoutArtifactLowHighState7_AllChr_01_213_143839A.rds", ClusterReorder=c(5,3,4,7,6))
CountCell_EachCluster<- function(HierarchialClusterInfoFile, ClusterReorder) {
      #### Finding cell IDs by cluster number.
      # RowSplitHeatmapDraw_NoRmvCell <- readRDS(HierarchialClusterInfoFile)
      # RowIndexExtract <-ComplexHeatmap::row_order(RowSplitHeatmapDraw_NoRmvCell); sum(length(unlist(RowIndexExtract)))  # "6" "7" "5" "3" "4" "1" "2"
      RowIndexExtract <- readRDS(HierarchialClusterInfoFile)


      IndexCell_ReorderClusterAll <- c();
      for (EachClusterNumb in ClusterReorder) {
          # EachClusterNumb <- ClusterReorder[1]; print(EachClusterNumb)
          IndexCell_EachCluster <- unlist(RowIndexExtract[names(RowIndexExtract)==as.character(EachClusterNumb)])
          LengthCell_EachCluster <- length(IndexCell_EachCluster)
          IndexCell_ReorderClusterAll <- c(IndexCell_ReorderClusterAll, LengthCell_EachCluster)
      }
      length(IndexCell_ReorderClusterAll) # [1] 118  31 135  92  11
      return(IndexCell_ReorderClusterAll)
}


#' Title
#'
#' @param HierarchialClusterInfoFile    cluster numbers of each cell
#' @param CNVStateData_AfterRmvHighGiniMeanFile   CNV data after removing cells of high Gini score and high mean.
#' @param MySampleID    Manually typed in
#' @param ClusterNumbToRmv    cluster number to remove
#' @param ClusterReorder      cluster number to reorder
#'
#' @return make CompHmapCNV_RmvByGiniMeanLowHighState_LowSt_Chr1only_ pdf file
#' @export
#'
#' @examples RemoveNormalClusterReorderHeatmap(HierarchialClusterInfoFile="RowSplitHeatmapDraw_WithoutArtifactLowHighState6_AllChr_01_213_143839A.rds",CNVStateData_AfterRmvHighGiniMeanFile="HMMcopy_StateIdeal_RmvCellByGiniMean_WithoutBinRegion_01_213_143939A.rds", MySampleID="01_213_143839A",ClusterNumbToRmv=c(1,2), ClusterReorder=c(3,6,4,5)
RemoveNormalClusterReorderHeatmap <- function (HierarchialClusterInfoFile="RowSplitHeatmapDraw_WithoutArtifactLowHighState6_AllChr_01_213_143839A.rds",
                                      CNVStateData_AfterRmvHighGiniMeanFile="HMMcopy_StateIdeal_RmvCellByGiniMean_WithoutBinRegion_01_213_143939A.rds",
                                       MySampleID,ClusterNumbToRmv=c(1,2), ClusterReorder=c(3,6,4,5)) {
      # HMMcopy_StateIdeal_NoRmvArtifactCell_ColNameProc <- readRDS(HMMcopyFile_BeforeRmvArtifact) ;
      # print(paste0("CNV state after removing by high Gini and Mean state: ", dim(HMMcopy_StateIdeal_NoRmvArtifactCell_ColNameProc))) # 881  5259

      if(is.null(CNVStateData_AfterRmvHighGiniMeanFile)) {
            HMMcopy_StateIdeal_AfterRmvArtifactCell_ColNameProc <- NULL
      } else {
            HMMcopy_StateIdeal_AfterRmvArtifactCell_ColNameProc <- readRDS(CNVStateData_AfterRmvHighGiniMeanFile)
      }
      dim(HMMcopy_StateIdeal_AfterRmvArtifactCell_ColNameProc) # Pt#3, 01-206, 704 4922;  Pt#8 09_025, 536 5105 # Pt9 BM_CD138, 584 5051

      ########## ++++++++++ Remove artifacts by cluster number of cells, and reorder cells by cluster number
      HMMcopy_StateIdeal_AfterRmvArtifactCell_RmvReorderCluster <-  RemoveReorderCell_ByClusterNumb(HMMcopy_StateIdeal_AfterRmvArtifactCell_ColNameProc,
                                     HierarchialClusterInfoFile, ClusterNumbToRmv, ClusterReorder, ClusterNumbToPick=NULL, MySampleID, ForPloidyPlot=FALSE)
      dim(HMMcopy_StateIdeal_AfterRmvArtifactCell_RmvReorderCluster) # Pt#3, 01-206, 651 4922;  Pt#8 09_025, 497 5105   # Pt9 BM_CD138, 475 5051
      saveRDS(HMMcopy_StateIdeal_AfterRmvArtifactCell_RmvReorderCluster,
              file=paste0("HMMcopy_StateIdeal_AfterRmvArtifactCell_RmvReorderCluster_",MySampleID,".rds"))

      ##### Count the number of bins in each chromosome
      Position_chrStart <- grep("chr", colnames(HMMcopy_StateIdeal_AfterRmvArtifactCell_RmvReorderCluster)); length(Position_chrStart) ; # print(Position_chrStart) ## 22
      # [1]    1  416  877 1262 1629 1969 2298 2597 2834 3050 3300 3555 3808 3992 4163 4315 4457 4603 4747 4854 4969 5027 5089 5364
      EndChr1Position <- (Position_chrStart[2]-1); print(EndChr1Position) # 415
      EndChr1PositionPercent <- (Position_chrStart[2]-1)/ncol(HMMcopy_StateIdeal_AfterRmvArtifactCell_RmvReorderCluster); print(EndChr1PositionPercent)

      EndChrYPosition <- Position_chrStart[22]-1; print(EndChrYPosition) # 5363
      EndChrYPositionPercent <- (Position_chrStart[22]-1)/ncol(HMMcopy_StateIdeal_AfterRmvArtifactCell_RmvReorderCluster); print(EndChrYPositionPercent)

      ##### Count the number of cells in each cluster.
      NumberCell_EachCluster <- CountCell_EachCluster(HierarchialClusterInfoFile, ClusterReorder) # [1] 118  31 135  92  11
      Numb_Cell <- sum(NumberCell_EachCluster)

      ## Split cells in hierarchical clustering
      if(length(ClusterReorder)>=10) {
            MySplit<-c(paste('clone0', c( rep(1:9, times=NumberCell_EachCluster[1:9])), sep="" ),
                       paste('clone',c(rep(10:length(ClusterReorder), times=NumberCell_EachCluster[10:length(ClusterReorder)] )), sep="" ) )
      } else {
            MySplit<-paste('clone', c( rep(1:length(ClusterReorder), times=NumberCell_EachCluster)) )
      }

      col_fun_CNV = circlize::colorRamp2(c(0:10),
                                         c( "#2EA3DE","#72C5EF","#D5E3EB","#FFCC99","lightsalmon","#FF9933","darkorange","red1","red2","#CC0066", "#99004C"))

      table(rowSums(HMMcopy_StateIdeal_AfterRmvArtifactCell_RmvReorderCluster) == "0")
      # sort(rowSums(HMMcopy_StateIdeal_AfterRmvArtifactCell_RmvReorderCluster))[1:5]

      ### Save heatmap to PDF file
      pdf(paste0("CompHmapCNV_RmvByGiniMeanLowHighState_LowSt_AllChr_", MySampleID,"_",Numb_Cell, "cells_RordClone.pdf"), width=10,height=10)  # height=3
              ##### Split heatmap cells by cluster number in horizontal, and by chromosome in vertical
              RowSplitHeatmap_ReorderClone <- ComplexHeatmap::Heatmap(as.matrix(HMMcopy_StateIdeal_AfterRmvArtifactCell_RmvReorderCluster), name="CloneHeatmap",
                                                                      col=col_fun_CNV, show_row_names=FALSE,show_column_names=FALSE, cluster_columns=FALSE, cluster_rows=FALSE,
                                                                      split=MySplit, use_raster=TRUE)
              # split=paste('clone', c( rep(1:length(ClusterReorder), times=NumberCell_EachCluster)) )
              RowSplitHeatmapDraw_ReorderClone <- ComplexHeatmap::draw(RowSplitHeatmap_ReorderClone)
              ComplexHeatmap::decorate_heatmap_body("CloneHeatmap", {
                # for (i in c(4,7)){  # This works.
                #     # grid::grid.lines(c(4/10, 4/10), c(-7.9, 1), gp = grid::gpar(col = "black", lwd = 1))  # (-7.9, 1) is end to end in vertical
                #     # grid::grid.lines(c(7/10, 7/10), c(0, 1), gp = grid::gpar(col = "red", lwd = 2))
                #     grid::grid.lines(c(i/10, i/10), c(-7.9, 1), gp = grid::gpar(col = "black", lwd = 1))
                # }
                for (i in 1:(length(Position_chrStart) +1) ) {
                  if (i==1) {
                      # EndChrPositionPercent <- 1/ncol(HMMcopy_StateIdeal_AfterRmvArtifactCell_RmvReorderCluster)
                      next
                  } else if (i==22) {
                      # EndChrPositionPercent <-ncol(HMMcopy_StateIdeal_AfterRmvArtifactCell_RmvReorderCluster)/ncol(HMMcopy_StateIdeal_AfterRmvArtifactCell_RmvReorderCluster)
                      next
                  } else {
                      EndChrPositionPercent <- (Position_chrStart[i]-1)/ncol(HMMcopy_StateIdeal_AfterRmvArtifactCell_RmvReorderCluster); # print(EndChrPositionPercent)
                  }
                  # grid::grid.lines(c( EndChrPositionPercent, EndChrPositionPercent), c(-nrow(HMMcopy_StateIdeal_AfterRmvArtifactCell_RmvReorderCluster)/4, 1), gp = grid::gpar(col = "white", lwd = 1)) # c(-20, 1) <= c(bot, top)
                   grid::grid.lines(x=c( EndChrPositionPercent, EndChrPositionPercent), y=c(-45, 1), gp = grid::gpar(col = "white", lwd = 1)) # c(-20, 1) <= c(bot, top)
                }
              } )

      dev.off()

      # RowSplitHeatmap_ReorderClone <- ComplexHeatmap::Heatmap(as.matrix(HMMcopy_StateIdeal_AfterRmvArtifactCell_RmvReorderCluster),
      #                                                         col=col_fun_CNV, show_row_names=FALSE,show_column_names=TRUE, cluster_columns=FALSE, cluster_rows=FALSE,
      #                                                         split=MySplit, use_raster=TRUE)
      # RowSplitHeatmapDraw_ReorderClone <- ComplexHeatmap::draw(RowSplitHeatmap_ReorderClone)

      #########################      ##########################      ##########################
      # options(expressions = 500000)
      # saveRDS(RowSplitHeatmapDraw_ReorderClone, file=paste0("RowSplitHeatmapDraw_WithoutArtifactLowHighState_AllChr_", MySampleID,"_RordClone.rds"))
      # Error: C stack usage  7956768 is too close to the limit
      ##########################      ##########################      ##########################

      ### Subset only chr1 and make a heatmap
      # RowIndexExtract <-ComplexHeatmap::row_order(RowSplitHeatmapDraw_ReorderClone); sum(length(unlist(RowIndexExtract))) # 387 # clone1, 2, 3,.. 5

      ###### =========== ########### ============== ############## ================ ###################### ===========
      ## ++++++++= Select only Chromosome1 bins and make a heatmap
      BinPosition_Chr02 <- grep("chr02", colnames(HMMcopy_StateIdeal_AfterRmvArtifactCell_RmvReorderCluster)); print(BinPosition_Chr02) # 371
      HMMcopy_StateIdeal_AfterRmvArtifactCell_RmvReorderCluster_Chr01 <- HMMcopy_StateIdeal_AfterRmvArtifactCell_RmvReorderCluster[, 1:(BinPosition_Chr02-1)];
      dim(HMMcopy_StateIdeal_AfterRmvArtifactCell_RmvReorderCluster_Chr01) # 387 370

      RowSplitHeatmap_ReorderClone_Chr1 <- ComplexHeatmap::Heatmap(as.matrix(HMMcopy_StateIdeal_AfterRmvArtifactCell_RmvReorderCluster_Chr01), col=col_fun_CNV,
                                                              show_row_names=FALSE,show_column_names=FALSE, cluster_columns=FALSE, cluster_rows=FALSE,
                                                              split=MySplit,
                                                              use_raster=TRUE)

      pdf(paste0("CompHmapCNV_RmvByGiniMeanLowHighState_LowSt_Chr1only_", MySampleID,"_",Numb_Cell, "cells_RordClone.pdf"), width=10,height=10)  # height=3
            RowSplitHeatmapDraw_ReorderClone <- ComplexHeatmap::draw(RowSplitHeatmap_ReorderClone_Chr1)
      dev.off()
}




#' Title Chr1Centromere_Heatmap
#'
#' @param HierarchialClusterInfoFile   HMMcopy data after remove artifact cells
#' @param HMMcopy_StateIdeal_RmvArtifactCNV_ByMeanGiniStateLowHighFile    HMMcopy data after remove artifact cells and after reordered
#' @param MySampleID   Manually typed in.
#' @param ClusterNumbToRmv Remove n oise clusters
#' @param ClusterReorder Reorder clusters by clone order.
#'
#' @return CompHeatmap_Chr1CentromereRegion_ pdf file
#' @export
#'
#' @examples Chr1Centromere_Heatmap(HierarchialClusterInfoFile="HMMcopy_StateIdeal_RmvArtifactCNV_ByMeanGini_01_213_143839A.rds", HMMcopy_StateIdeal_RmvArtifactCNV_ByMeanGiniStateLowHighFile="HMMcopy_StateIdeal_AfterRmvArtifactCell_RmvReorderCluster_01_213_143839A.rds", MySampleID="01_213_143839A")
Chr1Centromere_Heatmap <- function (HierarchialClusterInfoFile, HMMcopy_AfterRmvCell_WithoutNABinRegionFile, HMMcopy_StateIdeal_RmvArtifactCNV_ByMeanGiniStateLowHighFile, MySampleID) {

      HierarchialClusterInfo <- readRDS(HierarchialClusterInfoFile)
      class(HierarchialClusterInfo) #list

      HMMcopy_StateIdeal_RmvArtifactCNV_ByMeanGini_StateLowHighReorder <- readRDS(HMMcopy_StateIdeal_RmvArtifactCNV_ByMeanGiniStateLowHighFile)
      HMMcopy_RmvArtifactCNV_ByMeanGini_StateLowHighReorder_CellIDColumn <- HMMcopy_StateIdeal_RmvArtifactCNV_ByMeanGini_StateLowHighReorder %>% tibble::rownames_to_column("CellID")
      dim(HMMcopy_RmvArtifactCNV_ByMeanGini_StateLowHighReorder_CellIDColumn) #Pt3, 01-206, 704 4923 # Pt#5.  314 5259 #Pt8 09_025;  536 5105

      ## HMMcopyData After low quality cluster cell and sorted by clone order. - This cell order is determined by K-Means clustering
      HMMcopy_AfterRmvCell_WithoutNABinRegion<- readRDS(HMMcopy_AfterRmvCell_WithoutNABinRegionFile); dim(HMMcopy_AfterRmvCell_WithoutNABinRegion) # 651 4922
      HMMcopy_AfterRmvCell_WithoutNABinRegion_CellIDColumn <- HMMcopy_AfterRmvCell_WithoutNABinRegion %>% data.frame %>% tibble::rownames_to_column("CellID") %>% dplyr::mutate(FakeCol="FakeCol")
      dim(HMMcopy_AfterRmvCell_WithoutNABinRegion_CellIDColumn) # 651 4923

      # ########## ++++++++++ Remove artifacts by cluster number of cells, and reorder cells by cluster number
      HMMcopy_StateIdeal_RmvbyGiniMeanStateLowHigh_Reorder <- dplyr::inner_join(HMMcopy_AfterRmvCell_WithoutNABinRegion_CellIDColumn[, c("CellID","FakeCol")],
                                                          HMMcopy_RmvArtifactCNV_ByMeanGini_StateLowHighReorder_CellIDColumn) %>% dplyr::select(-FakeCol) %>% tibble::column_to_rownames("CellID")
      dim(HMMcopy_StateIdeal_RmvbyGiniMeanStateLowHigh_Reorder) # 651 4922


      ### Checking the reorderd cell IDs
      # identical(rownames(HMMcopy_StateIdeal_RmvbyGiniMeanStateLowHigh_Reorder), rownames()) # TRUE


      #### +++++++++++ #### +++++++++++ #### +++++++++++ #### +++++++++++ #### +++++++++++ #### +++++++++++ #### +++++++++++ #### +++++++++++
      #### +++++++++++ Find "233" "234" "235" "236” bins in Chr1 for #8, 09_025_143855A ++++++++++++ #### ++++++++++++ #### ++++++++++++
      ##### Old style
      FindChr1Bin <- table(grepl("chr01", colnames(HMMcopy_StateIdeal_RmvbyGiniMeanStateLowHigh_Reorder))); print(FindChr1Bin) #Pt4. TRUE: 370 #Pt5, 415  ### Pt8. TRUE: 426,  Pt9: 348
      Chr1BinNumb <- FindChr1Bin[names(FindChr1Bin)=="TRUE"]; print(Chr1BinNumb) # Pt4. 370 #Pt5. 415  # Pt8. 426

      ColMeanChr1 <- colMeans(HMMcopy_StateIdeal_RmvbyGiniMeanStateLowHigh_Reorder[,1:Chr1BinNumb])
      names(ColMeanChr1) <- colnames(HMMcopy_StateIdeal_RmvbyGiniMeanStateLowHigh_Reorder)[1:Chr1BinNumb]
      # plot(ColMeanChr1); abline(h=1.4, col="red")

      names(ColMeanChr1)[ColMeanChr1<2]; # Pt4, "chr01_103500001_104000000" "chr01_119500001_120000000"
      ## Pt8.  ColMeanChr1<1.5. [1] "chr01_145500001_146000000" "chr01_146500001_147000000" "chr01_147000001_147500000" "chr01_147500001_148000000"
      which(colnames(HMMcopy_StateIdeal_RmvbyGiniMeanStateLowHigh_Reorder)=="chr01_145500001_146000000")
      ## start point of 1q  #Pt4. "chr01_119500001_120000000", 232 #Pt8. chr01_145500001_146000000, 233 $ Pt9: 234

      NumbState1_AllBin <- c();
      for (i in 1:Chr1BinNumb) {
          # i <- 1
            TableState <- table(HMMcopy_StateIdeal_RmvbyGiniMeanStateLowHigh_Reorder[, i])
            if(any(names(TableState)==1)) {
              NumbState1 <- TableState[names(TableState)==1]
          } else {
              NumbState1 <- 0
          }
          NumbState1_AllBin <- c(NumbState1_AllBin, NumbState1)
      }
      names(NumbState1_AllBin) <- colnames(HMMcopy_StateIdeal_RmvbyGiniMeanStateLowHigh_Reorder)[1:Chr1BinNumb]
      sort(NumbState1_AllBin, decreasing=TRUE)[1:5]
      # chr01_145500001_146000000 chr01_146500001_147000000 chr01_119500001_120000000 chr01_119000001_119500000 chr01_118500001_119000000  # Pt4.
      # 118                        86                        67                        47                        42
      StartPoint1q <- which(colnames(HMMcopy_StateIdeal_RmvbyGiniMeanStateLowHigh_Reorder)=="chr01_145500001_146000000") ## start point of 1q # Pt8: 233
      ColName_226_247 <- colnames(HMMcopy_StateIdeal_RmvbyGiniMeanStateLowHigh_Reorder)[(StartPoint1q-7):(StartPoint1q+14)]; print(ColName_226_247) # [1] "chr01_119000001_119500000" "chr01_119500001_120000000" "chr01_145500001_146000000" "chr01_146500001_147000000"
      # plot(NumbState1_AllBin); abline(h=500, col="red")

      ## Let's see heatmap
      Chr1_HMMcopy_StateIdeal <- HMMcopy_StateIdeal_RmvbyGiniMeanStateLowHigh_Reorder[1:Chr1BinNumb]; dim(Chr1_HMMcopy_StateIdeal) # 637 370
      colnames(Chr1_HMMcopy_StateIdeal) <- c(1:(StartPoint1q-8), ColName_226_247, (StartPoint1q+15):280)

      col_fun_CNV = circlize::colorRamp2(c(0:10),
                                         c( "#2EA3DE","#72C5EF","#D5E3EB","#FFCC99","lightsalmon","#FF9933","darkorange","red1","red2","#CC0066", "#99004C"))
      RowSplitHeatmap_Chr1Only <- ComplexHeatmap::Heatmap(as.matrix(Chr1_HMMcopy_StateIdeal[, 225:250]), col=col_fun_CNV,
                                                          show_row_names=FALSE,show_column_names=TRUE, cluster_columns=FALSE,
                                                          cluster_rows=FALSE, use_raster=TRUE,  column_names_rot=90)
      pdf(paste0("CompHeatmap_Chr1CentromereRegion_", MySampleID,"_",nrow(Chr1_HMMcopy_StateIdeal),"cells.pdf"), width=6,height=8)  # height=3
              RowSplitHeatmapDraw_SubsetCluster <- ComplexHeatmap::draw(RowSplitHeatmap_Chr1Only)
      dev.off()
}


################################################################### Ploidy plot #####################################################################

#' Title
#'
#' @param karyoplot karyoplot obejct
#' @param snps karyoplot tiles
#' @param lrr.column lrr
#' @param labels LRR
#' @param ymin -4
#' @param ymax 2
#' @param out.of.range  points
#' @param out.of.range.col white
#' @param density.height 0.05
#' @param density.window 1e+05
#' @param line.at.0   TRUE
#' @param line.at.0.col   blue
#' @param r0  0
#' @param r1  1
#' @param points.cex  0.3
#' @param points.col  #333333
#' @param points.pch  16
#' @param label.cex   1.5
#' @param label.srt   90
#' @param label.margin  0.03
#' @param add.axis    TRUE
#' @param axis.cex    1.2
#' @param track.margin  0.1
#' @param data.panel  1
#' @param verbose   Boolean
#'
#' @returns NULL
#' @export
#'
#' @examples   plotLRR(karyoplot, snps, lrr.column, labels, ymin, ymax, out.of.range, out.of.range.col,density.height, density.window, line.at.0, line.at.0.col, r0, r1, points.cex, points.col, points.pch, label.cex, label.srt, label.margin, add.axis, axis.cex,track.margin, data.panel, verbose)
plotLRR <- function (karyoplot, snps, lrr.column = "lrr", labels = "LRR", ymin = -4, ymax = 2, out.of.range = "points", out.of.range.col = "white",
                     density.height = 0.05, density.window = 1e+05, line.at.0 = TRUE, line.at.0.col = "blue", r0 = 0, r1 = 1, points.cex = 0.3,
                     points.col = "#333333", points.pch = 16, label.cex = 1.5, label.srt = 90, label.margin = 0.03, add.axis = TRUE, axis.cex = 1.2,
                     track.margin = 0.1, data.panel = 1, verbose = FALSE)  {
      if (!methods::is(karyoplot, "KaryoPlot"))
        stop("karyoplot must be a KaryoPlot object")

      out.of.range <- match.arg(out.of.range, c("points", "density"))

      if (is.character(lrr.column)) {
        if (!any(names(GenomicRanges::mcols(snps)) == lrr.column)) {
          stop("The lrr.column (", lrr.column, ") has not been found in the data")
        }
      } else {
        if (!(is.numeric(lrr.column) && as.integer(lrr.column) ==
              lrr.column && lrr.column > 0 && lrr.column <= ncol(mcols(snps)))) {
          stop("lrr.column must be either the name of a column or an integer between 1 and the number of metadata columns in snps")
        }
      }
      snps <- removeNAs(snps, lrr.na = TRUE, baf.na = FALSE, id.na = FALSE, verbose = FALSE)


      below.min <- GenomicRanges::mcols(snps)[, lrr.column] < ymin
      above.max <- GenomicRanges::mcols(snps)[, lrr.column] > ymax
      if (any(!(below.min | above.max))) { # This step makes the main plot.  <<<<============
        karyoploteR::kpPoints(karyoplot, data = snps[!(below.min |
                                                         above.max)], y = GenomicRanges::mcols(snps)[, lrr.column][!(below.min | above.max)], ymin = ymin, ymax = ymax, r0 = r0, r1 = r1,
                              col = points.col, pch = points.pch, cex = points.cex)
      }

      if (line.at.0 == TRUE) { # TRUE  ## This makes x-axis blue line.
        karyoploteR::kpAbline(karyoplot, h = 0, ymin = ymin,  ymax = ymax, r0 = r0, r1 = r1, col = line.at.0.col)
        karyoploteR::kpAbline(karyoplot, h = 0, ymin = 0.5,   ymax = ymax, r0 = r0, r1 = 0.5, col = "red")
      }
      invisible(karyoplot)
}



#' Title
#'
#' @param CNVSeg_StateIdealSelect_RmvByGiniMeanState_RmvByCluster CNA state data after removing outliers by High Gini, high mean and cluster number
#' @param TotalSumThreshold 30
#'
#' @returns NULL
#' @export
#'
#' @examples  CNVStatePloidyLine(CNVSeg_StateIdealSelect_RmvByGiniMeanState_RmvByCluster, TotalSumThreshold=30)
CNVStatePloidyLine <- function (CNVSeg_StateIdealSelect_RmvByGiniMeanState_RmvByCluster, TotalSumThreshold=30 ) {
    AllChr <- unique(CNVSeg_StateIdealSelect_RmvByGiniMeanState_RmvByCluster$chr); print(AllChr); length(AllChr)
    ## [1] "chr1"  "chr10" "chr11" "chr12" "chr13" "chr14" "chr1

    for (EachChr in AllChr) {
        # EachChr <- AllChr[1]; print(EachChr)
        # if(EachChr=="chrY") next;
        print(paste0("current chromosome: ", EachChr))
        CNVSeg_RmvByGiniMeanState_Subset <- dplyr::filter(CNVSeg_StateIdealSelect_RmvByGiniMeanState_RmvByCluster, chr==EachChr);
        dim(CNVSeg_RmvByGiniMeanState_Subset); CNVSeg_RmvByGiniMeanState_Subset[1,] # 3091    6
        #         chr start     end state   median                 cell_id
        #       <char> <int>   <int> <int>    <num>                  <char>
        #   1:   chr1     1 2000000     2 2.068097 X01_213_143839A_R04_C60

        ## Find major ploidy number
        FindPloidy <-  table(CNVSeg_RmvByGiniMeanState_Subset$state); print(FindPloidy) # Major ploidy should have cells > 600
        # 1    2    3    4    5    6
        # 684 1052  688  303    9    5

        # MajorPloidy <- names(FindPloidy)[FindPloidy > 300] ; print(paste0("Major Ploidy: ", MajorPloidy)) #"1", "2", "3"
        # BotPloidy <-  min(MajorPloidy); TopPloidy <- max(MajorPloidy);
        # print(paste0("bottom Ploidy: ", BotPloidy)); print(paste0("Top Ploidy: ", TopPloidy)) # "1", "3"
        # TopPloidy <- names(FindPloidy)[FindPloidy==max(FindPloidy[FindPloidy>600])] ; print(paste0("Top Ploidy: ", TopPloidy)) # "2"
        # BotPloidy <- names(FindPloidy)[FindPloidy==min(FindPloidy[FindPloidy>600])] ; print(paste0("bottom Ploidy: ", BotPloidy)) # "1"

        ### Indicate ploidy for each segment with only major ploidy
        # CNVSeg_RmvByGiniMeanState_Subset$Ploidy <- ifelse(CNVSeg_RmvByGiniMeanState_Subset$state %in%
        #                                               c(BotPloidy:TopPloidy),CNVSeg_RmvByGiniMeanState_Subset$state, 1 )

        CNVSeg_RmvByGiniMeanState_Subset[1,]
        #       chr start     end state   median                 cell_id Ploidy
        #     <char> <int>   <int> <int>    <num>                  <char>  <num>
        # 1:   chr1     1 2000000     2 2.068097 X01_213_143839A_R04_C60      2

        ### New column of chr_start_end
        CNVSeg_RmvByGiniMeanState_Subset_chrstrend <- dplyr::mutate(CNVSeg_RmvByGiniMeanState_Subset,chrstrend=paste0(chr,"_",start,"_",end)); CNVSeg_RmvByGiniMeanState_Subset_chrstrend[1:2,]

        ## Find unique start positions
        UniqueStartPosition <- unique(CNVSeg_RmvByGiniMeanState_Subset_chrstrend$start); #  print(UniqueStartPosition)
        length(UniqueStartPosition); tail(UniqueStartPosition) # 434 # [1] 244500001 245000001 245500001 246000001 247000001 247500001

        ## Calculate multiplication of unique(ploidy) in each segment.  1p end 123400000
        MajorMinorPloidy <- data.frame(); MajorMinorPloidy_All <- data.frame()
        for(EachStartPosition in UniqueStartPosition) {
          # EachStartPosition <- UniqueStartPosition[150]; 24000001 # 176500001 #  150000001 # UniqueStartPosition[10] # 119500001, this is chr1p state=1 is main starting

          if(EachChr=="chr1" & EachStartPosition==1) next;

          StartPos_Major<-as.integer(); EndPos_Major<-as.integer(); StartPos_Minor<-as.integer(); EndPos_Minor<-as.integer();
          MajorPloidy<-0; MinorPloidy<-0; MinorPloidy<-0; EachMinorPloidy<-0; MajorPloidyDF <- data.frame(); MinorPloidyDF <- data.frame()

          ## subset segments with the same start position
          CNVSeg_RmvByGiniMeanState_Subset_chrstrend_Subset <- dplyr::filter(CNVSeg_RmvByGiniMeanState_Subset_chrstrend, state %in% c(1:4), start==EachStartPosition )

          ## Find the major ploidy in the segment
          PloidyTable <- table(CNVSeg_RmvByGiniMeanState_Subset_chrstrend_Subset$state); #  print(paste0("current ploidy: ", PloidyTable))
          # print(PloidyTable); print(paste0("sum of ploidy values: ", sum(PloidyTable)))

          if(sum(PloidyTable) <= TotalSumThreshold) next;

          if(length(PloidyTable)==1 & PloidyTable[1] >= 30) {  # & PloidyTable[1] < 2
              MajorPloidy <- as.integer(names(PloidyTable[1]))
          } else if(length(PloidyTable)==2 & PloidyTable[1]==PloidyTable[2] ) {
              MajorPloidy <- as.numeric(names(PloidyTable[1]))
          } else if( (length(PloidyTable)==3 & PloidyTable[1]==PloidyTable[2]) | (length(PloidyTable)==3 & PloidyTable[2]==PloidyTable[3]) |
                     (length(PloidyTable)==3 & PloidyTable[1]==PloidyTable[3])  )  {
            MajorPloidy <- as.integer(min(names(PloidyTable)[ PloidyTable == max(PloidyTable)])  )
            if(min(PloidyTable) >= 30 ) {
                MinorPloidy <- as.integer(min(names(PloidyTable)[ PloidyTable == min(PloidyTable)])  )
            } else{
                MinorPloidy <- 0
            }
          } else if( (length(PloidyTable)==4 & PloidyTable[1]==PloidyTable[2]) | (length(PloidyTable)==4 & PloidyTable[2]==PloidyTable[3]) |
                     (length(PloidyTable)==4 & PloidyTable[3]==PloidyTable[4]) | (length(PloidyTable)==4 & PloidyTable[1]==PloidyTable[3])    ) {
              MajorPloidy <- as.integer(min(names(PloidyTable)[ PloidyTable == max(PloidyTable)])  )
          } else  {
              MajorPloidy <- as.integer(names(PloidyTable)[ PloidyTable == max(PloidyTable) ]  ); print(paste0("Major ploidy: ", MajorPloidy))
              CNVSeg_RmvByGiniMeanState_Subset_chrstrend_MajorPloidySubset <- dplyr::filter(CNVSeg_RmvByGiniMeanState_Subset_chrstrend_Subset, state==MajorPloidy)
              # MinorPloidy <- (median( CNVSeg_RmvByGiniMeanState_Subset_chrstrend_MajorPloidySubset$median ) )

            ## To get MinorPloidy, only when there is MinorPloidy.
            if(length(PloidyTable) > 1 & sort(PloidyTable, decreasing=TRUE)[2] >= 20 ) {  ## I used 10 for first patient.
                MinorPloidy <- as.integer(names(PloidyTable)[ PloidyTable== sort(PloidyTable, decreasing=TRUE)[2]  ]  ); print(paste0("Minor ploidy: ", MinorPloidy))
                CNVSeg_RmvByGiniMeanState_Subset_chrstrend_MinorPloidySubset <- dplyr::filter(CNVSeg_RmvByGiniMeanState_Subset_chrstrend_Subset, state %in% MinorPloidy)
            } else if( length(PloidyTable) <= 1 | sort(PloidyTable, decreasing=TRUE)[2] < 20  ) {  ## I used 10 for first patient.
                MinorPloidy<-0
            }
          }

          print(paste0("major ploidy in the segment: ", MajorPloidy))
          print(paste0("minor ploidy in the segment: ", MinorPloidy))

          if(  (is.na(MajorPloidy[1]) | MajorPloidy==0) & (MinorPloidy==0 )  ) {
              break;
          }

          if(sum(PloidyTable) > 50) {
              print(paste0("Ploidy 1, start: ", StartPos_Major))   # 150000001
              print(paste0("Ploidy 1, end: ", EndPos_Major))
              # break;
          };

          if(MajorPloidy > 0) {
              # if(MinorPloidy ==4) {
              #       print(paste0("StartPos_Major:", StartPos_Major))
              #       print(paste0("EndPos_Major:", EndPos_Major))
              #       print("MinorPloidy 4 break"); #break;
              # }

              ## MajorPloidy
              StartPos_Major <- EachStartPosition; print(StartPos_Major) # as.numeric(sapply(strsplit(EachSegment, split="_"), "[", 2)); print(StartPos)

              CNVSeg_RmvByGiniMeanState_Subset_chrstrend_MajorPloidySubset <- dplyr::filter(CNVSeg_RmvByGiniMeanState_Subset_chrstrend_Subset, state==MajorPloidy)
              EndPos_Major <-  CNVSeg_RmvByGiniMeanState_Subset_chrstrend_MajorPloidySubset$end[nrow(CNVSeg_RmvByGiniMeanState_Subset_chrstrend_MajorPloidySubset)];  print(EndPos_Major)

              r0 = 0; r1 = 1;
              karyoplot<-kp
              ymin=0; ymax=5
              chr <- EachChr
              x0 <- StartPos_Major; #   start(karyoplot$genome[chr]) # start
              x1 <- EndPos_Major # end(karyoplot$genome[chr])   # end
              h_major <- MajorPloidy; print(h_major)
              kpSegments(karyoplot = karyoplot, chr = chr, x0 = x0,   ### This make a horizontal line. I need X0, X1, and h
                         x1 = x1, y0=h_major, y1=h_major, ymin=ymin, ymax=ymax,
                         r0 = r0, r1 = r1, data.panel = data.panel, clipping = clipping, col="green", lwd=7)

              ### collect MajorPloidy information.
              MajorPloidyDF <- data.frame(t( c(EachChr, StartPos_Major,EndPos_Major, h_major, "Major" ))) ;
              rownames(MajorPloidyDF) <- paste0("Major_",StartPos_Major,"_",EndPos_Major); colnames(MajorPloidyDF) <- c("ChrNumb", "Start","End","Ploidy","MajorMinor")

          }

          # if(is.null(MinorPloidy)) {
          #         print(paste0("StartPos_Minor:", StartPos_Major))
          #         print(paste0("EndPos_Minor:", EndPos_Major))
          #         print("MinorPloidy NULL break");
          #    next;
          # }

          if(MinorPloidy[1]>0) {
              ## MinorPloidy
              StartPos_Minor <- EachStartPosition; print(StartPos_Minor) # as.numeric(sapply(strsplit(EachSegment, split="_"), "[", 2)); print(StartPos)
              EndPos_Minor <-  CNVSeg_RmvByGiniMeanState_Subset_chrstrend_MinorPloidySubset$end[nrow(CNVSeg_RmvByGiniMeanState_Subset_chrstrend_MinorPloidySubset)];  print(EndPos_Minor)

              for(EachMinorPloidy in  MinorPloidy) {
                  h_minor <- EachMinorPloidy; print(h_minor)
                  kpSegments(karyoplot = karyoplot, chr = chr, x0 = StartPos_Minor,   ### This make a horizontal line. I need X0, X1, and h
                           x1 = EndPos_Minor, y0=h_minor, y1=h_minor, ymin=ymin, ymax=ymax,
                           r0 = r0, r1 = r1, data.panel = data.panel, clipping = clipping, col="orange", lwd=7)
              }

              ### collect MinorPloidy information.
              MinorPloidyDF <- data.frame(t(c(EachChr, StartPos_Minor,EndPos_Minor, h_minor, "Minor" )));
              rownames(MinorPloidyDF) <- paste0("Minor_",StartPos_Minor,"_",EndPos_Minor); colnames(MinorPloidyDF) <- c("ChrNumb", "Start","End","Ploidy","MajorMinor")
          }

          ### I was working on here on Tuesday Oct. 29th
          if(ncol(MajorPloidyDF)==5 & ncol(MinorPloidyDF)==5) {
                MajorMinorPloidy <- rbind(MajorPloidyDF, MinorPloidyDF)
                MajorMinorPloidy_All <- rbind(MajorMinorPloidy_All, MajorMinorPloidy)
          } else if(ncol(MajorPloidyDF)==5 & ncol(MinorPloidyDF)==0) {
                MajorMinorPloidy_All <- rbind(MajorMinorPloidy_All, MajorPloidyDF)
          }

          # if(EachMinorPloidy==0) {
          #       print(paste0("StartPos_Minor:", StartPos_Major))
          #       print(paste0("EndPos_Minor:", EndPos_Major))
          #       print("EachMinorPloidy 0 break"); break;
          # }

        } # end of the 2nd for loop

        #### Make horizontal line by the number of cells.

        if(EachChr=="chr1") { ## 1p : 1~123400000
              # CNVSeg_RmvByGiniMeanState_Subset_chrstrend  This is subset by only chr1.
              CNVSeg_RmvByGiniMeanState_Subset_chrstrend_1p <- dplyr::filter(CNVSeg_RmvByGiniMeanState_Subset_chrstrend, end < 123400000);
              dim(CNVSeg_RmvByGiniMeanState_Subset_chrstrend_1p) # 2303 7

              CNVSeg_RmvByGiniMeanState_Subset_chrstrend_1q <- dplyr::filter(CNVSeg_RmvByGiniMeanState_Subset_chrstrend, start > 150600001) # ,  end > 150600001);
              dim(CNVSeg_RmvByGiniMeanState_Subset_chrstrend_1q) # 2613 7

              ## Find the 2nd major ploidy in 1p or 1q, and make horizontal line.
              ## Chr1p
              PloidyTable_1p <- table(CNVSeg_RmvByGiniMeanState_Subset_chrstrend_1p$state); #
              print(PloidyTable_1p);
              # Main Ploidy
              MajorPloidy_ByChr1p <- NULL;
              if(sort(PloidyTable_1p, decreasing=TRUE)[1] >= 400) {
                    MajorPloidy_ByChr1p <- as.integer(names(PloidyTable_1p)[ PloidyTable_1p== sort(PloidyTable_1p, decreasing=TRUE)[1]  ]  );
                    ymin=0; ymax=5;r0=0; r1=1;
                    kpSegments(karyoplot=kp, chr=EachChr, x0=1,   ### This make a horizontal line. I need X0, X1, and h
                               x1=123400000, y0=MajorPloidy_ByChr1p, y1=MajorPloidy_ByChr1p, ymin=ymin, ymax=ymax,
                               r0=r0, r1=r1, data.panel=data.panel, clipping = clipping, col="green", lwd=7)
              }

              ## Minor Ploidy
              MinorPloidy_ByChr1p <- NULL;
              if(sort(PloidyTable_1p, decreasing=TRUE)[2] >= 400) {
                    MinorPloidy_ByChr1p <- as.integer(names(PloidyTable_1p)[ PloidyTable_1p== sort(PloidyTable_1p, decreasing=TRUE)[2]  ]  );
                    ymin=0; ymax=5; r0=0; r1=1;
                    kpSegments(karyoplot=kp, chr=EachChr, x0=1,   ### This make a horizontal line. I need X0, X1, and h
                               x1=123400000, y0=MinorPloidy_ByChr1p, y1=MinorPloidy_ByChr1p, ymin=ymin, ymax=ymax,
                               r0=r0, r1=r1, data.panel=data.panel, clipping = clipping, col="orange", lwd=7)
              }

              ## Chr1q
              PloidyTable_1q <- table(CNVSeg_RmvByGiniMeanState_Subset_chrstrend_1q$state); #
              print(PloidyTable_1q);
            } else {  # end of if-phrase to draw horizontal line by the number of cells.
                  PloidyTable_ByCellNumb <- table(CNVSeg_RmvByGiniMeanState_Subset_chrstrend$state); #
                  print(PloidyTable_ByCellNumb);

                  if(sort(PloidyTable_ByCellNumb, decreasing=TRUE)[1] >= 1200) {
                        MajorPloidy_ByCellNumb <- as.integer(names(PloidyTable_ByCellNumb)[ PloidyTable_ByCellNumb== sort(PloidyTable_ByCellNumb, decreasing=TRUE)[1]  ]  );
                        ymin=0; ymax=5;r0=0; r1=1;
                        # kpSegments(karyoplot=kp, chr=EachChr, x0=1,   ### This make a horizontal line. I need X0, X1, and h
                        #            x1=as.numeric(MajorMinorPloidy_All$End[nrow(MajorMinorPloidy_All)]), y0=MajorPloidy_ByCellNumb, y1=MajorPloidy_ByCellNumb, ymin=ymin, ymax=ymax,
                        #            r0=r0, r1=r1, data.panel=data.panel, clipping = clipping, col="green", lwd=7)
                  }
            }

    } # end of the first for loop
}



#' Title
#'
#' @param CNVSegmentFile  original CNV segment file
#' @param HierarchialClusterInfoFile cluster number of cells
#' @param HierarchialClusterInfoAfRmvClusterFile  cluster number of cells after removing outliers
#' @param HierarchialClusterInfoAfRmvClusterChr1File  cluster number of cells in chr 1
#' @param CNVStateData_AfterRmvHighGiniMeanFile   CNV state data after removing cells of high gini and high mean
#' @param CNVStateData_AfterRmvHighGiniMeanClusterFile  CNV state data after removing cells of high gini and high mean, and specific clusters
#' @param ClusterNumbToRmv  cluster number to remove, NULL
#' @param ClusterNumbToPick  cluster number to make ploidy plot
#' @param ClusterNumbToPick_InChr1  cluster number to make ploidy plots in chr1
#' @param MySampleID  manually typed in
#'
#' @returns NULL
#' @export
#'
#' @examples   Karyoplote_BySegment(CNVSegmentFile,HierarchialClusterInfoFile,HierarchialClusterInfoAfRmvClusterFile, HierarchialClusterInfoAfRmvClusterChr1File,CNVStateData_AfterRmvHighGiniMeanFile,CNVStateData_AfterRmvHighGiniMeanClusterFile, ClusterNumbToRmv,  ClusterNumbToPick,ClusterNumbToPick_InChr1, MySampleID)
Karyoplote_BySegment <- function(CNVSegmentFile="/Users/lees130/Library/CloudStorage/OneDrive-NYULangoneHealth/N07_DLPplus_scWGS/04_DLPplus_1qGain_8ptData/01-213_143839A/hmmcopy_segments.csv.gz",
                                 # HierarchialClusterInfoFile="/Users/lees130/Library/CloudStorage/OneDrive-NYULangoneHealth/N07_DLPplus_scWGS/09b_HMMCopyHierarchicalCluster_Rpackage/02_OutHMMcopyClusteringHeatmap_RmvBinCellByGiniMean/RowSplitHeatmapDraw_WithoutArtifact_AllChr_01_213_143839A.rds",
                                 HierarchialClusterInfoAfRmvClusterFile="/Users/lees130/Library/CloudStorage/OneDrive-NYULangoneHealth/N07_DLPplus_scWGS/09b_HMMCopyHierarchicalCluster_Rpackage/02_OutHMMcopyClusteringHeatmap_RmvBinCellByGiniMean/RowSplitHeatmapDraw_WithoutArtifactByCluster_AllChr_01_213_143839A.rds",
                                 HierarchialClusterInfoAfRmvClusterChr1File=NULL,
                                 CNVStateData_AfterRmvHighGiniMeanFile="/Users/lees130/Library/CloudStorage/OneDrive-NYULangoneHealth/N07_DLPplus_scWGS/09b_HMMCopyHierarchicalCluster_Rpackage/02_OutHMMcopyClusteringHeatmap_RmvBinCellByGiniMean/HMMcopy_StateIdeal_RmvCellByGiniMean_WithoutBinRegion_01_213_143839A.rds",
                                 CNVStateData_AfterRmvHighGiniMeanClusterFile=NULL,
                                 ClusterNumbToRmv=c(7,8,9),  ClusterNumbToPick=c(6,8),ClusterNumbToPick_InChr1=c(1,2),  MySampleID="01_213_143839A") {

    print(paste0("CNVSegmentFile: ", CNVSegmentFile))
    CNVSegment <- data.table::fread(CNVSegmentFile, header=TRUE, stringsAsFactors=FALSE); CNVSegment[1:2,];
    table(CNVSegment$state); table(CNVSegment$multiplier); CNVSegment[1:2,] ; length(unique(CNVSegment$cell_id)) #  0~11, # 1~6  # multiplier: median copy number value of segment # 1216
    #         chr   start       end state   median multiplier                cell_id
    #       <char>   <int>     <int> <int>    <num>      <int>                 <char>
    #   1:   chr1       1   7500000     5 4.564473          4 01-213-143839A-R03-C39
    # plot(CNVSegment$median, CNVSegment$state) # state and median have a good positive correlation.

    ## I need 'StartEndChrRead", "state", or "cell_id" columns.   ### Also, delete chrX and chrY
    # CNVSeg_StateIdealSelect <- dplyr::mutate(CNVSegment, StartEndChr=paste0(chr,"_", start,"_",end)) # %>%    # "_",seq(1:nrow(CNVSeg_StateIdeal
    CNVSeg_StateIdealSelect <- dplyr::select(CNVSegment, c(chr, start,end, state,median, cell_id)) %>% dplyr::filter(!chr %in% c("chrX","chrY"));
    CNVSeg_StateIdealSelect[1:2,]; dim(CNVSeg_StateIdealSelect) #  405339      5

    CNVSeg_StateIdealSelect_Delchr <- dplyr::mutate(CNVSeg_StateIdealSelect, chr=gsub("chr","", CNVSeg_StateIdealSelect$chr));
    CNVSeg_StateIdealSelect_Delchr[1:2,]; dim(CNVSeg_StateIdealSelect_Delchr) # 387836      5

    ## +++++++++  After sorting.
    # sort(CNVSeg_StateIdealSelect_AllCellID_Proc$Start[CNVSeg_StateIdealSelect_AllCellID_Proc$Chr=="chr01"])[1:5] # [1] 1500001 2000001 3000001 3500001 4000001 This is goal.
    CNVSeg_StateIdealSelect_Sort <- dplyr::arrange(CNVSeg_StateIdealSelect_Delchr, chr, as.numeric(start), as.numeric(end),cell_id);

    ## remove "NA" rows in 'state' data.
    table(is.na(CNVSeg_StateIdealSelect_Sort$state)) #All FALSE: 387836
    CNVSeg_StateIdealSelect_NonNA <- CNVSeg_StateIdealSelect_Sort[!is.na(CNVSeg_StateIdealSelect_Sort$state),]
    CNVSeg_StateIdealSelect_NonNA <- CNVSeg_StateIdealSelect_Sort[!is.na(CNVSeg_StateIdealSelect_Sort$median),]
    dim(CNVSeg_StateIdealSelect_NonNA); CNVSeg_StateIdealSelect_NonNA[1:2,]; table(is.na(CNVSeg_StateIdealSelect_NonNA$state)) # 387836       5

    CNVSeg_StateIdealSelect_Proc <- CNVSeg_StateIdealSelect_NonNA

    CNVSeg_StateIdealSelect_Proc$cell_id <- paste0("X", CNVSeg_StateIdealSelect_Proc$cell_id)
    CNVSeg_StateIdealSelect_Proc$cell_id <- gsub("-","_", CNVSeg_StateIdealSelect_Proc$cell_id)
    # CNVSeg_StateIdealSelect_Proc$median <- log2(as.numeric(CNVSeg_StateIdealSelect_Proc$median)) # This is not good.

    #CNVSeg_StateIdealSelect_ProcState <- dplyr::mutate(state=as.numeric(10^(-1*as.numeric(CNVSeg_StateIdealSelect_Sort$state))))

    CNVSeg_StateIdealSelect_Proc[1:2,]; dim(CNVSeg_StateIdealSelect_Proc) ; length(unique(CNVSeg_StateIdealSelect_Proc$cell_id))# 09_025, 618248      6
    #           chr start     end state   median                 cell_id
    #        <char> <int>   <int> <int>    <num>                  <char>
    #   1:      1     1     1500000  8 7.888887 X01_206_143839A_R06_C36

    ##########################################################################################################################
    ##### +++++++++++++++++ Compare Cell IDs between after removing by MeanState + Geni Coef  and By Cluster 6789
    if(is.null(CNVStateData_AfterRmvHighGiniMeanFile)) {
        CNVSeg_StateIdeal_AfterRmvArtifactCell_ColNameProc <- NULL
    } else {
        CNVSeg_StateIdeal_AfterRmvArtifactCell_ColNameProc <- readRDS(CNVStateData_AfterRmvHighGiniMeanFile)
    }
    dim(CNVSeg_StateIdeal_AfterRmvArtifactCell_ColNameProc); length(unique(rownames(CNVSeg_StateIdeal_AfterRmvArtifactCell_ColNameProc))) # 09_025, 497  5105
    #                         chr01   .1 .2 .3
    # X01_206_143839A_R19_C31     2 2  2  2  2

    # print(paste0("Number of cells after removing artifact by clusters: ", nrow(CNVSeg_StateIdeal_RmvCluster6789) ))

    if(is.null(CNVSeg_StateIdeal_AfterRmvArtifactCell_ColNameProc)) {
    } else {
        print(paste0("Number of cells after removing artifact by high Gini and Mean: ", nrow(CNVSeg_StateIdeal_AfterRmvArtifactCell_ColNameProc) )) # 656

        ## Remove cells of high Gini coef and High Mean State
        CNVSeg_StateIdealSelect_RmvByGiniMeanState <-dplyr::filter(CNVSeg_StateIdealSelect_Proc,
                                                                   (CNVSeg_StateIdealSelect_Proc$cell_id %in% rownames(CNVSeg_StateIdeal_AfterRmvArtifactCell_ColNameProc)))
        dim(CNVSeg_StateIdealSelect_RmvByGiniMeanState); length(unique(CNVSeg_StateIdealSelect_RmvByGiniMeanState$cell_id)) # 09_025, 30798  6 # 09_025, 497

    }
    NumbOfCell <- length(unique(CNVSeg_StateIdealSelect_RmvByGiniMeanState$cell_id)); print(paste0("Numb of cells: ",NumbOfCell)) # 497
    #     chromosome   start     end probes  segmean     sample       # I don't need 'probes' column
    # 1          1       1  880000     13 2.021323 H_OM-15920

    #### +++++++++++  Convert dta.frame to GRanges.  ++++++++++++++ ###
    CNVSeg_StateIdealSelect_RmvByGiniMeanState <- dplyr::mutate(CNVSeg_StateIdealSelect_RmvByGiniMeanState, chr=paste0("chr", chr))
    CNVSeg_StateIdealSelect_RmvByGiniMeanState[1:2,]
    #       chr start     end state   median                 cell_id
    #     <char> <int>   <int> <int>    <num>                  <char>
    # 1: chr1     1 1500000     4 3.999999 X01_213_143839A_R08_C55
    # 2: chr1     1 1500000     2 2.000005 X01_213_143839A_R16_C40

    CNVSeg_StateIdealSelect_GRanges <- makeGRangesFromDataFrame(CNVSeg_StateIdealSelect_RmvByGiniMeanState); #CNVSeg_StateIdealSelect_GRanges[100000:100040,]
    CNVSeg_StateIdealSelect_GRanges[1:5,]
    #       chromosome start     end state  segmean                  sample
    #           <int> <int>   <int> <int>    <num>                  <char>
    # 1:          1     1 1500000     4 3.999999 X01_213_143839A_R08_C55
    # 2:          1     1 1500000     2 2.000005 X01_213_143839A_R16_C40

    ### Karyoplote rr
    CNVSeg_StateIdealSelect_rr <- CNVSeg_StateIdealSelect_GRanges
    CNVSeg_StateIdealSelect_rr$cn <- 0

    ### Karyoplote Tiles
    CNVSeg_StateIdealSelect_Tiles <- CNVSeg_StateIdealSelect_GRanges
    CNVSeg_StateIdealSelect_Tiles$lrr <- CNVSeg_StateIdealSelect_RmvByGiniMeanState$median; CNVSeg_StateIdealSelect_Tiles[1:5,]

    # CNVSeg_Chr1 <- CNVSeg_StateIdealSelect_RmvByGiniMeanState %>% data.frame %>% dplyr::filter(chr=="chr1"); dim(CNVSeg_Chr1) # 3244 6
    CNVSeg_Chr1 <- dplyr::filter(data.frame(CNVSeg_StateIdealSelect_RmvByGiniMeanState), chr=="chr1"); dim(CNVSeg_Chr1) # 3244 6
    table(CNVSeg_Chr1$state)
    #### +++++++++  KaryoploteR   ++++++++++ ######
    #### +++++++++++++++++ #### +++++++++++++++++ #### +++++++++++++++++ #### +++++++++++++++++ #### +++++++++++++++++ #### +++++++++++++++++
    # pdf(paste0("OutKaryoplote_AllCell_BeforeRmvHighGiniMean_", MySampleID, ".pdf"), height=6,width=12)
    #       kp <- plotKaryotype("hg38", plot.type = 4, main="Random Data, 2Mb resolution", cex=2) #### Kayrotype plot background.  # plot.type is upto 7. "4" is the best
    #       plotLRR(kp, CNVSeg_StateIdealSelect_Tiles, ymin=0, ymax=16, add.axis = FALSE, labels=NA)
    #       plotCopyNumberCallsAsLines(kp, cn.calls = CNVSeg_StateIdealSelect_rr, ymin=0, ymax=12, labels="ploidy")
    # dev.off()

    ##########################################################################################################################
    ## remove Cells by Hierarchical cluster numbers
    #### +++++++++++++++++ #### +++++++++++++++++ #### +++++++++++++++++ #### +++++++++++++++++ #### +++++++++++++++++ #### +++++++++++++++++
    # RowSplitHeatmapDraw_NoRmvCell <- readRDS(HierarchialClusterInfoFile)
    # RowIndexExtract <-ComplexHeatmap::row_order(RowSplitHeatmapDraw_NoRmvCell); sum(length(unlist(RowIndexExtract)))  # 637

    ### Find Cell IDs after remove artifacts by cluster number of cells.
    # CNVSeg_StateIdeal_RmvCluster6789 <-  RemoveReorderCell_ByClusterNumb(HMMcopy_StateIdeal_AfterRmvArtifactCell_ColNameProc=CNVSeg_StateIdeal_AfterRmvArtifactCell_ColNameProc,
    #                                                                      HierarchialClusterInfoFile, ClusterNumbToRmv, ClusterNumbToPick,ClusterReorder=NULL,MySampleID,ForPloidyPlot=TRUE)
    # dim(CNVSeg_StateIdeal_RmvCluster6789) # 553 5259.   ##### Problem is here.
    CNVSeg_StateIdeal_RmvCluster6789 <- CNVSeg_StateIdeal_AfterRmvArtifactCell_ColNameProc; dim(CNVSeg_StateIdeal_RmvCluster6789)


    ### Filter in cell id that were obtained after removing bad cluster CellIDs
    CNVSeg_StateIdealSelect_RmvByGiniMeanState_RmvByCluster <- dplyr::filter(CNVSeg_StateIdealSelect_RmvByGiniMeanState,
                                                                             CNVSeg_StateIdealSelect_RmvByGiniMeanState$cell_id %in% rownames(CNVSeg_StateIdeal_RmvCluster6789))
    dim(CNVSeg_StateIdealSelect_RmvByGiniMeanState_RmvByCluster) # 50589       5
    # Checking one more time.
    table(names(table(CNVSeg_StateIdealSelect_RmvByGiniMeanState_RmvByCluster$cell_id)) %in% rownames(CNVSeg_StateIdeal_RmvCluster6789)) ## Only TRUE 427

    #### +++++++++++  Convert data.frame to GRanges.  ++++++++++++++ ###
    if(!unique(grepl("chr", CNVSeg_StateIdealSelect_RmvByGiniMeanState_RmvByCluster$chr))) {
      CNVSeg_StateIdealSelect_RmvByGiniMeanState_RmvByCluster <- dplyr::mutate(CNVSeg_StateIdealSelect_RmvByGiniMeanState_RmvByCluster, chr=paste0("chr", chr))
    }
    CNVSeg_StateIdealSelect_RmvByGiniMeanState_RmvByCluster[1:2,]; length(unique(CNVSeg_StateIdealSelect_RmvByGiniMeanState_RmvByCluster$cell_id))
    dim(CNVSeg_StateIdealSelect_RmvByGiniMeanState_RmvByCluster)#  91446     6 # [1] 45873     6

    #       chr start     end state   median                 cell_id
    #     <char> <int>   <int> <int>    <num>                  <char>
    # 1: chr1     1 1500000     4 3.999999 X01_213_143839A_R08_C55
    # 2: chr1     1 1500000     2 2.000005 X01_213_143839A_R16_C40

    CNVSeg_StateIdealSelect_RmvByCluster_GRanges <- makeGRangesFromDataFrame(CNVSeg_StateIdealSelect_RmvByGiniMeanState_RmvByCluster);
    CNVSeg_StateIdealSelect_RmvByCluster_GRanges[1:5,]
    #         GRanges object with 5 ranges and 0 metadata columns:
    #       seqnames    ranges strand
    #         <Rle> <IRanges>  <Rle>
    #   [1]     chr1 1-1500000      *

    # ### Karyoplote rr
    CNVSeg_StateIdealSelect_RmvByCluster_rr <- CNVSeg_StateIdealSelect_RmvByCluster_GRanges
    CNVSeg_StateIdealSelect_RmvByCluster_rr$cn <- 0

    ### Karyoplote Tiles
    CNVSeg_StateIdealSelect_RmvByCluster_Tiles <- CNVSeg_StateIdealSelect_RmvByCluster_GRanges
    CNVSeg_StateIdealSelect_RmvByCluster_Tiles$lrr <- CNVSeg_StateIdealSelect_RmvByGiniMeanState_RmvByCluster$median; CNVSeg_StateIdealSelect_RmvByCluster_Tiles[1:5,]

    #### +++++++++  KaryoploteR for all cells filtered by High Gini Coef and High mean and by cluster #  ++++++++++ ######
    #### +++++++++++++++++ #### +++++++++++++++++ #### +++++++++++++++++ #### +++++++++++++++++ #### +++++++++++++++++ #### +++++++++++++++++
    ClusterNumbConcat <- paste(ClusterNumbToRmv, collapse=""); print(ClusterNumbConcat)
    pdf(paste0("OutKaryoplote_AllCell_BeforeRmvHighGiniMean_BySegMedian_ByCluster",ClusterNumbConcat,"_", MySampleID, "_PloidyLine",NumbOfCell,"c.pdf"), height=6,width=12)
        kp <- karyoploteR::plotKaryotype("hg38", plot.type = 4, main="", cex=2) #### Kayrotype plot background.  # plot.type is upto 7. "4" is the best
        plotLRR(kp, snps=CNVSeg_StateIdealSelect_RmvByCluster_Tiles, ymin=0, ymax=5, add.axis=FALSE, labels=NA)
        plotCopyNumberCallsAsLines(kp, cn.calls = CNVSeg_StateIdealSelect_RmvByCluster_rr, ymin=0, ymax=5, labels="")  ## y-axis lines.
        ### Go to line 629 to draw horizontal line.
        CNVStatePloidyLine(CNVSeg_StateIdealSelect_RmvByGiniMeanState_RmvByCluster=CNVSeg_StateIdealSelect_RmvByGiniMeanState_RmvByCluster, TotalSumThreshold=30)
    dev.off()

    #### +++++++++++++++++ #### +++++++++++++++++ #### +++++++++++++++++ #### +++++++++++++++++ #### +++++++++++++++++  #### +++++++++++++++++  #### +++++++++++++++++  #### +++++++++++++++++
    ##########################################################################################################################
    # CNVSeg_StateIdealSelect_OneCell <- dplyr::filter(CNVSeg_StateIdealSelect_Proc, grepl("01-213-143839A-R16-C62", cell_id))
    # dim(CNVSeg_StateIdealSelect_OneCell) # 354 4
    # manhattan(x=CNVSeg_StateIdealSelect_OneCell, chr="chr", bp="start", snp="cell_id", p="state", cex=0.3 )

    ##########################################################################################################################
    #### +++++++++++++++++ ### Pick up just 2 cluster cells and make Karyotype plot.  +++++++++++++ ###########
    ##########################################################################################################################
    ### Find Cell IDs of my interestcluster number of cells.
    CNVSeg_StateIdeal_PickCluster79 <-  RemoveReorderCell_ByClusterNumb(HMMcopy_StateIdeal_AfterRmvArtifactCell_ColNameProc= CNVSeg_StateIdeal_RmvCluster6789,
                                                                 HierarchialClusterInfoFile=HierarchialClusterInfoAfRmvClusterFile, ClusterNumbToRmv=NULL,  ClusterNumbToPick=NULL,ClusterReorder=NULL,MySampleID,ForPloidyPlot=TRUE)
    dim(CNVSeg_StateIdeal_PickCluster79) # 291  5359

    for (EachClusterNumber in ClusterNumbToPick) {
        # EachClusterNumber <- ClusterNumbToPick[2]; print(paste0("Number of cluster to pick: ", EachClusterNumber))
        CNVSeg_StateIdeal_PickCluster7Only <- CNVSeg_StateIdeal_PickCluster79[grepl(paste0("Cluster",EachClusterNumber), rownames(CNVSeg_StateIdeal_PickCluster79)), ]; dim(CNVSeg_StateIdeal_PickCluster7Only) # 44 5259

        # ### Filter in cell id that were picked up with interesting clusters - Cluster 9, 7
        # CNVSeg_StateIdealSelect_RmvByGiniMeanState_PickByCluster <- dplyr::filter(CNVSeg_StateIdealSelect_RmvByGiniMeanState,
        #                                                         CNVSeg_StateIdealSelect_RmvByGiniMeanState$cell_id %in% gsub("_Cluster(.*)", "", rownames(CNVSeg_StateIdeal_PickCluster79)))
        # dim(CNVSeg_StateIdealSelect_RmvByGiniMeanState_PickByCluster) # 35981       5

        ## ++++ #########################  Filter in cell id that were picked up with Cluster7       ## ++++ #########################
        CNVSeg_StateIdealSelect_RmvByGiniMeanState_PickByCluster7 <- dplyr::filter(CNVSeg_StateIdealSelect_RmvByGiniMeanState,
                                                                                   CNVSeg_StateIdealSelect_RmvByGiniMeanState$cell_id %in% gsub("_Cluster(.*)", "", rownames(CNVSeg_StateIdeal_PickCluster7Only)))
        #### +++++++++++  Convert dta.frame to GRanges.  ++++++++++++++ ###
        # CNVSeg_StateIdealSelect_RmvByGiniMeanState_PickByCluster7 <- dplyr::mutate(CNVSeg_StateIdealSelect_RmvByGiniMeanState_PickByCluster7, chr=paste0("chr", chr))
        CNVSeg_StateIdealSelect_RmvByGiniMeanState_PickByCluster7[1:2,]
        #       chr start     end state   median                 cell_id
        #     <char> <int>   <int> <int>    <num>                  <char>
        # 1:  chr1     1 2500000     1 1.150863 X09_025_143855A_R15_C51

        CNVSeg_StateIdealSelect_PickByCluster7_GRanges <- makeGRangesFromDataFrame(CNVSeg_StateIdealSelect_RmvByGiniMeanState_PickByCluster7);
        CNVSeg_StateIdealSelect_PickByCluster7_GRanges[1:5,]
        #       chromosome start     end state  segmean                  sample
        #           <int> <int>   <int> <int>    <num>                  <char>
        # 1:          1     1 1500000     4 3.999999 X01_213_143839A_R08_C55

        ### Karyoplote rr
        CNVSeg_StateIdealSelect_PickByCluster7_rr <- CNVSeg_StateIdealSelect_PickByCluster7_GRanges
        CNVSeg_StateIdealSelect_PickByCluster7_rr$cn <- 0

        ### Karyoplote Tiles
        CNVSeg_StateIdealSelect_PickByCluster7_Tiles <- CNVSeg_StateIdealSelect_PickByCluster7_GRanges
        CNVSeg_StateIdealSelect_PickByCluster7_Tiles$lrr <- CNVSeg_StateIdealSelect_RmvByGiniMeanState_PickByCluster7$median; CNVSeg_StateIdealSelect_PickByCluster7_Tiles[1:5,]

        #### +++++++++  KaryoploteR   ++++++++++ ######
        #### +++++++++++++++++ #### +++++++++++++++++ #### +++++++++++++++++ #### +++++++++++++++++ #### +++++++++++++++++ #### +++++++++++++++++
        # pdf(paste0("OutKaryoplote_AllCell_BeforeRmvHighGiniMean_PickByCluster",EachClusterNumber, "_", MySampleID, "_PloidyLine.pdf"), height=6,width=12)
        #       kp <- plotKaryotype("hg38", plot.type = 4, main="",  cex=2) #### Kayrotype plot background.  # plot.type is upto 7. "4" is the best
        #       plotLRR(kp, CNVSeg_StateIdealSelect_PickByCluster7_Tiles, ymin=0, ymax=5, add.axis = FALSE, labels=NA)
        #       plotCopyNumberCallsAsLines(kp, cn.calls = CNVSeg_StateIdealSelect_PickByCluster7_rr, ymin=0, ymax=5, labels="")
        #       CNVStatePloidyLine(CNVSeg_StateIdealSelect_RmvByGiniMeanState_RmvByCluster=CNVSeg_StateIdealSelect_RmvByGiniMeanState_PickByCluster7,  TotalSumThreshold=30)
        # dev.off()
    } # end of for loop

    #### +++++++++++++++++ #### +++++++++++++++++ #### +++++++++++++++++ #### +++++++++++++++++ #### +++++++++++++++++  #### +++++++++++++++++  #### +++++++++++++++++  #### +++++++++++++++++
    #### +++++++++++++++++ #### +++++++++++++++++ #### +++++++++++++++++ #### +++++++++++++++++ #### +++++++++++++++++  #### +++++++++++++++++  #### +++++++++++++++++  #### +++++++++++++++++
    #### Just chromosome 1 segments and make Karyoplot.  ###############
    if(!is.null(HierarchialClusterInfoAfRmvClusterChr1File)) {
        CNVSeg_StateIdeal_RmvCluster6789_ReorderBy <- readRDS(CNVStateData_AfterRmvHighGiniMeanClusterFile); dim(CNVSeg_StateIdeal_RmvCluster6789_ReorderBy) # 587 426

        CNVSeg_StateIdeal_PickCluster79_Chr1 <-RemoveCell_ByClusterNumb(HMMcopy_StateIdeal_AfterRmvArtifactCell_ColNameProc=CNVSeg_StateIdeal_RmvCluster6789_ReorderBy,
                                                                        HierarchialClusterInfoFile=HierarchialClusterInfoAfRmvClusterChr1File, ClusterNumbToRmv=NULL, ClusterNumbToPick=ClusterNumbToPick_InChr1)
        dim(CNVSeg_StateIdeal_PickCluster79_Chr1); CNVSeg_StateIdeal_PickCluster79_Chr1[1:2,1:5] # 587 426

        for (EachClusterNumber in ClusterNumbToPick_InChr1) {
            # EachClusterNumber <- ClusterNumbToPick_InChr1[3]; print(paste0("Number of cluster to pick: ", EachClusterNumber))
            CNVSeg_StateIdeal_PickCluster7Only_Chr1 <- CNVSeg_StateIdeal_PickCluster79_Chr1[grepl(paste0("Cluster",EachClusterNumber), rownames(CNVSeg_StateIdeal_PickCluster79_Chr1)), ];
            dim(CNVSeg_StateIdeal_PickCluster7Only_Chr1) # 73 426

            ## ++++ #########################  Filter in cell id that were picked up with Cluster7       ## ++++ #########################
            CNVSeg_StateIdealSelect_RmvByGiniMeanState_PickByCluster7_Chr1 <- dplyr::filter(CNVSeg_StateIdealSelect_RmvByGiniMeanState,
                                                                                            CNVSeg_StateIdealSelect_RmvByGiniMeanState$cell_id %in% gsub("_Cluster(.*)", "", rownames(CNVSeg_StateIdeal_PickCluster7Only_Chr1)),
                                                                                            CNVSeg_StateIdealSelect_RmvByGiniMeanState$chr == "chr1")
            #### +++++++++++  Convert dta.frame to GRanges.  ++++++++++++++ ###
            CNVSeg_StateIdealSelect_RmvByGiniMeanState_PickByCluster7_Chr1[1:2,]
            #       chr start     end state   median                 cell_id
            #     <char> <int>   <int> <int>    <num>                  <char>
            # 1: chr1     1 1500000     4 3.999999 X01_213_143839A_R08_C55

            ## Test
            # CNVSeg_StateIdealSelect_RmvByGiniMeanState_PickByCluster7_Chr1 <- CNVSeg_StateIdealSelect_RmvByGiniMeanState_PickByCluster7_Chr1[1:20,]
            # CNVSeg_StateIdealSelect_RmvByGiniMeanState_PickByCluster7_Chr1$start <- 2000000

            CNVSeg_StateIdealSelect_PickByCluster7_Chr1_GRanges <- makeGRangesFromDataFrame(CNVSeg_StateIdealSelect_RmvByGiniMeanState_PickByCluster7_Chr1);
            CNVSeg_StateIdealSelect_PickByCluster7_Chr1_GRanges[1:5,]
            #       chromosome start     end state  segmean                  sample
            #           <int> <int>   <int> <int>    <num>                  <char>
            # 1:          1     1 1500000     4 3.999999 X01_213_143839A_R08_C55

            ### Karyoplote rr
            CNVSeg_StateIdealSelect_PickByCluster7_Chr1_rr <- CNVSeg_StateIdealSelect_PickByCluster7_Chr1_GRanges
            CNVSeg_StateIdealSelect_PickByCluster7_Chr1_rr$cn <- 0

            ### Karyoplote Tiles
            CNVSeg_StateIdealSelect_PickByCluster7_Chr1_Tiles <- CNVSeg_StateIdealSelect_PickByCluster7_Chr1_GRanges
            CNVSeg_StateIdealSelect_PickByCluster7_Chr1_Tiles$lrr <- CNVSeg_StateIdealSelect_RmvByGiniMeanState_PickByCluster7_Chr1$median;
            CNVSeg_StateIdealSelect_PickByCluster7_Chr1_Tiles[1:5,]

            #### +++++++++  KaryoploteR   ++++++++++ ######
            #### +++++++++++++++++ #### +++++++++++++++++ #### +++++++++++++++++ #### +++++++++++++++++ #### +++++++++++++++++ #### +++++++++++++++++
            # pdf(paste0("OutKaryoplote_AllCell_BeforeRmvHighGiniMean_OnlyChr1_PickByCluster",EachClusterNumber, "_", MySampleID, "_PloidyLine",NumbOfCell,"c.pdf"), height=6,width=12)
            #     kp <- plotKaryotype("hg38", plot.type = 4, chromosomes="chr1",  main="", cex=2) #### Kayrotype plot background.  # plot.type is upto 7. "4" is the best
            #     ∂plotLRR(kp, CNVSeg_StateIdealSelect_PickByCluster7_Chr1_Tiles, ymin=0, ymax=5, add.axis = FALSE, labels=NA)
            #     plotCopyNumberCallsAsLines(kp, cn.calls = CNVSeg_StateIdealSelect_PickByCluster7_Chr1_rr, ymin=0, ymax=5, labels="")
            #     CNVStatePloidyLine_Chr1(CNVSeg_StateIdealSelect_RmvByGiniMeanState_RmvByCluster=CNVSeg_StateIdealSelect_RmvByGiniMeanState_PickByCluster7_Chr1,  TotalSumThreshold=30)
            # dev.off()
        }
    }

}


########## ########## ########## ########## ########## ########## ########## ########## ########## ##########
# PhylogeneticTree function.
########## ########## ########## ########## ########## ########## ########## ########## ########## ##########

#' Title PhylogeneticTree
#'
#' @param df data.frame of tree structure
#' @param num_cells_in data.frame of cell numbers
#' @param MySampleID my sample ID
#'
#' @returns
#' @export
#'
#' @examples PhylogeneticTree(df=TreeData, num_cells_in=NumbCellData, MySampleID=MySampleID)
PhylogeneticTree <- function(df=TreeData, num_cells_in=NumbCellData, MySampleID=MySampleID)  {
      num_cells <- num_cells_in[which(num_cells_in$sample==MySampleID),][,!names(num_cells_in) %in% c("sample")]

      colors = c("Normal" = "black", "Clone1" = "#7b0033", "Clone2"="#008EA0FF", "Clone3"="#8A4198FF", "Clone4"="#CFC44AFF", "Clone5"="#1A5354FF", "Clone6"="#CF4A5BFF", "Clone7"="#B87333",
                 "Clone8"="#660099FF","Clone9"="#003399FF", "Clone10"="#D60047FF",  "Clone11"="#9C964A", "Clone12"="#00D68FFF","Clone13"= "#0A47FFFF", "Clone14"="#5050FF99")

      sample <- MySampleID

      real_Clones <- grepl('^CNV', num_cells$name) | grepl('^N', num_cells$name) | grepl('^spCNV', num_cells$name) | grepl('^Clone', num_cells$name)
      real_Clones_names = num_cells[which(real_Clones),]$name
      #set size of fake Clones to mean of the rest (if any are present)
      if(length(real_Clones_names) != length(num_cells$name)) {
        num_cells[which(!real_Clones),]$counts <- mean(num_cells[which(real_Clones),]$counts) # fake Clone cell number.
        fake_Clones_names = num_cells[which(!real_Clones),]$name
      }
      if(length(real_Clones_names) == length(num_cells$name)){message("no unobserved Clones detected")}

      df_from_to <- df[,c("from", "to")]
      ig_x <- graph_from_data_frame(df_from_to, directed=TRUE, vertices=num_cells$name) #https://igraph.org/r/doc/graph_from_data_frame.html
      ig_x_phylo <- as.phylo(ig_x) #sending to ggtree uses as.phylo() method
      ig_x_phylo$edge.length= 1 #df$score- hold off for now/ for simplified plots, make edge lengths constant
      #ig_x_phylo$edge.length=rescale(ig_x_phylo$edge.length, to=c(1,5)) #should the score also be normalized incase they're very large/ wide range?

      #add line type to num_cells since already linked/ in order so can tell where x is
      #https://yulab-smu.top/treedata-book/faq.html#change-colors-or-line-types-of-arbitrarily-selected-branches
      lty <- rep(1, dim(num_cells)[1]) #make all 1 = normal
      #take from fake names to account for more than one fake Clone
      lty[which(num_cells$name %in% unique(df[which(df$from %in% fake_Clones_names),]$to))] <- 3 #child of fakes
      lty[which(num_cells$name %in% unique(df[which(df$from %in% fake_Clones_names),]$from))] <- 3 #fakes
      num_cells$lty <- lty

      for(i in 1:nrow(num_cells)) {
        if(num_cells$name[i] %in% fake_Clones_names){num_cells$cols[i] = "gray95"}
        else(num_cells$cols[i] = colors[num_cells$name[i]])
      }

      man_cols <- num_cells %>% pull(cols, name) #named vector for colors

      parent_child <- as.data.frame(cbind(ig_x_phylo$edge, df_from_to))
      colnames(parent_child) <- c("parent", "child", "from", "to")
      parent_child

      gg <- ggtree(ig_x_phylo, layout="slanted", size=3) %<+% num_cells
      gg_rot <- gg #switch x and y to rotate the plot
      gg_rot$data$x <- gg$data$y
      gg_rot$data$y <- -1 * gg$data$x #negate so the root is at the top

      PtNumb <- paste0("Pt",LoopNumb)

      message(glue("plotting {MySampleID}"))
      gg1 <- gg_rot + aes(linetype=I(lty)) +
        geom_point(aes(size=(counts*1.2))) + #slightly larger than points, acts as a border
        geom_point(aes(size=counts, color = factor(label))) + guides(color="none") + #just show legend for size nodes are labeled (don't need color legend)
        scale_color_manual(values = man_cols, na.value = "black") + #color by labels added to num_cells so they're in order
        geom_text(aes(label=label), position = position_nudge(0)) +
        ggtitle(glue("{MySampleID} Clonal Evolution")) +
        scale_size_continuous(range=c(10,30), name="Number of Cells") #use this instead of scaling the counts data, range sets smallest/ largest sizes (default they were very small)

      OutFileName <- glue("OutPhylogeneticTree_{MySampleID}.pdf"); print(OutFileName)
      # ggsave(plot=gg1, file=glue(paste0("{out_dir}/",PtNumb,"_{s}_inferCNVphylo.pdf")), width=10, height=8, units="in")
      ggsave(plot=gg1, file=OutFileName, width=10, height=8, units="in")

      print(paste0("this s is done: ", MySampleID))
      return("PhylogeneticTree is done")
}




