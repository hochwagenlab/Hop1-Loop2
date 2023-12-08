
library(hwglabr2)
library(GenomicRanges)
library(IRanges)
library(regioneR)
library(hwglabr2)
library(EnrichedHeatmap)
library(ggplot2)
library(rtracklayer)
library(BSgenome)
library(data.table)
library(stringr)
library(readr)
library(patchwork)
library('Repitools')


## Set Ggplot theme 
ggplot2_theme <- theme_classic() +
  theme(plot.title=element_text(hjust=0.5, size=10),
        plot.subtitle=element_text(hjust=0.5, size=10),
        axis.text=element_text(colour='black', type="Arial"),
        axis.ticks=element_line(colour='black'))

theme_set(ggplot2_theme)


setwd('/Users/darmokandjalad/Documents/Hop1PHD-Loop2/RAnalysis')

WtFile <- paste0('/Users/darmokandjalad/Documents/Hop1PHD-Loop2/RAnalysis/Hop1-wildtype-334-340-32-177-Reps-SK1Yue-PM_B3W3_MACS2_FE.bdg.gz')
Wtsig <- import_bedGraph(WtFile)
Wtsig <- sort(Wtsig)

loopFile <- paste0('/Users/darmokandjalad/Documents/Hop1PHD-Loop2/RAnalysis/11644-NHantiHop1-20220331-20220125-20221004-Reps-SK1Yue-PM-PE_B4_W3_MACS2_FE.bdg.gz')
loopSig <- import_bedGraph(loopFile)
loopSig <- sort(loopSig)

loopPchFile <- paste0('/Users/darmokandjalad/Documents/Hop1PHD-Loop2/RAnalysis/11757-hop1NH-20220125-20220721-20221004-Reps-SK1Yue-PM-PE_B4_W3_MACS2_FE.bdg.gz')
lpPchSig <- import_bedGraph(loopPchFile)
lpPchSig <- sort(lpPchSig)

pchFile <- paste0('/Users/darmokandjalad/Documents/Hop1PHD-Loop2/RAnalysis/11758-hop1NH-20220125-20220721-20221004-Reps-SK1Yue-PM-PE_B4_W3_MACS2_FE.bdg.gz')
pchSig <- import_bedGraph(pchFile)
pchSig <- sort(pchSig)

Hop1rec8 <- paste0('rec8D-Hop1-61-91-Reps-SK1Yue-PM_B3W3_MACS2_FE.bdg.gz')
rec8Sig <- import_bedGraph(Hop1rec8)
rec8Sig <- sort(rec8Sig)

rec8loop <- import.bedGraph('11689-antiHop1-20220125_20220331-Reps-SK1Yue-PM-PE_B4_W3_MACS2_FE.bdg.gz')
rec8loop <- sort(rec8loop)







#Spike-in normalization Function
############################################################################################
#' Compute spike-in normalization factor from total read counts
#'
#' Computes spike-in normalization factor between two spiked-in samples using
#' total counts of aligned reads. Inputs paths to text files containing counts
#' of aligned reads per chromosome of a hybrid SK1:S288C genome.
#' @param ref_chip_counts Either a single or a list of paths to reference ChIP
#' samples' read counts file. No default.
#' @param ref_input_counts Either a single or a list of paths to reference input
#' samples' read counts file. No default.
#' @param test_chip_counts Either a single or a list of paths to test ChIP
#' samples' read counts file. No default.
#' @param test_input_counts Either a single or a list of paths to test input
#' samples' read counts file. No default.
#' @param return_counts Logical indicating whether to return the computed read
#' counts instead of the normalization factor. Defaults to \code{FALSE}.
#' @return Numeric normalization factor.
#' @examples
#' \dontrun{
#' spikein_normalization_factor_from_counts(
#'      ref_chip_counts='Counts_AH119_chip.txt',
#'      ref_input_counts='Counts_AH119_input.txt',
#'      test_chip_counts='Counts_AH8104_chip.txt',
#'      test_input_counts='Counts_AH8104_input.txt')
#'
#' spikein_normalization_factor_from_counts(
#'     ref_chip_counts=list('Counts_AH119_chip_1.txt',
#'                          'Counts_AH119_chip_2.txt',
#'                          'Counts_AH119_chip_3.txt'),
#'     ref_input_counts=list('Counts_AH119_inp_1.txt',
#'                           'Counts_AH119_inp_2.txt',
#'                           'Counts_AH119_inp_3.txt'),
#'     test_chip_counts='Counts_AH8104_chip.txt',
#'     test_input_counts='Counts_AH8104_input.txt')
#' }
#' @export
spikein_normalization_factor_from_counts <- function(
  ref_chip_counts, ref_input_counts, test_chip_counts, test_input_counts,
  return_counts=FALSE) {
  
  # Put paths in list
  files <- list(ref_chip=ref_chip_counts, ref_input=ref_input_counts,
                test_chip=test_chip_counts, test_input=test_input_counts)
  
  # Convert each element into list, if not one already
  for (i in seq_along(files)) {
    if (!is.list(files[[i]])) files[[i]] <- list(files[[i]])
  }
  
  # Print files to read to console
  message('>>> Read alignment count files:')
  for (i in seq_along(files)) {
    for (file in files[[i]]) {
      message('   ', basename(file))
    }
  }    
  
  message()
  # Read files into tibble in list
  tables <- list()
  for (i in seq_along(files)) {
    tables[[i]] <- sapply(files[[i]], FUN=read_tsv, col_names=F,
                          simplify=FALSE, USE.NAMES=TRUE)
  }
  names(tables) <- names(files)
  
  message()
  # Get read counts per chromosome
  message('>>> Count reads per genome:')
  counts <- list()
  for (i in seq_along(tables)) {
    counts[[i]] <- sapply(tables[[i]], FUN=sum_per_genome,
                          simplify=FALSE, USE.NAMES=TRUE)
  }
  names(counts) <- names(tables)
  
  # Add-up counts for replicates (results in nested lists)
  for (i in seq_along(counts)) {
    if (length(counts[[i]]) > 1) {
      total <- counts[[i]][[1]]
      for (j in 2:length(counts[[i]])) {
        total <- total + counts[[i]][[j]]
      }
      counts[[i]] <- total
    } else counts[[i]] <- unlist(counts[[i]])
  }
  
  if (return_counts) {
    message('---')
    message('Done!')
    return(counts)
  }
  
  # Compute normalization factor
  result <- normalization_factor(ctrl_input=counts$ref_input,
                                 ctrl_chip=counts$ref_chip,
                                 test_input=counts$test_input,
                                 test_chip=counts$test_chip)
  
  message('---')
  message('Done!')
  
  return(result)
}

# Helper functions
sum_per_genome <- function(df) {
  # Compute sum of reads aligned to each genome
  S288C <- sum(
    df[apply(df, 1, function(x) str_detect(x[1],'_S288C')), 2])
  SK1 <- sum(
    df[apply(df, 1, function(x) str_detect(x[1], '_SK1')), 2])
  
  # Print result to console
  message('  S288C: ', formatC(S288C, big.mark=",",
                               drop0trailing=TRUE, format="f"))
  message('  SK1: ', formatC(SK1, big.mark=",",
                             drop0trailing=TRUE, format="f"))
  message('      ', round(S288C * 100 / (SK1 + S288C), 1), '% spike-in reads')
  
  # Return result as named vector
  c('S288C'=S288C, 'SK1'=SK1)
}


normalization_factor <- function(ctrl_input, ctrl_chip,
                                 test_input, test_chip) {
  # Compute Q values
  Q_ctrl_input <- ctrl_input['S288C'] / ctrl_input['SK1']
  Q_ctrl_chip <- ctrl_chip['S288C'] / ctrl_chip['SK1']
  
  Q_test_input <- test_input['S288C'] / test_input['SK1']
  Q_test_chip <- test_chip['S288C'] / test_chip['SK1']
  
  # Compute normalization factors
  a_ctrl <- Q_ctrl_input / Q_ctrl_chip
  a_test <- Q_test_input / Q_test_chip
  
  # Return reference strain-centric normalization factor
  a_test/ a_ctrl
}
############################

#put in all test and input counts to determine the amount of protein pulled down
############################################################################################




setwd('/Users/darmokandjalad/Documents/Hop1PHD-Loop2/RAnalysis/SpikeValue')
###########################

#Normalize to Spike in for hop1-loop2, hop1-loop2 pch2, pch2 * save these as bedgraphs
########################################################################

#Loop2_Norm
read_counts_loop2_1 <- data.frame(
  Condition=c('WT', 'loop2'),
  NF=c(1,
       spikein_normalization_factor_from_counts(
         ref_chip_counts='stats_H2CMNAFX5_n01_Oct_sp7797hop1_S288c_SK1_Yue-PM.txt',
         ref_input_counts='stats_H2CMNAFX5_n01_Oct_sp7797in_S288c_SK1_Yue-PM.txt',
         test_chip_counts='stats_H2CMNAFX5_n01_Oct_sp11644hop1_S288c_SK1_Yue-PM.txt',
         test_input_counts='stats_H2CMNAFX5_n01_Oct_sp11644in_S288c_SK1_Yue-PM_sorted.txt')
  )
)
loopNF_1 <- read_counts_loop2_1[2,2]
#[1] 0.6041732
read_counts_loop2_2 <- data.frame(
  Condition=c('WT', 'loop2'),
  NF=c(1,
       spikein_normalization_factor_from_counts(
         ref_chip_counts='counts_HVV2HAFX3_n01_nov2022_7797-Spk-hop1_S288c_SK1_Yue-PM.txt',
         ref_input_counts='counts_HVV2HAFX3_n01_nov2022_7797-Spk-in_S288c_SK1_Yue-PM.txt',
         test_chip_counts='counts_HVV2HAFX3_n01_nov2022_11644-Spk-hop1_S288c_SK1_Yue-PM.txt',
         test_input_counts='counts_HVV2HAFX3_n01_nov2022_11644-Spk-in_S288c_SK1_Yue-PM.txt')
  )
)
loopNF_2 <- read_counts_loop2_2[2,2]
#[1] 0.5963675
lpNF_list <- c(loopNF_1, loopNF_2)
lpNF <- mean(lpNF_list)
#[1] 0.6002704

lpN <- loopSig$score*lpNF
lpSigdf <- toDataframe(loopSig)
lpSigdf <- cbind(lpSigdf, lpN)
lpSigdf<- subset(lpSigdf, select = -c(score))
names(lpSigdf)[names(lpSigdf) == "lpN"] <- "score"
lpSigN <- makeGRangesFromDataFrame(lpSigdf, keep.extra.columns = TRUE)

sum(lpSigN$score)
#[1] 5260857
library(rtracklayer)
export.bed(lpSigN,con='Hop1_loop2_Reps_Norm.bed')


#Loop2-pch2_Norm
read_counts_loop2pch2_1 <- data.frame(
  Condition=c('WT', 'loop2pch2'),
  NF=c(1,
       spikein_normalization_factor_from_counts(
         ref_chip_counts='stats_H2CMNAFX5_n01_Oct_sp7797hop1_S288c_SK1_Yue-PM.txt',
         ref_input_counts='stats_H2CMNAFX5_n01_Oct_sp7797in_S288c_SK1_Yue-PM.txt',
         test_chip_counts='stats_H2CMNAFX5_n01_Oct_sp11757hop1_S288c_SK1_Yue-PM.txt',
         test_input_counts='stats_H2CMNAFX5_n01_Oct_sp11757in_S288c_SK1_Yue-PM.txt')
  )
)

read_counts_loop2pch2_2 <- data.frame(
  Condition=c('WT', 'loop2pch2'),
  NF=c(1,
       spikein_normalization_factor_from_counts(
         ref_chip_counts='stats_HT5W2AFX3_n01_7797e-spike_S288c_SK1_Yue-PM.txt',
         ref_input_counts='stats_HT5W2AFX3_n01_7797in-spike_S288c_SK1_Yue-PM.txt',
         test_chip_counts='stats_HT5W2AFX3_n01_11757e-spike_S288c_SK1_Yue-PM.txt',
         test_input_counts='stats_HT5W2AFX3_n01_11757in-spike_S288c_SK1_Yue-PM.txt')
  )
)

lpPchNF_1 <- read_counts_loop2pch2_1[2,2]
lpPchNF_2 <- read_counts_loop2pch2_2[2,2]
lpPchNF_list <- c(lpPchNF_1, lpPchNF_2)
lpPchNF <- mean(lpPchNF_list)
#[1] 0.688359
lpPchN <- lpPchSig$score*lpPchNF
lpPchSigdf <- toDataframe(lpPchSig)
lpPchSigdf <- cbind(lpPchSigdf, lpPchN)
lpPchSigdf<- subset(lpPchSigdf, select = -c(score))
names(lpPchSigdf)[names(lpPchSigdf) == "lpPchN"] <- "score"
lpPchSigN <- makeGRangesFromDataFrame(lpPchSigdf, keep.extra.columns = TRUE)
export.bed(lpPchSigN,con='Hop1_loop2pch2_Reps_Norm.bed')


#Pch2_norm
read_counts_pch2_1 <- data.frame(
  Condition=c('WT', 'pch2'),
  NF=c(1,
       spikein_normalization_factor_from_counts(
         ref_chip_counts='stats_H2CMNAFX5_n01_Oct_sp7797hop1_S288c_SK1_Yue-PM.txt',
         ref_input_counts='stats_H2CMNAFX5_n01_Oct_sp7797in_S288c_SK1_Yue-PM.txt',
         test_chip_counts='stats_H2CMNAFX5_n01_Oct_sp11758hop1_S288c_SK1_Yue-PM.txt',
         test_input_counts='stats_H2CMNAFX5_n01_Oct_sp11758in_S288c_SK1_Yue-PM.txt')
  )
)

read_counts_pch2_2 <- data.frame(
  Condition=c('WT', 'pch2'),
  NF=c(1,
       spikein_normalization_factor_from_counts(
         ref_chip_counts='stats_HT5W2AFX3_n01_7797e-spike_S288c_SK1_Yue-PM.txt',
         ref_input_counts='stats_HT5W2AFX3_n01_7797in-spike_S288c_SK1_Yue-PM.txt',
         test_chip_counts='stats_HT5W2AFX3_n01_11758e-spike_S288c_SK1_Yue-PM.txt',
         test_input_counts='stats_HT5W2AFX3_n01_11758in-spike_S288c_SK1_Yue-PM.txt')
  )
)

PchNF_1 <- read_counts_pch2_1[2,2]
PchNF_2 <- read_counts_pch2_2[2,2]
PchNF_list <- c(PchNF_1, PchNF_2)
PchNF <- mean(PchNF_list)
#[1] 0.9388861
PchN <- pchSig$score*PchNF
PchSigdf <- toDataframe(pchSig)
PchSigdf <- cbind(PchSigdf, PchN)
PchSigdf<- subset(PchSigdf, select = -c(score))
names(PchSigdf)[names(PchSigdf) == "PchN"] <- "score"
PchSigN <- makeGRangesFromDataFrame(PchSigdf, keep.extra.columns = TRUE)
library(rtracklayer)
export.bed(PchSigN,con='Hop1_pch2_reps_norm.bed')

Hop1read_counts_rec8_1 <- data.frame(
  Condition=c('WT', 'rec8'),
  NF=c(1,
       spikein_normalization_factor_from_counts(
         ref_chip_counts='counts_HNWF2AFX2_n01_AH7797-3h-anti-Hop1-SNP_S288c_SK1_Yue-PM.txt',
         ref_input_counts='counts_HNWF2AFX2_n01_AH7797-3h-input-SNP_S288c_SK1_Yue-PM.txt',
         test_chip_counts='counts_HNWF2AFX2_n01_AH5187-3h-anti-Hop1-SNP_S288c_SK1_Yue-PM.txt',
         test_input_counts='counts_HNWF2AFX2_n01_AH5187-3h-input-SNP_S288c_SK1_Yue-PM.txt')
  )
)
Hop1read_counts_rec8_2 <- data.frame(
  Condition=c('WT', 'rec8'),
  NF=c(1,
       spikein_normalization_factor_from_counts(
         ref_chip_counts='stats_H3GCHAFX3_n01_7797chp-hop1-NH-SPKIN_S288c_SK1_Yue-PM.txt',
         ref_input_counts='stats_H3GCHAFX3_n01_7797in-hop1-NH-SPKIN_S288c_SK1_Yue-PM.txt',
         test_chip_counts='stats_H3GCHAFX3_n01_5187chp-hop1-NH-SPKIN_S288c_SK1_Yue-PM.txt',
         test_input_counts='stats_H3GCHAFX3_n01_5187in-hop1-NH-SPKIN_S288c_SK1_Yue-PM.txt')
  )
)



rec8NF_1 <- Hop1read_counts_rec8_1[2,2]
rec8NF_2 <- Hop1read_counts_rec8_2[2,2]
rec8NF <- ((rec8NF_1 + rec8NF_2)/ 2)
rec8NF <- rec8Sig$score*rec8NF
rec8Sigdf <- toDataframe(rec8Sig)
rec8Sigdf <- cbind(rec8Sigdf, rec8NF)
rec8Sigdf <- subset(rec8Sigdf, select = -c(score))
names(rec8Sigdf)[names(rec8Sigdf) == "rec8NF"] <- "score"
rec8SigN <- makeGRangesFromDataFrame(rec8Sigdf, keep.extra.columns = TRUE)
export.bed(rec8SigN,con='Hop1sig_rec8D_norm.bed')


read_counts_rec8Phd <- data.frame(
  Condition=c('WT', 'rec8phd'),
  NF=c(1,
       spikein_normalization_factor_from_counts(
         ref_chip_counts='stats_H3GCHAFX3_n01_7797chp-hop1-NH-SPKIN_S288c_SK1_Yue-PM.txt',
         ref_input_counts='stats_H3GCHAFX3_n01_7797in-hop1-NH-SPKIN_S288c_SK1_Yue-PM.txt',
         test_chip_counts='stats_H3GCHAFX3_n01_10517chp-hop1-NH-SPKIN_S288c_SK1_Yue-PM.txt',
         test_input_counts='stats_H3GCHAFX3_n01_10517in-hop1-NH-SPKIN_S288c_SK1_Yue-PM.txt')
  )
)
rec8Phd_NF_1 <- read_counts_rec8Phd[2,2]
rec8loopNF <- rec8loop$score*rec8Phd_NF_1
rec8loopDF <- toDataframe(rec8loop)
rec8loopDF <- cbind(rec8loopDF, rec8loopNF)
rec8loopDF <- subset(rec8loopDF, select = -c(score))
names(rec8loopDF)[names(rec8loopDF) == "rec8loopNF"] <- "score"
rec8loopSigN <- makeGRangesFromDataFrame(rec8loopDF, keep.extra.columns = TRUE)
export.bed(rec8loopSigN,con='Hop1sig_rec8loopMut_norm2phd.bed')




data <- data.frame( name=c("1Wildtype","2hop1-loop2","2hop1-loop2","3pch2", "3pch2","4hop1-loop2 pch2", "4hop1-loop2 pch2"), 
                    value=c(100,60.43,59.63675,96.83214,90.94508,68.70668,68.96511))
#save these spike in normalized bedgraphs
ggplot2_theme <- theme_classic()
theme_set(ggplot2_theme)

p <- ggplot(data, aes(x=name, y=value))
p <- p + geom_bar(fill = '#0099CD', stat ="summary")
p <- p + geom_point(color='black',shape=6, size=5)
p













