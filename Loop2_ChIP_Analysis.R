# Everything For Hop1-loop2 Manuscript
#working folder is (On cal's computer) /Users/darmokandjalad/Documents/Hop1PHD-Loop2/RAnalysis
#Load Packages
############################
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
library(ggpmisc)

## Set Ggplot theme 
ggplot2_theme <- theme_classic() +
  theme(plot.title=element_text(hjust=0.5, size=10),
        plot.subtitle=element_text(hjust=0.5, size=10),
        axis.text=element_text(colour='black'),
        axis.ticks=element_line(colour='black'))

theme_set(ggplot2_theme)
#####################

#Set working Directory and Load Bedgraph files, these are not yet spike in normalized
###################################################
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

rec8File <- paste0('/Users/darmokandjalad/Documents/Hop1PHD-Loop2/RAnalysis/Rec8-wildtype-50-75-Reps-SK1Yue-PM_B3W3_MACS2_FE.bdg.gz')
rec8Sig <- import_bedGraph(rec8File)
rec8Sig <- sort(rec8Sig)


###################

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
read_counts_loop2 <- data.frame(
  Condition=c('WT', 'loop2'),
  NF=c(1,
       spikein_normalization_factor_from_counts(
         ref_chip_counts='stats_H2CMNAFX5_n01_Oct_sp7797hop1_S288c_SK1_Yue-PM.txt',
         ref_input_counts='stats_H2CMNAFX5_n01_Oct_sp7797in_S288c_SK1_Yue-PM.txt',
         test_chip_counts='stats_H2CMNAFX5_n01_Oct_sp11644hop1_S288c_SK1_Yue-PM.txt',
         test_input_counts='stats_H2CMNAFX5_n01_Oct_sp11644in_S288c_SK1_Yue-PM_sorted.txt')
  )
)
loopNF_1 <- read_counts_loop2[2,2]
lpN <- loopSig$score*loopNF_1
lpSigdf <- toDataframe(loopSig)
lpSigdf <- cbind(lpSigdf, lpN)
lpSigdf<- subset(lpSigdf, select = -c(score))
names(lpSigdf)[names(lpSigdf) == "lpN"] <- "score"
lpSigN <- makeGRangesFromDataFrame(lpSigdf, keep.extra.columns = TRUE)

sum(lpSigN$score)
#[1] 5295062



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
lpPchN <- lpPchSig$score*lpPchNF
lpPchSigdf <- toDataframe(lpPchSig)
lpPchSigdf <- cbind(lpPchSigdf, lpPchN)
lpPchSigdf<- subset(lpPchSigdf, select = -c(score))
names(lpPchSigdf)[names(lpPchSigdf) == "lpPchN"] <- "score"
lpPchSigN <- makeGRangesFromDataFrame(lpPchSigdf, keep.extra.columns = TRUE)



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
PchNF <- pchSig$score*PchNF
PchSigdf <- toDataframe(pchSig)
PchSigdf <- cbind(PchSigdf, PchNF)
PchSigdf<- subset(PchSigdf, select = -c(score))
names(PchSigdf)[names(PchSigdf) == "PchNF"] <- "score"
PchSigN <- makeGRangesFromDataFrame(PchSigdf, keep.extra.columns = TRUE)


data <- data.frame( name=c("1Wildtype","2hop1-loop2","3pch2", "3pch2","4hop1-loop2 pch2", "4hop1-loop2 pch2"), value=c(100,60.43,96.83214,90.94508,68.70668,68.96511))
#save these spike in normalized bedgraphs

p <- ggplot(data, aes(x=name, y=value))

p <- p + geom_bar(fill = 'gray', stat ="summary")
p <- p + geom_point(color='red',shape=6, size=5)
p

via <- read.csv('/Users/darmokandjalad/Documents/Hop1PHD-Loop2/SporeViability.csv')
p <- ggplot(via, aes(x=strain, y=mean)) + geom_bar(fill="grey", stat = "summary")
p 
p <- p + geom_errorbar(aes(ymin=mean - standDev, ymax=value + standDev)) 
p

via <- read.csv('/Users/darmokandjalad/Documents/Hop1PHD-Loop2/SporeViability_summary.csv')
p <- ggplot(via, aes(x=strain, y=value)) + geom_bar(fill="grey", stat = "summary")
p 
p <- p + geom_errorbar(aes(ymin=mean - standDev, ymax=value + standDev)) 
p


df<-data.frame(Mean=c(94.33333333, 89.56, 96.14, 95.68, 41.62, 74.81666667),
               sd=c(4.077172223,2.707951255,2.07557221,2.735324478,13.06797612,5.482852056),
               Category=c("5tel1", "6tel1 hop1 loop2", "1Wild type", "3pch2", "4pch2 hop1 loop2", "2hop1 loop2")
)

# Load ggplot2
library(ggplot2)

ggplot(df, aes(x=Category, y=Mean, fill="blue")) +
  geom_bar(position=position_dodge(), stat="identity",
           colour='black') +
  geom_errorbar(aes(ymin=Mean-sd, ymax=Mean+sd), width=.2)



p

lpSigN
PchSigN
lpPchSigN

#Loop2NF = 0.6041732
#Loop2 pch2 NF = 0.688359
#pch2 NF = 0.9388861



############################


#Plot signal along chromosomes in one plot

#Load and Align Normalized signal along chromosomes *Normalized bedgraphs to compare Y axes
###########################################################################
#Normalized Grange Objects lpSigN = loop2, PchSigN = pch2, lpPchSigN = loop2 pch2

cen10kb = rtracklayer::import.bed('/Users/darmokandjalad/Documents/Hop1PHD-Loop2/RAnalysis/Cen10KbFlankRegions.bed')
#set Parameters
signal <- lpPchSigN
chrnum = "chrV"
antibodyTarget = "Hop1 Signal"
genotype = "hop1-loop2 pch2  (Nancy Hollingsworth Anti-Hop1)"

#Get Genome and centromeres
genome_info <- hwglabr2::get_chr_coordinates(genome = 'SK1Yue')
centromeres <- hwglabr2::get_chr_coordinates(genome='SK1Yue', as_df=FALSE)
cen = centromeres[seqnames(centromeres) == chrnum]
cen <- round((start(cen) + end(cen))/2)
cen <- data.frame(cen = cen, y = 0)
cen <- cen/1000


clusters <- rtracklayer::import.bed('/Users/darmokandjalad/Documents/Hop1PHD-Loop2/RAnalysis/clusters_joined.bed')
clsters <- clusters[seqnames(clusters) == chrnum]
clstersDF <- data.frame(clsters)

# Sort sequences and levels to make sure they match
signal <- sort(GenomeInfoDb::sortSeqlevels(signal))
genome_info <- GenomeInfoDb::sortSeqlevels(genome_info)

# Add info to signal object
GenomeInfoDb::seqlengths(signal) <- GenomeInfoDb::seqlengths(genome_info)
# Compute 100-bp tiling windows
bins <- GenomicRanges::tileGenome(GenomeInfoDb::seqlengths(signal),
                                  tilewidth=100,
                                  cut.last.tile.in.chrom=TRUE)
# Get signal as "RleList"; the signal is stored in the "score" metadata column
score <- GenomicRanges::coverage(signal, weight="score")
# Compute average signal per tile
bins <- GenomicRanges::binnedAverage(bins, score, "binned_score")
# Keep only chr V
bins <- GenomeInfoDb::keepSeqlevels(bins, chrnum, pruning.mode = "coarse")

# Get positions as the midpoints of the intervals
positions <- bins@ranges@start + floor(bins@ranges@width / 2)

# Make data frame (convert positions to Kb; signal is the binned score)
df <- data.frame(position=positions / 1000, signal=bins$binned_score)



#plot(df, xlab='Position on chrV (Kb)', ylab='Average signal', type='l', lwd=2, main="rec8D 5187", ylim = c(0,8))
p <- ggplot(df, aes(x=position, y=signal)) + ylim(0,8) + geom_line()
p <- p + geom_point(cen, mapping=aes(cen,0),
                    size = 2.0, colour = 'red') + ggtitle(genotype) + xlab(chrnum) + ylab(antibodyTarget)
p <- p + geom_segment(clstersDF, mapping=aes(x=start / 1000, xend=end /1000, y=8, yend=8),colour = 'green')

p



#save Plot


wtSigN_Xp <- p
lpSigN_Xp <- p
PchSigN_Xp <- p
lpPchSigN_Xp <- p

ChrX_All <- wtSigN_Xp / lpSigN_Xp / PchSigN_Xp / lpPchSigN_Xp

wtSigN_Vp <- p
lpSigN_Vp <- p
PchSigN_Vp <- p
lpPchSigN_Vp <- p

chrV_All <- wtSigN_Vp / lpSigN_Vp / PchSigN_Vp / lpPchSigN_Vp

ChrX_All <- wtSigN_Xp / lpSigN_Xp / pch2SigN_Xp / lppch2SigN_Xp
ChrX_All

ChrV_All <- wtSigN_Vp / lpSigN_Vp / PchSigN_Vp / lpPchSigN_Vp
ChrV_All

ChrXV_All <- wtSigN_XVp / lpSigN_XVp / PchSigN_XVp / lpPchSigN_XVp
ChrXV_All

#####################################################
#Quantify signal in islands vs non-islands
lpSigN
Wtsig
setwd('/Users/darmokandjalad/Documents/Hop1PHD-Loop2/RAnalysis')
clusters = rtracklayer::import.bed('clusters_joined.bed')
deserts = rtracklayer::import.bed('deserts_joined.bed')

axis = hwglabr2::get_Red1_summits("SK1Yue")

 lpSigN
 mean(lpSigN$score)
 #[1] 0.5730923
 sum(lpSigN$score)
 #[1] 5295062

hits = findOverlaps(clusters,lpSigN)
lpSigN_cluster <- lpSigN[subjectHits(hits)]
rm(hits)
mean(lpSigN_cluster$score)
  #[1] 0.636009
sum(lpSigN_cluster$score)
#[1] 974707.3


hits = overlapsAny(deserts,lpSigN, type = "within")
lpSigN_desert <- subsetByOverlaps(lpSigN, deserts)
rm(hits)
  mean(lpSigN_desert$score)
  #[1] 0.5605792
  sum(lpSigN_desert$score)
 # [1] 4320265

Wtsig
mean(Wtsig$score)
  #[1] 0.9354605
sum(Wtsig$score)
  #[1] 10266761

hits = findOverlaps(clusters,Wtsig)
Wtsig_cluster <- Wtsig[subjectHits(hits)]
rm(hits)
mean(Wtsig_cluster$score)
  #[1] 1.315606
sum(Wtsig_cluster$score)
  #[1] 2615032

Wtsig_desert <- subsetByOverlaps(Wtsig, deserts)
mean(lpSigN_desert$score)
  #0.5605792
sum(Wtsig_desert$score)
  #[1] 7651484

#OK but clusters are only a fraction of the genome, how large is the fraction
wdithClusters <- width(clusters)
head(wdithClusters)
#[1]  4999 64999 29999 24999  4999 14999
sum(wdithClusters)
#[1] 2149874

wdithDesert <- width(deserts)
head(wdithDesert)
sum(wdithDesert)
#[1] 9913145








loopSig
mean(loopSig$score)
      #[1] 0.9485562
hits = findOverlaps(clusters,loopSig)
loopSig_cluster <- loopSig[subjectHits(hits)]
rm(hits)
mean(loopSig_cluster$score)
  #1.052693
loopSig_desert <- subsetByOverlaps(loopSig, deserts)
mean(loopSig_desert$score)
  #0.9278451

#Qunatify enrichment of signal in island or non-islands* double check this. Signal compared after normalizing each to genome average
################################################################################################
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

rec8File <- paste0('/Users/darmokandjalad/Documents/Hop1PHD-Loop2/RAnalysis/Rec8-wildtype-50-75-Reps-SK1Yue-PM_B3W3_MACS2_FE.bdg.gz')
rec8Sig <- import_bedGraph(rec8File)
rec8Sig <- sort(rec8Sig)

signalfile <- rec8Sig


clusters = rtracklayer::import.bed('clusters_joined.bed')
deserts = rtracklayer::import.bed('deserts_joined.bed')
axis = hwglabr2::get_Red1_summits("SK1Yue")

hits = findOverlaps(clusters,axis)
axis_cluster <- axis[subjectHits(hits)]
rm(hits)
hits = findOverlaps(deserts,axis)
axis_desert <- axis[subjectHits(hits)]
rm(hits)
subset(axis_cluster, (name %in% axis_desert$name))
mcols(axis_desert)['class'] = 'desert'
mcols(axis_cluster)['class'] = 'cluster'


genAvg <- hwglabr2::average_chr_signal(signalfile)$genome_avrg
signalfile$score <- signalfile$score/genAvg

mat1 <- normalizeToMatrix(signalfile, axis_cluster, value_column = "score",
                          extend = 1000, mean_mode = "weighted", w = 5,empty_value=NA)

mat1_avrg <- hwglabr2::signal_mean_and_ci(mat1,
                                          ci=0.95, rep_bootstrap=1000,
                                          na_rm=TRUE)

mat2 <- normalizeToMatrix(signalfile, axis_desert, value_column = "score",
                          extend = 1000, mean_mode = "weighted", w = 5,empty_value=NA)

mat2_avrg <- hwglabr2::signal_mean_and_ci(mat2,
                                          ci=0.95, rep_bootstrap=1000,
                                          na_rm=TRUE)

mat1_avrg_df <- data.frame(Position=seq(-199, 200), sample = 'cluster',mat1_avrg)
mat2_avrg_df <- data.frame(Position=seq(-199, 200), sample = 'desert',mat2_avrg)
alldata = rbind(mat1_avrg_df,mat2_avrg_df)

p <- ggplot(alldata, aes(Position, Mean, group = sample, colour=sample, fill=sample)) +
  labs(title = "Signal at Red1 summits",
       x = "Distance to axis (bp)", y = "Average\nChIP-seq signal") +
  geom_vline(xintercept = 0, lty = 3) +
  scale_x_continuous(breaks = c(-199, 0, 200),
                     labels = c("-1 kb", "summit", "1 kb"))

# Add confidence interval as a ribbon
p <- p + geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha=0.3,colour=NA)
p <- p + geom_line()
p

pWt_enrichIslands <- p
pLoop2_enrichIslands <- p
pPch2_enrichIslands <- p
ploop2_Pch2_enrichIslands <- p

#########################################################


#Qunatify enrichment of signal in Centromere or centromere
################################################################################################
setwd('/Users/darmokandjalad/Documents/Hop1PHD-Loop2/RAnalysis')
#Normalized Grange Objects lpSigN = loop2, PchSigN = pch2, lpPchSigN = loop2 pch2
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

rec8File <- paste0('/Users/darmokandjalad/Documents/Hop1PHD-Loop2/RAnalysis/Rec8-wildtype-50-75-Reps-SK1Yue-PM_B3W3_MACS2_FE.bdg.gz')
rec8Sig <- import_bedGraph(rec8File)
rec8Sig <- sort(rec8Sig)

#Put signal over genome average
#pick signal to analyzie
X = rec8Sig
Cen_flank = 10000

signal <- X
genAvg <- hwglabr2::average_chr_signal(signal)$genome_avrg
signal$score <- signal$score/genAvg

#get genome
genome_info <- hwglabr2::get_chr_coordinates(genome = 'SK1Yue')

# Sort sequences and levels to make sure they match
signal <- sort(GenomeInfoDb::sortSeqlevels(signal))
genome_info <- GenomeInfoDb::sortSeqlevels(genome_info)

# isolate centromere regions
#how far out from cen


binsc <- toDataframe(genome_info)
cenStartminus <- binsc$start - Cen_flank
cenEndplus <- binsc$end + Cen_flank
chr <- binsc$chr
Cens_bpFlank <- data.frame(first_column = chr, second_column= cenStartminus, third_column = cenEndplus)
Cens_bpFlank <- toGRanges(Cens_bpFlank)
df <- toDataframe(Cens_bpFlank)
write.table(df, file="Cen10KbFlankRegions.bed", quote=F, sep="\t", row.names=F, col.names=F)



axis = hwglabr2::get_Red1_summits("SK1Yue")

hits = findOverlaps(Cens_bpFlank,axis)
axis_cens <- axis[subjectHits(hits)]
rm(hits)
axis_offcen <- axis[axis %outside% axis_cens]




mat1 <- normalizeToMatrix(signal, axis_cens, value_column = "score",
                          extend = 1000, mean_mode = "weighted", w = 5,empty_value=NA)

mat1_avrg <- hwglabr2::signal_mean_and_ci(mat1,
                                          ci=0.95, rep_bootstrap=1000,
                                          na_rm=TRUE)

mat2 <- normalizeToMatrix(signal, axis_offcen, value_column = "score",
                          extend = 1000, mean_mode = "weighted", w = 5,empty_value=NA)

mat2_avrg <- hwglabr2::signal_mean_and_ci(mat2,
                                          ci=0.95, rep_bootstrap=1000,
                                          na_rm=TRUE)

mat1_avrg_df <- data.frame(mat1_avrg)
row = c(-199:200)
mat1_avrg_df['Position'] = row
mat1_avrg_df['sample'] = 'On_cen'


mat2_avrg_df <- data.frame(mat2_avrg)
mat2_avrg_df['Position'] = row
mat2_avrg_df['sample'] = 'Off_cen'
alldata = rbind(mat1_avrg_df,mat2_avrg_df)


p <- ggplot(alldata, aes(Position, Mean, group = sample, colour=sample, fill=sample)) +
  labs(title = "Signal at Centromere or chr Arm",
       x = "Distance to axis site (bp)", y = "Average\nChIP-seq signal") +
  geom_vline(xintercept = 0, lty = 3) +
  scale_x_continuous(breaks = c(-199, 0, 200),
                     labels = c("-1 kb", "summit", "1 kb"))

# Add confidence interval as a ribbon
p <- p + geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha=0.3,colour=NA)
p <- p + geom_line()
p

pWt_Centromere_10000kb <- p
pLoop2_enrichCentromere_10000kb <- p
pPch2_enrichCentromere_10000kb <- p
ploop2_Pch2_enrichCentromere_10000kb <- p
pLoop2_enrichCentromere_10000kb_rep <- p
ploop2_Pch2_enrichCentromere_10000kb_rep <- p
Pch2_enrichCentromere_10000kb_rep <- p



########################################

#plot average signal per chromsomes to get chromosomes size effects
#######################################
setwd('/Users/darmokandjalad/Documents/Hop1PHD-Loop2/RAnalysis')
#Normalized Grange Objects lpSigN = loop2, PchSigN = pch2, lpPchSigN = loop2 pch2
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

rec8File <- paste0('/Users/darmokandjalad/Documents/Hop1PHD-Loop2/RAnalysis/Rec8-wildtype-50-75-Reps-SK1Yue-PM_B3W3_MACS2_FE.bdg.gz')
rec8Sig <- import_bedGraph(rec8File)
rec8Sig <- sort(rec8Sig)

genome_info <- hwglabr2::get_chr_coordinates(genome = 'SK1Yue')
chrSize <- read.csv('/Users/darmokandjalad/Documents/Hop1PHD-Loop2/RAnalysis/chrSize.csv')
#Set file to look at
signal <- rec8Sig

# Sort sequences and levels to make sure they match
signal <- sort(GenomeInfoDb::sortSeqlevels(signal))
genome_info <- GenomeInfoDb::sortSeqlevels(genome_info)

# Add info to signal object
GenomeInfoDb::seqlengths(signal) <- GenomeInfoDb::seqlengths(genome_info)
signal <- sort(GenomeInfoDb::sortSeqlevels(signal))

#calculate the signal per chromsomes
chrSig <- hwglabr2::average_chr_signal(signal)
chrSig_df <- toDataframe(chrSig$seq_avrg)

chrLength <- chrSize$length
chrSig_df <- cbind(chrSig_df, chrLength)

colnames(chrSig_df) <- c('chr', 'wt_avg_sig', 'chrlength')


p <- ggplot(chrSig_df, aes(x=chrLength, y = wt_avg_sig))
p <- p + geom_point()
p


#I just saved the csv because it felt easier than figuring out how to get the data frame together, although i'm sure its not actually that hard. 
#import file 
Chr_sig_log2 <- read.csv('/Users/darmokandjalad/Documents/Hop1PHD-Loop2/RAnalysis/chrLeng_Sig_log2.csv')


p <- ggplot(Chr_sig_log2, aes(x = length, y = sig, group = group, col = group ))
p <- p + geom_point()
p <- p + stat_poly_line(se = FALSE) 
p <- p + stat_poly_eq()
p
p <- ggplot(Chr_sig_log2, aes(x = length, y = sig, group = group, col = group ))
p <- p + geom_point()
p <- p + stat_poly_line(se = FALSE) 
p <- p + stat_poly_eq(aes(label = paste(after_stat(eq.label),
                                        after_stat(rr.label), sep = "*\", \"*")))
p



wtChrSig <- subset(Chr_sig_log2, group == 'wt')
cor(wtChrSig$length, wtChrSig$sig  , method = "pearson", use = "complete.obs")
  #[1] -0.8966651
lpChrSig <- subset(Chr_sig_log2, group == 'loop2')
cor(lpChrSig$length, lpChrSig$sig  , method = "pearson", use = "complete.obs")
  #[1] -0.5828771
pchChrSig <- subset(Chr_sig_log2, group == 'pch2')
cor(pchChrSig$length, pchChrSig$sig  , method = "pearson", use = "complete.obs")
  #[1] -0.8411006

lpPchChrSig <- subset(Chr_sig_log2, group == 'loop2pch2')
cor(lpPchChrSig$length, lpPchChrSig$sig  , method = "pearson", use = "complete.obs")
  #[1] -0.611651

rec8ChrSig <- subset(Chr_sig_log2, group == 'rec8')
cor(rec8ChrSig$length, rec8ChrSig$sig  , method = "pearson", use = "complete.obs")
  #[1] -0.8594557

library("ggpubr")
ggscatter(wtChrSig, x = "length", y = "sig", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "chrsize", ylab = "hop1sig")
#R = -0.9, p = 2.6e-06
ggscatter(lpChrSig, x = "length", y = "sig", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "chrsize", ylab = "hop1sig")
#R = -0.51, p = 0.018
ggscatter(pchChrSig, x = "length", y = "sig", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "chrsize", ylab = "hop1sig")
#R = -0.84, p = 4.5e-05
ggscatter(lpPchChrSig, x = "length", y = "sig", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "chrsize", ylab = "hop1sig")
#R = -0.61, p = 0.012
ggscatter(rec8ChrSig, x = "length", y = "sig", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "chrsize", ylab = "hop1sig")
#R = -0.86, p = 2e-05
ggscatter(Chr_sig_log2, x="length", y="sig", combine = TRUE, color = "group",
          add = "reg.line", cof.int = TRUE,
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Chromosome Size", ylab = "Hop1 Signal")

############################################################




SK1 = GRanges(seqnames=c("chrI","chrII","chrIII","chrIV","chrV","chrVI","chrVII","chrVIII","chrIX","chrX","chrXI","chrXII","chrXIII","chrXIV","chrXV","chrXVI"),
              ranges = IRanges(start=c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1),
                               end=c(203894,794509,342719,1490683,602515,284457,1067527,544539,435586,719295,687261,1008249,908608,812466,1054034,921189)
              ))
cens = GRanges(seqnames = c("chrI","chrII","chrIII","chrIV","chrV","chrVI","chrVII","chrVIII","chrIX","chrX","chrXI","chrXII","chrXIII","chrXIV","chrXV","chrXVI"),
               ranges = IRanges(start=c),
               end=c(154629,251816,108709,460753,171137,170911,501252,
               ))













