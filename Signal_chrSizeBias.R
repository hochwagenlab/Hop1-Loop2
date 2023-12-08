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
library("ggpubr")
library("ggpmisc")






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

Hop1_rec8DFile <- paste0('/Users/darmokandjalad/Documents/Hop1PHD-Loop2/RAnalysis/Hop1-rec8D-61-91-Reps-SK1Yue-PM_B3W3_MACS2_FE.bdg.gz')
rec8Sig <- import_bedGraph(Hop1_rec8DFile)
rec8Sig <- sort(rec8Sig)

Red1Sig_rec8mut <- paste0('/Users/darmokandjalad/Documents/HISEQ-Illumina_SequencingResults/Average_Reps/Red1-rec8D-39-62-193-90-Reps-SK1Yue-PM_B3W3_MACS2_FE.bdg.gz')
Red1Sig_rec8mut <- import_bedGraph(Red1Sig_rec8mut)
Red1Sig_rec8mut <- sort(Red1Sig_rec8mut)

genome_info <- hwglabr2::get_chr_coordinates(genome = 'SK1Yue')
chrSize <- read.csv('/Users/darmokandjalad/Documents/Hop1PHD-Loop2/RAnalysis/chrSize.csv')

#Set file to look at
signal <- Wtsig

genAvg <- hwglabr2::average_chr_signal(signal)$genome_avrg
signal$score <- signal$score/genAvg

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
Chr_sig_log2 <- read.csv('/Users/darmokandjalad/Documents/Hop1PHD-Loop2/RAnalysis/chrLeng_Sig_log2_.csv')


#p <- ggplot(Chr_sig_log2, aes(x = length, y = sig, group = group, col = group ))
#p <- p + geom_point()
#p <- p + stat_poly_line(se = FALSE) 
#p <- p + stat_poly_eq()
#p
#p <- ggplot(Chr_sig_log2, aes(x = length, y = sig, group = group, col = group ))
#p <- p + geom_point()
#p <- p + stat_poly_line(se = FALSE) 
#p <- p + stat_poly_eq(aes(label = paste(after_stat(eq.label),
                                        after_stat(rr.label), sep = "*\", \"*")))
#p

p <- ggplot(data = Chr_sig_log2, aes(x = length, y = sig, group = group, col = group)) + 
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_color_manual(values = c("#0099CD","#B72467","#F68B1F","#CBDB2A","#F68B1F","#CBDB2A"))
p

#The stats are below


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

