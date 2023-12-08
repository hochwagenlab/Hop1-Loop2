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

ggplot2_blanktheme <- theme_classic() +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank(),
        title = NULL)
theme_set(ggplot2_blanktheme)




###############################
#11570_A1
dsbSig_a1 <- read.table('/Users/darmokandjalad/Desktop/Loop2_Neale/Hochwagencollab2023/HpM/HpM.ASM205788v1Spike-RA77_11570_rad50S_6h_A1.txt', header = TRUE)
both <- dsbSig_a1$Watson + dsbSig_a1$Crick
dsbSig_a1$both <- both


#start <- dsbSig$Pos - 249
#dsbSig$start <- start
#dsbSig[dsbSig == "Chr"] <- "seqnames"
dsbSig_a1["Chr"][dsbSig_a1["Chr"] == 1] <- "chrI"
dsbSig_a1["Chr"][dsbSig_a1["Chr"] == 2] <- "chrII"
dsbSig_a1["Chr"][dsbSig_a1["Chr"] == 3] <- "chrIII"
dsbSig_a1["Chr"][dsbSig_a1["Chr"] == 4] <- "chrIV"
dsbSig_a1["Chr"][dsbSig_a1["Chr"] == 5] <- "chrV"
dsbSig_a1["Chr"][dsbSig_a1["Chr"] == 6] <- "chrVI"
dsbSig_a1["Chr"][dsbSig_a1["Chr"] == 7] <- "chrVII"
dsbSig_a1["Chr"][dsbSig_a1["Chr"] == 8] <- "chrVIII"
dsbSig_a1["Chr"][dsbSig_a1["Chr"] == 9] <- "chrIX"
dsbSig_a1["Chr"][dsbSig_a1["Chr"] == 10] <- "chrX"
dsbSig_a1["Chr"][dsbSig_a1["Chr"] == 11] <- "chrXI"
dsbSig_a1["Chr"][dsbSig_a1["Chr"] == 12] <- "chrXII"
dsbSig_a1["Chr"][dsbSig_a1["Chr"] == 13] <- "chrXIII"
dsbSig_a1["Chr"][dsbSig_a1["Chr"] == 14] <- "chrXIV"
dsbSig_a1["Chr"][dsbSig_a1["Chr"] == 15] <- "chrXV"
dsbSig_a1["Chr"][dsbSig_a1["Chr"] == 16] <- "chrXVI"


colnames(dsbSig_a1) <- c('seqnames','Pos','Watson','Crick','score')              
dsbSig_a1 <- makeGRangesFromDataFrame(dsbSig_a1, keep.extra.columns = TRUE, seqnames.field = "seqnames", start.field = "Pos", end.field = "Pos")              
dsbSig_a1 <- dropSeqlevels(dsbSig_a1, "chrIII", pruning.mode = "coarse")

setwd('/Users/darmokandjalad/Documents/Hop1PHD-Loop2/RAnalysis')
export.bed(dsbSig_a1,con='ccSeq_dsbsig_11570_A1.bed')

#bin 100bp
genome_info <- hwglabr2::get_chr_coordinates(genome = 'SK1Yue')
signal <- dsbSig_a1
# Sort sequences and levels to make sure they match
signal <- sort(GenomeInfoDb::sortSeqlevels(signal))
genome_info <- GenomeInfoDb::sortSeqlevels(genome_info)
genome_info <- dropSeqlevels(genome_info, "chrIII", pruning.mode = "coarse")
# Add info to signal object
GenomeInfoDb::seqlengths(signal) <- GenomeInfoDb::seqlengths(genome_info)
# Compute 100-bp tiling windows
bins <- GenomicRanges::tileGenome(GenomeInfoDb::seqlengths(signal),
                                  tilewidth=100,
                                  cut.last.tile.in.chrom=TRUE)
# Get signal as "RleList"; the signal is stored in the "score" metadata column
score <- GenomicRanges::coverage(signal, weight="score")
# Compute average signal per tile
dsbSig_a1 <- GenomicRanges::binnedAverage(bins, score, "binned_score")

###########11570_A2
dsbSig_a2 <- read.table('/Users/darmokandjalad/Desktop/Loop2_Neale/Hochwagencollab2023/HpM/HpM.ASM205788v1Spike-RA78_11570_rad50S_6h_B1.txt', header = TRUE)
both <- dsbSig_a2$Watson + dsbSig_a2$Crick
dsbSig_a2$both <- both


#start <- dsbSig_a2$Pos - 249
#dsbSig_a2$start <- start
#dsbSig_a2[dsbSig_a2 == "Chr"] <- "seqnames"
dsbSig_a2["Chr"][dsbSig_a2["Chr"] == 1] <- "chrI"
dsbSig_a2["Chr"][dsbSig_a2["Chr"] == 2] <- "chrII"
dsbSig_a2["Chr"][dsbSig_a2["Chr"] == 3] <- "chrIII"
dsbSig_a2["Chr"][dsbSig_a2["Chr"] == 4] <- "chrIV"
dsbSig_a2["Chr"][dsbSig_a2["Chr"] == 5] <- "chrV"
dsbSig_a2["Chr"][dsbSig_a2["Chr"] == 6] <- "chrVI"
dsbSig_a2["Chr"][dsbSig_a2["Chr"] == 7] <- "chrVII"
dsbSig_a2["Chr"][dsbSig_a2["Chr"] == 8] <- "chrVIII"
dsbSig_a2["Chr"][dsbSig_a2["Chr"] == 9] <- "chrIX"
dsbSig_a2["Chr"][dsbSig_a2["Chr"] == 10] <- "chrX"
dsbSig_a2["Chr"][dsbSig_a2["Chr"] == 11] <- "chrXI"
dsbSig_a2["Chr"][dsbSig_a2["Chr"] == 12] <- "chrXII"
dsbSig_a2["Chr"][dsbSig_a2["Chr"] == 13] <- "chrXIII"
dsbSig_a2["Chr"][dsbSig_a2["Chr"] == 14] <- "chrXIV"
dsbSig_a2["Chr"][dsbSig_a2["Chr"] == 15] <- "chrXV"
dsbSig_a2["Chr"][dsbSig_a2["Chr"] == 16] <- "chrXVI"


colnames(dsbSig_a2) <- c('seqnames','Pos','Watson','Crick','score')              
dsbSig_a2 <- makeGRangesFromDataFrame(dsbSig_a2, keep.extra.columns = TRUE, seqnames.field = "seqnames", start.field = "Pos", end.field = "Pos")              
dsbSig_a2 <- dropSeqlevels(dsbSig_a2, "chrIII", pruning.mode = "coarse")

setwd('/Users/darmokandjalad/Documents/Hop1PHD-Loop2/RAnalysis')
export.bed(dsbSig_a2,con='ccSeq_dsbsig_11570_A2.bed')

#bin 100bp
genome_info <- hwglabr2::get_chr_coordinates(genome = 'SK1Yue')
signal <- dsbSig_a2
# Sort sequences and levels to make sure they match
signal <- sort(GenomeInfoDb::sortSeqlevels(signal))
genome_info <- GenomeInfoDb::sortSeqlevels(genome_info)
genome_info <- dropSeqlevels(genome_info, "chrIII", pruning.mode = "coarse")
# Add info to signal object
GenomeInfoDb::seqlengths(signal) <- GenomeInfoDb::seqlengths(genome_info)
# Compute 100-bp tiling windows
bins <- GenomicRanges::tileGenome(GenomeInfoDb::seqlengths(signal),
                                  tilewidth=100,
                                  cut.last.tile.in.chrom=TRUE)
# Get signal as "RleList"; the signal is stored in the "score" metadata column
score <- GenomicRanges::coverage(signal, weight="score")
# Compute average signal per tile
dsbSig_a2 <- GenomicRanges::binnedAverage(bins, score, "binned_score")

########combine 11570A1 and A2
a1 <- dsbSig_a1$binned_score
a2 <- dsbSig_a2$binned_score
both <- a1+a2
dsbSig_a1$score <- both 
dsbSig11570 <- dsbSig_a1
export.bed(dsbSig11570,con='ccSeq_dsbsig_11570_reps.bed')

#####################################
#####11687
#11687_A1
dsbSig_a1 <- read.table('/Users/darmokandjalad/Desktop/Loop2_Neale/Hochwagencollab2023/HpM/HpM.ASM205788v1Spike-RA79_11687_rad50Shop1-loop2_6h_A1.txt', header = TRUE)
both <- dsbSig_a1$Watson + dsbSig_a1$Crick
dsbSig_a1$both <- both


#start <- dsbSig$Pos - 249
#dsbSig$start <- start
#dsbSig[dsbSig == "Chr"] <- "seqnames"
dsbSig_a1["Chr"][dsbSig_a1["Chr"] == 1] <- "chrI"
dsbSig_a1["Chr"][dsbSig_a1["Chr"] == 2] <- "chrII"
dsbSig_a1["Chr"][dsbSig_a1["Chr"] == 3] <- "chrIII"
dsbSig_a1["Chr"][dsbSig_a1["Chr"] == 4] <- "chrIV"
dsbSig_a1["Chr"][dsbSig_a1["Chr"] == 5] <- "chrV"
dsbSig_a1["Chr"][dsbSig_a1["Chr"] == 6] <- "chrVI"
dsbSig_a1["Chr"][dsbSig_a1["Chr"] == 7] <- "chrVII"
dsbSig_a1["Chr"][dsbSig_a1["Chr"] == 8] <- "chrVIII"
dsbSig_a1["Chr"][dsbSig_a1["Chr"] == 9] <- "chrIX"
dsbSig_a1["Chr"][dsbSig_a1["Chr"] == 10] <- "chrX"
dsbSig_a1["Chr"][dsbSig_a1["Chr"] == 11] <- "chrXI"
dsbSig_a1["Chr"][dsbSig_a1["Chr"] == 12] <- "chrXII"
dsbSig_a1["Chr"][dsbSig_a1["Chr"] == 13] <- "chrXIII"
dsbSig_a1["Chr"][dsbSig_a1["Chr"] == 14] <- "chrXIV"
dsbSig_a1["Chr"][dsbSig_a1["Chr"] == 15] <- "chrXV"
dsbSig_a1["Chr"][dsbSig_a1["Chr"] == 16] <- "chrXVI"


colnames(dsbSig_a1) <- c('seqnames','Pos','Watson','Crick','score')              
dsbSig_a1 <- makeGRangesFromDataFrame(dsbSig_a1, keep.extra.columns = TRUE, seqnames.field = "seqnames", start.field = "Pos", end.field = "Pos")              
dsbSig_a1 <- dropSeqlevels(dsbSig_a1, "chrIII", pruning.mode = "coarse")

setwd('/Users/darmokandjalad/Documents/Hop1PHD-Loop2/RAnalysis')
export.bed(dsbSig_a1,con='ccSeq_dsbsig_11687_A1.bed')

#bin 100bp
genome_info <- hwglabr2::get_chr_coordinates(genome = 'SK1Yue')
signal <- dsbSig_a1
# Sort sequences and levels to make sure they match
signal <- sort(GenomeInfoDb::sortSeqlevels(signal))
genome_info <- GenomeInfoDb::sortSeqlevels(genome_info)
genome_info <- dropSeqlevels(genome_info, "chrIII", pruning.mode = "coarse")
# Add info to signal object
GenomeInfoDb::seqlengths(signal) <- GenomeInfoDb::seqlengths(genome_info)
# Compute 100-bp tiling windows
bins <- GenomicRanges::tileGenome(GenomeInfoDb::seqlengths(signal),
                                  tilewidth=100,
                                  cut.last.tile.in.chrom=TRUE)
# Get signal as "RleList"; the signal is stored in the "score" metadata column
score <- GenomicRanges::coverage(signal, weight="score")
# Compute average signal per tile
dsbSig_a1 <- GenomicRanges::binnedAverage(bins, score, "binned_score")

###########11687_A2
dsbSig_a2 <- read.table('/Users/darmokandjalad/Desktop/Loop2_Neale/Hochwagencollab2023/HpM/HpM.ASM205788v1Spike-RA80_11687_rad50Shop1-loop2_6h_B1.txt', header = TRUE)
both <- dsbSig_a2$Watson + dsbSig_a2$Crick
dsbSig_a2$both <- both


#start <- dsbSig_a2$Pos - 249
#dsbSig_a2$start <- start
#dsbSig_a2[dsbSig_a2 == "Chr"] <- "seqnames"
dsbSig_a2["Chr"][dsbSig_a2["Chr"] == 1] <- "chrI"
dsbSig_a2["Chr"][dsbSig_a2["Chr"] == 2] <- "chrII"
dsbSig_a2["Chr"][dsbSig_a2["Chr"] == 3] <- "chrIII"
dsbSig_a2["Chr"][dsbSig_a2["Chr"] == 4] <- "chrIV"
dsbSig_a2["Chr"][dsbSig_a2["Chr"] == 5] <- "chrV"
dsbSig_a2["Chr"][dsbSig_a2["Chr"] == 6] <- "chrVI"
dsbSig_a2["Chr"][dsbSig_a2["Chr"] == 7] <- "chrVII"
dsbSig_a2["Chr"][dsbSig_a2["Chr"] == 8] <- "chrVIII"
dsbSig_a2["Chr"][dsbSig_a2["Chr"] == 9] <- "chrIX"
dsbSig_a2["Chr"][dsbSig_a2["Chr"] == 10] <- "chrX"
dsbSig_a2["Chr"][dsbSig_a2["Chr"] == 11] <- "chrXI"
dsbSig_a2["Chr"][dsbSig_a2["Chr"] == 12] <- "chrXII"
dsbSig_a2["Chr"][dsbSig_a2["Chr"] == 13] <- "chrXIII"
dsbSig_a2["Chr"][dsbSig_a2["Chr"] == 14] <- "chrXIV"
dsbSig_a2["Chr"][dsbSig_a2["Chr"] == 15] <- "chrXV"
dsbSig_a2["Chr"][dsbSig_a2["Chr"] == 16] <- "chrXVI"


colnames(dsbSig_a2) <- c('seqnames','Pos','Watson','Crick','score')              
dsbSig_a2 <- makeGRangesFromDataFrame(dsbSig_a2, keep.extra.columns = TRUE, seqnames.field = "seqnames", start.field = "Pos", end.field = "Pos")              
dsbSig_a2 <- dropSeqlevels(dsbSig_a2, "chrIII", pruning.mode = "coarse")

setwd('/Users/darmokandjalad/Documents/Hop1PHD-Loop2/RAnalysis')
export.bed(dsbSig_a2,con='ccSeq_dsbsig_11687_A2.bed')

#bin 100bp
genome_info <- hwglabr2::get_chr_coordinates(genome = 'SK1Yue')
signal <- dsbSig_a2
# Sort sequences and levels to make sure they match
signal <- sort(GenomeInfoDb::sortSeqlevels(signal))
genome_info <- GenomeInfoDb::sortSeqlevels(genome_info)
genome_info <- dropSeqlevels(genome_info, "chrIII", pruning.mode = "coarse")
# Add info to signal object
GenomeInfoDb::seqlengths(signal) <- GenomeInfoDb::seqlengths(genome_info)
# Compute 100-bp tiling windows
bins <- GenomicRanges::tileGenome(GenomeInfoDb::seqlengths(signal),
                                  tilewidth=100,
                                  cut.last.tile.in.chrom=TRUE)
# Get signal as "RleList"; the signal is stored in the "score" metadata column
score <- GenomicRanges::coverage(signal, weight="score")
# Compute average signal per tile
dsbSig_a2 <- GenomicRanges::binnedAverage(bins, score, "binned_score")

########combine 11570A1 and A2
a1 <- dsbSig_a1$binned_score
a2 <- dsbSig_a2$binned_score
both <- a1+a2
dsbSig_a1$score <- both 
dsbSig11570 <- dsbSig_a1
export.bed(dsbSig11570,con='ccSeq_dsbsig_11687_reps.bed')