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
library(magick)
library(EnrichedHeatmap)
library(circlize)
library(GenomicRanges)


#----------------------------------------------------------------#
# Figure 1A (for all chromosomes)                                #
#----------------------------------------------------------------#
# ggplot theme
ggplot2_blanktheme <- theme_classic() +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_text(angle = 0),
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())
theme_set(ggplot2_blanktheme)



WtFile <- paste0('/Users/darmokandjalad/Documents/Hop1PHD-Loop2/RAnalysis/Hop1-wildtype-334-340-32-177-Reps-SK1Yue-PM_B3W3_MACS2_FE.bdg.gz')
Wtsig <- import_bedGraph(WtFile)
Wtsig <- sort(Wtsig)
genAvg <- hwglabr2::average_chr_signal(Wtsig)$genome_avrg
Wtsig$score <- Wtsig$score/genAvg


loopFile <- paste0('/Users/darmokandjalad/Documents/Hop1PHD-Loop2/RAnalysis/11644-NHantiHop1-20220331-20220125-20221004-Reps-SK1Yue-PM-PE_B4_W3_MACS2_FE.bdg.gz')
loopSig <- import_bedGraph(loopFile)
loopSig <- sort(loopSig)
genAvg <- hwglabr2::average_chr_signal(loopSig)$genome_avrg
loopSig$score <- loopSig$score/genAvg

loopPchFile <- paste0('/Users/darmokandjalad/Documents/Hop1PHD-Loop2/RAnalysis/11757-hop1NH-20220125-20220721-20221004-Reps-SK1Yue-PM-PE_B4_W3_MACS2_FE.bdg.gz')
lpPchSig <- import_bedGraph(loopPchFile)
lpPchSig <- sort(lpPchSig)
genAvg <- hwglabr2::average_chr_signal(lpPchSig)$genome_avrg
lpPchSig$score <- lpPchSig$score/genAvg

pchFile <- paste0('/Users/darmokandjalad/Documents/Hop1PHD-Loop2/RAnalysis/11758-hop1NH-20220125-20220721-20221004-Reps-SK1Yue-PM-PE_B4_W3_MACS2_FE.bdg.gz')
pchSig <- import_bedGraph(pchFile)
pchSig <- sort(pchSig)
genAvg <- hwglabr2::average_chr_signal(pchSig)$genome_avrg
pchSig$score <- pchSig$score/genAvg

rec8File <- paste0('/Users/darmokandjalad/Documents/Hop1PHD-Loop2/RAnalysis/Hop1-rec8D-61-91-Reps-SK1Yue-PM_B3W3_MACS2_FE.bdg.gz')
rec8Sig <- import_bedGraph(rec8File)
rec8Sig <- sort(rec8Sig)
genAvg <- hwglabr2::average_chr_signal(rec8Sig)$genome_avrg
rec8Sig$score <- rec8Sig$score/genAvg


#Wtsig$group <- "WT"
#loopSig$group <- "loop2"
#all <- c(Wtsig, loopSig)
all <-lpPchSig

Flank = 20000

genome_info <- hwglabr2::get_chr_coordinates(genome = 'SK1Yue')
genome_infodf <- toDataframe(genome_info)
cen_mid <- genome_info
cen_mid <- round((start(cen_mid) + end(cen_mid))/2)
cen_plus20 <- cen_mid + Flank
cen_min20 <- cen_mid - Flank
chr <- genome_infodf$chr
flank_cen <- data.frame(first_column = genome_infodf$chr, second_column = cen_min20, third_column = cen_plus20)
colnames(flank_cen) <- c("seqnames","start","end")
flank_cen <- toGRanges(flank_cen)
cen_mid <- data.frame(cen_midpt = cen_mid, y = 0)

hits = findOverlaps(flank_cen,all)
all_cens <- all[subjectHits(hits)]
rm(hits)

sample <- all_cens
tile = 10
positions_ex <- seq(from = 1 , to = 40001 , by = 10)
# function for plotting signal across chromosomes
plot_signal <- function(sample,chrnum,tile) {
  genome_info <- hwglabr2::get_chr_coordinates(genome = 'SK1Yue')
  sample_ex <- sort(GenomeInfoDb::sortSeqlevels(sample))
  GenomeInfoDb::seqlengths(sample_ex) <- GenomeInfoDb::seqlengths(genome_info)
  bins_ex <- GenomicRanges::tileGenome(GenomeInfoDb::seqlengths(sample_ex),
                                       tilewidth=tile,
                                       cut.last.tile.in.chrom=TRUE)
  score_ex <- GenomicRanges::coverage(sample_ex, weight="score")
  bins_ex <- GenomicRanges::binnedAverage(bins_ex, score_ex, "binned_score")
  #bins_ex <- GenomeInfoDb::keepSeqlevels(bins_ex, "chrII", pruning.mode="coarse")
  
  bins_ex <- GenomeInfoDb::keepSeqlevels(bins_ex, paste0("chr",chrnum),pruning.mode="coarse")
  hits = findOverlaps(flank_cen,bins_ex)
  bins_ex <- bins_ex[subjectHits(hits)]
  rm(hits)


  positions_ex <- seq(from = 1 , to = 40001 , by = 10)
  df_ex <- data.frame(seqnames=paste0('chr',chrnum),position=positions_ex, signal=bins_ex$binned_score)
  #df_ex <- data.frame(seqnames="chrII",position=positions_ex, signal=bins_ex$binned_score)
  
  
  #df_exS <- df_ex$position - (df_ex$position/2)
  #df_exE <- df_ex$position + (df_ex$position/2)
  #df_ex$start <- df_exS
  #df_ex$end <- df_exE
  #df_ex$width <- tile
  
  #df_ex <- makeGRangesFromDataFrame(df_ex,
                                    #start.field="start",
                                    #end.field=c("end"))
  #en_ex <- GenomeInfoDb::keepSeqlevels(genome_info, "chrII", pruning.mode="coarse")
  #cen_mid_ex <- cen_ex
  #cen_mid_ex <- round((start(cen_mid_ex) + end(cen_mid_ex))/2)
  #cen_ex <- toDataframe(cen_ex)
  #cen_ex$start <- cen_ex$start - 20000
  #cen_ex$end <- cenr_ex$end + 20000
  #cen_ex <- makeGRangesFromDataFrame(cen_ex)
  #hits = findOverlaps(cen_ex,df_ex)
  #df_ex <- df_ex[subjectHits(hits)]
  #rm(hits)
}


max_score=max(sample$signal)
#genome_info <- hwglabr2::get_chr_coordinates(genome = 'SK1Yue')
#chrs <- list('I','VI','III','IX','VIII','V','XI','X','XIV','II','XIII','XVI','XII','VII','XV','IV')
chrs <- list('I','II','III','IV','V','VI','VII','VIII','IX','X','XI','XII','XIII','XIV','XV','XVI')

#chrs <- list('I','VI', 'III')
p=list()
for(i in chrs) {
  print(i)
  p[[i]] = plot_signal(sample,i,tile)
}

#clusters <- rtracklayer::import.bed('/Users/darmokandjalad/Documents/Hop1PHD-Loop2/RAnalysis/clusters_joined.bed')

plist=list()
for(i in chrs){
  a = data.frame()
  temp=data.frame()
  temp = p[[i]]
  cen_mid = 10001
  a <- ggplot(temp,aes(position,signal)) +
    geom_line() +
    ylab(i) + ylim(-0.65,11) + xlim(0,40001)
  #ylim is 10 for all except pch2, where ylim is 12
  a <- a + geom_vline(xintercept = 20000, alpha = 0.2)
  #a <- a + geom_point(cen_mid, mapping=aes(cen_midpt,-0.15),
                     # size = 1.6, colour = 'green',shape=20)
  #clusterregion = data.frame(clusters[seqnames(clusters)==paste0("chr",i)])
  #a <- a + geom_segment(clusterregion, size = 6,alpha = 0.6,
                        #mapping=aes(x = start, y = max_score/2, xend = end, yend = max_score/2, colour = "segment"))
  plist[[i]] <- a
}

#wrap_plots(plist,nrow = 16)

p <- wrap_plots(plist,nrow = 16)
p

#export pdf, portrait 4X8inches

####################################

#Put signal over genome average
#pick signal to analyzie
ggplot2_blanktheme <- theme_classic() +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank(),
        title = NULL)
theme_set(ggplot2_blanktheme)

            #For Wildtype
            signal = Wtsig
            Cen_flank = 10000
            
            
            genAvg <- hwglabr2::average_chr_signal(signal)$genome_avrg
            signal$score <- signal$score/genAvg
            
            #get genome
            genome_info <- hwglabr2::get_chr_coordinates(genome = 'SK1Yue')
            
            # Sort sequences and levels to make sure they match
            signal <- sort(GenomeInfoDb::sortSeqlevels(signal))
            genome_info <- GenomeInfoDb::sortSeqlevels(genome_info)
            
            # isolate centromere regions
            #how far out from cen
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
              
              geom_vline(xintercept = 0, lty = 3) +
              scale_x_continuous(breaks = c(-199, 0, 200),
                                 labels = c("-10 kb", "summit", "10 kb"))+
              scale_color_manual(values = c("#662D91", "#00A14B"))+
              theme(legend.position="none")
            
            # Add confidence interval as a ribbon
            p <- p + geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha=0.3,colour=NA) 
            p <- p + geom_line()
            p





#for loop2
            signal = loopSig
            Cen_flank = 10000
            
            
            genAvg <- hwglabr2::average_chr_signal(signal)$genome_avrg
            signal$score <- signal$score/genAvg
            
            #get genome
            genome_info <- hwglabr2::get_chr_coordinates(genome = 'SK1Yue')
            
            # Sort sequences and levels to make sure they match
            signal <- sort(GenomeInfoDb::sortSeqlevels(signal))
            genome_info <- GenomeInfoDb::sortSeqlevels(genome_info)
            
            # isolate centromere regions
            #how far out from cen
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
              labs(x = "Distance to axis site (bp)", y = "Average\nChIP-seq signal") +
              geom_vline(xintercept = 0, lty = 3) +
              scale_x_continuous(breaks = c(-199, 0, 200),
                                 labels = c("-10 kb", "summit", "10 kb"))+
              scale_color_manual(values = c("#662D91", "#00A14B"))+
              theme(legend.position="none")
            
            # Add confidence interval as a ribbon
            p <- p + geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha=0.3,colour=NA) 
            p <- p + geom_line()
            p
            
            


#for loop2Pch2
signal = lpPchSig
Cen_flank = 10000


genAvg <- hwglabr2::average_chr_signal(signal)$genome_avrg
signal$score <- signal$score/genAvg

#get genome
genome_info <- hwglabr2::get_chr_coordinates(genome = 'SK1Yue')

# Sort sequences and levels to make sure they match
signal <- sort(GenomeInfoDb::sortSeqlevels(signal))
genome_info <- GenomeInfoDb::sortSeqlevels(genome_info)

# isolate centromere regions
#how far out from cen
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
  labs(x = "Distance to axis site (bp)", y = "Average\nChIP-seq signal") +
  geom_vline(xintercept = 0, lty = 3) +
  scale_x_continuous(breaks = c(-199, 0, 200),
                     labels = c("-10 kb", "summit", "10 kb"))+
  scale_color_manual(values = c("#662D91", "#00A14B"))+
  theme(legend.position="none")

# Add confidence interval as a ribbon
p <- p + geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha=0.3,colour=NA) 
p <- p + geom_line()
p


###For Pch2

signal = pchSig
Cen_flank = 10000


genAvg <- hwglabr2::average_chr_signal(signal)$genome_avrg
signal$score <- signal$score/genAvg

#get genome
genome_info <- hwglabr2::get_chr_coordinates(genome = 'SK1Yue')

# Sort sequences and levels to make sure they match
signal <- sort(GenomeInfoDb::sortSeqlevels(signal))
genome_info <- GenomeInfoDb::sortSeqlevels(genome_info)

# isolate centromere regions
#how far out from cen
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
  labs(x = "Distance to axis site (bp)", y = "Average\nChIP-seq signal") +
  geom_vline(xintercept = 0, lty = 3) +
  scale_x_continuous(breaks = c(-199, 0, 200),
                     labels = c("-10 kb", "summit", "10 kb"))+
  scale_color_manual(values = c("#662D91", "#00A14B"))+
  theme(legend.position="none")

# Add confidence interval as a ribbon
p <- p + geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha=0.3,colour=NA) 
p <- p + geom_line()
p


###for rec8

signal = rec8Sig
Cen_flank = 10000


genAvg <- hwglabr2::average_chr_signal(signal)$genome_avrg
signal$score <- signal$score/genAvg

#get genome
genome_info <- hwglabr2::get_chr_coordinates(genome = 'SK1Yue')

# Sort sequences and levels to make sure they match
signal <- sort(GenomeInfoDb::sortSeqlevels(signal))
genome_info <- GenomeInfoDb::sortSeqlevels(genome_info)

# isolate centromere regions
#how far out from cen
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
  labs(x = "Distance to axis site (bp)", y = "Average\nChIP-seq signal") +
  geom_vline(xintercept = 0, lty = 3) +
  scale_x_continuous(breaks = c(-199, 0, 200),
                     labels = c("-10 kb", "summit", "10 kb"))+
  scale_color_manual(values = c("#662D91", "#00A14B"))+
  theme(legend.position="none")

# Add confidence interval as a ribbon
p <- p + geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha=0.3,colour=NA) 
p <- p + geom_line()
p  


###################################################       