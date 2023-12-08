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
library(ggpubr)

ggplot2_blanktheme <- theme_classic() +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank(),
        title = NULL)
theme_set(ggplot2_blanktheme)
ggplot2_nearblankTheme <- theme_classic() +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank(),
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank(),
        title = NULL)
theme_set(ggplot2_nearblankTheme)


lpN = rtracklayer::import.bed('/Users/darmokandjalad/Documents/Hop1PHD-Loop2/RAnalysis/SpikeValue/Hop1_loop2_Reps_Norm.bed')
lpPchN = rtracklayer::import.bed('/Users/darmokandjalad/Documents/Hop1PHD-Loop2/RAnalysis/SpikeValue/Hop1_loop2pch2_Reps_Norm.bed')
pchN = rtracklayer::import.bed('/Users/darmokandjalad/Documents/Hop1PHD-Loop2/RAnalysis/SpikeValue/Hop1_pch2_reps_norm.bed')
Wtsig <- hwglabr2::import_bedGraph('/Users/darmokandjalad/Documents/Hop1PHD-Loop2/RAnalysis/Hop1-wildtype-334-340-32-177-Reps-SK1Yue-PM_B3W3_MACS2_FE.bdg.gz')
genome_info <- hwglabr2::get_chr_coordinates(genome = 'SK1Yue')


signal <- Wtsig
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
        Wtsig <- bins
        
signal <- lpN      
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
          lpN <- bins
          
LeftOver <-  Wtsig$binned_score - lpN$binned_score   
Wtsig$score <- LeftOver
mcols(Wtsig)$binned_score <- NULL
signal <- Wtsig

        
              plot_signal <- function(sample,chrnum,tile) {
                genome_info <- hwglabr2::get_chr_coordinates(genome = 'SK1Yue')
                sample_ex <- sort(GenomeInfoDb::sortSeqlevels(sample))
                GenomeInfoDb::seqlengths(sample_ex) <- GenomeInfoDb::seqlengths(genome_info)
                bins_ex <- GenomicRanges::tileGenome(GenomeInfoDb::seqlengths(sample_ex),
                                                     tilewidth=tile,
                                                     cut.last.tile.in.chrom=TRUE)
                score_ex <- GenomicRanges::coverage(sample_ex, weight="score")
                bins_ex <- GenomicRanges::binnedAverage(bins_ex, score_ex, "binned_score")
                bins_ex <- GenomeInfoDb::keepSeqlevels(bins_ex, paste0("chr",chrnum),pruning.mode="coarse")
                positions_ex <- bins_ex@ranges@start + floor(bins_ex@ranges@width / 2)
                df_ex <- data.frame(seqnames=paste0('chr',chrnum),position=positions_ex, signal=bins_ex$binned_score)
              }
              
              
              setwd('/Users/darmokandjalad/Documents/Hop1PHD-Loop2/RAnalysis')
              max_score=max(signal$score)
              min_score=min(signal$score)
              genome_info <- hwglabr2::get_chr_coordinates(genome = 'SK1Yue')
              chrs <- list('I','VI','III','IX','VIII','V','XI','X','XIV','II','XIII','XVI','XII','VII','XV','IV')
              p=list()
              for(i in chrs) {
                print(i)
                p[[i]] = plot_signal(signal,i,1000)
              }
              
              
              clusters = rtracklayer::import.bed('clusters_joined.bed')
              
              plist=list()
              for(i in chrs){
                a = data.frame()
                temp=data.frame()
                temp = p[[i]]
                if (exists("cen_mid")) {
                  rm(cen_mid) }
                cen_mid = genome_info[seqnames(genome_info)==paste0("chr",i)]
                cen_mid <- round((start(cen_mid) + end(cen_mid))/2)
                cen_mid <- data.frame(cen_midpt = cen_mid, y = 0)
                a <- ggplot(temp,aes(position,signal)) +
                  geom_line(position='identity', size=0.25, colour = "#B72467") +
                  ylab(i) + ylim(min_score,max_score) + xlim(0,1531933) 
                a <- a + geom_hline(yintercept = 0, linetype="dashed", alpha =0.5, size =0.25)
                
                a <- a + geom_point(cen_mid, mapping=aes(cen_midpt,-0.15),
                                    colour = 'black',shape=2, size=2)
                clusterregion = data.frame(clusters[seqnames(clusters)==paste0("chr",i)])
                a <- a + geom_segment(clusterregion, mapping=aes(x=start, xend=end, y=0, yend=0),colour = '#0099CD', alpha=0.5) 
                a <- a + theme(legend.position = "none")
                plist[[i]] <- a
              }
              
              wrap_plots(plist,nrow = 16)
              
              pdf("Wt-minus-Loop2_all-.pdf", height = 16, width = 5)
              print(wrap_plots(plist,nrow = 16))
              dev.off()

              
Wtsig <- hwglabr2::import_bedGraph('/Users/darmokandjalad/Documents/Hop1PHD-Loop2/RAnalysis/Hop1-wildtype-334-340-32-177-Reps-SK1Yue-PM_B3W3_MACS2_FE.bdg.gz')
signal <- Wtsig
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
              Wtsig <- bins
              
signal <- lpPchN      
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
              lpN <- bins
              
              LeftOver <-  Wtsig$binned_score - lpN$binned_score   
              Wtsig$score <- LeftOver
              mcols(Wtsig)$binned_score <- NULL
              signal <- Wtsig              
              
              plot_signal <- function(sample,chrnum,tile) {
                genome_info <- hwglabr2::get_chr_coordinates(genome = 'SK1Yue')
                sample_ex <- sort(GenomeInfoDb::sortSeqlevels(sample))
                GenomeInfoDb::seqlengths(sample_ex) <- GenomeInfoDb::seqlengths(genome_info)
                bins_ex <- GenomicRanges::tileGenome(GenomeInfoDb::seqlengths(sample_ex),
                                                     tilewidth=tile,
                                                     cut.last.tile.in.chrom=TRUE)
                score_ex <- GenomicRanges::coverage(sample_ex, weight="score")
                bins_ex <- GenomicRanges::binnedAverage(bins_ex, score_ex, "binned_score")
                bins_ex <- GenomeInfoDb::keepSeqlevels(bins_ex, paste0("chr",chrnum),pruning.mode="coarse")
                positions_ex <- bins_ex@ranges@start + floor(bins_ex@ranges@width / 2)
                df_ex <- data.frame(seqnames=paste0('chr',chrnum),position=positions_ex, signal=bins_ex$binned_score)
              }
              
              
              setwd('/Users/darmokandjalad/Documents/Hop1PHD-Loop2/RAnalysis')
              max_score=max(signal$score)
              min_score=min(signal$score)
              genome_info <- hwglabr2::get_chr_coordinates(genome = 'SK1Yue')
              chrs <- list('I','VI','III','IX','VIII','V','XI','X','XIV','II','XIII','XVI','XII','VII','XV','IV')
              p=list()
              for(i in chrs) {
                print(i)
                p[[i]] = plot_signal(signal,i,1000)
              }
              
              
              clusters = rtracklayer::import.bed('clusters_joined.bed')
              
              plist=list()
              for(i in chrs){
                a = data.frame()
                temp=data.frame()
                temp = p[[i]]
                if (exists("cen_mid")) {
                  rm(cen_mid) }
                cen_mid = genome_info[seqnames(genome_info)==paste0("chr",i)]
                cen_mid <- round((start(cen_mid) + end(cen_mid))/2)
                cen_mid <- data.frame(cen_midpt = cen_mid, y = 0)
                a <- ggplot(temp,aes(position,signal)) +
                  geom_line(position='identity', size=0.25, colour = "#B72467") +
                  ylab(i) + ylim(min_score,max_score) + xlim(0,1531933) 
                a <- a + geom_hline(yintercept = 0, linetype="dashed", alpha =0.5, size =0.25)
                
                a <- a + geom_point(cen_mid, mapping=aes(cen_midpt,-0.15),
                                    colour = 'black',shape=2, size=2)
                clusterregion = data.frame(clusters[seqnames(clusters)==paste0("chr",i)])
                a <- a + geom_segment(clusterregion, mapping=aes(x=start, xend=end, y=0, yend=0),colour = '#0099CD', alpha=0.5) 
                a <- a + theme(legend.position = "none")
                plist[[i]] <- a
              }
              
              wrap_plots(plist,nrow = 16)
              
              pdf("Wt-minus-Loop2Pch2_all--.pdf", height = 16, width = 5)
              print(wrap_plots(plist,nrow = 16))
              dev.off()
              


