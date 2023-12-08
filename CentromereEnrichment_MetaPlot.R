#Quantify Signal in Centromere
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

ggplot2_theme <- theme_classic()
theme_set(ggplot2_theme)





bedgraphs <- c("loop2-Hop1_loop2_Reps_Norm.bed",
               "loop2pch2-Hop1_loop2pch2_Reps_Norm.bed",
               "pch2-Hop1_pch2_reps_norm.bed",
               "wildtype-Hop1-334-340-32-177-Reps-SK1Yue-PM_B3W3_MACS2_FE.bed")
#bedgraphs <- c(
#"wildtype-Hop1-334-340-32-177-Reps-SK1Yue-PM_B3W3_MACS2_FE.bed")

genome_info <- hwglabr2::get_chr_coordinates(genome = 'SK1Yue')
binsc <- toDataframe(genome_info)
Cen_flank = 10000
cenStartminus <- binsc$start - Cen_flank
cenEndplus <- binsc$end + Cen_flank
chr <- binsc$chr
Cens_bpFlank <- data.frame(first_column = chr, second_column= cenStartminus, third_column = cenEndplus)
Cens_bpFlank <- toGRanges(Cens_bpFlank)
#df <- toDataframe(Cens_bpFlank)
#write.table(df, file="Cen10KbFlankRegions.bed", quote=F, sep="\t", row.names=F, col.names=F)
axis = hwglabr2::get_Red1_summits("SK1Yue")
hits = findOverlaps(Cens_bpFlank,axis)
axis_cens <- axis[subjectHits(hits)]
rm(hits)
axis_offcen <- axis[axis %outside% axis_cens]


for (i in 1:length(bedgraphs)) {
  bedgraph_file <- bedgraphs[i]
  name <- strsplit(bedgraph_file,split = "-")
  name <- sapply(name, "[[",1)
  Sig <- rtracklayer::import.bed(bedgraph_file)
  #Sig <- hwglabr2::import_bedGraph(bedgraph_file, local_copy=T)
  genAvg <- hwglabr2::average_chr_signal(Sig)$genome_avrg
  Sig$score <- Sig$score/genAvg
  
  mat1 <- normalizeToMatrix(Sig, axis_cens, value_column = "score",
                            extend = 1000, mean_mode = "weighted", w = 5,empty_value=NA)
  
  mat1_avrg <- hwglabr2::signal_mean_and_ci(mat1,
                                            ci=0.95, rep_bootstrap=1000,
                                            na_rm=TRUE)
  
  mat2 <- normalizeToMatrix(Sig, axis_offcen, value_column = "score",
                            extend = 1000, mean_mode = "weighted", w = 5,empty_value=NA)
  
  mat2_avrg <- hwglabr2::signal_mean_and_ci(mat2,
                                            ci=0.95, rep_bootstrap=1000,
                                            na_rm=TRUE)
  
  mat1_avrg_df <- data.frame(Position=seq(-199, 200), sample = 'cluster',mat1_avrg)
  mat2_avrg_df <- data.frame(Position=seq(-199, 200), sample = 'desert',mat2_avrg)
  alldata = rbind(mat1_avrg_df,mat2_avrg_df)
  pdf(paste0(name[1],"sigOn_axisSites_onVoffCENS_BootMetaPlot.pdf"), height = 5, width = 5)
  print(ggplot(alldata, aes(Position, Mean, group = sample, colour=sample, fill=sample)) +
          labs(title = element_blank())+
          ylim(1,5)+
          geom_line() +
          geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha=0.3,colour=NA) +
          geom_vline(xintercept = 0, lty = 3) +
          scale_x_continuous(breaks = c(-199, 0, 200),labels = c("-1 kb", "Axis summit", "1 kb")) +
          scale_color_manual(values = c("#00A14B","#662D91"))+
          theme(legend.position = "none"))
  dev.off()
}







#not necessary to work with spike-in files because this analysis is an internal test.. regions are compared
#within a single chip, not between chips.. 
setwd('/Users/darmokandjalad/Documents/Hop1PHD-Loop2/RAnalysis')

#WtFile <- paste0('/Users/darmokandjalad/Documents/Hop1PHD-Loop2/RAnalysis/Hop1-wildtype-334-340-32-177-Reps-SK1Yue-PM_B3W3_MACS2_FE.bdg.gz')
#Wtsig <- import_bedGraph(WtFile)
#Wtsig <- sort(Wtsig)

Wtsig <- rtracklayer::import.bed('wildtype-Hop1-334-340-32-177-Reps-SK1Yue-PM_B3W3_MACS2_FE.bed')
Wtsig <- sort(Wtsig)

#loopFile <- paste0('/Users/darmokandjalad/Documents/Hop1PHD-Loop2/RAnalysis/11644-NHantiHop1-20220331-20220125-20221004-Reps-SK1Yue-PM-PE_B4_W3_MACS2_FE.bdg.gz')
#loopSig <- import_bedGraph(loopFile)
#loopSig <- sort(loopSig)

loopSig <- rtracklayer::import.bed('loop2-Hop1_loop2_Reps_Norm.bed')
loopSig <- sort(loopSig)

#loopPchFile <- paste0('/Users/darmokandjalad/Documents/Hop1PHD-Loop2/RAnalysis/11757-hop1NH-20220125-20220721-20221004-Reps-SK1Yue-PM-PE_B4_W3_MACS2_FE.bdg.gz')
#lpPchSig <- import_bedGraph(loopPchFile)
#lpPchSig <- sort(lpPchSig)

lpPchSig <- rtracklayer::import.bed('loop2pch2-Hop1_loop2pch2_Reps_Norm.bed')
lpPchSig <- sort(lpPchSig)

#pchFile <- paste0('/Users/darmokandjalad/Documents/Hop1PHD-Loop2/RAnalysis/11758-hop1NH-20220125-20220721-20221004-Reps-SK1Yue-PM-PE_B4_W3_MACS2_FE.bdg.gz')
#pchSig <- import_bedGraph(pchFile)
#pchSig <- sort(pchSig)

pchSig <- rtracklayer::import.bed('pch2-Hop1_pch2_reps_norm.bed')
pchSig <- sort(pchSig)

rec8File <- paste0('/Users/darmokandjalad/Documents/Hop1PHD-Loop2/RAnalysis/Hop1-rec8D-61-91-Reps-SK1Yue-PM_B3W3_MACS2_FE.bdg.gz')
rec8Sig <- import_bedGraph(rec8File)
rec8Sig <- sort(rec8Sig)




genome_info <- hwglabr2::get_chr_coordinates(genome = 'SK1Yue')
binsc <- toDataframe(genome_info)
Cen_flank = 10000
cenStartminus <- binsc$start - Cen_flank
cenEndplus <- binsc$end + Cen_flank
chr <- binsc$chr
Cens_bpFlank <- data.frame(first_column = chr, second_column= cenStartminus, third_column = cenEndplus)
Cens_bpFlank <- toGRanges(Cens_bpFlank)
#df <- toDataframe(Cens_bpFlank)
#write.table(df, file="Cen10KbFlankRegions.bed", quote=F, sep="\t", row.names=F, col.names=F)


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
                #Avg on cen 1.670775
                A <- mean(mat1_avrg_df$Mean)
                mat2_avrg_df <- data.frame(mat2_avrg)
                mat2_avrg_df['Position'] = row
                mat2_avrg_df['sample'] = 'Off_cen'
                #Avg off cen 1.63729
                B <- mean(mat2_avrg_df$Mean)
                #A on cen, B off cen A/B is fold diff on cen 
                Wtfold <- A/B
                #1.020452
                alldata = rbind(mat1_avrg_df,mat2_avrg_df)
                pdf("Wt_Norm_sigOn_axisSites_BootMetaPlot.pdf", height = 5, width = 5)
                print(alldata, aes(Position, Mean, group = sample, colour=sample, fill=sample)) +
                  ylim(1,5)+
                  geom_vline(xintercept = 0, lty = 3) +
                  scale_x_continuous(breaks = c(-199, 0, 200),
                                     labels = c("-10 kb", "summit", "10 kb"))+
                  scale_color_manual(values = c("#662D91", "#00A14B"))+
                  theme(legend.position="none")+
                  geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha=0.3,colour=NA)+
                  geom_line()
                dev.off()
                
                p <- ggplot(alldata, aes(Position, Mean, group = sample, colour=sample, fill=sample)) +
                  ylim(1,5)+
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
                  ylim(1,5)+
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
                  ylim(1,5)+
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
                  ylim(1,5)+
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
                  ylim(1,5)+
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




































#Heat Maps


# prepare centromere data for Loop2
                
                setwd('/Users/darmokandjalad/Documents/Hop1PHD-Loop2/RAnalysis/Centromere_Enrichment')
                centromeres <- hwglabr2::get_chr_coordinates(genome='SK1Yue', as_df=FALSE)
                
                midpoint <- floor(width(centromeres) / 2)
                start(centromeres) <- start(centromeres) + midpoint
                end(centromeres) <- start(centromeres)
#for Loop2                
                bedgraph_file <- "/Users/darmokandjalad/Documents/Hop1PHD-Loop2/RAnalysis/11644-NHantiHop1-20220331-20220125-20221004-Reps-SK1Yue-PM-PE_B4_W3_MACS2_FE.bdg.gz"
                naming <- "hop1_loop2"
                
                Hop1_bg <- hwglabr2::import_bedGraph(bedgraph_file)
                genAvg <- hwglabr2::average_chr_signal(Hop1_bg)$genome_avrg
                Hop1_bg$score <- Hop1_bg$score/genAvg
                extend <- 10000
                windowL <- 10
                mat1 <- normalizeToMatrix(Hop1_bg, centromeres, value_column = "score",
                                          extend = extend, mean_mode = "w0", w = windowL)
                col_fun1 <- colorRamp2(c(0,quantile(mat1, c( 0.5, 0.99))), c("#662D91", "white", "#00A14B"))
                col_fun <- col_fun1
                col_fun2 <- c("#CBDB2A")
                pdf("heatmap_hop1-loop2.pdf", width = 4, height = 8)
                print(EnrichedHeatmap(mat1, col = col_fun, name = "Signal", show_row_names = TRUE,
                                      top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:6),
                                                                                             show_error = TRUE)),
                                      
                                      row_order = 1:length(centromeres),
                                     axis_name = c(paste0("-",extend/1000,"kb"), "Centromeres", paste0(extend/1000,"kb")),
                                      column_title = naming))
               
                dev.off()

#for Wildtype
     
                bedgraph_file <- "/Users/darmokandjalad/Documents/Hop1PHD-Loop2/RAnalysis/Hop1-wildtype-334-340-32-177-Reps-SK1Yue-PM_B3W3_MACS2_FE.bdg.gz"
                naming <- "wildtype"
                
                Hop1_bg <- hwglabr2::import_bedGraph(bedgraph_file)
                genAvg <- hwglabr2::average_chr_signal(Hop1_bg)$genome_avrg
                Hop1_bg$score <- Hop1_bg$score/genAvg
                extend <- 10000
                windowL <- 10
                mat1 <- normalizeToMatrix(Hop1_bg, centromeres, value_column = "score",
                                          extend = extend, mean_mode = "w0", w = windowL)
                col_fun1 <- colorRamp2(c(0,quantile(mat1, c( 0.5, 0.99))), c("#662D91", "white", "#00A14B"))
                col_fun <- col_fun1
                pdf("heatmap_wildtype.pdf", width = 4, height = 8)
                print(EnrichedHeatmap(mat1, col = col_fun, name = "Signal",
                                      top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:6),
                                                                                               show_error = TRUE)),
                                      
                                      row_order = 1:length(centromeres),
                                      axis_name = c(paste0("-",extend/1000,"kb"), "Centromeres", paste0(extend/1000,"kb")),
                                      column_title = naming))
                dev.off()

                
                
#for loop2pch2
                
                bedgraph_file <- "/Users/darmokandjalad/Documents/Hop1PHD-Loop2/RAnalysis/11757-hop1NH-20220125-20220721-20221004-Reps-SK1Yue-PM-PE_B4_W3_MACS2_FE.bdg.gz"
                naming <- "loop2 pch2"
                
                Hop1_bg <- hwglabr2::import_bedGraph(bedgraph_file)
                genAvg <- hwglabr2::average_chr_signal(Hop1_bg)$genome_avrg
                Hop1_bg$score <- Hop1_bg$score/genAvg
                extend <- 10000
                windowL <- 10
                mat1 <- normalizeToMatrix(Hop1_bg, centromeres, value_column = "score",
                                          extend = extend, mean_mode = "w0", w = windowL)
                col_fun1 <- colorRamp2(c(0,quantile(mat1, c( 0.5, 0.99))), c("#662D91", "white", "#00A14B"))
                col_fun <- col_fun1
                pdf("heatmap_loop2pch2.pdf", width = 4, height = 8)
                print(EnrichedHeatmap(mat1, col = col_fun, name = "Signal",
                                      top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:6),
                                                                                               show_error = TRUE)),
                                      
                                      row_order = 1:length(centromeres),
                                      axis_name = c(paste0("-",extend/1000,"kb"), "Centromeres", paste0(extend/1000,"kb")),
                                      column_title = naming))
                dev.off()                

#for pch2
                
                bedgraph_file <- "/Users/darmokandjalad/Documents/Hop1PHD-Loop2/RAnalysis/11758-hop1NH-20220125-20220721-20221004-Reps-SK1Yue-PM-PE_B4_W3_MACS2_FE.bdg.gz"
                naming <- "pch2"
                
                Hop1_bg <- hwglabr2::import_bedGraph(bedgraph_file)
                genAvg <- hwglabr2::average_chr_signal(Hop1_bg)$genome_avrg
                Hop1_bg$score <- Hop1_bg$score/genAvg
                extend <- 10000
                windowL <- 10
                mat1 <- normalizeToMatrix(Hop1_bg, centromeres, value_column = "score",
                                          extend = extend, mean_mode = "w0", w = windowL)
                col_fun1 <- colorRamp2(c(0,quantile(mat1, c( 0.5, 0.99))), c("#662D91", "white", "#00A14B"))
                col_fun <- col_fun1
                pdf("heatmap_pch2.pdf", width = 4, height = 8)
                print(EnrichedHeatmap(mat1, col = col_fun, name = "Signal",
                                      top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:6),
                                                                                               show_error = TRUE)),
                                      
                                      row_order = 1:length(centromeres),
                                      axis_name = c(paste0("-",extend/1000,"kb"), "Centromeres", paste0(extend/1000,"kb")),
                                      column_title = naming))
                dev.off()  


#for hop1sig in rec8D
                
                bedgraph_file <- "/Users/darmokandjalad/Documents/Hop1PHD-Loop2/RAnalysis/Hop1-rec8D-61-91-Reps-SK1Yue-PM_B3W3_MACS2_FE.bdg.gz"
                naming <- "rec8"
                
                Hop1_bg <- hwglabr2::import_bedGraph(bedgraph_file)
                genAvg <- hwglabr2::average_chr_signal(Hop1_bg)$genome_avrg
                Hop1_bg$score <- Hop1_bg$score/genAvg
                extend <- 10000
                windowL <- 10
                mat1 <- normalizeToMatrix(Hop1_bg, centromeres, value_column = "score",
                                          extend = extend, mean_mode = "w0", w = windowL)
                col_fun1 <- colorRamp2(c(0,quantile(mat1, c( 0.5, 0.99))), c("#00A14B", "white", "#B72467"))
                col_fun <- col_fun1
                pdf("heatmap_Hop1sig-rec8D.pdf")
                print(EnrichedHeatmap(mat1, col = col_fun, name = "Signal",
                                      top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:6),
                                                                                               show_error = TRUE)),
                                      
                                      row_order = 1:length(centromeres),
                                      axis_name = c(paste0("-",extend/1000,"kb"), "Centromeres", paste0(extend/1000,"kb")),
                                      column_title = naming))
                dev.off()  












