#Qunatify enrichment of signal in island or non-islands* double check this. Signal compared after normalizing each to genome average
################################################################################################
setwd('/Users/darmokandjalad/Documents/Hop1PHD-Loop2/RAnalysis')
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




bedgraphs <- c("wildtype-Hop1-334-340-32-177-Reps-SK1Yue-PM_B3W3_MACS2_FE.bdg copy.gz",
               "11644-NHantiHop1-20220331-20220125-20221004-Reps-SK1Yue-PM-PE_B4_W3_MACS2_FE.bdg.gz",
               "11757-hop1NH-20220125-20220721-20221004-Reps-SK1Yue-PM-PE_B4_W3_MACS2_FE.bdg.gz",
               "11758-hop1NH-20220125-20220721-20221004-Reps-SK1Yue-PM-PE_B4_W3_MACS2_FE.bdg.gz",
               "rec8D-Hop1-61-91-Reps-SK1Yue-PM_B3W3_MACS2_FE.bdg.gz")
bedgraphs <- c("loop2-Hop1_loop2_Reps_Norm.bed",
               "loop2pch2-Hop1_loop2pch2_Reps_Norm.bed",
               "pch2-Hop1_pch2_reps_norm.bed",
               "wildtype-Hop1-334-340-32-177-Reps-SK1Yue-PM_B3W3_MACS2_FE.bed")
#bedgraphs <- c(
               "wildtype-Hop1-334-340-32-177-Reps-SK1Yue-PM_B3W3_MACS2_FE.bed")

ggplot2_blanktheme <- theme_classic() +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank(),
        title = NULL)
theme_set(ggplot2_blanktheme)

ggplot2_theme <- theme_classic()
theme_set(ggplot2_theme)

# function to make heatmap of signal at cluster and desert hotspots
for (i in 1:length(bedgraphs)) {
  bedgraph_file <- bedgraphs[i]
  name <- strsplit(bedgraph_file,split = "-")
  name <- sapply(name, "[[",1)
  Sig <- rtracklayer::import.bed(bedgraph_file)
  #Sig <- hwglabr2::import_bedGraph(bedgraph_file, local_copy=T)
  genAvg <- hwglabr2::average_chr_signal(Sig)$genome_avrg
  Sig$score <- Sig$score/genAvg

  mat1 <- normalizeToMatrix(Sig, axis_cluster, value_column = "score",
                            extend = 1000, mean_mode = "weighted", w = 5,empty_value=NA)
  
  mat1_avrg <- hwglabr2::signal_mean_and_ci(mat1,
                                            ci=0.95, rep_bootstrap=1000,
                                            na_rm=TRUE)
  
  mat2 <- normalizeToMatrix(Sig, axis_desert, value_column = "score",
                            extend = 1000, mean_mode = "weighted", w = 5,empty_value=NA)
  
  mat2_avrg <- hwglabr2::signal_mean_and_ci(mat2,
                                            ci=0.95, rep_bootstrap=1000,
                                            na_rm=TRUE)
  
  mat1_avrg_df <- data.frame(Position=seq(-199, 200), sample = 'cluster',mat1_avrg)
  mat2_avrg_df <- data.frame(Position=seq(-199, 200), sample = 'desert',mat2_avrg)
  alldata = rbind(mat1_avrg_df,mat2_avrg_df)
  pdf(paste0(name[1],"sigOn_axisSites_BootMetaPlot.pdf"), height = 5, width = 5)
  print(ggplot(alldata, aes(Position, Mean, group = sample, colour=sample, fill=sample)) +
          labs(title = element_blank())+
          ylim(1,3.5)+
        geom_line() +
        geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha=0.3,colour=NA) +
        geom_vline(xintercept = 0, lty = 3) +
        scale_x_continuous(breaks = c(-199, 0, 200),labels = c("-1 kb", "Axis summit", "1 kb")) +
        scale_color_manual(values = c("#00A14B","#662D91"))+
        theme(legend.position = "none"))
  dev.off()
}


###not used--- use total
#WtClustAvg <- mean(mat1_avrg_df$Mean)
#1.967167
#WtDestAvg <- mean(mat2_avrg_df$Mean)
#1.541556
#WtFold <- WtClustAvg/WtDestAvg
#1.276092

#loop2ClustAvg <- mean(mat1_avrg_df$Mean)
#1.560487
#loop2DestAvg <- mean(mat2_avrg_df$Mean)
#1.748088
#loop2Fold <- loop2ClustAvg/loop2DestAvg
#0.8926819