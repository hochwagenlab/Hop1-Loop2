ggplot2_blanktheme <- theme_classic() +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank(),
        title = NULL)
theme_set(ggplot2_blanktheme)



bedgraphs <- c("wildtype-Hop1-334-340-32-177-Reps-SK1Yue-PM_B3W3_MACS2_FE.bdg copy.gz",
               "11644-NHantiHop1-20220331-20220125-20221004-Reps-SK1Yue-PM-PE_B4_W3_MACS2_FE.bdg.gz",
               "11757-hop1NH-20220125-20220721-20221004-Reps-SK1Yue-PM-PE_B4_W3_MACS2_FE.bdg.gz",
               "11758-hop1NH-20220125-20220721-20221004-Reps-SK1Yue-PM-PE_B4_W3_MACS2_FE.bdg.gz",
               "rec8D-Hop1-61-91-Reps-SK1Yue-PM_B3W3_MACS2_FE.bdg.gz")
#bedgraphs <- c("wildtype-Hop1-334-340-32-177-Reps-SK1Yue-PM_B3W3_MACS2_FE.bdg copy.gz")

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




# function to make heatmap of signal at cluster and desert hotspots
for (i in 1:length(bedgraphs)) {
  bedgraph_file <- bedgraphs[i]
  name <- strsplit(bedgraph_file,split = "-")
  name <- sapply(name, "[[",1)
  Sig <- hwglabr2::import_bedGraph(bedgraph_file, local_copy=T)
  genAvg <- hwglabr2::average_chr_signal(Sig)$genome_avrg
  Sig$score <- Sig$score/genAvg


max_score=max(Sig$score)

genome_info <- hwglabr2::get_chr_coordinates(genome = "SK1Yue")
chrs <- list('I','VI','III','IX','VIII','V','XI','X','XIV','II','XIII','XVI','XII','VII','XV','IV')
#chrs <- list('I')
p=list()
for(i in chrs) {
  print(i)
  p[[i]] = plot_signal(Sig,i,1000)
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
    geom_line(position='identity', colour = "#0099CD") +
    ylab(i) + ylim(-0.65,max_score) + xlim(0,1531933)
    
  a <- a + geom_point(cen_mid, mapping=aes(cen_midpt,-0.15),
                      colour = 'black',shape=2, size=2)
  clusterregion = data.frame(clusters[seqnames(clusters)==paste0("chr",i)])
  a <- a + geom_segment(clusterregion, mapping=aes(x=start, xend=end, y=0, yend=0),colour = 'black') 
  a <- a + theme(legend.position = "none")
  plist[[i]] <- a
}
pdf(paste0(name[1],"All_chr.pdf"), height = 20, width = 10)
print(wrap_plots(plist,nrow = 16))
dev.off()
}