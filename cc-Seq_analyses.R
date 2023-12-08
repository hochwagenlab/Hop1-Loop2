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

theme_set(theme_classic())

setwd('/Users/darmokandjalad/Documents/Hop1PHD-Loop2/RAnalysis')

rad50 <- rtracklayer::import.bed('ccSeq_dsbsig_11570_reps.bed')
lprad50 <- rtracklayer::import.bed('ccSeq_dsbsig_11687_reps.bed')
genome_info <- hwglabr2::get_chr_coordinates(genome = 'SK1Yue')
chrSize <- read.csv('/Users/darmokandjalad/Documents/Hop1PHD-Loop2/RAnalysis/chrSize.csv')

signal <- lprad50
genAvg <- hwglabr2::average_chr_signal(signal)$genome_avrg
signal$score <- signal$score/genAvg

signal <- sort(GenomeInfoDb::sortSeqlevels(signal))
genome_info <- GenomeInfoDb::sortSeqlevels(genome_info)
genome_info <- dropSeqlevels(genome_info, "chrIII", pruning.mode = "coarse")

# Add info to signal object
GenomeInfoDb::seqlengths(signal) <- GenomeInfoDb::seqlengths(genome_info)
signal <- sort(GenomeInfoDb::sortSeqlevels(signal))

#calculate the signal per chromsomes
chrSig <- hwglabr2::average_chr_signal(signal)
chrSig_df <- toDataframe(chrSig$seq_avrg)

chrSize %>% filter(chr=='chr3')
chrSize <- chrSize[-3,]

chrLength <- chrSize$length

chrSig_df <- cbind(chrSig_df, chrLength)

colnames(chrSig_df) <- c('chr', 'wt_avg_sig', 'chrlength')


p <- ggplot(chrSig_df, aes(x=chrLength, y = wt_avg_sig))
p <- p + geom_point()
p


Chr_sig_log2 <- read.csv('/Users/darmokandjalad/Documents/Hop1PHD-Loop2/RAnalysis/chrLeng_Sig_log2_AND_dsb.csv')

p <- ggplot(data = Chr_sig_log2, aes(x = length, y = sig, group = group, col = group)) + 
        geom_point(size = 3) +
        geom_smooth(method = "lm", se = FALSE) +
        scale_color_manual(values = c("#0099CD","#B72467","#F68B1F","#CBDB2A","#000000","#6D6D6D"))
p



chr_sig <- read.csv('/Users/darmokandjalad/Documents/Hop1PHD-Loop2/RAnalysis/chrLeng_log2_dsb.csv')

p <- ggplot(data = chr_sig, aes(x = length, y = sig, group = group, col = group)) + 
        geom_point(size = 3) +
        geom_smooth(method = "lm", se = FALSE) +
        scale_color_manual(values = c("#0099CD","#B72467"))
p


p <- ggplot(data = chr_sig, aes(x = length, y = Actual_sig_genomeAveNorm, group = group, col = group)) + 
        geom_point(size = 3) +
        geom_smooth(method = "lm", se = FALSE) +
        scale_color_manual(values = c("#0099CD","#B72467"))
p

p <- ggplot(data = chr_sig, aes(x = Chr, y = overRad50, fill = group, col = group)) + 
        geom_bar(stat = "identity",position = "dodge") +
        scale_fill_manual(values = c("#0099CD","#B72467"))+
        scale_color_manual(values = c("#0099CD","#B72467"))

p














ggplot2_blanktheme <- theme_classic() +
        theme(axis.line=element_blank(),axis.text.x=element_blank(),
              axis.text.y=element_blank(),axis.ticks=element_blank(),
              axis.title.x=element_blank(),
              axis.title.y=element_blank(),
              panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),plot.background=element_blank(),
              title = NULL)
theme_set(ggplot2_blanktheme)

######################
#dsb plot signal on whole genome
####################
setwd('/Users/darmokandjalad/Documents/Hop1PHD-Loop2/RAnalysis')
bedgraphs <- c("ccSeq_dsbsig_11687_reps_nozero.bed",
               "ccSeq_dsbsig_11570_reps_nozero.bed",
               "ccSeq_dsbsig_11687_reps.bed",
               "ccSeq_dsbsig_11570_reps.bed")

plot_signal <- function(sample,chrnum,tile) {
        genome_info <- hwglabr2::get_chr_coordinates(genome = 'SK1Yue')
        genome_info <- dropSeqlevels(genome_info, "chrIII", pruning.mode = "coarse")
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
        Sig <- import.bed(bedgraph_file)
        
        max_score=max(Sig$score)
        
        genome_info <- hwglabr2::get_chr_coordinates(genome = "SK1Yue")
        genome_info <- dropSeqlevels(genome_info, "chrIII", pruning.mode = "coarse")
        chrs <- list('I','VI','IX','VIII','V','XI','X','XIV','II','XIII','XVI','XII','VII','XV','IV')
        #chrs <- list('I')
        p=list()
        for(i in chrs) {
                print(i)
                p[[i]] = plot_signal(Sig,i,1000)
        }
        
        clusters = rtracklayer::import.bed('clusters_joined.bed')
        clusters <- dropSeqlevels(clusters, "chrIII", pruning.mode = "coarse")
        
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



#################################
##IslandVdeserts ----on hotspots
########################################
# function to make heatmap of signal at cluster and desert hotspots
setwd('/Users/darmokandjalad/Documents/Hop1PHD-Loop2/RAnalysis')
bedgraphs <- c("ccSeq_dsbsig_11687_reps.bed",
               "ccSeq_dsbsig_11570_reps.bed")
##for11570    
#11687
        signal <- import.bed("ccSeq_dsbsig_11687_reps.bed")
        signal <- signal[signal$score > 0]
        
        setwd('/Users/darmokandjalad/Documents/Hop1PHD-Loop2/RAnalysis')
        clusters = rtracklayer::import.bed('clusters_joined.bed')
        deserts = rtracklayer::import.bed('deserts_joined.bed')
        clusters <- dropSeqlevels(clusters, "chrIII", pruning.mode = "coarse")
        deserts <- dropSeqlevels(deserts, "chrIII", pruning.mode = "coarse")
        hotspots = hwglabr2::get_dsb_hotspots("SK1Yue")
        hotspots <- dropSeqlevels(hotspots, "chrIII", pruning.mode = "coarse")
        #axis = hwglabr2::get_Red1_summits("SK1Yue")
        
        #hits = findOverlaps(axis,signal)
        #signal_axis <- signal[subjectHits(hits)]
        #2302
        #rm(hits)
        
        #hits = findOverlaps(hotspots, signal)
        #signal_hotspot <- signal[subjectHits(hits)]
        #rm(hits)
        #10002
        
        #hotspots=3191
        #clusters=126
        hits = findOverlaps(clusters, hotspots)
        clusters_hotspot <- hotspots[subjectHits(hits)]
        rm(hits)
        #594
        hits = findOverlaps(deserts, hotspots)
        deserts_hotspot <- hotspots[subjectHits(hits)]
        rm(hits)
        #2613
        
        loopIlsand <- mean(clusters_hotspot$score)
        loopDesert <- mean(deserts_hotspot$score)
        #hits = findOverlaps(clusters,signal_hotspot)
        #signal_hotspot_cluster <- signal_hotspot[subjectHits(hits)]
        #rm(hits)
        #1559
        #hits = findOverlaps(deserts,signal_hotspot)
        #signal_hotspot_desert <- signal_hotspot[subjectHits(hits)]
        #rm(hits)
        #8443
        #mcols(signal_hotspot_desert)['class'] = 'desert'
        #mcols(signal_hotspot_cluster)['class'] = 'cluster'
        
        
        #ok fine now if we change extend to only say 500, is this better? 
        #it is supposed to be 1000
        mat1 <- normalizeToMatrix(signal, clusters_hotspot, value_column = "score",
                                  extend = 1000, mean_mode = "weighted", w = 5,empty_value=NA)
        #rows are 1:441(for 11570)
        #rows are 1:446 for 11687
        mat1_avrg <- hwglabr2::signal_mean_and_ci(mat1,
                                                  ci=0.95, rep_bootstrap=1000,
                                                  na_rm=TRUE)
        #u1 <- mat1_avrg[grep("u", rownames(mat1_avrg)),]
        #d1 <- mat1_avrg[grep("d", rownames(mat1_avrg)),]
        #mat1_avrg <- rbind(u1, d1)
        
        mat2 <- normalizeToMatrix(signal, deserts_hotspot, value_column = "score",
                                  extend = 1000, mean_mode = "weighted", w = 5,empty_value=NA)
        #rows are 1:446 for 11570 and 11687
        mat2_avrg <- hwglabr2::signal_mean_and_ci(mat2,
                                                  ci=0.95, rep_bootstrap=1000,
                                                  na_rm=TRUE)
        #t2 <- grep("t", rownames(mat2_avrg), value = TRUE)
        #mat2_avrg <- mat2_avrg[rownames(mat2_avrg)!=t2,]
        #u2 <- mat2_avrg[grep("u", rownames(mat2_avrg)),]
        #d2 <- mat2_avrg[grep("d", rownames(mat2_avrg)),]
        #mat2_avrg <- rbind(u2, d2)
        mat1_avrg_df <- data.frame(Position=seq(-219, 221), sample = 'cluster',mat1_avrg)
        mat2_avrg_df <- data.frame(Position=seq(-222, 223), sample = 'desert',mat2_avrg)
        
        #mat1_avrg_df <- data.frame(Position=seq(-199, 200), sample = 'cluster',mat1_avrg)
        #mat2_avrg_df <- data.frame(Position=seq(-199, 200), sample = 'desert',mat2_avrg)
        alldata = rbind(mat1_avrg_df,mat2_avrg_df)
        pdf("11687_1kbsig_100bin_hotspot_islandVsDesertBootMetaPlot.pdf", height = 5, width = 5)
        print(ggplot(alldata, aes(Position, Mean, group = sample, colour=sample,fill=sample)) +
                             labs(title = element_blank())+
                             geom_line() +
                             geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha=0.3,colour=NA) +
                             geom_vline(xintercept = 0, lty = 3) +
                             scale_x_continuous(breaks = c(-220, 0, 220),labels = c("-1 kb", "Axis summit", "1 kb")) +
                             scale_color_manual(values = c("#00A14B","#662D91"))+
                             theme(legend.position = "none"))
        dev.off()
        
#now repeat for 11570     

signal <- import.bed("ccSeq_dsbsig_11570_reps.bed")
signal <- signal[signal$score > 0]

setwd('/Users/darmokandjalad/Documents/Hop1PHD-Loop2/RAnalysis')
clusters = rtracklayer::import.bed('clusters_joined.bed')
deserts = rtracklayer::import.bed('deserts_joined.bed')
clusters <- dropSeqlevels(clusters, "chrIII", pruning.mode = "coarse")
deserts <- dropSeqlevels(deserts, "chrIII", pruning.mode = "coarse")
hotspots = hwglabr2::get_dsb_hotspots("SK1Yue")
hotspots <- dropSeqlevels(hotspots, "chrIII", pruning.mode = "coarse")
#axis = hwglabr2::get_Red1_summits("SK1Yue")

#hits = findOverlaps(axis,signal)
#signal_axis <- signal[subjectHits(hits)]
#2302
#rm(hits)

hits = findOverlaps(hotspots, signal)
signal_hotspot <- signal[subjectHits(hits)]
#rm(hits)
#10002

#hotspots=3191
#clusters=126
hits = findOverlaps(clusters, signal_hotspot)
clusters_hotspot <- signal_hotspot[subjectHits(hits)]
rm(hits)
#594
hits = findOverlaps(deserts, signal_hotspot)
deserts_hotspot <- signal_hotspot[subjectHits(hits)]
rm(hits)
#2613
Wt_cluser <- mean(clusters_hotspot$score)
Wt_desert <- mean(deserts_hotspot$score)

Wenrich <- Wt_cluser/Wt_desert

loop_cluster <- mean(clusters_hotspot$score)
loop_desert <- mean(deserts_hotspot$score)

loopenrich <- loop_cluster / loop_desert

#hits = findOverlaps(clusters,signal_hotspot)
#signal_hotspot_cluster <- signal_hotspot[subjectHits(hits)]
#rm(hits)
#1559
#hits = findOverlaps(deserts,signal_hotspot)
#signal_hotspot_desert <- signal_hotspot[subjectHits(hits)]
#rm(hits)
#8443
#mcols(signal_hotspot_desert)['class'] = 'desert'
#mcols(signal_hotspot_cluster)['class'] = 'cluster'


#ok fine now if we change extend to only say 500, is this better? 
#it is supposed to be 1000
mat1 <- normalizeToMatrix(signal, clusters_hotspot, value_column = "score",
                          extend = 1000, mean_mode = "weighted", w = 5,empty_value=NA)
#rows are 1:441(for 11570)
#rows are 1:446 for 11687
mat1_avrg <- hwglabr2::signal_mean_and_ci(mat1,
                                          ci=0.95, rep_bootstrap=1000,
                                          na_rm=TRUE)
#u1 <- mat1_avrg[grep("u", rownames(mat1_avrg)),]
#d1 <- mat1_avrg[grep("d", rownames(mat1_avrg)),]
#mat1_avrg <- rbind(u1, d1)

mat2 <- normalizeToMatrix(signal, deserts_hotspot, value_column = "score",
                          extend = 1000, mean_mode = "weighted", w = 5,empty_value=NA)
#rows are 1:446 for 11570 and 11687
mat2_avrg <- hwglabr2::signal_mean_and_ci(mat2,
                                          ci=0.95, rep_bootstrap=1000,
                                          na_rm=TRUE)
#t2 <- grep("t", rownames(mat2_avrg), value = TRUE)
#mat2_avrg <- mat2_avrg[rownames(mat2_avrg)!=t2,]
#u2 <- mat2_avrg[grep("u", rownames(mat2_avrg)),]
#d2 <- mat2_avrg[grep("d", rownames(mat2_avrg)),]
#mat2_avrg <- rbind(u2, d2)
mat1_avrg_df <- data.frame(Position=seq(-220, 220), sample = 'cluster',mat1_avrg)
mat2_avrg_df <- data.frame(Position=seq(-222, 223), sample = 'desert',mat2_avrg)

#mat1_avrg_df <- data.frame(Position=seq(-199, 200), sample = 'cluster',mat1_avrg)
#mat2_avrg_df <- data.frame(Position=seq(-199, 200), sample = 'desert',mat2_avrg)
alldata = rbind(mat1_avrg_df,mat2_avrg_df)
pdf("11570_1kbsig_100bin_hotspot_islandVsDesertBootMetaPlot.pdf", height = 5, width = 5)
print(ggplot(alldata, aes(Position, Mean, group = sample, colour=sample,fill=sample)) +
              labs(title = element_blank())+
              geom_line() +
              geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha=0.3,colour=NA) +
              geom_vline(xintercept = 0, lty = 3) +
              scale_x_continuous(breaks = c(-220, 0, 220),labels = c("-1 kb", "Axis summit", "1 kb")) +
              scale_color_manual(values = c("#00A14B","#662D91"))+
              theme(legend.position = "none"))
dev.off()  
        
######################################
#now for centromeres 
#use 20kb cen flank for comparision in paper
#################################################
genome_info <- hwglabr2::get_chr_coordinates(genome = 'SK1Yue')
binsc <- toDataframe(genome_info)
Cen_flank = 20000
cenStartminus <- binsc$start - Cen_flank
cenEndplus <- binsc$end + Cen_flank
chr <- binsc$chr
Cens_bpFlank <- data.frame(first_column = chr, second_column= cenStartminus, third_column = cenEndplus)
Cens_bpFlank <- toGRanges(Cens_bpFlank)
df <- toDataframe(Cens_bpFlank)
write.table(df, file="Cen10KbFlankRegions.bed", quote=F, sep="\t", row.names=F, col.names=F)
Cens_bpFlank <- dropSeqlevels(Cens_bpFlank, "chrIII", pruning.mode = "coarse")


signal <- import.bed("ccSeq_dsbsig_11687_reps.bed")
signal <- signal[signal$score > 0]
Cen_flank = 20000

#get genome
genome_info <- hwglabr2::get_chr_coordinates(genome = 'SK1Yue')
genome_info <- dropSeqlevels(genome_info, "chrIII", pruning.mode = "coarse")
signal <- sort(GenomeInfoDb::sortSeqlevels(signal))
genome_info <- GenomeInfoDb::sortSeqlevels(genome_info)

hotspots = hwglabr2::get_dsb_hotspots("SK1Yue")
hotspots <- dropSeqlevels(hotspots, "chrIII", pruning.mode = "coarse")


hits = findOverlaps(Cens_bpFlank,hotspots)
hotspots_cens <- hotspots[subjectHits(hits)]
rm(hits)
hotspots_offcen <- hotspots[hotspots %outside% hotspots_cens]




mat1 <- normalizeToMatrix(signal, hotspots_cens, value_column = "score",
                          extend = 1000, mean_mode = "weighted", w = 5,empty_value=NA)

mat1_avrg <- hwglabr2::signal_mean_and_ci(mat1,
                                          ci=0.95, rep_bootstrap=1000,
                                          na_rm=TRUE)

mat2 <- normalizeToMatrix(signal, hotspots_offcen, value_column = "score",
                          extend = 1000, mean_mode = "weighted", w = 5,empty_value=NA)

mat2_avrg <- hwglabr2::signal_mean_and_ci(mat2,
                                          ci=0.95, rep_bootstrap=1000,
                                          na_rm=TRUE)
#20Kb
mat1_avrg_df <- data.frame(Position=seq(-223, 221), sample = 'OnCen',mat1_avrg)
mat2_avrg_df <- data.frame(Position=seq(-223, 221), sample = 'OffCen',mat2_avrg)
#10kb
mat1_avrg_df <- data.frame(Position=seq(-220, 220), sample = 'OnCen',mat1_avrg)
mat2_avrg_df <- data.frame(Position=seq(-223, 221), sample = 'OffCen',mat2_avrg)


alldata = rbind(mat1_avrg_df,mat2_avrg_df)

pdf("11687_1kbsig_100bin_hotspot_OnVoffCen20kBootMetaPlot.pdf", height = 5, width = 5)
print(ggplot(alldata, aes(Position, Mean, group = sample, colour=sample,fill=sample)) +
              labs(title = element_blank())+
              geom_line() +
              geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha=0.3,colour=NA) +
              geom_vline(xintercept = 0, lty = 3) +
              scale_x_continuous(breaks = c(-220, 0, 220),labels = c("-1 kb", "Axis summit", "1 kb")) +
              scale_color_manual(values = c("#00A14B","#662D91"))+
              theme(legend.position = "none"))
dev.off()
p <- ggplot(alldata, aes(Position, Mean, group = sample, colour=sample,fill=sample)) +
        labs(title = element_blank())+
        geom_line() +
        geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha=0.3,colour=NA) +
        geom_vline(xintercept = 0, lty = 3) +
        scale_x_continuous(breaks = c(-220, 0, 220),labels = c("-1 kb", "Axis summit", "1 kb")) +
        scale_color_manual(values = c("#00A14B","#662D91"))+
        theme(legend.position = "none")
p


signal <- import.bed("ccSeq_dsbsig_11570_reps.bed")
signal <- signal[signal$score > 0]
Cen_flank = 20000

#get genome
genome_info <- hwglabr2::get_chr_coordinates(genome = 'SK1Yue')
genome_info <- dropSeqlevels(genome_info, "chrIII", pruning.mode = "coarse")
signal <- sort(GenomeInfoDb::sortSeqlevels(signal))
genome_info <- GenomeInfoDb::sortSeqlevels(genome_info)

hotspots = hwglabr2::get_dsb_hotspots("SK1Yue")
hotspots <- dropSeqlevels(hotspots, "chrIII", pruning.mode = "coarse")


hits = findOverlaps(Cens_bpFlank,hotspots)
hotspots_cens <- hotspots[subjectHits(hits)]
rm(hits)
hotspots_offcen <- hotspots[hotspots %outside% hotspots_cens]




mat1 <- normalizeToMatrix(signal, hotspots_cens, value_column = "score",
                          extend = 1000, mean_mode = "weighted", w = 5,empty_value=NA)

mat1_avrg <- hwglabr2::signal_mean_and_ci(mat1,
                                          ci=0.95, rep_bootstrap=1000,
                                          na_rm=TRUE)

mat2 <- normalizeToMatrix(signal, hotspots_offcen, value_column = "score",
                          extend = 1000, mean_mode = "weighted", w = 5,empty_value=NA)

mat2_avrg <- hwglabr2::signal_mean_and_ci(mat2,
                                          ci=0.95, rep_bootstrap=1000,
                                          na_rm=TRUE)

mat1_avrg_df <- data.frame(Position=seq(-223, 221), sample = 'OnCen',mat1_avrg)
mat2_avrg_df <- data.frame(Position=seq(-223, 221), sample = 'OffCen',mat2_avrg)


alldata = rbind(mat1_avrg_df,mat2_avrg_df)

pdf("11570_1kbsig_100bin_hotspot_OnVoffCen20kBootMetaPlot.pdf", height = 5, width = 5)
print(ggplot(alldata, aes(Position, Mean, group = sample, colour=sample,fill=sample)) +
              labs(title = element_blank())+
              geom_line() +
              geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha=0.3,colour=NA) +
              geom_vline(xintercept = 0, lty = 3) +
              scale_x_continuous(breaks = c(-220, 0, 220),labels = c("-1 kb", "Axis summit", "1 kb")) +
              scale_color_manual(values = c("#00A14B","#662D91"))+
              theme(legend.position = "none"))
dev.off()

####now for wt


        
        





















##########################################
##cc-seq plot signal with zoom and wildtype or hop1-loop2 signal
#############################
##for11687for
#################################

#set Parameters... i just did this for all genotypes
setwd('/Users/darmokandjalad/Documents/Hop1PHD-Loop2/RAnalysis')
chrnum = "chrXIII"
antibodyTarget = "ccSeq"
genotype = "11687"
lpN = rtracklayer::import.bed('/Users/darmokandjalad/Documents/Hop1PHD-Loop2/RAnalysis/SpikeValue/Hop1_loop2_Reps_Norm.bed')
        signal <- lpN
        signal <- dropSeqlevels(signal, "chrIII", pruning.mode = "coarse")
        genome_info <- hwglabr2::get_chr_coordinates(genome = 'SK1Yue')
        genome_info <- dropSeqlevels(genome_info, "chrIII", pruning.mode = "coarse")
        centromeres <- hwglabr2::get_chr_coordinates(genome='SK1Yue', as_df=FALSE)
        centromeres <- dropSeqlevels(centromeres, "chrIII", pruning.mode = "coarse")
        cen = centromeres[seqnames(centromeres) == chrnum]
        cen <- round((start(cen) + end(cen))/2)
        cen <- data.frame(cen = cen, y = 0)
        cen <- cen/1000
        
        
        clusters <- rtracklayer::import.bed('/Users/darmokandjalad/Documents/Hop1PHD-Loop2/RAnalysis/clusters_joined.bed')
        clsters <- clusters[seqnames(clusters) == chrnum]
        clstersDF <- data.frame(clsters)
        Start <- clstersDF$start/1000
        End <- clstersDF$end/1000
        clstersDF <- subset(clstersDF, select = -start)
        clstersDF <- subset(clstersDF, select = -end)
        clstersDF$start <- Start
        clstersDF$end <- End
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
        lp2CHPdf <- data.frame(position=positions / 1000, signal=bins$binned_score)


signal <- import.bed("ccSeq_dsbsig_11687_reps.bed")
chrnum = "chrXIII"
antibodyTarget = "ccSeq"
genotype = "11687"

                genome_info <- hwglabr2::get_chr_coordinates(genome = 'SK1Yue')
                genome_info <- dropSeqlevels(genome_info, "chrIII", pruning.mode = "coarse")
                centromeres <- hwglabr2::get_chr_coordinates(genome='SK1Yue', as_df=FALSE)
                centromeres <- dropSeqlevels(centromeres, "chrIII", pruning.mode = "coarse")
                cen = centromeres[seqnames(centromeres) == chrnum]
                cen <- round((start(cen) + end(cen))/2)
                cen <- data.frame(cen = cen, y = 0)
                cen <- cen/1000
                
                
                clusters <- rtracklayer::import.bed('/Users/darmokandjalad/Documents/Hop1PHD-Loop2/RAnalysis/clusters_joined.bed')
                clsters <- clusters[seqnames(clusters) == chrnum]
                clstersDF <- data.frame(clsters)
                Start <- clstersDF$start/1000
                End <- clstersDF$end/1000
                clstersDF <- subset(clstersDF, select = -start)
                clstersDF <- subset(clstersDF, select = -end)
                clstersDF$start <- Start
                clstersDF$end <- End

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
                        loop2df <- data.frame(position=positions / 1000, signal=bins$binned_score)
                    mxMatch <- lp2CHPdf$signal*4
                    lp2CHPdf$signal <- mxMatch
                        loop2df$group <- "2loopCCseq"
                        lp2CHPdf$group <- "1loopCHIP"
                        theme_set(theme_classic())                       
                        All_loop_df <- rbind(lp2CHPdf,loop2df)
                        p <- ggplot(All_loop_df, aes(x=position, y=signal, colour = group)) +
                                geom_line()+
                                scale_color_manual(values = c("#0099CD","#B72467"))
                        p
                        p <- p + geom_point(cen, mapping=aes(cen,0),
                                            size = 5.0, colour = 'black', shape=17) + xlab(chrnum) + ylab(antibodyTarget)
                        p <- p + geom_segment(clstersDF, mapping=aes(x=start, xend=end, y=0, yend=0),colour = 'black')
                        p
                        p <- p + geom_vline(xintercept = 625, linetype="dashed")
                        p <- p + geom_vline(xintercept = 725, linetype="dashed")
                        
                        p
                        
                        #The plot with a zoom
                        p2 <- ggplot(All_loop_df, aes(x=position, y=signal, colour = group)) +
                                scale_color_manual(values = c("#0099CD","#B72467")) 
                        p2 <- p2 + ylim(-1,35) + xlim(625,725) + geom_line(aes(x=position, y=signal), alpha=.9)
                        p2 <- p2 + geom_point(cen, mapping=aes(cen,0),
                                              size = 5.0, colour = 'black', shape = 17) + xlab(chrnum) + ylab(antibodyTarget)
                        p2 <- p2 + geom_segment(clstersDF, mapping=aes(x=start, xend=end, y=0, yend=0),colour = 'black')
                        p2 <- p2  + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank())
                        p2
                        
                        
########for11687for                        
#######################################  
##for11570                       
#################################                        
##
signal <- import.bed("ccSeq_dsbsig_11570_reps.bed")
chrnum = "chrXIII"
antibodyTarget = "ccSeq"
genotype = "11570"
genome_info <- hwglabr2::get_chr_coordinates(genome = 'SK1Yue')
genome_info <- dropSeqlevels(genome_info, "chrIII", pruning.mode = "coarse")
centromeres <- hwglabr2::get_chr_coordinates(genome='SK1Yue', as_df=FALSE)
centromeres <- dropSeqlevels(centromeres, "chrIII", pruning.mode = "coarse")
cen = centromeres[seqnames(centromeres) == chrnum]
cen <- round((start(cen) + end(cen))/2)
cen <- data.frame(cen = cen, y = 0)
cen <- cen/1000


clusters <- rtracklayer::import.bed('/Users/darmokandjalad/Documents/Hop1PHD-Loop2/RAnalysis/clusters_joined.bed')
clsters <- clusters[seqnames(clusters) == chrnum]
clstersDF <- data.frame(clsters)
Start <- clstersDF$start/1000
End <- clstersDF$end/1000
clstersDF <- subset(clstersDF, select = -start)
clstersDF <- subset(clstersDF, select = -end)
clstersDF$start <- Start
clstersDF$end <- End
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
                        Wtdf <- data.frame(position=positions / 1000, signal=bins$binned_score)
 
                        
                                               
Wtsig <- hwglabr2::import_bedGraph('/Users/darmokandjalad/Documents/Hop1PHD-Loop2/RAnalysis/Hop1-wildtype-334-340-32-177-Reps-SK1Yue-PM_B3W3_MACS2_FE.bdg.gz')
Wtsig <- dropSeqlevels(Wtsig, "chrIII", pruning.mode = "coarse")
signal <- Wtsig
chrnum = "chrXIII"
antibodyTarget = "ccSeq"
genotype = "11570"
genome_info <- hwglabr2::get_chr_coordinates(genome = 'SK1Yue')
genome_info <- dropSeqlevels(genome_info, "chrIII", pruning.mode = "coarse")
centromeres <- hwglabr2::get_chr_coordinates(genome='SK1Yue', as_df=FALSE)
centromeres <- dropSeqlevels(centromeres, "chrIII", pruning.mode = "coarse")
cen = centromeres[seqnames(centromeres) == chrnum]
cen <- round((start(cen) + end(cen))/2)
cen <- data.frame(cen = cen, y = 0)
cen <- cen/1000


clusters <- rtracklayer::import.bed('/Users/darmokandjalad/Documents/Hop1PHD-Loop2/RAnalysis/clusters_joined.bed')
clsters <- clusters[seqnames(clusters) == chrnum]
clstersDF <- data.frame(clsters)
Start <- clstersDF$start/1000
End <- clstersDF$end/1000
clstersDF <- subset(clstersDF, select = -start)
clstersDF <- subset(clstersDF, select = -end)
clstersDF$start <- Start
clstersDF$end <- End
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
WtCHIPdf <- data.frame(position=positions / 1000, signal=bins$binned_score)

FixMax <- WtCHIPdf$signal*4
WtCHIPdf$signal <- FixMax
                        
                        
                        
WtCHIPdf$group <- "1wtCHIP"
Wtdf$group <- "2wtCCseq"

wtdf <- rbind(Wtdf,WtCHIPdf)


## Set Ggplot theme 
theme_set(theme_classic())

#All, too much on plot, cant see it. 
#plot(df, xlab='Position on chrV (Kb)', ylab='Average signal', type='l', lwd=2, main="rec8D 5187", ylim = c(0,8))
p <- ggplot(wtdf, aes(x=position, y=signal, colour = group)) +
        geom_line()+
        scale_color_manual(values = c("#0099CD","#B72467"))
p
p <- p + geom_point(cen, mapping=aes(cen,0),
                    size = 5.0, colour = 'black', shape=17) + xlab(chrnum) + ylab(antibodyTarget)
p <- p + geom_segment(clstersDF, mapping=aes(x=start, xend=end, y=0, yend=0),colour = 'black')
p
p <- p + geom_vline(xintercept = 625, linetype="dashed")
p <- p + geom_vline(xintercept = 725, linetype="dashed")

p
575,725
#The plot with a zoom
p2 <- ggplot(wtdf, aes(x=position, y=signal, colour = group)) +
        scale_color_manual(values = c("#0099CD","#B72467")) 
p2 <- p2 + ylim(-1,35) + xlim(625,725) + geom_line(aes(x=position, y=signal), alpha=.9)
p2 <- p2 + geom_point(cen, mapping=aes(cen,0),
                      size = 5.0, colour = 'black', shape = 17) + xlab(chrnum) + ylab(antibodyTarget)
p2 <- p2 + geom_segment(clstersDF, mapping=aes(x=start, xend=end, y=0, yend=0),colour = 'black')
p2 <- p2  + theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank())
p2







########for11570for 
############################################
########whats next?
################################




q<- ggplot(alldata, aes(Position, Mean, group = sample, colour=sample,fill=sample)) +
        labs(title = element_blank())+
        geom_line() +
        geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha=0.3,colour=NA) +
        geom_vline(xintercept = 0, lty = 3) +
        scale_x_continuous(breaks = c(-220, 0, 220),labels = c("-1 kb", "Axis summit", "1 kb")) +
        scale_color_manual(values = c("#00A14B","#662D91"))
q




p<-ggplot(alldata, aes(Position, Mean, group = sample, colour=sample, fill=sample)) +
        labs(title = element_blank())+
        geom_line() +
        geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha=0.3,colour=NA) +
        geom_vline(xintercept = 0, lty = 3) +
        scale_x_continuous(breaks = c(-199, 0, 200),labels = c("-1 kb", "Axis summit", "1 kb")) +
        scale_color_manual(values = c("#00A14B","#662D91"))+
        theme(legend.position = "none")
p

p<-ggplot(alldata, aes(Position, Mean, group = sample)) +
        labs(title = element_blank())+
        geom_line()
       
p

#seq in the middle of the Hotspot will be depleted
#I am better off looking next to hotspots.
#not sure how to do this. 
#also why is this so jagged am I supposed to smooth? 








setwd('/Users/darmokandjalad/Documents/Hop1PHD-Loop2/RAnalysis')
clusters = rtracklayer::import.bed('clusters_joined.bed')
deserts = rtracklayer::import.bed('deserts_joined.bed')
#hotspots = hwglabr2::get_dsb_hotspots("SK1Yue")
axis = hwglabr2::get_Red1_summits("SK1Yue")

hits = findOverlaps(axis,dsbSig)
dsbSig_cluster <- dsbSig[subjectHits(hits)]
rm(hits)




hits = findOverlaps(clusters,dsbSig)
dsbSig_cluster <- dsbSig[subjectHits(hits)]
rm(hits)
hits = findOverlaps(deserts,dsbSig)
dsbSig_desert <- dsbSig[subjectHits(hits)]
rm(hits)
subset(axis_cluster, (name %in% axis_desert$name))
mcols(dsbSig_desert)['class'] = 'desert'
mcols(dsbSig_cluster)['class'] = 'cluster'

###
clusters = rtracklayer::import.bed('clusters_joined.bed')
deserts = rtracklayer::import.bed('deserts_joined.bed')

              



###









mat1 <- normalizeToMatrix(dsbSig, dsbSig_cluster, value_column = "score",
                          extend = 1000, mean_mode = "weighted", w = 5,empty_value=NA)

mat1_avrg <- hwglabr2::signal_mean_and_ci(mat1,
                                          ci=0.95, rep_bootstrap=1000,
                                          na_rm=TRUE)

mat2 <- normalizeToMatrix(dsbSig, dsbSig_desert, value_column = "score",
                          extend = 1000, mean_mode = "weighted", w = 5,empty_value=NA)

mat2_avrg <- hwglabr2::signal_mean_and_ci(mat2,
                                          ci=0.95, rep_bootstrap=1000,
                                          na_rm=TRUE)

mat1_avrg_df <- data.frame(Position=seq(-199, 200), sample = 'cluster',mat1_avrg)
mat2_avrg_df <- data.frame(Position=seq(-199, 200), sample = 'desert',mat2_avrg)
alldata = rbind(mat1_avrg_df,mat2_avrg_df)
pdf("All_dsbSig_BootMetaPlot.pdf", height = 5, width = 5)
print(ggplot(alldata, aes(Position, Mean, group = sample, colour=sample, fill=sample)) +
        labs(title = element_blank())+
        geom_line() +
        geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha=0.3,colour=NA) +
        geom_vline(xintercept = 0, lty = 3) +
        scale_x_continuous(breaks = c(-199, 0, 200),labels = c("-1 kb", "Axis summit", "1 kb")) +
        scale_color_manual(values = c("#00A14B","#662D91"))+
        theme(legend.position = "none"))
dev.off()


p<-ggplot(alldata, aes(Position, Mean, group = sample, colour=sample, fill=sample)) +
        labs(title = element_blank())+
        geom_line() +
        geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha=0.3,colour=NA) +
        geom_vline(xintercept = 0, lty = 3) +
        scale_x_continuous(breaks = c(-199, 0, 200),labels = c("-1 kb", "Axis summit", "1 kb")) +
        scale_color_manual(values = c("#00A14B","#662D91"))
p

