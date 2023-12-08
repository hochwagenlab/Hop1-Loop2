#---------------------------------------------------------------#
# Axis to adapt to enrichement of axis in islands vs deserts  (deserts = non-islands)                 #
# Figure 2 code                                                 #
#---------------------------------------------------------------#
# load packages
library(GenomicRanges)
library(regioneR)
library(hwglabr2)
library(EnrichedHeatmap)
library(ggplot2)
library(patchwork)

# Create working folder with necessary files 
setwd('/Users/darmokandjalad/Documents/Hop1PHD-Loop2/RAnalysis')

#----------------------------------------------------------------#
# Figure 2A  This is the heatmap clusters and deserts                                                    #
#----------------------------------------------------------------#




clusters = rtracklayer::import.bed('clusters_joined.bed')
deserts = rtracklayer::import.bed('deserts_joined.bed')
axisSites = hwglabr2::get_Red1_summits(genome = "SK1Yue")

hits = findOverlaps(clusters,axisSites)
overlaps <- pintersect(clusters[queryHits(hits)], axisSites[subjectHits(hits)])
percentOverlap <- width(overlaps) / width(axisSites[subjectHits(hits)])
hits <- hits[percentOverlap > 0.5]
axisSites_cluster <- axisSites[subjectHits(hits)]
rm(hits);rm(overlaps)

hits = findOverlaps(deserts,axisSites)
overlaps <- pintersect(deserts[queryHits(hits)], axisSites[subjectHits(hits)])
percentOverlap <- width(overlaps) / width(axisSites[subjectHits(hits)])
hits <- hits[percentOverlap > 0.5]
axisSites_desert <- axisSites[subjectHits(hits)]

mcols(axisSites_desert)['class'] = 'desert'
mcols(axisSites_cluster)['class'] = 'cluster'

axisSites_desert_sort = axisSites_desert[order(width(axisSites_desert),decreasing = T)]
axisSites_desert_sort = axisSites_desert[order(score(axisSites_desert),decreasing = T)]
axisSites_cluster_sort = axisSites_cluster[order(width(axisSites_cluster),decreasing = T)]
axisSites_cluster_sort = axisSites_cluster[order(score(axisSites_cluster),decreasing = T)]

axisSites_all = c(axisSites_desert_sort,axisSites_cluster_sort)
midpoint <- floor(width(axisSites_all) / 2)
start(axisSites_all) <- start(axisSites_all) + midpoint
end(axisSites_all) <- start(axisSites_all)

#dir.create("axisSites_pdf")
bedgraphs <- c("wildtype-Hop1-334-340-32-177-Reps-SK1Yue-PM_B3W3_MACS2_FE.bdg copy.gz",
               "11644-NHantiHop1-20220331-20220125-20221004-Reps-SK1Yue-PM-PE_B4_W3_MACS2_FE.bdg.gz",
               "11757-hop1NH-20220125-20220721-20221004-Reps-SK1Yue-PM-PE_B4_W3_MACS2_FE.bdg.gz",
               "11758-hop1NH-20220125-20220721-20221004-Reps-SK1Yue-PM-PE_B4_W3_MACS2_FE.bdg.gz",
               "rec8D-Hop1-61-91-Reps-SK1Yue-PM_B3W3_MACS2_FE.bdg.gz")
#bedgraphs <- c("wildtype-Hop1-334-340-32-177-Reps-SK1Yue-PM_B3W3_MACS2_FE.bdg copy.gz")

# function to make heatmap of signal at cluster and desert hotspots
for (i in 1:length(bedgraphs)) {
  bedgraph_file <- bedgraphs[i]
  name <- strsplit(bedgraph_file,split = "-")
  name <- sapply(name, "[[",1)
  Red1_bg <- hwglabr2::import_bedGraph(bedgraph_file, local_copy=T)
  genAvg <- hwglabr2::average_chr_signal(Red1_bg)$genome_avrg
  Red1_bg$score <- Red1_bg$score/genAvg
  
  mat1 <- normalizeToMatrix(Red1_bg, axisSites_all, value_column = "score",
                            extend = 1000, mean_mode = "weighted", w = 5,empty_value=NA)
  
  col_fun <- circlize::colorRamp2(quantile(mat1, c( 0.01,0.25, 0.5, 0.75, 0.95),na.rm=T), c("#662D91","#9C77BC", "white","#87D8A9", "#00A14B"))
  
  pdf(paste0(name[1],"sigOn_axisSites_heatmap.pdf"), height = 8, width = 4)
  print(EnrichedHeatmap(mat1, col = col_fun, name = "Signal",pos_line=FALSE, use_raster=FALSE,
                        top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = c("#00A14B","#662D91")),
                                                                                 show_error = TRUE,pos_line=FALSE)),
                        row_title_rot = 0,
                        axis_name = c("-1 kb", "axisSites", "1 kb"),
                        row_order = 1:length(axisSites_all),
                        split=axisSites_all$class,
                        column_title ="Hop1 signal"))
  dev.off()
}

