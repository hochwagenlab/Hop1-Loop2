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




via <- read.csv('/Users/darmokandjalad/Documents/Hop1PHD-Loop2/RAnalysis/SporeViability/SporeViability_summary.csv')

##This is spore viaiblity 
df<-data.frame(Mean=c(92.5, 89.56, 96.14, 95.68, 41.62, 74.81666667),
               sd=c(3.908324449,2.707951255,2.07557221,2.735324478,13.06797612,5.482852056),
               Category=c("5tel1", "6tel1 hop1 loop2", "1Wild type", "3pch2", "4pch2 hop1 loop2", "2hop1 loop2")
)

# Load ggplot2
library(ggplot2)

p <- ggplot(df, aes(x=Category, y=Mean)) +
  geom_bar(position=position_dodge(), stat="identity",
           fill = '#0099CD') +
  geom_errorbar(aes(ymin=Mean-sd, ymax=Mean+sd), width=.2)
p










