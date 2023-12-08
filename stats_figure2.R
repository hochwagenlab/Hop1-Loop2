## check for statistical significance 
#is data normally distributed?

library(dplyr)
library(ggpubr)


#Calculate Statistics for Viability data
#######################################
Viability_data <- read.csv('/Users/darmokandjalad/Documents/Hop1PHD-Loop2/RAnalysis/sporeViability.csv')

#We start by displaying a random sample of 10 rows using the function sample_n()[in dplyr package].

#Show 10 random rows:
  
set.seed(1234)
dplyr::sample_n(Viability_data, 10)

ggdensity(Viability_data$Viability, main = "Density plot of yeast Viaiblity ", xlab="Viability")
ggqqplot(Viability_data$Viability)

shapiro.test(Viability_data$Viability)
#Shapiro-Wilk normality test

#data:  Viability_data$Viability
#W = 0.83121, p-value = 0.009517
#From the output, the p-value > 0.05 implying that the distribution of the data are 
#not significantly different from normal distribution. 
#In other words, we can assume the normality.

Wt_viab <- subset(Viability_data, Viability_data$Strain == "wt")
x <- Wt_viab$Viability
lp_viab <- subset(Viability_data, Viability_data$Strain =='loop2')
y <- lp_viab$Viability
lpPch_viab <- subset(Viability_data, Viability_data$Strain == "pch2loop2")
y <- lpPch_viab$Viability


wilcox.test(x,y, paired = FALSE, exact = FALSE)
  #Wilcoxon rank sum test with continuity correction

  #data:  x and y
  #W = 8.5, p-value = 0.1212
  #alternative hypothesis: true location shift is not equal to 0
t.test(x,y, paired = FALSE)
  #Welch Two Sample t-test

  #data:  x and y
  #t = 2.6857, df = 3.9875, p-value = 0.05508
  #alternative hypothesis: true difference in means is not equal to 0
  #95 percent confidence interval:
  #  -0.5447029 31.6113695
  #sample estimates:
  #  mean of x mean of y 
  #92.43333  76.90000 


wilcox.test(x, y, paired = FALSE, alternative = "less")
  #Wilcoxon rank sum test with continuity correction

  #data:  x and y
  #W = 8.5, p-value = 0.9768
  #alternative hypothesis: true location shift is less than 0

lpPch_viab <- subset(Viability_data, Viability_data$Strain =='pch2loop2')
y <- lpPch_viab$Viability
wilcox.test(x,y, paired = FALSE, exact = FALSE)

  #Wilcoxon rank sum test with continuity correction

    #W = 9, p-value = 0.08086
  #alternative hypothesis: true location shift is not equal to 0
wilcox.test(x, y, paired = FALSE, alternative = "less")
  #Wilcoxon rank sum exact test

  #data:  x and y
  #W = 9, p-value = 1
  #alternative hypothesis: true location shift is less than 0




#IF
##########################################################
integrated_sig <- read.csv('/Users/darmokandjalad/Documents/Hop1PHD-Loop2/Manuscript/Hop1_Loop2_figures/Figure2/Edited_Intensity_Data-allstrains.csv')

Hop1_intSig <- subset(integrated_sig, integrated_sig$Channel == "Hop1")
ggdensity(Hop1_intSig$IntDen)
ggqqplot(Hop1_intSig$IntDen)
shapiro.test(Hop1_intSig$IntDen)

  #Shapiro-Wilk normality test

  #data:  Hop1_intSig$IntDen
  #W = 0.89988, p-value = 3.997e-12
#This is not a normal distribution of data


hop1H2 <- subset(Hop1_intSig, Hop1_intSig$Hour == "H2")
Wthop1_H2 <- subset(hop1H2, hop1H2$Label == "WT")
loop2Hop1_H2 <- subset(hop1H2, hop1H2$Label == "hop1-loop2")
wilcox.test(x= Wthop1_H2$IntDen, y = loop2Hop1_H2$IntDen, paired = FALSE)

  #Wilcoxon rank sum exact test

  #data:  Wthop1_H2$IntDen and loop2Hop1_H2$IntDen
  #W = 392, p-value = 9.721e-10
  #alternative hypothesis: true location shift is not equal to 0
 
lp2Pch2Hop1H2 <- subset(hop1H2, hop1H2$Label == "hop1-loop2-pch2Δ") 
wilcox.test(x=Wthop1_H2$IntDen, y = lp2Pch2Hop1H2$IntDen, paired = FALSE)
  #Wilcoxon rank sum exact test

  #data:  Wthop1_H2$IntDen and lp2Pch2Hop1H2$IntDen
  #W = 366, p-value = 9.249e-07
  #alternative hypothesis: true location shift is not equal to 0

pch2Hop1H2 <- subset(hop1H2, hop1H2$Label == "pch2Δ")
wilcox.test(x=Wthop1_H2$IntDen, y = pch2Hop1H2$IntDen, paired = FALSE)
  #Wilcoxon rank sum exact test

  #data:  Wthop1_H2$IntDen and pch2Hop1H2$IntDen
  #W = 253, p-value = 0.7251
  #alternative hypothesis: true location shift is not equal to 0




hop1H3 <- subset(Hop1_intSig, Hop1_intSig$Hour == "H3")
Wthop1_H3 <- subset(hop1H3, hop1H3$Label == "WT")
loop2Hop1_H3 <- subset(hop1H3, hop1H3$Label == "hop1-loop2")
wilcox.test(x= Wthop1_H3$IntDen, y = loop2Hop1_H3$IntDen)
  #Wilcoxon rank sum exact test

  #data:  Wthop1_H3$IntDen and loop2Hop1_H3$IntDen
  #W = 195, p-value = 0.9042
  #alternative hypothesis: true location shift is not equal to 0

lp2Pch2Hop1H3 <- subset(hop1H3, hop1H3$Label == "hop1-loop2-pch2Δ") 
wilcox.test(x=Wthop1_H3$IntDen, y = lp2Pch2Hop1H3$IntDen)
  #Wilcoxon rank sum exact test

  #data:  Wthop1_H3$IntDen and lp2Pch2Hop1H3$IntDen
  #W = 91, p-value = 0.002643
  #alternative hypothesis: true location shift is not equal to 0

pch2Hop1H3 <- subset(hop1H3, hop1H3$Label == "pch2Δ")
wilcox.test(x=Wthop1_H3$IntDen, y = pch2Hop1H3$IntDen)
  #Wilcoxon rank sum exact test

  #data:  Wthop1_H3$IntDen and pch2Hop1H3$IntDen
  #W = 57, p-value = 6.963e-06
  #alternative hypothesis: true location shift is not equal to 0


hop1H4 <- subset(Hop1_intSig, Hop1_intSig$Hour == "H4")
Wthop1_H4 <- subset(hop1H4, hop1H4$Label == "WT")
loop2Hop1_H4 <- subset(hop1H4, hop1H4$Label == "hop1-loop2")
wilcox.test(x= Wthop1_H4$IntDen, y = loop2Hop1_H4$IntDen)
  #Wilcoxon rank sum exact test

  #data:  Wthop1_H4$IntDen and loop2Hop1_H4$IntDen
  #W = 109, p-value = 0.01319
  #alternative hypothesis: true location shift is not equal to 0

lp2Pch2Hop1H4 <- subset(hop1H4, hop1H4$Label == "hop1-loop2-pch2Δ") 
wilcox.test(x=Wthop1_H4$IntDen, y = lp2Pch2Hop1H4$IntDen)
  #Wilcoxon rank sum exact test

  #data:  Wthop1_H4$IntDen and lp2Pch2Hop1H4$IntDen
  #W = 249, p-value = 0.01736
  #alternative hypothesis: true location shift is not equal to 0

pch2Hop1H4 <- subset(hop1H4, hop1H4$Label == "pch2Δ")
wilcox.test(x=Wthop1_H4$IntDen, y = pch2Hop1H4$IntDen)
  #Wilcoxon rank sum exact test

  #data:  Wthop1_H4$IntDen and pch2Hop1H4$IntDen
  #W = 11, p-value = 1.161e-08
  #alternative hypothesis: true location shift is not equal to 0





###
Zip1_intSig <- subset(integrated_sig, integrated_sig$Channel == "Zip1")
ggdensity(Zip1_intSig$IntDen)
ggqqplot(Zip1_intSig$IntDen)
shapiro.test(Zip1_intSig$IntDen)


  #Shapiro-Wilk normality test

  #data:  Zip1_intSig$IntDen
  #W = 0.96517, p-value = 4.519e-07

Zip1H2 <- subset(Zip1_intSig, Zip1_intSig$Hour == "H2")
Wtzip1_H2 <- subset(Zip1H2, Zip1H2$Label == "WT")
loop2zip1_H2 <- subset(Zip1H2, Zip1H2$Label == "hop1-loop2")
wilcox.test(x= Wtzip1_H2$IntDen, y = loop2zip1_H2$IntDen)

  #Wilcoxon rank sum exact test

  #data:  Wtzip1_H2$IntDen and loop2zip1_H2$IntDen
  #W = 305, p-value = 0.003885
  #alternative hypothesis: true location shift is not equal to 0


lp2Pch2zip1H2 <- subset(Zip1H2, Zip1H2$Label == "hop1-loop2-pch2Δ") 
wilcox.test(x=Wtzip1_H2$IntDen, y = lp2Pch2zip1H2$IntDen)
  #Wilcoxon rank sum exact test

  #data:  Wtzip1_H2$IntDen and lp2Pch2zip1H2$IntDen
  #W = 243, p-value = 0.2534
  #alternative hypothesis: true location shift is not equal to 0


pch2zip1H2 <- subset(Zip1H2, Zip1H2$Label == "pch2Δ")
wilcox.test(x=Wtzip1_H2$IntDen, y = pch2zip1H2$IntDen)
  #Wilcoxon rank sum exact test

  #data:  Wtzip1_H2$IntDen and pch2zip1H2$IntDen
  #W = 345, p-value = 0.1095
  #alternative hypothesis: true location shift is not equal to 0

Zip1H3 <- subset(Zip1_intSig, Zip1_intSig$Hour == "H3")
Wtzip1_H3 <- subset(Zip1H3, Zip1H3$Label == "WT")
loop2zip1_H3 <- subset(Zip1H3, Zip1H3$Label == "hop1-loop2")
wilcox.test(x= Wtzip1_H3$IntDen, y = loop2zip1_H3$IntDen)

  #Wilcoxon rank sum exact test

  #data:  Wtzip1_H3$IntDen and loop2zip1_H3$IntDen
  #W = 368, p-value = 6.178e-07
  #alternative hypothesis: true location shift is not equal to 0


lp2Pch2zip1H3 <- subset(Zip1H3, Zip1H3$Label == "hop1-loop2-pch2Δ") 
wilcox.test(x=Wtzip1_H3$IntDen, y = lp2Pch2zip1H3$IntDen)
  #Wilcoxon rank sum exact test

  #data:  Wtzip1_H3$IntDen and lp2Pch2zip1H3$IntDen
  #W = 169, p-value = 0.4135
  #alternative hypothesis: true location shift is not equal to 0


pch2zip1H3 <- subset(Zip1H3, Zip1H3$Label == "pch2Δ")
wilcox.test(x=Wtzip1_H3$IntDen, y = pch2zip1H3$IntDen)
  #Wilcoxon rank sum exact test

  #data:  Wtzip1_H3$IntDen and pch2zip1H3$IntDen
  #W = 314, p-value = 0.04117
  #alternative hypothesis: true location shift is not equal to 0

Zip1H4 <- subset(Zip1_intSig, Zip1_intSig$Hour == "H4")
Wtzip1_H4 <- subset(Zip1H4, Zip1H4$Label == "WT")
loop2zip1_H4 <- subset(Zip1H4, Zip1H4$Label == "hop1-loop2")
wilcox.test(x= Wtzip1_H4$IntDen, y = loop2zip1_H4$IntDen)

  #Wilcoxon rank sum exact test

  #data:  Wtzip1_H4$IntDen and loop2zip1_H4$IntDen
  #W = 367, p-value = 7.571e-07
  #alternative hypothesis: true location shift is not equal to 0



lp2Pch2zip1H4 <- subset(Zip1H4, Zip1H4$Label == "hop1-loop2-pch2Δ") 
wilcox.test(x=Wtzip1_H4$IntDen, y = lp2Pch2zip1H4$IntDen)
  #Wilcoxon rank sum exact test

  #data:  Wtzip1_H4$IntDen and lp2Pch2zip1H4$IntDen
  #W = 646, p-value = 6.039e-05
  #alternative hypothesis: true location shift is not equal to 0

pch2zip1H4 <- subset(Zip1H4, Zip1H4$Label == "pch2Δ")
wilcox.test(x=Wtzip1_H4$IntDen, y = pch2zip1H4$IntDen)
  #Wilcoxon rank sum exact test

  #data:  Wtzip1_H4$IntDen and pch2zip1H4$IntDen
  #W = 202, p-value = 0.534
  #alternative hypothesis: true location shift is not equal to 0




















Wt_hop1_intSig <- subset(Hop1_intSig, Hop1_intSig$Label == "WT")
loop2_hop1_intSig <- subset(Hop1_intSig, Hop1_intSig$Label == "hop1-loop2")
wilcox.test(x= Wt_hop1_intSig$IntDen, y = loop2_hop1_intSig$IntDen)
  #Wilcoxon rank sum test with continuity correction

  #data:  Wt_hop1_intSig$IntDen and loop2_hop1_intSig$IntDen
  #W = 2319, p-value = 0.0065
  #alternative hypothesis: true location shift is not equal to 0
Loop2pch2_hop1_intSig <- subset(Hop1_intSig, Hop1_intSig$Label == "hop1-loop2-pch2Δ")
wilcox.test(x= Wt_hop1_intSig$IntDen, y = Loop2pch2_hop1_intSig$IntDen)
  #Wilcoxon rank sum test with continuity correction

  #data:  Wt_hop1_intSig$IntDen and Loop2pch2_hop1_intSig$IntDen
  #W = 2502, p-value = 0.6691
  #alternative hypothesis: true location shift is not equal to 0
pch2_hop1_intSig <- subset(Hop1_intSig, Hop1_intSig$Label == "pch2Δ")
wilcox.test(x= Wt_hop1_intSig$IntDen, y = pch2_hop1_intSig$IntDen)
  #Wilcoxon rank sum test with continuity correction

  #data:  Wt_hop1_intSig$IntDen and pch2_hop1_intSig$IntDen
  #W = 845, p-value = 1.173e-08
  #alternative hypothesis: true location shift is not equal to 0




Zip1_intSig <- subset(integrated_sig, integrated_sig$Channel == "Zip1")
ggdensity(Zip1_intSig$IntDen)
ggqqplot(Zip1_intSig$IntDen)
shapiro.test(Zip1_intSig$IntDen)

  #Shapiro-Wilk normality test

  #data:  Zip1_intSig$IntDen
  #W = 0.96517, p-value = 4.519e-07
Wt_zip1_intSig <- subset(Zip1_intSig, Zip1_intSig$Label == "WT")
loop2_zip1_intSig <- subset(Zip1_intSig, Zip1_intSig$Label == "hop1-loop2")
wilcox.test(x= Wt_zip1_intSig$IntDen, y = loop2_zip1_intSig$IntDen)
  #Wilcoxon rank sum test with continuity correction

  #data:  Wt_zip1_intSig$IntDen and loop2_zip1_intSig$IntDen
  #W = 3178, p-value = 4.83e-13
  #alternative hypothesis: true location shift is not equal to 0
Loop2pch2_zip1_intSig <- subset(Zip1_intSig, Zip1_intSig$Label == "hop1-loop2-pch2Δ")
wilcox.test(x= Wt_zip1_intSig$IntDen, y = Loop2pch2_zip1_intSig$IntDen)
  #Wilcoxon rank sum test with continuity correction

  #data:  Wt_zip1_intSig$IntDen and Loop2pch2_zip1_intSig$IntDen
  #W = 3111, p-value = 0.002774
  #alternative hypothesis: true location shift is not equal to 0
pch2_zip1_intSig <- subset(Zip1_intSig, Zip1_intSig$Label == "pch2Δ")
wilcox.test(x= Wt_zip1_intSig$IntDen, y = pch2_zip1_intSig$IntDen)
  #Wilcoxon rank sum test with continuity correction

  #data:  Wt_zip1_intSig$IntDen and pch2_zip1_intSig$IntDen
  #W = 2561, p-value = 0.01294
  #alternative hypothesis: true location shift is not equal to 0








###
all_zip1_length <- read.csv('/Users/darmokandjalad/Documents/Hop1PHD-Loop2/Manuscript/Hop1_Loop2_figures/Figure2/pch2 + Zip1TC_loop2/EditedZip1Length-allstrains.csv')
ggdensity(all_zip1_length$Total.Branch.Length)
ggqqplot(all_zip1_length$Total.Branch.Length)
shapiro.test(all_zip1_length$Total.Branch.Length)
AllTotalLeng <- all_zip1_length$Total.Branch.Length
uMAllTotalLeng <- AllTotalLeng*0.0799
all_zip1_length <- cbind(all_zip1_length, uMAllTotalLeng)

Wt_zip1_leng <- subset(all_zip1_length, all_zip1_length$Strain == "7797")
H2WtLeng <- subset(Wt_zip1_leng, Wt_zip1_leng$Hour == "H2")
lop2Zip2_leng <- subset(all_zip1_length, all_zip1_length$Strain == "11644")
H2lp2_zip1Leng <- subset(lop2Zip2_leng, lop2Zip2_leng$Hour == "H2")
wilcox.test(x=H2WtLeng$Total.Branch.Length, y= H2lp2_zip1Leng$Total.Branch.Length)
  #Wilcoxon rank sum test with continuity correction

  #data:  H2WtLeng$Total.Branch.Length and H2lp2_zip1Leng$Total.Branch.Length
  #W = 154.5, p-value = 0.2235
  #alternative hypothesis: true location shift is not equal to 0

pch2_zip1_leng <- subset(all_zip1_length, all_zip1_length$Strain == "11758")
H2pch2_zip1leng <- subset(pch2_zip1_leng, pch2_zip1_leng$Hour == "H2")
wilcox.test(x=H2WtLeng$Total.Branch.Length, y= H2pch2_zip1leng$Total.Branch.Length)
  #Wilcoxon rank sum exact test

  #data:  H2WtLeng$Total.Branch.Length and H2pch2_zip1leng$Total.Branch.Length
  #W = 80, p-value = 1.595e-05
  #alternative hypothesis: true location shift is not equal to 0

loop2pch2_zip1_leng <- subset(all_zip1_length, all_zip1_length$Strain == "11757")
H2lp2Pch2_zip1_leng <- subset(loop2pch2_zip1_leng, loop2pch2_zip1_leng$Hour == "H2")
wilcox.test(x=H2WtLeng$Total.Branch.Length, y= H2lp2Pch2_zip1_leng$Total.Branch.Length)
  #Wilcoxon rank sum exact test

  #data:  H2WtLeng$Total.Branch.Length and H2lp2Pch2_zip1_leng$Total.Branch.Length
  #W = 174, p-value = 0.4945
  #alternative hypothesis: true location shift is not equal to 0


H3WtLeng <- subset(Wt_zip1_leng, Wt_zip1_leng$Hour == "H3")
H3lp2_zip1Leng <- subset(lop2Zip2_leng, lop2Zip2_leng$Hour == "H3")
wilcox.test(x=H3WtLeng$Total.Branch.Length, y= H3lp2_zip1Leng$Total.Branch.Length)
  #Wilcoxon rank sum test with continuity correction

  #data:  H3WtLeng$Total.Branch.Length and H3lp2_zip1Leng$Total.Branch.Length
  #W = 346, p-value = 8.286e-05
  #alternative hypothesis: true location shift is not equal to 0

H3pch2_zip1leng <- subset(pch2_zip1_leng, pch2_zip1_leng$Hour == "H3")
wilcox.test(x=H3WtLeng$Total.Branch.Length, y= H3pch2_zip1leng$Total.Branch.Length)
  #Wilcoxon rank sum exact test

  #data:  H3WtLeng$Total.Branch.Length and H3pch2_zip1leng$Total.Branch.Length
  #W = 130, p-value = 0.0143
  #alternative hypothesis: true location shift is not equal to 0

H3lp2Pch2_zip1_leng <- subset(loop2pch2_zip1_leng, loop2pch2_zip1_leng$Hour == "H3")
wilcox.test(x=H3WtLeng$Total.Branch.Length, y= H3lp2Pch2_zip1_leng$Total.Branch.Length)
  #Wilcoxon rank sum exact test

  #data:  H3WtLeng$Total.Branch.Length and H3lp2Pch2_zip1_leng$Total.Branch.Length
  #W = 261, p-value = 0.1022
  #alternative hypothesis: true location shift is not equal to 0


H4WtLeng <- subset(Wt_zip1_leng, Wt_zip1_leng$Hour == "H4")
H4lp2_zip1Leng <- subset(lop2Zip2_leng, lop2Zip2_leng$Hour == "H4")
wilcox.test(x=H4WtLeng$Total.Branch.Length, y= H4lp2_zip1Leng$Total.Branch.Length)

  #Wilcoxon rank sum exact test

  #data:  H4WtLeng$Total.Branch.Length and H4lp2_zip1Leng$Total.Branch.Length
  #W = 395, p-value = 2.757e-10
  #alternative hypothesis: true location shift is not equal to 0
H4pch2_zip1leng <- subset(pch2_zip1_leng, pch2_zip1_leng$Hour == "H4")
wilcox.test(x=H4WtLeng$Total.Branch.Length, y= H4pch2_zip1leng$Total.Branch.Length)
  #Wilcoxon rank sum exact test

  #data:  H4WtLeng$Total.Branch.Length and H4pch2_zip1leng$Total.Branch.Length
  #W = 202, p-value = 0.534
  #alternative hypothesis: true location shift is not equal to 0
H4lp2Pch2_zip1_leng <- subset(loop2pch2_zip1_leng, loop2pch2_zip1_leng$Hour == "H4")
wilcox.test(x=H4WtLeng$Total.Branch.Length, y= H4lp2Pch2_zip1_leng$Total.Branch.Length)
  #Wilcoxon rank sum test with continuity correction

  #data:  H4WtLeng$Total.Branch.Length and H4lp2Pch2_zip1_leng$Total.Branch.Length
  #W = 583, p-value = 0.004212
  #alternative hypothesis: true location shift is not equal to 0



















Wt_zip1_leng <- subset(all_zip1_length, all_zip1_length$Strain == "7797")
loop2_zip1_leng <- subset(all_zip1_length, all_zip1_length$Strain == "11644")
wilcox.test(x=Wt_zip1_leng$Total.Branch.Length, y=loop2_zip1_leng$Total.Branch.Length)

  #Wilcoxon rank sum test with continuity correction

  #data:  Wt_zip1_leng$Total.Branch.Length and loop2_zip1_leng$Total.Branch.Length
  #W = 2832.5, p-value = 6.073e-08
  #alternative hypothesis: true location shift is not equal to 0

pch2_zip1_leng <- subset(all_zip1_length, all_zip1_length$Strain == "11758")
wilcox.test(x=Wt_zip1_leng$Total.Branch.Length, y=pch2_zip1_leng$Total.Branch.Length)
  #Wilcoxon rank sum test with continuity correction

  #data:  Wt_zip1_leng$Total.Branch.Length and pch2_zip1_leng$Total.Branch.Length
  #W = 1427.5, p-value = 0.003475
  #alternative hypothesis: true location shift is not equal to 0

loop2pch2_zip1_leng <- subset(all_zip1_length, all_zip1_length$Strain == "11757")
wilcox.test(x=Wt_zip1_leng$Total.Branch.Length, y=loop2pch2_zip1_leng$Total.Branch.Length)
  #Wilcoxon rank sum test with continuity correction

  #data:  Wt_zip1_leng$Total.Branch.Length and loop2pch2_zip1_leng$Total.Branch.Length
  #W = 2698, p-value = 0.2103
  #alternative hypothesis: true location shift is not equal to 0





hop1FociArea <- read.csv('/Users/darmokandjalad/Documents/Hop1PHD-Loop2/IF_results/IF_Analysis_CSV/Combined_DapiAreaHop1Area.csv')

ggdensity(hop1FociArea$Hop1.Total.Area)
ggqqplot(hop1FociArea$Hop1.Total.Area)
shapiro.test(hop1FociArea$Hop1.Total.Area)

wtHop1Area <- subset(hop1FociArea, hop1FociArea$Strain =="WT")
Wt_2_hop1Area <- subset(wtHop1Area, wtHop1Area$Hour == "H2")
loop2Hop1Area <- subset(hop1FociArea, hop1FociArea$Strain == "hop1-loop2")
lp_2_hop1Are  <- subset(loop2Hop1Area, loop2Hop1Area$Hour == "H2")
wilcox.test(x=wtHop1Area$Percent.Area.Hop1OverDapi, y=loop2Hop1Area$Percent.Area.Hop1OverDapi)
  #Wilcoxon rank sum test with continuity correction

  #data:  wtHop1Area$Percent.Area.Hop1OverDapi and loop2Hop1Area$Percent.Area.Hop1OverDapi
  #W = 1990, p-value = 0.3199
  #alternative hypothesis: true location shift is not equal to 0
wilcox.test(x=Wt_2_hop1Area$Percent.Area.Hop1OverDapi, y=lp_2_hop1Are$Percent.Area.Hop1OverDapi)
  #Wilcoxon rank sum exact test

  #data:  Wt_2_hop1Area$Percent.Area.Hop1OverDapi and lp_2_hop1Are$Percent.Area.Hop1OverDapi
  #W = 383, p-value = 1.758e-08
  #alternative hypothesis: true location shift is not equal to 0


pch2Hop1Area <- subset(hop1FociArea, hop1FociArea$Strain =="pch2∆")
wilcox.test(x=wtHop1Area$Percent.Area.Hop1OverDapi, y=pch2Hop1Area$Percent.Area.Hop1OverDapi)
  #Wilcoxon rank sum test with continuity correction

  #data:  wtHop1Area$Percent.Area.Hop1OverDapi and pch2Hop1Area$Percent.Area.Hop1OverDapi
  #W = 30, p-value < 2.2e-16
  #alternative hypothesis: true location shift is not equal to 0
pch_2_hopArea <- subset(pch2Hop1Area, pch2Hop1Area$Hour == "H2")
wilcox.test(x=Wt_2_hop1Area$Percent.Area.Hop1OverDapi, y=pch_2_hopArea$Percent.Area.Hop1OverDapi)
  #Wilcoxon rank sum exact test

  #data:  Wt_2_hop1Area$Percent.Area.Hop1OverDapi and pch_2_hopArea$Percent.Area.Hop1OverDapi
  #W = 17, p-value = 2.483e-10
  #alternative hypothesis: true location shift is not equal to 0



pchloopArea <- subset(hop1FociArea, hop1FociArea$Strain =="hop1-loop2-pch2∆")
wilcox.test(x=wtHop1Area$Percent.Area.Hop1OverDapi, y=pchloopArea$Percent.Area.Hop1OverDapi)

  #Wilcoxon rank sum test with continuity correction

  #data:  wtHop1Area$Percent.Area.Hop1OverDapi and pchloopArea$Percent.Area.Hop1OverDapi
  #W = 2143, p-value = 0.2801
  #alternative hypothesis: true location shift is not equal to 0
pchloop_2_Area <- subset(pchloopArea, pchloopArea$Hour =="H2")
wilcox.test(x=Wt_2_hop1Area$Percent.Area.Hop1OverDapi, y=pchloop_2_Area$Percent.Area.Hop1OverDapi)

  #Wilcoxon rank sum exact test

  #data:  Wt_2_hop1Area$Percent.Area.Hop1OverDapi and pchloop_2_Area$Percent.Area.Hop1OverDapi
  #W = 360, p-value = 2.884e-06
  #alternative hypothesis: true location shift is not equal to 0


Wt_3_hop1Area <- subset(wtHop1Area, wtHop1Area$Hour == "H3")
lp_3_hop1Are  <- subset(loop2Hop1Area, loop2Hop1Area$Hour == "H3")
wilcox.test(x=Wt_3_hop1Area$Percent.Area.Hop1OverDapi, y=lp_3_hop1Are$Percent.Area.Hop1OverDapi)
  #Wilcoxon rank sum exact test

  #data:  Wt_3_hop1Area$Percent.Area.Hop1OverDapi and lp_3_hop1Are$Percent.Area.Hop1OverDapi
  #W = 127, p-value = 0.04909
  #alternative hypothesis: true location shift is not equal to 0

pch_3_hopArea <- subset(pch2Hop1Area, pch2Hop1Area$Hour == "H3")
wilcox.test(x=Wt_3_hop1Area$Percent.Area.Hop1OverDapi, y=pch_3_hopArea$Percent.Area.Hop1OverDapi)

  #Wilcoxon rank sum exact test

  #data:  Wt_3_hop1Area$Percent.Area.Hop1OverDapi and pch_3_hopArea$Percent.Area.Hop1OverDapi
  #W = 1, p-value = 4.164e-12
  #alternative hypothesis: true location shift is not equal to 0

pchloop_3_Area <- subset(pchloopArea, pchloopArea$Hour =="H3")
wilcox.test(x=Wt_3_hop1Area$Percent.Area.Hop1OverDapi, y=pchloop_3_Area$Percent.Area.Hop1OverDapi)
  #Wilcoxon rank sum exact test

  #data:  Wt_3_hop1Area$Percent.Area.Hop1OverDapi and pchloop_3_Area$Percent.Area.Hop1OverDapi
  #W = 143, p-value = 0.1274
  #alternative hypothesis: true location shift is not equal to 0
Wt_4_hop1Area <- subset(wtHop1Area, wtHop1Area$Hour == "H4")
lp_4_hop1Are  <- subset(loop2Hop1Area, loop2Hop1Area$Hour == "H4")
wilcox.test(x=Wt_4_hop1Area$Percent.Area.Hop1OverDapi, y=lp_4_hop1Are$Percent.Area.Hop1OverDapi)

  #Wilcoxon rank sum exact test

  #data:  Wt_4_hop1Area$Percent.Area.Hop1OverDapi and lp_4_hop1Are$Percent.Area.Hop1OverDapi
  #W = 76, p-value = 0.0005305
  #alternative hypothesis: true location shift is not equal to 0
pch_4_hopArea <- subset(pch2Hop1Area, pch2Hop1Area$Hour == "H4")
wilcox.test(x=Wt_4_hop1Area$Percent.Area.Hop1OverDapi, y=pch_4_hopArea$Percent.Area.Hop1OverDapi)

  #Wilcoxon rank sum exact test

  #data:  Wt_4_hop1Area$Percent.Area.Hop1OverDapi and pch_4_hopArea$Percent.Area.Hop1OverDapi
  #   W = 0, p-value = 5.956e-11
  #alternative hypothesis: true location shift is not equal to 0
pchloop_4_Area <- subset(pchloopArea, pchloopArea$Hour =="H4")
wilcox.test(x=Wt_4_hop1Area$Percent.Area.Hop1OverDapi, y=pchloop_4_Area$Percent.Area.Hop1OverDapi)

  #Wilcoxon rank sum exact test

  #data:  Wt_4_hop1Area$Percent.Area.Hop1OverDapi and pchloop_4_Area$Percent.Area.Hop1OverDapi
  #W = 109, p-value = 1.102e-06
  #alternative hypothesis: true location shift is not equal to 0





