rm(list = ls())
library(dplyr)
setwd("~/Desktop/Stuff/IremResults_OverlapComparison/Coverages")


operonstart <- 2120000
operonend <- 2121500
operonsize <- operonend - operonstart + 1
#########################################################################################

                                          #CHLORAMPHENICOL#

##########################################################################################
ChlCoverageNDCT0 <- read.delim(file = "~/Desktop/Stuff/IremResults_OverlapComparison/Coverages/Chlndc_T0.coverage", header = FALSE)
colnames(ChlCoverageNDCT0) <- c("chromosome", "locus", "NDCT0.Cov")

RPM_scale <- sum(ChlCoverageNDCT0$NDCT0.Cov)/1000000

x <- ChlCoverageNDCT0$NDCT0.Cov

for(i in 1:length(x)){
  x[i] <- x[i]/RPM_scale
}
ChlCoverageNDCT0$RPM.NDCT0.Cov <- x
##############################################
ChlCoverageNDCT60 <- read.delim(file = "~/Desktop/Stuff/IremResults_OverlapComparison/Coverages/Chlndc_T60.coverage", header = FALSE)
colnames(ChlCoverageNDCT60) <- c("chromosome", "Locus", "NDCT60.Cov")

RPM_scale <- sum(ChlCoverageNDCT60$NDCT60.Cov)/1000000

x <- ChlCoverageNDCT60$NDCT60.Cov

for(i in 1:length(x)){
  x[i] <- x[i]/RPM_scale
}
ChlCoverageNDCT60$RPM.NDCT60.Cov <- x
#############################################

ChlCoverage1Q <- read.delim(file = "~/Desktop/Stuff/IremResults_OverlapComparison/Coverages/Chl1Q_T60.coverage", header = FALSE)
colnames(ChlCoverage1Q) <- c("chromosome", "Locus", "1Q.Cov")

RPM_scale <- sum(ChlCoverage1Q$`1Q.Cov`)/1000000

x <- ChlCoverage1Q$`1Q.Cov`

for(i in 1:length(x)){
  x[i] <- x[i]/RPM_scale
}
ChlCoverage1Q$RPM.1QT60.Cov <- x
#############################################

ChlCoverage3Q <- read.delim(file = "~/Desktop/Stuff/IremResults_OverlapComparison/Coverages/Chl3Q_T60.coverage", header = FALSE)
colnames(ChlCoverage3Q) <- c("Chromosome", "Locus", "3Q.Cov")

RPM_scale <- sum(ChlCoverage3Q$`3Q.Cov`)/1000000

x <- ChlCoverage3Q$`3Q.Cov`

for(i in 1:length(x)){
  x[i] <- x[i]/RPM_scale
}
ChlCoverage3Q$RPM.3QT60.Cov <- x

#################################################################################################

                          ################  Chl Graph  #############

#################################################################################################

SP2199Operon <- as.data.frame(matrix(nrow= operonsize, ncol=5))

colnames(SP2199Operon) <- c("Locus", "RPM.NDC.T0", "RPM.NDC.T60", "RPM.1Q", "RPM.3Q")
SP2199Operon$Locus <- ChlCoverageNDCT0[operonstart:operonend, 2]/1000
SP2199Operon$RPM.NDC.T0 <- ChlCoverageNDCT0[operonstart:operonend, 4]
SP2199Operon$RPM.NDC.T60 <- ChlCoverageNDCT60[operonstart:operonend, 4]
SP2199Operon$RPM.1Q <- ChlCoverage1Q[operonstart:operonend, 4]
SP2199Operon$RPM.3Q <- ChlCoverage3Q[operonstart:operonend, 4]


SP2199Operon.log10 <- as.data.frame(matrix(nrow= operonsize, ncol=5))
colnames(SP2199Operon.log10) <- colnames(SP2199Operon)
SP2199Operon.log10$Locus <- SP2199Operon$Locus
SP2199Operon.log10$RPM.NDC.T0 <- log10(SP2199Operon$RPM.NDC.T0)
SP2199Operon.log10$RPM.NDC.T60 <- log10(SP2199Operon$RPM.NDC.T60)
SP2199Operon.log10$RPM.1Q <- log10(SP2199Operon$RPM.1Q)
SP2199Operon.log10$RPM.3Q <- log10(SP2199Operon$RPM.3Q)


plot.ts(SP2199Operon.log10)

plot(SP2199Operon.log10$Locus, SP2199Operon.log10$RPM.NDC.T60, type = "l", col = "blue")
lines(SP2199Operon.log10$Locus, SP2199Operon.log10$RPM.1Q, type = "l", col = "yellow")
lines(SP2199Operon.log10$Locus, SP2199Operon.log10$RPM.3Q, type = "l", col = "red")
write.csv(SP2199Operon, file = "SP2199Operon_Chl_RNAseq_Coverage.csv", quote = FALSE, sep = "", row.names = FALSE, col.names = FALSE)

#########################################################################################

                      #KASUGAMYCIN#

##########################################################################################

KsgCoverageNDCT0 <- read.delim(file = "~/Desktop/Stuff/IremResults_OverlapComparison/Coverages/Ksgndc_T0.coverage", header = FALSE)
colnames(KsgCoverageNDCT0) <- c("chromosome", "locus", "NDCT0.Cov")

RPM_scale <- sum(KsgCoverageNDCT0$NDCT0.Cov)/1000000

x <- KsgCoverageNDCT0$NDCT0.Cov

for(i in 1:length(x)){
  x[i] <- x[i]/RPM_scale
}
KsgCoverageNDCT0$RPM.NDCT0.Cov <- x
##############################################
KsgCoverageNDCT60 <- read.delim(file = "~/Desktop/Stuff/IremResults_OverlapComparison/Coverages/Ksgndc_T60.coverage", header = FALSE)
colnames(KsgCoverageNDCT60) <- c("chromosome", "Locus", "NDCT60.Cov")

RPM_scale <- sum(KsgCoverageNDCT60$NDCT60.Cov)/1000000

x <- KsgCoverageNDCT60$NDCT60.Cov

for(i in 1:length(x)){
  x[i] <- x[i]/RPM_scale
}
KsgCoverageNDCT60$RPM.NDCT60.Cov <- x
#############################################

KsgCoverage1Q <- read.delim(file = "~/Desktop/Stuff/IremResults_OverlapComparison/Coverages/Ksg1Q_T60.coverage", header = FALSE)
colnames(KsgCoverage1Q) <- c("chromosome", "Locus", "1Q.Cov")

RPM_scale <- sum(KsgCoverage1Q$`1Q.Cov`)/1000000

x <- KsgCoverage1Q$`1Q.Cov`

for(i in 1:length(x)){
  x[i] <- x[i]/RPM_scale
}
KsgCoverage1Q$RPM.1QT60.Cov <- x
#############################################

KsgCoverage3Q <- read.delim(file = "~/Desktop/Stuff/IremResults_OverlapComparison/Coverages/Ksg3Q_T60.coverage", header = FALSE)
colnames(KsgCoverage3Q) <- c("Chromosome", "Locus", "3Q.Cov")

RPM_scale <- sum(KsgCoverage3Q$`3Q.Cov`)/1000000

x <- KsgCoverage3Q$`3Q.Cov`

for(i in 1:length(x)){
  x[i] <- x[i]/RPM_scale
}
KsgCoverage3Q$RPM.3QT60.Cov <- x

#################################################################################################

                          ################  Ksg Graph  #############

#################################################################################################

SP2199Operon <- as.data.frame(matrix(nrow= operonsize, ncol=5))

colnames(SP2199Operon) <- c("Locus", "RPM.NDC.T0", "RPM.NDC.T60", "RPM.1Q", "RPM.3Q")
SP2199Operon$Locus <- KsgCoverageNDCT0[operonstart:operonend, 2]/1000
SP2199Operon$RPM.NDC.T0 <- KsgCoverageNDCT0[operonstart:operonend, 4]
SP2199Operon$RPM.NDC.T60 <- KsgCoverageNDCT60[operonstart:operonend, 4]
SP2199Operon$RPM.1Q <- KsgCoverage1Q[operonstart:operonend, 4]
SP2199Operon$RPM.3Q <- KsgCoverage3Q[operonstart:operonend, 4]


SP2199Operon.log10 <- as.data.frame(matrix(nrow= operonsize, ncol=5))
colnames(SP2199Operon.log10) <- colnames(SP2199Operon)
SP2199Operon.log10$Locus <- SP2199Operon$Locus
SP2199Operon.log10$RPM.NDC.T0 <- log10(SP2199Operon$RPM.NDC.T0)
SP2199Operon.log10$RPM.NDC.T60 <- log10(SP2199Operon$RPM.NDC.T60)
SP2199Operon.log10$RPM.1Q <- log10(SP2199Operon$RPM.1Q)
SP2199Operon.log10$RPM.3Q <- log10(SP2199Operon$RPM.3Q)


plot.ts(SP2199Operon.log10)

plot(SP2199Operon.log10$Locus, SP2199Operon.log10$RPM.NDC.T0, type = "l", col = "green")
lines(SP2199Operon.log10$Locus, SP2199Operon.log10$RPM.NDC.T60, type = "l", col = "blue")
lines(SP2199Operon.log10$Locus, SP2199Operon.log10$RPM.1Q, type = "l", col = "yellow")
lines(SP2199Operon.log10$Locus, SP2199Operon.log10$RPM.3Q, type = "l", col = "red")
write.csv(SP2199Operon, file = "SP2199Operon_Ksg_RNAseq_Coverage.csv", quote = FALSE, sep = "", row.names = FALSE, col.names = FALSE)
