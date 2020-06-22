files <- list.files(pattern=".txt")
# wdth(20); colnames(read.table(files[1], sep="\t", header=TRUE, as.is=TRUE, nrow=3)); wdth(100)
#  [1] "Protein.Name"
#  [2] "Peptide.Modified.Sequence"
#  [3] "PeptideLabel"
#  [4] "avg"
#  [5] "LCL57_ATMi_0Gy_1"
#  [6] "LCL57_ATMi_5Gy_0.25"
#  [7] "LCL57_ATMi_5Gy_1"
#  [8] "LCL57_ATMi_5Gy_6"
#  [9] "LCL57_ATMi_5Gy_24"
# [10] "avg.1"
# [11] "LCL57_DMSO_0Gy_1"
# [12] "LCL57_DMSO_5Gy_0.25"
# [13] "LCL57_DMSO_5Gy_1"
# [14] "LCL57_DMSO_5Gy_6"
# [15] "LCL57_DMSO_5Gy_24"
# [16] "stdev"
# [17] "LCL57_ATMi_0Gy_1.1"
# [18] "LCL57_ATMi_5Gy_0.25.1"
# [19] "LCL57_ATMi_5Gy_1.1"
# [20] "LCL57_ATMi_5Gy_6.1"
# [21] "LCL57_ATMi_5Gy_24.1"
# [22] "stdev.1"
# [23] "LCL57_DMSO_0Gy_1.1"
# [24] "LCL57_DMSO_5Gy_0.25.1"
# [25] "LCL57_DMSO_5Gy_1.1"
# [26] "LCL57_DMSO_5Gy_6.1"
# [27] "LCL57_DMSO_5Gy_24.1"

tbl.orig <- read.table(files[1], sep="\t", header=TRUE, as.is=TRUE, nrow=-1)
dim(tbl.orig)  #37 27
save(tbl.orig, file="tbl.orig.RData")

tbl.atmi.avg <- tbl.orig[, 5:9]
tbl.atmi.sd  <- tbl.orig[, 17:21]

tbl.dmso.avg <- tbl.orig[, 11:15]
tbl.dmso.sd  <- tbl.orig[, 23:27]
rownames(tbl.atmi.avg) <- tbl.orig[,3]
rownames(tbl.dmso.avg) <- tbl.orig[,3]

rownames(tbl.atmi.sd) <- tbl.orig[,3]
rownames(tbl.dmso.sd) <- tbl.orig[,3]

save(tbl.atmi.avg, tbl.atmi.sd, tbl.dmso.avg, tbl.dmso.sd, file="LCL57_1.RData")

# tbl.avg <- cbind(tbl.atmi.avg, tbl.dmso.avg)
# tbl.sd <- cbind(tbl.atmi.sd, tbl.dmso.sd)
# colnames(tbl.sd) <- colnames(tbl.avg)
# save(tbl.avg, tbl.sd, file="LCL57_1.RData")

library(ggplot2)
time <- c(0, 0.25, 1, 6, 24)
atmi <- as.numeric(tbl.atmi.avg[1,])
atmi.sd <-    as.numeric(tbl.atmi.sd[1,])
dmso <- as.numeric(tbl.dmso.avg[1,])
dmso.sd <-  as.numeric(tbl.dmso.sd[1,])
tbl.tmp <- data.frame(time=as.factor(time), atmi=atmi, atmi.sd=atmi.sd, dmso=dmso, dmso.sd=dmso.sd)


avg <- with(tbl.tmp, c(atmi, dmso))
time <- as.numeric(rep(tbl.tmp$time,2))
sd <- with(tbl.tmp, c(atmi.sd, dmso.sd))
experiment <- c(rep("ATMi", 5), rep("DMSO", 5))

length(experiment)
length(time)

tbl <- data.frame(experiment=experiment, time=time, area=avg, sd=sd)
lapply(tbl, class)

p <- ggplot(tbl) +
     geom_line(aes(y = area, x = time, colour = experiment), stat="identity", size=1.5) +
     geom_errorbar(aes(ymin=area-sd, ymax=area+sd, x=time, y=area, colour=experiment), width=0.2) +
     ggtitle("Time Course for individual analytes") +
     theme(plot.title = element_text(size = 18, face = "bold")) +
     labs(x="Time", y="Ratio") +
     scale_x_continuous(breaks=1:5,
                  labels=c("Mock", "0.25", "1.0", "6", "24"))

p


#----------------------------------------------------------------------------------------------------
toPlotFriendlyTable <- function(tbl.in)
{
   analyte <- tbl.in$PeptideLabel
   tbl.atmi.avg <- tbl.in[, 5:9]
   tbl.atmi.sd  <- tbl.in[, 17:21]

   tbl.dmso.avg <- tbl.in[, 11:15]
   tbl.dmso.sd  <- tbl.in[, 23:27]

   time <- c(0, 0.25, 1, 6, 24)
   atmi <- as.numeric(tbl.atmi.avg[1,])
   atmi.sd <-    as.numeric(tbl.atmi.sd[1,])
   dmso <- as.numeric(tbl.dmso.avg[1,])
   dmso.sd <-  as.numeric(tbl.dmso.sd[1,])
   tbl.tmp <- data.frame(time=as.factor(time), atmi=atmi, atmi.sd=atmi.sd, dmso=dmso, dmso.sd=dmso.sd)


   avg <- with(tbl.tmp, c(atmi, dmso))
   time <- as.numeric(rep(tbl.tmp$time,2))
   sd <- with(tbl.tmp, c(atmi.sd, dmso.sd))
   experiment <- c(rep("ATMi", 5), rep("DMSO", 5))

   length(experiment)
   length(time)

   data.frame(analyte=rep(analyte, length(experiment)), experiment=experiment, time=time, area=avg, sd=sd)

} # toPlotFriendlyTable
#----------------------------------------------------------------------------------------------------
test_toPlotFriendlyTable <- function()
{
   if(!exists("tbl.orig"))
     tbl.orig <- get(load("tbl.orig.RData"))

   tbl.good <- toPlotFriendlyTable(tbl.orig[1,, drop=FALSE])

} # test_toPlotFriendlyTable
#----------------------------------------------------------------------------------------------------
plotAnalyte <- function(tbl.analyte)
{
   analyteName <- tbl.analyte$analyte[1]

   p <- ggplot(tbl.analyte) +
         geom_line(aes(y = area, x = time, colour = experiment), stat="identity", size=1.5) +
         geom_errorbar(aes(ymin=area-sd, ymax=area+sd, x=time, y=area, colour=experiment), width=0.2) +
         ggtitle(sprintf("%s Time Course", analyteName)) +
         theme(plot.title = element_text(size = 18, face = "bold")) +
         labs(x="Time", y="Ratio") +
         scale_x_continuous(breaks=1:5, labels=c("Mock", "0.25", "1.0", "6", "24"))

   p

} # plotAnalyte
#----------------------------------------------------------------------------------------------------
