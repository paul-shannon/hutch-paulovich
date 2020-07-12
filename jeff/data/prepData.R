library(ggplot2)
library(RUnit)
#----------------------------------------------------------------------------------------------------
files <- list.files(pattern=".txt")
#----------------------------------------------------------------------------------------------------
# LCL - Lymphoblastoid Cell lines: a Continuous in Vitro Source of
# Cells to Study Carcinogen Sensitivity and DNA Repair
#----------------------------------------------------------------------------------------------------
runTests <- function ()
{
  test_toPlotFriendlyTable.01()
  test_toPlotFriendlyTable.02()
  test_toPlotFriendlyTable.03()

} # runTests
#----------------------------------------------------------------------------------------------------
toPlotFriendlyTable.01 <- function(tbl.in)
{
   analyte <- tbl.in$PeptideLabel
   tbl.atmi.avg <- tbl.in[, 5:9]
   tbl.atmi.sd  <- tbl.in[, 17:21]

   tbl.dmso.avg <- tbl.in[, 11:15]
   tbl.dmso.sd  <- tbl.in[, 23:27]

   time <- c(0, 0.25, 1, 6, 24)
   radiation <- c(0, 5, 5, 5, 5)
   atmi <- as.numeric(tbl.atmi.avg[1,])
   atmi.sd <-    as.numeric(tbl.atmi.sd[1,])
   dmso <- as.numeric(tbl.dmso.avg[1,])
   dmso.sd <-  as.numeric(tbl.dmso.sd[1,])

   tbl.tmp <- data.frame(time=as.factor(time), radiation=radiation, atmi=atmi, atmi.sd=atmi.sd, dmso=dmso, dmso.sd=dmso.sd)


   avg <- with(tbl.tmp, c(atmi, dmso))
   time <- as.numeric(rep(tbl.tmp$time,2))
   sd <- with(tbl.tmp, c(atmi.sd, dmso.sd))
   experiment <- c(rep("ATMi", 5), rep("DMSO", 5))
   radiation <- as.numeric(rep(tbl.tmp$radiation, 2))
   length(experiment)
   length(time)
   length(radiation)

   tbl.out <- data.frame(analyte=rep(analyte, length(experiment)),
                         experiment=as.character(experiment),
                         time=time,
                         radiation=radiation,
                         area=avg,
                         sd=sd)
   tbl.out$analyte <- as.character(tbl.out$analyte)
   tbl.out$experiment <- as.character(tbl.out$experiment)

   tbl.out

} # toPlotFriendlyTable.01
#----------------------------------------------------------------------------------------------------
test_toPlotFriendlyTable.01 <- function()
{
   printf("--- test_toPlotFriendlyTable.01")

   tbl.orig <- read.table(files[1], sep="\t", header=TRUE, as.is=TRUE, nrow=-1)

   tbl.good <- toPlotFriendlyTable.01(tbl.orig[1,, drop=FALSE])
   checkEquals(colnames(tbl.good), c("analyte", "experiment", "time", "radiation", "area", "sd"))
   checkEquals(dim(tbl.good), c(10, 6))
   checkEquals(sort(unique(tbl.good$experiment)), c("ATMi", "DMSO"))
   checkTrue(all(tbl.good$analyte %in% all.peptides))

   tbls <- lapply(1:nrow(tbl.orig), function(row) toPlotFriendlyTable.01(tbl.orig[row,, drop=FALSE]))
   checkEquals(length(tbls), 37)
   tbl.01 <- do.call(rbind, tbls)
   checkEquals(dim(tbl.01), c(370, 6))
   checkEquals(length(unique(tbl.01$analyte)), 37)
   checkTrue(all(tbl.01$analyte %in% all.peptides))

   save(tbl.01, file="tbl.01.RData")

} # test_toPlotFriendlyTable.01
#----------------------------------------------------------------------------------------------------
toPlotFriendlyTable.02 <- function(tbl.in)
{
   analyte <- tbl.in$Peptide.Label
   tbl.ir.avg <- as.numeric(tbl.in[, 5:9])
   tbl.ir.sd  <- as.numeric(tbl.in[, 11:15])

   time <- c(0, 0.25, 1, 6, 24)
   radiation <- c(0, 1, 2, 5, 10)
   #atmi <- as.numeric(tbl.atmi.avg[1,])
   #atmi.sd <-    as.numeric(tbl.atmi.sd[1,])
   #dmso <- as.numeric(tbl.dmso.avg[1,])
   #dmso.sd <-  as.numeric(tbl.dmso.sd[1,])
   tbl.out <- data.frame(analyte=rep(analyte, length(time)),
                         experiment=rep("ir", length(time)),
                         time=as.factor(time),
                         radiation=radiation,
                         area=tbl.ir.avg,
                         sd=tbl.ir.sd,
                         stringsAsFactors=FALSE)

   tbl.out$time <- as.factor(tbl.out$time)
   tbl.out

} # toPlotFriendlyTable.02
#----------------------------------------------------------------------------------------------------
test_toPlotFriendlyTable.02 <- function()
{
   printf("--- test_toPlotFriendlyTable.02")

   tbl.orig <- read.table(files[2], sep="\t", header=TRUE, as.is=TRUE, nrow=-1)

   tbl.good <- toPlotFriendlyTable.02(tbl.orig[1,, drop=FALSE])
   checkEquals(colnames(tbl.good), c("analyte", "experiment", "time", "radiation", "area", "sd"))
   checkEquals(dim(tbl.good), c(5, 6))
   checkEquals(sort(unique(tbl.good$experiment)), "ir")
   checkTrue(all(tbl.good$analyte %in% all.peptides))

   tbls <- lapply(1:nrow(tbl.orig), function(row) toPlotFriendlyTable.02(tbl.orig[row,, drop=FALSE]))
   checkEquals(length(tbls), 39)
   tbl.02 <- do.call(rbind, tbls)
   checkEquals(dim(tbl.02), c(195, 6))
   checkEquals(length(unique(tbl.02$analyte)), 39)
   checkTrue(all(tbl.02$analyte %in% all.peptides))

   save(tbl.02, file="tbl.02.RData")

} # test_toPlotFriendlyTable.02
#----------------------------------------------------------------------------------------------------
toPlotFriendlyTable.03 <- function(tbl.in)
{
   analyte <- tbl.in$PeptideLabel
      #----------------------------------------
      # PBMC + and - SEB_DMSO
      #----------------------------------------
      # grep("minus.SEB_DMSO", colnames(tbl.in))  # avg: 13 14 sd: 30 31
      # grep("plus.SEB_DMSO", colnames(tbl.in))  # avg: 5 6  sd: 22 23

   tbl.plus.seb.dmso.avg <- as.numeric(tbl.in[, c(5, 6)])
   tbl.plus.seb.dmso.sd <- as.numeric(tbl.in[, c(22, 23)])

   tbl.minus.seb.dmso.avg <- as.numeric(tbl.in[, c(13, 14)])
   tbl.minus.seb.dmso.sd <- as.numeric(tbl.in[, c(30, 31)])


      #----------------------------------------
      # PBMC + and - SEB_ATMi
      #----------------------------------------
      # grep("plus.SEB_ATMi", colnames(tbl.in))  # avg: 7 8  sd: 24 25
      # grep("minus.SEB_ATMi", colnames(tbl.in))  # avg: 15 16 sd: 32 33

   tbl.minus.seb.atmi.avg <- as.numeric(tbl.in[, c(13, 14)])
   tbl.minus.seb.atmi.sd <- as.numeric(tbl.in[, c(30, 31)])

      #----------------------------------------
      # PBMC + and - SEB_DNAPKi
      #----------------------------------------
      # grep("plus.SEB_DNAPKi", colnames(tbl.in))  # avg:  9 10 sd: 26 27
      # grep("minus.SEB_DNAPKi", colnames(tbl.in)) # avg: 17 18 sd: 34 35

   tbl.plus.seb.dnapki.avg <- as.numeric(tbl.in[, c(9, 10)])
   tbl.plus.seb.dnapki.sd <- as.numeric(tbl.in[, c(26, 27)])

   tbl.minus.seb.dnapki.avg <- as.numeric(tbl.in[, c(17, 18)])
   tbl.minus.seb.dnapki.sd <- as.numeric(tbl.in[, c(34, 35)])

      #----------------------------------------
      # PBMC + and - SEB_ATRi
      #----------------------------------------
      # grep("plus.SEB_ATRi", colnames(tbl.in))  #  avg: 11 12 sd: 28 29
      # grep("minus.SEB_ATRi", colnames(tbl.in))  # avg: 19 20 sd: 36 37

   tbl.plus.seb.atri.avg <- as.numeric(tbl.in[,  c(11, 12)])
   tbl.minus.seb.atri.avg <- as.numeric(tbl.in[, c(19, 20)])
   tbl.plus.seb.atri.sd <- as.numeric(tbl.in[,  c(28, 29)])
   tbl.minus.seb.atri.sd <- as.numeric(tbl.in[, c(36, 37)])

      #---------------------------------------
      # conditions summarized
      # no time course
      # radiation: 0 5Gy
      #  SEB_DMSO: +/-
      #  SEB_ATMi: +/-
      #  SEB_DNAPKi: +/-
      #  SEB_ATRi: +/-
      #---------------------------------------
      # DNA-PKcs: DNA-dependent protein kinase, catalytic subunit, also known as DNA-PKcs
      #  repair doubled-stranded breaks
      # SEB? The Supplemented Enzymatic Buffer (SEB) is a proprietary substrate stabilizing
      #   component that may be used with KinEASE TK assay to optimize the enzymatic reaction
      #   for some tyrosine kinases while maintaining inhibitor properties such as IC50.

   time <- rep(1, 16)
   radiation <- rep(c(0, 5), 8)
   experiment <- c("dmso.plus", "dmso.plus", "dmso.minus", "dmso.minus",
                   "atmi.plus", "atmi.plus", "atmi.minus", "atmi.minus",
                   "dnapki.plus", "dnapki.plus", "dnapki.minus", "dnapki.minus",
                   "atri.plus", "atri.plus", "atri.minus", "atri.minus")
   area <- as.numeric(c(tbl.in[, c(grep("PBMC.plus.SEB_DMSO_0Gy_1$", colnames(tbl.in)),
                                   grep("PBMC.plus.SEB_DMSO_5Gy_1$", colnames(tbl.in)),
                                   grep("PBMC.minus.SEB_DMSO_0Gy_1$", colnames(tbl.in)),
                                   grep("PBMC.minus.SEB_DMSO_5Gy_1$", colnames(tbl.in)),

                                   grep("PBMC.plus.SEB_ATMi_0Gy_1$", colnames(tbl.in)),
                                   grep("PBMC.plus.SEB_ATMi_5Gy_1$", colnames(tbl.in)),
                                   grep("PBMC.minus.SEB_ATMi_0Gy_1$", colnames(tbl.in)),
                                   grep("PBMC.minus.SEB_ATMi_5Gy_1$", colnames(tbl.in)),

                                   grep("PBMC.plus.SEB_DNAPKi_0Gy_1$", colnames(tbl.in)),
                                   grep("PBMC.plus.SEB_DNAPKi_5Gy_1$", colnames(tbl.in)),
                                   grep("PBMC.minus.SEB_DNAPKi_0Gy_1$", colnames(tbl.in)),
                                   grep("PBMC.minus.SEB_DNAPKi_5Gy_1$", colnames(tbl.in)),

                                   grep("PBMC.plus.SEB_ATRi_0Gy_1$", colnames(tbl.in)),
                                   grep("PBMC.plus.SEB_ATRi_5Gy_1$", colnames(tbl.in)),
                                   grep("PBMC.minus.SEB_ATRi_0Gy_1$", colnames(tbl.in)),
                                   grep("PBMC.minus.SEB_ATRi_5Gy_1$", colnames(tbl.in)))]))

  sd <-  as.numeric(c(tbl.in[, c(grep("PBMC.plus.SEB_DMSO_0Gy_1.1$", colnames(tbl.in)),
                                   grep("PBMC.plus.SEB_DMSO_5Gy_1.1$", colnames(tbl.in)),
                                   grep("PBMC.minus.SEB_DMSO_0Gy_1.1$", colnames(tbl.in)),
                                   grep("PBMC.minus.SEB_DMSO_5Gy_1.1$", colnames(tbl.in)),

                                   grep("PBMC.plus.SEB_ATMi_0Gy_1.1$", colnames(tbl.in)),
                                   grep("PBMC.plus.SEB_ATMi_5Gy_1.1$", colnames(tbl.in)),
                                   grep("PBMC.minus.SEB_ATMi_0Gy_1.1$", colnames(tbl.in)),
                                   grep("PBMC.minus.SEB_ATMi_5Gy_1.1$", colnames(tbl.in)),

                                   grep("PBMC.plus.SEB_DNAPKi_0Gy_1.1$", colnames(tbl.in)),
                                   grep("PBMC.plus.SEB_DNAPKi_5Gy_1.1$", colnames(tbl.in)),
                                   grep("PBMC.minus.SEB_DNAPKi_0Gy_1.1$", colnames(tbl.in)),
                                   grep("PBMC.minus.SEB_DNAPKi_5Gy_1.1$", colnames(tbl.in)),

                                   grep("PBMC.plus.SEB_ATRi_0Gy_1.1$", colnames(tbl.in)),
                                   grep("PBMC.plus.SEB_ATRi_5Gy_1.1$", colnames(tbl.in)),
                                   grep("PBMC.minus.SEB_ATRi_0Gy_1.1$", colnames(tbl.in)),
                                   grep("PBMC.minus.SEB_ATRi_5Gy_1.1$", colnames(tbl.in)))]))


   tbl.out <- data.frame(analyte=rep(analyte, length(time)),
                         experiment=experiment,
                         time=as.factor(time),
                         radiation=radiation,
                         area=area,
                         sd=sd,
                         stringsAsFactors=FALSE)

  tbl.out$time <- as.factor(tbl.out$time)
  tbl.out

} # toPlotFriendlyTable.03
#----------------------------------------------------------------------------------------------------
test_toPlotFriendlyTable.03 <- function()
{
   printf("--- test_toPlotFriendlyTable.03")

     # fix, eg, PBMC+ and PBMC- to PBMC.plus.  and PBMC.minus. in colnames, by
     # manual editing the .txt file.  R does not permit arithmetic operators in colnames
   tbl.orig <- read.table(files[3], sep="\t", header=TRUE, as.is=TRUE, nrow=-1)

   tbl.good <- toPlotFriendlyTable.03(tbl.orig[1,, drop=FALSE])
   checkEquals(colnames(tbl.good), c("analyte", "experiment", "time", "radiation", "area", "sd"))
   checkEquals(dim(tbl.good), c(16, 6))

     # eight treatments x 2 radiation levels (0 and 5Gy)
   checkEquals(sort(unique(tbl.good$experiment)),
              c("atmi.minus", "atmi.plus", "atri.minus", "atri.plus",
                "dmso.minus", "dmso.plus", "dnapki.minus", "dnapki.plus"))

   checkTrue(all(tbl.good$analyte %in% all.peptides))

   tbls <- lapply(1:nrow(tbl.orig), function(row) toPlotFriendlyTable.03(tbl.orig[row,, drop=FALSE]))
   checkEquals(length(tbls), 37)
   tbl.03 <- do.call(rbind, tbls)
   checkEquals(dim(tbl.03), c(592, 6))
   checkEquals(length(unique(tbl.03$analyte)), 37)
   checkTrue(all(tbl.03$analyte %in% all.peptides))

   save(tbl.03, file="tbl.03.RData")

} # test_toPlotFriendlyTable.03
#----------------------------------------------------------------------------------------------------
plotAnalyte.01 <- function(tbl.analyte)
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

} # plotAnalyte.01
#----------------------------------------------------------------------------------------------------
plotAnalyte.01 <- function(tbl.analyte)
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

} # plotAnalyte.01
#----------------------------------------------------------------------------------------------------
plotAnalyte.03 <- function(tbl.analyte)
{
   analyteName <- tbl.analyte$analyte[1]
   tbl.analyte$condition <- with(tbl.analyte, sprintf("%s:%s", experiment, radiation))
   browser()

   p <- ggplot(tbl.analyte) +
         geom_line(aes(y = area, x = condition, colour = experiment), stat="identity", size=1.5) +
         geom_errorbar(aes(ymin=area-sd, ymax=area+sd, x=condition, y=area, colour=experiment), width=0.7) +
         ggtitle(sprintf("%s", analyteName)) +
         theme(plot.title = element_text(size = 18, face = "bold")) +
         labs(x="condition", y="Ratio")# +
         #scale_x_continuous(breaks=1:5, labels=c("Mock", "0.25", "1.0", "6", "24"))

   p

} # plotAnalyte.01
#----------------------------------------------------------------------------------------------------
# the first file distributes colnames across treatment/control, and time
# does this work for all files?
#
#   tbl.atmi.avg <- tbl.in[, 5:9]
#   tbl.atmi.sd  <- tbl.in[, 17:21]
#
#   tbl.dmso.avg <- tbl.in[, 11:15]
#   tbl.dmso.sd  <- tbl.in[, 23:27]
#
#   time <- c(0, 0.25, 1, 6, 24)
#   atmi <- as.numeric(tbl.atmi.avg[1,])
#   atmi.sd <-    as.numeric(tbl.atmi.sd[1,])
#   dmso <- as.numeric(tbl.dmso.avg[1,])
#   dmso.sd <-  as.numeric(tbl.dmso.sd[1,])
#   tbl.tmp <- data.frame(time=as.factor(time), atmi=atmi, atmi.sd=atmi.sd, dmso=dmso, dmso.sd=dmso.sd)
#
examineColNames <- function()
{
  files <- list.files(pattern=".txt")
  length(files)

  tbl.01 <- read.table(files[1], sep="\t", header=TRUE, as.is=TRUE, nrow=-1)
  tbl.02 <- read.table(files[2], sep="\t", header=TRUE, as.is=TRUE, nrow=-1)
  tbl.03 <- read.table(files[3], sep="\t", header=TRUE, as.is=TRUE, nrow=-1)
  tbl.04 <- read.table(files[4], sep="\t", header=TRUE, as.is=TRUE, nrow=-1)

  peptides.01 <- sort(tbl.01$PeptideLabel)
  peptides.02 <- sort(tbl.02$Peptide.Label)
  peptides.03 <- sort(tbl.03$PeptideLabel)
  peptides.04 <- sort(tbl.04$PeptideLabel)

  head(peptides.01)
  head(peptides.02)
  head(peptides.03)
  head(peptides.04)

  length(peptides.01)   # 37
  length(peptides.02)   # 39
  length(peptides.03)   # 37
  length(peptides.04)   # 33

  all.peptides <- sort(unique(c(peptides.01, peptides.02, peptides.03, peptides.04)))
  length(all.peptides)  # 39
  save(all.peptides, file="all.peptides.RData")

  wdth(20);colnames(tbl.01)
  colnames(tbl.02)
  colnames(tbl.03)
  colnames(tbl.04)

  return(all.peptides)

} # examineColNames
#----------------------------------------------------------------------------------------------------
createWebAppReadyDataFrame <- function()
{
   tbl.orig <- read.table(files[1], sep="\t", header=TRUE, as.is=TRUE, nrow=-1)
   tbls.01 <- lapply(1:nrow(tbl.orig), function(row) toPlotFriendlyTable.01(tbl.orig[row,, drop=FALSE]))
   length(tbls.01)  # 39
   tbl.01 <- do.call(rbind, tbls.01)
   dim(tbl.01)  # 370 6

   tbl.orig <- read.table(files[2], sep="\t", header=TRUE, as.is=TRUE, nrow=-1)
   tbls.02 <- lapply(1:nrow(tbl.orig), function(row) toPlotFriendlyTable.02(tbl.orig[row,, drop=FALSE]))
   length(tbls.02)  # 39
   tbl.02 <- do.call(rbind, tbls.02)
   dim(tbl.02)  # 195 6



   tbl.orig <- read.table(files[3], sep="\t", header=TRUE, as.is=TRUE, nrow=-1)
   tbls.03 <- lapply(1:nrow(tbl.orig), function(row) toPlotFriendlyTable.03(tbl.orig[row,, drop=FALSE]))
   length(tbls.03)  # 37
   tbl.03 <- do.call(rbind, tbls.03)
   dim(tbl.03)  # 592 6

   tbl.123 <- rbind(tbl.01, tbl.02, tbl.03)
   dim(tbl.123)   # 1157  6
   head(tbl.123)
   length(unique(tbl.123$analyte))  # 39
   tbl <- tbl.123
   group <- rep("none", nrow(tbl))
   atmi.pure <- grep("ATMi", tbl$experiment)
   atmi.mixed <- sort(c(grep("atmi.minus", tbl$experiment), grep("atmi.plus", tbl$experiment)))
   atri.mixed <- sort(c(grep("atri.minus", tbl$experiment), grep("atri.plus", tbl$experiment)))
   dmso.pure <- grep("DMSO", tbl$experiment)
   dmso.mixed <- sort(c(grep("dmso.minus", tbl$experiment), grep("dmso.plus", tbl$experiment)))
   dnapki.mixed <- sort(c(grep("dnapki.minus", tbl$experiment), grep("dnapki.plus", tbl$experiment)))
   ir.pure <- grep("ir", tbl$experiment)

   group[atmi.pure]  <- "ATMi"
   group[atmi.mixed] <- "ATMi +/-"
   group[atri.mixed] <- "ATRi +/-"
   group[dmso.pure]  <- "DMSO"
   group[dmso.mixed] <- "DMSO +/-"
   group[dnapki.mixed] <- "DNAPKI +/-"
   group[ir.pure] <- "IR"

   stopifnot(length(which(group == "none")) == 0)
   tbl$group <- group
   na.area <- which(is.na(tbl$area))
   na.sd <- which(is.na(tbl$sd))
   if(length(na.area) > 0)
       tbl$area[na.area] <- 0
   if(length(na.sd) > 0)
       tbl$sd[na.sd] <- 0

   stopifnot(length(which(is.na(tbl$area))) == 0)
   stopifnot(length(which(is.na(tbl$sd))) == 0)

   save(tbl, file="tbl.39analytes.1157x6.RData")

} # createWebAppReadyDataFrame
#----------------------------------------------------------------------------------------------------
all.peptides <- examineColNames()
