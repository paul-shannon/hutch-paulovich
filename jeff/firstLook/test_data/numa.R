library(RUnit)
files <- dir(pattern="prelim")
files
tbl.02 <- read.csv(files[2], as.is=TRUE)
dim(tbl.02)  # 42924 53
peptide.rows <- grep("NUMA", tbl.02[, "Protein.Name"])
length(peptide.rows) # 686
tbl.numa <- tbl.02[peptide.rows,]
dim(tbl.numa)
# tbl.numa[, "SampleVariable1"]
# coi.area <- grep("area", colnames(tbl.numa), ignore.case=TRUE,v=T)       # "Area"          "Area.Ratio" "Total.Area"
# coi.ratio <- grep("Ratio", colnames(tbl.numa), ignore.case=FALSE,v=T)     #  "Area.Ratio"  "Ratio.To.Standard"

coi <- c("Protein.Name", "Peptide.Modified.Sequence", "Transition",
         "Replicate.Name", "SampleVariable1",
         "Area.Ratio", "Ratio.To.Standard", "Total.Area",
         "Isotope.Label.Type", "Retention.Time", "DotProductLightToHeavy", "Area", "Max.Height")
tbl.numa <- tbl.numa[, coi]
dim(tbl.numa)  # 686 13
full.name <- tbl.numa$Protein.Name
geneSymbol <- sub("_HUMAN", "", unlist(lapply(strsplit(full.name, "\\|"), "[", 3)))
tbl.numa$geneSymbol  <- geneSymbol
dim(tbl.numa)  # 686 14

#----------------------------------------------------------------------------------------------------
# 2,5, 10 Gy: radiation measure, the gray is the unit of absorbed dose and has replaced the rad.
# 1 gray = 1 Joule/kilogram and also equals 100 rad.
parseSampleName <- function(s)
{
   tokens <- strsplit(s, "_")[[1]]
   stopifnot(length(tokens) %in% 3:5)
   stopifnot(tokens[1] %in% c("RefQC", "LCL57"))

   if(length(tokens) == 3){
      #printf("tokens 3")
      return(data.frame(cellType=tokens[1],
                        treatment=tokens[2],
                        radiation=NA,
                        time=NA,
                        rep=as.integer(tokens[3]),
                        stringsAsFactors=FALSE))
      }

   if(length(tokens) == 4){
      #printf("tokens 4")
      return(data.frame(cellType=tokens[1],
                        treatment=NA,
                        radiation=tokens[2],
                        time=60 * as.numeric(sub("Hr", "", tokens[3])),
                        rep=as.integer(sub("Rep", "", tokens[4])),
                        stringsAsFactors=FALSE))
      }

   if(length(tokens) == 5){
      #printf("tokens 5")
      return(data.frame(cellType=tokens[1],
                        treatment=tokens[2],
                        radiation=tokens[3],
                        time=60 * as.numeric(sub("Hr", "", tokens[4])),
                        rep=as.integer(sub("Rep", "", tokens[5])),
                        stringsAsFactors=FALSE))

      }

} # parseSampleName
#----------------------------------------------------------------------------------------------------
test_parseSampleName <- function()
{
   x <- tbl.numa$SampleVariable1
   set.seed(17)
   checkEquals(parseSampleName(x[1]),
               data.frame(cellType="RefQC", treatment="01", radiation=NA, time=NA, rep=1,
                          stringsAsFactors=FALSE))
   checkEquals(parseSampleName(x[5]),
               data.frame(cellType="LCL57", treatment=NA, radiation="5Gy", time=60, rep=3,
                          stringsAsFactors=FALSE))


   checkEquals(parseSampleName(x[6]),
               data.frame(cellType="LCL57", treatment="ATMi", radiation="5Gy", time=15, rep=1,
                          stringsAsFactors=FALSE))




} # test_parseSampleName
#----------------------------------------------------------------------------------------------------
dim(tbl.numa) # 686 14
sampleNames <- tbl.numa$SampleVariable1

tbl.parsed <- do.call(rbind, lapply(sampleNames, parseSampleName))

dim(tbl.numa)
dim(tbl.parsed)
tbl.x <- cbind(tbl.numa, tbl.parsed)

coi <-  c(
    "geneSymbol",
    "treatment",
    "radiation",
    "time",
    "Transition",
    "rep",
    "Peptide.Modified.Sequence",
    "Area.Ratio",
    "Ratio.To.Standard",
    "Total.Area",
    "Isotope.Label.Type",
    "Retention.Time",
    "DotProductLightToHeavy",
    "Area",
    "Max.Height",
    "Protein.Name",
    "Replicate.Name",
    "SampleVariable1",
    "cellType")

tbl.numa <- tbl.x[, coi]
dim(tbl.numa) # 686 14

tbl.numaLight <- subset(tbl.numa, Isotope.Label.Type=="light")
tbl.numaLight <- subset(tbl.numaLight, cellType == "LCL57")
rownames(tbl.numaLight) <- NULL
dim(tbl.numaLight)  # 343 19

sampleNames <- c("LCL57_1Gy",
                 "LCL57_2Gy",
                 "LCL57_5Gy",
                 "LCL57_10Gy",
                 "LCL57_ATMi_mock",
                 "LCL57_ATMi_5Gy",
                 "LCL57_DMSO_mock",
                 "LCL57_DMSO_5Gy",
                 "LCL57_mock")

groupName <- rep("untreated", nrow(tbl.numaLight))
groupName[grep("LCL57_1Gy", tbl.numaLight$SampleVariable1)] <- "untreated_1Gy"
groupName[grep("LCL57_2Gy", tbl.numaLight$SampleVariable1)] <- "untreated_2Gy"
groupName[grep("LCL57_5Gy", tbl.numaLight$SampleVariable1)] <- "untreated_5Gy"
groupName[grep("LCL57_10Gy", tbl.numaLight$SampleVariable1)] <- "untreated_10Gy"

groupName[grep("LCL57_ATMi_mock", tbl.numaLight$SampleVariable1)] <- "ATMi_mock"
groupName[grep("LCL57_ATMi_5Gy", tbl.numaLight$SampleVariable1)] <- "ATMi_5Gy"
groupName[grep("LCL57_DMSO_mock", tbl.numaLight$SampleVariable1)] <- "DMSO_mock"
groupName[grep("LCL57_DMSO_5Gy", tbl.numaLight$SampleVariable1)] <- "DMSO_5Gy"
groupName[grep("LCL57_mock", tbl.numaLight$SampleVariable1)] <- "mock"

tbl.numaLight$group <- groupName
dim(tbl.numaLight) # 315 20

wdth(20); colnames(tbl.numaLight); wdth(120)

coi <- c(
"geneSymbol",
"group",
"treatment",
"radiation",
"time",
"Transition",
"rep",
"Total.Area",
"Area",
"Max.Height",
"Peptide.Modified.Sequence",
"Area.Ratio",
"Ratio.To.Standard",
"Isotope.Label.Type",
"Retention.Time",
"DotProductLightToHeavy",
"Protein.Name",
"Replicate.Name",
"SampleVariable1",
"cellType")

tbl.numaLight$Area.Ratio <- as.numeric(tbl.numaLight$Area.Ratio)

save(tbl.numaLight, file="tbl.numa.RData")







# preferred.column.order <- c("transition", "treatment", "radiation", "time", "rep",
#                             "totalArea", "retentionTime",
#                             "gene", "protein", "sequence", "cellType", "label",
#                             "ratioToStandard", "sampleName", "group")
#
# tbl.numaLight <- tbl.numaLight[, preferred.column.order]
# save(tbl.numa, tbl.numaLight, file="tbl.numa.RData")
#
#
# varNames <- grep("LCL57", tbl.numaLight$sampleName, v=T)
# varNames.2 <- unique(sub("Hr_Rep.*$", "", varNames))
#
# varNames <- grep("LCL57", tbl.numa$SampleVariable1, v=T)
# length(varNames.2)
# # remove trailing time values: 1, 0.25, 6, 24
#
#
#
#
# sort(unique(c(
# "LCL57_5Gy",
# "LCL57_ATMi_5Gy",
# "LCL57_DMSO_mock",
# "LCL57_ATMi_5Gy",
# "LCL57_ATMi_mock",
# "LCL57_10Gy",
# "LCL57_DMSO_5Gy",
# "LCL57_2Gy",
# "LCL57_DMSO_5Gy",
# "LCL57_1Gy",
# "LCL57_ATMi_5Gy",
# "LCL57_ATMi_5Gy",
# "LCL57_mock",
# "LCL57_DMSO_5Gy",
# "LCL57_DMSO_5Gy")))

