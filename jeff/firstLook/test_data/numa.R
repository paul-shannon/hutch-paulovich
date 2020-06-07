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

coi <- c("Protein.Name", "SampleVariable1", "Area.Ratio", "Ratio.To.Standard", "Replicate", "Total.Area", "SampleVariable2")
tbl.numa <- tbl.numa[, coi]
full.name <- tbl.numa$Protein.Name
geneSymbol <- sub("_HUMAN", "", unlist(lapply(strsplit(full.name, "\\|"), "[", 3)))
tbl.numa$geneSymbol  <- geneSymbol

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
dim(tbl.numa)
sampleNames <- tbl.numa$SampleVariable1

tbl.parsed <- do.call(rbind, lapply(sampleNames, parseSampleName))

tbl.x <- cbind(tbl.numa[, c(1,3,4,6,8,2)], tbl.parsed)

preferred.column.order <- c(
 "geneSymbol",
 "Protein.Name",
 "cellType",
 "treatment",
 "radiation",
 "time",
 "rep",
 "Area.Ratio",
 "Ratio.To.Standard",
 "Total.Area",
 "SampleVariable1")
tbl.numa <- tbl.x[, preferred.column.order]
save(tbl.numa, file="tbl.numa.RData")


