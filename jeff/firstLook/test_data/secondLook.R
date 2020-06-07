files <- dir(pattern="prelim")
files
tbl.01 <- read.table(files[1], sep=",", as.is=TRUE, nrow=-1, fill=TRUE, quote="", header=TRUE)
dim(tbl.01)  # 29694  71

tbl.02 <- read.table(files[2], sep=",", as.is=TRUE, nrow=-1, fill=FALSE, quote="", header=TRUE)
tbl.02 <- read.csv(files[2], as.is=TRUE)
dim(tbl.02)  # 42924 53
peptide.rows <- grep("NUMA", tbl.02[, "Protein.Name"])
length(peptide.rows) # 686
tbl.numa <- tbl.02[peptide.rows,]
dim(tbl.numa)
tbl.numa[, "SampleVariable1"]

#counts <- lapply(tbl.01, function(col) length(grep("QLQDNPPQEK", col)))
#counts[counts > 0]
  # $Peptide.Sequence [1] 686
  # $Peptide.Modified.Sequence [1] 686

tbl <- tbl.01[grep("QLQDNPPQEK", tbl.01$Peptide.Sequence),]
dim(tbl)
wdth(40); colnames(tbl); wdth(120)
wdth(10000);head(tbl)

as.data.frame(t(tbl[1:2,]))

#--------------------------------------------------------------------------------
# try to find the useful columns in tbl.01
#--------------------------------------------------------------------------------
empty.columns <- which(as.logical(lapply(tbl.01, function(col) all(is.na(col)))))
if(length(empty.columns) > 0)
    tbl.01 <- tbl.01[, -empty.columns]
dim(tbl.01)   # 29694 14
tbl.numa.01 <- tbl.01[grep("NUMA", tbl.01$Protein.Name),]
dim(tbl.numa.01)   # 686 14
fivenum(tbl.numa.01$Ratio.To.Standard) # 0.0100 0.0275 0.1591 1.2782 2.6117

#--------------------------------------------------------------------------------
# try to find the useful columns in tbl.02
#--------------------------------------------------------------------------------
dim(tbl.02)
empty.columns <- which(as.logical(lapply(tbl.02, function(col) all(is.na(col)))))
length(empty.columns) # 21
if(length(empty.columns) > 0)
    tbl.02 <- tbl.02[, -empty.columns]
dim(tbl.02)   # 29694 32
tbl.numa.02 <- tbl.02[grep("NUMA", tbl.01$Protein.Name),]
dim(tbl.numa.02)   # 686 32
fivenum(tbl.numa.02$Ratio.To.Standard) # 0.0100 0.0275 0.1591 1.2782 2.6117
set.seed(17)
some.random.rows <- sort(sample(1:686, 3))
as.data.frame(t(tbl.numa.02[some.random.rows,]))


wdth(20)
colnames(tbl.numa.01)
wdth(2000)


tbl.01.coi.area <- grep("area", colnames(tbl.01), ignore.case=TRUE,v=T)
tbl.01.coi.ratio <- grep("ratio", colnames(tbl.01), ignore.case=TRUE,v=T)
tbl.01.coi.ratio.max <- lapply(tbl.01.coi.ratio, function(colname) max(tbl.01[, colname], na.rm=TRUE))
names(tbl.01.coi.ratio.max) <- tbl.01.coi.ratio
fivenum(tbl.01$Ratio.To.Standard) # [1]   0.0017   0.0722   0.3523   2.2450 377.0333
tbl <- tbl.01[grep("QLQDNPPQEK", tbl.01$Peptide.Sequence),]
fivenum(tbl$Ratio.To.Standard)   # 0.0100 0.0275 0.1591 1.2782 2.6117

names(tbl.01.coi.max) <- tbl.01.coi.area
tbl.01.coi.max




tbl.02.coi <- grep("area", colnames(tbl.02), ignore.case=TRUE,v=T)

tbl.01.coi.max <- lapply(tbl.01.coi.area, function(colname) max(tbl.01[, colname], na.rm=TRUE))
names(tbl.01.coi.max) <- tbl.01.coi.area
tbl.01.coi.max



tbl.02.coi.max <- lapply(tbl.02.coi, function(colname) max(tbl.02[, colname], na.rm=TRUE))
names(tbl.02.coi.max) <- tbl.02.coi
tbl.02.coi.max
