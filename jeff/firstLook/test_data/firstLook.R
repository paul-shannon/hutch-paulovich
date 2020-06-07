files <- dir(pattern="prelim")
files

tbl.01 <- read.table(files[1], sep=",", as.is=TRUE, nrow=-1, fill=TRUE, quote="", header=TRUE)
dim(tbl.01)  # 29694  71
wdth(40); colnames(tbl.01); wdth(500)

tbl.02 <- read.table(files[2], sep=",", as.is=TRUE, nrow=-1, fill=TRUE, quote="", header=TRUE)
dim(tbl.02)  # 42924  53
wdth(40); colnames(tbl.02); wdth(500)


length(unique(tbl.02$Protein.Name))  #32
as.data.frame(sort(table(tbl.02$Protein.Name), decreasing=TRUE))
#
#                     Var1 Freq
# 1  sp|P20700|LMNB1_HUMAN 3136
# 2  sp|Q12888|TP53B_HUMAN 2842
# 3   sp|O96017|CHK2_HUMAN 2548
# 4   sp|P46013|KI67_HUMAN 2156
# 5  sp|P30305|MPIP2_HUMAN 2058
# 6    sp|P04637|P53_HUMAN 1568
# 7   sp|O96013|PAK4_HUMAN 1470
# 8   sp|Q00987|MDM2_HUMAN 1470
# 9   sp|Q14566|MCM6_HUMAN 1470
# 10  sp|O14757|CHK1_HUMAN 1372
# 11   sp|P04406|G3P_HUMAN 1372
# 12  sp|P06493|CDK1_HUMAN 1372
# 13 sp|P30307|MPIP3_HUMAN 1372
# 14 sp|P49959|MRE11_HUMAN 1372
# 15 sp|Q02223|TNR17_HUMAN 1372
# 16 sp|Q6IBW4|CNDH2_HUMAN 1372
# 17 sp|Q8NG31|CASC5_HUMAN 1372
# 18  sp|Q92541|RTF1_HUMAN 1372
# 19 sp|Q96ER3|SAAL1_HUMAN 1372
# 20 sp|Q9BVJ6|UT14A_HUMAN 1372
# 21 sp|Q9H400|LIME1_HUMAN 1372
# 22 sp|Q99638|RAD9A_HUMAN 1274
# 23   sp|O43561|LAT_HUMAN  686
# 24 sp|O75417|DPOLQ_HUMAN  686
# 25 sp|P09874|PARP1_HUMAN  686
# 26 sp|P42574|CASP3_HUMAN  686
# 27  sp|P50613|CDK7_HUMAN  686
# 28 sp|Q14980|NUMA1_HUMAN  686
# 29 sp|Q53HL2|BOREA_HUMAN  686
# 30 sp|Q92878|RAD50_HUMAN  686
# 31  sp|P07437|TBB5_HUMAN  490
# 32  sp|P68133|ACTS_HUMAN  490
