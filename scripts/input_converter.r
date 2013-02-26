### Convert HapMap ###
cat("Converting HapMap download data\n")

source("C:/Documents and Settings/msalm/My Documents/scripts_0110/Hapmap_to_LINKAGE.r")

ceu.tran<-h2l(X="C:/Documents and Settings/msalm/My Documents/HapMapIII_r28/CEU_r27")
jpt.tran<-h2l(X="C:/Documents and Settings/msalm/My Documents/HapMapIII_r28/JPT")
yri.tran<-h2l(X="C:/Documents and Settings/msalm/My Documents/HapMapIII_r28/YRI")
tsi.tran<-h2l(X="C:/Documents and Settings/msalm/My Documents/HapMapIII_r28/TSI")
chb.tran<-h2l(X="C:/Documents and Settings/msalm/My Documents/HapMapIII_r28/CHB")

setwd("C:/Documents and Settings/msalm/My Documents/scripts_0110/PFIDO/input")

temp<-ceu.tran$ped[,c(2,2,3:ncol(ceu.tran$ped))]
write.table(temp,file="CEU.ped",col.names=F, row.names=F, quote=F, sep="\t")
write.table(ceu.tran$map,file="CEU.map",col.names=F, row.names=F, quote=F, sep="\t")

temp<-tsi.tran$ped[,c(2,2,3:ncol(tsi.tran$ped))]
write.table(temp,file="TSI.ped",col.names=F, row.names=F, quote=F, sep="\t")
write.table(tsi.tran$map,file="TSI.map",col.names=F, row.names=F, quote=F, sep="\t")

temp<-chb.tran$ped[,c(2,2,3:ncol(chb.tran$ped))]
write.table(temp,file="CHB.ped",col.names=F, row.names=F, quote=F, sep="\t")
write.table(chb.tran$map,file="CHB.map",col.names=F, row.names=F, quote=F, sep="\t")


