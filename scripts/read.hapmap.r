  upload<-function(testurl=testurl, SAVE=T)
  {
  if(SAVE){YRI<- read.HapMap.data(testurl, save=paste("genotypes_chr8_",ref.in,"_r27_nr.b36_fwd.txt.gz", sep=""))}else{YRI<-read.HapMap.data(testurl)};
  use<-which(YRI$snp.support[,4] > 7562060 & YRI$snp.support[,4] < 12099352);
  YRI$snp.data<-YRI$snp.data[,use];
  YRI$snp.support<-YRI$snp.support[use,];
  return(YRI)
  }
ref.in<-param[which(param[,1] == "REFERENCE"),2];
setwd(paste(root, "/tmp",sep=""))
test<-which(dir() == paste("genotypes_chr8_",ref.in,"_r27_nr.b36_fwd.txt", sep=""))
test.2<-file.info(paste("genotypes_chr8_",ref.in,"_r27_nr.b36_fwd.txt", sep=""))
if(is.na(test.2[1])){test.2<-0}
 if(length(test) == 1 & test.2[1] > 0)
 {
 testurl<-paste("file://genotypes_chr8_",ref.in,"_r27_nr.b36_fwd.txt", sep="")
 hapmap<-upload(testurl, SAVE=F);
 }else{
 testurl <- paste("ftp://ftp.ncbi.nlm.nih.gov/hapmap/genotypes/latest/forward/non-redundant/genotypes_chr8_",ref.in,"_r27_nr.b36_fwd.txt.gz", sep="")
hapmap<-upload(testurl, SAVE=T);
  }
full.hapmap<-hapmap 
hapmap<-hapmap$snp.data;

setwd(paste(root,"/scripts",sep=""))
# Founder filter
if(1 == 1)
{
  h3<-as.matrix(read.delim("relationships_w_pops_121708.txt", header=TRUE))
  h3<-h3[which(h3[,3] == "0" & h3[,3] == "0"),]
  hit<-which(rownames(hapmap) %in% h3[,2])
  if(length(hit) > 0){hapmap<-hapmap[hit,]}else{cat("HapMap II filter failed\n")}
}

# Filter to HapMap II
if(1 == 1 & ref.in != "YRI")
{
  h2<-as.matrix(read.delim("hapmap2samples.txt", header=FALSE))
  hit<-which(rownames(hapmap) %in% h2)
  if(length(hit) > 0){hapmap<-hapmap[hit,]}else{cat("HapMap II filter failed\n")}
}
