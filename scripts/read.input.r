setwd(paste(root, "/input",sep=""))
ped.in<-param[which(param[,1] == "PED"),2]
map.in<-param[which(param[,1] == "MAP"),2]

if(ped.in == "")
  {
  LOG<-c(LOG, "No input files specified. Continuing with reference-only analysis")
  cat("No input files specified. Continuing with reference-only analysis\n")
  }else{
cat("Loading specified input files\n")
ped<-as.matrix(read.delim(ped.in, header=FALSE))
map<-as.matrix(read.delim(map.in, header=FALSE))


if(nrow(ped) > 200)
{
  #confirmDialog(LABEL="Warning",message="samples > 200; you are strongly advised to split the analysis into groups of < 200 samples")
}


# Check formatting

test<-suppressWarnings(as.numeric(ped[,1]))
test<-which(is.na(test))
if(length(test) > 0){stop("\nNon-numeric characters found in family identifiers: Only numeric identifiers are accepted.\n")}
test<-suppressWarnings(as.numeric(ped[,2]))
test<-which(is.na(test))
if(length(test) > 0){stop("\nNon-numeric characters found in sample identifiers: Only numeric identifiers are accepted.\n")}

test<-nchar(ped[1,])

snp.count.a<-length(which(nchar(ped[1,(3:ncol(ped))]) == 3))
if(snp.count.a != nrow(map)){stop("Number of SNPs in ped file doesn't match those in map file\n")}

test.pheno<-unique(nchar(ped[,6]))
if(test.pheno == 3)
{
 cat("No phenotype column detected. Setting all samples to controls\n")
 ped<-cbind(ped[,1:5], rep(1,nrow(ped)),ped[,6:ncol(ped)])
}
colnames(ped)<-c("FID","IID","dad","mum","sex","pheno",map[,4])

test<-median(test)
if(test > 3)
{
 stop("Genotype data contains too many characters\n")
}
if(test == 2)
{
 cat("Genotypes are not delimitted. Adding space\n")
 ped.a<-ped[,1:6]
 ped.b<-ped[, 7:ncol(ped)]
  ped.b<-apply(X=ped.b,FUN=sub,MARGIN=1, pattern="A", replacement="A ")
  ped.b<-apply(X=ped.b,FUN=sub,MARGIN=1, pattern="T", replacement="T ")
  ped.b<-apply(X=ped.b,FUN=sub,MARGIN=1, pattern="G", replacement="G ")
  ped.b<-apply(X=ped.b,FUN=sub,MARGIN=1, pattern="C", replacement="C ")
  ped.b<-substr(ped.b,1,3)
  blank<-grep("N", ped.b)
  if(length(blank) > 0){ped.b<-gsub(ped.b[blank[1]],"0 0",ped.b)}
  ped<-cbind(ped.a, ped.b)
}
if(test == 3)
{
  test.2<-substr(ped[1,10],2,2)
  if(test.2 != " ")
  {
    cat("Changing genotype delimitter from ", test.2, " to space\n")
    ped.a<-ped[,1:6]
    ped.b<-ped[, 7:ncol(ped)]
    ped.b<-gsub(test.2," ",ped.b)
    ped<-cbind(ped.a, ped.b)
  }
}

test<-ncol(map)
if(test != 4)
{
 stop("Incorrect number of columns in map input file. Expecting 4 columns\n")
}

### Filter to genomic interval
# SNP
  keep<-which(as.numeric(map[,4]) >= 7562060 & as.numeric(map[,4]) <= 12099352)
  map<-map[keep,]
  sub.ped<-ped[,(keep+6)]
     test<-as.numeric(colnames(sub.ped))
     test<-length(which(test >= 7562060 &  test <= 12099352))/ncol(sub.ped)
     if(test != 1){stop("Genomic interval range filter failed\n")}
  ped<-cbind(ped[,1:6], sub.ped)

}# first "ifelse" toggle