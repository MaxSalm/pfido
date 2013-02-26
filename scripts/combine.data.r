
  new.snps<-new.ped$snp.data

  keep<-which(colnames(hapmap) %in% colnames(new.snps))
  if(length(keep) > 0){hapmap<-hapmap[,keep]}else{stop("No matching SNP identifiers found between reference and test dataset - aborting\n")}
  keep<-which(colnames(new.snps) %in% colnames(hapmap))
  if(length(keep) > 0){new.snps<-new.snps[,keep]}else{stop("No matching SNP identifiers found between reference and test dataset - aborting\n")}

  re.order<-order(colnames(hapmap))
  hapmap<-hapmap[,re.order]
  re.order<-order(colnames(new.snps))
  new.snps<-new.snps[,re.order]

  # Check against original genotype calls, and flip alleles if necessary
  get.me<-colnames(hapmap)
  get.me<-full.hapmap$snp.support[which(rownames(full.hapmap$snp.support) %in% get.me),]
  get.me.2<-colnames(new.snps)
  get.me.2<-new.ped$snp.support[which(names(new.ped$snp.support) %in% get.me.2)]

  get.me<-cbind(get.me, get.me.2)
  one<-as.vector(get.me[,"Assignment"])
  two<-as.vector(get.me[,"get.me.2"])

  # Refined version
    list.snp<-cbind(one, two)
    flip.me<-vector()
    for(zz in 1:nrow(list.snp))
    {
      geno.a<-unlist(strsplit(list.snp[zz,1],"/"))
      geno.b<-unlist(strsplit(list.snp[zz,2],"/"))
      coll<-c(geno.a, geno.b)
        rem<-which(coll == "-")
        if(length(rem) > 0){coll<-coll[-rem]}
        coll<-unique(coll)
      if(length(coll) > 2){flip.me<-c(flip.me, zz)}
    }
    if(length(flip.me) > 0)
    {
    cat("Flipping ",length(flip.me)," SNPs\n")
  get.me<-rownames(get.me[flip.me,])
  get.me<-which(colnames(new.snps) %in% get.me)
  new.snps<-switch.alleles(x=new.snps, snps=get.me)
    }

  temp<-snp.rbind(hapmap, new.snps)
  FA<-rownames(temp)
  rem<-grep("NA", FA)
  FA[rem]<-2
  FA[-rem]<-1
  FA<-as.numeric(FA)
  flips<-test.allele.switch(snps=temp, split=FA)
  caught<-which(flips > 50)  # Remove those with a decibit > 5; http://en.wikipedia.org/wiki/Bayes_factor
  if(length(caught) > 0){cat("Removing ",length(caught)," SNPs with unresolved strandness issues\n");temp<-temp[,-caught]}

combined.data<-temp