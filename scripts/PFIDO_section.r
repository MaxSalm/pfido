LOG<-c(LOG, "")
LOG<-c(LOG, "PFIDO stage")
LOG<-c(LOG, "-----------")
genotypes<-combined.data

  ##############
  ### SNP QC ###
  ##############
# Restrict analyses to SNPs specified in "snp.ids"
test<-param[which(param[,1] == "restrict"),2]
if(test=="Y")
{
  setwd(paste(root,"/scripts",sep=""))                      # SNP name list (without header) is expected here
  get.snp.names<-param[which(param[,1] == "snp.ids"),2]   # Retrieve user supplied filename
  snp.ids<-as.matrix(read.delim(get.snp.names, header=TRUE))
  snp.ids=snp.ids[,ref.in]
  snp.ids=snp.ids[which(snp.ids!="")]

  hit<-which(colnames(genotypes) %in% snp.ids)
  if(length(hit) > 0)
  {
  LOG<-c(LOG, paste("Found ",((length(hit))/(length(snp.ids)))*100,"% requested SNPs", sep=""))
   genotypes<-genotypes[,hit]
  }else{stop("Specified dbSNP identifiers not found\n")}
}

# Remove SNPs with low genotyping (N.B. Doing this first improves sample retention rate)
# and individual samples are more informative than individuals SNPs
qc.snp<-col.summary(genotypes)
call.rate<-as.numeric(param[which(param[,1] == "call.rate"),2])
out.snp<-which(qc.snp[,"Call.rate"] < call.rate)
	if(length(out.snp)>0){genotypes<-genotypes[,-out.snp]}

# Remove individulas with low genotyping
qc.ind<-row.summary(genotypes)
mind<-as.numeric(param[which(param[,1] == "mind"),2])
out.ind<-which(qc.ind[,"Call.rate"] < mind)
	if(length(out.ind)>0){genotypes<-genotypes[-out.ind,]}

# Remove monomorphic alleles
qc.snp<-col.summary(genotypes)
mono<-which(qc.snp[,"MAF"] == 0)
	if(length(mono)>0){genotypes<-genotypes[,-mono]}

# Report on the QC stage
LOG<-c(LOG, "After filtering:")
  test<-combined.data
  test<-dim(test)
  ref<-dim(genotypes)
   LOG<-c(LOG, paste("At geno= ", call.rate ," -SNPs removed: ",test[2]-ref[2] ," or ",((test[2]-ref[2])/test[2])*100,"%",sep=""))
   LOG<-c(LOG, paste("At mind= ", mind ," -Samples removed: ",test[1] - ref[1] ," or ",((test[1]-ref[1])/test[1])*100,"%",sep=""))
   LOG<-c(LOG, paste(ncol(genotypes),"SNPs in analysis"))
   LOG<-c(LOG,paste(nrow(genotypes)," samples in analysis"))



  #######################
  ### IBS calculation ###
  #######################
count<-ibsCount(genotypes)    # Count the total number of IBS alleles

test.count<-rep(NA, nrow(count))
for(zz in 1:nrow(count))
{
  test.count[zz]<-length(which(count[zz,] == 0))
}
rem<-which(test.count > nrow(count)/4)
  if(length(rem) > 0)
  {
    cat("Removing ",length(rem)," more samples: uninformative in ibsCount stage\n")
    rem<-which(rownames(genotypes) %in% rownames(count)[rem])
    genotypes<-genotypes[-rem, ]
    count<-ibsCount(genotypes)
  }

  ### Multidimensional scaling ###
ibs.dist<-ibsDist(count)          # Convert a matrix of IBS counts into a distance matrix

    test.ibs<-which(is.na(ibs.dist), arr.ind=TRUE)
    if(length(test.ibs) > 0){stop("NaNs found in distance matrix: increase SNP filtering stringency\n")}
DIM<-as.numeric(param[which(param[,1] == "DIM"),2])
mds<-cmdscale(d=ibs.dist, k=DIM, eig=TRUE, x.ret = TRUE)  # classical MDS
LOG<-c(LOG,"MDS analysis complete")

	###########################
	### Shapiro-Wilk's Test ###
	###########################
eigener<-cbind(mds$eig[1:DIM], rep(NA,DIM))
	for(i in 1:DIM)
	{
	eigener[i,2]<-shapiro.test(mds$points[,i])$p
	}
sig.dim<-which(eigener[,2] == min(eigener[,2]))
  test.dim<-which(eigener[,2] < 0.05)
  if(length(test.dim) == 0){stop("A Gaussian distribution found in all dimensions tested...try increasing DIM\n")}
if(sig.dim != 1){sig.dim<-1;cat("Dimension forced: ",sig.dim,"\n")}

  #########################
  ### Outlier detection ###
  #########################
  # Optional
  OUTLIER<-param[which(param[,1] == "OUTLIER"),2]
  if(OUTLIER == "Y")
  {
    setwd(paste(root,"/output",sep=""))
    pdf("outlier_boxplot.pdf")
    outlier<-boxplot(as.numeric(mds$points[,sig.dim]), plot=TRUE, main="Outlier detection")
    dev.off()
    box.outliers<-outlier$out
    box.outliers<-which(mds$points[,sig.dim] %in% box.outliers)
    box.outliers<-rownames(mds$points)[box.outliers]           # Sample IDs of boxplot outliers
    outlier<-getOutliers(y=as.numeric(mds$points[,sig.dim]), method="II", alpha=c(0.05, 0.05), distribution = "normal")
    outlier<-c(outlier$iRight, outlier$iLeft)
      if(length(outlier) > 0)
      {
        LOG<-c(LOG,paste(length(outlier), " possible sample outliers found...removing"))
        # Repeat IBS-MDS stage, without outliers identified by getOutliers
    non.outlier<-rownames(mds$points)[-outlier]
    count<-ibsCount(genotypes[non.outlier])
    ibs.dist<-ibsDist(count)
    mds<-cmdscale(d=ibs.dist, k=DIM, eig=TRUE, x.ret = TRUE)
      }else{LOG<-c(LOG,"No outliers detected")}
  }


  #########################
  ### Data distribution ###
  #########################
MODEL<-param[which(param[,1] == "MODEL"),2]
MODEL<-unlist(strsplit(MODEL,""))
MODEL<-MODEL[which(MODEL != ",")]
  da<-agostino.test(mds$points[,sig.dim])
  test.again<-length(which(MODEL == "V"))/length(MODEL)
  if(da$p.value < 0.05 & test.again == 0.5){warning("Data is significantly skewed: consider using variable variance only in mclust\n")}

  ##############
	### mclust ###
	##############
# Perform a cluster analysis combining model-based hierarchical clustering,
# EM for Gaussian mixture models and BIC. Single dimension only (univariate)
test.1<-try(Mclust(mds$points[,sig.dim], G=c(1:9), modelNames=MODEL, control=emControl()))
den<-test.1

    # Explore BIC to make sure the correct model has been called
    test.me<-which(MODEL %in% c("E","V"))
    if(length(test.me) == 2)
    {
    explore<-test.1$BIC
     E.score<-which.max(explore[,1])
     V.score<-which.max(explore[,2])
     explore<-round(explore - max(explore))
     if(E.score == V.score){LOG<-c(LOG,"Same mode chosen whether variance set to E or V")}
     iqr<-IQR(mds$points[,sig.dim])
     }

modes<-den$G
OBJ<-as.matrix(mds$points[,sig.dim])
cat("mclust first pass: ",modes,"\n")

  ############################
  ### Validate clustering  ###
  ############################
VALIDATE<-param[which(param[,1] == "VALIDATE"),2]
  if(VALIDATE == "Y")
  {
  mode.range<-c(2,9)    # Silhouette method requires K > 2
     validate<-suppressWarnings(try(clValid(obj=OBJ, nClust=mode.range[1]:mode.range[2], maxitems = 10000, clMethods = c("model"), validation=c("internal"))))
     test.5<- class(validate)
     test.5<-grep("clValid", test.5)
     if(length(test.5) > 0)
     {
      setwd(paste(root,"/output",sep=""))
      
        out.validate<-measures(validate)
      #out.validate<-try(optimalScores(validate))  # De-preciated method
      #no.err<-grep("Error", out.validate)
      #if(length(no.err) == 0){
        write.table(out.validate, file="Cluster_validation_metrics.txt", quote=FALSE, sep="\t")
        LOG<-c(LOG, "Cluster validation metrics written to output folder")
         # }else{LOG<-c(LOG, "Cluster validation metrics NOT written to output folder: optimalScores error")}
   
   
          conn.temp<- measures(validate, measures="Connectivity")
        conn.temp<-dimnames(conn.temp)[[2]][which.min(conn.temp)]
       dunn.temp<-measures(validate, measures="Dunn")
         dunn.temp<-dimnames(dunn.temp)[[2]][which.max(dunn.temp)]       
      sil.temp<- measures(validate, measures="Silhouette")
         sil.temp<-dimnames(sil.temp)[[2]][which.max(sil.temp)]  
      validate<-measures(validate)
      tempa<-as.numeric(c(conn.temp,dunn.temp,sil.temp))
      best.cluster<-length(which(modes == tempa))    
        LOG<-c(LOG, paste(best.cluster, " out of 3 validation metrics concur with BIC"))
      
      if(best.cluster == 0)
      {
        override.BIC <- param[which(param[,1] == "override.BIC"),2]
        if(override.BIC == "Y")
        {
        forced.G<-median(tempa)
        LOG<-c(LOG, paste("In this event, you selected to over-ride BIC component selection, and use validation metrics instead: ",forced.G," components fit",sep=""))
         test.1<-try(Mclust(mds$points[,sig.dim], G=forced.G, modelNames=MODEL, control=emControl()))
         den<-test.1
         OBJ<-as.matrix(mds$points[,sig.dim])
        }else{cat("Continuing with original clustering\n")}
      }
     }else{cat("mclust validation step failed\n");validate<-NA}
  }else{validate<-NA}

  ###  Also provide singular comparisons
  # e.g. Dunn index only
  # This addition is crucial to JPT sample clustering assesment
  # N.B. This reports specific metrics for the mclust stage solution (i.e. "modes")
  dunn.index<-dunn(Data=OBJ,clusters=den$classification)
  conn.index<-connectivity(Data=OBJ,clusters=den$classification)
  silh.width<-silhouette(x=den$classification, dist=dist(OBJ))
  silh.width<-summary(silh.width)$avg.width


 ### Retrieve mclust results
calls<-den$classification
p.value<-rep(NA, length(calls))
	for(i in 1:length(calls))
	{
	p.value[i]<-den$z[i,calls[i]]
	}
p.value<-1-p.value
THRESHOLD<-as.numeric(param[which(param[,1] == "THRESHOLD"),2])
false.call<-which(p.value >= THRESHOLD)
called<-((length(calls) - length(false.call)) / length(calls)) *100
	if(length(false.call)>0)
	{
	calls[false.call]<-paste(calls[false.call],NA, sep="-")
	}

  cat("\tClustering complete\n")
  
  # Inversion predictions
PED<-names(calls)
PED<-cbind(PED, calls)
PED<-cbind(PED, p.value)
PED<-cbind(PED, mds$points[,sig.dim])
PED<-cbind(PED, mds$points[,(sig.dim+1)])
colnames(PED)[4:5]<-c("sig.dim", "sig.dim+1")

cat("PFIDO complete.\n\n\n")
result<-PED
setwd(paste(root, "/scripts",sep=""))

  ### Validate internal positive controls
FISH<-as.matrix(read.delim("FISHres_8pinv.txt", header=TRUE))
  add.inv<-rep(NA, nrow(result))
     for(B in 1:nrow(result))
     {
      hit<-which(FISH[,1] == PED[B,1])
         if(length(hit) > 0){add.inv[B]<-FISH[hit,2]}
     }
    rosetta<-cbind(result[,"calls"],add.inv)
    rosetta<-rosetta[!(is.na(rosetta[,2])),]
      rem<-grep("NA", rosetta[,1])
      if(length(rem) > 0){rosetta<-rosetta[-rem,]}
      one<-unique(rosetta[which(rosetta[,1] == "1"),2])
      two<-unique(rosetta[which(rosetta[,1] == "2"),2])
      three<-unique(rosetta[which(rosetta[,1] == "3"),2])
        tester<-length(one)+length(two)+length(three)
        if(tester != 3 & ref.in == "CEU"){stop("Not all CEU have been called correctly, aborting\n")}
        if(tester != 2 & ref.in == "JPT"){stop("Not all JPT have been called correctly, aborting\n")}

    # Re-format inversion-type calls, from cluster # to genotype
    temp<-split(rosetta[,2], f=rosetta[,1])
    temp<-lapply(X=temp, FUN=unique)
    for(TT in 1:length(temp))
    {
     result[,"calls"]<-gsub(labels(temp)[TT],temp[TT],result[,"calls"])
    }
    
    #result[,"calls"]<-gsub("1",one,result[,"calls"])   # Depreciated section
    #result[,"calls"]<-gsub("2",two,result[,"calls"])
    #result[,"calls"]<-gsub("3",three,result[,"calls"])

  # Re-format identifiers
     targets<-grep("Family",result[,1])
     if(length(targets) > 0)
     {
     temp<-matrix(unlist(strsplit(result[targets,1],"\\.")), ncol=4, byrow=TRUE)
     result[targets,1]<-temp[,4]
     }
     
setwd(paste(root, "/output",sep=""))
out.nom<-param[which(param[,1] == "out.nom"),2]
write.table(result, file=out.nom, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
     