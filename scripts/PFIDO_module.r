### PFIDO modules for snpMatrix ###
# 23/02/2010
# "X" are the SNPs from a snp.matrix object, with SNPs only mapping to the inversion

### Change Log
# 26/8/2010   - Added outlier detection section

### Options explained:
# X         = genotype data
# file.type = Is it the "snp.data" entry of the snp.matrix (0) or not (1)?
# call.rate = the per SNP succesful genotyping rate threshold
# mind      = the per sample genotyping rate threshold 
# DIM       = The number of dimensions to be calculated during MDS
# THRESHOLD = The PFIDO p-value threshold
# plots     = Plot MDS values (boolean)
# NOM       = Designate output file name for MDS plot
# out.nom   = Designate output file name for PFIDO results 
# restrict  = Are only a subset of SNPs meant to be analysed (boolean)?
# snp.ids   = SNP identifiers for the SNPs meant to be analysed
# GFORCE    = Would you like to restrict the model-fitting stage to 1:3 components(boolean)?
# MODEL     = Specify mclust model (default is "E" & "V")
# OUTLIER   = Would you like to identify & remove potential outliers (boolean)?
# VALIDATE   = Would you like to perform the clValid stage (boolean)?

### Short-cut for de-bugging
# X=SNP$snp.data; file.type=0; call.rate = 0.9; mind = 0.9; DIM = 4; THRESHOLD = 0.05;plots = FALSE; out.nom="PFIDO_output.txt";restrict=FALSE; snp.ids=NA; GFORCE=FALSE; OUTLIER=FALSE; MODEL=c("E","V"); VALIDATE=TRUE


PFIDO<-function(X, file.type=1, call.rate = 0.9, mind = 0.9, DIM = 4, THRESHOLD = 0.05,plots = TRUE, NOM="PFIDO_diagnostics.pdf", out.nom="PFIDO_output.txt",restrict=FALSE, snp.ids=NA, GFORCE=FALSE, OUTLIER=FALSE, MODEL=c("E","V"), VALIDATE=TRUE)
{
### Ensure that all dependencies are loaded, and of a correct version
temp<-require(clValid)
      if(!(temp)){stop("Please install the clValid package\n")}
temp<-require(mclust)
      if(!(temp)){stop("Please install the mclust package\n")}
temp<-require(extremevalues)
      if(!(temp)){stop("Please install the extremevalues package\n")}
temp<-require(moments)
      if(!(temp)){stop("Please install the moments package\n")}      
ver<-Biobase::package.version("snpMatrix")
if(ver != "1.8.0"){warning("PFIDO was tested with snpMatrix v.1.8.0\n")}
ver<-Biobase::package.version("mclust")
if(ver != "3.4.6"){warning("PFIDO was tested with mclust v.3.4.6\n")}

      
flush.console()

 cat("PFIDO initiated...\n")
	#############
	### Input ###
	#############
  CLASS<-class(X)
  if(CLASS[length(CLASS)] != "snp.matrix"){stop("Input data format not recognised\n")}

  if(file.type==1)
  {
    genotypes<-X$snp.data
  }else{genotypes<-X}
cat("\tFile loaded\n")

  ##############
  ### SNP QC ###
  ##############
# Restrict analyses to SNPs specified in "snp.ids"
if(restrict==TRUE)
{
  hit<-which(colnames(genotypes) %in% snp.ids)    
  if(length(hit) > 0)
  {
  cat("Found ",((length(hit))/(length(snp.ids)))*100,"% requested SNPs\n")
   genotypes<-genotypes[,hit]
  }else{warning("Specified dbSNP identifiers not found\n")}
}

# Remove SNPs with low genotyping (N.B. Doing this first improves sample retention rate)
# and individual samples are more informative than individuals SNPs
qc.snp<-col.summary(genotypes)
out.snp<-which(qc.snp[,"Call.rate"] < call.rate)
	if(length(out.snp)>0){genotypes<-genotypes[,-out.snp]}

# Remove individulas with low genotyping
qc.ind<-row.summary(genotypes)
out.ind<-which(qc.ind[,"Call.rate"] < mind)
	if(length(out.ind)>0){genotypes<-genotypes[-out.ind,]}

# Remove monomorphic alleles
qc.snp<-col.summary(genotypes)
mono<-which(qc.snp[,"MAF"] == 0)
	if(length(mono)>0){genotypes<-genotypes[,-mono]}

# Report on the QC stage
cat("\tAfter filtering:\n")
  if(file.type == 1)
  {
  test<-X$snp.data
  test<-dim(test)
  ref<-dim(genotypes)
    cat("\t\tAt geno= ", call.rate ," -SNPs removed: ",test[2]-ref[2] ," or ",(test[2]-ref[2])/test[2],"\n")
    cat("\t\tAt mind= ",mind," -Samples removed: ",test[1] - ref[1] ," or ",(test[1]-ref[1])/test[1],"\n") 
  }else{
  test<-X
  test<-dim(test)
  ref<-dim(genotypes)
    cat("\t\tAt geno= ", call.rate ," -SNPs removed: ",test[2]-ref[2] ," or ",(test[2]-ref[2])/test[2],"\n")
    cat("\t\tAt mind= ", mind ," -Samples removed: ",test[1] - ref[1] ," or ",(test[1]-ref[1])/test[1],"\n") 
    cat(ncol(genotypes)," SNPs in analysis\n")
    cat(nrow(genotypes)," samples in analysis\n")
  }


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

mds<-cmdscale(d=ibs.dist, k=DIM, eig=TRUE, x.ret = TRUE)  # classical MDS 
cat("\tMDS analysis complete\n")

	###########################
	### Shapiro-Wilk's Test ###
	###########################
eigener<-cbind(mds$eig, rep(NA,length(mds$eig)))
	for(i in 1:nrow(eigener))
	{
	eigener[i,2]<-shapiro.test(mds$points[,i])$p
	}
sig.dim<-which(eigener[,2] == min(eigener[,2]))
  test.dim<-which(eigener[,2] < 0.05)
  if(length(test.dim) == 0){stop("A Gaussian distribution found in all dimensions tested...try increasing DIM\n")}
cat("Dimension: ",sig.dim,"\n")
if(sig.dim != 1){sig.dim<-1;cat("Dimension forced: ",sig.dim,"\n")}

  #########################
  ### Outlier detection ###
  #########################
  # Optional
  if(OUTLIER)
  {
    outlier<-boxplot(as.numeric(mds$points[,sig.dim]), plot=TRUE, main="Outlier detection")
    box.outliers<-outlier$out
    box.outliers<-which(mds$points[,sig.dim] %in% box.outliers)
    box.outliers<-rownames(mds$points)[box.outliers]           # Sample IDs of boxplot outliers
    outlier<-getOutliers(y=as.numeric(mds$points[,sig.dim]), rho = 0.1, pval = c(0.05, 0.85), method = "normal")
    outlier<-outlier$iOut
      if(length(which(outlier) > 0))
      {
        cat(length(which(outlier)), " possible outliers found...removing\n")
        # Repeat IBS-MDS stage, without outliers 
    non.outlier<-rownames(mds$points)[!outlier]
    count<-ibsCount(genotypes[non.outlier])
    ibs.dist<-ibsDist(count) 
    mds<-cmdscale(d=ibs.dist, k=DIM, eig=TRUE, x.ret = TRUE)
      }else{cat("No outliers detected\n")} 
  }
  
  #########################
  ### Data distribution ###
  #########################
  da<-agostino.test(mds$points[,sig.dim])
  test.again<-length(which(MODEL == "V"))/length(MODEL)
  if(da$p.value < 0.05 & test.again == 0.5){warning("Data is significantly skewed: consider using variable variance only in mclust\n")}
  
	##############
	### mclust ###
	##############
# Perform a cluster analysis combining model-based hierarchical clustering,
# EM for Gaussian mixture models and BIC. Single dimension only (univariate: for efficacy with YRI 8p23) 

test.1<-try(Mclust(mds$points[,sig.dim], G=c(1:9), modelNames=MODEL, control=emControl()))     

if(GFORCE){test.1<-Mclust(mds$points[,sig.dim], G=c(1:3), modelNames=MODEL)}
den<-test.1

    # Explore BIC to make sure the correct model has been called
    test.me<-which(MODEL %in% c("E","V"))
    if(length(test.me) == 2)
    {
    explore<-test.1$BIC
     E.score<-which.max(explore[,1])
     V.score<-which.max(explore[,2])
     explore<-round(explore - max(explore))
     if(E.score == V.score){cat("Same mode chosen whether variance set to E or V \n")}
     iqr<-IQR(mds$points[,sig.dim])
     }
      
modes<-den$G
OBJ<-as.matrix(mds$points[,sig.dim])
cat("mclust first pass: ",modes,"\n")
  
  ############################
  ### Validate clustering  ###
  ############################
  # See http://cran.r-project.org/web/views/Cluster.html
  # See also Handl et al (2005):"Computational cluster validation in post-genomic data analysis"
  if(VALIDATE)
  {
  mode.range<-c(modes-2, modes+2)
    if(mode.range[1] < 2){mode.range[1]<-2}    # Silhouette method requires K > 2
    if(mode.range[2] < mode.range[1]){mode.range[2]<-2}
    validate<-try(clValid(obj=OBJ, nClust=mode.range[1]:mode.range[2], maxitems = 3000, clMethods = c("model"), validation=c("internal")))
     test.5<- class(validate)
     test.5<-grep("clValid", test.5)   
     if(length(test.5) > 0)
     { 
      summary(validate)
      #report.me<-optimalScores(validate)
      #report.me<-as.numeric(as.matrix(report.me)[,3])
      #report.me<-length(which(report.me == modes))
      #cat(report.me, " out of 3 cluster validation methods concur with mclust clustering\t")
      validate<-measures(validate)
     }else{cat("core mclust validation step aborted\n");validate<-NA}
  }else{validate<-NA} 
   
  ###  Also provide singular comparisons
  # e.g. Dunn index only
  # This addition is crucial to JPT sample clustering assesment
  # N.B. This reports specific metrics for the mclust stage solution (i.e. "modes")
  dunn.index<-dunn(Data=OBJ,clusters=den$classification)
  conn.index<-connectivity(Data=OBJ,clusters=den$classification)
  silh.width<-silhouette(x=den$classification, dist=dist(OBJ)) 
  silh.width<-summary(silh.width)$avg.width 
  
  ### Compare to PAM (Partitioning Around Mediods, a robust k-means)
  # Vestigial section
  #PAM<-clValid(obj=OBJ, nClust=(modes-1):(modes+1), clMethods = c("pam"), validation=c("internal"))
  #PAM<-measures(PAM)
  
  ### Compare to hclust, using "complete" agglomeration method
  # Vestigial section
  if(1 == 0)
  {
  PAM<-try(clValid(obj=OBJ, nClust=(modes-1):(modes+1), clMethods = c("hierarchical"), validation=c("internal"), method = "complete",))
    test.5<- class(PAM)
    test.5<-grep("clValid", test.5)   
   if(length(test.5) > 0)
   {
    PAM<-optimalScores(PAM)
    temp<-as.numeric(PAM[[3]][2])
    if(temp == modes){cat("hclust/Dunn concur with clustering\n")}
    temp<-as.numeric(PAM[[3]][3])
    if(temp == modes){cat("hclust/Silhouette concur with clustering\n")}
   }else{cat("Comparison to hclust failed\n")}
  }

  ### Retrieve original mclust results    
calls<-den$classification
p.value<-rep(NA, length(calls))
	for(i in 1:length(calls))
	{
	p.value[i]<-den$z[i,calls[i]]
	}
p.value<-1-p.value
false.call<-which(p.value >= THRESHOLD)
called<-((length(calls) - length(false.call)) / length(calls)) *100
	if(length(false.call)>0)
	{
	calls[false.call]<-paste(calls[false.call],NA, sep="-")
	}

  cat("\tClustering complete\n")

 	##############
	### Output ###
	##############
# Diagnostic plots
if(plots == TRUE)
{
	pdf(NOM)
	plot(mds$points[,sig.dim:(sig.dim+1)], pch=20, main="MDS plot", type="n", xlab=paste("Eigenvector ",sig.dim, sep=""), ylab=paste("Eigenvector ",sig.dim+1, sep=""))
	N<-length(unique(calls))
	types<-unique(calls)
	shape<-c(18,15,17,4)
	mycol<-c("brown4", "dodgerblue3", "cyan4", "black")
		for(i in 1:N)
		{
		draw<-which(calls == types[i])
			if(is.na(types[i])){draw<-which(is.na(calls))}
		points(x=mds$points[draw,sig.dim],y=mds$points[draw,(sig.dim+1)], pch=shape[i], col=mycol[i])
		}
	legend(x="bottomright", legend=types, col=mycol, pch=shape)

	plot(density(mds$points[,1]), main=paste("Distribution along axis", sig.dim))
	dev.off()
}

# Inversion predictions
PED<-names(calls)
PED<-cbind(PED, calls)
PED<-cbind(PED, p.value)
PED<-cbind(PED, mds$points[,sig.dim])
PED<-cbind(PED, mds$points[,(sig.dim+1)])
colnames(PED)[4:5]<-c("sig.dim", "sig.dim+1")
write.table(PED, file=out.nom, sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
cat("PFIDO complete.\n\n\n")

pfido.out<-PED
cluster.validate<-validate

return(list(pfido.out=pfido.out, cluster.validate=cluster.validate, dunn.index=dunn.index,conn=conn.index,silh.width=silh.width,snps.used=colnames(genotypes)))
}


##############################################################
# A function characterising simple 1-D cluster variance for hets, on a PFIDO object
cl.spread<-function(X=keep, K=2)
{
  hit<-X[which(as.numeric(X[,"calls"]) == K),ncol(X)]
  out<-as.numeric(sd(hit, na.rm=TRUE))
return(out)
}

##############################################################
# A function to draw MDS output relative to FISHED samples
fish.verify<-function(X=RES1, Z=FISH, LEG=TRUE)
{
  temp<-X
  temp[,"p.value"]<-NA
    for(i in 1:nrow(temp))
    {
     hit<-which(Z[,1] == temp[i,1])
      if(length(hit) > 0)
      {
       temp[i,"p.value"]<-Z[hit,2]
      }
    }
rem<-which(is.na(temp[,"calls"]))
temp[rem,"calls"]<-"NS"    
    
# Draw PFIDO-typed
plot(x=temp[,"sig.dim"],y=temp[,"sig.dim+1"], pch=20, main="MDS plot", type="n")
	N<-length(unique(temp[,"calls"]))
	types<-unique(temp[,"calls"])
	types<-types[order(types, decreasing=FALSE)]
	shape<-c(18,15,17,4)
	mycol<-c("brown4", "dodgerblue3", "cyan4", "black")
		for(i in 1:N)
		{
		draw<-which(temp[,"calls"] == types[i])
		points(x=temp[draw,"sig.dim"],y=temp[draw,"sig.dim+1"], pch=shape[i], col=mycol[i])
		}
	if(LEG){legend(x="topright", legend=types, col=mycol, pch=shape)}
  
# Draw FISH-typed
    text(x=as.numeric(temp[,"sig.dim"]),y=as.numeric(temp[,"sig.dim+1"]),labels=temp[,"p.value"])
}

##############################################################
# A function to report on IBS/category

uIBS<-function(X=YRI.inv, pfido.res=keep.3, mind=0.585, call.rate=0.9)
{
genotypes<-X

# Remove individulas with low genotyping
qc.ind<-row.summary(genotypes)
out.ind<-which(qc.ind[,"Call.rate"] < mind)
	if(length(out.ind)>0){genotypes<-genotypes[-out.ind,]}

# Remove SNPs with low genotyping (N.B. Doing this first improves sample retention rate)
qc.snp<-col.summary(genotypes)
out.snp<-which(qc.snp[,"Call.rate"] < call.rate)
	if(length(out.snp)>0){genotypes<-genotypes[,-out.snp]}

# Remove monomorphic alleles
qc.snp<-col.summary(genotypes)
mono<-which(qc.snp[,"MAF"] == 0)
	if(length(mono)>0){genotypes<-genotypes[,-mono]}

# Remove NI/NS individuals
rem<-pfido.res[which(pfido.res[,"calls"] == "2"),1]
rem<-which(rownames(genotypes) %in% rem)
genotypes<-genotypes[-rem,]
rem<-pfido.res[which(is.na(pfido.res[,"calls"])),1]
rem<-which(rownames(genotypes) %in% rem)
genotypes<-genotypes[-rem,]

count<-ibsCount(genotypes)
ibs.dist<-ibsDist(count)
explore<-as.matrix(ibs.dist)

II<-pfido.res[which(pfido.res[,"calls"] == "1"),1]
II<-which(rownames(explore) %in% II)
NN<-pfido.res[which(pfido.res[,"calls"] == "3"),1]
NN<-which(rownames(explore) %in% NN)
 # NA out self-self comparisons, leaving only IBS calculations between II and NN samples
  for(i in 1:length(II))
  {
   ROW<-II[i]
   explore[ROW, II]<-NA 
  }
   for(i in 1:length(NN))
  {
   ROW<-NN[i]
   explore[ROW, NN]<-NA 
  }
explore<-as.vector(explore)
rem<-which(is.na(explore))
explore<-explore[-rem]
rem<-which(duplicated(explore))
explore<-explore[-rem]

MEAN<-mean(explore, na.rm=TRUE)
SD<-sd(explore, na.rm=TRUE)
return(list(MEAN=MEAN, SD=SD))
} 


##############################################################
# A function to report on individual Gaussian distribution parameters.

Gdist<-function(X, file.type=1, call.rate = 0.3, mind = 0.585, DIM = 4, THRESHOLD = 0.05,plots = TRUE, NOM="PFIDO_diagnostics.pdf", out.nom="PFIDO_output.txt",restrict=FALSE, snp.ids=NA)
{
flush.console()

if(file.type==1)
{
genotypes<-X$snp.data
}else{genotypes<-X}
cat("\tFile loaded\n")

# Restrict to specified SNPs
if(restrict==TRUE)
{
  hit<-which(colnames(genotypes) %in% snp.ids)    
  if(length(hit) > 0)
  {
  cat("Found ",((length(hit))/(length(snp.ids)))*100,"% requested SNPs\n")
   genotypes<-genotypes[,hit]
  }
}

# Remove individulas with low genotyping
qc.ind<-row.summary(genotypes)
out.ind<-which(qc.ind[,"Call.rate"] < mind)
	if(length(out.ind)>0){genotypes<-genotypes[-out.ind,]}

# Remove SNPs with low genotyping (N.B. Doing this first improves sample retention rate)
qc.snp<-col.summary(genotypes)
out.snp<-which(qc.snp[,"Call.rate"] < call.rate)
	if(length(out.snp)>0){genotypes<-genotypes[,-out.snp]}

# Remove monomorphic alleles
qc.snp<-col.summary(genotypes)
mono<-which(qc.snp[,"MAF"] == 0)
	if(length(mono)>0){genotypes<-genotypes[,-mono]}


cat("\tAfter filtering:\n")
  if(file.type == 1)
  {
  test<-X$snp.data
  test<-dim(test)
  ref<-dim(genotypes)
    cat("\t\tSNPs removed: ",test[2]-ref[2] ,"\n")
    cat("\t\tSamples removed: ",test[1] - ref[1] ,"\n") 
  }else{
  test<-X
  test<-dim(test)
  ref<-dim(genotypes)
    cat("\t\tSNPs removed: ",test[2]-ref[2] ,"\n")
    cat("\t\tSamples removed: ",test[1] - ref[1] ,"\n") 
  }

count<-ibsCount(genotypes)
ibs.dist<-ibsDist(count)
mds<-cmdscale(d=ibs.dist, k=DIM, eig=TRUE, x.ret = TRUE)

cat("\tMDS analysis complete\n")

	###########################
	### Shapiro-Wilk's Test ###
	###########################
eigener<-cbind(mds$eig, rep(NA,length(mds$eig)))
	for(i in 1:nrow(eigener))
	{
	eigener[i,2]<-shapiro.test(mds$points[,i])$p
	}
sig.dim<-which(eigener[,2] == min(eigener[,2]))

	##############
	### mclust ###
	##############
# Perform a cluster analysis combining model-based hierarchical clustering,
# EM for Gaussian mixture models and BIC. Single dimension only (for efficacy with YRI 8p23) 
 # test.1<-Mclust(mds$points[,sig.dim], prior= priorControl())    # With prior
test.1<-Mclust(mds$points[,sig.dim])                              # Without prior
den<-test.1

NA18912<-den$classification
NA18912<-NA18912[grep("18912",names(NA18912))]

CT<-summary(factor(den$classification))
CT<-rbind(CT,den$parameters$mean)
SD<-sqrt(den$parameters$variance$sigmasq)
if(length(SD) != ncol(CT)){SD<-rep(SD,3)}
CT<-rbind(CT, SD)
return(list(CT, NA18912))
}

