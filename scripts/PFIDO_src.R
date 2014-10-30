      ###########################
      ### PFIDO master script ###
      ###########################

#rm(list=ls())     # Remove all existing objects in R workspace
LOG<-date()       # Create an empty LOG file

### Load (and download) required libraries
temp<-require(snpMatrix)
      if(!(temp)){
       source("http://bioconductor.org/biocLite.R")
       biocLite("snpMatrix")
      }
temp<-require(clValid)
      if(!(temp)){install.packages("clValid", dependencies=TRUE, repos="http://cran.r-project.org")}
temp<-require(extremevalues)
      if(!(temp)){install.packages("extremevalues", dependencies=TRUE, repos="http://cran.r-project.org")}
temp<-require(moments)
      if(!(temp)){install.packages("moments", dependencies=TRUE, repos="http://cran.r-project.org")}
temp<-require(ade4)
      if(!(temp)){install.packages("ade4", dependencies=TRUE, repos="http://cran.r-project.org")}
temp<-require(mclust)
      if(!(temp)){install.packages("mclust", dependencies=TRUE, repos="http://cran.r-project.org")}

### Check package versions
ver<-Biobase::package.version("snpMatrix")
if(ver != "1.8.0"){warning("PFIDO was tested with snpMatrix v.1.8.0\n")}
ver<-Biobase::package.version("mclust")
if(ver != "3.4.8"){warning("PFIDO was tested with mclust v.3.4.8\n")}
if(ver == "3.4.1"){warning("Some windows binaries of this mclust version (v.3.4.1) are not stable. Please download an alternative\n")}
ver<-Biobase::package.version("clValid")
if(ver != "0.6-1"){warning("PFIDO was tested with clValid v.0.5-7\n")}
ver<-Biobase::package.version("extremevalues")
if(ver != "2.1"){warning("PFIDO was tested with extremevalues v.1.0\n")}
ver<-Biobase::package.version("moments")
if(ver != "0.11"){warning("PFIDO was tested with moments v.0.11\n")}

### Directory structure
# "C:/Documents and Settings/msalm/My Documents/scripts_0110/PFIDO"
root<-getwd();
root<-unlist(strsplit(root, "/"));
root<-paste(root[1:which(root == "PFIDO_gui")], sep="/", collapse="/");
setwd(root);

### Check that mclust functions as expected
setwd(paste(root, "/scripts", sep=""));
source("mclust.test.r");

### Retrieve analysis parameters
setwd(paste(root,"/tmp",collapse="", sep=""));
param<-as.matrix(read.delim("parameters.txt",header=FALSE));
  LOG<-c(LOG, "Parameters selected:");
  LOG<-c(LOG, param);
param<-matrix(unlist(strsplit(param,"=")), ncol=2, byrow=TRUE);


###########################
### Retreive input data ###
###########################
setwd(paste(root, "/scripts",sep=""));
source("read.input.r");

################################
### Create snp.matrix object ###
################################
setwd(paste(root,"/tmp",sep=""));
ref.only<-ls()
ref.only<-length(which(ref.only == "ped"))
if(ref.only > 0)
  {
  write.table(ped,file="temp.ped", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE);
  write.table(map,file="temp.map", sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE);
  
  rm(list=c("ped", "map"))
  new.map<-read.pedfile.map("temp.map");
  new.ped<-read.snps.pedfile(file="temp.ped", snp.names = new.map[,1], X = FALSE, sep = " ", low.mem = FALSE);
  }

###################
### HapMap data ###
###################
setwd(paste(root,"/scripts",sep=""));
source("read.hapmap.r");

########################
### Combine datasets ###
########################
setwd(paste(root,"/scripts",sep=""));
if(ref.only > 0)
  {
  source("combine.data.r");
  }else{
  combined.data<-hapmap
  }
##################
### Run PFIDO  ###
##################
setwd(paste(root,"/scripts",sep=""));
source("PFIDO_section.r");

################################
### Create diagnostic figure ###
################################

test<-param[which(param[,1] == "plot"),2]
if(test == "Y")
{
  setwd(paste(root,"/scripts",sep=""));
  source("PFIDOtoMDSplot.r");
  setwd(paste(root,"/output",sep=""));
  file.nom<-param[which(param[,1] == "plot.name"),2]
  pdf(file.nom)
    PFIDOtoMDSplot(X=result, add.text=F, add.inset=T, LEG=T, add.prediction=T)
  dev.off()
  LOG<-c(LOG, "Annotated MDS plot saved to output folder")
}

LOG<-c(LOG, "...Finished")
LOG<-c(LOG, date())

setwd(paste(root,"/output" , sep=""))
write.table(LOG, file="log.txt", sep="\n", col.names=FALSE, row.names=FALSE, quote=FALSE)

cat("...Finished\n")
setwd(root)

### GUI to confirm program has run
confirmDialog <- function(message)
	{
	window <- gwindow("PFIDO")
	group <- ggroup(container = window)
	gimage("info", dirname="stock", container=group)
  ## A group for the message and buttons
	inner.group <- ggroup(horizontal=FALSE, container = group)
	glabel(message, container=inner.group, expand=TRUE)
  ## A group to organize the buttons
	button.group <- ggroup(container = inner.group)
  ## Push buttons to right
	addSpring(button.group)
	gbutton("ok", handler=function(h,...) dispose(window), container=button.group)
	return()
	}
confirmDialog(message="PFIDO run completed.\nPlease check the output folder for results")
setwd(root)


