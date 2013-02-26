###############################################
### Accesory function to plot PFIDO results ###
###############################################
# version 2 (24/10/2010)
 
### Parameters
# ALPHA=0.7*255               # Transparency  (0 < alpha < 1)
# add.prediction=TRUE         # Colour in PFIDO predictions
# add.inset=TRUE              # Add a density estimate plot of the first dimension, as an inset panel
# add.text.label=TRUE         # Add inversion-type labels to the plot
# ASP=2.8/3                   # Plot aspect ratio (y:x)
# LEGEND=TRUE                 # Add a plot legend or not.
# STEP.X=2                    # Inter-tick ratio on x-axis
# STEP.Y=2                    # Inter-tick ratio on y-axis

PFIDOtoMDSplot<-function(X=x, ASP=3/3, ALPHA=0.7*255,  add.prediction=FALSE, add.inset=TRUE, add.text.label=TRUE, LEGEND=TRUE, STEP.X=2, STEP.Y=2)
{

  # Reference data
  REF<-as.matrix(read.delim(file=paste(root, "/scripts/FISHres_8pinv.txt", sep=""), header=TRUE))
     
   ### Annotate input data
   # Expecting no header in row 1
   test<-grep("sig.dim", X[1,])
   if(length(test)>0){stop("Header in first row, please remove\n")}
   add<-rep(NA, nrow(X))
   for(i in 1:length(add))
   {
    hit<-grep(substr(X[i,1],1,7),REF[,1])
     if(length(hit)==1){add[i]<-REF[hit,2]}
   }
   rem<-grep("NA", X[,2])
   if(length(rem) > 0){X[rem,2]<-NA}
   
   ### Create colour scheme
   require(RColorBrewer)
   my.palette<-brewer.pal(n=12,"Paired")
my.palette<-col2rgb(my.palette)
  col.a<-1  # Light blue
col.a<-rgb(red=my.palette[1,col.a], green=my.palette[2,col.a], blue=my.palette[3,col.a],alpha=ALPHA,maxColorValue =255)
  col.b<-2  # Dark blue
col.b<-rgb(red=my.palette[1,col.b], green=my.palette[2,col.b], blue=my.palette[3,col.b],alpha=ALPHA,maxColorValue =255)
  col.c<-7  # Light orange
col.c<-rgb(red=my.palette[1,col.c], green=my.palette[2,col.c], blue=my.palette[3,col.c],alpha=ALPHA,maxColorValue =255)
  col.d<-8  # Dark orange
col.d<-rgb(red=my.palette[1,col.d], green=my.palette[2,col.d], blue=my.palette[3,col.d],alpha=ALPHA,maxColorValue =255)
  join.line<-col2rgb(brewer.pal(n=8,"Set2")[8])
join.line<-rgb(red=join.line[1], green=join.line[2], blue=join.line[3],maxColorValue =255)
  highlight<-6  # Dark red
highlight<-rgb(red=my.palette[1,highlight], green=my.palette[2,highlight], blue=my.palette[3,highlight],alpha=ALPHA,maxColorValue =255)
  outline<-brewer.pal(n=8,"Set2")[8]   # Grey
   purple<-10   # Dark purple
purple<-rgb(red=my.palette[1,purple], green=my.palette[2,purple], blue=my.palette[3,purple],alpha=ALPHA,maxColorValue =255)
  line.3<-11   # Yellow
line.3<-rgb(red=my.palette[1,line.3], green=my.palette[2,line.3], blue=my.palette[3,line.3],maxColorValue =255)
  light.green<-3  #Light green
light.green<-rgb(red=my.palette[1,light.green], green=my.palette[2,light.green], blue=my.palette[3,light.green],alpha=ALPHA,maxColorValue =255)
  dark.green<-4   # Dark green
dark.green<-rgb(red=my.palette[1,dark.green], green=my.palette[2,dark.green], blue=my.palette[3,dark.green],alpha=ALPHA,maxColorValue =255)
   light.red<-5   # Light red
light.red<-rgb(red=my.palette[1,light.red], green=my.palette[2,light.red], blue=my.palette[3,light.red],alpha=ALPHA,maxColorValue =255)
darker.grey<-brewer.pal(n=8,"Accent")[8]


   ### Create empty plot
   par(col="black",col.axis=darker.grey)
   plot(x=X[,4], y=X[,5] ,type="n", xlab=NA, ylab=NA, axes=FALSE, asp=ASP)
   axis(pretty(x=as.numeric(X[,4]), n=STEP.X), side=1, pos=0, cex=0.3, col.axis=join.line)
   axis(pretty(x=as.numeric(X[,5]), n=STEP.Y), side=2, pos=0, cex=0.3, col.axis=join.line)

   # Add legend
   if(add.prediction & LEGEND)
   {
   LEG<-c("II (FISH)","IN (FISH)","NN (FISH)","","II (PFIDO)", "IN (PFIDO)","NN (PFIDO)","No cluster")
   COL<-c(highlight,col.b, dark.green, "white",light.red, col.a, light.green, purple)
   SYM<-c(18,15,17,17,18,15,17,4)
   legend(x="bottomright", legend=LEG, col = COL, pch=SYM, cex=0.8, box.col=join.line, ncol=2)
   }
   if(!(add.prediction) & LEGEND)
   {
   LEG<-c("II","IN","NN","-")
   COL<-c(highlight,col.b, dark.green, purple)
   SYM<-c(18,15,17,4)
   legend(x="bottomright", legend=LEG, col = COL, pch=SYM, cex=0.8, box.col=join.line)
   }

   # Annotate PFIDO predictions (internal)
   if(add.prediction)
   {
     temp<-cbind(X[,2], add)
     hit<-which(is.na(add))
     addition<-paste(X[hit,2], "p", sep="")
     add[hit]<-addition
      rem<-grep("NA", add)
      if(length(rem)>0){add[rem]<-NA}
   }
   add<-gsub(" ","",add)

    # Add sample points without predictions called
      draw<-which(is.na(add))
   if(length(draw) > 0){
    if(length(draw) == 1){points(x=X[draw,4],y=X[draw,5], pch=4, col=purple, cex=0.7)}else{
   points(X[draw,4:5], pch=4, col=purple, cex=0.7)}}

   # Add data points - FISH first, then overlay with PFIDO prediction
      draw<-which(add == "II")
   points(X[draw,4:5], pch=18, col=highlight, cex=1.1)
      draw<-which(add == "NI")
   points(X[draw,4:5], pch=15, col=col.b, cex=1.1)
      draw<-which(add == "NN")
   points(X[draw,4:5], pch=17, col=dark.green, cex=1.1)
     if(add.prediction)
     {
          draw<-which(add == "IIp")
         points(X[draw,4:5], pch=18, col=light.red, cex=0.9)
          draw<-which(add == "NIp")
         points(X[draw,4:5], pch=15, col=col.a, cex=0.8)
           draw<-which(add == "NNp")
         points(X[draw,4:5], pch=17, col=light.green, cex=0.8)
     }
   # Mark FISH results called as NA by PFIDO
    draw<-which(!is.na(add) & is.na(X[,2]))
    #draw<-which(is.na(X[,2]))
    if(length(draw) > 0){points(X[draw,4:5], pch=4, col=purple, cex=0.7)}
 
   

   if(add.text.label)
   {
     hit<-which(add %in% c("II","NI","NN"))
     text(x=X[hit,4:5], pos=3, labels=add[hit])
   }

   if(add.inset)
   {
     require(ade4)
     f1 <- function(a){
        opar=par("mar","xaxt","yaxt","plt","col")
        on.exit(par(opar))
        par(mar=rep(.1,4),xaxt="n",yaxt="n",plt=par("plt"), col=join.line)
        plot(density(as.numeric(a)), axes=TRUE, main=NA,xlab=NA, ylab=NA, col="darkgrey")
        #axis(pretty(x=as.numeric(X[,4]), n=4), side=1, pos=0)
        #axis(pretty(x=density(as.numeric(X[,4])), n=2), side=2, pos=0)

        }
     add.scatter(f1(X[,4]),posi="topright", inset=c(0,0),bg.col="transparent")
   }

} # End of function


