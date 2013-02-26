### KNOWN BUGS
# Have to do getwd(PFIDO folder) for anything to work...pain that needs to be resolved.
# If you click start again, it won't work. Have to close GUI, and reload source.
# WARNING sample size window; program continues uninterrupted - remove buttons.


### GUI libraries
require(gWidgets)
require(gWidgetstcltk)
options("guiToolkit"="tcltk")

### Common widget
confirmDialog <- function(LABEL="Confirm",message, handler=function(h,...)quit(save = "no" ))
{
 window <- gwindow(LABEL)
 group <- ggroup(container = window)

 ## A group for the message and buttons
 inner.group <- ggroup(horizontal=FALSE, container = group)
 glabel(message, container=inner.group, expand=TRUE)

 ## A group to organize the buttons
 button.group <- ggroup(container = inner.group)
 ## Push buttons to right
 addSpring(button.group)
 gbutton("Abort", handler=handler, container=button.group)
 gbutton("Continue", handler = function(h,...) dispose(window), container=button.group)
return()
}



### Directory structure
# "C:/Documents and Settings/msalm/My Documents/scripts_0110/PFIDO"
root<-getwd();
root<-unlist(strsplit(root, "/"));
root<-paste(root[1:which(root == "PFIDO_gui")], sep="/", collapse="/");
setwd(root);

# The starting paramter set...
temp=c("PED", "","MAP","","REFERENCE","CEU","HWE","0.001","MAF","0.01","call.rate","0.9","mind","0.9","DIM","4","THRESHOLD","1","plot","Y","plot.name","mds_plot.pdf","out.nom","MP1_output.txt","restrict","N","snp.ids","NA","MODEL","E,V","OUTLIER","Y","VALIDATE","Y","override.BIC","Y")
params=matrix(temp, ncol=2, byrow=T)
rownames(params)=params[,1]


##################
### Containers ###
##################
# The main window
win <- gwindow("PFIDO")
nb = gnotebook(container=win, expand=F)# Notebook would be ideal
MAIN = ggroup(horizontal=F, container=nb, label="Main")
ADV = ggroup(container=nb, label="Options")

# The "sub-windows" (containers within containers within containers...)
pedgrp = ggroup(horizontal = T, spacing = 10, container = MAIN)
mapgrp = ggroup(horizontal = T, spacing = 10, container=MAIN)
hapgrp = ggroup(horizontal = F, spacing = 10, container =MAIN)
outgrp = ggroup(horizontal =T, spacing =10, container=MAIN)
grp_name <- ggroup(horizontal = F, spacing = 10, container=MAIN)
qc = ggroup(horizontal=F, spacing =10, container=ADV, anchor=c(-1,-1)) # N.B. order of groups matters
  new_label= ggroup(horizontal =T, spacing =10, container=qc, anchor=c(-1,-1))
  gmaf= ggroup(horizontal =T, spacing =10, container=qc, anchor=c(-1,-1))
  ghwe= ggroup(horizontal =T, spacing =10, container=qc, anchor=c(-1,-1))
  gmind=ggroup(horizontal =T, spacing =10, container=qc, anchor=c(-1,-1))
  gcall=ggroup(horizontal =T, spacing =10, container=qc, anchor=c(-1,-1))
  new_label2= ggroup(horizontal =T, spacing =10, container=qc, anchor=c(-1,-1))
  gthresh = ggroup(horizontal =T, spacing =10, container=qc, anchor=c(-1,-1))
  ckbox = ggroup(horizontal=F, container=qc, anchor=c(-1,1))
  
##################################
### Fill the sub-windows: MAIN ###
##################################
ped_name <- glabel("PED file: ", container = pedgrp)
ped_in = gfilebrowse(text = "Select a PED file...", type = "open", quote = TRUE, container = pedgrp)
map_name <- glabel("MAP file: ", container = mapgrp)
map_in = gfilebrowse (text = "Select a MAP file...", type = "open", quote = TRUE, container = mapgrp)
hapmap_pops <- c("CEU","YRI","JPT")
glabel("HapMap Reference population:", cont=hapgrp)
user_pop <- gradio(hapmap_pops, container=hapgrp)
glabel("Output filename:", cont=outgrp)
out_nom = gedit("output.txt", container=outgrp)

#################################
### Fill the sub-windows: ADV ###
#################################
glabel("Data inclusion thresholds:", container = new_label)
glabel("Exclude samples missing (%) <", container=gmind)
umind=gedit("90", container=gmind)
glabel("Exclude SNPs missing (%) <", container=gcall)
ucall=gedit("90", container=gcall)
glabel("MAF >", container=gmaf)
umaf=gedit("0.01", container=gmaf)
glabel("HWE, p<", container=ghwe)
uhwe=gedit("0.001", container=ghwe)
glabel("PFIDO run parameters:", container = new_label2)
glabel("Call threshold, p <", container=gthresh)
uthresh=gedit("0.05",container=gthresh)
uoutlier=gcheckbox("Remove putative outliers?", checked=T, container=ckbox)
uoptsnp=gcheckbox("Use optimised SNP subsets?", checked=F, container=ckbox)

######################
### Action buttons ###
######################
out=gbutton("Start", cont = grp_name, handler = function(h,...)
                                              {
                                                params["PED",2]=svalue(ped_in)
                                                params["MAP",2]=svalue(map_in)
                                                params["REFERENCE",2]=svalue(user_pop)
                                                params["out.nom",2]=svalue(out_nom)
                                                params["HWE",2]=svalue(uhwe)
                                                params["MAF",2]=svalue(umaf)
                                                if(svalue(uoutlier))
                                                {
                                                  params["OUTLIER",2]="Y"
                                                }
                                                if(svalue(uoptsnp))
                                                {
                                                  params["restrict",2]="Y"
                                                  params["snp.ids",2]="snpsFILE.txt"
                                                }
                                                params["mind",2]=as.numeric(svalue(umind))/100
                                                params["call.rate",2]=as.numeric(svalue(ucall))/100
                                                params["THRESHOLD",2] = as.numeric(svalue(uthresh))
                                                write.table(params,file="tmp/parameters.txt", col.names=F, row.names=F, quote=F, sep="=")
                                                print(params)
                                                source(paste(root,"/scripts/PFIDO_src.R", collapse="", sep=""))
                                              })

### A status bar...OR a window that updates and identifies program stage

