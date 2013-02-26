irisMclust <- try(Mclust(iris[,-5]))
test<-grep("Error", irisMclust)
if(length(test) > 0){stop("A quick test on iris dataset indicates that mclust is not working\n")}
