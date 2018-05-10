library(Hmisc)
setwd("/Users/jtfield/git-repos/biostats_data_DO_NOT_PUSH")

d = read.table("chr21-subset.vcf", 
               sep="\t", 
               fill=FALSE, 
               strip.white=TRUE,stringsAsFactors=FALSE)
head(d[1:15],n=15)
str(d)


new.d=d[c(-1:-9)]
str(new.d)

new.d[sapply(new.d, is.factor)] = lapply(new.d[sapply(new.d, is.factor)], as.character)

new.d[new.d=="0|0"]<-0
new.d[new.d=="0|1"]<-1
new.d[new.d=="1|0"]<-1
new.d[new.d=="1|1"]<-2

indx <- sapply(new.d, as.numeric)
# at this point, the dataframe is read in, the values have been replaced by numbers
# and the dataframe has been converted to a matrix of numbers (as.numeric)



model = prcomp(indx, scale=TRUE)

model$rotation
model$scale

par(mfrow=c(2,2))
plot(model$x[,1], col=indx[,5])
plot(model$x[,2], col=indx[,5])
plot(model$x[,3], col=indx[,5])
plot(model$x[,4], col=indx[,5])




head(new.d[1:15],n=15)

