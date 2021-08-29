#MGrofit v2.0
#JV 15.10.2014


library("grofit")
library("cluster")
library("gclus")
library("genefilter")
library(xlsx)

# error bar function
error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
  if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}

MGrofit=function(fileformat="standard",
                 file,
                 fit="s",
                 model=c("logistic", "richards","gompertz", "gompertz.exp"),
                 Interactive=F,
                 writeTable=T,
                 exclude=0,
                 makeMean=T,
                 makePlot=T){

# get the data
  if(fileformat == "standard"){
    data=read.csv(file)
  }
  if(fileformat == "tecan"){
    t <- readline("export tecan out output xlsx into csv file [ENTER]")
    data=read.csv(file, skip=34)
  }  

# formatting the data
  if(fileformat == "standard"){
    time=matrix(ncol=length(data[1,]), nrow=length(data[,1])-1)
    
    for (x in 2:length(data[,1])){
     time[x-1,]=as.numeric(data[1,])
     }
    
    time=time[,-c(1:3)]
    data=data[-1,]
  }
  if(fileformat == "tecan"){
    data <- cbind(Content=data[,1], Well.Col=data[,1], Well.Row=data[,1], data[,2:ncol(data)])
    data <- data[-2,]
    data <- data[!is.na(data[,6]),]
    
    time=matrix(ncol=length(data[1,]), nrow=length(data[,1])-1)
    
    for (x in 2:length(data[,1])){
      time[x-1,]=as.numeric(data[1,])
    }
    
    time=time[,-c(1:3)]
    data=data[-1,]
    
    levels(data$Well.Col) <- sub("\\D", "", levels(data$Well.Col), perl=T)
    for(i in 1:10){
    levels(data$Well.Row) <- sub("\\d", "", levels(data$Well.Row), perl=T)
    } 
  }  

# exclude initial data points... sometimes allows then to fit data
  if (exclude>0){
    data=data[,-c(4:(3+exclude))]
    time=time[,-c(1:(0+exclude))]
  }

# get rid of curves with less than 7 datapoints
  for (i in length(data[,1]):1){
    if (sum(!is.na(data[i,]))<10){
      data=data[-i,]
      time=time[-i,]
    }
  }

#change header type to character  
  data[,1] <- as.character(data[,1])
  data[,2] <- as.character(data[,2])
  data[,3] <- as.character(data[,3])
  
# chose the fit and the model and start grofit
  if (fit=="s"){
    b <- gcFit(time, data, grofit.control(fit.opt=fit, interactive=Interactive))
    a=summary(gcFit(time, data, grofit.control(fit.opt=fit, interactive=Interactive))) 
     }
  if (fit=="b"){
    b <- gcFit(time, data, grofit.control(fit.opt=fit,model.type=c(model), interactive=Interactive))
    a=summary(gcFit(time, data, grofit.control(fit.opt=fit,model.type=c(model), interactive=Interactive))) 
  }
  if (fit=="m"){
    b <- gcFit(time, data, grofit.control(fit.opt="m",model.type=c(model), interactive=Interactive))
    a=summary(gcFit(time, data, grofit.control(fit.opt="m",model.type=c(model), interactive=Interactive)))
    }
colnames(a)[1:3]=c("Sample", "Column", "Row")

# plot individual fits into pdf
  if(makePlot==T){
    filepdf <- sub(".csv", ".pdf", file, fixed=T)
    if(fit=="s"){
      output <- paste("Grofit", "spline", "plot", filepdf, sep="_")
    }
    if (fit=="b"){ 
      if(length(model)>1){
        output <- paste("Grofit", "splineNmodel", filepdf, sep="_") }
      else{ output <- paste("Grofit","splineNmodel", model, file, sep="_")
      }  
    }
    if (fit=="m"){
      if(length(model)>1){
        output <- paste("Grofit", "model", file, sep="_")}
      else{ output <- paste("Grofit","model", model, file, sep="_")  
      }
    }  

    pdf(output)
    plot(b, opt=fit, raw=T, slope=T, pch=1:7, cex=1, colSpline=1, colData=1:7)
    dev.off()
  }
  
  
# get rid of columns with only NA  
  ind <- colSums(is.na(a)) != nrow(a)
  a=a[, ind]

# get rid of columns that you dont trust
  for (u in length(a[,1]):1){
    if (a[u,4]==FALSE){
      a=a[-u,]
    }  
  }  

# means with spline fit
  if (fit=="s"& makeMean==T){
    a2=a[,-c(1:7)]
    Data=matrix(NA, nrow=length(a2[1,]), ncol=length(a2[,1]))
    
    for (z in 1:length(a2[1,])){
      Data[z,]=a2[,z]
    }
    colnames(Data)=a[,1]
    row.names(Data)=colnames(a2)
    
    Mean=sapply(unique(colnames(Data)), function(x) rowMeans( Data[ , which(x==colnames(Data)), drop=FALSE], na.rm=TRUE) )
    Sds=sapply(unique(colnames(Data)),function(x) rowSds( Data[ , which(x==colnames(Data)), drop=FALSE], na.rm=TRUE) )
    
    Result=matrix(NA, nrow=length(Mean[1,]), ncol=1+2*(length(Mean[,1])))
    Result[,1]=colnames(Mean)
    
    coln=c("sample")
    for (x in 1:length(Mean[,1])){
      Result[,2*x]=Mean[x,]
      Result[,2*x+1]=Sds[x,]
      coln=c(coln, paste("Mean",row.names(Mean)[x], sep="_"), paste("Stdev", row.names(Sds)[x], sep="_"))
    }
    colnames(Result)=coln
    
    a=Result
    
}
  
# means for all other fits/models...  take care no propagation of the error from model
  if (fit!="s" & makeMean==T){
    
    a2=a[,-c(1:8)]
    Data=matrix(NA, nrow=length(a2[1,]), ncol=length(a2[,1]))
    
    for (z in 1:length(a2[1,])){
      Data[z,]=a2[,z]
    }
    colnames(Data)=a[,1]
    row.names(Data)=colnames(a2)
    
    Mean=sapply(unique(colnames(Data)), function(x) rowMeans( Data[ , which(x==colnames(Data)), drop=FALSE], na.rm=TRUE) )
    Sds=sapply(unique(colnames(Data)),function(x) rowSds( Data[ , which(x==colnames(Data)), drop=FALSE], na.rm=TRUE) )
    
    Result=matrix(NA, nrow=length(Mean[1,]), ncol=1+2*(length(Mean[,1])))
    Result[,1]=colnames(Mean)
    
    coln=c("sample")
    for (x in 1:length(Mean[,1])){
      Result[,2*x]=Mean[x,]
      Result[,2*x+1]=Sds[x,]
      coln=c(coln, paste("Mean",row.names(Mean)[x], sep="_"), paste("Stdev", row.names(Sds)[x], sep="_"))
    }
    colnames(Result)=coln
    
    a=Result
  }


  View(a)
  
# write the table
if (writeTable==T){
  if (fit=="s"){write.table(a,file=paste("Grofit", "spline", file, sep="_"), sep=",",row.names=FALSE)
    }
  if (fit=="b"){ 
      if(length(model)>1){
        write.table(a,file=paste("Grofit", "splineNmodel", file, sep="_"), sep=",",row.names=FALSE)}
      else{ write.table(a,file=paste("Grofit","splineNmodel", model, file, sep="_"), sep=",",row.names=FALSE)  
      }  
    }
  if (fit=="m"){
    if(length(model)>1){
      write.table(a,file=paste("Grofit", "model", file, sep="_"), sep=",",row.names=FALSE)}
    else{ write.table(a,file=paste("Grofit","model", model, file, sep="_"), sep=",",row.names=FALSE)  
      }
    }  

# plot µ, ð and A
  par(las=2) # make label text perpendicular to axis
#  par(mar=c(16,6,4,4)) # increase y-axis margin.
  
  
  barx <- barplot(as.numeric(a[,4]), names.arg=a[,1],col="red", axis.lty=1, xpd=F, ylab="h", space = 0.4, cex.names=0.8)
  error.bar(barx,as.numeric(a[,4]), as.numeric(a[,5]))
  title(main="lag phase")
  abline(h=mean(as.numeric(a[,4])), col="red")
  abline(h=c(mean(as.numeric(a[,4]))+sd(as.numeric(a[,4])), mean(as.numeric(a[,4]))-sd(as.numeric(a[,4]))), col="brown")
  
  
  barx <- barplot(as.numeric(a[,2]), names.arg=a[,1],col="blue", axis.lty=1, xpd=F, ylab=expression(Delta * OD[595]/h), space = 0.4, cex.names=0.8)
  error.bar(barx,as.numeric(a[,2]), as.numeric(a[,3]))
  title(main="growth rate")
  abline(h=mean(as.numeric(a[,2])), col="red")
  abline(h=c(mean(as.numeric(a[,2]))+sd(as.numeric(a[,2])), mean(as.numeric(a[,2]))-sd(as.numeric(a[,2]))), col="brown")
  
  barx <- barplot(as.numeric(a[,6]), names.arg=a[,1],col="yellow", axis.lty=1, xpd=F, ylab=expression(OD[595]), space = 0.4, cex.names=0.8)
  error.bar(barx,as.numeric(a[,6]), as.numeric(a[,7]))
  title(main="final OD595")
  abline(h=mean(as.numeric(a[,6])), col="red")
  abline(h=c(mean(as.numeric(a[,6]))+sd(as.numeric(a[,6])), mean(as.numeric(a[,6]))-sd(as.numeric(a[,6]))), col="brown")
  
  }
#par(las=0,mar=c(5.1, 4.1, 4.1, 2.1))
}