#
vieuxlines2keep <- function(firstline,lastline){

  lengthentry=21
  numberline=lastline-firstline
  iterations=numberline/lengthentry+1

  codons=c(firstline,firstline+5)
  aa=c(firstline,firstline+10)
  classif=c(firstline,firstline+15)

  i=2
  while (i<=iterations){
    codons=c(codons,(i-1)*lengthentry+firstline,(i-1)*lengthentry+firstline+5)
    aa=c(aa,(i-1)*lengthentry+firstline,(i-1)*lengthentry+firstline+10)
    classif=c(classif,(i-1)*lengthentry+firstline,(i-1)*lengthentry+firstline+15)
    i=i+1
  }

  result=rbind(codons,aa,classif)
  return(result)

}

readpvalbias <- function(filename,firstline,lastline,type){

  lines2read=lines2keep(firstline,lastline)

  table=c()
  for (i in 1:length(lines2read[1,])){
    if (i%%2==1){
      title=read.table(filename,skip=lines2read[type,i]-1,nrows=1)
    }

    if (i%%2==0){
      line=cbind(title,read.table(filename,skip=lines2read[type,i]-1,nrows=1))
      table=rbind(table,line)
    }
  }


  result=c()

  for (j in 5:length(table[1,])){
    result=c(result,0)
    for (i in 1:length(table[,2])){
      if (table[i,j]>0.95){
        result[length(result)]=result[length(result)]+1
      }
      if (table[i,j]<0.05){
        result[length(result)]=result[length(result)]-1
      }
    }
  }


  if (type==1){
    labbels=read.table(filename,skip=6,nrows=1)
  }
  if (type==2){
    labbels=read.table(filename,skip=11,nrows=1)
  }
  if (type==3){
    labbels=read.table(filename,skip=16,nrows=1)
  }


  resultf=result


  labbels=unlist(labbels,use.names=FALSE)
  resultf=replace(resultf, resultf==0, NA)

  names(resultf)=labbels

  return(resultf)

}


readcountbias <- function(filename,firstline,lastline,type){

  lines2read=lines2keep(firstline,lastline)-1


  table=c()
  for (i in 1:length(lines2read[1,])){
    if (i%%2==1){
      title=read.table(filename,skip=lines2read[type,i],nrows=1)
    }
    if (i%%2==0){
      line=cbind(title,read.table(filename,skip=lines2read[type,i]-1,nrows=1))
      table=rbind(table,line)
    }
  }


  result=c()

  for (i in 4:length(table[1,])){
    result=c(result,mean(table[,i]))
  }


  if (type==1){
    labbels=read.table(filename,skip=6,nrows=1)
  }
  if (type==2){
    labbels=read.table(filename,skip=11,nrows=1)
  }
  if (type==3){
    labbels=read.table(filename,skip=16,nrows=1)
  }


  resultf=result



  labbels=unlist(labbels,use.names=FALSE)
  resultf=replace(resultf, resultf==0, NA)

  names(resultf)=labbels

  return(resultf)

}

filename,firstline,lastline,type

library(ggplot2)
library(scales)

iterations=((lastline-firstline)/21)+1

if (type=="codons"){
  type=1
}
if (type=="aa"){
  type=2
}
if (type=="classif"){
  type=3
}

value=readpvalbias(filename,firstline,lastline,type)
count=readcountbias(filename,firstline,lastline,type)*100

namesofcols=names(count)

colv=c()
rowv=c()
l=length(value)
circlefill=c()
circlesize=c()

for (j in 1:length(value)){

  circlefill=c(circlefill,value[j])
  circlesize=c(circlesize,count[j])

  colv=c(colv,namesofcols[j])
  rowv=c(rowv,1)

}

result=cbind(rowv, colv, circlesize, circlefill)
circlefill2=circlefill
circlesize2=circlesize

result=data.frame(result)

require(ggplot2)
heatmap <-  ggplot(result, aes(y=factor(rowv),  x=factor(colv))) +
    xlab(" ") +
    ylab(" ") +
    theme_bw() +
    theme(axis.ticks=element_blank(), axis.text.x=element_text(size=10, angle=270, hjust=0, colour="black"),axis.text.y=element_text(size=0, angle=270, hjust=0, colour="black"))

heatmap + geom_point(aes(colour=circlefill2, size=circlesize2), na.rm=TRUE) + 
    scale_size("relative abundance", range=c(65/l, 260/l)) +
    scale_color_gradientn("significant tests", colours=c("blue","cyan","grey80", "yellow","red"), values=rescale(c(-iterations,-1,0,1,iterations)), limits=c(-iterations,iterations), na.value="grey80")