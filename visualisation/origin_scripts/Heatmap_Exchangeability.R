#
lines2keepbias <- function(firstline,lastline){

  length1entry=18
  numberline=lastline-firstline
  iterations=numberline/length1entry+1

  codons=firstline+5
  aa=firstline+10
  classif=firstline+15

  i=2
  while (i<=iterations){
    codons=c(codons,tail(codons,n=1)+length1entry)
    aa=c(aa,tail(aa,n=1)+length1entry)
    classif=c(classif,tail(classif,n=1)+length1entry)
  i=i+1
  }

  result=rbind(codons,aa,classif)
  return(result)

}

readcountbias <- function(filename,line,type){

  line=line-1

  table=read.table(filename,skip=line-1,nrows=1)

  result=table
  names(result)=NA


  if (type==1){
    number=61
    labbels=read.table(filename,skip=6,nrows=1)
  }
  if (type==2){
    number=20
    labbels=read.table(filename,skip=11,nrows=1)
  }
  if (type==3){
    number=4
    labbels=read.table(filename,skip=16,nrows=1)
  }


  resultf=c()


  for (i in 1:number){
    a=(i-1)*number+1
    b=i*number
    resultf=rbind(resultf,result[a:b])
  }

  labbels=unlist(labbels,use.names=FALSE)
  rownames(resultf)=labbels
  colnames(resultf)=labbels

  return(resultf)

}

lines2keepcount <- function(firstline,lastline){

  lengthentry=21
  numberline=lastline-firstline
  iterations=numberline/lengthentry+1

  codons=c(firstline+5)
  aa=c(firstline+10)
  classif=c(firstline+15)

  i=2
  while (i<=iterations){
    codons=c(codons,(i-1)*lengthentry+firstline+5)
    aa=c(aa,(i-1)*lengthentry+firstline+10)
    classif=c(classif,(i-1)*lengthentry+firstline+15)
    i=i+1
  }

  result=rbind(codons,aa,classif)
  return(result)

}

readcount <- function(filename,line,type){

  line=line-1

  table=read.table(filename,skip=line-1,nrows=1)

  result=table

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

  names(resultf)=labbels

  return(resultf)

}

filename,firstlinebias, lastlinebias, firstlinecount, lastlinecount,type

table=c()
lines2readcount=lines2keepcount(firstlinecount,lastlinecount)
lines2readbias=lines2keepbias(firstlinebias,lastlinebias)
iterations=length(lines2readbias[1,])/length(lines2readcount[1,])

for (i in 1:iterations){
  for (j in 1:length(lines2readcount[1,])){
    count=readcount(filename,lines2readcount[type,j],type)
    bias=readcountbias(filename,lines2readbias[type,j+(i-1)*length(lines2readcount[1,])],type)

    for (k in 1:length(bias[1,])){
      for (l in 1:length(bias[,1])){
        bias[l,k]=bias[l,k]*count[k]
      }
    }
    gain=colSums(bias)
    loss=rowSums(bias)
    result=(gain-loss)*1000
    table=rbind(table,result)
  }
}

tablemean=c()
tablepairs=c()
for (i in 1:length(table[1,])){
  tablemean=cbind(tablemean,mean(table[,i]))
  valuecolumn=0
  for (j in 1:length(table[,1])){
    if (table[j,i]>0){
      valuecolumn=valuecolumn+1
    }
  }
  tablepairs=cbind(tablepairs,valuecolumn*100/length(table[,1]))
}
tablemean=as.vector(tablemean)
tablepairs=as.vector(tablepairs)
names(tablepairs)=colnames(table)
names(tablemean)=colnames(table)

library(ggplot2)
library(scales)

if (type=="codons"){
  type=1
}
if (type=="aa"){
  type=2
}
if (type=="classif"){
  type=3
}

namesofcols=names(tablemean)

colv=c()
rowv=c()
l=length(tablemean)
circlefill=c()
circlesize=c()

for (j in 1:length(tablemean)){

  circlefill=c(circlefill,tablepairs[j])
  circlesize=c(circlesize,tablemean[j])

  colv=c(colv,namesofcols[j])
  rowv=c(rowv,1)

}

result=cbind(rowv, colv, circlesize, circlefill)
circlefill2=circlefill
circlesize2=circlesize

result=data.frame(result)

require(ggplot2)
heatmap <- ggplot(result, aes(y=factor(rowv), x=factor(colv))) +
    xlab(" ") + 
    ylab(" ") + 
    theme_bw() +
    theme(axis.ticks=element_blank(), axis.text.x=element_text(size=10, angle=270, hjust=0, colour="black"), axis.text.y=element_text(size=0, angle=270, hjust=0, colour="black"))

heatmap + geom_point(aes(colour=circlefill2, size=abs(circlesize2)), na.rm=TRUE) + 
    scale_size("absolute gain or loss", range = c(1,10)) + 
    scale_color_gradientn("propotion of gaining/losing pairs", colours=c("blue","cyan","grey80", "yellow","red"), values=rescale(c(0,49,50,51,100)), limits=c(0,100), na.value="grey80")