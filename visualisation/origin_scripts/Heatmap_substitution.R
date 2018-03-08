# Lines to keep for heatmaps of substitutions
lines2keep <- function(firstline,lastline){

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

readpvalbias <- function(filename,firstline,lastline,type){

  lines2read=lines2keep(firstline,lastline)

  table=c()
  for (i in 1:length(lines2read[1,])){
    table=rbind(table,read.table(filename,skip=lines2read[type,i]-1,nrows=1))
  }

  result=c()

  for (j in 2:length(table[1,])){
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
  resultf=replace(resultf, resultf==0, NA)
  rownames(resultf)=labbels
  colnames(resultf)=labbels

  return(resultf)

}


readcountbias <- function(filename,firstline,lastline,type){

  lines2read=lines2keep(firstline,lastline)-1


  table=c()
  for (i in 1:length(lines2read[1,])){
    table=rbind(table,read.table(filename,skip=lines2read[type,i]-1,nrows=1))
  }

  result=c()

  for (i in 1:length(table[1,])){
    result=c(result,mean(table[,i]))
  }


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
  resultf=replace(resultf, resultf==0, NA)
  rownames(resultf)=labbels
  colnames(resultf)=labbels

  return(resultf)

}

filename,firstline,lastline,type

library(ggplot2)
library(scales)

iterations=((lastline-firstline)/18)+1

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

namesofrows=rownames(count)
namesofcols=colnames(count)

rowv=c()
colv=c()

l=length(value[1,])
circlefill=c()
circlesize=c()

for (i in 1:length(value[1,])){

  for (j in 1:length(value[,1])){

    circlefill=c(circlefill,value[j,i])
    circlesize=c(circlesize,count[j,i])
    rowv=c(rowv,namesofrows[j])
    colv=c(colv,namesofcols[i])

  }
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
    theme(axis.ticks=element_blank(), axis.text.x=element_text(size=10, angle=270, hjust=0, colour="black"))

heatmap + geom_point(aes(colour=circlefill2, size =circlesize2), na.rm=TRUE) +
    scale_size("proportion of mutation", range = c(80/l, 240/l)) +
    scale_color_gradientn("significant tests", colours=c("blue","cyan", "grey80", "yellow","red"), values=rescale(c(-iterations,-1,0,1,iterations)), limits=c(-iterations,iterations), na.value="grey80")