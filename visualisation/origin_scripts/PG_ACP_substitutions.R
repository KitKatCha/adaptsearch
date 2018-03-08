# Lines to keep for ACP of substitutions
lines2keep <- function(firstline,lastline){

  lengthentry=15
  numberline=lastline-firstline
  iterations=numberline/lengthentry+1

  codons=firstline+5
  aa=firstline+9
  classif=firstline+13

  i=2
  while (i<=iterations){
    codons=c(codons,tail(codons,n=1)+lengthentry)
    aa=c(aa,tail(aa,n=1)+lengthentry)
    classif=c(classif,tail(classif,n=1)+lengthentry)
    i=i+1
  }

  result=rbind(codons,aa,classif)
  return(result)

}

# Create data table for ACP of subsitutions
readcountbias <- function(filename,firstline,lastline,type){

  lines2read=lines2keep(firstline,lastline)-1

  table=c()
  for (i in 1:length(lines2read[1,])){
    table=rbind(table,read.table(filename,skip=lines2read[type,i]-1,nrows=1))
  }


  if (type==1){
    number=61 ## AA
    labbels=read.table(filename,skip=firstline+2,nrows=1)
  }
  if (type==2){
    number=20  # codons
    labbels=read.table(filename,skip=firstline+6,nrows=1)
  }
  if (type==3){
    number=4 # AA types
    labbels=read.table(filename,skip=firstline+10,nrows=1)
  }


  result=c()

  for (k in 1:length(table[,1])){
    line=c()
    for (i in 1:number){
      for (j in 1:number){
        if (j<=i){
          line=c(line,NA)
        }
        if (j>i){
          line=c(line,table[k,(i-1)*number+j])
        }
      }
    }
    result=rbind(result,line)
  }

  table=result

  labbels=unlist(labbels,use.names=FALSE)
  colnames(table)=labbels
  table <- table[,colSums(is.na(table))<nrow(table)]
  return(table)

}

# function to create a circle
circle <- function(center = c(0, 0), npoints = 100) {
    r = 1
    tt = seq(0, 2 * pi, length = npoints)
    xx = center[1] + r * cos(tt)
    yy = center[1] + r * sin(tt)
    return(data.frame(x = xx, y = yy))
}

filename,firstline,lastline,type

species=c("fi","se","cl","wi","pl","br","tr","ex","ja","re","cr","fu","sp")
coast=c(0,0,1,0,0,0,0,1,0,0,1,1,0)
temp=c(2,2,1,1,1,1,1,1,0,0,0,0,0)
depth=c()
colours=c()
namesf=c()

# Color code
for (i in 1:(length(species)-1)){
  for (j in (i+1):length(species)){
    namesf=c(namesf,paste(species[j],species[i],sep=">"))
    if (temp[i]==temp[j]){
      colours=c(colours,"grey50")
    }
    if (temp[i]>temp[j]){
      colours=c(colours,"red")
    }
    if (coast[i]>coast[j]){
      depth=c(depth,"SeaGreen")
    }
    if (coast[i]<coast[j]){
      depth=c(depth,"DeepSkyBlue")
    }
    if (coast[i]==coast[j]){
      depth=c(depth,"grey50")
    }
  }
}

library(ggplot2)
library(ade4)

table=readcountbias(filename,firstline,lastline,type)
rownames(table)=namesf
pca=dudi.pca(table,nf=5,scannf=FALSE,center=TRUE,scale=TRUE)
# create data frame with scores
scores = as.data.frame(pca$li)

corcir = circle(c(0, 0), npoints = 100)

# create data frame with correlations between variables and PCs
correlations = as.data.frame(cor(table, pca$li))

# data frame with arrows coordinates
arrows = data.frame(x1=rep(0,length(correlations$Axis1)), y1=rep(0,length(correlations$Axis1)), x2=correlations$Axis1, y2=correlations$Axis2, l=rownames(correlations))
rownames(arrows)=rownames(correlations)
arrows = arrows[rowSums(is.na(arrows))<(ncol(arrows)-3),]
arrowsx=arrows[order(abs(arrows$x2)),]
arrowsy=arrows[order(abs(arrows$y2)),]
arrowsx=tail(arrowsx,5)
arrowsy=tail(arrowsy,5)
arrows=rbind(arrowsx,arrowsy)

# geom_path will do open circles 
ggplot() +geom_point(data = scores, aes(x = Axis1, y = Axis2, label = rownames(scores)), colour=depth, size=4) + geom_text(data= scores, aes(x= Axis1, y=Axis2, label=rownames(scores)), colour=colours, hjust=-0.15, vjust=-0.1)+ 
  geom_hline(yintercept = 0, colour = "gray65") +
  geom_vline(xintercept = 0, colour = "gray65") +
  geom_text(colour = "tomato", alpha = 0.8, size = 4) +
  ggtitle("PCA plot of Abundance") +
geom_path(data = corcir, aes(x = x, y = y), colour = "gray65") + 
    geom_segment(data = arrows, aes(x = x1, y = y1, xend = x2, yend = y2), colour = "gray65") + 
    geom_text(data = arrows, aes(x = x2, y = y2, label = arrows$l)) + theme_bw()
