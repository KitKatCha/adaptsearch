library(rbokeh)

setwd("~/Documents/Fork_AdaptSearch/adaptsearch/galaxy_wrappers/07_MutCount/test-data/OUT_concat")
data <- read.table("aatypes_freqs.csv", header=TRUE, dec=".", sep=",")
data_1 <- data[seq(2, nrow(data), 3),]
data_2 <- data[seq(3, nrow(data), 3),]

i <- figure(width=600, heigth=400) %>%
  ly_bar(x=Species, y=aromatics, data=data_1, hover=TRUE) %>%
  ly_bar(x=Species, y=charged, data=data_1, hover=TRUE) %>%
  ly_bar(x=Species, y=polar, data=data_1, hover=TRUE) %>%
  ly_bar(x=Species, y=unpolar, data=data_1, hover=TRUE) %>%
  theme_axis("x", major_label_orientation = 90) %>%
  y_axis(label="Frequencies")

f <- figure(width=300, heigth=200) %>%
  ly_bar(x=colnames(data_1)[2:5], y=as.numeric(data_1[1,][2:5]), 
         data=data_1, hover=TRUE) %>%
  ly_points(x=colnames(data_2)[2:5], y=as.numeric(data_2[1,][2:5]), 
            data=data_2, color="red") %>%
  y_axis(label="Frequencies and pvalues") %>%
  x_axis(label="Amino-acides types")

data_aa <- read.table("aa_freqs.csv", header=TRUE, dec=".", sep=",")
data_aa_1 <- data_aa[seq(2, nrow(data), 3),]
data_aa_2 <- data_aa[seq(3, nrow(data), 3),]

g <- figure(width=1000, heigth=200) %>%
  ly_bar(x=colnames(data_aa_1)[2:21], y=as.numeric(data_aa_1[1,][2:21]), 
         data=data_aa_1, hover=TRUE) %>%
  ly_points(x=colnames(data_aa_2)[2:21], y=as.numeric(data_aa_2[1,][2:21]), 
             data=data_aa_2) %>%
  y_axis(label="Frequencies and pvalues") %>%
  x_axis(label="Amino-acids")
