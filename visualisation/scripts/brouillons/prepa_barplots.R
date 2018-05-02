dat <- data.frame(
  group = rep(c("Above", "Below"), each=10),
  x = rep(1:10, 2),
  y = c(runif(10, 0, 1), runif(10, -1, 0))
)

library(ggplot2)
ggplot(dat, aes(x=x, y=y, fill=group)) + 
  geom_bar(stat="identity", position="identity")

dat <- data.frame(
  transi = c("aromatic->aromatic","polar->aromatics","unpolar->aromatics","charged->aromatics",
            "aromatics->polar","polar->polar","unpolar->polar","charged->polar","aromatics->unpolar",
            "polar->unpolar","unpolar->unpolar","charged->unpolar","aromatics->charged",
            "polar->charged","unpolar->charged","charged->charged"),
  freq = c(0.000000,0.047619,0.021739,0.000000,0.823223,0.000000,0.021739,0.045455,0.990157,-1.000000,
      0.000000,0.000000,0.000000,-1.000000,0.000000,0.000000),
  pval = c(0.495050,0.663366,0.445545,0.495050,0.613861,0.495050,0.198020,0.287129,0.623762,0.128713,
           0.495050,0.247525,0.495050,0.247525,0.613861,0.495050)
)
newdat <- dat[!dat$freq == 0,]

ggplot(newdat, aes(x=reorder(transi, -freq), y=freq)) + 
  geom_bar(stat="identity", width = 0.25, color = "black", fill="steelblue") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# ----------------------------------------------------------------------------------------------------------

file <- read.table("forbarplot.txt",sep=" ", dec=".", header=TRUE)
t <- t(file)
ggplot(newdat, aes(x=t[,0] y=t[,1])) + 
  geom_bar(stat="identity", width = 0.25, color = "black", fill="steelblue") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
