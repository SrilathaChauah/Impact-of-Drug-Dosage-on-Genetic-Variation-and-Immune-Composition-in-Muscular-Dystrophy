if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Biobase", version = "3.8")
library(Biobase)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("GEOquery", version = "3.8")
library(GEOquery)

install.packages("reshape2")
library(reshape2)

install.packages("ggplot2")
library(ggplot2)

install.packages("ggrepel")
library(ggrepel)



data <- read.table("GSE95682_RPKM_-_all_genes.csv", header = T, sep = ',',stringsAsFactors = F)
data <- na.omit(data)

###########################################---------PCA--------##########################################

install.packages("devtools")
library(devtools)
install_github("vqv/ggbiplot")
library(ggbiplot)
install.packages("ggplot2")
library(ggplot2)
data1 <- read.table("GSE95682.csv", header = T, sep = ',',stringsAsFactors = F)
data <- data1
data$GeneSymbol <- NULL
data2 <- log(data+1)

data3 <- data.frame(GeneSymbols = data1$GeneSymbol, data2)

data3.pca <- prcomp(t(data2), center = TRUE,scale. = TRUE)

library("GGally")
df <- data.frame(data3.pca$x[,1:4])
df$id  <- row.names(df)
df$tre  <- sapply(df$id,FUN = function(x){substr(x,1,5)})

ggpairs(data=df,columns = 1:4,aes(color=tre))

ggbiplot(data3.pca, labels=rownames(data3))

#===========
ggpairs(data=df,columns = 1:5,
        aes(color=tre, shape =tre))

ggbiplot(data3.pca, labels=rownames(data3))


barplot(data3.pca$sdev * 100/sum(data3.pca$sdev),
        las=1,
        xlab="Dimensions", 
        ylab="Proportion of explained variance (%)", y.axis=NULL,
        col="darkgrey")

####################################################################################################################


## Density plots for RNAseq data

library(meltt)
melted_data = melt(data1, id.vars = "GeneSymbol")

melted_data$group <- "NA"
melted_data$group[melted_data$variable %in% c("Vehicle1", 
                                              "Vehicle2", 
                                              "Vehicle3",
                                              "Vehicle4",
                                              "Vehicle5")] <- "untreated"

melted_data$group[melted_data$variable %in% c("Daily1", 
                                              "Daily2", 
                                              "Daily3",
                                              "Daily4",
                                              "Daily5")] <- "daily"

melted_data$group[melted_data$variable %in% c("Weekly1", 
                                              "Weekly2", 
                                              "Weekly3",
                                              "Weekly4",
                                              "Weekly5")] <- "weekly"
daily <- melted_data[melted_data$group == 'daily',]
weekly <- melted_data[melted_data$group == 'weekly',]
untreated <- melted_data[melted_data$group == 'untreated',]

nountreated <- rbind(daily,weekly)
nodaily <- rbind(weekly, untreated)

melted_data$group <- factor(melted_data$group, levels = c("untreated", "daily", "weekly"), ordered = TRUE)
library(ggplot2)
ggplot(subset(melted_data, melted_data$GeneSymbol %in% "Klf15"),
       aes(x=value,fill=group)) + geom_density() + ggtitle(label = "Melted plot for Quattrocelli's dataset", 
      subtitle = "Gene Symbol: Klf15") + theme_bw() + theme(plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)) 

ggplot(subset(melted_data, melted_data$GeneSymbol %in% "Klf15"),
       aes(x=group,y=value, fill=group)) + geom_boxplot() + geom_jitter(width = 0.25) + ggtitle(label = "Box plot for Quattrocelli's dataset", 
        subtitle = "Gene Symbol: Klf15") + theme_bw() + theme(plot.title = element_text(hjust = 0.5),
         plot.subtitle = element_text(hjust = 0.5)) 


ggplot(subset(melted_data, melted_data$GeneSymbol %in% "Dmd"),
       aes(x=value,fill=group)) + geom_density() + ggtitle(label = "Melted plot for Quattrocelli's dataset", 
       subtitle = "Gene Symbol: Dmd") + theme_bw() + theme(plot.title = element_text(hjust = 0.5),
       plot.subtitle = element_text(hjust = 0.5)) 


ggplot(subset(nountreated, nountreated$GeneSymbol %in% "Mef2a"),
       aes(x=value,fill=group)) + geom_density() + ggtitle(label = "Melted plot for Quattrocelli's dataset", 
       subtitle = "Gene Symbol: Mef2a") + theme_bw() + theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) 

ggplot(subset(nountreated, nountreated$GeneSymbol %in% "Anxa6"),
       aes(x=value,fill=group)) + geom_density() + ggtitle(label = "Melted plot for Quattrocelli's dataset", 
       subtitle = "Gene Symbol: Anxa6") + theme_bw() + theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) 

ggplot(subset(nodaily, nodaily$GeneSymbol %in% "Fbxo32"),
       aes(x=value,fill=group)) + geom_density() + ggtitle(label = "Melted plot for Quattrocelli's dataset", 
       subtitle = "Gene Symbol: Fbxo32") + theme_bw() + theme(plot.title = element_text(hjust = 0.5),
       plot.subtitle = element_text(hjust = 0.5)) 

ggplot(subset(melted_data, melted_data$GeneSymbol %in% "Cxcl12"),
       aes(x=value,fill=group)) + geom_density() + ggtitle(label = "Melted plot for Quattrocelli's dataset", 
       subtitle = "Gene Symbol: Cxcl12") + theme_bw() + theme(plot.title = element_text(hjust = 0.5),
       plot.subtitle = element_text(hjust = 0.5)) 

ggplot(subset(nountreated, nountreated$GeneSymbol %in% "Ccl7"),
       aes(x=value,fill=group)) + geom_bin2d() + ggtitle(label = "Melted plot for Quattrocelli's dataset", 
       subtitle = "Gene Symbol: Ccl7") + theme_bw() + theme(plot.title = element_text(hjust = 0.5),
       plot.subtitle = element_text(hjust = 0.5)) 

counts <- table(melted_data$GeneSymbol %in% "Ccl7")
barplot(counts, main="melteddata.csv", horiz=TRUE, names.arg=c("Daily1", "Weekly1", "Vehicle1"), cex.names=0.8)
#========================================================================================
##box plot

ggplot(subset(melted_data, melted_data$GeneSymbol %in% "Klf15"),
       aes(x=group,y=value, fill=group)) + geom_boxplot() + geom_jitter(width = 0.25) + ggtitle(label = "Box plot for Quattrocelli's dataset", 
                                                                                                subtitle = "Gene Symbol: Klf15") + theme_bw() + theme(plot.title = element_text(hjust = 0.5),
                                                                                                                                                      plot.subtitle = element_text(hjust = 0.5)) 


ggplot(subset(melted_data, melted_data$GeneSymbol %in% "Dmd"),
       aes(x=group,y=value, fill=group)) + geom_boxplot() + geom_jitter(width = 0.25) + ggtitle(label = "Box plot for Quattrocelli's dataset", 
                                                                                                subtitle = "Gene Symbol: Dmd") + theme_bw() + theme(plot.title = element_text(hjust = 0.5),
                                                                                                                                                      plot.subtitle = element_text(hjust = 0.5)) 


ggplot(subset(nountreated, nountreated$GeneSymbol %in% "Mef2a"),
       aes(x=group,y=value, fill=group)) + geom_boxplot() + geom_jitter(width = 0.25) + ggtitle(label = "Box plot for Quattrocelli's dataset", 
                                                                                                subtitle = "Gene Symbol: Mef2a") + theme_bw() + theme(plot.title = element_text(hjust = 0.5),
                                                                                                                                                    plot.subtitle = element_text(hjust = 0.5)) 

ggplot(subset(melted_data, melted_data$GeneSymbol %in% "Anxa6"),
       aes(x=group,y=value, fill=group)) + geom_boxplot() + geom_jitter(width = 0.25) + ggtitle(label = "Box plot for Quattrocelli's dataset", 
                                                                                                subtitle = "Gene Symbol: Anxa6") + theme_bw() + theme(plot.title = element_text(hjust = 0.5),
                                                                                                                                                      plot.subtitle = element_text(hjust = 0.5)) 

ggplot(subset(nodaily, nodaily$GeneSymbol %in% "Fbxo32"),
       aes(x=group,y=value, fill=group)) + geom_boxplot() + geom_jitter(width = 0.25) + ggtitle(label = "Box plot for Quattrocelli's dataset", 
                                                                                                subtitle = "Gene Symbol: Fbxo32") + theme_bw() + theme(plot.title = element_text(hjust = 0.5),
                                                                                                                                                      plot.subtitle = element_text(hjust = 0.5)) 

ggplot(subset(nountreated, nountreated$GeneSymbol %in% "Cxcl12"),
       aes(x=group,y=value, fill=group)) + geom_boxplot() + geom_jitter(width = 0.25) + ggtitle(label = "Box plot for Quattrocelli's dataset", 
                                                                                                subtitle = "Gene Symbol: Cxcl12") + theme_bw() + theme(plot.title = element_text(hjust = 0.5),
                                                                                                                                                      plot.subtitle = element_text(hjust = 0.5)) 

ggplot(subset(nountreated, nountreated$GeneSymbol %in% "Ccl7"),
       aes(x=group,y=value, fill=group)) + geom_boxplot() + geom_jitter(width = 0.25) + ggtitle(label = "Box plot for Quattrocelli's dataset", 
                                                                                                subtitle = "Gene Symbol: Ccl7") + theme_bw() + theme(plot.title = element_text(hjust = 0.5),
                                                                                                                                                       plot.subtitle = element_text(hjust = 0.5)) 

#========================================================================================
ggplot(subset(melted_data, melted_data$GeneSymbol %in% "A730049H05Rik"),
       aes(x=value,fill=group)) + geom_density() + ggtitle(label = "Melted plot for Quattrocelli's dataset", 
                                                           subtitle = "Gene Symbol: A730049H05Rik") + theme_bw() + theme(plot.title = element_text(hjust = 0.5),
                                                                                                                 plot.subtitle = element_text(hjust = 0.5)) 

ggplot(subset(melted_data, melted_data$GeneSymbol %in% "Gm12840"),
       aes(x=value,fill=group)) + geom_density() + ggtitle(label = "Melted plot for Quattrocelli's dataset", 
                                                           subtitle = "Gene Symbol: Gm12840") + theme_bw() + theme(plot.title = element_text(hjust = 0.5),
                                                                                                                         plot.subtitle = element_text(hjust = 0.5)) 

ggplot(subset(melted_data, melted_data$GeneSymbol %in% "Fosl2"),
       aes(x=value,fill=group)) + geom_density() + ggtitle(label = "Melted plot for Quattrocelli's dataset", 
                                                           subtitle = "Gene Symbol: Fosl2") + theme_bw() + theme(plot.title = element_text(hjust = 0.5),
                                                                                                                         plot.subtitle = element_text(hjust = 0.5)) 

ggplot(subset(melted_data, melted_data$GeneSymbol %in% "Zfp36"),
       aes(x=value,fill=group)) + geom_density() + ggtitle(label = "Melted plot for Quattrocelli's dataset", 
                                                           subtitle = "Gene Symbol: Zfp36") + theme_bw() + theme(plot.title = element_text(hjust = 0.5),
                                                                                                                         plot.subtitle = element_text(hjust = 0.5)) 

ggplot(subset(melted_data, melted_data$GeneSymbol %in% "Csrnp1"),
       aes(x=value,fill=group)) + geom_density() + ggtitle(label = "Melted plot for Quattrocelli's dataset", 
                                                           subtitle = "Gene Symbol: Csrnp1") + theme_bw() + theme(plot.title = element_text(hjust = 0.5),
                                                                                                                         plot.subtitle = element_text(hjust = 0.5)) 

ggplot(subset(melted_data, melted_data$GeneSymbol %in% "Pde4b"),
       aes(x=value,fill=group)) + geom_density() + ggtitle(label = "Melted plot for Quattrocelli's dataset", 
                                                           subtitle = "Gene Symbol: Pde4b") + theme_bw() + theme(plot.title = element_text(hjust = 0.5),
                                                                                                                         plot.subtitle = element_text(hjust = 0.5)) 
ggplot(subset(melted_data, melted_data$GeneSymbol %in% "Tmem252"),
       aes(x=value,fill=group)) + geom_density() + ggtitle(label = "Melted plot for Quattrocelli's dataset", 
                                                           subtitle = "Gene Symbol: Tmem252") + theme_bw() + theme(plot.title = element_text(hjust = 0.5),
                                                                                                                 plot.subtitle = element_text(hjust = 0.5)) 

ggplot(subset(melted_data, melted_data$GeneSymbol %in% "Rapgef5"),
       aes(x=value,fill=group)) + geom_density() + ggtitle(label = "Melted plot for Quattrocelli's dataset", 
                                                           subtitle = "Gene Symbol: Rapgef5") + theme_bw() + theme(plot.title = element_text(hjust = 0.5),
                                                                                                                 plot.subtitle = element_text(hjust = 0.5)) 

ggplot(subset(melted_data, melted_data$GeneSymbol %in% "Mnda"),
       aes(x=value,fill=group)) + geom_density() + ggtitle(label = "Melted plot for Quattrocelli's dataset", 
                                                           subtitle = "Gene Symbol: Mnda") + theme_bw() + theme(plot.title = element_text(hjust = 0.5),
                                                                                                                 plot.subtitle = element_text(hjust = 0.5)) 

ggplot(subset(melted_data, melted_data$GeneSymbol %in% "Ret"),
       aes(x=value,fill=group)) + geom_density() + ggtitle(label = "Melted plot for Quattrocelli's dataset", 
                                                           subtitle = "Gene Symbol: Ret") + theme_bw() + theme(plot.title = element_text(hjust = 0.5),
                                                                                                                 plot.subtitle = element_text(hjust = 0.5)) 

ggplot(subset(melted_data, melted_data$GeneSymbol %in% "Dusp5"),
       aes(x=value,fill=group)) + geom_density() + ggtitle(label = "Melted plot for Quattrocelli's dataset", 
                                                           subtitle = "Gene Symbol: Dusp5") + theme_bw() + theme(plot.title = element_text(hjust = 0.5),
                                                                                                                 plot.subtitle = element_text(hjust = 0.5)) 

ggplot(subset(melted_data, melted_data$GeneSymbol %in% "Tmem2"),
       aes(x=value,fill=group)) + geom_density() + ggtitle(label = "Melted plot for Quattrocelli's dataset", 
                                                           subtitle = "Gene Symbol: Tmem2") + theme_bw() + theme(plot.title = element_text(hjust = 0.5),
                                                                                                                 plot.subtitle = element_text(hjust = 0.5)) 

ggplot(subset(melted_data, melted_data$GeneSymbol %in% "RP24-196H4.1"),
       aes(x=value,fill=group)) + geom_density() + ggtitle(label = "Melted plot for Quattrocelli's dataset", 
                                                           subtitle = "Gene Symbol: RP24-196H4.1") + theme_bw() + theme(plot.title = element_text(hjust = 0.5),
                                                                                                                 plot.subtitle = element_text(hjust = 0.5)) 

ggplot(subset(melted_data, melted_data$GeneSymbol %in% "Stat2"),
       aes(x=value,fill=group)) + geom_density() + ggtitle(label = "Melted plot for Quattrocelli's dataset", 
                                                           subtitle = "Gene Symbol: Stat2") + theme_bw() + theme(plot.title = element_text(hjust = 0.5),
                                                                                                                 plot.subtitle = element_text(hjust = 0.5)) 




#===========================================================================================
ggplot(subset(melted_data, melted_data$gene.y %in% "Snx9"),
       aes(x=value,fill=group)) + geom_density() + ggtitle(label = "Melted plot for Quattrocelli's dataset", 
       subtitle = "Gene Symbol: Snx9") + theme_bw() + theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) 


ggplot(subset(melted_data, melted_data$gene.y %in% "Dbp"),
       aes(x=value,fill=group)) + geom_density() + ggtitle(label = "Melted plot for Quattrocelli's dataset", 
       subtitle = "Gene Symbol: Dbp") + theme_bw() + theme(plot.title = element_text(hjust = 0.5),
       plot.subtitle = element_text(hjust = 0.5)) 

ggplot(subset(melted_data, melted_data$gene.y %in% "U6"),
       aes(x=value,fill=group)) + geom_density() + ggtitle(label = "Melted plot for Quattrocelli's dataset", 
       subtitle = "Gene Symbol: U6") + theme_bw() + theme(plot.title = element_text(hjust = 0.5),
       plot.subtitle = element_text(hjust = 0.5)) 

ggplot(subset(melted_data, melted_data$gene.y %in% "Abcc8"),
       aes(x=value,fill=group)) + geom_density() + ggtitle(label = "Melted plot for Quattrocelli's dataset", 
        subtitle = "Gene Symbol: Abcc8") + theme_bw() + theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) 

ggplot(subset(melted_data, melted_data$gene.y %in% "Arntl"),
       aes(x=value,fill=group)) + geom_density() + ggtitle(label = "Melted plot for Quattrocelli's dataset", 
       subtitle = "Gene Symbol: Arntl") + theme_bw() + theme(plot.title = element_text(hjust = 0.5),
       plot.subtitle = element_text(hjust = 0.5)) 

ggplot(subset(melted_data, melted_data$gene.y %in% "Car14"),
       aes(x=value,fill=group)) + geom_density() + ggtitle(label = "Melted plot for Quattrocelli's dataset", 
       subtitle = "Gene Symbol: Car14") + theme_bw() + theme(plot.title = element_text(hjust = 0.5),
       plot.subtitle = element_text(hjust = 0.5)) 

ggplot(subset(melted_data, melted_data$gene.y %in% "C7"),
       aes(x=value,fill=group)) + geom_density() + ggtitle(label = "Melted plot for Quattrocelli's dataset", 
       subtitle = "Gene Symbol: C7") + theme_bw() + theme(plot.title = element_text(hjust = 0.5),
       plot.subtitle = element_text(hjust = 0.5)) 

ggplot(subset(melted_data, melted_data$gene.y %in% "Agt"),
       aes(x=value,fill=group)) + geom_density() + ggtitle(label = "Melted plot for Quattrocelli's dataset", 
       subtitle = "Gene Symbol: Agt") + theme_bw() + theme(plot.title = element_text(hjust = 0.5),
       plot.subtitle = element_text(hjust = 0.5))

ggplot(subset(melted_data, melted_data$gene.y %in% "Npas2"),
       aes(x=value,fill=group)) + geom_density() + ggtitle(label = "Melted plot for Quattrocelli's dataset", 
       subtitle = "Gene Symbol: Npas2") + theme_bw() + theme(plot.title = element_text(hjust = 0.5),
       plot.subtitle = element_text(hjust = 0.5))

ggplot(subset(melted_data, melted_data$gene.y %in% "Odf3l2"),
       aes(x=value,fill=group)) + geom_density() + ggtitle(label = "Melted plot for Quattrocelli's dataset", 
       subtitle = "Gene Symbol: Odf3l2") + theme_bw() + theme(plot.title = element_text(hjust = 0.5),
       plot.subtitle = element_text(hjust = 0.5))

#########################################  CIBERSORTx plots  ############################################################
install.packages("ggplot2")
install.packages("ggpairs")
install.packages("GGally")
library(GGally)
library(ggpairs)
library(ggplot2)
require(GGally)
cibersortFractions <- read.table("CIBERSORTxGEP_Job9_Fractions.csv", header = T, sep = ',',stringsAsFactors = F)

cibersortFractions1 <- read.table("CIBERSORTxGEP_Job12_Fractions.csv", header = T, sep = ',',stringsAsFactors = F)

ggpairs(data=cibersortFractions, # data.frame with variables
        columns=2:5, # columns to plot, default to all.
        title="Fraction data", # title of the plot
        aes(color = sapply(cibersortFractions$Mixture,FUN = function(x){strsplit(x=x,split=".",fixed=TRUE)[[1]][1]}))) # aesthetics, ggplot2 style

ggpairs(data=cibersortFractions, # data.frame with variables
        columns=6:9, # columns to plot, default to all.
        title="Fraction data", # title of the plot
        aes(color = sapply(cibersortFractions$Mixture,FUN = function(x){strsplit(x=x,split=".",fixed=TRUE)[[1]][1]}))) # aesthetics, ggplot2 style

ggpairs(data=cibersortFractions, # data.frame with variables
        columns=10:13, # columns to plot, default to all.
        title="Fraction data", # title of the plot
        aes(color = sapply(cibersortFractions$Mixture,FUN = function(x){strsplit(x=x,split=".",fixed=TRUE)[[1]][1]}))) # aesthetics, ggplot2 style

ggpairs(data=cibersortFractions, # data.frame with variables
        columns=18:20, # columns to plot, default to all.
        title="Fraction data", # title of the plot
        aes(color = sapply(cibersortFractions$Mixture,FUN = function(x){strsplit(x=x,split=".",fixed=TRUE)[[1]][1]}))) # aesthetics, ggplot2 style

ggpairs(data=cibersortFractions, # data.frame with variables
        columns=21:23, # columns to plot, default to all.
        title="Fraction data", # title of the plot
        aes(color = sapply(cibersortFractions$Mixture,FUN = function(x){strsplit(x=x,split=".",fixed=TRUE)[[1]][1]}))) # aesthetics, ggplot2 style

#sapply(Mixture,FUN = function(x){split(x,".")[1]})

ggpairs(data=cibersortFractions1, # data.frame with variables
        columns=2:5, # columns to plot, default to all.
        title="Fraction data", # title of the plot
        aes(color = sapply(cibersortFractions1$Mixture,FUN = function(x){strsplit(x=x,split=".",fixed=TRUE)[[1]][1]}))) # aesthetics, ggplot2 style

ggpairs(data=cibersortFractions1, # data.frame with variables
        columns=6:9, # columns to plot, default to all.
        title="Fraction data", # title of the plot
        aes(color = sapply(cibersortFractions1$Mixture,FUN = function(x){strsplit(x=x,split=".",fixed=TRUE)[[1]][1]}))) # aesthetics, ggplot2 style

ggpairs(data=cibersortFractions1, # data.frame with variables
        columns=10:13, # columns to plot, default to all.
        title="Fraction data", # title of the plot
        aes(color = sapply(cibersortFractions1$Mixture,FUN = function(x){strsplit(x=x,split=".",fixed=TRUE)[[1]][1]}))) # aesthetics, ggplot2 style

ggpairs(data=cibersortFractions1, # data.frame with variables
        columns=14:17, # columns to plot, default to all.
        title="Fraction data", # title of the plot
        aes(color = sapply(cibersortFractions1$Mixture,FUN = function(x){strsplit(x=x,split=".",fixed=TRUE)[[1]][1]}))) # aesthetics, ggplot2 style

ggpairs(data=cibersortFractions1, # data.frame with variables
        columns=18:20, # columns to plot, default to all.
        title="Fraction data", # title of the plot
        aes(color = sapply(cibersortFractions1$Mixture,FUN = function(x){strsplit(x=x,split=".",fixed=TRUE)[[1]][1]}))) # aesthetics, ggplot2 style

ggpairs(data=cibersortFractions1, # data.frame with variables
        columns=21:23, # columns to plot, default to all.
        title="Fraction data", # title of the plot
        aes(color = sapply(cibersortFractions1$Mixture,FUN = function(x){strsplit(x=x,split=".",fixed=TRUE)[[1]][1]}))) # aesthetics, ggplot2 style



combinedFractions <- read.table("fraction.csv", header = T, sep = ',',stringsAsFactors = F)
ggpairs(data=combinedFractions, # data.frame with variables
        columns=2:4, # columns to plot, default to all.
        title="Fraction data", # title of the plot
        aes(color = sapply(combinedFractions$Mixture,FUN = function(x){strsplit(x=x,split=".",fixed=TRUE)[[1]][1]}))) # aesthetics, ggplot2 style

ggpairs(data=combinedFractions, # data.frame with variables
        columns=5:7, # columns to plot, default to all.
        title="Fraction data", # title of the plot
        aes(color = sapply(combinedFractions$Mixture,FUN = function(x){strsplit(x=x,split=".",fixed=TRUE)[[1]][1]}))) # aesthetics, ggplot2 style

ggpairs(data=combinedFractions, # data.frame with variables
        columns=8:10, # columns to plot, default to all.
        title="Fraction data", # title of the plot
        aes(color = sapply(combinedFractions$Mixture,FUN = function(x){strsplit(x=x,split=".",fixed=TRUE)[[1]][1]}))) # aesthetics, ggplot2 style

ggpairs(data=combinedFractions, # data.frame with variables
        columns=11:13, # columns to plot, default to all.
        title="Fraction data", # title of the plot
        aes(color = sapply(combinedFractions$Mixture,FUN = function(x){strsplit(x=x,split=".",fixed=TRUE)[[1]][1]}))) # aesthetics, ggplot2 style

ggpairs(data=combinedFractions, # data.frame with variables
        columns=14:16, # columns to plot, default to all.
        title="Fraction data", # title of the plot
        aes(color = sapply(combinedFractions$Mixture,FUN = function(x){strsplit(x=x,split=".",fixed=TRUE)[[1]][1]}))) # aesthetics, ggplot2 style

ggpairs(data=combinedFractions, # data.frame with variables
        columns=17:19, # columns to plot, default to all.
        title="Fraction data", # title of the plot
        aes(color = sapply(combinedFractions$Mixture,FUN = function(x){strsplit(x=x,split=".",fixed=TRUE)[[1]][1]}))) # aesthetics, ggplot2 style


ModifiedFractions <- read.table("modified_fractions.csv", header = T, sep = ',',stringsAsFactors = F)
ggpairs(data=ModifiedFractions, # data.frame with variables
        columns=6:6, # columns to plot, default to all.
        title="Fraction data - Daily and Weekly", # title of the plot
        aes(color = sapply(ModifiedFractions$Mixture,FUN = function(x){strsplit(x=x,split=".",fixed=TRUE)[[1]][1]}))) # aesthetics, ggplot2 style

ggpairs(data=ModifiedFractions, # data.frame with variables
        columns=5:5, # columns to plot, default to all.
        title="Fraction data - Daily and Weekly", # title of the plot
        aes(color = sapply(ModifiedFractions$Mixture,FUN = function(x){strsplit(x=x,split=".",fixed=TRUE)[[1]][1]}))) # aesthetics, ggplot2 style

ggpairs(data=ModifiedFractions, # data.frame with variables
        columns=2:2, # columns to plot, default to all.
        title="Fraction data - Daily and Weekly", # title of the plot
        aes(color = sapply(ModifiedFractions$Mixture,FUN = function(x){strsplit(x=x,split=".",fixed=TRUE)[[1]][1]}))) # aesthetics, ggplot2 style

ggpairs(data=ModifiedFractions, # data.frame with variables
        columns=4:4, # columns to plot, default to all.
        title="Fraction data - Daily and Weekly", # title of the plot
        aes(color = sapply(ModifiedFractions$Mixture,FUN = function(x){strsplit(x=x,split=".",fixed=TRUE)[[1]][1]}))) # aesthetics, ggplot2 style

ggpairs(data=ModifiedFractions, # data.frame with variables
        columns=9:9, # columns to plot, default to all.
        title="Fraction data - Daily and Weekly", # title of the plot
        aes(color = sapply(ModifiedFractions$Mixture,FUN = function(x){strsplit(x=x,split=".",fixed=TRUE)[[1]][1]}))) # aesthetics, ggplot2 style

ggpairs(data=ModifiedFractions, # data.frame with variables
        columns=16:16, # columns to plot, default to all.
        title="Fraction data - Daily and Weekly", # title of the plot
        aes(color = sapply(ModifiedFractions$Mixture,FUN = function(x){strsplit(x=x,split=".",fixed=TRUE)[[1]][1]}))) # aesthetics, ggplot2 style

ggpairs(data=ModifiedFractions, # data.frame with variables
        columns=11:11, # columns to plot, default to all.
        title="Fraction data - Daily and Weekly", # title of the plot
        aes(color = sapply(ModifiedFractions$Mixture,FUN = function(x){strsplit(x=x,split=".",fixed=TRUE)[[1]][1]}))) # aesthetics, ggplot2 style


######################################################################################################
ggpairs(cibersortFractions, mapping = NULL, columns = 1:ncol(cibersortFractions), title = NULL,
        upper = list(continuous = "cor", combo = "box_no_facet", discrete ="facetbar", na = "na"), 
        lower = list(continuous = "points", combo ="facethist", discrete = "facetbar", na = "na"), 
        diag = list(continuous ="densityDiag", discrete = "barDiag", na = "naDiag"), params = NULL,
        xlab = NULL, ylab = NULL, axisLabels = c("show", "internal", "none"),
        columnLabels = colnames(cibersortFractions[columns]), 
        labeller = "label_value",switch = NULL, showStrips = NULL, legend = NULL,cardinality_threshold = 15, 
        progress = NULL)

#############################  Differential Expression Analysis for RNAseq data  ###########################################

row.names(data) <- data$gene.y

data$gene.y <- NULL

data2 <- as.data.frame(t(data))

data.daily <- data2[-c(6, 7, 8, 9,10,11,12,13,14,15), ]
data.weekly <- data2[-c(1, 2, 3, 4,5,11,12,13,14,15), ]
data.vehicle <- data2[-c(1, 2, 3, 4,5,6,7,8,9,10), ]
data.daily_weekly <- data2[-c(11,12,13,14,15), ]

#write.table(x = data.daily, file = "C:/My_Work/R/Projects/data.dailyGS.txt", sep = "\t", row.names = F, col.names = T)

#write.table(x = data.weekly, file = "C:/My_Work/R/Projects/data.weeklyGS.txt", sep = "\t", row.names = F, col.names = T)

#data.deltatype <- data.frame(Delta$cell_type1)
data.daily$gene.y <- NULL
data.weekly$gene.y <- NULL
data.vehicle$gene.y <- NULL
data.daily_weekly$gene.y <- NULL

data.daily["Mean.daily",] <- colMeans(data.daily)
data.weekly["Mean.weekly",] <- colMeans(data.weekly)
data.vehicle["Mean.vehicle",] <- colMeans(data.vehicle)
data.daily_weekly["Mean.daily_weekly",] <- colMeans(data.daily_weekly)

data.daily.2 <- as.data.frame(t(data.daily))
data.weekly.2 <- as.data.frame(t(data.weekly))
data.vehicle.2 <- as.data.frame(t(data.vehicle))
data.daily_weekly.2 <- as.data.frame(t(data.daily_weekly))
######################################################################################################################
##DE analysis for vehicle vs daily

vehicle.V.daily <- data.frame(gene.y = rownames(data.vehicle.2),
                              Mean.vehicle = data.vehicle.2$Mean.vehicle,
                              Mean.daily = data.daily.2$Mean.daily)


vehicle.V.daily$log2FC <- with(vehicle.V.daily, log2(Mean.vehicle/Mean.daily))

data.DF1 <- data2
data.DF1$gene.y <- NULL
data.DF1 <- as.data.frame(t(data.DF1))


## Describe a function to calculate the p-value using student's t-test using the data.DF4 data frame.
calc.pval <- function(x){
  C <- x[1:5]
  D <- x[11:15]
  p <- t.test(C,D, alternative = "two.sided", paired = F, var.equal = FALSE)$p.value
  
  return(p)
}

## Applying the p-value calculation function across each row of the DF3 dataframe
vehicle.V.daily$pval <- apply(X = data.DF1, MARGIN = 1, FUN = calc.pval)

## calculating negative log10 value of the pvalue.
vehicle.V.daily$log10.pval <- -log10(vehicle.V.daily$pval)

## Multiple Testing Correction using "Benjamini & Hochberg" method.
vehicle.V.daily$p.adj <- p.adjust(vehicle.V.daily$pval, method = "fdr")
vehicle.V.daily$log10.padj <- -log10(vehicle.V.daily$p.adj)

## Volcano Plot using ggplot

# Set the thresholds for p-value and fold change
p.cutoff <- -log10(0.05)     # p-value cutoff of p = 0.05
fc.cutoff <- log2(1.15)      # Log2FC cutoff of 15%

# With corrected p-values.
# Initiating the ggplot with log2FC on X-axis and log10.pval on Y-axis.
ggplot(vehicle.V.daily, aes(log2FC, log10.pval, label = gene.y)) +    
  geom_point(col = "black") +   # Baseline points as black
  xlim(c(-3,3))+           # Set the X-axis limits such the plot is symmetrical
  ylim(c(0,8))+
  geom_point(data = subset(vehicle.V.daily, log10.pval > p.cutoff), col = "darkgray") +   
  geom_point(data = subset(vehicle.V.daily, log10.pval > p.cutoff & log2FC < -fc.cutoff), col = "green") +  
  geom_point(data = subset(vehicle.V.daily, log10.pval > p.cutoff & log2FC > fc.cutoff), col = "red") + 
  geom_hline(yintercept = p.cutoff, col = "black") +   # Horizontal Line for p-value threshold
  geom_vline(xintercept = -fc.cutoff, col = "green") +  # Vertical Line for negative Fold Change Threshold
  geom_vline(xintercept = fc.cutoff, col = "red") +     # Vertical Line for positive Fold Change Threshold
  #geom_text_repel(data = subset(data.plot, log10.pval > p.cutoff & abs(log2FC) > fc.cutoff)) +    # Add labels to significant probes. Only use this if the number of significant probes are small
  xlab("log2(Vehicle/Daily)") + ylab("-log10(p-value)") + # Adding appropriate labels for X and Y axes.
  ggtitle(label = "Differential Expression Analysis of Vehicle/Daily", 
          subtitle = "With FDR corrected p-values") +     # Plot Title
  theme_bw() +  # ggplot usually puts a grey background. This theme puts a simple white background with nice box outlines
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) 

vehicle_daily_up <- subset(vehicle.V.daily [order(vehicle.V.daily$pval),], log10.pval > p.cutoff 
                       & log2FC > fc.cutoff)
write.table(x = vehicle_daily_up, file = "C:/My_Work/R/Projects/vehicle_daily_up.txt", sep = "\t",
            row.names = F, col.names = T)

vehicle_daily_down <- subset(vehicle.V.daily [order(vehicle.V.daily$pval),], log10.pval > p.cutoff 
                           & log2FC < -fc.cutoff)
write.table(x = vehicle_daily_down, file = "C:/My_Work/R/Projects/vehicle_daily_down.txt", sep = "\t",
            row.names = F, col.names = T)
######################################################################################################################
##DE analysis for vehicle vs weekly

vehicle.V.weekly <- data.frame(gene.y = rownames(data.vehicle.2),
                              Mean.vehicle = data.vehicle.2$Mean.vehicle,
                              Mean.weekly = data.weekly.2$Mean.weekly)


vehicle.V.weekly$log2FC <- with(vehicle.V.weekly, log2(Mean.vehicle/Mean.weekly))

data.DF1 <- data2
data.DF1$gene.y <- NULL
data.DF1 <- as.data.frame(t(data.DF1))


## Describe a function to calculate the p-value using student's t-test using the data.DF4 data frame.
calc.pval <- function(x){
  C <- x[6:10]
  D <- x[11:15]
  p <- t.test(C,D, alternative = "two.sided", paired = F, var.equal = FALSE)$p.value
  
  return(p)
}

## Applying the p-value calculation function across each row of the DF3 dataframe
vehicle.V.weekly$pval <- apply(X = data.DF1, MARGIN = 1, FUN = calc.pval)

## calculating negative log10 value of the pvalue.
vehicle.V.weekly$log10.pval <- -log10(vehicle.V.weekly$pval)

## Multiple Testing Correction using "Benjamini & Hochberg" method.
vehicle.V.weekly$p.adj <- p.adjust(vehicle.V.weekly$pval, method = "fdr")
vehicle.V.weekly$log10.padj <- -log10(vehicle.V.weekly$p.adj)

## Volcano Plot using ggplot

# Set the thresholds for p-value and fold change
p.cutoff <- -log10(0.05)     # p-value cutoff of p = 0.05
fc.cutoff <- log2(1.15)      # Log2FC cutoff of 15%

# With corrected p-values.
# Initiating the ggplot with log2FC on X-axis and log10.pval on Y-axis.
ggplot(vehicle.V.weekly, aes(log2FC, log10.pval, label = gene.y)) +    
  geom_point(col = "black") +   # Baseline points as black
  xlim(c(-3,3))+           # Set the X-axis limits such the plot is symmetrical
  ylim(c(0,8))+
  geom_point(data = subset(vehicle.V.weekly, log10.pval > p.cutoff), col = "darkgray") +   
  geom_point(data = subset(vehicle.V.weekly, log10.pval > p.cutoff & log2FC < -fc.cutoff), col = "green") +  
  geom_point(data = subset(vehicle.V.weekly, log10.pval > p.cutoff & log2FC > fc.cutoff), col = "red") + 
  geom_hline(yintercept = p.cutoff, col = "black") +   # Horizontal Line for p-value threshold
  geom_vline(xintercept = -fc.cutoff, col = "green") +  # Vertical Line for negative Fold Change Threshold
  geom_vline(xintercept = fc.cutoff, col = "red") +     # Vertical Line for positive Fold Change Threshold
  #geom_text_repel(data = subset(data.plot, log10.pval > p.cutoff & abs(log2FC) > fc.cutoff)) +    # Add labels to significant probes. Only use this if the number of significant probes are small
  xlab("log2(Vehicle/weekly)") + ylab("-log10(p-value)") + # Adding appropriate labels for X and Y axes.
  ggtitle(label = "Differential Expression Analysis of Vehicle/weekly", 
          subtitle = "With FDR corrected p-values") +     # Plot Title
  theme_bw() +  # ggplot usually puts a grey background. This theme puts a simple white background with nice box outlines
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) 

vehicle_weekly_up <- subset(vehicle.V.weekly [order(vehicle.V.weekly$pval),], log10.pval > p.cutoff 
                           & log2FC > fc.cutoff)
write.table(x = vehicle_weekly_up, file = "C:/My_Work/R/Projects/vehicle_weekly_up.txt", sep = "\t",
            row.names = F, col.names = T)

vehicle_weekly_down <- subset(vehicle.V.weekly [order(vehicle.V.weekly$pval),], log10.pval > p.cutoff 
                             & log2FC < -fc.cutoff)
write.table(x = vehicle_weekly_down, file = "C:/My_Work/R/Projects/vehicle_weekly_down.txt", sep = "\t",
            row.names = F, col.names = T)
########################################################################################################################

########################################################################################################################
##DE analysis for vehicle vs Daily+weekly

vehicle.V.daily_weekly <- data.frame(gene.y = rownames(data.vehicle.2),
                               Mean.vehicle = data.vehicle.2$Mean.vehicle,
                               Mean.daily_weekly = data.daily_weekly.2$Mean.daily_weekly)


vehicle.V.daily_weekly$log2FC <- with(vehicle.V.daily_weekly, log2(Mean.vehicle/Mean.daily_weekly))

data.DF1 <- data2
data.DF1$gene.y <- NULL
data.DF1 <- as.data.frame(t(data.DF1))


## Describe a function to calculate the p-value using student's t-test using the data.DF4 data frame.
calc.pval <- function(x){
  C <- x[1:10]
  D <- x[11:15]
  p <- t.test(D, C, alternative = "two.sided", paired = F, var.equal = FALSE)$p.value
  
  return(p)
}

## Applying the p-value calculation function across each row of the DF3 dataframe
vehicle.V.daily_weekly$pval <- apply(X = data.DF1, MARGIN = 1, FUN = calc.pval)

## calculating negative log10 value of the pvalue.
vehicle.V.daily_weekly$log10.pval <- -log10(vehicle.V.daily_weekly$pval)

## Multiple Testing Correction using "Benjamini & Hochberg" method.
vehicle.V.daily_weekly$p.adj <- p.adjust(vehicle.V.daily_weekly$pval, method = "fdr")
vehicle.V.daily_weekly$log10.padj <- -log10(vehicle.V.daily_weekly$p.adj)

## Volcano Plot using ggplot

# Set the thresholds for p-value and fold change
p.cutoff <- -log10(0.05)     # p-value cutoff of p = 0.05
fc.cutoff <- log2(1.15)      # Log2FC cutoff of 15%

# With corrected p-values.
# Initiating the ggplot with log2FC on X-axis and log10.pval on Y-axis.
ggplot(vehicle.V.daily_weekly, aes(log2FC, log10.pval, label = gene.y)) +    
  geom_point(col = "black") +   # Baseline points as black
  xlim(c(-3,3))+           # Set the X-axis limits such the plot is symmetrical
  ylim(c(0,8))+
  geom_point(data = subset(vehicle.V.daily_weekly, log10.pval > p.cutoff), col = "darkgray") +   
  geom_point(data = subset(vehicle.V.daily_weekly, log10.pval > p.cutoff & log2FC < -fc.cutoff), col = "green") +  
  geom_point(data = subset(vehicle.V.daily_weekly, log10.pval > p.cutoff & log2FC > fc.cutoff), col = "red") + 
  geom_hline(yintercept = p.cutoff, col = "black") +   # Horizontal Line for p-value threshold
  geom_vline(xintercept = -fc.cutoff, col = "green") +  # Vertical Line for negative Fold Change Threshold
  geom_vline(xintercept = fc.cutoff, col = "red") +     # Vertical Line for positive Fold Change Threshold
  #geom_text_repel(data = subset(data.plot, log10.pval > p.cutoff & abs(log2FC) > fc.cutoff)) +    # Add labels to significant probes. Only use this if the number of significant probes are small
  xlab("log2(Vehicle/daily+weekly)") + ylab("-log10(p-value)") + # Adding appropriate labels for X and Y axes.
  ggtitle(label = "Differential Expression Analysis of Vehicle/daily+weekly", 
          subtitle = "With FDR corrected p-values") +     # Plot Title
  theme_bw() +  # ggplot usually puts a grey background. This theme puts a simple white background with nice box outlines
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) 

vehicle_daily_weekly_up <- subset(vehicle.V.daily_weekly [order(vehicle.V.daily_weekly$pval),], 
                                  log10.pval > p.cutoff & log2FC > fc.cutoff)
write.table(x = vehicle_daily_weekly_up, file = "C:/My_Work/R/Projects/vehicle_daily_weekly_up.txt", sep = "\t",
            row.names = F, col.names = T)

vehicle_daily_weekly_down <- subset(vehicle.V.daily_weekly [order(vehicle.V.daily_weekly$pval),], 
                                    log10.pval > p.cutoff & log2FC < -fc.cutoff)
write.table(x = vehicle_daily_weekly_down, file = "C:/My_Work/R/Projects/vehicle_daily_weekly_down.txt", sep = "\t",
            row.names = F, col.names = T)
########################################################################################################################

########################################################################################################################
Daily.vs.Weekly <- data.frame(gene.y = rownames(data.daily.2),
                           Mean.daily = data.daily.2$Mean.daily,
                           Mean.weekly = data.weekly.2$Mean.weekly)


Daily.vs.Weekly$log2FC <- with(Daily.vs.Weekly, log2(Mean.daily/Mean.weekly))

data.DF1 <- data2
data.DF1$gene.y <- NULL
data.DF1 <- as.data.frame(t(data.DF1))


## Describe a function to calculate the p-value using student's t-test using the data.DF4 data frame.
calc.pval <- function(x){
  C <- x[1:5]
  D <- x[6:10]
  p <- t.test(C,D, alternative = "two.sided", paired = F, var.equal = FALSE)$p.value
  
  return(p)
}

## Applying the p-value calculation function across each row of the DF3 dataframe
Daily.vs.Weekly$pval <- apply(X = data.DF1, MARGIN = 1, FUN = calc.pval)

## calculating negative log10 value of the pvalue.
Daily.vs.Weekly$log10.pval <- -log10(Daily.vs.Weekly$pval)

## Multiple Testing Correction using "Benjamini & Hochberg" method.
Daily.vs.Weekly$p.adj <- p.adjust(Daily.vs.Weekly$pval, method = "fdr")
Daily.vs.Weekly$log10.padj <- -log10(Daily.vs.Weekly$p.adj)

## Volcano Plot using ggplot

# Set the thresholds for p-value and fold change
p.cutoff <- -log10(0.05)     # p-value cutoff of p = 0.05
fc.cutoff <- log2(1.15)      # Log2FC cutoff of 15%

# With corrected p-values.
# Initiating the ggplot with log2FC on X-axis and log10.pval on Y-axis.
ggplot(Daily.vs.Weekly, aes(log2FC, log10.pval, label = gene.y)) +    
  geom_point(col = "black") +   # Baseline points as black
  xlim(c(-3,2))+           # Set the X-axis limits such the plot is symmetrical
  ylim(c(0,5))+
  geom_point(data = subset(Daily.vs.Weekly, log10.pval > p.cutoff), col = "darkgray") +   
  geom_point(data = subset(Daily.vs.Weekly, log10.pval > p.cutoff & log2FC < -fc.cutoff), col = "green") +  
  geom_point(data = subset(Daily.vs.Weekly, log10.pval > p.cutoff & log2FC > fc.cutoff), col = "red") + 
  geom_hline(yintercept = p.cutoff, col = "black") +   # Horizontal Line for p-value threshold
  geom_vline(xintercept = -fc.cutoff, col = "green") +  # Vertical Line for negative Fold Change Threshold
  geom_vline(xintercept = fc.cutoff, col = "red") +     # Vertical Line for positive Fold Change Threshold
  #geom_text_repel(data = subset(data.plot, log10.pval > p.cutoff & abs(log2FC) > fc.cutoff)) +    # Add labels to significant probes. Only use this if the number of significant probes are small
  xlab("log2(Daily/Weekly)") + ylab("-log10(p-value)") + # Adding appropriate labels for X and Y axes.
  ggtitle(label = "Differential Expression Analysis of Daily/Weekly", 
          subtitle = "With FDR corrected p-values") +     # Plot Title
  theme_bw() +  # ggplot usually puts a grey background. This theme puts a simple white background with nice box outlines
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) 

daily_weekly_up <- subset(Daily.vs.Weekly [order(Daily.vs.Weekly$pval),], 
                                  log10.pval > p.cutoff & log2FC > fc.cutoff)
write.table(x = daily_weekly_up, file = "C:/My_Work/R/Projects/daily_weekly_up.txt", sep = "\t",
            row.names = F, col.names = T)

daily_weekly_down <- subset(Daily.vs.Weekly [order(Daily.vs.Weekly$pval),], 
                                    log10.pval > p.cutoff & log2FC < -fc.cutoff)
write.table(x = daily_weekly_down, file = "C:/My_Work/R/Projects/daily_weekly_down.txt", sep = "\t",
            row.names = F, col.names = T)