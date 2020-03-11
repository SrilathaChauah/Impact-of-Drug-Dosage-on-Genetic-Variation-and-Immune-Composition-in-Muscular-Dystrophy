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

## load Gene microarray data
gset <- getGEO("GSE84876", GSEMatrix =TRUE, getGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL20258", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

## Extract the expression data
edata <- exprs(gset)

## Extract the phenotype information
pdata <- pData(gset)

## Extract the features (probes) information
fdata <- fData(gset)

## Transposing the expression data such that the rows are the samples and columns are the probes
edata2 <- t(edata)

## Adding the Sample type information to be able to group all the control and DM2 samples
data.DF3 <- data.frame(Sample.Type = pdata$source_name_ch1, edata2)
data.DF3.2 <- as.data.frame(t(data.DF3))
write.table(x = data.DF3.2, file = "C:/My_Work/R/Projects/data.DF3.2.txt", sep = "\t", row.names = F, col.names = T)

## Separating the control and DM2 samples into separate data frames.
data.untreated <- subset(data.DF3, Sample.Type %in% "quadriceps, utrn+/-;mdx, untreated")
data.spir_lisi <- subset(data.DF3, Sample.Type %in% "quadriceps, utrn+/-;mdx, spironolactone plus lisinopril-treated")
data.eple_lisi <- subset(data.DF3, Sample.Type %in% "quadriceps, utrn+/-;mdx, eplerenone plus lisinopril-treated")
data.prednisolone <- subset(data.DF3, Sample.Type %in% "quadriceps, utrn+/-;mdx, prednisolone-treated")

library(reshape2)
library(meltt)
DATA <- read.table("melteddata.csv", header = T, stringsAsFactors = F, sep = ',')
######################################################################################################################
##Density Plots

melted_data1 = melt(DATA, id.vars = "gene.y")

melted_data1$group <- "NA"
melted_data1$group[melted_data1$variable %in% c("untreated", 
                                              "untreated1", 
                                              "untreated2")] <- "untreated"

melted_data1$group[melted_data1$variable %in% c("prednisolone.treated", 
                                              "prednisolone.treated1", 
                                              "prednisolone.treated2")] <- "treated"

library(ggplot2)

ggplot(subset(melted_data1, melted_data1$gene.y %in% "Dbp"),
       aes(x=value,fill=group)) + geom_density() + ggtitle(label = "Melted plot for Chadwick's dataset", 
       subtitle = "Gene Symbol: Dbp") + theme_bw() + theme(plot.title = element_text(hjust = 0.5),
       plot.subtitle = element_text(hjust = 0.5))

ggplot(subset(melted_data1, melted_data1$gene.y %in% "U6"),
       aes(x=value,fill=group)) + geom_density() + ggtitle(label = "Melted plot for Chadwick's dataset", 
       subtitle = "Gene Symbol: U6") + theme_bw() + theme(plot.title = element_text(hjust = 0.5),
       plot.subtitle = element_text(hjust = 0.5))


ggplot(subset(melted_data1, melted_data1$gene.y %in% "Abcc8"),
       aes(x=value,fill=group)) + geom_density() + ggtitle(label = "Melted plot for Chadwick's dataset", 
       subtitle = "Gene Symbol: Abcc8") + theme_bw() + theme(plot.title = element_text(hjust = 0.5),
       plot.subtitle = element_text(hjust = 0.5))

ggplot(subset(melted_data1, melted_data1$gene.y %in% "Arntl"),
       aes(x=value,fill=group)) + geom_density() + ggtitle(label = "Melted plot for Chadwick's dataset", 
       subtitle = "Gene Symbol: Arntl") + theme_bw() + theme(plot.title = element_text(hjust = 0.5),
       plot.subtitle = element_text(hjust = 0.5))

ggplot(subset(melted_data1, melted_data1$gene.y %in% "Car14"),
       aes(x=value,fill=group)) + geom_density() + ggtitle(label = "Melted plot for Chadwick's dataset", 
       subtitle = "Gene Symbol: Car14") + theme_bw() + theme(plot.title = element_text(hjust = 0.5),
       plot.subtitle = element_text(hjust = 0.5))

ggplot(subset(melted_data1, melted_data1$gene.y %in% "C7"),
       aes(x=value,fill=group)) + geom_density() + ggtitle(label = "Melted plot for Chadwick's dataset", 
       subtitle = "Gene Symbol: C7") + theme_bw() + theme(plot.title = element_text(hjust = 0.5),
       plot.subtitle = element_text(hjust = 0.5))

ggplot(subset(melted_data1, melted_data1$gene.y %in% "Agt"),
       aes(x=value,fill=group)) + geom_density() + ggtitle(label = "Melted plot for Chadwick's dataset", 
       subtitle = "Gene Symbol: Agt") + theme_bw() + theme(plot.title = element_text(hjust = 0.5),
       plot.subtitle = element_text(hjust = 0.5))

ggplot(subset(melted_data1, melted_data1$gene.y %in% "Npas2"),
       aes(x=value,fill=group)) + geom_density() + ggtitle(label = "Melted plot for Chadwick's dataset", 
       subtitle = "Gene Symbol: Npas2") + theme_bw() + theme(plot.title = element_text(hjust = 0.5),
       plot.subtitle = element_text(hjust = 0.5))

ggplot(subset(melted_data1, melted_data1$gene.y %in% "Odf3l2"),
       aes(x=value,fill=group)) + geom_density() + ggtitle(label = "Melted plot for Chadwick's dataset", 
       subtitle = "Gene Symbol: Odf3l2") + theme_bw() + theme(plot.title = element_text(hjust = 0.5),
       plot.subtitle = element_text(hjust = 0.5))
#========================================================================================
ggplot(subset(melted_data1, melted_data1$gene.y %in% "A730049H05Rik"),
       aes(x=value,fill=group)) + geom_density() + ggtitle(label = "Melted plot for Chadwick's dataset", 
                                                           subtitle = "Gene Symbol: A730049H05Rik") + theme_bw() + theme(plot.title = element_text(hjust = 0.5),
                                                                                                                         plot.subtitle = element_text(hjust = 0.5)) 

ggplot(subset(melted_data1, melted_data1$gene.y %in% "Gm12840"),
       aes(x=value,fill=group)) + geom_density() + ggtitle(label = "Melted plot for Chadwick's dataset", 
                                                           subtitle = "Gene Symbol: Gm12840") + theme_bw() + theme(plot.title = element_text(hjust = 0.5),
                                                                                                                   plot.subtitle = element_text(hjust = 0.5)) 

ggplot(subset(melted_data1, melted_data1$gene.y %in% "Fosl2"),
       aes(x=value,fill=group)) + geom_density() + ggtitle(label = "Melted plot for Chadwick's dataset", 
                                                           subtitle = "Gene Symbol: Fosl2") + theme_bw() + theme(plot.title = element_text(hjust = 0.5),
                                                                                                                 plot.subtitle = element_text(hjust = 0.5)) 

ggplot(subset(melted_data1, melted_data1$gene.y %in% "Zfp36"),
       aes(x=value,fill=group)) + geom_density() + ggtitle(label = "Melted plot for Chadwick's dataset", 
                                                           subtitle = "Gene Symbol: Zfp36") + theme_bw() + theme(plot.title = element_text(hjust = 0.5),
                                                                                                                 plot.subtitle = element_text(hjust = 0.5)) 

ggplot(subset(melted_data1, melted_data1$gene.y %in% "Csrnp1"),
       aes(x=value,fill=group)) + geom_density() + ggtitle(label = "Melted plot for Chadwick's dataset", 
                                                           subtitle = "Gene Symbol: Csrnp1") + theme_bw() + theme(plot.title = element_text(hjust = 0.5),
                                                                                                                  plot.subtitle = element_text(hjust = 0.5)) 

ggplot(subset(melted_data1, melted_data1$gene.y %in% "Pde4b"),
       aes(x=value,fill=group)) + geom_density() + ggtitle(label = "Melted plot for Chadwick's dataset", 
                                                           subtitle = "Gene Symbol: Pde4b") + theme_bw() + theme(plot.title = element_text(hjust = 0.5),
                                                                                                                 plot.subtitle = element_text(hjust = 0.5)) 
ggplot(subset(melted_data1, melted_data1$gene.y %in% "Tmem252"),
       aes(x=value,fill=group)) + geom_density() + ggtitle(label = "Melted plot for Chadwick's dataset", 
                                                           subtitle = "Gene Symbol: Tmem252") + theme_bw() + theme(plot.title = element_text(hjust = 0.5),
                                                                                                                   plot.subtitle = element_text(hjust = 0.5)) 

ggplot(subset(melted_data1, melted_data1$gene.y %in% "Rapgef5"),
       aes(x=value,fill=group)) + geom_density() + ggtitle(label = "Melted plot for Chadwick's dataset", 
                                                           subtitle = "Gene Symbol: Rapgef5") + theme_bw() + theme(plot.title = element_text(hjust = 0.5),
                                                                                                                   plot.subtitle = element_text(hjust = 0.5)) 

ggplot(subset(melted_data1, melted_data1$gene.y %in% "Mnda"),
       aes(x=value,fill=group)) + geom_density() + ggtitle(label = "Melted plot for Chadwick's dataset", 
                                                           subtitle = "Gene Symbol: Mnda") + theme_bw() + theme(plot.title = element_text(hjust = 0.5),
                                                                                                                plot.subtitle = element_text(hjust = 0.5)) 

ggplot(subset(melted_data1, melted_data1$gene.y %in% "Ret"),
       aes(x=value,fill=group)) + geom_density() + ggtitle(label = "Melted plot for Chadwick's dataset", 
                                                           subtitle = "Gene Symbol: Ret") + theme_bw() + theme(plot.title = element_text(hjust = 0.5),
                                                                                                               plot.subtitle = element_text(hjust = 0.5)) 

ggplot(subset(melted_data1, melted_data1$gene.y %in% "Dusp5"),
       aes(x=value,fill=group)) + geom_density() + ggtitle(label = "Melted plot for Chadwick's dataset", 
                                                           subtitle = "Gene Symbol: Dusp5") + theme_bw() + theme(plot.title = element_text(hjust = 0.5),
                                                                                                                 plot.subtitle = element_text(hjust = 0.5)) 

ggplot(subset(melted_data1, melted_data1$gene.y %in% "Tmem2"),
       aes(x=value,fill=group)) + geom_density() + ggtitle(label = "Melted plot for Chadwick's dataset", 
                                                           subtitle = "Gene Symbol: Tmem2") + theme_bw() + theme(plot.title = element_text(hjust = 0.5),
                                                                                                                 plot.subtitle = element_text(hjust = 0.5)) 

ggplot(subset(melted_data1, melted_data1$gene.y %in% "RP24-196H4.1"),
       aes(x=value,fill=group)) + geom_density() + ggtitle(label = "Melted plot for Chadwick's dataset", 
                                                           subtitle = "Gene Symbol: RP24-196H4.1") + theme_bw() + theme(plot.title = element_text(hjust = 0.5),
                                                                                                                        plot.subtitle = element_text(hjust = 0.5)) 

ggplot(subset(melted_data1, melted_data1$gene.y %in% "Stat2"),
       aes(x=value,fill=group)) + geom_density() + ggtitle(label = "Melted plot for Chadwick's dataset", 
                                                           subtitle = "Gene Symbol: Stat2") + theme_bw() + theme(plot.title = element_text(hjust = 0.5),
                                                                                                                 plot.subtitle = element_text(hjust = 0.5)) 


#===========================================================================================
## Differential Expression Analysis for Microarray data

## Removing the sample type information so as to convert the entire data frame to numerical values
data.untreated$Sample.Type <- NULL
data.prednisolone$Sample.Type <- NULL

## Calculating the mean expression values for each gene across the control and DM2 samples
data.untreated["Mean.untreated",] <- colMeans(data.untreated)
data.prednisolone["Mean.prednisolone",] <- colMeans(data.prednisolone)

## Transposing the previous data frames so that it is easier to access the mean values as a column
data.untreated.2 <- as.data.frame(t(data.untreated))
data.prednisolone.2 <- as.data.frame(t(data.prednisolone))

#### DE Analysis for prednisolone treated vs untreated

## Creating a data frame for plotting which would include the individual mean values, log2 fold ratios and the p-values. 
data.plot <- data.frame(Probe.ID = rownames(data.prednisolone.2),
                        Mean.prednisolone = data.prednisolone.2$Mean.prednisolone,
                        Mean.untreated = data.untreated.2$Mean.untreated)

## Log2FC calculation
data.plot$log2FC <- with(data.plot, log2(Mean.prednisolone/Mean.untreated))

####
myData <- data.DF3[-c(4, 5, 6, 7, 8, 9), ]

## Transpose the data.DF3 after removing the Sample type information to remove the non-numerical values,
## so that the probes are across the rows

data.DF4 <- myData
data.DF4$Sample.Type <- NULL
data.DF4 <- as.data.frame(t(data.DF4))


## Describe a function to calculate the p-value using student's t-test using the data.DF4 data frame.
calc.pval <- function(x){
  C <- x[1:3]
  D <- x[4:6]
  p <- t.test(C,D, alternative = "two.sided", paired = F, var.equal = FALSE)$p.value
  
  return(p)
}

## Applying the p-value calculation function across each row of the DF3 dataframe
data.plot$pval <- apply(X = data.DF4, MARGIN = 1, FUN = calc.pval)

## calculating negative log10 value of the pvalue.
data.plot$log10.pval <- -log10(data.plot$pval)

## Multiple Testing Correction using "Benjamini & Hochberg" method.
data.plot$p.adj <- p.adjust(data.plot$pval, method = "fdr")
data.plot$log10.padj <- -log10(data.plot$p.adj)

## Volcano Plot using ggplot

# Set the thresholds for p-value and fold change
p.cutoff <- -log10(0.05)     # p-value cutoff of p = 0.05
fc.cutoff <- log2(1.15)      # Log2FC cutoff of 15%

# With corrected p-values.
# Initiating the ggplot with log2FC on X-axis and log10.pval on Y-axis.
ggplot(data.plot, aes(log2FC, log10.pval, label = Probe.ID)) +    
  geom_point(col = "black") +   # Baseline points as black
  xlim(c(-2,2))+           # Set the X-axis limits such the plot is symmetrical
  ylim(c(0,6))+
  geom_point(data = subset(data.plot, log10.pval > p.cutoff), col = "darkgray") +   
  geom_point(data = subset(data.plot, log10.pval > p.cutoff & log2FC < -fc.cutoff), col = "green") +  
  geom_point(data = subset(data.plot, log10.pval > p.cutoff & log2FC > fc.cutoff), col = "red") + 
  geom_hline(yintercept = p.cutoff, col = "black") +   # Horizontal Line for p-value threshold
  geom_vline(xintercept = -fc.cutoff, col = "green") +  # Vertical Line for negative Fold Change Threshold
  geom_vline(xintercept = fc.cutoff, col = "red") +     # Vertical Line for positive Fold Change Threshold
  #geom_text_repel(data = subset(data.plot, log10.pval > p.cutoff & abs(log2FC) > fc.cutoff)) +    # Add labels to significant probes. Only use this if the number of significant probes are small
  xlab("log2(Prednisolone/Untreated)") + ylab("-log10(p-value)") + # Adding appropriate labels for X and Y axes.
  ggtitle(label = "Differential Expression Analysis of Prednisolone/Untreated", 
          subtitle = "With FDR corrected p-values") +     # Plot Title
  theme_bw() +  # ggplot usually puts a grey background. This theme puts a simple white background with nice box outlines
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) 

#### DE Analysis for prednisolone untreated vs treated

## Creating a data frame for plotting which would include the individual mean values, log2 fold ratios and the p-values. 
data.plot1 <- data.frame(Probe.ID = rownames(data.untreated.2),
                        Mean.untreated = data.untreated.2$Mean.untreated,
                        Mean.prednisolone = data.prednisolone.2$Mean.prednisolone)

## Log2FC calculation
data.plot1$log2FC <- with(data.plot1, log2(Mean.untreated/Mean.prednisolone))

####
myData1 <- data.DF3[-c(4, 5, 6, 7, 8, 9), ]

## Transpose the data.DF3 after removing the Sample type information to remove the non-numerical values,
## so that the probes are across the rows

data.DF4 <- myData1
data.DF4$Sample.Type <- NULL
data.DF4 <- as.data.frame(t(data.DF4))


## Describe a function to calculate the p-value using student's t-test using the data.DF4 data frame.
calc.pval <- function(x){
  C <- x[1:3]
  D <- x[4:6]
  p <- t.test(C,D, alternative = "two.sided", paired = F, var.equal = FALSE)$p.value
  
  return(p)
}

## Applying the p-value calculation function across each row of the DF3 dataframe
data.plot1$pval <- apply(X = data.DF4, MARGIN = 1, FUN = calc.pval)

## calculating negative log10 value of the pvalue.
data.plot1$log10.pval <- -log10(data.plot1$pval)

## Multiple Testing Correction using "Benjamini & Hochberg" method.
data.plot1$p.adj <- p.adjust(data.plot1$pval, method = "fdr")
data.plot1$log10.padj <- -log10(data.plot1$p.adj)

## Volcano Plot using ggplot

# Set the thresholds for p-value and fold change
p.cutoff <- -log10(0.05)     # p-value cutoff of p = 0.05
fc.cutoff <- log2(1.15)      # Log2FC cutoff of 15%

# With corrected p-values.
# Initiating the ggplot with log2FC on X-axis and log10.pval on Y-axis.
ggplot(data.plot1, aes(log2FC, log10.pval, label = Probe.ID)) +    
  geom_point(col = "black") +   # Baseline points as black
  xlim(c(-2,2))+           # Set the X-axis limits such the plot is symmetrical
  ylim(c(0,6))+
  geom_point(data = subset(data.plot1, log10.pval > p.cutoff), col = "darkgray") +   
  geom_point(data = subset(data.plot1, log10.pval > p.cutoff & log2FC < -fc.cutoff), col = "green") +  
  geom_point(data = subset(data.plot1, log10.pval > p.cutoff & log2FC > fc.cutoff), col = "red") + 
  geom_hline(yintercept = p.cutoff, col = "black") +   # Horizontal Line for p-value threshold
  geom_vline(xintercept = -fc.cutoff, col = "green") +  # Vertical Line for negative Fold Change Threshold
  geom_vline(xintercept = fc.cutoff, col = "red") +     # Vertical Line for positive Fold Change Threshold
  #geom_text_repel(data = subset(data.plot, log10.pval > p.cutoff & abs(log2FC) > fc.cutoff)) +    # Add labels to significant probes. Only use this if the number of significant probes are small
  xlab("log2(Untreated/Prednisolone)") + ylab("-log10(p-value)") + # Adding appropriate labels for X and Y axes.
  ggtitle(label = "Differential Expression Analysis of Untreated/Prednisolone", 
          subtitle = "With FDR corrected p-values") +     # Plot Title
  theme_bw() +  # ggplot usually puts a grey background. This theme puts a simple white background with nice box outlines
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5)) 

untreated_treated_up <- subset(data.plot1 [order(data.plot1$pval),], 
                                  log10.pval > p.cutoff & log2FC > fc.cutoff)
write.table(x = untreated_treated_up, file = "C:/My_Work/R/Projects/untreated_prednisolone_up.txt", sep = "\t",
            row.names = F, col.names = T)

untreated_treated_down <- subset(data.plot1 [order(data.plot1$pval),], 
                                    log10.pval > p.cutoff & log2FC < -fc.cutoff)
write.table(x = untreated_treated_down, file = "C:/My_Work/R/Projects/untreated_prednisolone_down.txt", sep = "\t",
            row.names = F, col.names = T)


##Separating gene symbols from gene assignment column from fdata since we had probe ID's instead of gene symbols

data.DF5 <- data.frame(Gene.symbol = fdata$gene_assignment, Transcript.ID = fdata$transcript_cluster_id)

write.table(x = data.DF5, file = "C:/My_Work/R/Projects/gene.list.txt", sep = "\t", row.names = F, col.names = T)

All.GeneSymbols <- read.table("All.gene.csv", header = T, stringsAsFactors = F)






