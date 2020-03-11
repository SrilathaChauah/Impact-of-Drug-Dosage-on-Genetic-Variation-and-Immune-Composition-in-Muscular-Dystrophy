install.packages(VennDiagram)
library(VennDiagram)

treated_up <- read.table("untreated_prednisolone_up.csv", header = F, stringsAsFactors = F)
A <- treated_up$V1

Prednisolone_down <- read.table("untreated_prednisolone_down.csv", header = F, stringsAsFactors = F)
B <- Prednisolone_down$V1

daily_vehicle_down <- read.table("vehicle_daily_up.csv", header = F, stringsAsFactors = F)
C <- daily_vehicle_down$V1

daily_up <- read.table("vehicle_daily_down.csv", header = F, stringsAsFactors = F)
D <- daily_up$V1

weekly_vehicle_down <- read.table("vehicle_weekly_up.csv", header = F, stringsAsFactors = F)
E <- weekly_vehicle_down$V1

weekly_up <- read.table("vehicle_weekly_down.csv", header = F, stringsAsFactors = F)
G <- weekly_up$V1

daily_weekly_vehicle_down <- read.table("vehicle_daily_weekly_up.csv", header = F, stringsAsFactors = F)
H <- daily_weekly_vehicle_down$V1

daily_weekly_vehicle_up <- read.table("vehicle_daily_weekly_down.csv", header = F, stringsAsFactors = F)
I <- daily_weekly_vehicle_up$V1

weekly_daily_down <- read.table("daily_weekly_up.csv", header = F, stringsAsFactors = F)
J <- weekly_daily_down$V1

weekly_daily_up <- read.table("daily_weekly_down.csv", header = F, stringsAsFactors = F)
K <- weekly_daily_up$V1


venn.diagram(list("treated_up"=A, "daily_up"=D, "weekly_up"=G), 
             fill = c("yellow", "green","cyan"), cex = 1.5, filename = "venn.diagram_UUP.png")

venn.diagram(list("Prednisolone_down"=B, "daily_vehicle_down"=C, "weekly_vehicle_down"=E), 
             fill = c("yellow", "green","cyan"), cex = 1.5, filename = "venn.diagram_DOWN.png")

venn.diagram(list("daily_weekly_up "=J, "vehicle_daily_up"=C, "vehicle_weekly_up"=E), 
             fill = c("yellow", "green","cyan"), cex = 1.5, filename = "Vehicle_Daily_Weekly_UP.png")

venn.diagram(list("daily_weekly_down"=K, "vehicle_daily_down"=D, "vehicle_weekly_down"=G), 
             fill = c("yellow", "green","cyan"), cex = 1.5, filename = "Vehicle_Daily_Weekly_DOWN.png")

venn.diagram(list("vehicle_daily_weekly_up"=H, "vehicle_daily_weekly_down"=I), 
             fill = c("yellow", "green"), cex = 1.5, filename = "venn.diagram_vehicle_daily_weekly.png")

venn.diagram(list("Prednisolone_up"=A, "Daily_up"=C, "Weekly_up"=E, "Prednisolone_down"=B, "Daily_down"=D, 
                  "Weekly_down"=G),fill = c("red", "pink","orange","yellow", "green","cyan"), cex = 1.5, 
             filename = "venn.diagram_UP&DOWNgenes.png")
#####################################################################################################################

##Kegg Pathway Test and Enrichment analysis

library(clusterProfiler)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("
                     ")
library(org.Mm.eg.db)

data2 <- read.table("GSE95682_RPKM_-_all_genes.csv", header = T, sep = ',',stringsAsFactors = F)
All.GeneSymbols2 <- data.frame(Gene.symbol = data2$gene.y)

universe.df <- bitr(All.GeneSymbols2$Gene.symbol, 
                    fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Mm.eg.db)


DE_genes <- read.table("ALL_DE_genes1.csv", header = T, sep = ',',stringsAsFactors = F)
###Kegg Pathway Test
##daily up
daily_vehicle_up.df <- bitr(DE_genes$Daily_up,
                fromType = "SYMBOL", toType = c( "ENTREZID"), OrgDb = org.Mm.eg.db)

kk <- enrichKEGG(gene = daily_vehicle_up.df$ENTREZID, organism = 'mmu', universe=universe.df$ENTREZID, 
                 maxGSSize = 2000, minGSSize = , pvalueCutoff = 0.05)

barplot(kk, showCategory=20)

##daily down
daily_vehicle_down.df <- bitr(DE_genes$Daily_down,
                            fromType = "SYMBOL", toType = c( "ENTREZID"), OrgDb = org.Mm.eg.db)

kk <- enrichKEGG(gene = daily_vehicle_down.df$ENTREZID, organism = 'mmu', universe=universe.df$ENTREZID, 
                 maxGSSize = 2000, minGSSize = , pvalueCutoff = 0.05)

barplot(kk, showCategory=20)

##Weekly up
weekly_vehicle_up.df <- bitr(DE_genes$Weekly_up,
                            fromType = "SYMBOL", toType = c( "ENTREZID"), OrgDb = org.Mm.eg.db)

kk1 <- enrichKEGG(gene = weekly_vehicle_up.df$ENTREZID, organism = 'mmu', universe=universe.df$ENTREZID, 
                 maxGSSize = 2000, minGSSize = , pvalueCutoff = 0.05)

barplot(kk1, showCategory=20)

##Daily_weekly_vehicle_up
Daily_weekly_vehicle_up.df <- bitr(DE_genes$DW_Vehicle_up,
                             fromType = "SYMBOL", toType = c( "ENTREZID"), OrgDb = org.Mm.eg.db)

kk <- enrichKEGG(gene = Daily_weekly_vehicle_up.df$ENTREZID, organism = 'mmu', universe=universe.df$ENTREZID, 
                 maxGSSize = 2000, minGSSize = , pvalueCutoff = 0.05)

barplot(kk, showCategory=20)

##weekly_Daily_up
weekly_Daily_up.df <- bitr(DE_genes$Weekly_Daily_up,
                                   fromType = "SYMBOL", toType = c( "ENTREZID"), OrgDb = org.Mm.eg.db)

kk <- enrichKEGG(gene = weekly_Daily_up.df$ENTREZID, organism = 'mmu', universe=universe.df$ENTREZID, 
                 maxGSSize = 2000, minGSSize = , pvalueCutoff = 0.05)

barplot(kk, showCategory=20)

##weekly_Daily_down
weekly_Daily_down.df <- bitr(DE_genes$Weekly_Daily_down,
                           fromType = "SYMBOL", toType = c( "ENTREZID"), OrgDb = org.Mm.eg.db)

kk <- enrichKEGG(gene = weekly_Daily_down.df$ENTREZID, organism = 'mmu', universe=universe.df$ENTREZID, 
                 maxGSSize = 2000, minGSSize = , pvalueCutoff = 0.05)

barplot(kk, showCategory=20)

####
vehicle_daily_down.df <- bitr(vehicle_daily_down$V1,
                fromType = "SYMBOL", toType = c( "ENTREZID"), OrgDb = org.Mm.eg.db)

kk <- enrichKEGG(gene = vehicle_daily_down.df$ENTREZID, organism = 'mmu', universe=universe.df$ENTREZID, 
                 maxGSSize = 2000, minGSSize = , pvalueCutoff = 0.05)

barplot(kk, showCategory=20)

##GO Enrichment, Replace ont = "BP" with "MF" or "CC"
ego <- enrichGO(gene = vehicle_daily_down.df$ENTREZID, universe=universe.df$ENTREZID, ont = "CC", OrgDb = org.Mm.eg.db,
                maxGSSize = 2000, pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05)

barplot(ego, showCategory=30)
#####################################################################################################################
###Kegg Pathway Test
##Vehicle over Weekly
vehicle_weekly_up.df <- bitr(vehicle_weekly_up$V1,
                            fromType = "SYMBOL", toType = c( "ENTREZID"), OrgDb = org.Mm.eg.db)

kk <- enrichKEGG(gene = vehicle_weekly_up.df$ENTREZID, organism = 'mmu', universe=universe.df$ENTREZID, 
                 maxGSSize = 2000, minGSSize = , pvalueCutoff = 0.05)

barplot(kk, showCategory=20)

##GO Enrichment, Replace ont = "BP" with "MF" or "CC"
ego <- enrichGO(gene = vehicle_weekly_up.df$ENTREZID, universe=universe.df$ENTREZID, ont = "BP", OrgDb = org.Mm.eg.db,
                maxGSSize = 2000, pAdjustMethod = "BH", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05)

barplot(ego, showCategory=30)

####
vehicle_weekly_down.df <- bitr(vehicle_weekly_down$V1,
                              fromType = "SYMBOL", toType = c( "ENTREZID"), OrgDb = org.Mm.eg.db)

kk <- enrichKEGG(gene = vehicle_weekly_down.df$ENTREZID, organism = 'mmu', universe=universe.df$ENTREZID, 
                 maxGSSize = 2000, minGSSize = , pvalueCutoff = 0.05)

barplot(kk, showCategory=20)
#####################################################################################################################
###Kegg Pathway Test
##Vehicle over Daily&Weekly
vehicle_daily_weekly_up.df <- bitr(vehicle_daily_weekly_up$V1,
                             fromType = "SYMBOL", toType = c( "ENTREZID"), OrgDb = org.Mm.eg.db)

kk <- enrichKEGG(gene = vehicle_daily_weekly_up.df$ENTREZID, organism = 'mmu', universe=universe.df$ENTREZID, 
                 maxGSSize = 2000, minGSSize = , pvalueCutoff = 0.05)

barplot(kk, showCategory=20)
####
vehicle_daily_weekly_down.df <- bitr(vehicle_daily_weekly_down$V1,
                               fromType = "SYMBOL", toType = c( "ENTREZID"), OrgDb = org.Mm.eg.db)

kk <- enrichKEGG(gene = vehicle_daily_weekly_down.df$ENTREZID, organism = 'mmu', universe=universe.df$ENTREZID, 
                 maxGSSize = 2000, minGSSize = , pvalueCutoff = 0.05)

barplot(kk, showCategory=20)
#####################################################################################################################
###Kegg Pathway Test
##Daily overWeekly
daily_weekly_up.df <- bitr(daily_weekly_up$V1,
                                   fromType = "SYMBOL", toType = c( "ENTREZID"), OrgDb = org.Mm.eg.db)

kk <- enrichKEGG(gene = daily_weekly_up.df$ENTREZID, organism = 'mmu', universe=universe.df$ENTREZID, 
                 maxGSSize = 2000, minGSSize = , pvalueCutoff = 0.05)

barplot(kk, showCategory=20)
####
daily_weekly_down.df <- bitr(daily_weekly_down$V1,
                                     fromType = "SYMBOL", toType = c( "ENTREZID"), OrgDb = org.Mm.eg.db)

kk <- enrichKEGG(gene = daily_weekly_down.df$ENTREZID, organism = 'mmu', universe=universe.df$ENTREZID, 
                 maxGSSize = 2000, minGSSize = , pvalueCutoff = 0.05)

barplot(kk, showCategory=20)
#####################################################################################################################
