###Genes in renal Cancer - Dataset GSE66272)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(version = "3.11")
BiocManager::install(c('affy','affyPLM', 'sva', 'AnnotationDbi', 'hgu133plus2.db', 'simpleaffy', 'arrayQualityMetrics','affyQCReport', 'gcrma', 'ggplot2'))
install.packages("readtable")


library(sva)
library(readxl)
library(affy)
library(simpleaffy)
library("arrayQualityMetrics")
library(affyQCReport)
library(preprocessCore)
library(affyPLM)
library(ggplot2)
library(pheatmap)

#read cel files from the computer, create an affybatch
dataset <- ReadAffy(celfile.path="/Users/sebi/Desktop/Final Project/data/GSE66272_RAW/")

#express data as a data frame 
frame <- as.data.frame(exprs(dataset))

#create QC Plot, each array being represented by a separate line in the figure.
plot(qc(dataset))

#create arrayQualityMetrics data sets to examine the info from dataset.
arrayQualityMetrics(expressionset = dataset, outdir = "Report_for_dataset", force = TRUE, do.logtransform = TRUE)

#create affyQCReport
QCReport(dataset,file="QCReport.pdf")

fset<- fitPLM(dataset) ## fitthedata
png("/Users/sebi/Desktop/Final Project/RLE_combined_datasets.png", width = 12, height = 9, units = 'in', res = 300)
RLE(fset, main = "RLE for CSE65272 dataset" )
dev.off()
## generating a png file for outputing the result of RLE

RLE_Stat = RLE(fset, type = "stats" )
png("/Users/sebi/Desktop/Final Project/RLE_combined_histogram.png", width = 12, height = 9, units = 'in', res = 300)
hist(RLE_Stat, main = "Histogram of relative log expression - Median", xlab = "Median", cex.main = 2, cex.axis = 1.5, cex.lab = 1.5, col = "green")
dev.off()
## generating the histogram for RLE for the combined data and saving it as a png file


png("/Users/sebi/Desktop/Final Project/NUSE_dataset.png", width = 12, height = 9, units = 'in', res = 300)
NUSE(fset, main = "Normalized unscaled standard error for GSE65272 datasets")
dev.off()

png("/Users/sebi/Desktop/Final Project//NUSE_histogram.png", width = 12, height = 9, units = 'in', res = 300)
NUSE_Stat = NUSE(fset, type = "stats")
NUSE(fset, main="NUSE boxplot PLM", las=2, cex.names=0.1)
hist(NUSE_Stat[1,], breaks=20, main = "Histogram of normalized Unscaled standard error - Median",
     xlab = "NUSE values", col = "red")
dev.off()

NUSE_statistics = NUSE(fset,type="stats")
Nuse_S <- data.frame(Median = NUSE_statistics[1, ]ggplot(Nuse_S, aes(Median))+ geom_histogram()+ xlab("Median")+ ylab("Frequency")+ ggtitle("NUSE Histogram")+ theme_bw()
                     par(mai=c(3.5,1,1,1))  NUSE(pset, main="NUSE", las=2)

##### Data normalization using mas5 package
## visulaize the data before normalization
png("/Users/sebi/Desktop/Final Project/boxplots_before_norm.png", width = 12, height = 9, units = 'in', res = 300)
x_labs<- factor(rownames(pData(dataset))) 
x_labs<- gsub("\\_.+", "", x_labs)
par(mar=c(15,2,1,1))
cols<- c("blue", "green", "red", "brown")
boxplot(dataset, xaxt = "n", col =cols, main = "Gene Expression data before normalization")
text(seq_along(x_labs), par("usr")[3] - 0.5, labels = x_labs, srt = 90, adj = 1, xpd = TRUE)
dev.off()

## data normalization using mas5 function
normdata<- mas5(dataset)

## visualise data after normalization 
png("/Users/sebi/Desktop/Final Project/boxplots_after_norm.png", width = 12, height = 9, units = 'in', res = 300)
x_labs<- factor(rownames(pData(normdata))) ## extracting the annotation data
x_labs<- gsub("\\_.+", "", x_labs)## removing the parts not necessary
par(mar=c(15,2,1,1)) ## specifying the ploting area
cols<- c("blue", "green", "red", "brown") ## chosing the colour 
boxplot(normdata, xaxt = "n", col =cols, main = "Gene Expression data after normalization")
text(seq_along(x_labs), par("usr")[3] - 0.5, labels = x_labs, srt = 90, adj = 1, xpd = TRUE)
dev.off()

# probing the normalized data
normdata <- log2(exprs(normdata))
write.csv(normdata, file="Mas5data.csv")
boxplot(normdata, main = "Gene expression level after Mas5 normalization", las=2)

#perform another QC analysis after normalization
arrayQualityMetrics(expressionset = normdata, 
                    outdir = "Users/sebi/Desktop/Final Project/report_norm_dataset",
                    force = T, do.logtransform = T)


### Batch correction - we use combat just when we have multiple sets of data
#meta <- read.csv("/Users/sebi/Desktop/Final Project/Metadata.csv")
#design <- model.matrix(~CN, meta)
#batch = meta$BATCH
#comdata <- ComBat(normdata, meta$BATCH, design)
write.csv(normdata,"/Users/sebi/Desktop/Final Project/normdata.csv")


mas5 <- read.csv("/Users/sebi/Desktop/Final Project/normdata.csv")

row.names(mas5) <- mas5$X
rma <- mas5[-1]


### ANNOTATION
library(WGCNA)
library(limma)
library(EnhancedVolcano)
library(ggplot2)
library(pheatmap)
annot <- function(d) {
  symbols <- AnnotationDbi::select(hgu133plus2.db, keys=rownames(d), columns=c("SYMBOL"), keytype="PROBEID")
  symbols <- symbols[!duplicated(symbols$PROBEID),]
  row.names(symbols) <- symbols$PROBEID
  d <- merge(symbols, d, by="row.names")
  row.names(d) <- d$Row.names
  d <- d[-c(1,2)]
  d <- na.omit(d)
  d <- collapseRows(d[-1], rowGroup=d$SYMBOL, rowID=rownames(d))
  return(d)
}

annotRma <- annot(rma)


### GENE FILTERING
filt <- function(d) {
  means <- rowMeans(d)
  perc2 <- as.numeric(quantile(means, probs=0.02, na.rm=T))
  temp <- d[which(means > perc2),]
  return(temp)
}

filtRma <- filt(annotRma$datETcollapsed)


### LIMMA
type <- factor(meta$CN, levels = c('normal','ccRCC'), ordered = F)
row.names(meta) <- meta$Sample_geo_accession
design <- model.matrix(~0+type, meta)

max <- 50000

limmaFit <- function(d) {
  lm <- lmFit(d, design)
  contrast <- makeContrasts(typeccRCC-typenormal, levels=design)
  colnames(contrast) <- c("Cancer-Control")
  fit <- contrasts.fit(lm, contrast)
  fit <- eBayes(fit)
  return(fit)
}

geneRma <- rownames(filtRma)
fitRma <- limmaFit(filtRma)

tTRma <- topTable(fitRma, coef=1, p.value=0.05, sort="P", genelist=geneRma, number=max)
top100 <- topTable(fitRma, coef=1, p.value=0.05, sort="P", genelist=geneRma, number=100)
write.table(tTRma, "/Users/sebi/Desktop/Final Project/fitdata.txt")
write.csv(tTRma, "/Users/sebi/Desktop/Final Project/fitdata.csv")

length(rownames(tTRma[which(tTRma$adj.P.Val < 0.05),]))
fitdata <- tTRma

EnhancedVolcano(fitdata, lab=row.names(fitdata), x="logFC", y="adj.P.Val", 
                pointSize=1, legendLabSize=10, labSize=3.0,
                title="Enhanced Volcano Plot", subtitle="Genes in Kidney Cancer")


## HEATMAP
hm <- function(dFit, dFilt, gl, t) {
  t50 <- topTable(dFit, genelist=gl, adjust.method="fdr", sort.by="p.value", number=50)
  input <- dFilt[(rownames(dFilt) %in% rownames(t50)),]
  group <- data.frame(sample=factor(meta$CN))
  row.names(group) <- colnames(dFilt)
  par(2,2,2,2)
  pheatmap(input, main=t, cluster_rows=F, annotation_col=group)
}

hm(fitRma, filtRma, geneRma, "Heatmap after Limma Analysis")


#Functional Analysis
library(topGO)
library(clusterProfiler)
library(pathview)
library(magrittr)
library(tidyr)
library(msigdbr)
library(enrichplot)

#We are working with our normalized data from mas5 file.
### FUNCTIONAL ANALYSIS
## DEG VECTOR
genelist <- function(d) {
  gene <- d$logFC
  names(gene) <- d$ID
  DEG.gene <- gene[gene > 1.5]
  DEG.gene <- sort(DEG.gene, decreasing=T)
  return(DEG.gene)
}

rGene <- genelist(tTRma)
DEG.r <- AnnotationDbi::select(org.Hs.eg.db, keys=names(rGene), columns=c("ENTREZID"), keytype="SYMBOL")


# GENE ONTOLOGY
eGoD <- function(d, o, f) {
  ego <- enrichGO(d$ENTREZID, org.Hs.eg.db, ont=o, readable=T)
  write.csv(ego, f)
}

eGoD(DEG.r, "ALL", "/Users/sebi/Desktop/Final Project/mas5CC.csv")
eGoD(DEG.r, "BP", "./FA/data/GO/rmaBP.csv")
eGoD(DEG.r, "MF", "./FA/data/GO/rmaMF.csv")

egoPlot  <- function(d, o, t) {
  ego <- enrichGO(d$ENTREZID, org.Hs.eg.db, ont=o, readable=F)
  egoR <- setReadable(ego, org.Hs.eg.db, keyType="ENTREZID")
  barplot(egoR, title=t, drop = T,showCategory=20)
}

egoPlot(DEG.r, "MF", "EnrichmentGO MF")
egoPlot(DEG.r, "BP", "EnrichmentGo BP")
egoPlot(DEG.r, "CC", "EnrichmentGo CC")

# OPTIONAL***
ggoPlot <- function(d, o, l, t) {
  ggo <- groupGO(d$ENTREZID, org.Hs.eg.db, keyType="ENTREZID", ont=o, level=l, readable=F)
  barplot(ggo, title=t, showCategory=10)
}

ggoPlot(DEG.r, "MF", 10, "Mas5 MF")
### KEGG ANALYSIS
keggPlot <- function(d, gene, c, t) {
  ek <- enrichKEGG(d$ENTREZID)
  dotplot(ek, title=t)
}

keggPlot(DEG.r, rGene, T, "KEGG pathway enrichment analysis")


### GENE CONCEPT NETWORK
cnetPlot <- function(d, gene, c) {
  # ed <- enrichDGN(d$ENTREZID)
  ek <- enrichKEGG(d$ENTREZID)
  edR <- setReadable(ek, org.Hs.eg.db, keyType="ENTREZID")
  cnetplot(edR, foldChange=gene, categorySize="pvalue", colorEdge=T, vertex.label.font=c, showCategory =7)
}

cnetPlot(DEG.r, rGene, 6)
cnetPlot(DEG.r, rGene, T)


## TRANSCRIPTIONAL FACTOR ANALYSIS
tfa <- function(d, gene) {
  c3 <- read.gmt("./FA/data/c3.gmt")
  e <- enricher(DEG.m$ENTREZID, TERM2GENE=c3)
  temp <- setReadable(e, org.Hs.eg.db, keyType="ENTREZID")
  cnetplot(temp, foldChange=gene, categorySize="pvalue", colorEdge=T)
}

tfa(DEG.r, rGene)

###GSEA ANALYSIS
gsealist <- function(d) {
  gene <- d$logFC
  names(gene) <- d$ID
  
  DEG <- AnnotationDbi::select(org.Hs.eg.db, keys=names(gene), columns=c("ENTREZID"), keytype="SYMBOL")
  DEG <- DEG[!duplicated(DEG$SYMBOL),]
  row.names(DEG) <- DEG$SYMBOL
  DEG <- merge(DEG, gene, by="row.names")
  
  genelist <- DEG$y
  names(genelist) <- DEG$ENTREZID
  genelist <- sort(genelist, decreasing=T)
  
  return(genelist)
}

mGene <- gsealist(fitdata)
rGene <- gsealist(rma)
gGene <- gsealist(gcrma)

msig <- msigdbr(species="Homo sapiens", category="H")
msigS <- msig %>% select(gs_name, entrez_gene)

gseaM <- GSEA(mGene, TERM2GENE=msigS)


gseaplot2(gseaM, geneSetID=1:6, pvalue=gseaM$pvalue)
install.packages("ridgeplot")
ridgeplot(gseaM, scale=15) + labs(x = "enrichment distribution")

e <- enricher(names(mGene), TERM2GENE=msigS)
eM <- setReadable(e, org.Hs.eg.db, keyType="ENTREZID")
heatplot(eM, foldChange=mGene)

#Create treemap with revigo

revigo.names <- c("term_ID","description","freqInDbPercent","abslog10pvalue","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0001906","cell killing",0.629,10.1107,0.994,0.000,"cell killing"),
                     c("GO:0001909","leukocyte mediated cytotoxicity",0.473,10.1226,0.960,0.000,"leukocyte mediated cytotoxicity"),
                     c("GO:0006968","cellular defense response",0.340,19.0731,0.922,0.000,"cellular defense response"),
                     c("GO:0070482","response to oxygen levels",2.089,6.2062,0.946,0.646,"cellular defense response"),
                     c("GO:0001666","response to hypoxia",1.945,6.7986,0.923,0.253,"cellular defense response"),
                     c("GO:0042110","T cell activation",2.637,35.3799,0.562,0.000,"T cell activation"),
                     c("GO:0002467","germinal center formation",0.087,4.7167,0.735,0.627,"T cell activation"),
                     c("GO:0002460","adaptive immune response based on somatic recombination of immune receptors built from immunoglobulin superfamily domains",1.881,13.6536,0.673,0.396,"T cell activation"),
                     c("GO:0050900","leukocyte migration",2.649,24.8794,0.656,0.415,"T cell activation"),
                     c("GO:0002440","production of molecular mediator of immune response",1.333,4.9706,0.733,0.379,"T cell activation"),
                     c("GO:0048002","antigen processing and presentation of peptide antigen",1.062,12.2765,0.719,0.368,"T cell activation"),
                     c("GO:0034341","response to interferon-gamma",0.935,13.0400,0.667,0.613,"T cell activation"),
                     c("GO:0002532","production of molecular mediator involved in inflammatory response",0.283,5.0991,0.865,0.633,"T cell activation"),
                     c("GO:0042092","type 2 immune response",0.196,5.0000,0.749,0.305,"T cell activation"),
                     c("GO:1903706","regulation of hemopoiesis",1.823,12.1141,0.562,0.580,"T cell activation"),
                     c("GO:0042116","macrophage activation",0.329,5.0410,0.694,0.694,"T cell activation"),
                     c("GO:0050728","negative regulation of inflammatory response",0.600,4.5544,0.785,0.685,"T cell activation"),
                     c("GO:0050727","regulation of inflammatory response",1.685,12.4547,0.802,0.495,"T cell activation"),
                     c("GO:0007204","positive regulation of cytosolic calcium ion concentration",1.529,8.6925,0.911,0.505,"T cell activation"),
                     c("GO:1990869","cellular response to chemokine",0.012,10.1945,0.932,0.438,"T cell activation"),
                     c("GO:0019882","antigen processing and presentation",1.339,14.0910,0.733,0.379,"T cell activation"),
                     c("GO:0002828","regulation of type 2 immune response",0.156,4.6289,0.679,0.442,"T cell activation"),
                     c("GO:0007162","negative regulation of cell adhesion",1.321,11.7496,0.801,0.566,"T cell activation"),
                     c("GO:0071346","cellular response to interferon-gamma",0.808,12.9136,0.666,0.428,"T cell activation"),
                     c("GO:0036336","dendritic cell migration",0.173,5.0009,0.707,0.651,"T cell activation"),
                     c("GO:0072678","T cell migration",0.237,9.0482,0.689,0.670,"T cell activation"),
                     c("GO:0033627","cell adhesion mediated by integrin",0.317,5.3862,0.908,0.480,"T cell activation"),
                     c("GO:0097028","dendritic cell differentiation",0.242,5.1662,0.688,0.611,"T cell activation"),
                     c("GO:0046629","gamma-delta T cell activation",0.081,4.9245,0.673,0.658,"T cell activation"),
                     c("GO:0002683","negative regulation of immune system process",2.152,21.1918,0.612,0.403,"T cell activation"),
                     c("GO:0034113","heterotypic cell-cell adhesion",0.289,5.1518,0.873,0.640,"T cell activation"),
                     c("GO:0002703","regulation of leukocyte mediated immunity",0.912,14.0872,0.625,0.696,"T cell activation"),
                     c("GO:0002228","natural killer cell mediated immunity",0.335,6.5072,0.677,0.536,"T cell activation"),
                     c("GO:0002697","regulation of immune effector process",1.922,14.8416,0.622,0.584,"T cell activation"),
                     c("GO:0002695","negative regulation of leukocyte activation",0.837,19.8996,0.565,0.587,"T cell activation"),
                     c("GO:0031589","cell-substrate adhesion",1.783,5.6198,0.893,0.588,"T cell activation"),
                     c("GO:0045576","mast cell activation",0.340,6.8794,0.693,0.697,"T cell activation"),
                     c("GO:0002291","T cell activation via T cell receptor contact with antigen bound to MHC molecule on antigen presenting cell",0.046,4.5702,0.649,0.626,"T cell activation"),
                     c("GO:0048872","homeostasis of number of cells",1.339,9.0367,0.940,0.446,"T cell activation"),
                     c("GO:0002757","immune response-activating signal transduction",3.405,9.6799,0.543,0.631,"T cell activation"),
                     c("GO:0002761","regulation of myeloid leukocyte differentiation",0.641,5.7570,0.579,0.676,"T cell activation"),
                     c("GO:0043374","CD8-positive, alpha-beta T cell differentiation",0.063,4.9245,0.626,0.644,"T cell activation"),
                     c("GO:0001774","microglial cell activation",0.104,6.3325,0.721,0.470,"T cell activation"),
                     c("GO:0001776","leukocyte homeostasis",0.490,10.4101,0.723,0.337,"T cell activation"),
                     c("GO:0034612","response to tumor necrosis factor",1.673,5.9747,0.914,0.664,"T cell activation"),
                     c("GO:0006959","humoral immune response",1.702,5.9957,0.696,0.467,"T cell activation"),
                     c("GO:0070098","chemokine-mediated signaling pathway",0.467,10.2062,0.871,0.598,"T cell activation"),
                     c("GO:0043312","neutrophil degranulation",2.799,14.9747,0.533,0.685,"T cell activation"),
                     c("GO:0030595","leukocyte chemotaxis",1.108,14.9547,0.613,0.589,"T cell activation"),
                     c("GO:0046718","viral entry into host cell",0.589,5.5421,0.983,0.002,"viral entry into host cell"),
                     c("GO:0019058","viral life cycle",2.366,4.6596,0.984,0.691,"viral entry into host cell"),
                     c("GO:0000280","nuclear division",3.508,11.6003,0.903,0.003,"nuclear division"),
                     c("GO:0045931","positive regulation of mitotic cell cycle",0.756,5.3279,0.799,0.493,"nuclear division"),
                     c("GO:0098883","synapse disassembly",0.006,8.8697,0.955,0.101,"nuclear division"),
                     c("GO:0043062","extracellular structure organization",1.835,11.5376,0.934,0.180,"nuclear division"),
                     c("GO:0030198","extracellular matrix organization",1.829,11.5817,0.934,0.180,"nuclear division"),
                     c("GO:0034508","centromere complex assembly",0.231,4.7645,0.922,0.497,"nuclear division"),
                     c("GO:0051383","kinetochore organization",0.104,5.9431,0.926,0.259,"nuclear division"),
                     c("GO:0048285","organelle fission",3.774,10.0114,0.941,0.399,"nuclear division"),
                     c("GO:0007077","mitotic nuclear envelope disassembly",0.254,5.4045,0.868,0.431,"nuclear division"),
                     c("GO:0045730","respiratory burst",0.162,7.7235,0.977,0.011,"respiratory burst"),
                     c("GO:0070661","leukocyte proliferation",1.570,26.3478,0.939,0.014,"leukocyte proliferation"),
                     c("GO:0006909","phagocytosis",1.858,9.9031,0.921,0.014,"phagocytosis"),
                     c("GO:0051235","maintenance of location",1.702,5.0079,0.913,0.162,"phagocytosis"),
                     c("GO:0050764","regulation of phagocytosis",0.381,8.4685,0.850,0.635,"phagocytosis"),
                     c("GO:0051924","regulation of calcium ion transport",1.206,5.5114,0.886,0.447,"phagocytosis"),
                     c("GO:0002791","regulation of peptide secretion",1.166,5.5287,0.869,0.681,"phagocytosis"),
                     c("GO:0007599","hemostasis",1.973,4.7282,0.943,0.319,"phagocytosis"),
                     c("GO:1903531","negative regulation of secretion by cell",1.004,6.0287,0.820,0.314,"phagocytosis"),
                     c("GO:0001819","positive regulation of cytokine production",2.256,14.1884,0.726,0.014,"positive regulation of cytokine production"),
                     c("GO:0035987","endodermal cell differentiation",0.260,4.8794,0.896,0.527,"positive regulation of cytokine production"),
                     c("GO:0048771","tissue remodeling",0.837,6.9031,0.928,0.158,"positive regulation of cytokine production"),
                     c("GO:0032611","interleukin-1 beta production",0.329,9.0164,0.824,0.665,"positive regulation of cytokine production"),
                     c("GO:0032612","interleukin-1 production",0.398,8.6904,0.836,0.679,"positive regulation of cytokine production"),
                     c("GO:0032623","interleukin-2 production",0.329,7.5560,0.838,0.665,"positive regulation of cytokine production"),
                     c("GO:0032615","interleukin-12 production",0.294,5.5622,0.840,0.657,"positive regulation of cytokine production"),
                     c("GO:0032613","interleukin-10 production",0.254,9.3045,0.841,0.647,"positive regulation of cytokine production"),
                     c("GO:0045765","regulation of angiogenesis",1.264,5.2596,0.833,0.390,"positive regulation of cytokine production"),
                     c("GO:0032637","interleukin-8 production",0.392,5.1487,0.836,0.678,"positive regulation of cytokine production"),
                     c("GO:0032653","regulation of interleukin-10 production",0.242,8.7033,0.802,0.666,"positive regulation of cytokine production"),
                     c("GO:0032655","regulation of interleukin-12 production",0.283,5.7399,0.802,0.677,"positive regulation of cytokine production"),
                     c("GO:0032663","regulation of interleukin-2 production",0.289,6.6021,0.802,0.679,"positive regulation of cytokine production"),
                     c("GO:0032677","regulation of interleukin-8 production",0.358,5.6968,0.796,0.694,"positive regulation of cytokine production"),
                     c("GO:0031100","animal organ regeneration",0.450,5.0306,0.911,0.148,"positive regulation of cytokine production"),
                     c("GO:0018108","peptidyl-tyrosine phosphorylation",2.135,6.6737,0.975,0.030,"peptidyl-tyrosine phosphorylation"),
                     c("GO:0018212","peptidyl-tyrosine modification",2.147,6.5817,0.978,0.647,"peptidyl-tyrosine phosphorylation"),
                     c("GO:0042107","cytokine metabolic process",0.635,6.0150,0.986,0.174,"peptidyl-tyrosine phosphorylation"),
                     c("GO:0070371","ERK1 and ERK2 cascade",1.564,4.5406,0.878,0.653,"peptidyl-tyrosine phosphorylation"),
                     c("GO:0070374","positive regulation of ERK1 and ERK2 cascade",1.068,5.9101,0.809,0.472,"peptidyl-tyrosine phosphorylation"),
                     c("GO:0051304","chromosome separation",0.404,10.8861,0.851,0.038,"chromosome separation"),
                     c("GO:0071173","spindle assembly checkpoint",0.173,8.2020,0.832,0.349,"chromosome separation"),
                     c("GO:0000075","cell cycle checkpoint",1.310,4.8697,0.812,0.616,"chromosome separation"),
                     c("GO:0031577","spindle checkpoint",0.225,8.2020,0.831,0.687,"chromosome separation"),
                     c("GO:0007051","spindle organization",0.912,5.5850,0.861,0.464,"chromosome separation"),
                     c("GO:2000106","regulation of leukocyte apoptotic process",0.473,8.9666,0.932,0.039,"regulation of leukocyte apoptotic process"),
                     c("GO:0071887","leukocyte apoptotic process",0.589,7.2541,0.953,0.384,"regulation of leukocyte apoptotic process"),
                     c("GO:0033032","regulation of myeloid cell apoptotic process",0.150,4.7471,0.940,0.388,"regulation of leukocyte apoptotic process"),
                     c("GO:0006801","superoxide metabolic process",0.329,5.0306,0.986,0.044,"superoxide metabolism"),
                     c("GO:0042554","superoxide anion generation",0.138,4.8894,0.986,0.688,"superoxide metabolism"),
                     c("GO:0007059","chromosome segregation",1.939,9.1451,0.955,0.046,"chromosome segregation"),
                     c("GO:0071222","cellular response to lipopolysaccharide",0.871,15.5969,0.870,0.069,"cellular response to lipopolysaccharide"),
                     c("GO:0019722","calcium-mediated signaling",0.877,7.4260,0.903,0.355,"cellular response to lipopolysaccharide"),
                     c("GO:0071216","cellular response to biotic stimulus",1.033,16.3224,0.912,0.606,"cellular response to lipopolysaccharide"),
                     c("GO:0019932","second-messenger-mediated signaling",1.339,7.5186,0.899,0.120,"cellular response to lipopolysaccharide"),
                     c("GO:0009595","detection of biotic stimulus",0.115,7.7670,0.934,0.491,"cellular response to lipopolysaccharide"),
                     c("GO:0007229","integrin-mediated signaling pathway",0.548,5.9066,0.905,0.152,"cellular response to lipopolysaccharide"),
                     c("GO:0007249","I-kappaB kinase/NF-kappaB signaling",1.535,5.0334,0.898,0.378,"cellular response to lipopolysaccharide"),
                     c("GO:0009615","response to virus",1.771,11.8928,0.884,0.656,"cellular response to lipopolysaccharide"),
                     c("GO:0002831","regulation of response to biotic stimulus",0.802,8.7645,0.862,0.590,"cellular response to lipopolysaccharide"),
                     c("GO:0002833","positive regulation of response to biotic stimulus",0.231,5.2027,0.842,0.522,"cellular response to lipopolysaccharide"),
                     c("GO:1990868","response to chemokine",0.012,10.1945,0.941,0.218,"cellular response to lipopolysaccharide"),
                     c("GO:0098581","detection of external biotic stimulus",0.098,6.3391,0.910,0.493,"cellular response to lipopolysaccharide"),
                     c("GO:0051607","defense response to virus",1.327,12.6308,0.626,0.689,"cellular response to lipopolysaccharide"),
                     c("GO:0035590","purinergic nucleotide receptor signaling pathway",0.127,5.3507,0.916,0.261,"cellular response to lipopolysaccharide"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$abslog10pvalue <- as.numeric( as.character(stuff$abslog10pvalue) );
stuff$freqInDbPercent <- as.numeric( as.character(stuff$freqInDbPercent) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="revigo_treemap2.pdf", width=16, height=9 ) # width and height are in inches

# check the tmPlot command documentation for all possible parameters - there are a lot more
treemap(
  stuff,
  index = c("representative","description"),
  vSize = "abslog10pvalue",
  type = "categorical",
  vColor = "representative",
  title = "REVIGO Gene Ontology treemap",
  inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
  lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
  bg.labels = "#CCCCCCAA",     # define background color of group labels
  # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
  position.legend = "none"
)

dev.off()



