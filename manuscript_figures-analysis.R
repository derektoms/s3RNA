# __________________________________________________/\\\\\\_____/\\\\\\_________________/\\\\\\\\\______/\\\\\_____/\\\_____/\\\\\\\\\____
#  _________________________________________________\////\\\____\////\\\_______________/\\\///////\\\___\/\\\\\\___\/\\\___/\\\\\\\\\\\\\__
#   ____________________________________________________\/\\\_______\/\\\______________\/\\\_____\/\\\___\/\\\/\\\__\/\\\__/\\\/////////\\\_
#    __/\\\\\\\\\\____/\\\\\__/\\\\\____/\\\\\\\\\_______\/\\\_______\/\\\______________\/\\\\\\\\\\\/____\/\\\//\\\_\/\\\_\/\\\_______\/\\\_
#     _\/\\\//////___/\\\///\\\\\///\\\_\////////\\\______\/\\\_______\/\\\______________\/\\\//////\\\____\/\\\\//\\\\/\\\_\/\\\\\\\\\\\\\\\_
#      _\/\\\\\\\\\\_\/\\\_\//\\\__\/\\\___/\\\\\\\\\\_____\/\\\_______\/\\\______________\/\\\____\//\\\___\/\\\_\//\\\/\\\_\/\\\/////////\\\_
#       _\////////\\\_\/\\\__\/\\\__\/\\\__/\\\/////\\\_____\/\\\_______\/\\\______________\/\\\_____\//\\\__\/\\\__\//\\\\\\_\/\\\_______\/\\\_
#        __/\\\\\\\\\\_\/\\\__\/\\\__\/\\\_\//\\\\\\\\/\\__/\\\\\\\\\__/\\\\\\\\\___________\/\\\______\//\\\_\/\\\___\//\\\\\_\/\\\_______\/\\\_
#         _\//////////__\///___\///___\///___\////////\//__\/////////__\/////////____________\///________\///__\///_____\/////__\///________\///__

## 2020-01-11
## This will serve as the final code / dataset for the manuscript
## Read data loaded from s3RNA package
## Last update: 2020-01-30

#   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-
#  / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \
# `-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `
#
# devtools::install_github('kevinblighe/EnhancedVolcano')
    
# Libraries
library(devtools)
library(dplyr)
library(tidyr)
library(stringr)
library(hexbin)
library(rnaseqGene) # loads 'ggplot2' and 'RColorBrewer'
library(DESeq2)
library(vsn)
library(data.table)
library(VennDiagram)
library(ggplots2)
library(BioCircos)
library(RColorBrewer)
library(rmarkdown)
library(EnhancedVolcano)
library(s3RNA)
library(treemap)
library(gridExtra)
library(cowplot)

#   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-
#  / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \
# `-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `

setwd("~/Documents/Manuscripts/Nuclear microRNA/")
PATH <- "Figures 2020-01-31/R output/"

# Load data and define two-factor model (#1)
miRReads$col.dat$Batch<-factor(miRReads$col.dat$Batch) ## this should be fixed in the package!
mir <- DESeqDataSetFromMatrix(countData = miRReads$raw.count, colData = miRReads$col.dat, design=~size*subcell+Batch)
sno <- DESeqDataSetFromMatrix(countData = snoRReads$raw.count, colData = snoRReads$col.dat, design=~size*subcell+batch)

mirD1 <- DESeq(mir)
snoD1 <- DESeq(sno)

## Group-based model (#2)
### As recommended online, use size and subcellular localization to make groups for contrasts
mir$group <- factor(paste0(mir$size,mir$subcell))
sno$group <- factor(paste0(sno$size,sno$subcell))

design(mir) <- ~group+Batch
design(sno) <- ~group+batch

mirD2 <- DESeq(mir)
snoD2 <- DESeq(sno)

#   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-
#  / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \
# `-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `
#
# Colour space
set1<-brewer.pal(4,"Set1")
lite1<-brewer.pal(4,"Pastel1")

ann_col <- list(size=c(small=set1[1],large=set1[2]),subcell=c(nucleus=set1[3],cytosol=set1[4]))
grad_col <- colorRampPalette(brewer.pal(9, "YlOrBr"))(100)

desat = function(cols, sat=0.5) {
    X <- diag(c(1, sat, 1)) %*% rgb2hsv(col2rgb(cols))
    hsv(X[1,], X[2,], X[3,])
}


class_col <- desat(c("#b7c572", "#d1bc46", "#f99575", "#b4ca48", "#69d289", "#e7a961", "#eca642", "#89cc63", "#cbc05e", "#4dd8bd", "#cfad28", "#4eb4f3", "#dd9aef", "#f392c4"),0.6) ## palette created with https://medialab.github.io/iwanthue/ using colour blind-friendly colours
## they line up nice with the groups and result in all-black text over top of each box.
#
#   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-
#  / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \
# `-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `
#
# Overview read analysis
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# fig2.dist <- read.csv(file = "2020-01-21 submission/Fig 2 data-dist.csv", header = TRUE, stringsAsFactors = FALSE)
# fig2.spec <- read.csv(file = "2020-01-21 submission/Fig 2 data-species.csv", header = TRUE, stringsAsFactors = FALSE)
# fig2 <- list("read_dist"=fig2.dist,"species"=fig2.spec)
# save(fig2, "2020-01-22_Figure2_data.rds")

load("2020-01-22_Figure2_data.rds")

# Figure 2A
## Size distribution
read_dist <- fig2$read_dist %>% gather(group,percent,SGC.nuc:LGC.cyt, factor_key=TRUE) %>% separate(group, into=c("size","subcell"))
read_dist$size <- factor(read_dist$size,levels=c("SGC","LGC")) # order for plotting
read_dist$subcell <- factor(read_dist$subcell, levels=c("nuc","cyt")) # order for plotting
ggplot(read_dist,aes(x=length,y=percent)) + geom_histogram(stat="identity",aes(fill=subcell,alpha=size),position=position_dodge()) + theme_bw(base_size=14) + theme(legend.position="bottom") + scale_y_continuous(trans="log1p") + scale_colour_manual(values=paste(ann_col$subcell),aesthetics="fill") + scale_alpha_manual(values=c(0.6,1.0)) + theme(legend.position=c(0.6,0.8),legend.box="horizontal",legend.direction="horizontal")
ggsave("seq_size_dist.pdf", width=9, height=3, units="in", path=PATH)

# Figure 2B
## Tree maps of RNA species
read_spec <- fig2$species
read_spec$col <- 1:14 # assign specific colours to each RNA species
read_spec$species[1]<-"lncRNA"
read_spec$species[2]<-"misc. RNA\nspecies"
read_spec$species[3]<-"rRNA\n(mitochondrial)"
read_spec$species[4]<-"tRNA\n(mitochondrial)"
read_spec$species[13]<-"miRNA hairpin\nprecursors"
read_spec$species[14]<-"mature miRNA"

pdf(paste0(PATH,"cyt_treemap.pdf"), width=6,height=10)
treemap(read_spec, index=c("species"), vSize="pc.CYT", title="Cytoplasmic small RNA composition (%)", vColor="col", type="value", palette=class_col, lowerbound.cex.labels=0) # cytoplamic
dev.off()

pdf(paste0(PATH,"nuc_treemap.pdf"), width=6,height=10)
treemap(read_spec, index=c("species"), vSize="pc.NUC",  title="Nuclear small RNA composition (%)", vColor="col", type="value", palette=class_col, lowerbound.cex.labels=0) # nuclear
dev.off()

#
#   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-
#  / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \
# `-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `
#
# MicroRNA analysis
## Full linear model, two factor
sizes<-results(mirD1, contrast=c("size","small","large")) # 0
locals<-results(mirD1, contrast=c("subcell","cytosol","nucleus")) # 83

## Specific groupings (comparable to Ram's analysis)
nuc.sizes<-results(mirD2, contrast=c("group","largenucleus","smallnucleus")) # 7; we got 378 back!
cyt.sizes<-results(mirD2, contrast=c("group","largecytosol","smallcytosol")) # 0
lg.local<-results(mirD2, contrast=c("group","largenucleus","largecytosol")) # 83
sm.local<-results(mirD2, contrast=c("group","smallnucleus","smallcytosol")) # 101
#
sizes[which(sizes$padj<0.1),] # there are no significant values! most are .99
locals[which(locals$padj<0.1),] # but plenty of significant differences here! 83 with padj<0.1
#
# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#
# Figure 3A
## Clustered heatmap
#
mirRLOG <-rlog(mirD1,blind=FALSE)
anno <- as.data.frame(colData(mirRLOG)[, c("size","subcell")])

topvarmiR<-head(order(rowVars(assay(mirRLOG)),decreasing=TRUE),50)
mat <- assay(mirRLOG)[topvarmiR, ]
mat <- mat - rowMeans(mat)
mirRow <- str_remove(rownames(mat),"ssc-")
pdf(paste0(PATH,"miRNA_heatmap.pdf"), width=8,height=10)
pheatmap(mat, annotation_col = anno, annotation_colors = ann_col, color=grad_col, show_colnames=FALSE, fontsize_row=8, labels_row = mirRow, main="Figure 3A",legend_labels=c("","nucleus","cytoplasm","","small","large"))
dev.off()

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#
# Figure 3B
## Boxplot of DE nuclear miR between sizes
#
nucsizediff<-rownames(nuc.sizes)[which(nuc.sizes$padj<0.1)]
mycounts <- t(log2((counts(mirD1[nucsizediff, ], normalized=TRUE, replaced=FALSE)+.5))) %>%
  merge(colData(mirD1), ., by="row.names") %>%
  gather(gene, expression, (ncol(.)-length(nucsizediff)+1):ncol(.))

mycounts$subcell <- factor(mycounts$subcell,levels=c("nucleus","cytosol"))

ggplot(mycounts, aes(size,expression)) + geom_boxplot(aes(colour=subcell), outlier.shape=NA) + geom_point(aes(colour=subcell), position=position_jitterdodge(), cex=2,alpha=1) + facet_wrap(~gene, scales="free_y") + scale_color_manual(name="subcellular\nlocalization", values=ann_col$subcell) + scale_x_discrete(limits=c("small","large"),labels=c("SGC","LGC")) + scale_y_continuous() + labs(x="", y="Expression (log normalized counts)", colour="subcellular \nlocalization", title="Figure 3B") + theme_bw(base_size=14) + theme(legend.direction="horizontal",legend.position=c(0.6,0.1))

ggsave("sm-v-lg_nuc-miRNA.pdf", width=8, height=6, units="in", path=PATH)
#
#   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-
#  / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \
# `-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `
#
# Digital PCR analysis
validation.data <- read.csv("~/Documents/Manuscripts/Nuclear microRNA/2019-12-06 Data from Eddie/ddPCR summary/Final data-Table 1.csv")
val.pcr <- validation.data %>% filter(type=='ddPCR')
val.seq <- validation.data %>% filter(type=='seq')

pcr.dat <- val.pcr %>% group_by(miRNA,GC) %>% summarize(pcr.exp = mean(ratio),pcr.min=pcr.exp-stdev,pcr.max=pcr.exp+stdev)
seq.dat <- val.seq %>% group_by(miRNA,GC) %>% summarize(seq.exp = mean(ratio),seq.min=seq.exp-stdev,seq.max=seq.exp+stdev)
left_join(pcr.dat,seq.dat)->fulldat

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Correlation
## Test the correlation between all data
full.x<-pull(fulldat,pcr.exp) # use 'pull()' to get a vector from a tbl
full.y<-pull(fulldat,seq.exp)
cor.test(full.x,full.y)

# data:  x and y
# t = 2.7064, df = 22, p-value = 0.01289
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.1207244 0.7516401
# sample estimates:
#       cor
# 0.4997818
#

cor.test(full.x,full.y,method="spearman") # ranked correlation, much stronger
## decided to go with the Pearson correlation because these are already transformed values (expression ratio) so I feel looking at rank-order is too far given the additional manipulation / distance from the original data
#     Spearman's rank correlation rho
#
# data:  full.x and full.y
# S = 844, p-value = 0.001161
# alternative hypothesis: true rho is not equal to 0
# sample estimates:
#       rho
# 0.6330435


## then check the correlation separately in nuclear and cytoplasmic miRNAs
val.pcr <- validation.data %>% filter(type=='ddPCR') %>% filter(enrich=="nuc")
val.seq <- validation.data %>% filter(type=='seq') %>% filter(enrich=="nuc")

pcr.dat <- val.pcr %>% group_by(miRNA,GC) %>% summarize(pcr.exp = mean(ratio),pcr.min=pcr.exp-stdev,pcr.max=pcr.exp+stdev)
seq.dat <- val.seq %>% group_by(miRNA,GC) %>% summarize(seq.exp = mean(ratio),seq.min=seq.exp-stdev,seq.max=seq.exp+stdev)
left_join(pcr.dat,seq.dat)->fulldat

nuc.x<-pull(fulldat,pcr.exp) # use 'pull()' to get a vector from a tbl
nuc.y<-pull(fulldat,seq.exp)

val.pcr <- validation.data %>% filter(type=='ddPCR') %>% filter(enrich=="cyto")
val.seq <- validation.data %>% filter(type=='seq') %>% filter(enrich=="cyto")

pcr.dat <- val.pcr %>% group_by(miRNA,GC) %>% summarize(pcr.exp = mean(ratio),pcr.min=pcr.exp-stdev,pcr.max=pcr.exp+stdev)
seq.dat <- val.seq %>% group_by(miRNA,GC) %>% summarize(seq.exp = mean(ratio),seq.min=seq.exp-stdev,seq.max=seq.exp+stdev)
left_join(pcr.dat,seq.dat)->fulldat

cyto.x<-pull(fulldat,pcr.exp) # use 'pull()' to get a vector from a tbl
cyto.y<-pull(fulldat,seq.exp)


cor.test(nuc.x,nuc.y)
#     Pearson's product-moment correlation
#
# data:  nuc.x and nuc.y
# t = 1.4304, df = 12, p-value = 0.1781
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  -0.1867403  0.7586154
# sample estimates:
#       cor
# 0.3816537

cor.test(cyto.x,cyto.y)
#     Pearson's product-moment correlation
#
# data:  cyto.x and cyto.y
# t = 2.4413, df = 8, p-value = 0.04048
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.04040279 0.90904874
# sample estimates:
#       cor
# 0.6534073

## Plot
# Pearson (linear) correlation
ggplot(fulldat,aes(x=pcr.exp,y=seq.exp))+ geom_smooth(method="lm",colour="darkgrey") + geom_point(aes(colour=miRNA,shape=GC),cex=4)+geom_errorbar(aes(ymin=seq.min,ymax=seq.max,colour=miRNA))+geom_errorbarh(aes(xmin=pcr.min,xmax=pcr.max,colour=miRNA))+ylab("Nuclear enrichment ratio (seq)")+xlab("Nuclear enrichment ratio (qPCR)") + theme_bw() + scale_colour_manual(name="microRNA", values=class_col[2:14]) + annotate("text",label="R = 0.500\nP = 0.013", x=7,y=12)

# Spearman's rank correlation (need a monotonic function to have geom_smooth meaningful)
ggplot(fulldat,aes(x=pcr.exp,y=seq.exp)) + geom_point(aes(colour=miRNA,shape=GC),cex=4) + geom_errorbar(aes(ymin=seq.min,ymax=seq.max,colour=miRNA)) + geom_errorbarh(aes(xmin=pcr.min,xmax=pcr.max,colour=miRNA)) + ylab("Nuclear enrichment ratio (seq)") + xlab("Nuclear enrichment ratio (qPCR)") + theme_bw() + scale_colour_manual(name="microRNA", values=class_col[2:14]) + annotate("text",label=expression(rho~"= 0.633, P = 0.001"), x=7,y=12)

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#
# Figure 4
## nuclear-enriched microRNA
nuc.dat <- validation.data %>% filter(enrich=="nuc")
nuc.dat$group <- factor(paste0(nuc.dat$GC,"-",nuc.dat$type))
nuc.dat$group <- factor(nuc.dat$group,levels=c("S-ddPCR","S-seq","L-ddPCR","L-seq")) # order factor for plotting
nuc.dat$miRNA_lab <- str_remove(nuc.dat$miRNA,"miR-")

nuc.plot <- ggplot(nuc.dat,aes(x=miRNA,y=ratio,group=group))+geom_col(aes(fill=GC,alpha=type),position='dodge',width=0.8) + scale_alpha_manual(values=c(1,0.5),name="assay",labels=c("ddPCR","seq.")) + scale_fill_manual(values=rev(paste(ann_col$size)),name="",labels=c("LGC","SGC"))+ scale_x_discrete(labels=nuc.dat$miRNA_lab)+ geom_errorbar(aes(ymin=ratio-stdev,ymax=ratio+stdev),position=position_dodge(0.8),width=0.2) + labs(x="",y="expression ratio\n(nuclear/cytoplasmic)") + theme_bw(base_size=18) + theme(axis.ticks.x = element_blank(),legend.position="bottom")+panel_border(color=paste(ann_col$subcell[1]), size=4)
#
## cytoplasmic-enriched microRNA
cyto.dat <- validation.data %>% filter(enrich=="cyto")
cyto.dat$group <- factor(paste0(cyto.dat$GC,"-",cyto.dat$type))
cyto.dat$group<-factor(cyto.dat$group,levels=c("S-ddPCR","S-seq","L-ddPCR","L-seq")) # order factor for plotting
cyto.dat$miRNA_lab <- str_remove(cyto.dat$miRNA,"miR-")
cyto.dat$GC <- factor(cyto.dat$GC,levels=c("S","L"))

cyto.plot <- ggplot(cyto.dat,aes(x=miRNA,y=ratio,group=group))+geom_col(aes(fill=GC,alpha=type),position='dodge',width=0.8) + scale_alpha_manual(values=c(1,0.5),name="",labels=c("ddPCR","seq.")) + scale_fill_manual(values=paste(ann_col$size),name="",labels=c("SGC","LGC"))+ scale_x_discrete(labels=cyto.dat$miRNA_lab)+ geom_errorbar(aes(ymin=ratio-stdev,ymax=ratio+stdev),position=position_dodge(0.8),width=0.2) + labs(x="",y="expression ratio\n(cytoplasmic/nuclear)") + theme_bw(base_size=18) + theme(axis.ticks.x = element_blank(),legend.position="bottom")+panel_border(color=paste(ann_col$subcell[2]), size=4)

# grid.arrange(nuc.plot,cyto.plot,nrow=2,top="Figure 4")

nuc.plot
ggsave("ddPCR_val_nuc.pdf", width=8, height=6, units="in", path=PATH)
cyto.plot
ggsave("ddPCR_val_cyt.pdf", width=8, height=6, units="in", path=PATH)

#   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-
#  / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \
# `-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `
#
# snoRNA analysis
## Full linear model, two factor
sizes<-results(snoD1, contrast=c("size","small","large")) # 0
locals<-results(snoD1, contrast=c("subcell","cytosol","nucleus")) # 66

## Specific groupings (comparable to Ram's analysis)
nuc.sizes<-results(snoD2, contrast=c("group","largenucleus","smallnucleus")) # 0
cyt.sizes<-results(snoD2, contrast=c("group","largecytosol","smallcytosol")) # 0
lg.local<-results(snoD2, contrast=c("group","largenucleus","largecytosol")) # 66
sm.local<-results(snoD2, contrast=c("group","smallnucleus","smallcytosol")) # 71
#
## Venn diagram
s<-sm.local
s<-s[which(s$padj<0.1),]
sno.smV<-row.names(s)

s<-lg.local
s<-s[which(s$padj<0.1),]
sno.lgV<-row.names(s)

venn.diagram(x=list("SGC"=sno.smV,"LGC" = sno.lgV),fill=ann_col$size,filename=NULL) -> venny
# grid.draw(venny)

## Significantly different between SM and LG nuclei
sm.loc <- anti_join(tbl_df(sno.smV),tbl_df(sno.lgV))
lg.loc <- anti_join(tbl_df(sno.lgV),tbl_df(sno.smV))

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#
# Figure 6
## Volcano plot of differentially expressed snoRNA
#

locres1s <- lfcShrink(snoD1, coef="subcell_nucleus_vs_cytosol", type="apeglm")
keyvals <- rep("grey",nrow(locres1s))
names(keyvals) <- rep('similar expression',nrow(locres1s))
keyvals[which(locres1s$log2FoldChange > 1)] <- ann_col$subcell[1]
names(keyvals)[which(locres1s$log2FoldChange > 1)] <- 'nucleus-enriched'
keyvals[which(locres1s$log2FoldChange < -1)] <- ann_col$subcell[2]
names(keyvals)[which(locres1s$log2FoldChange < -1)] <- 'cytoplasm-enriched'
## optionally add difference of LGC and SGC (i.e. genes that only differ in their subcellular location in one group)
keyvals[which(rownames(locres1s) %in% pull(lg.loc))] <- ann_col$size[2]
names(keyvals)[which(rownames(locres1s) %in% pull(lg.loc))] <- 'size:LGC'
keyvals[which(rownames(locres1s) %in% pull(sm.loc))] <- ann_col$size[1]
names(keyvals)[which(rownames(locres1s) %in% pull(sm.loc))] <- 'size:SGC'

## Now adding shape information based on the snoRNA class
snoRAnnotate <- read.csv(file = "piRNA_snoRNA/snoRNA-RFAM.csv", header = TRUE, stringsAsFactors = FALSE)
sno.names <- substring(snoRAnnotate$x,3)
CD <- snoRAnnotate$x[which(snoRAnnotate$class=="C/D")]
HACA <- snoRAnnotate$x[which(snoRAnnotate$class=="H/ACA")]

keyvals.shape <- rep(16,nrow(locres1s))
names(keyvals.shape) <- rep('uncategorized',nrow(locres1s))

keyvals.shape[which(rownames(locres1s) %in% CD)] <- 17
names(keyvals.shape)[which(rownames(locres1s) %in% CD)] <- 'C/D'

keyvals.shape[which(rownames(locres1s) %in% HACA)] <- 15
names(keyvals.shape)[which(rownames(locres1s) %in% HACA)] <- 'H/ACA'

EnhancedVolcano(locres1s, lab=sno.names, x='log2FoldChange', y='padj', ylab = bquote(~-Log[10]~italic(P[adj])), xlim=c(-3,3), ylim=c(0,8), FCcutoff = 1, pCutoff = 0.1, colCustom=keyvals, shapeCustom=keyvals.shape, title = "Figure 6", transcriptLabSize = 4.0, transcriptPointSize =2, drawConnectors=TRUE, colConnectors="grey50", legendPosition='right', legendLabSize = 10, legendIconSize = 3.0,subtitle="")

ggsave("volcano_snoR-localization.pdf", width=8, height=8, units="in", path=PATH)
#
#   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-.   .-.-
#  / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \ \ / / \
# `-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `-`-'   `