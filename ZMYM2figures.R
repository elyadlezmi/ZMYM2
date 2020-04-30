library("ggrepel")
library("ggplot2")
library("ggpubr")
library("gplots")
library("edgeR")
library("reshape2")
library("RColorBrewer")
library("ggsignif")
library("stringr")
library("metagene2")
library("ggfortify")
setwd('./tables/')

######################################## F1A ######################################################

df <- read.csv('ESCessentialomeCS.csv')  # logFC off all genes in CRISPR screen of Yilmaz et al. 2018
epifactors <- read.csv('EpiGenes_main.csv') 

df$isepi <- 0
df[df$Gene.Symbol %in% epifactors$HGNC_symbol,]$isepi <- 1
df$group <- 'Other'
df[df$p.value < 0.05 & df$CRISPR.score > mean(df$CRISPR.score),]$group <- 'Growth restricting'
df[df$p.value < 0.05 & df$CRISPR.score < mean(df$CRISPR.score),]$group <- 'Essential'

ggplot() + geom_point(data=df[df$isepi == 0,],aes(x=CRISPR.score,y=-log(p.value)),alpha=0.05, size=1) + 
  theme_classic() + 
  geom_point(data=df[df$isepi == 1,],aes(x=CRISPR.score,y=-log(p.value),color=group),alpha=0.9, size=1) + 
  scale_color_manual(values=c("dodgerblue4","red2","#ffefd5")) + 
  geom_text_repel(data=df[df$Gene.Symbol == 'ZMYM2',], aes(x=CRISPR.score,y=-log(p.value),
                                                               label=Gene.Symbol),size=3.5,col='black', fontface = "italic",min.segment.length = 0) +
  theme(axis.title=element_text(size=12),axis.text = element_text(size=10),legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlim(-3,2.5) + xlab('CRISPR Score') + scale_y_continuous(expand = c(0, 0),limits = c(0,37))


##### F1Apie ########
epi <- df[df$isepi == 1,]  
epi <- data.frame(
  group = c('Growth restricting', 'Essential', 'Other'),
  value = c(dim(epi[epi$group=='Growth restricting',])[1], dim(epi[epi$group=='Essential',])[1], dim(epi[epi$group=='Other',])[1]))

ggdonutchart(epi, "value", label = "group",lab.pos = 'out',fill=c("red2","dodgerblue4","#ffefd5"),lab.font = c(12, "bold", "black"))

###########################################  F1B  #######################################################

df <- read.csv('ESCessentialomeCS.csv')
epifactors <- read.csv('EpiGenes_main.csv')
complexes <- read.csv('EpiGenes_complexes.csv') # manually added KDM1A to LSD-CoREST

comp <- gsub(' ','',as.character(complexes$UniProt_ID))
comp <- strsplit(comp,',')
names(comp) <- complexes$Complex_name

df$isepi <- 0
df[df$Gene.Symbol %in% epifactors$HGNC_symbol,]$isepi <- 1
df$group <- 'Other'
df[df$p.value < 0.05 & df$CRISPR.score > mean(df$CRISPR.score),]$group <- 'Growth restricting'
df[df$p.value < 0.05 & df$CRISPR.score < mean(df$CRISPR.score),]$group <- 'Essential'

cols <- c()
for( i in names(comp)) {
  epi <- epifactors[epifactors$UniProt_ID %in% comp[[i]],]$HGNC_symbol
  epi <- df[df$Gene.Symbol %in% epi,]
  cols <- rbind(cols, c(i,dim(epi[epi$group=='Growth restricting',])[1],'Growth restricting',as.numeric(dim(epi)[1])))
  cols <- rbind(cols, c(i,dim(epi[epi$group=='Essential',])[1],'Essential',as.numeric(dim(epi)[1])))
  cols <- rbind(cols, c(i,dim(epi[epi$group=='Other',])[1],'Other',as.numeric(dim(epi)[1]))) }

epi <- data.frame(cols)
names(epi) <- c('complex','no_of_genes', 'classification','total_in_complex')
epi$total_in_complex <- as.numeric(as.character(epi$total_in_complex))
epi <- epi[epi$total_in_complex > 2 ,]  # filter for complexes with more than 3 genes
epi$no_of_genes <- as.numeric(as.character(epi$no_of_genes))

# filter out complexes who are mostly neutral
outs1 <- epi[epi$no_of_genes <= 1 & epi$classification == 'Growth restricting','complex'] 
outs2 <- epi[epi$no_of_genes <= 1 & epi$classification == 'Essential','complex']
outs <- intersect(outs1,outs2)
epi <- epi[!epi$complex %in% outs,]

epi$no_of_genes <- as.numeric(as.character(epi$no_of_genes))/epi$total_in_complex*100 # calc percentage

epi <- epi[rev(order(epi$total_in_complex)),]
epi$classification = factor(epi$classification, levels = c('Essential','Other','Growth restricting'))

sorted <- epi[epi$classification == 'Growth restricting',]
sorted <- sorted[order(sorted$no_of_genes),]
growth <- as.character(sorted[sorted$no_of_genes !=0,'complex'])

sorted <- epi[epi$classification == 'Essential',]
sorted <- sorted[!sorted$complex %in% growth,]
sorted <- sorted[order(-sorted$no_of_genes),]

growth <- c(as.character(sorted$complex),growth)

epi$complex = factor(epi$complex, levels = growth)

ggplot() + geom_bar(data=epi,aes(x=complex,y=no_of_genes,fill=classification),
                    stat="summary",color='black',size=0.2,width=0.8) + 
  scale_fill_manual(values=c("dodgerblue4","#ffefd5","red2")) + theme_classic() + 
  theme(axis.title=element_text(size=12),
        axis.text.x = element_text(size=10), axis.text.y = element_text(size=8),
        legend.position=c(0.5,0.1),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.title = element_blank(),legend.text = element_text(size = 9)) +
  xlab('Epigenetic complex') + ylab('Genes in complex (%)') + coord_flip()

###########################################  F1C  #######################################################

df <- read.csv('ESCessentialomeCS.csv')
epifactors <- read.csv('EpiGenes_main.csv')
complexes <- read.csv('EpiGenes_complexes.csv') # manually added KDM1A to LSD-CoREST

comp <- gsub(' ','',as.character(complexes$UniProt_ID))
comp <- strsplit(comp,',')
names(comp) <- complexes$Complex_name
LCH <- intersect(epifactors[epifactors$UniProt_ID %in% comp[['LSD-CoREST']],]$HGNC_symbol, epifactors[epifactors$UniProt_ID %in% comp[['BHC']],]$HGNC_symbol)
#LCH <- c('RCOR1','ZMYM2','HDAC1','HDAC2','KDM1A')

gtex <- read.csv('gtex_and_ESC.csv')

gtex <- gtex[gtex$gene_name %in% LCH,]

df <-  df[df$Gene.Symbol %in% LCH,]
df <- df[order(df$CRISPR.score,decreasing = T),]

gtex <- merge(df[,c('Gene.Symbol','CRISPR.score')],gtex,by.x = 'Gene.Symbol',by.y = 'gene_name')
gtex <- gtex[order(gtex$CRISPR.score,decreasing = T),]

gtex$X <- NULL
gtex$ensemblID <- NULL
gtex$CRISPR.score <- NULL
row.names(gtex) <- gtex$Gene.Symbol
gtex$Gene.Symbol <- NULL

#gtex <- gtex[,seq(1,11)]

gtex <- t(gtex)
gtex <- scale(gtex)
gtex <- t(gtex)

pellet <- colorRampPalette(c("dodgerblue4", "white", "red2"))
col_breaks = c(seq(-4,-1,length=10),  # for red
               seq(-0.9,0.9,length=10),           # for yellow
               seq(1,4,length=10))  

heatmap.2(gtex,trace="none",Rowv = T,col=bluered(256),density.info="none",dendrogram="both",colsep = c(31),srtCol=45)

########################################### F1D ###########################################

load('Allcounts.RData')

samples <- c("esc","eb","teratoma")
df <- df[,index %in% c('ensemblID','gene_name',samples)]
index <- index[index %in% c('ensembleID','gene_name',samples)]

# TMM normalization
y <- DGEList(df[-c(1,2)], genes=df[c(1,2)])
y <- calcNormFactors(y,method = 'TMM')
df <- data.frame(cbind(y$genes,cpm(y)))

expression <- c()
for(i in samples) {
  p <- df[df$gene_name=='ZMYM2',index == i]
  p <- as.numeric(p)
  p <- p[!p %in% boxplot(p,plot=F)$out]
  p <- cbind(p,rep(i,length(p)))
  expression <- rbind(expression,p) }

expression <- data.frame(expression=as.numeric(expression[,1]),cell_type=expression[,2])
expression <- na.omit(expression)
expression$cell_type = factor(expression$cell_type, levels = c("esc","eb","teratoma"))

#normalize (relative) to esc
expression$expression <- expression$expression / mean(expression[expression$cell_type=='esc','expression'])

ggplot(data=expression,aes(x=cell_type,y=expression,fill=cell_type)) + geom_bar(stat="summary",color='black') +
  scale_color_manual(values = c("#999999", "#A6D49F", "#C59849")) + 
  theme_classic() + scale_y_continuous(expand = c(0, 0)) +
  theme(plot.title=element_text(size=14), axis.title.y=element_text(size=12),
    axis.title.x=element_blank(),axis.text.x = element_text(size=12),legend.position="none",
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    legend.text = element_text(size = 10)) +
    geom_signif(comparisons = list(c("esc","eb"),c("esc","teratoma")),test = t.test) +
    stat_summary(fun.data = mean_se, geom = "errorbar",size=0.25, width=.2)

t.test(expression[expression$cell_type == 'esc','expression'],expression[expression$cell_type == 'eb','expression'], var.equal = TRUE)
t.test(expression[expression$cell_type == 'esc','expression'],expression[expression$cell_type == 'teratoma','expression'], var.equal = TRUE  )


########################################### F1E ###########################################

studies <- list(c("Guo_2017_Naive","Guo_2017_Primed"),c("Messmer_2019_Naive","Messmer_2019_Primed"),c("Pastor_2016_Naive","Pastor_2016_Primed"),c("Theunissen_2016_Naive","Theunissen_2016_Primed"))
expression <- c()

for(study in studies) {
  load('Allcounts.RData')
  df <- df[,index %in% c('ensembleID','gene_name',study)]
  index <- index[index %in% c('ensembleID','gene_name',study)]
  
  # TMM normalization
  y <- DGEList(df[-c(1,2)], genes=df[c(1,2)])
  keep <- rowSums(cpm(y)>1) > (length(df)-2)/length(study)
  y <- y[keep, , keep.lib.sizes=FALSE]
  y <- calcNormFactors(y,method = 'TMM')
  df <- data.frame(cbind(y$genes,cpm(y)))
  
  for(i in study) {
    p <- df[df$gene_name=='ZMYM2',index == i]
    p <- as.numeric(p)
    p <- p[!p %in% boxplot(p,plot=F)$out]
    p <- cbind(p,rep(str_sub(i,1,6),length(p)),rep(str_sub(i,-6,-1),length(p)))
    expression <- rbind(expression,p) }   }

expression <- data.frame(expression=as.numeric(expression[,1]),Study=expression[,2],Pluripotency=expression[,3])
expression <- na.omit(expression)

for(i in unique(expression$Study)) {
  expression[expression$Study==i,'expression'] <- expression[expression$Study==i,'expression'] / 
    mean(expression[expression$Study==i & expression$Pluripotency=='Primed', 'expression']) }

expression$Pluripotency = factor(expression$Pluripotency, levels = c("Primed","_Naive"))

ggplot(data=expression,aes(x=Study,y=expression,fill=Pluripotency)) +  
  geom_bar(stat="summary",position = position_dodge(width=0.9),color='black') +
  theme_classic() + scale_y_continuous(expand = c(0, 0)) +
  stat_summary(fun.data = mean_se, geom = "errorbar",size=0.25, width=.2,position = position_dodge(width=0.9))

for(i in unique(expression$Study)) {
  t <- expression[expression$Study==i,]
  print(i)
  print(t.test(t$expression~t$Pluripotency, var.equal =TRUE)) }


########################################### F1F ###########################################'

studies <- list(c("Warrier_2018_blastocyst","Warrier_2018_ESCs"),c("Yan_2013_late_blastocyst","Yan_2013_ESC_p10"))
expression <- c()

for(study in studies) {
  load('Allcounts.RData')
  n <- 1
  df <- df[,index %in% c('ensembleID','gene_name',study)]
  index <- index[index %in% c('ensembleID','gene_name',study)]
  
  # TMM normalization
  y <- DGEList(df[-c(1,2)], genes=df[c(1,2)])
  keep <- rowSums(cpm(y)>1) > (length(df)-2)/length(study)
  y <- y[keep, , keep.lib.sizes=FALSE]
  y <- calcNormFactors(y,method = 'TMM')
  df <- data.frame(cbind(y$genes,cpm(y)))
  
  for(i in study) {
    p <- df[df$gene_name=='ZMYM2',index == i]
    p <- as.numeric(p)
    p <- p[!p %in% boxplot(p,plot=F)$out]
    p <- cbind(p,rep(str_sub(i,1,6),length(p)),rep(n,length(p)))
    expression <- rbind(expression,p)
    n <- n+1}
}

expression <- data.frame(expression=as.numeric(expression[,1]),Study=expression[,2],Pluripotency=expression[,3])
expression <- na.omit(expression)

for(i in unique(expression$Study)) {
  expression[expression$Study==i,'expression'] <- expression[expression$Study==i,'expression'] / 
    mean(expression[expression$Study==i & expression$Pluripotency==2, 'expression']) }

ggplot(data=expression,aes(x=Study,y=expression,fill=Pluripotency)) +  
  geom_bar(stat="summary",position = position_dodge(width=0.9),color='black') +
  theme_classic() + scale_y_continuous(expand = c(0, 0)) +
  #  geom_signif(test = t.test,comparisons = c("naive",'primed')) +
  stat_summary(fun.data = mean_se, geom = "errorbar",size=0.25, width=.2,position = position_dodge(width=0.9))

for(i in unique(expression$Study)) {
  t <- expression[expression$Study==i,]
  print(i)
  print(t.test(t$expression~t$Pluripotency, var.equal =TRUE)) }

########################################### F2C ###########################################

df <- read.csv("./tables/annexin.csv")
df$condition <- 'CTRL'
df[df$sample %in% c('z6','z15','z22'),'condition'] <- 'ZMYM2-null'
df$experiment <- 'apoptosis'
df <- df[df$Name=="Q2",]
df[,c('Name',"X.Cells","sample")] <- NULL
all <- df

df <- read.csv('./tables/CellsCycle16.csv')
df$condition <- 'CTRL'
df[df$Name %in% c('z6','z15','z22'),'condition'] <- 'ZMYM2-null'
df$experiment <- 'proliferation'
df[,c('Name',"X.Cells")] <- NULL
all <- rbind(all,df)

ggplot(data=all,aes(x=experiment,y=Statistic,fill=condition)) + geom_bar(stat="summary",width = 0.6,color='black',position=position_dodge(0.6)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", size=0.25, width=.1,position=position_dodge(0.6)) + theme_classic() + scale_fill_manual(values=c("#3172C0","#FFAF30")) +
  scale_y_continuous(expand = c(0, 0))

########################################### F2D  ###########################################

DE <- read.csv('koDE.csv')[-1]
DE <- DE[order(DE$FDR, decreasing = F),]

df <- read.csv('pluripotency_markers.csv')
markers <- df[df$Marker == 'core_pluripotency',]$gene
naive_markers <- df[df$Marker == 'naive_pluripotency',]$gene
primed_markers <- df[df$Marker == 'primed_pluripotency',]$gene

naive_genes <- read.csv('NaivePrimed6studiesDE.csv')
naive_genes <- naive_genes[order(naive_genes$FDR,-abs(naive_genes$logFC)),]
naive_genes <- naive_genes[naive_genes$FDR < 0.05,]

primed_genes <- naive_genes[naive_genes$logFC < 0,]$gene_name
naive_genes <-  naive_genes[naive_genes$logFC > 0,]$gene_name

write_names <- c(as.character(markers),as.character(naive_markers),as.character(primed_markers))

### all genes
ggplot(data=DE) + theme_classic() +
  geom_point(data=DE,aes(x=logFC, y=-log(FDR)), alpha=0.5, size=2,colour = "#A0A0A0") +
  geom_point(data=DE[DE$gene_name %in% naive_markers & DE$FDR < 0.05,],aes(x=logFC, y=-log(FDR)), alpha=1, size=3.5,colour = "#990000") +
  geom_point(data=DE[DE$gene_name %in% primed_markers & DE$FDR < 0.05,],aes(x=logFC, y=-log(FDR)), alpha=1, size=3.5,colour = "#004C99") +
  geom_point(data=DE[DE$gene_name %in% markers & DE$FDR < 0.05,],aes(x=logFC, y=-log(FDR)), alpha=1, size=3.5,colour = "black") +
  xlab("logFC ZMYM2 mutant vs. WT") + ylab("-logFDR") + labs(color = "Gene expression profile") +
  theme(plot.title=element_text(size=14), axis.title.y=element_text(size=12),
        axis.text.x = element_text(size=12),legend.position=c(0.3,0.85),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.text = element_text(size = 12),legend.title =element_text(size = 12) ) + 
  geom_hline(yintercept = -log(0.05), col='red') + scale_y_continuous(expand = c(0, 0)) + 
  geom_text_repel(data=DE[DE$gene_name %in% write_names & DE$FDR < 0.05,], aes(x=logFC, y=-log(FDR),label=gene_name),size=3,fontface='italic',col='black')


################  F2Dpie  ############

library("ggpubr")
df <- data.frame(
  group = c('Upregulated', 'Downregulated'),
  value = c(dim(DE[DE$logFC > 0 & DE$FDR <= 0.05,])[1], dim(DE[DE$logFC < 0 & DE$FDR <= 0.05,])[1]))

ggdonutchart(df, "value", label = "group",lab.pos = 'out',fill=c("red2","dodgerblue4"),lab.font = c(12, "bold", "black"))


############################################ F2E  ###########################################

DE <- read.csv('koDE.csv')[-1]
DE <- DE[order(DE$FDR, decreasing = F),]

naive_genes <- read.csv('NaivePrimed6studiesDE.csv')
naive_genes <- naive_genes[order(naive_genes$FDR,-abs(naive_genes$logFC)),]
naive_genes <- naive_genes[naive_genes$FDR < 0.05,]

primed_genes <- naive_genes[naive_genes$logFC < 0,]$gene_name
naive_genes <-  naive_genes[naive_genes$logFC > 0,]$gene_name
ggplot(data=DE) + theme_classic() +
  geom_point(data=DE,aes(x=logFC, y=-log(FDR)), alpha=0.5, size=2,colour = "#A0A0A0") +
  geom_point(data=DE[DE$gene_name %in% primed_genes & DE$FDR < 0.05,],aes(x=logFC, y=-log(FDR)), alpha=0.4, size=3.5,colour = "#004C99",fill=NA,stroke=1,pch=21) +
  geom_point(data=DE[DE$gene_name %in% naive_genes & DE$FDR < 0.05,],aes(x=logFC, y=-log(FDR)), alpha=0.4, size=3.5,colour = "#990000",fill=NA,stroke=1.5,pch=21) +
  xlab("logFC ZMYM2 mutant vs. WT") + ylab("-logFDR") + labs(color = "Gene expression profile") +
  theme(plot.title=element_text(size=14), axis.title.y=element_text(size=12),
        axis.text.x = element_text(size=12),legend.position=c(0.3,0.85),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.text = element_text(size = 12),legend.title =element_text(size = 12) ) + 
  geom_hline(yintercept = -log(0.05), col='red') + scale_y_continuous(expand = c(0, 0))
  
#### chi-square

naive_genes <- read.csv('NaivePrimed6studiesDE.csv')
up_in_naive <- naive_genes[naive_genes$FDR < 0.05 & naive_genes$logFC > 0,]$ensemblID
down_in_naive <- naive_genes[naive_genes$FDR < 0.05 & naive_genes$logFC < 0,]$ensemblID
ns_in_naive <- naive_genes[naive_genes$FDR > 0.05,]$ensemblID

DE <- read.csv('koDE.csv')[-1]
up_in_de <- DE[DE$FDR < 0.05 & DE$logFC > 0,]$ensemblID
down_in_de <- DE[DE$FDR < 0.05 & DE$logFC < 0,]$ensemblID
ns_in_de <- DE[DE$FDR > 0.05,]$ensemblID

col1 <- c(length(intersect(up_in_naive,up_in_de)),length(intersect(up_in_naive,down_in_de)),length(intersect(up_in_naive,ns_in_de)))
col2 <- c(length(intersect(down_in_naive,up_in_de)),length(intersect(down_in_naive,down_in_de)),length(intersect(down_in_naive,ns_in_de)))
col3 <- c(length(intersect(ns_in_naive,up_in_de)),length(intersect(ns_in_naive,down_in_de)),length(intersect(ns_in_naive,ns_in_de)))

chi <- cbind(col1,col2,col3)
chisq.test(chi)


############################################ F2F  ###########################################

df <- read.csv('koDE.csv')
df <- df[order(df$FDR),]
topde <- df$ensemblID[1:100]

df <-  read.csv('normed_expression_medians.csv',row.names = 1)
samples <- c("ensemblID","gene_name","Guo_2017_Naive","Guo_2017_Primed","Messmer_2019_Naive","Messmer_2019_Primed","Pastor_2016_Naive","Pastor_2016_Primed" ,
             "Theunissen_2016_Naive","Theunissen_2016_Primed","Warrier_2018_blastocyst","Warrier_2018_ESCs","Yan_2013_ESC_p10","Yan_2013_late_blastocyst","control","ZMYM2.mutant")
df <- df[samples]
df <- df[df$ensemblID %in% topde,]

row.names(df) <- df$gene_name
df$ensemblID <- NULL
df$gene_name <- NULL
df <- df[complete.cases(df),]

df <- t(df)
df <- scale(df)
df <- t(df)
heatmap.2(as.matrix(df),trace="none",Rowv = T,col=bluered(256),density.info="none",dendrogram="both",srtCol=20)


#################################### F3A #######################################

library(DiffBind)

ac <- dba(sampleSheet="samplesheetH3Ac.csv")
plot(ac)

#################################### F3B #######################################

igg <- "H3Ac"
feature <- "./metaPromoter/metaPromoter.csv"
meta <- data.frame(matrix(ncol = 3,nrow = 0))
colnames(meta) <- c("mean_CPM","Distance_from_TSS","sample")

for(i in c("emp1","emp2","emp3","z15","z22")) {
  
  df <- read.csv(paste(i,igg,feature,sep = "."),row.names = 1)
  df[,c("X.4950","X.5000","X5000","X4950")]  <- NULL
  p <- colMeans(df[,-c(1:8)])
  p <- stack(p)
  names(p) <- c("mean_CPM","Distance_from_TSS")
  p[,2] <- gsub('X','',p[,2])
  p[,2] <- gsub('\\.','-',p[,2])
  p[,2] <- as.numeric(p[,2])
  p$sample <- i 
  
  meta <- rbind(meta,p)  }

meta$condition <- 'Control'
meta[meta$sample %in% c("z6","z15","z22"),"condition"] <- "ZMYM2-null"

a <- aggregate(meta[, c(1,2)], list(meta$Distance_from_TSS,meta$condition),mean)

ggplot() + geom_path(data = a,aes(x=Distance_from_TSS,y=mean_CPM,color=Group.2)) + theme_classic() + 
  scale_color_manual(values=c("#3172C0","#FFAF30"))

#################################### F3C #######################################

DE <- read.csv('differential_H3Ac.csv',row.names = 1)
names(DE)[c(1,12)] <- c('gene_name','FDR')
DE <- DE[order(DE$FDR, decreasing = F),]

### all genes
ggplot(data=DE) + theme_classic() +
  geom_point(data=DE[!(DE$upreg | DE$downreg),],aes(x=Fold, y=-log(FDR)), alpha=0.5, size=2,colour = "#A0A0A0") +
  geom_point(data=DE[DE$upreg,],aes(x=Fold, y=-log(FDR)), alpha=0.5, size=2,colour = "red2") +
  geom_point(data=DE[DE$downreg,],aes(x=Fold, y=-log(FDR)), alpha=0.5, size=2,colour = "dodgerblue4") +
  xlab("logFC ZMYM2 mutant vs. WT") + ylab("-logFDR") + labs(color = "Gene expression profile") +
  theme_classic() + 
  geom_hline(yintercept = -log(0.05), col='red') + scale_y_continuous(expand = c(0, 0)) 

### pie
df <- read.csv('differential_H3Ac.csv',row.names = 1)
library("ggpubr")
DE <- df
df <- data.frame(
  group = c('Upregulated', 'Downregulated'),
  value = c(dim(DE[DE$Fold > 0.5 & DE$FDR.x <= 0.05,])[1], dim(DE[DE$Fold < -0.5 & DE$FDR.x <= 0.05,])[1]) )

ggdonutchart(df, "value", label = "group",lab.pos = 'out',fill=c("red2","dodgerblue4"),lab.font = c(12, "bold", "black"))

### GSEA

df <- read.csv('GSEA.csv', col.names = c('term','FDR','group'))
df <- df[df$group=='DOWN',]
df$term <- factor(df$term, levels = c('GO_EMBRYONIC_MORPHOGENESIS','GO_CENTRAL_NERVOUS_SYSTEM_DEVELOPMENT',
                                      'GO_DIGESTIVE_TRACT_MORPHOGENESIS','GO_EMBRYONIC_ORGAN_DEVELOPMENT'))

ggplot(df[df$group=='DOWN',]) + geom_bar(aes(x=term,y=-log(FDR)),stat = 'identity',fill="dodgerblue4") + 
  theme_classic() + coord_flip()

df <- read.csv('GSEA.csv',header = F,col.names = c('term','FDR','group'))
df <- df[df$group=='UP',]
df$term <- factor(df$term, levels = c('GO_REPRODUCTION','GO_GAMETE_GENERATION','GO_GERM_CELL_DEVELOPMENT',
                                      'GO_CELLULAR_PROCESS_INVOLVED_IN_REPRODUCTION_IN_MULTICELLULAR_ORGANISM'))

ggplot(df[df$group=='UP',]) +
  geom_bar(aes(x=term,y=-log(FDR)),stat = 'identity',fill="red2") +
  theme_classic() + coord_flip()

#################################### F3D  #######################################

DE <- read.csv('differential_H3Ac.csv',row.names = 1)
names(DE)[c(1,12)] <- c('gene_name','FDR')
DE <- DE[order(DE$FDR, decreasing = F),]

df <- read.csv('./tables/pluripotency_markers.csv')
markers <- df[df$Marker == 'core_pluripotency',]$gene
naive_markers <- df[df$Marker == 'naive_pluripotency',]$gene
primed_markers <- df[df$Marker == 'primed_pluripotency',]$gene

write_names <- c(as.character(markers),as.character(naive_markers),as.character(primed_markers))

DE <- DE[DE$FDR <= 0.05,]

### all genes
ggplot(data=DE) + theme_classic() +
  geom_smooth(data = DE,aes(x=Fold, y=logFC),method='lm', color='grey') +
  geom_point(data=DE,aes(x=Fold, y=logFC), alpha=0.4, size=2,colour = "#A0A0A0") +
  geom_point(data=DE[DE$gene_name %in% naive_markers & DE$FDR < 0.05,],aes(x=Fold, y=logFC), alpha=1, size=3.5,colour = "#990000") +
  geom_point(data=DE[DE$gene_name %in% primed_markers & DE$FDR < 0.05,],aes(x=Fold, y=logFC), alpha=1, size=3.5,colour = "#004C99") +
  geom_point(data=DE[DE$gene_name %in% markers & DE$FDR < 0.05,],aes(x=Fold, y=logFC), alpha=1, size=3.5,colour = "black") +
  xlab("logFC ZMYM2 mutant vs. WT") + ylab("-logFDR") + labs(color = "Gene expression profile") +
  theme_classic() + 
  geom_text_repel(data=DE[DE$gene_name %in% write_names & DE$FDR < 0.05,], aes(x=Fold, y=logFC,label=gene_name),size=3,fontface='italic',col='black')


DE <- DE[!is.na(DE$logFC),]
cor.test(DE$Fold,DE$logFC, method='pearson')

##################################### F4A ##########################################

df <- read.csv('4diff.csv')
df <- melt(df)
df$zmym <- 'Control'
df[df$variable %in% c('Z6','Z15','Z22'),'zmym'] <- 'ZMYM2-null'
df$Condition = factor(df$Condition, levels = c('mTeSR','mTeSR+RA','EB','EB+RA'))

t <- df[df$Gene == 'POU5F1',]

ggplot(data=t,aes(x=Condition,y=value,fill=zmym)) +  
  geom_bar(stat="summary",position = position_dodge(width=0.9),color='black') +
  theme_classic() + scale_y_continuous(expand = c(0, 0)) + scale_fill_manual(values=c("#66B2FF","#FFB266")) +
  stat_summary(fun.data = mean_se, geom = "errorbar",size=0.25, width=.2,position = position_dodge(width=0.9))

for(i in unique(t$Condition)[-1]) {
  a <- t[t$Condition == i,]
  print(i)
  print(t.test(a$value~a$zmym, var.equal=T))  }

t <- df[df$Gene == 'NANOG',]

ggplot(data=t,aes(x=Condition,y=value,fill=zmym)) +  
  geom_bar(stat="summary",position = position_dodge(width=0.9),color='black') +
  theme_classic() + scale_y_continuous(expand = c(0, 0)) + scale_fill_manual(values=c("#66B2FF","#FFB266")) +
  stat_summary(fun.data = mean_se, geom = "errorbar",size=0.25, width=.2,position = position_dodge(width=0.9))

for(i in unique(t$Condition)[-1]) {
  a <- t[t$Condition == i,]
  print(i)
  print(t.test(a$value~a$zmym, var.equal=T))  }


##################################### F4B ##########################################

df <- read.csv('nanogRAdiff.csv')
df$zmym <- 'Control'
df[df$line %in% c('z6','z15','z22'),'zmym'] <- 'ZMYM2-null'
df$day <- as.factor(df$day)

ggplot(data=df,aes(x=day,y=expression,fill=zmym)) + 
  geom_bar(stat="summary",position = position_dodge(width=0.9),color='black') +
  theme_classic() + scale_y_continuous(expand = c(0, 0)) + scale_fill_manual(values=c("#66B2FF","#FFB266")) +
  stat_summary(fun.data = mean_se, geom = "errorbar",size=0.25, width=.2,position = position_dodge(width=0.9))

for(i in unique(df$day)[-1]) {
  a <- df[df$day == i,]
  print(i)
  print(t.test(a$expression~a$zmym, var.equal=T))  }

df <- read.csv('oct4RAdiff.csv')
df$zmym <- 'Control'
df[df$line %in% c('z6','z15','z22'),'zmym'] <- 'ZMYM2-null'
df$day <- as.factor(df$day)

ggplot(data=df,aes(x=day,y=expression,fill=zmym)) + 
  geom_bar(stat="summary",position = position_dodge(width=0.9),color='black') +
  theme_classic() + scale_y_continuous(expand = c(0, 0)) + scale_fill_manual(values=c("#66B2FF","#FFB266")) +
  stat_summary(fun.data = mean_se, geom = "errorbar",size=0.25, width=.2,position = position_dodge(width=0.9))

for(i in unique(df$day)[-1]) {
  a <- df[df$day == i,]
  print(i)
  print(t.test(a$expression~a$zmym, var.equal=T))  }


#################################### F4C ###########################################

df <- read.csv('Teratoma_size.csv')  

ggplot(data=df,aes(x=condition,y=Size)) + theme_classic() + stat_boxplot(geom = "errorbar", width = 0.5) + geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color=condition),position=position_jitter(0.2),size=5,alpha=0.8) +
  scale_color_manual(values=c("#66B2FF","#FFB266")) +  expand_limits(x = 0, y = 0) + ylab("Teratoma weight (g)") + 
  theme(legend.position=c(0.8,0.85), axis.title.y = element_text(size=14),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.text.x = element_text(size=14),axis.text.y = element_text(size=12), 
        axis.ticks.x = element_blank())

t.test(df$Size~df$condition, var.equal = T)


#################################### F4E ###########################################

df <- read.csv('pluripotency_markers.csv')
pluri <- as.character(df[df$Marker == 'core_pluripotency',]$gene)

DE <- read.csv('zmym_teratomaDE.csv')

ggplot(data=DE) + theme_classic() +
  geom_point(data=DE,aes(x=logFC, y=-log(PValue)), alpha=0.5, size=2,colour = "#A0A0A0") +
  geom_point(data=DE[DE$gene_name %in% pluri  & DE$PValue < 0.05,],aes(x=logFC, y=-log(PValue)), alpha=1, size=3.5,colour = 'black') +
  xlab("logFC ZMYM2-null vs. Control") + ylab("-log(P-Value)") +
  geom_hline(yintercept = -log(0.05), col='red') + scale_y_continuous(expand = c(0, 0)) + 
  geom_text_repel(data=DE[DE$gene_name %in% pluri & DE$PValue < 0.05,], aes(x=logFC, y=-log(PValue),label=gene_name),size=3,fontface='italic',col='black')

#################################### FS1D ###########################################

# eKAR.R

#################################### FS1E ###########################################

df <- read.csv('zmym2_koTPM.csv', row.names = 2)
df <- df[,-c(1,3)]
df <- melt(df)
df$zmym <- 'Control'
df[df$variable %in% c('z6','z15','z22'),'zmym'] <- 'ZMYM2-null'

df <- df[df$gene_name %in% c("POU5F1","NANOG","LIN28A", "DNMT3B"),]

ggplot(data=df,aes(x=gene_name,y=value,fill=zmym)) + 
  geom_bar(stat="summary",position = position_dodge(width=0.9),color='black') +
  theme_classic() + scale_y_continuous(expand = c(0, 0)) + scale_fill_manual(values=c("#66B2FF","#FFB266")) +
  stat_summary(fun.data = mean_se, geom = "errorbar",size=0.25, width=.2,position = position_dodge(width=0.9))

#################################### FS1F ###########################################

y <- read.csv('zmym2_koTPM.csv', row.names = 2)
y <- y[-c(1,2,3)]
g <- read.csv('BJandMoretpm.csv', row.names = 2)
#g <- g[-c(1,2,3)]
g <- g[-c(1,2,3,11,12,13,14)]

y <- cbind(y,g)

y <- y[rowSums(y) != 0,]
y <- t(y)
y <- prcomp(y,scale. = T, center = T)

autoplot(y, label.repel=T) + theme_classic()


#################################### FS2A ###########################################

igg <- "H3K4me1"
feature <- "metaPromoter.csv"
meta <- data.frame(matrix(ncol = 3,nrow = 0))
colnames(meta) <- c("mean_CPM","Distance_from_TSS","sample")

for(i in c("emp1","emp2","emp3","z15","z22")) {
  
  df <- read.csv(paste(i,igg,feature,sep = "."),row.names = 1)
  df[,c("X.4950","X.5000","X5000","X4950")]  <- NULL
  p <- colMeans(df[,-c(1:8)])
  p <- stack(p)
  names(p) <- c("mean_CPM","Distance_from_TSS")
  p[,2] <- gsub('X','',p[,2])
  p[,2] <- gsub('\\.','-',p[,2])
  p[,2] <- as.numeric(p[,2])
  p$sample <- i 
  
  meta <- rbind(meta,p)  }

meta$condition <- 'Control'
meta[meta$sample %in% c("z6","z15","z22"),"condition"] <- "ZMYM2-null"

a <- aggregate(meta[, c(1,2)], list(meta$Distance_from_TSS,meta$condition),mean)

ggplot() + geom_path(data = a,aes(x=Distance_from_TSS,y=mean_CPM,color=Group.2)) + theme_classic() + 
  scale_color_manual(values=c("#3172C0","#FFAF30"))

#################################### FS2B ###########################################

df <- read.csv('zmym2_koTPM.csv')
df$wt <- rowMeans(df[,c("emp1","emp2","emp3")])

groups <- list("TPM<1"=as.character(df[df$wt > 0 & df$wt <= 1 ,"gene_name"]),"TPM<5"=as.character(df[df$wt > 1 & df$wt <= 5 ,"gene_name"]),
               "TPM<25"=as.character(df[df$wt > 5 & df$wt <= 25 ,"gene_name"]),"TPM<100"=as.character(df[df$wt > 25 & df$wt <= 100 ,"gene_name"]),
               "TPM>100"=as.character(df[df$wt > 100,"gene_name"]))

meta <- data.frame(matrix(ncol = 3,nrow = 0))
colnames(meta) <- c("mean_CPM","Distance_from_TSS","expression")

for(g in names(groups)) {
  df <- rbind(read.csv("emp1.H3Ac.metaPromoter.csv",row.names = 1),read.csv("emp2.H3Ac.metaPromoter.csv",row.names = 1),read.csv("emp3.H3Ac.metaPromoter.csv",row.names = 1))
  df <- df[df$gene_name %in% groups[g][[1]],]
  df[,c("X.4950","X.5000","X5000","X4950")]  <- NULL
  p <- colMeans(df[,-c(1:8)])
  p <- stack(p)
  names(p) <- c("mean_CPM","Distance_from_TSS")
  p[,2] <- gsub('X','',p[,2])
  p[,2] <- gsub('\\.','-',p[,2])
  p[,2] <- as.numeric(p[,2])
  p$expression <- g
  
  meta <- rbind(meta,p)} 

meta$expression <- factor(meta$expression, levels=c("TPM<1","TPM<5","TPM<25","TPM<100", "TPM>100"))

ggplot() + geom_path(data = meta,aes(x=Distance_from_TSS,y=mean_CPM,color=expression)) + theme_classic() +
  scale_color_manual(values=c("#0F5EBD","#5389C9","#3172C0","#094995","#063A78")) + scale_y_continuous(expand = c(0,0),limits = c(0,0.9) ) +
  scale_x_continuous(expand = c(0,0),limits = c(-5000,5000) )


#################################### FS2C ###########################################

library(ggplot2)

igg <- "H3Ac"
feature <- "metaPromoter.csv"

df <- read.csv('zmym2_koTPM.csv')
df$ko <- rowMeans(df[,c("z6","z22","z15")])

groups <- list("TPM<1"=as.character(df[df$ko > 0 & df$ko <= 1 ,"gene_name"]),"TPM<5"=as.character(df[df$ko > 1 & df$ko <= 5 ,"gene_name"]),
               "TPM<25"=as.character(df[df$ko > 5 & df$ko <= 25 ,"gene_name"]),"TPM<100"=as.character(df[df$ko > 25 & df$ko <= 100 ,"gene_name"]),
               "TPM>100"=as.character(df[df$ko > 100,"gene_name"]))

meta <- data.frame(matrix(ncol = 3,nrow = 0))
colnames(meta) <- c("mean_CPM","Distance_from_TSS","expression")

for(g in names(groups)) {
  df <- rbind(read.csv("z6.H3Ac.metaPromoter.csv",row.names = 1),read.csv("z15.H3Ac.metaPromoter.csv",row.names = 1),read.csv("z22.H3Ac.metaPromoter.csv",row.names = 1))
  df <- df[df$gene_name %in% groups[g][[1]],]
  df[,c("X.4950","X.5000","X5000","X4950")]  <- NULL
  p <- colMeans(df[,-c(1:8)])
  p <- stack(p)
  names(p) <- c("mean_CPM","Distance_from_TSS")
  p[,2] <- gsub('X','',p[,2])
  p[,2] <- gsub('\\.','-',p[,2])
  p[,2] <- as.numeric(p[,2])
  p$expression <- g
  
  meta <- rbind(meta,p)} 

meta$expression <- factor(meta$expression, levels=c("TPM<1","TPM<5","TPM<25","TPM<100", "TPM>100"))


ggplot() + geom_path(data = meta,aes(x=Distance_from_TSS,y=mean_CPM,color=expression)) + theme_classic() +
  scale_color_manual(values=c("#FF9C00","#FFC05B","#FFAF30","#E48C00","#B77000")) + scale_y_continuous(expand = c(0,0),limits = c(0,0.9) ) +
  scale_x_continuous(expand = c(0,0),limits = c(-5000,5000) )
