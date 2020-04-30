
############################## filter GTEx table ###########################################

gtex <- read.csv(gzfile('GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_median_tpm.gct.gz','rt'),sep = '\t',skip = 2)

names(gtex)

gtex$Adipose <- apply(gtex[,3:4],1,mean,na.rm=TRUE)
gtex$Artery <- apply(gtex[,6:8],1,mean,na.rm=TRUE)
gtex$Brain <- apply(gtex[,10:22],1,mean,na.rm=TRUE)
gtex$Cervix <- apply(gtex[,26:27],1,mean,na.rm=TRUE)
gtex$Colon <- apply(gtex[,28:29],1,mean,na.rm=TRUE)
gtex$Esophagus <- apply(gtex[,30:32],1,mean,na.rm=TRUE)
gtex$Heart <- apply(gtex[,34:35],1,mean,na.rm=TRUE)
gtec$Skin <- apply(gtex[,46:47],1,mean,na.rm=TRUE) 

gtex[,c(3,4,6,7,8,10:22,26,27,28,29,30,31,32,34,35,46,47)] <- NULL

names(gtex)[1] <- 'ensemblID'
names(gtex)[2] <- 'gene_name'
names(gtex) <- gsub("\\.\\.\\."," ",names(gtex))
names(gtex) <- gsub("\\."," ",names(gtex))

write.csv(gtex,'GTEx_2016-01-15_v7_SimilarCombined.csv')


##### add ESCs to GTEx

gtex <- read.csv('GTEx_2016-01-15_v7_SimilarCombined.csv',row.names = 1)
gtex$ensemblID <- gsub("\\..*","",gtex$ensemblID)
gtex$gene_name <- NULL 

df <- read.csv('./expression/escTPM.csv',row.names = 1)
df$ESCs <- apply(df[,-c(1,2,3)],1,median)
df$ensemblID <- gsub("\\..*","",df$ensemblID)
df <- df[,c('ensemblID','gene_name','ESCs')]

gtex <- merge(df,gtex,by='ensemblID')
write.csv(gtex,'gtex_and_ESC.csv')

############################### median normalized expression of ALL samples ################################

library("edgeR")
cell_types <- system('ls /media/elyad/HDD/EXPRESSION/', intern = T)

df <-  read.csv(paste('/media/elyad/HDD/EXPRESSION/',cell_types[1],'/counts.csv',sep=''))[c(2,3)]             
index <- c('ensemblID','gene_name')

for(i in cell_types) {
  y <- read.csv(paste('/media/elyad/HDD/EXPRESSION/',i,'/counts.csv',sep=''))[-c(1,3)]
  y <- y[colSums(!is.na(y)) > 0]
  index <- c(index,rep(i,length(y)-1))
  df <- merge(df,y,by='ensemblID')  }

y <- read.csv('zmym2ReadCounts.csv')[-c(1,3,4)]
y$emp2 <- NULL #outlier
index <- c(index,c('control','control','ZMYM2-mutant','ZMYM2-mutant','ZMYM2-mutant'))
df <- merge(df,y,by='ensemblID')  

y <- read.csv('/home/elyad/reference_genome/data/escCounts.csv')[-c(1,3,4)]
index <- c(index,rep('esc',length(y)-1))
df <- merge(df,y,by='ensemblID')  

index <- index[c(T,T,colSums(df[-c(1,2)]) != 0)]
df <- df[,c(T,T,colSums(df[-c(1,2)]) != 0)]

save(df,index,file = 'Allcounts.RData') # Save all counts of all samples

y <- DGEList(df[-c(1,2)], genes=df[c(1,2)])

# TMM normalization
y <- calcNormFactors(y,method = 'TMM')

y <- data.frame(cbind(y$genes,cpm(y)))

df <- y[c(1,2)]
for(i in c(cell_types,'control','ZMYM2-mutant','esc')) {
  df[i] <- apply(y[,index == i],1,median) }

write.csv(df,'normed_expression_medians.csv') # median expression of each condition of all experiments


######################## differential expression ZMYM2 KO VS WT ESCs #######################

library("edgeR")

group <- c(1,1,1,2,2,2)
counts <- read.csv('./tables/zmym2ReadCounts.csv')[-c(1,4)]

cpm <- apply(counts[-c(1,2)],2, function(x) (x/sum(x))*1000000)
keep <- rowSums(cpm > 1) > 1
# DE analysis 
y <- DGEList(counts[c(-1,-2)], group, genes=counts[c(1,2)])
# remove unexpressed genes

# TMM normalization
y <- calcNormFactors(y,method = 'TMM')

# MDS plot 
plotMDS(y)                      

## start over after removing emp2 (outlier)
group <- c(1,1,2,2,2)
counts <- read.csv('zmym2ReadCounts.csv')[-c(1,4)]
counts$emp2 <- NULL
# DE analysis 
y <- DGEList(counts[c(-1,-2)], group, genes=counts[c(1,2)])
# remove unexpressed genes
#keep <- rowSums(cpm(y)>1) > 1
y <- y[keep, , keep.lib.sizes=FALSE]
# TMM normalization
y <- calcNormFactors(y,method = 'TMM')
# MDS plot 
plotMDS(y)                      

# fit model and calculate DE
design <- model.matrix(~group)
y <- estimateDisp(y, design)
y <- glmFit(y, design)
y <- glmLRT(y)
# list of ALL genes with logFC and FDR
DE <- data.frame(topTags(y,n = Inf,sort.by='p.value'))  

write.csv(DE,'koDE.csv')

############################ differential expression naive VS primed ############################

library("edgeR")
load('Allcounts.RData')

primed <- c("Guo_2017_Primed","Messmer_2019_Primed", "Pastor_2016_Primed","Theunissen_2016_Primed","Warrier_2018_ESCs","Yan_2013_ESC_p10")
naive <- c("Guo_2017_Naive","Messmer_2019_Naive","Pastor_2016_Naive","Theunissen_2016_Naive","Warrier_2018_blastocyst","Yan_2013_late_blastocyst")

for( i in c(1:6)) {

  p <- df[,index %in% c('ensemblID','gene_name',primed[i])]
  primed_count <- length(names(p)) - 2
  n <- df[,index %in% c(naive[i])]
  p <- cbind(p,n)
  naive_count <- length(names(p)) - 2 - primed_count
  group <- c(rep(1,primed_count),rep(2,naive_count))
  
  # DE analysis 
  y <- DGEList(p[c(-1,-2)], group, genes=p[c(1,2)])
  # remove unexpressed genes
  keep <- rowSums(cpm(y)>1) > 1
  y <- y[keep, , keep.lib.sizes=FALSE]
  # TMM normalization
  y <- calcNormFactors(y)
  # MDS plot ### figure 4 ###
  #plotMDS(y)                      
  
  # fit model and calculate DE
  design <- model.matrix(~group)
  y <- estimateDisp(y, design)
  y <- glmFit(y, design)
  y <- glmLRT(y)
  # list of ALL genes with logFC and FDR
  DE <- data.frame(topTags(y,n = Inf,sort.by='p.value'))  
  
  write.csv(DE,paste(naive[i],'VsPrimedDE.csv',sep='')) }

####### combine all 6 studies
df <- read.csv("Guo_2017_NaiveVsPrimedDE.csv",row.names = 1)[c(1,2)]

for(i in c("Guo_2017_NaiveVsPrimedDE.csv","Messmer_2019_NaiveVsPrimedDE.csv","Pastor_2016_NaiveVsPrimedDE.csv",
           "Theunissen_2016_NaiveVsPrimedDE.csv","Warrier_2018_blastocystVsPrimedDE.csv","Yan_2013_late_blastocystVsPrimedDE.csv")) {
  temp <- read.csv(i)
  temp[paste(i,'logFC',sep="")] <- temp$logFC
  temp[paste(i,'FDR',sep="")] <- temp$FDR
  df <- merge(df,temp[,c('ensemblID',paste(i,'logFC',sep=""),paste(i,'FDR',sep=""))],by='ensemblID')}

# sig in at least 5 studies
sig_index <- rowSums(df[c(4,6,8,10,12,14)] < 0.05) >= 5
df$FDR <- 1
df[sig_index,]$FDR <- apply(df[sig_index, c(4,6,8,10,12,14)],1,median)

df$logFC <- apply(df[c(3,5,7,9,11,13)],1,median)


write.csv(df[c('ensemblID','gene_name','logFC','FDR')],'NaivePrimed6studiesDE.csv') 


#################################### metagene analysis from bedgraph ####################################

library(dplyr)

findcpm <- function(coord,bedgraph) {
  if(length(filter(bedgraph,start<coord & end>coord)$cpm) != 1) {return(0)}
  else {return(filter(bedgraph,start<coord & end>coord)$cpm) } }

igg <- 'H3Ac'
for(i in c("emp1","emp2","emp3","z6","z15","z22")) {
  system(paste('bedtools intersect -a ',i,'_',igg,'.bed', ' -b /home/elyad/LabSync/ChipReference/gencode.v33.promoters.bed > ', i,'.',igg,'.promoters.bed',sep = ''))}
# "z6","z15","z22"
beds <- c()
for(i in c("emp1","emp2","emp3","z6","z15","z22")) {beds <- c(beds,paste( i,'.',igg,'.promoters.bed',sep = ''))}

anno <- read.csv('/home/elyad/LabSync/ChipReference/genecodeV33table.csv')

Bins <- 100
Genome_range <- 5000

#####
anno$chr <- as.character(anno$chr)

for(i in seq(-5000,5000,Genome_range/Bins)) {
  anno[,as.character(i)] <- anno$TSS + i
  anno[anno$strand == "-",as.character(i)] <- anno[anno$strand == "-","TSS"] - i}

###

for(bedfile in beds) {
  
  meta <- data.frame(matrix(ncol = length(anno),nrow = 0))
  colnames(meta) <- names(anno)
  
  chroms <- unique(anno$chr)
  for(chr in chroms) {
    
    bed <- read.csv(bedfile,sep='\t',col.names = c('chr','start','end','cpm'))
    bed <- bed[bed$chr == chr,]
    bed <- bed[!duplicated(bed),]
    
    annoT <- anno[anno$chr == chr,]
    
    meta <- rbind(meta,cbind(annoT[,1:8],apply(annoT[,-c(1:8)],1:2,function(x) {findcpm(coord = x,bedgraph = bed)})))
    
    print(chr)  }
  
  for(i in seq(-5000,5000,Genome_range/Bins)) {meta[,as.character(i)] <- as.numeric(unlist(meta[,as.character(i)]))}
  
  write.csv(meta, gsub("promoters.bed","metaPromoter.csv",bedfile))
  
  print(bedfile)}

####################################   diffbind analysis    ####################################

#BiocManager::install("DiffBind")
library(DiffBind)
library(data.table)

split_doubles <- function(row) {
  rows <- row[0,]
  for(g in seq(1:length(row$close_genes))) {
    for(i in seq(1:length(row[g,"close_genes"][[1]]))) {
      p <- row[g,]
      p$close_genes <- row[g,"close_genes"][[1]][i]
      p$gene_type <- row[g,"gene_type"][[1]][i]
      rows <- rbind(rows,p) }  }
  return(rows)}

############ diffbind analysis ##### 

tamoxifen <- dba(sampleSheet="samplesheetH3Ac.csv")
plot(tamoxifen)

tamoxifen <- dba.count(tamoxifen, minOverlap = 3)
tamoxifen <- dba.contrast(tamoxifen, categories=DBA_CONDITION)
tamoxifen <- dba.analyze(tamoxifen,method=DBA_ALL_METHODS)
tamoxifen.DB <- dba.report(tamoxifen,bCalled=TRUE,th=1,method=DBA_EDGER)


write.csv(tamoxifen.DB,'EDGERdifferentially_bound_sites.csv')


#### annotate peaks close to enhancers and TSS
annotation <- read.csv('/home/elyad/LabSync/ChipReference/genecodeV33table.csv',row.names=1)
annotation <- annotation[,c("gene_name","chr","start","end","TSS","gene_type")]

enhancers <- read.csv('/home/elyad/LabSync/ChipReference/fantom5_enhancers_GRCh38.bed',sep='\t',col.names = c("chr","start","end","gene_type"))
enhancers$TSS <- (enhancers$start+enhancers$end)/2
gene_name <- c()
for(i in 1:length(enhancers$chr)) {gene_name <- c(gene_name,paste0('enhancer',as.character(i)))}
enhancers$gene_name <- gene_name

annotation <- rbind(enhancers,annotation)
#
peaks <- read.csv('EDGERdifferentially_bound_sites.csv',row.names=1)
promoters <- c()
types <- c()
#
for(i in seq(1,nrow(peaks))) {
  genes <- annotation[annotation$chr == as.character(peaks[i,'seqnames']),]
  genes <- genes[genes$TSS > peaks[i,'start']-4000 & genes$TSS < peaks[i,'end']+4000,]
  promoters <- c(promoters,paste(genes$gene_name,collapse = ' '))
  types <- c(types,paste(genes$gene_type,collapse = ' ')) }

peaks$close_genes <- promoters
peaks$gene_type <- types

peaks$close_genes <- strsplit(peaks$close_genes,' ')
peaks$gene_type <- strsplit(peaks$gene_type,' ')

temp <- peaks[as.logical(lapply(peaks$close_genes, function(x) {length(x) >1 })),]
temp <- split_doubles(temp)

peaks <- peaks[!as.logical(lapply(peaks$close_genes, function(x) {length(x) >1 })),]
peaks <- rbind(peaks,temp) 

peaks$close_genes <- as.character(peaks$close_genes)
peaks$gene_type <- as.character(peaks$gene_type)

peaks[unlist(lapply(peaks$close_genes,function(x) {identical(peaks$close_genes[2], x)})),c("close_genes","gene_type")]  <- ""

peaks <- peaks[order(peaks$Fold,decreasing = T),]
peaks <- peaks[!(duplicated(peaks$close_genes) & peaks$close_genes != ""),]
write.csv(peaks,'EDGERannotated_differentially_bound_sites.csv')


###### data analysis #######

df <- read.csv('EDGERannotated_differentially_bound_sites.csv',row.names = 1)
df <- df[df$gene_type != "" & df$gene_type != "Enhancer",]

df$upreg <- df$FDR<=0.05 & df$Fold > 0
df$downreg <- df$FDR<=0.05 & df$Fold < -0
df$novel <- df$Called1 >= 2 & df$Called2 == 0 & df$FDR <= 0.05
df$uniqueWT <- df$Called1 == 0 & df$Called2 >= 2 & df$FDR <= 0.05
df$unchanged <- df$FDR > 0.9 & df$Called1 >= 2 & df$Called2 >= 2 & abs(df$fold) < 0.
df$close_genes <- as.character(df$close_genes)
#df[df$close_genes== "", "close_genes"] <- "AnnoLess"

df <- merge(df,read.csv('./tables/koDE.csv'),by.x="close_genes",by.y="gene_name",all.x = T)

df$onlyHyperAC <- df$Fold > 0.5 & df$FDR.x < 0.05 & abs(df$logFC) < 0.5
df[is.na(df$onlyHyperAC),'onlyHyperAC'] <- FALSE
df$correlatedHyper <- df$Fold > 0.5 & df$logFC > 0.5 & df$FDR.x < 0.05
df[is.na(df$correlatedHyper),'correlatedHyper'] <- FALSE
df$onlyHypoAC <- df$Fold < -0.5 & df$FDR.x < 0.05 & abs(df$logFC) < 0.5
df[is.na(df$onlyHypoAC),'onlyHypoAC'] <- FALSE
df$correlatedHypo <- df$Fold < -0.5 & df$logFC < -0.5 & df$FDR.x < 0.05
df[is.na(df$correlatedHypo),'correlatedHypo'] <- FALSE

df <- df[order(df$Fold,decreasing = T),]
df <- df[!duplicated(df$close_genes),]
summary(df)
write.csv(df,"differential_H3Ac.csv")


######### differential binding analysis on diffbind output ##########
df <- read.csv('EDGERannotated_differentially_bound_sites.csv',row.names = 1)
df <- df[df$gene_type != "" & df$gene_type != "Enhancer",]

df$upreg <- df$FDR<=0.05 & df$Fold > 0.5
df$downreg <- df$FDR<=0.05 & df$Fold < -0.5
df$novel <- df$Called1 >= 2 & df$Called2 == 0 & df$FDR <= 0.05
df$uniqueWT <- df$Called1 == 0 & df$Called2 >= 2 & df$FDR <= 0.05
df$unchanged <- df$FDR > 0.9 & df$Called1 >= 2 & df$Called2 >= 2 & abs(df$Fold) < 0.2
df$close_genes <- as.character(df$close_genes)
#df[df$close_genes== "", "close_genes"] <- "AnnoLess"

df <- merge(df,read.csv('./tables/koDE.csv'),by.x="close_genes",by.y="gene_name",all.x = T)

df <- df[order(df$logCPM,decreasing=T),]
df <- df[!(duplicated(df$end) & duplicated(df$start) & duplicated(df$width) & duplicated(df$Conc)),]
df <- df[!duplicated(df$close_genes),]

df$onlyHyperAC <- df$Fold > 0.5 & df$FDR.x < 0.05 & abs(df$logFC) < 0.5
df[is.na(df$onlyHyperAC),'onlyHyperAC'] <- FALSE
df$correlatedHyper <- df$Fold > 0.5 & df$logFC > 0.5 & df$FDR.x < 0.05
df[is.na(df$correlatedHyper),'correlatedHyper'] <- FALSE
df$onlyHypoAC <- df$Fold < -0.5 & df$FDR.x < 0.05 & abs(df$logFC) < 0.5
df[is.na(df$onlyHypoAC),'onlyHypoAC'] <- FALSE
df$correlatedHypo <- df$Fold < -0.5 & df$logFC < -0.5 & df$FDR.x < 0.05
df[is.na(df$correlatedHypo),'correlatedHypo'] <- FALSE


summary(df)
write.csv(df,"differential_H3Ac.csv")


########################differential expression ZMYM2 KO VS WT Teratomas ########################

library("edgeR")

group <- c(1,1,2,2)
counts <- read.csv('zmym2TeratomasReadCounts.csv')[-c(1,3)]

# DE analysis 
y <- DGEList(counts[-1], group, genes=counts[1])
# remove unexpressed genes
keep <- rowSums(cpm(y)>1) > 1
y <- y[keep, , keep.lib.sizes=FALSE]
# TMM normalization
y <- calcNormFactors(y,method = 'TMM')

# MDS plot ### figure 4 ###
plotMDS(y)                      

# fit model and calculate DE
design <- model.matrix(~group)
y <- estimateDisp(y, design)
y <- glmFit(y, design)
y <- glmLRT(y)
# list of ALL genes with logFC and FDR
DE <- data.frame(topTags(y,n = Inf,sort.by='p.value'))  
DE <- merge(read.csv('/home/elyad/reference_genome/data/gencode.v29primary_names_type_length.csv',row.names=1),DE,by='ensemblID')

write.csv(DE,'zmymm2KOteratomasDE.csv')

