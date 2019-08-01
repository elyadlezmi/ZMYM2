
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
counts <- read.csv('zmym2ReadCounts.csv')[-c(1,4)]

# DE analysis 
y <- DGEList(counts[c(-1,-2)], group, genes=counts[c(1,2)])
# remove unexpressed genes
keep <- rowSums(cpm(y)>1) > 1
y <- y[keep, , keep.lib.sizes=FALSE]
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
keep <- rowSums(cpm(y)>1) > 1
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


### create TSS table 

df <- read.csv('gencode.v29.annotation.gtf.gz',sep='\t',skip = 5, col.names = c('chr','x','feature','start','end','y','strand','z','w'))
df <- df[df$feature == 'gene',]
df$gene_name <- sub('gene_id .+; gene_name (.+); level.+','\\1',df[,'w'])
df$ensemblID <- sub('gene_id (.+\\.\\d+); gene_type.+','\\1',df[,'w'])
df <- df[,c('ensemblID','gene_name','chr','start','end','strand')]

df$TSS <- df$start
df[df$strand=='-','TSS'] <- df[df$strand=='-','end']
write.csv(df,'gencodeV29table.csv')

####### create bed for metagene

df <- read.csv('gencodeV29table.csv')
bed <- df[c('chr')]
bed$start <- df$TSS-5000
bed$end <- df$TSS+5000
bed$name <- df$ensemblID
bed$zero <- 0
bed$strand <- df$strand
options(scipen=999)
write.table(bed,'gencodeV29promoters_ALLgenes.bed',sep='\t',col.names = F,row.names = F, quote=F)

#####################################################################################################