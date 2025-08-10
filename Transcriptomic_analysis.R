library(tidyverse)
library(Seurat)
library(edgeR)
library(patchwork)
library(DistributionUtils)
library(ggplot2); theme_set(theme_linedraw()) #theme_set(theme_bw()) #theme_set(theme_minimal()) #theme_set(theme_light())
library(MetaCycle)
library(org.Dm.eg.db)
library(cowplot)
library(gridExtra)
library(grid)
library(plotROC)
library(flextable)
library(svglite)
library(VennDiagram)
library(ggpubr)
library(spatstat)
library(zoo)
library(ComplexHeatmap)
library(shadowtext)
library(ggprism)
library(fitdistrplus)
library(openxlsx)
library(data.table)
library(parallel)
library(circlize)
library(harmonicmeanp)
library(gt)
library(scales)
library(ggrepel)


library(AnnotationDbi)
library(GO.db)
library(memes)
library(org.Dm.eg.db)
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
library(BSgenome.Dmelanogaster.UCSC.dm6)
library(GenomicRanges)
library(GenomicFeatures)
library(MotifDb)
library(universalmotif)

library("ggseqlogo")


library(magrittr)
library(biomaRt)
#library(BiocFileCache)

#for enrichment analysis
#library(Biostrings); library(seqLogo); library(BCRANK)

'%!in%' <- function(x,y)!('%in%'(x,y))

circ_genes <- c('tim', 'per', 'Clk', 'cwo', 'vri', 'Pdp1', 'cry')
WD = "/Users/Teddy/Desktop/CCG_distribution"

color_list = c(C3_LD = '#FFCC88', C3_DD='#FF6644', C2_LD='#44BBFF', C2_DD='#4455FF')
color_list_dark = c(C3_LD = '#BB8844', C3_DD='#BB2200', C2_LD='#0077BB', C2_DD='#0011BB')

color_list_alt = c(C3_LD = '#FF8844', C3_DD='#FF6666', C2_LD='#44BBFF', C2_DD='#4455FF')
color_list_dark_alt = c(C3_LD = '#BB4400', C3_DD='#AA2222', C2_LD='#0077BB', C2_DD='#0011BB')


# 1. Create data_C2C3 from data3 (pcat 2-4; biological concat = group2)
################################################################################
# hi <- readRDS('data3.rds)
# hi$Type <- mutate(hi[['seurat_clusters']], Type = ifelse(seurat_clusters %in% c(0, 2), 'C3', ifelse(seurat_clusters %in% c(1, 3, 5), 'C2', '')))[['Type']]

data_C2C3 <- readRDS('data3.rds')

data_C2C3[['sample']] <- data_C2C3@active.ident
group2 <- c('CT3', 'ZT23', 'ZT3', 'ZT11', 'CT11', 'CT19', 'ZT15', 'CT7', 'CT23', 
            'ZT7', 'ZT31', 'CT27', 'CT35', 'ZT19', 'ZT35', 'CT43', 'ZT27', 'ZT43', 
            'CT15', 'ZT47', 'CT39', 'CT47', 'CT31', 'ZT39')
hold <- data_C2C3@meta.data[['sample1']]
for (i in 1:length(group2)) {
  hold = replace(hold, hold==unique(hold)[i], group2[i])
}
data_C2C3@meta.data[['group2']] <- hold
rm(hold, i, group2)


data_C2C3$mode1 <- factor(data_C2C3$mode1, levels = c('LD', 'DD'))
data_C2C3$sample <- factor(data_C2C3$sample, levels = c('C3', 'C2'))
data_C2C3$group2 <- factor(data_C2C3$group2, levels = c(paste0('CT', seq(3, 47, 4)), paste0('ZT', seq(3, 47, 4))))

data_C2C3@meta.data[['sample_mode']] <- data_C2C3[[c('sample', 'mode1')]] %>% 
  mutate(sample_mode <- paste(sample, mode1, sep='_'), .keep='unused') %>% unlist %>%
  factor(levels=c('C3_LD', 'C2_LD', 'C3_DD', 'C2_DD'))

data_C2C3@meta.data[['sample_mode_sample1']] <- data_C2C3[[c('sample', 'mode1', 'sample1')]] %>% 
  mutate(sample_mode_sample1 <- paste(sample, mode1, sample1, sep='_'), .keep='unused') %>% unlist

data_C2C3@meta.data[['sample_mode_group2']] <- data_C2C3[[c('sample', 'mode1', 'group2')]] %>% 
  mutate(sample_mode_group2 <- paste(sample, mode1, group2, sep='_'), .keep='unused') %>% unlist

data_C2C3@meta.data[['time4']] <- data_C2C3[[c('mode1', 'group2')]] %>% 
  mutate(time4 = as.numeric(as.numeric(gsub("[a-zA-Z]", "", group2))+ifelse(mode1=='LD', 0, 48))) %>%
  mutate(time4 = factor(time4, sort(unique(time4)))) %>% {.[, 'time4']}

data_C2C3@meta.data[['time2']] <- data_C2C3[[c('group2')]] %>% 
  mutate(time2 = as.numeric(gsub("[a-zA-Z]", "", group2))) %>%
  mutate(time2 = factor(time2, sort(unique(time2)))) %>% {.[, 'time2']}

saveRDS(data_C2C3, 'data_C2C3.rds')


subset_one_cluster_timepoint = function(data, Samp_Mode, timepoint, subsets) {
  data[[c('sample_mode', 'time1')]] %>% subset(sample_mode == Samp_Mode & time1 == timepoint) %>%
    rownames %>%
    sample %>% 
    {split(., ceiling(seq_along(.)/length(.)*subsets))} %>%
    {lapply(names(.), function(name, lst=.) data.frame(cells = lst[[name]], subset = name))} %>%
    do.call(what='rbind') %>%
    cbind(Sample_Mode = Samp_Mode, timepoint = timepoint)
}
add_cell_subsets = function(data, subsets) {
  res = data.frame(matrix(ncol=4, nrow=0)) %>%
    `colnames<-`(c('cells', 'subset', 'Cluster', 'timepoint'))
  
  for (s_m in unique(data$sample_mode)) {
    for (timepoint in unique(data$time1)) {
      res = rbind(res, subset_one_cluster_timepoint(data, s_m, timepoint, subsets))
    }
  }
  
  res %>%
    mutate(timepoint = as.integer(timepoint)+24*(as.integer(subset)-1)) %>%
    {inner_join(data.frame(cells = colnames(data)), .)} %>%
    return
}



# added_C2C3
{set.seed(0)
added_C2C3 <- data_C2C3
added_C2C3@meta.data[['time3']] <- added_C2C3[[c('time2')]] %>%
  mutate(time3 = (as.numeric(as.character(time2)))) %>% {.$time3}
  
add_cells <- added_C2C3[[c('sample_mode', 'time2')]] %>%
  {split(., list(.$sample_mode, .$time2))} %>% 
  lapply(function(DF) {
    N = nrow(DF)
    DF[sample(N, floor(N/2) + sample(c(0, 1), 1)*(N%%2), replace=FALSE), ]
  }) %>%
  do.call(what=rbind) %>%
  rownames %>% gsub('.*\\.', '', .) %>% 
  {subset(added_C2C3, cells = .)}

add_cells@meta.data[['time3']] <- add_cells[[c('time3')]] %>% 
  mutate(time3 = time3 %% 24 + 48) %>% {.$time3}

added_C2C3 <- merge(added_C2C3, add_cells)
rm(add_cells)

added_C2C3$sample <- factor(added_C2C3$sample, levels = c('C3', 'C2'))
added_C2C3$mode1 <- factor(added_C2C3$mode1, levels = c('LD', 'DD'))

added_C2C3@meta.data[['time6']] <- added_C2C3[[c('mode1', 'time3')]] %>% 
  mutate(time6 = time3+ifelse(mode1=='LD', 0, 72)) %>%
  mutate(time6 = factor(time6, sort(unique(time6)))) %>% {.[, 'time6']}

added_C2C3@meta.data[['time3']] <- added_C2C3[[c('time3')]] %>%
  mutate(time3 = factor(time3, levels=sort(unique(time3)))) %>% {.[, 'time3']}

added_C2C3@meta.data[['sample_mode_time3']] <- added_C2C3[[c('sample', 'mode1', 'time3')]] %>% 
  mutate(sample_mode_time3 <- paste(sample, mode1, time3, sep='_'), .keep='unused') %>% {.$sample_mode_time3}

added_C2C3_AGG <- AggregateExpression(added_C2C3, return.seurat=TRUE, group.by='sample_mode_time3')
added_C2C3_AGG@meta.data[['sample']] <- added_C2C3_AGG[[c('sample_mode_time3')]] %>% 
  mutate(sample = gsub('-.*', '', sample_mode_time3), .keep='unused') %>% unlist %>%
  factor(levels=c('C3', 'C2'))
added_C2C3_AGG@meta.data[['mode1']] <- added_C2C3_AGG[[c('sample_mode_time3')]] %>% 
  mutate(sample = sub('D-.*', 'D', sample_mode_time3) %>% sub('.*-', '', .), .keep='unused') %>% unlist %>%
  factor(levels=c('LD', 'DD'))
added_C2C3_AGG@meta.data[['time3']] <- added_C2C3_AGG[[c('sample_mode_time3')]] %>% 
  mutate(sample = gsub('.*-', '', sample_mode_time3), .keep='unused') %>% unlist %>%
  {factor(., levels=sort(unique(.)))}
added_C2C3_AGG@meta.data[['sample_mode']] <- added_C2C3_AGG[[c('sample', 'mode1')]] %>% 
  mutate(sample_mode <- paste(sample, mode1, sep='_'), .keep='unused') %>% unlist %>%
  factor(levels=c('C3_LD', 'C2_LD', 'C3_DD', 'C2_DD'))
}

saveRDS(added_C2C3, 'added_C2C3.rds')
saveRDS(added_C2C3_AGG, 'added_C2C3_AGG.rds')




# added_C2C3_pool3
{set.seed(0)
added_C2C3_pool3 <- data_C2C3
added_C2C3_pool3@meta.data[['time3']] <- added_C2C3_pool3[[c('time2')]] %>%
  mutate(time3 = (as.numeric(as.character(time2)))) %>% {.$time3}
  
add_cells <- added_C2C3_pool3[[c('sample_mode', 'time1')]] %>%
  {split(., list(.$sample_mode, .$time1))} %>% 
  lapply(function(DF) {
    N = nrow(DF)
    DF[sample(N, floor(N/2) + sample(c(0, 1), 1)*(N%%2), replace=FALSE), ]
  }) %>%
  do.call(what=rbind) %>%
  rownames %>% gsub('.*\\.', '', .) %>% 
  {subset(added_C2C3_pool3, cells = .)}

add_cells@meta.data[['time3']] <- add_cells[[c('time3')]] %>% 
  mutate(time3 = time3 %% 24 + 48) %>% {.$time3}

added_C2C3_pool3 <- merge(added_C2C3_pool3, add_cells)
rm(add_cells)

added_C2C3_pool3$sample <- factor(added_C2C3_pool3$sample, levels = c('C3', 'C2'))
added_C2C3_pool3$mode1 <- factor(added_C2C3_pool3$mode1, levels = c('LD', 'DD'))

added_C2C3_pool3@meta.data[['time6']] <- added_C2C3_pool3[[c('mode1', 'time3')]] %>% 
  mutate(time6 = time3+ifelse(mode1=='LD', 0, 72)) %>%
  mutate(time6 = factor(time6, sort(unique(time6)))) %>% {.[, 'time6']}

added_C2C3_pool3@meta.data[['time3']] <- added_C2C3_pool3[[c('time3')]] %>%
  mutate(time3 = factor(time3, levels=sort(unique(time3)))) %>% {.[, 'time3']}

added_C2C3_pool3@meta.data[['sample_mode_time3']] <- added_C2C3_pool3[[c('sample', 'mode1', 'time3')]] %>% 
  mutate(sample_mode_time3 <- paste(sample, mode1, time3, sep='_'), .keep='unused') %>% {.$sample_mode_time3}

added_C2C3_pool3_AGG <- AggregateExpression(added_C2C3_pool3, return.seurat=TRUE, group.by='sample_mode_time3')
added_C2C3_pool3_AGG@meta.data[['sample']] <- added_C2C3_pool3_AGG[[c('sample_mode_time3')]] %>% 
  mutate(sample = gsub('-.*', '', sample_mode_time3), .keep='unused') %>% unlist %>%
  factor(levels=c('C3', 'C2'))
added_C2C3_pool3_AGG@meta.data[['mode1']] <- added_C2C3_pool3_AGG[[c('sample_mode_time3')]] %>% 
  mutate(sample = sub('D-.*', 'D', sample_mode_time3) %>% sub('.*-', '', .), .keep='unused') %>% unlist %>%
  factor(levels=c('LD', 'DD'))
added_C2C3_pool3_AGG@meta.data[['time3']] <- added_C2C3_pool3_AGG[[c('sample_mode_time3')]] %>% 
  mutate(sample = gsub('.*-', '', sample_mode_time3), .keep='unused') %>% unlist %>%
  {factor(., levels=sort(unique(.)))}
added_C2C3_pool3_AGG@meta.data[['sample_mode']] <- added_C2C3_pool3_AGG[[c('sample', 'mode1')]] %>% 
  mutate(sample_mode <- paste(sample, mode1, sep='_'), .keep='unused') %>% unlist %>%
  factor(levels=c('C3_LD', 'C2_LD', 'C3_DD', 'C2_DD'))
}

saveRDS(added_C2C3_pool3, 'added_C2C3_pool3.rds')
saveRDS(added_C2C3_pool3_AGG, 'added_C2C3_pool3_AGG.rds')


# added_C2C3_pool
{set.seed(0)
added_C2C3_pool <- data_C2C3
added_C2C3_pool@meta.data[['time3']] <- added_C2C3_pool[[c('sample_mode', 'time1')]] %>%
  {split(., list(.$sample_mode, .$time1))} %>% 
  lapply(function(DF) {
    data.frame(DF, rep = sample(seq(1:nrow(DF))%%3, replace=FALSE)) %>%
      mutate(time3 = as.numeric(as.character(time1))+24*rep)
  }) %>%
  do.call(what=rbind) %>%
  as.data.frame %>%
  {`rownames<-`(., gsub('.*\\.', '', rownames(.)))} %>%
  {.[Cells(added_C2C3_pool), ]} %>%
  {.$time3}



added_C2C3_pool$group2 <- NULL
added_C2C3_pool$sample_mode_group2 <- NULL
added_C2C3_pool$time4 <- NULL
added_C2C3_pool$time2 <- NULL

added_C2C3_pool@meta.data[['time6']] <- added_C2C3_pool[[c('mode1', 'time3')]] %>% 
  mutate(time6 = time3+ifelse(mode1=='LD', 0, 72)) %>%
  mutate(time6 = factor(time6, sort(unique(time6)))) %>% {.[, 'time6']}

added_C2C3_pool@meta.data[['time3']] <- added_C2C3_pool[[c('time3')]] %>%
  mutate(time3 = factor(time3, levels=sort(unique(time3)))) %>% {.[, 'time3']}

added_C2C3_pool@meta.data[['sample_mode_time3']] <- added_C2C3_pool[[c('sample', 'mode1', 'time3')]] %>% 
  mutate(sample_mode_time3 <- paste(sample, mode1, time3, sep='_'), .keep='unused') %>% {.$sample_mode_time3}

added_C2C3_pool_AGG <- AggregateExpression(added_C2C3_pool, return.seurat=TRUE, group.by='sample_mode_time3')
added_C2C3_pool_AGG@meta.data[['sample']] <- added_C2C3_pool_AGG[[c('sample_mode_time3')]] %>% 
  mutate(sample = gsub('-.*', '', sample_mode_time3), .keep='unused') %>% unlist %>%
  factor(levels=c('C3', 'C2'))
added_C2C3_pool_AGG@meta.data[['mode1']] <- added_C2C3_pool_AGG[[c('sample_mode_time3')]] %>% 
  mutate(sample = sub('D-.*', 'D', sample_mode_time3) %>% sub('.*-', '', .), .keep='unused') %>% unlist %>%
  factor(levels=c('LD', 'DD'))
added_C2C3_pool_AGG@meta.data[['time3']] <- added_C2C3_pool_AGG[[c('sample_mode_time3')]] %>% 
  mutate(sample = gsub('.*-', '', sample_mode_time3), .keep='unused') %>% unlist %>%
  {factor(., levels=sort(unique(.)))}
added_C2C3_pool_AGG@meta.data[['sample_mode']] <- added_C2C3_pool_AGG[[c('sample', 'mode1')]] %>% 
  mutate(sample_mode <- paste(sample, mode1, sep='_'), .keep='unused') %>% unlist %>%
  factor(levels=c('C3_LD', 'C2_LD', 'C3_DD', 'C2_DD'))
}

saveRDS(added_C2C3_pool, 'added_C2C3_pool.rds')
saveRDS(added_C2C3_pool_AGG, 'added_C2C3_pool_AGG.rds')



# Repeating added_C2C3_pool
reps_adding_pool <- lapply(c(0:2), function(s) {
  {set.seed(s)
    added_C2C3_pool = data_C2C3
    added_C2C3_pool@meta.data[['time3']] = added_C2C3_pool[[c('sample_mode', 'time1')]] %>%
      {split(., list(.$sample_mode, .$time1))} %>% 
      lapply(function(DF) {
        data.frame(DF, rep = sample(seq(1:nrow(DF))%%3, replace=FALSE)) %>%
          mutate(time3 = as.numeric(as.character(time1))+24*rep)
      }) %>%
      do.call(what=rbind) %>%
      as.data.frame %>%
      {`rownames<-`(., gsub('.*\\.', '', rownames(.)))} %>%
      {.[Cells(added_C2C3_pool), ]} %>%
      {.$time3}
    
    
    added_C2C3_pool$group2 = NULL
    added_C2C3_pool$sample_mode_group2 = NULL
    added_C2C3_pool$time4 = NULL
    added_C2C3_pool$time2 = NULL
    
    added_C2C3_pool@meta.data[['time6']] = added_C2C3_pool[[c('mode1', 'time3')]] %>% 
      mutate(time6 = time3+ifelse(mode1=='LD', 0, 72)) %>%
      mutate(time6 = factor(time6, sort(unique(time6)))) %>% {.[, 'time6']}
    
    added_C2C3_pool@meta.data[['time3']] = added_C2C3_pool[[c('time3')]] %>%
      mutate(time3 = factor(time3, levels=sort(unique(time3)))) %>% {.[, 'time3']}
    
    added_C2C3_pool@meta.data[['sample_mode_time3']] = added_C2C3_pool[[c('sample', 'mode1', 'time3')]] %>% 
      mutate(sample_mode_time3 = paste(sample, mode1, time3, sep='_'), .keep='unused') %>% {.$sample_mode_time3}
    
    added_C2C3_pool_AGG = AggregateExpression(added_C2C3_pool, return.seurat=TRUE, group.by='sample_mode_time3')
    added_C2C3_pool_AGG@meta.data[['sample']] = added_C2C3_pool_AGG[[c('sample_mode_time3')]] %>% 
      mutate(sample = gsub('-.*', '', sample_mode_time3), .keep='unused') %>% unlist %>%
      factor(levels=c('C3', 'C2'))
    added_C2C3_pool_AGG@meta.data[['mode1']] = added_C2C3_pool_AGG[[c('sample_mode_time3')]] %>% 
      mutate(sample = sub('D-.*', 'D', sample_mode_time3) %>% sub('.*-', '', .), .keep='unused') %>% unlist %>%
      factor(levels=c('LD', 'DD'))
    added_C2C3_pool_AGG@meta.data[['time3']] = added_C2C3_pool_AGG[[c('sample_mode_time3')]] %>% 
      mutate(sample = gsub('.*-', '', sample_mode_time3), .keep='unused') %>% unlist %>%
      {factor(., levels=sort(unique(.)))}
    added_C2C3_pool_AGG@meta.data[['sample_mode']] = added_C2C3_pool_AGG[[c('sample', 'mode1')]] %>% 
      mutate(sample_mode = paste(sample, mode1, sep='_'), .keep='unused') %>% unlist %>%
      factor(levels=c('C3_LD', 'C2_LD', 'C3_DD', 'C2_DD'))
    added_C2C3_pool_AGG@meta.data[['time6']] <- added_C2C3_pool_AGG[[c('mode1', 'time3')]] %>% 
      mutate(time6 = as.numeric(as.numeric(as.character(time3))+ifelse(mode1=='LD', 0, 72))) %>%
      mutate(time6 = factor(time6, sort(unique(time6)))) %>% {.[, 'time6']}
    
    return(list(avg = added_C2C3_pool, AGG = added_C2C3_pool_AGG))
  }
}) %>%
  `names<-`(paste0('R', c(0:2)))

#saveRDS(reps_adding_pool, 'reps_adding_pool.rds')
reps_adding_pool.rds <- readRDS('reps_adding_pool.rds')

avg_reps_adding_pool_AGG <- lapply(names(reps_adding_pool), function(name_L2) {
  reps_adding_pool[[name_L2]]$AGG@assays$RNA$counts %>%
    as.data.frame %>% 
    {mutate(., Gene = rownames(.))} %>%
    pivot_longer(cols=-Gene, names_to='cond', values_to=name_L2)
}) %>%
  Reduce(f=merge) %>%
  mutate(mean_value = rowMeans(subset(., select=c(R0, R1, R2))), 
         R0=NULL, R1=NULL, R2=NULL) %>%
  pivot_wider(names_from=cond, values_from=mean_value) %>% 
  as.data.frame %>%
  {`rownames<-`(., .$Gene)} %>%
  mutate(Gene = NULL) %>%
  CreateSeuratObject(meta.data=reps_adding_pool$R0$AGG@meta.data) %>%
  NormalizeData()

saveRDS(avg_reps_adding_pool_AGG, 'avg_reps_adding_pool_AGG.rds')


{avg_reps_adding_pool_AGG_01_1d <- avg_reps_adding_pool_AGG@assays$RNA$data %>% 
  expm1 %>%
  as.data.frame %>%
  {mutate(., Gene = rownames(.))} %>%
  pivot_longer(-Gene, values_to='expr') %>%
  mutate(`Sample-Mode` = gsub('D-.*', 'D', name), 
         timepoint = gsub('.*-', '', name), .keep='unused') %>% 
  mutate(day = as.numeric(as.character(timepoint)) %/% 24) %>%
  {split(., list(.$`Sample-Mode`, .$Gene, .$day))} %>% 
  lapply(function(df) {
    mutate(df, expr = (expr - min(expr))) %>% 
      mutate(expr = expr/max(expr))
  }) %>%
  do.call(what='rbind') %>%
  as.data.frame %>%
  mutate(cond = paste(`Sample-Mode`, timepoint, sep='_'),
         day=NULL, .keep='unused') %>%
  pivot_wider(names_from = 'cond', values_from = 'expr') %>%
  as.matrix %>%
  {`rownames<-`(., .[,'Gene'])} %>%
  {.[, -which(colnames(.)=='Gene')]} %>%
  CreateSeuratObject

avg_reps_adding_pool_AGG_01_1d@meta.data[['sample-mode_time3']] <- Cells(avg_reps_adding_pool_AGG_01_1d)
avg_reps_adding_pool_AGG_01_1d@meta.data[['sample']] <- avg_reps_adding_pool_AGG_01_1d[[c('sample-mode_time3')]] %>% 
  mutate(sample = gsub('-.*', '', `sample-mode_time3`), .keep='unused') %>% unlist %>%
  factor(levels=c('C3', 'C2'))
avg_reps_adding_pool_AGG_01_1d@meta.data[['mode1']] <- avg_reps_adding_pool_AGG_01_1d[[c('sample-mode_time3')]] %>% 
  mutate(sample = sub('_.*', '', `sample-mode_time3`) %>% sub('.*-', '', .), .keep='unused') %>% unlist %>%
  factor(levels=c('LD', 'DD'))
avg_reps_adding_pool_AGG_01_1d@meta.data[['time3']] <- avg_reps_adding_pool_AGG_01_1d[[c('sample-mode_time3')]] %>% 
  mutate(sample = gsub('.*_', '', `sample-mode_time3`), .keep='unused') %>% unlist %>%
  factor(levels = seq(3, 71, 4))
avg_reps_adding_pool_AGG_01_1d@meta.data[['sample_mode']] <- avg_reps_adding_pool_AGG_01_1d[[c('sample', 'mode1')]] %>% 
  mutate(sample_mode <- paste(sample, mode1, sep='_'), .keep='unused') %>% unlist %>%
  factor(levels=c('C3_LD', 'C2_LD', 'C3_DD', 'C2_DD'))
avg_reps_adding_pool_AGG_01_1d@meta.data[['time6']] <- avg_reps_adding_pool_AGG_01_1d[[c('mode1', 'time3')]] %>% 
  mutate(time6 = as.numeric(as.numeric(as.character(time3))+ifelse(mode1=='LD', 0, 72))) %>%
  mutate(time6 = factor(time6, sort(unique(time6)))) %>% {.[, 'time6']}

saveRDS(avg_reps_adding_pool_AGG_01_1d, 'avg_reps_adding_pool_AGG_01_1d.rds')
}

{avg_reps_adding_pool_AGG_1_3d <- avg_reps_adding_pool_AGG@assays$RNA$data %>% 
    expm1 %>%
    as.data.frame %>%
    {mutate(., Gene = rownames(.))} %>%
    pivot_longer(-Gene, values_to='expr') %>%
    mutate(`Sample-Mode` = gsub('D-.*', 'D', name), 
           timepoint = gsub('.*-', '', name), .keep='unused') %>% 
    {split(., list(.$`Sample-Mode`, .$Gene))} %>% 
    lapply(function(df) {
      mutate(df, expr = expr/max(expr))
    }) %>%
    do.call(what='rbind') %>%
    as.data.frame %>%
    mutate(cond = paste(`Sample-Mode`, timepoint, sep='_'),
           day=NULL, .keep='unused') %>%
    pivot_wider(names_from = 'cond', values_from = 'expr') %>%
    as.matrix %>%
    {`rownames<-`(., .[,'Gene'])} %>%
    {.[, -which(colnames(.)=='Gene')]} %>%
    CreateSeuratObject
  
  avg_reps_adding_pool_AGG_1_3d@meta.data[['sample-mode_time3']] <- Cells(avg_reps_adding_pool_AGG_1_3d)
  avg_reps_adding_pool_AGG_1_3d@meta.data[['sample']] <- avg_reps_adding_pool_AGG_1_3d[[c('sample-mode_time3')]] %>% 
    mutate(sample = gsub('-.*', '', `sample-mode_time3`), .keep='unused') %>% unlist %>%
    factor(levels=c('C3', 'C2'))
  avg_reps_adding_pool_AGG_1_3d@meta.data[['mode1']] <- avg_reps_adding_pool_AGG_1_3d[[c('sample-mode_time3')]] %>% 
    mutate(sample = sub('_.*', '', `sample-mode_time3`) %>% sub('.*-', '', .), .keep='unused') %>% unlist %>%
    factor(levels=c('LD', 'DD'))
  avg_reps_adding_pool_AGG_1_3d@meta.data[['time3']] <- avg_reps_adding_pool_AGG_1_3d[[c('sample-mode_time3')]] %>% 
    mutate(sample = gsub('.*_', '', `sample-mode_time3`), .keep='unused') %>% unlist %>%
    factor(levels = seq(3, 71, 4))
  avg_reps_adding_pool_AGG_1_3d@meta.data[['sample_mode']] <- avg_reps_adding_pool_AGG_1_3d[[c('sample', 'mode1')]] %>% 
    mutate(sample_mode <- paste(sample, mode1, sep='_'), .keep='unused') %>% unlist %>%
    factor(levels=c('C3_LD', 'C2_LD', 'C3_DD', 'C2_DD'))
  avg_reps_adding_pool_AGG_1_3d@meta.data[['time6']] <- avg_reps_adding_pool_AGG_1_3d[[c('mode1', 'time3')]] %>% 
    mutate(time6 = as.numeric(as.numeric(as.character(time3))+ifelse(mode1=='LD', 0, 72))) %>%
    mutate(time6 = factor(time6, sort(unique(time6)))) %>% {.[, 'time6']}
  
  saveRDS(avg_reps_adding_pool_AGG_1_3d, 'avg_reps_adding_pool_AGG_1_3d.rds')
}

avg_data_adding_pool_AGG_01_1d <- AggregateExpression(avg_reps_adding_pool_AGG_01_1d, group.by = 'sample-mode_time3')$RNA %>% 
  as.data.frame %>%
  {mutate(., Gene = rownames(.))} %>%
  pivot_longer(cols=-'Gene', names_to='sample_mode_time', values_to='avg') %>%
  mutate(Sample = gsub('-.*', '', sample_mode_time), 
         Mode = gsub('D-.*', 'D', sample_mode_time) %>% gsub('.*-', '', .), 
         timepoint = gsub('.*-', '', sample_mode_time),
         sample_mode_time = gsub('_', '-', sample_mode_time))

saveRDS(avg_data_adding_pool_AGG_01_1d, 'avg_data_adding_pool_AGG_01_1d.rds')

data_C2C3_AGG <- AggregateExpression(data_C2C3, return.seurat=TRUE, group.by='sample_mode_group2') %>%
  NormalizeData
data_C2C3_AGG@meta.data[['sample']] <- data_C2C3_AGG[[c('sample_mode_group2')]] %>% 
  mutate(sample = gsub('-.*', '', sample_mode_group2), .keep='unused') %>% unlist %>%
  factor(levels=c('C3', 'C2'))
data_C2C3_AGG@meta.data[['mode1']] <- data_C2C3_AGG[[c('sample_mode_group2')]] %>% 
  mutate(sample = sub('D-.*', 'D', sample_mode_group2) %>% sub('.*-', '', .), .keep='unused') %>% unlist %>%
  factor(levels=c('LD', 'DD'))
data_C2C3_AGG@meta.data[['group2']] <- data_C2C3_AGG[[c('sample_mode_group2')]] %>% 
  mutate(sample = gsub('.*-', '', sample_mode_group2), .keep='unused') %>% unlist %>%
  factor(levels = c(paste0('CT', seq(3, 47, 4)), paste0('ZT', seq(3, 47, 4))))
data_C2C3_AGG@meta.data[['sample_mode']] <- data_C2C3_AGG[[c('sample', 'mode1')]] %>% 
  mutate(sample_mode <- paste(sample, mode1, sep='_'), .keep='unused') %>% unlist %>%
  factor(levels=c('C3_LD', 'C2_LD', 'C3_DD', 'C2_DD'))
data_C2C3_AGG@meta.data[['time4']] <- data_C2C3_AGG[[c('mode1', 'group2')]] %>% 
  mutate(time4 = as.numeric(as.numeric(gsub("[a-zA-Z]", "", group2))+ifelse(mode1=='LD', 0, 48))) %>%
  mutate(time4 = factor(time4, sort(unique(time4)))) %>% {.[, 'time4']}



data_C2C3_AGG_CLRnorm <- AggregateExpression(data_C2C3, return.seurat=TRUE, group.by='sample_mode_group2', normalization.method = 'CLR')
data_C2C3_AGG_CLRnorm@meta.data[['sample']] <- data_C2C3_AGG_CLRnorm[[c('sample_mode_group2')]] %>% 
  mutate(sample = gsub('-.*', '', sample_mode_group2), .keep='unused') %>% unlist %>%
  factor(levels=c('C3', 'C2'))
data_C2C3_AGG_CLRnorm@meta.data[['mode1']] <- data_C2C3_AGG_CLRnorm[[c('sample_mode_group2')]] %>% 
  mutate(sample = sub('D-.*', 'D', sample_mode_group2) %>% sub('.*-', '', .), .keep='unused') %>% unlist %>%
  factor(levels=c('LD', 'DD'))
data_C2C3_AGG_CLRnorm@meta.data[['group2']] <- data_C2C3_AGG_CLRnorm[[c('sample_mode_group2')]] %>% 
  mutate(sample = gsub('.*-', '', sample_mode_group2), .keep='unused') %>% unlist %>%
  factor(levels = c(paste0('CT', seq(3, 47, 4)), paste0('ZT', seq(3, 47, 4))))
data_C2C3_AGG_CLRnorm@meta.data[['sample_mode']] <- data_C2C3_AGG_CLRnorm[[c('sample', 'mode1')]] %>% 
  mutate(sample_mode <- paste(sample, mode1, sep='_'), .keep='unused') %>% unlist %>%
  factor(levels=c('C3_LD', 'C2_LD', 'C3_DD', 'C2_DD'))


saveRDS(data_C2C3_AGG, 'data_C2C3_AGG.rds')
saveRDS(data_C2C3_AGG_CLRnorm, 'data_C2C3_AGG_CLRnorm.rds')

# All genes scale between 0-1 over 2 days
data_C2C3_AGG_01_2d <- data_C2C3_AGG@assays$RNA$data %>% 
  expm1 %>%
  as.data.frame %>%
  {mutate(., Gene = rownames(.))} %>%
  pivot_longer(-Gene, values_to='expr') %>%
  mutate(`Sample-Mode` = gsub('D-.*', 'D', name), 
         timepoint = gsub('.*-', '', name), .keep='unused') %>%
  {split(., list(.$`Sample-Mode`, .$Gene))} %>% 
  lapply(function(df) {
    mutate(df, expr = (expr - min(expr))) %>% 
      mutate(expr = expr/max(expr))
  }) %>%
  do.call(what='rbind') %>%
  as.data.frame %>%
  mutate(cond = paste(`Sample-Mode`, timepoint, sep='_'), .keep='unused') %>%
  pivot_wider(names_from = 'cond', values_from = 'expr') %>%
  as.matrix %>%
  {`rownames<-`(., .[,'Gene'])} %>%
  {.[, -which(colnames(.)=='Gene')]} %>%
  CreateSeuratObject

data_C2C3_AGG_01_2d@meta.data[['sample-mode_group2']] <- Cells(data_C2C3_AGG_01_2d)
data_C2C3_AGG_01_2d@meta.data[['sample']] <- data_C2C3_AGG_01_2d[[c('sample-mode_group2')]] %>% 
  mutate(sample = gsub('-.*', '', `sample-mode_group2`), .keep='unused') %>% unlist %>%
  factor(levels=c('C3', 'C2'))
data_C2C3_AGG_01_2d@meta.data[['mode1']] <- data_C2C3_AGG_01_2d[[c('sample-mode_group2')]] %>% 
  mutate(sample = sub('_.*', '', `sample-mode_group2`) %>% sub('.*-', '', .), .keep='unused') %>% unlist %>%
  factor(levels=c('LD', 'DD'))
data_C2C3_AGG_01_2d@meta.data[['group2']] <- data_C2C3_AGG_01_2d[[c('sample-mode_group2')]] %>% 
  mutate(sample = gsub('.*_', '', `sample-mode_group2`), .keep='unused') %>% unlist %>%
  factor(levels = c(paste0('CT', seq(3, 47, 4)), paste0('ZT', seq(3, 47, 4))))
data_C2C3_AGG_01_2d@meta.data[['sample_mode']] <- data_C2C3_AGG_01_2d[[c('sample', 'mode1')]] %>% 
  mutate(sample_mode <- paste(sample, mode1, sep='_'), .keep='unused') %>% unlist %>%
  factor(levels=c('C3_LD', 'C2_LD', 'C3_DD', 'C2_DD'))
data_C2C3_AGG_01_2d@meta.data[['time4']] <- data_C2C3_AGG_01_2d[[c('mode1', 'group2')]] %>% 
  mutate(time4 = as.numeric(as.numeric(gsub("[a-zA-Z]", "", group2))+ifelse(mode1=='LD', 0, 48))) %>%
  mutate(time4 = factor(time4, sort(unique(time4)))) %>% {.[, 'time4']}

saveRDS(data_C2C3_AGG_01_2d, 'data_C2C3_AGG_01_2d.rds')


# All genes scale between -0.5-0.5 each day
data_C2C3_AGG_01_1d <- data_C2C3_AGG@assays$RNA$data %>% 
  expm1 %>%
  as.data.frame %>%
  {mutate(., Gene = rownames(.))} %>%
  pivot_longer(-Gene, values_to='expr') %>%
  mutate(`Sample-Mode` = gsub('D-.*', 'D', name), 
         timepoint = gsub('.*-', '', name), .keep='unused') %>% 
  mutate(time = as.numeric(gsub("\\D", "", timepoint))) %>% 
  mutate(rep = ceiling(time/24), .keep='unused') %>%
  {split(., list(.$`Sample-Mode`, .$Gene, .$rep))} %>% 
  lapply(function(df) {
    mutate(df, expr = (expr - min(expr))) %>% 
      mutate(expr = expr/max(expr))
  }) %>%
  do.call(what='rbind') %>%
  as.data.frame %>%
  mutate(cond = paste(`Sample-Mode`, timepoint, sep='_'), 
         rep=NULL, .keep='unused') %>%
  pivot_wider(names_from = 'cond', values_from = 'expr') %>%
  as.matrix %>%
  {`rownames<-`(., .[,'Gene'])} %>%
  {.[, -which(colnames(.)=='Gene')]} %>%
  CreateSeuratObject

data_C2C3_AGG_01_1d@meta.data[['sample-mode_group2']] <- Cells(data_C2C3_AGG_01_1d)
data_C2C3_AGG_01_1d@meta.data[['sample']] <- data_C2C3_AGG_01_1d[[c('sample-mode_group2')]] %>% 
  mutate(sample = gsub('-.*', '', `sample-mode_group2`), .keep='unused') %>% unlist %>%
  factor(levels=c('C3', 'C2'))
data_C2C3_AGG_01_1d@meta.data[['mode1']] <- data_C2C3_AGG_01_1d[[c('sample-mode_group2')]] %>% 
  mutate(sample = sub('_.*', '', `sample-mode_group2`) %>% sub('.*-', '', .), .keep='unused') %>% unlist %>%
  factor(levels=c('LD', 'DD'))
data_C2C3_AGG_01_1d@meta.data[['group2']] <- data_C2C3_AGG_01_1d[[c('sample-mode_group2')]] %>% 
  mutate(sample = gsub('.*_', '', `sample-mode_group2`), .keep='unused') %>% unlist %>%
  factor(levels = c(paste0('CT', seq(3, 47, 4)), paste0('ZT', seq(3, 47, 4))))
data_C2C3_AGG_01_1d@meta.data[['sample_mode']] <- data_C2C3_AGG_01_1d[[c('sample', 'mode1')]] %>% 
  mutate(sample_mode <- paste(sample, mode1, sep='_'), .keep='unused') %>% unlist %>%
  factor(levels=c('C3_LD', 'C2_LD', 'C3_DD', 'C2_DD'))
data_C2C3_AGG_01_1d@meta.data[['time4']] <- data_C2C3_AGG_01_1d[[c('mode1', 'group2')]] %>% 
  mutate(time4 = as.numeric(as.numeric(gsub("[a-zA-Z]", "", group2))+ifelse(mode1=='LD', 0, 48))) %>%
  mutate(time4 = factor(time4, sort(unique(time4)))) %>% {.[, 'time4']}

saveRDS(data_C2C3_AGG_01_1d, 'data_C2C3_AGG_01_1d.rds')

# All genes scale to max of 1 over 2 days
data_C2C3_AGG_max2 <- data_C2C3_AGG@assays$RNA$data %>% 
  expm1 %>%
  as.data.frame %>%
  {mutate(., Gene = rownames(.))} %>%
  pivot_longer(-Gene, values_to='expr') %>%
  mutate(`Sample-Mode` = gsub('D-.*', 'D', name), 
         timepoint = gsub('.*-', '', name), .keep='unused') %>%
  {split(., list(.$`Sample-Mode`, .$Gene))} %>% 
  lapply(function(df) {
    mutate(df, expr = expr/max(expr))
  }) %>%
  do.call(what='rbind') %>%
  as.data.frame %>%
  mutate(cond = paste(`Sample-Mode`, timepoint, sep='_'), .keep='unused') %>%
  pivot_wider(names_from = 'cond', values_from = 'expr') %>%
  as.matrix %>%
  {`rownames<-`(., .[,'Gene'])} %>%
  {.[, -which(colnames(.)=='Gene')]} %>%
  CreateSeuratObject

data_C2C3_AGG_max2@meta.data[['sample-mode_group2']] <- Cells(data_C2C3_AGG_max2)
data_C2C3_AGG_max2@meta.data[['sample']] <- data_C2C3_AGG_max2[[c('sample-mode_group2')]] %>% 
  mutate(sample = gsub('-.*', '', `sample-mode_group2`), .keep='unused') %>% unlist %>%
  factor(levels=c('C3', 'C2'))
data_C2C3_AGG_max2@meta.data[['mode1']] <- data_C2C3_AGG_max2[[c('sample-mode_group2')]] %>% 
  mutate(sample = sub('_.*', '', `sample-mode_group2`) %>% sub('.*-', '', .), .keep='unused') %>% unlist %>%
  factor(levels=c('LD', 'DD'))
data_C2C3_AGG_max2@meta.data[['group2']] <- data_C2C3_AGG_max2[[c('sample-mode_group2')]] %>% 
  mutate(sample = gsub('.*_', '', `sample-mode_group2`), .keep='unused') %>% unlist %>%
  factor(levels = c(paste0('CT', seq(3, 47, 4)), paste0('ZT', seq(3, 47, 4))))
data_C2C3_AGG_max2@meta.data[['sample_mode']] <- data_C2C3_AGG_max2[[c('sample', 'mode1')]] %>% 
  mutate(sample_mode <- paste(sample, mode1, sep='_'), .keep='unused') %>% unlist %>%
  factor(levels=c('C3_LD', 'C2_LD', 'C3_DD', 'C2_DD'))

saveRDS(data_C2C3_AGG_max2, 'data_C2C3_AGG_max2.rds')


data_C2C3_AGG_max1 <- data_C2C3_AGG@assays$RNA$data %>% 
  expm1 %>% 
  as.data.frame %>% 
  {mutate(., Gene = rownames(.))} %>% 
  pivot_longer(-Gene, values_to='expr') %>%
  mutate(`Sample-Mode` = gsub('D-.*', 'D', name), 
         timepoint = gsub('.*-', '', name), .keep='unused') %>%
  mutate(rep = gsub('.*T', '', timepoint) %>%
           as.integer %>%
           {.>24} %>% as.integer %>% {.+1}) %>% 
  {split(., list(.$`Sample-Mode`, .$rep, .$Gene))} %>%
  lapply(function(df) {
    mutate(df, expr = expr/max(expr))
  }) %>%
  do.call(what='rbind') %>%
  as.data.frame %>%
  mutate(cond = paste(`Sample-Mode`, timepoint, sep='_'),
         rep=NULL, .keep='unused') %>%
  pivot_wider(names_from = 'cond', values_from = 'expr') %>% 
  as.data.frame %>%
  {`rownames<-`(., .[,'Gene'])} %>% 
  mutate(across(everything(), .fns = ~replace_na(.,0))) %>%
  as.matrix %>%
  {.[, -which(colnames(.)=='Gene')]} %>%
  CreateSeuratObject

data_C2C3_AGG_max1@meta.data[['sample-mode_group2']] <- Cells(data_C2C3_AGG_max1)
data_C2C3_AGG_max1@meta.data[['sample']] <- data_C2C3_AGG_max1[[c('sample-mode_group2')]] %>% 
  mutate(sample = gsub('-.*', '', `sample-mode_group2`), .keep='unused') %>% unlist %>%
  factor(levels=c('C3', 'C2'))
data_C2C3_AGG_max1@meta.data[['mode1']] <- data_C2C3_AGG_max1[[c('sample-mode_group2')]] %>% 
  mutate(sample = sub('_.*', '', `sample-mode_group2`) %>% sub('.*-', '', .), .keep='unused') %>% unlist %>%
  factor(levels=c('LD', 'DD'))
data_C2C3_AGG_max1@meta.data[['group2']] <- data_C2C3_AGG_max1[[c('sample-mode_group2')]] %>% 
  mutate(sample = gsub('.*_', '', `sample-mode_group2`), .keep='unused') %>% unlist %>%
  factor(levels = c(paste0('CT', seq(3, 47, 4)), paste0('ZT', seq(3, 47, 4))))
data_C2C3_AGG_max1@meta.data[['sample_mode']] <- data_C2C3_AGG_max1[[c('sample', 'mode1')]] %>% 
  mutate(sample_mode <- paste(sample, mode1, sep='_'), .keep='unused') %>% unlist %>%
  factor(levels=c('C3_LD', 'C2_LD', 'C3_DD', 'C2_DD'))

avg_plotter_concat_mode(circ_genes, data=data_C2C3_AGG_max1, COUNTS=TRUE)


saveRDS(data_C2C3_AGG_max1, 'data_C2C3_AGG_max1')

####    ####    ####    ####    ####    ####    ####    ####    ####    ####   


# 2. START FROM FILES (from above code)
################################################################################

data_C2C3 <- readRDS('data_C2C3.rds')
data_C2C3_pcts <- CreateSeuratObject(counts = 1*(data_C2C3[["RNA"]]$counts !=0), meta.data = data_C2C3@meta.data)

genes_C2C3 <- rownames(data_C2C3)[rowSums(data_C2C3[["RNA"]]$counts)!=0]

genes_atLeast0.02pct <- AverageExpression(data_C2C3_pcts, features=genes_C2C3, assays='RNA', layer='counts')[[1]] %>%
  as.data.frame %>% 
  {rownames(.)[. > 0.02]}

data_C2C3_normGroup1 <- readRDS('data_C2C3_normGroup1.rds')
data_C2C3_normGroup2 <- readRDS('data_C2C3_normGroup2.rds')

data_C2C3_AGG <- readRDS('data_C2C3_AGG.rds')
#data_C2C3_AGG_CLRnorm <- readRDS('data_C2C3_AGG_CLRnorm.rds')

data_C2C3_AGG_max2 <- readRDS('data_C2C3_AGG_max2.rds')
data_C2C3_AGG_01_2d <- readRDS('data_C2C3_AGG_01_2d.rds')
data_C2C3_AGG_01_1d <- readRDS('data_C2C3_AGG_01_1d.rds')

added_C2C3 <- readRDS('added_C2C3.rds')
added_C2C3_AGG <- readRDS('added_C2C3_AGG.rds')
added_C2C3_pool3 <- readRDS('added_C2C3_pool3.rds')
added_C2C3_pool3_AGG <- readRDS('added_C2C3_pool3_AGG.rds')
added_C2C3_pool <- readRDS('added_C2C3_pool.rds')
added_C2C3_pool_AGG <- readRDS('added_C2C3_pool_AGG.rds')

reps_adding_pool <- readRDS('reps_adding_pool.rds')
avg_reps_adding_pool_AGG <- readRDS('avg_reps_adding_pool_AGG.rds')
avg_reps_adding_pool_AGG_01_1d <- readRDS('avg_reps_adding_pool_AGG_01_1d.rds')
avg_reps_adding_pool_AGG_1_3d <- readRDS('avg_reps_adding_pool_AGG_1_3d.rds')

avg_data_adding_pool_AGG_01_1d <- readRDS('avg_data_adding_pool_AGG_01_1d.rds')

mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "dmelanogaster_gene_ensembl", host = 'ensembl.org')

symbol_ensembl <- getBM(attributes = c("external_gene_name", "ensembl_gene_id"),
                        filters = "external_gene_name",
                        values = rownames(data_C2C3),
                        mart = mart) %>% 
  `colnames<-`(c('SYMBOL', 'ENSEMBL'))

merge(reconstituting_DFs_10$reconstituted, symbol_ensembl, by.x='Gene', by.y='SYMBOL', all.y=FALSE) %>% {split(., list(.$Sample, .$Mode))} %>% lapply(function(x) dim(x))
lapply(res10, function(res) dim(res$candidates))
lapply(res10, function(res) dim(subset(res$candidates, !is.na(ENSEMBL))))

symbol_ensembl[symbol_ensembl$SYMBOL %in% res10$C2_LD$candidates$SYMBOL,]

Cyclers <- merge(reconstituting_DFs_10$reconstituted, symbol_ensembl, by.x='Gene', by.y='SYMBOL', all.x=TRUE) %>%
  subset(select = c(Gene, ENSEMBL)) %>% 
  unique
subset(Cyclers, is.na(ENSEMBL))
{data.frame(c('312', 'FBgn0029514'),
           c('ABCA', 'FBgn0031170'),
           c('AdenoK', 'FBgn0036337'), 
           c('CG14812', 'FBgn0026090'), 
           c('CG1646', 'FBgn0039600'), 
           c('CG17556', 'FBgn0038462'), 
           c('CG17565', 'FBgn0038424'), 
           c('CG17802', 'FBgn0038549'), 
           c('CG2277', 'FBgn0035204'), 
           c('CG30069', 'FBgn0050069'), 
           c('CG31030', 'FBgn0051030'),
           c('CG32506', 'FBgn0052506'),
           c('CG32549', 'FBgn0052549'),
           c('CG3678', 'FBgn0038461'),
           c('CG42797', 'FBgn0261931'), 
           c('CG43778', 'FBgn0264308'),
           c('CG4538', 'FBgn0038745'),
           c('CG6145', 'FBgn0033853'),
           c('CG7028', 'FBgn0027587'),
           c('CG9253', 'FBgn0032919'),
           c('CG9384', 'FBgn0036446'),
           c('CG9425', 'FBgn0036451'),
           c('Dhap-at', 'FBgn0040212'),
           c('lid', 'FBgn0031759'),
           c('Mes4', 'FBgn0034726'),
           c('Nc73EF', 'FBgn0010352'),
           c('Nuak1', 'FBgn0262617'),
           c('prominin-like', 'FBgn0026189'),
           c('Rpb11', 'FBgn0032634'),
           c('SRm160', 'FBgn0036340'), 
           c('Thd1', 'FBgn0026869'),
           c('Updo', 'FBgn0033428')
)} %>% t %>% `colnames<-`(c('Gene', 'ENSEMBL_2')) %>%
  merge(Cyclers, all.y=TRUE) %>% 
  mutate(ENSEMBL = ifelse(is.na(ENSEMBL), ENSEMBL_2, ENSEMBL), 
         .keep='unused') -> Cyclers

passQC <- merge(subset(to_filter, second_pct>0.1), symbol_ensembl, by.x='Gene', by.y='SYMBOL', all.x=TRUE) %>%
  subset(select = c(Gene, ENSEMBL)) %>% 
  unique
subset(passQC, is.na(ENSEMBL))

{data.frame(c('312', 'FBgn0029514'),
           c('ABCA', 'FBgn0031170'),
           c('AdenoK', 'FBgn0036337'), 
           c('CG14812', 'FBgn0026090'), 
           c('CG1646', 'FBgn0039600'), 
           c('CG17556', 'FBgn0038462'), 
           c('CG17565', 'FBgn0038424'), 
           c('CG17802', 'FBgn0038549'), 
           c('CG2277', 'FBgn0035204'), 
           c('CG30069', 'FBgn0050069'), 
           c('CG31030', 'FBgn0051030'),
           c('CG32506', 'FBgn0052506'),
           c('CG32549', 'FBgn0052549'),
           c('CG3678', 'FBgn0038461'),
           c('CG42797', 'FBgn0261931'), 
           c('CG43778', 'FBgn0264308'),
           c('CG4538', 'FBgn0038745'),
           c('CG6145', 'FBgn0033853'),
           c('CG7028', 'FBgn0027587'),
           c('CG9253', 'FBgn0032919'),
           c('CG9384', 'FBgn0036446'),
           c('CG9425', 'FBgn0036451'),
           c('Dhap-at', 'FBgn0040212'),
           c('lid', 'FBgn0031759'),
           c('Mes4', 'FBgn0034726'),
           c('Nc73EF', 'FBgn0010352'),
           c('Nuak1', 'FBgn0262617'),
           c('prominin-like', 'FBgn0026189'),
           c('Rpb11', 'FBgn0032634'),
           c('SRm160', 'FBgn0036340'), 
           c('Thd1', 'FBgn0026869'),
           c('Updo', 'FBgn0033428'),
           c("A16","FBgn0028965"),
           c("ABCD","FBgn0039890"),
           c("ADPS","FBgn0033983"),
           c("Ady43A","FBgn0026602"),
           c("Akt1","FBgn0010379"),
           c("Arf102F","FBgn0013749"),
           c("Arf51F","FBgn0013750"),
           c("Arf79F","FBgn0010348"),
           c("Argk","FBgn0000116"),
           c("betaggt-II","FBgn0028970"),
           c("bip2","FBgn0026262"),
           c("Bsg25D","FBgn0000228"),
           c("Btk29A","FBgn0003502"),
           c("CG10082","FBgn0034644"),
           c("CG10103","FBgn0035715"),
           c("CG10188","FBgn0032796"),
           c("CG10672","FBgn0035588"),
           c("CG10960","FBgn0036316"),
           c("CG11070","FBgn0028467"),
           c("CG11137","FBgn0037199"),
           c("CG11319","FBgn0031835"),
           c("CG11399","FBgn0037021"),
           c("CG11486","FBgn0035397"),
           c("CG11753","FBgn0037603"),
           c("CG11781","FBgn0039259"),
           c("CG11873","FBgn0039633"),
           c("CG12004","FBgn0035236"),
           c("CG12163","FBgn0260462"),
           c("CG12688","FBgn0029707"),
           c("CG12913","FBgn0033500"),
           c("CG13366","FBgn0025633"),
           c("CG14184","FBgn0036932"),
           c("CG14407","FBgn0030584"),
           c("CG14806","FBgn0029594"),
           c("CG14977","FBgn0035469"),
           c("CG15168","FBgn0032732"),
           c("CG15309","FBgn0030183"),
           c("CG16791","FBgn0038881"),
           c("CG17454","FBgn0039977"),
           c("CG17716","FBgn0000633"),
           c("CG1814","FBgn0033426"),
           c("CG18815","FBgn0042138"),
           c("CG1882","FBgn0033226"),
           c("CG30424","FBgn0050424"),
           c("CG31324","FBgn0051324"),
           c("CG31357","FBgn0051357"),
           c("CG31694","FBgn0051694"),
           c("CG32039","FBgn0052039"),
           c("CG32441","FBgn0052441"),
           c("CG32521","FBgn0052521"),
           c("CG32581","FBgn0052581"),
           c("CG32850","FBgn0052850"),
           c("CG33977","FBgn0053977"),
           c("CG34351","FBgn0085380"),
           c("CG34401","FBgn0085430"),
           c("CG3530","FBgn0028497"),
           c("CG3542","FBgn0031492"),
           c("CG42668","FBgn0261550"),
           c("CG43222","FBgn0262858"),
           c("CG4552","FBgn0031304"),
           c("CG4587","FBgn0028863"),
           c("CG4603","FBgn0035593"),
           c("CG46338","FBgn0288856"),
           c("CG4825","FBgn0037010"),
           c("CG4908","FBgn0032195"),
           c("CG5027","FBgn0036579"),
           c("CG5028","FBgn0039358"),
           c("CG5037","FBgn0032222"),
           c("CG5044","FBgn0038326"),
           c("CG5059","FBgn0037007"),
           c("CG5110","FBgn0032642"),
           c("CG5522","FBgn0034158"),
           c("CG5830","FBgn0036556"),
           c("CG5903","FBgn0038400"),
           c("CG6024","FBgn0036202"),
           c("CG6123","FBgn0030913"),
           c("CG6218","FBgn0038321"),
           c("CG6379","FBgn0029693"),
           c("CG6512","FBgn0036702"),
           c("CG6617","FBgn0030944"),
           c("CG6767","FBgn0036030"),
           c("CG7154","FBgn0031947"),
           c("CG7324","FBgn0037074"),
           c("CG7414","FBgn0037135"),
           c("CG7946","FBgn0039743"),
           c("CG7956","FBgn0038890"),
           c("CG7971","FBgn0035253"),
           c("CG8009","FBgn0036090"),
           c("CG8108","FBgn0027567"),
           c("CG8177","FBgn0036043"),
           c("CG8397","FBgn0034066"),
           c("CG8485","FBgn0033915"),
           c("CG8878","FBgn0027504"),
           c("CG8974","FBgn0030693"),
           c("CG9065","FBgn0030610"),
           c("CG9601","FBgn0037578"),
           c("CG9911","FBgn0030734"),
           c("CG9934","FBgn0032467"),
           c("CG9941","FBgn0030514"),
           c("CtsB1","FBgn0030521"),
           c("Dak1","FBgn0028833"),
           c("Ect4","FBgn0262579"),
           c("Eip71CD","FBgn0000565"),
           c("Gtp-bp","FBgn0010391"),
           c("l(1)G0007","FBgn0026713"),
           c("l(1)G0193","FBgn0027280"),
           c("l(2)37Cc","FBgn0002031"),
           c("l(3)72Ab","FBgn0263599"),
           c("metl","FBgn0035247"),
           c("mtRNApol","FBgn0261938"),
           c("olf186-F","FBgn0041585"),
           c("PI4KIIIalpha","FBgn0267350"),
           c("Pmp70","FBgn0031069"),
           c("Rpb10","FBgn0039218"),
           c("Rpb12","FBgn0262954"),
           c("Rpb8","FBgn0037121"),
           c("RpII18","FBgn0003275"),
           c("RpII215","FBgn0003277"),
           c("Sep1","FBgn0011710"),
           c("Sep2","FBgn0014029"),
           c("Sep4","FBgn0259923"),
           c("smt3","FBgn0264922"),
           c("snRNA:7SK","FBgn0065099"),
           c("Sra-1","FBgn0038320"),
           c("TMS1","FBgn0028399"),
           c("Ufd1-like","FBgn0036136"),
           c("Yippee","FBgn0026749")
)} %>% t %>% `colnames<-`(c('Gene', 'ENSEMBL_2')) %>%
  merge(passQC, all.y=TRUE) %>% 
  mutate(ENSEMBL = ifelse(is.na(ENSEMBL), ENSEMBL_2, ENSEMBL), 
         .keep='unused') -> passQC

####    ####    ####    ####    ####    ####    ####    ####    ####    ####   


# 3. FUNCTIONS TO INITIALIZE
################################################################################
# Split a list into 3 parts
split_list = function(my_list, N) {split(my_list, ceiling(seq_along(my_list) / ceiling(length(my_list) / N)))}


# FILTERING/PREPARING:
PrctCellExpringGene <- function(object, genes, group.by = "all"){
  if(group.by == "all"){
    prct = unlist(lapply(genes,calc_helper, object=object))
    result = data.frame(Markers = genes, Cell_proportion = prct)
    return(result)
  }
  
  else{        
    list = SplitObject(object, group.by)
    factors = names(list)
    
    results = lapply(list, PrctCellExpringGene, genes=genes)
    for(i in 1:length(factors)){
      results[[i]]$Feature = factors[i]
    }
    combined = do.call("rbind", results)
    return(combined)
  }
}
calc_helper <- function(object,genes){
  counts = object[['RNA']]@counts
  ncells = ncol(counts)
  sum(counts[genes,]>0)/ncells
}

PrctCellExpringGene_presplit = function(object_list, genes) {
  # make object_list by splitting: SplitObject(data_C2C3, 'sample') %>% lapply(SplitObject, split.by='mode1') %>% unlist()
  prct = unlist(lapply(genes,calc_helper, object=object_list))
  result = data.frame(Markers = genes, Cell_proportion = prct)
}

append_kruskal_scores_from_above = function(a_df, COUNTS=TRUE, NO_ZEROS=FALSE, SEB_FILT=FALSE, zt_name='group1', data=data_C2C3) {
  
  name = paste("KW_scores", ifelse(COUNTS, "_cts", "_norm"),  sep='')
  a_df = a_df[which(a_df$Gene!=' '),]
  
  apply(a_df, MARGIN = 1, FUN = function(x, data=data, zt_name=zt_name, COUNTS=COUNTS) {
    Kruskal_one_gene_sample(Sample = x['Sample'], Gene = x['Gene'], zt_name=zt_name, COUNTS=COUNTS, data=data)
  }, data=data, zt_name=zt_name, COUNTS=COUNTS) %>% unlist %>% 
    {data.frame(Gene=a_df['Gene'], Sample=a_df['Sample'], KW_scores=.)} %>%
    rename_with(~name, "KW_scores") %>% return
}

Kruskal_one_gene_sample = function(Sample, Gene, COUNTS=TRUE, data=data_C2C3, zt_name='group1'){
  
  FetchData(data, vars=c("sample", zt_name, Gene), layer = ifelse(COUNTS,"counts","data")) %>%
    subset(sample == Sample) %>% 
    pivot_longer(cols = -c(1,2), names_to = "Gene", values_to = "expr") %>% arrange(desc(expr)) %>%
    `colnames<-`(c('sample', 'group1', 'Gene', 'expr')) %>%
    {if (!COUNTS) mutate(., expr=expm1(expr)) else .} %>% 
    {if (length(which(!.$expr==0))<=2) {
      NaN
    } else {
      kruskal.test(expr ~ group1, data=.) %>% 
        .$p.value %>% return()
    }} %>% return()
}

filterer = function(to_filter, expression_vec) {
  expression_vec = paste0('to_filter$', expression_vec)
  for (expr in expression_vec) {
    to_filter = na.omit(to_filter[eval(parse(text = expr)), ])
  }
  return(to_filter)
}


# PLOTTING:

subplotter_concat = function(names, data=data_C2C3, samples=levels(data$sample), 
                             COUNTS=FALSE, N_days = 2, rib_alpha=0.4, point_size=2, 
                             BLACKOUT=FALSE, NO_X_TEXT=TRUE, NO_GRIDLINES=TRUE, ANNOTATE=FALSE,
                             func = 'avg_plotter_concat', more_args_list=NULL,
                             n_plots=ceiling(length(names)/16), n_row=1) {

  len = length(names)
  more_args_list = c(more_args_list, data = data, COUNTS=COUNTS, N_days=N_days, rib_alpha=rib_alpha, point_size=point_size,
                     ANNOTATE=ANNOTATE, BLACKOUT=BLACKOUT, NO_GRIDLINES=NO_GRIDLINES, NO_X_TEXT=NO_X_TEXT) %>% as.list
  
  lapply(1:n_plots, function(n, names, len, samples, func, more_args_list){
    list(gene = names[(1+floor((n-1)/n_plots*len)):floor(n/n_plots*len)], smpl = samples) %>%
      {if(length(more_args_list)) append(., unlist(more_args_list)) else .} %>% 
      {do.call(func, .)}
  }, names=unlist(names), len=len, samples=samples, func=func, more_args_list=more_args_list) %>%
    {grid.arrange(grobs = ., nrow=n_row)}
}


avg_plotter_concat_mode = function(gene, smpl=levels(data$sample), dgrp_colname='group2', 
                                   COUNTS=FALSE, quicknorm_object=NULL, include_SE=TRUE,include_SE_bars=FALSE, ZERO_Y=TRUE, 
                                   data=data_C2C3, point_size=2, rib_alpha=1, NO_X_TEXT=TRUE, BLACKOUT=FALSE, N_days = 2, 
                                   FACET_2D=FALSE, SWAP_FACET=FALSE, NO_GRIDLINES=TRUE, GENEMODE=FALSE) {
  if (!GENEMODE) {
    plot_grid(nrow=1, 
              subset(data, mode1 == 'LD') %>%
                avg_plotter_concat(gene=gene, smpl=smpl, dgrp_colname=dgrp_colname,
                                   COUNTS=COUNTS, quicknorm_object=quicknorm_object,
                                   include_SE=include_SE, include_SE_bars=include_SE_bars,
                                   ZERO_Y=ZERO_Y, point_size=point_size, rib_alpha=rib_alpha,
                                   NO_X_TEXT=NO_X_TEXT, BLACKOUT=BLACKOUT, N_days=N_days, FACET_2D=FACET_2D, SWAP_FACET=SWAP_FACET, NO_GRIDLINES=NO_GRIDLINES), 
              subset(data, mode1 == 'DD') %>%
                avg_plotter_concat(gene=gene, smpl=smpl, dgrp_colname=dgrp_colname,
                                   COUNTS=COUNTS, quicknorm_object=quicknorm_object,
                                   include_SE=include_SE, include_SE_bars=include_SE_bars,
                                   ZERO_Y=ZERO_Y, point_size=point_size, rib_alpha=rib_alpha,
                                   NO_X_TEXT=NO_X_TEXT, BLACKOUT=BLACKOUT, N_days=N_days, FACET_2D=FACET_2D, SWAP_FACET=SWAP_FACET, NO_GRIDLINES=NO_GRIDLINES)) %>% print
  } else {
    plot_grid(nrow=2, plotlist = lapply(levels(data$sample_mode), function(s_m) {
      avg_plotter_concat(data = subset(data, sample_mode == s_m), 
                         gene=gene, smpl=smpl, dgrp_colname=dgrp_colname,
                         COUNTS=COUNTS, quicknorm_object=quicknorm_object,
                         include_SE=include_SE, include_SE_bars=include_SE_bars,
                         ZERO_Y=ZERO_Y, point_size=point_size, rib_alpha=rib_alpha,
                         NO_X_TEXT=NO_X_TEXT, BLACKOUT=BLACKOUT, N_days=N_days, FACET_2D=FACET_2D, SWAP_FACET=TRUE, NO_GRIDLINES=NO_GRIDLINES)
    }))
  }
  
  #SplitObject(data_C2C3, 'mode1') %>%
  #  lapply(function(df) avg_plotter_concat(data=df,
  #                                         gene=gene, smpl=smpl, dgrp_colname=dgrp_colname,
  #                                         COUNTS=COUNTS, quicknorm_object=quicknorm_object,
  #                                         include_SE=include_SE, include_SE_bars=include_SE_bars,
  #                                         ZERO_Y=ZERO_Y, point_size=point_size, rib_alpha=rib_alpha,
  #                                         NO_X_TEXT=NO_X_TEXT, N_days=N_days)) %>%
  #  plot_grid(nrow = 1, plotlist = .)
}



avg_plotter_concat = function(gene, smpl=unique(data$sample), dgrp_colname='group2', 
                              COUNTS=FALSE, quicknorm_object=NULL, include_SE=TRUE,include_SE_bars=FALSE, ZERO_Y=TRUE, 
                              data=data_C2C3, point_size=2, line_size=0.5, rib_alpha=1, NO_X_TEXT=FALSE, NO_Y_TEXT=FALSE, BLACKOUT=FALSE, N_days = NULL, 
                              ANNOTATE=FALSE, FACET_2D=FALSE, SWAP_FACET=FALSE, NO_TAB_X=FALSE, NO_TAB_Y=FALSE, NO_GRIDLINES=FALSE, LEGEND=FALSE, BRIGHTEN_DD_DAY=TRUE, 
                              strip_x_manual=NULL, strip_y_manual=NULL) {
  
  Mode=unique(data$mode1)
  LD = Mode=='LD'
  
  p = FetchData(data, vars = c("sample", dgrp_colname, gene), layer = ifelse(COUNTS,"counts","data")) %>% 
    subset(sample %in% smpl) %>%
    pivot_longer(cols = -c(1,2), names_to = "Gene", values_to = "expr") %>%
    `colnames<-`(c('sample', 'ZT_concat', 'Gene', 'expr')) %>% 
    mutate(ZT_concat = gsub('.*T', '', ZT_concat) %>% as.integer) %>%
    {if (!is.null(quicknorm_object)) {
      merge(., cbind(quicknorm_object %>% dplyr::select(Gene, sample), 
                     ZTmax = quicknorm_object %>% subset(select = -c(Gene, sample)) %>% as.matrix %>% 
                       rowMax) %>% 
              mutate(sample=sample, .keep='unused')) %>%
        mutate(expr = expr/ZTmax, .keep='unused')
    } else .} %>% 
    mutate(Gene = factor(Gene, levels = gene %>% unique), 
           sample = factor(sample, levels = smpl %>% unique)) %>% 
    {if (!COUNTS) mutate(., expr = expm1(expr)) else .} %>% 
    ggplot(aes(x = ZT_concat, y = expr, group = Gene, col = Gene, fill = Gene)) +
    stat_summary(geom = "line", fun = "mean", linewidth=line_size) +
    {if (point_size) stat_summary(geom = "point", fun = "mean", size=point_size)} +
    #{if (include_SE) {
    #  stat_summary(fun = mean, geom = "ribbon", alpha=rib_alpha, color = NA,
    #               fun.max = function(x) mean(x) + sd(x) / sqrt(length(x)),
    #               fun.min = function(x) mean(x) - sd(x) / sqrt(length(x))) }} +
    #{if (include_SE_bars) stat_summary(fun.data = "mean_cl_boot")} +
    {if (ZERO_Y) expand_limits(y=0)} +
    {if (FACET_2D) facet_grid(. ~ factor(sample, levels = smpl), scales = "free_y")
      else if (SWAP_FACET) facet_grid(factor(sample, levels = smpl) ~ Gene, scales = "free_y")
      else facet_grid(Gene ~ factor(sample, levels = smpl), scales = "free_y")} +
    {if (SWAP_FACET) theme(strip.text.x = element_text(face = "italic"))
      else if (!FACET_2D) theme(strip.text.y = element_text(face = "italic"))} +
    {if (NO_GRIDLINES) theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())} +
    theme(axis.text.x = element_text(size = 4, angle=45, hjust=0.5, vjust=1), axis.title = element_blank()) +
    {if (!LEGEND) theme(legend.position = "none")} +
    {if (NO_X_TEXT) theme(axis.text.x = element_blank(), axis.ticks.x=element_blank())} +
    {if (NO_Y_TEXT) theme(axis.text.y = element_blank(), axis.ticks.y=element_blank())} +
    {if (BLACKOUT) {scale_color_manual(values=rep('black', length(unique(gene))))}} +
    {if (BLACKOUT) {scale_fill_manual(values=rep('black', length(unique(gene))))}} +
    {if (NO_TAB_X) theme(strip.text.x = element_blank())} +
    {if (NO_TAB_Y) theme(strip.text.y = element_blank())} +
    {if (!is.null(strip_x_manual)) theme(strip.text.x = element_text(margin = margin(strip_x_manual,0,strip_x_manual,0, "in")))} +
    {if (!is.null(strip_y_manual)) theme(strip.text.y = element_text(margin = margin(0,strip_y_manual,0,strip_y_manual, "in")))}
  
  
  if (!is.null(N_days)) {
    for (day in 1:N_days) {
      if (LD) {
      p = p + 
        geom_rect(data=data.frame(xmin=0+(24*(day-1)), xmax=12+(24*(day-1)), ymin=-Inf, ymax=Inf), 
                  aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="yellow", colour='black', linewidth=0.1, alpha=0.3, inherit.aes = FALSE) +
        geom_rect(data=data.frame(xmin=12+(24*(day-1)), xmax=24+(24*(day-1)), ymin=-Inf, ymax=Inf), 
                  aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey35", colour='black', linewidth=0.1, alpha=0.3, inherit.aes = FALSE)
      } else {
      p = p + 
        geom_rect(data=data.frame(xmin=0+(24*(day-1)), xmax=12+(24*(day-1)), ymin=-Inf, ymax=Inf), 
                  aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill=ifelse(BRIGHTEN_DD_DAY, "grey65", "grey35"), colour='black', linetype='dashed', linewidth=0.1, alpha=0.3, inherit.aes = FALSE) +
        geom_rect(data=data.frame(xmin=12+(24*(day-1)), xmax=24+(24*(day-1)), ymin=-Inf, ymax=Inf), 
                  aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey35", colour='black', linetype='dashed', linewidth=0.1, alpha=0.3, inherit.aes = FALSE)
      }
    }
    p = p +
      scale_x_continuous(breaks = data[[dgrp_colname]][[1]] %>% unique %>% gsub('.*T', '', .) %>% as.integer %>% sort, limits = c(0, N_days*24)) +
      scale_x_continuous(expand = c(0, 0))
    
  } else  p = p +
      scale_x_continuous(breaks = data[[dgrp_colname]][[1]] %>% unique %>% gsub('.*T', '', .) %>% as.integer %>% sort, limits = c(NA, NA)) +
      coord_cartesian(xlim=c(1,23))
  
  if (ANNOTATE) {
    
    p = lapply(gene, plot_annotation_table, sample = smpl, mode = Mode) %>% 
      plot_grid(ncol = 1, plotlist = .) %>%
      plot_grid(nrow=1, p, ., rel_widths = c(1.5, 1))
  }
  
  return(p)
}

add_rectangles_to_plots = function(p, N_days=2, LD=TRUE, BRIGHTEN_DD_DAY=TRUE, Alpha=0.3) {
  if (!is.null(N_days)) {
    for (day in 1:N_days) {
      if (LD) {
        p = p + 
          geom_rect(data=data.frame(xmin=0+(24*(day-1)), xmax=12+(24*(day-1)), ymin=-Inf, ymax=Inf), 
                    aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="yellow", colour='black', linewidth=0.1, alpha=Alpha, inherit.aes = FALSE) +
          geom_rect(data=data.frame(xmin=12+(24*(day-1)), xmax=24+(24*(day-1)), ymin=-Inf, ymax=Inf), 
                    aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey35", colour='black', linewidth=0.1, alpha=Alpha, inherit.aes = FALSE)
      } else {
        p = p + 
          geom_rect(data=data.frame(xmin=0+(24*(day-1)), xmax=12+(24*(day-1)), ymin=-Inf, ymax=Inf), 
                    aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill=ifelse(BRIGHTEN_DD_DAY, "grey65", "grey35"), colour='black', linetype='dashed', linewidth=0.1, alpha=Alpha, inherit.aes = FALSE) +
          geom_rect(data=data.frame(xmin=12+(24*(day-1)), xmax=24+(24*(day-1)), ymin=-Inf, ymax=Inf), 
                    aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey35", colour='black', linetype='dashed', linewidth=0.1, alpha=Alpha, inherit.aes = FALSE)
      }
    }
  }
  p = p +
    scale_x_continuous(breaks = seq(3, 47, 4), limits = c(0, N_days*24)) +
    scale_x_continuous(expand = c(0, 0))
  return(p)
}


avg_plotter_repDot = function(gene, data=data_C2C3, COUNTS=FALSE, rib_alpha=0.4, 
                              ZERO_Y=TRUE, NO_X_TEXT=FALSE, NO_DD=FALSE, NO_LD=FALSE, sample_remove=NULL, NO_GRIDLINES=FALSE) {
  df = AverageExpression(object=data, features = gene, group.by = c('group2', 'sample_mode'), layer = ifelse(COUNTS,"counts","data")) %>%
    .[['RNA']] %>% as.data.frame %>% 
    {mutate(., Gene = rownames(.))} %>%
    pivot_longer(cols = -Gene, names_to = 'condition') %>%
    mutate(Sample_Mode = gsub('.*_', '', condition) %>% gsub('-', '_', .) %>%
             factor(levels = c('C3_LD', 'C2_LD', 'C3_DD', 'C2_DD')),
           timepoint = gsub('.*T', '', condition) %>% gsub('_.*', '', .) %>% as.integer, .keep='unused') %>%
    mutate(rep = ifelse(timepoint > 24, 2, 1), 
           timepoint = timepoint %% 24, 
           isDD = ifelse(gsub('.*_', '', Sample_Mode)=='DD', TRUE, FALSE), 
           Gene = factor(Gene, levels = gene %>% unique)) %>%
    subset(gsub('_.*', '', Sample_Mode) %!in% sample_remove)
  
  ldPlot = subset(df, isDD==FALSE) %>%
    {ggplot(., aes(x=timepoint, y=value, group=rep, col = Gene, fill = Gene)) +
    geom_point() +
    geom_line() +
    #geom_ribbon(mapping=aes(x=timepoint, y=value, fill=Gene), inherit.aes = FALSE,
    #            alpha = 0.4, stat = "summary", fun.max = max, fun.min = min) +
    geom_ribbon(data=pivot_wider(., names_from = 'rep', values_from = 'value'),
                mapping=aes(x=timepoint, ymin=`1`, ymax=`2`, fill=Gene), inherit.aes=FALSE, alpha=rib_alpha) +
    {if (ZERO_Y) expand_limits(y=0)} +
    facet_grid(Gene ~ Sample_Mode, scales = 'free_y') +
    {if (NO_X_TEXT) theme(axis.text.x = element_blank(), axis.ticks.x=element_blank())} +
    theme(axis.text.x = element_text(size = 6), axis.title = element_blank(), legend.position = "none") +
    {if (NO_GRIDLINES) theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())} +
    scale_x_continuous(breaks = seq(3, 23, 4)) +
    coord_cartesian(xlim=c(1,23)) +
    geom_rect(data=data.frame(xmin=0, xmax=12, ymin=-Inf, ymax=Inf), 
              aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="yellow", alpha=0.2, inherit.aes = FALSE) +
    geom_rect(data=data.frame(xmin=12, xmax=24, ymin=-Inf, ymax=Inf), 
              aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey", alpha=0.2, inherit.aes = FALSE)}
  
  ddPlot = subset(df, isDD==TRUE) %>%
    {ggplot(., aes(x=timepoint, y=value, group=rep, col = Gene, fill = Gene)) +
    geom_point() +
    geom_line() +
    #geom_ribbon(mapping=aes(x=timepoint, y=value, fill=Gene), inherit.aes = FALSE,
    #            alpha = 0.4, stat = "summary", fun.max = max, fun.min = min) +
    geom_ribbon(data=pivot_wider(., names_from = 'rep', values_from = 'value'),
                mapping=aes(x=timepoint, ymin=`1`, ymax=`2`, fill=Gene), inherit.aes=FALSE, alpha=rib_alpha) +
    {if (ZERO_Y) expand_limits(y=0)} +
    facet_grid(Gene ~ Sample_Mode, scales = 'free_y') +
    {if (NO_X_TEXT) theme(axis.text.x = element_blank(), axis.ticks.x=element_blank())} +
    theme(axis.text.x = element_text(size = 6), axis.title = element_blank(), legend.position = "none") +
    {if (NO_GRIDLINES) theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())} +
    scale_x_continuous(breaks = seq(3, 23, 4)) +
    coord_cartesian(xlim=c(1,23)) +
    geom_rect(data=data.frame(xmin=0, xmax=12, ymin=-Inf, ymax=Inf), 
              aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey", alpha=0.2, inherit.aes = FALSE) +
    geom_rect(data=data.frame(xmin=12, xmax=24, ymin=-Inf, ymax=Inf), 
              aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey", alpha=0.2, inherit.aes = FALSE)}
  
  if(!NO_DD*NO_LD) {
    if(NO_DD) ldPlot
    else if(NO_LD) ddPlot
    else plot_grid(nrow=1, ldPlot, ddPlot)
  }
}

plot_annotation_table = function(gene, sample, mode, to_filt=to_filter) {
  newTable = to_filt %>% 
    subset(Sample == sample & Mode == mode & Gene == gene, 
           select = c(Gene, Cell_proportion, most_pct, second_pct, ratio, ampl, meta2d_pvalue, meta2d_pvalue_AGG_norm)) %>% 
    mutate_at(vars(-Gene), as.numeric) %>%
    mutate_at(vars(-Gene), round, digits=3) %>% 
    t %>% as.data.frame %>% 
    {mutate(., Gene = rownames(.))} %>% 
    {`colnames<-`(., .[1,])} %>% .[-1, c(2,1)] %>% 
    flextable()
  
  newTable = ggplot() + theme_void() +
    annotation_custom(gen_grob(newTable),
                      xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)
}



plotter_geneCombine_modes = function(Genes_, data_ = data_C2C3_AGG_01_1d, dgrp_colname_='group2', Sample_=c('C3', 'C2'), 
                                     include_SE_=TRUE, include_SE_bars_=FALSE, COUNTS_=TRUE, ZERO_Y_=TRUE, NO_X_TEXT_=FALSE, 
                                     rib_alpha_=0.4, point_size_=0, N_days_=2, ALIGN_=TRUE, offset_=-0.5, NO_GRIDLINES_=TRUE,
                                     PRINT_MEAN_=TRUE, PRINT_GENES_=FALSE, LEGEND_=FALSE) {
  
  plot_grid(nrow=1, 
            plotter_geneCombine(genes=Genes_, data=data_, Mode='LD', Sample=Sample_, dgrp_colname=dgrp_colname_,
                                include_SE=include_SE_, include_SE_bars=include_SE_bars_, 
                                COUNTS=COUNTS_, ZERO_Y=ZERO_Y_, NO_X_TEXT=NO_X_TEXT_, 
                                rib_alpha=rib_alpha_, point_size=point_size_, N_days=N_days_,
                                ALIGN=ALIGN_, offset=offset_, NO_GRIDLINES=NO_GRIDLINES_,
                                PRINT_MEAN=PRINT_MEAN_, PRINT_GENES=PRINT_GENES_, LEGEND=LEGEND_), 
            plotter_geneCombine(genes=Genes_, data=data_, Mode='DD', Sample=Sample_, dgrp_colname=dgrp_colname_,
                                include_SE=include_SE_, include_SE_bars=include_SE_bars_, 
                                COUNTS=COUNTS_, ZERO_Y=ZERO_Y_, NO_X_TEXT=NO_X_TEXT_, 
                                rib_alpha=rib_alpha_, point_size=point_size_, N_days=N_days_,
                                ALIGN=ALIGN_, offset=offset_, NO_GRIDLINES=NO_GRIDLINES_, 
                                PRINT_MEAN=PRINT_MEAN_, PRINT_GENES=PRINT_GENES_, LEGEND=LEGEND_))
}



plotter_geneCombine = function(genes, data = data_C2C3_AGG_01_1d, Mode='LD', Sample=c('C3', 'C2'), dgrp_colname='group2', 
                               include_SE=TRUE, include_SE_bars=FALSE, COUNTS=TRUE, ZERO_Y=TRUE, NO_X_TEXT=FALSE, NO_Y_TEXT=FALSE,
                               rib_alpha=0.75, gene_alpha=1, point_size=0, mean_linesize=0.5, gene_linesize=0.5,
                               errorWidth=0, errorThick=1, ErrorFun="mean_se",
                               N_days=2, ALIGN=FALSE, offset=-0.5, NO_GRIDLINES=FALSE, 
                               PRINT_MEAN=TRUE, PRINT_GENES=FALSE, LEGEND=FALSE, Y_LIM_10_OFFSET=TRUE, PRISM=FALSE, CLASSIC=FALSE, GENECOLOR=TRUE) {
  
  Sample_Mode = do.call(paste, c(expand.grid(Sample, Mode), sep='_'))
  LD = all(Mode=='LD')
  
  p = FetchData(data, vars = c('sample_mode', dgrp_colname, genes), layer=ifelse(COUNTS,'counts','data')) %>%
    subset(sample_mode %in% Sample_Mode) %>%
    pivot_longer(cols = -c(1,2), names_to = "Gene", values_to = "expr") %>%
    `colnames<-`(c('sample_mode', 'ZT_concat', 'Gene', 'expr')) %>% 
    mutate(ZT_concat = gsub('.*T', '', ZT_concat) %>% as.integer) %>%
    mutate(Gene = factor(Gene, levels = genes %>% unique), 
           sample_mode = factor(sample_mode, levels = Sample_Mode)) %>%
    mutate(Sample = factor(gsub('_.*', '', sample_mode), levels = c('C3', 'C2'))) %>% 
    {if (!COUNTS) mutate(., expr = expm1(expr)) else .} %>% 
    {if (ALIGN) {
      {split(., list(.$Gene, .$sample_mode))} %>%
        lapply(function(df) {
          mutate(df, expr = (expr - min(expr))) %>% 
            mutate(expr = expr/max(expr))
        }) %>%
        do.call(what='rbind')
    } else .} %>% 
    mutate(expr = expr + offset) %>%
    ggplot(aes(x = ZT_concat, y = expr))
  
  if (is.null(N_days)) N_days = max(as.numeric(gsub('[^0-9]', '', data[[dgrp_colname]])))%/%24 + 1
  
  p = allRect(p, N_days, LD)
  
  p = p +
    {if (PRINT_GENES) {
      if(GENECOLOR) geom_line(aes(color=Gene), alpha=gene_alpha, linewidth=gene_linesize)
      else geom_line(aes(group=Gene), alpha=gene_alpha, linewidth=gene_linesize)}} +
    {if (point_size) stat_summary(geom = "point", fun = "mean", size=point_size)} +
    {if (include_SE) {
      stat_summary(fun = mean, geom = "ribbon", alpha=rib_alpha, color = NA, 
                   fun.max = function(x) mean(x) + sd(x) / sqrt(length(x)),
                   fun.min = function(x) mean(x) - sd(x) / sqrt(length(x))) }
    } +
    {if (include_SE_bars) stat_summary(fun.data = ErrorFun, geom = "errorbar", width = errorWidth, size=errorThick)} +
    scale_alpha(range = c(0, 1)) +
    {if (PRINT_MEAN) stat_summary(geom = "line", fun = "mean", linewidth=mean_linesize)} +
    theme_linedraw() + theme(axis.text=element_text(colour="black")) +
    {if (PRISM) theme_prism()} +
    {if (CLASSIC) theme_classic()} +
    {if (ZERO_Y) expand_limits(y=0)} +
    {if (Y_LIM_10_OFFSET) ylim(0+offset,1+offset)} +
    facet_grid(. ~ Sample, scales = "free_y") +
    {if (NO_GRIDLINES) theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())} +
    theme(axis.text.x = element_text(hjust=0.5, vjust=1), axis.title = element_blank(), legend.text = element_text(face = "italic")) +
    {if (!LEGEND) theme(legend.position = 'none')} +
    {if (NO_X_TEXT) theme(axis.text.x = element_blank(), axis.ticks.x=element_blank())} +
    {if (NO_Y_TEXT) theme(axis.text.y = element_blank(), axis.ticks.y=element_blank())} + 
    scale_x_continuous(expand = c(0, 0))
  
  #if (!is.null(N_days)) {
  #  for (day in 1:N_days) {
  #    if (LD) {
  #      p = p + 
  #        geom_rect(data=data.frame(xmin=0+(24*(day-1)), xmax=12+(24*(day-1)), ymin=-Inf, ymax=Inf), 
  #                  aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="yellow", colour='black', linewidth=0.1, alpha=0.3, inherit.aes = FALSE) +
  #        geom_rect(data=data.frame(xmin=12+(24*(day-1)), xmax=24+(24*(day-1)), ymin=-Inf, ymax=Inf), 
  #                  aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey35", colour='black', linewidth=0.1, alpha=0.3, inherit.aes = FALSE)
  #    } else {
  #      p = p + 
  #        geom_rect(data=data.frame(xmin=0+(24*(day-1)), xmax=12+(24*(day-1)), ymin=-Inf, ymax=Inf), 
  #                  aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill=ifelse(1, "grey65", "grey35"), colour='black', linetype='dashed', linewidth=0.1, alpha=0.3, inherit.aes = FALSE) +
  #        geom_rect(data=data.frame(xmin=12+(24*(day-1)), xmax=24+(24*(day-1)), ymin=-Inf, ymax=Inf), 
  #                  aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey35", colour='black', linetype='dashed', linewidth=0.1, alpha=0.3, inherit.aes = FALSE)
  #    }
  #  }
  #  
  #  p = p +
  #    scale_x_continuous(breaks = data[[dgrp_colname]][[1]] %>% unique %>% gsub('.*T', '', .) %>% as.integer %>% sort, limits = c(0, N_days*24)) +
  #    scale_x_continuous(expand = c(0, 0))
  #  
  #} else {
  #  p = p +
  #    scale_x_continuous(breaks = data[[dgrp_colname]][[1]] %>% unique %>% gsub('.*T', '', .) %>% as.integer %>% sort, limits = c(NA, NA)) +
  #    coord_cartesian(xlim=c(1,23))
  #}
  
  return(p)
}


repplotter_ofSummary = function(DF, Gene0=NULL, Sample0=NULL, Mode0=NULL, Timecol='Time', Repcol='Rep', Valcol='value', 
                                include_SE=TRUE, include_SE_bars=FALSE, ErrorFun="mean_se",
                                PRINT_MEAN=TRUE, mean_linesize=0.5, mean_point_size=0, rib_alpha=0.4, errorWidth=1, errorThick=1,
                                PRISM=FALSE, CLASSIC=FALSE, ZERO_Y=TRUE, Y_LIM_10_OFFSET=FALSE, 
                                NO_GRIDLINES=TRUE, LEGEND=FALSE, NO_X_TEXT=TRUE, NO_Y_TEXT=FALSE,
                                FACET_STRIPS=FALSE, hard_ys=NULL, expand_ys=c(0.05, 0.05)) {
  if (!is.null(Gene0)) {DF = subset(DF, Gene==Gene0)}
  if (!is.null(Sample0)) {DF = subset(DF, Sample==Sample0)}
  if (!is.null(Mode0)) {DF = subset(DF, Mode==Mode0)}
  
  LD = all(DF$Mode=='LD')
  
  DF = DF[c(Timecol, Repcol, Valcol, 'Sample', 'Mode')] %>%
    `colnames<-`(c('Time', 'Rep', 'value', 'Sample', 'Mode')) %>%
    mutate(across(c(value, Time), as.numeric))
  p = DF %>%
    ggplot(aes(x = Time, y = value))
  
  N_days = max(as.numeric(gsub('[^0-9]', '', DF$Time)))%/%24 + 1
  
  p = allRect(p, N_days, LD)
  p = p +
    {if (mean_point_size) stat_summary(geom = "point", fun = "mean", size=mean_point_size)} +
    {if (PRINT_MEAN) stat_summary(geom = "line", fun = "mean", linewidth=mean_linesize)} +
    {if (include_SE) {
      stat_summary(fun = mean, geom = "ribbon", alpha=rib_alpha, color = NA, 
                   fun.max = function(x) mean(x) + sd(x) / sqrt(length(x)),
                   fun.min = function(x) mean(x) - sd(x) / sqrt(length(x))) }
    } +
    {if (include_SE_bars) stat_summary(fun.data = ErrorFun, geom = "errorbar", width = errorWidth, size=errorThick)} +
    scale_alpha(range = c(0, 1)) +
    theme_linedraw() + theme(axis.text=element_text(colour="black")) +
    {if (PRISM) theme_prism()} +
    {if (CLASSIC) theme_classic()} +
    {if (ZERO_Y) expand_limits(y=0)} +
    {if (Y_LIM_10_OFFSET) ylim(0+offset,1+offset)} +
    facet_grid(. ~ Sample, scales = "free_y") +
    {if (NO_GRIDLINES) theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())} +
    theme(axis.text.x = element_text(hjust=0.5, vjust=1), axis.title = element_blank(), legend.text = element_text(face = "italic")) +
    {if (!LEGEND) theme(legend.position = 'none')} +
    {if (NO_X_TEXT) theme(axis.text.x = element_blank(), axis.ticks.x=element_blank())} +
    {if (NO_Y_TEXT) theme(axis.text.y = element_blank(), axis.ticks.y=element_blank())} + 
    {if (!FACET_STRIPS) theme(strip.text = element_blank())} + 
    {if (!is.null(hard_ys)) coord_cartesian(ylim = hard_ys)} + 
    {if (!is.null(hard_ys)) scale_y_continuous(breaks=seq(0, 1, 0.25))
      else scale_y_continuous(expand = expand_scale(expand_ys, 0))} +
    scale_x_continuous(expand = c(0, 0))
  
  return(p)
}


avg_plotter = function(gene, smpl=unique(data$sample),
                       COUNTS=FALSE, SLOWNORM=FALSE,
                       include_SE=TRUE, rib_alpha = 1,
                       include_SE_bars=FALSE, 
                       ZERO_Y=TRUE, data=data_C2C3, no_Zeros = FALSE, 
                       SEB_FILT = FALSE, NORMALIZE_MEGACELL = FALSE,
                       point_size=0.1, NO_X_TEXT=FALSE) {
  if(NORMALIZE_MEGACELL) {Megacell = AggregateExpression(data, group.by = c('sample', 'group1')) %>%
    .$RNA %>% as.data.frame %>% 
    {cbind(Smpl = names(.) %>% gsub(pattern = '_.*', replacement = '') %>%
             gsub(pattern='g', replacement=''), 
           timepoint = names(.) %>% gsub(pattern = '.*ZT', replacement = ''), 
           Value = colSums(.))} %>% as.data.frame %>% 
    subset(Smpl %in% smpl) %>%
    merge(table(data$sample, data$group1) %>% as.data.frame %>%
            mutate(Smpl = as.character(Var1), 
                   timepoint = as.character(Var2 %>% gsub(pattern = 'ZT', replacement = '')), 
                   .keep='unused')) %>%
    mutate_all(as.numeric) %>%
    mutate(Val_normCell = Value/Freq) %>% 
    mutate(c_ZT_factor_norm = Val_normCell/mean(Val_normCell))}
  
  FetchData(data, vars = c("sample", "group1", gene), layer = ifelse(COUNTS,"counts","data")) %>%
    subset(sample %in% smpl) %>%
    pivot_longer(cols = -c(1,2), names_to = "Gene", values_to = "expr") %>%
    {if (no_Zeros) .[which(.$expr != 0),] else .} %>%
    {if (SEB_FILT) Sebastian_filter_per_gene_smpl(., Smpl = smpl, gene = Gene, 
                                                  pct_expressions = data_pcts) else .} %>% 
    {if (NORMALIZE_MEGACELL) Normalize_megacell(., megacells = Megacell) else .} %>% 
    {if (SLOWNORM) {
      merge(., 
            lapply(levels(factor(.$Gene)), function(Gene_) 
              lapply(levels(factor(.$sample)), function(smpl_) 
                lapply(levels(factor(.$group1)), function (time_) {
                  subset(., subset = Gene==Gene_ & group1==time_ & sample==smpl_) %>%
                    .$expr %>% mean() 
                }) %>% unlist() %>% max() %>% c(Gene_, smpl_)
              )
            ) %>% 
              as.data.frame() %>% t() %>% as.data.frame() %>%
              `colnames<-`(c('ZTmax', 'Gene', 'sample')) %>%
              mutate(ZTmax = as.numeric(ZTmax)), 
            by = c('Gene', 'sample')) %>% 
        mutate(expr = expr/ZTmax)
    } else .} %>%
    mutate(Gene = factor(Gene, levels = gene %>% unique), 
           sample = factor(sample, levels = smpl %>% unique)) %>%
    {if (!COUNTS) mutate(., expr = expm1(expr)) else .} %>% 
    ggplot(aes(x = group1, y = expr, group = Gene, col = Gene, fill = Gene)) +
    stat_summary(geom = "line", fun = "mean") +
    {if (include_SE) {
      stat_summary(fun = mean, geom = "ribbon", alpha=rib_alpha, color = NA,
                   fun.max = function(x) mean(x) + sd(x) / sqrt(length(x)),
                   fun.min = function(x) mean(x) - sd(x) / sqrt(length(x)))
    }} +
    {if (include_SE_bars) {
      stat_summary(fun.data = "mean_cl_boot", size = point_size)
    }} +
    {if (ZERO_Y) expand_limits(y=0)} +
    facet_grid(Gene ~ factor(sample, levels = smpl), scales = "free_y") +
    theme(axis.text.x = element_text(size = 6, angle = 45), axis.title = element_blank(), legend.position = "none",
          axis.text.y = element_text(size = 6)) + 
    {if (NO_X_TEXT) theme(axis.text.x = element_blank(), axis.ticks.x=element_blank())}
}


daynorm = function(DF, whichNorm='10', whichCols=NULL) {
  
  if (is.null(whichCols)) whichCols=c(colnames(DF)[colnames(DF) %!in% c('value', 'Time')], 'rep')
  
  if (whichNorm=='10' | whichNorm=='01') {
    DF = DF %>%
      mutate(rep = Time %/% 24) %>%
      group_by(across(all_of(whichCols))) %>%
      reframe(value = (value-min(value))/max(value-min(value)), Time=Time) %>% 
      as.data.frame %>%
      mutate(value = ifelse(is.na(value), 0, value))
  } else if (grep('1', whichNorm)) {
    DF = DF %>%
      mutate(rep = Time %/% 24) %>%
      group_by(across(all_of(whichCols))) %>%
      reframe(value = value/max(value), Time=Time) %>% 
      as.data.frame %>%
      mutate(value = ifelse(is.na(value), 0, value))
  }
  
  return(DF)
}

make_plot = function(Genes = c('Clk', 'cry', 'tim', 'per', 'vri', 'Pdp1', 'cwo', 'Gclc'), 
                     S=data_C2C3, Group='time4', Mode0=NULL, Type0=NULL, ZERO_Y=TRUE, 
                     WhichNorm=NULL, TIME_INT=TRUE, X_TEXT=FALSE, TITLE=NULL, COUNTS=FALSE,
                     PRISM=FALSE, CLASSIC=FALSE, MINIMAL=FALSE,
                     GRIDLESS=TRUE, line_size=1, point_size=2,
                     FACET_STRIPS=FALSE, expand_ys=c(0.05, 0.05)) {
  
  if(is.null(Mode0)) Mode0=levels(S$mode1)
  if(is.null(Type0)) Type0=levels(S$sample)
  
  if(is.factor(Genes)) GeneLevels = levels(Genes)
  else GeneLevels = unique(Genes)
  
  if (!TIME_INT) x_levels = levels(S[[Group]][[1]]) %>%
      gsub('_', ' ', .)
  
  S = subset(S, mode1 %in% Mode0 &
               sample %in% Type0)

  S = AverageExpression(S, features=Genes, group.by = Group, layer = ifelse(COUNTS, 'counts', 'data'))$RNA %>%
    as.data.frame %>%
    {mutate(., Gene = rownames(.))} %>% 
    pivot_longer(cols=-Gene, names_to='Time') %>% 
    mutate(Time = gsub('^g', '', Time) %>%
             gsub('-', ' ', .), 
           Gene = factor(Gene, levels=GeneLevels)) %>% 
    {if (TIME_INT) mutate(., Time = gsub('[^0-9]', '', Time) %>% as.numeric)
      else mutate(., Time = factor(Time, levels=x_levels))} %>% 
    {if (TIME_INT && !is.null(WhichNorm)) daynorm(., whichNorm=WhichNorm)
      else . }
  
  
  p = S %>%
    ggplot(aes(Time, value, group=Gene))
  
  N_days = max(as.numeric(gsub('[^0-9]', '', S$Time)))%/%24 + 1
  
  LD = all(Mode0=='LD')
  
  p = allRect(p, N_days, LD)
  
  p = p +
    geom_point(size = point_size) +
    geom_line(linewidth = line_size) +
    facet_grid(Gene~., scales = 'free') +
    theme_linedraw() + 
    theme(axis.text=element_text(colour="black", size=8)) +
    {if (MINIMAL) theme_minimal()} +
    {if (PRISM) theme_prism()} +
    {if (CLASSIC) theme_classic()} +
    theme(axis.title.y=element_blank()) +
    {if (!is.null(TITLE)) ggtitle(TITLE)} +
    {if (!X_TEXT) theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.x=element_blank())} +
    {if (ZERO_Y) expand_limits(y=0)} +
    {if (GRIDLESS) theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())} +
    {if (!FACET_STRIPS) theme(strip.text = element_blank())} + 
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = expand_scale(expand_ys, 0))
  
  return(p)
}

oneRect = function(P, xMin=0, ygb='yellow', LW=0.1, Alpha=0.3, xMax=NULL) {
  if (is.null(xMax)) xMax = xMin + 12
  
  P = P + geom_rect(data=data.frame(xmin=xMin, xmax=xMax, ymin=-Inf, ymax=Inf), 
                    aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill=ygb, colour='black', linewidth=LW, alpha=Alpha, inherit.aes = FALSE)
}

allRect = function(P, n_days, ld) {
  
  if (n_days == 6) {
    P = P %>% 
      oneRect(xMin=0, ygb='yellow', LW=0.1) %>%
      oneRect(xMin=12, ygb='grey35', LW=0.1) %>%
      oneRect(xMin=24, ygb='yellow', LW=0.1) %>%
      oneRect(xMin=36, ygb='grey35', LW=0.1) %>%
      oneRect(xMin=48, ygb='yellow', LW=0.1) %>%
      oneRect(xMin=60, ygb='grey35', LW=0.1) %>%
      oneRect(xMin=72, ygb='grey35', LW=0.1) %>%
      oneRect(xMin=84, ygb='grey35', LW=0.1) %>%
      oneRect(xMin=96, ygb='grey35', LW=0.1) %>%
      oneRect(xMin=108, ygb='grey35', LW=0.1) %>%
      oneRect(xMin=120, ygb='grey35', LW=0.1) %>%
      oneRect(xMin=132, ygb='grey35', LW=0.1) %>%
      oneRect(xMin=0, xMax=24, Alpha=0, LW=0.4) %>%
      oneRect(xMin=24, xMax=48, Alpha=0, LW=0.4) %>%
      oneRect(xMin=48, xMax=72, Alpha=0, LW=0.4) %>%
      oneRect(xMin=72, xMax=96, Alpha=0, LW=0.4) %>%
      oneRect(xMin=96, xMax=120, Alpha=0, LW=0.4) %>%
      oneRect(xMin=120, xMax=144, Alpha=0, LW=0.4)
  }
  
  if (n_days == 4) {
    P = P %>% 
      oneRect(xMin=0, ygb='yellow', LW=0.1) %>%
      oneRect(xMin=12, ygb='grey35', LW=0.1) %>%
      oneRect(xMin=24, ygb='yellow', LW=0.1) %>%
      oneRect(xMin=36, ygb='grey35', LW=0.1) %>%
      oneRect(xMin=48, ygb='grey35', LW=0.1) %>%
      oneRect(xMin=60, ygb='grey35', LW=0.1) %>%
      oneRect(xMin=72, ygb='grey35', LW=0.1) %>%
      oneRect(xMin=84, ygb='grey35', LW=0.1) %>%
      oneRect(xMin=0, xMax=24, Alpha=0, LW=0.4) %>%
      oneRect(xMin=24, xMax=48, Alpha=0, LW=0.4) %>%
      oneRect(xMin=48, xMax=72, Alpha=0, LW=0.4) %>%
      oneRect(xMin=72, xMax=96, Alpha=0, LW=0.4)
  }
  
  if (n_days == 3) {
    P = P %>%
      oneRect(xMin=0, ygb=ifelse(!ld, 'grey35', 'yellow'), LW=0) %>%
      oneRect(xMin=12, ygb='grey35', LW=0.1) %>%
      oneRect(xMin=24, ygb=ifelse(!ld, 'grey35', 'yellow'), LW=0) %>%
      oneRect(xMin=36, ygb='grey35', LW=0.1) %>%
      oneRect(xMin=48, ygb=ifelse(!ld, 'grey35', 'yellow'), LW=0) %>%
      oneRect(xMin=60, ygb='grey35', LW=0.1) %>%
      oneRect(xMin=0, xMax=24, Alpha=0, LW=0.4) %>%
      oneRect(xMin=24, xMax=48, Alpha=0, LW=0.4) %>%
      oneRect(xMin=48, xMax=72, Alpha=0, LW=0.4)
  }
  
  if (n_days == 2) {
    P = P %>%
      oneRect(xMin=0, ygb=ifelse(!ld, 'grey35', 'yellow1'), LW=0) %>%
      oneRect(xMin=12, ygb='grey35', LW=0.1) %>%
      oneRect(xMin=24, ygb=ifelse(!ld, 'grey35', 'yellow'), LW=0) %>%
      oneRect(xMin=36, ygb='grey35', LW=0.1) %>%
      oneRect(xMin=0, xMax=24, Alpha=0, LW=0.4) %>%
      oneRect(xMin=24, xMax=48, Alpha=0, LW=0.4)
  }
  
  return(P)
}


make_pdf_from_gene_list = function(gene_list, Samples=unique(data$sample), filename, 
                                   COUNTS=FALSE, rib_alpha=1, ZERO_Y=TRUE, include_SE=TRUE, include_SE_bars=FALSE, 
                                   data=data_C2C3, dgrp_colname = 'group2', N_days = 2) {
  
  pdf(file = filename, height=length(gene_list))

  avg_plotter_concat_mode(gene=gene_list, smpl=Samples, dgrp_colname=dgrp_colname, 
                          COUNTS=COUNTS, include_SE=include_SE, include_SE_bars=include_SE_bars, 
                          ZERO_Y=ZERO_Y, data=data, rib_alpha=rib_alpha, N_days=N_days) %>% print
  dev.off()
}
pdf_gene_lists = function(list, Samples=unique(data$sample), 
                          subdivs=ceiling(length(list)/40), name_base, 
                          COUNTS=FALSE, data=data_C2C3, rib_alpha=1, 
                          dgrp_colname='group2', N_days=2) {
  len = length(list)
  for (i in 1:subdivs) {
    filename = paste0(name_base, '_', i, '.pdf')
    print(filename)
    make_pdf_from_gene_list(list[(floor((i-1)*len/subdivs)+1) : floor((i)*len/subdivs)], 
                            filename=filename, COUNTS=COUNTS, rib_alpha=rib_alpha, 
                            Samples=Samples, dgrp_colname=dgrp_colname, N_days=N_days,
                            data=data)
  }
}




make_pdf_filtered = function(to_filter, filter_list=NULL, meta2d_threshold=NULL, name_distinguisher=NULL, ANNOTATE=TRUE,
                             data=data_C2C3, data_cond=data_C2C3_normGroup2, data_AGG=data_C2C3_AGG, data_AGG_CLRnorm=data_C2C3_AGG_CLRnorm,
                              CELL=TRUE, CELL_NO_ERROR=FALSE, COND=FALSE, COND_NO_ERROR=FALSE, AGG_NORM=FALSE, AGG_CLRNORM=FALSE,
                              BIO_RIB_CELL=FALSE, BIO_RIB_COND=FALSE) {
  
  if (!is.null(filter_list)) {
    for (filt in filter_list) {
    to_filter = filterer(to_filter, filt)
    }
  }
  
  to_filter = mutate(to_filter, Sample_Mode = paste(Sample, Mode, sep='_'))
  
  to_filter = lapply(unique(to_filter$Sample_Mode), 
                     function(S_M, df=to_filter) subset(df, Sample_Mode == S_M))
  
  if (CELL) {
    for (df in to_filter) {
      if (!is.null(meta2d_threshold)) df = subset(df, meta2d_pvalue<meta2d_threshold)
      
      name = ifelse(is.null(name_distinguisher), 
                    paste0('plot_cell_', unique(df$Sample_Mode), '_', paste(sapply(filter_list, paste, collapse=""), collapse='_'), '.pdf'),
                    paste0('plot_cell_', unique(df$Sample_Mode), '_', name_distinguisher, '.pdf')
      )
      pdf(name, width=12, height=nrow(df)/8/(1+as.integer(!ANNOTATE)*0.5))
      df$Gene %>% unique %>%
        subplotter_concat(N_days=2, n_plots = 6*(1+as.integer(!ANNOTATE)*0.5), data=subset(data, mode1 == unique(df$Mode) & sample == unique(df$Sample)), ANNOTATE=ANNOTATE, point_size = 1) %>% print
      dev.off()
    }
  }
  if (CELL_NO_ERROR) {
    for (df in to_filter) {
      if (!is.null(meta2d_threshold)) df = subset(df, meta2d_pvalue<meta2d_threshold)
      
      name = ifelse(is.null(name_distinguisher), 
                    paste0('plot_cellNo_', unique(df$Sample_Mode), '_', paste(sapply(filter_list, paste, collapse=""), collapse='_'), '.pdf'),
                    paste0('plot_cellNo_', unique(df$Sample_Mode), '_', name_distinguisher, '.pdf')
      )
      pdf(name, width=12, height=nrow(df)/8/(1+as.integer(!ANNOTATE)*0.5))
      df$Gene %>% unique %>%
        subplotter_concat(N_days=2, n_plots = 6*(1+as.integer(!ANNOTATE)*0.5), rib_alpha=0, data=subset(data, mode1 == unique(df$Mode) & sample == unique(df$Sample)), ANNOTATE=ANNOTATE, point_size = 1) %>% print
      dev.off()
    }
  }
  
  if (COND) {
    for (df in to_filter) {
      if (!is.null(meta2d_threshold)) df = subset(df, meta2d_pvalue_AGG_norm<meta2d_threshold)
      
      name = ifelse(is.null(name_distinguisher), 
                    paste0('plot_cond_', unique(df$Sample_Mode), '_', paste(sapply(filter_list, paste, collapse=""), collapse='_'), '.pdf'),
                    paste0('plot_cond_', unique(df$Sample_Mode), '_', name_distinguisher, '.pdf')
      )
      pdf(name, width=12, height=nrow(df)/8/(1+as.integer(!ANNOTATE)*0.5))
      df$Gene %>% unique %>%
        subplotter_concat(N_days=2, n_plots = 6*(1+as.integer(!ANNOTATE)*0.5), data=subset(data_cond, mode1 == unique(df$Mode) & sample == unique(df$Sample)), ANNOTATE=ANNOTATE, COUNTS=TRUE, point_size = 1) %>% print
      dev.off()
    }
  }
  if (COND_NO_ERROR) {
    for (df in to_filter) {
      if (!is.null(meta2d_threshold)) df = subset(df, meta2d_pvalue_AGG_norm<meta2d_threshold)
      
      name = ifelse(is.null(name_distinguisher), 
                    paste0('plot_condNo_', unique(df$Sample_Mode), '_', paste(sapply(filter_list, paste, collapse=""), collapse='_'), '.pdf'),
                    paste0('plot_condNo_', unique(df$Sample_Mode), '_', name_distinguisher, '.pdf')
      )
      pdf(name, width=12, height=nrow(df)/8/(1+as.integer(!ANNOTATE)*0.5))
      df$Gene %>% unique %>%
        subplotter_concat(N_days=2, n_plots = 6*(1+as.integer(!ANNOTATE)*0.5), rib_alpha=0, data=subset(data_cond, mode1 == unique(df$Mode) & sample == unique(df$Sample)), ANNOTATE=ANNOTATE, COUNTS=TRUE, point_size = 1) %>% print
      dev.off()
    }
  }
  
  if (AGG_NORM) {
    for (df in to_filter) {
      if (!is.null(meta2d_threshold)) df = subset(df, meta2d_pvalue_AGG_norm<meta2d_threshold)
      
      name = ifelse(is.null(name_distinguisher), 
                    paste0('plot_AGGnorm_', unique(df$Sample_Mode), '_', paste(sapply(filter_list, paste, collapse=""), collapse='_'), '.pdf'),
                    paste0('plot_AGGnorm_', unique(df$Sample_Mode), '_', name_distinguisher, '.pdf')
      )
      pdf(name, width=12, height=nrow(df)/8/(1+as.integer(!ANNOTATE)*0.5))
      df$Gene %>% unique %>%
        subplotter_concat(N_days=2, n_plots = 6*(1+as.integer(!ANNOTATE)*0.5), data=subset(data_AGG, mode1 == unique(df$Mode) & sample == unique(df$Sample)), ANNOTATE=ANNOTATE, point_size = 1) %>% print
      dev.off()
    }
  }
  if (AGG_CLRNORM) {
    for (df in to_filter) {
      if (!is.null(meta2d_threshold)) df = subset(df, meta2d_pvalue_AGG_CLRnorm<meta2d_threshold)
      
      name = ifelse(is.null(name_distinguisher), 
                    paste0('plot_AGG_CLRnorm_', unique(df$Sample_Mode), '_', paste(sapply(filter_list, paste, collapse=""), collapse='_'), '.pdf'),
                    paste0('plot_AGG_CLRnorm_', unique(df$Sample_Mode), '_', name_distinguisher, '.pdf')
      )
      pdf(name, width=12, height=nrow(df)/8/(1+as.integer(!ANNOTATE)*0.5))
      df$Gene %>% unique %>%
        subplotter_concat(N_days=2, n_plots = 6*(1+as.integer(!ANNOTATE)*0.5), data=subset(data_AGG_CLRnorm, mode1 == unique(df$Mode) & sample == unique(df$Sample)), ANNOTATE=ANNOTATE, point_size = 1) %>% print
      dev.off()
    }
  }
  
  # if (BIO_RIB_CELL*0) {
  #   for (df in to_filter) {
  #     if (!is.null(meta2d_threshold)) df = subset(df, meta2d_pvalue<meta2d_threshold)
  # 
  #     name = ifelse(is.null(name_distinguisher),
  #                   paste0('plot_bioCell_', unique(df$Sample_Mode), '_', paste(sapply(filter_list, paste, collapse=""), collapse='_'), '.pdf'),
  #                   paste0('plot_bioCell_', unique(df$Sample_Mode), '_', name_distinguisher, '.pdf')
  #     )
  #     pdf(name, width=12, height=nrow(df)/16*(1+ANNOTATE))
  #     df$Gene %>% unique %>%
  #       subplotter_concat(func = 'avg_plotter_repDot', n_plots = 6*(1+!ANNOTATE), data=subset(data, mode1 == unique(df$Mode)),
  #                         more_args_list=list(NO_DD=unique(df$Mode)=='LD', NO_LD=unique(df$Mode)=='DD')) %>% print
  #     dev.off()
  #   }
  # }
  # if (BIO_RIB_COND*0) {
  #   for (df in to_filter) {
  #     if (!is.null(meta2d_threshold)) df = subset(df, meta2d_pvalue_AGG_norm<meta2d_threshold)
  # 
  #     name = ifelse(is.null(name_distinguisher),
  #                   paste0('plot_bioCond_', unique(df$Sample_Mode), '_', paste(sapply(filter_list, paste, collapse=""), collapse='_'), '.pdf'),
  #                   paste0('plot_bioCond_', unique(df$Sample_Mode), '_', name_distinguisher, '.pdf')
  #     )
  #     pdf(name, width=12, height=nrow(df)/16*(1+ANNOTATE))
  #     df$Gene %>% unique %>%
  #       subplotter_concat(func = 'avg_plotter_repDot', n_plots = 6*(1+!ANNOTATE), data=subset(data_cond, mode1 == unique(df$Mode)), COUNTS=TRUE,
  #                         more_args_list=list(NO_DD=unique(df$Mode)=='LD', NO_LD=unique(df$Mode)=='DD')) %>% print
  #     dev.off()
  #   }
  # }
  
  return(to_filter)
}


pdf_combine_and_delete_N = function(name_base, N) {
  qpdf::pdf_combine(input = paste0(name_base, '_', 1:N, '.pdf'),
                    output = paste0(name_base, '.pdf'))
  lapply(paste0(name_base, '_', 1:N, '.pdf'), file.remove)
}
pdf_combine_and_delete = function(name_base, N=NULL) { #N doesn't matter
  files = list.files(pattern=paste0("^", name_base))
  qpdf::pdf_combine(input = files,
                    output = paste0(name_base, '.pdf'))
  lapply(files, file.remove)
}




hmap_phase = function(df_genes=filtered_preliminary, 
                      avg_data_times=avg_data_times2_01, 
                      sample_remove=NULL, mode_remove=NULL, TITLES=FALSE, background_fill=NULL, 
                      Colors=color_list_alt, col_low='white', COLOR_GENES=FALSE, SQUARE_COLORS=FALSE) {
  
  avg_data_time = avg_data_times %>% 
    mutate(sample_mode_time = NULL) %>%
    pivot_wider(names_from = timepoint, values_from = avg)
  
  avg_data_time = avg_data_time %>%
    as.data.frame %>%
    mutate(Mode = factor(Mode, levels=c('LD', 'DD'))) %>%
    subset(Sample %!in% sample_remove & Mode %!in% mode_remove) %>% 
    mutate(Sample_Mode = paste(Sample, Mode, sep='_') %>%
             factor(levels = c('C3_LD', 'C3_DD', 'C2_LD', 'C2_DD')), .keep='unused')
  
  df_genes = df_genes %>%
    as.data.frame %>%
    mutate(Mode = factor(Mode, levels=c('LD', 'DD'))) %>%
    subset(Sample %!in% sample_remove & Mode %!in% mode_remove) %>% 
    mutate(Sample_Mode = paste(Sample, Mode, sep='_') %>%
             factor(levels = c('C3_LD', 'C3_DD', 'C2_LD', 'C2_DD')), .keep='unused')
  
  Sample_Modes = unique(avg_data_time$Sample_Mode)
  
  res = list()
  for (S_M in Sample_Modes) {
    res$Genes[[S_M]] = df_genes[df_genes$Sample_Mode %in% S_M, 'Gene']

    #Color = hcl(h = seq(15, 375, length = length(Sample_Modes)+1), l = 65, c = 100)[which(Sample_Modes == S_M)]
    Color = Colors[which(Sample_Modes == S_M)]
    
    print(S_M)
    
    df = avg_data_time %>% 
      subset(Gene %in% res$Genes[[S_M]] & Sample_Mode %in% S_M)
    
    df = df %>% 
      {`rownames<-`(., .$Gene)} %>% 
      subset(select = c(-Gene, -Sample_Mode)) %>%
      #t %>% scale %>% t %>% as.data.frame %>%
      {mutate(., Gene = rownames(.))} %>%
      mutate(Gene = factor(Gene, levels = arrange(subset(df_genes, Sample_Mode %in% S_M), 
                                                  desc(meta2d_phase))$Gene %>% unique)) %>%
      arrange(Gene) 
    
    if (COLOR_GENES) {
      df = cbind(df, Color = hcl(h = seq(15, 375, length = nrow(df)+1), l = 65, c = 100) %>% {.[1:nrow(df)]} %>% rev )
      
      p = df %>%
        pivot_longer(cols = as.character(seq(3,max(as.numeric(colnames(df)), na.rm=TRUE),4))) %>% 
        mutate(timepoint = factor(name, levels = unique(name)), 
               value = ifelse(is.nan(value), 0, value), .keep='unused') %>%
        ggplot(aes(timepoint, Gene)) +
        {if(SQUARE_COLORS) geom_tile(aes(alpha = value^2, fill=Gene)) 
          else geom_tile(aes(alpha = value, fill=Gene))} +
        scale_fill_manual(values=df$Color)
      
    } else {
      p = df %>%
        pivot_longer(cols = as.character(seq(3,max(as.numeric(colnames(df)), na.rm=TRUE),4))) %>% 
        mutate(timepoint = factor(name, levels = unique(name)), 
               value = ifelse(is.nan(value), 0, value), .keep='unused') %>%
        ggplot(aes(timepoint, Gene)) +
        geom_tile(aes(fill=value)) +
        scale_fill_gradient(low=col_low, high=Color) 
    }
    
    res$plots[[S_M]] = p + 
      theme_bw() + 
      theme(axis.line=element_blank(),
            axis.text.x=element_blank(),
            axis.text.y=element_blank(),axis.ticks=element_blank(),
            axis.title.x=element_blank(),
            axis.title.y=element_blank(),legend.position="none", 
            panel.border=element_blank(), 
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
      {if(is.null(background_fill)) theme(panel.background=element_blank()) 
        else theme(panel.background=element_rect(fill = background_fill))} +
      {if(TITLES) ggtitle(S_M)}
  }

  plot_grid(plotlist=res$plots, ncol=1, rel_heights=unlist(lapply(res$Genes, length)))
}


hist_phase(Filtereds$DFs$reps_adding_pool , 
           CIRCULAR=FALSE, 
           TRIM_TO_24=TRUE, 
           binWidth=2, 
           n_W=0, 
           Colors=rep('red', 4) %>% `names<-`(c('C3_LD', 'C3_DD', 'C2_LD', 'C2_DD')), 
           HIST_FRONT=TRUE, LD=FALSE, to_select='meta2d_phase')

plot_grid(
Filtereds$DFs$reps_adding_pool %>%
  ggplot(aes(meta2d_phase)) +
  geom_histogram(bins=12) +
  facet_grid(Sample~Mode),
Filtereds$DFs$reps_adding_pool %>%
  mutate(meta2d_phase = meta2d_phase %% 24) %>% 
  ggplot(aes(meta2d_phase)) +
  geom_histogram(bins=12) +
  facet_grid(Sample~Mode)
)
   
plot_grid(
  to_filter %>%
    ggplot(aes(meta2d_phase)) +
    geom_histogram(bins=12) +
    facet_grid(Sample~Mode),
  to_filter %>%
    mutate(meta2d_phase = meta2d_phase %% 24) %>% 
    ggplot(aes(meta2d_phase)) +
    geom_histogram(bins=12) +
    facet_grid(Sample~Mode)
)



hist_phase = function(df_genes=filtered_preliminary, to_select = 'meta2d_phase_AGG_norm', CIRCULAR=FALSE, TRIM_TO_24=FALSE, CHOP_TO_24=TRUE, binWidth=2, 
                      kernel_adjust=1, n_W=0, W_hr=0.5, TITLE=FALSE, YSCALE=FALSE, FREEY=TRUE, Colors = color_list, LD=TRUE, HIST_FRONT=TRUE) {
  period_max = 28
  
  df_genes = df_genes[, c('Gene', 'Sample', 'Mode', to_select)] %>%
    `colnames<-`(c('Gene', 'Sample', 'Mode', 'meta2d_phase'))
  
  df_genes = mutate(df_genes, Sample_Mode = paste0(Sample, '_', Mode) %>% 
                      gsub('_', ' ', .) %>%
                      factor(levels = c('C3 LD', 'C3 DD', 'C2 LD', 'C2 DD')))
  
  names(Colors) = gsub('_', ' ', names(Colors))
  
  if (TRIM_TO_24 | CHOP_TO_24) {
    if (CHOP_TO_24) {
      df_genes = subset(df_genes, meta2d_phase <= 24)
    } else {
      df_genes = mutate(df_genes, meta2d_phase = ifelse(meta2d_phase >= 24, 24, meta2d_phase))
    }
    period_max = 24
    
    if (n_W) {
      df_genes = split(df_genes, df_genes$Sample_Mode) %>% 
        lapply(function(df) {
          {cut(df$meta2d_phase, breaks=seq(0, 24, W_hr), labels=1:(24/W_hr)*W_hr-(W_hr/2))} %>%
            table %>%
            as.data.frame %>%
            `colnames<-`(c('val', 'Freq')) %>% 
            {zoo::zoo(.$Freq, .$val)} %>%
            #rollapply(., n_W, mean) %>%
            rollmean(., k=n_W, fill = NA, align = 'center') %>%
            na.locf() %>%
            as.data.frame %>% 
            `colnames<-`('count') %>%
            {mutate(., time = rownames(.))} %>%
            apply(1, function(L) {rep(L[2], L[1])}) %>%
            do.call(what='c') %>%
            {data.frame(meta2d_phase = as.integer(.), Sample_Mode = df[1, 'Sample_Mode'])}
        }) %>% 
        do.call(what='rbind')
    }
  }
  
  
  p = df_genes %>%
    ggplot(aes(x = meta2d_phase)) +
    {if(YSCALE) geom_histogram(aes(fill = Sample_Mode), breaks=seq(0,24,binWidth), color='black')
      else geom_histogram(aes(y = after_stat(density), fill = Sample_Mode), breaks=seq(0,24,binWidth), color='black')} +
    #facet_grid(Mode ~ Sample) +
    #facet_grid(Sample_Mode ~ .) +
    facet_wrap('Sample_Mode', ncol=1, scales=ifelse(FREEY, 'free', 'fixed')) +
    theme_minimal() +
    theme(axis.title = element_blank()) +
    guides(fill = "none", colour='none') +
    scale_x_continuous(breaks = seq(0, period_max, 2)) +
    {if (!YSCALE) theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())} +
    scale_fill_manual(values = Colors)
  
  day=1
  BRIGHTEN_DD_DAY=TRUE
    if (LD) {
      p = p +
        geom_rect(data=data.frame(xmin=0+(24*(day-1)), xmax=12+(24*(day-1)), ymin=-0, ymax=Inf), 
                  aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="yellow", colour='black', linewidth=0.4, alpha=0.3, inherit.aes = FALSE) +
        geom_rect(data=data.frame(xmin=12+(24*(day-1)), xmax=24+(24*(day-1)), ymin=-0, ymax=Inf), 
                  aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey35", colour='black', linewidth=0.4, alpha=0.3, inherit.aes = FALSE)
    } else {
      p = p + 
        geom_rect(data=data.frame(xmin=0+(24*(day-1)), xmax=12+(24*(day-1)), ymin=-0, ymax=Inf), 
                  aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill=ifelse(BRIGHTEN_DD_DAY, "grey75", "grey35"), colour='black', linewidth=0.4, alpha=0.3, inherit.aes = FALSE) +
        geom_rect(data=data.frame(xmin=12+(24*(day-1)), xmax=24+(24*(day-1)), ymin=-0, ymax=Inf), 
                  aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey35", colour='black', linewidth=0.4, alpha=0.3, inherit.aes = FALSE)
    }
  
  if (CIRCULAR) {
    p = p +
      coord_radial(theta = "x", expand = FALSE)
  } else if (!YSCALE) {
    p = p +
      geom_density(adjust = kernel_adjust)
  }
  if (!CIRCULAR) p = p + 
    ylim(0, NA) +
    scale_x_continuous(expand = c(0, 0), breaks = seq(0, period_max, 6)) +
    scale_y_continuous(expand = c(0, 0))
    
  if(TITLE) p = p + ggtitle(n_W*W_hr)
  
  if(HIST_FRONT) p$layers = rev(p$layers)
  return(p)
}


counter = function(DF) {
  if (nrow(DF)) ret = data.frame(table(DF$Sample, DF$Mode)) 
  else ret = data.frame(c('C2', 'C3', 'C2', 'C3'), c('DD', 'DD', 'LD', 'LD'), rep(0, 4))
  
  colnames(ret) <- c('Sample', 'Mode', 'count')
  ret = ret %>%
    mutate(Sample_Mode = paste(Sample, Mode, sep='_'), .keep='unused') %>%
    pivot_wider(names_from=Sample_Mode, values_from = 'count')
  return(ret)
}


BH_procedure = function(DF, Q, Qcol='BH.Q', Pcol='ADJ.P', Qname='BH.Q', CUT=TRUE) {
  # ASSUMING NA's ARE REMOVED PRIOR IF EXCLUSION INTENDED
  
  if (is.null(Qcol)) {
    DF = cbind(DF, BH.Q = p.adjust(DF[[Pcol]], method='BH'))
    Qcol = 'BH.Q'
  }
  
  Pmax = max(DF[ DF[[Qcol]] < Q, Pcol ], na.rm = TRUE)
  
  if (CUT) DF = DF[which(DF[[Pcol]] <= Pmax), ]
  else DF = cbind(DF, Qpass = DF[[Pcol]] <= Pmax)
  
  colnames(DF) <- gsub('BH.Q', Qname, colnames(DF))
  return(DF)
}

substatter = function(Filtered, Subsets=Subset_Stats$subsets, R_NOT_P = TRUE) {
  
  if (R_NOT_P) {
    lapply(Subsets, function(Set) {
      lapply(levels(data_C2C3$sample_mode), function(S_M) {
        S = gsub('_.*', '', S_M)
        M = gsub('.*_', '', S_M)
        data.frame(Gene = Set) %>%
          make_corrplot(GO_anno = passQC, filtered=Filtered, sample=S, mode=M) %>%
          from_corrplot_meanmaker
      }) %>%
        do.call(what='cbind') %>%
        `colnames<-`(levels(data_C2C3$sample_mode)) 
    }) %>%
      do.call(what='rbind') %>%
      `rownames<-`(names(Subsets))
  } else {
    Subset_Stats$p_vals_pool3 <- lapply(Subsets, function(Set) {
      lapply(levels(data_C2C3$sample_mode), function(S_M) {
        S = gsub('_.*', '', S_M)
        M = gsub('.*_', '', S_M)
        run_fisher_general(B = subset(Filtered, Sample == S & Mode == M)$Gene,
                           A = subset(to_filter, Gene %in% Set & second_pct > 0.10 & Sample == S & Mode == M)$Gene, 
                           U = subset(to_filter, second_pct > 0.10 & Sample == S & Mode == M)$Gene)
      }) %>%
        do.call(what='cbind') %>%
        `colnames<-`(levels(data_C2C3$sample_mode))
    }) %>%
      do.call(what='rbind') %>%
      `rownames<-`(names(Subsets))
  }
}


####    ####    ####    ####    ####    ####    ####    ####    ####    ####   


# 4. PREPARING KW, meta2d
################################################################################
KW_scores_C2C3_LD <- merge(rownames(data_C2C3), unique(data_C2C3$sample)) %>%
  setNames(c('Gene', 'Sample')) %>%
  append_kruskal_scores_from_above(data=subset(data_C2C3, mode1 == 'LD'), zt_name='group1', COUNTS = FALSE) %>%
  mutate(KW_scores_norm = replace(KW_scores_norm, is.na(KW_scores_norm), 1))
#write.csv(KW_scores_C2C3_LD, 'KW_scores_C2C3_LD.csv', row.names=FALSE)
KW_scores_C2C3_LD <- read.csv("KW_scores_C2C3_LD.csv")

KW_scores_C2C3_DD <- merge(rownames(data_C2C3), unique(data_C2C3$sample)) %>%
  setNames(c('Gene', 'Sample')) %>% 
  append_kruskal_scores_from_above(data=subset(data_C2C3, mode1 == 'DD'), zt_name='group2', COUNTS = FALSE) %>%
  mutate(KW_scores_norm = replace(KW_scores_norm, is.na(KW_scores_norm), 1))
#write.csv(KW_scores_C2C3_DD, 'KW_scores_C2C3_DD.csv', row.names=FALSE)
KW_scores_C2C3_DD <- read.csv("KW_scores_C2C3_DD.csv")

KW_scores_C2C3 <- rbind(cbind(KW_scores_C2C3_LD, Mode = 'LD'), cbind(KW_scores_C2C3_DD, Mode = 'DD')); rm(KW_scores_C2C3_LD, KW_scores_C2C3_DD)

to_meta_C2C3_group2 <- AverageExpression(data_C2C3, group.by = c("group2", "sample", "mode1"), layer="data")$RNA %>% as.data.frame %>%
  {mutate(., Gene = rownames(.))} %>%
  pivot_longer(cols=-'Gene', names_to='TP_sample_mode') %>%
  mutate(Timepoint = TP_sample_mode %>% gsub('_.*', '', .) %>% gsub('.*T', '', .), 
         Gene_SM = paste(Gene, TP_sample_mode %>% gsub('.*_C', 'C', .), sep = '_'), .keep='unused') %>%
  pivot_wider(names_from = 'Timepoint', values_from = 'value')
  
write.csv(to_meta_C2C3_group2, file = "to_meta_C2C3_group2.csv", row.names=FALSE)

meta2d(infile = "to_meta_C2C3_group2.csv", filestyle = "csv", timepoints = "line1", minper=20, 
       maxper=28, outdir=getwd(), parallelize = TRUE)

meta_C2C3_group2 <- read.csv('meta2d_to_meta_C2C3_group2.csv') %>%
  mutate(Gene = gsub('_.*', "", CycID), 
         Mode = gsub('.*_', "", CycID), 
         Sample = gsub('.*_C', 'C', CycID) %>% sub('_.*', '', .),
         .keep='unused')
#write.csv(meta_C2C3_group2, file = "meta_C2C3_group2.csv", row.names=FALSE)
meta_C2C3_group2 <- read.csv('meta_C2C3_group2.csv')

# group2 divides group2
{N = 2
  for (i in 1:N) {
    subCol = 'group2'
    tag = paste0('meta_C2C3_group2_g', i)
    times = seq(3, 23, 4)+24*(i-1)
    infile = paste0('to_', tag, '.csv')
    fromfile = paste0('meta2d_to_', tag, '.csv')
    outfile = paste0(tag, '.csv')
    
    to_meta = data_C2C3[, gsub('.*T', '', data_C2C3@meta.data[[subCol]]) %in% times] %>%
      AverageExpression(group.by = c("group2", "sample", "mode1"), layer="data") %>%
      .$RNA %>% as.data.frame %>%
      {mutate(., Gene = rownames(.))} %>%
      pivot_longer(cols=-'Gene', names_to='TP_sample_mode') %>%
      mutate(Timepoint = TP_sample_mode %>% gsub('_.*', '', .) %>% gsub('.*T', '', .), 
             Gene_SM = paste(Gene, TP_sample_mode %>% gsub('.*_C', 'C', .), sep = '_'), .keep='unused') %>%
      pivot_wider(names_from = 'Timepoint', values_from = 'value')
    
    write.csv(to_meta, file = infile, row.names=FALSE)
    
    meta2d(infile = infile, filestyle = "csv", timepoints = "line1", minper=20, 
           maxper=28, outdir=getwd(), parallelize = TRUE)
    
    meta <- read.csv(fromfile) %>%
      mutate(Gene = gsub('_.*', "", CycID), 
             Mode = gsub('.*_', "", CycID), 
             Sample = gsub('.*_C', 'C', CycID) %>% sub('_.*', '', .),
             .keep='unused')
    write.csv(meta, file = outfile, row.names=FALSE)
  }
  rm(N, i, to_meta, meta)
}

{N = 2 # to import:
  for (i in 1:N) {
    outfile = paste0('meta_C2C3_group2_g', i, '.csv')
    tag = paste0('meta_C2C3_group2_g', i)
    
    meta_C2C3_multi[['group2']][[i]] = read.csv(outfile) %>%
      cbind(i = i)
  }
rm(N, i, outfile, tag)
meta_C2C3_multi[['group2']] <- meta_C2C3_multi[['group2']] %>%
  do.call(what='rbind')
}



# group2-replicates of pcat2
{N = 2
  for (i in 1:N) {
    subCol = paste0('group', N)
    tag = paste0('meta_C2C3_pcat2_2s', i)
    times = seq(3, 23, 4)+24*(i-1)
    infile = paste0('to_', tag, '.csv')
    fromfile = paste0('meta2d_to_', tag, '.csv')
    outfile = paste0(tag, '.csv')
    
    to_meta = data_C2C3[, gsub('.*T', '', data_C2C3@meta.data[[subCol]]) %in% times] %>%
      AverageExpression(group.by = c("pcat2", "sample", "mode1"), layer="data") %>%
      .$RNA %>% as.data.frame %>%
      {mutate(., Gene = rownames(.))} %>%
      pivot_longer(cols=-'Gene', names_to='TP_sample_mode') %>%
      mutate(Timepoint = TP_sample_mode %>% gsub('_.*', '', .) %>% gsub('.*T', '', .) %>% gsub('g', '', .), 
             Gene_SM = paste(Gene, TP_sample_mode %>% gsub('.*_C', 'C', .), sep = '_'), .keep='unused') %>%
      pivot_wider(names_from = 'Timepoint', values_from = 'value')
    
    write.csv(to_meta, file = infile, row.names=FALSE)
    
    meta2d(infile = infile, filestyle = "csv", timepoints = "line1", minper=20, 
           maxper=28, outdir=getwd(), parallelize = TRUE)
    
    meta <- read.csv(fromfile) %>%
      mutate(Gene = gsub('_.*', "", CycID), 
             Mode = gsub('.*_', "", CycID), 
             Sample = gsub('.*_C', 'C', CycID) %>% sub('_.*', '', .),
             .keep='unused')
    write.csv(meta, file = outfile, row.names=FALSE)
  }
  rm(N, i, to_meta, meta)
}

{N = 2 # to import:
  for (i in 1:N) {
    outfile = paste0('meta_C2C3_pcat2_2s', i, '.csv')
    tag = paste0('meta_C2C3_pcat2_2s', i)
    
    meta_C2C3_multi[['group2-pcat2']][[i]] = read.csv(outfile) %>%
      cbind(i = i)
  }
rm(N, i, outfile, tag)
meta_C2C3_multi[['group2-pcat2']] <- meta_C2C3_multi[['group2-pcat2']] %>%
  do.call(what='rbind')
}

# pcat-2 replicates of group2
{N = 2
  for (i in 1:N) {
    subCol = paste0('pcat', N)
    tag = paste0('meta_C2C3_group2_2s', i)
    times = seq(3, 23, 4)+24*(i-1)
    infile = paste0('to_', tag, '.csv')
    fromfile = paste0('meta2d_to_', tag, '.csv')
    outfile = paste0(tag, '.csv')
    
    to_meta = data_C2C3[, data_C2C3@meta.data[[subCol]] %in% times] %>%
      AverageExpression(group.by = c("group2", "sample", "mode1"), layer="data") %>%
      .$RNA %>% as.data.frame %>%
      {mutate(., Gene = rownames(.))} %>%
      pivot_longer(cols=-'Gene', names_to='TP_sample_mode') %>%
      mutate(Timepoint = TP_sample_mode %>% gsub('_.*', '', .) %>% gsub('.*T', '', .), 
             Gene_SM = paste(Gene, TP_sample_mode %>% gsub('.*_C', 'C', .), sep = '_'), .keep='unused') %>%
      pivot_wider(names_from = 'Timepoint', values_from = 'value')
    
    write.csv(to_meta, file = infile, row.names=FALSE)
    
    meta2d(infile = infile, filestyle = "csv", timepoints = "line1", minper=20, 
           maxper=28, outdir=getwd(), parallelize = TRUE)
    
    meta <- read.csv(fromfile) %>%
      mutate(Gene = gsub('_.*', "", CycID), 
             Mode = gsub('.*_', "", CycID), 
             Sample = gsub('.*_C', 'C', CycID) %>% sub('_.*', '', .),
             .keep='unused')
    write.csv(meta, file = outfile, row.names=FALSE)
  }
  rm(N, i, to_meta, meta)
}

{N = 2 # to import:
  for (i in 1:N) {
    outfile = paste0('meta_C2C3_group2_2s', i, '.csv')
    tag = paste0('meta_C2C3_group2_2s', i)
    
    meta_C2C3_multi[['pcat2-group2']][[i]] = read.csv(outfile) %>%
      cbind(i = i)
  }
rm(N, i, outfile, tag)
meta_C2C3_multi[['pcat2-group2']] <- meta_C2C3_multi[['pcat2-group2']] %>%
  do.call(what='rbind')
}

meta_C2C3_multi$`group2` %>%
  .[, c('Gene', 'Mode', 'Sample', 'JTK_pvalue', 'JTK_adjphase', 'meta2d_pvalue', 'meta2d_phase', 'i')] %>%
  subset(Gene == 'Hr38') # NO: no sr pass below 0.4, Hr38 is lost, cwo-0.5, Pdp1 C3 is lost - can't do this

meta_C2C3_multi$`pcat2-group2` %>%
  .[, c('Gene', 'Mode', 'Sample', 'JTK_pvalue', 'JTK_adjphase', 'meta2d_pvalue', 'meta2d_phase', 'i')] %>%
  subset(Gene == 'cwo') # 0.5 needed to pass sr LD, cwo fails C3-LD

meta_C2C3_multi$`group2-pcat2` %>%
  .[, c('Gene', 'Mode', 'Sample', 'JTK_pvalue', 'JTK_adjphase', 'meta2d_pvalue', 'meta2d_phase')] %>% 
  subset(Gene == 'cwo') # cwo fails all


# NORMALIZED AGGREGATES - cell-norm
to_meta_C2C3_AGG_group2 <- AverageExpression(data_C2C3_AGG, group.by = c("group2", "sample", "mode1"), layer="data")$RNA %>% as.data.frame %>%
  {mutate(., Gene = rownames(.))} %>%
  pivot_longer(cols=-'Gene', names_to='TP_sample_mode') %>%
  mutate(Timepoint = TP_sample_mode %>% gsub('_.*', '', .) %>% gsub('.*T', '', .), 
         Gene_SM = paste(Gene, TP_sample_mode %>% gsub('.*_C', 'C', .), sep = '_'), .keep='unused') %>%
  pivot_wider(names_from = 'Timepoint', values_from = 'value')

write.csv(to_meta_C2C3_AGG_group2, file = "to_meta_C2C3_AGG_group2.csv", row.names=FALSE)

meta2d(infile = "to_meta_C2C3_AGG_group2.csv", filestyle = "csv", timepoints = "line1", minper=20, 
       maxper=28, outdir=getwd(), parallelize = TRUE)

meta_C2C3_AGG_group2 <- read.csv('meta2d_to_meta_C2C3_AGG_group2.csv') %>%
  mutate(Gene = gsub('_.*', "", CycID), 
         Mode = gsub('.*_', "", CycID), 
         Sample = gsub('.*_C', 'C', CycID) %>% sub('_.*', '', .),
         .keep='unused')
write.csv(meta_C2C3_AGG_group2, file = "meta_C2C3_AGG_group2.csv", row.names=FALSE)
meta_C2C3_AGG_group2 <- read.csv('meta_C2C3_AGG_group2.csv')

# FOR MEGACELL_NORM (p-values equivalent to condition-norm)
MC_norm %>%
  mutate(Gene_SM = paste0(Gene, '_', Sample, '__', Mode), .keep='unused') %>%
  pivot_wider(names_from = "timepoint", values_from = 'exp') %>%
  write.csv('to_meta_MC_norm_group2.csv', row.names=FALSE)

meta2d(infile = "to_meta_MC_norm_group2.csv", filestyle = "csv", timepoints = "line1", minper=20, 
       maxper=28, outdir=getwd(), parallelize = TRUE)

meta_MC_norm_group2 <- read.csv('meta2d_to_meta_MC_norm_group2.csv') %>%
  mutate(Gene = gsub('_.*', "", CycID), 
         Mode = gsub('.*_', "", CycID), 
         Sample = gsub('.*_C', 'C', CycID) %>% sub('_.*', '', .),
         .keep='unused')
#write.csv(meta_MC_norm_group2, file = "meta_MC_norm_group2.csv", row.names=FALSE)
meta_MC_norm_group2 <- read.csv('meta_MC_norm_group2.csv')


# NORMALIZED AGGREGATES - CLR-norm
to_meta_C2C3_AGG_CLRnorm_group2 <- AverageExpression(data_C2C3_AGG_CLRnorm, group.by = c("group2", "sample", "mode1"), layer="data")$RNA %>% as.data.frame %>%
  {mutate(., Gene = rownames(.))} %>%
  pivot_longer(cols=-'Gene', names_to='TP_sample_mode') %>%
  mutate(Timepoint = TP_sample_mode %>% gsub('_.*', '', .) %>% gsub('.*T', '', .), 
         Gene_SM = paste(Gene, TP_sample_mode %>% gsub('.*_C', 'C', .), sep = '_'), .keep='unused') %>%
  pivot_wider(names_from = 'Timepoint', values_from = 'value')

write.csv(to_meta_C2C3_AGG_CLRnorm_group2, file = "to_meta_C2C3_AGG_CLRnorm_group2.csv", row.names=FALSE)

meta2d(infile = "to_meta_C2C3_AGG_CLRnorm_group2.csv", filestyle = "csv", timepoints = "line1", minper=20, 
       maxper=28, outdir=getwd(), parallelize = TRUE)

meta_C2C3_AGG_CLRnorm_group2 <- read.csv('meta2d_to_meta_C2C3_AGG_CLRnorm_group2.csv') %>%
  mutate(Gene = gsub('_.*', "", CycID), 
         Mode = gsub('.*_', "", CycID), 
         Sample = gsub('.*_C', 'C', CycID) %>% sub('_.*', '', .),
         .keep='unused')
write.csv(meta_C2C3_AGG_CLRnorm_group2, file = "meta_C2C3_AGG_CLRnorm_group2.csv", row.names=FALSE)
meta_C2C3_AGG_CLRnorm_group2 <- read.csv('meta_C2C3_AGG_CLRnorm_group2.csv')




# NORMALIZED AGGREGATES - 1-day-max-norm
to_meta_C2C3_AGG_max1 <- AverageExpression(data_C2C3_AGG_max1, group.by = c("group2", "sample", "mode1"), layer="counts")$RNA %>% as.data.frame %>%
  {mutate(., Gene = rownames(.))} %>%
  pivot_longer(cols=-'Gene', names_to='TP_sample_mode') %>%
  mutate(Timepoint = TP_sample_mode %>% gsub('_.*', '', .) %>% gsub('.*T', '', .), 
         Gene_SM = paste(Gene, TP_sample_mode %>% gsub('.*_C', 'C', .), sep = '_'), .keep='unused') %>%
  pivot_wider(names_from = 'Timepoint', values_from = 'value')

write.csv(to_meta_C2C3_AGG_max1, file = "to_meta_C2C3_AGG_max1.csv", row.names=FALSE)

meta2d(infile = "to_meta_C2C3_AGG_max1.csv", filestyle = "csv", timepoints = "line1", minper=20, 
       maxper=28, outdir=getwd(), parallelize = TRUE)

meta_C2C3_AGG_max1 <- read.csv('meta2d_to_meta_C2C3_AGG_max1.csv') %>%
  mutate(Gene = gsub('_.*', "", CycID), 
         Mode = gsub('.*_', "", CycID), 
         Sample = gsub('.*_C', 'C', CycID) %>% sub('_.*', '', .),
         .keep='unused')
write.csv(meta_C2C3_AGG_max1, file = "meta_C2C3_AGG_max1.csv", row.names=FALSE)
meta_C2C3_AGG_max1 <- read.csv('meta_C2C3_AGG_max1.csv')



to_meta_C2C3_AGG_max01 <- AverageExpression(data_C2C3_AGG_01_1d, group.by = c("group2", "sample", "mode1"), layer="counts")$RNA %>% as.data.frame %>%
  {mutate(., Gene = rownames(.))} %>%
  pivot_longer(cols=-'Gene', names_to='TP_sample_mode') %>%
  mutate(Timepoint = TP_sample_mode %>% gsub('_.*', '', .) %>% gsub('.*T', '', .), 
         Gene_SM = paste(Gene, TP_sample_mode %>% gsub('.*_C', 'C', .), sep = '_'), .keep='unused') %>%
  mutate(value = ifelse(is.na(value), 0, value)) %>%
  pivot_wider(names_from = 'Timepoint', values_from = 'value')

write.csv(to_meta_C2C3_AGG_max01, file = "to_meta_C2C3_AGG_max01.csv", row.names=FALSE)

meta2d(infile = "to_meta_C2C3_AGG_max01.csv", filestyle = "csv", timepoints = "line1", minper=20, 
       maxper=28, outdir=getwd(), parallelize = TRUE)

meta_C2C3_AGG_max01 <- read.csv('meta2d_to_meta_C2C3_AGG_max01.csv') %>%
  mutate(Gene = gsub('_.*', "", CycID), 
         Mode = gsub('.*_', "", CycID), 
         Sample = gsub('.*_C', 'C', CycID) %>% sub('_.*', '', .),
         .keep='unused')
write.csv(meta_C2C3_AGG_max01, file = "meta_C2C3_AGG_max01.csv", row.names=FALSE)
meta_C2C3_AGG_max01 <- read.csv('meta_C2C3_AGG_max01.csv')



# time4

to_meta_C2C3_time4 <- AverageExpression(data_C2C3, group.by = c("time4", "sample"), layer="data")$RNA %>% as.data.frame %>%
  {mutate(., Gene = rownames(.))} %>%
  pivot_longer(cols=-'Gene', names_to='TP_sample') %>%
  mutate(Timepoint = TP_sample %>% gsub('_.*', '', .) %>% gsub('.*T', '', .) %>% gsub('g', '', .) , 
         Gene_S = paste(Gene, TP_sample %>% gsub('.*_C', 'C', .), sep = '_'), .keep='unused') %>%
  pivot_wider(names_from = 'Timepoint', values_from = 'value')

write.csv(to_meta_C2C3_time4, file = "to_meta_C2C3_time4.csv", row.names=FALSE)

meta_C2C3_time4 <- meta2d(infile = "to_meta_C2C3_time4.csv", outputFile = FALSE, filestyle = "csv", timepoints = "line1", minper=20, nCores = detectCores()-1,
       maxper=28, parallelize = TRUE)$meta %>%
  mutate(Gene = gsub('_.*', "", CycID), 
         Sample = gsub('.*_', "", CycID), 
         .keep='unused')
write.csv(meta_C2C3_time4, file = "meta_C2C3_time4.csv", row.names=FALSE)
meta_C2C3_time4 <- read.csv('meta_C2C3_time4.csv')


to_meta_C2C3_AGG_time4 <- AverageExpression(data_C2C3_AGG, group.by = c("time4", "sample"), layer="data")$RNA %>% as.data.frame %>%
  {mutate(., Gene = rownames(.))} %>%
  pivot_longer(cols=-'Gene', names_to='TP_sample') %>%
  mutate(Timepoint = TP_sample %>% gsub('_.*', '', .) %>% gsub('.*T', '', .) %>% gsub('g', '', .) , 
         Gene_S = paste(Gene, TP_sample %>% gsub('.*_C', 'C', .), sep = '_'), .keep='unused') %>%
  pivot_wider(names_from = 'Timepoint', values_from = 'value')

write.csv(to_meta_C2C3_AGG_time4, file = "to_meta_C2C3_AGG_time4.csv", row.names=FALSE)

meta_C2C3_AGG_time4 <- meta2d(infile = "to_meta_C2C3_AGG_time4.csv", outputFile = FALSE, filestyle = "csv", timepoints = "line1", minper=20, nCores = detectCores()-1,
                          maxper=28, parallelize = TRUE)$meta %>%
  mutate(Gene = gsub('_.*', "", CycID), 
         Sample = gsub('.*_', "", CycID), 
         .keep='unused')
#write.csv(meta_C2C3_time4, file = "meta_C2C3_time4.csv", row.names=FALSE)

meta_C2C3_AGG_time4 <- read.csv('meta_C2C3_AGG_time4.csv')



to_meta_C2C3_AGG_max1_4 <- AverageExpression(data_C2C3_AGG_max1, group.by = c("time4", "sample"), layer="data")$RNA %>% as.data.frame %>%
  {mutate(., Gene = rownames(.))} %>%
  pivot_longer(cols=-'Gene', names_to='TP_sample') %>%
  mutate(Timepoint = TP_sample %>% gsub('_.*', '', .) %>% gsub('.*T', '', .) %>% gsub('g', '', .) , 
         Gene_S = paste(Gene, TP_sample %>% gsub('.*_C', 'C', .), sep = '_'), .keep='unused') %>%
  pivot_wider(names_from = 'Timepoint', values_from = 'value')

write.csv(to_meta_C2C3_AGG_max1_4, file = "to_meta_C2C3_AGG_max1_4.csv", row.names=FALSE)

meta2d(infile = "to_meta_C2C3_AGG_max1_4.csv", filestyle = "csv", timepoints = "line1", minper=20, 
       maxper=28, outdir=getwd(), parallelize = TRUE)

meta_C2C3_AGG_max1_4 <- read.csv('meta2d_to_meta_C2C3_AGG_max1_4.csv') %>%
  mutate(Gene = gsub('_.*', "", CycID), 
         Mode = gsub('.*_', "", CycID), 
         Sample = gsub('.*_C', 'C', CycID) %>% sub('_.*', '', .),
         .keep='unused')
write.csv(meta_C2C3_AGG_max1_4, file = "meta_C2C3_AGG_max1_4.csv", row.names=FALSE)
meta_C2C3_AGG_max1_4 <- read.csv('meta_C2C3_AGG_max1_4.csv')



to_meta_C2C3_AGG_max01_4 <- AverageExpression(data_C2C3_AGG_01_1d, group.by = c("time4", "sample", "mode1"), layer="counts")$RNA %>% as.data.frame %>%
  {mutate(., Gene = rownames(.))} %>%
  pivot_longer(cols=-'Gene', names_to='TP_sample_mode') %>%
  mutate(Timepoint = TP_sample_mode %>% gsub('_.*', '', .) %>% gsub('.*T', '', .) %>% gsub('g', '', .), 
         Gene_SM = paste(Gene, TP_sample_mode %>% gsub('.*_C', 'C', .), sep = '_'), .keep='unused') %>%
  mutate(value = ifelse(is.na(value), 0, value)) %>%
  pivot_wider(names_from = 'Timepoint', values_from = 'value')

write.csv(to_meta_C2C3_AGG_max01_4, file = "to_meta_C2C3_AGG_max01_4.csv", row.names=FALSE)

meta2d(infile = "to_meta_C2C3_AGG_max01_4.csv", filestyle = "csv", timepoints = "line1", minper=20, 
       maxper=28, outdir=getwd(), parallelize = TRUE)

meta_C2C3_AGG_max01_4 <- read.csv('meta2d_to_meta_C2C3_AGG_max01_4.csv') %>%
  mutate(Gene = gsub('_.*', "", CycID), 
         Mode = gsub('.*_', "", CycID), 
         Sample = gsub('.*_C', 'C', CycID) %>% sub('_.*', '', .),
         .keep='unused')
write.csv(meta_C2C3_AGG_max01_4, file = "meta_C2C3_AGG_max01_4.csv", row.names=FALSE)
meta_C2C3_AGG_max01_4 <- read.csv('meta_C2C3_AGG_max01_4.csv')




{to_meta_C2C3_time3 <- AverageExpression(added_C2C3, group.by = c("time3", "sample", "mode1"), layer="data")$RNA %>% as.data.frame %>%
  {mutate(., Gene = rownames(.))} %>%
  pivot_longer(cols=-'Gene', names_to='TP_sample_mode') %>%
  mutate(Timepoint = TP_sample_mode %>% gsub('_.*', '', .) %>% gsub('.*T', '', .) %>% gsub('g', '', .), 
         Gene_SM = paste(Gene, TP_sample_mode %>% gsub('.*_C', 'C', .), sep = '_'), .keep='unused') %>%
  pivot_wider(names_from = 'Timepoint', values_from = 'value')

write.csv(to_meta_C2C3_time3, file = "to_meta_C2C3_time3.csv", row.names=FALSE)

meta2d(infile = "to_meta_C2C3_time3.csv", filestyle = "csv", timepoints = "line1", minper=20, 
       maxper=28, outdir=getwd(), parallelize = TRUE)

meta_C2C3_time3 <- read.csv('meta2d_to_meta_C2C3_time3.csv') %>%
  mutate(Gene = gsub('_.*', "", CycID), 
         Mode = gsub('.*_', "", CycID), 
         Sample = gsub('.*_C', 'C', CycID) %>% sub('_.*', '', .),
         .keep='unused')
write.csv(meta_C2C3_time3, file = "meta_C2C3_time3.csv", row.names=FALSE)
meta_C2C3_time3 <- read.csv('meta_C2C3_time3.csv')
}

{to_meta_C2C3_AGG_time3 <- AverageExpression(added_C2C3_AGG, group.by = c("time3", "sample", "mode1"), layer="data")$RNA %>% as.data.frame %>%
  {mutate(., Gene = rownames(.))} %>%
  pivot_longer(cols=-'Gene', names_to='TP_sample_mode') %>%
  mutate(Timepoint = TP_sample_mode %>% gsub('_.*', '', .) %>% gsub('.*T', '', .) %>% gsub('g', '', .), 
         Gene_SM = paste(Gene, TP_sample_mode %>% gsub('.*_C', 'C', .), sep = '_'), .keep='unused') %>%
  pivot_wider(names_from = 'Timepoint', values_from = 'value')

write.csv(to_meta_C2C3_AGG_time3, file = "to_meta_C2C3_AGG_time3.csv", row.names=FALSE)

meta2d(infile = "to_meta_C2C3_AGG_time3.csv", filestyle = "csv", timepoints = "line1", minper=20, 
       maxper=28, outdir=getwd(), parallelize = TRUE)

meta_C2C3_AGG_time3 <- read.csv('meta2d_to_meta_C2C3_AGG_time3.csv') %>%
  mutate(Gene = gsub('_.*', "", CycID), 
         Mode = gsub('.*_', "", CycID), 
         Sample = gsub('.*_C', 'C', CycID) %>% sub('_.*', '', .),
         .keep='unused')
write.csv(meta_C2C3_AGG_time3, file = "meta_C2C3_AGG_time3.csv", row.names=FALSE)
}
meta_C2C3_AGG_time3 <- read.csv('meta_C2C3_AGG_time3.csv')


{to_meta_C2C3_pool3 <- AverageExpression(added_C2C3_pool3, group.by = c("time3", "sample", "mode1"), layer="data")$RNA %>% as.data.frame %>%
  {mutate(., Gene = rownames(.))} %>%
  pivot_longer(cols=-'Gene', names_to='TP_sample_mode') %>%
  mutate(Timepoint = TP_sample_mode %>% gsub('_.*', '', .) %>% gsub('.*T', '', .) %>% gsub('g', '', .), 
         Gene_SM = paste(Gene, TP_sample_mode %>% gsub('.*_C', 'C', .), sep = '_'), .keep='unused') %>%
  pivot_wider(names_from = 'Timepoint', values_from = 'value')

write.csv(to_meta_C2C3_pool3, file = "to_meta_C2C3_pool3.csv", row.names=FALSE)

meta2d(infile = "to_meta_C2C3_pool3.csv", filestyle = "csv", timepoints = "line1", minper=20, 
       maxper=28, outdir=getwd(), parallelize = TRUE)

meta_C2C3_pool3 <- read.csv('meta2d_to_meta_C2C3_pool3.csv') %>%
  mutate(Gene = gsub('_.*', "", CycID), 
         Mode = gsub('.*_', "", CycID), 
         Sample = gsub('.*_C', 'C', CycID) %>% sub('_.*', '', .),
         .keep='unused')
write.csv(meta_C2C3_pool3, file = "meta_C2C3_pool3.csv", row.names=FALSE)
}
meta_C2C3_pool3 <- read.csv('meta_C2C3_pool3.csv')


{to_meta_C2C3_AGG_pool3 <- AverageExpression(added_C2C3_pool3_AGG, group.by = c("time3", "sample", "mode1"), layer="data")$RNA %>% as.data.frame %>%
  {mutate(., Gene = rownames(.))} %>%
  pivot_longer(cols=-'Gene', names_to='TP_sample_mode') %>%
  mutate(Timepoint = TP_sample_mode %>% gsub('_.*', '', .) %>% gsub('.*T', '', .) %>% gsub('g', '', .), 
         Gene_SM = paste(Gene, TP_sample_mode %>% gsub('.*_C', 'C', .), sep = '_'), .keep='unused') %>%
  pivot_wider(names_from = 'Timepoint', values_from = 'value')

write.csv(to_meta_C2C3_AGG_pool3, file = "to_meta_C2C3_AGG_pool3.csv", row.names=FALSE)

meta2d(infile = "to_meta_C2C3_AGG_pool3.csv", filestyle = "csv", timepoints = "line1", minper=20, 
       maxper=28, outdir=getwd(), parallelize = TRUE)

meta_C2C3_AGG_pool3 <- read.csv('meta2d_to_meta_C2C3_AGG_pool3.csv') %>%
  mutate(Gene = gsub('_.*', "", CycID), 
         Mode = gsub('.*_', "", CycID), 
         Sample = gsub('.*_C', 'C', CycID) %>% sub('_.*', '', .),
         .keep='unused')
write.csv(meta_C2C3_AGG_pool3, file = "meta_C2C3_AGG_pool3.csv", row.names=FALSE)
}
meta_C2C3_AGG_pool3 <- read.csv('meta_C2C3_AGG_pool3.csv')



{to_meta_C2C3_pool <- AverageExpression(added_C2C3_pool, group.by = c("time3", "sample", "mode1"), layer="data")$RNA %>% as.data.frame %>%
  {mutate(., Gene = rownames(.))} %>%
  pivot_longer(cols=-'Gene', names_to='TP_sample_mode') %>%
  mutate(Timepoint = TP_sample_mode %>% gsub('_.*', '', .) %>% gsub('.*T', '', .) %>% gsub('g', '', .), 
         Gene_SM = paste(Gene, TP_sample_mode %>% gsub('.*_C', 'C', .), sep = '_'), .keep='unused') %>%
  pivot_wider(names_from = 'Timepoint', values_from = 'value')

write.csv(to_meta_C2C3_pool, file = "to_meta_C2C3_pool.csv", row.names=FALSE)

meta2d(infile = "to_meta_C2C3_pool.csv", filestyle = "csv", timepoints = "line1", minper=20, 
       maxper=28, outdir=getwd(), parallelize = TRUE)

meta_C2C3_pool <- read.csv('meta2d_to_meta_C2C3_pool.csv') %>%
  mutate(Gene = gsub('_.*', "", CycID), 
         Mode = gsub('.*_', "", CycID), 
         Sample = gsub('.*_C', 'C', CycID) %>% sub('_.*', '', .),
         .keep='unused')
write.csv(meta_C2C3_pool, file = "meta_C2C3_pool.csv", row.names=FALSE)
}
meta_C2C3_pool <- read.csv('meta_C2C3_pool.csv')


{to_meta_C2C3_AGG_pool <- AverageExpression(added_C2C3_pool_AGG, group.by = c("time3", "sample", "mode1"), layer="data")$RNA %>% as.data.frame %>%
  {mutate(., Gene = rownames(.))} %>%
  pivot_longer(cols=-'Gene', names_to='TP_sample_mode') %>%
  mutate(Timepoint = TP_sample_mode %>% gsub('_.*', '', .) %>% gsub('.*T', '', .) %>% gsub('g', '', .), 
         Gene_SM = paste(Gene, TP_sample_mode %>% gsub('.*_C', 'C', .), sep = '_'), .keep='unused') %>%
  pivot_wider(names_from = 'Timepoint', values_from = 'value')

write.csv(to_meta_C2C3_AGG_pool, file = "to_meta_C2C3_AGG_pool.csv", row.names=FALSE)

meta2d(infile = "to_meta_C2C3_AGG_pool.csv", filestyle = "csv", timepoints = "line1", minper=20, 
       maxper=28, outdir=getwd(), parallelize = TRUE)

meta_C2C3_AGG_pool <- read.csv('meta2d_to_meta_C2C3_AGG_pool.csv') %>%
  mutate(Gene = gsub('_.*', "", CycID), 
         Mode = gsub('.*_', "", CycID), 
         Sample = gsub('.*_C', 'C', CycID) %>% sub('_.*', '', .),
         .keep='unused')
write.csv(meta_C2C3_AGG_pool, file = "meta_C2C3_AGG_pool.csv", row.names=FALSE)
}
meta_C2C3_AGG_pool <- read.csv('meta_C2C3_AGG_pool.csv')



reps_adding_pool$res$avg <- lapply(paste0('R', c(0:2)), function(R) {
  hold = AverageExpression(reps_adding_pool[[R]]$avg, group.by = c("time3", "sample", "mode1"), layer="data")$RNA %>% as.data.frame %>%
    {mutate(., Gene = rownames(.))} %>%
    pivot_longer(cols=-'Gene', names_to='TP_sample_mode') %>% 
    mutate(Timepoint = TP_sample_mode %>% gsub('_.*', '', .) %>% gsub('.*T', '', .) %>% gsub('g', '', .), 
           Gene__SM = paste(Gene, TP_sample_mode %>% gsub('.*_C', 'C', .), sep = '__'), .keep='unused')
  
  Summary = hold %>% group_by(Gene__SM) %>% summarise(Ratio = max(value)/min(value)) %>%
    mutate(Ratio = replace(Ratio, is.na(Ratio), 0), 
           CycID = Gene__SM, .keep='unused')
  
  hold %>% 
    pivot_wider(names_from = 'Timepoint', values_from = 'value') %>% 
    {.[ which(apply(., 1, max, na.rm=TRUE) != 0), order(as.numeric(colnames(.)), na.last = FALSE)]} %>% 
    write.csv(file = "reps_adding_pool.csv", row.names=FALSE)
  
  meta2d(infile = "reps_adding_pool.csv", filestyle = "csv", timepoints = "line1",
         outdir=getwd(), parallelize = TRUE, outputFile = FALSE, nCores = detectCores())$meta %>%
    as.data.frame %>%
    merge(Summary) %>%
    mutate(Gene = gsub('_.*', '', CycID), 
           Sample = gsub('.*__', '', CycID) %>% gsub('_.*', '', .),
           Mode = gsub('.*__', '', CycID) %>% gsub('.*_', '', .),
           run = R,
           .keep='unused') %>%
    return()
}) %>% `names<-`(paste0('R', c(0:2)))

reps_adding_pool$res$AGG <- lapply(paste0('R', c(0:2)), function(R) {
  hold = AverageExpression(reps_adding_pool[[R]]$AGG, group.by = c("time3", "sample", "mode1"), layer="data")$RNA %>% as.data.frame %>%
    {mutate(., Gene = rownames(.))} %>%
    pivot_longer(cols=-'Gene', names_to='TP_sample_mode') %>% 
    mutate(Timepoint = TP_sample_mode %>% gsub('_.*', '', .) %>% gsub('.*T', '', .) %>% gsub('g', '', .), 
           Gene__SM = paste(Gene, TP_sample_mode %>% gsub('.*_C', 'C', .), sep = '__'), .keep='unused')
  
  Summary = hold %>% group_by(Gene__SM) %>% summarise(Ratio = max(value)/min(value)) %>%
    mutate(Ratio = replace(Ratio, is.na(Ratio), 0), 
           CycID = Gene__SM, .keep='unused')
  
  hold %>% 
    pivot_wider(names_from = 'Timepoint', values_from = 'value') %>% 
    {.[ which(apply(., 1, max, na.rm=TRUE) != 0), order(as.numeric(colnames(.)), na.last = FALSE)]} %>% 
    write.csv(file = "reps_adding_pool.csv", row.names=FALSE)
  
  meta2d(infile = "reps_adding_pool.csv", filestyle = "csv", timepoints = "line1",
         outdir=getwd(), parallelize = TRUE, outputFile = FALSE, nCores = detectCores())$meta %>%
    as.data.frame %>%
    merge(Summary) %>%
    mutate(Gene = gsub('_.*', '', CycID), 
           Sample = gsub('.*__', '', CycID) %>% gsub('_.*', '', .),
           Mode = gsub('.*__', '', CycID) %>% gsub('.*_', '', .),
           run = R,
           .keep='unused') %>%
    return()
}) %>% `names<-`(paste0('R', c(0:2)))

library(DESeq2)

dds <- AggregateExpression(added_C2C3, group.by = c('sample_mode_time3'), return.seurat = TRUE)

dds <- DESeqDataSetFromMatrix(countData = dds@assays$RNA$counts, colData =dds@meta.data, design = ~ sample_mode_time3)
dds <- estimateSizeFactors(dds)

counts(dds, normalized=TRUE) %>%
  as.data.frame %>%
  mutate(Gene = rownames(.)) %>%
  pivot_longer(cols=-Gene, names_to='cond') %>% 
  mutate(Type = gsub('D-.*', 'D', cond) %>% gsub('-', '_', .),
         time = as.numeric(gsub('.*-', '', cond)), .keep='unused') %>%
  mutate(Gene_Type = paste(Gene, Type)) %>%
  pivot_wider(names_from = time, values_from = value) %>% 
  as.data.frame %>%
  {`rownames<-`(., .$Gene_Type)} %>%
  mutate(Gene_Type=NULL, Gene=NULL, Type=NULL) %>% 
  {.[ which(apply(., 1, max, na.rm=TRUE) != 0), order(as.numeric(colnames(.)), na.last = FALSE)]} %>% 
  write.csv(file = "files_cycling/dds_SMtime3.csv", row.names=TRUE)

meta2d(infile = "dds_SMtime3.csv", filestyle = "csv", timepoints = "line1", minper=20, 
       maxper=28, outdir=getwd(), parallelize = TRUE)

meta_dds_SMtime3 <- read.csv('files_cycling/meta2d_dds_SMtime3.csv') %>%
  mutate(Gene = gsub(' .*', "", CycID), 
         Mode = gsub('.*_', "", CycID), 
         Sample = gsub('.* C', 'C', CycID) %>% sub('_.*', '', .),
         .keep='unused')
write.csv(meta_dds_SMtime3, file = "files_cycling/meta_dds_SMtime3.csv", row.names=FALSE)



DDS <- CreateSeuratObject(counts(dds, normalized=TRUE))

DDS@meta.data[['sample']] <- Cells(DDS) %>% 
  gsub('-.*', '', .) %>% 
  factor(levels=c('C3', 'C2'))
DDS@meta.data[['mode1']] <- Cells(DDS) %>%
  sub('D-.*', 'D', .) %>% 
  sub('.*-', '', .) %>%
  factor(levels=c('LD', 'DD'))
DDS@meta.data[['time3']] <- Cells(DDS) %>%
  gsub('.*-', '', .) %>% 
  factor(levels = seq(3, 71, 4))
DDS@meta.data[['sample_mode']] <- DDS[[c('sample', 'mode1')]] %>% 
  mutate(sample_mode <- paste(sample, mode1, sep='_'), .keep='unused') %>% unlist %>%
  factor(levels=c('C3_LD', 'C2_LD', 'C3_DD', 'C2_DD'))
DDS@meta.data[['time6']] <- DDS[[c('mode1', 'time3')]] %>% 
  mutate(time6 = as.numeric(as.numeric(as.character(time3))+ifelse(mode1=='LD', 0, 48))) %>%
  mutate(time6 = factor(time6, sort(unique(time6)))) %>% {.[, 'time6']}




####    ####    ####    ####    ####    ####    ####    ####    ####    ####   


# 5. PREPARING AMPLITUDES, EdgeR
################################################################################
avg_data_times2_max1 <- AverageExpression(data_C2C3_AGG_max1, group.by = c('sample_mode', 'group2'), layer="counts")$RNA %>% 
  as.data.frame %>%
  {mutate(., Gene = rownames(.))} %>%
  pivot_longer(cols=-'Gene', names_to='sample_mode_time', values_to='avg') %>%
  mutate(Mode = gsub('_.*', '', sample_mode_time) %>% gsub('.*-', '', .), 
         Sample = gsub('-.*', '', sample_mode_time), 
         timepoint = gsub('.*_', '', sample_mode_time) %>% gsub('.*T', '', .),
         sample_mode_time = gsub('_', '-', sample_mode_time))

avg_data_times2_01 <- AverageExpression(data_C2C3_AGG_01_1d, group.by = c('sample_mode', 'group2'), layer="counts")$RNA %>% 
  as.data.frame %>%
  {mutate(., Gene = rownames(.))} %>%
  pivot_longer(cols=-'Gene', names_to='sample_mode_time', values_to='avg') %>%
  mutate(Mode = gsub('_.*', '', sample_mode_time) %>% gsub('.*-', '', .), 
         Sample = gsub('-.*', '', sample_mode_time), 
         timepoint = gsub('.*_', '', sample_mode_time) %>% gsub('.*T', '', .),
         sample_mode_time = gsub('_', '-', sample_mode_time)) 

avg_data_times2_01_2d <- AverageExpression(data_C2C3_AGG_01_2d, group.by = c('sample_mode', 'group2'), layer="counts")$RNA %>% 
  as.data.frame %>%
  {mutate(., Gene = rownames(.))} %>%
  pivot_longer(cols=-'Gene', names_to='sample_mode_time', values_to='avg') %>%
  mutate(Mode = gsub('_.*', '', sample_mode_time) %>% gsub('.*-', '', .), 
         Sample = gsub('-.*', '', sample_mode_time), 
         timepoint = gsub('.*_', '', sample_mode_time) %>% gsub('.*T', '', .),
         sample_mode_time = gsub('_', '-', sample_mode_time))

avg_data_times2 <- AverageExpression(data_C2C3, group.by = c('sample_mode', 'group2'), layer="data")$RNA %>% as.data.frame %>%
  {mutate(., Gene = rownames(.))} %>%
  pivot_longer(cols=-'Gene', names_to='sample_mode_time', values_to='avg') %>%
  mutate(Mode = gsub('_.*', '', sample_mode_time) %>% gsub('.*-', '', .), 
         Sample = gsub('-.*', '', sample_mode_time), 
         timepoint = gsub('.*_', '', sample_mode_time) %>% gsub('.*T', '', .),
         sample_mode_time = gsub('_', '-', sample_mode_time))

avg_data_times2 %>%
  pivot_wider(names_from = 'timepoint', values_from = 'avg')

maxMin_data_times2 <- avg_data_times2 %>%
  {split(., list(.$Gene, .$Sample, .$Mode))} %>%
  lapply(function(df) rbind(cbind(df[which.max(df$avg),], tag = 'max'), 
                            cbind(df[which.min(df$avg),], tag='min'))) %>%
  do.call(what='rbind')

maxMin_data_params2 <- maxMin_data_times2 %>%
  {split(., .$Gene)} %>%
  lapply(function(df) {cbind(Gene = unique(df$Gene), 
                             max = df[df$tag == 'max', 'avg'],
                             min = df[df$tag == 'min', 'avg'])}) %>%
  do.call(what='rbind') %>%
  as.data.frame %>%
  mutate_at(c('max', 'min'), as.numeric) %>%
  mutate(ratio = max/min,
         ampl = max-min,
         rel_ampl = 2*(max-min)/(max+min))



pcts_data_times2 <- SplitObject(data_C2C3, 'sample_mode') %>% lapply(SplitObject, split.by='group2') %>% unlist %>% 
  {lapply(names(.), function(name) PrctCellExpringGene_presplit(.[[name]], genes=rownames(.[[name]])) %>%
            cbind(cond = name) %>%
            mutate(Gene = Markers,
                   Sample_Mode = gsub("\\..*", '', cond),
                   group2 = sub(".*\\.", '', cond),
                   .keep='unused'))} %>%
  do.call(what='rbind') %>%
  as.data.frame

pcts_data_vals2 <- pcts_data_times2 %>%
  {split(., list(.$Gene, .$Sample_Mode))} %>%
  lapply(function(df) {
    df = arrange(df, desc(Cell_proportion))
    cbind(Gene = unique(df$Gene), Sample_Mode = unique(df$Sample_Mode), 
          most_pct = df[1, 'Cell_proportion'], most_time = df[1, 'group2'],
          second_pct = df[2, 'Cell_proportion'], second_time = df[2, 'group2'])
  }) %>%
  do.call(what='rbind') %>%
  as.data.frame %>%
  mutate(Sample = gsub('_.*', '', Sample_Mode), 
         Mode = gsub('.*_', '', Sample_Mode), .keep='unused')

write.csv(pcts_data_vals2, 'pcts_data_vals2.csv', row.names=FALSE)
pcts_data_vals2 <- read.csv('pcts_data_vals2.csv')


{ # RUN EdgeR on genes_atLeast0.02pct, giving edgeR_lists and edgeR_results - scratch that
  bulks <- AggregateExpression(data_C2C3, group.by='sample_mode_group2', features=rownames(data_C2C3), assays='RNA')[[1]]
  bulk_groups = names(as.data.frame(bulks)) %>% gsub('-line.*', '', .)
  d <- DGEList(counts=bulks, group=factor(bulk_groups)) %>% 
    calcNormFactors()
  
  design.mat <- model.matrix(~ 0 + d$samples$group)
  colnames(design.mat) <- levels(d$samples$group)
  d <- estimateGLMCommonDisp(d, design.mat)
  d <- estimateGLMTagwiseDisp(d, design.mat)
  
  {edgeR_tests = list()
    for (t1 in seq(3, 23, 4)) {
      for (t2 in seq(3, 23, 4)) {
        for (s in c('C2', 'C3')) {
          for (c in c('-LD-ZT', '-DD-CT')) {
            cond1 = paste0(s, c, t1)
            cond2 = paste0(s, c, t2)
            edgeR_tests[[cond1]][[cond2]] = exactTest(d, pair=c(which(levels(d$samples$group) == cond1),
                                                                which(levels(d$samples$group) == cond2)))
          }
        }
      }
    }
    rm(t1, t2, s)
  }
  
  edgeR_results <- lapply(names(edgeR_tests), function(c1) {
    lapply(names(edgeR_tests[[c1]]), function(c2) {
      cbind(edgeR_tests[[c1]][[c2]][['table']], c1=c1, c2=c2) %>%
        {cbind(., Gene = rownames(.))}
    }) %>% do.call(what='rbind')
  }) %>% do.call(what='rbind') %>%
    as.data.frame
} # RUN EdgeR on genes_atLeast0.02pct, giving edgeR_lists and edgeR_results


avg_data_times2 <- AverageExpression(data_C2C3, group.by = c('sample_mode', 'group2'), layer="data")$RNA %>% as.data.frame %>%
  {mutate(., Gene = rownames(.))} %>%
  pivot_longer(cols=-'Gene', names_to='sample_mode_time', values_to='avg') %>%
  mutate(Mode = gsub('_.*', '', sample_mode_time) %>% gsub('.*-', '', .), 
         Sample = gsub('-.*', '', sample_mode_time), 
         timepoint = gsub('.*_', '', sample_mode_time) %>% gsub('.*T', '', .),
         sample_mode_time = gsub('_', '-', sample_mode_time))

maxMin_data_times2 <- avg_data_times2 %>%
  {split(., list(.$Gene, .$Sample, .$Mode))} %>%
  lapply(function(df) rbind(cbind(df[which.max(df$avg),], tag = 'max'), 
                            cbind(df[which.min(df$avg),], tag='min'))) %>%
  do.call(what='rbind')

maxMin_data_params2 <- maxMin_data_times2 %>%
  {split(., list(.$Gene, .$Sample, .$Mode))} %>%
  lapply(function(df) {cbind(Gene = unique(df$Gene), 
                             Sample = unique(df$Sample),
                             Mode = unique(df$Mode), 
                             max = df[df$tag == 'max', 'avg'],
                             min = df[df$tag == 'min', 'avg'])}) %>%
  do.call(what='rbind') %>%
  as.data.frame %>%
  mutate_at(c('max', 'min'), as.numeric) %>%
  mutate(ratio = max/min,
         ampl = max-min,
         rel_ampl = 2*(max-min)/(max+min))

#maxMin_data_params2 <- maxMin_data_times2 %>%
#  {split(., list(.$Gene, .$Sample, .$Mode))} %>%
#  lapply(function(df) cbind(Gene = unique(df$Gene), 
#                            c1 = df$sample_mode_time[1], 
#                            c2 = df$sample_mode_time[2])) %>%
#  do.call(what='rbind') %>%
#  as.data.frame %>%
#  merge(edgeR_results) %>%
#  mutate(Sample = gsub('-.*', '', c1), Mode = gsub('D-.*', 'D', c1) %>% gsub('.*-', '', .)) %>%
#  merge(maxMin_data_params2) %>%
#  mutate(FDR_EdgeR = p.adjust(PValue, method='BH'))


avg_data <- AverageExpression(data_C2C3, group.by = c("sample", "mode1"), layer="data")$RNA %>% as.data.frame %>%
  {mutate(., Gene = rownames(.))} %>%
  pivot_longer(cols=-'Gene', names_to='sample_mode', values_to='avg') %>%
  mutate(Sample = gsub('_.*', '', sample_mode),
         Mode = gsub('.*_', '', sample_mode), .keep='unused')

pcts_data <- SplitObject(data_C2C3, 'sample') %>% lapply(SplitObject, split.by='mode1') %>% unlist %>% 
  {lapply(names(.), function(name) PrctCellExpringGene_presplit(.[[name]], genes=rownames(.[[name]])) %>%
            cbind(cond = name) %>%
            mutate(Gene = Markers,
                   Sample = gsub("\\..*", '', cond),
                   Mode = sub(".*\\.", '', cond),
                   .keep='unused'))} %>%
  do.call(what='rbind')

expr_data <- merge(avg_data, pcts_data) %>%
  merge(maxMin_data_params2) %>%
  merge(pcts_data_vals2)




{pcts_data_time2 <- SplitObject(data_C2C3, split.by='sample') %>% 
    lapply(SplitObject, split.by='mode1') %>% unlist %>%
    lapply(SplitObject, split.by='group2') %>% unlist %>%
    {lapply(names(.), function(name) PrctCellExpringGene_presplit(.[[name]], genes=rownames(.[[name]])) %>%
              cbind(cond = name) %>%
              mutate(Gene = Markers,
                     group2 = gsub('.*\\.', '', cond), 
                     Sample = gsub('\\..*', '', cond),
                     Mode = sub('D\\..*', 'D', cond) %>% sub('.*\\.', '', .),
                     .keep='unused'))} %>%
    do.call(what='rbind')} # to make pcts_data_time2
#write.csv(pcts_data_time2, 'pcts_data_time2.csv', row.names=FALSE)
pcts_data_time2 <- read.csv('pcts_data_time2.csv')


# FOR NORMALIZED AGGREGATES

avg_data_times2_AGG <- AverageExpression(data_C2C3_AGG, group.by = c('sample_mode', 'group2'), layer="data")$RNA %>% as.data.frame %>%
  {mutate(., Gene = rownames(.))} %>%
  pivot_longer(cols=-'Gene', names_to='sample_mode_time', values_to='avg') %>%
  mutate(Mode = gsub('_.*', '', sample_mode_time) %>% gsub('.*-', '', .), 
         Sample = gsub('-.*', '', sample_mode_time), 
         timepoint = gsub('.*_', '', sample_mode_time) %>% gsub('.*T', '', .),
         sample_mode_time = gsub('_', '-', sample_mode_time))


maxMin_data_times2_AGG <- avg_data_times2_AGG %>%
  {split(., list(.$Gene, .$Sample, .$Mode))} %>%
  lapply(function(df) rbind(cbind(df[which.max(df$avg),], tag = 'max'), 
                            cbind(df[which.min(df$avg),], tag='min'))) %>%
  do.call(what='rbind')

maxMin_data_params2_AGG <- maxMin_data_times2_AGG %>%
  {split(., list(.$Gene, .$Sample, .$Mode))} %>%
  lapply(function(df) {cbind(Gene = unique(df$Gene), 
                             Sample = unique(df$Sample),
                             Mode = unique(df$Mode), 
                             max = df[df$tag == 'max', 'avg'],
                             min = df[df$tag == 'min', 'avg'])}) %>%
  do.call(what='rbind') %>%
  as.data.frame %>%
  mutate_at(c('max', 'min'), as.numeric) %>%
  mutate(ratio = max/min,
         ampl = max-min,
         rel_ampl = 2*(max-min)/(max+min)) %>%
  {`colnames<-`(., c("Gene", "Sample", "Mode", 
                     paste0(names(.)[which(names(.) %!in% c("Gene", "Sample", "Mode"))], '_AGG')))}

avg_data_times2_AGG %>% 
  mutate(tag = paste0(Gene, '-', Sample, '_', Mode), 
         sample_mode_time = NULL, .keep='unused') %>%
  pivot_wider(names_from='timepoint', values_from='avg') %>%
  mutate(Gene = gsub('-.*', '', tag), 
         Sample_Mode = gsub('.*-', '', tag), .keep='unused') %>%
  mutate(Mode = gsub('.*_', '', Sample_Mode), 
         Sample = gsub('_.*', '', Sample_Mode)) %>% 
  left_join(reconstituting_DFs_10$reconstituted, .) %>% 
  write.csv(file='C2C3_filtered_wAggNormExprs.csv')

avg_data_times2_AGG %>% 
  mutate(tag = paste0(Gene, '-', Sample, '_', Mode), 
         sample_mode_time = NULL, .keep='unused') %>%
  pivot_wider(names_from='timepoint', values_from='avg') %>%
  mutate(Gene = gsub('-.*', '', tag), 
         Sample_Mode = gsub('.*-', '', tag), .keep='unused') %>%
  mutate(Mode = gsub('.*_', '', Sample_Mode), 
         Sample = gsub('_.*', '', Sample_Mode)) %>% 
  left_join(subset(to_filter, second_pct>0.1), .) %>% 
  write.csv(file='C2C3_passQC_wAggNormExprs.csv')

avg_data_ALL_adding_pool <- lapply(names(reps_adding_pool[1:3]), function(name_L2) {
  reps_adding_pool[[name_L2]]$AGG@assays$RNA$data %>%
    expm1 %>%
    as.data.frame %>% 
    {mutate(., Gene = rownames(.))} %>%
    pivot_longer(cols=-Gene, names_to='cond', values_to=name_L2)
}) %>%
  Reduce(f=merge) %>%
  mutate(
    Time = gsub('.*-', '', cond), 
    Sample = gsub('-.*', '', cond),
    Mode = sub('D-.*', 'D', cond) %>% sub('.*-', '', .),
    .keep='unused') %>%
  pivot_longer(cols=colnames(.)[gsub("[0-9]", "", colnames(.)) == 'R'], 
               names_to='Rep')

avg_data_ALL_adding_pool_1d_01 <- avg_data_ALL_adding_pool %>%
  mutate(Time = as.numeric(Time)) %>%
  mutate(Day = Time%/%24) %>%
  group_by(Sample, Mode, Gene, Time) %>%
  mutate(TimeMean = mean(value)) %>%
  group_by(Day, Sample, Mode, Gene) %>%
  reframe(value = (value-min(TimeMean))/max(TimeMean-min(TimeMean)), Time=Time, Rep=Rep)

avg_data_ALL_adding_pool_3d_01 <- avg_data_ALL_adding_pool %>%
  mutate(Time = as.numeric(Time)) %>%
  group_by(Sample, Mode, Gene, Time) %>%
  mutate(TimeMean = mean(value)) %>%
  group_by(Sample, Mode, Gene) %>%
  reframe(value = (value-min(TimeMean))/max(TimeMean-min(TimeMean)), Time=Time, Rep=Rep)
  
avg_data_ALL_adding_pool_3d_1 <- avg_data_ALL_adding_pool %>%
  mutate(Time = as.numeric(Time)) %>%
  group_by(Sample, Mode, Gene, Time) %>%
  mutate(TimeMean = mean(value)) %>%
  group_by(Sample, Mode, Gene) %>%
  reframe(value = value/max(TimeMean), Time=Time, Rep=Rep)


avg_data_ALL_adding_pool_6d_01 <- avg_data_ALL_adding_pool %>%
  mutate(Time = as.numeric(Time)+72*as.numeric(Mode=='DD')) %>%
  group_by(Sample, Mode, Gene, Time) %>%
  mutate(TimeMean = mean(value)) %>%
  group_by(Sample, Gene) %>%
  reframe(value = (value-min(TimeMean))/max(TimeMean-min(TimeMean)), Time=Time, Rep=Rep, Mode=Mode)

avg_data_ALL_adding_pool_6d_1 <- avg_data_ALL_adding_pool %>%
  mutate(Time = as.numeric(Time)+72*as.numeric(Mode=='DD')) %>%
  group_by(Sample, Mode, Gene, Time) %>%
  mutate(TimeMean = mean(value)) %>%
  group_by(Sample, Gene) %>%
  reframe(value = value/max(TimeMean), Time=Time, Rep=Rep, Mode=Mode)



avg_data_ALL_adding_pool_1dx9_01 <- avg_data_ALL_adding_pool %>%
  mutate(Time = as.numeric(Time)) %>%
  mutate(Day = Time%/%24, 
         Time = Time%%24) %>%
  group_by(Day, Rep, Sample, Mode, Gene) %>%
  mutate(value = (value-min(value))/(max(value)-min(value)))

saveRDS(avg_data_ALL_adding_pool_1d_01, 'avg_data_ALL_adding_pool_1d_01.rds')
saveRDS(avg_data_ALL_adding_pool_3d_01, 'avg_data_ALL_adding_pool_3d_01.rds')
saveRDS(avg_data_ALL_adding_pool_3d_1, 'avg_data_ALL_adding_pool_3d_1.rds')
saveRDS(avg_data_ALL_adding_pool_6d_01, 'avg_data_ALL_adding_pool_6d_01.rds')
saveRDS(avg_data_ALL_adding_pool_6d_1, 'avg_data_ALL_adding_pool_6d_1.rds')
saveRDS(avg_data_ALL_adding_pool_1dx9_01, 'avg_data_ALL_adding_pool_1dx9_01.rds')


####    ####    ####    ####    ####    ####    ####    ####    ####    ####   

avg_data_ALL_adding_pool <- readRDS('avg_data_ALL_adding_pool.rds')
avg_data_ALL_adding_pool_1d_01 <- readRDS('avg_data_ALL_adding_pool_1d_01.rds')
avg_data_ALL_adding_pool_3d_01 <- readRDS('avg_data_ALL_adding_pool_3d_01.rds')
avg_data_ALL_adding_pool_3d_1 <- readRDS('avg_data_ALL_adding_pool_3d_1.rds')
avg_data_ALL_adding_pool_6d_01 <- readRDS('avg_data_ALL_adding_pool_6d_01.rds')
avg_data_ALL_adding_pool_6d_1 <- readRDS('avg_data_ALL_adding_pool_6d_1.rds')
avg_data_ALL_adding_pool_1dx9_01 <- readRDS('avg_data_ALL_adding_pool_1dx9_01.rds')

#reps_adding_pool$res = readRDS('res.rds')

reps_adding_pool <- readRDS('reps_adding_pool.rds')

meta_C2C3_time3 <- read.csv('meta_C2C3_time3.csv')
meta_C2C3_AGG_time3 <- read.csv('meta_C2C3_AGG_time3.csv')
meta_C2C3_pool3 <- read.csv('meta_C2C3_pool3.csv')
meta_C2C3_AGG_pool3 <- read.csv('meta_C2C3_AGG_pool3.csv')
meta_C2C3_pool <- read.csv('meta_C2C3_pool.csv')
meta_C2C3_AGG_pool <- read.csv('meta_C2C3_AGG_pool.csv')

meta_dds_SMtime3 <- read.csv('files_cycling/meta_dds_SMtime3.csv')


to_filter <- merge(KW_scores_C2C3, expr_data[, c('Gene', 'Sample', 'Mode', 'max', 'ampl', 'ratio', 'Cell_proportion')]) %>% 
  merge(pcts_data_vals2[, c('Gene', 'Sample', 'Mode', 'most_pct', 'second_pct')]) %>%
  merge(maxMin_data_params2_AGG[, c('Gene', 'Sample', 'Mode', 'max_AGG', 'ampl_AGG', 'ratio_AGG')]) %>%
  merge(meta_C2C3_group2[, c('Gene', 'Sample', 'Mode', 'JTK_pvalue', 
                             'JTK_adjphase', 'meta2d_pvalue', 'meta2d_phase')]) %>% 
  merge(meta_C2C3_AGG_group2[, c('Gene', 'Sample', 'Mode', 'JTK_pvalue', 
                                         'JTK_adjphase', 'meta2d_pvalue', 'meta2d_phase')] %>% 
          rename('JTK_pvalue_AGG_norm' = 'JTK_pvalue') %>%
          rename('JTK_adjphase_AGG_norm' = 'JTK_adjphase') %>%
          rename('meta2d_pvalue_AGG_norm' = 'meta2d_pvalue') %>%
          rename('meta2d_phase_AGG_norm' = 'meta2d_phase')) %>% 
  # merge(meta_C2C3_AGG_CLRnorm_group2[, c('Gene', 'Sample', 'Mode', 'JTK_pvalue', 
  #                                        'JTK_adjphase', 'meta2d_pvalue', 'meta2d_phase')] %>%
  #         rename('JTK_pvalue_AGG_CLRnorm' = 'JTK_pvalue') %>%
  #         rename('JTK_adjphase_AGG_CLRnorm' = 'JTK_adjphase') %>%
  #         rename('meta2d_pvalue_AGG_CLRnorm' = 'meta2d_pvalue') %>%
  #         rename('meta2d_phase_AGG_CLRnorm' = 'meta2d_phase')) %>% 
  # merge(meta_MC_norm_group2[, c('Gene', 'Sample', 'Mode', 'JTK_pvalue', 
  #                               'JTK_adjphase', 'meta2d_pvalue', 'meta2d_phase')] %>%
  #         rename('JTK_pvalue_MCNorm' = 'JTK_pvalue') %>%
  #         rename('JTK_adjphase_MCNorm' = 'JTK_adjphase') %>%
  #         rename('meta2d_pvalue_MCNorm' = 'meta2d_pvalue') %>%
  #         rename('meta2d_phase_MCNorm' = 'meta2d_phase')) %>% 
  merge(meta_C2C3_AGG_max1[, c('Gene', 'Sample', 'Mode', 'JTK_pvalue',
                                'JTK_adjphase', 'meta2d_pvalue', 'meta2d_phase')] %>%
          rename('JTK_pvalue_norm1d' = 'JTK_pvalue') %>%
          rename('JTK_adjphase_norm1d' = 'JTK_adjphase') %>%
          rename('meta2d_pvalue_norm1d' = 'meta2d_pvalue') %>%
          rename('meta2d_phase_norm1d' = 'meta2d_phase')) %>%
  #left_join(meta_C2C3_time4[, c('Gene', 'Sample', 'ARS_pvalue', 'JTK_pvalue', 'LS_pvalue', 'meta2d_pvalue')] %>% 
  #            {`colnames<-`(., gsub('pvalue', 'pvalue_avg_4', colnames(.)))}) %>%
  #left_join(meta_C2C3_AGG_time4[, c('Gene', 'Sample', 'ARS_pvalue', 'JTK_pvalue', 'LS_pvalue', 'meta2d_pvalue')] %>% 
  #            {`colnames<-`(., gsub('pvalue', 'pvalue_AGG_4', colnames(.)))}) %>% 
  
  merge(meta_C2C3_time3[, c('Gene', 'Sample', 'Mode', 'ARS_pvalue', 'JTK_pvalue', 'LS_pvalue', 'meta2d_pvalue')] %>% 
              {`colnames<-`(., gsub('pvalue', 'pvalue_avg_3', colnames(.)))}) %>% 
  merge(meta_C2C3_AGG_time3[, c('Gene', 'Sample', 'Mode', 'ARS_pvalue', 'JTK_pvalue', 'LS_pvalue', 'meta2d_pvalue')] %>% 
              {`colnames<-`(., gsub('pvalue', 'pvalue_AGG_3', colnames(.)))}) %>%
  
  merge(meta_C2C3_pool3[, c('Gene', 'Sample', 'Mode', 'ARS_pvalue', 'JTK_pvalue', 'LS_pvalue', 'meta2d_pvalue')] %>% 
          {`colnames<-`(., gsub('pvalue', 'pvalue_avg_pool3', colnames(.)))}) %>% 
  merge(meta_C2C3_AGG_pool3[, c('Gene', 'Sample', 'Mode', 'ARS_pvalue', 'JTK_pvalue', 'LS_pvalue', 'meta2d_pvalue')] %>% 
          {`colnames<-`(., gsub('pvalue', 'pvalue_AGG_pool3', colnames(.)))}) %>%
  
  merge(meta_C2C3_pool[, c('Gene', 'Sample', 'Mode', 'ARS_pvalue', 'JTK_pvalue', 'LS_pvalue', 'meta2d_pvalue')] %>% 
          {`colnames<-`(., gsub('pvalue', 'pvalue_avg_pool', colnames(.)))}) %>% 
  merge(meta_C2C3_AGG_pool[, c('Gene', 'Sample', 'Mode', 'ARS_pvalue', 'JTK_pvalue', 'LS_pvalue', 'meta2d_pvalue')] %>% 
          {`colnames<-`(., gsub('pvalue', 'pvalue_AGG_pool', colnames(.)))}) %>%
  
  #merge(meta_dds_SMtime3[, c('Gene', 'Sample', 'Mode', 'ARS_pvalue', 'JTK_pvalue', 'LS_pvalue', 'meta2d_pvalue')] %>% 
  #        {`colnames<-`(., gsub('pvalue', 'pvalue_avg_DESeq_time3', colnames(.)))}) %>% 
  
  {mutate_at(., colnames(.[1, ])[!grepl("[[:alpha:]]", .[1, ]) | grepl("Inf", .[1, ])], as.numeric)}





# 6. New filterings
################################################################################

new_filter <- merge(KW_scores_C2C3, expr_data[, c('Gene', 'Sample', 'Mode', 'max', 'min', 'ampl', 'ratio', 'Cell_proportion')]) %>% 
  merge(pcts_data_vals2[, c('Gene', 'Sample', 'Mode', 'most_pct', 'second_pct')]) %>%
  merge(maxMin_data_params2_AGG[, c('Gene', 'Sample', 'Mode', 'max_AGG', 'min_AGG', 'ampl_AGG', 'ratio_AGG')]) %>%
  merge(meta_C2C3_group2[, c('ARS_pvalue', 'JTK_pvalue', 'LS_pvalue', 'meta2d_pvalue', 'Gene', 'Mode', 'Sample')] %>% 
          {`colnames<-`(., gsub('pvalue', 'pvalue_avg', colnames(.)))}) %>%
  merge(meta_C2C3_AGG_group2[, c('ARS_pvalue', 'JTK_pvalue', 'LS_pvalue', 'meta2d_pvalue', 'Gene', 'Mode', 'Sample')] %>% 
          {`colnames<-`(., gsub('pvalue', 'pvalue_AGG', colnames(.)))}) %>%
  left_join(meta_C2C3_time4[, c('ARS_pvalue', 'JTK_pvalue', 'LS_pvalue', 'meta2d_pvalue', 'Gene', 'Sample')] %>% 
          {`colnames<-`(., gsub('pvalue', 'pvalue_avg_4', colnames(.)))}) %>%
  left_join(meta_C2C3_AGG_time4[, c('ARS_pvalue', 'JTK_pvalue', 'LS_pvalue', 'meta2d_pvalue', 'Gene', 'Sample')] %>% 
          {`colnames<-`(., gsub('pvalue', 'pvalue_AGG_4', colnames(.)))}) %>%
  {mutate_at(., colnames(.[1, ])[!grepl("[[:alpha:]]", .[1, ]) | grepl("Inf", .[1, ])], as.numeric)}



lapply(c(paste0(c('ARS_pvalue', 'JTK_pvalue', 'LS_pvalue', 'meta2d_pvalue'), '_avg'), 
       paste0(c('ARS_pvalue', 'JTK_pvalue', 'LS_pvalue', 'meta2d_pvalue'), '_AGG')), 
       function(Col) BH_procedure(subset(new_filter, second_pct>0.1), Q=0.1, Pcol=Col, Qcol=NULL, Qname=paste0(Col, '_BH.Q')) %>% counter) %>%
  do.call(what=rbind) %>%
  mutate(method = c(paste0(c('ARS_pvalue', 'JTK_pvalue', 'LS_pvalue', 'meta2d_pvalue'), '_avg'), 
                    paste0(c('ARS_pvalue', 'JTK_pvalue', 'LS_pvalue', 'meta2d_pvalue'), '_AGG')))

lapply(c(paste0(c('ARS_pvalue', 'JTK_pvalue', 'LS_pvalue', 'meta2d_pvalue'), '_avg'), 
         paste0(c('ARS_pvalue', 'JTK_pvalue', 'LS_pvalue', 'meta2d_pvalue'), '_AGG')), 
       function(Col) BH_procedure(subset(new_filter, second_pct>0.1 & ratio>1.5), Q=0.1, Pcol=Col, Qcol=NULL, Qname=paste0(Col, '_BH.Q')) %>% counter) %>%
  do.call(what=rbind) %>%
  mutate(method = c(paste0(c('ARS_pvalue', 'JTK_pvalue', 'LS_pvalue', 'meta2d_pvalue'), '_avg'), 
                    paste0(c('ARS_pvalue', 'JTK_pvalue', 'LS_pvalue', 'meta2d_pvalue'), '_AGG')))


lapply(c(paste0(c('ARS_pvalue', 'JTK_pvalue', 'LS_pvalue', 'meta2d_pvalue'), '_avg'), 
         paste0(c('ARS_pvalue', 'JTK_pvalue', 'LS_pvalue', 'meta2d_pvalue'), '_AGG')), 
       function(Col) BH_procedure(subset(new_filter, second_pct>0.1 & KW_scores_norm<0.05), Q=0.1, Pcol=Col, Qcol=NULL, Qname=paste0(Col, '_BH.Q')) %>% counter) %>%
  rbindlist(fill=TRUE) %>%
  mutate(method = c(paste0(c('ARS_pvalue', 'JTK_pvalue', 'LS_pvalue', 'meta2d_pvalue'), '_avg'), 
                    paste0(c('ARS_pvalue', 'JTK_pvalue', 'LS_pvalue', 'meta2d_pvalue'), '_AGG')))

lapply(c(paste0(c('ARS_pvalue', 'JTK_pvalue', 'LS_pvalue', 'meta2d_pvalue'), '_avg'), 
         paste0(c('ARS_pvalue', 'JTK_pvalue', 'LS_pvalue', 'meta2d_pvalue'), '_AGG')), 
       function(Col) BH_procedure(subset(new_filter, KW_scores_norm<0.05), Q=0.1, Pcol=Col, Qcol=NULL, Qname=paste0(Col, '_BH.Q')) %>% counter) %>%
  rbindlist(fill=TRUE) %>%
  mutate(method = c(paste0(c('ARS_pvalue', 'JTK_pvalue', 'LS_pvalue', 'meta2d_pvalue'), '_avg'), 
                    paste0(c('ARS_pvalue', 'JTK_pvalue', 'LS_pvalue', 'meta2d_pvalue'), '_AGG')))

o <- rbind(
  cbind(Method='avg', 
        BH_procedure(subset(new_filter, second_pct>0.1 & ratio>1.5), 
                     Q=0.1, 
                     Pcol='meta2d_pvalue_avg', 
                     Qcol=NULL, 
                     Qname=paste0('meta2d_pvalue_avg_BH.Q')) %>% 
          subset(select=c(Gene, Sample, Mode))), 
  cbind(Method='AGG', 
        BH_procedure(subset(new_filter, second_pct>0.1 & ratio_AGG>1.5), 
                     Q=0.1, 
                     Pcol='meta2d_pvalue_AGG', 
                     Qcol=NULL, 
                     Qname=paste0('meta2d_pvalue_AGG_BH.Q')) %>% 
          subset(select=c(Gene, Sample, Mode))
  )
)

subset(o, select=c(Gene, Sample, Mode)) %>%
  unique %>%
  {table(.$Sample, .$Mode)}



filtered_new <- subset(new_filter, second_pct>0.1 & ratio_AGG>1.5, select=c(Gene, Sample, Mode, meta2d_pvalue_avg, meta2d_pvalue_AGG)) %>%
  pivot_longer(cols = c(meta2d_pvalue_avg, meta2d_pvalue_AGG), names_to='Method', values_to='pvalue') %>%
  BH_procedure(Q=0.1, Pcol='pvalue', Qcol=NULL)
filtered_new %>%
  {table(.$Sample, .$Mode)}
filtered_new %>%
  subset(select=c(Gene, Sample, Mode)) %>%
  {table(.$Sample, .$Mode)}

  group_by(Gene, Sample, Mode) %>% summarise(N = n()) %>%
  subset(N == 2) %>%
  {table(.$Sample, .$Mode)}


lapply(levels(data_C2C3$sample_mode), function(S_M) {
  S = gsub('_.*', '', S_M) 
  M = gsub('.*_', '', S_M)
  
  G = subset(filtered_new, Sample == S & Mode == M)$Gene %>% unique
  subplotter_concat(names = G, 
                    data=data_C2C3_AGG[, data_C2C3_AGG$sample_mode == S_M], 
                    N_days=2, n_plots=4, ANNOTATE=TRUE) %>%
    ggsave(filename=paste0('refilt_', S_M, '.pdf'), width=12, height=length(G)/3)
}
)



HMQ <- subset(new_filter, second_pct>0.1 & ratio_AGG>1.5, select=c(Gene, Sample, Mode, meta2d_pvalue_avg, meta2d_pvalue_AGG, ampl, max)) %>%
  cbind(., HMP = apply(., 1, function(L) hmp.stat(c(L[['meta2d_pvalue_avg']], L[['meta2d_pvalue_AGG']], w=NULL)))) %>% 
  BH_procedure(Q=0.1, Pcol='HMP', Qcol=NULL)

  #subset(select=c(Gene, Sample, Mode)) %>%
  #unique %>%
  #{table(.$Sample, .$Mode)}

lapply(levels(data_C2C3$sample_mode), function(S_M) {
  S = gsub('_.*', '', S_M) 
  M = gsub('.*_', '', S_M)
  
  G = subset(HMQ, Sample == S & Mode == M)$Gene %>% unique
  subplotter_concat(names = G, 
                    data=data_C2C3_AGG[, data_C2C3_AGG$sample_mode == S_M], 
                    N_days=2, n_plots=4, ANNOTATE=TRUE) %>%
    ggsave(filename=paste0('HMQ_', S_M, '.pdf'), width=12, height=length(G)/3)
}
)



lapply(levels(data_C2C3$sample_mode), function(S_M) {
  S = gsub('_.*', '', S_M) 
  M = gsub('.*_', '', S_M)
  
  Gn = subset(HMQ, Sample == S & Mode == M)$Gene %>% unique
  Go = subset(reconstituting_DFs_10$reconstituted, Sample == S & Mode == M)$Gene %>% unique
  Gc = Go[Go %!in% Gn]
  
  subplotter_concat(names = Gc, 
                    data=data_C2C3_AGG[, data_C2C3_AGG$sample_mode == S_M], 
                    N_days=2, n_plots=4, ANNOTATE=TRUE) %>%
    ggsave(filename=paste0('Cutting_', S_M, '.pdf'), width=12, height=length(Gc)/3, limitsize=FALSE)
}
)

library(plotROC)

subset(new_filter, second_pct>0.1 & ratio_AGG>1.5, select=c(Gene, Sample, Mode, meta2d_pvalue_avg, meta2d_pvalue_AGG, ampl, max, min, ratio, ampl_AGG, max_AGG, min_AGG, ratio_AGG, KW_scores_norm)) %>%
  cbind(., HMP = apply(., 1, function(L) hmp.stat(c(L[['meta2d_pvalue_avg']], L[['meta2d_pvalue_AGG']], w=NULL)))) %>% 
  BH_procedure(Q=0.1, Pcol='HMP', Qcol=NULL, CUT = FALSE) %>% 
  pivot_longer(cols = c(ampl, max, min, ratio, ampl_AGG, max_AGG, min_AGG, ratio_AGG, KW_scores_norm), names_to='Method', values_to='value') %>%
  ggplot(aes(d=Qpass, m=value, color=Method)) +
  geom_roc()

hi <- subset(new_filter, second_pct>0.1 & ratio_AGG>1.5, select=c(Gene, Sample, Mode, meta2d_pvalue_avg, meta2d_pvalue_AGG, ampl, max, min, ratio, ampl_AGG, max_AGG, min_AGG, ratio_AGG, KW_scores_norm)) %>%
  cbind(., HMP = apply(., 1, function(L) hmp.stat(c(L[['meta2d_pvalue_avg']], L[['meta2d_pvalue_AGG']], w=NULL)))) %>% 
  BH_procedure(Q=0.1, Pcol='HMP', Qcol=NULL, CUT = TRUE)


lapply(c(1, 0.5, 0.1, 0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001), 
       function(L) subset(new_filter, second_pct>0.1 & ratio_AGG>1.5 & KW_scores_norm<=L) %>%
         cbind(., HMP = apply(., 1, function(L) hmp.stat(c(L[['meta2d_pvalue_avg']], L[['meta2d_pvalue_AGG']], w=NULL)))) %>% 
         BH_procedure(Q=0.1, Pcol='HMP', Qcol=NULL, CUT = TRUE) %>%
         counter) %>%
  do.call(what=rbind) %>%
  as.data.frame %>%
  `rownames<-`(c(1, 0.5, 0.1, 0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001))

lapply(seq(0, 1, 0.05), 
       function(L) subset(new_filter, second_pct>0.1 & ratio_AGG>1.5 & min_AGG>=L) %>%
         cbind(., HMP = apply(., 1, function(L) hmp.stat(c(L[['meta2d_pvalue_avg']], L[['meta2d_pvalue_AGG']], w=NULL)))) %>% 
         BH_procedure(Q=0.1, Pcol='HMP', Qcol=NULL, CUT = TRUE) %>%
         counter) %>%
  do.call(what=rbind) %>%
  as.data.frame %>%
  `rownames<-`(seq(0, 1, 0.05))

lapply(seq(0, 1, 0.05), 
       function(L) subset(new_filter, second_pct>0.1 & ratio_AGG>1.5 & min_AGG>=L) %>%
                      cbind(., HMP = apply(., 1, function(L) hmp.stat(c(L[['meta2d_pvalue_avg']], L[['meta2d_pvalue_AGG']], w=NULL)))) %>% 
                      BH_procedure(Q=0.1, Pcol='HMP', Qcol=NULL, CUT = TRUE) %>% 
                      {.[grepl('Rp', .$Gene), ]} %>% counter) %>%
         do.call(what=rbind) %>%
         as.data.frame %>%
         `rownames<-`(seq(0, 1, 0.05))

lapply(levels(data_C2C3$sample), function(Sa) {
  
  G = subset(meta_C2C3_time4, Sample == Sa & meta2d_BH.Q < 0.1)$Gene %>% unique %>% sort %>%
    {factor(., levels=.)}
  make_plot(Genes=G, S=data_C2C3, Type0 = Sa, Group = 'time4') %>%
    ggsave(filename=paste0('refilt_', Sa, '.pdf'), width=4, height=length(G)*0.8, limitsize=FALSE)
}
)

HMQ_4 <- subset(new_filter, second_pct>0.1 & ratio_AGG>1.5, select=c(Gene, Sample, meta2d_pvalue_avg_4, meta2d_pvalue_AGG_4)) %>% unique %>% 
  cbind(., HMP = apply(., 1, function(L) hmp.stat(c(L[['meta2d_pvalue_avg_4']], L[['meta2d_pvalue_AGG_4']], w=NULL)))) %>% 
  BH_procedure(Q=0.1, Pcol='HMP', Qcol=NULL, CUT=FALSE)

lapply(levels(data_C2C3$sample), function(Sa) {
  
  G = subset(HMQ_4, Sample==Sa & Qpass)$Gene %>% unique %>% sort %>%
    {factor(., levels=.)}
  make_plot(Genes=G, S=data_C2C3, Type0 = Sa, Group = 'time4') %>%
    ggsave(filename=paste0('HMQ4', Sa, '.pdf'), width=4, height=length(G)*0.8, limitsize=FALSE)
}
)

subset(HMQ_4, Qpass)$Sample %>% table

HMQ_4_comp <- HMQ_4[, c('Gene', 'Sample', 'HMP', 'BH.Q')] %>%
  pivot_wider(names_from=Sample, values_from=BH.Q, 1) %>% 
  mutate(C3divC2 = C3/C2, 
         C3xC2 = C3*C2)
arrange(HMQ_4_comp, C3divC2) %>% head
  
plot_grid(ncol=2, rel_heights=c(1, 1, 2),
          make_plot(head(arrange(HMQ_4_comp, C3divC2), 3)$Gene, S = data_C2C3_AGG, Type0='C3', ZERO_Y=TRUE), 
          make_plot(head(arrange(HMQ_4_comp, C3divC2), 3)$Gene, S = data_C2C3_AGG, Type0='C2', ZERO_Y=TRUE), 
          
          make_plot(head(arrange(HMQ_4_comp, desc(C3divC2)), 3)$Gene, S = data_C2C3_AGG, Type0='C3', ZERO_Y=TRUE), 
          make_plot(head(arrange(HMQ_4_comp, desc(C3divC2)), 3)$Gene, S = data_C2C3_AGG, Type0='C2', ZERO_Y=TRUE), 
          
          make_plot(head(arrange(HMQ_4_comp, C3xC2), 6)$Gene, S = data_C2C3_AGG, Type0='C3', ZERO_Y=TRUE), 
          make_plot(head(arrange(HMQ_4_comp, C3xC2), 6)$Gene, S = data_C2C3_AGG, Type0='C2', ZERO_Y=TRUE))


HMQ_4[, c('Gene', 'Sample', 'HMP', 'BH.Q')] %>%
  pivot_wider(names_from=Sample, values_from=BH.Q, 1) %>% 
  ggplot(aes(C3, C2)) +
  geom_point() +
  scale_x_continuous(trans='log10') +
  scale_y_continuous(trans='log10')

subset(meta_C2C3_time3, meta2d_BH.Q < 0.05) %>% {table(.$Sample, .$Mode)}

subset(meta_C2C3_time3, JTK_BH.Q < 0.5) %>% {table(.$Sample, .$Mode)}

subset(meta_C2C3_time3, meta2d_BH.Q < 0.05) %>% {.[grep('^Rp|^mRp', .$Gene), ]} %>%  
  mutate(sample_mode = paste(Sample, Mode, sep='_')) %>%
  {table(.$sample_mode, .$Gene)}

subset(meta_C2C3_time3, meta2d_BH.Q < 0.05 & Mode == 'DD') %>% {.[grep('^Rp|^mRp', .$Gene), ]} %>%  
  {table(.$Sample)} 


# adding 3rd day distribution
dist = rbind(cbind(sample = 'three', 
                   data = subset(added_C2C3, time3 %in% seq(51, 71, 4))@assays$RNA@data %>%
                    as.data.frame %>%
                    pivot_longer(cols=everything(), names_to='sample', values_to='data') %>% {.$data}),
             cbind(sample = 'two', 
                   data = data_C2C3@assays$RNA@data %>%
                     as.data.frame %>%
                     pivot_longer(cols=everything(), names_to='sample', values_to='data') %>% {.$data})) %>%
  as.data.frame %>%
  mutate(data = as.numeric(data))

plot_grid(ncol=1, 
          ggplot(dist, aes(data, fill=sample)) +
            geom_density(alpha=0.2) +
            scale_x_log10(),
          ggplot(dist, aes(data, fill=sample)) +
            geom_density(alpha=0.2) +
            scale_x_log10() +
            facet_wrap(~sample)) %>%
  ggsave(filename='time3/dist_t2t3.pdf')

rm(dist)


####    ####    ####    ####    ####    ####    ####    ####    ####    ####    


# 7. ASSESSING RBPs+TFs BY EYE - initialize lists
################################################################################
all_TFs_COs_w_child_gomf <- select(org.Dm.eg.db,
                                   keys=c("GO:0003712", GOMFOFFSPRING[["GO:0003712"]],
                                          "GO:0003700", GOMFOFFSPRING[["GO:0003700"]]),
                                   columns=c("SYMBOL", "ENSEMBL"),
                                   keytype="GO") #725, per, Clk, vri, cwo, Pdp1 (no tim, cry)

all_RBPs_w_child_gomf <- select(org.Dm.eg.db,
                                keys=c("GO:0003723", GOMFOFFSPRING[["GO:0003723"]]),                                
                                columns=c("SYMBOL", "ENSEMBL"),
                                keytype="GO") #1098

found_TFRBPs <- c(all_TFs_COs_w_child_gomf$SYMBOL, all_RBPs_w_child_gomf$SYMBOL) %>%
  {.[. %in% rownames(data_C2C3)[rowSums(data_C2C3[["RNA"]]$counts)!=0]]}

subset(to_filter, KW_scores_norm != 1 & Gene %in% all_TFs_COs_w_child_gomf$SYMBOL)$Gene %>% unique %>%
  pdf_gene_lists(name_base = 'all_TFs_COs', rib_alpha=0.4)

subset(to_filter, KW_scores_norm != 1 & Gene %in% all_RBPs_w_child_gomf$SYMBOL)$Gene %>% unique %>%
  pdf_gene_lists(name_base = 'all_RBPs', rib_alpha=0.4)
pdf_combine_and_delete("all_RBPs")

pdf_combine_and_delete("all_TFs_COs")
pdf_combine_and_delete("all_RBPs")

vetted <-  list()
vetted$LD$C3$plausible$TFs <- c('ab', 'alien', 'arm', 'Pdp1', 'pho', 'Pur-alpha', 'REPTOR-BP', 'Smr', 'SoxN', 'sr', 'twi', 'Bx', 'cbt', 
                         'CG31388', 'CG8944', 'cic', 'crc', 'CrebA', 'CrebB', 'cwo', 'CycC', 'Eip75B', 'FoxP', 'fru', 'Hr38', 'kay')
vetted$LD$C3$plausible$RBPs <- c('Pur-alpha', 'RpL10Ab', 'RpL13', 'RpL13A', 'RpL19', 'RpL23', 'RpL23A', 'RpL26', 'RpL30', 'RpL9',
                                 'RpS10b', 'RpS13', 'RpS16', 'RpS27', 'RpS4', 'RpS5a', 'RpS9', 'smg', 'Smg5', 'Spf45', 'Zcchc7',
                                 'CG10077', 'CG14100', 'Cpsf6', 'fne', 'Larp7', 'mRpS18C', 'mRpS6', 'msi', 'nonA')

subplotter_concat(vetted$LD$C3$plausible$TFs, N_days=2, func = 'avg_plotter_concat_mode')
subplotter_concat(vetted$LD$C3$plausible$RBPs, N_days=2, func = 'avg_plotter_concat_mode')

vetted$LD$C3$strong$TFs <- c('ab', 'Pdp1', 'pho', 'Pur-alpha', 'SoxN', 'sr',
                             'twi', 'cbt', 'cic', 'crc', 'CrebA', 'cwo', 'fru', 'Hr38')
vetted$LD$C3$strong$RBPs <- c('Pur-alpha', 'RpL10Ab', 'RpL13', 'RpL13A', 'RpL19', 'RpL23', 
                              'RpL26', 'RpL30','RpS10b', 'RpS13', 'RpS4', 'RpS5a', 'RpS9')

subplotter_concat(c(vetted$LD$C3$strong$RBPs, vetted$LD$C3$strong$TFs), N_days=2, func = 'avg_plotter_concat_mode')
####    ####    ####    ####    ####    ####    ####    ####    ####    ####   

# 8. ROC
################################################################################

{param_test <- list()
param_test$df <- to_filter %>%
  subset(Mode == 'LD' & Sample == 'C3') %>%
  mutate(IS_STRONG = as.integer(Gene %in% c(vetted$LD$C3$strong$TFs, vetted$LD$C3$strong$RBPs))) %>%
  mutate(IS_PLAUSIBLE = as.integer(Gene %in% c(vetted$LD$C3$plausible$TFs, vetted$LD$C3$plausible$RBPs))) %>%
  {dplyr::select(., -c(Mode, Sample, meta2d_phase, JTK_adjphase))} %>%
  {`colnames<-`(., colnames(.) %>%
                  gsub(pattern='_', replacement=' ') %>%
                  gsub(pattern='pvalue', replacement='p-value') %>%
                  gsub(pattern='scores norm', replacement='p-value'))}

param_test$strong <- param_test$df %>%
  subset(!xor(`IS STRONG`, `IS PLAUSIBLE`)) %>%
  pivot_longer(cols = -c('Gene', 'IS STRONG', 'IS PLAUSIBLE'), names_to = 'Filter:', values_to = 'M')

param_test$plausible <- param_test$df %>%
  pivot_longer(cols = -c('Gene', 'IS STRONG', 'IS PLAUSIBLE'), names_to = 'Filter:', values_to = 'M')
}

#filters <- c('avg', 'Cell proportion', 'KW p-value', 'JTK p-value', 'meta2d p-value')
#filters_notinvert <- c('avg', 'Cell proportion')
#filters_invert <- setdiff(filters, filters_notinvert)


filters <- names(param_test$df) %>% {.[. %!in% c('Gene', 'IS STRONG', 'IS PLAUSIBLE', 'logCPM', 'min', 'avg', 'rel ampl', 
                                                 'PValue', 'KW p-value')]}
filters_invert <- c('FDR EdgeR','KW FDR', "JTK p-value", "meta2d p-value")
filters_notinvert <- setdiff(filters, filters_invert)

# NO LABELS
ggplot(param_test$strong, aes(d = `IS STRONG`, m = M, color = `Filter:`)) + 
  geom_roc(data = subset(param_test$strong, `Filter:` %in% filters_invert), 
           increasing = FALSE, n.cuts = 0) +
  geom_roc(data = subset(param_test$strong, `Filter:` %in% filters_notinvert), 
           increasing = TRUE, n.cuts = 0) + 
  style_roc()

# WITH LABELS
ggplot(param_test$strong, aes(d = `IS STRONG`, m = M, color = `Filter:`)) + 
  geom_roc(data = subset(param_test$strong, `Filter:` %in% filters_invert), 
           increasing = FALSE, labelround=5, labelsize = 4, cutoffs.at = c(1e-2, 0.1, 0.05, 0.5),
           cutoff.labels = c(1e-2, 0.1, 0.05, 0.5)) +
  geom_roc(data = subset(param_test$strong, `Filter:` %in% filters_notinvert), 
           increasing = TRUE, labelround=5, labelsize = 4, cutoffs.at = c(0.05, 0.5),
           cutoff.labels = c(0.05, 0.5)) + 
  style_roc()

filters_invert <- c('FDR EdgeR','PValue', 'KW FDR', 'KW p-value', "meta2d p-value")

ggplot(param_test$strong, aes(d = `IS STRONG`, m = M, color = `Filter:`)) + 
  geom_roc(data = subset(param_test$strong, `Filter:` %in% filters_invert), 
           increasing = FALSE, labelround=5, labelsize = 4, cutoffs.at = c(1e-2, 0.1, 0.2, 0.05, 0.5),
           cutoff.labels = c(1e-2, 0.1, 0.2, 0.05, 0.5))

####    ####    ####    ####    ####    ####    ####    ####    ####    ####   


# 9. TESTING FILTERING APPROACHES; PRINTING
################################################################################

to_filter_approaches <- to_filter
  

# INITIAL NOTES
###############

# optimal at meta2d 0.2 (preserves C3 Pdp1)
to_filter %>%
  filterer("second_pct>0.1") %>%
  filterer("ratio>1.5") %>% 
  {lapply(paste0('meta2d_pvalue<', c(1e-4, 1e-3, seq(0.01, 0.1, 0.01), seq(0.1, 0.5, 0.1), 1)),  # meta2d_pvalue, meta2d_pvalue_AGG_norm, meta2d_pvalue_AGG_CLRnorm, meta2d_pvalue_MCNorm
          function(expr_vec, df=.) {
            df = filterer(df, expr_vec) %>%
              mutate(Sample_Mode = paste0(Sample, '_', Mode), .keep='unused')
            #df = mutate(df, KW_FDR = p.adjust(KW_scores_norm, method='BH'),
            #            FDR_EdgeR = p.adjust(PValue, method='BH'), 
            #            FDR_meta2d = p.adjust(meta2d_pvalue, method='BH')) 
            
            per_in = df[df$Gene =='per', 'Sample_Mode'] %>% paste(collapse = ', ')
            Pdp1_in = df[df$Gene =='Pdp1', 'Sample_Mode'] %>% paste(collapse = ', ')
            cwo_in = df[df$Gene =='cwo', 'Sample_Mode'] %>% paste(collapse = ', ')
            Hr38_in = df[df$Gene =='Hr38', 'Sample_Mode'] %>% paste(collapse = ', ')
            
            #df = filterer(df, "FDR_EdgeR<0.05")
            #df = filterer(df, "KW_FDR<0.05")
            
            
            #df = mutate(df, KW_FDR = p.adjust(KW_scores_norm, method='BH'),
            #            FDR_EdgeR = p.adjust(PValue, method='BH'), 
            #            FDR_meta2d = p.adjust(meta2d_pvalue, method='BH'))
            
            #Pdp1_FDR_meta2d = subset(df, Gene == 'Pdp1' & Sample == 'C3' & Mode == 'LD')$FDR_meta2d
            #cwo_FDR_meta2d = subset(df, Gene == 'cwo' & Sample == 'C3' & Mode == 'LD')$FDR_meta2d
            
            nGenes = length(unique(df$Gene))
          
            nGenesCond = split(df, list(df$Sample_Mode)) %>%
              lapply(nrow) %>%
              paste(collapse = ',')
            
            df = subset(df, Sample_Mode='C3_LD')
            
            df = subset(df, Gene %in% found_TFRBPs)

            g = unique(df$Gene)

            cbind(Filter = expr_vec, 
                  `retained-sens` = length(which(c(vetted$LD$C3$plausible$TFs, vetted$LD$C3$plausible$RBPs) %in% g))/length(c(vetted$LD$C3$plausible$TFs, vetted$LD$C3$plausible$RBPs)),
                  `included-spef` = length(which(g %in% c(vetted$LD$C3$plausible$TFs, vetted$LD$C3$plausible$RBPs)))/length(g),
                  nGenes = nGenes, `C2_DD,LD,C3_DD,LD`=nGenesCond, 
                  per_in=per_in, Pdp1_in=Pdp1_in, cwo_in=cwo_in, Hr38_in=Hr38_in
            )}
  )} %>% 
  do.call(what='rbind') %>%
  as.data.frame %>% 
  mutate(FOM = sqrt(as.numeric(`retained-sens`))*sqrt(as.numeric(`included-spef`)))%>%
  .[,c(1:3, 10, 4:9)] %>% print

# meta2d-0.2; amp-0.5; ratio-1.3; then FDR_EdgeR<0.05
to_filter %>%
  filterer("meta2d_pvalue<0.2") %>%                        # ADJUST THIS LINE
  filterer("amp>0.5") %>%                                  # ADJUST THIS LINE
  filterer('ratio>1.3') %>%
  {lapply(paste0('ratio>', seq(1, 2, 0.1)), 
          function(expr_vec, df=.) {
            df = filterer(df, expr_vec)
            df = mutate(df, KW_FDR = p.adjust(KW_scores_norm, method='BH'),
                        FDR_EdgeR = p.adjust(PValue, method='BH')) %>%
              filterer("FDR_EdgeR<0.05")
            
            nGenes = length(unique(df$Genes))
            
            df = subset(df, Sample == 'C3' & Mode == 'LD' & Gene %in% found_TFRBPs)
            g = unique(df$Gene)
            
            cbind(Filter = expr_vec, 
                  `retained-sens` = length(which(c(vetted$LD$C3$plausible$TFs, vetted$LD$C3$plausible$RBPs) %in% g))/length(c(vetted$LD$C3$plausible$TFs, vetted$LD$C3$plausible$RBPs)),
                  `included-spef` = length(which(g %in% c(vetted$LD$C3$plausible$TFs, vetted$LD$C3$plausible$RBPs)))/length(g),
                  nGenes = nGenes
            )}
  )} %>%
  do.call(what='rbind') %>%
  as.data.frame %>%
  mutate(FOM = sqrt(as.numeric(`retained-sens`))*sqrt(as.numeric(`included-spef`)))
# For meta2d 0.2: PValue 0.02 to 0.5, 0.58 at 0.05
# For meta2d 0.05: PValue 0.02 to 0.5, 0.61 at 0.05

# For meta2d 0.05: amp 0.4, ratio 1.3 - 0.63
# For meta2d 0.2: amp 0.6 or 0.7, ratio 1.3 - 0.62


####    

# RATIO
#######

# RATIO
#######
subset(to_filter, second_pct > 0.05) %>%
  {table(.$Sample, .$Mode)}
subset(to_filter, second_pct > 0.1) %>%
  {table(.$Sample, .$Mode)}
subset(to_filter, second_pct > 0.15) %>%
  {table(.$Sample, .$Mode)}

{run_name = 'Ratio'
setwd(WD)
dir.create(paste0('./', run_name), showWarnings = FALSE)
setwd(paste0('./', run_name))

name = "QCavg>0.05_meta2d_pvalue<0.05_ratio>1.5"
dir.create(paste0('./', name), showWarnings = FALSE)
setwd(paste0('./', name))

filtered <- to_filter %>%
  subset(ratio > 1.5 & Cell_proportion > 0.05)

make_pdf_filtered(filtered, meta2d_threshold = 0.05, name_distinguisher = name, 
                  CELL = TRUE, CELL_NO_ERROR = TRUE, COND = TRUE, COND_NO_ERROR = TRUE, AGG_NORM = TRUE, AGG_CLRNORM = TRUE,
                  BIO_RIB_CELL = TRUE, BIO_RIB_COND = TRUE)

filtered <- to_filter %>%
  subset(meta2d_pvalue<0.05 & ratio > 1.5 & Cell_proportion > 0.05)


png(filename = paste0('PHASE', name, '.png'))
plot_grid(nrow=1, hist_phase(filtered), hmap_phase(filtered))
dev.off()

name = "QC2most>0.1_meta2d_pvalue<0.05_ratio>1.5"
dir.create(paste0('./', name), showWarnings = FALSE)
setwd(paste0('./', name))


filtered <- to_filter %>%
  subset(ratio > 1.5 & second_pct > 0.1)

make_pdf_filtered(filtered, meta2d_threshold = 0.05, name_distinguisher = name, 
                  CELL = TRUE, CELL_NO_ERROR = TRUE, COND = TRUE, COND_NO_ERROR = TRUE, AGG_NORM = TRUE, AGG_CLRNORM = TRUE,
                  BIO_RIB_CELL = TRUE, BIO_RIB_COND = TRUE)

filtered <- to_filter %>%
  subset(meta2d_pvalue<0.05 & ratio > 1.5 & second_pct > 0.1)

png(filename = paste0('PHASE', name, '.png'))
plot_grid(nrow=1, hist_phase(filtered), hmap_phase(filtered))
dev.off()
setwd(WD)
rm(distinguisher, filtered, name)
}

# FILTERING - RECONSTITUTING




# strict + relaxing = relaxed
# 

#######

# 2d-p, RATIO
#############

{to_filter_approaches %>%
  filterer("second_pct>0.15") %>%
  filterer("meta2d_pvalue<0.05") %>% 
  {lapply(paste0('ratio>', seq(1, 2, 0.1)), 
          function(expr_vec, df=.) {
            df = filterer(df, expr_vec) %>%
              mutate(Sample_Mode = paste0(Sample, '_', Mode), .keep='unused')
            
            per_in = df[df$Gene =='per', 'Sample_Mode'] %>% paste(collapse = ', ')
            Pdp1_in = df[df$Gene =='Pdp1', 'Sample_Mode'] %>% paste(collapse = ', ')
            cwo_in = df[df$Gene =='cwo', 'Sample_Mode'] %>% paste(collapse = ', ')
            Hr38_in = df[df$Gene =='Hr38', 'Sample_Mode'] %>% paste(collapse = ', ')
            
            nGenes = length(unique(df$Genes))
            nGenesCond = split(df, df$Sample_Mode) %>%
              lapply(nrow) %>%
              paste(collapse = ',')
            
            
            df = subset(df, Sample_Mode == 'C3_LD')
            
            df = subset(df, Gene %in% found_TFRBPs)
            g = unique(df$Gene)
            
            cbind(Filter = expr_vec, 
                  `retained-sens` = length(which(c(vetted$LD$C3$plausible$TFs, vetted$LD$C3$plausible$RBPs) %in% g))/length(c(vetted$LD$C3$plausible$TFs, vetted$LD$C3$plausible$RBPs)),
                  `included-spef` = length(which(g %in% c(vetted$LD$C3$plausible$TFs, vetted$LD$C3$plausible$RBPs)))/length(g),
                  nGenes = nGenes, `C2_DD,LD,C3_DD,LD`=nGenesCond, 
                  per_in=per_in, Pdp1_in=Pdp1_in, cwo_in=cwo_in, Hr38_in=Hr38_in
            )}
  )} %>%
  do.call(what='rbind') %>%
  as.data.frame%>%
  mutate(FOM = sqrt(as.numeric(`retained-sens`))*sqrt(as.numeric(`included-spef`)))%>%
  .[,c(1:3, 10, 4:9)] %>% print
  
run_name = 'Ratio'
dir.create(paste0('./', run_name), showWarnings = FALSE)
setwd(paste0('./', run_name))

filtered <- to_filter_approaches %>%
  subset(meta2d_pvalue<0.1 & ratio > 1.5 & Cell_proportion > 0.05)

make_pdf_filtered(filtered, list('meta2d_pvalue<0.1', 'ratio>1.5', 'Cell_proportion>0.05'), 
                  CELL = FALSE, CELL_NO_ERROR = TRUE, COND = TRUE, COND_NO_ERROR = TRUE, 
                  BIO_RIB_CELL = TRUE, BIO_RIB_COND = TRUE)


png(filename = 'PHASE_meta2d_pvalue<0.1_ratio>1.5.png_Cell_proportion>0.05')
plot_grid(nrow=1, hist_phase(filtered), hmap_phase(filtered))
dev.off()



filtered <- to_filter_approaches %>%
  subset(meta2d_pvalue<0.05 & ratio > 1.5 & Cell_proportion > 0.05)

make_pdf_filtered(filtered, list('meta2d_pvalue<0.05', 'ratio>1.5', 'Cell_proportion>0.05'), 
                  CELL = FALSE, CELL_NO_ERROR = TRUE, COND = TRUE, COND_NO_ERROR = TRUE,
                  BIO_RIB_CELL = TRUE, BIO_RIB_COND = TRUE)


png(filename = 'PHASE_meta2d_pvalue<0.1_ratio>1.5_Cell_proportion>0.05.png')
plot_grid(nrow=1, hist_phase(filtered), hmap_phase(filtered))
dev.off()


filtered <- to_filter_approaches %>%
  subset(meta2d_pvalue<0.1 & ratio > 1.5 & Cell_proportion > 0.1)

make_pdf_filtered(filtered, list('meta2d_pvalue<0.1', 'ratio>1.5', 'Cell_proportion>0.1'), 
                  CELL = FALSE, CELL_NO_ERROR = TRUE, COND = TRUE, COND_NO_ERROR = TRUE,
                  BIO_RIB_CELL = TRUE, BIO_RIB_COND = TRUE)


png(filename = 'PHASE_meta2d_pvalue<0.1_ratio>1.5.png_Cell_proportion>0.1')
plot_grid(nrow=1, hist_phase(filtered), hmap_phase(filtered))
dev.off()


filtered <- to_filter_approaches %>%
  subset(meta2d_pvalue<0.05 & ratio > 1.5 & Cell_proportion > 0.1)

make_pdf_filtered(filtered, list('meta2d_pvalue<0.05', 'ratio>1.5', 'Cell_proportion>0.1'), 
                  CELL = FALSE, CELL_NO_ERROR = TRUE, COND = TRUE, COND_NO_ERROR = TRUE,
                  BIO_RIB_CELL = TRUE, BIO_RIB_COND = TRUE)


png(filename = 'PHASE_meta2d_pvalue<0.1_ratio>1.5_Cell_proportion>0.1.png')
plot_grid(nrow=1, hist_phase(filtered), hmap_phase(filtered))
dev.off()



setwd('/Users/teddyrashkover/Documents/GitHub/CCG_distribution')
rm(run_name, filtered)
}
#######

# 2d-p, RATIO+AMPL
##################

{run_name = 'Ratio_Ampl'
dir.create(paste0('./', run_name), showWarnings = FALSE)
setwd(paste0('./', run_name))

to_filter_approaches %>%
  filterer("second_pct>0.10") %>%
  filterer("meta2d_pvalue<0.05") %>% 
  filterer("ratio>1.3") %>% 
  {lapply(paste0('ampl>', seq(0, 1.5, 0.1)), 
          function(expr_vec, df=.) {
            df = filterer(df, expr_vec) %>%
              mutate(Sample_Mode = paste0(Sample, '_', Mode), .keep='unused')
            
            per_in = df[df$Gene =='per', 'Sample_Mode'] %>% paste(collapse = ', ')
            Pdp1_in = df[df$Gene =='Pdp1', 'Sample_Mode'] %>% paste(collapse = ', ')
            cwo_in = df[df$Gene =='cwo', 'Sample_Mode'] %>% paste(collapse = ', ')
            Hr38_in = df[df$Gene =='Hr38', 'Sample_Mode'] %>% paste(collapse = ', ')
            
            nGenes = length(unique(df$Genes))
            nGenesCond = split(df, df$Sample_Mode) %>%
              lapply(nrow) %>%
              paste(collapse = ',')
            

            df = subset(df, Sample_Mode == 'C3_LD')
            
            df = subset(df, Gene %in% found_TFRBPs)
            g = unique(df$Gene)
            
            cbind(Filter = expr_vec, 
                  `retained-sens` = length(which(c(vetted$LD$C3$plausible$TFs, vetted$LD$C3$plausible$RBPs) %in% g))/length(c(vetted$LD$C3$plausible$TFs, vetted$LD$C3$plausible$RBPs)),
                  `included-spef` = length(which(g %in% c(vetted$LD$C3$plausible$TFs, vetted$LD$C3$plausible$RBPs)))/length(g),
                  nGenes = nGenes, `C2_DD,LD,C3_DD,LD`=nGenesCond, 
                  per_in=per_in, Pdp1_in=Pdp1_in, cwo_in=cwo_in, Hr38_in=Hr38_in
            )}
  )} %>%
  do.call(what='rbind') %>%
  as.data.frame %>%
  mutate(FOM = sqrt(as.numeric(`retained-sens`))*sqrt(as.numeric(`included-spef`))) %>%
  .[,c(1:3, 10, 4:9)] %>% print

dir.create('./lenient', showWarnings = FALSE)
setwd('./lenient/')

filtered <- to_filter_approaches %>%
  subset(meta2d_pvalue<0.05 & ratio > 1.5 & ampl > 0.3)

make_pdf_filtered(filtered, list('meta2d_pvalue<0.05', 'ratio>1.5', 'ampl>0.3'), 
                  CELL = TRUE, CELL_NO_ERROR = TRUE, COND = TRUE, COND_NO_ERROR = TRUE,
                  BIO_RIB_CELL = TRUE, BIO_RIB_COND = TRUE)

png(filename = 'PHASE_meta2d_pvalue<0.05_ratio>1.5_ampl>0.3.png')
plot_grid(nrow=1, hist_phase(filtered), hmap_phase(filtered))
dev.off()

setwd(paste0('/Users/teddyrashkover/Documents/GitHub/CCG_distribution/', run_name))
dir.create('./QCstringent', showWarnings = FALSE)
setwd('./QCstringent')



filtered <- to_filter_approaches %>%
  subset(meta2d_pvalue<0.05 & ratio > 1.5 & ampl > 0.8)

make_pdf_filtered(filtered, list('meta2d_pvalue<0.05', 'ratio>1.5', 'ampl>0.8'), 
                  CELL = TRUE, CELL_NO_ERROR = TRUE, COND = TRUE, COND_NO_ERROR = TRUE,
                  BIO_RIB_CELL = TRUE, BIO_RIB_COND = TRUE)

png(filename = 'PHASE_meta2d_pvalue<0.05_ratio>1.5_ampl>0.8.png')
plot_grid(nrow=1, hist_phase(filtered), hmap_phase(filtered))
dev.off()



setwd(paste0('/Users/teddyrashkover/Documents/GitHub/CCG_distribution/', run_name))
dir.create('./QCstringent_lowRatio', showWarnings = FALSE)
setwd('./QCstringent_lowRatio')

filtered <- to_filter_approaches %>%
  subset(meta2d_pvalue<0.05 & ratio > 1.3 & ampl > 0.8)

make_pdf_filtered(filtered, list('meta2d_pvalue<0.05', 'ratio>1.3', 'ampl>0.8'), 
                  CELL = TRUE, CELL_NO_ERROR = TRUE, COND = TRUE, COND_NO_ERROR = TRUE,
                  BIO_RIB_CELL = TRUE, BIO_RIB_COND = TRUE)

png(filename = 'PHASE_meta2d_pvalue<0.05_ratio>1.3_ampl>0.8.png')
plot_grid(nrow=1, hist_phase(filtered), hmap_phase(filtered))
dev.off()

setwd('/Users/teddyrashkover/Documents/GitHub/CCG_distribution')
rm(run_name, filtered)
}

#######

# 2d-p, KW-p
############

{run_name = 'KW-p'
dir.create(paste0('./', run_name), showWarnings = FALSE)
setwd(paste0('./', run_name))

to_filter_approaches %>%
  filterer("second_pct>0.15") %>%
  filterer("meta2d_pvalue<0.05") %>% 
  {lapply(paste0('KW_scores_norm<', c(1e-4, 1e-3, seq(0.01, 0.1, 0.01), seq(0.1, 0.5, 0.1), 1)), 
          function(expr_vec, df=.) {
            
            df = filterer(df, expr_vec) %>%
              mutate(Sample_Mode = paste0(Sample, '_', Mode), .keep='unused')
            
            per_in = df[df$Gene =='per', 'Sample_Mode'] %>% paste(collapse = ', ')
            sr_in = df[df$Gene =='sr', 'Sample_Mode'] %>% paste(collapse = ', ')
            CrebA_in = df[df$Gene =='CrebA', 'Sample_Mode'] %>% paste(collapse = ', ')
            fru_in = df[df$Gene =='fru', 'Sample_Mode'] %>% paste(collapse = ', ')
            Pka_in = df[df$Gene =='Pka-C1', 'Sample_Mode'] %>% paste(collapse = ', ')
            Hr38_in = df[df$Gene =='Hr38', 'Sample_Mode'] %>% paste(collapse = ', ')
            
            nGenes = length(unique(df$Genes))
            nGenesCond = split(df, df$Sample_Mode) %>%
              lapply(nrow) %>%
              paste(collapse = ',')
            
            
            df = subset(df, Sample_Mode == 'C3_LD')
            
            df = subset(df, Gene %in% found_TFRBPs)
            g = unique(df$Gene)
            
            cbind(Filter = expr_vec, 
                  `retained-sens` = length(which(c(vetted$LD$C3$plausible$TFs, vetted$LD$C3$plausible$RBPs) %in% g))/length(c(vetted$LD$C3$plausible$TFs, vetted$LD$C3$plausible$RBPs)),
                  `included-spef` = length(which(g %in% c(vetted$LD$C3$plausible$TFs, vetted$LD$C3$plausible$RBPs)))/length(g),
                  nGenes = nGenes, `C2_DD,LD,C3_DD,LD`=nGenesCond, 
                  per_in=per_in, sr_in=sr_in, fru_in=fru_in, Pka_in=Pka_in, Hr38_in=Hr38_in, CrebA_in=CrebA_in
            )}
  )} %>%
  do.call(what='rbind') %>%
  as.data.frame%>%
  mutate(FOM = sqrt(as.numeric(`retained-sens`))*sqrt(as.numeric(`included-spef`)))%>%
  .[,c(1:3, 12, 4:11)] %>% print

filtered <- to_filter_approaches %>%
  subset(meta2d_pvalue<0.05 & KW_scores_norm < 0.05)

make_pdf_filtered(filtered, list('meta2d_pvalue<0.05', 'KW_scores_norm < 0.05'), 
                  CELL = TRUE, CELL_NO_ERROR = TRUE, COND = TRUE, COND_NO_ERROR = TRUE,
                  BIO_RIB_CELL = TRUE, BIO_RIB_COND = TRUE)


png(filename = 'PHASE_meta2d_pvalue<0.05_KW_scores_norm<0.05.png')
plot_grid(nrow=1, hist_phase(filtered), hmap_phase(filtered))
dev.off()

setwd('/Users/teddyrashkover/Documents/GitHub/CCG_distribution')
rm(run_name, filtered)
}

#######

# 2d-q, RATIO - unpicked
#############

{run_name = '2Dq_Ratio'
dir.create(paste0('./', run_name), showWarnings = FALSE)
setwd(paste0('./', run_name))

to_filter_approaches %>%
  #filterer("second_pct>0.10") %>%
  filterer("ratio>1.5") %>%
  #{lapply(paste0('ratio>', seq(1, 2, 0.1)), 
  {lapply(paste0('second_pct>', seq(0,0.25,0.025)), 
          function(expr_vec, df=.) {
            df = filterer(df, expr_vec) %>%
              mutate(Sample_Mode = paste0(Sample, '_', Mode), .keep='unused')
            
            #df = mutate(df, FDR_meta2d = p.adjust(meta2d_pvalue, method='BH'))
            #df = mutate(df, FDR_meta2d_AGG = p.adjust(meta2d_pvalue_AGG_norm, method='BH'))
            df = mutate(df, FDR_meta2d_1d = p.adjust(meta2d_pvalue_norm1d, method='BH'))
            #df = filterer(df, 'FDR_meta2d<0.1')
            #df = filterer(df, 'FDR_meta2d_AGG<0.1')
            df = filterer(df, 'FDR_meta2d_1d<0.1')
            
            
            per_in = df[df$Gene =='per', 'Sample_Mode'] %>% paste(collapse = ', ')
            sr_in = df[df$Gene =='sr', 'Sample_Mode'] %>% paste(collapse = ', ')
            CrebA_in = df[df$Gene =='CrebA', 'Sample_Mode'] %>% paste(collapse = ', ')
            fru_in = df[df$Gene =='fru', 'Sample_Mode'] %>% paste(collapse = ', ')
            Pka_in = df[df$Gene =='Pka-C1', 'Sample_Mode'] %>% paste(collapse = ', ')
            Hr38_in = df[df$Gene =='Hr38', 'Sample_Mode'] %>% paste(collapse = ', ')
            
            nGenes = length(unique(df$Genes))
            nGenesCond = split(df, df$Sample_Mode) %>%
              lapply(nrow) %>%
              paste(collapse = ',')
            
            
            df = subset(df, Sample_Mode == 'C3_LD')
            
            df = subset(df, Gene %in% found_TFRBPs)
            g = unique(df$Gene)
            
            cbind(Filter = expr_vec, 
                  `retained-sens` = length(which(c(vetted$LD$C3$plausible$TFs, vetted$LD$C3$plausible$RBPs) %in% g))/length(c(vetted$LD$C3$plausible$TFs, vetted$LD$C3$plausible$RBPs)),
                  `included-spef` = length(which(g %in% c(vetted$LD$C3$plausible$TFs, vetted$LD$C3$plausible$RBPs)))/length(g),
                  nGenes = nGenes, `C2_DD,LD,C3_DD,LD`=nGenesCond, 
                  per_in=per_in, sr_in=sr_in, fru_in=fru_in, Pka_in=Pka_in, Hr38_in=Hr38_in, CrebA_in=CrebA_in
            )}
  )} %>%
  do.call(what='rbind') %>%
  as.data.frame%>%
  mutate(FOM = sqrt(as.numeric(`retained-sens`))*sqrt(as.numeric(`included-spef`)))%>%
  .[,c(1:3, 12, 4:11)] %>% print

filtered <- to_filter_approaches %>%
  subset(meta2d_pvalue<0.05 & ratio > 1.5)

make_pdf_filtered(filtered, list('meta2d_pvalue<0.05', 'ratio>1.5'), 
                  CELL = FALSE, CELL_NO_ERROR = TRUE, COND = TRUE, COND_NO_ERROR = TRUE,
                  BIO_RIB_CELL = TRUE, BIO_RIB_COND = TRUE)


png(filename = 'PHASE_meta2d_pvalue<0.05_ratio>1.5.png')
plot_grid(nrow=1, hist_phase(filtered), hmap_phase(filtered))
dev.off()

setwd('/Users/teddyrashkover/Documents/GitHub/CCG_distribution')
rm(run_name, filtered)
}
#######

# 2d-q, RATIO+AMPL - unpicked
##################

{run_name = '2Dq_Ratio_Ampl'
dir.create(paste0('./', run_name), showWarnings = FALSE)
setwd(paste0('./', run_name))

to_filter_approaches %>%
  filterer("second_pct>0.15") %>%
  filterer("ratio>1.5") %>% 
  {lapply(paste0('ampl>', seq(0, 1.5, 0.1)), 
          function(expr_vec, df=.) {
            df = filterer(df, expr_vec) %>%
              mutate(Sample_Mode = paste0(Sample, '_', Mode), .keep='unused')
            
            df = mutate(df, FDR_meta2d = p.adjust(meta2d_pvalue, method='BH'))
            df = filterer(df, 'FDR_meta2d<0.1')
            
            per_in = df[df$Gene =='per', 'Sample_Mode'] %>% paste(collapse = ', ')
            Pdp1_in = df[df$Gene =='Pdp1', 'Sample_Mode'] %>% paste(collapse = ', ')
            cwo_in = df[df$Gene =='cwo', 'Sample_Mode'] %>% paste(collapse = ', ')
            Hr38_in = df[df$Gene =='Hr38', 'Sample_Mode'] %>% paste(collapse = ', ')
            
            nGenes = length(unique(df$Genes))
            nGenesCond = split(df, df$Sample_Mode) %>%
              lapply(nrow) %>%
              paste(collapse = ',')
            
            
            df = subset(df, Sample_Mode == 'C3_LD')
            
            df = subset(df, Gene %in% found_TFRBPs)
            g = unique(df$Gene)
            
            cbind(Filter = expr_vec, 
                  `retained-sens` = length(which(c(vetted$LD$C3$plausible$TFs, vetted$LD$C3$plausible$RBPs) %in% g))/length(c(vetted$LD$C3$plausible$TFs, vetted$LD$C3$plausible$RBPs)),
                  `included-spef` = length(which(g %in% c(vetted$LD$C3$plausible$TFs, vetted$LD$C3$plausible$RBPs)))/length(g),
                  nGenes = nGenes, `C2_DD,LD,C3_DD,LD`=nGenesCond, 
                  per_in=per_in, Pdp1_in=Pdp1_in, cwo_in=cwo_in, Hr38_in=Hr38_in
            )}
  )} %>%
  do.call(what='rbind') %>%
  as.data.frame %>%
  mutate(FOM = sqrt(as.numeric(`retained-sens`))*sqrt(as.numeric(`included-spef`))) %>%
  .[,c(1:3, 10, 4:9)] %>% print

filtered <- to_filter_approaches %>%
  subset(meta2d_pvalue<0.05 & ratio > 1.5 & ampl > 1.5)

make_pdf_filtered(filtered, list('meta2d_pvalue<0.05', 'ratio>1.5', 'ampl>1.5'), 
                  CELL = TRUE, CELL_NO_ERROR = TRUE, COND = TRUE, COND_NO_ERROR = TRUE,
                  BIO_RIB_CELL = TRUE, BIO_RIB_COND = TRUE)

png(filename = 'PHASE_meta2d_pvalue<0.05_ratio>1.5_ampl>1.5.png')
plot_grid(nrow=1, hist_phase(filtered), hmap_phase(filtered))
dev.off()

setwd('/Users/teddyrashkover/Documents/GitHub/CCG_distribution')
rm(run_name, filtered)
}

#######


# 2d-q, KW-p - unpicked
############

{run_name = '2Dq_KW-p'
dir.create(paste0('./', run_name), showWarnings = FALSE)
setwd(paste0('./', run_name))

to_filter_approaches %>%
  filterer("second_pct>0.10") %>%
  filterer('ratio>1.5') %>%
  #filterer('ampl>0.4') %>%
  {lapply(paste0('KW_scores_norm<', c(1e-4, 1e-3, seq(0.01, 0.1, 0.01), seq(0.1, 0.5, 0.1), 1)), 
          function(expr_vec, df=.) {
            df = filterer(df, expr_vec) %>%
              mutate(Sample_Mode = paste0(Sample, '_', Mode), .keep='unused')
            
            df = mutate(df, FDR_meta2d = p.adjust(meta2d_pvalue, method='BH'))
            df = filterer(df, 'FDR_meta2d<0.1')
            
            per_in = df[df$Gene =='per', 'Sample_Mode'] %>% paste(collapse = ', ')
            sr_in = df[df$Gene =='sr', 'Sample_Mode'] %>% paste(collapse = ', ')
            CrebA_in = df[df$Gene =='CrebA', 'Sample_Mode'] %>% paste(collapse = ', ')
            fru_in = df[df$Gene =='fru', 'Sample_Mode'] %>% paste(collapse = ', ')
            Pka_in = df[df$Gene =='Pka-C1', 'Sample_Mode'] %>% paste(collapse = ', ')
            Hr38_in = df[df$Gene =='Hr38', 'Sample_Mode'] %>% paste(collapse = ', ')
            
            nGenes = length(unique(df$Genes))
            nGenesCond = split(df, df$Sample_Mode) %>%
              lapply(nrow) %>%
              paste(collapse = ',')
            
            
            df = subset(df, Sample_Mode == 'C3_LD')
            
            df = subset(df, Gene %in% found_TFRBPs)
            g = unique(df$Gene)
            
            cbind(Filter = expr_vec, 
                  `retained-sens` = length(which(c(vetted$LD$C3$plausible$TFs, vetted$LD$C3$plausible$RBPs) %in% g))/length(c(vetted$LD$C3$plausible$TFs, vetted$LD$C3$plausible$RBPs)),
                  `included-spef` = length(which(g %in% c(vetted$LD$C3$plausible$TFs, vetted$LD$C3$plausible$RBPs)))/length(g),
                  nGenes = nGenes, `C2_DD,LD,C3_DD,LD`=nGenesCond, 
                  per_in=per_in, sr_in=sr_in, fru_in=fru_in, Pka_in=Pka_in, Hr38_in=Hr38_in, CrebA_in=CrebA_in
            )}
  )} %>%
  do.call(what='rbind') %>%
  as.data.frame%>%
  mutate(FOM = sqrt(as.numeric(`retained-sens`))*sqrt(as.numeric(`included-spef`)))%>%
  .[,c(1:3, 12, 4:11)] %>% print

{pdf('onepage_noC3LD_KW_q_meta2d_q.pdf', height=60, width=10)
my_goodness = function(L) {
  plot_grid(ncol = 2, 
    avg_plotter_concat_mode(L[['C3.LD']]$Gene, data=data_C2C3, N_days=2, rib_alpha=0.4),
    plot_grid(ncol=1, rel_heights=c(5,3),
              avg_plotter_concat_mode(L[['C2.LD']]$Gene, data=data_C2C3, N_days=2, rib_alpha=0.4),
              plot_grid(ncol=1, rel_heights=c(5,3), 
                        avg_plotter_concat_mode(L[['C3.DD']]$Gene, data=data_C2C3, N_days=2, rib_alpha=0.4),
                        avg_plotter_concat_mode(L[['C2.DD']]$Gene, data=data_C2C3, N_days=2, rib_alpha=0.4))))
}
to_filter_approaches %>%
  filterer("second_pct>0.10") %>%
  #filterer('ratio>1.5') %>%
  mutate(FDR_KW = p.adjust(KW_scores_norm, method='BH')) %>%
  filterer('FDR_KW<0.1') %>%
  mutate(FDR_meta2d = p.adjust(meta2d_pvalue, method='BH')) %>%
  filterer('FDR_meta2d<0.1') %>%
  subset(select = c(Gene, Sample, Mode)) %>%
  {split(., list(.$Sample, .$Mode))} %>%
  my_goodness %>% print
rm(my_goodness)
dev.off()
}



filtered_KWq_2Dq_ratio <- to_filter_approaches %>%
  filterer("second_pct>0.10") %>%
  filterer('ratio>1.5') %>%
  mutate(FDR_KW = p.adjust(KW_scores_norm, method='BH')) %>%
  filterer('FDR_KW<0.1') %>%
  mutate(FDR_meta2d = p.adjust(meta2d_pvalue, method='BH')) %>%
  filterer('FDR_meta2d<0.1') %>%
  subset(select = c(Gene, Sample, Mode))

# # # # #
to_filter_approaches %>%
  subset(meta2d_pvalue<0.05 & KW_scores_norm < 0.05) %>% 
  make_pdf_filtered(list('meta2d_pvalue<0.05', 'KW_scores_norm < 0.05'), 
                  CELL = TRUE, CELL_NO_ERROR = TRUE, COND = TRUE, COND_NO_ERROR = TRUE,
                  BIO_RIB_CELL = TRUE, BIO_RIB_COND = TRUE)


png(filename = 'PHASE_meta2d_pvalue<0.05_KW_scores_norm<0.05.png')
plot_grid(nrow=1, hist_phase(filtered), hmap_phase(filtered))
dev.off()

setwd('/Users/teddyrashkover/Documents/GitHub/CCG_distribution')
rm(run_name, filtered)
}



ggplot(to_filter, aes(y=KW_scores_norm, x=meta2d_pvalue)) + 
  geom_point() + 
  geom_smooth(method="lm", col="black") + 
  annotate("text", x=20, y=4.5, label=paste0("r = ", round(cor(mtcars$wt, mtcars$mpg), 2)), hjust=0) +
  annotate("text", x=20, y=4.25, label=paste0("p = ", round(cor.test(mtcars$wt, mtcars$mpg)$p.value, 3)), hjust=0) +
  theme_classic()


######

# 2d-q, KW=p+RATIO+AMPL - unpicked
#######################

{run_name = '2Dq_KW-p_RATIO_AMPL'
dir.create(paste0('./', run_name), showWarnings = FALSE)
setwd(paste0('./', run_name))

to_filter_approaches %>%
  filterer("Cell_proportion>0.15") %>%
  {lapply(paste0('KW_scores_norm<', c(1e-4, 1e-3, seq(0.01, 0.1, 0.01), seq(0.1, 0.5, 0.1), 1)), 
          function(expr_vec, df=.) {
            df = filterer(df, expr_vec) %>%
              mutate(Sample_Mode = paste0(Sample, '_', Mode), .keep='unused')
            
            df = mutate(df, FDR_meta2d = p.adjust(meta2d_pvalue, method='BH'))
            df = filterer(df, 'FDR_meta2d<0.1')
            
            per_in = df[df$Gene =='per', 'Sample_Mode'] %>% paste(collapse = ', ')
            Pdp1_in = df[df$Gene =='Pdp1', 'Sample_Mode'] %>% paste(collapse = ', ')
            cwo_in = df[df$Gene =='cwo', 'Sample_Mode'] %>% paste(collapse = ', ')
            Hr38_in = df[df$Gene =='Hr38', 'Sample_Mode'] %>% paste(collapse = ', ')
            
            nGenes = length(unique(df$Genes))
            nGenesCond = split(df, df$Sample_Mode) %>%
              lapply(nrow) %>%
              paste(collapse = ',')
            
            
            df = subset(df, Sample_Mode == 'C3_LD')
            
            df = subset(df, Gene %in% found_TFRBPs)
            g = unique(df$Gene)
            
            cbind(Filter = expr_vec, 
                  `retained-sens` = length(which(c(vetted$LD$C3$plausible$TFs, vetted$LD$C3$plausible$RBPs) %in% g))/length(c(vetted$LD$C3$plausible$TFs, vetted$LD$C3$plausible$RBPs)),
                  `included-spef` = length(which(g %in% c(vetted$LD$C3$plausible$TFs, vetted$LD$C3$plausible$RBPs)))/length(g),
                  nGenes = nGenes, `C2_DD,LD,C3_DD,LD`=nGenesCond, 
                  per_in=per_in, Pdp1_in=Pdp1_in, cwo_in=cwo_in, Hr38_in=Hr38_in
            )}
  )} %>%
  do.call(what='rbind') %>%
  as.data.frame %>%
  mutate(FOM = sqrt(as.numeric(`retained-sens`))*sqrt(as.numeric(`included-spef`))) %>%
  .[,c(1:3, 10, 4:9)] %>% print

filtered <- to_filter_approaches %>%
  subset(meta2d_pvalue<0.05 & KW_scores_norm < 0.05)

make_pdf_filtered(filtered, list('meta2d_pvalue<0.05', 'KW_scores_norm < 0.05'), 
                  CELL = TRUE, CELL_NO_ERROR = TRUE, COND = TRUE, COND_NO_ERROR = TRUE,
                  BIO_RIB_CELL = TRUE, BIO_RIB_COND = TRUE)


png(filename = 'PHASE_meta2d_pvalue<0.05_KW_scores_norm<0.05.png')
plot_grid(nrow=1, hist_phase(filtered), hmap_phase(filtered))
dev.off()

setwd('/Users/teddyrashkover/Documents/GitHub/CCG_distribution')
rm(run_name, filtered)
}
####    ####    ####    ####    ####    ####    ####    ####    ####    ####   


# 10. FINAL FILTERING: RECONSTITUTING (filtered_plots)
################################################################################
genes_activity = c('sr', 'CG14186', 'Hr38')
genes_CREB_immediate = c('CrebA', 'CrebB', 'Pka-C1', 'Pka-R1', 'Pka-C3')
genes_CREB_involved = c('Rala', 'dnc', 'kay', 'crc')
#genes_MAPK_involved = c('cic')


avg_plotter_concat_mode(c(genes_activity, genes_CREB_immediate, genes_CREB_involved, genes_MAPK_involved, 'cwo'), rib_alpha = 0.4, N_days = 2)
avg_plotter_concat_mode(c(genes_activity, genes_CREB_immediate, genes_CREB_involved, genes_MAPK_involved, 'cwo'), rib_alpha = 0.4, N_days = 2, data = data_C2C3_AGG)
avg_plotter_concat_mode(c(genes_activity, genes_CREB_immediate, genes_CREB_involved, genes_MAPK_involved, 'cwo'), rib_alpha = 0.4, N_days = 2, data = data_C2C3_normGroup2, COUNTS=TRUE)


reconstituting_DFs_15 <- list(in_both = subset(to_filter, second_pct > 0.15 & ratio > 1.5 & ratio_AGG > 1.5 & meta2d_pvalue < 0.05 & meta2d_pvalue_AGG_norm < 0.05),
                              fail_cell = subset(to_filter, second_pct > 0.15 & ratio_AGG > 1.5 & meta2d_pvalue_AGG_norm < 0.05 & (meta2d_pvalue < 0.05 | ratio < 1.5)),
                              fail_AGG = subset(to_filter, second_pct > 0.15 & ratio > 1.5 & meta2d_pvalue < 0.05 & (meta2d_pvalue_AGG_norm < 0.05 | ratio_AGG < 1.5)),
                              fails_returned = subset(to_filter, second_pct > 0.15 & (meta2d_pvalue > 0.05 | ratio < 1.5) & meta2d_pvalue_AGG_norm < 0.05 & ratio_AGG > 2),
                              fails_not_returned = subset(to_filter, second_pct > 0.15 & ratio > 1.5 & meta2d_pvalue_AGG_norm < 0.05 & meta2d_pvalue > 0.05 & ratio_AGG < 2),
                              reconstituted = subset(to_filter, second_pct > 0.15 & ((ratio > 1.5 & meta2d_pvalue < 0.05 & ratio_AGG > 1.5 & meta2d_pvalue_AGG_norm < 0.05) | 
                                                                                       (meta2d_pvalue_AGG_norm < 0.05 & ratio_AGG > 2))),
                              QC_removed = subset(to_filter, second_pct < 0.15))

reconstituting_DFs_10 <- list(in_both = subset(to_filter, second_pct > 0.10 & ratio > 1.5 & ratio_AGG > 1.5 & meta2d_pvalue < 0.05 & meta2d_pvalue_AGG_norm < 0.05),
                              fail_cell = subset(to_filter, second_pct > 0.10 & ratio_AGG > 1.5 & meta2d_pvalue_AGG_norm < 0.05 & (meta2d_pvalue < 0.05 | ratio < 1.5)),
                              fail_AGG = subset(to_filter, second_pct > 0.10 & ratio > 1.5 & meta2d_pvalue < 0.05 & (meta2d_pvalue_AGG_norm < 0.05 | ratio_AGG < 1.5)),
                              fails_returned = subset(to_filter, second_pct > 0.10 & (meta2d_pvalue > 0.05 | ratio < 1.5) & meta2d_pvalue_AGG_norm < 0.05 & ratio_AGG > 2),
                              fails_not_returned = subset(to_filter, second_pct > 0.10 & ratio > 1.5 & meta2d_pvalue_AGG_norm < 0.05 & meta2d_pvalue > 0.05 & ratio_AGG < 2),
                              reconstituted = subset(to_filter, second_pct > 0.10 & ((ratio > 1.5 & meta2d_pvalue < 0.05 & ratio_AGG > 1.5 & meta2d_pvalue_AGG_norm < 0.05) | 
                                                                                       (meta2d_pvalue_AGG_norm < 0.05 & ratio_AGG > 2))),
                              QC_removed = subset(to_filter, second_pct < 0.10))


# RETENTION OF IMPORTANT GENES
{reconstituting_DFs = reconstituting_DFs_10
  #reconstituting_DFs = reconstituting_DFs_15
  reconstituting_genes <- c(circ_genes, 'cyc', genes_activity, genes_CREB_immediate, genes_CREB_involved)
  #reconstituting_genes = 'tyf'
  reconstituting_genes <- lapply(unique(paste(to_filter$Sample, to_filter$Mode, sep='_')), function(S_M) {
    lapply(reconstituting_DFs, function(df) {
      df = mutate(df, Sample_Mode = paste(df$Sample, df$Mode, sep='_')) %>%
        subset(Sample_Mode == S_M)
      lapply(reconstituting_genes, function(gene, DF=df) {
        assign(gene, gene %in% DF$Gene)
      }) %>% 
        do.call(what='cbind') %>%
        `colnames<-`(reconstituting_genes)
    }) %>% 
      do.call(what='rbind') %>%
      `rownames<-`(names(reconstituting_DFs))
  }) %>%
    `names<-`(unique(paste(to_filter$Sample, to_filter$Mode, sep='_')))
  reconstituting_genes %>% print
  rm(reconstituting_DFs)
}

make_pdf_filtered(reconstituting_DFs_10$reconstituted, name_distinguisher = 'reconstituted', 
                  CELL = FALSE, CELL_NO_ERROR = FALSE, COND = FALSE, COND_NO_ERROR = FALSE, AGG_NORM = TRUE, AGG_CLRNORM = FALSE,
                  BIO_RIB_CELL = FALSE, BIO_RIB_COND = FALSE, ANNOTATE=TRUE)
make_pdf_filtered(reconstituting_DFs_10$reconstituted, name_distinguisher = 'reconstituted_noAnno', 
                  CELL = FALSE, CELL_NO_ERROR = FALSE, COND = FALSE, COND_NO_ERROR = FALSE, AGG_NORM = TRUE, AGG_CLRNORM = FALSE,
                  BIO_RIB_CELL = FALSE, BIO_RIB_COND = FALSE, ANNOTATE=FALSE)


subset(to_filter, second_pct > 0.15 & ratio > 1.5 & ratio_AGG > 1.5 & meta2d_pvalue < 0.05 & meta2d_pvalue_AGG_norm < 0.05) %>%
  make_pdf_filtered(name_distinguisher = 'in_both', 
                  CELL = FALSE, CELL_NO_ERROR = FALSE, COND = FALSE, COND_NO_ERROR = FALSE, AGG_NORM = TRUE, AGG_CLRNORM = FALSE,
                  BIO_RIB_CELL = FALSE, BIO_RIB_COND = FALSE, ANNOTATE=TRUE)

subset(to_filter, second_pct > 0.15 & ratio_AGG > 1.5 & meta2d_pvalue_AGG_norm < 0.05 & (ratio < 1.5 | meta2d_pvalue > 0.05)) %>%
  make_pdf_filtered(name_distinguisher = 'inAGG_notCell', 
                    CELL = FALSE, CELL_NO_ERROR = FALSE, COND = FALSE, COND_NO_ERROR = FALSE, AGG_NORM = TRUE, AGG_CLRNORM = FALSE,
                    BIO_RIB_CELL = FALSE, BIO_RIB_COND = FALSE, ANNOTATE=TRUE)

subset(to_filter, second_pct > 0.15 & ratio > 1.5 & meta2d_pvalue < 0.05 & (ratio_AGG < 1.5 | meta2d_pvalue_AGG_norm > 0.05)) %>%
  make_pdf_filtered(name_distinguisher = 'inCell_notAGG', 
                    CELL = FALSE, CELL_NO_ERROR = FALSE, COND = FALSE, COND_NO_ERROR = FALSE, AGG_NORM = TRUE, AGG_CLRNORM = FALSE,
                    BIO_RIB_CELL = FALSE, BIO_RIB_COND = FALSE, ANNOTATE=TRUE)


plot_grid(nrow=1, hist_phase(reconstituting_DFs_10$reconstituted, CIRCULAR=FALSE, TRIM_TO_24=TRUE, binWidth=2), 
          hist_phase(reconstituting_DFs_10$reconstituted, CIRCULAR=TRUE, TRIM_TO_24=TRUE, binWidth=2), 
          hmap_phase(reconstituting_DFs_10$reconstituted))

{setwd(WD)
dir.create('./filtered_plots', showWarnings=FALSE)
setwd('filtered_plots')
dir.create('./filtered_plots_small_noTabY', showWarnings = FALSE)
setwd('filtered_plots_small_noTabY')
lapply(names(color_list), dir.create, showWarnings=FALSE)
reconstituting_DFs_10$reconstituted %>%
  {split(., list(.$Sample, .$Mode))} %>%
  lapply(function(df) {
    S_M = paste(df[1, 'Sample'], df[1, 'Mode'], sep='_')
    setwd(paste(WD, 'filtered_plots/filtered_plots_small_noTabY', S_M, sep='/'))
    lapply(unique(df$Gene), function(G) {
      avg_plotter_concat(G, data=subset(data_C2C3_AGG, sample_mode == S_M), 
                         BLACKOUT=TRUE, N_days=2, SWAP_FACET=TRUE, NO_GRIDLINES=TRUE, NO_X_TEXT=TRUE, NO_TAB_Y=TRUE) %>%
      ggsave(file=paste0(G, '.svg'), plot=., width=1, height=1)
    })
  })
setwd(WD)
}
{setwd(WD)
  dir.create('./filtered_plots', showWarnings=FALSE)
  setwd('filtered_plots')
  dir.create('./filtered_plots_medium_noTabY', showWarnings = FALSE)
  setwd('filtered_plots_medium_noTabY')
  lapply(names(color_list), dir.create, showWarnings=FALSE)
  reconstituting_DFs_10$reconstituted %>%
    {split(., list(.$Sample, .$Mode))} %>%
    lapply(function(df) {
      S_M = paste(df[1, 'Sample'], df[1, 'Mode'], sep='_')
      setwd(paste(WD, 'filtered_plots/filtered_plots_medium_noTabY', S_M, sep='/'))
      lapply(unique(df$Gene), function(G) {
        avg_plotter_concat(G, data=subset(data_C2C3_AGG, sample_mode == S_M), 
                           BLACKOUT=TRUE, N_days=2, SWAP_FACET=TRUE, NO_GRIDLINES=TRUE, NO_X_TEXT=TRUE, NO_TAB_Y=TRUE) %>%
          ggsave(file=paste0(G, '.svg'), plot=., width=2, height=2)
      })
    })
  setwd(WD)
}
{setwd(WD)
  dir.create('./filtered_plots', showWarnings=FALSE)
  setwd('filtered_plots')
  dir.create('./filtered_plots_large_noTabY', showWarnings = FALSE)
  setwd('filtered_plots_large_noTabY')
  lapply(names(color_list), dir.create, showWarnings=FALSE)
  reconstituting_DFs_10$reconstituted %>%
    {split(., list(.$Sample, .$Mode))} %>%
    lapply(function(df) {
      S_M = paste(df[1, 'Sample'], df[1, 'Mode'], sep='_')
      setwd(paste(WD, 'filtered_plots/filtered_plots_large_noTabY', S_M, sep='/'))
      lapply(unique(df$Gene), function(G) {
        avg_plotter_concat(G, data=subset(data_C2C3_AGG, sample_mode == S_M), 
                           BLACKOUT=TRUE, N_days=2, SWAP_FACET=TRUE, NO_GRIDLINES=TRUE, NO_X_TEXT=TRUE, NO_TAB_Y=TRUE) %>%
          ggsave(file=paste0(G, '.svg'), plot=., width=5, height=5)
      })
    })
  setwd(WD)
}
####    ####    ####    ####    ####    ####    ####    ####    ####    ####   


# KW-p, 2D-q
############
filtered_KWp_2Dq <- to_filter %>%
  subset(second_pct > 0.1 & ratio > 1.5) %>% 
  subset(KW_scores_norm < 0.05) %>%     
  mutate(FDR_meta2d = p.adjust(meta2d_pvalue, method='BH')) %>% 
  subset(FDR_meta2d < 0.1)
# filtered_KWp_2Dq %>% 
#make_FIMO_run(add_to_name='KWp_2dQ', WRITE_PROMOTERS=TRUE,RUN=FALSE)
#make_panel(filename='shuffle_2cells_KWp_2dCellQ.pdf')
#hist_phase(to_select = 'meta2d_phase')
#{split(., list(.$Sample, .$Mode))} %>% 
#lapply(function(df) hist_phase(to_select = 'meta2d_phase')) %>%
#plot_grid(plotlist = .)
#lapply(function(df) df$Gene) %>%
#####


# SHUFFLE
################################################################################
# 2-DAY QUICK COLUMN SHUFFLE
{set.seed(1)
  shuffle_2cq <- list()
  shuffle_2cq$prefilt <- subset(to_filter, second_pct > 0.10 & ratio > 1.5 & ratio_AGG > 1.5, select=c(Gene, Sample, Mode, second_pct, ratio, ratio_AGG, meta2d_pvalue_AGG_norm)) %>%
    mutate(Gene_SM = paste(Gene, Sample, Mode, sep='_'))
  
  shuffle_2cq$AGGnorm$prefiltPreMeta <- AverageExpression(data_C2C3_AGG, group.by = c("group2", "sample", "mode1"), layer="data")$RNA %>% as.data.frame %>%
    {mutate(., Gene = rownames(.))} %>%
    pivot_longer(cols=-'Gene', names_to='TP_sample_mode') %>%
    mutate(Timepoint = TP_sample_mode %>% gsub('_.*', '', .) %>% gsub('.*T', '', .), 
           Gene_SM = paste(Gene, TP_sample_mode %>% gsub('.*_C', 'C', .), sep = '_'), .keep='unused') %>%
    subset(Gene_SM %in% shuffle_2cq$prefilt$Gene_SM) %>%
    pivot_wider(names_from = 'Timepoint', values_from = 'value')
  
  i = 1
  while(1) {
    print(i)
    shuffle_2cq$AGGnorm$prefiltPreMeta %>%
      {`colnames<-`(., colnames(.)[c(1, sample(2:13, replace = FALSE))])} %>%
      write.csv(file = "shuffle_2cq.csv", row.names=FALSE)
    
    meta2d(infile = "shuffle_2cq.csv", filestyle = "csv", timepoints = "line1", minper=20, 
           maxper=28, outdir=getwd(), parallelize = TRUE)
    
    shuffle_2cq$AGGnorm$runPasses[[i]] <- read.csv('meta2d_shuffle_2cq.csv') %>%
      mutate(Gene = gsub('_.*', "", CycID), 
             Mode = gsub('.*_', "", CycID), 
             Sample = gsub('.*_C', 'C', CycID) %>% sub('_.*', '', .),
             .keep='unused') %>%
      subset(select = c(Gene, Mode, Sample, meta2d_pvalue), meta2d_pvalue<0.05) %>% 
      {table(.$Sample, .$Mode)} %>%
      as.data.frame %>%
      `colnames<-`(c('Sample', 'Mode', 'count'))
    i = i+1
  }
  rm(i)
  shuffle_2cq$AGGnorm$counts <- lapply(seq_along(shuffle_2cq$AGGnorm$runPasses), function(i) cbind(shuffle_2cq$AGGnorm$runPasses[[i]], run=i)) %>% 
    c(list(data = cbind(subset(shuffle_2cq$prefilt, meta2d_pvalue_AGG_norm<0.05) %>% 
                          {table(.$Sample, .$Mode)} %>% 
                          as.data.frame %>% 
                          `colnames<-`(c('Sample', 'Mode', 'count')), run='data'))) %>% 
    do.call(what='rbind') %>%
    mutate(Sample = factor(Sample, levels=c('C3', 'C2')), 
           Mode = factor(Mode, levels=c('LD', 'DD')))
  
  shuffle_2cq$AGGnorm$percentiles <- shuffle_2cq$AGGnorm$counts %>%
    {split(., list(.$Sample, .$Mode))} %>%
    lapply(function(df) cbind(df, Rank = rank(df$count, ties.method='average')) %>%
             mutate(percentile = Rank/max(Rank)) %>%
             subset(run == 'data'))
}

shuffle_2cq$AGGnorm$counts %>%
  {ggplot() +
      geom_histogram(data=subset(., run != 'data'), aes(x = count), binwidth=1) +
      geom_vline(data=subset(., run == 'data'), aes(xintercept = count))} +
  facet_grid(Mode ~ Sample)

# 2-DAY PER-ROW SHUFFLE, FULL CYCLING ANALYSIS
{set.seed(1)
  shuffle_2rf <- list()
  shuffle_2rf$prefilt <- subset(to_filter, second_pct > 0.10 & ratio > 1.5 & ratio_AGG > 1.5, select=c(Gene, Sample, Mode, second_pct, ratio, ratio_AGG, meta2d_pvalue, meta2d_pvalue_AGG_norm)) %>%
    mutate(Gene_SM = paste(Gene, Sample, Mode, sep='_'))
  
  shuffle_2rf$cellnorm$prefiltPreMeta <- AverageExpression(data_C2C3, group.by = c("group2", "sample", "mode1"), layer="data")$RNA %>% as.data.frame %>%
    {mutate(., Gene = rownames(.))} %>%
    pivot_longer(cols=-'Gene', names_to='TP_sample_mode') %>%
    mutate(Timepoint = TP_sample_mode %>% gsub('_.*', '', .) %>% gsub('.*T', '', .), 
           Gene_SM = paste(Gene, TP_sample_mode %>% gsub('.*_C', 'C', .), sep = '_'), .keep='unused') %>%
    subset(Gene_SM %in% shuffle_2rf$prefilt$Gene_SM) %>%
    pivot_wider(names_from = 'Timepoint', values_from = 'value')
  
  shuffle_2rf$AGGnorm$prefiltPreMeta <- AverageExpression(data_C2C3_AGG, group.by = c("group2", "sample", "mode1"), layer="data")$RNA %>% as.data.frame %>%
    {mutate(., Gene = rownames(.))} %>%
    pivot_longer(cols=-'Gene', names_to='TP_sample_mode') %>%
    mutate(Timepoint = TP_sample_mode %>% gsub('_.*', '', .) %>% gsub('.*T', '', .), 
           Gene_SM = paste(Gene, TP_sample_mode %>% gsub('.*_C', 'C', .), sep = '_'), .keep='unused') %>%
    subset(Gene_SM %in% shuffle_2rf$prefilt$Gene_SM) %>%
    pivot_wider(names_from = 'Timepoint', values_from = 'value')
  
  i = 1
  while(i<=500) {
    print(i)
    shuffle_indices = lapply(1:nrow(shuffle_2rf$prefilt), function(c) {
      sample(2:13, replace=FALSE)
    })
    
    lapply(seq_along(shuffle_indices), function(x) {
      shuffle_2rf$cellnorm$prefiltPreMeta[x, c(1, shuffle_indices[[x]])] %>%
        `colnames<-`(c('Gene_SM', seq(3,47,4)))}) %>% do.call(what='rbind') %>%
      write.csv(file = "shuffle_2rf_cell.csv", row.names=FALSE)
    
    meta2d(infile = "shuffle_2rf_cell.csv", filestyle = "csv", timepoints = "line1", minper=20, 
           maxper=28, outdir=getwd(), parallelize = TRUE)
    
    resCell = read.csv('meta2d_shuffle_2rf_cell.csv') %>%
      mutate(Gene = gsub('_.*', "", CycID), 
             Mode = gsub('.*_', "", CycID), 
             Sample = gsub('.*_C', 'C', CycID) %>% sub('_.*', '', .),
             p_cell = meta2d_pvalue,
             .keep='unused') %>%
      subset(select = c(Gene, Mode, Sample, p_cell))
    
    
    lapply(seq_along(shuffle_indices), function(x) {
      shuffle_2rf$AGGnorm$prefiltPreMeta[x, c(1, shuffle_indices[[x]])] %>%
        `colnames<-`(c('Gene_SM', seq(3,47,4)))}) %>% do.call(what='rbind') %>%
      write.csv(file = "shuffle_2rf_AGG.csv", row.names=FALSE)
    
    meta2d(infile = "shuffle_2rf_AGG.csv", filestyle = "csv", timepoints = "line1", minper=20, 
           maxper=28, outdir=getwd(), parallelize = TRUE)
    
    resAGG = read.csv('meta2d_shuffle_2rf_AGG.csv') %>%
      mutate(Gene = gsub('_.*', "", CycID), 
             Mode = gsub('.*_', "", CycID), 
             Sample = gsub('.*_C', 'C', CycID) %>% sub('_.*', '', .),
             p_AGG = meta2d_pvalue,
             .keep='unused') %>%
      subset(select = c(Gene, Mode, Sample, p_AGG))
    
    shuffle_2rf$retention[[i]] = merge(shuffle_2rf$prefilt, resCell) %>%
      merge(resAGG) %>% 
      {data.frame(.[,c('Gene', 'Sample', 'Mode', 'p_cell', 'p_AGG')], pass = ((.$ratio > 1.5 & .$p_cell < 0.05 & .$ratio_AGG > 1.5 & .$p_AGG < 0.05) | 
                                                                                (.$p_AGG < 0.05 & .$ratio_AGG > 2)))}
    
    shuffle_2rf$passes[[i]] = subset(shuffle_2rf$retention[[i]], pass==TRUE) %>% 
      {table(.$Sample, .$Mode)} %>%
      as.data.frame %>%
      `colnames<-`(c('Sample', 'Mode', 'count'))
    
    i = i+1
  }
  rm(i, shuffle_indices, resCell, resAGG)
  
  shuffle_2rf$counts <- lapply(seq_along(shuffle_2rf$passes), function(i) cbind(shuffle_2rf$passes[[i]], run=i)) %>% 
    c(list(data = cbind(subset(shuffle_2rf$prefilt, ((ratio > 1.5 & meta2d_pvalue < 0.05 & ratio_AGG > 1.5 & meta2d_pvalue_AGG_norm < 0.05) | 
                                                       (meta2d_pvalue_AGG_norm < 0.05 & ratio_AGG > 2))) %>% 
                          {table(.$Sample, .$Mode)} %>% 
                          as.data.frame %>% 
                          `colnames<-`(c('Sample', 'Mode', 'count')), run='data'))) %>% 
    do.call(what='rbind') %>%
    mutate(Sample = factor(Sample, levels=c('C3', 'C2')), 
           Mode = factor(Mode, levels=c('LD', 'DD')))
  
  shuffle_2rf$percentiles <- shuffle_2rf$counts %>%
    {split(., list(.$Sample, .$Mode))} %>%
    lapply(function(df) cbind(df, Rank = rank(df$count, ties.method='average')) %>%
             mutate(percentile = Rank/max(Rank)) %>%
             subset(run == 'data'))
  }

saveRDS(shuffle_2rf, 'shuffle_2rf.rds')

shuffle_2rf$plot = shuffle_2rf$counts %>%
  {ggplot() +
      geom_histogram(data=subset(., run != 'data'), aes(x = count), binwidth=1) +
      geom_vline(data=subset(., run == 'data'), aes(xintercept = count))} +
  labs(x = 'Cyclers', y = 'Frequency') +   xlim(0, NA) + 
  facet_grid(Mode ~ Sample)
ggsave(file='shuffle_2rf_500.svg', plot=shuffle_2rf$plot, width=5, height=5)


shuffle_2rf$percentiles

# derive p-values from shuffle
shuffle_2rf$meta2d_percentiles <- lapply(seq_along(shuffle_2rf$retention), function(i) {cbind(shuffle_2rf$retention[[i]], run=i)}) %>% 
  do.call(what='rbind') %>% 
  mutate(meta2d_pvalue = p_cell, 
         meta2d_pvalue_AGG_norm=p_AGG, 
         pass=NULL, 
         .keep='unused') %>% 
  {rbind(., cbind(shuffle_2rf$prefilt, run='data')[, colnames(.)] )}
  

oi <- shuffle_2rf$meta2d_percentiles %>%
  mutate(G_S_M = paste(Gene, Sample, Mode, sep='_')) %>%
  {split(., .$G_S_M)} %>% 
  lapply(function(df) {
    cbind(Gene=df$Gene[1], Sample=df$Sample[1], Mode=df$Mode[1],
          pct_cell = rank(df$meta2d_pvalue, ties.method='average') %>% {.[length(.)]/length(.)},
          pct_AGG = rank(df$meta2d_pvalue_AGG_norm, ties.method='average') %>% {.[length(.)]/length(.)},
          meta2d_pvalue = df[nrow(df), 'meta2d_pvalue'], 
          meta2d_pvalue_AGG_norm = df[nrow(df), 'meta2d_pvalue_AGG_norm']
    )
  }) %>%
  do.call(what='rbind') %>%
  as.data.frame

subset(oi, meta2d_pvalue_AGG_norm < pct_AGG)$meta2d_pvalue_AGG_norm %>% as.numeric %>% {hist(., breaks = length(.))}
subset(oi, meta2d_pvalue_AGG_norm > pct_AGG)$meta2d_pvalue_AGG_norm %>% as.numeric %>% {hist(., breaks = length(.))}

subset(oi, meta2d_pvalue_AGG_norm < pct_AGG) %>% dim
dim(oi)
# 90% of p-values increased
subset(oi, meta2d_pvalue_AGG_norm < 0.05) %>% dim
subset(oi, pct_AGG < 0.05) %>% dim

shuffle_2rf$meta2d_percentiles %>%
  merge(subset(shuffle_2rf$meta2d_percentiles, run=='data' & meta2d_pvalue_AGG_norm<0.05, select=c(Gene, Sample, Mode))) %>%
  #subset(Gene == 'mAChR-A') %>%
  {ggplot() + 
      geom_histogram(data=subset(., run != 'data'), aes(x = meta2d_pvalue_AGG_norm, y = after_stat(density)), binwidth=0.01) +
      geom_density(data=subset(., run == 'data'), aes(x = meta2d_pvalue_AGG_norm)) } +
  #geom_vline(data=subset(., run == 'data'), aes(xintercept = meta2d_pvalue_AGG_norm))} +
  facet_grid(Mode ~ Sample)


shuffle_2rf$meta2d_percentiles %>%
  {ggplot() + 
      geom_histogram(data=subset(., run != 'data'), aes(x = meta2d_pvalue), binwidth=0.01)} +
  facet_grid(Mode ~ Sample)

shuffle_2rf$meta2d_percentiles %>%
  merge(reconstituting_DFs_10$reconstituted[, c('Gene', 'Sample', 'Mode')]) %>%
  {ggplot() + 
      geom_vline(data=subset(., run == 'data'), aes(xintercept = meta2d_pvalue))} +
  facet_grid(Mode ~ Sample)

# plotting new distributions using ranks for all
hey <- shuffle_2rf$meta2d_percentiles %>%
  mutate(G_S_M = paste(Gene, Sample, Mode, sep='_')) %>%
  {split(., .$G_S_M)} %>% 
  lapply(function(df) {
    cbind(df, pct_cell = rank(df$meta2d_pvalue, ties.method='average')/nrow(df),
          pct_AGG = rank(df$meta2d_pvalue_AGG_norm, ties.method='average')/nrow(df))
  }) %>%
  do.call(what='rbind') %>%
  {split(., .$run)} %>% 
  lapply(function(df) {
    merge(df, shuffle_2rf$prefilt, by=c('Gene', 'Sample', 'Mode')) %>% 
      subset(((ratio > 1.5 & pct_cell < 0.05 & ratio_AGG > 1.5 & pct_AGG < 0.05) | 
                (pct_AGG < 0.05 & ratio_AGG > 2))) %>%
      {table(.$Sample, .$Mode)} %>%
      as.data.frame %>%
      `colnames<-`(c('Sample', 'Mode', 'count'))
  })
hey_plot <- hey %>%
  {lapply(names(.), function(name) cbind(.[[name]], run=name))} %>%
  do.call(what='rbind') %>%
  mutate(Sample = factor(Sample, levels=c('C3', 'C2')), 
         Mode = factor(Mode, levels=c('LD', 'DD'))) %>%
  {ggplot() +
      geom_histogram(data=subset(., run != 'data'), aes(x = count), binwidth=1) +
      geom_vline(data=subset(., run == 'data'), aes(xintercept = count))} +
  labs(x = 'Cyclers', y = 'Frequency') +   xlim(0, NA) + 
  facet_grid(Mode ~ Sample)

ggsave(file='shuffle_2rf_shuffleP.svg', plot=hey_plot, width=5, height=5)

hey_plot_2 <- shuffle_2rf$meta2d_percentiles %>%
  mutate(G_S_M = paste(Gene, Sample, Mode, sep='_')) %>%
  {split(., .$G_S_M)} %>% 
  lapply(function(df) {
    cbind(df, pct_cell = rank(df$meta2d_pvalue, ties.method='average')/nrow(df),
          pct_AGG = rank(df$meta2d_pvalue_AGG_norm, ties.method='average')/nrow(df))
  }) %>%
  do.call(what='rbind') %>%
  {split(., .$run)} %>% 
  lapply(function(df) {
    merge(df, shuffle_2rf$prefilt, by=c('Gene', 'Sample', 'Mode')) %>% 
      subset(ratio_AGG > 2 & pct_AGG < 0.05) %>%
      {table(.$Sample, .$Mode)} %>%
      as.data.frame %>%
      `colnames<-`(c('Sample', 'Mode', 'count'))
  }) %>%
  {lapply(names(.), function(name) cbind(.[[name]], run=name))} %>%
  do.call(what='rbind') %>%
  mutate(Sample = factor(Sample, levels=c('C3', 'C2')), 
         Mode = factor(Mode, levels=c('LD', 'DD'))) %>%
  {ggplot() +
      geom_histogram(data=subset(., run != 'data'), aes(x = count), binwidth=1) +
      geom_vline(data=subset(., run == 'data'), aes(xintercept = count))} +
  labs(x = 'Cyclers', y = 'Frequency') +   xlim(0, NA) + 
  facet_grid(Mode ~ Sample)
# ^ for only pct_AGG<0.05 but ratio_AGG>2
ggsave(file='shuffle_2rf_shuffleP_simple.svg', plot=hey_plot_2, width=5, height=5)


# make distribution plot using q-values (just using AGG)
subset(to_filter, second_pct>0.1 & ratio_AGG>1.5) %>%
  mutate(FDR_meta2d_AGG = p.adjust(meta2d_pvalue_AGG_norm, method='BH')) %>%
  subset(FDR_meta2d_AGG < 0.1) %>%
  {plot_grid(nrow=1, hist_phase(., CIRCULAR=FALSE, TRIM_TO_24=TRUE, binWidth=2, n_W=3), 
             hist_phase(., CIRCULAR=TRUE, TRIM_TO_24=TRUE, binWidth=2, n_W=3), 
             hmap_phase(.))}

subset(to_filter, second_pct>0.1 & ratio_AGG>1.5) %>%
  mutate(FDR_meta2d_AGG = p.adjust(meta2d_pvalue_AGG_norm, method='BH')) %>%
  subset(FDR_meta2d_AGG < 0.1) %>%
  subset(Sample == 'C3' & Mode == 'LD') %>%
  .$Gene %>%
  {plot_grid(nrow=1,
             avg_plotter_concat(.[1:25], data=subset(data_C2C3_AGG, sample=='C3' & mode1=='LD'), N_days = 2), 
             avg_plotter_concat(.[26:50], data=subset(data_C2C3_AGG, sample=='C3' & mode1=='LD'), N_days = 2), 
             avg_plotter_concat(.[51:77], data=subset(data_C2C3_AGG, sample=='C3' & mode1=='LD'), N_days = 2)
  )}
subset(to_filter, second_pct>0.1 & ratio_AGG>1.5) %>%
  mutate(FDR_meta2d_AGG = p.adjust(meta2d_pvalue_AGG_norm, method='BH')) %>%
  subset(FDR_meta2d_AGG < 0.1) %>%
  subset(Sample == 'C2' & Mode == 'DD') %>%
  .$Gene %>%
  avg_plotter_concat(data=subset(data_C2C3_AGG, sample=='C2' & mode1=='DD'), N_days = 2)


shuffle_2rf$AGGQ$counts <- lapply(shuffle_2rf$retention, function(df) {
  merge(df, shuffle_2rf$prefilt) %>%
    subset(ratio_AGG > 1.5) %>% 
    mutate(FDR_meta2d_AGG = p.adjust(p_AGG, method='BH')) %>% 
    subset(FDR_meta2d_AGG < 0.1) %>% 
    {table(.$Sample, .$Mode)} %>% 
    as.data.frame %>% 
    `colnames<-`(c('Sample', 'Mode', 'count'))
}) %>%
  {lapply(seq_along(.), function(i) cbind(.[[i]], run=i))} %>%
  c(list(data = cbind(subset(shuffle_2rf$prefilt, ratio_AGG > 1.5) %>%
                        mutate(FDR_meta2d_AGG = p.adjust(meta2d_pvalue_AGG_norm, method='BH')) %>%
                        subset(FDR_meta2d_AGG < 0.1) %>% 
                        {table(.$Sample, .$Mode)} %>% 
                        as.data.frame %>% 
                        `colnames<-`(c('Sample', 'Mode', 'count')), run='data'))) %>% 
  do.call(what='rbind') %>%
  mutate(Sample = factor(Sample, levels=c('C3', 'C2')), 
         Mode = factor(Mode, levels=c('LD', 'DD'))) %>%
  as.data.frame

shuffle_2rf$AGGQ$plot = shuffle_2rf$AGGQ$counts %>% as.data.frame %>%
  {ggplot() +
      geom_histogram(data=subset(., run != 'data'), aes(x = count), binwidth=1) +
      geom_vline(data=subset(., run == 'data'), aes(xintercept = count))} +
  labs(x = 'Cyclers', y = 'Frequency') +   xlim(0, NA) + 
  facet_grid(Mode ~ Sample)

ggsave(file='shuffle_2rf_500_Q1.5.svg', plot=shuffle_2rf$AGGQ$plot, width=5, height=5)


shuffle_2rf$AGGQ$percentiles = shuffle_2rf$AGGQ$counts %>%
  {split(., list(.$Sample, .$Mode))} %>%
  lapply(function(df) cbind(df, Rank = rank(df$count, ties.method='average')) %>%
           mutate(percentile = Rank/max(Rank)) %>%
           subset(run == 'data'))



# KW

shuffle_2rf$KW$counts <- lapply(shuffle_2rf$retention, function(df) {
  merge(df, to_filter[,c('Gene', 'Sample', 'Mode', 'second_pct', 'ratio', 'ratio_AGG', 'meta2d_pvalue', 'meta2d_pvalue_AGG_norm', 'KW_scores_norm')]) %>%
    subset(second_pct > 0.1 & ratio > 1.5 & p_cell < 0.05 & KW_scores_norm < 0.05) %>% 
    {table(.$Sample, .$Mode)} %>% 
    as.data.frame %>% 
    `colnames<-`(c('Sample', 'Mode', 'count'))
}) %>%
  {lapply(seq_along(.), function(i) cbind(.[[i]], run=i))} %>%
  c(list(data = cbind(subset(to_filter, second_pct > 0.1 & ratio > 1.5 & meta2d_pvalue < 0.05 & KW_scores_norm < 0.05) %>%
                        {table(.$Sample, .$Mode)} %>% 
                        as.data.frame %>% 
                        `colnames<-`(c('Sample', 'Mode', 'count')), run='data'))) %>% 
  do.call(what='rbind') %>%
  mutate(Sample = factor(Sample, levels=c('C3', 'C2')), 
         Mode = factor(Mode, levels=c('LD', 'DD'))) %>%
  as.data.frame

shuffle_2rf$KW$plot = shuffle_2rf$KW$counts %>% as.data.frame %>%
  {ggplot() +
      geom_histogram(data=subset(., run != 'data'), aes(x = count), binwidth=1) +
      geom_vline(data=subset(., run == 'data'), aes(xintercept = count))} +
  labs(x = 'Cyclers', y = 'Frequency') +   xlim(0, NA) + 
  facet_grid(Mode ~ Sample)

ggsave(file='shuffle_2rf_500_KW.svg', plot=shuffle_2rf$KW$plot, width=5, height=5)

shuffle_2rf$KW$percentiles = shuffle_2rf$KW$counts %>%
  {split(., list(.$Sample, .$Mode))} %>%
  lapply(function(df) cbind(df, Rank = rank(df$count, ties.method='average')) %>%
           mutate(percentile = Rank/max(Rank)) %>%
           subset(run == 'data'))



# KW p w/ Q (cell)

shuffle_2rf$KWp_2dCellQ$counts <- lapply(shuffle_2rf$retention, function(df) {
  merge(df, to_filter[,c('Gene', 'Sample', 'Mode', 'second_pct', 'ratio', 'ratio_AGG', 'meta2d_pvalue', 'meta2d_pvalue_AGG_norm', 'KW_scores_norm')]) %>%
    subset(second_pct > 0.1 & ratio > 1.5 & KW_scores_norm < 0.05) %>% 
    mutate(FDR_meta2d = p.adjust(p_cell, method='BH')) %>% 
    subset(FDR_meta2d < 0.1) %>% 
    {if (!nrow(.)) data.frame(c('C2', 'C3', 'C2', 'C3'), c('DD', 'DD', 'LD', 'LD'), rep(0, 4))
      else table(.$Sample, .$Mode)} %>% 
    as.data.frame %>% 
    `colnames<-`(c('Sample', 'Mode', 'count'))
}) %>%
  {lapply(seq_along(.), function(i) cbind(.[[i]], run=i))} %>%
  c(list(data = cbind(subset(to_filter, second_pct > 0.1 & ratio > 1.5 & KW_scores_norm < 0.05) %>% 
                        mutate(FDR_meta2d = p.adjust(meta2d_pvalue, method='BH')) %>%
                        subset(FDR_meta2d < 0.1) %>% 
                        {table(.$Sample, .$Mode)} %>% 
                        as.data.frame %>% 
                        `colnames<-`(c('Sample', 'Mode', 'count')), run='data'))) %>% 
  do.call(what='rbind') %>%
  mutate(Sample = factor(Sample, levels=c('C3', 'C2')), 
         Mode = factor(Mode, levels=c('LD', 'DD'))) %>%
  as.data.frame

shuffle_2rf$KWp_2dCellQ$plot = shuffle_2rf$KWp_2dCellQ$counts %>% as.data.frame %>%
  {ggplot() +
      geom_histogram(data=subset(., run != 'data'), aes(x = count), binwidth=1) +
      geom_vline(data=subset(., run == 'data'), aes(xintercept = count))} +
  labs(x = 'Cyclers', y = 'Frequency') +   xlim(0, NA) + 
  facet_grid(Mode ~ Sample)

ggsave(file='shuffle_2rf_500_KWp_2dCellQ.svg', plot=shuffle_2rf$KWp_2dCellQ$plot, width=5, height=5)


shuffle_2rf$KWp_2dCellQ$percentiles = shuffle_2rf$KWp_2dCellQ$counts %>%
  {split(., list(.$Sample, .$Mode))} %>%
  lapply(function(df) cbind(df, Rank = rank(df$count, ties.method='average')) %>%
           mutate(percentile = Rank/max(Rank)) %>%
           subset(run == 'data'))

# KW Q w/ Q (cell)

shuffle_2rf$KWq_2dCellQ$counts <- lapply(shuffle_2rf$retention, function(df) {
  merge(df, to_filter[,c('Gene', 'Sample', 'Mode', 'second_pct', 'ratio', 'ratio_AGG', 'meta2d_pvalue', 'meta2d_pvalue_AGG_norm', 'KW_scores_norm')]) %>%
    subset(second_pct > 0.1 & ratio > 1.5) %>% 
    mutate(FDR_KW = p.adjust(KW_scores_norm, method='BH')) %>% 
    subset(FDR_KW < 0.1) %>%     
    mutate(FDR_meta2d = p.adjust(p_cell, method='BH')) %>% 
    subset(FDR_meta2d < 0.1) %>% 
    {if (!nrow(.)) data.frame(c('C2', 'C3', 'C2', 'C3'), c('DD', 'DD', 'LD', 'LD'), rep(0, 4))
      else table(.$Sample, .$Mode)} %>% 
    as.data.frame %>% 
    `colnames<-`(c('Sample', 'Mode', 'count'))
}) %>%
  {lapply(seq_along(.), function(i) cbind(.[[i]], run=i))} %>%
  c(list(data = cbind(subset(to_filter, second_pct > 0.1 & ratio > 1.5) %>% 
                        mutate(FDR_KW = p.adjust(KW_scores_norm, method='BH')) %>% 
                        subset(FDR_KW < 0.1) %>% 
                        mutate(FDR_meta2d = p.adjust(meta2d_pvalue, method='BH')) %>%
                        subset(FDR_meta2d < 0.1) %>% 
                        {table(.$Sample, .$Mode)} %>% 
                        as.data.frame %>% 
                        `colnames<-`(c('Sample', 'Mode', 'count')), run='data'))) %>% 
  do.call(what='rbind') %>%
  mutate(Sample = factor(Sample, levels=c('C3', 'C2')), 
         Mode = factor(Mode, levels=c('LD', 'DD'))) %>%
  as.data.frame

shuffle_2rf$KWq_2dCellQ$plot = shuffle_2rf$KWq_2dCellQ$counts %>% as.data.frame %>%
  {ggplot() +
      geom_histogram(data=subset(., run != 'data'), aes(x = count), binwidth=1) +
      geom_vline(data=subset(., run == 'data'), aes(xintercept = count))} +
  labs(x = 'Cyclers', y = 'Frequency') +   xlim(0, NA) + 
  facet_grid(Mode ~ Sample)

ggsave(file='shuffle_2rf_500_KWq_2dCellQ.svg', plot=shuffle_2rf$KWq_2dCellQ$plot, width=5, height=5)

shuffle_2rf$KWq_2dCellQ$percentiles = shuffle_2rf$KWq_2dCellQ$counts %>%
  {split(., list(.$Sample, .$Mode))} %>%
  lapply(function(df) cbind(df, Rank = rank(df$count, ties.method='average')) %>%
           mutate(percentile = Rank/max(Rank)) %>%
           subset(run == 'data'))
shuffle_2rf$KWq_2dCellQ$plot
shuffle_2rf$KWq_2dCellQ$percentiles
####    ####    ####    ####    ####    ####    ####    ####    ####    ####   


# SHUFFLE CELLS
################################################################################
subset(data_C2C3, sample == 'C3' & mode1 == 'LD', features = unique(res10$C3_LD$filtered$Gene))@assays$RNA$counts %>%
  write.table('counts_C3LD_cycling.txt', quote=FALSE)
subset(data_C2C3, sample == 'C2' & mode1 == 'LD', features = unique(res10$C2_LD$filtered$Gene))@assays$RNA$counts %>%
  write.table('counts_C2LD_cycling.txt', quote=FALSE) 
subset(data_C2C3, sample == 'C3' & mode1 == 'DD', features = unique(res10$C3_DD$filtered$Gene))@assays$RNA$counts %>%
  write.table('counts_C3DD_cycling.txt', quote=FALSE)
subset(data_C2C3, sample == 'C2' & mode1 == 'DD', features = unique(res10$C2_DD$filtered$Gene))@assays$RNA$counts %>%
  write.table('counts_C2DD_cycling.txt', quote=FALSE)


subset(data_C2C3, sample == 'C3' & mode1 == 'LD', 
       features = unique(subset(to_filter, Sample == 'C3' & Mode == 'LD' & second_pct > 0.1)$Gene))@assays$RNA$counts %>%
  write.table('counts_C3LD_expressed.txt', quote=FALSE)
subset(data_C2C3, sample == 'C2' & mode1 == 'LD',
       features = unique(subset(to_filter, Sample == 'C2' & Mode == 'LD' & second_pct > 0.1)$Gene))@assays$RNA$counts %>%
  write.table('counts_C2LD_expressed.txt', quote=FALSE) 
subset(data_C2C3, sample == 'C3' & mode1 == 'DD',
       features = unique(subset(to_filter, Sample == 'C3' & Mode == 'DD' & second_pct > 0.1)$Gene))@assays$RNA$counts %>%
  write.table('counts_C3DD_expressed.txt', quote=FALSE)
subset(data_C2C3, sample == 'C2' & mode1 == 'DD',
       features = unique(subset(to_filter, Sample == 'C2' & Mode == 'DD' & second_pct > 0.1)$Gene))@assays$RNA$counts %>%
  write.table('counts_C2DD_expressed.txt', quote=FALSE)

# 2-DAY PER-ROW SHUFFLE, FULL CYCLING ANALYSIS
{set.seed(1)
  #shuffle_2cells <- list()
  #shuffle_2cells$prefilt <- subset(to_filter, second_pct > 0.10, select=c(Gene, Sample, Mode, second_pct, ratio, ratio_AGG, meta2d_pvalue, meta2d_pvalue_AGG_norm)) %>%
  #  mutate(Gene_SM = paste(Gene, Sample, Mode, sep='_'))
  
  split_object <- SplitObject(data_C2C3, split.by = "sample_mode")
  data_lengths <- lapply(split_object, ncol)
  
  i = 300
  while(i<=300) {
    print(i)
    
    set.seed(i)
    indices <- lapply(data_lengths, function(x) sample(1:x, replace=FALSE))
    
    shuffled_data <- list()
    for (j in 1:length(split_object)) {
      Data <- split_object[[j]]
      indices <- sample(1:ncol(Data), replace=FALSE)
      Data@meta.data[['sample_mode_group2']] <- Data@meta.data[['sample_mode_group2']][indices]
      Data@meta.data[['group2']] <- Data[['sample_mode_group2']] %>% 
        mutate(sample = gsub('.*-', '', sample_mode_group2), .keep='unused') %>% unlist
      shuffled_data[[j]] <- Data
    }
    shuffled_data = merge(shuffled_data[[1]], shuffled_data[2:4])
    
    shuffled_data_AGG <- AggregateExpression(shuffled_data, return.seurat=TRUE, group.by='sample_mode_group2')
    shuffled_data_AGG@meta.data[['sample']] <- shuffled_data_AGG[[c('sample_mode_group2')]] %>% 
      mutate(sample = gsub('-.*', '', sample_mode_group2), .keep='unused') %>% unlist %>%
      factor(levels=c('C3', 'C2'))
    shuffled_data_AGG@meta.data[['mode1']] <- shuffled_data_AGG[[c('sample_mode_group2')]] %>% 
      mutate(sample = sub('D-.*', 'D', sample_mode_group2) %>% sub('.*-', '', .), .keep='unused') %>% unlist %>%
      factor(levels=c('LD', 'DD'))
    shuffled_data_AGG@meta.data[['group2']] <- shuffled_data_AGG[[c('sample_mode_group2')]] %>% 
      mutate(sample = gsub('.*-', '', sample_mode_group2), .keep='unused') %>% unlist %>%
      factor(levels = c(paste0('CT', seq(3, 47, 4)), paste0('ZT', seq(3, 47, 4))))
    
    preMeta_cell <- AverageExpression(shuffled_data, group.by = c("group2", "sample", "mode1"), layer="data")$RNA %>% as.data.frame %>%
      {mutate(., Gene = rownames(.))} %>%
      pivot_longer(cols=-'Gene', names_to='TP_sample_mode') %>%
      mutate(Timepoint = TP_sample_mode %>% gsub('_.*', '', .) %>% gsub('.*T', '', .), 
             Gene_SM = paste(Gene, TP_sample_mode %>% gsub('.*_C', 'C', .), sep = '_'), .keep='unused') %>%
      subset(Gene_SM %in% shuffle_2cells$prefilt$Gene_SM) %>%
      pivot_wider(names_from = 'Timepoint', values_from = 'value')
    
    preMeta_AGG <- AverageExpression(shuffled_data_AGG, group.by = c("group2", "sample", "mode1"), layer="data")$RNA %>% as.data.frame %>%
      {mutate(., Gene = rownames(.))} %>%
      pivot_longer(cols=-'Gene', names_to='TP_sample_mode') %>%
      mutate(Timepoint = TP_sample_mode %>% gsub('_.*', '', .) %>% gsub('.*T', '', .), 
             Gene_SM = paste(Gene, TP_sample_mode %>% gsub('.*_C', 'C', .), sep = '_'), .keep='unused') %>%
      subset(Gene_SM %in% shuffle_2cells$prefilt$Gene_SM) %>%
      pivot_wider(names_from = 'Timepoint', values_from = 'value')
    
    shuffle_indices = lapply(1:nrow(shuffle_2cells$prefilt), function(c) {
      sample(2:13, replace=FALSE)
    })
    
    lapply(seq_along(shuffle_indices), function(x) {
      preMeta_cell[x, c(1, shuffle_indices[[x]])] %>%
        `colnames<-`(c('Gene_SM', seq(3,47,4)))}) %>% do.call(what='rbind') %>%
      write.csv(file = "shuffle_2cells_cell.csv", row.names=FALSE)
    
    meta2d(infile = "shuffle_2cells_cell.csv", filestyle = "csv", timepoints = "line1", minper=20, 
           maxper=28, outdir=getwd(), parallelize = TRUE)
    
    resCell = merge(read.csv('meta2d_shuffle_2cells_cell.csv'), 
                    data.frame(CycID = preMeta_cell$Gene_SM, ratio = preMeta_cell[, 2:13] %>% apply(1, function(x) max(x)/min(x)))) %>% 
      mutate(Gene = gsub('_.*', "", CycID), 
             Mode = gsub('.*_', "", CycID), 
             Sample = gsub('.*_C', 'C', CycID) %>% sub('_.*', '', .),
             p_cell = meta2d_pvalue,
             .keep='unused') %>%
      subset(select = c(Gene, Mode, Sample, p_cell, ratio))
    
    
    lapply(seq_along(shuffle_indices), function(x) {
      preMeta_AGG[x, c(1, shuffle_indices[[x]])] %>%
        `colnames<-`(c('Gene_SM', seq(3,47,4)))}) %>% do.call(what='rbind') %>%
      write.csv(file = "shuffle_2cells_AGG.csv", row.names=FALSE)
    
    meta2d(infile = "shuffle_2cells_AGG.csv", filestyle = "csv", timepoints = "line1", minper=20, 
           maxper=28, outdir=getwd(), parallelize = TRUE)
    
    resAGG = merge(read.csv('meta2d_shuffle_2cells_AGG.csv'), 
                    data.frame(CycID = preMeta_AGG$Gene_SM, ratio_AGG = preMeta_AGG[, 2:13] %>% apply(1, function(x) max(x)/min(x)))) %>%
      mutate(Gene = gsub('_.*', "", CycID), 
             Mode = gsub('.*_', "", CycID), 
             Sample = gsub('.*_C', 'C', CycID) %>% sub('_.*', '', .),
             p_AGG= meta2d_pvalue,
             .keep='unused') %>%
      subset(select = c(Gene, Mode, Sample, p_AGG, ratio_AGG))
    
    shuffle_2cells$retention[[i]] = merge(resAGG, resCell) %>% 
      {data.frame(.[,c('Gene', 'Sample', 'Mode', 'p_cell', 'p_AGG', 'ratio', 'ratio_AGG')], 
                  pass = ((.$ratio > 1.5 & .$p_cell < 0.05 & .$ratio_AGG > 1.5 & .$p_AGG < 0.05) | 
                                                                                (.$p_AGG < 0.05 & .$ratio_AGG > 2)))}
    
    shuffle_2cells$passes[[i]] = subset(shuffle_2cells$retention[[i]], pass==TRUE) %>% 
      {table(.$Sample, .$Mode)} %>%
      as.data.frame %>%
      `colnames<-`(c('Sample', 'Mode', 'count'))
    
    i = i+1
  }
  rm(i, x, data_lengths, split_object, indices, Data, shuffle_indices, resCell, resAGG, shuffled_data, shuffled_data_AGG, preMeta_cell, preMeta_AGG)
  
  shuffle_2cells$counts <- lapply(seq_along(shuffle_2cells$passes), function(i) cbind(shuffle_2cells$passes[[i]], run=i)) %>% 
    c(list(data = cbind(subset(shuffle_2cells$prefilt, ((ratio > 1.5 & meta2d_pvalue < 0.05 & ratio_AGG > 1.5 & meta2d_pvalue_AGG_norm < 0.05) | 
                                                       (meta2d_pvalue_AGG_norm < 0.05 & ratio_AGG > 2))) %>% 
                          {table(.$Sample, .$Mode)} %>% 
                          as.data.frame %>% 
                          `colnames<-`(c('Sample', 'Mode', 'count')), run='data'))) %>% 
    do.call(what='rbind') %>%
    mutate(Sample = factor(Sample, levels=c('C3', 'C2')), 
           Mode = factor(Mode, levels=c('LD', 'DD')))
  
  shuffle_2cells$percentiles <- shuffle_2cells$counts %>%
    {split(., list(.$Sample, .$Mode))} %>%
    lapply(function(df) cbind(df, Rank = rank(df$count, ties.method='average')) %>%
             mutate(percentile = Rank/max(Rank)) %>%
             subset(run == 'data'))
}

saveRDS(shuffle_2cells, 'shuffle_2cells.rds')

shuffle_2cells <- readRDS('shuffle_2cells.rds')

shuffle_2cells$plot = shuffle_2cells$counts %>% as.data.frame %>%
  {merge(., group_by(subset(., run != 'data'), Sample, Mode) %>% summarise(Mean = mean(count)))} %>% 
  {ggplot(data=.) +
      geom_histogram(data=subset(., run != 'data'), aes(x = count), binwidth=1, fill="gray30",color="black") +
      geom_vline(data=subset(., run == 'data'), aes(xintercept = count)) +
      geom_vline(data=., aes(xintercept = Mean), linetype = "dashed", linewidth = 1) +
      geom_shadowtext(data=subset(., run==1), aes(label=round(Mean),y=Inf,x=Mean), vjust=4,hjust=1.25,col='white',size=4) + 
      geom_text(data=subset(., run == 'data'), aes(label=round(count),y=Inf,x=count),vjust=2,hjust=1.25,col='black',size=4)} +
  labs(x = 'Cyclers', y = 'Frequency') +   xlim(0, NA) + 
  facet_grid(Mode ~ Sample, scales = "free", axes = "all", axis.labels = "all")

ggsave(file='shuffle_2cells.svg', plot=shuffle_2cells$plot, width=5, height=5)
ggsave(file='shuffle_2cells.pdf', plot=shuffle_2cells$plot, width=5, height=5)
ggsave(file='shuffle_2cells.jpg', plot=shuffle_2cells$plot, width=5, height=5)

# FDR q-values
shuffle_2cells$AGGQ$counts <- lapply(shuffle_2cells$retention, function(df) {
  merge(df, subset(shuffle_2cells$prefilt, select = -c(ratio, ratio_AGG))) %>%
    subset(ratio_AGG > 1.5) %>% 
    mutate(FDR_meta2d_AGG = p.adjust(p_AGG, method='BH')) %>% 
    subset(FDR_meta2d_AGG < 0.1) %>% 
    {table(.$Sample, .$Mode)} %>% 
    as.data.frame %>% 
    `colnames<-`(c('Sample', 'Mode', 'count'))
}) %>% 
  {lapply(seq_along(.), function(i) cbind(.[[i]], run=i))} %>%
  c(list(data = cbind(subset(shuffle_2cells$prefilt, ratio_AGG > 1.5
                             ) %>%
                        mutate(FDR_meta2d_AGG = p.adjust(meta2d_pvalue_AGG_norm, method='BH')) %>%
                        subset(FDR_meta2d_AGG < 0.1) %>% 
                        {table(.$Sample, .$Mode)} %>% 
                        as.data.frame %>% 
                        `colnames<-`(c('Sample', 'Mode', 'count')), run='data'))) %>% 
  do.call(what='rbind') %>%
  mutate(Sample = factor(Sample, levels=c('C3', 'C2')), 
         Mode = factor(Mode, levels=c('LD', 'DD'))) %>%
  as.data.frame

shuffle_2cells$AGGQ$plot = shuffle_2cells$AGGQ$counts %>% 
  {merge(., group_by(subset(., run != 'data'), Sample, Mode) %>% summarise(Mean = mean(count)))} %>% 
  {ggplot(data=.) +
      geom_histogram(data=subset(., run != 'data'), aes(x = count), binwidth=1, fill="gray30",color="black") +
      geom_vline(data=subset(., run == 'data'), aes(xintercept = count)) +
      geom_vline(data=., aes(xintercept = Mean), linetype = "dashed", linewidth = 1) +
      geom_shadowtext(data=subset(., run==1), aes(label=round(Mean),y=Inf,x=Mean), vjust=4,hjust=1.25,col='white',size=4) + 
      geom_text(data=subset(., run == 'data'), aes(label=round(count),y=Inf,x=count),vjust=2,hjust=1.25,col='black',size=4)} +
  labs(x = 'Cyclers', y = 'Frequency') +   xlim(0, NA) + 
  facet_grid(Mode ~ Sample, scales = "free", axes = "all", axis.labels = "all")


ggsave(file='shuffle_2cells_Q1.5.svg', plot=shuffle_2cells$AGGQ$plot, width=5, height=5)
ggsave(file='shuffle_2cells_Q1.5.pdf', plot=shuffle_2cells$AGGQ$plot, width=5, height=5)
ggsave(file='shuffle_2cells_Q1.5.jpg', plot=shuffle_2cells$AGGQ$plot, width=5, height=5)



shuffle_2cells$AGGQ$percentiles = shuffle_2cells$AGGQ$counts %>%
  {split(., list(.$Sample, .$Mode))} %>%
  lapply(function(df) cbind(df, Rank = rank(df$count, ties.method='average')) %>%
           mutate(percentile = Rank/max(Rank)) %>%
           subset(run == 'data'))

# KW

shuffle_2cells$KW$counts <- lapply(shuffle_2cells$retention, function(df) {
  merge(df, to_filter[,c('Gene', 'Sample', 'Mode', 'second_pct', 'meta2d_pvalue', 'meta2d_pvalue_AGG_norm', 'KW_scores_norm')]) %>%
    subset(second_pct > 0.1 & p_cell < 0.05 & KW_scores_norm < 0.05 & ratio > 1.5
           ) %>% 
    {table(.$Sample, .$Mode)} %>% 
    as.data.frame %>% 
    `colnames<-`(c('Sample', 'Mode', 'count'))
}) %>%
  {lapply(seq_along(.), function(i) cbind(.[[i]], run=i))} %>%
  c(list(data = cbind(subset(to_filter, second_pct > 0.1 & meta2d_pvalue < 0.05 & KW_scores_norm < 0.05 & ratio > 1.5
           ) %>% 
                        {table(.$Sample, .$Mode)} %>% 
                        as.data.frame %>% 
                        `colnames<-`(c('Sample', 'Mode', 'count')), run='data'))) %>% 
  do.call(what='rbind') %>%
  mutate(Sample = factor(Sample, levels=c('C3', 'C2')), 
         Mode = factor(Mode, levels=c('LD', 'DD'))) %>%
  as.data.frame

shuffle_2cells$KW$plot = shuffle_2cells$KW$counts %>%
  {merge(., group_by(subset(., run != 'data'), Sample, Mode) %>% summarise(Mean = mean(count)))} %>% 
  {ggplot(data=.) +
      geom_histogram(data=subset(., run != 'data'), aes(x = count), binwidth=1, fill="gray30",color="black") +
      geom_vline(data=subset(., run == 'data'), aes(xintercept = count)) +
      geom_vline(data=., aes(xintercept = Mean), linetype = "dashed", linewidth = 1) +
      geom_shadowtext(data=subset(., run==1), aes(label=round(Mean),y=Inf,x=Mean), vjust=4,hjust=1.25,col='white',size=4) + 
      geom_text(data=subset(., run == 'data'), aes(label=round(count),y=Inf,x=count),vjust=2,hjust=1.25,col='black',size=4)} +
  labs(x = 'Cyclers', y = 'Frequency') +   xlim(0, NA) + 
  facet_grid(Mode ~ Sample, scales = "free", axes = "all", axis.labels = "all")

ggsave(file='shuffle_2cells_KW.svg', plot=shuffle_2cells$KW$plot, width=5, height=5)
ggsave(file='shuffle_2cells_KW.pdf', plot=shuffle_2cells$KW$plot, width=5, height=5)
ggsave(file='shuffle_2cells_KW.jpg', plot=shuffle_2cells$KW$plot, width=5, height=5)

shuffle_2cells$KW$percentiles = shuffle_2cells$KW$counts %>%
  {split(., list(.$Sample, .$Mode))} %>%
  lapply(function(df) cbind(df, Rank = rank(df$count, ties.method='average')) %>%
           mutate(percentile = Rank/max(Rank)) %>%
           subset(run == 'data'))



# KW p w/ Q (cell)

shuffle_2cells$KWp_2dCellQ$counts <- lapply(shuffle_2cells$retention, function(df) {
  merge(df, to_filter[,c('Gene', 'Sample', 'Mode', 'second_pct', 'meta2d_pvalue', 'meta2d_pvalue_AGG_norm', 'KW_scores_norm')]) %>%
    subset(second_pct > 0.1 & KW_scores_norm < 0.05 & ratio > 1.5 
           ) %>% 
    mutate(FDR_meta2d = p.adjust(p_cell, method='BH')) %>% 
    subset(FDR_meta2d < 0.1) %>% 
    {if (!nrow(.)) data.frame(c('C2', 'C3', 'C2', 'C3'), c('DD', 'DD', 'LD', 'LD'), rep(0, 4))
      else table(.$Sample, .$Mode)} %>% 
    as.data.frame %>% 
    `colnames<-`(c('Sample', 'Mode', 'count'))
}) %>%
  {lapply(seq_along(.), function(i) cbind(.[[i]], run=i))} %>%
  c(list(data = cbind(subset(to_filter, second_pct > 0.1 & KW_scores_norm < 0.05 & ratio > 1.5 
                             ) %>% 
                        mutate(FDR_meta2d = p.adjust(meta2d_pvalue, method='BH')) %>%
                        subset(FDR_meta2d < 0.1) %>% 
                        {table(.$Sample, .$Mode)} %>% 
                        as.data.frame %>% 
                        `colnames<-`(c('Sample', 'Mode', 'count')), run='data'))) %>% 
  do.call(what='rbind') %>%
  mutate(Sample = factor(Sample, levels=c('C3', 'C2')), 
         Mode = factor(Mode, levels=c('LD', 'DD'))) %>%
  as.data.frame


shuffle_2cells$KWp_2dCellQ$plot = shuffle_2cells$KWp_2dCellQ$counts %>%
  {merge(., group_by(subset(., run != 'data'), Sample, Mode) %>% summarise(Mean = mean(count)))} %>% 
  {ggplot(data=.) +
      geom_histogram(data=subset(., run != 'data'), aes(x = count), binwidth=1, fill="gray30",color="black") +
      geom_vline(data=subset(., run == 'data'), aes(xintercept = count)) +
      geom_vline(data=., aes(xintercept = Mean), linetype = "dashed", linewidth = 1) +
      geom_shadowtext(data=subset(., run==1), aes(label=round(Mean),y=Inf,x=Mean), vjust=4,hjust=1.25,col='white',size=4) + 
      geom_text(data=subset(., run == 'data'), aes(label=round(count),y=Inf,x=count),vjust=2,hjust=1.25,col='black',size=4)} +
  labs(x = 'Cyclers', y = 'Frequency') +   xlim(0, NA) + 
  facet_grid(Mode ~ Sample, scales = "free", axes = "all", axis.labels = "all")

ggsave(file='shuffle_2cells_KWp_2dCellQ.svg', plot=shuffle_2cells$KWp_2dCellQ$plot, width=5, height=5)
ggsave(file='shuffle_2cells_KWp_2dCellQ.pdf', plot=shuffle_2cells$KWp_2dCellQ$plot, width=5, height=5)
ggsave(file='shuffle_2cells_KWp_2dCellQ.jpg', plot=shuffle_2cells$KWp_2dCellQ$plot, width=5, height=5)


shuffle_2cells$KWp_2dCellQ$percentiles = shuffle_2cells$KWp_2dCellQ$counts %>%
  {split(., list(.$Sample, .$Mode))} %>%
  lapply(function(df) cbind(df, Rank = rank(df$count, ties.method='average')) %>%
           mutate(percentile = Rank/max(Rank)) %>%
           subset(run == 'data'))


subset(to_filter, second_pct > 0.1 & ratio > 1.5 & KW_scores_norm < 0.05) %>% 
  mutate(FDR_meta2d = p.adjust(meta2d_pvalue, method='BH')) %>%
  subset(FDR_meta2d < 0.1) %>% subset(Mode == 'DD' & Sample == 'C3') %>% {.$Gene}



# KW Q w/ p (cell) - best (DD 49-9, 30-7, 100 cyclers)

shuffle_2cells$KWQ_2dCellP$counts <- lapply(shuffle_2cells$retention, function(df) {
  merge(df, to_filter[,c('Gene', 'Sample', 'Mode', 'second_pct', 'meta2d_pvalue', 'meta2d_pvalue_AGG_norm', 'KW_scores_norm')]) %>%
    subset(second_pct > 0.1 & p_cell < 0.05 & ratio > 1.5 
    ) %>% 
    mutate(FDR_meta2d = p.adjust(KW_scores_norm, method='BH')) %>% 
    subset(FDR_meta2d < 0.1) %>% 
    {if (!nrow(.)) data.frame(c('C2', 'C3', 'C2', 'C3'), c('DD', 'DD', 'LD', 'LD'), rep(0, 4))
      else table(.$Sample, .$Mode)} %>% 
    as.data.frame %>% 
    `colnames<-`(c('Sample', 'Mode', 'count'))
}) %>%
  {lapply(seq_along(.), function(i) cbind(.[[i]], run=i))} %>%
  c(list(data = cbind(subset(to_filter, second_pct > 0.1 & meta2d_pvalue < 0.05 & ratio > 1.5 
  ) %>% 
    mutate(FDR_meta2d = p.adjust(KW_scores_norm, method='BH')) %>%
    subset(FDR_meta2d < 0.1) %>% 
    {table(.$Sample, .$Mode)} %>% 
    as.data.frame %>% 
    `colnames<-`(c('Sample', 'Mode', 'count')), run='data'))) %>% 
  do.call(what='rbind') %>%
  mutate(Sample = factor(Sample, levels=c('C3', 'C2')), 
         Mode = factor(Mode, levels=c('LD', 'DD'))) %>%
  as.data.frame

shuffle_2cells$KWQ_2dCellP$plot = shuffle_2cells$KWQ_2dCellP$counts %>%
  {merge(., group_by(subset(., run != 'data'), Sample, Mode) %>% summarise(Mean = mean(count)))} %>% 
  {ggplot(data=.) +
      geom_histogram(data=subset(., run != 'data'), aes(x = count), binwidth=1, fill="gray30",color="black") +
      geom_vline(data=subset(., run == 'data'), aes(xintercept = count)) +
      geom_vline(data=., aes(xintercept = Mean), linetype = "dashed", linewidth = 1) +
      geom_shadowtext(data=subset(., run==1), aes(label=round(Mean),y=Inf,x=Mean), vjust=4,hjust=1.25,col='white',size=4) + 
      geom_text(data=subset(., run == 'data'), aes(label=round(count),y=Inf,x=count),vjust=2,hjust=1.25,col='black',size=4)} +
  labs(x = 'Cyclers', y = 'Frequency') +   xlim(0, NA) + 
  facet_grid(Mode ~ Sample, scales = "free", axes = "all", axis.labels = "all")

my_filtered <- subset(to_filter, second_pct > 0.1 & meta2d_pvalue < 0.05 & ratio > 1.5) %>% 
  mutate(FDR_meta2d = p.adjust(KW_scores_norm, method='BH')) %>%
  subset(FDR_meta2d < 0.1)


# KW Q w/ Q (cell)

shuffle_2cells$KWq_2dCellQ$counts <- lapply(shuffle_2cells$retention, function(df) {
  merge(df, to_filter[,c('Gene', 'Sample', 'Mode', 'second_pct', 'meta2d_pvalue', 'meta2d_pvalue_AGG_norm', 'KW_scores_norm')]) %>%
    subset(second_pct > 0.1 & ratio > 1.5) %>% 
    mutate(FDR_KW = p.adjust(KW_scores_norm, method='BH')) %>% 
    subset(FDR_KW < 0.1) %>%     
    mutate(FDR_meta2d = p.adjust(p_cell, method='BH')) %>% 
    subset(FDR_meta2d < 0.1) %>% 
    {if (!nrow(.)) data.frame(c('C2', 'C3', 'C2', 'C3'), c('DD', 'DD', 'LD', 'LD'), rep(0, 4))
      else table(.$Sample, .$Mode)} %>% 
    as.data.frame %>% 
    `colnames<-`(c('Sample', 'Mode', 'count'))
}) %>%
  {lapply(seq_along(.), function(i) cbind(.[[i]], run=i))} %>%
  c(list(data = cbind(subset(to_filter, second_pct > 0.1 & ratio > 1.5) %>% 
                        mutate(FDR_KW = p.adjust(KW_scores_norm, method='BH')) %>% 
                        subset(FDR_KW < 0.1) %>% 
                        mutate(FDR_meta2d = p.adjust(meta2d_pvalue, method='BH')) %>%
                        subset(FDR_meta2d < 0.1) %>% 
                        {table(.$Sample, .$Mode)} %>% 
                        as.data.frame %>% 
                        `colnames<-`(c('Sample', 'Mode', 'count')), run='data'))) %>% 
  do.call(what='rbind') %>%
  mutate(Sample = factor(Sample, levels=c('C3', 'C2')), 
         Mode = factor(Mode, levels=c('LD', 'DD'))) %>%
  as.data.frame

shuffle_2cells$KWq_2dCellQ$plot = shuffle_2cells$KWq_2dCellQ$counts %>% 
  {merge(., group_by(subset(., run != 'data'), Sample, Mode) %>% summarise(Mean = mean(count)))} %>% 
  {ggplot(data=.) +
      geom_histogram(data=subset(., run != 'data'), aes(x = count), binwidth=1, fill="gray30",color="black") +
      geom_vline(data=subset(., run == 'data'), aes(xintercept = count)) +
      geom_vline(data=., aes(xintercept = Mean), linetype = "dashed", linewidth = 1) +
      geom_shadowtext(data=subset(., run==1), aes(label=round(Mean),y=Inf,x=Mean), vjust=4,hjust=1.25,col='white',size=4) + 
      geom_text(data=subset(., run == 'data'), aes(label=round(count),y=Inf,x=count),vjust=2,hjust=1.25,col='black',size=4)} +
  labs(x = 'Cyclers', y = 'Frequency') +   xlim(0, NA) + 
  facet_grid(Mode ~ Sample, scales = "free", axes = "all", axis.labels = "all")

ggsave(file='shuffle_2cells_KWq_2dCellQ.svg', plot=shuffle_2cells$KWq_2dCellQ$plot, width=5, height=5)
ggsave(file='shuffle_2cells_KWq_2dCellQ.pdf', plot=shuffle_2cells$KWq_2dCellQ$plot, width=5, height=5)
ggsave(file='shuffle_2cells_KWq_2dCellQ.jpg', plot=shuffle_2cells$KWq_2dCellQ$plot, width=5, height=5)

shuffle_2cells$KWq_2dCellQ$percentiles = shuffle_2cells$KWq_2dCellQ$counts %>%
  {split(., list(.$Sample, .$Mode))} %>%
  lapply(function(df) cbind(df, Rank = rank(df$count, ties.method='average')) %>%
           mutate(percentile = Rank/max(Rank)) %>%
           subset(run == 'data'))
shuffle_2cells$KWq_2dCellQ$plot
shuffle_2cells$KWq_2dCellQ$percentiles




shuffle_2cq$AGGnorm$counts %>%
  subset(run != 'data') %>%
  {split(., list(.$Sample, .$Mode))} %>%
  lapply(., function(x) {
    data.frame(Sample = x[1, 'Sample'],
               Mode = x[1, 'Mode'],
               SD = sd(x$count), 
               Mean = mean(x$count))})

shuffle_2cells$KW$counts %>%
  subset(run != 'data') %>%
  {split(., list(.$Sample, .$Mode))} %>%
  lapply(., function(x) {
    data.frame(Sample = x[1, 'Sample'],
               Mode = x[1, 'Mode'],
               SD = sd(x$count), 
               Mean = mean(x$count))})

shuffle_2cells$KWp_2dCellQ$counts %>%
  subset(run != 'data') %>%
  {split(., list(.$Sample, .$Mode))} %>%
  lapply(., function(x) {
    data.frame(Sample = x[1, 'Sample'],
               Mode = x[1, 'Mode'],
               SD = sd(x$count), 
               Mean = mean(x$count))})

shuffle_2cells$KWq_2dCellQ$counts %>%
  subset(run != 'data') %>%
  {split(., list(.$Sample, .$Mode))} %>%
  lapply(., function(x) {
    data.frame(Sample = x[1, 'Sample'],
               Mode = x[1, 'Mode'],
               SD = sd(x$count), 
               Mean = mean(x$count))})



lapply(unique(data_C2C3$sample_mode), function(S_M) {
  S = gsub('_.*', '', S_M)
  M = gsub('.*_', '', S_M)
  sm = subset(to_filter, second_pct > 0.1 & ratio > 1.5) %>%
    mutate(FDR_meta2d = p.adjust(meta2d_pvalue_AGG_norm, method='BH')) %>% 
    subset(FDR_meta2d < 0.1 & Sample==S & Mode==M, select=c(Gene, Sample, Mode))  
  
  lg = subset(reconstituting_DFs_10$reconstituted, Sample==S & Mode==M,
              select=c(Gene, Sample, Mode))
  return(data.frame(Sample = S, 
                    Mode = M, 
                    sm = nrow(sm), 
                    lg = nrow(lg), 
                    I = length(intersect(paste0(sm$Gene, sm$Sample, sm$Mode), paste0(lg$Gene, lg$Sample, lg$Mode))), 
                    U = length(union(paste0(sm$Gene, sm$Sample, sm$Mode), paste0(lg$Gene, lg$Sample, lg$Mode)))) %>%
           mutate(J = I/U))
}) %>%
  do.call(what=rbind)





####    ####    ####    ####    ####    ####    ####    ####    ####    ####   


# BOOT CELLS
################################################################################

{boot = list()
i=1

while(i<501) {
  print(i)
  set.seed(i)
  
  Bsh <- data_C2C3[[c('sample_mode', 'time2')]] %>%
    {split(., list(.$sample_mode, .$time2))} %>%
    #{.[lapply(., nrow) != 0]} %>% 
    lapply(function(DF) if(nrow(DF) >= 2) cbind(old=rownames(DF), new=sample(rownames(DF), replace=TRUE)) else cbind(old=rownames(DF), new=rownames(DF))) %>%
    do.call(what=rbind) %>%
    as.data.frame
  
  boot_data <- data_C2C3@assays$RNA$counts %>% {.[, match(Bsh$new, colnames(.))]} %>% 
    `colnames<-`(Bsh$old) %>% 
    CreateSeuratObject() %>%
    {AddMetaData(., data_C2C3@meta.data[, colnames(data_C2C3@meta.data) %!in% c('nCount_RNA', 'nFeature_RNA', 'orig.ident')])} %>%
    NormalizeData()
  
  
  
  boot_times <- AverageExpression(boot_data, group.by = c('sample_mode', 'time2')) %>% {.$RNA} %>%
    as.data.frame %>%
    {mutate(., Gene = rownames(.))} %>% 
    pivot_longer(cols=-Gene, names_to='cond') %>% 
    mutate(S_M = gsub('_.*', '', cond), 
           Time = gsub('.*_', '', cond) %>%
             as.numeric, 
           .keep='unused')
  
  boot_summary <- boot_times %>%
    mutate(rep = 1+(Time %/% 24)) %>%
    group_by(Gene, S_M, rep) %>%
    summarise(max = max(value), min = min(value)) %>%
    pivot_wider(names_from=rep, values_from=c(max, min))
  
  
  boot_times %>% 
    mutate(Gene_S_M = paste(Gene, S_M), .keep='unused') %>%
    as.data.frame %>% 
    pivot_wider(names_from = Time, values_from = value) %>% 
    as.data.frame %>%
    {`rownames<-`(., .$Gene_S_M)} %>%
    mutate(Gene_S_M = NULL) %>% 
    {.[ which(apply(., 1, max, na.rm=TRUE) != 0), order(as.numeric(colnames(.)), na.last = FALSE)]} %>%
    write.csv(file = "files_cycling/boot.csv", row.names=TRUE)
  
  boot$avg[[i]] <- meta2d(infile = "files_cycling/boot.csv", filestyle = "csv", timepoints = "line1",
                                      outdir=getwd(), parallelize = TRUE, outputFile = FALSE, nCores = detectCores())$meta %>%
    mutate(Gene = gsub(' .*', "", CycID), 
           S_M = gsub('.* ', "", CycID), .keep='unused') %>%
    merge(unique(subset(boot_times, select=c(Gene, S_M))), all=TRUE) %>%
    merge(boot_summary) %>%
    {.[, c('Gene', 'S_M', 'ARS_pvalue', 'JTK_pvalue', 'LS_pvalue', 'meta2d_pvalue', 'meta2d_phase', 'max_1', 'max_2', 'min_1', 'min_2')]}
  
  rm(boot_times, boot_summary)
  
  boot_times_AGG <- AggregateExpression(boot_data, group.by = 'sample_mode_group2', return.seurat = TRUE) %>% 
    NormalizeData %>%
    AverageExpression %>%
    {.$RNA} %>%
    as.data.frame %>%
    {mutate(., Gene = rownames(.))} %>% 
    pivot_longer(cols=-Gene, names_to='cond') %>% 
    mutate(S_M = gsub('D-.*', 'D', cond), 
           Time = gsub('.*-', '', cond) %>%
             gsub('[a-zA-Z]', '', .) %>% 
             as.numeric, 
           .keep='unused')
  
  boot_summary_AGG <- boot_times_AGG %>%
    mutate(rep = 1+(Time %/% 24)) %>%
    group_by(Gene, S_M, rep) %>%
    summarise(max = max(value), min = min(value)) %>%
    pivot_wider(names_from=rep, values_from=c(max, min))
  
  
  boot_times_AGG %>% 
    mutate(Gene_S_M = paste(Gene, S_M), .keep='unused') %>%
    as.data.frame %>% 
    pivot_wider(names_from = Time, values_from = value) %>% 
    as.data.frame %>%
    {`rownames<-`(., .$Gene_S_M)} %>%
    mutate(Gene_S_M = NULL) %>% 
    {.[ which(apply(., 1, max, na.rm=TRUE) != 0), order(as.numeric(colnames(.)), na.last = FALSE)]} %>%
    write.csv(file = "files_cycling/boot.csv", row.names=TRUE)
  
  boot$AGG[[i]] <- meta2d(infile = "files_cycling/boot.csv", filestyle = "csv", timepoints = "line1",
                      outdir=getwd(), parallelize = TRUE, outputFile = FALSE, nCores = detectCores())$meta %>%
    mutate(Gene = gsub(' .*', "", CycID), 
           S_M = gsub('.* ', "", CycID), .keep='unused') %>%
    merge(unique(subset(boot_times_AGG, select=c(Gene, S_M))), all=TRUE) %>%
    merge(boot_summary_AGG) %>%
    {.[, c('Gene', 'S_M', 'ARS_pvalue', 'JTK_pvalue', 'LS_pvalue', 'meta2d_pvalue', 'meta2d_phase', 'max_1', 'max_2', 'min_1', 'min_2')]}
  
  rm(boot_times_AGG, boot_data, boot_summary_AGG)
  
  i=i+1
  
}
rm(i)
}



####    ####    ####    ####    ####    ####    ####    ####    ####    ####   



# NO RATIO FILTERING
####################

shuffle_2cells$noRatio$retention <- shuffle_2cells$retention
shuffle_2cells$noRatio$prefilt <- shuffle_2cells$prefilt

# FDR q-values
shuffle_2cells$noRatio$AGGQ$counts <- lapply(shuffle_2cells$noRatio$retention, function(df) {
  merge(df, subset(shuffle_2cells$noRatio$prefilt, select = -c(ratio, ratio_AGG))) %>%
    mutate(FDR_meta2d_AGG = p.adjust(p_AGG, method='BH')) %>% 
    subset(FDR_meta2d_AGG < 0.1) %>% 
    {table(.$Sample, .$Mode)} %>% 
    as.data.frame %>% 
    `colnames<-`(c('Sample', 'Mode', 'count'))
}) %>% 
  {lapply(seq_along(.), function(i) cbind(.[[i]], run=i))} %>%
  c(list(data = cbind(subset(shuffle_2cells$noRatio$prefilt) %>%
                        mutate(FDR_meta2d_AGG = p.adjust(meta2d_pvalue_AGG_norm, method='BH')) %>%
                        subset(FDR_meta2d_AGG < 0.1) %>% 
                        {table(.$Sample, .$Mode)} %>% 
                        as.data.frame %>% 
                        `colnames<-`(c('Sample', 'Mode', 'count')), run='data'))) %>% 
  do.call(what='rbind') %>%
  mutate(Sample = factor(Sample, levels=c('C3', 'C2')), 
         Mode = factor(Mode, levels=c('LD', 'DD'))) %>%
  as.data.frame

shuffle_2cells$noRatio$AGGQ$plot = shuffle_2cells$noRatio$AGGQ$counts %>% 
  {merge(., group_by(subset(., run != 'data'), Sample, Mode) %>% summarise(Mean = mean(count)))} %>% 
  {ggplot(data=.) +
      geom_histogram(data=subset(., run != 'data'), aes(x = count), binwidth=1, fill="gray30",color="black") +
      geom_vline(data=subset(., run == 'data'), aes(xintercept = count)) +
      geom_vline(data=., aes(xintercept = Mean), linetype = "dashed", linewidth = 1) +
      geom_shadowtext(data=subset(., run==1), aes(label=round(Mean),y=Inf,x=Mean), vjust=4,hjust=1.25,col='white',size=4) + 
      geom_text(data=subset(., run == 'data'), aes(label=round(count),y=Inf,x=count),vjust=2,hjust=1.25,col='black',size=4)} +
  labs(x = 'Cyclers', y = 'Frequency') +   xlim(0, NA) + 
  facet_grid(Mode ~ Sample, scales = "free", axes = "all", axis.labels = "all")


ggsave(file='shuffle_2cells_noRatio_Q1.5.svg', plot=shuffle_2cells$noRatio$AGGQ$plot, width=5, height=5)
ggsave(file='shuffle_2cells_noRatio_Q1.5.pdf', plot=shuffle_2cells$noRatio$AGGQ$plot, width=5, height=5)
ggsave(file='shuffle_2cells_noRatio_Q1.5.jpg', plot=shuffle_2cells$noRatio$AGGQ$plot, width=5, height=5)



shuffle_2cells$noRatio$AGGQ$percentiles = shuffle_2cells$noRatio$AGGQ$counts %>%
  {split(., list(.$Sample, .$Mode))} %>%
  lapply(function(df) cbind(df, Rank = rank(df$count, ties.method='average')) %>%
           mutate(percentile = Rank/max(Rank)) %>%
           subset(run == 'data'))

# KW

shuffle_2cells$noRatio$KW$counts <- lapply(shuffle_2cells$noRatio$retention, function(df) {
  merge(df, to_filter[,c('Gene', 'Sample', 'Mode', 'second_pct', 'meta2d_pvalue', 'meta2d_pvalue_AGG_norm', 'KW_scores_norm')]) %>%
    subset(second_pct > 0.1 & p_cell < 0.05 & KW_scores_norm < 0.05) %>% 
    {table(.$Sample, .$Mode)} %>% 
    as.data.frame %>% 
    `colnames<-`(c('Sample', 'Mode', 'count'))
}) %>%
  {lapply(seq_along(.), function(i) cbind(.[[i]], run=i))} %>%
  c(list(data = cbind(subset(to_filter, second_pct > 0.1 & meta2d_pvalue < 0.05 & KW_scores_norm < 0.05) %>%
                        {table(.$Sample, .$Mode)} %>% 
                        as.data.frame %>% 
                        `colnames<-`(c('Sample', 'Mode', 'count')), run='data'))) %>% 
  do.call(what='rbind') %>%
  mutate(Sample = factor(Sample, levels=c('C3', 'C2')), 
         Mode = factor(Mode, levels=c('LD', 'DD'))) %>%
  as.data.frame

shuffle_2cells$noRatio$KW$plot = shuffle_2cells$noRatio$KW$counts %>%
  {merge(., group_by(subset(., run != 'data'), Sample, Mode) %>% summarise(Mean = mean(count)))} %>% 
  {ggplot(data=.) +
      geom_histogram(data=subset(., run != 'data'), aes(x = count), binwidth=1, fill="gray30",color="black") +
      geom_vline(data=subset(., run == 'data'), aes(xintercept = count)) +
      geom_vline(data=., aes(xintercept = Mean), linetype = "dashed", linewidth = 1) +
      geom_shadowtext(data=subset(., run==1), aes(label=round(Mean),y=Inf,x=Mean), vjust=4,hjust=1.25,col='white',size=4) + 
      geom_text(data=subset(., run == 'data'), aes(label=round(count),y=Inf,x=count),vjust=2,hjust=1.25,col='black',size=4)} +
  labs(x = 'Cyclers', y = 'Frequency') +   xlim(0, NA) + 
  facet_grid(Mode ~ Sample, scales = "free", axes = "all", axis.labels = "all")

ggsave(file='shuffle_2cells_noRatio_KW.svg', plot=shuffle_2cells$noRatio$KW$plot, width=5, height=5)
ggsave(file='shuffle_2cells_noRatio_KW.pdf', plot=shuffle_2cells$noRatio$KW$plot, width=5, height=5)
ggsave(file='shuffle_2cells_noRatio_KW.jpg', plot=shuffle_2cells$noRatio$KW$plot, width=5, height=5)

shuffle_2cells$noRatio$KW$percentiles = shuffle_2cells$noRatio$KW$counts %>%
  {split(., list(.$Sample, .$Mode))} %>%
  lapply(function(df) cbind(df, Rank = rank(df$count, ties.method='average')) %>%
           mutate(percentile = Rank/max(Rank)) %>%
           subset(run == 'data'))



# KW p w/ Q (cell)

shuffle_2cells$noRatio$KWp_2dCellQ$counts <- lapply(shuffle_2cells$noRatio$retention, function(df) {
  merge(df, to_filter[,c('Gene', 'Sample', 'Mode', 'second_pct', 'meta2d_pvalue', 'meta2d_pvalue_AGG_norm', 'KW_scores_norm')]) %>%
    subset(second_pct > 0.1 & KW_scores_norm < 0.05) %>% 
    mutate(FDR_meta2d = p.adjust(p_cell, method='BH')) %>% 
    subset(FDR_meta2d < 0.1) %>% 
    {if (!nrow(.)) data.frame(c('C2', 'C3', 'C2', 'C3'), c('DD', 'DD', 'LD', 'LD'), rep(0, 4))
      else table(.$Sample, .$Mode)} %>% 
    as.data.frame %>% 
    `colnames<-`(c('Sample', 'Mode', 'count'))
}) %>%
  {lapply(seq_along(.), function(i) cbind(.[[i]], run=i))} %>%
  c(list(data = cbind(subset(to_filter, second_pct > 0.1 & KW_scores_norm < 0.05) %>% 
                        mutate(FDR_meta2d = p.adjust(meta2d_pvalue, method='BH')) %>%
                        subset(FDR_meta2d < 0.1) %>% 
                        {table(.$Sample, .$Mode)} %>% 
                        as.data.frame %>% 
                        `colnames<-`(c('Sample', 'Mode', 'count')), run='data'))) %>% 
  do.call(what='rbind') %>%
  mutate(Sample = factor(Sample, levels=c('C3', 'C2')), 
         Mode = factor(Mode, levels=c('LD', 'DD'))) %>%
  as.data.frame

shuffle_2cells$noRatio$KWp_2dCellQ$plot = shuffle_2cells$noRatio$KWp_2dCellQ$counts %>%
  {merge(., group_by(subset(., run != 'data'), Sample, Mode) %>% summarise(Mean = mean(count)))} %>% 
  {ggplot(data=.) +
      geom_histogram(data=subset(., run != 'data'), aes(x = count), binwidth=1, fill="gray30",color="black") +
      geom_vline(data=subset(., run == 'data'), aes(xintercept = count)) +
      geom_vline(data=., aes(xintercept = Mean), linetype = "dashed", linewidth = 1) +
      geom_shadowtext(data=subset(., run==1), aes(label=round(Mean),y=Inf,x=Mean), vjust=4,hjust=1.25,col='white',size=4) + 
      geom_text(data=subset(., run == 'data'), aes(label=round(count),y=Inf,x=count),vjust=2,hjust=1.25,col='black',size=4)} +
  labs(x = 'Cyclers', y = 'Frequency') +   xlim(0, NA) + 
  facet_grid(Mode ~ Sample, scales = "free", axes = "all", axis.labels = "all")

ggsave(file='shuffle_2cells_noRatio_KWp_2dCellQ.svg', plot=shuffle_2cells$noRatio$KWp_2dCellQ$plot, width=5, height=5)
ggsave(file='shuffle_2cells_noRatio_KWp_2dCellQ.pdf', plot=shuffle_2cells$noRatio$KWp_2dCellQ$plot, width=5, height=5)
ggsave(file='shuffle_2cells_noRatio_KWp_2dCellQ.jpg', plot=shuffle_2cells$noRatio$KWp_2dCellQ$plot, width=5, height=5)


shuffle_2cells$noRatio$KWp_2dCellQ$percentiles = shuffle_2cells$noRatio$KWp_2dCellQ$counts %>%
  {split(., list(.$Sample, .$Mode))} %>%
  lapply(function(df) cbind(df, Rank = rank(df$count, ties.method='average')) %>%
           mutate(percentile = Rank/max(Rank)) %>%
           subset(run == 'data'))

# KW Q w/ Q (cell)

shuffle_2cells$noRatio$KWq_2dCellQ$counts <- lapply(shuffle_2cells$noRatio$retention, function(df) {
  merge(df, to_filter[,c('Gene', 'Sample', 'Mode', 'second_pct', 'meta2d_pvalue', 'meta2d_pvalue_AGG_norm', 'KW_scores_norm')]) %>%
    subset(second_pct > 0.1) %>% 
    mutate(FDR_KW = p.adjust(KW_scores_norm, method='BH')) %>% 
    subset(FDR_KW < 0.1) %>%     
    mutate(FDR_meta2d = p.adjust(p_cell, method='BH')) %>% 
    subset(FDR_meta2d < 0.1) %>% 
    {if (!nrow(.)) data.frame(c('C2', 'C3', 'C2', 'C3'), c('DD', 'DD', 'LD', 'LD'), rep(0, 4))
      else table(.$Sample, .$Mode)} %>% 
    as.data.frame %>% 
    `colnames<-`(c('Sample', 'Mode', 'count'))
}) %>%
  {lapply(seq_along(.), function(i) cbind(.[[i]], run=i))} %>%
  c(list(data = cbind(subset(to_filter, second_pct > 0.1) %>% 
                        mutate(FDR_KW = p.adjust(KW_scores_norm, method='BH')) %>% 
                        subset(FDR_KW < 0.1) %>% 
                        mutate(FDR_meta2d = p.adjust(meta2d_pvalue, method='BH')) %>%
                        subset(FDR_meta2d < 0.1) %>% 
                        {table(.$Sample, .$Mode)} %>% 
                        as.data.frame %>% 
                        `colnames<-`(c('Sample', 'Mode', 'count')), run='data'))) %>% 
  do.call(what='rbind') %>%
  mutate(Sample = factor(Sample, levels=c('C3', 'C2')), 
         Mode = factor(Mode, levels=c('LD', 'DD'))) %>%
  as.data.frame

shuffle_2cells$noRatio$KWq_2dCellQ$plot = shuffle_2cells$noRatio$KWq_2dCellQ$counts %>% 
  {merge(., group_by(subset(., run != 'data'), Sample, Mode) %>% summarise(Mean = mean(count)))} %>% 
  {ggplot(data=.) +
      geom_histogram(data=subset(., run != 'data'), aes(x = count), binwidth=1, fill="gray30",color="black") +
      geom_vline(data=subset(., run == 'data'), aes(xintercept = count)) +
      geom_vline(data=., aes(xintercept = Mean), linetype = "dashed", linewidth = 1) +
      geom_shadowtext(data=subset(., run==1), aes(label=round(Mean),y=Inf,x=Mean), vjust=4,hjust=1.25,col='white',size=4) + 
      geom_text(data=subset(., run == 'data'), aes(label=round(count),y=Inf,x=count),vjust=2,hjust=1.25,col='black',size=4)} +
  labs(x = 'Cyclers', y = 'Frequency') +   xlim(0, NA) + 
  facet_grid(Mode ~ Sample, scales = "free", axes = "all", axis.labels = "all")

ggsave(file='shuffle_2cells_noRatio_KWq_2dCellQ.svg', plot=shuffle_2cells$noRatio$KWq_2dCellQ$plot, width=5, height=5)
ggsave(file='shuffle_2cells_noRatio_KWq_2dCellQ.pdf', plot=shuffle_2cells$noRatio$KWq_2dCellQ$plot, width=5, height=5)
ggsave(file='shuffle_2cells_noRatio_KWq_2dCellQ.jpg', plot=shuffle_2cells$noRatio$KWq_2dCellQ$plot, width=5, height=5)

shuffle_2cells$noRatio$KWq_2dCellQ$percentiles = shuffle_2cells$noRatio$KWq_2dCellQ$counts %>%
  {split(., list(.$Sample, .$Mode))} %>%
  lapply(function(df) cbind(df, Rank = rank(df$count, ties.method='average')) %>%
           mutate(percentile = Rank/max(Rank)) %>%
           subset(run == 'data'))
shuffle_2cells$noRatio$KWq_2dCellQ$plot
shuffle_2cells$noRatio$KWq_2dCellQ$percentiles

####    ####    ####


# 11. HIERARCHICAL CLUSTERING
################################################################################

# Gene, ENSEMBL_ID, Type, Mode, Cycler, second_pct, ratio_avg_norm, ratio_bulk_norm, meta2d_pvalue_avg_norm, meta2d_pvalue_bulk_norm


{clusterings <- list()
N = 10

clusterings[[as.character(N)]] <- reconstituting_DFs_10$reconstituted %>% 
  subset(select = c(Gene, Sample, Mode)) %>%
  merge(avg_data_times2_01) %>%
  mutate(sample_mode_time = NULL) %>%
  pivot_wider(names_from='timepoint', values_from='avg') %>% 
  {split(., list(.$Sample, .$Mode))} %>%
  lapply(function(df) {
    cbind(df, clst = hclust(dist(df[,colnames(df) %in% seq(3, 47, 4)])) %>%
      cutree(N))
  }) %>%
  do.call(what='rbind')
clusterings[[as.character(N)]] %>% 
  {split(., list(.$clst, .$Sample, .$Mode))} %>% 
  lapply(function(df) plotter_geneCombine(df$Gene, data=data_C2C3_AGG_01_1d, COUNTS=TRUE, NO_X_TEXT=TRUE, Sample=df[1,'Sample'], Mode=df[1,'Mode'], PRINT_GENES=TRUE)) %>%
  plot_grid(plotlist=., ncol=N)
rm(N)
}
clusterings[['10']] %>%
  {split(., list(.$clst, .$Sample, .$Mode))} %>% 
  lapply(function(df) print(paste('Sample Mode:', df[1,'Sample'], df[1,'Mode'], '\n', 'Genes:', df$Gene)))
  
####    ####    ####    ####    ####    ####    ####    ####    ####    ####    


# 12. PRINTING SUPPLIMENTARY TABLES
################################################################################

# Gene, ENSEMBL_ID, Type, Mode, Cycler, second_pct, ratio_avg_norm, ratio_bulk_norm, meta2d_pvalue_avg_norm, meta2d_pvalue_bulk_norm

subset(to_filter, second_pct>0.1, select=c(Gene, Sample, Mode, ratio, ratio_AGG, 
                                           meta2d_pvalue, meta2d_pvalue_AGG_norm, meta2d_phase, meta2d_phase_AGG_norm)) %>%
  merge(passQC) %>% 
  mutate(PASS = ((ratio > 1.5 & meta2d_pvalue < 0.05 & ratio_AGG > 1.5 & meta2d_pvalue_AGG_norm < 0.05) | 
                                       (meta2d_pvalue_AGG_norm < 0.05 & ratio_AGG > 2))) %>%
  .[,c(1,10,2,3,11,4:9)] %>% 
  `names<-`(c('Gene', 'ENSEMBL_ID', 'Cell_type', 'Mode', 'PASS',
              'ratio_avg_norm', 'ratio_bulk_norm', 'meta2d_pvalue_avg_norm', 'meta2d_pvalue_bulk_norm',
              'meta2d_phase_norm', 'meta2d_phase_AGG_norm')) %>% 
  write.table(file='Supplementary1.txt', quote=FALSE, row.names=FALSE)



####    ####    ####    ####    ####    ####    ####    ####    ####    ####    


# 13. Variation
################################################################################
avg_data_down_times2_AGG <- SetIdent(data_C2C3, value="sample_mode") %>%
  subset(downsample = 900) %>%
  AverageExpression(group.by = c('sample_mode', 'group2'), layer="data")$RNA %>% as.data.frame %>%
  {mutate(., Gene = rownames(.))} %>%
  pivot_longer(cols=-'Gene', names_to='sample_mode_time', values_to='avg') %>%
  mutate(Mode = gsub('_.*', '', sample_mode_time) %>% gsub('.*-', '', .), 
         Sample = gsub('-.*', '', sample_mode_time), 
         timepoint = gsub('.*_', '', sample_mode_time) %>% gsub('.*T', '', .),
         sample_mode_time = gsub('_', '-', sample_mode_time))


var_df <- avg_data_times2_AGG %>%
  #avg_data_down_times2_AGG %>%
  mutate(sample_mode_time = NULL) %>%
  pivot_wider(names_from = 'timepoint', values_from = 'avg') %>%
  merge(subset(to_filter, second_pct > 0.1, select=c(Gene, Sample, Mode))) %>%
  {cbind(., SD = apply(.[, colnames(.) %!in% c('Gene', 'Mode', 'Sample')], 1, sd),
         Mean = apply(.[, colnames(.) %!in% c('Gene', 'Mode', 'Sample')], 1, mean))} %>%
  mutate(CV = SD/Mean) %>%
  mutate(Mode = factor(Mode, levels=c('LD', 'DD')),
         Sample = factor(Sample, levels=c('C3', 'C2')), 
         Sample_Mode = factor(paste0(Sample, Mode), levels=c('C3LD', 'C3DD', 'C2LD', 'C2DD')))

ggplot(var_df, aes(CV, group=Sample_Mode, fill=Sample_Mode)) +
  geom_density(alpha=0.2)

split(var_df, var_df$Sample_Mode) %>%
  {`names<-`(lapply(., function(x) mean(x$CV)), names(.))}

ggplot(var_df, aes(SD, group=Sample_Mode, fill=Sample_Mode)) +
  geom_density(alpha=0.2) +
  scale_x_continuous(trans='log10')

split(var_df, var_df$Sample_Mode) %>%
  {`names<-`(lapply(., function(x) mean(x$SD)), names(.))}

table(data_C2C3$sample_mode)

split(var_df, var_df$Sample_Mode) %>% 
  {ks.test(.$C3LD$CV, .$C2LD$CV)}

split(var_df, var_df$Sample_Mode) %>% 
  {ks.test(.$C3LD$CV, .$C3DD$CV)}
split(var_df, var_df$Sample_Mode) %>% 
  {ks.test(.$C2LD$CV, .$C2DD$CV)}

split(var_df, var_df$Sample_Mode) %>% 
  {ks.test(.$C3LD$CV, .$C2LD$CV)}

split(var_df, var_df$Sample_Mode) %>% 
  {wilcox.test(.$C3LD$CV, .$C2LD$CV, paired = FALSE)}




rm(data_down, avg_data_down_times2_AGG)

####    ####    ####    ####    ####    ####    ####    ####    ####    ####    


# 14. DOTPLOTS
################################################################################

{setwd(WD)
dir.create('DOTPLOTS', showWarnings = FALSE)
for (S_M in levels(data_C2C3$sample_mode)) {
  print(S_M)
  setwd(paste(WD, 'DOTPLOTS', sep='/'))
  dir.create(S_M, showWarnings = FALSE)
  setwd(S_M)
  
  DotPlot(subset(data_C2C3, sample_mode==S_M), features=circ_genes, group.by = 'time1') %>%
    ggsave(file=paste0('CCG_', S_M, '.svg'), plot=., width=6, height=6)
  
  #DotPlot(subset(data_C2C3, sample_mode==S_M), features=setdiff(circ_genes, 'Pdp1'), group.by = 'time1') %>%
  #  ggsave(file=paste0('CCG-Pdp1_', S_M, '.svg'), plot=., width=6, height=6)
  
  #DotPlot(subset(data_C2C3, sample_mode==S_M), features=c(circ_genes, 'Pdfr'), group.by = 'time1') %>%
  #  ggsave(file=paste0('CCG+Pdfr_', S_M, '.svg'), plot=., width=6, height=6)
  
  #DotPlot(subset(data_C2C3, sample_mode==S_M), features=c(setdiff(circ_genes, 'Pdp1'), 'Pdfr'), group.by = 'time1') %>%
  #  ggsave(file=paste0('CCG+Pdp1-Pdp1_', S_M, '.svg'), plot=., width=6, height=6)
}
setwd(WD)
}

{setwd(WD)
dir.create('DOTPLOTS', showWarnings = FALSE)
setwd('DOTPLOTS')
lapply(levels(data_C2C3$sample_mode), function(S_M) DotPlot(subset(data_C2C3, sample_mode==S_M), features=c(circ_genes, 'cyc'), group.by = 'sample_mode', scale.min=0, sac)) %>%
  ggarrange(nrow=4,ncol=1, plotlist = ., common.legend = TRUE, legend = 'right') %>% 
  ggsave(file=paste0('CCG+cyc_flat', '.svg'), plot=., width=6, height=6)

lapply(levels(data_C2C3$sample_mode), function(S_M) DotPlot(subset(data_C2C3, sample_mode==S_M), features=c(setdiff(circ_genes, 'Pdp1'), 'cyc'), group.by = 'sample_mode')) %>%
  ggarrange(nrow=4,ncol=1, plotlist = ., common.legend = TRUE, legend = 'right') %>% 
  ggsave(file=paste0('CCG+cyc-Pdp1_flat', '.svg'), plot=., width=6, height=6)

lapply(levels(data_C2C3$sample_mode), function(S_M) DotPlot(subset(data_C2C3, sample_mode==S_M), features=c(circ_genes, 'Pdfr', 'cyc'), group.by = 'sample_mode')) %>%
  ggarrange(nrow=4,ncol=1, plotlist = ., common.legend = TRUE, legend = 'right') %>% 
  ggsave(file=paste0('CCG+cyc+Pdfr_flat_lg', '.svg'), plot=., width=6, height=6)

lapply(levels(data_C2C3$sample_mode), function(S_M) DotPlot(subset(data_C2C3, sample_mode==S_M), features=c(setdiff(circ_genes, 'Pdp1'), 'Pdfr', 'cyc'), group.by = 'sample_mode')) %>%
  ggarrange(nrow=4,ncol=1, plotlist = ., common.legend = TRUE, legend = 'right') %>% 
  ggsave(file=paste0('CCG+cyc-Pdp1+Pdfr_flat', '.svg'), plot=., width=6, height=6)

setwd(WD)
}

{setwd(WD)
  dir.create('DOTPLOTS', showWarnings = FALSE)
  for (Sa in levels(data_C2C3$sample)) {
    print(Sa)
    setwd(paste(WD, 'DOTPLOTS', sep='/'))
    dir.create(Sa, showWarnings = FALSE)
    setwd(Sa)
    
    DotPlot(subset(data_C2C3, sample==Sa), features=circ_genes, group.by = 'time1') %>%
      ggsave(file=paste0('sample_CCG_', Sa, '.svg'), plot=., width=6, height=6)
    
    DotPlot(subset(data_C2C3, sample==Sa), features=setdiff(circ_genes, 'Pdp1'), group.by = 'time1') %>%
      ggsave(file=paste0('sample_CCG-Pdp1_', Sa, '.svg'), plot=., width=6, height=6)
    
    DotPlot(subset(data_C2C3, sample==Sa), features=c(circ_genes, 'Pdfr'), group.by = 'time1') %>%
      ggsave(file=paste0('sample_CCG+Pdfr_', Sa, '.svg'), plot=., width=6, height=6)
    
    DotPlot(subset(data_C2C3, sample==Sa), features=c(setdiff(circ_genes, 'Pdp1'), 'Pdfr'), group.by = 'time1') %>%
      ggsave(file=paste0('sample_CCG+Pdp1-Pdp1_', Sa, '.svg'), plot=., width=6, height=6)
  }
  setwd(WD)
}


{setwd(WD)
  dir.create('DOTPLOTS', showWarnings = FALSE)
  setwd('DOTPLOTS')
  lapply(levels(data_C2C3$sample), function(Sa) DotPlot(subset(data_C2C3, sample==Sa), features=c(circ_genes, 'cyc'), group.by = 'sample')) %>%
    ggarrange(nrow=2,ncol=1, plotlist = ., common.legend = TRUE, legend = 'right') %>% 
    {ggsave(file=paste0('sample_CCG+cyc_flat', '.pdf'), plot=., width=5, height=3)
      ggsave(file=paste0('sample_CCG+cyc_flat', '.svg'), plot=., width=5, height=3)}
  
  lapply(levels(data_C2C3$sample), function(Sa) DotPlot(subset(data_C2C3, sample==Sa), features=c(setdiff(circ_genes, 'Pdp1'), 'cyc'), group.by = 'sample')) %>%
    ggarrange(nrow=2,ncol=1, plotlist = ., common.legend = TRUE, legend = 'right') %>% 
    {ggsave(file=paste0('sample_CCG+cyc-Pdp1_flat', '.pdf'), plot=., width=5, height=3)
      ggsave(file=paste0('sample_CCG+cyc-Pdp1_flat', '.svg'), plot=., width=5, height=3)}
  
  lapply(levels(data_C2C3$sample), function(Sa) DotPlot(subset(data_C2C3, sample==Sa), features=c(circ_genes, 'Pdfr', 'cyc'), group.by = 'sample')) %>%
    ggarrange(nrow=2,ncol=1, plotlist = ., common.legend = TRUE, legend = 'right') %>% 
    {ggsave(file=paste0('sample_CCG+cyc+Pdfr_flat', '.pdf'), plot=., width=5, height=3)
      ggsave(file=paste0('sample_CCG+cyc+Pdfr_flat', '.svg'), plot=., width=5, height=3)}
  
  lapply(levels(data_C2C3$sample), function(Sa) DotPlot(subset(data_C2C3, sample==Sa), features=c(setdiff(circ_genes, 'Pdp1'), 'Pdfr', 'cyc'), group.by = 'sample')) %>%
    ggarrange(nrow=2,ncol=1, plotlist = ., common.legend = TRUE, legend = 'right') %>% 
    {ggsave(file=paste0('sample_CCG+cyc-Pdp1+Pdfr_flat', '.pdf'), plot=., width=5, height=3)
      ggsave(file=paste0('sample_CCG+cyc-Pdp1+Pdfr_flat', '.svg'), plot=., width=5, height=3)}
  
  setwd(WD)
}

{setwd(WD)
  dir.create('DOTPLOTS', showWarnings = FALSE)
  setwd('DOTPLOTS')
  lapply(levels(data_C2C3$sample), function(Sa) DotPlot(subset(data_C2C3, sample==Sa), features=c(circ_genes, 'cyc'), group.by = 'sample') + 
    theme(axis.text = element_text(size = 6), text = element_text(size = 6))) %>%
  ggarrange(nrow=2,ncol=1, plotlist = ., common.legend = TRUE, legend = 'right') %>% 
  {ggsave(file=paste0('square_sample_CCG+cyc_flat', '.pdf'), plot=., width=3, height=3)
    ggsave(file=paste0('square_sample_CCG+cyc_flat', '.svg'), plot=., width=3, height=3)}

  lapply(levels(data_C2C3$sample), function(Sa) DotPlot(subset(data_C2C3, sample==Sa), features=c(setdiff(circ_genes, 'Pdp1'), 'cyc'), group.by = 'sample') + 
           theme(axis.text = element_text(size = 6), text = element_text(size = 6))) %>%
    ggarrange(nrow=2,ncol=1, plotlist = ., common.legend = TRUE, legend = 'right') %>% 
    {ggsave(file=paste0('square_sample_CCG+cyc-Pdp1_flat', '.pdf'), plot=., width=3, height=3)
      ggsave(file=paste0('square_sample_CCG+cyc-Pdp1_flat', '.svg'), plot=., width=3, height=3)}
  
  lapply(levels(data_C2C3$sample), function(Sa) DotPlot(subset(data_C2C3, sample==Sa), features=c(circ_genes, 'Pdfr', 'cyc'), group.by = 'sample') + 
           theme(axis.text = element_text(size = 6), text = element_text(size = 6))) %>%
    ggarrange(nrow=2,ncol=1, plotlist = ., common.legend = TRUE, legend = 'right') %>% 
    {ggsave(file=paste0('square_sample_CCG+cyc+Pdfr_flat', '.pdf'), plot=., width=3, height=3)
      ggsave(file=paste0('square_sample_CCG+cyc+Pdfr_flat', '.svg'), plot=., width=3, height=3)}
  
  lapply(levels(data_C2C3$sample), function(Sa) DotPlot(subset(data_C2C3, sample==Sa), features=c(setdiff(circ_genes, 'Pdp1'), 'Pdfr', 'cyc'), group.by = 'sample') + 
           theme(axis.text = element_text(size = 6), text = element_text(size = 6))) %>%
    ggarrange(nrow=2,ncol=1, plotlist = ., common.legend = TRUE, legend = 'right') %>% 
    {ggsave(file=paste0('square_sample_CCG+cyc-Pdp1+Pdfr_flat', '.pdf'), plot=., width=3, height=3)
      ggsave(file=paste0('square_sample_CCG+cyc-Pdp1+Pdfr_flat', '.svg'), plot=., width=3, height=3)}
  setwd(WD)
}

setwd(paste0(WD, '/DOTPLOTS'))
lapply(levels(data_C2C3$sample), function(Sa) DotPlot(subset(data_C2C3, sample==Sa), features=c(circ_genes, 'Pdfr', 'cyc'), group.by = 'sample', scale.min=0, scale.max=100, scale.by='radius') + 
         theme(axis.text = element_text(size = 6), text = element_text(size = 6))) %>%
  ggarrange(nrow=2,ncol=1, plotlist = ., common.legend = TRUE, legend = 'right') %>% 
  {lapply(c('.pdf', '.svg'), function(ext) ggsave(file=paste0('square_radius_CCG+cyc+Pdfr_flat', ext), plot=., width=3, height=3))}

lapply(levels(data_C2C3$sample), function(Sa) DotPlot(subset(data_C2C3, sample==Sa), features=c(circ_genes, 'Pdfr', 'cyc'), group.by = 'sample', scale.min=0, scale.max=100, scale.by='size') + 
         theme(axis.text = element_text(size = 6), text = element_text(size = 6))) %>%
  ggarrange(nrow=2,ncol=1, plotlist = ., common.legend = TRUE, legend = 'right') %>% 
  {lapply(c('.pdf', '.svg'), function(ext) ggsave(file=paste0('square_size_CCG+cyc+Pdfr_flat', ext), plot=., width=3, height=3))}

lapply(levels(data_C2C3$sample), function(Sa) DotPlot(subset(data_C2C3, sample==Sa), features=c(circ_genes, 'Pdfr', 'cyc'), group.by = 'sample', scale.min=0, scale.max=100, scale.by='radius') + 
         theme(axis.text = element_text(size = 8), text = element_text(size = 8))) %>%  ggarrange(nrow=2,ncol=1, plotlist = ., common.legend = TRUE, legend = 'right') %>% 
  {lapply(c('.pdf', '.svg'), function(ext) ggsave(file=paste0('radius_CCG+cyc+Pdfr_flat', ext), plot=., width=5, height=3))}

lapply(levels(data_C2C3$sample), function(Sa) DotPlot(subset(data_C2C3, sample==Sa), features=c(circ_genes, 'Pdfr', 'cyc'), group.by = 'sample', scale.min=0, scale.max=100, scale.by='size') +
         theme(axis.text = element_text(size = 8), text = element_text(size = 8))) %>%  ggarrange(nrow=2,ncol=1, plotlist = ., common.legend = TRUE, legend = 'right') %>% 
  {lapply(c('.pdf', '.svg'), function(ext) ggsave(file=paste0('size_CCG+cyc+Pdfr_flat', ext), plot=., width=5, height=3))}
setwd(WD)


setwd(paste0(WD, '/DOTPLOTS'))
lapply(levels(data_C2C3$sample_mode), function(S_M) {
  DotPlot(subset(data_C2C3, sample_mode==S_M), features=circ_genes, group.by = 'time1', scale.min=0, scale.max=100, scale.by='size', scale=FALSE) %>%
    ggsave(file=paste0('CCG_size_', S_M, '.svg'), plot=., width=6, height=6)
  DotPlot(subset(data_C2C3, sample_mode==S_M), features=circ_genes, group.by = 'time1', scale.min=0, scale.max=100, scale.by='radius', scale=FALSE) %>%
    ggsave(file=paste0('CCG_radius_', S_M, '.svg'), plot=., width=6, height=6)
  
  DotPlot(subset(data_C2C3, sample_mode==S_M), features=setdiff(circ_genes, 'Pdp1'), group.by = 'time1', scale.by='size', scale=FALSE) %>%
    ggsave(file=paste0('CCG_size_noPdp1', S_M, '.svg'), plot=., width=6, height=6)
  DotPlot(subset(data_C2C3, sample_mode==S_M), features=setdiff(circ_genes, 'Pdp1'), group.by = 'time1', scale.by='radius', scale=FALSE) %>%
    ggsave(file=paste0('CCG_radius_noPdp1_', S_M, '.svg'), plot=., width=6, height=6)
})


lapply(levels(data_C2C3$sample_mode), function(S_M) DotPlot(subset(data_C2C3, sample_mode==S_M), features=c('hdc', 'ort', 'HisCl1'), group.by = 'sample_mode')) %>%
  ggarrange(nrow=4,ncol=1, plotlist = ., common.legend = TRUE, legend = 'right')


subset(reconstituting_DFs_10$reconstituted, Gene %in% c('hdc', 'ort', 'HisCl1', 'Ace', 'ChAT', rownames(data_C2C3) %>% {.[grep('AChR', .)]} %>% sort))
avg_plotter_concat_mode(c('hdc', 'ort', 'HisCl1', 'Ace', 'ChAT', rownames(data_C2C3) %>% {.[grep('AChR', .)]} %>% sort), data=data_C2C3_AGG, N_days=2) %>% ggsave(file='GO_C3DD_receptors.svg', plot=., width=6, height=6)

####    ####    ####    ####    ####    ####    ####    ####    ####    ####    


# 15. VENN DIAGRAMS
################################################################################

dir.create('./Venn', showWarnings = FALSE)

P = avg_plotter_concat('tim', data=subset(data_C2C3_AGG, sample_mode == 'C3_LD'), BLACKOUT=TRUE, N_days=2, SWAP_FACET=TRUE, NO_GRIDLINES=TRUE, NO_X_TEXT=TRUE)
ggsave(file="test.svg", plot=P, width=10, height=8)

reconstituting_DFs_10$reconstituted %>%
  {split(., list(.$Sample, .$Mode))} %>%
  {`names<-`(., gsub('\\.', ' ', names(.)))} %>%
  .[names(color_list) %>% gsub('_', ' ', .)] %>% 
  lapply(function(df) df$Gene) %>%
  venn.diagram(fill = color_list, filename=NULL, cat.fontface='bold') %>%
  plot_grid %>%
  ggsave(file="venn_cycle_all.svg", plot=., width=5, height=5)

reconstituting_DFs_10$reconstituted %>%
  subset(meta2d_phase_AGG_norm > 16 & meta2d_phase_AGG_norm < 24) %>%
  {split(., list(.$Sample, .$Mode))} %>%
  {`names<-`(., gsub('\\.', ' ', names(.)))} %>%
  .[names(color_list) %>% gsub('_', ' ', .)] %>% 
  lapply(function(df) df$Gene) %>%
  venn.diagram(fill = color_list, filename=NULL, cat.fontface='bold') %>%
  plot_grid

reconstituting_DFs_10$reconstituted %>%
  subset(meta2d_phase_AGG_norm > 16 & meta2d_phase_AGG_norm < 24) %>%
  {split(., list(.$Sample, .$Mode))} %>% 
  lapply(function(df) df$Gene) %>%
  {cat(paste('LD:', paste(intersect(.$C2.LD, .$C3.LD), collapse=', '), 
             '\nC3:', paste(intersect(.$C3.LD, .$C3.DD), collapse=', '),
             '\nC2:', paste(intersect(.$C2.LD, .$C2.DD), collapse=', ')))}
             
subset(to_filter, second_pct>0.10) %>%
  {split(., list(.$Sample, .$Mode))} %>%
  {`names<-`(., gsub('\\.', ' ', names(.)))} %>%
  .[names(color_list) %>% gsub('_', ' ', .)] %>% 
  lapply(function(df) df$Gene) %>%
  venn.diagram(fill = color_list, filename=NULL, cat.fontface='bold') %>%
  plot_grid %>%
  ggsave(file="venn_QC_all.svg", plot=., width=5, height=5)


subset(to_filter, second_pct>0.10) %>%
  {split(., list(.$Sample, .$Mode))} %>%
  {`names<-`(., gsub('\\.', ' ', names(.)))} %>%
  .[names(color_list) %>% gsub('_', ' ', .)] %>% 
  lapply(function(df) df$Gene) %>%
  {lapply(., function(c1) {
    lapply(., function(c2) length(unique(c(c1, c2)))) %>%
      do.call(what='cbind')
  })} %>% do.call(what='rbind') %>%
  {`rownames<-`(., colnames(.))}


reconstituting_DFs_10$reconstituted %>%
  {split(., list(.$Sample, .$Mode))} %>%
  {`names<-`(., gsub('\\.', ' ', names(.)))} %>%
  .[names(color_list) %>% gsub('_', ' ', .)] %>% 
  {.[startsWith(names(.), 'C3')]} %>%    # MODIFY THIS LINE
  #{.[endsWith(names(.), 'LD')]} %>%    # MODIFY THIS LINE
  lapply(function(df) df$Gene) %>%
  {venn.diagram(., fill = color_list[gsub(' ', '_', names(.))], filename=NULL, cat.fontface='bold', fontfamily = 'Arial', cat.fontfamily='Arial') %>%
      plot_grid(labels = 'C3: fisher-p = 0.09531; \nchisq-p = 0.1152; \nE(chance)=12')
  } %>%
  ggsave(file="venn_cycle_C3_LAB.svg", plot=., width=5, height=5)
fisher.test(cbind(c(21,273), c(180, 17824-273-21-180))) # Using all genes as background: 1.93e-11
fisher.test(cbind(c(21,273), c(180, 3989-273-21-180))) # Using genes passing QC in both as background: 0.09531
chisq.test(cbind(c(21,273), c(180, 3989-273-21-180))) # Using genes passing QC in both as background: 0.1152
round(273/3989*180) # 12


reconstituting_DFs_10$reconstituted %>%
  {split(., list(.$Sample, .$Mode))} %>%
  {`names<-`(., gsub('\\.', ' ', names(.)))} %>%
  .[names(color_list) %>% gsub('_', ' ', .)] %>% 
  {.[startsWith(names(.), 'C2')]} %>%    # MODIFY THIS LINE
  #{.[endsWith(names(.), 'LD')]} %>%    # MODIFY THIS LINE
  lapply(function(df) df$Gene) %>%
  {venn.diagram(., fill = color_list[gsub(' ', '_', names(.))], filename=NULL, cat.fontface='bold', fontfamily = 'Arial', cat.fontfamily='Arial') %>%
      plot_grid(labels = 'C2: fisher-p = 0.02771; \nchisq-p = 0.02956, \nE(chance)=8')
  } %>%
  ggsave(file="venn_cycle_C2_LAB.svg", plot=., width=5, height=5)
fisher.test(cbind(c(17,215), c(127, 17824-215-17-215))) # Using all genes as background: 6.356e-12
fisher.test(cbind(c(17,215), c(127, 3427-215-17-215))) # Using genes passing QC in both as background: 0.02771
chisq.test(cbind(c(17,215), c(127, 3427-215-17-215))) # Using genes passing QC in both as background: 0.02956
round(127/3427*215) # 8


reconstituting_DFs_10$reconstituted %>%
  {split(., list(.$Sample, .$Mode))} %>%
  {`names<-`(., gsub('\\.', ' ', names(.)))} %>%
  .[names(color_list) %>% gsub('_', ' ', .)] %>% 
  #{.[startsWith(names(.), 'C2')]} %>%    # MODIFY THIS LINE
  {.[endsWith(names(.), 'LD')]} %>%    # MODIFY THIS LINE
  lapply(function(df) df$Gene) %>%
  {venn.diagram(., fill = color_list[gsub(' ', '_', names(.))], filename=NULL, cat.fontface='bold', fontfamily = 'Arial', cat.fontfamily='Arial') %>%
      plot_grid(labels = 'LD: fisher-p = 1.105e-05; \nchisq-p = 1.575e-06; \nE(chance)=17')
  } %>%
  ggsave(file="venn_cycle_LD_LAB.svg", plot=., width=5, height=5)
fisher.test(cbind(c(36,258), c(196, 17824-258-36-196))) # Using all genes as background: < 2.2e-16
fisher.test(cbind(c(36,258), c(196, 4012-258-36-196))) # Using genes passing QC in both as background: 1.105e-05
chisq.test(cbind(c(36,258), c(196, 4012-258-36-196))) # Using genes passing QC in both as background: 1.575e-06
round(258/4012*258) # 17


reconstituting_DFs_10$reconstituted %>%
  {split(., list(.$Sample, .$Mode))} %>%
  {`names<-`(., gsub('\\.', ' ', names(.)))} %>%
  .[names(color_list) %>% gsub('_', ' ', .)] %>% 
  #{.[startsWith(names(.), 'C2')]} %>%    # MODIFY THIS LINE
  {.[endsWith(names(.), 'DD')]} %>%    # MODIFY THIS LINE
  lapply(function(df) df$Gene) %>%
  {venn.diagram(., fill = color_list[gsub(' ', '_', names(.))], filename=NULL, cat.fontface='bold', fontfamily = 'Arial', cat.fontfamily='Arial') %>%
      plot_grid(labels = 'DD: fisher-p = 0.1228; \nchisq-p = 0.1343; \nE(chance)=7')
  } %>%
  ggsave(file="venn_cycle_DD_LAB.svg", plot=., width=5, height=5)
fisher.test(cbind(c(12,189), c(132, 17824-189-12-132))) # Using all genes as background: 8.591e-08
fisher.test(cbind(c(12,189), c(132, 3825-189-12-132))) # Using genes passing QC in both as background: 0.1228
chisq.test(cbind(c(12,189), c(132, 3825-189-12-132))) # Using genes passing QC in both as background: 0.1343
round(189/3825*132) # 7

{plot_grid(nrow=2, 
          reconstituting_DFs_10$reconstituted %>%
            {split(., list(.$Sample, .$Mode))} %>%
            {`names<-`(., gsub('\\.', ' ', names(.)))} %>%
            .[names(color_list) %>% gsub('_', ' ', .)] %>% 
            {.[startsWith(names(.), 'C3')]} %>%    # MODIFY THIS LINE
            #{.[endsWith(names(.), 'LD')]} %>%    # MODIFY THIS LINE
            lapply(function(df) df$Gene) %>%
            {venn.diagram(., fill = color_list[gsub(' ', '_', names(.))], filename=NULL, 
                          cat.fontface='bold', fontfamily = 'Arial', cat.fontfamily='Arial',
                          cat.just=rep(list(c(0.5,-2)), length(names(.)))) %>%
                plot_grid(#labels = 'C3: fisher-p = 0.09531; chisq-p = 0.1152; E(chance)=12'
                          )
            },
          reconstituting_DFs_10$reconstituted %>%
            {split(., list(.$Sample, .$Mode))} %>%
            {`names<-`(., gsub('\\.', ' ', names(.)))} %>%
            .[names(color_list) %>% gsub('_', ' ', .)] %>% 
            {.[startsWith(names(.), 'C2')]} %>%    # MODIFY THIS LINE
            #{.[endsWith(names(.), 'LD')]} %>%    # MODIFY THIS LINE
            lapply(function(df) df$Gene) %>%
            {venn.diagram(., fill = color_list[gsub(' ', '_', names(.))], filename=NULL, 
                          cat.fontface='bold', fontfamily = 'Arial', cat.fontfamily='Arial',
                          cat.just=rep(list(c(0.5,-2)), length(names(.)))) %>%
                plot_grid(#labels = 'C2: fisher-p = 0.02771; chisq-p = 0.02956, E(chance)=8'
                          )
            },
          reconstituting_DFs_10$reconstituted %>%
            {split(., list(.$Sample, .$Mode))} %>%
            {`names<-`(., gsub('\\.', ' ', names(.)))} %>%
            .[names(color_list) %>% gsub('_', ' ', .)] %>% 
            #{.[startsWith(names(.), 'C2')]} %>%    # MODIFY THIS LINE
            {.[endsWith(names(.), 'LD')]} %>%    # MODIFY THIS LINE
            lapply(function(df) df$Gene) %>%
            {venn.diagram(., fill = color_list[gsub(' ', '_', names(.))], filename=NULL, 
                          cat.fontface='bold', fontfamily = 'Arial', cat.fontfamily='Arial',
                          cat.just=rep(list(c(0.5,-2)), length(names(.)))) %>%
                plot_grid(#labels = 'LD: fisher-p = 1.105e-05; chisq-p = 1.575e-06; E(chance)=17'
                          )
            },
          reconstituting_DFs_10$reconstituted %>%
            {split(., list(.$Sample, .$Mode))} %>%
            {`names<-`(., gsub('\\.', ' ', names(.)))} %>%
            .[names(color_list) %>% gsub('_', ' ', .)] %>% 
            #{.[startsWith(names(.), 'C2')]} %>%    # MODIFY THIS LINE
            {.[endsWith(names(.), 'DD')]} %>%    # MODIFY THIS LINE
            lapply(function(df) df$Gene) %>%
            {venn.diagram(., fill = color_list[gsub(' ', '_', names(.))], filename=NULL, 
                          cat.fontface='bold', fontfamily = 'Arial', cat.fontfamily='Arial', 
                          cat.just=rep(list(c(0.5,-2)), length(names(.)))) %>%
                plot_grid(#labels = 'DD: fisher-p = 0.1228; chisq-p = 0.1343; E(chance)=7'
                          )
            }
) %>% ggsave(file="venn_cycle_pairs_noLAB.svg", plot=., width=10, height=10)
}




####    ####    ####    ####    ####    ####    ####    ####    ####    ####   


# 16. PHASE HEATMAPS
################################################################################
plot_grid(nrow=1, hist_phase(reconstituting_DFs_10$reconstituted, CIRCULAR=FALSE, TRIM_TO_24=TRUE, binWidth=2, n_W=0, YSCALE=TRUE, Colors = color_list), 
          hist_phase(reconstituting_DFs_10$reconstituted, CIRCULAR=TRUE, TRIM_TO_24=TRUE, binWidth=2, n_W=0, Colors=color_list), 
          hmap_phase(reconstituting_DFs_10$reconstituted, Colors = color_list_dark)) %>%
  {lapply(c('.svg', '.pdf', '.png'), function(ext)  ggsave(file=paste0('phaseHmap_nW0_dark', ext), plot=., width=5, height=5))}

plot_grid(nrow=1, hist_phase(reconstituting_DFs_10$reconstituted, CIRCULAR=FALSE, TRIM_TO_24=TRUE, binWidth=2, n_W=0, YSCALE=TRUE, Colors = color_list), 
          hist_phase(reconstituting_DFs_10$reconstituted, CIRCULAR=TRUE, TRIM_TO_24=TRUE, binWidth=2, n_W=0, Colors=color_list), 
          hmap_phase(reconstituting_DFs_10$reconstituted, Colors = list(C3_LD="#945CB0", C3_DD="#945CB0", C2_LD="#945CB0", C2_DD="#945CB0"), col_low = "#308663")) %>%
  {lapply(c('.svg', '.pdf', '.png'), function(ext)  ggsave(file=paste0('phaseHmap_nW0_InesCols', ext), plot=., width=5, height=5))}

plot_grid(nrow=1, hist_phase(reconstituting_DFs_10$reconstituted, CIRCULAR=FALSE, TRIM_TO_24=TRUE, binWidth=2, n_W=0, YSCALE=TRUE, Colors=rep('blue', 4) %>% `names<-`(c('C3_LD', 'C3_DD', 'C2_LD', 'C2_DD')), HIST_FRONT=TRUE, LD=TRUE), 
          hist_phase(reconstituting_DFs_10$reconstituted, CIRCULAR=TRUE, TRIM_TO_24=TRUE, binWidth=2, n_W=0, Colors=rep('blue', 4) %>% `names<-`(c('C3_LD', 'C3_DD', 'C2_LD', 'C2_DD')), HIST_FRONT=TRUE, LD=TRUE), 
          hmap_phase(reconstituting_DFs_10$reconstituted, Colors = list(C3_LD="blue", C3_DD="blue", C2_LD="blue", C2_DD="blue"), col_low = "red")) %>%
  {lapply(c('.pdf'), function(ext)  ggsave(file=paste0('phaseHmap_nW0_JanV1', ext), plot=., width=5, height=5))}



plot_grid(nrow=1, hist_phase(reconstituting_DFs_10$reconstituted, CIRCULAR=FALSE, TRIM_TO_24=TRUE, binWidth=2, n_W=0, YSCALE=TRUE, Colors = color_list_alt), 
          hist_phase(reconstituting_DFs_10$reconstituted, CIRCULAR=TRUE, TRIM_TO_24=TRUE, binWidth=2, n_W=0, Colors=color_list_alt), 
          hmap_phase(reconstituting_DFs_10$reconstituted, Colors = color_list_alt)) %>%
  {lapply(c('.svg', '.pdf', '.png'), function(ext)  ggsave(file=paste0('phaseHmap_nW0_alt', ext), plot=., width=5, height=5))}


plot_grid(nrow=1, hist_phase(reconstituting_DFs_10$reconstituted, CIRCULAR=FALSE, TRIM_TO_24=TRUE, binWidth=2, n_W=0, YSCALE=TRUE, Colors = color_list_alt), 
          hist_phase(reconstituting_DFs_10$reconstituted, CIRCULAR=TRUE, TRIM_TO_24=TRUE, binWidth=2, n_W=0, Colors=color_list_alt), 
          hmap_phase(reconstituting_DFs_10$reconstituted, Colors = color_list_alt, col_low='black')) %>% 
  {lapply(c('.svg', '.pdf', '.png'), function(ext)  ggsave(file=paste0('phaseHmap_nW0_alt_black', ext), plot=., width=5, height=5))}


plot_grid(nrow=1, hist_phase(reconstituting_DFs_10$reconstituted, CIRCULAR=FALSE, TRIM_TO_24=TRUE, binWidth=1, n_W=0, YSCALE=TRUE, Colors = color_list_alt), 
          hist_phase(reconstituting_DFs_10$reconstituted, CIRCULAR=TRUE, TRIM_TO_24=TRUE, binWidth=1, n_W=0, Colors=color_list_alt), 
          hmap_phase(reconstituting_DFs_10$reconstituted, Colors = color_list_alt)) %>%
  {lapply(c('.svg', '.pdf', '.png'), function(ext)  ggsave(file=paste0('phaseHmap_nW0_1h_alt', ext), plot=., width=5, height=5))}

plot_grid(nrow=1, hist_phase(reconstituting_DFs_10$reconstituted, CIRCULAR=FALSE, TRIM_TO_24=TRUE, binWidth=2, n_W=0, YSCALE=TRUE, Colors = color_list_alt, LD=FALSE), 
          hist_phase(reconstituting_DFs_10$reconstituted, CIRCULAR=TRUE, TRIM_TO_24=TRUE, binWidth=2, n_W=0, Colors=color_list_alt, LD=FALSE), 
          hmap_phase(reconstituting_DFs_10$reconstituted, Colors = color_list_alt, )) %>%
  {lapply(c('.svg', '.pdf', '.png'), function(ext)  ggsave(file=paste0('phaseHmap_nW0_alt_DD', ext), plot=., width=5, height=5))}
  
plot_grid(nrow=1, hist_phase(reconstituting_DFs_10$reconstituted, CIRCULAR=FALSE, TRIM_TO_24=TRUE, binWidth=2, n_W=0, YSCALE=FALSE), 
          hist_phase(reconstituting_DFs_10$reconstituted, CIRCULAR=TRUE, TRIM_TO_24=TRUE, binWidth=2, n_W=0), 
          hmap_phase(reconstituting_DFs_10$reconstituted)) %>%
  ggsave(file='phaseHmap_nW0_dens.svg', plot=., width=5, height=5)
plot_grid(nrow=1, hist_phase(reconstituting_DFs_10$reconstituted, CIRCULAR=FALSE, TRIM_TO_24=TRUE, binWidth=2, n_W=3, YSCALE=TRUE), 
          hist_phase(reconstituting_DFs_10$reconstituted, CIRCULAR=TRUE, TRIM_TO_24=TRUE, binWidth=2, n_W=3), 
          hmap_phase(reconstituting_DFs_10$reconstituted)) %>%
  ggsave(file='phaseHmap_nW3.svg', plot=., width=5, height=5)
plot_grid(nrow=1, hist_phase(reconstituting_DFs_10$reconstituted, CIRCULAR=FALSE, TRIM_TO_24=TRUE, binWidth=2, n_W=3, YSCALE=FALSE), 
          hist_phase(reconstituting_DFs_10$reconstituted, CIRCULAR=TRUE, TRIM_TO_24=TRUE, binWidth=2, n_W=3), 
          hmap_phase(reconstituting_DFs_10$reconstituted)) %>%
  ggsave(file='phaseHmap_nW3_dens.svg', plot=., width=5, height=5)

hmap_phase(reconstituting_DFs_10$reconstituted, sample_remove = 'C2', mode_remove = 'DD', Colors = color_list) %>%
  ggsave(file='hmap_C3LD.svg', plot=., width=5, height=5)
hmap_phase(reconstituting_DFs_10$reconstituted, sample_remove = 'C3', mode_remove = 'DD', Colors = color_list) %>%
  ggsave(file='hmap_C2LD.svg', plot=., width=5, height=5)
hmap_phase(reconstituting_DFs_10$reconstituted, sample_remove = 'C2', mode_remove = 'LD', Colors = color_list) %>%
  ggsave(file='hmap_C3DD.svg', plot=., width=5, height=5)
hmap_phase(reconstituting_DFs_10$reconstituted, sample_remove = 'C3', mode_remove = 'LD', Colors = color_list) %>%
  ggsave(file='hmap_C2DD.svg', plot=., width=5, height=5)

plot_grid(nrow=1, 
          plot_grid(nrow=1, hist_phase(reconstituting_DFs_10$reconstituted, CIRCULAR=FALSE, TRIM_TO_24=TRUE, binWidth=2), 
                    hist_phase(reconstituting_DFs_10$reconstituted, CIRCULAR=TRUE, TRIM_TO_24=TRUE, binWidth=2), 
                    ggplot()+theme_minimal()), 
          plot_grid(nrow=1, hist_phase(reconstituting_DFs_10$reconstituted, CIRCULAR=FALSE, TRIM_TO_24=TRUE, binWidth=2, n_W=2), 
                    hist_phase(reconstituting_DFs_10$reconstituted, CIRCULAR=TRUE, TRIM_TO_24=TRUE, binWidth=2, n_W=2), 
                    ggplot()+theme_minimal()), 
          plot_grid(nrow=1, hist_phase(reconstituting_DFs_10$reconstituted, CIRCULAR=FALSE, TRIM_TO_24=TRUE, binWidth=2, n_W=4), 
                    hist_phase(reconstituting_DFs_10$reconstituted, CIRCULAR=TRUE, TRIM_TO_24=TRUE, binWidth=2, n_W=4), 
                    ggplot()+theme_minimal()), 
          plot_grid(nrow=1, hist_phase(reconstituting_DFs_10$reconstituted, CIRCULAR=FALSE, TRIM_TO_24=TRUE, binWidth=2, n_W=6), 
                    hist_phase(reconstituting_DFs_10$reconstituted, CIRCULAR=TRUE, TRIM_TO_24=TRUE, binWidth=2, n_W=6), 
                    ggplot()+theme_minimal()), 
          plot_grid(nrow=1, hist_phase(reconstituting_DFs_10$reconstituted, CIRCULAR=FALSE, TRIM_TO_24=TRUE, binWidth=2, n_W=8), 
                    hist_phase(reconstituting_DFs_10$reconstituted, CIRCULAR=TRUE, TRIM_TO_24=TRUE, binWidth=2, n_W=8), 
                    ggplot()+theme_minimal()), 
          plot_grid(nrow=1, hist_phase(reconstituting_DFs_10$reconstituted, CIRCULAR=FALSE, TRIM_TO_24=TRUE, binWidth=2, n_W=10), 
                    hist_phase(reconstituting_DFs_10$reconstituted, CIRCULAR=TRUE, TRIM_TO_24=TRUE, binWidth=2, n_W=10), 
                    ggplot()+theme_minimal()), 
          plot_grid(nrow=1, hist_phase(reconstituting_DFs_10$reconstituted, CIRCULAR=FALSE, TRIM_TO_24=TRUE, binWidth=2, n_W=12), 
                    hist_phase(reconstituting_DFs_10$reconstituted, CIRCULAR=TRUE, TRIM_TO_24=TRUE, binWidth=2, n_W=12), 
                    ggplot()+theme_minimal()), 
          plot_grid(nrow=1, hist_phase(reconstituting_DFs_10$reconstituted, CIRCULAR=FALSE, TRIM_TO_24=TRUE, binWidth=2, n_W=14), 
                    hist_phase(reconstituting_DFs_10$reconstituted, CIRCULAR=TRUE, TRIM_TO_24=TRUE, binWidth=2, n_W=14), 
                    ggplot()+theme_minimal()))



####    ####    ####    ####    ####    ####    ####    ####    ####    ####   


# 17. plotting random stuff
################################################################################

# cwo: C3-LD passes meta2d 0.01 and KW 1e-3; C2-DD passes meta2d 0.05
# Pdp1: C3-LD passes meta2d 0.05, KW 1e-6; C2-LD passes meta2d 0.05, KW 0.05
# per: C3-DD passes meta2d 5e-5 (KW 0.8)
# tyf: both DD pass meta2d 0.05, both fail KW


# ARE DD PASSES GOOD? MANY, YES
to_filter %>% 
  subset(meta2d_pvalue < 0.05 & KW_scores_norm < 0.05 & Mode == 'DD') %>% 
  .[['Gene']] %>% 
  subplotter_concat(N_days=2, func = 'avg_plotter_concat_mode', n_plots = 8) # 123, image saved (DD_0.05_0.05)

to_filter %>% 
  subset(meta2d_pvalue < 0.05 & KW_scores_norm < 0.05 & Mode == 'DD') %>% 
  .[['Gene']] %>% 
  {.[startsWith(., 'Rp')]} %>%
  subplotter_concat(N_days=2, func = 'avg_plotter_concat_mode') # 21

to_filter %>% 
  subset(meta2d_pvalue < 0.05 & KW_scores_norm < 1e-3 & Mode == 'DD') %>% 
  .[['Gene']] %>% 
  subplotter_concat(N_days=2, func = 'avg_plotter_concat_mode') # 22, image saved


# QUICK RETENTION COUNTS - JTK sweep, PValue
to_filter %>%
  {lapply(c(1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6), function(m, k=0.05, df=.) {
    df = subset(df, JTK_pvalue < m & PValue < k)
    c(C3_LD = length(df[df$Mode == 'LD' & df$Sample == 'C3', 'Gene']), 
      C2_LD = length(df[df$Mode == 'LD' & df$Sample == 'C2', 'Gene']), 
      C3_DD = length(df[df$Mode == 'DD' & df$Sample == 'C3', 'Gene']), 
      C2_DD = length(df[df$Mode == 'DD' & df$Sample == 'C2', 'Gene']),
      meta2d = m, KW = k)
  })} %>% 
  do.call(what='rbind')
# QUICK RETENTION COUNTS - meta2d sweep
to_filter %>%
    {lapply(c(1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6), function(m, k=0.05, df=.) {
      df = subset(df, meta2d_pvalue < m & KW_scores_norm < k)
      c(C3_LD = length(df[df$Mode == 'LD' & df$Sample == 'C3', 'Gene']), 
        C2_LD = length(df[df$Mode == 'LD' & df$Sample == 'C2', 'Gene']), 
        C3_DD = length(df[df$Mode == 'DD' & df$Sample == 'C3', 'Gene']), 
        C2_DD = length(df[df$Mode == 'DD' & df$Sample == 'C2', 'Gene']),
        meta2d = m, KW = k)
    })} %>% 
  do.call(what='rbind')
# QUICK RETENTION COUNTS - KW sweep
to_filter %>%
  {lapply(c(1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6), function(k, m=0.2, df=.) {
    df = subset(df, meta2d_pvalue < m & KW_scores_norm < k)
    c(C3_LD = length(df[df$Mode == 'LD' & df$Sample == 'C3', 'Gene'] %>% unique), 
      C2_LD = length(df[df$Mode == 'LD' & df$Sample == 'C2', 'Gene'] %>% unique), 
      C3_DD = length(df[df$Mode == 'DD' & df$Sample == 'C3', 'Gene'] %>% unique), 
      C2_DD = length(df[df$Mode == 'DD' & df$Sample == 'C2', 'Gene'] %>% unique),
      meta2d = m, KW = k)
  })} %>% 
  do.call(what='rbind')

# QUICK RETENTION COUNTS - PValue sweep
to_filter %>%
  {lapply(c(1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6), function(k, m=0.05, df=.) {
    df = subset(df, meta2d_pvalue < m & PValue < k)
    c(C3_LD = length(df[df$Mode == 'LD' & df$Sample == 'C3', 'Gene'] %>% unique), 
      C2_LD = length(df[df$Mode == 'LD' & df$Sample == 'C2', 'Gene'] %>% unique), 
      C3_DD = length(df[df$Mode == 'DD' & df$Sample == 'C3', 'Gene'] %>% unique), 
      C2_DD = length(df[df$Mode == 'DD' & df$Sample == 'C2', 'Gene'] %>% unique),
      meta2d = m, PValue = k)
  })} %>% 
  do.call(what='rbind')



to_filter %>% 
  subset(JTK_pvalue < 0.05 & PValue < 0.05 & Mode == 'DD') %>% 
  .[['Gene']] %>% unique %>% 
  subplotter_concat(N_days=2)

to_filter %>% 
  subset(JTK_pvalue < 0.05 & PValue < 0.05 & Sample == 'C3' & Mode == 'LD') %>% 
  .[['Gene']] %>% unique %>% 
  subplotter_concat(N_days=2, samples='C3')

# BH correct p-values, then plot ROC (with both) and note comparable values
# Capture screenshots
# Make rainbow plots to show optimal regimes and make optimal amplitude filters
# Give Sebastian options: Using x filters I got a score of X, using y filters...

to_filter %>% 
  subset(meta2d_pvalue < 0.05 & KW_scores_norm < 0.05 & Mode == 'DD') %>% 
  .[['Gene']] %>% 
  subplotter_concat(N_days=2, func = 'avg_plotter_concat_mode', n_plots = 8) # 123, image saved (DD_0.05_0.05)

  

filtered_preliminary <- to_filter %>% 
  subset(meta2d_pvalue < 0.05 & KW_scores_norm < 0.05)

table(filtered_preliminary$Mode, filtered_preliminary$Sample) 
# for meta2d 0.01, KW 0.05: C3-LD: 132, C2-LD: 71, C3-DD: 36, C2-DD: 25
# for meta2d 0.05, KW 0.05: C3-LD: 212, C2-LD: 150, C3-DD: 78, C2-DD: 45

filtered_preliminary %>% 
  subset(Mode == 'DD') %>%
  .[['Gene']] %>% 
  subplotter_concat(N_days=2, func = 'avg_plotter_concat_mode') # 61, none in common, a few don't bad for 0.01

filtered_preliminary[filtered_preliminary$Mode == 'DD', 'Gene'] %in% filtered_preliminary[filtered_preliminary$Mode == 'LD', 'Gene'] %>%
  {length(which(.))/length(.)} # 13% of DDs are in LD; Kr-h1 is in C2 DD (for 0.01 at least)

filtered_preliminary[filtered_preliminary$Mode == 'DD', 'Gene'] %>%
  .[.%in% filtered_preliminary[filtered_preliminary$Mode == 'LD', 'Gene']] %>%
  avg_plotter_concat_mode(N_days=2, rib_alpha = 0.4) # only 8, but some look pretty great (RpL26, RpS25)

genes_preliminary <- filtered_preliminary %>% 
  {.[.$Mode == 'LD', 'Gene']} %>% unique

FindMarkers(data_C2C3, ident.1='C2', group.by = 'sample')
####    ####    ####    ####    ####    ####    ####    ####    ####    ####   


# 18. DISTRIBUTIONS (USED)
################################################################################
WD_Yerbol <- '/Users/Teddy/Desktop/Yerbol_project/'

markers_G_C <- read.table(paste0(WD_Yerbol, 'G_C_markers_0.5_0.05.txt'))[[1]]
markers_R_C <- read.table(paste0(WD_Yerbol, 'R_C_markers_0.5_0.05.txt'))[[1]]

markers_Gany_C <- read.table(paste0(WD_Yerbol, 'G_C_markers_0.05.txt'))[[1]]
markers_Rany_C <- read.table(paste0(WD_Yerbol, 'R_C_markers_0.05.txt'))[[1]]

markers_Gany_C01 <- read.table(paste0(WD_Yerbol, 'G_C_markers_0.01.txt'))[[1]]
markers_Rany_C01 <- read.table(paste0(WD_Yerbol, 'R_C_markers_0.01.txt'))[[1]]

markers_G_N <- read.table(paste0(WD_Yerbol, 'G_N_markers_0.5_0.05.txt'))[[1]]
markers_R_N <- read.table(paste0(WD_Yerbol, 'R_N_markers_0.5_0.05.txt'))[[1]]

manual_R_N <- c('ninaE', paste0('Rh', 2:6),
                'trp', 'trpl', 'ninaA', # from aging single cell atlas marker genes
                'HisCl1' # from histamine KD paper
                )

manual_G_N <- c('repo', 'alrm', 'wrapper', 'Hml', 'moody', 'Indy', 'Vmat', # from fly cell atlas/shegal
                'e', #https://onlinelibrary.wiley.com/doi/10.1002/cne.10360
                'pnt', 'ttk' #https://doi.org/10.1242/dev.00468)
) 
#               'HisT', #https://www.science.org/doi/10.1126/sciadv.abq1780
#               'bdl', #https://doi.org/10.1016/j.ydbio.2016.04.020

contaminants <- list(markers_G_C = markers_G_C,
                     markers_R_C = markers_R_C,
                     Comb_GO_ribo_mitoG_N = markers_G_N,
                     markers_R_N = markers_R_N, 
                     markers_Gany_C = markers_Gany_C,
                     markers_Rany_C = markers_Rany_C,
                     markers_Gany_C01 = markers_Gany_C01,
                     markers_Rany_C01 = markers_Rany_C01, 
                     manual_R_N = manual_R_N,
                     manual_G_N = manual_G_N)


plot_gene_expression_distributions = function(Genes = circ_genes, text_genes=circ_genes, Expr_data=expr_data, kernel_genes = NULL,
                                              AVG_NOT_PROP=FALSE, LOGSWAP=FALSE, LEGEND=TRUE, ymax=NA, xmax=NA, add_line=0, kernel_adjust=1, circ_adjust=0.1, add_circ=NULL) {
  if (!is.null(add_circ)) {
    Expr_data = rbind(cbind(Expr_data, sepCirc = FALSE), 
                      cbind(merge(Expr_data, subset(add_circ, select=c(Gene, Sample, Mode))), sepCirc = TRUE))
  } else Expr_data = cbind(Expr_data, sepCirc=FALSE)
  
  Expr_data[, c('Gene', 'Sample', 'Mode', ifelse(AVG_NOT_PROP, 'avg', 'Cell_proportion'), 'sepCirc')] %>%
    `colnames<-`(., c('Gene', 'Sample', 'Mode', 'value', 'sepCirc')) %>% 
    mutate(value = as.numeric(value),
           Sample = factor(Sample, levels=c('C3', 'C2')), 
           Mode = factor(Mode, levels=c('LD', 'DD'))) %>%
    arrange(desc(value)) %>%
    mutate(Gene = factor(Gene, levels=unique(Gene))) %>%
    {ggplot(., aes(value)) +
        geom_density(data=subset(., sepCirc==FALSE)) +
        #geom_density(data=subset(., sepCirc==TRUE), linetype='dashed', color='gray', adjust=circ_adjust) +
        geom_histogram(binwidth = 500/nrow(.), data=subset(., sepCirc==TRUE), aes(y=after_stat(density)), fill="transparent", color='gray') +
        {if(!is.null(kernel_genes)) geom_density(data = subset(., Gene %in% kernel_genes & sepCirc==FALSE), color = 'gray45', linewidth=0.8, adjust=kernel_adjust)} +
        #{if(!is.null(kernel_genes)) geom_histogram(binwidth = 500/nrow(.), data = subset(., Gene %in% kernel_genes & sepCirc==FALSE), color = 'black')} +
        scale_y_continuous(transform = ifelse(LOGSWAP, 'log1p', 'identity'), limits = c(NA, ymax)) +
        scale_x_continuous(transform = ifelse(LOGSWAP, 'identity', 'log10'), limits = c(NA, xmax)) +
        geom_vline(data = subset(., Gene %in% Genes & sepCirc==FALSE), mapping = aes(xintercept = value, color=Gene, group=Gene), linetype='solid', linewidth = 0.3, key_glyph='rect') +
        {if (add_line) geom_vline(xintercept=add_line, linewidth=0.5)} +
        {if(!LEGEND) geom_text(data = subset(., Gene %in% text_genes & sepCirc==FALSE), aes(x = value, y = -0.12, label = Gene, color=Gene, group=Gene), angle=45)} +
        #ggrepel::geom_text_repel(data = subset(., Gene %in% text_genes), aes(x = value, y = 0, label = Gene, colour = Gene), angle = -90, vjust = -1.1, direction = 'x') +
        ylab('Density') +
        xlab(ifelse(AVG_NOT_PROP, 'Average Expression', 'Proportion of Cells Expressing')) +
        facet_grid(Mode ~ Sample) + theme_minimal() +
        {if(!LEGEND) theme(legend.position="none")} +
        scale_color_manual(values = scales::hue_pal()(length(Genes))) 
    }
}

plot_distributions_gamma = function(genes_highlight = circ_genes, markers=contaminants$markers_R_C, 
                                    Expr_data=expr_data, add_circ_df=NULL, samp, mode, to_use='second_pct', LEGEND=TRUE) {
  if (!is.null(add_circ_df)) {
    Expr_data = rbind(cbind(Expr_data, sepCirc = FALSE), 
                      cbind(merge(Expr_data, subset(add_circ_df, select=c(Gene, Sample, Mode))), sepCirc = TRUE))
  } else Expr_data = cbind(Expr_data, sepCirc=FALSE)
  
  Expr_data = Expr_data[, c('Gene', 'Sample', 'Mode', to_use, 'sepCirc')] %>%
    `colnames<-`(., c('Gene', 'Sample', 'Mode', 'value', 'sepCirc')) %>% 
    mutate(value = as.numeric(value)) %>%
    subset(Sample %in% samp & Mode %in% mode) %>%
    arrange(desc(value)) %>%
    mutate(Gene = factor(Gene, levels=unique(Gene)))
  
  Data = subset(Expr_data, Sample %in% samp & Mode %in% mode)
  
  Data_contaminants = subset(Data, Gene %in% markers)
  
  D_cont = Data_contaminants$value %>%
    as.numeric %>% 
    density %>%
    {data.frame(x=.$x, y=.$y)}
  
  D = Data$value %>%
    as.numeric %>% 
    density %>%
    {data.frame(x=.$x, y=.$y)}
  
  fit = fitdistrplus::fitdist(Data_contaminants$value, 'gamma', lower=c(0,0), method = "mme")
  
  D = cbind(D, ygamma = dgamma(D$x, fit$estimate['shape'], fit$estimate['rate']))
  
  ggplot(D, aes(x, y)) +
    geom_line() +
    geom_line(aes(x=x, y=ygamma)) +
    geom_histogram(data=subset(Data, sepCirc==TRUE), aes(x=value, y=after_stat(density)), 
                   binwidth = 250/nrow(Data), fill="transparent", color='gray', alpha=0.1) +
    geom_ribbon(aes(ymax=y, ymin=ygamma), fill='#63E5D3', alpha=0.5) +
    geom_ribbon(aes(ymax=ygamma, ymin=0), fill='#E66372', alpha=0.5) +
    geom_vline(xintercept = 0.1) + 
    scale_x_continuous(limits=c(0, 1)) + 
    scale_y_continuous(transform='log1p', limits=c(0, NA)) +
    geom_segment(data = subset(Data, Gene %in% genes_highlight & sepCirc==FALSE), 
               mapping = aes(x = value, xend=value, y=0, yend=Inf, color=Gene, group=Gene), 
               linetype='dashed', linewidth = 0.5) +
    {if(!LEGEND) theme(legend.position="none")} +
    scale_color_manual(values = scales::hue_pal()(length(genes_highlight))) +
    xlab(ifelse(to_use=='second_pct', 'second greatest percentage expression', to_use)) +

    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
    theme(axis.title.x=element_blank())
}
plot_grid(nrow=1, 
          plot_distributions_gamma(samp='C3', mode='LD', markers=contaminants$markers_R_C, add_circ_df=reconstituting_DFs_10$reconstituted, to_use = 'Cell_proportion'),
          plot_distributions_gamma(samp='C3', mode='LD', markers=contaminants$markers_R_C, add_circ_df=reconstituting_DFs_10$reconstituted))


# USED IN SUPP
# manual optic lobe genes
lapply(c('C3', 'C2'), function(S) {
  plot_distributions_gamma(samp=S, mode=c('LD', 'DD'), markers=manual_G_N, add_circ_df=reconstituting_DFs_10$reconstituted) %>%
    {lapply(c('.pdf'), function(extension) ggsave(paste0('G_N_manual_', S, extension), plot=., width=3,height=2))}
})
lapply(c('C3', 'C2'), function(S) {
  plot_distributions_gamma(samp=S, mode=c('LD', 'DD'), markers=manual_R_N, add_circ_df=reconstituting_DFs_10$reconstituted) %>%
    {lapply(c('.pdf'), function(extension) ggsave(paste0('R_N_manual_', S, extension), plot=., width=3,height=2))}
})


 
# not used
lapply(names(contaminants), function(name) {
  lapply(c('C3_LD', 'C3_DD', 'C2_LD', 'C2_DD'), function(S_M) {
    S = gsub('_.*', '', S_M)
    M = gsub('.*_', '', S_M)
    plot_distributions_gamma(samp=S, mode=M, markers=contaminants[[name]], add_circ_df=reconstituting_DFs_10$reconstituted) +
      ggtitle(S_M)
  }) %>%
    ggarrange(plotlist=., common.legend=TRUE, legend='right') %>%
    {lapply(c('.svg', '.pdf', '.jpg'), function(extension) ggsave(paste0(name, extension), plot=., width=6,height=6))}
})

pcts_data_vals2 %>% mutate(Cell_proportion = second_pct) %>%
  plot_gene_expression_distributions(Genes=c(circ_genes, 'cyc'), 
                                     AVG_NOT_PROP=FALSE, LOGSWAP=TRUE, kernel_genes = markers_R_C,
                                     Expr_data = ., xmax=1, add_line = 0.1, kernel_adjust=5, add_circ = reconstituting_DFs_10$reconstituted)
%>%
  {lapply(c('.svg', '.pdf', '.jpg'), function(extension) ggsave(paste0('distribution_G_C', extension), plot=., width=5.5,height=5))}


pcts_data_vals2 %>% mutate(Cell_proportion = second_pct) %>%
  plot_gene_expression_distributions(Genes=c(circ_genes, 'cyc'), 
                                     AVG_NOT_PROP=FALSE, LOGSWAP=TRUE, kernel_genes = markers_R_C,
                                     Expr_data = ., xmax=1, add_line = 0.1, kernel_adjust=3, add_circ = reconstituting_DFs_10$reconstituted)
%>%
  {lapply(c('.svg', '.pdf', '.jpg'), function(extension) ggsave(paste0('distribution_R_C', extension), plot=., width=5.5,height=5))}

####    ####    ####    ####    ####    ####    ####    ####    ####    ####   


# 19. DISTRIBUTIONS
################################################################################
genes_contaminantsOL_57 <- read.csv('genes_contaminantsOL_57.csv')[['x']]
optic_contaminantsOL_30 <- read.csv('optic_contaminants_30.csv')[['x']]
R_N_markers_45 <- read.csv('R_N_markers_0.5_0.05.txt')[[1]]

genes_exprDist <- c(circ_genes, 'Pdfr', 
                    'repo', 'ninaE', 'Rh2', 'Rh3', 'Rh4', 'Rh5', 'Rh6', 
                    genes_contaminantsOL_57, optic_contaminantsOL_30)


pcts_data_vals2 %>% mutate(Cell_proportion = second_pct) %>%
  plot_gene_expression_distributions(Genes=c(circ_genes, 'cyc'), text_genes = c('cwo', 'per', 'Pdp1', 'vri', 'cyc'),
                                     AVG_NOT_PROP=FALSE, LOGSWAP=TRUE, kernel_genes = c(genes_contaminantsOL_57),
                                     Expr_data = ., xmax=1, add_line = 0.1, kernel_adjust=5, add_circ = reconstituting_DFs_10$reconstituted)
  

plot_grid(nrow=1, 
          pcts_data_vals2 %>% mutate(Cell_proportion = second_pct) %>%
            plot_gene_expression_distributions(Genes=c(circ_genes, 'cyc'), text_genes = c('cwo', 'per', 'Pdp1', 'vri', 'cyc'),
                                               AVG_NOT_PROP=FALSE, LOGSWAP=TRUE, kernel_genes = genes_contaminantsOL_57,
                                               Expr_data = ., xmax=1, add_line = 0.1, kernel_adjust=5), 
          pcts_data_vals2 %>% mutate(Cell_proportion = second_pct) %>%
            plot_gene_expression_distributions(Genes=c(circ_genes, 'cyc'), text_genes = c('cwo', 'per', 'Pdp1', 'vri', 'cyc'),
                                               AVG_NOT_PROP=FALSE, LOGSWAP=TRUE, kernel_genes = optic_contaminantsOL_30,
                                               Expr_data = ., xmax=1, add_line = 0.1, kernel_adjust=5), 
          pcts_data_vals2 %>% mutate(Cell_proportion = second_pct) %>%
            plot_gene_expression_distributions(Genes=c(circ_genes, 'cyc'), text_genes = c('cwo', 'per', 'Pdp1', 'vri', 'cyc'),
                                               AVG_NOT_PROP=FALSE, LOGSWAP=TRUE, kernel_genes = R_N_markers_45,
                                               Expr_data = ., xmax=1, add_line = 0.1, kernel_adjust=5))

pcts_data_vals2 %>% mutate(Cell_proportion = second_pct) %>%
  plot_gene_expression_distributions(Genes=c(circ_genes, 'cyc'), text_genes = c('cwo', 'per', 'Pdp1', 'vri', 'cyc'),
                                     AVG_NOT_PROP=FALSE, LOGSWAP=TRUE, kernel_genes = genes_contaminantsOL_57,
                                     Expr_data = ., xmax=1, add_line = 0.1, kernel_adjust=5) %>%
  ggsave(file="glial_neur_contaminants_adj5_lg.svg", plot=., width=10, height=10)

pcts_data_vals2 %>% mutate(Cell_proportion = second_pct) %>%
  plot_gene_expression_distributions(Genes=c(circ_genes, 'cyc'), text_genes = c(circ_genes, 'cyc'),
                                     AVG_NOT_PROP=FALSE, LOGSWAP=TRUE, kernel_genes = genes_contaminantsOL_57,
                                     Expr_data = ., xmax=1, add_line = 0.1, kernel_adjust=5) %>%
  ggsave(file="glial_neur_contaminants_adj5_lg_circ.svg", plot=., width=10, height=10)

avg_plotter_repDot(c('cwo', 'per', 'vri'), data = data_C2C3_AGG, NO_GRIDLINES = TRUE) %>%
  ggsave(file="CCG_repdot.svg", plot=., width=5, height=5)
avg_plotter_concat_mode(c('cwo', 'per', 'vri'), data = data_C2C3_AGG, NO_GRIDLINES = TRUE) %>%
  ggsave(file="CCG_concat_mode.svg", plot=., width=5, height=5)

pcts_data_vals2 %>% mutate(Cell_proportion = second_pct) %>%
  plot_gene_expression_distributions(Genes=c(circ_genes, 'cyc'), text_genes = c('cwo', 'per', 'Pdp1', 'vri', 'cyc'),
                                     AVG_NOT_PROP=FALSE, LOGSWAP=TRUE, kernel_genes = optic_contaminantsOL_30,
                                     Expr_data = ., xmax=1, add_line = 0.1, kernel_adjust=5) %>%
  ggsave(file="optic_neur_contaminants_adj5_lg.svg", plot=., width=10, height=10)


pcts_data_vals2 %>% mutate(Cell_proportion = second_pct) %>%
  plot_gene_expression_distributions(Genes=c(circ_genes, 'cyc'), text_genes = c(circ_genes, 'cyc'),
                                     AVG_NOT_PROP=FALSE, LOGSWAP=TRUE, kernel_genes = optic_contaminantsOL_30,
                                     Expr_data = ., xmax=1, add_line = 0.1, kernel_adjust=5) %>%
  ggsave(file="optic_neur_contaminants_adj5_all_lg.svg", plot=., width=10, height=10)

plot_grid(nrow=1, 
          pcts_data_vals2 %>% mutate(Cell_proportion = second_pct) %>%
            plot_gene_expression_distributions(Genes=circ_genes, text_genes = c('cwo', 'per', 'Pdp1', 'vri'),
                                               AVG_NOT_PROP=FALSE, LOGSWAP=TRUE, kernel_genes = genes_contaminantsOL_57,
                                               Expr_data = ., xmax=1, add_line = 0.1, kernel_adjust=10),
          pcts_data_vals2 %>% mutate(Cell_proportion = second_pct) %>%
            plot_gene_expression_distributions(Genes=circ_genes, text_genes = c('cwo', 'per', 'Pdp1', 'vri'),
                                               AVG_NOT_PROP=FALSE, LOGSWAP=TRUE, kernel_genes = genes_contaminantsOL_57,
                                               Expr_data = ., xmax=1, add_line = 0.1, kernel_adjust=5),
          pcts_data_vals2 %>% mutate(Cell_proportion = second_pct) %>%
            plot_gene_expression_distributions(Genes=circ_genes, text_genes = c('cwo', 'per', 'Pdp1', 'vri'),
                                               AVG_NOT_PROP=FALSE, LOGSWAP=TRUE, kernel_genes = genes_contaminantsOL_57,
                                               Expr_data = ., xmax=1, add_line = 0.1, kernel_adjust = 1))



plot_grid(nrow=1, plot_gene_expression_distributions(genes_exprDist, AVG_NOT_PROP=FALSE, 
                                                     Expr_data = merge(avg_data, pcts_data)), 
          plot_gene_expression_distributions(genes_exprDist, AVG_NOT_PROP=TRUE, 
                                             Expr_data = merge(avg_data, pcts_data)))



plot_grid(nrow=1, plot_gene_expression_distributions(Genes=genes_exprDist, text_genes = c('cwo', 'per', 'vri', 'Rh4', 'Rh5'), 
                                                     AVG_NOT_PROP=FALSE, LOGSWAP=TRUE, 
                                                    Expr_data = merge(avg_data, pcts_data), xmax=0.5), 
          plot_gene_expression_distributions(Genes=genes_exprDist, text_genes = c('cwo', 'per', 'vri', 'Rh4', 'Rh5'),
                                             AVG_NOT_PROP=TRUE, LOGSWAP=TRUE, 
                                             Expr_data = merge(avg_data, pcts_data), xmax=2))

plot_grid(nrow=1, 
          pcts_data_vals2 %>%
            subset(Gene %in% optic_contaminantsOL_30) %>%
            mutate(val = as.numeric(most_pct)) %>%
            {split(., list(.$Sample, .$Mode))} %>% 
            lapply(function(x) cbind(x, twoSD = quantile.density(density(x$val, adjust=1), probs = 0.954))) %>%
            do.call(what='rbind') %>%
            {ggplot(., aes(val)) +
                geom_density() +
                geom_vline(aes(xintercept = twoSD)) +
                facet_grid(Sample ~ Mode) + 
                ggtitle('max %')}, 
          pcts_data_vals2 %>%
            subset(Gene %in% optic_contaminantsOL_30) %>%
            mutate(val = as.numeric(second_pct)) %>%
            {split(., list(.$Sample, .$Mode))} %>% 
            lapply(function(x) cbind(x, twoSD = quantile.density(density(x$val, adjust=1), probs = 0.954))) %>%
            do.call(what='rbind') %>%
            {ggplot(., aes(val)) +
                geom_density() +
                geom_vline(aes(xintercept = twoSD)) +
                facet_grid(Sample ~ Mode) + 
                ggtitle('second max %')})



plot_grid(nrow=1, 
          pcts_data %>% 
            plot_gene_expression_distributions(Genes=circ_genes, text_genes = c('cwo', 'per', 'Pdp1', 'vri'),
                                               AVG_NOT_PROP=FALSE, LOGSWAP=TRUE, 
                                               Expr_data = ., xmax=1) + ggtitle('avg %'),
          pcts_data_vals2 %>% mutate(Cell_proportion = most_pct) %>%
            plot_gene_expression_distributions(Genes=circ_genes, text_genes = c('cwo', 'per', 'Pdp1', 'vri'), kernel_genes = genes_contaminantsOL_57,
                                               AVG_NOT_PROP=FALSE, LOGSWAP=TRUE, 
                                               Expr_data = ., xmax=1) + ggtitle('max %'),
          pcts_data_vals2 %>% mutate(Cell_proportion = second_pct) %>%
            plot_gene_expression_distributions(Genes=circ_genes, text_genes = c('cwo', 'per', 'Pdp1', 'vri'),
                                               AVG_NOT_PROP=FALSE, LOGSWAP=TRUE, 
                                               Expr_data = ., xmax=1, solid_line = 0.1) + ggtitle('second-most %'))



take_percentiles = function(values, percentiles = c(0.954), Kernel_adjust=1, x_min=0) {
  df_density = density(values, adjust=Kernel_adjust) %>%
    {data.frame(x=.$x, y=.$y)}
  df_density = subset(df_density, x>x_min)
  quantiles = quantile(df_density$y, prob=percentiles)
}


pcts_data %>%
  subset(Gene %in% optic_contaminantsOL_30) %>%
  {split(., list(.$Sample, .$Mode))} %>% 
  lapply(function(x) cbind(x, take_percentiles(x$Cell_proportion, Kernel_adjust = 1/1, x_min=0.000))) %>%
  do.call(what='rbind') %>% 
  `colnames<-`(c('Cell_proportion', 'Gene', 'Sample', 'Mode', 'x')) %>%
  {ggplot(., aes(Cell_proportion)) +
  geom_density() +
  geom_vline(aes(xintercept = x)) +
  facet_grid(Sample ~ Mode)}


pcts_data %>%
  subset(Gene %in% optic_contaminantsOL_30) %>%
  {split(., list(.$Sample, .$Mode))} %>% 
  lapply(function(df) {
    cbind(df, rank = rank(df$Cell_proportion)) %>%
      mutate(percentile = rank/max(rank)) %>%
      arrange(percentile) %>%
      {.[which.min(abs(.$percentile-0.954)),]}
  })

plot_grid(nrow = 1,
  pcts_data %>%
    subset(Gene %in% optic_contaminantsOL_30) %>%
    {split(., list(.$Sample, .$Mode))} %>% 
    lapply(function(x) cbind(x, twoSD = quantile.density(density(x$Cell_proportion, adjust=1), probs = 0.954))) %>%
    do.call(what='rbind') %>%
    {ggplot(., aes(Cell_proportion)) +
        geom_density() +
        geom_vline(aes(xintercept = twoSD)) +
        facet_grid(Sample ~ Mode)}, 
  pcts_data %>%
    subset(Gene %in% optic_contaminantsOL_30) %>%
    {split(., list(.$Sample, .$Mode))} %>% 
    lapply(function(x) cbind(x, quantile(fitdist(x$Cell_proportion,  distr = "gamma", method = "mme"), probs = 0.954)[[1]])) %>%
    do.call(what='rbind') %>%
    {ggplot(., aes(Cell_proportion)) +
        geom_density() +
        geom_vline(aes(xintercept = `p=0.954`)) +
        facet_grid(Sample ~ Mode)})

plot_grid(nrow = 1,
          avg_data %>%
            subset(Gene %in% optic_contaminantsOL_30) %>%
            {split(., list(.$Sample, .$Mode))} %>% 
            lapply(function(x) cbind(x, twoSD = quantile.density(density(x$avg, adjust=1), probs = 0.954))) %>%
            do.call(what='rbind') %>%
            {ggplot(., aes(avg)) +
                geom_density() +
                geom_vline(aes(xintercept = twoSD)) +
                facet_grid(Sample ~ Mode)}, 
          avg_data %>%
            subset(Gene %in% optic_contaminantsOL_30) %>%
            {split(., list(.$Sample, .$Mode))} %>% 
            lapply(function(x) cbind(x, quantile(fitdist(x$avg,  distr = "gamma", method = "mme"), probs = 0.954)[[1]])) %>%
            do.call(what='rbind') %>%
            {ggplot(., aes(avg)) +
                geom_density() +
                geom_vline(aes(xintercept = `p=0.954`)) +
                facet_grid(Sample ~ Mode)})

hi$Cell_proportion %>% ecdf %>% plot
library(fitdistrplus)
fit.gamma = fitdist(hi$Cell_proportion, distr = "gamma", method = "mme")
plot(fit.gamma)
quantile(fit.gamma, probs=0.954)[[1]]

cdfcomp(fit.gamma)
abline(h=0.954, col="blue")

{par(mfrow = c(2, 4))
pcts_data %>%
  subset(Gene %in% optic_contaminantsOL_30) %>%
  {split(., list(.$Sample, .$Mode))} %>% 
  for (x in .) {
    fit.gamma = fitdist(x$Cell_proportion, distr = "gamma", method = "mme")
    cdfcomp(fit.gamma)
    abline(h=0.954, col="blue")
  }
pcts_data %>%
  subset(Gene %in% optic_contaminantsOL_30) %>%
  {split(., list(.$Sample, .$Mode))} %>% 
  for (x in .) {
    fit.gamma = fitdist(x$Cell_proportion, distr = "gamma", method = "mme")
    denscomp(fit.gamma, xlim=c(0, 0.5), ylim=c(0,20))
    abline(v=quantile(fit.gamma, probs=0.954)[[1]], col="blue")
  }
}


{par(mfrow = c(2, 4))
pcts_data %>%
  subset(Gene %in% genes_contaminantsOL_57) %>%
  {split(., list(.$Sample, .$Mode))} %>% 
  for (x in .) {
    fit.gamma = fitdist(x$Cell_proportion, distr = "gamma", method = "mme")
    cdfcomp(fit.gamma)
    abline(h=0.954, col="blue")
  }
pcts_data %>%
  subset(Gene %in% genes_contaminantsOL_57) %>%
  {split(., list(.$Sample, .$Mode))} %>% 
  for (x in .) {
    fit.gamma = fitdist(x$Cell_proportion, distr = "gamma", method = "mme")
    denscomp(fit.gamma, xlim=c(0, 0.5), ylim=c(0,20))
    abline(v=quantile(fit.gamma, probs=0.954)[[1]], col="blue")
  }
}





plot_grid(nrow=1, 
          pcts_data %>%
            subset(Gene %in% optic_contaminantsOL_30) %>%
            {split(., list(.$Sample, .$Mode))} %>% 
            lapply(function(x) cbind(x, take_percentiles(x$Cell_proportion))) %>%
            do.call(what='rbind') %>%
            {ggplot(., aes(Cell_proportion)) +
                geom_density() +
                #geom_vline(aes(xintercept = x)) +
                facet_grid(Sample ~ Mode)},
          avg_data %>%
            subset(Gene %in% optic_contaminantsOL_30) %>%
            {split(., list(.$Sample, .$Mode))} %>% 
            lapply(function(x) cbind(x, take_percentiles(x$avg))) %>%
            do.call(what='rbind') %>% 
            {ggplot(., aes(avg)) +
                geom_density() +
                #geom_vline(aes(xintercept = x)) +
                facet_grid(Sample ~ Mode)})


####    ####    ####    ####    ####    ####    ####    ####    ####    ####   


# 20. Differential Expression Testing
################################################################################
# Differential Expression Testing - DESeq2 (NEEDS >2 samples)
pbulk_C2C3 <- AggregateExpression(data_C2C3, assays = "RNA", return.seurat = T, group.by = c("group1", 'sample1', "mode1", "sample"))
pbulk_C2C3$ID <- paste(pbulk_C2C3$group1, pbulk_C2C3$mode1, pbulk_C2C3$sample, sep='_')
Idents(pbulk_C2C3) <- 'ID'

#DE = list()
DE$C3$CT3CT15 <- FindMarkers(object = pbulk_C2C3, 
                            ident.1 = "ZT3_LD_C3", 
                            ident.2 = "ZT15_LD_C3",
                            test.use = "DESeq2")
#DE_cellbased <- list()
DE_cellbased$C3$ZT3ZT15 <- data_C2C3 %>%
  subset(sample == 'C3') %>%
  FindMarkers(ident.1 = 'ZT3', ident.2 = 'ZT15', group.by = 'group1')
DE_cellbased$C2$ZT3ZT15 <- data_C2C3 %>%
  subset(sample == 'C2') %>%
  FindMarkers(ident.1 = 'ZT3', ident.2 = 'ZT15', group.by = 'group1')

DE_cellbased$C3$ZT3ZT15['per', ] # 0.906 p-value, 1 p_val_adj
DE_cellbased$C2$ZT3ZT15['per', ] # 0.49 p-value, 1 p_val_adj

to_filter %>% subset(Gene %in% circ_genes)

####    ####    ####    ####    ####    ####    ####    ####    ####    ####   


# 21. Promoter analysis
################################################################################
genes <- genes(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
names(genes) <- genes$gene_id

# FOR Pdp1, meta2d < 0.2; KW < 0.01 (LD: C3-237, C2-148; DD:C3-73, C2-39)



FIMO_promoters <- promoters(genes, upstream = 1000, downstream=0)[
  Cyclers$ENSEMBL[Cyclers$ENSEMBL %in% names(genes)]] %>%
  get_sequence(BSgenome.Dmelanogaster.UCSC.dm6) %>% 
  do.call(what=xscat) %>%
  DNAStringSet %>% 
  `names<-`(1) %>%
  runFimo(motifs = '/Users/Teddy/Desktop/motif_databases/FLY/fly_combined.meme',
        outdir = "./",meme_path = "/opt/local/bin/")

hi <- promoters(genes, upstream = 1000, downstream=0)[
  Cyclers$ENSEMBL[Cyclers$ENSEMBL %in% names(genes)]] %>%
  get_sequence(BSgenome.Dmelanogaster.UCSC.dm6)

makeGRangesFromDataFrame

hi[[1]] 


hi %>% as.data.frame %>% ggplot(aes(x=log10(pvalue))) + geom_histogram() + scale_y_continuous(trans='log10')

hi$motif_id %>% unique %>% length

make_FIMO_run = function(to_filt=to_filter, samples = 'all', modes = 'all', add_to_name=NULL, filter_expr_list=NULL, 
                         WRITE_PROMOTERS=FALSE, RUN=TRUE, motifPath='/Users/Teddy/Desktop/motif_databases/FLY/fly_combined.meme',
                         promoter_up = 1000, promoter_down=0, cyclers=NULL, df_Gene_ENS=NULL, QNOTP=FALSE, thresh=1e-4, CONCAT=FALSE) {
  print(motifPath)
  res = list()
  res$filters = c(filter_expr_list, samples, modes)
  res$promoter_def = c(start = -promoter_up, end=promoter_down)
  res$filtered = to_filt
  
  if(samples != 'all') res$filtered = subset(res$filtered, Sample %in% samples)
  if(modes != 'all') res$filtered = subset(res$filtered, Mode %in% modes)
  
  for (filt_expr in filter_expr_list) {
    res$filtered = filterer(res$filtered, filt_expr)
  }
  
  if (is.null(cyclers)) res$candidates = res$filtered$Gene %>% unique %>%
    {AnnotationDbi::select(org.Dm.eg.db, keys=., columns=c('ENSEMBL'), keytype='SYMBOL')}
  else res$candidates = Cyclers[Cyclers$Gene %in% res$filtered$Gene, ] %>% 
    mutate(SYMBOL = Gene, .keep='unused')
  
  genes = genes(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
  names(genes) = genes$gene_id
  res$promoters = promoters(genes, upstream = promoter_up, downstream=promoter_down)[
    res$candidates$ENSEMBL[res$candidates$ENSEMBL %in% names(genes)]]
  
  res$sequences = res$promoters %>%
    get_sequence(BSgenome.Dmelanogaster.UCSC.dm6) 
  
  if(CONCAT) res$sequences = do.call(xscat, res$sequences) %>% DNAStringSet %>% `names<-`(1)
  
  if(WRITE_PROMOTERS) {
    c('promoters', samples, modes, unlist(filter_expr_list), add_to_name) %>% 
      list %>%
      sapply(paste, collapse='_') %>%
      paste0('.fa') %>% 
      writeXStringSet(res$sequences, filepath = .)
  }
  
  if (RUN) {
    res$fimo <- runFimo(res$sequences,
                        motifs = motifPath,
                        outdir = "./",meme_path = "/opt/local/bin/",
                        qv_thresh = QNOTP, thresh=thresh, text=!QNOTP, no_qvalue=!QNOTP)
    
    #print(head(res$fimo))
    #print(head(res$promoters))
    
    res$results = mergeByOverlaps(query = res$fimo, subject = res$promoters) %>%
      {data.frame(data.frame(.$`res$fimo`), gene_id=.$gene_id)} %>% 
      mutate(motif_id = gsub("_.*", "", motif_id)) %>%
      {if(is.null(df_Gene_ENS)) merge(., 
                                      AnnotationDbi::select(org.Dm.eg.db, keys=.$motif_id, columns='SYMBOL', keytype='ENSEMBL') %>% unique,
                                      by.x = 'motif_id', by.y='ENSEMBL', all.x=TRUE) %>% 
          mutate(Gene_from_alt_id = sub("[^[:alnum:] ].*", "", motif_alt_id), .keep='unused') %>% 
          mutate(TF = ifelse(is.na(SYMBOL), Gene_from_alt_id, SYMBOL), .keep='unused')
          else merge(., mutate(TF = Gene, .keep='unused'), df_Gene_ENS, by.x='motif_id', by.y='ENSEMBL') 
      } %>% 
      merge(res$candidates %>% mutate(gene_id = ENSEMBL)) %>%
      mutate(Target = SYMBOL, target_id = gene_id, .keep='unused') 
      
    res$ResinCandidates = subset(res$results, TF %in% res$candidates$SYMBOL)
  }
  
  names(res) = gsub('\\.', '_', names(res))
  
  return(res)
}


promoters(genes, 1000, 0) %>% head %>% get_sequence(BSgenome.Dmelanogaster.UCSC.dm6) %>% {xscat(.[[1]], .[[2]])}
promoters(genes, 1000, 0) %>% head %>% get_sequence(BSgenome.Dmelanogaster.UCSC.dm6) %>% do.call(what='xscat')


WD_split = '/Users/Teddy/Desktop/motif_databases/FLY/split_files'
res10_QC_Q10 <- list()
res10_QC_Q10 <- reconstituting_DFs_10$reconstituted %>%
  {split(., list(.$Sample, .$Mode))} %>%
  {lapply(names(.), function(name) {
    df = .[[name]]
    make_FIMO_run(df, samples = unique(df$Sample), modes = unique(df$Mode), 
                  motifPath = paste0(WD_split, '/onlyExpr_', name, '.meme'), cyclers = Cyclers, QNOTP = TRUE, thres=0.1)
  })}
names(res10_QC_Q10) <- names(res10)

res10_QC_Q05 <- list()
res10_QC_Q05 <- reconstituting_DFs_10$reconstituted %>%
  {split(., list(.$Sample, .$Mode))} %>%
  {lapply(names(.), function(name) {
    df = .[[name]]
    make_FIMO_run(df, samples = unique(df$Sample), modes = unique(df$Mode), 
                  motifPath = paste0(WD_split, '/onlyExpr_', name, '.meme'), cyclers = Cyclers, QNOTP = TRUE, thres=0.05)
  })}
names(res10_QC_Q05) <- names(res10)



res10_anno_comb <- list()
res10_anno_comb <- reconstituting_DFs_10$reconstituted %>%
  {split(., list(.$Sample, .$Mode))} %>%
  lapply(function(df) make_FIMO_run(df, samples = unique(df$Sample), modes = unique(df$Mode), WRITE_PROMOTERS=TRUE, cyclers = Cyclers))
lapply(res10_anno, function(res) {res$results$Target %>% unique %>% length})

res10_anno_comb_Q_2 <- list()
res10_anno_comb_Q_2 <- reconstituting_DFs_10$reconstituted %>%
  {split(., list(.$Sample, .$Mode))} %>%
  lapply(function(df) make_FIMO_run(df, samples = unique(df$Sample), modes = unique(df$Mode), cyclers = Cyclers, QNOTP = TRUE, thres=0.1))
lapply(res10_anno_comb_Q_2, function(res) {res$results$TF %>% unique %>% length})
saveRDS(res10_anno_comb_Q_2, 'res10_anno_comb_Q_2.rds')

saveRDS(res10_anno_comb_Q, 'res10_anno_comb_Q.rds')

res10_anno_comb_Q_2 <- readRDS('res10_anno_comb_Q_2.rds')

res10_anno_comb_Q_2$C3.LD$results$TF %>% 
  {.[. %in% subset(to_filter, second_pct<0.1 & Sample == 'C3' & Mode == 'LD')$Gene]} %>% unique %>% length
res10_anno_comb_Q_2$C3.DD$results$TF %>% 
  {.[. %in% subset(to_filter, second_pct<0.1 & Sample == 'C3' & Mode == 'DD')$Gene]} %>% unique %>% length
res10_anno_comb_Q_2$C2.LD$results$TF %>% 
  {.[. %in% subset(to_filter, second_pct<0.1 & Sample == 'C3' & Mode == 'LD')$Gene]} %>% unique %>% length
res10_anno_comb_Q_2$C2.DD$results$TF %>% 
  {.[. %in% subset(to_filter, second_pct<0.1 & Sample == 'C2' & Mode == 'DD')$Gene]} %>% unique %>% length


# SET SHUFFLE TO RUN, THEN MAKE ENRICHMENT TESTER


res10_all <- list()
res10_all <- reconstituting_DFs_10$reconstituted %>%
  {split(., list(.$Sample, .$Mode))} %>%
  lapply(function(df) make_FIMO_run(df, samples = unique(df$Sample), modes = unique(df$Mode), cyclers = Cyclers, QNOTP = FALSE, thres=2))

lapply(res10_all, function(res) {res$results$TF %>% unique %>% length})

# from C3LD: 245; C2LD: 207;  C3DD: 227; C2DD: 203
#      C3LD: 115; C2LD: 95;   C3DD: 110; C2DD: 97

res10_combined <- list()
res10_combined <- reconstituting_DFs_10$reconstituted %>%
  {split(., list(.$Sample, .$Mode))} %>%
  lapply(function(df) make_FIMO_run(df, samples = unique(df$Sample), modes = unique(df$Mode), 
                                    WRITE_PROMOTERS=TRUE#, motifPath='/Users/Teddy/Desktop/motif_databases/FLY/flyreg.v2.meme'
  )
  )

length(which(res10$C2.DD$results$TF %>%print %in% res10$C2.DD$candidates$SYMBOL %>% print))

res15 <- list()
res10 <- list()
res15 <- reconstituting_DFs_15$reconstituted %>%
  {split(., list(.$Sample, .$Mode))} %>%
  lapply(function(df) make_FIMO_run(df, samples = unique(df$Sample), modes = unique(df$Mode), WRITE_PROMOTERS=FALSE))

res10 <- reconstituting_DFs_10$reconstituted %>%
  {split(., list(.$Sample, .$Mode))} %>%
  lapply(function(df) make_FIMO_run(df, samples = unique(df$Sample), modes = unique(df$Mode), WRITE_PROMOTERS=FALSE))

res10_concat <- reconstituting_DFs_10$reconstituted %>%
  {split(., list(.$Sample, .$Mode))} %>%
  lapply(function(df) make_FIMO_run(df, samples = unique(df$Sample), modes = unique(df$Mode), WRITE_PROMOTERS=FALSE, CONCAT=TRUE))


res10_combined <- list()
res10_combined <- reconstituting_DFs_10$reconstituted %>%
  {split(., list(.$Sample, .$Mode))} %>%
  lapply(function(df) make_FIMO_run(df, samples = unique(df$Sample), modes = unique(df$Mode), WRITE_PROMOTERS=TRUE, 
                                    motifPath='/Users/Teddy/Desktop/motif_databases/FLY/fly_combined.meme'))

res10 <- res10 %>%
  lapply(function(df) sort(unique(df$ResinCandidates$TF))) %>%
  {`names<-`(lapply(names(.), function(df_name)
    if (length(.[[df_name]] >= 1)) {
      plot_cyc <- avg_plotter_concat(.[[df_name]], data=subset(data_C2C3_AGG, sample == gsub('\\..*', '', df_name) &
                                                                 mode1 == gsub('.*\\.', '', df_name)), N_days=2)
      append(res10[[df_name]], list(plot_cyc=plot_cyc))
    } else append(res10[[df_name]], list(plot_cyc=NULL))
    
  ), names(.))
  }

res10_QC_Q10 <- lapply(res10_QC_Q10, function(df) {
  RBPs = subset(df$filtered, Gene %in% merge(GO_RBPs, Cyclers)$Gene)$Gene %>% unique
  TFs = subset(df$filtered, Gene %in% merge(GO_TRs, Cyclers)$Gene)$Gene %>% unique
  print(c(unique(df$filtered$Sample), unique(df$filtered$Mode)))
  print('RBPs:')
  print(RBPs)
  print('TFs:')
  print(TFs)
  print('nTFs:')
  print(length(TFs))
  cat("\n")
  append(df, list(RBPs=RBPs, TFs=TFs))
})

#saveRDS(res10, 'res10.rds')
res10 <- readRDS('res10.rds')

#short <- list()
short$FIMO <- subset(to_filter, Gene %in% c(genes_activity, 'CrebA', 'CrebB', 'Pka-C1',
                                            'Pka-R1', 'Pka-C3', 'Rala', 'dnc', 'kay', 'crc')) %>%
  {split(., list(.$Sample, .$Mode))} %>%
  lapply(function(df) make_FIMO_run(df, samples = unique(df$Sample), modes = unique(df$Mode), WRITE_PROMOTERS=TRUE))


lapply(names(res10), function(name) {
  print(name)
  res = res10[[name]] %>% rename(Gene = 'Target')
  
  namebase = paste0(name, '_passQC10')
  make_pdf_targets(dfRes = res$results, dfFiltered = res$filtered, 
                   TFs_in=subset(to_filter, Sample == gsub('\\..*', '', name) & Mode == gsub('.*\\.', '', name) & Cell_proportion>0.1)$Gene,
                   namebase=namebase, splitInto = 5)
  pdf_combine_and_delete(namebase)
})

lapply(names(res10), function(name) {
  res = res10[[name]]
  namebase = paste0(name, '_inData')
  make_pdf_targets(dfRes = res$results, dfFiltered = res$filtered, 
                   TFs_in=rownames(data_C2C3),
                   namebase=namebase, splitInto = 20)
  pdf_combine_and_delete(namebase)
})

####    ####    ####    ####    ####    ####    ####    ####    ####    ####   


# 22. TO MAKE BACKGROUND FASTAs
################################################################################
to_filter %>%
  subset(second_pct > 0.1) %>%
  {split(., list(.$Sample, .$Mode))} %>%
  lapply(function(df) make_FIMO_run(df, samples = unique(df$Sample), modes = unique(df$Mode), add_to_name='background_10',
                                    WRITE_PROMOTERS=TRUE, RUN=FALSE))

to_filter %>%
  subset(second_pct > 0.15) %>%
  {split(., list(.$Sample, .$Mode))} %>%
  lapply(function(df) make_FIMO_run(df, samples = unique(df$Sample), modes = unique(df$Mode), add_to_name='background_15',
                                    WRITE_PROMOTERS=TRUE, RUN=FALSE))

{reconstituting_DFs_10$reconstituted %>%
  {split(., list(.$Sample, .$Mode))} %>%
  {lapply(names(.), function(df_name) {
    Genes = .[[df_name]][['Gene']] %>% unique
    df_name = gsub('\\.', '_', df_name) %>% unique
    write.table(Genes, file=paste0('Genes_10_', df_name, '.csv'), row.names=FALSE, col.names=FALSE, quote=FALSE)
    IDs = Cyclers[Cyclers$Gene %in% Genes, 'ENSEMBL'] %>% unique
    #write.table(IDs, file=paste0('IDs_10_', df_name, '.csv'), row.names=FALSE, col.names=FALSE, quote=FALSE)
  })}
} # GENES 2x 0.10

{to_filter %>%
    subset(second_pct > 0.1) %>%
    {split(., list(.$Sample, .$Mode))} %>%
    {lapply(names(.), function(df_name) {
      Genes = .[[df_name]][['Gene']] %>% unique
      df_name = gsub('\\.', '_', df_name) %>% unique
      write.table(Genes, file=paste0('Background_Genes_10_', df_name, '.csv'), row.names=FALSE, col.names=FALSE, quote=FALSE)
      IDs = passQC[passQC$Gene %in% Genes, 'ENSEMBL'] %>% unique
      #write.table(IDs, file=paste0('Background_IDs_10_', df_name, '.csv'), row.names=FALSE, col.names=FALSE, quote=FALSE)
    })}
} # BACKGROUND 2x 0.10

{to_filter %>%
  subset(second_pct > 0.1 & ratio > 1.5) %>% 
  subset(KW_scores_norm < 0.05) %>%     
  mutate(FDR_meta2d = p.adjust(meta2d_pvalue, method='BH')) %>% 
  subset(FDR_meta2d < 0.1) %>% 
    {split(., list(.$Sample, .$Mode))} %>%
    {lapply(names(.), function(df_name) {
      Genes = .[[df_name]][['Gene']] %>% unique
      df_name = gsub('\\.', '_', df_name) %>% unique
      #write.table(Genes, file=paste0('Genes_KWp005_2Dq01_', df_name, '.csv'), row.names=FALSE, col.names=FALSE, quote=FALSE)
      IDs = Cyclers[Cyclers$Gene %in% Genes, 'ENSEMBL'] %>% unique
      write.table(IDs, file=paste0('IDs_KWp005_2Dq01_', df_name, '.csv'), row.names=FALSE, col.names=FALSE, quote=FALSE)
    })}
} # STRINGENT 2x 0.10

to_filter %>%
  subset(second_pct > 0.1 & KW_scores_norm < 0.05 & meta2d_pvalue_AGG_norm < 0.05) %>%
  {split(., list(.$Sample, .$Mode))} %>%
  {lapply(names(.), function(df_name) {
    Genes = .[[df_name]][['Gene']] %>% unique
    df_name = gsub('\\.', '_', df_name) %>% unique
    IDs = Cyclers[Gene %in% Genes, 'ENSEMBL'] %>% unique
    write.table(IDs, file=paste0('myFilt10', df_name, '.csv'), row.names=FALSE, col.names=FALSE, quote=FALSE)
  })}

to_filter %>%
  subset(second_pct > 0.1 & ratio_AGG>1.5 & ampl_AGG>0.5 & meta2d_pvalue_AGG_norm < 0.05) %>%
  {split(., list(.$Sample, .$Mode))} %>%
  {lapply(names(.), function(df_name) {
    Genes = .[[df_name]][['Gene']] %>% unique
    df_name = gsub('\\.', '_', df_name) %>% unique
    IDs = select(org.Dm.eg.db, keys=Genes, columns='ENSEMBL', keytype='SYMBOL')[['ENSEMBL']] %>% unique
    write.table(IDs, file=paste0('myFilt_r1.5_a0.5', df_name, '.csv'), row.names=FALSE, col.names=FALSE, quote=FALSE)
  })}


to_filter %>%
  subset(second_pct > 0.10 & ((ratio > 1.5 & ampl > 0.5 & meta2d_pvalue < 0.05 & 
                                 ratio_AGG > 1.5 & ampl_AGG > 0.5 & meta2d_pvalue_AGG_norm < 0.05) | 
                                (meta2d_pvalue_AGG_norm < 0.05 & ratio_AGG > 2 & ampl_AGG > 0.5))) %>%
  {split(., list(.$Sample, .$Mode))} %>%
  {lapply(names(.), function(df_name) {
    Genes = .[[df_name]][['Gene']] %>% unique
    df_name = gsub('\\.', '_', df_name) %>% unique
    IDs = select(org.Dm.eg.db, keys=Genes, columns='ENSEMBL', keytype='SYMBOL')[['ENSEMBL']] %>% unique
    write.table(IDs, file=paste0('reconstituted_a0.5', df_name, '.csv'), row.names=FALSE, col.names=FALSE, quote=FALSE)
  })}


####    ####    ####    ####    ####    ####    ####    ####    ####    ####   


# 23. PA Analysis
################################################################################

res10_anno_comb_Q_2$C3.LD$results$TF %>% unique
res10$C3_LD$results %>% subset(TF == 'sr') %>% .$Target %>%
  {subset(avg_data_times2_01, Gene %in% ., select = c(Gene, Sample, Mode, timepoint, avg))} %>%
  pivot_wider(names_from=timepoint, values_from=avg) %>%
  {split(., list(.$Sample, .$Mode))} %>%
  lapply(function(obj) {
    S = obj[1, 'Sample']
    M = obj[1, 'Mode']
    subset(obj, select=-c(Sample, Mode)) %>%
      as.data.frame %>%
      {`rownames<-`(., .$Gene)} %>%
      mutate(Gene = NULL) %>%
      as.matrix %>% 
      dist(method='euclidean') %>%
      hclust(method='average') %>%
      cutree(hc, k=4) %>% 
      {data.frame(Gene = names(.), cluster = .)} %>%
      {split(., .$cluster)} %>%
      lapply(function(df) {
        plotter_geneCombine(df$Gene, data=data_C2C3_AGG_01_1d, Mode=M, Sample=S, 
                            NO_X_TEXT=TRUE, NO_GRIDLINES = TRUE, PRINT_GENES = TRUE, LEGEND=TRUE, offset = -0.5, ALIGN=TRUE) 
      })
  }) %>%
  unlist(recursive=FALSE) %>%
  plot_grid(plotlist=., nrow=4)

gene_hclusterer = function(G, S, M, cutoff, K_NOT_H_CUT=TRUE, gene_expr=avg_data_times2_01, times=seq(3,47,4)) {
  G=unique(G)
  hc = subset(gene_expr, Gene %in% G & Sample == S & Mode == M, select=c(Gene, avg, timepoint)) %>% 
    pivot_wider(names_from=timepoint, values_from=avg) %>% 
    as.data.frame %>%
    {`rownames<-`(., .$Gene)} %>% 
    mutate(Gene = NULL) %>%
    as.matrix %>% 
    dist(method='euclidean') %>%
    hclust(method='average')
  res_df = hc %>%
    {if (K_NOT_H_CUT) cutree(., k=cutoff) else cutree(., h=cutoff)} %>%
    {data.frame(Gene = names(.), cluster = .)}
  plots = res_df %>% 
    {split(., .$cluster)} %>%
    lapply(function(df) {
      plotter_geneCombine(df$Gene, data=data_C2C3_AGG_01_1d, Mode=M, Sample=S, 
                          NO_X_TEXT=TRUE, NO_GRIDLINES = TRUE, PRINT_GENES = TRUE, LEGEND=FALSE, offset = -0.5, ALIGN=TRUE) 
    })
  return(list(res_df=res_df, hc=hc, p=plot_grid(plotlist=plots, nrow=1)))
}

sqz <- gene_hclusterer(subset(res10$C3_LD$results, TF == 'sqz')$Target, S='C3', M='LD', cutoff=4)
pho <- gene_hclusterer(subset(res10$C3_DD$results, TF == 'pho')$Target, S='C3', M='DD', cutoff=4)
sr <- gene_hclusterer(subset(res10$C3_LD$results, TF == 'sr')$Target, S='C3', M='LD', cutoff=4)
BtbVII <- gene_hclusterer(subset(res10$C3_LD$results, TF == 'BtbVII')$Target, S='C3', M='LD', cutoff=4)
CrebA <- gene_hclusterer(subset(res10$C3_LD$results, TF == 'CrebA')$Target, S='C3', M='LD', cutoff=4)

res10 <- lapply(res10, function(res) {
  Sa = res$filtered[1, 'Sample']
  Mo = res$filtered[1, 'Mode']
  
  paste("- - - -")
  print(paste(Sa, Mo))
  
  res = c(res, hclustered = list(subset(res$results, TF %in% subset(to_filter, second_pct>0.1 & Sample==Sa & Mode==Mo)$Gene) %>%
    {split(., .$TF)} %>% 
    lapply(function(df) ifelse(length(unique(df$Target)) < 5, 
                               NA, 
                               list(unique(df$Target)))) %>% 
    {.[!is.na(.)]} %>% 
    lapply(function(Ge) gene_hclusterer(G = Ge[[1]], S=Sa, M=Mo, cutoff=4, K_NOT_H_CUT=TRUE)))) %>%
    return
})


a_res10 <- lapply(res10, function(res) {
  res$hclustered %>%
    lapply(function(df_TF) {
      df_TF$res_df %>%
        group_by(cluster) %>%
        filter(n()>=4) %>% ungroup %>%
        {split(., .$cluster)} %>%
        lapply(function(df) list(targets = df$Gene))
    }) %>% 
    Filter(f=length)
})

a_res10 <- lapply(names(a_res10), function(name_L2) {
  Sa = gsub('_.*', '', name_L2)
  Mo = gsub('.*_', '', name_L2)
  lapply(names(a_res10[[name_L2]]), function(name_L3) {
    
    if (name_L3 != 'df') {
      L3 = a_res10[[name_L2]][[name_L3]]
      #print(names(L3))
      
      lapply(names(L3), function(name_L4) {
        L4 = L3[[name_L4]]
        
        print(length(L4$targets))
        
        phases = subset(to_filter, Gene %in% L4$targets & Sample==Sa & Mode==Mo)$meta2d_phase_AGG_norm %>%
          floor %>% unique %>% sort
        df_cyc = subset(reconstituting_DFs_10$reconstituted, Sample==Sa & Mode==Mo, select=c(Gene, meta2d_phase_AGG_norm))
        corrs = make_corrs(L4$targets, Sa, Mo)
        p_fisher = run_fisher_general(A = subset(df_cyc, floor(meta2d_phase_AGG_norm) %in% phases)$Gene,
                                      B = unique(subset(res10[[name_L2]]$results, TF==name_L3)$Target),
                                      U = unique(df_cyc$Gene))
        r_80 = subset(corrs, x_gene != y_gene) %>%
          subset(r > 0.8) %>%
          count(x_gene) %>%
          arrange(desc(n))
        
        S = gsub('_.*', '', name_L2)
        M = gsub('.*_', '', name_L2)
        
        amplitude = subset(avg_data_times2_01, Sample==S & Mode==M & Gene %in% L4$targets) %>%
          group_by(timepoint) %>% 
          summarise(avg = mean(avg)) %>% 
          {(max(.$avg)-min(.$avg))/2}
        
        L4 = list(df = data.frame(TF = name_L3, 
                                  cluster = name_L4,
                                  n_targets = length(L4$targets),
                                  amplitude = amplitude,
                                  r_mean = from_corrplot_meanmaker(corrs),
                                  p_fisher = p_fisher, 
                                  r_80_gene = r_80[1, 'x_gene'],
                                  r_80_n = r_80[1, 'n'] %>% 
                                    {ifelse(is.na(.), 0, .)}, 
                                  cond = name_L2, 
                                  n_phases = length(phases)),
                  cond = name_L2,
                  TF = name_L3, 
                  cluster = name_L4,
                  targets = L4$targets,
                  phases = phases,
                  corrs = corrs,
                  r_80 = r_80)
      }) %>% `names<-`(names(a_res10[[name_L2]][[name_L3]])) 
    }
    }) %>% `names<-`(names(a_res10[[name_L2]]))
  }) %>% `names<-`(names(a_res10)) %>%
  lapply(function(L2) L2 = list(L2, df=list(
    lapply(L2[names(L2) != 'df'], function(L3) lapply(L3, function(L4) L4$df) %>%
             do.call(what='rbind')) %>%
      do.call(what='rbind') 
  )) %>%
    unlist(recursive=FALSE) 
  )

# DO ANY CYCLING TFs PASS FISHER

filt_res_df %>% mutate(Sample = gsub('_.*', '', cond), 
                       Mode = gsub('.*_', '', cond)) %>%
  merge(subset(reconstituting_DFs_10$reconstituted, select=c(Gene, Sample, Mode)) %>%
          mutate(TF = Gene, .keep='unused'))

lapply(a_res10, function(L2) {
  DF_cyc = L2$df %>%
    mutate(Sample = gsub('_.*', '', cond), 
           Mode = gsub('.*_', '', cond)) %>%
    merge(subset(reconstituting_DFs_10$reconstituted, select=c(Gene, Sample, Mode)) %>%
            mutate(TF = Gene, .keep='unused'))
  print(nrow(DF_cyc))
  DF_cyc %>%
    subset(amplitude > 0.35) %>%
    mutate(FDR_fisher = p.adjust(as.numeric(p_fisher), method='BH')) %>%
    subset(FDR_fisher < 0.1) %>% print %>%
    {if (!nrow(.)) ggplot() + theme_void()
      else 
        apply(., 1, function(df) {
          G = df[['TF']][[1]]
          C = df[['cluster']][[1]]
          cond = df[['cond']][[1]]
          S = gsub('_.*', '', cond)
          M = gsub('.*_', '', cond)
          targets = a_res10[[cond]][[G]][[C]]$targets
          
          print(paste(G, C))
          a_res10[[cond]][[G]][[C]]$phases %>% print
          
          make_pdf_targets(subset(res10[[cond]]$results, TF == G), res10[[cond]]$filtered, namebase=paste('TF', G, paste0(S, M), sep='_'), Cols=color_list_alt)
          
          make_pdf_targets(subset(res10[[cond]]$results, TF == G & Target %in% targets), res10[[cond]]$filtered, 
                           namebase=paste('TF', G, paste0(S, M), 'clst', C, sep='_'), Cols=color_list_alt)
          
          plot_grid(nrow=1, ncol=2,
                    avg_plotter_concat(gene=G, data=subset(data_C2C3_AGG, sample==S[[1]] & mode1==M[[1]]), N_days = 2, NO_X_TEXT = TRUE, BLACKOUT = TRUE),
                    subset(res10[[cond]]$hclustered[[G]]$res_df, cluster == C)$Gene %>%
                      plotter_geneCombine(data=data_C2C3_AGG_01_1d, Mode=M, Sample=S, 
                                          NO_X_TEXT=TRUE, NO_GRIDLINES = TRUE, PRINT_GENES = TRUE, LEGEND=FALSE, offset = -0.5, ALIGN=TRUE)
          )
        }) %>%
        plot_grid(plotlist=., nrow=5, ncol=1)
    }
}) %>%
  plot_grid(plotlist=., nrow=1, ncol=4)

res10$C3_LD$hclustered$twi$p

res10$C3_LD$hclustered$pho$p
res10$C3_DD$hclustered$pho$p

res10$C3_DD$hclustered$ab$p
res10$C3_DD$hclustered$ab$res_df %>% arrange(cluster)


lapply(a_res10, function(L2) {
  L2$df %>%
    mutate(Sample = gsub('_.*', '', cond), 
           Mode = gsub('.*_', '', cond)) %>%
    merge(subset(reconstituting_DFs_10$reconstituted, select=c(Gene, Sample, Mode)) %>%
            mutate(TF = Gene, .keep='unused')) %>%
    #subset(amplitude > 0.35) %>%
    mutate(FDR_fisher = p.adjust(as.numeric(p_fisher), method='BH')) %>%
    subset(FDR_fisher < 0.05)
})

a_res10_df <- lapply(a_res10, function(L2) {
  L2 = L2[names(L2) != 'df']
  
  lapply(L2, function(L3) {
    lapply(L3, function(L4) {
      data.frame(cond=L4$cond,
                 TF=L4$TF,
                 cluster=L4$cluster,
                 targets=paste(L4$targets, collapse=','),
                 phases=paste(L4$phases, collapse=','))
    }) %>% do.call(what=rbind)
  }) %>% do.call(what=rbind)
}) %>% do.call(what=rbind) %>% 
  merge(
    lapply(a_res10, function(L2) {
      subset(L2$df, select=c(cond, TF, cluster, n_targets, n_phases, amplitude, r_mean, p_fisher)) %>%
        mutate(CYCLE = paste(TF, cond, sep='_') %in% (reconstituting_DFs_10$reconstituted %>%
                                                        {paste(.$Gene, .$Sample, .$Mode, sep='_')}) ) %>%
        {left_join(., 
                   subset(., CYCLE) %>%
                     mutate(FDR_fisher_cycle = p.adjust(as.numeric(p_fisher), method='BH')) %>%
                     mutate(PASS = FDR_fisher_cycle < 0.05) 
        )}
    }) %>%
      do.call(what=rbind)
  )




trim_TFclusters = function(df, th_n_targets=5, th_r_80_n=1, th_r_mean=0.5, th_p_fisher=0.05, th_FDR_fisher=NULL, PRINT=TRUE) {
  th_n_targets = as.numeric(th_n_targets)
  th_r_80_n = as.numeric(th_r_80_n)
  th_r_mean = as.numeric(th_r_mean)
  th_p_fisher = as.numeric(th_p_fisher)
  th_FDR_fisher = as.numeric(th_FDR_fisher)
  
  if (PRINT) cat(paste('Init:', nrow(df), '\n'))
  
  df = subset(df, n_targets >= th_n_targets)
  if (PRINT) cat(paste0('After n_targets >= ', th_n_targets, ': ', nrow(df), '\n'))

  df = subset(df, r_80_n >= th_r_80_n)
  if (PRINT) cat(paste0('After r_80_n >= ', th_r_80_n, ': ', nrow(df), '\n'))
  
  df = subset(df, r_mean >= th_r_mean)
  if (PRINT) cat(paste0('After r_mean >= ', th_r_mean, ': ', nrow(df), '\n'))
  
  
  if (is.null(th_FDR_fisher)) df = subset(df, p_fisher < th_p_fisher)
  if (PRINT) cat(paste0('pval_fisher <= ', th_FDR_fisher, ': ', nrow(
    subset(df, p_fisher < th_FDR_fisher)), '\n'
  ))
  
  df = mutate(df, FDR_fisher = p.adjust(as.numeric(p_fisher), method='BH')) %>%
    subset(FDR_fisher <= th_FDR_fisher)
  if (PRINT) cat(paste0('FDR_fisher <= ', th_FDR_fisher, ': ', nrow(df), '\n'))
  
  return(df)
}


lapply(a_res10, function(L2) L2$df) %>%
  do.call(what='rbind') %>%
  pivot_longer(cols = -c(TF, cond, r_80_gene, cluster), names_to='filter') %>%
  ggplot(aes(value)) +
  geom_density() +
  scale_x_continuous(trans='log10') +
  facet_wrap(~filter, scales = 'free')

lapply(a_res10, function(L2) {
  L2$df
}) %>%
  do.call(what='rbind') %>%
  mutate(PASS = p_fisher<0.05) %>%
  ggplot(aes(amplitude, n_phases, color=cond, shape=PASS)) +
  geom_point() +
  facet_wrap(~cond)

lapply(a_res10, function(L2) {
  L2$df
}) %>%
  do.call(what='rbind') %>%
  mutate(PASS = p_fisher<0.05) %>%
  ggplot(aes(n_phases, p_fisher, color=cond)) +
  geom_point() +
  facet_wrap(~cond)

filt_res_df <- lapply(names(a_res10), function(name_L2) {
  L2 = a_res10[[name_L2]]
  print(name_L2)
  trim_TFclusters(L2$df,
                  th_n_targets = 10,
                  th_r_80_n = 0, 
                  th_r_mean = 0.5, 
                  th_FDR_fisher=0.1)
}) %>%
  do.call(what='rbind')

filt_res_df_2 <- lapply(names(a_res10), function(name_L2) {
  L2 = a_res10[[name_L2]]
  print(name_L2)
  df = subset(L2$df, amplitude>0.4) %>%
    trim_TFclusters(th_n_targets = 0, th_r_80_n = 0, th_r_mean=0.5, th_FDR_fisher=0.1, PRINT=TRUE)
  print(nrow(df))
  return(df)
}) %>%
  do.call(what='rbind')
table(filt_res_df_2$cond)

filt_res_df %>%
  {split(., .$cond)} %>%
  lapply(function(DF) {
    apply(DF, 1, function(df) {
      TF = df[['TF']]
      C = df[['cluster']]
      cond = df[['cond']]
      S = gsub('_.*', '', cond)
      M = gsub('.*_', '', cond)
      subset(res10[[cond]]$hclustered[[TF]]$res_df, cluster == C)$Gene %>%
        plotter_geneCombine(data=data_C2C3_AGG_01_1d, Mode=M, Sample=S, 
                            NO_X_TEXT=TRUE, NO_GRIDLINES = TRUE, PRINT_GENES = TRUE, LEGEND=FALSE, offset = -0.5, ALIGN=TRUE) 
    }) %>%
      plot_grid(plotlist=., nrow=1, ncol=8)
  }) %>%
  plot_grid(plotlist=., nrow=4, ncol=1) 



expand.grid(Amplitude = seq(0, 0.5, 0.05), 
            Th_r_mean = c(seq(0, 0.6, 0.1)), 
            cond = names(a_res10)) %>% 
  apply(1, function(df) {
    DF = a_res10[[df[['cond']]]]$df %>%
      subset(amplitude > df[['Amplitude']] & r_mean>df[['Th_r_mean']])
    df = cbind(df, N = nrow(trim_TFclusters(df = DF, 
                                            th_r_mean = 0,
                                            th_n_targets = 0,
                                            th_r_80_n=0,
                                            th_FDR_fisher=0.1,
                                            PRINT=TRUE)))
    #df = cbind(df, N = nrow(DF))
  }) %>%
  t %>%
  as.data.frame %>% 
  .[1:4] %>%
  `colnames<-`(c('amplitude', 'th_r_mean', 'cond', 'N')) %>%
  mutate_at(c('amplitude', 'th_r_mean', 'N'), as.numeric) %>% 
  {ggplot(., aes(x=amplitude, y=th_r_mean, color=cond, size=N)) +
      geom_point() +
      scale_size(limits=c(1, NA)) +
      facet_wrap(~cond, ncol=2)}


expand.grid(Amplitude = seq(0, 0.5, 0.05), 
            N_phases = 0:16, 
            cond = names(a_res10)) %>% 
  apply(1, function(df) {
    DF = a_res10[[df[['cond']]]]$df %>%
      subset(amplitude > df[['Amplitude']] & n_phases<df[['N_phases']])
     df = cbind(df, N = nrow(trim_TFclusters(df = DF, 
                                             th_r_mean = 0,
                                             th_n_targets = 0,
                                             th_r_80_n=0,
                                             th_FDR_fisher=0.1,
                                             PRINT=TRUE)))
    #df = cbind(df, N = nrow(DF))
  }) %>%
  t %>%
  as.data.frame %>% 
  .[1:4] %>%
  `colnames<-`(c('amplitude', 'n_phases', 'cond', 'N')) %>%
  mutate_at(c('amplitude', 'n_phases', 'N'), as.numeric) %>% 
  {ggplot(., aes(x=amplitude, y=n_phases, color=cond, size=N)) +
      geom_point() +
      scale_size(limits=c(1, NA)) +
      facet_wrap(~cond, ncol=2)}


# th_r_mean and th_n_targets
expand.grid(th_r_mean = c(seq(0, 0.6, 0.1)), 
            th_n_targets = 5:15, 
            cond = names(a_res10)) %>% 
  apply(1, function(df) {
    df = cbind(df, N = nrow(trim_TFclusters(df = a_res10[[df[['cond']]]]$df, 
                                            th_r_mean = df[['th_r_mean']],
                                            th_n_targets = df[['th_n_targets']], 
                                            th_r_80_n=0,
                                            PRINT=FALSE)))
  }) %>%
  t %>%
  as.data.frame %>% 
  .[1:4] %>%
  `colnames<-`(c('th_r_mean', 'th_n_targets', 'cond', 'N')) %>%
  mutate_at(c('th_r_mean', 'th_n_targets', 'N'), as.numeric) %>% 
  {ggplot(., aes(x=th_r_mean, y=th_n_targets, color=cond, size=N)) +
      geom_point() +
      scale_size(limits=c(1, NA)) +
      facet_wrap(~cond, ncol=2)}


# th_r_mean & th_r_80_n
expand.grid(th_r_mean = c(seq(0, 0.6, 0.1)), 
            th_r_80_n = c(0:6), 
            cond = names(a_res10)) %>%
  apply(1, function(df) {
    df = cbind(df, N = nrow(trim_TFclusters(df = a_res10[[df[['cond']]]]$df, 
                                            th_r_mean = df[['th_r_mean']],
                                            th_r_80_n = df[['th_r_80_n']], 
                                            th_n_targets = 5,
                                            PRINT=FALSE)))
  }) %>%
  t %>%
  as.data.frame %>% 
  .[1:4] %>%
  `colnames<-`(c('th_r_mean', 'th_r_80_n', 'cond', 'N')) %>%
  mutate_at(c('th_r_mean', 'th_r_80_n', 'N'), as.numeric) %>%  
  ggplot(aes(x=th_r_mean, y=th_r_80_n, color=cond, size=N)) +
      geom_point() +
      scale_size(limits=c(1, NA)) +
      facet_wrap(~cond, ncol=2)


# th_r_80_n and th_n_targets
expand.grid(th_n_targets = 0:15, 
            th_r_80_n = c(0:6), 
            cond = names(a_res10)) %>%
  apply(1, function(df) {
    df = cbind(df, N = nrow(trim_TFclusters(df = a_res10[[df[['cond']]]]$df, 
                                            th_n_targets = df[['th_n_targets']],
                                            th_r_80_n = df[['th_r_80_n']], 
                                            th_r_mean = 0,
                                            PRINT=FALSE)))
  }) %>%
  t %>%
  as.data.frame %>% 
  .[1:4] %>%
  `colnames<-`(c('th_n_targets', 'th_r_80_n', 'cond', 'N')) %>%
  mutate_at(c('th_n_targets', 'th_r_80_n', 'N'), as.numeric) %>% 
  {ggplot(., aes(x=th_n_targets, y=th_r_80_n, color=cond, size=N)) +
      geom_point() +
      scale_size(limits=c(1, NA)) +
      facet_wrap(~cond, ncol=2)}

lapply(names(a_res10), function(name_L2) {
  L2 = a_res10[[name_L2]]
  nrow(trim_TFclusters(L2$df,
                       th_n_targets = 10,
                       th_r_80_n = 0, 
                       th_r_mean = 0.5))
})

  
  
{hi <- lapply(a_res10, function(L2) L2$df) %>%
    do.call(what='rbind')
categories = c('r_mean', 'p_fisher', 'r_80_n', 'n_targets')
corrs = as.data.frame(matrix(ncol=length(categories), nrow=0))
for (a in 1:(length(categories)-1)) {
  for (b in (a+1):length(categories)) {
    A = categories[a]
    B = categories[b]
    print(A)
    print(B)
    corrs = rbind(corrs, data.frame(A, B, cor(x=hi[[A]], y=hi[[B]])))
  }
}
print(corrs)
rm(categories, corrs, a, A, b, B, hi)
}


# to functionalize - replace all res10


# 1: TF is expressed
# 2: Targets >= 5
# 3: Targets per cluster > 5
# 4: agreement > 80%

# bring make_corrs above code chunk
# bring run_fisher_general above code chunk (or move analysis below)






make_pdf_targets = function(dfRes, dfFiltered, TFs_in=rownames(data_C2C3), namebase=NULL, splitInto=1, Cols=color_list) {
  # dfRes must have TFs, Targets called Gene, dfFiltered must have genes called Gene

  dfRes = subset(dfRes, TF %in% TFs_in, select = c(TF, Target)) %>% 
    rename(Gene = Target) %>%
    merge(dfFiltered) %>% 
    {split(., .$TF)}
  print(length(dfRes))
  
  LD = 'LD' %in% unique(dfRes[[1]]$Mode)
  
  if (length(dfRes) != 0) {
    for (i in 1:splitInto) {
      dfRes_split = dfRes[ (1+floor((i-1)*length(dfRes)/splitInto)):(floor(i*length(dfRes)/splitInto)) ]
      
      p = lapply(dfRes_split, function(df_split) {
        tf = unique(df_split$TF)
        Genes = unique(df_split$Gene) 
        nTargets = length(unique(df_split$Gene))
        print(tf)
        plot_grid(nrow=1, 
                  ggplot()+theme_minimal(),
                  avg_plotter_concat_mode(tf, N_days=2, rib_alpha=0.4, FACET_2D=TRUE, data=data_C2C3_AGG, NO_GRIDLINES = TRUE, BLACKOUT = TRUE), 
                  hist_phase(df_split, binWidth=24/nTargets, YSCALE=TRUE, LD=LD, Colors = Cols),
                  hist_phase(df_split, CIRCULAR=TRUE, YSCALE=FALSE, LD=LD, Colors = Cols), 
                  plotter_geneCombine_modes(Genes, NO_GRIDLINES = TRUE, NO_X_TEXT_ = TRUE),
                  rel_widths = c(0.2, 1.5, 1, 0.4, 1.5), labels=tf)
      }
      ) %>%
        plot_grid(plotlist=., ncol=3)
      if (!is.null(namebase)) ggsave(file=paste0(namebase, '_', i, '.pdf'), plot=p,  height=length(dfRes_split)/3+0.85, width=40)
      else print(p)
    }
  }
  return(p)
}

make_pdf_targets(subset(res10$C3_LD$results, TF == 'fru'), res10$C3_LD$filtered, namebase='TF_fru_C3LD')
make_pdf_targets(subset(res10$C3_LD$results, TF == 'sr'), res10$C3_LD$filtered, namebase='TF_sr_C3LD')

make_pdf_targets(subset(res10$C3_DD$results, TF == 'BtbVII'), res10$C3_DD$filtered, namebase='TF_BtbVII_C3DD')
make_pdf_targets(subset(res10$C3_DD$results, TF == 'pho'), res10$C3_DD$filtered, namebase='TF_pho_C3DD')


plot_grid(nrow=1, 
          ggplot()+theme_minimal(),
          avg_plotter_concat_mode(tf, N_days=2, rib_alpha=0.4, FACET_2D=TRUE, data=data_C2C3_AGG, NO_GRIDLINES = TRUE), 
          hist_phase(df_split, binWidth=24/nTargets, YSCALE=TRUE, LD=LD),
          hist_phase(df_split, CIRCULAR=TRUE, YSCALE=FALSE, LD=LD), 
          plotter_geneCombine_modes(Genes, NO_GRIDLINES = TRUE),
          rel_widths = c(0.2, 1.5, 1, 0.4, 1.5), labels=tf)



global_print_type='svg'

subset(res10$C3_LD$results, TF == 'fru') %>%
  mutate(Gene = Target, .keep='unused') %>%
  merge(res10$C3_LD$filtered) %>%
  hist_phase(CIRCULAR=TRUE, YSCALE=FALSE) %>%
  ggsave(file=paste0('phase_fru_C3LD.', global_print_type), plot=., width=1.3, height=1.3)

subset(res10$C3_LD$results, TF == 'sr') %>%
  mutate(Gene = Target, .keep='unused') %>%
  merge(res10$C3_LD$filtered) %>%
  hist_phase(CIRCULAR=TRUE, YSCALE=FALSE) %>%
  ggsave(file=paste0('phase_sr_C3LD.', global_print_type), plot=., width=1.3, height=1.3)

subset(res10$C3_DD$results, TF == 'pho') %>%
  mutate(Gene = Target, .keep='unused') %>%
  merge(res10$C3_DD$filtered) %>%
  hist_phase(CIRCULAR=TRUE, YSCALE=FALSE, LD=FALSE) %>%
  ggsave(file=paste0('phase_pho_C3DD.', global_print_type), plot=., width=1.3, height=1.3)

subset(res10$C3_DD$results, TF == 'BtbVII') %>%
  mutate(Gene = Target, .keep='unused') %>%
  merge(res10$C3_DD$filtered) %>%
  hist_phase(CIRCULAR=TRUE, YSCALE=FALSE, LD=FALSE) %>%
  ggsave(file=paste0('phase_BtbVII_C3LD.', global_print_type), plot=., width=1.3, height=1.3)


res10$C3_LD$filtered %>%
  hist_phase(CIRCULAR=TRUE, YSCALE=FALSE) %>%
  ggsave(file='phase_C3LD.pdf', plot=., width=1.3, height=1.3)
res10$C3_DD$filtered %>%
  hist_phase(CIRCULAR=TRUE, YSCALE=FALSE, LD=FALSE) %>%
  ggsave(file='phase_C3DD.pdf', plot=., width=1.3, height=1.3)
res10$C2_LD$filtered %>%
  hist_phase(CIRCULAR=TRUE, YSCALE=FALSE) %>%
  ggsave(file='phase_C2LD.pdf', plot=., width=1.3, height=1.3)
res10$C2_DD$filtered %>%
  hist_phase(CIRCULAR=TRUE, YSCALE=FALSE, LD=FALSE) %>%
  ggsave(file='phase_C2DD.pdf', plot=., width=1.3, height=1.3)


lapply(names(res10_QC_Q10), function(name) {
  res = res10_QC_Q10[[name]]
  namebase = paste0(name, '_Q10_cycle')
  make_pdf_targets(dfRes = res$results, dfFiltered = res$filtered, 
                   TFs_in=res$candidates$SYMBOL,
                   namebase=namebase)
  pdf_combine_and_delete(namebase)
})

lapply(names(res10_QC_Q05), function(name) {
  res = res10_QC_Q05[[name]]
  namebase = paste0(name, '_Q05_cycle')
  make_pdf_targets(dfRes = res$results, dfFiltered = res$filtered, 
                   TFs_in=res$candidates$SYMBOL,
                   namebase=namebase)
  pdf_combine_and_delete(namebase)
})


lapply(names(res10_QC_Q10), function(name) {
  print(name)
  res = res10_QC_Q10[[name]]
  
  namebase = paste0(name, '_Q10_passQC10')
  make_pdf_targets(dfRes = res$results, dfFiltered = res$filtered, 
                   TFs_in=subset(to_filter, Sample == res$filters[[1]] & Mode == res$filters[[2]] & Cell_proportion>0.1)$Gene,
                   namebase=namebase, splitInto = 5)
  pdf_combine_and_delete(namebase)
})


lapply(names(res10_anno_comb_Q_2), function(name) {
  res = res10_anno_comb_Q_2[[name]]
  namebase = paste0(name, '_
                    ')
  make_pdf_targets(dfRes = res$results, dfFiltered = res$filtered, 
                   TFs_in=res$candidates$SYMBOL,
                   namebase=namebase)
  pdf_combine_and_delete(namebase)
})

lapply(names(res10), function(name) {
  res = res10[[name]]
  namebase = paste0(name, '_cycle')
  make_pdf_targets(dfRes = res$results, dfFiltered = res$filtered, 
                   TFs_in=res$candidates$SYMBOL,
                   namebase=namebase)
  pdf_combine_and_delete(namebase)
})

lapply(names(res10), function(name) {
  print(name)
  res = res10[[name]] %>% rename(Gene = 'Target')

  namebase = paste0(name, '_passQC10')
  make_pdf_targets(dfRes = res$results, dfFiltered = res$filtered, 
                   TFs_in=subset(to_filter, Sample == gsub('\\..*', '', name) & Mode == gsub('.*\\.', '', name) & Cell_proportion>0.1)$Gene,
                   namebase=namebase, splitInto = 5)
  pdf_combine_and_delete(namebase)
})

lapply(names(res10), function(name) {
  res = res10[[name]]
  namebase = paste0(name, '_inData')
  make_pdf_targets(dfRes = res$results, dfFiltered = res$filtered, 
                   TFs_in=rownames(data_C2C3),
                   namebase=namebase, splitInto = 20)
  pdf_combine_and_delete(namebase)
})


lapply(names(res10), function(name) {
  print(name)
  res = res10[[name]] %>% rename(Gene = 'Target')
  
  namebase = paste0(name, '_passQC10')
  make_pdf_targets(dfRes = res$results, dfFiltered = res$filtered, 
                   TFs_in=subset(to_filter, Sample == gsub('\\..*', '', name) & Mode == gsub('.*\\.', '', name) & Cell_proportion>0.1)$Gene,
                   namebase=namebase, splitInto = 5)
  pdf_combine_and_delete(namebase)
})

lapply(names(res10), function(name) {
  res = res10[[name]]
  namebase = paste0(name, '_inData')
  make_pdf_targets(dfRes = res$results, dfFiltered = res$filtered, 
                   TFs_in=rownames(data_C2C3),
                   namebase=namebase, splitInto = 20)
  pdf_combine_and_delete(namebase)
})


subset(to_filter, second_pct > 0.15 & ratio > 1.5 & meta2d_pvalue < 0.05 & (ratio_AGG < 1.5 | meta2d_pvalue_AGG_norm > 0.05)) %>%
  make_pdf_filtered(name_distinguisher = 'inCell_notAGG', 
                    CELL = FALSE, CELL_NO_ERROR = FALSE, COND = FALSE, COND_NO_ERROR = FALSE, AGG_NORM = TRUE, AGG_CLRNORM = FALSE,
                    BIO_RIB_CELL = FALSE, BIO_RIB_COND = FALSE, ANNOTATE=TRUE)

#print TFs
lapply(res10, function(res) {
  subset(res$filtered, Gene %in% res$TFs) %>%
    make_pdf_filtered(name_distinguisher = 'TFs', 
                      CELL = FALSE, CELL_NO_ERROR = FALSE, COND = FALSE, COND_NO_ERROR = FALSE, AGG_NORM = TRUE, AGG_CLRNORM = FALSE,
                      BIO_RIB_CELL = FALSE, BIO_RIB_COND = FALSE, ANNOTATE=FALSE)
})



####    ####    ####    ####    ####    ####    ####    ####    ####    ####   


# 24. BiSS
################################################################################

BiSS_WD <- './Promoter_analysis_2x10/BiSS_QC10'

pattern_to_gene_helper = Vectorize(function(pattern) {
  pattern = sub('^[^_]*_', '', pattern)
  if (str_length(pattern) >= 2) {
    if (substr(pattern, 1, 2) == 'l_') {
      pattern = pattern %>% sub('_', '(', .) %>% sub('_', ')', .)
    }
  }
  pattern = gsub('_.*', '', pattern)
  return(pattern)
})

load_BiSS_results = function(wd, filetags=c('C3LD', 'C2LD', 'C3DD', 'C2DD')) {
  genes = genes(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
  names(genes) = genes$gene_id
  
  lapply(filetags, function(filetag) {
    list.files(wd) %>% {.[grep(filetag, .)]} %>%
      lapply(function(File) read.table(file.path(wd, File))) %>%
      do.call(what='rbind') %>%
      `colnames<-`(c('pattern', 'coordinates', 'motif_start', 'motif_end', 'strand', 'score', 'p-value', 'q-value', 'motif')) %>%
      mutate(genePattern = pattern_to_gene_helper(pattern),
             ENSEMBL = gsub('.*_FBgn', 'FBgn', pattern) %>% gsub('_.*', '', .)) %>% unique %>%
      {merge(., select(org.Dm.eg.db,
                       keys=.$ENSEMBL, 
                       columns=c("SYMBOL", "ENSEMBL"),
                       keytype="ENSEMBL"), all=TRUE)} %>% unique %>%
      mutate(Target = ifelse(is.na(SYMBOL), genePattern, SYMBOL), 
             geneEnsembl = SYMBOL, 
             SYMBOL = NULL) %>%
      mutate(seqnames = gsub(':.*', '', coordinates), 
             coordinates = gsub('.*:', '', coordinates), .keep='unused') %>%
      mutate(start = as.integer(gsub('-.*', '', coordinates)), 
             end = as.integer(gsub('.*-', '', coordinates))) %>% 
      merge(subset(as.data.frame(promoters(genes, 1000, 0)), select=-strand)) %>%
      rename('start' = 'promoter_start') %>%
      rename('end' = 'promoter_end') %>%
      {merge(., mutate(select(org.Dm.eg.db,
                              keys=.$gene_id, 
                              columns=c("SYMBOL", "ENSEMBL"),
                              keytype="ENSEMBL"),
                       TF=SYMBOL, gene_id=ENSEMBL, .keep='unused'),
             all=TRUE)} %>% unique
  }) %>%
    `names<-`(gsub('\\.', '_', filetags))
}

BiSS = BiSS_2 %>%
  {`names<-`(.,  paste0(substr(names(BiSS_2), 1,2), '_', substr(names(BiSS_2), 3,4)))}


lapply(names(BiSS_2), function(name) cbind(BiSS_2[[name]], Sample_Mode=name)) %>%
  do.call(what='rbind') %>%
  as.data.frame %>%
  .$`p-value` %>% min
qvalue::qvalue(.$`p-value`, lambda=0)
{cbind(., newQ = qvalue::qvalue(.$`p-value`, lambda=0))}



BiSS_df <- load_BiSS_results(BiSS_WD)
#saveRDS(BiSS_df, 'BiSS_df.rds')

BiSS_df <- readRDS('BiSS_df.rds')


load_BiSS_results_2 = function(wd, filetags=c('C3LD', 'C2LD', 'C3DD', 'C2DD')) {
  genes = genes(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
  names(genes) = genes$gene_id
  
  lapply(filetags, function(filetag) {
    a = Sys.time()
    hold = list.files(wd) %>% {.[grep(filetag, .)]} %>% 
      lapply(function(File) read.table(file.path(wd, File))) %>% 
      do.call(what='rbind') %>%
      `colnames<-`(c('pattern', 'coordinates', 'motif_start', 'motif_end', 'strand', 'score', 'p-value', 'q-value', 'motif')) %>%
      mutate(genePattern = pattern_to_gene_helper(pattern),
             ENSEMBL = gsub('.*_FBgn', 'FBgn', pattern) %>% gsub('_.*', '', .)) %>% unique
    print(paste0('1: ', Sys.time()-a))
    a = Sys.time()
    
    hold = merge(hold, select(org.Dm.eg.db,
                              keys=hold$ENSEMBL, 
                              columns=c("SYMBOL", "ENSEMBL"),
                              keytype="ENSEMBL"), all=TRUE) %>% unique
    print(paste0('2: ', Sys.time()-a))
    a = Sys.time()
    
    hold = hold %>%
      mutate(Target = ifelse(is.na(SYMBOL), genePattern, SYMBOL), 
             geneEnsembl = SYMBOL, 
             SYMBOL = NULL) %>%
      mutate(seqnames = gsub(':.*', '', coordinates), 
             coordinates = gsub('.*:', '', coordinates), .keep='unused') %>%
      mutate(start = as.integer(gsub('-.*', '', coordinates)), 
             end = as.integer(gsub('.*-', '', coordinates)))
    print(paste0('3: ', Sys.time()-a))
    a = Sys.time()
    
    hold = merge(hold, subset(as.data.frame(promoters(genes, 1000, 0)), select=-strand))
    
    print(paste0('4: ', Sys.time()-a))
    a = Sys.time()
    
    hold = hold %>% 
      rename('start' = 'promoter_start') %>%
      rename('end' = 'promoter_end')
    
    print(paste0('5: ', Sys.time()-a))
    a = Sys.time()
    
    hold = merge(hold, mutate(select(org.Dm.eg.db,
                                     keys=hold$gene_id, 
                                     columns=c("SYMBOL", "ENSEMBL"),
                                     keytype="ENSEMBL"),
                              TF=SYMBOL, gene_id=ENSEMBL, .keep='unused'),
                 all=TRUE) %>% unique
  }) %>%
    `names<-`(gsub('\\.', '_', filetags))
}


BiSS_2 <- load_BiSS_results(BiSS_WD)
saveRDS(BiSS_2, 'BiSS_2.rds')
BiSS_2 <- readRDS('BiSS_2.rds')


BiSS = BiSS_2 %>%
  lapply(function(df) mutate(df, Target = Gene, .keep='unused'))



lapply(names(res10), function(name) {
  res = res10[[name]]
  BiSS = BiSS_2[[name]]
  make_pdf_targets(dfRes = BiSS, dfFiltered = res$filtered, TFs_in=unique(res$filtered$Gene), namebase=paste0(name, '_BiSS_cycle'))
  print(name)
})

lapply(names(res10), function(name) {
  res = res10[[name]]
  BiSS = BiSS_df[[name]]
  make_pdf_targets(dfRes = BiSS, dfFiltered = res$filtered, TFs_in=unique(res$filtered$Gene), namebase=paste0(name, '_BiSS_cycle'))
  print(name)
})

####    ####    ####    ####    ####    ####    ####    ####    ####    ####   


# 25. Activity regulated genes (ARGs)
#################################
{ARG = list()
ARG$ChR2_dTrpA1 <- c('CG7995', 'CG7218', 'CG15745', 'CG13999', 'CG42856', 'CG11221', 'Atf-2',
                     'CG6201', 'CG17734', 'CG10960', 'CG1607', 'sty', 'CG6051', 'Eip75B', 
                     'CG8177', 'Gap1', 'cv-c', 'raw', 'BRWD3', 'Hr78', 'CG7510', 'Rph', 'qtc')
ARG$dTrpA1_KCL <- c('Jra')
ARG$ChR2_KCL <- c('Hsromega', 'Dref', 'CG42261', 'Hop', 'CG4577', 'Su(z)2', 
                  'Hsp23', 'CG9328', 'tipE', 'Sln', 'phyl', 'Ubi-p63E')
ARG$ALL <- c('Hr38', 'sr', 'CG14186', 'CG30497', 'CG13055', 'l(1)G0148', 
             'CG17778', 'CG8910', 'CG14024', 'CG42708', 'Mctp','CG13868')
ARG$min2 <- c(ARG$ChR2_dTrpA1, ARG$dTrpA1_KCL, ARG$ChR2_KCL, ARG$ALL)

ARG$converting$converter <- data.frame(paper=ARG$min2, hold=ifelse(ARG$min2 %!in% rownames(data_C2C3), ARG$min2, NA)) %>%
  merge(data.frame(c("CG7995", "Gk2"), 
                   c("CG15745", "IP3K2"), 
                   c("CG42856", "Sik3"), 
                   c("CG11221", "meng"), 
                   c("Gap1", "RasGAP1"), 
                   c("CG7510", "anchor"), 
                   c("Hsromega", "lncRNA:Hsromega"), 
                   c("CG42261", "fid"), 
                   c("Hop", "Stip1"), 
                   c("CG30497", "CG46385"), 
                   c("l(1)G0148", "Cdc7"),
                   c("CG42708", "GLS")) %>% t %>%
          `colnames<-`(c("hold", 'data')), all=TRUE) %>%
  mutate(data = ifelse(is.na(data), paper, data), hold=NULL)

ARG$converted$min2 <- ARG$converting$converter$data
}


ARGs_check <- readxl::read_excel(file.path(getwd(), 'Gene families/AEG/Joined.xlsx'), na = 'NA') %>%
  `colnames<-`(c('ChR2', 'dTrpA1', 'KCl', 'del', 'ChR2+dTrpA1', 'dTrpA1+KCl', 'ChR2+KCl', 'All')) %>%
  subset(select=-del)

ARGs_check_min2 <- c(ARGs_check$ChR2, ARGs_check$dTrpA1, ARGs_check$KCl) %>%
  na.omit %>% 
  table %>% 
  {.[. >= 2]}

rm(ARGs_check, ARGs_check_min2)

# SAVE ARGs
lapply(levels(data_C2C3$sample_mode), function(S_M) {
  S = gsub('_.*', '', S_M)
  M = gsub('.*_', '', S_M)
  subset(to_filter, Sample==S & Mode==M & Gene %in% ARG$converted$min2 & second_pct>0.1)$Gene %>% print %>%
    unique %>% print %>%
  write.csv(file = paste0('ARGs_', S, M, '.csv'), row.names=FALSE, col.names=FALSE)
})


# day vs night
avg_data_times2_max1 %>%
  subset(Gene %in% ARG$converting$converter$data) %>%
  subset(Gene %in% passQC$Gene) %>%
  mutate(timepoint = as.numeric(timepoint)) %>%
  {split(., list(.$Sample, .$Mode))} %>%
  lapply(function(df) {
    df = mutate(df, day = ifelse(floor(timepoint/12)%%2, 'Night', 'Day'), 
                phase = floor(timepoint/12)) 
    #ggboxplot(df, x='phase', y='avg') +
    ggboxplot(df, x='timepoint', y='avg') +
      #ggboxplot(df, x='day', y='avg') +
      stat_anova_test()
  }) %>%
  rev %>%
  {plot_grid(plotlist=., labels = names(.))}

# BOXPLOT OR BOOTPLOT OF ARGs
avg_data_times2_01 %>%
  subset(Gene %in% ARG$converting$converter$data) %>%
  subset(Gene %in% passQC$Gene) %>%
  mutate(avg = avg-0.5) %>%
  mutate(timepoint = as.numeric(timepoint)) %>%
  {split(., list(.$Sample, .$Mode))} %>%
  lapply(function(df) {
    #ggboxplot(df, x='timepoint', y='avg') +
    #  stat_anova_test()
    df = mutate(df, day = ifelse(floor(timepoint/12)%%2, 'Night', 'Day'), 
                phase = floor(timepoint/12)) 
    p = ggplot(df, aes(x=timepoint, y=avg, group=timepoint)) +
      geom_rect(data=data.frame(xmin=0, xmax=12, ymin=-Inf, ymax=Inf), 
                aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill=ifelse(df[1, 'Mode'] == 'DD', "grey65", "yellow"), colour='black', linetype='dashed', linewidth=0.1, alpha=0.3, inherit.aes = FALSE) +
      geom_rect(data=data.frame(xmin=12, xmax=24, ymin=-Inf, ymax=Inf), 
                aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey35", colour='black', linetype='dashed', linewidth=0.1, alpha=0.3, inherit.aes = FALSE) +
      geom_rect(data=data.frame(xmin=24, xmax=36, ymin=-Inf, ymax=Inf), 
                aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill=ifelse(df[1, 'Mode'] == 'DD', "grey65", "yellow"), colour='black', linetype='dashed', linewidth=0.1, alpha=0.3, inherit.aes = FALSE) +
      geom_rect(data=data.frame(xmin=36, xmax=48, ymin=-Inf, ymax=Inf), 
                aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey35", colour='black', linetype='dashed', linewidth=0.1, alpha=0.3, inherit.aes = FALSE) +
      #stat_summary(fun.data = "mean_cl_boot", geom='errorbar') + #stat_summary(geom = "line", fun = "mean") + 
      stat_summary(fun.data = mean_se, geom = "errorbar", width=1) + stat_summary(geom = 'point', fun = "mean") +
      #geom_boxplot(width=2) +
      coord_cartesian(xlim=c(1,23+24), ylim=c(-0.5, 0.5)) +
      theme(axis.title = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), 
            panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      facet_wrap(~Sample)
    if(TRUE) {
      p = p + annotate("text", x = 24, y=0.45, size=2, label = paste('p =', signif(kruskal.test(data=df, avg ~ timepoint)$p.value, digits=3)))
    }
  }) %>% 
  rev %>% 
  #{lapply(names(.), function(name) {
  #  lapply(c('.svg', '.pdf', '.jpg'), function(extension) ggsave(paste0('ARG_boot_', name, extension), plot=.[[name]], width=2.3,height=2))
  #})}
  {plot_grid(plotlist=.#, labels = names(.)
  )} %>%
  {lapply(c('.svg', '.pdf', '.jpg'), function(extension) ggsave(paste0('ARG_boot', extension), plot=., width=5.5,height=5))}

# COMBINED PLOT OF CYCLING ARGs
avg_data_times2_01 %>%
  subset(Gene %in% ARG$converting$converter$data) %>% .$Gene %>%
  {subset(reconstituting_DFs_10$reconstituted, Gene %in% .)} %>%
  {split(., list(.$Sample, .$Mode))} %>%
  lapply(function(x) plotter_geneCombine(x$Gene, data=data_C2C3_AGG_01_1d, Mode=x[1,'Mode'], Sample=x[1,'Sample'], 
                                         NO_X_TEXT=TRUE, NO_GRIDLINES = TRUE, PRINT_GENES = TRUE, LEGEND=TRUE, offset = -0.5, ALIGN=TRUE)) %>%
  rev %>%
  {lapply(names(.), function(name) {
    lapply(c('.svg', '.pdf', '.jpg'), function(extension) ggsave(paste0('ARG_plots_', name, extension), plot=.[[name]], width=4.2,height=3))
  })}
#plot_grid(plotlist=.) %>%
#{lapply(c('.svg', '.pdf', '.jpg'), function(extension) ggsave(paste0('ARG_plots', extension), plot=., width=8,height=5))}

# INDIVIDUAL GENE DAY/NIGHT
{setwd('/Users/Teddy/Desktop/CCG_distribution/ARG plots/DayNight')
  avg_data_times2_01 %>%
    subset(Gene %in% ARG$converting$converter$data) %>%
    subset(Gene %in% passQC$Gene) %>%
    #subset(Gene %in% c('sr', 'Hr38', 'CG14186')) %>%
    mutate(avg = avg-0.5) %>%
    mutate(timepoint = as.numeric(timepoint)) %>%
    {split(., list(.$Sample, .$Mode))} %>%
    lapply(function(df) {
      #ggboxplot(df, x='timepoint', y='avg') +
      #  stat_anova_test()
      df = mutate(df, day = ifelse(floor(timepoint/12)%%2, 'Night', 'Day'), 
                  phase = floor(timepoint/12)) 
      
      print(df)
      
      bartop = max(plotrix::std.error(subset(df, day=='Day')$avg)+mean(subset(df, day=='Day')$avg), 
                   plotrix::std.error(subset(df, day=='Night')$avg)+mean(subset(df, day=='Night')$avg))
      
      p = ggplot(df, aes(x=day, y=avg, group=day)) +
        geom_rect(data=data.frame(xmin=-Inf, xmax=1.5, ymin=-Inf, ymax=Inf), 
                  aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill=ifelse(df[1, 'Mode'] == 'DD', "grey65", "yellow"), colour='black', linetype='dashed', linewidth=0.1, alpha=0.3, inherit.aes = FALSE) +
        geom_rect(data=data.frame(xmin=1.5, xmax=Inf, ymin=-Inf, ymax=Inf), 
                  aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="grey35", colour='black', linetype='dashed', linewidth=0.1, alpha=0.3, inherit.aes = FALSE) +
        #stat_summary(fun.data = "mean_cl_boot", geom='errorbar') + #stat_summary(geom = "line", fun = "mean") + 
        #stat_summary(fun.data = mean_se, geom = "errorbar") +
        geom_boxplot(width=0.25) +
        theme(axis.title = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), 
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text = element_text(face = "italic")) +
        coord_cartesian(ylim=c(-0.5, 0.5))
      #coord_cartesian(ylim(-0.5, 0.5)) +
      #coord_cartesian(xlim=c(1,23+24))
      
      print(p)
      
      if(TRUE) {
        p = p + add_pvalue(data.frame(label = paste('p =', signif(wilcox.test(data=df, avg ~ day, alternative='greater')$p.value, digits=3)),
                                      group1='Day', group2='Night', y.position =bartop+0.3), label.size=2)
        
      }
    }) %>% 
    {lapply(names(.), function(name) { 
      lapply(c('.jpg', '.svg', '.pdf'), function(extension) ggsave(paste0('DAYNIGHT_', name, '_wilcox_onesided', extension), plot=.[[name]], width=2.3,height=2))
    })}
  setwd(WD)
}



ARG$min2$DayNight$AGGall <- avg_data_times2_AGG %>%
  subset(Gene %in% ARG$converting$converter$data) %>%
  subset(Gene %in% passQC$Gene) %>%
  #subset(Gene %in% c('sr', 'Hr38', 'CG14186', 'RasGAP1')) %>%
  mutate(avg = avg-0.5) %>%
  mutate(timepoint = as.numeric(timepoint)) %>%
  {split(., list(.$Sample, .$Mode, .$Gene))} %>%
  lapply(function(df) {
    #ggboxplot(df, x='timepoint', y='avg') +
    #  stat_anova_test()
    df = mutate(df, day = ifelse(floor(timepoint/12)%%2, 'Night', 'Day'), 
                phase = floor(timepoint/12)) 
    data.frame(Gene = df[1, 'Gene'],
               Sample = df[1, 'Sample'],
               Mode = df[1, 'Mode'],
               p_wilcoxon_dayUp = signif(wilcox.test(data=df, avg ~ day, alternative='greater')$p.value, digits=3), 
               p_wilcoxon = signif(wilcox.test(data=df, avg ~ day)$p.value, digits=3))
    
  }) %>%
  do.call(what='rbind') %>% 
  #anti_join(reconstituting_DFs_10$reconstituted) %>% 
  mutate(FDR = p.adjust(as.numeric(p_wilcoxon_dayUp), method='BH'), 
         FDR_two = p.adjust(as.numeric(p_wilcoxon), method='BH')) %>%
  arrange(p_wilcoxon_dayUp) %>%
  mutate(N = FDR/p_wilcoxon_dayUp)

{split(., list(.$Sample, .$Mode))} %>%
  lapply(function(df) df %>%
           mutate(FDR = p.adjust(p_wilcoxon_dayUp, method='BH'), 
                  FDR_two = p.adjust(p_wilcoxon, method='BH'))
  ) %>% do.call(what='rbind') %>%
  arrange(FDR)

lapply(c('AGGall'), function(name) {
  #lapply(names(ARG$min2$DayNight), function(name) {
  cat('\n- - - - - - - - - - - - - - - - - - \n')
  print(name)
  pass = ARG$min2$DayNight[[name]] %>%
    mutate(PASS = FDR < 0.1) %>%
    mutate(pPass = p_wilcoxon_dayUp < 0.05) %>%
    mutate(PASS_two = FDR_two < 0.1) %>%
    mutate(pPass_two = p_wilcoxon < 0.05)
  
  print(subset(pass, pPass==TRUE) %>% arrange(Mode) %>% arrange(Sample))
  print('FDR one')
  print(subset(pass, PASS==TRUE) %>%
          {table(.$Sample, .$Mode)})
  print('p one')
  print(subset(pass, pPass==TRUE) %>%
          {table(.$Sample, .$Mode)})
  print('FDR two')
  print(subset(pass, PASS_two==TRUE) %>%
          {table(.$Sample, .$Mode)})
  print('p two')
  print(subset(pass, pPass_two==TRUE) %>%
          {table(.$Sample, .$Mode)})
  return('')
})

{setwd('/Users/Teddy/Desktop/CCG_distribution/ARG plots/plots_3Genes')
  for (s in c('C3', 'C2')) {
    for (m in c('LD', 'DD')) {
      plotter_geneCombine(genes_activity, data=data_C2C3_AGG_01_1d, COUNTS=TRUE, Sample=s, Mode=m, 
                          PRINT_GENES=TRUE, offset=-0.5, ALIGN=TRUE, NO_GRIDLINES=TRUE, LEGEND=TRUE, NO_X_TEXT=TRUE, PRINT_MEAN=FALSE, point_size=0) %>%
        {lapply(c('.jpg', '.svg', '.pdf'), function(ext)   ggsave(file=paste0('ARG_comb_', s,m, ext), plot=., width=2.9, height=1.25*1.25))}
      plotter_geneCombine(genes_activity, data=data_C2C3_AGG_01_1d, COUNTS=TRUE, Sample=s, Mode=m, 
                          PRINT_GENES=TRUE, offset=-0.5, ALIGN=TRUE, NO_GRIDLINES=TRUE, LEGEND=FALSE, NO_X_TEXT=TRUE, PRINT_MEAN=FALSE, point_size=0) %>%
        {lapply(c('.jpg', '.svg', '.pdf'), function(ext)   ggsave(file=paste0('ARG_comb_noLegend_', s,m, ext), plot=., width=1.25*1.25+0.25, height=1.25*1.25))}
    }
  } 
  setwd(WD)
  rm(s, m)
}

{setwd('/Users/Teddy/Desktop/CCG_distribution/ARG plots/plots_collective')
  for (s in c('C3', 'C2')) {
    for (m in c('LD', 'DD')) {
      subset(reconstituting_DFs_10$reconstituted, Sample == s & Mode == m & Gene %in% ARG$converted$min2)$Gene %>%
        plotter_geneCombine(data=data_C2C3_AGG_01_1d, COUNTS=TRUE, Sample=s, Mode=m, 
                            PRINT_GENES=TRUE, offset=-0.5, ALIGN=TRUE, NO_GRIDLINES=TRUE, LEGEND=TRUE, NO_X_TEXT=TRUE, PRINT_MEAN=FALSE, point_size=0) %>%
        {lapply(c('.jpg', '.svg', '.pdf'), function(ext)   ggsave(file=paste0('ARG_plots_', s,m, ext), plot=., width=2.9, height=1.25*1.25))}
      subset(reconstituting_DFs_10$reconstituted, Sample == s & Mode == m & Gene %in% ARG$converted$min2)$Gene %>%
        plotter_geneCombine(data=data_C2C3_AGG_01_1d, COUNTS=TRUE, Sample=s, Mode=m, 
                            PRINT_GENES=TRUE, offset=-0.5, ALIGN=TRUE, NO_GRIDLINES=TRUE, LEGEND=FALSE, NO_X_TEXT=TRUE, PRINT_MEAN=FALSE, point_size=0) %>%
        {lapply(c('.jpg', '.svg', '.pdf'), function(ext)   ggsave(file=paste0('ARG_plots_noLegend_', s,m, ext), plot=., width=1.25*1.25+0.25, height=1.25*1.25))}
    }
  } 
  setwd(WD)
  rm(s, m)
}


lapply(c('C3_LD', 'C2_LD', 'C3_DD', 'C2_DD'), function(S_M) {
  S = gsub('_.*', '', S_M)
  M = gsub('.*_', '', S_M)
  
  D = subset(data_C2C3_AGG, sample_mode == S_M)
  
  subset(to_filter, Gene %in% ARG$converted$min2 & Sample==S & Mode==M & second_pct>0.1)$Gene %>% 
    avg_plotter_concat(data = D, N_days=2, NO_X_TEXT = TRUE, NO_GRIDLINES = TRUE)
}) %>%
  plot_grid(nrow=1, plotlist=.) %>% print %>%
  ggsave(file=paste0('expressed_ARGs.pdf'), plot=., width=10, height=40)



####    ####


# 26. GENES OF INTEREST/GENE ENRICHMENT
################################################################################

# TF notes:
lapply(res10, function(res) {
  length(which(unique(res$results$TF) %in% unique(res$filtered$Gene)))
})
lapply(names(res10), function(name) {
  res = res10[[name]]
  length(which(unique(res$results$TF) %in% unique(subset(to_filter, Sample == gsub('_.*', '', name) & Mode == gsub('.*_', '', name) & Cell_proportion>0.1)$Gene)))
})

lapply(names(res10), function(name) {
  res = res10[[name]]
  length(which(unique(res$results$TF) %in% unique(subset(to_filter, Sample == gsub('_.*', '', name) & Mode == gsub('.*_', '', name) & to_filter$max_AGG>0)$Gene)))
})

unique(res10$C3_LD$ResinCandidates$TF) %>% {.[. %in% unique(res10$C3_DD$ResinCandidates$TF)]} # pho and cwo
unique(res10$C3_LD$ResinCandidates$TF) %>% {.[. %in% unique(res10$C2_LD$ResinCandidates$TF)]} # sr and fru

make_panel(data.frame(Gene = ARG$converted$min2), pass_df = subset(to_filter, second_pct<0.1), cyclers=subset(to_filter, second_pct<0.1), PLOT_NESTING=FALSE)

make_panel = function(ID_df, filename=NULL, pass_df=reconstituting_DFs_10$reconstituted, cyclers=passQC, PLOT_NESTING=TRUE) {
  L_p = subset(pass_df, Gene %in% merge(ID_df, cyclers)$Gene, 
               select=c(Gene, Sample, Mode)) %>%
    {split(., list(.$Sample, .$Mode))} %>% 
    lapply(function(df) {
      if (nrow(df)) {
        avg_plotter_concat(gene = unique(df$Gene), COUNTS=TRUE,
                           data=subset(data_C2C3_AGG_max1, sample == df[1,'Sample'] & mode1 == df[1,'Mode']), 
                           N_days=2, SWAP_FACET = TRUE, NO_GRIDLINES=TRUE, NO_X_TEXT=TRUE)
      } else (ggplot()+theme_minimal())
    })
  
  print(names(L_p))
  
  if (PLOT_NESTING) {
    p = plot_grid(ncol=1, L_p[[4]], rel_heights = c(2+nrow(L_p[[3]]$data)+nrow(L_p[[2]]$data)+nrow(L_p[[1]]$data), 2+nrow(L_p[[4]]$data)),
                  plot_grid(nrow=1, L_p[[3]], rel_widths = c(2+nrow(L_p[[3]]$data), nrow(L_p[[2]]$data)+2+nrow(L_p[[1]]$data)),
                            plot_grid(nrow=1, L_p[[2]], L_p[[1]], rel_widths = c(2+nrow(L_p[[2]]$data), 2+nrow(L_p[[1]]$data)))))
    if (!is_null(filename)) ggsave(file=filename, plot=p, width=3*length(unique(L_p[[4]]$data$Gene)), height=4.5, limitsize = FALSE)
  } else {
    p = plot_grid(ncol=1, plotlist=L_p) 
    if (!is_null(filename)) ggsave(file=filename, plot=p, width=2.5*length(unique(L_p[[4]]$data$Gene)), height=6, limitsize = FALSE)
  }
  return(p)
}

make_phase_from_GOdf = function(GO_df, n_W=0, binWidth=2, cyclers=Cyclers, pass_df=reconstituting_DFs_10$reconstituted) {
  merge(merge(GO_df, cyclers), pass_df) %>% 
    {plot_grid(nrow=1, hist_phase(., CIRCULAR=FALSE, TRIM_TO_24=TRUE, binWidth=2, n_W=n_W, YSCALE=TRUE), 
               hist_phase(., CIRCULAR=TRUE, TRIM_TO_24=TRUE, binWidth=2, n_W=n_W), 
               hmap_phase(.))}
}
  


####    ####


# 27. Gene rankings, Pdfr
#############
expr_rankings <- subset(expr_data, select=c(Gene, Sample, Mode, avg, Cell_proportion), Cell_proportion>0.1) %>%
  group_by(Gene, Sample) %>%
  summarise(avg = mean(avg), 
            Cell_proportion = mean(Cell_proportion)) %>%
  ungroup %>%
  mutate(avg_percentile = rank(avg) %>% {./max(.)}, 
         prop_percentile = rank(Cell_proportion) %>% {./max(.)})

subset(expr_data, select=c(Gene, Sample, Mode, avg, Cell_proportion)) %>%
  group_by(Gene, Sample) %>%
  summarise(avg = mean(avg), 
            Cell_proportion = mean(Cell_proportion)) %>%
  ungroup %>%
  mutate(avg_percentile = rank(avg) %>% {./max(.)}, 
         prop_percentile = rank(Cell_proportion) %>% {./max(.)}) %>%
  subset(Gene == 'per')


subset(expr_rankings, Gene == 'Pdfr')
####     ####


# 28. GO GROUPS, MAKING PLOTS
#########################



# MITOCHONDRIAL GENES:
rownames(data_C2C3) %>% {.[startsWith(., 'mt:')]} # length 38
rownames(data_C2C3) %>% {.[startsWith(., 'mt:t')]} %>% # length 22
  avg_plotter_concat_mode(data=data_C2C3_AGG)
rownames(data_C2C3) %>% {.[startsWith(., 'mt:')]} %>% {.[!startsWith(., 'mt:t')]} %>%
  avg_plotter_concat_mode(data=data_C2C3_AGG)

# INITIALIZING AND PRINTING GO TERMS
{GO_TRs <- select(org.Dm.eg.db,
                 keys=c("GO:0003712", GOMFOFFSPRING[["GO:0003712"]],
                        "GO:0003700", GOMFOFFSPRING[["GO:0003700"]]),
                 columns="ENSEMBL",
                 keytype="GO")['ENSEMBL'] %>% unique
make_panel(GO_TRs, "GO_TR_panel.svg")
make_panel(GO_TRs, filename='GO_TRs_panel_expressed.svg', pass_df=subset(to_filter, second_pct > 0.1), PLOT_NESTING=FALSE)
make_phase_from_GOdf(GO_TRs) %>% ggsave(file='GO_TR_hmap.svg', plot=., width=5, height=5)

GO_RBPs <- select(org.Dm.eg.db,
                 keys=c("GO:0003723", GOMFOFFSPRING[["GO:0003723"]]),
                 columns="ENSEMBL",
                 keytype="GO")['ENSEMBL'] %>% unique
make_panel(GO_RBPs, "GO_RBP_panel.svg")
make_panel(GO_RBPs, filename='GO_RBPs_panel_expressed.svg', pass_df=subset(to_filter, second_pct > 0.1), PLOT_NESTING=FALSE)
make_phase_from_GOdf(GO_RBPs) %>% ggsave(file='GO_RBP_hmap.pdf', plot=., width=5, height=5)

GO_translation <- select(org.Dm.eg.db,
                         keys=c("GO:0006412", GOBPOFFSPRING[["GO:0006412"]]),
                         columns="ENSEMBL",
                         keytype="GO")['ENSEMBL'] %>% unique
make_panel(GO_translation, "GO_translation_panel.svg")
make_panel(GO_translation, filename='GO_translation_panel_expressed.svg', pass_df=subset(to_filter, second_pct > 0.1), PLOT_NESTING=FALSE)
make_phase_from_GOdf(GO_translation) %>% ggsave(file='GO_translation_hmap.svg', plot=., width=5, height=5)

GO_trans_init <- select(org.Dm.eg.db,
                        keys=c("GO:0003743", GOMFOFFSPRING[["GO:0003743"]]),
                        columns="ENSEMBL",
                        keytype="GO")['ENSEMBL'] %>% unique
make_panel(GO_trans_init, "GO_trans_init_panel.svg")
make_panel(GO_trans_init, filename='GO_trans_init_panel_expressed.svg', pass_df=subset(to_filter, second_pct > 0.1), PLOT_NESTING=FALSE)
make_phase_from_GOdf(GO_trans_init) %>% ggsave(file='GO_trans_init_hmap.svg', plot=., width=5, height=5)
  
GO_ribo <- AnnotationDbi::select(org.Dm.eg.db,
                                 keys=c("GO:0005840", GOCCOFFSPRING[["GO:0005840"]]),
                                 columns="ENSEMBL",
                                 keytype="GO")['ENSEMBL'] %>% unique

GO_ribo_cell <- select(org.Dm.eg.db,
                       keys=c("GO:0022626", GOCCOFFSPRING[["GO:0022626"]]),
                       columns="ENSEMBL",
                       keytype="GO")['ENSEMBL'] %>% unique
make_panel(GO_ribo_cell, "GO_ribo_cell_panel.svg")
make_panel(GO_ribo_cell, filename='GO_ribo_cell_panel_expressed.svg', pass_df=subset(to_filter, second_pct > 0.1), PLOT_NESTING=FALSE)
make_phase_from_GOdf(GO_ribo_cell) %>% ggsave(file='GO_ribo_cell_hmap.svg', plot=., width=5, height=5)

avg_plotter_concat_mode(merge(GO_ribo_cell, passQC)$Gene) %>%
  ggsave(filename='GO_ribo_cell_passQC.pdf', height=80, width=6, limitsize=FALSE)
avg_plotter_concat_mode(merge(GO_ribo_mito, passQC)$Gene) %>%
  ggsave(filename='GO_ribo_mito_passQC.pdf', height=80, width=6, limitsize=FALSE)

A = merge(GO_ribo_cell, passQC)$Gene
B = subset(Filtereds$DFs$reps_adding_pool, Sample == 'C3' & Mode == 'LD')$Gene
U = subset(to_filter, second_pct > 0.1 & Sample == 'C3' & Mode == 'LD')$Gene


data.frame(c(length(unique( A[A %in% B])), 
             length(unique( A[A %!in% B] ))),
           c(length(unique( B[B %!in% A] )), 
             length(unique( U[U %!in% c(A, B)] )))) %>%
  fisher.test()


GO_ribo_mito <- select(org.Dm.eg.db,
                       keys=c("GO:0005761", GOCCOFFSPRING[["GO:0005761"]]),
                       columns=c("ENSEMBL"),
                       keytype="GO")['ENSEMBL'] %>% unique
make_panel(GO_ribo_mito, "GO_ribo_mito_panel.svg")
make_panel(GO_ribo_mito, filename='GO_ribo_mito_panel_expressed.svg', pass_df=subset(to_filter, second_pct > 0.1), PLOT_NESTING=FALSE)
make_phase_from_GOdf(GO_ribo_mito) %>% ggsave(file='GO_ribo_mito_hmap.svg', plot=., width=5, height=5)


GO_mito <- select(org.Dm.eg.db,
                  keys=c("GO:0005739", GOCCOFFSPRING[["GO:0005739"]]),
                  columns=c("ENSEMBL"),
                  keytype="GO")['ENSEMBL'] %>% unique
make_panel(GO_mito, "GO_mito_panel.svg")
make_panel(GO_mito, filename='GO_mito_panel_expressed.svg', pass_df=subset(to_filter, second_pct > 0.1), PLOT_NESTING=FALSE)
make_phase_from_GOdf(GO_mito) %>% ggsave(file='GO_mito_hmap.svg', plot=., width=5, height=5)


GO_mito_notComplex <- select(org.Dm.eg.db,
                  keys=c("GO:0005759", GOCCOFFSPRING[["GO:0005759"]],
                         "GO:0044290", GOCCOFFSPRING[["GO:0044290"]],
                         "GO:0005740", GOCCOFFSPRING[["GO:0005740"]],
                         "GO:0020023", GOCCOFFSPRING[["GO:0020023"]]),
                  columns=c("ENSEMBL"),
                  keytype="GO")['ENSEMBL'] %>% unique
make_panel(GO_mito_notComplex, "GO_mito_notComplex_panel.svg")
make_panel(GO_mito_notComplex, filename='GO_mito_notComplex_panel_expressed.svg', pass_df=subset(to_filter, second_pct > 0.1), PLOT_NESTING=FALSE)
make_phase_from_GOdf(GO_mito_notComplex) %>% ggsave(file='GO_mito_notComplex_hmap.svg', plot=., width=5, height=5)

GO_channel <- select(org.Dm.eg.db,
                     keys=c("GO:0005216", GOBPOFFSPRING[["GO:0005216"]]),
                     columns="ENSEMBL",
                     keytype="GO")['ENSEMBL'] %>% unique
make_panel(GO_channel, "GO_channel_panel.svg")
make_panel(GO_channel, filename='GO_channel_panel_expressed.svg', pass_df=subset(to_filter, second_pct > 0.1), PLOT_NESTING=FALSE)
make_phase_from_GOdf(GO_channel) %>% ggsave(file='GO_channel_hmap.svg', plot=., width=5, height=5)

#pie = list()
pie_genes <- list(ARG = ARG$converted$min2,   GPCR = genes_GPCR, 
                  TRs = GO_TRs %>% merge(Cyclers) %>% .$Gene, 
               Channel = GO_channel %>% merge(Cyclers) %>% .$Gene, 
               cRibo = GO_ribo_cell %>% merge(Cyclers) %>% .$Gene, 
               mRibo = GO_ribo_mito %>% merge(Cyclers) %>% .$Gene,
               ETC_mt = GO_ETC %>% merge(Cyclers) %>% {.[startsWith(.$Gene, 'mt:'), 'Gene']},
               ETC_Nomt = GO_ETC %>% merge(Cyclers) %>% {.[!startsWith(.$Gene, 'mt:'), 'Gene']},
               mt = GO_mito_notComplex %>% merge(Cyclers) %>% {.[startsWith(.$Gene, 'mt:'), 'Gene']},
               Nomt = GO_mito_notComplex %>% merge(Cyclers) %>% {.[!startsWith(.$Gene, 'mt:'), 'Gene']}
)


library(webr)

make_pie_from_categoryGenes = function(categoryGeneList, cycGenes, DROP_OTHER=FALSE, nameReplacements=NULL) {
  plot_genes = categoryGeneList
  if (!is.null(nameReplacements)) names(plot_genes) = nameReplacements
  
  if (DROP_OTHER) Levs = plot_genes
  else Levs = c(names(plot_genes), 'other')
  
  lapply(names(plot_genes), function(catName) {
    catGenes = plot_genes[[catName]]
    data.frame(cycGenes %in% catGenes) %>%
      `names<-`(catName) %>%
      `rownames<-`(cycGenes)
  }) %>% 
    do.call(what='cbind') %>% 
    {apply(., 1, function(L) {
      seq_along(L) == min(which(L))
    })} %>% 
    {rbind(., !colSums(.))} %>%
    data.frame(category = c(names(plot_genes), 'other')) %>% 
    pivot_longer(cols = -category) %>% 
    subset(value==TRUE, select = -value) %>% 
    mutate(category = factor(category, levels=Levs)) %>% 
    PieDonut(aes(category), pieLabelSize = 2)
}

reconstituting_DFs_10$reconstituted %>%
  {split(., list(.$Sample, .$Mode))} %>% 
  {plot_grid(plotlist = lapply(., function(df) {
    make_pie_from_categoryGenes(pie_genes, unique(df$Gene))
  }), 
  labels = names(.))}

filtered_KWp_2Dq %>%
  {split(., list(.$Sample, .$Mode))} %>% 
  {plot_grid(plotlist = lapply(., function(df) {
    make_pie_from_categoryGenes(pie_genes, unique(df$Gene), DROP_OTHER = FALSE)
  }), 
  labels = names(.))}

# AMPK (with regards to mitochondrial cycling - in C2 LD and DD, (lowest mRNA levels) AMPKalpha peaks during the night)
c('AMPKalpha', 'SNF4Agamma', 'alc') %>%
  avg_plotter_concat_mode(data=data_C2C3_AGG)

GO_ETC <- select(org.Dm.eg.db,
                 keys=c("GO:0022900", GOBPOFFSPRING[["GO:0022900"]]),
                 columns="ENSEMBL",
                 keytype="GO")['ENSEMBL'] %>% unique
make_panel(GO_ETC, "GO_ETC_panel.svg")
make_panel(GO_ETC, filename='GO_ETC_panel_expressed.svg', pass_df=subset(to_filter, second_pct > 0.1), PLOT_NESTING=FALSE)
make_phase_from_GOdf(GO_ETC) %>% ggsave(file='GO_ETC_hmap.svg', plot=., width=5, height=5)

GO_ATP_SYN_COUP_ET <- select(org.Dm.eg.db,
                             keys=c("GO:0042775", GOBPOFFSPRING[["GO:0042775"]]),
                             columns="ENSEMBL",
                             keytype="GO")['ENSEMBL'] %>% unique
make_panel(GO_ATP_SYN_COUP_ET, "GO_ATPSYN_coupled_ET.svg")
make_panel(GO_ATP_SYN_COUP_ET, filename='GO_ATP_SYN_COUP_ET_panel_expressed.svg', pass_df=subset(to_filter, second_pct > 0.1), PLOT_NESTING=FALSE)
make_phase_from_GOdf(GO_ATP_SYN_COUP_ET) %>% ggsave(file='GO_ATP_SYN_COUP_ET_hmap.svg', plot=., width=5, height=5)

# write code ^ for single GO printing to work with res10 (ribosomes)
# while running code, manually annotate missing cyclers (or missing in rownames(data_C2C3))


GO_circ <- select(org.Dm.eg.db,
                  keys=c("GO:0042752", GOBPOFFSPRING[["GO:0042752"]]),
                  columns="ENSEMBL",
                  keytype="GO")['ENSEMBL'] %>% unique
make_panel(GO_circ, filename='GO_circ_panel.svg')
make_panel(GO_circ, filename='GO_circ_panel_expressed.svg', pass_df=subset(to_filter, second_pct > 0.1), PLOT_NESTING=FALSE)
make_phase_from_GOdf(GO_circ) %>% ggsave(file='GO_circ_hmap.svg', plot=., width=5, height=5)

GO_EIF3 <- select(org.Dm.eg.db,
                  keys=c("GO:0005852", GOCCOFFSPRING[["GO:0005852"]]),
                  columns="ENSEMBL",
                  keytype="GO")['ENSEMBL'] %>% unique
make_panel(GO_EIF3, filename='GO_EIF3_panel.svg')
make_panel(GO_EIF3, filename='GO_EIF3_panel_expressed.svg', pass_df=subset(to_filter, second_pct > 0.1), PLOT_NESTING=FALSE)
make_phase_from_GOdf(GO_EIF3) %>% ggsave(file='GO_EIF3_hmap.svg', plot=., width=5, height=5)

GO_nitro <- AnnotationDbi::select(org.Dm.eg.db,
                                  keys=c("GO:0044271", GOBPOFFSPRING[["GO:0044271"]]),
                                  columns="ENSEMBL",
                                  keytype="GO")['ENSEMBL'] %>% unique
GO_NPR <- AnnotationDbi::select(org.Dm.eg.db,
                                keys=c("GO:0008188", GOCCOFFSPRING[["GO:0008188"]]),
                                columns="ENSEMBL",
                                keytype="GO")['ENSEMBL'] %>% unique


GOs <- list(GO_TRs=GO_TRs, GO_RBPs=GO_RBPs, GO_translation=GO_translation, GO_EIF3=GO_EIF3, GO_trans_init=GO_trans_init,
            GO_ribo_cell=GO_ribo_cell, GO_ribo_mito=GO_ribo_mito, GO_ribo=GO_ribo, GO_nitro=GO_nitro, GO_NPR=GO_NPR, GO_ETC=GO_ETC, GO_circ=GO_circ)
}






lapply(names(GOs), function(GO_name) {
  Ps = subset(reconstituting_DFs_10$reconstituted, Gene %in% merge(GOs[[GO_name]], Cyclers)$Gene) %>%
  {split(., list(.$Sample, .$Mode))} %>% 
  lapply(function(x) {if(nrow(x)) {
    plotter_geneCombine(x$Gene, data=data_C2C3_AGG_01_1d, Mode=x[1,'Mode'], Sample=x[1,'Sample'], 
                                         NO_X_TEXT=TRUE, NO_GRIDLINES = TRUE, PRINT_GENES = TRUE, LEGEND=TRUE, offset = -0.5, ALIGN=TRUE)
  } else return(NULL)
  })
  
  Ps %>% rev %>%
  #{lapply(names(.), function(name) {
  #  lapply(c('.svg', '.pdf', '.jpg'), function(extension) ggsave(paste0('Comb_', GO_name, name, extension), plot=.[[name]], width=4.2,height=3))
  #})}
    plot_grid(plotlist=.) %>%
    {lapply(c('.svg', '.pdf', '.jpg'), function(extension) ggsave(paste0('Comb_', GO_name, extension), plot=., width=8,height=5))}
})


lapply(names(GOs), function(GO_name) {
  Ps = subset(to_filter, second_pct>0.1 & Gene %in% merge(GOs[[GO_name]], passQC)$Gene) %>%
    {split(., list(.$Sample, .$Mode))} %>% 
    lapply(function(x) {if(nrow(x)) {
      plotter_geneCombine(x$Gene, data=data_C2C3_AGG_01_1d, Mode=x[1,'Mode'], Sample=x[1,'Sample'], 
                          NO_X_TEXT=TRUE, NO_GRIDLINES = TRUE, PRINT_GENES = FALSE, LEGEND=FALSE, offset = -0.5, ALIGN=TRUE)
    } else return(NULL)
    })
  
  Ps %>% rev %>%
    #{lapply(names(.), function(name) {
    #  lapply(c('.svg', '.pdf', '.jpg'), function(extension) ggsave(paste0('Comb_', GO_name, name, extension), plot=.[[name]], width=4.2,height=3))
    #})}
    plot_grid(plotlist=.) %>%
    {lapply(c('.svg', '.pdf', '.jpg'), function(extension) ggsave(paste0('all_Comb_', GO_name, extension), plot=., width=8,height=5))}
})



lapply('GO_ribo', function(GO_name) {
  Ps = subset(reconstituting_DFs_10$reconstituted, Gene %in% merge(GOs[[GO_name]], Cyclers)$Gene) %>%
    {split(., list(.$Sample, .$Mode))} %>% 
    lapply(function(x) {if(nrow(x)) {
      plotter_geneCombine(x$Gene, data=data_C2C3_AGG_01_1d, Mode=x[1,'Mode'], Sample=x[1,'Sample'], 
                          NO_X_TEXT=TRUE, NO_GRIDLINES = TRUE, PRINT_GENES = TRUE, LEGEND=TRUE, offset = -0.5, ALIGN=TRUE)
    } else return(NULL)
    })
  
  Ps %>% rev %>%
    {lapply(names(.), function(name) {
      lapply(c('.svg', '.pdf', '.jpg'), function(extension) ggsave(paste0('Comb_', GO_name, name, extension), plot=.[[name]], width=4.2,height=3))
    })}
    #plot_grid(plotlist=.) %>%
    #{lapply(c('.svg', '.pdf', '.jpg'), function(extension) ggsave(paste0('Comb_', GO_name, extension), plot=., width=8,height=5))}
})


save.image()

compute_ccf = function(expr1, expr2, lagMax=11) {
  ccf(expr1, expr2, lag.max=lagMax) %>% 
    {data.frame(.$acf, .$lag)} %>%
    `colnames<-`(c('R', 'lag')) %>%
    mutate(hr_lag = lag*4, 
           mod_hr_lag = (lag*4)%%24, .keep='unused')
}

takeMax_from_ccfComp = function(df, hr_lag_min=-44, hr_lag_max=44) {
  df = subset(df, hr_lag > hr_lag_min & hr_lag < hr_lag_max)
  df = df[which.max(df$R), ]
}

make_crosscor_unfinished = function(GO_df, Sample='C2', Mode='LD', filtered=reconstituting_DFs_10$reconstituted, cyclers=Cyclers) {
  
  hi <- GO_df %>% 
    {subset(filtered, 
            Gene %in% merge(., cyclers)$Gene & Sample == Sample & Mode == Mode, select=c(Gene, Sample, Mode))} %>%
    merge(avg_data_times2_max1) %>%
    mutate(timepoint = as.numeric(timepoint)) %>%
    arrange(timepoint)
  
  ho <- data.frame()
  gene_vec = unique(hi$Gene)
  for (x in seq_along(gene_vec)) {
    for (y in x:length(gene_vec)) {
      x_gene = gene_vec[[x]]
      y_gene = gene_vec[[y]]
      ho = rbind(ho, cbind(x_gene = x_gene, y_gene=y_gene,
                           compute_ccf(subset(hi, Gene %in% x_gene)$avg, subset(hi, Gene %in% y_gene)$avg) %>% 
                             #takeMax_from_ccfComp
                             {.[.$hr_lag == 0, ]}
      ))
    }
  }
  
  p2 <- ho %>% 
    `colnames<-`(c('Var1', 'Var2', 'value', 'hr_lag', 'mod_hr_lag')) %>%
    arrange(Var1) %>%
    group_by(Var1) %>%
    filter(row_number() >= which(Var1 == Var2)) %>%
    ggplot(aes(x=Var1, y=Var2, fill=value)) + 
    geom_tile() +
    theme(axis.text.x = element_text(angle=90, hjust=TRUE)) +
    xlab("") + 
    ylab("") +
    #geom_text(aes(label=hr_lag), size=2, colour='white') +
    lims(fill = c(-1, 1))
  return(p2)
}

#subset(reconstituting_DFs_10$reconstituted, Gene %in% merge(GO_ETC, Cyclers)$Gene, select=c(Gene, Sample, Mode))


make_corrs = function(Genes, sample, mode, data_times = avg_data_times2_01) {
  df = subset(data_times, Gene %in% Genes & Sample == sample & Mode == mode) %>%
    mutate(timepoint = as.numeric(timepoint)) %>%
    arrange(timepoint)
  
  
  cors <- data.frame()
  gene_vec = unique(df$Gene)
  if(length(gene_vec) <= 1) return(FALSE)
  
  for (x in seq_along(gene_vec)) {
    for (y in seq_along(gene_vec)) {
      #for (y in x:length(gene_vec)) {
      x_gene = gene_vec[[x]]
      y_gene = gene_vec[[y]]
      
      subset(df, Gene %in% x_gene)
      subset(df, Gene %in% y_gene)
      
      cors = rbind(cors, cbind(x_gene = x_gene, y_gene=y_gene,
                               r = cor(subset(df, Gene %in% x_gene)$avg, subset(df, Gene %in% y_gene)$avg)
      ))
    }
  }
  return(cors)
}

make_corrplot = function(GO_df, sample='C2', mode='LD', TRI=FALSE, LABEL=FALSE, plotName=NULL, dendName=NULL,
                         filtered=reconstituting_DFs_10$reconstituted, GO_anno=Cyclers, avg_data_times=avg_data_times2_01) {
  df_cyclers <- GO_df %>% 
    {subset(filtered, 
            Gene %in% merge(., GO_anno)$Gene & sample == Sample & mode == Mode, select=c(Gene, Sample, Mode, meta2d_phase))} 
  df = df_cyclers %>%
    merge(avg_data_times) %>%
    mutate(timepoint = as.numeric(timepoint)) %>%
    arrange(timepoint) 
  
  cors <- data.frame()
  gene_vec = unique(df$Gene)
  if(length(gene_vec) <= 1) return(FALSE)
  
  for (x in seq_along(gene_vec)) {
    for (y in seq_along(gene_vec)) {
      #for (y in x:length(gene_vec)) {
      x_gene = gene_vec[[x]]
      y_gene = gene_vec[[y]]
      cors = rbind(cors, cbind(x_gene = x_gene, y_gene=y_gene,
                               r = cor(subset(df, Gene %in% x_gene)$avg, subset(df, Gene %in% y_gene)$avg)
      ))
    }
  }
  
  wide_tri <- mutate(cors, r = as.numeric(r)) %>%
    pivot_wider(names_from='y_gene', values_from='r') %>%
    as.data.frame %>%
    {`rownames<-`(., .$x_gene)} %>%
    mutate(x_gene = NULL) %>%
    .[gene_vec, gene_vec]
  
  wide <- wide_tri
  
  hc = hclust( dist(wide, method = "euclidean"), method = "complete" )
  ord = hc$order
  
  wide_tri <- wide_tri[gene_vec[ord], gene_vec[ord]]
  wide_tri[lower.tri(wide_tri)] = NA
  
  cors_tri <- wide_tri %>%
    {mutate(., x_gene = rownames(.))} %>%
    pivot_longer(cols=-x_gene, names_to = 'y_gene', values_to = 'r')
  
  if (TRI) df_clst = cors_tri
  else df_clst = cors
  df_clst = df_clst %>%
    mutate(x_gene = factor(x_gene, levels=gene_vec[ord]),
           y_gene = factor(y_gene, levels=rev(gene_vec[ord])),
           r = as.numeric(r)) %>% 
    merge(unique(df_cyclers) , by.x = 'y_gene', by.y='Gene') %>%
    mutate(meta2d_phase = round(meta2d_phase%%24, 1),
           r = round(r, 2))
  
  p = df_clst %>%
    ggplot(aes(x=x_gene, y=y_gene, fill=r)) +
    geom_tile() +
    #scale_x_discrete(position='top') +
    theme_minimal() +
    theme(axis.text.x = element_text(angle=90, hjust=TRUE)) +
    {if(!LABEL) theme(legend.position = "none")} +
    labs(x='', y='', fill='R value') +
    geom_text(aes(label=r), size=2, colour='white') +
    scale_fill_gradient2(limits=c(-1,1), low = 'orange', high='blue', na.value='white')
  
  if (!is.null(plotName)) {
    ggsave(file=plotName, plot=p, width=1+length(gene_vec)/4, height=1+length(gene_vec)/4, limitsize=FALSE)
  }
  if (!is.null(dendName)){
    ggsave(file=dendName, plot=grid.arrange(plot(hc)), width=5, height=5)
  }
  
  p
  
  return(cors)
  #return(list(p = p, dend = hc))
  
  if (0) {
    # ALSO PHASES
    L1 <- levels(df_clst$x_gene)
    L2 <- unique(df_clst$meta2d_phase)%%24 %>% round(1)
    p = df_clst %>% 
      ggplot(aes(x=as.numeric(x_gene), y=y_gene, fill=r)) + 
      geom_tile() +
      scale_x_continuous(breaks=1:length(L1), labels=L1,
                         sec.axis = sec_axis(~., breaks=1:length(L2), labels=L2)) +
      scale_fill_gradient2(limits=c(-1,1), low = 'orange', high='blue') +
      theme(axis.text.x = element_text(angle=90, hjust=FALSE), legend.title=element_blank()) +
      scale_x_discrete(position='top')
  }
}

from_corrplot_meanmaker = function(df) {
  if (is.null(ncol(df))) {
    return(NA)
  } else if (!ncol(df)) {
    return(NA)
  } else subset(df, x_gene != y_gene)$r %>% as.numeric %>% mean %>% return
}


run_fisher_general = function(A, B, U) {
  A = unique(A)
  B = unique(B)
  U = unique(U)
  
  if (! all(A %in% U)) stop(paste(paste(A[!A %in% U], collapse=','), "of A not in set U"))
  if (! all(B %in% U)) stop(paste(paste(B[!B %in% U], collapse=','), "of A not in set U"))
  
  
  # lists A and B contained in U
  data.frame(c(length(which(A %in% B)), 
               length(which(A %!in% B))), 
             c(length(which(B %!in% A)), 
               length(which(U %!in% c(A, B))))) %>%
    fisher.test %>%
    .$p %>%
    return
}


{Subset_Stats <- list()
Subset_Stats$subsets <- lapply(c('GO_ribo', "GO_nitro", "GO_NPR", 'GO_ribo_cell', 'GO_ribo_mito', 'GO_ETC'), function(x) merge(GOs[[x]], passQC)$Gene) %>%
  `names<-`(c('GO_ribo', "GO_nitro", "GO_NPR", 'GO_ribo_cell', 'GO_ribo_mito', 'GO_ETC')) %>%
  c(list(ARGs = ARG$converted$min2,
         ETC_mi = merge(GO_ETC, passQC)$Gene %>% {.[startsWith(., 'mt:')]}, 
         ETC_NOmi = merge(GO_ETC, passQC)$Gene %>% {.[!startsWith(., 'mt:')]}, 
         GPCR_4 = c('CrebA', 'Pka-C1', 'Rala', 'dnc'), 
         GPCR_5 = c('CrebA', 'Pka-C1', 'Rala', 'kay', 'dnc'),
         GPCR_6 = c('CrebA', 'Pka-C1', 'Rala', 'kay', 'CrebB', 'dnc'),
         GPCR_cyc = c('CrebA', 'Pka-C1', 'Rala', 'kay'),
         ARGs_main = c('CG14186', 'Hr38', 'sr'), 
         ARGs2FC2x = c('sr', 'CG14186', 'Hr38', 'CG46385'))) %>%
  {append(., list(GO_ribo_noPre = .$GO_ribo[grep('pre', .$GO_ribo, invert=TRUE)]))}


Subset_Stats$r_means <- lapply(Subset_Stats$subsets, function(Set) {
  lapply(levels(data_C2C3$sample_mode), function(S_M) {
    S = gsub('_.*', '', S_M)
    M = gsub('.*_', '', S_M)
    data.frame(Gene = Set) %>%
      make_corrplot(GO_anno = Cyclers, filtered=reconstituting_DFs_10$reconstituted, sample=S, mode=M) %>%
      from_corrplot_meanmaker
  }) %>%
    do.call(what='cbind') %>%
    `colnames<-`(levels(data_C2C3$sample_mode)) 
  }) %>%
  do.call(what='rbind') %>%
  `rownames<-`(names(Subset_Stats$subsets))

Subset_Stats$p_vals <- lapply(Subset_Stats$subsets, function(Set) {
  lapply(levels(data_C2C3$sample_mode), function(S_M) {
    S = gsub('_.*', '', S_M)
    M = gsub('.*_', '', S_M)
    run_fisher_general(B = subset(reconstituting_DFs_10$reconstituted, Sample == S & Mode == M)$Gene,
                       A = subset(to_filter, Gene %in% Set & second_pct > 0.10 & Sample == S & Mode == M)$Gene, 
                       U = subset(to_filter, second_pct > 0.10 & Sample == S & Mode == M)$Gene)
  }) %>%
    do.call(what='cbind') %>%
    `colnames<-`(levels(data_C2C3$sample_mode))
}) %>%
  do.call(what='rbind') %>%
  `rownames<-`(names(Subset_Stats$subsets))

#Subset_Stats$r_means_Teddy <- lapply(Subset_Stats$subsets, function(Set) {
#  lapply(levels(data_C2C3$sample_mode), function(S_M) {
#    S = gsub('_.*', '', S_M)
#    M = gsub('.*_', '', S_M)
#    data.frame(Gene = Set) %>%
#      make_corrplot(GO_anno = Cyclers, filtered=my_filtered, sample=S, mode=M) %>%
#      from_corrplot_meanmaker
#  }) %>%
#    do.call(what='cbind') %>%
#    `colnames<-`(levels(data_C2C3$sample_mode)) 
#}) %>%
#  do.call(what='rbind') %>%
#  `rownames<-`(names(Subset_Stats$subsets))

#Subset_Stats$p_vals_Teddy <- lapply(Subset_Stats$subsets, function(Set) {
#  lapply(levels(data_C2C3$sample_mode), function(S_M) {
#    S = gsub('_.*', '', S_M)
#    M = gsub('.*_', '', S_M)
#    run_fisher_general(B = subset(my_filtered, Sample == S & Mode == M)$Gene,
#                       A = subset(to_filter, Gene %in% Set & second_pct > 0.10 & Sample == S & Mode == M)$Gene, 
#                       U = subset(to_filter, second_pct > 0.10 & Sample == S & Mode == M)$Gene)
#  }) %>%
#    do.call(what='cbind') %>%
#    `colnames<-`(levels(data_C2C3$sample_mode))
#}) %>%
#  do.call(what='rbind') %>%
#  `rownames<-`(names(Subset_Stats$subsets))
}



lapply(levels(data_C2C3$sample_mode), function(S_M) {
  S = gsub('_.*', '', S_M)
  M = gsub('.*_', '', S_M)
  data.frame(Gene = c('CG14186', 'Hr38', 'sr')) %>%
    make_corrplot(GO_anno = passQC, filtered=to_filter, sample=S, mode=M) %>%
    from_corrplot_meanmaker
}) %>%
  do.call(what='cbind') %>%
  `colnames<-`(levels(data_C2C3$sample_mode))

c('CrebA', 'Pka-C1', 'Rala', 'kay', 'CrebB', 'dnc') %>%
  {subset(reconstituting_DFs_10$reconstituted, Gene %in% . & Sample=='C3' & Mode=='LD')}


lapply(Subset_Stats$subsets[c('GPCR_4', 'GPCR_5', 'GPCR_6')], function(Set) {
  lapply(levels(data_C2C3$sample_mode), function(S_M) {
    S = gsub('_.*', '', S_M)
    M = gsub('.*_', '', S_M)
    data.frame(Gene = Set) %>%
      make_corrplot(GO_anno=passQC, filtered=subset(to_filter, Sample==S & Mode==M & second_pct>0.1), sample=S, mode=M) %>%
      from_corrplot_meanmaker
  }) %>%
    do.call(what='cbind') %>%
    `colnames<-`(levels(data_C2C3$sample_mode)) 
}) %>%
  do.call(what='rbind') %>%
  `rownames<-`(c('GPCR_4', 'GPCR_5', 'GPCR_6'))

Subset_Stats$subsets$GO_ribo %>%
  {subset(reconstituting_DFs_10$reconstituted, Gene %in% . & Sample=='C3'&Mode=='LD')} %>%
  .$Gene %>%
  avg_plotter_concat(data=subset(data_C2C3_AGG, sample=='C3'&mode1=='LD'), N_days = 2)



####    ####


#### LESS IMPORTANT ####

hi <- data.frame(GO=c('GO:0044271', 'GO:0003735', 'GO:0008188'), 
           type=c('BP', 'MF', 'MF')) %>%
  apply(1, function(df) {
    GO = df[['GO']]
    print(GO)
    if (df[['type']] == 'BP') {
      the_GO = AnnotationDbi::select(org.Dm.eg.db,
                      keys=c(GO, GOBPOFFSPRING[[GO]]),
                      columns="ENSEMBL",
                      keytype="GO")['ENSEMBL'] %>% unique
    } else if (df[['type']] == 'MF') {
      the_GO = AnnotationDbi::select(org.Dm.eg.db,
                      keys=c(GO, GOMFOFFSPRING[[GO]]),
                      columns="ENSEMBL",
                      keytype="GO")['ENSEMBL'] %>% unique
    } else if (df[['type']] == 'CC') {
      the_GO = AnnotationDbi::select(org.Dm.eg.db,
                      keys=c(GO, GOCCOFFSPRING[[GO]]),
                      columns="ENSEMBL",
                      keytype="GO")['ENSEMBL'] %>% unique
    }
    data.frame(
      make_corrplot(the_GO, sample='C3', mode='LD', GO_anno = Cyclers, filtered=reconstituting_DFs_10$reconstituted) %>%
        from_corrplot_meanmaker,
      make_corrplot(the_GO, sample='C2', mode='LD', GO_anno = Cyclers, filtered=reconstituting_DFs_10$reconstituted) %>%
        from_corrplot_meanmaker,
      make_corrplot(the_GO, sample='C3', mode='DD', GO_anno = Cyclers, filtered=reconstituting_DFs_10$reconstituted) %>%
        from_corrplot_meanmaker,
      make_corrplot(the_GO, sample='C2', mode='DD', GO_anno = Cyclers, filtered=reconstituting_DFs_10$reconstituted) %>%
        from_corrplot_meanmaker) %>%
      `colnames<-`(levels(data_C2C3$sample_mode)) %>%
      return
  })
hi %>%
  do.call(what='rbind') %>%
  `rownames<-`(c('nitrogen biosynthesis', 'ribosome', 'neuropeptide receptor'))

apply(c('GO:0044271', 'GO:0003735', 'GO:0008188'), function(GO) {
  print(GO)
  the_GO = select(org.Dm.eg.db,
                   keys=c(GO, GOCCOFFSPRING[[GO]]),
                   columns="ENSEMBL",
                   keytype="GO")['ENSEMBL'] %>% unique
  
  data.frame(
    make_corrplot(the_GO, sample='C3', mode='LD', GO_anno = Cyclers, filtered=reconstituting_DFs_10$reconstituted) %>%
      from_corrplot_meanmaker,
    make_corrplot(the_GO, sample='C2', mode='LD', GO_anno = Cyclers, filtered=reconstituting_DFs_10$reconstituted) %>%
      from_corrplot_meanmaker,
    make_corrplot(the_GO, sample='C3', mode='DD', GO_anno = Cyclers, filtered=reconstituting_DFs_10$reconstituted) %>%
      from_corrplot_meanmaker,
    make_corrplot(the_GO, sample='C2', mode='DD', GO_anno = Cyclers, filtered=reconstituting_DFs_10$reconstituted) %>%
      from_corrplot_meanmaker) %>%
    `colnames<-`(levels(data_C2C3$sample_mode)) %>%
    return
}) %>%
  do.call(what='rbind') %>%
  `rownames<-`(c('nitrogen biosynthesis', 'ribosome', 'neuropeptide receptor'))


  
# FOR TOMORROW, YOU GET TO MAKE CORR CHECKER:
#   PER CONDITION
#     PULL GO TERMS FOR CYCLERS
#     TAKE UNIQUE GO TERMS, NOTE BP/MF/CC
#     FOR EACH GO TERM
#       PULL ALL TERM AND CHILDREN (BP/MF/CC) GENES IN merge(passQC, subset(to_filter, second_pct>0.1, Sample = ???, Mode = ???)$Gene)
#     SAVE THIS
#   PER CONDITION
#     REMOVE OVER 50, UNDER 5
#     CALCULATE R's, return R average for each
#   check out GO terms with high coordination, pull # cyclers per go term, maybe rate by % cyclers times r value


lapply(c('GO_ribo_cell', 'GO_ribo_mito', 'GO_ETC'), function(name) {
  make_corrplot(GOs[[name]], sample='C3', mode='LD', plotName = paste0('all_', name, '_C3LD.pdf'), GO_anno = passQC, filtered=subset(to_filter, second_pct>0.10))
  make_corrplot(GOs[[name]], sample='C2', mode='LD', plotName = paste0('all_', name, '_C2LD.pdf'), GO_anno = passQC, filtered=subset(to_filter, second_pct>0.10))
  make_corrplot(GOs[[name]], sample='C3', mode='DD', plotName = paste0('all_', name, '_C3DD.pdf'), GO_anno = passQC, filtered=subset(to_filter, second_pct>0.10))
  make_corrplot(GOs[[name]], sample='C2', mode='DD', plotName = paste0('all_', name, '_C2DD.pdf'), GO_anno = passQC, filtered=subset(to_filter, second_pct>0.10))
}) 

# corr_plots for all ETC_mi:
lapply(levels(data_C2C3$sample_mode), function(S_M) {
  S = gsub('_.*', '', S_M)
  M = gsub('.*_', '', S_M)
  merge(GO_ETC, passQC)$Gene %>% {.[startsWith(., 'mt:')]} %>% data.frame(Gene = .) %>%
    make_corrplot(GO_anno=passQC, filtered=subset(to_filter, second_pct>0.10), sample=S, mode=M, plotName=paste0('all_MT_passQC_mi_', S, M, '.pdf'))
})

# Collective plot for mt: only
lapply('GO_ETC', function(GO_name) {
  Ps = subset(reconstituting_DFs_10$reconstituted, Gene %in% (merge(GOs[[GO_name]], Cyclers)$Gene %>% {.[startsWith(., 'mt:')]})) %>%
    {split(., list(.$Sample, .$Mode))} %>% 
    lapply(function(x) {if(nrow(x)) {
      plotter_geneCombine(x$Gene, data=data_C2C3_AGG_01_1d, Mode=x[1,'Mode'], Sample=x[1,'Sample'], 
                          NO_X_TEXT=TRUE, NO_GRIDLINES = TRUE, PRINT_GENES = TRUE, LEGEND=TRUE, offset = -0.5, ALIGN=TRUE)
    } else return(NULL)
    })
  Ps
  Ps %>% rev %>%
    #{lapply(names(.), function(name) {
    #  lapply(c('.svg', '.pdf', '.jpg'), function(extension) ggsave(paste0('Comb_', GO_name, name, extension), plot=.[[name]], width=4.2,height=3))
    #})}
    plot_grid(plotlist=.) %>%
    {lapply(c('.svg', '.pdf', '.jpg'), function(extension) ggsave(paste0('Comb_', GO_name, '_mi', extension), plot=., width=8,height=5))}
})

# Collective plot for ETC-non-mt: only
lapply('GO_ETC', function(GO_name) {
  Ps = subset(reconstituting_DFs_10$reconstituted, Gene %in% (merge(GOs[[GO_name]], Cyclers)$Gene %>% {.[!startsWith(., 'mt:')]})) %>%
    {split(., list(.$Sample, .$Mode))} %>% 
    lapply(function(x) {if(nrow(x)) {
      plotter_geneCombine(x$Gene, data=data_C2C3_AGG_01_1d, Mode=x[1,'Mode'], Sample=x[1,'Sample'], 
                          NO_X_TEXT=TRUE, NO_GRIDLINES = TRUE, PRINT_GENES = TRUE, LEGEND=TRUE, offset = -0.5, ALIGN=TRUE)
    } else return(NULL)
    })
  Ps
  Ps %>% rev %>%
    #{lapply(names(.), function(name) {
    #  lapply(c('.svg', '.pdf', '.jpg'), function(extension) ggsave(paste0('Comb_', GO_name, name, extension), plot=.[[name]], width=4.2,height=3))
    #})}
    plot_grid(plotlist=.) %>%
    {lapply(c('.svg', '.pdf', '.jpg'), function(extension) ggsave(paste0('Comb_', GO_name, '_NOTmt', extension), plot=., width=8,height=5))}
})

# Collective plot for expressed mt: only
lapply('GO_ETC', function(GO_name) {
  Ps = subset(to_filter, Gene %in% (merge(GOs[[GO_name]], passQC)$Gene %>% {.[startsWith(., 'mt:')]})) %>%
    {split(., list(.$Sample, .$Mode))} %>% 
    lapply(function(x) {if(nrow(x)) {
      plotter_geneCombine(x$Gene, data=data_C2C3_AGG_01_1d, Mode=x[1,'Mode'], Sample=x[1,'Sample'], 
                          NO_X_TEXT=TRUE, NO_GRIDLINES = TRUE, PRINT_GENES = TRUE, LEGEND=TRUE, offset = -0.5, ALIGN=TRUE)
    } else return(NULL)
    })
  Ps
  Ps %>% rev %>%
    #{lapply(names(.), function(name) {
    #  lapply(c('.svg', '.pdf', '.jpg'), function(extension) ggsave(paste0('Comb_', GO_name, name, extension), plot=.[[name]], width=4.2,height=3))
    #})}
    plot_grid(plotlist=.) %>%
    {lapply(c('.svg', '.pdf', '.jpg'), function(extension) ggsave(paste0('all_Comb_', GO_name, '_mi', extension), plot=., width=8,height=5))}
})

# Collective plot for expressed ETC-non-mt: only
lapply('GO_ETC', function(GO_name) {
  Ps = subset(to_filter, Gene %in% (merge(GOs[[GO_name]], passQC)$Gene %>% {.[!startsWith(., 'mt:')]})) %>%
    {split(., list(.$Sample, .$Mode))} %>% 
    lapply(function(x) {if(nrow(x)) {
      plotter_geneCombine(x$Gene, data=data_C2C3_AGG_01_1d, Mode=x[1,'Mode'], Sample=x[1,'Sample'], 
                          NO_X_TEXT=TRUE, NO_GRIDLINES = TRUE, PRINT_GENES = TRUE, LEGEND=FALSE, offset = -0.5, ALIGN=TRUE)
    } else return(NULL)
    })
  Ps
  Ps %>% rev %>%
    #{lapply(names(.), function(name) {
    #  lapply(c('.svg', '.pdf', '.jpg'), function(extension) ggsave(paste0('Comb_', GO_name, name, extension), plot=.[[name]], width=4.2,height=3))
    #})}
    plot_grid(plotlist=.) %>%
    {lapply(c('.svg', '.pdf', '.jpg'), function(extension) ggsave(paste0('all_Comb_', GO_name, '_NOTmt', extension), plot=., width=8,height=5))}
})



# PKA TRANSDUCTION:
c('CrebA', 'CrebB', 'Pka-C1', 'Rala', 'dnc') %>%
  {subset(reconstituting_DFs_10$reconstituted, Gene %in% .)}
c('CrebA', 'Pka-C1', 'Rala') %>%
  {subset(to_filter, Gene %in% .)} %>%
  {split(., list(.$Sample, .$Mode))} %>% print %>%
  lapply(function(x) {if(nrow(x)) {
    plotter_geneCombine(x$Gene, data=data_C2C3_AGG_01_1d, Mode=x[1,'Mode'], Sample=x[1,'Sample'], 
                        NO_X_TEXT=TRUE, NO_GRIDLINES = TRUE, PRINT_GENES = TRUE, LEGEND=TRUE, offset = -0.5, ALIGN=TRUE)
  } else return(NULL)
  }
  ) %>% 
  rev %>%
  #{lapply(names(.), function(name) {
  #  lapply(c('.svg', '.pdf', '.jpg'), function(extension) ggsave(paste0('all_Comb_', GO_name, name, extension), plot=.[[name]], width=4.2,height=3))
  #})}
  plot_grid(plotlist=.) %>%
  {lapply(c('.svg', '.pdf', '.jpg'), function(extension) ggsave(paste0('all_Comb_', 'PKA_transduction', extension), plot=., width=8,height=5))}


# mitochondrial genome (in passQC)
lapply(levels(data_C2C3$sample_mode), function(S_M) {
  S = gsub('_.*', '', S_M)
  M = gsub('.*_', '', S_M)
  rownames(data_C2C3) %>% {.[startsWith(., 'mt:')]} %>% data.frame(Gene = .) %>%
    make_corrplot(GO_anno=passQC, filtered=subset(to_filter, second_pct>0.10), sample=S, mode=M, plotName=paste0('all_MT_passQC_', S, M, '.pdf'))
})

hi <- rownames(data_C2C3) %>% {.[startsWith(., 'mt:')]} %>% data.frame(Gene = .) %>%
  make_corrplot(GO_anno=passQC, filtered=subset(to_filter, second_pct>0.10), sample='C3', mode='LD')

subset(hi, x_gene != y_gene)$r %>% as.numeric %>% mean

{x = subset(passQC, Gene %in% ARG$converted$min2)['ENSEMBL']
  make_corrplot(x, sample='C3', mode='LD', plotName = paste0('all_ARG2x_C3LD.png'), GO_anno = passQC, filtered=subset(to_filter, second_pct>0.10&Sample=='C3'&Mode=='LD'))
  make_corrplot(x, sample='C2', mode='LD', plotName = paste0('all_ARG2x_C2LD.png'), GO_anno = passQC, filtered=subset(to_filter, second_pct>0.10&Sample=='C2'&Mode=='LD'))
  make_corrplot(x, sample='C3', mode='DD', plotName = paste0('all_ARG2x_C3DD.png'), GO_anno = passQC, filtered=subset(to_filter, second_pct>0.10&Sample=='C3'&Mode=='DD'))
  make_corrplot(x, sample='C2', mode='DD', plotName = paste0('all_ARG2x_C2DD.png'), GO_anno = passQC, filtered=subset(to_filter, second_pct>0.10&Sample=='C2'&Mode=='DD'))
  rm(x)
}








####    ####    ####    ####    ####    ####    ####    ####    ####    ####   


# 29. PEPTIDE/NT PLOTS
##################

# ^ run on avg_max1; scale up; inerpolate

genes_GPCR_but_smaller <- c("Dh44-R2", "CCAP-R", "CCKLR-17D3", "FMRFaR", "CG13579", "CG33639", 
                    "CG30340", "TkR86C", "CCHa1-R", "RYa-R", "CNMaR", "rk", "hec", 
                    "CCHa2-R", "AstA-R1", "AstA-R2", "CapaR", "MsR1", "Lkr", "CG13575", 
                    "TrissinR", "Proc-R", "PK2-R1", "MsR2", "AkhR", "SIFaR", "PK2-R2", 
                    "CCKLR-17D1", "CG32547", "Lgr1", "CG13995", "PK1-R", "Lgr3", 
                    "CG13229", "Dh31-R", "AstC-R1", "Pdfr", "Dh44-R1", "TkR99D", 
                    "NPFR", "ETHR", "Lgr4", "sNPF-R", "CrzR", "AstC-R2") %>% unique %>% sort

genes_GPCR <- c("5-HT1A", "5-HT1B", "5-HT2A", "5-HT2B", "5-HT7", "AdoR", "AkhR", "AstA-R1",
             "AstA-R2", "AstC-R1", "AstC-R2", "boss", "CapaR", "CCAP-R", "CCHa1-R", 
             "CCHa2-R", "CCKLR-14D1", "CCKLR-17D3", "CG11318", "CG12290", "mAChR-C", 
             "CG13229", "CG13575", "CG13579", "CG13995", "CG15556", "CG15614", "CG15744", 
             "Octalpha2R", "CG30340", "CG31760", "CG32447", "CG32547", "CG33639", "CG4313", 
             "CG44153", "CG7497", "Cirl", "CNMaR", "CrzR", "Dh31-R", "Dh44-R1", "Dh44-R2", 
             "Dop1R1", "Dop1R2", "Dop2R", "DopEcR", "ETHR", "FMRFaR", "fz", "fz2", "fz3", 
             "fz4", "GABA-B-R1", "GABA-B-R2", "GABA-B-R3", "hec", "Lgr1", "Lgr3", "Lgr4", 
             "Lkr", "mAChR-A", "mAChR-B", "mAChR-C", "mGluR", "moody", "MsR1", "MsR2", 
             "mth", "mthl1", "mthl10", "mthl11", "mthl12", "mthl13", "mthl14", "mthl15", 
             "mthl2", "mthl3", "mthl4", "mthl5", "mthl6", "mthl7", "mthl8", "mthl9", 
             "mtt", "ninaE", "NPFR", "Oamb", "Octbeta1R", "Octbeta2R", "Octbeta3R", 
             "Oct-TyrR", "Pdfr", "PK1-R", "PK2-R1", "PK2-R2", "Proc-R", "Rh2", "Rh3", 
             "Rh4", "Rh5", "Rh6", "Rh7", "rho", "rho-4", "rho-5", "rho-6", "rho-7", 
             "rk", "ru", "RYa-R", "smo", "smog", "sNPF-R", "SPR", "stan", "stet", 
             "TkR86C", "TkR99D", "Tre1", "TrissinR", "TyrR", "TyrRII") %>% unique %>% sort

genes_neuropeptides <- c("amn", "Akh", "CCHa2", "Lk", "SIFa", "Crz", "Dsk", "Nplp2", "CCHa1", 
                   "NPF", "Proc", "Hug", "ITP", "AstA", "Tk", "Nplp4", "Ms", "ETH", "Burs", 
                   "Nplp1", "Nplp3", "Eh", "Pdf", "Pburs", "Capa", "sNPF", "CNMa", "Mip", "Ptth") %>% unique %>% sort

make_heatmap = function(Genes, avgs = mutate(avg_data, Type = paste(Sample, Mode, sep='')), GENE_SCALE=FALSE, CELL_SCALE=FALSE, LOG1P=FALSE, CLUST_ROWS=FALSE, CLUST_COLS=FALSE) {
  if (GENE_SCALE*LOG1P + CELL_SCALE*LOG1P) warning("SCALEs override LOG1P")
  if (GENE_SCALE*CELL_SCALE) title = 'Cell+Gene-scaled\nAverage Expression'
  else if (GENE_SCALE) title = 'Gene-scaled\nAverage Expression' 
  else if (CELL_SCALE) title = 'Cell-scaled\nAverage Expression' 
  else if (LOG1P) title = 'log1p Average\nExpression'
  else title = 'Average\nExpression'
  
  
  
  if (!all(Genes %in% avgs$Gene)) {
    warning(paste0(Genes[Genes %!in% avgs$Gene], " not in avg_data", collapse=', '))
    Genes = subset(avgs, Gene %in% Genes)$Gene
  }
  
  avgs %>%
    subset(Gene %in% Genes, select=c(Gene, avg, Type)) %>% 
    #group_by(Gene, Type) %>%
    #summarise_at(vars(avg), list(avg = mean)) %>%
    {if (LOG1P) mutate(., avg = log1p(avg)) else .} %>% 
    pivot_wider(names_from = 'Gene', values_from = 'avg') %>% 
    as.data.frame %>% 
    #{.[, c('Type', Genes)]} %>% 
    {`rownames<-`(., .$Type)} %>% 
    subset(select = -Type) %>% 
    {if (CELL_SCALE) t(scale(t(.))) else .} %>% 
    {if (GENE_SCALE) scale(as.matrix(.)) else .} %>% 
    as.matrix %>% print %>%
    Heatmap(cluster_rows=CLUST_ROWS, cluster_columns=CLUST_COLS, 
            column_names_gp = grid::gpar(fontsize = 8),
            row_names_gp = grid::gpar(fontsize = 8), heatmap_legend_param = list(title=title), row_names_side='left')
}


{setwd(WD)
data_avg <- avg_data %>% group_by(Sample, Gene) %>%
  summarise(avg = mean(avg)) %>%
  mutate(Type = Sample, .keep='unused')

  dir.create('heatmaps', showWarnings = FALSE)
  setwd('./heatmaps')
  make_heatmap(genes_neuropeptides, LOG1P=TRUE, avgs=data_avg) %>% 
    draw %>% grid.grabExpr %>% 
    {lapply(c('.svg', '.pdf', '.jpg'), function(extension) ggsave(paste0('NP', extension), plot=., width=length(genes_neuropeptides)/6+1,height=2/6+0.55))}
  
  make_heatmap(genes_GPCR, LOG1P=TRUE, avgs=data_avg) %>% 
    draw %>% grid.grabExpr %>% 
    {lapply(c('.svg', '.pdf', '.jpg'), function(extension) ggsave(paste0('NPR', extension), plot=., width=length(genes_GPCR)/6+1,height=2/6+0.9))}
  
  G_NP = subset(data_avg, Gene %in% genes_neuropeptides) %>% group_by(Gene) %>% summarize(max_avg = max(avg)) %>% subset(max_avg > 0.1) %>% .$Gene
  G_NPR = subset(data_avg, Gene %in% genes_GPCR) %>% group_by(Gene) %>% summarize(max_avg = max(avg)) %>% subset(max_avg > 0.1) %>% .$Gene
  
  make_heatmap(G_NP, LOG1P=TRUE, avgs=data_avg) %>% 
    draw %>% grid.grabExpr %>% 
    {lapply(c('.svg', '.pdf', '.jpg'), function(extension) ggsave(paste0('NP_min0.1', extension), plot=., width=length(G_NP)/6+1,height=2/6+0.55))}
  
  make_heatmap(G_NPR, LOG1P=TRUE, avgs=data_avg) %>% 
    draw %>% grid.grabExpr %>% 
    {lapply(c('.svg', '.pdf', '.jpg'), function(extension) ggsave(paste0('NPR_min0.1', extension), plot=., width=length(G_NPR)/6+1,height=2/6+0.9))}
  
  
  make_heatmap(genes_neuropeptides, LOG1P=TRUE, avgs=data_avg) %>% 
    draw %>% grid.grabExpr %>% 
    {lapply(c('.svg', '.pdf', '.jpg'), function(extension) ggsave(paste0('NP_forLegend', extension), plot=., width=length(genes_neuropeptides)/6+1,height=2))}
  make_heatmap(genes_GPCR, LOG1P=TRUE, avgs=data_avg) %>% 
    draw %>% grid.grabExpr %>% 
    {lapply(c('.svg', '.pdf', '.jpg'), function(extension) ggsave(paste0('NPR_forLegend', extension), plot=., width=length(genes_GPCR)/6+1,height=2))}
  make_heatmap(G_NP, LOG1P=TRUE, avgs=data_avg) %>% 
    draw %>% grid.grabExpr %>% 
    {lapply(c('.svg', '.pdf', '.jpg'), function(extension) ggsave(paste0('NP_min0.1_forLegend', extension), plot=., width=length(G_NP)/6+1,height=2))}
  make_heatmap(G_NPR, LOG1P=TRUE, avgs=data_avg) %>% 
    draw %>% grid.grabExpr %>% 
    {lapply(c('.svg', '.pdf', '.jpg'), function(extension) ggsave(paste0('NPR_min0.1_forLegend', extension), plot=., width=length(G_NPR)/6+1,height=2))}
  
  rm(G_NP, G_NPR, data_avg)
  setwd(WD)
}


make_heatmap(genes_GPCR) %>% 
  draw %>% grid.grabExpr %>%
  ggsave(file='hmap_NPR.svg', plot=., width=length(genes_GPCR)/6, height=2, limitsize=FALSE)
make_heatmap(genes_GPCR, LOG1P=TRUE) %>% 
  draw %>% grid.grabExpr %>%
  ggsave(file='hmap_NPR_log1p.svg', width=length(genes_GPCR)/6, height=2, limitsize=FALSE)
make_heatmap(genes_GPCR, GENE_SCALE=TRUE) %>% 
  draw %>% grid.grabExpr %>%
  ggsave(file='hmap_NPR_genescale.svg', width=length(genes_GPCR)/6, height=2, limitsize=FALSE)
make_heatmap(genes_GPCR, CELL_SCALE=TRUE) %>% 
  draw %>% grid.grabExpr %>%
  ggsave(file='hmap_NPR_cellscale.svg', width=length(genes_GPCR)/6, height=2, limitsize=FALSE)
make_heatmap(genes_GPCR, GENE_SCALE=TRUE, CELL_SCALE=TRUE) %>% 
  draw %>% grid.grabExpr %>%
  ggsave(file='hmap_NPR_bothscale.svg', width=length(genes_GPCR)/6, height=2, limitsize=FALSE)

make_heatmap(genes_neuropeptides) %>% 
  draw %>% grid.grabExpr %>%
  ggsave(file='hmap_NP.svg', plot=., width=length(genes_neuropeptides)/6, height=2, limitsize=FALSE)
make_heatmap(genes_neuropeptides, LOG1P=TRUE) %>% 
  draw %>% grid.grabExpr %>%
  ggsave(file='hmap_NP_log1p.svg', plot=., width=length(genes_neuropeptides)/6, height=2, limitsize=FALSE)
make_heatmap(genes_neuropeptides, GENE_SCALE=TRUE) %>% 
  draw %>% grid.grabExpr %>%
  ggsave(file='hmap_NP_genescale.svg', plot=., width=length(genes_neuropeptides)/6, height=2, limitsize=FALSE)
make_heatmap(genes_neuropeptides, CELL_SCALE=TRUE) %>% 
  draw %>% grid.grabExpr %>%
  ggsave(file='hmap_NP_cellscale.svg', plot=., width=length(genes_neuropeptides)/6, height=2, limitsize=FALSE)
make_heatmap(genes_neuropeptides, GENE_SCALE=TRUE, CELL_SCALE=TRUE) %>% 
  draw %>% grid.grabExpr %>%
  ggsave(file='hmap_NP_bothscale.svg', plot=., width=length(genes_neuropeptides)/6, height=2, limitsize=FALSE)

####    ####    ####    ####    ####    ####    ####    ####    ####    ####   


# 30. OTHER
#######

avg_plotter_concat_mode(c('REPTOR', 'RapGAP1', 'RasGAP1'), data=data_C2C3_AGG)
avg_plotter_concat_mode(c('REPTOR', 'Akt1', 's6k', 's6kII', 'Thor', 'sgg', 'InR', 'chico',  'RapGAP1', 'RasGAP1', 'foxo'), data=data_C2C3_AGG)

c('REPTOR', 'Akt', 'S6K', 'S6KII', 'Thor', 'sgg')

select(org.Dm.eg.db,
       keys=c("GO:0005942", GOCCOFFSPRING[["GO:0005942"]]),
       columns=c("ENSEMBL", "SYMBOL"),
       keytype="GO")$SYMBOL %>% unique %>%
  avg_plotter_concat_mode(data=data_C2C3_AGG)

subset(symbol_ensembl, ENSEMBL %in% c('FBgn0283472', 'FBgn0010379', 'FBgn0262866', 'FBgn0024248'))

GO_PI3K <- select(org.Dm.eg.db,
                  keys=c("GO:0005942", GOCCOFFSPRING[["GO:0005942"]]),
                  columns="ENSEMBL",
                  keytype="GO")['ENSEMBL'] %>% unique

GO_AC <- select(org.Dm.eg.db,
                           keys=c("GO:0004016", GOMFOFFSPRING[["GO:0004016"]]),
                           columns="ENSEMBL",
                           keytype="GO")['ENSEMBL'] %>% unique
make_panel(GO_AC, filename='GO_AC_panel.svg')
make_panel(GO_AC, filename='GO_AC_panel_expressed.svg', pass_df=subset(to_filter, second_pct > 0.1), PLOT_NESTING=FALSE)
make_phase_from_GOdf(GO_AC) %>% ggsave(file='GO_AC_hmap.svg', plot=., width=5, height=5)

DotPlot(data_C2C3, features=merge(GO_AC, Cyclers)$Gene, scale=FALSE)
DotPlot(data_C2C3, features=merge(GO_AC, symbol_ensembl)$SYMBOL, scale=FALSE)

subset(reconstituting_DFs_10$reconstituted, Gene %in% merge(GO_AC, Cyclers)$Gene)
# ^ Ac3 C3DD (also looks cycling in C2DD); CG10738, CG43373
# CONCLUSION: perhaps Ac3 cycling is driven by activity; both express Ac3 (circadian, Gs) and 



GO_Gprot <- select(org.Dm.eg.db,
                keys=c("GO:0005834", GOCCOFFSPRING[["GO:0005834"]]),
                columns="ENSEMBL",
                keytype="GO")['ENSEMBL'] %>% unique

DotPlot(data_C2C3, features=merge(GO_Gprot, Cyclers)$Gene)
DotPlot(data_C2C3, features=merge(GO_Gprot, symbol_ensembl)$SYMBOL, scale=FALSE)

subset(reconstituting_DFs_10$reconstituted, Gene %in% merge(GO_Gprot, Cyclers)$Gene)

avg_plotter_concat_mode('Galphas', rib_alpha=0.4)

DotPlot(data_C2C3, features=rownames(data_C2C3) %>% {c(.[startsWith(., 'Galpha')], 'CG17760', 'CG30054', 'cta')}, scale=FALSE)
# ^ alpha q, alpha o, alpha s <- C2 don't lack alpha s; 


GO_PLC <- select(org.Dm.eg.db,
                   keys=c("GO:0004629", GOMFOFFSPRING[["GO:0004629"]]),
                   columns="ENSEMBL",
                   keytype="GO")['ENSEMBL'] %>% unique
DotPlot(data_C2C3, features=merge(GO_PLC, symbol_ensembl)$SYMBOL, scale=FALSE)
# ^ YES Plc21C and norpA <- expecting IP3 signaling


# RAS
all_TFs_COs_w_child_gomf <- select(org.Dm.eg.db,
                                   keys=c("GO:GO:0005096", GOMFOFFSPRING[["GO:0003712"]],
                                          "GO:0003700", GOMFOFFSPRING[["GO:0003700"]]),
                                   columns=c("SYMBOL", "ENSEMBL"),
                                   keytype="GO") #725, per, Clk, vri, cwo, Pdp1 (no tim, cry)

#GO = list()
GO$Ras <- select(org.Dm.eg.db,
       keys=c("GO:0005096"), 
       columns=c("SYMBOL", "ENSEMBL"),
       keytype="GO")

GO$RTK <- select(org.Dm.eg.db,
                       keys=c("GO:0004714"), 
                       columns=c("SYMBOL", "ENSEMBL"),
                       keytype="GO")


select(org.Dm.eg.db,
       keys=c("FBgn0265276"), 
       columns=c("SYMBOL", "ENSEMBL"),
       keytype="ENSEMBL")



# Pdfr-related genes of interest:

c('CaMKII', 'rl', 'Dsor1', 'Erk7', 'nmo', 'bsk', 'p38a', 'p38b', 'p38c', 'Raf', 'Ras85D', 'Pp1-87B', 'Pp1-Y1', 'Pp1-Y2', 'Pp1-13C', 'dnc', 'rut', 'Jra', 'Atf3', 'kay', 'CG7786', 'crc') %>%
  avg_plotter_concat_mode(rib_alpha = 0.4, N_days = 2)

# Cyclers from above:
c('rl', 'dnc', 'kay', 'crc') # ATF crc, phosphodiesterase dnc, ERK rl (questionable), Jun homologue Fos

c('Crtc', 'CanA1', 'CanB', 'CanB2', 'Pp2B-14D') %>% 
  avg_plotter_concat_mode(rib_alpha = 0.4, N_days = 2)



ChIP_Hardin$CWO$`Gene Name` %>% {.[. %in% res10$C3_LD$candidates$SYMBOL]} %>% unique %>%
  avg_plotter_concat_mode(data=data_C2C3_AGG)


res_ribo <- list()
res_ribo <- subset(reconstituting_DFs_10$reconstituted, Gene %in% merge(rbind(GO_ribo_cell, GO_ribo_mito), Cyclers)$Gene) %>%
  {split(., list(.$Sample, .$Mode))} %>%
  lapply(function(df) make_FIMO_run(df, samples = unique(df$Sample), modes = unique(df$Mode), cyclers = Cyclers, QNOTP = FALSE, WRITE_PROMOTERS=TRUE)) %>%
  {`names<-`(., gsub('\\.', '_', names(.)))}
lapply(res_ribo, function(res) {res$results$TF %>% unique %>% length})

# GET CYCLING TFs
lapply(names(res_ribo), function(name) {
  res = res_ribo[[name]]
  resAll = res10[[name]]
  hi = resAll$candidates$SYMBOL
  res$results$TF %>% 
    {.[. %in% resAll$candidates$SYMBOL]} %>%
    table %>% sort
})


# thresholding overall:
lapply(names(res_ribo), function(name) {
  res = res_ribo[[name]]
  res$results$TF %>% 
    table %>% sort %>%
    c %>% {.[. > 10]}
})

lapply(res_ribo, function(res) res$candidates$SYMBOL %>% length)

lapply(names(res_ribo), function(name) {
  res = res_ribo[[name]]
  resAll = res10[[name]]
  namebase = paste0(name, '_ribo_p_cycle')
  make_pdf_targets(dfRes = res$results, dfFiltered = resAll$filtered, 
                   TFs_in=resAll$candidates$SYMBOL,
                   namebase=namebase)
  pdf_combine_and_delete(namebase)
})

res_ribo$C3_DD$results$TF %>% table %>% sort

plot(0:100, dbinom(0:100, size=25000, prob=4^-5), type='h')
plot(0:100, dbinom(0:100, size=25000, prob=4^-6), type='h')
plot(0:100, dbinom(0:100, size=25000, prob=4^-7), type='h')
# search for enrichment of RBPs in hits

plot(0:100, dbinom(0:100, size=100, prob=.3), type='h')

####    ####    ####    ####    ####    ####    ####    ####    ####    ####   



# 31. Rosbash agreement
################################################################################

pacemaker <- readxl::read_excel('rosbash_cyclers.xlsx', sheet=1) %>% 
  {rbind(., readxl::read_excel('rosbash_cyclers.xlsx', sheet=2, skip=1, col_names = names(.)))}

lapply(res10, function(res) {
  subset(pacemaker, Condition == 'DD' & Gene %in% res$filtered$Gene) %>%
    {table(.$Gene)} %>% sort %>%
    as.data.frame %>%
    {ggplot(., aes(Freq)) +
        geom_bar() +
        scale_x_binned(breaks = 1:max(.$Freq))}
}) %>% {plot_grid(plotlist=., labels=names(.))}

reconstituting_DFs_10$reconstituted %>%
  subset(select = c(Gene, Sample, Mode)) %>%
  merge(subset(pacemaker, select = c(Gene, Condition, cluster), Condition == 'LD')) %>%
  {table(.$Gene, .$Sample, .$Mode)} %>%
  as.data.frame %>% 
  `names<-`(c('Gene', 'Sample', 'Mode', 'Freq')) %>%
  arrange(desc(Freq)) %>%
  subset(Freq != 0)

lapply(res10, function(res) {
  subset(pacemaker, Gene %in% res$filtered$Gene & Condition == 'DD') %>%
    {table(.$Gene)} %>% sort})

# Trpm constitutively permissive calcium channel
'Trpm' %in% unique(c(ARG$ChR2_dTrpA1, ARG$dTrpA1_KCL, ARG$ChR2_KCL)) # NOPE
reconstituting_DFs_10$reconstituted %>% subset(Gene == 'Trpm') # ONLY C3
subset(pacemaker, Condition == 'DD' & Gene == 'Trpm') # COULD BE ALL

# cyclers in pacemaker clusters
subset(pacemaker, Gene %in% Cyclers$Gene) %>%
  {table(.$cluster)} %>% as.data.frame %>%
  `colnames<-`(c('cluster', 'Cyclers')) %>%
  arrange(desc(Cyclers)) %>%
  mutate(cluster = factor(cluster, levels=unique(cluster))) %>%
  ggplot(aes(y=Cyclers, x=cluster, fill=cluster)) +
  geom_bar(stat="identity") +
  theme(legend.position="none") +
  xlab('Pacemaker neuron cluster')

pacemaker %>%
  mutate(inC2C3 = Gene %in% subset(reconstituting_DFs_10$reconstituted, Mode == 'DD')$Gene)  %>% 
  {table(.$cluster, .$inC2C3)} %>% 
  as.data.frame %>%
  pivot_wider(names_from = 'Var2', values_from = 'Freq') %>% 
  `names<-`(c('cluster', 'NOT_IN', 'IN')) %>%
  mutate(pct_in = IN/(IN+NOT_IN)) %>%
  arrange(desc(pct_in)) %>%
  mutate(cluster = factor(cluster, levels=unique(cluster))) %>%
  ggplot(aes(y=pct_in, x=cluster, fill=cluster)) +
  geom_bar(stat="identity") +
  theme(legend.position="none") +
  xlab('Pacemaker neuron cluster') +
  ylab('Percentage C-cyclers in pacemaker, DD')

G = subset(reconstituting_DFs_10$reconstituted, Mode == 'DD')$Gene %>% unique
#G = res10$C2_LD$filtered$Gene
G = Cyclers$Gene %>% unique
pacemaker %>%
  mutate(inC2C3 = Gene %in% G)  %>% 
  {table(.$cluster, .$inC2C3)} %>% 
  as.data.frame %>%
  pivot_wider(names_from = 'Var2', values_from = 'Freq') %>% 
  `names<-`(c('cluster', 'NOT_IN', 'IN')) %>%
  mutate(pct_in = IN/(IN+NOT_IN)) %>%
  mutate(clkNO_cYES = length(G)-IN) %>%
  mutate(clkNO_cNO = length(rownames(data_C2C3)) - (IN+NOT_IN+clkNO_cYES)) %>%
  {data.frame(., fisher_pval = apply(., 1, function(L) {
    data.frame(c(L[['IN']], L[['NOT_IN']]), c(L[['clkNO_cYES']], L[['clkNO_cNO']])) %>% 
      mutate_all(as.numeric) %>% 
      fisher.test %>% .$p.value
  }))} %>% 
  mutate(fisher_q = p.adjust(fisher_pval, method='BH')) %>% 
  arrange(fisher_q) %>% 
  mutate(cluster = factor(cluster, levels=unique(cluster))) %>%
  mutate(neg_log_q = -log10(fisher_q)) %>% print %>%
  ggplot(aes(y=neg_log_q, x=NOT_IN+IN, label=cluster)) +
  geom_point() +
  geom_text()
  
  ggplot(aes(y=neg_log_q, x=cluster, fill=cluster)) +
  geom_bar(stat="identity") +
  theme(legend.position="none") +
  xlab('Pacemaker neuron cluster') +
  ylab('Enrichment (-log(q-value))') +
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed')


####    ####    ####    ####    ####    ####    ####    ####    ####    ####   



# 32. OTHER ENRICHMENT
################################################################################

# CLK/CWO ChIP-seq
WD_ChIP <- "/Users/Teddy/Desktop/CCG_distribution/cwo_ChIP.xlsx"

ChIP_Hardin <- list(CWO = rbind(readxl::read_excel(WD_ChIP, 
                                                   sheet = 1, col_names = TRUE) %>% cbind(ZT = 'ZT2'),
                                readxl::read_excel(WD_ChIP, 
                                                   sheet = 2, col_names = TRUE) %>% cbind(ZT = 'ZT14')),
                    CLK = rbind(readxl::read_excel(WD_ChIP, 
                                                   sheet = 3, col_names = TRUE) %>% cbind(ZT = 'ZT2'), 
                                readxl::read_excel(WD_ChIP, 
                                                   sheet = 4, col_names = TRUE) %>% cbind(ZT = 'ZT14'))) %>% # NOTE, RESULTS ARE GENE-MAPPED, AFTER REMOVAL OF INTERGENIC MAPPINGS; adult heads
  {c(., list(CWO_common_genes = .$CWO %>% 
               group_by(`Gene Name`) %>% filter(n()>1) %>% .[['Gene Name']] %>% unique, 
             CLK_common_genes = .$CLK %>%
               group_by(`Gene Name`) %>% filter(n()>1) %>% .[['Gene Name']] %>% unique,
             CWO_either_genes = .$CWO[['Gene Name']] %>% unique,
             CLK_either_genes = .$CLK[['Gene Name']] %>% unique, 
             both_either = intersect(.$CWO[['Gene Name']] %>% unique, .$CLK[['Gene Name']] %>% unique)))}

FIMO$Pdp1_C3LD$ChIP <- list()
FIMO$Pdp1_C3LD$filtered$Gene %>% {.[. %in% ChIP_Hardin$CLK_either_genes]} # CWO: 27 either, 18 common, CLK: 5 either, 0 common

sample(genes_C2C3, 235) %>% {.[. %in% ChIP_Hardin$CWO_either_genes]} %>% length # CWO: 8 either, 5 common, CLK: 2 either, 0 common

ChIP_shuffle_C3LD <- list()

ChIP_shuffle_C3LD$CWO_common <- lapply(1:50, function(x) sample(genes_C2C3, 235) %>% {.[. %in% ChIP_Hardin$CWO_common_genes]} %>% length) %>% unlist
ChIP_shuffle_C3LD$CWO_either <- lapply(1:50, function(x) sample(genes_C2C3, 235) %>% {.[. %in% ChIP_Hardin$CWO_either_genes]} %>% length) %>% unlist
ChIP_shuffle_C3LD$CLK_common <- lapply(1:50, function(x) sample(genes_C2C3, 235) %>% {.[. %in% ChIP_Hardin$CLK_common_genes]} %>% length) %>% unlist
ChIP_shuffle_C3LD$CLK_either <- lapply(1:50, function(x) sample(genes_C2C3, 235) %>% {.[. %in% ChIP_Hardin$CLK_either_genes]} %>% length) %>% unlist

lapply(ChIP_shuffle_C3LD, mean)


background = data_C2C3 %>% subset(sample_mode == 'C3_LD') %>%
  {rownames(.)[rowSums(.[["RNA"]]$counts)!=0]}

candidates = FIMO$Pdp1_C3LD$filtered$Gene %>% unique
non_candidates = background[background %!in% candidates]

process = ChIP_Hardin$CWO_common_genes %>% {.[. %in% background]}

to_fisher <- data.frame(in_list = c(which(candidates %in% process) %>% length, 
                                      which(candidates %!in% process) %>% length), 
                        not_in_list = c(which(non_candidates %in% process) %>% length,
                                        which(non_candidates %!in% process) %>% length), 
                        row.names = c('in_process', 'not_in_process'))

chisq.test(to_fisher)$expected # one is below 5
fisher.test(to_fisher)

hi <- run_fisher(candidates = FIMO$Pdp1_C3LD$filtered$Gene %>% unique, 
           background = data_C2C3 %>% subset(sample_mode == 'C3_LD') %>%
             {rownames(.)[rowSums(.[["RNA"]]$counts)!=0]})

lapply(levels(data_C2C3$sample_mode), function(S_M) {
  S = gsub('_.*', '', S_M)
  M = gsub('.*_', '', S_M)
  run_fisher_general(A = subset(reconstituting_DFs_10$reconstituted, Sample==S & Mode==M)$Gene, 
                     B = subset(to_filter, second_pct>0.1 & Gene %in% ChIP_Hardin$CWO_either_genes & Sample==S & Mode==M)$Gene, 
                     U = subset(to_filter, second_pct>0.1 & Sample==S & Mode==M)$Gene)
}) %>%
  `names<-`(levels(data_C2C3$sample_mode))

fisher_runs_amp2 <- lapply(filtered_amp2, function(df) {
  Sample_Mode = unique(df$Sample_Mode)
  Background = subset(data_C2C3, sample_mode == Sample_Mode) %>%
    {rownames(.)[rowSums(.[["RNA"]]$counts)!=0]}
  Candidates = unique(df$Gene)
  c(run_fisher(candidates=Candidates, background=Background), Sample_Mode=Sample_Mode)
})

run_fisher = function(candidates, background, process_list = ChIP_Hardin[c('CWO_common_genes', 'CWO_either_genes', 
                                                                'CLK_common_genes', 'CLK_either_genes')]) {
  res = list()
  non_candidates = background[background %!in% candidates]
  
  for (i in 1:length(process_list)) {
    process = process_list[i][[1]]
    res[[names(process_list[i])]]$df = data.frame(in_list = c(which(candidates %in% process) %>% length, 
                                                            which(candidates %!in% process) %>% length), 
                                                not_in_list = c(which(non_candidates %in% process) %>% length,
                                                                which(non_candidates %!in% process) %>% length), 
                                                row.names = c('in_process', 'not_in_process'))
    res[[names(process_list[i])]]$fisher = fisher.test(res[[names(process_list[i])]]$df)
    res[[names(process_list[i])]]$p = res[[names(process_list[i])]]$fisher$p.value
  }
  return(res)
}




####    ####    ####    ####    ####    ####    ####    ####    ####    ####   



# 33. MEGACELL BUSINESS
################################################################################
data_C2C3_normGroup1 <- readRDS('data_C2C3_normGroup1.rds')
data_C2C3_normGroup2 <- readRDS('data_C2C3_normGroup2.rds')

MC_norm_group2 <- AggregateExpression(data_C2C3, assays='RNA', group.by=c('sample_mode', 'group2'))[[1]] %>%
  {t(.)/colSums(.)} %>% t %>% as.data.frame

AggregateExpression(data_C2C3, assays='RNA', group.by=c('sample_mode', 'group2'))[[1]] %>%
  colSums

FetchData(data_C2C3, vars = c('per', 'group2', 'sample_mode'), layer = 'data') %>% 
  as.data.frame %>% 
  subset(group2 == 'ZT3' & sample_mode == 'C3_LD') %>%
  .$per %>% expm1 %>% mean 

AverageExpression(data_C2C3, features = 'per', group.by = c('sample_mode', 'group2'), layer = 'data') %>%
  .[[1]] %>% as.data.frame %>%
  .[,'C3-LD_ZT3']

avgs_normGroup2['per', 'C3-LD_ZT3']

avgs_normGroup2 <- readRDS('data_C2C3_normGroup2.rds') %>% 
  AverageExpression(group.by=c('sample_mode', 'group2'), layer = 'counts') %>% 
  .[[1]] %>% as.data.frame

MC_norm_group2 <- MC_norm_group2[rownames(avgs_normGroup2), colnames(avgs_normGroup2)]

avgs_normGroup2/MC_norm_group2

write.csv(avgs_normGroup2, 'avgs_normGroup2.csv')
write.csv(MC_norm_group2, 'MC_norm_group2.csv')


genes_DD_both <- to_filter %>%
  subset(meta2d_pvalue<0.2 & ratio>1.3 & ampl>0.7) %>%
  subset(Mode == 'DD') %>%
  group_by(Gene) %>%
  filter(n()>=2) %>%
  .$Gene %>% unique

#genes_display <- c(circ_genes, 'cyc', 'sr', 'Hr38')
genes_display <- genes_DD_both

plot_grid(nrow=1, rel_widths = c(0.55, 0.55, 0.2, 1, 0.2, 1, 0.2, 1),
          MC_norm_group2 %>% 
            {mutate(., Gene = rownames(.))} %>%
            pivot_longer(cols = -Gene, names_to = 'condition') %>%
            mutate(Sample_Mode = gsub('_.*', '', condition) %>% gsub('-', '_', .) %>%
                     factor(levels = c('C3_LD', 'C2_LD', 'C3_DD', 'C2_DD')),
                   timepoint = gsub('.*T', '', condition) %>% as.integer, .keep='unused') %>%
            subset(Gene %in% genes_display) %>%
            mutate(Gene = factor(Gene, levels=)) %>%
            subset(Sample_Mode %in% c('C3_LD', 'C2_LD')) %>%
            ggplot(aes(timepoint, value, color=Gene)) +
            geom_line() + 
            ylim(0, NA) +
            facet_grid(Gene ~ Sample_Mode, scales = 'free_y') +
            theme(axis.text.x = element_text(size = 6), axis.title = element_blank(), legend.position = "none"), 
          
          MC_norm_group2 %>% 
            {mutate(., Gene = rownames(.))} %>%
            pivot_longer(cols = -Gene, names_to = 'condition') %>%
            mutate(Sample_Mode = gsub('_.*', '', condition) %>% gsub('-', '_', .) %>%
                     factor(levels = c('C3_LD', 'C2_LD', 'C3_DD', 'C2_DD')),
                   timepoint = gsub('.*T', '', condition) %>% as.integer, .keep='unused') %>%
            subset(Gene %in% genes_display) %>%
            mutate(Gene = factor(Gene, levels=genes_display)) %>%
            subset(Sample_Mode %in% c('C3_DD', 'C2_DD')) %>%
            ggplot(aes(timepoint, value, color=Gene)) +
            geom_line() + 
            ylim(0, NA) +
            facet_grid(Gene ~ Sample_Mode, scales = 'free_y') +
            theme(axis.text.x = element_text(size = 6), axis.title = element_blank(), legend.position = "none"), 
          
          ggplot()+theme_minimal(),
          
          avg_plotter_concat_mode(genes_display, N_days = 2, rib_alpha = 0.4, data=data_C2C3_normGroup2, COUNTS = TRUE) +
            ggtitle('Gene-cell normalized averages'), 
          
          ggplot()+theme_minimal(),
            
          avg_plotter_concat_mode(genes_display, N_days = 2, rib_alpha = 0.4) +
            ggtitle('Cell-normalized averages'), 
          
          ggplot()+theme_minimal(),
            
          avg_plotter_repDot(genes_display))







####    ####    ####    ####    ####    ####    ####    ####    ####    ####   



# 34. Plotting good candidates, standard_plotter
################################################################################

# move this below
standard_plotter = function(Genes, filetag, Sample, Mode, data=data_C2C3_AGG, COUNTS=FALSE, HORIZ=TRUE, NO_Y_TEXT=FALSE, NO_TAB_Y=FALSE, BRIGHTEN_DD_DAY=TRUE, filetype='pdf', scale_factor=1.25) {
  
  D1 = scale_factor*length(Genes)+0.25
  D2 = scale_factor*1.25
  W = ifelse(HORIZ, D1+0.25*as.integer(!NO_Y_TEXT), D2+as.integer(!NO_Y_TEXT)*0.25)
  H = ifelse(HORIZ, D2, D1)
  
  avg_plotter_concat(Genes, data=subset(data, sample==Sample & mode1==Mode), COUNTS=COUNTS, 
                     NO_Y_TEXT=NO_Y_TEXT, SWAP_FACET=HORIZ, NO_X_TEXT=TRUE, N_days=2, 
                     NO_GRIDLINES=TRUE, NO_TAB_Y=NO_TAB_Y, BRIGHTEN_DD_DAY=BRIGHTEN_DD_DAY, BLACKOUT=TRUE) %>%
    ggsave(file=paste(filetag, filetype, sep='.'), plot=., width=W, height=H)
}



standard_plotter_2 = function(Gene, Sample, Mode, filetype='pdf', p_size=1, Point_size=1, Line_size=0.5, Margin_p_size_scaler=1/72,
                              P_DESC=TRUE, P_NODESC=TRUE, P_NOTHING=FALSE, data=data_C2C3_AGG, COUNTS=FALSE, HORIZ=FALSE, 
                              BRIGHTEN_DD_DAY=TRUE, filestart='') {
  
  Margin = ifelse( is.null(Margin_p_size_scaler), NULL, p_size * Margin_p_size_scaler)
  
  
  if (P_DESC) {
    avg_plotter_concat(Gene, data=subset(data, sample==Sample & mode1==Mode), COUNTS=COUNTS, 
                       NO_Y_TEXT=FALSE, SWAP_FACET=HORIZ, NO_X_TEXT=TRUE, N_days=2, 
                       NO_GRIDLINES=TRUE, NO_TAB_Y=FALSE, BRIGHTEN_DD_DAY=BRIGHTEN_DD_DAY, 
                       BLACKOUT=TRUE, point_size=Point_size, line_size=Line_size, 
                       strip_x_manual=Margin, strip_y_manual=Margin) %>%
      ggsave(file=paste0(filestart, Gene, '_', Sample, '_', Mode, '_', p_size, '_DESC.', filetype), plot=., width=p_size, height=p_size)
  }
  
  if (P_NODESC) {
    avg_plotter_concat(Gene, data=subset(data, sample==Sample & mode1==Mode), COUNTS=COUNTS, 
                       NO_Y_TEXT=TRUE, SWAP_FACET=HORIZ, NO_X_TEXT=TRUE, N_days=2, 
                       NO_GRIDLINES=TRUE, NO_TAB_Y=TRUE, BRIGHTEN_DD_DAY=BRIGHTEN_DD_DAY, 
                       BLACKOUT=TRUE, point_size=Point_size, line_size=Line_size, 
                       strip_x_manual=Margin) %>%
      ggsave(file=paste0(filestart, Gene, '_', Sample, '_', Mode, '_', p_size, '_NODESC.', filetype), plot=., width=p_size, height=p_size)
  }
  
  if (P_NOTHING) {
    avg_plotter_concat(Gene, data=subset(data, sample==Sample & mode1==Mode), COUNTS=COUNTS, 
                       NO_Y_TEXT=TRUE, SWAP_FACET=HORIZ, NO_X_TEXT=TRUE, N_days=2, 
                       NO_GRIDLINES=TRUE, NO_TAB_Y=TRUE, NO_TAB_X=TRUE, 
                       BRIGHTEN_DD_DAY=BRIGHTEN_DD_DAY, BLACKOUT=TRUE, point_size=Point_size, line_size=Line_size) %>%
      ggsave(file=paste0(filestart, Gene, '_', Sample, '_', Mode, '_', p_size, '_NOTHING.', filetype), plot=., width=p_size, height=p_size)
  }
}

standard_plotter_2('sr', 'C3', 'LD', P_NODESC=TRUE, P_NOTHING=TRUE, Point_size=1, p_size=1, Margin_p_size_scaler=1/72)
avg_plotter_concat('sr', data=subset(data_C2C3_AGG, sample=='C3' & mode1=='LD')) +
  theme(strip.text.x = element_text(margin = margin(1/6,0,1/6,0, "in")), 
        strip.text.y = element_text(margin = margin(0,1/6,0,1/6, "in")))


lapply(levels(data_C2C3$sample_mode), function(S_M) {
  S = gsub('_.*', '', S_M)
  M = gsub('.*_', '', S_M)
  standard_plotter_2('cwo', S, M, P_NODESC=TRUE, P_NOTHING=FALSE)
}
)

merge(GO_ETC, Cyclers) %>%
  merge(subset(reconstituting_DFs_10$reconstituted, Sample == 'C2' & Mode == 'LD')) %>%
  {.$Gene} %>%
  avg_plotter_concat(data=subset(data_C2C3_AGG, sample=='C2' & mode1=='LD'), N_days = 2)

lapply(c(
  'RpL13A','RpS11','RpL31','RpL6',
  'mt:CoI','mt:CoII','mt:Cyt-b','mt:ND6',
  'ND-30','ND-51','ND-B14.7','NP15.6'
), function(G) {
  standard_plotter_2(G, 'C2', 'LD', P_NODESC=TRUE, P_NOTHING=FALSE, P_DESC=FALSE)
}
)





standard_plotter(genes_activity, filetag='ARG_C3LD', Sample='C3', Mode='LD', HORIZ=FALSE, NO_Y_TEXT=FALSE)
standard_plotter(genes_activity, filetag='ARG_C3LD', Sample='C3', Mode='LD', HORIZ=FALSE, NO_Y_TEXT=FALSE, filetype='svg')
standard_plotter(genes_activity, filetag='ARG_C2LD', Sample='C2', Mode='LD', HORIZ=FALSE, NO_Y_TEXT=FALSE)
standard_plotter(genes_activity, filetag='ARG_C2LD', Sample='C2', Mode='LD', HORIZ=FALSE, NO_Y_TEXT=FALSE, filetype='svg')
standard_plotter(genes_activity, filetag='ARG_C3DD', Sample='C3', Mode='DD', HORIZ=FALSE, NO_Y_TEXT=FALSE)
standard_plotter(genes_activity, filetag='ARG_C3DD', Sample='C3', Mode='DD', HORIZ=FALSE, NO_Y_TEXT=FALSE, filetype='svg')

standard_plotter(c('shn', 'REPTOR-BP', 'CG6843'), filetag='TFs_C3LD', Sample='C3', Mode='LD', HORIZ=FALSE, NO_Y_TEXT=FALSE)
standard_plotter(c('shn', 'REPTOR-BP', 'CG6843'), filetag='TFs_C3LD', Sample='C3', Mode='LD', HORIZ=FALSE, NO_Y_TEXT=FALSE, filetype='svg')

standard_plotter(c('Atf3', 'ab', 'CG9932'), filetag='TFs_C3DD', Sample='C3', Mode='DD', HORIZ=FALSE, NO_Y_TEXT=FALSE)
standard_plotter(c('Atf3', 'ab', 'CG9932'), filetag='TFs_C3DD', Sample='C3', Mode='DD', HORIZ=FALSE, NO_Y_TEXT=FALSE, filetype='svg')

standard_plotter(genes_activity, filetag='TFs_C3DD', Sample='C3', Mode='DD', HORIZ=FALSE, NO_Y_TEXT=FALSE)
standard_plotter(genes_activity, filetag='TFs_C3DD', Sample='C3', Mode='DD', HORIZ=FALSE, NO_Y_TEXT=FALSE, filetype='svg')

standard_plotter(c('eIF3k', 'eIF3e', 'RpL30', 'mRpL40', 'ND-B14.5A', 'mt:CoII'), 
                 filetag='ofInterest_C3LD', Sample='C3', Mode='LD', HORIZ=FALSE, NO_Y_TEXT=FALSE)
standard_plotter(c('eIF3k', 'eIF3e', 'RpL30', 'mRpL40', 'ND-B14.5A', 'mt:CoII')[1:3], 
                 filetag='ofInterest_C3LD_1', Sample='C3', Mode='LD', HORIZ=FALSE, NO_Y_TEXT=FALSE)
standard_plotter(c('eIF3k', 'eIF3e', 'RpL30', 'mRpL40', 'ND-B14.5A', 'mt:CoII')[4:6], 
                 filetag='ofInterest_C3LD_2', Sample='C3', Mode='LD', HORIZ=FALSE, NO_Y_TEXT=FALSE)
standard_plotter(c('eIF3k', 'eIF3e', 'RpL30', 'mRpL40', 'ND-B14.5A', 'mt:CoII'), 
                 filetag='ofInterest_C3LD', Sample='C3', Mode='LD', HORIZ=FALSE, NO_Y_TEXT=FALSE, filetype='svg')
standard_plotter(c('eIF3k', 'eIF3e', 'RpL30', 'mRpL40', 'ND-B14.5A', 'mt:CoII')[1:3], 
                 filetag='ofInterest_C3LD_1', Sample='C3', Mode='LD', HORIZ=FALSE, NO_Y_TEXT=FALSE, filetype='svg')
standard_plotter(c('eIF3k', 'eIF3e', 'RpL30', 'mRpL40', 'ND-B14.5A', 'mt:CoII')[4:6], 
                 filetag='ofInterest_C3LD_2', Sample='C3', Mode='LD', HORIZ=FALSE, NO_Y_TEXT=FALSE, filetype='svg')

standard_plotter(c('eIF3a', 'eIF3g1', 'RpL31', 'mRpS35', 'ND-51', 'mt:CoI'), 
                 filetag='ofInterest_C2LD', Sample='C2', Mode='LD', HORIZ=FALSE, NO_Y_TEXT=FALSE)
standard_plotter(c('eIF3a', 'eIF3g1', 'RpL31', 'mRpS35', 'ND-51', 'mt:CoI')[1:3],
                 filetag='ofInterest_C2LD_1', Sample='C2', Mode='LD', HORIZ=FALSE, NO_Y_TEXT=FALSE)
standard_plotter(c('eIF3a', 'eIF3g1', 'RpL31', 'mRpS35', 'ND-51', 'mt:CoI')[4:6],
                 filetag='ofInterest_C2LD_2', Sample='C2', Mode='LD', HORIZ=FALSE, NO_Y_TEXT=FALSE)
standard_plotter(c('eIF3a', 'eIF3g1', 'RpL31', 'mRpS35', 'ND-51', 'mt:CoI'),
                 filetag='ofInterest_C2LD', Sample='C2', Mode='LD', HORIZ=FALSE, NO_Y_TEXT=FALSE, filetype='svg')
standard_plotter(c('eIF3a', 'eIF3g1', 'RpL31', 'mRpS35', 'ND-51', 'mt:CoI')[1:3], 
                 filetag='ofInterest_C2LD_1', Sample='C2', Mode='LD', HORIZ=FALSE, NO_Y_TEXT=FALSE, filetype='svg')
standard_plotter(c('eIF3a', 'eIF3g1', 'RpL31', 'mRpS35', 'ND-51', 'mt:CoI')[4:6],
                 filetag='ofInterest_C2LD_2', Sample='C2', Mode='LD', HORIZ=FALSE, NO_Y_TEXT=FALSE, filetype='svg')

standard_plotter(c('eIF3c', 'eIF3d1', 'RpL29', 'mRpL49', 'mRpS18C', 'SdhC'), 
                 filetag='ofInterest_C3DD', Sample='C3', Mode='DD', HORIZ=FALSE, NO_Y_TEXT=FALSE)
standard_plotter(c('eIF3c', 'eIF3d1', 'RpL29', 'mRpL49', 'mRpS18C', 'SdhC')[1:3],
                 filetag='ofInterest_C3DD_1', Sample='C3', Mode='DD', HORIZ=FALSE, NO_Y_TEXT=FALSE)
standard_plotter(c('eIF3c', 'eIF3d1', 'RpL29', 'mRpL49', 'mRpS18C', 'SdhC')[4:6], 
                 filetag='ofInterest_C3DD_2', Sample='C3', Mode='DD', HORIZ=FALSE, NO_Y_TEXT=FALSE)
standard_plotter(c('eIF3c', 'eIF3d1', 'RpL29', 'mRpL49', 'mRpS18C', 'SdhC'),
                 filetag='ofInterest_C3DD', Sample='C3', Mode='DD', HORIZ=FALSE, NO_Y_TEXT=FALSE, filetype='svg')
standard_plotter(c('eIF3c', 'eIF3d1', 'RpL29', 'mRpL49', 'mRpS18C', 'SdhC')[1:3],
                 filetag='ofInterest_C3DD_1', Sample='C3', Mode='DD', HORIZ=FALSE, NO_Y_TEXT=FALSE, filetype='svg')
standard_plotter(c('eIF3c', 'eIF3d1', 'RpL29', 'mRpL49', 'mRpS18C', 'SdhC')[4:6], 
                 filetag='ofInterest_C3DD_2', Sample='C3', Mode='DD', HORIZ=FALSE, NO_Y_TEXT=FALSE, filetype='svg')

global_print_type='pdf'

lapply(c('sr', 'CrebA', 'kay', 'Pka-C1', 'cwo', 'Rala', 'fru', 'mt:Cyt-b', genes_activity), function(G) {
  standard_plotter(G, filetag=paste0(G, '_C3LD_noYtext'), Sample='C3', Mode='LD', HORIZ=TRUE, NO_Y_TEXT=TRUE, NO_TAB_Y=TRUE, filetype=global_print_type, data=data_C2C3_AGG_max2, COUNTS=TRUE, scale_factor=1)})

lapply(c('sr', 'CrebA'), function(G) {
  standard_plotter(G, filetag=paste0(G, '_C3LD_YLAB'), Sample='C3', Mode='LD', HORIZ=TRUE, NO_Y_TEXT=FALSE, NO_TAB_Y=TRUE, filetype=global_print_type, data=data_C2C3_AGG_max2, COUNTS=TRUE, scale_factor=1)})

lapply(c('fru', 'sr'), function(G) {
  standard_plotter(G, filetag=paste0(G, '_C3LD_YLAB'), Sample='C3', Mode='LD', HORIZ=TRUE, NO_Y_TEXT=FALSE, NO_TAB_Y=TRUE, filetype=global_print_type, scale_factor=1)})

lapply(c('eIF3k', 'eIF3e', 'RpL30', 'mRpL40', 'ND-B14.5A', 'mt:CoII'), function(G) {
  standard_plotter(G, filetag=paste0(G, '_C3LD_noYtext'), Sample='C3', Mode='LD', HORIZ=TRUE, NO_Y_TEXT=TRUE, NO_TAB_Y=TRUE, filetype=global_print_type, data=data_C2C3_AGG_max2, COUNTS=TRUE, scale_factor=1)})

lapply(c('mRpL49', 'CG9932', 'eIF3c', 'ab', 'CG34250', 'Atf3', 'pho', 'BtbVII'), function(G) {
  standard_plotter(G, filetag=paste0(G, '_C3DD_noYtext'), Sample='C3', Mode='DD', HORIZ=TRUE, NO_Y_TEXT=TRUE, NO_TAB_Y=TRUE, filetype=global_print_type, data=data_C2C3_AGG_max2, COUNTS=TRUE, scale_factor=1)})

lapply(c('mRpL49', 'CG9932'), function(G) {
  standard_plotter(G, filetag=paste0(G, '_C3DD_YLAB'), Sample='C3', Mode='DD', HORIZ=TRUE, NO_Y_TEXT=FALSE, NO_TAB_Y=TRUE, filetype=global_print_type, data=data_C2C3_AGG_max2, COUNTS=TRUE, scale_factor=1)})
lapply(c('pho', 'BtbVII'), function(G) {
  standard_plotter(G, filetag=paste0(G, '_C3DD_YLAB'), Sample='C3', Mode='DD', HORIZ=TRUE, NO_Y_TEXT=FALSE, NO_TAB_Y=TRUE, filetype=global_print_type, scale_factor=1)})


{for (S in c('C2', 'C3')) {
  for (M in c('LD', 'DD')) {
    for (ext in c('pdf', 'svg', 'png')) {
      standard_plotter('RasGAP1', filetag=paste('RasGAP1', S, M, sep='_'), Sample=S, Mode=M, HORIZ=TRUE, NO_Y_TEXT=FALSE, NO_TAB_Y=FALSE, filetype=ext, data=data_C2C3_AGG, COUNTS=FALSE, scale_factor=1.25)
    }
  }
}
  rm(S, M, ext)
}

merge(GO_ETC, Cyclers)$Gene %>% {.[. %in% res10$C3_LD$filtered$Gene]} %>%
  avg_plotter_concat(data=subset(data_C2C3_AGG_max1, sample=='C3' & mode1=='LD'), COUNTS=TRUE, NO_Y_TEXT=FALSE, FACET_2D=TRUE, NO_X_TEXT=TRUE, N_days=2, NO_GRIDLINES=TRUE, LEGEND=TRUE)
merge(GO_ribo_cell, Cyclers)$Gene %>% {.[. %in% res10$C3_LD$filtered$Gene]} %>%
  plotter_geneCombine(Mode='LD', Sample='C3', NO_X_TEXT=TRUE, offset=0, data=data_C2C3_AGG_max1, ALIGN=TRUE)

merge(GO_ETC, passQC)$Gene %>% {.[. %in% res10$C3_LD$filtered$Gene]} %>%
  avg_plotter_concat(data=subset(data_C2C3_AGG, sample=='C3' & mode1=='LD'), COUNTS=FALSE, NO_Y_TEXT=FALSE, FACET_2D=FALSE, NO_X_TEXT=TRUE, N_days=2, NO_GRIDLINES=TRUE, LEGEND=TRUE)

####    ####    ####    ####    ####    ####    ####    ####    ####    ####    



# 35. Paper plots:
################################################################################

# 5C LD
c('CG14186', 'Hr38', 'cwo') %>%
  lapply(function(G) {
    standard_plotter_2(G, 'C3', 'LD', P_NODESC=TRUE, P_NOTHING=FALSE)
  })

# 5C DD
c('CG34250', 'BtbVII', 'eIF3c', 'pho', 'mRpL49') %>%
  lapply(function(G) {
    standard_plotter_2(G, 'C3', 'DD', P_NODESC=TRUE, P_NOTHING=FALSE, Point_size=1.5, Line_size = 2)
  })


# 5C LD alt
c('CrebA', 'mRpL48', 'RpL27', 'mt:CoIII') %>%
  lapply(function(G) {
    standard_plotter_2(G, 'C3', 'LD', P_NODESC=TRUE, P_NOTHING=FALSE)
  })

# ADDED OPTIONS, LD:
c("Adk1","alc","alpha-Man-Ib","CG14657","CG17839","CG1998","CG31510","CG31717","CG32700","CG3652","CG43066",
  "CG4577","CG46385","CG7083","CrebA","Dys","Fer1HCH","Fer2LCH","fru","hig","Ilk","Kul","lncRNA:CR42861","lncRNA:CR45054",
  "meng","mRpL48","mt:Cyt-b","mt:CoIII","Nuak1","Oscillin","path","Rala","RapGAP1","Snap24","sr","Trpm") %>%
  lapply(function(G) {
    standard_plotter_2(G, 'C3', 'LD', P_NODESC=TRUE, P_NOTHING=FALSE)
  })
# ADDED OPTIONS, DD:
c("bsf","CG30392","CG34250","CG41099","CG4660","CG4658","CG7139","DCTN2-p50","Fis1","hng3","Nna1","Rab4","Rrp47","Sema2a","sick","Synj","TER94") %>%
  lapply(function(G) {
    standard_plotter_2(G, 'C3', 'DD', P_NODESC=TRUE, P_NOTHING=FALSE)
  })

# PER REQUEST
c('CG14186', 'Hr38', 'sr') %>%
  lapply(function(G) {
    standard_plotter_2(G, 'C3', 'LD', P_NODESC=TRUE, P_NOTHING=FALSE)
  })


plot_grid(
plotter_geneCombine(subset(reconstituting_DFs_10$reconstituted, Sample=='C3' & Mode=='DD')$Gene %>% {.[grep('^RpL', .)]},
                    data=data_C2C3_AGG_01_1d, Mode='DD', Sample='C3', 
                    NO_X_TEXT=TRUE, NO_Y_TEXT=FALSE, NO_GRIDLINES = TRUE, PRINT_GENES = FALSE, LEGEND=TRUE, offset = -0.5, ALIGN=TRUE, 
                    gene_linesize = 0.5, gene_alpha = 0.7, mean_linesize=0.5, rib_alpha=0.8) +
  #theme(strip.text.x = element_text(margin = margin(1/72,0,1/72,0, "in"))) +
  theme(strip.text.y = element_text(margin = margin(0,1/72,0,1/72, "in"))),

plotter_geneCombine(subset(reconstituting_DFs_10$reconstituted, Sample=='C3' & Mode=='DD')$Gene %>% {.[grep('^RpL', .)]},
                    data=data_C2C3_AGG_01_1d, Mode='DD', Sample='C3', 
                    NO_X_TEXT=TRUE, NO_Y_TEXT=FALSE, NO_GRIDLINES = TRUE, PRINT_GENES = FALSE, LEGEND=TRUE, offset = -0.5, ALIGN=TRUE, 
                    gene_linesize = 0.5, gene_alpha = 0.7, mean_linesize=0.5, rib_alpha=0.8, PRISM=TRUE) +
  #theme(strip.text.x = element_text(margin = margin(1/72,0,1/72,0, "in"))) +
  theme(strip.text.y = element_text(margin = margin(0,1/72,0,1/72, "in")))
)

# 5C COMB PLOTS
Subset_Stats$subsets[c('GO_ribo_noPre', 'GO_ribo_cell', 'GO_ribo_mito', 'ETC_mi', 'ETC_NOmi', 'ARGs', 'ARGs_main', 'GPCR_cyc', 'GCPR_4', 'GPCR_5', 'GPCR_6')] %>%
  {lapply(names(.), function(name) {
    Ps = subset(reconstituting_DFs_10$reconstituted, Gene %in% .[[name]]) %>% 
      {split(., list(.$Sample, .$Mode))} %>% 
      lapply(function(x) {if(nrow(x)) {
        plotter_geneCombine(x$Gene, data=data_C2C3_AGG_01_1d, Mode=x[1,'Mode'], Sample=x[1,'Sample'], 
                            NO_X_TEXT=TRUE, NO_Y_TEXT=TRUE, NO_GRIDLINES = TRUE, PRINT_GENES = TRUE, LEGEND=FALSE, offset = -0.5, ALIGN=TRUE, 
                            gene_linesize = 0.5, gene_alpha = 0.7, mean_linesize=0.5, rib_alpha=0.8) +
          #theme(strip.text.x = element_text(margin = margin(1/72,0,1/72,0, "in"))) +
          theme(strip.text.y = element_text(margin = margin(0,1/72,0,1/72, "in")))
      } else return(NULL)
      }) 
    
    Ps %>% rev %>%
      {lapply(names(.), function(S_M_name) {
        #lapply(c('.svg', '.pdf', '.jpg'), function(extension) ggsave(paste0('Comb_', name, '_', S_M_name, '_CLEAN', extension), plot=.[[S_M_name]], width=1.6,height=1.6))
        #ggsave(paste0('Comb_', name, '_', S_M_name, '_CLEAN.pdf'), plot=.[[S_M_name]], width=1.6,height=1.6)
        ggsave(paste0('Comb_', name, '_', S_M_name, '_CLEAN_small.pdf'), plot=.[[S_M_name]], width=1.3,height=1.3)
      }
      )}
  })}

Subset_Stats$subsets[c('GO_ribo_noPre', 'GO_ribo_cell', 'GO_ribo_mito', 'ETC_mi', 'ETC_NOmi', 'ARGs', 'ARGs_main', 'GPCR_cyc', 'GCPR_4', 'GPCR_5', 'GPCR_6')] %>%
  {lapply(names(.), function(name) {
    Pa = subset(reconstituting_DFs_10$reconstituted, Gene %in% .[[name]]) %>% 
      {split(., list(.$Sample, .$Mode))} %>% 
      lapply(function(x) {if(nrow(x)) {
        plotter_geneCombine(x$Gene, data=data_C2C3_AGG_01_1d, Mode=x[1,'Mode'], Sample=x[1,'Sample'], 
                            NO_X_TEXT=TRUE, NO_Y_TEXT=FALSE, NO_GRIDLINES = TRUE, PRINT_GENES = TRUE, LEGEND=TRUE, offset = -0.5, ALIGN=TRUE, 
                            gene_linesize = 0.5, gene_alpha = 0.7, mean_linesize=0.5, rib_alpha=0.8) +
          #theme(strip.text.x = element_text(margin = margin(1/72,0,1/72,0, "in"))) +
          theme(strip.text.y = element_text(margin = margin(0,1/72,0,1/72, "in")))
      } else return(NULL)
      })
    
    Pa %>% rev %>%
      {lapply(names(.), function(S_M_name) {
        #lapply(c('.svg', '.pdf', '.jpg'), function(extension) ggsave(paste0('Comb_', name, '_', S_M_name, '_CLEAN', extension), plot=.[[S_M_name]], width=1.6,height=1.6))
        #ggsave(paste0('Comb_', name, '_', S_M_name, '_NOTES.pdf'), plot=.[[S_M_name]], width=1.6,height=1.6)
        ggsave(paste0('Comb_', name, '_', S_M_name, '_NOTES_small.pdf'), plot=.[[S_M_name]], width=1.3,height=1.3)
      }
      )}
  })}

subset(to_filter, Gene %in% c('Hr38', 'sr', 'CG14186')) %>% 
  {split(., list(.$Sample, .$Mode))} %>% 
  lapply(function(x) {if(nrow(x)) {
    plotter_geneCombine(x$Gene, data=data_C2C3_AGG_01_1d, Mode=x[1,'Mode'], Sample=x[1,'Sample'], 
                        NO_X_TEXT=TRUE, NO_Y_TEXT=TRUE, NO_GRIDLINES = TRUE, PRINT_GENES = TRUE, LEGEND=FALSE, offset = -0.5, ALIGN=TRUE,
                        gene_linesize = 0.5, gene_alpha = 0.7, mean_linesize=0.5, rib_alpha=0.8) +
      #theme(strip.text.x = element_text(margin = margin(1/72,0,1/72,0, "in"))) +
      theme(strip.text.y = element_text(margin = margin(0,1/72,0,1/72, "in")))
  } else return(NULL)
  }) %>%
  {lapply(names(.), function(S_M_name) {
    ggsave(paste0('Comb_3ARGs_NoGenes', S_M_name, '.pdf'), plot=.[[S_M_name]], width=1.6,height=1.6)
  }
  )}



# 
lapply(a_res10, function(L2) {
  DF_cyc = L2$df %>%
    mutate(Sample = gsub('_.*', '', cond), 
           Mode = gsub('.*_', '', cond)) %>%
    merge(subset(reconstituting_DFs_10$reconstituted, select=c(Gene, Sample, Mode)) %>%
            mutate(TF = Gene, .keep='unused'))
  print(nrow(DF_cyc))
  DF_cyc %>%
    mutate(FDR_fisher = p.adjust(as.numeric(p_fisher), method='BH')) %>%
    subset(FDR_fisher < 0.1) %>% print %>%
    {if (!nrow(.)) ggplot() + theme_void()
      else 
        apply(., 1, function(df) {
          G = df[['TF']][[1]]
          C = df[['cluster']][[1]]
          cond = df[['cond']][[1]]
          S = gsub('_.*', '', cond)
          M = gsub('.*_', '', cond)
          targets = a_res10[[cond]][[G]][[C]]$targets
          
          print(paste(G, C))
          a_res10[[cond]][[G]][[C]]$phases %>% print
          
          #standard_plotter_2(G, S, M, P_DESC=TRUE, P_NODESC=TRUE, filetype='pdf')
          
          plotter_geneCombine(targets, data=data_C2C3_AGG_01_1d, Mode=M, Sample=S, 
                              NO_X_TEXT=TRUE, NO_Y_TEXT=TRUE, NO_GRIDLINES = TRUE, PRINT_GENES = FALSE, LEGEND=FALSE, offset = -0.5, ALIGN=TRUE,
                              gene_linesize = 0.5, gene_alpha = 0.7, mean_linesize=0.5, rib_alpha=0.8) %>%
            #ggsave(filename=paste0('TargetsComb_', G, '_', S, M, '.pdf'), plot=., width=1.6,height=1.6)
            ggsave(filename=paste0('TargetsComb_', G, '_', S, M, '_geneSize.pdf'), plot=., width=1,height=1)
          #ggsave(filename=paste0('TargetsComb_', G, '_', S, M, '_small.pdf'), plot=., width=1.3,height=1.3)
          
          #make_pdf_targets(subset(res10[[cond]]$results, TF == G), res10[[cond]]$filtered, namebase=paste('TF', G, paste0(S, M), sep='_'), Cols=color_list_alt)
          
          #make_pdf_targets(subset(res10[[cond]]$results, TF == G & Target %in% targets), res10[[cond]]$filtered, namebase=paste('TF', G, paste0(S, M), 'clst', C, sep='_'), Cols=color_list_alt)
          
        }) 
    }
})

# ab

res10$C3_DD$hclustered$ab$p
res10$C3_DD$hclustered$ab$res_df %>%
  {split(., .$cluster)} %>%
  lapply(function(df) {
    C = df[1, 'cluster']
    targets = df$Gene
    
    plotter_geneCombine(targets, data=data_C2C3_AGG_01_1d, Mode='DD', Sample='C3', 
                        NO_X_TEXT=TRUE, NO_Y_TEXT=TRUE, NO_GRIDLINES = TRUE, PRINT_GENES = TRUE, LEGEND=FALSE, offset = -0.5, ALIGN=TRUE,
                        gene_linesize = 0.5, gene_alpha = 0.7, mean_linesize=0.5, rib_alpha=0.8) %>%
      ggsave(filename=paste0('TargetsCombCluster_', C, '_ab_C3_DD.pdf'), plot=., width=1.6,height=1.6)
  })

res10$C3_DD$hclustered$ab$res_df %>%
  list %>%
  lapply(function(df) {
    targets = df$Gene
    
    plotter_geneCombine(targets, data=data_C2C3_AGG_01_1d, Mode='DD', Sample='C3', 
                        NO_X_TEXT=TRUE, NO_Y_TEXT=TRUE, NO_GRIDLINES = TRUE, PRINT_GENES = TRUE, LEGEND=FALSE, offset = -0.5, ALIGN=TRUE,
                        gene_linesize = 0.5, gene_alpha = 0.7, mean_linesize=0.5, rib_alpha=0.8) %>%
      #ggsave(filename='TargetsCombCluster_ab_C3_DD.pdf', plot=., width=1.6,height=1.6)
      ggsave(filename='TargetsCombCluster_ab_C3_DD_small.pdf', plot=., width=1.3,height=1.3)
    #ggsave(filename='TargetsCombCluster_ab_C3_DD_genesize.pdf', plot=., width=1,height=1)
  })

hist_phase(subset(reconstituting_DFs_10$reconstituted, Sample=='C3'&Mode=='DD'& Gene %in% res10$C3_DD$hclustered$ab$res_df$Gene), 
           CIRCULAR=FALSE, TRIM_TO_24=TRUE, binWidth=1, n_W=0, Colors=color_list_alt, LD = FALSE, kernel_adjust = 0) %>%
  ggsave(filename='phase_targets_ab_C3DD_w1.pdf', plot=., width=2,height=1)

hist_phase(subset(reconstituting_DFs_10$reconstituted, Sample=='C3'&Mode=='DD'), 
           CIRCULAR=FALSE, TRIM_TO_24=TRUE, binWidth=1, n_W=0, Colors=color_list_alt, LD = FALSE, kernel_adjust = 0) %>%
  ggsave(filename='phase_C3DD_w1.pdf', plot=., width=2,height=1)



# twi
a_res10$C3_LD$df %>% subset(TF == 'twi')

hist_phase(subset(reconstituting_DFs_10$reconstituted, Sample=='C3'&Mode=='LD'& Gene %in% res10$C3_LD$hclustered$twi$res_df$Gene), 
           CIRCULAR=FALSE, TRIM_TO_24=TRUE, binWidth=1, n_W=0, Colors=color_list_alt, LD = FALSE, kernel_adjust = 0) %>%
  ggsave(filename='phase_targets_twi_C3LD_w1.pdf', plot=., width=2,height=1)
hist_phase(subset(reconstituting_DFs_10$reconstituted, Sample=='C3'&Mode=='LD'), 
           CIRCULAR=FALSE, TRIM_TO_24=TRUE, binWidth=1, n_W=0, Colors=color_list_alt, LD = FALSE, kernel_adjust = 0) %>%
  ggsave(filename='phase_C3LD_w01.pdf', plot=., width=2,height=1)

hist_phase(subset(reconstituting_DFs_10$reconstituted, Sample=='C3'&Mode=='LD'& Gene %in% res10$C3_LD$hclustered$twi$res_df$Gene), 
           CIRCULAR=FALSE, TRIM_TO_24=TRUE, binWidth=0.5, n_W=0, Colors=color_list_alt, LD = FALSE, kernel_adjust = 0) %>%
  ggsave(filename='phase_targets_twi_C3LD_w05.pdf', plot=., width=2,height=1)
hist_phase(subset(reconstituting_DFs_10$reconstituted, Sample=='C3'&Mode=='LD'), 
           CIRCULAR=FALSE, TRIM_TO_24=TRUE, binWidth=0.5, n_W=0, Colors=color_list_alt, LD = FALSE, kernel_adjust = 0) %>%
  ggsave(filename='phase_C3LD_w05.pdf', plot=., width=2,height=1)


res10$C3_LD$hclustered$twi$res_df %>%
  {split(., .$cluster)} %>%
  lapply(function(df) {
    C = df[1, 'cluster']
    targets = df$Gene
    
    plotter_geneCombine(targets, data=data_C2C3_AGG_01_1d, Mode='LD', Sample='C3', 
                        NO_X_TEXT=TRUE, NO_Y_TEXT=TRUE, NO_GRIDLINES = TRUE, PRINT_GENES = TRUE, LEGEND=FALSE, offset = -0.5, ALIGN=TRUE,
                        gene_linesize = 0.5, gene_alpha = 0.7, mean_linesize=0.5, rib_alpha=0.8) %>%
      ggsave(filename=paste0('TargetsCombCluster_', C, '_twi_C3_LD_small.pdf'), plot=., width=1,height=1)
  })

res10$C3_LD$hclustered$twi$res_df %>%
  list %>%
  lapply(function(df) {
    targets = df$Gene
    
    plotter_geneCombine(targets, data=data_C2C3_AGG_01_1d, Mode='LD', Sample='C3', 
                        NO_X_TEXT=TRUE, NO_Y_TEXT=TRUE, NO_GRIDLINES = TRUE, PRINT_GENES = TRUE, LEGEND=FALSE, offset = -0.5, ALIGN=TRUE,
                        gene_linesize = 0.5, gene_alpha = 0.7, mean_linesize=0.5, rib_alpha=0.8) %>%
      #ggsave(filename='TargetsCombCluster_twi_C3_LD.pdf', plot=., width=1.6,height=1.6)
      ggsave(filename='TargetsCombCluster_twi_C3_LD_small.pdf', plot=., width=1.3,height=1.3)
    #ggsave(filename='TargetsCombCluster_twi_C3_LD_genesize.pdf', plot=., width=1,height=1)
  })


# ALTERNATIVE PLOTTING STYLES not used
plotter_geneCombine(merge(subset(reconstituting_DFs_10$reconstituted, Sample=='C3'&Mode=='DD'), merge(GOs$GO_ribo_cell, Cyclers))$Gene, data=data_C2C3_AGG_01_1d, Mode='DD', Sample='C3',
                    NO_X_TEXT=TRUE, NO_GRIDLINES = TRUE, PRINT_GENES = TRUE, LEGEND=FALSE, NO_Y_TEXT=TRUE, offset = -0.5, ALIGN=TRUE, gene_linesize = 0.5, mean_linesize=0, rib_alpha=0, point_size=0)

hmap_phase(df_genes = subset(to_filter, Gene %in% merge(Cyclers, GO_ribo_cell)$Gene & Sample=='C3' & Mode=='LD'))

Subset_Stats$subsets[c('GO_ribo_noPre', 'GO_ribo_cell', 'GO_ribo_mito', 'ETC_mi', 'ETC_NOmi', 'ARGs', 'ARGs_main', 'GPCR_cyc')] %>%
  lapply(function(Genes) {
    plotter_geneCombine(subset(reconstituting_DFs_10$reconstituted, Sample=='C3' & Mode=='LD' & Gene %in% Genes)$Gene %>% print, data=data_C2C3_AGG_01_1d, Mode='LD', Sample='C3',
                        NO_X_TEXT=TRUE, NO_GRIDLINES = TRUE, PRINT_GENES = TRUE, LEGEND=FALSE, NO_Y_TEXT=TRUE, offset = -0.5, ALIGN=TRUE, gene_linesize = 0.5, gene_alpha = 0.6, mean_linesize=0.5, rib_alpha=1, point_size=0)
  }) %>%
  plot_grid(plotlist=.)

Subset_Stats$subsets[c('GO_ribo_cell', 'GO_ribo_mito', 'ETC_mi', 'ETC_NOmi', 'ARGs', 'ARGs_main', 'GPCR_cyc')] %>%
  lapply(function(Genes) {
    hmap_phase(subset(reconstituting_DFs_10$reconstituted, Sample=='C3' & Mode=='LD' & Gene %in% Genes)$Gene)
  }) %>%
  plot_grid(plotlist=.)

Subset_Stats$subsets[c('GO_ribo_cell', 'GO_ribo_mito', 'ETC_mi', 'ARGs', 'GPCR_cyc')] %>%
  lapply(function(Genes) {
    df = subset(reconstituting_DFs_10$reconstituted, Sample=='C3' & Mode=='LD' & Gene %in% Genes)
    list(p = hmap_phase(df, sample_remove='C2', mode_remove='DD', col_low = 'white', COLOR_GENES = TRUE, background_fill='white'), 
         l = nrow(df))
  }) %>% 
  {plot_grid(plotlist= lapply(., function(x) x$p), ncol=1, rel_heights= unlist(lapply(., function(x) x$l)))}


Subset_Stats$subsets[c('GO_ribo_cell', 'GO_ribo_mito', 'ETC_mi', 'ARGs', 'GPCR_cyc')] %>%
  lapply(function(Genes) {
    df = subset(reconstituting_DFs_10$reconstituted, Sample=='C3' & Mode=='LD' & Gene %in% Genes)
    list(p = hmap_phase(df, sample_remove='C2', mode_remove='DD', col_low = 'white', COLOR_GENES = FALSE), 
         l = nrow(df))
  }) %>% 
  {plot_grid(plotlist= lapply(., function(x) x$p), ncol=1, rel_heights= unlist(lapply(., function(x) x$l)))}



tempfun=function(M) {
  plot_grid(ncol=2,
            plotter_geneCombine(subset(reconstituting_DFs_10$reconstituted, Sample=='C3' & Mode==M)$Gene %>% {.[grep('^Rp', .)]},
                                data=data_C2C3_AGG_01_1d, Mode=M, Sample='C3', 
                                NO_X_TEXT=TRUE, NO_Y_TEXT=FALSE, NO_GRIDLINES = TRUE, PRINT_GENES = TRUE, LEGEND=FALSE, offset = -0.5, ALIGN=TRUE, 
                                gene_linesize = 0.5, gene_alpha = 0.7, mean_linesize=0.5, rib_alpha=0.8) +
              #theme(strip.text.x = element_text(margin = margin(1/72,0,1/72,0, "in"))) +
              theme(strip.text.y = element_text(margin = margin(0,1/72,0,1/72, "in"))),
            plotter_geneCombine(subset(reconstituting_DFs_10$reconstituted, Sample=='C3' & Mode==M)$Gene %>% {.[grep('^Rp', .)]},
                                data=data_C2C3_AGG_01_1d, Mode=M, Sample='C3', 
                                NO_X_TEXT=TRUE, NO_Y_TEXT=FALSE, NO_GRIDLINES = TRUE, PRINT_GENES = TRUE, LEGEND=FALSE, offset = -0.5, ALIGN=TRUE, 
                                gene_linesize = 0.5, gene_alpha = 0.7, mean_linesize=0.5, rib_alpha=0.8, PRISM=TRUE, GENECOLOR=FALSE) +
              #theme(strip.text.x = element_text(margin = margin(1/72,0,1/72,0, "in"))) +
              theme(strip.text.y = element_text(margin = margin(0,1/72,0,1/72, "in"))), 
            
            
            plotter_geneCombine(subset(reconstituting_DFs_10$reconstituted, Sample=='C3' & Mode==M)$Gene %>% {.[grep('^Rp', .)]},
                                data=data_C2C3_AGG_01_1d, Mode=M, Sample='C3', 
                                NO_X_TEXT=TRUE, NO_Y_TEXT=FALSE, NO_GRIDLINES = TRUE, PRINT_GENES = TRUE, LEGEND=FALSE, offset = -0.5, ALIGN=TRUE, 
                                gene_linesize = 0.2, gene_alpha = 0.7, mean_linesize=0.5, rib_alpha=0.8) +
              #theme(strip.text.x = element_text(margin = margin(1/72,0,1/72,0, "in"))) +
              theme(strip.text.y = element_text(margin = margin(0,1/72,0,1/72, "in"))),
            plotter_geneCombine(subset(reconstituting_DFs_10$reconstituted, Sample=='C3' & Mode==M)$Gene %>% {.[grep('^Rp', .)]},
                                data=data_C2C3_AGG_01_1d, Mode=M, Sample='C3', 
                                NO_X_TEXT=TRUE, NO_Y_TEXT=FALSE, NO_GRIDLINES = TRUE, PRINT_GENES = TRUE, LEGEND=FALSE, offset = -0.5, ALIGN=TRUE, 
                                gene_linesize = 0.2, gene_alpha = 0.7, mean_linesize=0.5, rib_alpha=0.8, PRISM=TRUE, GENECOLOR=FALSE) +
              #theme(strip.text.x = element_text(margin = margin(1/72,0,1/72,0, "in"))) +
              theme(strip.text.y = element_text(margin = margin(0,1/72,0,1/72, "in"))), 
            
            plotter_geneCombine(subset(reconstituting_DFs_10$reconstituted, Sample=='C3' & Mode==M)$Gene %>% {.[grep('^Rp', .)]},
                                data=data_C2C3_AGG_01_1d, Mode=M, Sample='C3', 
                                NO_X_TEXT=TRUE, NO_Y_TEXT=FALSE, NO_GRIDLINES = TRUE, PRINT_GENES = FALSE, LEGEND=FALSE, offset = -0.5, ALIGN=TRUE, 
                                gene_linesize = 0.5, gene_alpha = 0.7, mean_linesize=0.5, rib_alpha=0.8) +
              #theme(strip.text.x = element_text(margin = margin(1/72,0,1/72,0, "in"))) +
              theme(strip.text.y = element_text(margin = margin(0,1/72,0,1/72, "in"))),
            plotter_geneCombine(subset(reconstituting_DFs_10$reconstituted, Sample=='C3' & Mode==M)$Gene %>% {.[grep('^Rp', .)]},
                                data=data_C2C3_AGG_01_1d, Mode=M, Sample='C3', 
                                NO_X_TEXT=TRUE, NO_Y_TEXT=FALSE, NO_GRIDLINES = TRUE, PRINT_GENES = FALSE, LEGEND=FALSE, offset = -0.5, ALIGN=TRUE, 
                                gene_linesize = 0.5, gene_alpha = 0.7, mean_linesize=0.5, rib_alpha=0.8, PRISM=TRUE, GENECOLOR=FALSE) +
              #theme(strip.text.x = element_text(margin = margin(1/72,0,1/72,0, "in"))) +
              theme(strip.text.y = element_text(margin = margin(0,1/72,0,1/72, "in")))
  )
}

plot_grid(nrow=1, rel_widths=c(1, 0.2, 1), 
          tempfun('LD'), 
          ggplot() + theme_void(),
          tempfun('DD'))

####    ####    ####    ####    ####    ####    ####    ####    ####    ####   


# 36. PRINTING TABLES
################################################################################

# 4 tabs – passing genes 			      include expression threshold
createWorkbook() %T>% 
  addWorksheet('C2_LD') %T>% 
  writeData('C2_LD', rowNames = FALSE,
            subset(to_filter, Sample=='C2' & Mode=='LD' & 
                     second_pct>0.1, select=c(Gene, Sample, Mode, second_pct))
            ) %T>% 
  addWorksheet('C3_LD') %T>% 
  writeData('C3_LD', rowNames = FALSE,
            subset(to_filter, Sample=='C3' & Mode=='LD' & 
                     second_pct>0.1, select=c(Gene, Sample, Mode, second_pct))
  ) %T>% 
  addWorksheet('C2_DD') %T>% 
  writeData('C2_DD', rowNames = FALSE,
            subset(to_filter, Sample=='C2' & Mode=='DD' & 
                     second_pct>0.1, select=c(Gene, Sample, Mode, second_pct))
  ) %T>% 
  addWorksheet('C3_DD') %T>% 
  writeData('C3_DD', rowNames = FALSE,
            subset(to_filter, Sample=='C3' & Mode=='DD' & 
                     second_pct>0.1, select=c(Gene, Sample, Mode, second_pct))
  ) %>% 
  saveWorkbook('genesExpressed.xlsx', overwrite=TRUE)


# 4 tabs – cycling genes			      include cycling parameters
createWorkbook() %T>% 
  addWorksheet('C2_LD') %T>% 
  writeData('C2_LD', rowNames = FALSE,
            subset(reconstituting_DFs_10$reconstituted, Sample=='C2' & Mode=='LD',
                   select=c(Gene, Sample, Mode, second_pct, 
                            ratio, ratio_AGG, meta2d_pvalue, meta2d_phase, meta2d_pvalue_AGG_norm, meta2d_phase_AGG_norm))
  ) %T>% 
  addWorksheet('C3_LD') %T>% 
  writeData('C3_LD', rowNames = FALSE,
            subset(reconstituting_DFs_10$reconstituted, Sample=='C3' & Mode=='LD',
                   select=c(Gene, Sample, Mode, second_pct, 
                            ratio, ratio_AGG, meta2d_pvalue, meta2d_phase, meta2d_pvalue_AGG_norm, meta2d_phase_AGG_norm))
  ) %T>% 
  addWorksheet('C2_DD') %T>% 
  writeData('C2_DD', rowNames = FALSE,
            subset(reconstituting_DFs_10$reconstituted, Sample=='C2' & Mode=='DD',
                   select=c(Gene, Sample, Mode, second_pct, 
                            ratio, ratio_AGG, meta2d_pvalue, meta2d_phase, meta2d_pvalue_AGG_norm, meta2d_phase_AGG_norm))
  ) %T>% 
  addWorksheet('C3_DD') %T>% 
  writeData('C3_DD', rowNames = FALSE,
            subset(reconstituting_DFs_10$reconstituted, Sample=='C3' & Mode=='DD',
                   select=c(Gene, Sample, Mode, second_pct, 
                            ratio, ratio_AGG, meta2d_pvalue, meta2d_phase, meta2d_pvalue_AGG_norm, meta2d_phase_AGG_norm))
  ) %>% 
  saveWorkbook('genesCycling.xlsx', overwrite=TRUE)


# 4 tabs – FIMO					            FIMO output
createWorkbook() %T>% 
  addWorksheet('C2_LD') %T>% 
  writeData('C2_LD', rowNames = FALSE,
            res10$C2_LD$results
  ) %T>% 
  addWorksheet('C3_LD') %T>% 
  writeData('C3_LD', rowNames = FALSE,
            res10$C3_LD$results
  ) %T>% 
  addWorksheet('C2_DD') %T>% 
  writeData('C2_DD', rowNames = FALSE,
            res10$C2_DD$results
  ) %T>% 
  addWorksheet('C3_DD') %T>% 
  writeData('C3_DD', rowNames = FALSE,
            res10$C3_DD$results
  ) %>% 
  saveWorkbook('FIMO_results.xlsx', overwrite=TRUE)

res10$C2_DD$results %>% head
as.data.frame(res10$C2_DD$promoters) %>% head

hi <- mergeByOverlaps(query = res10$C2_DD$fimo, subject = res10$C2_DD$promoters) %>%
  {data.frame(.$`res10$C2_DD$fimo`, .$gene_id)}


# 4 tabs – FIMO analysis					  Cycling TFs/in Fimo/pass
createWorkbook() %T>% 
  addWorksheet('C2_LD') %T>% 
  writeData('C2_LD', rowNames = FALSE,
            subset(a_res10_df, cond == 'C2_LD')
  ) %T>% 
  addWorksheet('C3_LD') %T>% 
  writeData('C3_LD', rowNames = FALSE,
            subset(a_res10_df, cond == 'C3_LD')
  ) %T>% 
  addWorksheet('C2_DD') %T>% 
  writeData('C2_DD', rowNames = FALSE,
            subset(a_res10_df, cond == 'C2_DD')
  ) %T>% 
  addWorksheet('C3_DD') %T>% 
  writeData('C3_DD', rowNames = FALSE,
            subset(a_res10_df, cond == 'C3_DD')
  ) %>% 
  saveWorkbook('FIMO_analysis.xlsx', overwrite=TRUE)


# 9 – 4 tabs – GO					        Gene/Category-mother-term expressed/cycle
GOdf_expr <- Subset_Stats$subsets %>%
  {.[c('ARGs', 'GO_ribo_cell', 'GO_ribo_mito', 'ETC_mi', 'ETC_NOmi')]} %>%
  stack %>%
  `colnames<-`(c('Gene', 'Group')) %>%
  merge(data.frame(Group = c('ARGs', 'GO_ribo_cell', 'GO_ribo_mito', 'ETC_mi', 'ETC_NOmi'),
                   Origin = c('Chen,2016', 'GOCC:0022626', 'GOCC:0005761', 'GOBP:0022900', 'GOBP:0022900'))) %>%
  merge(subset(to_filter, select=c('Gene', 'Sample', 'Mode', 'second_pct'))) %>%
  subset(second_pct > 0.1, select=-second_pct) %>%
  mutate(Cycle = (paste0(Gene, Sample, Mode) %in% paste0(reconstituting_DFs_10$reconstituted$Gene,
                                                        reconstituting_DFs_10$reconstituted$Sample,
                                                        reconstituting_DFs_10$reconstituted$Mode) ))

createWorkbook() %T>% 
  addWorksheet('C2_LD') %T>% 
  writeData('C2_LD', rowNames = FALSE,
            subset(GOdf_expr, Sample == 'C2' & Mode == 'LD')
  ) %T>% 
  addWorksheet('C3_LD') %T>% 
  writeData('C3_LD', rowNames = FALSE,
            subset(GOdf_expr, Sample == 'C2' & Mode == 'LD')
  ) %T>% 
  addWorksheet('C2_DD') %T>% 
  writeData('C2_DD', rowNames = FALSE,
            subset(GOdf_expr, Sample == 'C2' & Mode == 'LD')
  ) %T>% 
  addWorksheet('C3_DD') %T>% 
  writeData('C3_DD', rowNames = FALSE,
            subset(GOdf_expr, Sample == 'C2' & Mode == 'LD')
  ) %>% 
  saveWorkbook('FIMO_analysis.xlsx', overwrite=TRUE)




# normalizedData
#   dataAGG - Aggregate expression of all genes per sample/timepoint (normalized)
#   dataAVG - Average expression of all genes per sample/timepoint (normalized)
createWorkbook() %T>% 
  addWorksheet('Aggregated, normalized') %T>% 
  writeData('Aggregated, normalized', rowNames=FALSE, 
            avg_data_times2 %>% subset(select=-sample_mode_time) %>%
              mutate(timepoint = as.integer(timepoint)) %>%
              arrange(timepoint) %>% pivot_wider(names_from='timepoint', values_from='avg')
  ) %T>%
  addWorksheet('Averaged, normalized') %T>% 
  writeData('Averaged, normalized', rowNames=FALSE, 
            avg_data_times2_AGG %>% subset(select=-sample_mode_time) %>%
              mutate(timepoint = as.integer(timepoint)) %>%
              arrange(timepoint) %>% pivot_wider(names_from='timepoint', values_from='avg')
  ) %>% 
  saveWorkbook('normalizedData.xlsx', overwrite=TRUE)


# metadata
createWorkbook() %T>% 
  addWorksheet('metadata') %T>% 
  writeData('metadata', 
              subset(data_C2C3@meta.data, 
                     select = c('nCount_RNA', 'nFeature_RNA', 'Perc_MT', 'is_sng', 'dgrp1', 'mode1',
                                'group1', 'group2', 'time1', 'seurat_clusters', 'sample')), rowNames = TRUE) %>% 
  saveWorkbook('metadata.xlsx', overwrite=TRUE)

# counts MTX
createWorkbook() %T>% 
  addWorksheet('counts') %T>% 
  writeData('counts', data_C2C3@assays$RNA$counts, rowNames = TRUE) %>% 
  saveWorkbook('counts.xlsx', overwrite=TRUE)




detach("package:openxlsx", unload=TRUE)

library(openxlsx)

####    ####    ####    ####    ####    ####    ####    ####    ####    ####   


plotter_geneCombine(merge(GO_ribo_cell, passQC)$Gene %>% {. [. %in% subset(HMQ, Sample=='C3')$Gene]},
                    data=data_C2C3_AGG_01_1d, Sample='C3', dgrp_colname='time4', Mode=c('LD', 'DD'), N_days=4, COUNTS=TRUE,
                    NO_X_TEXT=TRUE, NO_Y_TEXT=FALSE, NO_GRIDLINES = TRUE, PRINT_GENES = FALSE, LEGEND=FALSE, offset = -0.5, ALIGN=FALSE, 
                    gene_linesize = 0.5, gene_alpha = 0.7, mean_linesize=0.5, rib_alpha=0.8, PRISM=TRUE, GENECOLOR=FALSE) +
  #theme(strip.text.x = element_text(margin = margin(1/72,0,1/72,0, "in"))) +
  theme(strip.text.y = element_text(margin = margin(0,1/72,0,1/72, "in")))



plotter_geneCombine( unlist(ARG$min2[1:48])[unlist(ARG$min2[1:48]) %in% subset(HMQ, Sample=='C3')$Gene],
                    data=data_C2C3_AGG_01_1d, Sample='C3', dgrp_colname='time4', Mode=c('LD', 'DD'), N_days=4, COUNTS=TRUE,
                    NO_X_TEXT=TRUE, NO_Y_TEXT=FALSE, NO_GRIDLINES = TRUE, PRINT_GENES = FALSE, LEGEND=FALSE, offset = -0.5, ALIGN=FALSE, 
                    gene_linesize = 0.5, gene_alpha = 0.7, mean_linesize=0.5, rib_alpha=0.8, PRISM=TRUE, GENECOLOR=FALSE) +
  #theme(strip.text.x = element_text(margin = margin(1/72,0,1/72,0, "in"))) +
  theme(strip.text.y = element_text(margin = margin(0,1/72,0,1/72, "in")))


plotter_geneCombine(merge(GO_ribo_cell, passQC)$Gene %>% {. [. %in% subset(HMQ_4, Sample=='C3')$Gene]},
                    data=data_C2C3_AGG_01_1d, Sample='C3', dgrp_colname='time4', Mode=c('LD', 'DD'), N_days=4, COUNTS=TRUE,
                    NO_X_TEXT=TRUE, NO_Y_TEXT=FALSE, NO_GRIDLINES = TRUE, PRINT_GENES = FALSE, LEGEND=FALSE, offset = -0.5, ALIGN=FALSE, 
                    gene_linesize = 0.5, gene_alpha = 0.7, mean_linesize=0.5, rib_alpha=0.8, PRISM=TRUE, GENECOLOR=FALSE) +
  #theme(strip.text.x = element_text(margin = margin(1/72,0,1/72,0, "in"))) +
  theme(strip.text.y = element_text(margin = margin(0,1/72,0,1/72, "in")))



plotter_geneCombine( unlist(ARG$min2[1:48])[unlist(ARG$min2[1:48]) %in% subset(Filtereds$DFs$reps_adding_pool, Sample=='C3')$Gene],
                     data=avg_reps_adding_pool_AGG_01_1d, Sample='C3', dgrp_colname='time3', Mode=c('LD', 'DD'), N_days=3, COUNTS=TRUE,
                     NO_X_TEXT=TRUE, NO_Y_TEXT=FALSE, NO_GRIDLINES = TRUE, PRINT_GENES = FALSE, LEGEND=FALSE, offset = -0.5, ALIGN=FALSE, 
                     gene_linesize = 0.5, gene_alpha = 0.7, mean_linesize=0.5, rib_alpha=0.8, PRISM=TRUE, GENECOLOR=FALSE) +
  #theme(strip.text.x = element_text(margin = margin(1/72,0,1/72,0, "in"))) +
  theme(strip.text.y = element_text(margin = margin(0,1/72,0,1/72, "in")))


# 37. HMQ PLOTS
################################################################################
Subset_Stats$r_means_HMQ <- lapply(Subset_Stats$subsets, function(Set) {
  lapply(levels(data_C2C3$sample_mode), function(S_M) {
    S = gsub('_.*', '', S_M)
    M = gsub('.*_', '', S_M)
    data.frame(Gene = Set) %>%
      make_corrplot(GO_anno = passQC, filtered=merge(to_filter, HMQ[, c('Gene', 'Sample', 'Mode')]), sample=S, mode=M) %>%
      from_corrplot_meanmaker
  }) %>%
    do.call(what='cbind') %>%
    `colnames<-`(levels(data_C2C3$sample_mode)) 
}) %>%
  do.call(what='rbind') %>%
  `rownames<-`(names(Subset_Stats$subsets))

Subset_Stats$p_vals_HMQ <- lapply(Subset_Stats$subsets, function(Set) {
  lapply(levels(data_C2C3$sample_mode), function(S_M) {
    S = gsub('_.*', '', S_M)
    M = gsub('.*_', '', S_M)
    run_fisher_general(B = subset(HMQ, Sample == S & Mode == M)$Gene,
                       A = subset(to_filter, Gene %in% Set & second_pct > 0.10 & Sample == S & Mode == M)$Gene, 
                       U = subset(to_filter, second_pct > 0.10 & Sample == S & Mode == M)$Gene)
  }) %>%
    do.call(what='cbind') %>%
    `colnames<-`(levels(data_C2C3$sample_mode))
}) %>%
  do.call(what='rbind') %>%
  `rownames<-`(names(Subset_Stats$subsets))

# 5D
c('CrebA', 'Rala', 'mt:CoIII', 'mRpL48', 'fru', 'lncRNA:CR42861') %>% {.[. %!in% subset(HMQ, Sample == 'C3')$Gene]}
c('Hr38', 'sr', 'CG14186') %>% {.[. %!in% subset(HMQ, Sample == 'C3')$Gene]}
c('CG34250', 'CG41099', 'DGTN2-p50', 'Rrp47', 'hng3', 'mRpL49', 'CG7139', 'Fis1') %>% {.[. %!in% subset(HMQ, Sample == 'C3')$Gene]}


# 5C COMB PLOTS
Subset_Stats$subsets[c('GO_ribo_noPre', 'GO_ribo_cell', 'ARGs')] %>%
  {lapply(names(.), function(name) {
    Ps = subset(HMQ, Gene %in% .[[name]] & Sample == 'C3') %>% 
      {split(., list(.$Mode))} %>% 
      lapply(function(x) {if(nrow(x)) {
        plotter_geneCombine(x$Gene, data=data_C2C3_AGG_01_1d, Mode=x[1,'Mode'], Sample='C3', 
                            NO_X_TEXT=TRUE, NO_Y_TEXT=TRUE, NO_GRIDLINES = TRUE, PRINT_GENES = TRUE, LEGEND=FALSE, offset = -0.5, ALIGN=TRUE, 
                            gene_linesize = 0.2, gene_alpha = 0.7, mean_linesize=0.5, rib_alpha=0.8, PRISM=TRUE, GENECOLOR=FALSE) +
          #theme(strip.text.x = element_text(margin = margin(1/72,0,1/72,0, "in"))) +
          theme(strip.text.y = element_text(margin = margin(0,1/72,0,1/72, "in")))
      } else return(NULL)
      }) 
    
    Ps %>% rev %>%
      {lapply(names(.), function(S_M_name) {
        #lapply(c('.svg', '.pdf', '.jpg'), function(extension) ggsave(paste0('Comb_', name, '_', S_M_name, '_CLEAN', extension), plot=.[[S_M_name]], width=1.6,height=1.6))
        #ggsave(paste0('Comb_', name, '_', S_M_name, '_CLEAN.pdf'), plot=.[[S_M_name]], width=1.6,height=1.6)
        ggsave(paste0('Comb_', name, '_', S_M_name, '_CLEAN_small.pdf'), plot=.[[S_M_name]], width=1.3,height=1.45)
      }
      )}
  })}

Subset_Stats$subsets[c('GO_ribo_noPre', 'GO_ribo_cell', 'ARGs')] %>%
  {lapply(names(.), function(name) {
    Pa = subset(HMQ, Gene %in% .[[name]] & Sample == 'C3') %>% 
      {split(., list(.$Sample, .$Mode))} %>% 
      lapply(function(x) {if(nrow(x)) {
        plotter_geneCombine(x$Gene, data=data_C2C3_AGG_01_1d, Mode=x[1,'Mode'], Sample='C3', 
                            NO_X_TEXT=TRUE, NO_Y_TEXT=FALSE, NO_GRIDLINES = TRUE, PRINT_GENES = TRUE, LEGEND=FALSE, offset = -0.5, ALIGN=TRUE, 
                            gene_linesize = 0.2, gene_alpha = 0.7, mean_linesize=0.5, rib_alpha=0.8, PRISM=TRUE, GENECOLOR=FALSE) +
          #theme(strip.text.x = element_text(margin = margin(1/72,0,1/72,0, "in"))) +
          theme(strip.text.y = element_text(margin = margin(0,1/72,0,1/72, "in")))
      } else return(NULL)
      })
    
    Pa %>% rev %>%
      {lapply(names(.), function(S_M_name) {
        #lapply(c('.svg', '.pdf', '.jpg'), function(extension) ggsave(paste0('Comb_', name, '_', S_M_name, '_CLEAN', extension), plot=.[[S_M_name]], width=1.6,height=1.6))
        #ggsave(paste0('Comb_', name, '_', S_M_name, '_NOTES.pdf'), plot=.[[S_M_name]], width=1.6,height=1.6)
        ggsave(paste0('Comb_', name, '_', S_M_name, '_NOTES_small.pdf'), plot=.[[S_M_name]], width=1.8,height=1.5)
      }
      )}
  })}


# JUST LD
Subset_Stats$subsets[c('GO_ribo_mito', 'ETC_mi', 'ETC_NOmi')] %>%
  {lapply(names(.), function(name) {
    Ps = subset(HMQ, Gene %in% .[[name]] & Sample == 'C3' & Mode == 'LD') %>% 
      {split(., list(.$Sample, .$Mode))} %>% 
      lapply(function(x) {if(nrow(x)) {
        plotter_geneCombine(x$Gene, data=data_C2C3_AGG_01_1d, Mode=x[1,'Mode'], Sample=x[1,'Sample'], 
                            NO_X_TEXT=TRUE, NO_Y_TEXT=TRUE, NO_GRIDLINES = TRUE, PRINT_GENES = TRUE, LEGEND=FALSE, offset = -0.5, ALIGN=TRUE, 
                            gene_linesize = 0.2, gene_alpha = 0.7, mean_linesize=0.5, rib_alpha=0.8, PRISM=TRUE, GENECOLOR=FALSE) +
          #theme(strip.text.x = element_text(margin = margin(1/72,0,1/72,0, "in"))) +
          theme(strip.text.y = element_text(margin = margin(0,1/72,0,1/72, "in")))
      } else return(NULL)
      }) 
    
    Ps %>% rev %>%
      {lapply(names(.), function(S_M_name) {
        #lapply(c('.svg', '.pdf', '.jpg'), function(extension) ggsave(paste0('Comb_', name, '_', S_M_name, '_CLEAN', extension), plot=.[[S_M_name]], width=1.6,height=1.6))
        #ggsave(paste0('Comb_', name, '_', S_M_name, '_CLEAN.pdf'), plot=.[[S_M_name]], width=1.6,height=1.6)
        ggsave(paste0('Comb_', name, '_', S_M_name, '_CLEAN_small.pdf'), plot=.[[S_M_name]], width=1.3,height=1.45)
      }
      )}
  })}

Subset_Stats$subsets[c('GO_ribo_mito', 'ETC_mi', 'ETC_NOmi')] %>%
  {lapply(names(.), function(name) {
    Pa = subset(HMQ, Gene %in% .[[name]] & Sample == 'C3' & Mode == 'LD') %>% 
      {split(., list(.$Sample, .$Mode))} %>% 
      lapply(function(x) {if(nrow(x)) {
        plotter_geneCombine(x$Gene, data=data_C2C3_AGG_01_1d, Mode=x[1,'Mode'], Sample=x[1,'Sample'], 
                            NO_X_TEXT=TRUE, NO_Y_TEXT=FALSE, NO_GRIDLINES = TRUE, PRINT_GENES = TRUE, LEGEND=FALSE, offset = -0.5, ALIGN=TRUE, 
                            gene_linesize = 0.2, gene_alpha = 0.7, mean_linesize=0.5, rib_alpha=0.8, PRISM=TRUE, GENECOLOR=FALSE) +
          #theme(strip.text.x = element_text(margin = margin(1/72,0,1/72,0, "in"))) +
          theme(strip.text.y = element_text(margin = margin(0,1/72,0,1/72, "in")))
      } else return(NULL)
      })
    
    Pa %>% rev %>%
      {lapply(names(.), function(S_M_name) {
        #lapply(c('.svg', '.pdf', '.jpg'), function(extension) ggsave(paste0('Comb_', name, '_', S_M_name, '_CLEAN', extension), plot=.[[S_M_name]], width=1.6,height=1.6))
        #ggsave(paste0('Comb_', name, '_', S_M_name, '_NOTES.pdf'), plot=.[[S_M_name]], width=1.6,height=1.6)
        ggsave(paste0('Comb_', name, '_', S_M_name, '_NOTES_small.pdf'), plot=.[[S_M_name]], width=1.75,height=1.45)
      }
      )}
  })}




plot_grid(nrow=1, 
          hist_phase(merge(to_filter, HMQ[, c('Gene', 'Sample', 'Mode')]) , CIRCULAR=FALSE, TRIM_TO_24=TRUE, binWidth=2, n_W=0, YSCALE=TRUE, Colors=rep('blue', 4) %>% `names<-`(c('C3_LD', 'C3_DD', 'C2_LD', 'C2_DD')), HIST_FRONT=TRUE, LD=TRUE), 
          hist_phase(merge(to_filter, HMQ[, c('Gene', 'Sample', 'Mode')]) , CIRCULAR=TRUE, TRIM_TO_24=TRUE, binWidth=2, n_W=0, Colors=rep('blue', 4) %>% `names<-`(c('C3_LD', 'C3_DD', 'C2_LD', 'C2_DD')), HIST_FRONT=TRUE, LD=TRUE), 
          hmap_phase(merge(to_filter, HMQ[, c('Gene', 'Sample', 'Mode')]) , Colors = list(C3_LD="blue", C3_DD="blue", C2_LD="blue", C2_DD="blue"), col_low = "red")) %>%
  {lapply(c('.pdf'), function(ext)  ggsave(file=paste0('HMQ_phases_LD', ext), plot=., width=5, height=5))}

plot_grid(nrow=1, 
          hist_phase(merge(to_filter, HMQ[, c('Gene', 'Sample', 'Mode')]) , CIRCULAR=FALSE, TRIM_TO_24=TRUE, binWidth=2, n_W=0, YSCALE=TRUE, Colors=rep('blue', 4) %>% `names<-`(c('C3_LD', 'C3_DD', 'C2_LD', 'C2_DD')), HIST_FRONT=TRUE, LD=FALSE), 
          hist_phase(merge(to_filter, HMQ[, c('Gene', 'Sample', 'Mode')]) , CIRCULAR=TRUE, TRIM_TO_24=TRUE, binWidth=2, n_W=0, Colors=rep('blue', 4) %>% `names<-`(c('C3_LD', 'C3_DD', 'C2_LD', 'C2_DD')), HIST_FRONT=TRUE, LD=FALSE), 
          hmap_phase(merge(to_filter, HMQ[, c('Gene', 'Sample', 'Mode')]) , Colors = list(C3_LD="blue", C3_DD="blue", C2_LD="blue", C2_DD="blue"), col_low = "red")) %>%
  {lapply(c('.pdf'), function(ext)  ggsave(file=paste0('HMQ_phases_DD', ext), plot=., width=5, height=5))}
####    ####    ####    ####    ####    ####    ####    ####    ####    ####   


# 38. time3 PLOTS
################################################################################


BH_procedure(subset(to_filter, second_pct > 0.1 & ratio>1.5), Qcol=NULL, Pcol='meta2d_pvalue_avg_3', Qname='avg3_BH.Q', Q=0.05) %>% {table(.$Sample, .$Mode)}


Subset_Stats$r_means_time3 <- lapply(Subset_Stats$subsets, function(Set) {
  lapply(levels(data_C2C3$sample_mode), function(S_M) {
    S = gsub('_.*', '', S_M)
    M = gsub('.*_', '', S_M)
    data.frame(Gene = Set) %>%
      make_corrplot(GO_anno = passQC, filtered=BH_procedure(subset(to_filter, second_pct > 0.1 & ratio>1.5), Qcol=NULL, Pcol='meta2d_pvalue_avg_3', Qname='avg3_BH.Q', Q=0.05), sample=S, mode=M) %>%
      from_corrplot_meanmaker
  }) %>%
    do.call(what='cbind') %>%
    `colnames<-`(levels(data_C2C3$sample_mode)) 
}) %>%
  do.call(what='rbind') %>%
  `rownames<-`(names(Subset_Stats$subsets))

Subset_Stats$p_vals_time3 <- lapply(Subset_Stats$subsets, function(Set) {
  lapply(levels(data_C2C3$sample_mode), function(S_M) {
    S = gsub('_.*', '', S_M)
    M = gsub('.*_', '', S_M)
    run_fisher_general(B = subset(BH_procedure(subset(to_filter, second_pct > 0.1 & ratio>1.5), Qcol=NULL, Pcol='meta2d_pvalue_avg_3', Qname='avg3_BH.Q', Q=0.05), Sample == S & Mode == M)$Gene,
                       A = subset(to_filter, Gene %in% Set & second_pct > 0.10 & Sample == S & Mode == M)$Gene, 
                       U = subset(to_filter, second_pct > 0.10 & Sample == S & Mode == M)$Gene)
  }) %>%
    do.call(what='cbind') %>%
    `colnames<-`(levels(data_C2C3$sample_mode))
}) %>%
  do.call(what='rbind') %>%
  `rownames<-`(names(Subset_Stats$subsets))


BH_procedure(subset(to_filter, second_pct > 0.1 & ratio>1.5), Qcol=NULL, Pcol='meta2d_pvalue_avg_3', Qname='avg3_BH.Q', Q=0.01) %>%
  {table(.$Sample, .$Mode)}

BH_procedure(subset(to_filter, second_pct > 0.1 & ratio>1.5), Qcol=NULL, Pcol='meta2d_pvalue_avg_3', Qname='avg3_BH.Q', Q=0.05) %>%
  {split(., list(.$Sample, .$Mode))} %>% 
  lapply(function(DF) {
    S_M = paste0(DF$Sample[1], '_', DF$Mode[1])
    Tp = gsub('_.*', '', S_M)
    M = gsub('.*_', '', S_M)
    
    split_list(DF$Gene, 10) %>%
      lapply(function(G) make_plot(G, Group='time2', Type0 = Tp, Mode0 = M)) %>%
      plot_grid(plotlist = ., nrow=1) %>%
      ggsave(filename = paste0('time3/Cyclers_Q05_avg3_', S_M, '.pdf'), width=30, height=length(DF$Gene)/10, limitsize=FALSE)
  } )


filtered_t3Q5 <- BH_procedure(subset(to_filter, second_pct > 0.1 & ratio>1.5), Qcol=NULL, Pcol='meta2d_pvalue_avg_3', Qname='avg3_BH.Q', Q=0.05)


    

# set printing to run; make barplots, then print all Rp's and grade

# phase plots filtered_t3Q5
plot_grid(nrow=1, 
          hist_phase(merge(to_filter, filtered_t3Q5[, c('Gene', 'Sample', 'Mode')]) , CIRCULAR=FALSE, TRIM_TO_24=TRUE, binWidth=2, n_W=0, YSCALE=TRUE, Colors=rep('blue', 4) %>% `names<-`(c('C3_LD', 'C3_DD', 'C2_LD', 'C2_DD')), HIST_FRONT=TRUE, LD=TRUE), 
          hist_phase(merge(to_filter, filtered_t3Q5[, c('Gene', 'Sample', 'Mode')]) , CIRCULAR=TRUE, TRIM_TO_24=TRUE, binWidth=2, n_W=0, Colors=rep('blue', 4) %>% `names<-`(c('C3_LD', 'C3_DD', 'C2_LD', 'C2_DD')), HIST_FRONT=TRUE, LD=TRUE), 
          hmap_phase(merge(to_filter, filtered_t3Q5[, c('Gene', 'Sample', 'Mode')]) , Colors = list(C3_LD="blue", C3_DD="blue", C2_LD="blue", C2_DD="blue"), col_low = "red")) %>%
  {lapply(c('.pdf'), function(ext)  ggsave(file=paste0('HMQ_phases_LD', ext), plot=., width=5, height=5))}

plot_grid(nrow=1, 
          hist_phase(merge(to_filter, filtered_t3Q5[, c('Gene', 'Sample', 'Mode')]) , CIRCULAR=FALSE, TRIM_TO_24=TRUE, binWidth=2, n_W=0, YSCALE=TRUE, Colors=rep('blue', 4) %>% `names<-`(c('C3_LD', 'C3_DD', 'C2_LD', 'C2_DD')), HIST_FRONT=TRUE, LD=FALSE), 
          hist_phase(merge(to_filter, filtered_t3Q5[, c('Gene', 'Sample', 'Mode')]) , CIRCULAR=TRUE, TRIM_TO_24=TRUE, binWidth=2, n_W=0, Colors=rep('blue', 4) %>% `names<-`(c('C3_LD', 'C3_DD', 'C2_LD', 'C2_DD')), HIST_FRONT=TRUE, LD=FALSE), 
          hmap_phase(merge(to_filter, filtered_t3Q5[, c('Gene', 'Sample', 'Mode')]) , Colors = list(C3_LD="blue", C3_DD="blue", C2_LD="blue", C2_DD="blue"), col_low = "red")) %>%
  {lapply(c('.pdf'), function(ext)  ggsave(file=paste0('HMQ_phases_DD', ext), plot=., width=5, height=5))}




####    ####    ####    ####    ####    ####    ####    ####    ####    ####   


# 40. TALK PLOTS, making Filtereds
################################################################################
#Filtereds <- list()

Filtereds$DFs <- list(
  
  expressed = subset(to_filter, second_pct > 0.1, select=c(Gene, Sample, Mode, meta2d_pvalue, meta2d_phase)),
  
  init = reconstituting_DFs_10$reconstituted, 
  
  avgP05 = subset(to_filter, second_pct>0.1 & ratio>1.5 & meta2d_pvalue < 0.05, select=c(Gene, Sample, Mode, meta2d_pvalue, meta2d_phase)), 
  
  AGGP05 = subset(to_filter, second_pct>0.1 & ratio_AGG>1.5 & meta2d_pvalue_AGG_norm < 0.05, select=c(Gene, Sample, Mode, meta2d_pvalue_AGG_norm, meta2d_phase_AGG_norm, meta2d_phase)),
  
  avgQ1 = subset(to_filter, second_pct>0.1 & ratio>1.5, select=c(Gene, Sample, Mode, meta2d_pvalue, meta2d_phase)) %>%
    BH_procedure(Qcol=NULL, Pcol='meta2d_pvalue', Q=0.1),
  
  AGGQ1 = subset(to_filter, second_pct>0.1 & ratio_AGG>1.5, select=c(Gene, Sample, Mode, meta2d_pvalue_AGG_norm, meta2d_phase_AGG_norm, meta2d_phase)) %>%
    BH_procedure(Qcol=NULL, Pcol='meta2d_pvalue_AGG_norm', Q=0.1), 
  
  HMP = subset(to_filter, second_pct>0.1 & ratio>1.5 & ratio_AGG>1.5, select=c(Gene, Sample, Mode, meta2d_pvalue, meta2d_pvalue_AGG_norm, meta2d_phase)) %>%
    cbind(., HMP = apply(., 1, function(L) hmp.stat(c(L[['meta2d_pvalue']], L[['meta2d_pvalue_AGG_norm']], w=NULL)))) %>% 
    subset(HMP<0.05), 
  
  HMQ1 = subset(to_filter, second_pct>0.1 & ratio>1.5 & ratio_AGG>1.5, select=c(Gene, Sample, Mode, meta2d_pvalue, meta2d_pvalue_AGG_norm, meta2d_phase)) %>%
    cbind(., HMP = apply(., 1, function(L) hmp.stat(c(L[['meta2d_pvalue']], L[['meta2d_pvalue_AGG_norm']], w=NULL)))) %>% 
    BH_procedure(Qcol=NULL, Pcol='HMP', Qname='HMQ', Q=0.1), 
  
  time3_avg_P05 = subset(to_filter, second_pct > 0.1 & ratio>1.5 & meta2d_pvalue_avg_3 < 0.05),
  time3_AGG_P05 = subset(to_filter, second_pct > 0.1 & ratio_AGG>1.5 & meta2d_pvalue_AGG_3 < 0.05),
  
  time3_HMP05 = subset(to_filter, second_pct>0.1 & ratio>1.5 & ratio_AGG>1.5, select=c(Gene, Sample, Mode, meta2d_pvalue_avg_3, meta2d_pvalue_AGG_3, meta2d_phase)) %>%
    cbind(., HMP = apply(., 1, function(L) hmp.stat(c(L[['meta2d_pvalue_avg_3']], L[['meta2d_pvalue_AGG_3']], w=NULL)))) %>% 
    subset(HMP<0.05), 
  time3_HMQ05 = subset(to_filter, second_pct>0.1 & ratio>1.5 & ratio_AGG>1.5, select=c(Gene, Sample, Mode, meta2d_pvalue_avg_3, meta2d_pvalue_AGG_3, meta2d_phase)) %>%
    cbind(., HMP = apply(., 1, function(L) hmp.stat(c(L[['meta2d_pvalue_avg_3']], L[['meta2d_pvalue_AGG_3']], w=NULL)))) %>% 
    BH_procedure(Qcol=NULL, Pcol='HMP', Qname='HMQ', Q=0.1), 
  
  time3_avg_Q05 = BH_procedure(subset(to_filter, second_pct > 0.1 & ratio>1.5), Qcol=NULL, Pcol='meta2d_pvalue_avg_3', Qname='avg3_BH.Q', Q=0.05), 
  time3_AGG_Q05 = BH_procedure(subset(to_filter, second_pct > 0.1 & ratio_AGG>1.5), Qcol=NULL, Pcol='meta2d_pvalue_AGG_3', Qname='AGG3_BH.Q', Q=0.05),
  
  pool3_avg = BH_procedure(subset(to_filter, second_pct > 0.1 & ratio>1.5), Qcol=NULL, Pcol='meta2d_pvalue_avg_pool3', Qname='avg_pool3_BH.Q', Q=0.05), 
  pool3_AGG = BH_procedure(subset(to_filter, second_pct > 0.1 & ratio_AGG>1.5), Qcol=NULL, Pcol='meta2d_pvalue_AGG_pool3', Qname='AGG_pool3_BH.Q', Q=0.05),
  
  pool_avg = BH_procedure(subset(to_filter, second_pct > 0.1 & ratio>1.5), Qcol=NULL, Pcol='meta2d_pvalue_avg_pool', Qname='avg_pool_BH.Q', Q=0.05), 
  pool_AGG = BH_procedure(subset(to_filter, second_pct > 0.1 & ratio_AGG>1.5), Qcol=NULL, Pcol='meta2d_pvalue_AGG_pool', Qname='AGG_pool_BH.Q', Q=0.05),
  pool_HMQ05 = subset(to_filter, second_pct>0.1 & ratio>1.5 & ratio_AGG>1.5, select=c(Gene, Sample, Mode, meta2d_pvalue_avg_pool, meta2d_pvalue_AGG_pool, meta2d_phase)) %>%
    cbind(., HMP = apply(., 1, function(L) hmp.stat(c(L[['meta2d_pvalue_avg_pool']], L[['meta2d_pvalue_AGG_pool']], w=NULL)))) %>% 
    BH_procedure(Qcol=NULL, Pcol='HMP', Qname='HMQ', Q=0.1),
  
  #dds_time3_avg_P05 = subset(to_filter, second_pct > 0.1 & ratio > 1.5 & meta2d_pvalue_avg_DESeq_time3 < 0.05),
  #dds_time3_avg_Q05 = BH_procedure(subset(to_filter, second_pct > 0.1 & ratio > 1.5), Qcol=NULL, Pcol='meta2d_pvalue_avg_DESeq_time3', Q=0.05), 
  
  reps_adding_pool = reps_adding_pool$res$AGG %>% do.call(what=rbind) %>%
    group_by(Gene, Sample, Mode) %>%
    summarise(HMP = hmp.stat(meta2d_pvalue), MinRatio = max(Ratio), meta2d_phase=Arg(mean(complex(argument = (2*pi*meta2d_phase)/24)))*24/(2*pi)) %>%
    mutate(meta2d_phase = meta2d_phase %% 24) %>%
    merge(., subset(to_filter, second_pct>0.1, select=c(Gene, Sample, Mode, second_pct, max, ampl))) %>% 
    subset(MinRatio > 1.5) %>%
    BH_procedure(Qcol=NULL, Pcol='HMP', Q=0.05) %>%
    mutate(Mean = (2*max-ampl)/2)
  
  #reps_adding_pool = Filtereds$DFs$reps_adding_pool
)



Filtereds$p_vals <- lapply(Filtereds$DFs, substatter, R_NOT_P=FALSE) %>%
  `names<-`(names(Filtereds$DFs))

Filtereds$r_means <- lapply(Filtereds$DFs, substatter, R_NOT_P=TRUE) %>%
  `names<-`(names(Filtereds$DFs))

Filtereds$p_val_DF <- lapply(rownames(Filtereds$p_vals[[1]]), function(RowName) {
  lapply(Filtereds$p_vals, function(DF) {
    DF = data.frame(value = DF[RowName, ])
    mutate(DF, Sample_Mode = rownames(DF), 
           Sample = gsub('_.*', '', rownames(DF)), 
           Mode = gsub('.*_', '', rownames(DF)), 
           stat = paste0('p_val_', RowName), 
    )}) %>%
    rbindlist(idcol='method')
}) %>% do.call(what=rbind) %>%
  as.data.frame

Filtereds$r_means_DF <- lapply(rownames(Filtereds$r_means[[1]]), function(RowName) {
  lapply(Filtereds$r_means, function(DF) {
    DF = data.frame(value = DF[RowName, ])
    mutate(DF, Sample_Mode = rownames(DF), 
           Sample = gsub('_.*', '', rownames(DF)), 
           Mode = gsub('.*_', '', rownames(DF)), 
           stat = paste0('r_means_', RowName), 
    )}) %>%
    rbindlist(idcol='method')
}) %>% do.call(what=rbind) %>%
  as.data.frame

Filtereds$tables <- lapply(Filtereds$DFs, function(DF) {
  Tb = table(DF$Sample, DF$Mode) %>%
    data.frame 
  if (nrow(Tb)) {
    colnames(Tb) = c('Sample', 'Mode', 'N')
  } else Tb = data.frame(Sample = character(0), Mode = character(0), N = numeric(0))
  return(Tb) 
}) %>%
  `names<-`(names(Filtereds$DFs))
  
Filtereds$table_RiboCell <- lapply(Filtereds$DFs, function(DF) {
  if (nrow(DF)) {
    DF = subset(DF, Gene %in% merge(GO_ribo_cell, passQC)$Gene)
    Tb = table(DF$Sample, DF$Mode) %>%
      data.frame  
    if (nrow(Tb)) {
      colnames(Tb) = c('Sample', 'Mode', 'N')
    } else Tb = data.frame(Sample = character(0), Mode = character(0), N = numeric(0))}
  else Tb = data.frame(Sample = character(0), Mode = character(0), N = numeric(0))
  return(Tb) 
}) %>%
  `names<-`(names(Filtereds$DFs))

expressed_r_means <- substatter(Filtereds$DFs$expressed, R_NOT_P=TRUE, Subsets = Subset_Stats$subsets['GO_ribo_noPre', 'GO_ribo_mito', 'GO_ribo_cell', 'ARGs'])

Filtereds_stats <- merge(mutate(rbindlist(Filtereds$tables, idcol='method'), 
                                count = N, .keep='unused'), 
                         mutate(rbindlist(Filtereds$table_RiboCell, idcol='method'), 
                                count_RiboCell = N, .keep='unused')) %>%
  mutate(pct_ribo = count_RiboCell/count) %>%
  pivot_longer(cols=-c(method, Sample, Mode), names_to = 'stat') %>%
  mutate(Sample_Mode = paste(Sample, Mode, sep='_')) %>%
  rbind(subset(Filtereds$p_val_DF, stat=='p_val_GO_ribo_cell') %>%
          mutate(value = -log10(value),
                 stat = '-log10(p_val_GO_ribo_cell)'), 
        subset(Filtereds$r_means_DF, stat=='r_means_GO_ribo_cell')) %>%
  mutate(QNOTP = grepl('Q', method))
  
saveRDS(Filtereds, file = 'Filtereds.rds')
#Filtereds <- readRDS('Filtereds.rds')
saveRDS(Filtereds_stats, file = 'Filtereds_stats.rds')



subset(Filtereds_stats, method %in% c("time3_avg_Q05", "time3_AGG_Q05", "dds_time3_avg_Q05", "time3_HMQ05", "pool3_avg", "pool3_AGG", "pool_avg", "pool_AGG", 'pool_HMQ05', 'reps_adding_pool') & 
         stat %in% c('count', 'pct_ribo', 'r_means_GO_ribo_cell')) %>%
  mutate(method = factor(method, levels = c("time3_avg_P05", "time3_AGG_P05", "dds_time3_avg_P05", "time3_HMP05", "time3_avg_Q05", "time3_AGG_Q05", "dds_time3_avg_Q05", "time3_HMQ05", "pool3_avg", "pool3_AGG", "pool_avg", "pool_AGG", 'pool_HMQ05', 'reps_adding_pool')), 
         stat = factor(stat, levels=c('count', 'r_means_GO_ribo_cell', 'pct_ribo'))) %>%
  ggplot(aes(x=method, y=value, fill=Sample_Mode)) +
  geom_bar(position="stack", stat="identity") +
  facet_grid(stat~., scales='free')



  
ggplot(Filtereds_stats, aes(x=QNOTP, y=value, fill=method)) +
  geom_bar(position="dodge", stat="identity") +
  facet_grid(stat~., scales='free')

subset(Filtereds_stats, Sample_Mode == 'C3 DD') %>%
  ggplot(aes(x=QNOTP, y=value, fill=method)) +
  geom_bar(position="dodge", stat="identity") +
  facet_grid(stat~., scales='free')


ggplot(Filtereds_stats, aes(x=method, y=value, fill=Sample_Mode)) +
  geom_bar(position="stack", stat="identity") +
  facet_grid(stat~., scales='free')


# SLIDE: Q-values!
subset(Filtereds_stats, method %in% c('avgP05', 'avgQ1') & 
         stat %in% c('count', 'pct_ribo', '-log10(p_val_GO_ribo_cell)', 'r_means_GO_ribo_cell')) %>%
  mutate(method = factor(method, levels = c('avgP05', 'avgQ1')), 
         stat = factor(stat, levels=c('count', 'pct_ribo', 'r_means_GO_ribo_cell', '-log10(p_val_GO_ribo_cell)'))) %>%
  ggplot(aes(x=method, y=value, fill=Sample_Mode)) +
  geom_bar(position="stack", stat="identity") +
  facet_wrap(~stat, scales='free')

subset(Filtereds_stats, method %in% c('avgP05', 'avgQ1') & 
         stat %in% c('count', 'pct_ribo', '-log10(p_val_GO_ribo_cell)', 'r_means_GO_ribo_cell') &
         Sample_Mode == 'C3_DD') %>%
  mutate(method = factor(method, levels = c('avgP05', 'avgQ1')), 
         stat = factor(stat, levels=c('count', 'pct_ribo', 'r_means_GO_ribo_cell', '-log10(p_val_GO_ribo_cell)'))) %>%
  ggplot(aes(x=method, y=value, fill=Sample_Mode)) +
  geom_bar(position="stack", stat="identity") +
  facet_wrap(~stat, scales='free')

# SLIDE: Normalization?
subset(Filtereds_stats, method %in% c('avgP05', 'AGGP05', 'avgQ1', 'AGGQ1') & 
         stat %in% c('count', 'pct_ribo', '-log10(p_val_GO_ribo_cell)', 'r_means_GO_ribo_cell')) %>%
  mutate(method = factor(method, levels = c('avgP05', 'AGGP05', 'avgQ1', 'AGGQ1')), 
         stat = factor(stat, levels=c('count', 'pct_ribo', 'r_means_GO_ribo_cell', '-log10(p_val_GO_ribo_cell)'))) %>%
  ggplot(aes(x=method, y=value, fill=Sample_Mode)) +
  geom_bar(position="stack", stat="identity") +
  facet_wrap(~stat, scales='free')

subset(Filtereds_stats, method %in% c('avgP05', 'AGGP05', 'avgQ1', 'AGGQ1') & 
         stat %in% c('count', 'pct_ribo', '-log10(p_val_GO_ribo_cell)', 'r_means_GO_ribo_cell') &
         Sample_Mode == 'C3_DD') %>%
  mutate(method = factor(method, levels = c('avgP05', 'AGGP05', 'avgQ1', 'AGGQ1')), 
         stat = factor(stat, levels=c('count', 'pct_ribo', 'r_means_GO_ribo_cell', '-log10(p_val_GO_ribo_cell)'))) %>%
  ggplot(aes(x=method, y=value, fill=Sample_Mode)) +
  geom_bar(position="stack", stat="identity") +
  facet_wrap(~stat, scales='free')

# SLIDE: HMP?
subset(Filtereds_stats, method %in% c('avgP05', 'AGGP05', 'avgQ1', 'AGGQ1', 'HMP', 'HMQ1', 'init') & 
         stat %in% c('count', 'pct_ribo', '-log10(p_val_GO_ribo_cell)', 'r_means_GO_ribo_cell')) %>%
  mutate(method = factor(method, levels = c('avgP05', 'AGGP05', 'avgQ1', 'AGGQ1', 'HMP', 'HMQ1', 'init')), 
         stat = factor(stat, levels=c('count', 'pct_ribo', 'r_means_GO_ribo_cell', '-log10(p_val_GO_ribo_cell)'))) %>%
  ggplot(aes(x=method, y=value, fill=Sample_Mode)) +
  geom_bar(position="stack", stat="identity") +
  facet_wrap(~stat, scales='free')

subset(Filtereds_stats, method %in% c('avgP05', 'AGGP05', 'avgQ1', 'AGGQ1', 'HMP', 'HMQ1', 'init') & 
         stat %in% c('count', 'pct_ribo', '-log10(p_val_GO_ribo_cell)', 'r_means_GO_ribo_cell') &
         Sample_Mode == 'C3_DD') %>%
  mutate(method = factor(method, levels = c('avgP05', 'AGGP05', 'avgQ1', 'AGGQ1', 'HMP', 'HMQ1', 'init')), 
         stat = factor(stat, levels=c('count', 'pct_ribo', 'r_means_GO_ribo_cell', '-log10(p_val_GO_ribo_cell)'))) %>%
  ggplot(aes(x=method, y=value, fill=Sample_Mode)) +
  geom_bar(position="stack", stat="identity") +
  facet_wrap(~stat, scales='free')

# Replicates Limit Performance
subset(Filtereds_stats, method %in% c('avgP05', 'HMP', 'init', 'time3_avg_P05', 'time3_AGG_P05', 'avgQ1', 'HMQ1', 'time3_avg_Q05', 'time3_AGG_Q05') & 
         stat %in% c('count', 'pct_ribo', '-log10(p_val_GO_ribo_cell)', 'r_means_GO_ribo_cell')) %>%
  mutate(method = factor(method, levels = c('avgP05', 'HMP', 'init', 'time3_avg_P05', 'time3_AGG_P05', 'avgQ1', 'HMQ1', 'time3_avg_Q05', 'time3_AGG_Q05')), 
         stat = factor(stat, levels=c('count', 'r_means_GO_ribo_cell', 'pct_ribo', '-log10(p_val_GO_ribo_cell)'))) %>%
  ggplot(aes(x=method, y=value, fill=Sample_Mode)) +
  geom_bar(position="stack", stat="identity") +
  facet_grid(stat~., scales='free')

# Which shuffle
subset(Filtereds_stats, method %in% c("time3_avg_P05", "time3_AGG_P05", "dds_time3_avg_P05", "time3_HMP05", "time3_avg_Q05", "time3_AGG_Q05", "dds_time3_avg_Q05", "time3_HMQ05", "pool3_avg", "pool3_AGG", "pool_avg", "pool_AGG", 'pool_HMQ05') & 
         stat %in% c('count', 'pct_ribo', '-log10(p_val_GO_ribo_cell)', 'r_means_GO_ribo_cell')) %>%
  mutate(method = factor(method, levels = c("time3_avg_P05", "time3_AGG_P05", "dds_time3_avg_P05", "time3_HMP05", "time3_avg_Q05", "time3_AGG_Q05", "dds_time3_avg_Q05", "time3_HMQ05", "pool3_avg", "pool3_AGG", "pool_avg", "pool_AGG", 'pool_HMQ05')), 
         stat = factor(stat, levels=c('count', 'r_means_GO_ribo_cell', 'pct_ribo', '-log10(p_val_GO_ribo_cell)'))) %>%
  ggplot(aes(x=method, y=value, fill=Sample_Mode)) +
  geom_bar(position="stack", stat="identity") +
  facet_grid(stat~., scales='free')

subset(Filtereds_stats, method %in% c("time3_avg_P05", "time3_AGG_P05", "dds_time3_avg_P05", "time3_HMP05", "time3_avg_Q05", "time3_AGG_Q05", "dds_time3_avg_Q05", "time3_HMQ05", "pool3_avg", "pool3_AGG", "pool_avg", "pool_AGG", 'pool_HMQ05') & 
         stat %in% c('count', 'pct_ribo', 'r_means_GO_ribo_cell')) %>%
  mutate(method = factor(method, levels = c("time3_avg_P05", "time3_AGG_P05", "dds_time3_avg_P05", "time3_HMP05", "time3_avg_Q05", "time3_AGG_Q05", "dds_time3_avg_Q05", "time3_HMQ05", "pool3_avg", "pool3_AGG", "pool_avg", "pool_AGG", 'pool_HMQ05')), 
         stat = factor(stat, levels=c('count', 'r_means_GO_ribo_cell', 'pct_ribo'))) %>%
  ggplot(aes(x=method, y=value, fill=Sample_Mode)) +
  geom_bar(position="stack", stat="identity") +
  facet_grid(stat~., scales='free')


lapply(levels(data_C2C3$sample_mode), function(S_M) { d
  S = gsub('_.*', '', S_M)
  M = gsub('.*_', '', S_M)
  G = subset(Filtereds$DFs$time3_avg_Q05, Sample==S & Mode==M)$Gene %>% unique %>% sort %>%
    {factor(., levels=.)}
  make_plot(Genes=G, S=data_C2C3, Type0 = S, Mode0=M, Group = 'time2') %>%
    ggsave(filename=paste0('avg_Q05_', S_M, '.pdf'), width=3, height=length(G)*0.8, limitsize=FALSE)
})

lapply(levels(data_C2C3$sample_mode), function(S_M) {
  S = gsub('_.*', '', S_M)
  M = gsub('.*_', '', S_M)
  G = subset(Filtereds$DFs$HMQ1, Sample==S & Mode==M)$Gene %>% unique %>% sort %>%
    {factor(., levels=.)}
  make_plot(Genes=G, S=data_C2C3, Type0 = S, Mode0=M, Group = 'time2') %>%
    ggsave(filename=paste0('HMQ1_', S_M, '.pdf'), width=3, height=length(G)*0.8, limitsize=FALSE)
})

lapply(levels(data_C2C3$sample_mode), function(S_M) {
  S = gsub('_.*', '', S_M)
  M = gsub('.*_', '', S_M)
  G = subset(Filtereds$DFs$time3_HMQ05, Sample==S & Mode==M)$Gene %>% unique %>% sort %>%
    {factor(., levels=.)}
  make_plot(Genes=G, S=data_C2C3, Type0 = S, Mode0=M, Group = 'time2') %>%
    ggsave(filename=paste0('HMQ05_', S_M, '.pdf'), width=3, height=length(G)*0.8, limitsize=FALSE)
})

lapply(levels(data_C2C3$sample_mode), function(S_M) {
  S = gsub('_.*', '', S_M)
  M = gsub('.*_', '', S_M)
  G = subset(Filtereds$DFs$pool_avg, Sample==S & Mode==M)$Gene %>% unique %>% sort %>%
    {factor(., levels=.)}
  make_plot(Genes=G, S=data_C2C3, Type0 = S, Mode0=M, Group = 'time2') %>%
    ggsave(filename=paste0('PoolQ05_', S_M, '.pdf'), width=3, height=length(G)*0.8, limitsize=FALSE)
})

lapply(levels(data_C2C3$sample_mode), function(S_M) {
  S = gsub('_.*', '', S_M)
  M = gsub('.*_', '', S_M)
  G = subset(Filtereds$DFs$pool_HMQ05, Sample==S & Mode==M)$Gene %>% unique %>% sort %>%
    {factor(., levels=.)}
  make_plot(Genes=G, S=data_C2C3, Type0 = S, Mode0=M, Group = 'time2') %>%
    ggsave(filename=paste0('PoolHMQ05_', S_M, '.pdf'), width=3, height=length(G)*0.8, limitsize=FALSE)
})



lapply(levels(data_C2C3$sample_mode), function(S_M) { d
  S = gsub('_.*', '', S_M)
  M = gsub('.*_', '', S_M)
  G = subset(Filtereds$DFs$time3_avg_Q05, Sample==S & Mode==M)$Gene %>% unique %>% sort %>%
    {factor(., levels=.)}
  make_plot(Genes=G, S=added_C2C3, Type0 = S, Mode0=M, Group = 'time3') %>%
    ggsave(filename=paste0('avg_time3_Q05_3_', S_M, '.pdf'), width=3, height=length(G)*0.8, limitsize=FALSE)
})

lapply(levels(data_C2C3$sample_mode), function(S_M) { d
  S = gsub('_.*', '', S_M)
  M = gsub('.*_', '', S_M)
  G = subset(Filtereds$DFs$time3_AGG_Q05, Sample==S & Mode==M)$Gene %>% unique %>% sort %>%
    {factor(., levels=.)}
  make_plot(Genes=G, S=added_C2C3_AGG, Type0 = S, Mode0=M, Group = 'time3') %>%
    ggsave(filename=paste0('AGG_time3_Q05_3_', S_M, '.pdf'), width=3, height=length(G)*0.8, limitsize=FALSE)
})

lapply(levels(data_C2C3$sample_mode), function(S_M) {
  S = gsub('_.*', '', S_M)
  M = gsub('.*_', '', S_M)
  G = subset(Filtereds$DFs$dds_time3_avg_Q05, Sample==S & Mode==M)$Gene %>% unique %>% sort %>%
    {factor(., levels=.)}
  make_plot(Genes=G, S=DDS, Type0 = S, Mode0=M, Group = 'time3', COUNTS=TRUE) %>%
    ggsave(filename=paste0('DDS_time3_Q05_3_', S_M, '.pdf'), width=3, height=length(G)*0.8, limitsize=FALSE)
})

lapply(levels(data_C2C3$sample_mode), function(S_M) {
  S = gsub('_.*', '', S_M)
  M = gsub('.*_', '', S_M)
  G = subset(Filtereds$DFs$reps_adding_pool, Sample==S & Mode==M)$Gene %>% unique %>% sort %>%
    {factor(., levels=.)}
  make_plot(Genes=G, S=reps_adding_pool$R0$AGG, Type0 = S, Mode0=M, Group = 'time3', COUNTS=FALSE) %>%
    ggsave(filename=paste0('reps_adding_pool_', S_M, '.pdf'), width=3, height=length(G)*0.8, limitsize=FALSE)
})





Bulk_Ines <- read.table('/Users/Teddy/Desktop/to_GEO/Bulk/GSE233184_PolyAnormtomax_6tp.txt', header=TRUE) %>%
  pivot_longer(cols=-gene, names_to='cond') %>%
  mutate(time1 = gsub('w1118_ZT', '', cond) %>% 
           gsub('_.*', '', .) %>% as.numeric, 
         rep = gsub('_S.*', '', cond) %>%
           gsub('.*_', '', .) %>%
           as.numeric) %>%
  mutate(time2 = time1+24*(rep-1))

subset(Bulk_Ines, gene %in% merge(passQC, GO_ribo_mito)$Gene) %>%
  ggplot(aes(time2, value)) +
  geom_point() +
  geom_line() +
  facet_wrap(~gene, scales='free')

subset(Bulk_Ines, gene %in% merge(passQC, GO_ribo_cell)$Gene) %>%
  ggplot(aes(time2, value)) +
  geom_point() +
  geom_line() +
  facet_wrap(~gene, scales='free')



heads <- read.xlsx('heads_cyclers.xlsx', sheet=2) %>% #sheet2 is 25C
  as.data.frame %>%
  subset(AMP_min > 1.5 & JTK_pvalue < 0.05)

View(heads)
nrow(heads)
subset(heads, gene.name %in% merge(passQC, GO_ribo_cell))
heads[grep('mRp', heads$gene.name), grepl("[[:digit:]]", colnames(heads)) & grepl("zt", colnames(heads)) | colnames(heads)=='gene.name'] %>%
  pivot_longer(cols=-gene.name, names_to='time') %>% 
  mutate(time = factor(time, levels=unique(time))) %>%
  ggplot(aes(time, value, group=gene.name)) +
  geom_point() +
  geom_line() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), size=5) +
  facet_grid(gene.name~., scales='free')



Bulk_Litovchenko <- read.xlsx('/Users/Teddy/Desktop/Brain cycle/Bulk_Litovchenko.xlsx', sheet='S2', startRow=29)

Bulk_Litovchenko %>%
  mutate('Rp' = grepl('^RpS|^RpL|^mRpS|^mRpL', Gene)) %>%
  group_by(Tissue, Rp) %>%
  summarise(N = n()) %>%
  pivot_wider(names_from=Rp, values_from=N)

rownames(data_C2C3)
Bulk_Litovchenko %>%
  mutate('Rp' = grepl('^RpS|^RpL|^mRpS|^mRpL', Gene)) %>%
  {split(., .$Tissue)} %>% 
  lapply(function(DF) run_fisher_general(
    A = DF$Gene %>% {.[. %in% rownames(data_C2C3)]}, 
    B = rownames(data_C2C3) %>% {.[grepl('^RpS|^RpL|^mRpS|^mRpL', .)]}, 
    U = rownames(data_C2C3)
  ))

ggarrange(nrow=2, ncol=3, common.legend = TRUE, legend = 'right',
          Bulk_Litovchenko %>%
            mutate('Rp' = grepl('^RpS|^RpL|^mRpS|^mRpL', Gene)) %>%
            ggplot(aes(LAG, color=Rp)) +
            geom_density() + 
            ggtitle('in all') +
            xlab('Phase'), 
          
          Bulk_Litovchenko %>%
            subset(Tissue=='fat body') %>%
            mutate('Rp' = grepl('^RpS|^RpL|^mRpS|^mRpL', Gene)) %>%
            ggplot(aes(LAG, color=Rp)) +
            geom_density() + 
            ggtitle('in fat body') +
            xlab('Phase'),
          
          Bulk_Litovchenko %>%
            subset(Tissue=='brain') %>%
            mutate('Rp' = grepl('^RpS|^RpL|^mRpS|^mRpL', Gene)) %>%
            ggplot(aes(LAG, color=Rp)) +
            geom_density() + 
            ggtitle('in brain') +
            xlab('Phase'),
          
          Bulk_Litovchenko %>%
            subset(Tissue=='gut') %>%
            mutate('Rp' = grepl('^RpS|^RpL|^mRpS|^mRpL', Gene)) %>%
            ggplot(aes(LAG, color=Rp)) +
            geom_density() + 
            ggtitle('in gut') +
            xlab('Phase'), 
          
          Bulk_Litovchenko %>%
            subset(Tissue=='malpighian tubules') %>%
            mutate('Rp' = grepl('^RpS|^RpL|^mRpS|^mRpL', Gene)) %>%
            ggplot(aes(LAG, color=Rp)) +
            geom_density() + 
            ggtitle('in malpighian tubules') +
            xlab('Phase')
)

# BELOW DEPENDS ON #41
subset(resCell, p_cell < 0.05 & ratio > 1.5) %>% {table(.$Sample, .$Mode)}

subset(resCell, p_cell < 0.05 & ratio > 1.5 & Sample == 'C3' & Mode == 'LD')$Gene %>%
  split_list(10) %>%
  lapply(function(G) make_plot(G, Group='time2', Type0 = 'C3', Mode0 = 'LD', S=shuffled_data, COUNTS = FALSE)) %>%
  plot_grid(plotlist = ., nrow=1) %>%
  ggsave(filename = paste0('shuffle_cyclers.pdf'), width=30, height=200/10, limitsize=FALSE)


data_TRM_cluster <- read.csv('/Users/Teddy/Dropbox/repo2/data_TRM_cluster.csv', row.names=1)
data_TRM_primary <- read.csv('/Users/Teddy/Dropbox/repo2/data_TRM_primary.csv', row.names=1)

data_TRM_cluster %>% 
  group_by(Gene, Cluster) %>%
  summarise(Max = max(value), 
            Min = min(value), 
            Mean = mean(value)) %>%
  mutate(ratio = Max/Min, 
         amplitude = Max-Min, 
         nAmp = (Max-Min)*2/(Max+Min)) %>%
  subset(grepl('^RpS|^RpL|^mRpS|^mRpL', Gene)) %>% 
  mutate(Cluster = factor(Cluster, levels=unique(Cluster))) %>% 
  ggplot(aes(amplitude, fill=Cluster, group=Cluster)) +
  geom_density(alpha=0.2)



data_TRM_cluster %>% 
  group_by(Gene, Cluster) %>%
  summarise(Max = max(value), 
            Min = min(value), 
            Mean = mean(value)) %>%
  mutate(ratio = Max/Min, 
         amplitude = Max-Min, 
         nAmp = (Max-Min)*2/(Max+Min)) %>%
  subset(grepl('^RpS|^RpL|^mRpS|^mRpL', Gene)) %>% 
  mutate(Cluster = factor(Cluster, levels=unique(Cluster))) %>% 
  ggplot(aes(nAmp)) +
  geom_density()

data_TRM_primary %>% 
  group_by(Gene, TypePrimary) %>%
  summarise(Max = max(value), 
            Min = min(value), 
            Mean = mean(value)) %>%
  mutate(ratio = Max/Min, 
         amplitude = Max-Min, 
         nAmp = (Max-Min)*2/(Max+Min)) %>%
  subset(grepl('^RpS|^RpL|^mRpS|^mRpL', Gene)) %>% 
  mutate(TypePrimary = factor(TypePrimary, levels=unique(TypePrimary))) %>% 
  ggplot(aes(nAmp, fill=TypePrimary, group=TypePrimary)) +
  geom_density()

rbind(
  cbind(Data='C2C3', 
        expr_data %>% 
          subset(grepl('^RpS|^RpL|^mRpS|^mRpL', Gene), select = c(Gene, rel_ampl, ratio))) %>% 
    `colnames<-`(c('Data', 'Gene', 'nAmp', 'ratio')),
  
  cbind(Data='glia', 
        data_TRM_primary %>% 
          group_by(Gene, TypePrimary) %>%
          summarise(Max = max(value), 
                    Min = min(value), 
                    Mean = mean(value)) %>%
          mutate(ratio = Max/Min, 
                 amplitude = Max-Min, 
                 nAmp = (Max-Min)*2/(Max+Min)) %>%
          subset(grepl('^RpS|^RpL|^mRpS|^mRpL', Gene), select = c(Gene, nAmp, ratio)))) %>%
  ggplot(aes(nAmp, fill=Data, group=Data)) +
  geom_density()
  

# "~/Desktop/Rhythm_1_backup/ext_data/Dopp.xlsx"
Dopps <- read.xlsx("~/Desktop/Rhythm_1_backup/ext_data/Dopp.xlsx", sheet = 'S3', startRow = 3) %>%
  mutate(Type = gsub("[^a-zA-Z]", "", replicate)) %>%
  mutate(HAS_CLOCK = ifelse(grepl('CXG|ALG|PG|EG|clock', .$cluster), 'Clock cell', 'Non-clock')) %>%
  mutate(Rp = grepl('^Rp', gene))






Dopp_sum <- Dopps %>%
  group_by(cluster, HAS_CLOCK, Rp, .drop=FALSE) %>%
  summarise(N = n()) %>%
  pivot_wider(names_from='Rp', values_from=N) %>%
  `colnames<-`(c('cluster', 'HAS_CLOCK', 'Non-Rp', 'Rp')) %>%
  mutate(Rp = replace_na(Rp, 0)) %>%
  mutate(pct_Rp = Rp/(Rp + `Non-Rp`)) %>%
  subset(`Non-Rp`+Rp > 15) %>%
  arrange(desc(Rp))

subset(Dopps, Rp==TRUE)$JTK_adjphase %>%
  table


plot_grid(nrow=1,
          ggboxplot(Dopp_sum, x='HAS_CLOCK', y='pct_Rp', 
                    color='HAS_CLOCK', 
                    palette =c("#00AFBB", "#FC4E07"), 
                    add='jitter') +
            stat_compare_means() +
            ylab('Fraction of Cyclers in Rp') +
            theme_prism() +
            theme(legend.position='none', axis.title.x=element_blank()),
          
          ggboxplot(Dopp_sum, x='HAS_CLOCK', y='Rp', 
                    color='HAS_CLOCK', 
                    palette =c("#00AFBB", "#FC4E07"), 
                    add='jitter') +
            stat_compare_means() +
            ylab('Number of Cyclers in Rp') +
            theme_prism() +
            theme(legend.position='none', axis.title.x=element_blank())) %>%
  ggsave(filename='barplot.pdf', width=8, height=3)




# ALSO, for Seb:


WD_REPO2 <- '/Users/Teddy/Dropbox/repo2/files_cycling/JTK_TRM_type_AGG.csv'
WD_Litovchenko <- '/Users/Teddy/Desktop/Brain cycle/Bulk_Litovchenko.xlsx' # Download supplementary table from Litovchenko et al to this directory
WD_ANE <- "/Users/Teddy/Downloads/elife-44642-supp1.xlsx" # Download supplementary table from thermosensitive splicing, 25C used
WD_DOPP <- "~/Desktop/Rhythm_1_backup/ext_data/Dopp.xlsx" # Download supplementary table Dopp et al to this directory

DOPP_passes <- read.xlsx(WD_DOPP, sheet = 'S3', startRow = 3) %>%
  mutate(Type = gsub("[^a-zA-Z]", "", replicate)) %>%
  mutate(HAS_CLOCK = ifelse(grepl('CXG|ALG|PG|EG|clock', .$cluster), 'Clock cell', 'Non-clock'))

ANE_bulk <- read.xlsx(WD_ANE, sheet='25C', startRow = 1)

Litovchenko_bulk <- read.xlsx(WD_Litovchenko, sheet='S2', startRow=29)



ANE_passes <- subset(ANE_bulk, temp == '25C' & AMP_min > 1.5 & JTK_pvalue < 0.05)
Litovchenko_passes <- subset(Litovchenko_bulk, ADJ.P < 0.05 & Tissue == 'brain')

nonDopp_bulk <- unique(c(ANE_passes$gene.name, Litovchenko_passes$Gene) %>% {.[! . %in% DOPP_passes$gene]})
write.table(nonDopp_bulk, file='/Users/Teddy/Downloads/nonDopp_bulk.txt', col.names=FALSE, row.names=FALSE)


REPO2_passes <- read.csv(WD_REPO2) %>%
  subset(AMP < 1.5) %>%
  cbind(BH.Q = p.adjust(.[['ADJ.P']], method='BH')) %>%
  subset(BH.Q < 0.05)

nonRepo2_bulk <- unique(c(ANE_passes$gene.name, Litovchenko_passes$Gene) %>% {.[! . %in% REPO2_passes$Gene]})
write.table(nonRepo2_bulk, file='/Users/Teddy/Downloads/nonRepo2_bulk.txt', col.names=FALSE, row.names=FALSE)

####    ####    ####    ####    ####    ####    ####    ####    ####    ####   


# 41. Running one shuffle
################################################################################
{set.seed(1)
  #shuffle_2cells <- list()
  prefilt <- subset(to_filter, second_pct > 0.10, select=c(Gene, Sample, Mode, second_pct, ratio, ratio_AGG, meta2d_pvalue, meta2d_pvalue_AGG_norm)) %>%
    mutate(Gene_SM = paste(Gene, Sample, Mode, sep='_'))
  
  split_object <- SplitObject(data_C2C3, split.by = "sample_mode")
  data_lengths <- lapply(split_object, ncol)
  
  
  indices <- lapply(data_lengths, function(x) sample(1:x, replace=FALSE))
  
  shuffled_data <- list()
  for (j in 1:length(split_object)) {
    Data <- split_object[[j]]
    indices <- sample(1:ncol(Data), replace=FALSE)
    Data@meta.data[['sample_mode_group2']] <- Data@meta.data[['sample_mode_group2']][indices]
    Data@meta.data[['group2']] <- Data[['sample_mode_group2']] %>% 
      mutate(sample = gsub('.*-', '', sample_mode_group2), .keep='unused') %>% unlist
    shuffled_data[[j]] <- Data
  }
  shuffled_data = merge(shuffled_data[[1]], shuffled_data[2:4])
  
  shuffled_data@meta.data[['time2']] <- shuffled_data[[c('group2')]] %>% 
    mutate(group2 = gsub('.*_', '', as.character(group2))) %>% 
    mutate(time2 = as.numeric(gsub("[a-zA-Z]", "", group2))) %>% 
    mutate(time2 = factor(time2, sort(unique(time2)))) %>% {.[, 'time2']}
  
  preMeta_cell <- AverageExpression(shuffled_data, group.by = c("group2", "sample", "mode1"), layer="data")$RNA %>% as.data.frame %>%
    {mutate(., Gene = rownames(.))} %>%
    pivot_longer(cols=-'Gene', names_to='TP_sample_mode') %>%
    mutate(Timepoint = TP_sample_mode %>% gsub('_.*', '', .) %>% gsub('.*T', '', .), 
           Gene_SM = paste(Gene, TP_sample_mode %>% gsub('.*_C', 'C', .), sep = '_'), .keep='unused') %>%
    subset(Gene_SM %in% prefilt$Gene_SM) %>%
    pivot_wider(names_from = 'Timepoint', values_from = 'value')
  
  shuffle_indices = lapply(1:nrow(prefilt), function(c) {
    sample(2:13, replace=FALSE)
  })
  
  lapply(seq_along(shuffle_indices), function(x) {
    preMeta_cell[x, c(1, shuffle_indices[[x]])] %>%
      `colnames<-`(c('Gene_SM', seq(3,47,4)))}) %>% do.call(what='rbind') %>%
    write.csv(file = "shuffle_2cells_cell.csv", row.names=FALSE)
  
  meta2d(infile = "shuffle_2cells_cell.csv", filestyle = "csv", timepoints = "line1", minper=20, 
         maxper=28, outdir=getwd(), parallelize = TRUE)
  
  resCell = merge(read.csv('meta2d_shuffle_2cells_cell.csv'), 
                  data.frame(CycID = preMeta_cell$Gene_SM, ratio = preMeta_cell[, 2:13] %>% apply(1, function(x) max(x)/min(x)))) %>% 
    mutate(Gene = gsub('_.*', "", CycID), 
           Mode = gsub('.*_', "", CycID), 
           Sample = gsub('.*_C', 'C', CycID) %>% sub('_.*', '', .),
           p_cell = meta2d_pvalue,
           .keep='unused') %>%
    subset(select = c(Gene, Mode, Sample, p_cell, ratio))
}



####    ####    ####    ####    ####    ####    ####    ####    ####    ####   


# 42. Data2
################################################################################


####    ####    ####    ####    ####    ####    ####    ####    ####    ####   


# 2.1. Analysis of Filtereds$DFs$reps_adding_pool
################################################################################


Filtereds$DFs$reps_adding_pool

# Print genes
lapply(levels(data_C2C3$sample_mode), function(S_M) {
  S = gsub('_.*', '', S_M)
  M = gsub('.*_', '', S_M)
  G = subset(Filtereds$DFs$reps_adding_pool, Sample==S & Mode==M)$Gene %>% unique %>% sort %>%
    {factor(., levels=.)}
  make_plot(Genes=G, S=reps_adding_pool$R0$AGG, Type0 = S, Mode0=M, Group = 'time3', COUNTS=FALSE) %>%
    ggsave(filename=paste0('reps_adding_pool/reps_adding_pool_', S_M, '.pdf'), width=3, height=length(G)*0.8, limitsize=FALSE)
})


lapply(levels(data_C2C3$sample_mode), function(S_M) {
  S = gsub('_.*', '', S_M)
  M = gsub('.*_', '', S_M)
  repplotter_ofSummary(DF=avg_data_ALL_adding_pool, Gene0='Eip74EF', Sample0 = S, Mode0=M)
}) %>% {plot_grid(plotlist=.)}


####    ####    ####    ####    ####    ####    ####    ####    ####    ####   


# 2.2. Print genes for GO
################################################################################
merge(Filtereds$DFs$reps_adding_pool, passQC)

# Genes - EnsID
lapply(levels(data_C2C3$sample_mode), function(S_M) {
  S = gsub('_.*', '', S_M)
  M = gsub('.*_', '', S_M)
  subset(merge(Filtereds$DFs$reps_adding_pool, passQC), Sample==S & Mode==M)$ENSEMBL %>%
    write.table(file=paste0('reps_adding_pool/EnsID_', S_M, '.txt'), quote = FALSE, row.names = FALSE, col.names = FALSE)
})
# Genes 
lapply(levels(data_C2C3$sample_mode), function(S_M) {
  S = gsub('_.*', '', S_M)
  M = gsub('.*_', '', S_M)
  subset(Filtereds$DFs$reps_adding_pool, Sample==S & Mode==M)$Gene %>%
    write.table(file=paste0('reps_adding_pool/Gene_', S_M, '.txt'), quote = FALSE, row.names = FALSE, col.names = FALSE)
})
# Backgrounds - EnsID
lapply(levels(data_C2C3$sample_mode), function(S_M) {
  S = gsub('_.*', '', S_M)
  M = gsub('.*_', '', S_M)
  subset(merge(to_filter, passQC), Sample==S & Mode==M & second_pct > 0.1)$ENSEMBL %>%
    write.table(file=paste0('reps_adding_pool/EnsID_background_', S_M, '.txt'), quote = FALSE, row.names = FALSE, col.names = FALSE)
})
# Backgrounds
lapply(levels(data_C2C3$sample_mode), function(S_M) {
  S = gsub('_.*', '', S_M)
  M = gsub('.*_', '', S_M)
  subset(to_filter, Sample==S & Mode==M & second_pct > 0.1)$Gene %>%
    write.table(file=paste0('reps_adding_pool/Gene_background_', S_M, '.txt'), quote = FALSE, row.names = FALSE, col.names = FALSE)
})




####    ####    ####    ####    ####    ####    ####    ####    ####    ####   


# 2.3. Print gene types
################################################################################
# print TFs
print_overlaps = function(DF_overlapping, FileStart, PLOT=TRUE, TABLE=TRUE, 
                          Data=reps_adding_pool$R0$AGG, DataCol='time3') {
  # DF_overlapping contains at least Gene, Sample, Mode already filtered to overlap
  res = list()
  
  if (nrow(DF_overlapping) != nrow(filter(group_by(DF_overlapping, Gene), n() > 1))) {
    stop('DF contains unique genes')
  }
  
  if (TABLE) {
    DF_overlapping %>%
      subset(select = c(Gene, Sample, Mode)) %>%
      write.table(file=paste0(FileStart, '_overlaps.txt'), quote = FALSE, row.names = FALSE)
  }
  
  if (PLOT) {
    DF_overlapping %>%
      {split(., .$Gene)} %>% 
      {ggsave(filename=paste0(FileStart, '_overlaps.pdf'), width=4, height=4*4*length(.)/8, limitsize = FALSE,
        plot = plot_grid(rel_heights = 1+1/unlist(lapply(., nrow)), ncol=1,
          plotlist = lapply(., 
            function(DF) {
              plot_grid(nrow=1,
                plotlist = apply(DF, 1, FUN = function(DF2) {
                  S = DF2[['Sample']]
                  M = DF2[['Mode']]
                  G = DF2[['Gene']]
                  make_plot(Genes=G, Type0=S, Mode0=M, S=Data, Group=DataCol,
                            TITLE = paste(S, M, G), FACET_STRIPS = FALSE, X_TEXT = FALSE) +
                    theme(plot.title = element_text(size=10), 
                          axis.title.x=element_blank())
                })
              )
            }
          ) 
        ) %>% print
      )
      }
  }
}

GO_group_maker = function(Genes_DF, Path, NameStart) {
  dir.create(file.path(getwd(), Path), showWarnings = FALSE, recursive = TRUE)
  NameStart = paste(Path, NameStart, sep='/')
  
  lapply(levels(data_C2C3$sample_mode), function(S_M) {
    S = gsub('_.*', '', S_M)
    M = gsub('.*_', '', S_M)
    subset(merge(Genes_DF, merge(Filtereds$DFs$reps_adding_pool, passQC)), Sample==S & Mode==M, select = c(ENSEMBL, Gene, meta2d_phase, HMP, MinRatio)) %>%
      mutate_if(is.numeric, signif, digits=3) %>%
      write.table(file=paste0(NameStart, '_', S_M, '.txt'), quote = FALSE, row.names = FALSE)
  })
  
  subset(merge(Genes_DF, merge(Filtereds$DFs$reps_adding_pool, passQC)), select = c(Gene, Sample, Mode)) %>%
    group_by(Gene) %>%
    filter(n() > 1) %>%
    print_overlaps(FileStart = NameStart)
  
  lapply(levels(data_C2C3$sample_mode), function(S_M) {
    S = gsub('_.*', '', S_M)
    M = gsub('.*_', '', S_M)
    G = subset(merge(Genes_DF, merge(Filtereds$DFs$reps_adding_pool, passQC)), Sample==S & Mode==M)$Gene
    make_plot(Genes=G, S=reps_adding_pool$R0$AGG, Type0 = S, Group = 'time6', COUNTS=FALSE) %>%
      ggsave(filename=paste0(NameStart, '_', S_M, '.pdf'), width=4, height=length(G)*0.8, limitsize=FALSE)
  })
}

GO_group_maker(Genes_DF=GO_TRs, Path = 'reps_adding_pool/3_groupings/TFs', NameStart='TFs')

GO_group_maker(Genes_DF=GO_RBPs, Path = 'reps_adding_pool/3_groupings/RBPs', NameStart='RBPs')

GO_group_maker(Genes_DF=GO_circ, Path = 'reps_adding_pool/3_groupings/Circ', NameStart='Circ')

data.frame(Gene = ARG$converted$min2) %>%
 GO_group_maker(Path = 'reps_adding_pool/3_groupings/ARGs', NameStart='ARGs')



# reps_adding_pool/rmeans.pdf
Filtereds$r_means$reps_adding_pool %>%
  data.frame %>%
  gt(rownames_to_stub = T) %>%
  fmt_number(columns = everything(), n_sigfig = 3, suffixing = TRUE) %>%
  data_color(
    columns = everything(),
    colors = col_bin(palette = c("blue", "red"), domain=NULL, bins=c(-1, seq(0.1, 1, 0.1)))) %>%
  gtsave("reps_adding_pool/rmeans.png")

# reps_adding_pool/pvals.pdf
Filtereds$p_vals$reps_adding_pool %>%
  data.frame %>%
  gt(rownames_to_stub = T) %>%
  fmt_number(columns = everything(), n_sigfig = 3, suffixing = TRUE) %>%
  data_color(
    columns = everything(),
    colors = col_bin(palette = c("red", "blue"), domain=NULL, bins=10^(0:floor(log10(min(Filtereds$p_vals$reps_adding_pool)))))) %>%
  gtsave("reps_adding_pool/pvals.png")

library(readxl)
library(writexl)
# prune sheet 'Genes' entries
lapply(excel_sheets('reps_adding_pool/2_GO/GOrilla.xlsx'), function(sheet) {
  DF = read_excel('reps_adding_pool/2_GO/GOrilla.xlsx', sheet = sheet) 
  DF = cbind(subset(DF, select=-Genes),
             Genes = matrix(unlist(
               DF$Genes %>%
                 lapply(function(G) {
                   gsub('\\[|\\]', '', G) %>%
                     strsplit(', ') %>% 
                     lapply(function(x) gsub('  - .*', '', x)) %>% {.[[1]]} %>%
                     paste(collapse=',')
                 })))
  )
}) %>%
  `names<-`(excel_sheets('reps_adding_pool/2_GO/GOrilla.xlsx')) %>%
  write_xlsx('reps_adding_pool/2_GO/GOrilla_trimmed.xlsx')


####    ####    ####    ####    ####    ####    ####    ####    ####    ####   


# 2.6. ARG timescales - USED
################################################################################
ARGs_DF <- list.files(file.path(WD, 'ARGs/Timeseries')) %>%
  lapply(function(Dir) {
    print(Dir)
    Files = list.files(file.path(WD, 'ARGs/Timeseries', Dir)) 
    metadata = Files %>% 
      {.[grep('.csv', .)]} %>% 
      {read.csv(file.path(WD, 'ARGs/Timeseries', Dir, .), header = FALSE)}
    
    Data = lapply(metadata[[1]], function(GSM) {
      Files[grep(GSM, Files)] %>%
        {file.path(WD, 'ARGs/Timeseries', Dir, .)} %>%
        read.table(header=TRUE)
    }) %>%
      `names<-`(metadata[[1]]) %>%
      rbindlist(idcol='GSM')
    
    list(DF = Data, 
         metadata = metadata) %>%
      return()
  }) %>% 
  {list(DF = do.call(what=rbind, lapply(., function(L2) L2$DF )), 
        metadata = do.call(what=rbind, lapply(., function(L2) L2$metadata)) %>%
          `names<-`(c('GSM', 'Stimulus', 'Time', 'Rep')))} %>%
  {merge(.$DF, .$metadata)} %>%
  subset(select = c(c(GSM, gene_id, FPKM, Stimulus, Time, Rep))) %>%
  as.data.frame %>%
  merge(ARG$converting$converter %>% `colnames<-`(c('gene_id', 'Gene'))) %>%
  mutate()



ARGs_quickness_rankings <- ARGs_DF %>%
  merge(passQC) %>% 
  group_by(Gene, Stimulus, Rep) %>%
  left_join(mutate(filter(., Time == 0), baseline_FPKM = FPKM, Time=NULL,GSM=NULL, .keep='unused')) %>%
  mutate(FC = FPKM/baseline_FPKM, 
         delta_FPKM = FPKM-baseline_FPKM) %>%
  group_by(Gene, Stimulus, Rep) %>%
  mutate(normFC = FC/max(FC), 
         norm_delta_FPKM = delta_FPKM/max(delta_FPKM), 
         TimeScore = (max(Time)-Time)/max(Time)) %>%
  mutate(Score = max(0, norm_delta_FPKM*TimeScore) %>% replace_na(0)) %>% 
  group_by(Gene, Stimulus) %>%
  summarise(Score = sum(Score)) %>%
  group_by(Gene) %>% 
  mutate(meanScore = mean(Score)) %>%
  group_by(Stimulus) %>%
  mutate(Rank = rank(Score, ties.method = 'first'), 
         `mean score rank` = rank(meanScore, ties.method = 'first'), 
         .keep='unused') %>% 
  pivot_wider(names_from = Stimulus, values_from = Rank) %>%
  as.data.frame %>%
  arrange(`mean score rank`) %>%
  merge(reps_adding_pool$res$AGG %>% do.call(what=rbind) %>%
          group_by(Gene, Sample, Mode) %>%
          summarise(HMP = hmp.stat(meta2d_pvalue), MinRatio = max(Ratio), meta2d_phase=Arg(mean(complex(argument = (2*pi*meta2d_phase)/24)))*24/(2*pi)) %>%
          merge(., subset(to_filter, second_pct>0.1, select=c(Gene, Sample, Mode))))



# PLOT p-values as a function of quickness

ARGs_2FC <- ARGs_DF %>%
  merge(passQC) %>% 
  group_by(Gene, Stimulus, Rep) %>%
  left_join(mutate(filter(., Time == 0), baseline_FPKM = FPKM, Time=NULL,GSM=NULL, .keep='unused')) %>%
  group_by(Gene, Stimulus, Time) %>%
  mutate(baseline_FPKM = mean(baseline_FPKM), 
         FPKM=mean(FPKM)) %>%
  mutate(FC = FPKM/baseline_FPKM, 
         delta_FPKM = FPKM-baseline_FPKM) %>% 
  subset(FC > 2) %>%
  group_by(Gene, Stimulus) %>%
  summarise(Time = min(Time)) %>%
  arrange(Time) %>% 
  ungroup

ARGs_2FC %>% 
  {split(., .$Stimulus)}

ARGs_2FC %>%
  arrange(Time) %>%
  group_by(Gene) %>%
  filter(n() >= 3) %>%
  group_by(Gene) %>%
  summarise(Time = min(Time)) # CG14186, CG46385, Hr38, Hsp23







ARGs_strongest <- ARGs_2FC %>%
  arrange(Time) %>%
  group_by(Gene) %>%
  filter(n() >= 2) %>%
  group_by(Gene) %>%
  summarise(Time = min(Time)) # CG13055, CG14186, CG46385, Hr38, sr

ARGs_strongest_cyc <- subset(ARGs_strongest, Gene %in% Filtereds$DFs$reps_adding_pool$Gene)

merge(ARGs_strongest_cyc, Filtereds$DFs$reps_adding_pool) %>%
  {split(., list(.$Sample, .$Mode))}


plot_grid(ncol=1, # C3
          plotter_geneCombine( subset(to_filter, second_pct > 0.1 & Gene %in% ARGs_strongest_cyc$Gene & Sample=='C3')$Gene,
                               data=avg_reps_adding_pool_AGG_01_1d, Sample='C3', dgrp_colname='time6', Mode=c('LD', 'DD'), N_days=6, COUNTS=TRUE,
                               NO_X_TEXT=TRUE, NO_Y_TEXT=FALSE, NO_GRIDLINES = TRUE, PRINT_GENES = FALSE, LEGEND=FALSE, offset = -0.5, ALIGN=FALSE, 
                               gene_linesize = 0.5, gene_alpha = 0.7, mean_linesize=0.5, 
                               include_SE_bars = TRUE, rib_alpha=0.8, errorWidth=0, errorThick=0.5,
                               PRISM=TRUE, GENECOLOR=FALSE) +
            theme(strip.text.y = element_text(margin = margin(0,1/72,0,1/72, "in"))),
          # C2
          plotter_geneCombine( subset(to_filter, second_pct > 0.1 & Gene %in% ARGs_strongest_cyc$Gene & Sample=='C2')$Gene,
                               data=avg_reps_adding_pool_AGG_01_1d, Sample='C2', dgrp_colname='time6', Mode=c('LD', 'DD'), N_days=6, COUNTS=TRUE,
                               NO_X_TEXT=TRUE, NO_Y_TEXT=FALSE, NO_GRIDLINES = TRUE, PRINT_GENES = FALSE, LEGEND=FALSE, offset = -0.5, ALIGN=FALSE, 
                               gene_linesize = 0.5, gene_alpha = 0.7, mean_linesize=0.5, 
                               include_SE_bars = TRUE, rib_alpha=0.8, errorWidth=0, errorThick=0.5,
                               PRISM=TRUE, GENECOLOR=FALSE) +
            theme(strip.text.y = element_text(margin = margin(0,1/72,0,1/72, "in")))) %>%
  ggsave(filename='reps_adding_pool/4_plots/collective/ARGs_strongest_cyc.pdf', height=3, width=4)


####    ####    ####    ####    ####    ####    ####    ####    ####    ####   


# 2.7. ARG timescales - NOT USED
################################################################################

# Take FC against 0
#   subset for our ARGs (before converted)
#   using median score and/or least score, make heatmap across time of FC/LFC
# Choose threshold (% of max) to categorize genes by
#   > plot waveforms by subset of genes (by on time required)

merge(ARGs_DF, passQC) %>% 
  #mutate(Time = as.character(Time)) %>%
  ggplot(aes(Time, FPKM)) +
  stat_summary(geom='errorbar', fun.data='mean_se') +
  facet_grid(Gene~Stimulus , scales='free')

ARGs_DF %>%
  merge(passQC) %>% 
  group_by(Gene, Stimulus, Rep) %>%
  left_join(mutate(filter(., Time == 0), baseline_FPKM = FPKM, Time=NULL,GSM=NULL, .keep='unused')) %>%
  mutate(FC = FPKM/baseline_FPKM) %>% 
  group_by(Time, Gene, Stimulus) %>%
  summarise(FC = mean(FC)) %>% 
  group_by(Gene) %>%
  filter(FC >= 1.5) %>% 
  filter(Time == min(Time)) %>%
  subset(select = c(Time, Gene)) %>%
  unique %>% {table(.$Time)}

ARGs_DF %>%
  merge(passQC) %>% 
  group_by(Gene, Stimulus, Rep) %>%
  left_join(mutate(filter(., Time == 0), baseline_FPKM = FPKM, Time=NULL,GSM=NULL, .keep='unused')) %>%
  mutate(FC = FPKM/baseline_FPKM) %>% 
  group_by(Time, Gene, Stimulus) %>%
  summarise(FC = mean(FC)) %>% 
  ungroup %>%
  mutate(L2FC = log2(FC)) %>%
  group_by(Gene) %>%
  filter(L2FC >= 0.5*max(L2FC)) %>%
  filter(Time == min(Time)) %>% 
  subset(select = c(Time, Gene)) %>%
  unique %>% {table(.$Time)}


ARGs_hmap <- ARGs_DF %>%
  merge(passQC) %>% 
  group_by(Gene, Stimulus, Rep) %>%
  left_join(mutate(filter(., Time == 0), baseline_FPKM = FPKM, Time=NULL,GSM=NULL, .keep='unused')) %>%
  mutate(FC = FPKM/baseline_FPKM) %>% 
  group_by(Time, Gene, Stimulus) %>%
  summarise(FC = mean(FC)) %>% 
  group_by(Time, Gene) %>%
  summarise(FC = median(FC)) %>%
  pivot_wider(names_from = Time, values_from = FC) %>% 
  #{.[, colSums(!is.na(.)) > 0 & !(colnames(.) == "Stimulus")] %>% 
  subset(select=-`90`) %>% 
  as.data.frame %>%
  {`rownames<-`(., .$Gene)} %>%
  mutate(Gene = NULL) %>% 
  mutate_all(log2) %>% 
  as.matrix %>%
  Heatmap(cluster_columns = FALSE, show_row_dend = FALSE, name='Median',
          col=colorRamp2(c(0, 1, max(.)), c("white", "blue", 'red'))) %>%
  draw

ARGs_DF %>%
  merge(passQC) %>% 
  group_by(Gene, Stimulus, Rep) %>%
  left_join(mutate(filter(., Time == 0), baseline_FPKM = FPKM, Time=NULL,GSM=NULL, .keep='unused')) %>%
  mutate(FC = FPKM/baseline_FPKM) %>%
  group_by(Time, Gene, Stimulus) %>%
  summarise(FC = mean(FC)) %>%
  pivot_wider(names_from = Time, values_from = FC) %>%
  
  {split(., .$Stimulus)} %>%
  lapply(function(DF) {
    Name = unique(DF$Stimulus)
    DF[, colSums(!is.na(DF)) > 0 & !(colnames(DF) == "Stimulus")] %>% 
      as.data.frame %>%
      {`rownames<-`(., .$Gene)} %>%
      mutate(Gene = NULL) %>% 
      mutate_all(log2) %>%
      {.[row_order(ARGs_hmap), ]} %>%
      as.matrix %>%
      {Heatmap(., cluster_columns = FALSE, cluster_rows = FALSE, name=Name,
               col=colorRamp2(c(0, 1, max(.)), c("white", "blue", 'red')))} %>%
      {grid.grabExpr(draw(.))}
  }) %>% 
  append(list(grid.grabExpr(draw(ARGs_hmap)))) %>% 
  plot_grid(plotlist=.)




# scaled:
ARGs_hmap_scaled <- ARGs_DF %>%
  merge(passQC) %>% 
  group_by(Gene, Stimulus, Rep) %>%
  left_join(mutate(filter(., Time == 0), baseline_FPKM = FPKM, Time=NULL,GSM=NULL, .keep='unused')) %>%
  mutate(FC = FPKM/baseline_FPKM) %>% 
  group_by(Time, Gene, Stimulus) %>%
  summarise(FC = mean(FC)) %>% 
  group_by(Time, Gene) %>%
  summarise(FC = median(FC)) %>%
  pivot_wider(names_from = Time, values_from = FC) %>% 
  #{.[, colSums(!is.na(.)) > 0 & !(colnames(.) == "Stimulus")] %>% 
  subset(select=-`90`) %>% 
  as.data.frame %>%
  {`rownames<-`(., .$Gene)} %>%
  mutate(Gene = NULL) %>% 
  mutate_all(log2) %>% 
  as.matrix %>% 
  apply(1, function(x) x/max(x)) %>% t %>%
  {.[row_order(ARGs_hmap), ]} %>%
  Heatmap(cluster_columns = FALSE, show_row_dend=FALSE, name='Median',
          col=colorRamp2(c(0, 1/2, 1), c("blue", "yellow", 'red'))) %>%
  draw


ARGs_DF %>%
  merge(passQC) %>% 
  group_by(Gene, Stimulus, Rep) %>%
  left_join(mutate(filter(., Time == 0), baseline_FPKM = FPKM, Time=NULL,GSM=NULL, .keep='unused')) %>%
  mutate(FC = FPKM/baseline_FPKM) %>%
  group_by(Time, Gene, Stimulus) %>%
  summarise(FC = mean(FC)) %>%
  pivot_wider(names_from = Time, values_from = FC) %>%
  
  {split(., .$Stimulus)} %>%
  lapply(function(DF) {
    Name = unique(DF$Stimulus)
    DF[, colSums(!is.na(DF)) > 0 & !(colnames(DF) == "Stimulus")] %>% 
      as.data.frame %>%
      {`rownames<-`(., .$Gene)} %>%
      mutate(Gene = NULL) %>% 
      mutate_all(log2) %>%
      apply(1, function(x) x/max(x)) %>% t %>%
      {.[row_order(ARGs_hmap_scaled), ]} %>%
      as.matrix %>%
      {Heatmap(., cluster_columns = FALSE, cluster_rows = FALSE, name=Name,
               col=colorRamp2(c(0, 1/2, 1), c("blue", "yellow", 'red'))) %>%
          {grid.grabExpr(draw(.))}
      }}) %>% 
  append(list(grid.grabExpr(draw(ARGs_hmap_scaled)))) %>% 
  plot_grid(plotlist=.)


ARGs_DF %>%
  merge(passQC) %>% 
  group_by(Gene, Stimulus, Rep) %>%
  left_join(mutate(filter(., Time == 0), baseline_FPKM = FPKM, Time=NULL,GSM=NULL, .keep='unused')) %>%
  mutate(FC = FPKM/baseline_FPKM) %>%
  group_by(Time, Gene, Stimulus) %>%
  summarise(FC = mean(FC)) %>%
  pivot_wider(names_from = Time, values_from = FC) %>%
  
  {split(., .$Stimulus)} %>%
  lapply(function(DF) {
    Name = unique(DF$Stimulus)
    DF[, colSums(!is.na(DF)) > 0 & !(colnames(DF) == "Stimulus")] %>% 
      as.data.frame %>%
      {`rownames<-`(., .$Gene)} %>%
      mutate(Gene = NULL) %>% 
      mutate_all(log2) %>% 
      as.matrix %>% 
      {.[rowMax(.) > 0, ]} %>% 
      apply(1, function(x) x/max(x)) %>% t %>%
      {Heatmap(., cluster_columns = FALSE, show_row_dend = FALSE, name=Name,
               col=colorRamp2(c(0, 1/2, 1), c("blue", "yellow", 'red'))) %>%
          {grid.grabExpr(draw(.))}
      }}) %>%
  plot_grid(plotlist=.)



# Quickness:

plot_grid(
  ARGs_DF %>%
    merge(passQC) %>% 
    group_by(Gene, Stimulus, Rep) %>%
    left_join(mutate(filter(., Time == 0), baseline_FPKM = FPKM, Time=NULL,GSM=NULL, .keep='unused')) %>%
    mutate(FC = FPKM/baseline_FPKM, 
           delta_FPKM = FPKM-baseline_FPKM) %>%
    group_by(Gene, Stimulus, Rep) %>%
    mutate(normFC = FC/max(FC), 
           norm_delta_FPKM = delta_FPKM/max(delta_FPKM), 
           TimeScore = (max(Time)-Time)/max(Time)) %>%
    mutate(Score = max(0, norm_delta_FPKM*TimeScore) %>% replace_na(0)) %>% 
    group_by(Gene, Stimulus) %>%
    summarise(Score = sum(Score)) %>%
    group_by(Gene) %>% 
    mutate(meanScore = mean(Score)) %>%
    group_by(Stimulus) %>%
    mutate(Rank = rank(Score, ties.method = 'first'), 
           `mean score rank` = rank(meanScore, ties.method = 'first'), 
           .keep='unused') %>% 
    pivot_wider(names_from = Stimulus, values_from = Rank) %>%
    as.data.frame %>%
    `rownames<-`(., .$Gene) %>%
    mutate(Gene = NULL) %>%
    arrange(`mean score rank`) %>%
    as.matrix %>%
    Heatmap(cluster_columns = FALSE, cluster_rows=FALSE, name='Quickness Ranking') %>%
    {grid.grabExpr(draw(.))}, 
  ARGs_quickness_rankings %>%
    merge(Filtereds$DFs$reps_adding_pool) %>%
    ggplot(aes(`mean score rank`)) +
    geom_histogram() +
    facet_grid(Sample ~ Mode)
)

ARGs_quickness_rankings %>%
  mutate(nlog10HMP = -log10(HMP)) %>%
  ggplot(aes(`mean score rank`, nlog10HMP)) +
  geom_point() + 
  facet_grid(Sample~Mode) +
  geom_smooth()

####    ####    ####    ####    ####    ####    ####    ####    ####    ####   


# 2.4. Paper plots - 5
################################################################################
make_standard_small = function(Gene, Sample=c('C2', 'C3'), Mode=c('LD', 'DD'), filePath='', Data=avg_reps_adding_pool_AGG_1_3d, Group0='time3', COUNTS0=TRUE) {
  fileName = file.path(filePath, paste0(Gene, '_', paste(Sample, collapse=''), paste(Mode, collapse=''), '.pdf') )
  
  p = make_plot(Genes = Gene, S = avg_reps_adding_pool_AGG_1_3d, Mode0=Mode, Type0=Sample, Group = 'time3', 
            COUNTS=TRUE, point_size = 0.75, line_size=0.6, PRISM=FALSE, CLASSIC=FALSE)

  ggsave(filename=fileName, plot=p, width=1.05, height=0.8)
}

make_error_small = function(Gene, Sample=c('C2', 'C3'), Mode=c('LD', 'DD'), Y10_0=FALSE, Data=avg_data_ALL_adding_pool, filePath='', name_add='Error_', hardYs=NULL) {
  fileName = file.path(filePath, paste0(name_add, Gene, '_', paste(Sample, collapse=''), paste(Mode, collapse=''), '.pdf') )
  
  p = repplotter_ofSummary(DF=Data, Gene0=Gene, Sample0=Sample, Mode0=Mode, 
                           include_SE = TRUE, include_SE_bars = TRUE, ErrorFun="mean_se", hard_ys = hardYs,
                           mean_point_size = 0, mean_linesize = 0.65, errorWidth = 0, errorThick=0.35, rib_alpha=0.5)
  
  ggsave(filename=fileName, plot=p, width=1.05, height=0.8)
}


# 5A, dotplots of CCG/Pdfr
DotPlot(data_C2C3, features=c('tim', 'per', 'Clk', 'cwo', 'vri', 'Pdp1', 'cry', 'Pdfr', 'cyc'), cols = c('white', 'red'), scale = FALSE, scale.min=0, scavle.max=100)


# 5B, phase plots
# phase plots filtered_t3Q5
{plot_grid(nrow=1, 
          hist_phase(Filtereds$DFs$reps_adding_pool , CIRCULAR=FALSE, TRIM_TO_24=TRUE, binWidth=2, n_W=0, YSCALE=TRUE, Colors=rep('red', 4) %>% `names<-`(c('C3_LD', 'C3_DD', 'C2_LD', 'C2_DD')), HIST_FRONT=TRUE, LD=TRUE, to_select='meta2d_phase'), 
          hist_phase(Filtereds$DFs$reps_adding_pool , CIRCULAR=TRUE, TRIM_TO_24=TRUE, binWidth=2, n_W=0, Colors=rep('red', 4) %>% `names<-`(c('C3_LD', 'C3_DD', 'C2_LD', 'C2_DD')), HIST_FRONT=TRUE, LD=TRUE, to_select='meta2d_phase'), 
          hmap_phase(Filtereds$DFs$reps_adding_pool , Colors = as.list(rep('red', 4)) %>% `names<-`(rev(c('C3_LD', 'C3_DD', 'C2_LD', 'C2_DD'))), col_low = "black", avg_data_time = avg_data_adding_pool_AGG_01_1d)) %>%
    {ggsave(file='reps_adding_pool/4_plots/Phases_Light_alt.pdf', width=5, height=5)}
  plot_grid(nrow=1, 
          hist_phase(Filtereds$DFs$reps_adding_pool , CIRCULAR=FALSE, TRIM_TO_24=TRUE, binWidth=2, n_W=0, YSCALE=TRUE, Colors=rep('red', 4) %>% `names<-`(c('C3_LD', 'C3_DD', 'C2_LD', 'C2_DD')), HIST_FRONT=TRUE, LD=FALSE, to_select='meta2d_phase'), 
          hist_phase(Filtereds$DFs$reps_adding_pool , CIRCULAR=TRUE, TRIM_TO_24=TRUE, binWidth=2, n_W=0, Colors=rep('red', 4) %>% `names<-`(c('C3_LD', 'C3_DD', 'C2_LD', 'C2_DD')), HIST_FRONT=TRUE, LD=FALSE, to_select='meta2d_phase'), 
          hmap_phase(Filtereds$DFs$reps_adding_pool , Colors = as.list(rep('red', 4)) %>% `names<-`(rev(c('C3_LD', 'C3_DD', 'C2_LD', 'C2_DD'))), col_low = "black", avg_data_time = avg_data_adding_pool_AGG_01_1d)) %>%
    {ggsave(file='reps_adding_pool/4_plots/Phases_Dark_alt.pdf', width=5, height=5)}
}

{
a_plot = hist_phase(Filtereds$DFs$reps_adding_pool , CIRCULAR=FALSE, TRIM_TO_24=TRUE, binWidth=2, n_W=0, YSCALE=TRUE, 
            Colors=rep('red', 4) %>% `names<-`(c('C3_LD', 'C3_DD', 'C2_LD', 'C2_DD')), HIST_FRONT=TRUE, LD=TRUE, to_select='meta2d_phase') +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave(a_plot, file='reps_adding_pool/4_plots/Phases_Light_short.pdf', width=1.5, height=3)

a_plot = hist_phase(Filtereds$DFs$reps_adding_pool , CIRCULAR=FALSE, TRIM_TO_24=TRUE, binWidth=2, n_W=0, YSCALE=TRUE, 
           Colors=rep('red', 4) %>% `names<-`(c('C3_LD', 'C3_DD', 'C2_LD', 'C2_DD')), HIST_FRONT=TRUE, LD=FALSE, to_select='meta2d_phase') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave(a_plot, file='reps_adding_pool/4_plots/Phases_Dark_short.pdf', width=1.5, height=3)

rm(a_plot)
}

avg_reps_adding_pool_AGG_01_1d

# 5C - Groups cycling: Ribosomes cell, ribosomes mito, ribosomes all, ETC, ARGs
lapply(c('base', 'prism'), function(Format) {
  Subset_Stats$subsets[c('GO_ribo_noPre', 'GO_ribo_cell', 'GO_ribo_mito', 'ETC_mi', 'ETC_NOmi', 'ARGs', 'ARGs2FC2x', 'ARGs_main', 'GPCR_cyc', 'GCPR_4', 'GPCR_5', 'GPCR_6')] %>%
    {lapply(names(.), function(name) {
      Ps = subset(Filtereds$DFs$reps_adding_pool, Gene %in% .[[name]]) %>% 
        {split(., list(.$Sample, .$Mode))} %>% 
        lapply(function(x) {if(nrow(x)) {
          plotter_geneCombine(x$Gene, data=avg_reps_adding_pool_AGG_01_1d, Mode=x[1,'Mode'], Sample=x[1,'Sample'], 
                              dgrp_colname='time3', N_days=3, COUNTS=TRUE, NO_X_TEXT=TRUE, NO_Y_TEXT=FALSE, 
                              NO_GRIDLINES = TRUE, PRINT_GENES = FALSE, LEGEND=FALSE, offset = -0.5, ALIGN=FALSE, 
                              gene_linesize = 0.5, gene_alpha = 0.7, mean_linesize=0.5, rib_alpha=0.8, PRISM=ifelse(Format=='prism', TRUE, FALSE), GENECOLOR=FALSE) +
            #theme(strip.text.x = element_text(margin = margin(1/72,0,1/72,0, "in"))) +
            theme(strip.text.y = element_text(margin = margin(0,1/72,0,1/72, "in")))
        } else return(NULL)
        }) 
      
      Ps %>% rev %>%
        {lapply(names(.), function(S_M_name) {
          #lapply(c('.svg', '.pdf', '.jpg'), function(extension) ggsave(paste0('Comb_', name, '_', S_M_name, '_CLEAN', extension), plot=.[[S_M_name]], width=1.6,height=1.6))
          #ggsave(paste0('Comb_', name, '_', S_M_name, '_CLEAN.pdf'), plot=.[[S_M_name]], width=1.6,height=1.6)
          ggsave(paste0('reps_adding_pool/4_plots/Not_printing_genes/', Format, '/', name, '_', S_M_name, '.pdf'), plot=.[[S_M_name]], width=1.9,height=1.5)
        }
        )}
    })}
})

# 5C - Groups cycling: Ribosomes cell, ribosomes mito, ribosomes all, ETC, ARGs
lapply(c('base', 'prism'), function(Format) {
  Subset_Stats$subsets[c('GO_ribo_noPre', 'GO_ribo_cell', 'GO_ribo_mito', 'ETC_mi', 'ETC_NOmi', 'ARGs', 'ARGs2FC2x', 'ARGs_main', 'GPCR_cyc', 'GCPR_4', 'GPCR_5', 'GPCR_6')] %>%
    {lapply(names(.), function(name) {
      Ps = subset(Filtereds$DFs$reps_adding_pool, Gene %in% .[[name]]) %>% 
        {split(., list(.$Sample, .$Mode))} %>% 
        lapply(function(x) {if(nrow(x)) {
          plotter_geneCombine(x$Gene, data=avg_reps_adding_pool_AGG_01_1d, Mode=x[1,'Mode'], Sample=x[1,'Sample'], 
                              dgrp_colname='time3', N_days=3, COUNTS=TRUE, NO_X_TEXT=TRUE, NO_Y_TEXT=FALSE, 
                              NO_GRIDLINES = TRUE, PRINT_GENES = TRUE, LEGEND=FALSE, offset = -0.5, ALIGN=FALSE, 
                              gene_linesize = 0.2, gene_alpha = 0.2, mean_linesize=0.5, rib_alpha=0.8, PRISM=ifelse(Format=='prism', TRUE, FALSE), GENECOLOR=FALSE) +
            #theme(strip.text.x = element_text(margin = margin(1/72,0,1/72,0, "in"))) +
            theme(strip.text.y = element_text(margin = margin(0,1/72,0,1/72, "in")))
        } else return(NULL)
        }) 
      
      Ps %>% rev %>%
        {lapply(names(.), function(S_M_name) {
          #lapply(c('.svg', '.pdf', '.jpg'), function(extension) ggsave(paste0('Comb_', name, '_', S_M_name, '_CLEAN', extension), plot=.[[S_M_name]], width=1.6,height=1.6))
          #ggsave(paste0('Comb_', name, '_', S_M_name, '_CLEAN.pdf'), plot=.[[S_M_name]], width=1.6,height=1.6)
          ggsave(paste0('reps_adding_pool/4_plots/Printing_genes/', Format, '/', name, '_', S_M_name, '.pdf'), plot=.[[S_M_name]], width=1.9,height=1.5)
        }
        )}
    })}
})

# 5C - C3LD Individual
{
# top 3 ARGs
lapply(c('Hr38', 'sr', 'CG14186'), function(G) {
  make_standard_small(G, 'C3', 'LD', 'reps_adding_pool/4_plots/indivGenes/C3LD/noErr')
  make_error_small(G, 'C3', 'LD', Data=avg_data_ALL_adding_pool, 
                   filePath='reps_adding_pool/4_plots/indivGenes/C3LD/Err', name_add='Err_')
  make_error_small(G, 'C3', 'LD', Data=avg_data_ALL_adding_pool_1d_01, 
                   filePath='reps_adding_pool/4_plots/indivGenes/C3LD/ErrNorm1d01', name_add='ErrNorm1d01_', hardYs=c(0,1))
  make_error_small(G, 'C3', 'LD', Data=avg_data_ALL_adding_pool_3d_01, 
                   filePath='reps_adding_pool/4_plots/indivGenes/C3LD/ErrNorm3d01', name_add='ErrNorm3d01_', hardYs=c(0,1))
  make_error_small(G, 'C3', 'LD', Data=avg_data_ALL_adding_pool_3d_1, 
                   filePath='reps_adding_pool/4_plots/indivGenes/C3LD/ErrNorm3d1', name_add='ErrNorm3d1_', hardYs=c(0,1))
})

c('CrebA', 'Rala', 'mt:CoIII', 'mRpL48', 'fru', 'lncRNA:CR42861') %>% {.[. %!in% subset(Filtereds$DFs$reps_adding_pool, Sample=='C3' & Mode=='LD')$Gene ]}
# C3LD Genes:
lapply(c('CrebA', 'Rala', 'mt:CoIII', 'mRpL48', 'fru', 'MP1', 'Mes4', 'RpS5a', 'mt:ATPase6'), function(G) {
  make_standard_small(G, 'C3', 'LD', 'reps_adding_pool/4_plots/indivGenes/C3LD/noErr')
  make_error_small(G, 'C3', 'LD', Data=avg_data_ALL_adding_pool, 
                   filePath='reps_adding_pool/4_plots/indivGenes/C3LD/Err', name_add='Err_')
  make_error_small(G, 'C3', 'LD', Data=avg_data_ALL_adding_pool_1d_01, 
                   filePath='reps_adding_pool/4_plots/indivGenes/C3LD/ErrNorm1d01', name_add='ErrNorm1d01_', hardYs=c(0,1))
  make_error_small(G, 'C3', 'LD', Data=avg_data_ALL_adding_pool_3d_01, 
                   filePath='reps_adding_pool/4_plots/indivGenes/C3LD/ErrNorm3d01', name_add='ErrNorm3d10_', hardYs=c(0,1))
  make_error_small(G, 'C3', 'LD', Data=avg_data_ALL_adding_pool_3d_1, 
                   filePath='reps_adding_pool/4_plots/indivGenes/C3LD/ErrNorm3d1', name_add='ErrNorm3d1_', hardYs=c(0,1))
})

# CCG
lapply(c('cwo', 'Pdp1', 'per'), function(G) {
  lapply(c('C2', 'C3'), function(S) {
    lapply(c('LD', 'DD'), function(M) {
      make_standard_small(G, S, M, file.path('reps_adding_pool/4_plots/indivGenes', paste0(S, M), 'noErr'))
      make_error_small(G, S, M, Data=avg_data_ALL_adding_pool, 
                       filePath=file.path('reps_adding_pool/4_plots/indivGenes/CCG', paste0(S, M), 'Err'), name_add='Err_')
      make_error_small(G, S, M, Data=avg_data_ALL_adding_pool_1d_01, 
                       filePath=file.path('reps_adding_pool/4_plots/indivGenes/CCG', paste0(S, M), 'ErrNorm1d01'), name_add='ErrNorm1d01_', hardYs=c(0,1))
      make_error_small(G, S, M, Data=avg_data_ALL_adding_pool_3d_01, 
                       filePath=file.path('reps_adding_pool/4_plots/indivGenes/CCG', paste0(S, M), 'ErrNorm3d01_'), name_add='ErrNorm3d01_', hardYs=c(0,1))
      make_error_small(G, S, M, Data=avg_data_ALL_adding_pool_3d_1, 
                       filePath=file.path('reps_adding_pool/4_plots/indivGenes/CCG', paste0(S, M), 'ErrNorm3d1_'), name_add='ErrNorm3d1_', hardYs=c(0,1))
    })
  })
})



# C3DD Genes:
lapply(c('Hr38', 'sr', "ab", "CG10103", "CG14186", "CG42594", "pho", "CG7139", "CG41099", "Fis1", "hng3", "sty", "Xrp1", 'Eip74EF', 'REPTOR-BP', 'RpL34a'), function(G) {
  make_standard_small(G, 'C3', 'DD', 'reps_adding_pool/4_plots/indivGenes/C3DD/noErr')
  make_error_small(G, 'C3', 'DD', Data=avg_data_ALL_adding_pool, 
                   filePath='reps_adding_pool/4_plots/indivGenes/C3DD/Err', name_add='Err_')
  make_error_small(G, 'C3', 'DD', Data=avg_data_ALL_adding_pool_1d_01, 
                   filePath='reps_adding_pool/4_plots/indivGenes/C3DD/ErrNorm1d01', name_add='ErrNorm1d01_', hardYs=c(0,1))
  make_error_small(G, 'C3', 'DD', Data=avg_data_ALL_adding_pool_3d_01, 
                   filePath='reps_adding_pool/4_plots/indivGenes/C3DD/ErrNorm3d01', name_add='ErrNorm3d10_', hardYs=c(0,1))
  make_error_small(G, 'C3', 'DD', Data=avg_data_ALL_adding_pool_3d_1, 
                   filePath='reps_adding_pool/4_plots/indivGenes/C3DD/ErrNorm3d1', name_add='ErrNorm3d1_', hardYs=c(0,1))
})
}

# C2LD Genes:
lapply(c("Hr38", "sr", "CG14186", "mAChR-A", "Rh7", "spoon", "RpL40", "mRpL28", "mt:ATPase6", "mt:CoI", "mt:CoIII", "mt:Cyt-b"), function(G) {
  make_standard_small(G, 'C2', 'LD', 'reps_adding_pool/4_plots/indivGenes/C2LD/noErr')
  make_error_small(G, 'C2', 'LD', Data=avg_data_ALL_adding_pool, 
                   filePath='reps_adding_pool/4_plots/indivGenes/C2LD/Err', name_add='Err_')
  make_error_small(G, 'C2', 'LD', Data=avg_data_ALL_adding_pool_1d_01, 
                   filePath='reps_adding_pool/4_plots/indivGenes/C2LD/ErrNorm1d01', name_add='ErrNorm1d01_', hardYs=c(0,1))
  make_error_small(G, 'C2', 'LD', Data=avg_data_ALL_adding_pool_3d_01, 
                   filePath='reps_adding_pool/4_plots/indivGenes/C2LD/ErrNorm3d01', name_add='ErrNorm3d10_', hardYs=c(0,1))
  make_error_small(G, 'C2', 'LD', Data=avg_data_ALL_adding_pool_3d_1, 
                   filePath='reps_adding_pool/4_plots/indivGenes/C2LD/ErrNorm3d1', name_add='ErrNorm3d1_', hardYs=c(0,1))
})

# "Hr38", "sr", "CG14186", "mAChR-A", "Rh7", "msi", "RpL40", "mRpL28", "mt:ATPase", "mt:CoI", "mt:CoIII", "mt:Cyt-b"


####    ####    ####    ####    ####    ####    ####    ####    ####    ####   


# 2.5. Paper plots - 5G
################################################################################
Filtereds$DFs$reps_adding_pool %>% subset(Gene %in% Subset_Stats$subsets$ARGs) %>%
  {split(., list(.$Sample, .$Mode))}

subset(to_filter, second_pct>0.1 & Gene %in% ARG$converted$min2, select=c(Gene, Sample, Mode)) %>%
  {split(., list(.$Sample, .$Mode))} %>%
  lapply(nrow)


# GO_ribo_cell
plot_grid(ncol=1, # C3
          plotter_geneCombine( subset(to_filter, second_pct > 0.1 & Gene %in% merge(GO_ribo_cell, passQC)$Gene & Sample=='C3')$Gene,
                               data=avg_reps_adding_pool_AGG_01_1d, Sample='C3', dgrp_colname='time6', Mode=c('LD', 'DD'), N_days=6, COUNTS=TRUE,
                               NO_X_TEXT=TRUE, NO_Y_TEXT=FALSE, NO_GRIDLINES = TRUE, PRINT_GENES = FALSE, LEGEND=FALSE, offset = -0.5, ALIGN=FALSE, 
                               gene_linesize = 0.5, gene_alpha = 0.7, mean_linesize=0.5, 
                               include_SE_bars = TRUE, rib_alpha=0.8, errorWidth=0, errorThick=0.5,
                               PRISM=TRUE, GENECOLOR=FALSE) +
            theme(strip.text.y = element_text(margin = margin(0,1/72,0,1/72, "in"))),
          # C2
          plotter_geneCombine( subset(to_filter, second_pct > 0.1 & Gene %in% merge(GO_ribo_cell, passQC)$Gene & Sample=='C2')$Gene,
                               data=avg_reps_adding_pool_AGG_01_1d, Sample='C2', dgrp_colname='time6', Mode=c('LD', 'DD'), N_days=6, COUNTS=TRUE,
                               NO_X_TEXT=TRUE, NO_Y_TEXT=FALSE, NO_GRIDLINES = TRUE, PRINT_GENES = FALSE, LEGEND=FALSE, offset = -0.5, ALIGN=FALSE, 
                               gene_linesize = 0.5, gene_alpha = 0.7, mean_linesize=0.5, 
                               include_SE_bars = TRUE, rib_alpha=0.8, errorWidth=0, errorThick=0.5,
                               PRISM=TRUE, GENECOLOR=FALSE) +
            theme(strip.text.y = element_text(margin = margin(0,1/72,0,1/72, "in")))) %>%
  ggsave(filename='reps_adding_pool/4_plots/collective/GO_ribo.pdf', height=3, width=4)


# Ribo no pre
plot_grid(ncol=1, # C3
          plotter_geneCombine( subset(to_filter, second_pct > 0.1 & Gene %in% Subset_Stats$subsets$GO_ribo_noPre & Sample=='C3')$Gene,
                               data=avg_reps_adding_pool_AGG_01_1d, Sample='C3', dgrp_colname='time6', Mode=c('LD', 'DD'), N_days=6, COUNTS=TRUE,
                               NO_X_TEXT=TRUE, NO_Y_TEXT=FALSE, NO_GRIDLINES = TRUE, PRINT_GENES = FALSE, LEGEND=FALSE, offset = -0.5, ALIGN=FALSE, 
                               gene_linesize = 0.5, gene_alpha = 0.7, mean_linesize=0.5, 
                               include_SE_bars = TRUE, rib_alpha=0.8, errorWidth=0, errorThick=0.5,
                               PRISM=TRUE, GENECOLOR=FALSE) +
            theme(strip.text.y = element_text(margin = margin(0,1/72,0,1/72, "in"))),
          # C2
          plotter_geneCombine( subset(to_filter, second_pct > 0.1 & Gene %in% Subset_Stats$subsets$GO_ribo_noPre & Sample=='C2')$Gene,
                               data=avg_reps_adding_pool_AGG_01_1d, Sample='C2', dgrp_colname='time6', Mode=c('LD', 'DD'), N_days=6, COUNTS=TRUE,
                               NO_X_TEXT=TRUE, NO_Y_TEXT=FALSE, NO_GRIDLINES = TRUE, PRINT_GENES = FALSE, LEGEND=FALSE, offset = -0.5, ALIGN=FALSE, 
                               gene_linesize = 0.5, gene_alpha = 0.7, mean_linesize=0.5, 
                               include_SE_bars = TRUE, rib_alpha=0.8, errorWidth=0, errorThick=0.5,
                               PRISM=TRUE, GENECOLOR=FALSE) +
            theme(strip.text.y = element_text(margin = margin(0,1/72,0,1/72, "in")))) %>%
  ggsave(filename='reps_adding_pool/4_plots/collective/GO_ribo_noPre.pdf', height=3, width=4)

#ARGs
plot_grid(ncol=1, # C3
          plotter_geneCombine( subset(to_filter, second_pct > 0.1 & Gene %in% ARG$converted$min2 & Sample=='C3')$Gene,
                               data=avg_reps_adding_pool_AGG_01_1d, Sample='C3', dgrp_colname='time6', Mode=c('LD', 'DD'), N_days=6, COUNTS=TRUE,
                               NO_X_TEXT=TRUE, NO_Y_TEXT=FALSE, NO_GRIDLINES = TRUE, PRINT_GENES = FALSE, LEGEND=FALSE, offset = -0.5, ALIGN=FALSE, 
                               gene_linesize = 0.5, gene_alpha = 0.7, mean_linesize=0.5, 
                               include_SE_bars = TRUE, rib_alpha=0.8, errorWidth=0, errorThick=0.5,
                               PRISM=TRUE, GENECOLOR=FALSE) +
            theme(strip.text.y = element_text(margin = margin(0,1/72,0,1/72, "in"))),
          # C2
          plotter_geneCombine( subset(to_filter, second_pct > 0.1 & Gene %in% ARG$converted$min2 & Sample=='C2')$Gene,
                               data=avg_reps_adding_pool_AGG_01_1d, Sample='C2', dgrp_colname='time6', Mode=c('LD', 'DD'), N_days=6, COUNTS=TRUE,
                               NO_X_TEXT=TRUE, NO_Y_TEXT=FALSE, NO_GRIDLINES = TRUE, PRINT_GENES = FALSE, LEGEND=FALSE, offset = -0.5, ALIGN=FALSE, 
                               gene_linesize = 0.5, gene_alpha = 0.7, mean_linesize=0.5, 
                               include_SE_bars = TRUE, rib_alpha=0.8, errorWidth=0, errorThick=0.5,
                               PRISM=TRUE, GENECOLOR=FALSE) +
            theme(strip.text.y = element_text(margin = margin(0,1/72,0,1/72, "in")))) %>%
  ggsave(filename='reps_adding_pool/4_plots/collective/ARGs.pdf', height=3, width=4)


lapply(unique(ARGs_strongest_cyc$Gene) %>% {.[. %in% passQC$Gene]}, function(G) {
  repplotter_ofSummary(DF=avg_data_ALL_adding_pool, Gene0=G, Sample0='C2', Mode0='DD', 
                       include_SE = TRUE, include_SE_bars = TRUE, ErrorFun="mean_se", hard_ys = NULL,
                       mean_point_size = 0, mean_linesize = 0.65, errorWidth = 0, errorThick=0.35, rib_alpha=0.5) +
    ggtitle(G)
}) %>%
  {plot_grid(plotlist=.)}





# ARG_RP - generating collective stats
{ARG_RP <- list(waveforms = avg_data_ALL_adding_pool_1d_01 %>% 
                 merge(rbind(data.frame(Gene = ARG$converted$min2, Group = 'ARG'),
                             data.frame(Gene = merge(GO_ribo_cell, passQC)$Gene, Group = 'CellRibo')
                             #data.frame(Gene = unique(ARGs_strongest_cyc$Gene), Group = 'ARGs2FC2x')
                       )) %>%
                 merge(subset(to_filter, second_pct>0.1, select=c(Gene, Sample, Mode))) %>%
                 group_by(Sample, Mode, Group, Time) %>%
                 summarise(value = mean(value)))


# meta2d
ARG_RP$waveforms %>%
  pivot_wider(names_from=Time, values_from=value) %>%
  as.data.frame %>%
  {`rownames<-`(., paste0(.$Group, '_', .$Sample, '-', .$Mode))} %>%
  mutate(Sample=NULL, Mode=NULL, Group=NULL) %>%
  write.csv(file = "reps_adding_pool/comp_waveforms.csv", row.names=TRUE)

ARG_RP$waveforms_1dx9 <- avg_data_ALL_adding_pool_1dx9_01 %>%
  merge(rbind(data.frame(Gene = ARG$converted$min2, Group = 'ARG'),
              data.frame(Gene = merge(GO_ribo_cell, passQC)$Gene, Group = 'CellRibo')
              #data.frame(Gene = unique(ARGs_strongest_cyc$Gene), Group = 'ARGs2FC2x')
              )) %>%
  merge(subset(to_filter, second_pct>0.1, select=c(Gene, Sample, Mode)))

ARG_RP$meta_pvals <- meta2d(infile ='reps_adding_pool/comp_waveforms.csv', outputFile=FALSE, filestyle='csv', timepoints = "line1", minper=20, maxper=28)$meta %>%
  mutate(Group = gsub('_.*', '', CycID), 
         Sample = gsub('-.*', '', CycID) %>% gsub('.*_', '', .), 
         Mode = gsub('.*-', '', CycID), .keep='unused') %>%
  subset(select = c(Group, Sample, Mode, meta2d_pvalue, meta2d_BH.Q, meta2d_phase))

ARG_RP$KW_3 <- ARG_RP$waveforms %>%
  mutate(Time = Time%%24) %>%
  {split(., list(.$Group, .$Sample, .$Mode))} %>%
  lapply(function(x) {
    kruskal.test(value ~ Time, data=x) %>%
      {data.frame(Group = x[1,'Group'], Sample = x[1,'Sample'], Mode = x[1,'Mode'], KW_pvalue_3 = .$p.value)}
  }) %>%
  do.call(what='rbind') %>%
  mutate(KW_BH.Q_3 = p.adjust(KW_pvalue_3, method='BH'))
  
ARG_RP$KW_9 <- ARG_RP$waveforms_1dx9 %>%
  {split(., list(.$Group, .$Sample, .$Mode))} %>%
  lapply(function(x) {
    kruskal.test(value ~ Time, data=x) %>%
      {data.frame(Group = x[1,'Group'], Sample = x[1,'Sample'], Mode = x[1,'Mode'], KW_pvalue_9 = .$p.value)}
  }) %>%
  do.call(what='rbind') %>%
  mutate(KW_BH.Q_9 = p.adjust(KW_pvalue_9, method='BH'))

ARG_RP$DayNight_3 <- ARG_RP$waveforms %>%
  mutate(Time = Time%%24) %>%
  mutate(DayNight = Time%/%12) %>%
  {split(., list(.$Group, .$Sample, .$Mode))} %>%
  lapply(function(x) {
    wilcox.test(value ~ DayNight, data=x) %>%
      {data.frame(Group = x[1,'Group'], Sample = x[1,'Sample'], Mode = x[1,'Mode'], DayNight_pvalue_3 = .$p.value)}
  }) %>% do.call(what='rbind') %>%
  mutate(DayNight_BH.Q_3 = p.adjust(DayNight_pvalue_3, method='BH'))

ARG_RP$DayNight_9 <- ARG_RP$waveforms_1dx9 %>%
  mutate(DayNight = Time%/%12) %>%
  {split(., list(.$Group, .$Sample, .$Mode))} %>%
  lapply(function(x) {
    wilcox.test(value ~ DayNight, data=x) %>%
      {data.frame(Group = x[1,'Group'], Sample = x[1,'Sample'], Mode = x[1,'Mode'], DayNight_pvalue_9 = .$p.value)}
  }) %>% do.call(what='rbind') %>%
  mutate(DayNight_BH.Q_9 = p.adjust(DayNight_pvalue_9, method='BH'))

ARG_RP$mean_ratio <- ARG_RP$waveforms_1dx9 %>%
  group_by(Sample, Mode, Group, Time, Day, Rep) %>%
  summarise(value=mean(value)) %>%
  group_by(Sample, Mode, Group, Day, Rep) %>%
  summarise(ratio = max(value)/min(value)) %>%
  group_by(Sample, Mode, Group) %>%
  summarise(ratio = mean(ratio))

ARG_RP$Fisher_p <- Filtereds$p_vals$reps_adding_pool[c('ARGs', 'GO_ribo_cell'), ] %>%
  t() %>% as.data.frame %>%
  {`colnames<-`(., c('ARG', 'CellRibo'))} %>%
  {mutate(S_M = rownames(.), .)} %>%
  {`rownames<-`(., NULL)} %>%
  mutate(Mode = gsub('.*_', '', S_M), 
         Sample = gsub('_.*', '', S_M), .keep='unused') %>%
  pivot_longer(cols=c(ARG, CellRibo), names_to='Group', values_to='Fisher_p') %>%
  mutate(Fisher_BH.Q = p.adjust(Fisher_p, method='BH'))

ARG_RP$mean_r <- substatter(Filtereds$DFs$expressed, R_NOT_P=TRUE, Subsets = Subset_Stats$subsets[c('ARGs', 'GO_ribo_cell')]) %>%
  t() %>% as.data.frame %>%
  {`colnames<-`(., c('ARG', 'CellRibo'))} %>%
  {mutate(S_M = rownames(.), .)} %>%
  {`rownames<-`(., NULL)} %>%
  mutate(Mode = gsub('.*_', '', S_M), 
         Sample = gsub('_.*', '', S_M), .keep='unused') %>%
  pivot_longer(cols=c(ARG, CellRibo), names_to='Group', values_to='mean_r')

ARG_RP$DF = Reduce(merge, ARG_RP[!grepl('waveform|DF', names(ARG_RP))])

ARG_RP$DF %>%
  {`rownames<-`(., paste(.$Sample, .$Mode, .$Group))} %>%
  mutate(Group=NULL, Sample=NULL, Mode=NULL) %>%
  {.[, !grepl('9|pvalue|_p', colnames(.))]} %>%
  mutate(across(contains('BH'), ~-log10(.), .names = "nLog10_{.col}"), .keep='unused') %>%
  as.matrix %>%
  t %>% {./rowMax(.)} %>% 
  pheatmap(cluster_rows = FALSE,  # Disable clustering of rows
           cluster_cols = FALSE,  # Disable clustering of columns
           display_numbers = TRUE,  # Display numbers in cells
           number_format = "%.2f",  # Format numbers to 2 decimal places
           color = colorRampPalette(c("blue", "white", "red"))(100),  # Custom color scale
           main = "scaled\nvalue",  # Title with a line break
           show_rownames = TRUE,  # Show row names
           show_colnames = TRUE,  # Show column names
           fontsize_row = 10,  # Adjust the font size of row names
           fontsize_col = 10,  # Adjust the font size of column names
           number_color = "black" 
  )
}

ARG_RP$DF %>%
  {.[, !grepl('9|pvalue|_p', colnames(.))]} %>%
  {.[c( 1:4, 7, 9, 5, 6, 8)]} %>%
  {split(., .$Group)}

ggsave(filename='reps_adding_pool/4_plots/collective/ARGs_stats.pdf', height=5, width=11,
       ARG_RP$DF %>%
         {.[, !grepl('9|pvalue|_p|Fisher', colnames(.))]} %>%
         mutate(Group = replace(Group, Group=='CellRibo', 'Cell Ribosomes'),# %>% {replace(Group, Group == 'ARGs2FC2x', 'ARGs-2FC')}, 
                Mode = factor(Mode, levels=c("LD", 'DD')),
                Sample = factor(Sample, levels=c("C2", 'C3')),
                across(contains('BH'), ~-log10(.), .names = "nLog10_{.col}"), .keep='unused') %>%
         pivot_longer(cols=-c(Group, Sample, Mode), names_to='Stat') %>%
         merge(data.frame(Stat = c("ratio", "mean_r", "nLog10_meta2d_BH.Q", "nLog10_KW_BH.Q_3", "nLog10_DayNight_BH.Q_3"), 
                          BetterName = c('< ratio > (across day-replicates)', '< pairwise r > (across genes)', '-log10( meta2d BH-Q value )', '-log10( Kruskal-Wallis BH-Q value) ', '-log10( Day vs Night BH-Q value )'))) %>%
         mutate(Sample_Group=paste(Sample, Group)) %>%
         ggplot(aes(x=Mode, y=value, color=BetterName, group=BetterName)) +
         geom_line() + 
         geom_point() +
         facet_wrap(~Sample_Group, scales='free') +
         #facet_grid(Group~Sample) +
         theme_prism()
)



  
ARG_RP$DF %>%
  {.[, !grepl('9|pvalue|_p|Fisher', colnames(.))]} %>%
  mutate(Mode = factor(Mode, levels=c("LD", 'DD')),
         Sample = factor(Sample, levels=c("C2", 'C3')),
         across(contains('BH'), ~-log10(.), .names = "nLog10_{.col}"), .keep='unused') %>%
  pivot_longer(cols=-c(Group, Sample, Mode), names_to='Stat') %>%
  {split(., list(.$Group, .$Sample))} %>%
  lapply(function(x) {
    wilcox.test(value ~ Mode, data=x) %>%
      {data.frame(Group = x[1,'Group'], Sample = x[1,'Sample'], Stat = .$p.value)}
  }) %>% do.call(what='rbind')



# SCATTER PLOTS of individual (group) stats - ratios not helpful
ARGs_DF_scatter <- avg_data_ALL_adding_pool_1dx9_01 %>%
  group_by(Gene, Time, Sample, Mode) %>%
  summarise(value = mean(value)) %>%
  group_by(Gene, Sample, Mode) %>%
  summarise(amp2 = max(value)-min(value)) %>%
  merge(rbind(cbind(Gene = ARG$converted$min2, Group = 'ARGs'),
              cbind(Gene = Subset_Stats$subsets$GO_ribo_cell, Group = 'Cell ribosome')#,
              #data.frame(Gene = unique(ARGs_strongest_cyc$Gene), Group = 'ARGs_2FC_2x')
              )) %>%
  merge(  reps_adding_pool$res$AGG %>% do.call(what=rbind) %>%
            group_by(Gene, Sample, Mode) %>%
            summarise(HMP = hmp.stat(meta2d_pvalue), MinRatio = max(Ratio), meta2d_phase=Arg(mean(complex(argument = (2*pi*meta2d_phase)/24)))*24/(2*pi)) %>%
            merge(., subset(to_filter, Gene %in% passQC$Gene, select=c(Gene, Sample, Mode)))) %>%
  mutate(`-log10(HMP)` = -log10(HMP))


ggsave(
  ARGs_DF_scatter %>%
    group_by(Group) %>% mutate(mean_phase = mean(meta2d_phase)) %>% ungroup %>%
    mutate(phase_offset_score = (cos( (mean_phase - meta2d_phase) / 24 * 2*pi)+1)/2 ) %>% 
    group_by(Gene, Sample, Group) %>% 
    mutate(phase_offset_score = mean(phase_offset_score), 
           maxNLogHMP = max(`-log10(HMP)`), 
           meanAmp = mean(amp2)) %>% 
    ungroup %>%
    mutate(Group_Sample = paste(Group, Sample)) %>%
    pivot_wider(names_from=Mode, values_from=c(amp2, HMP, MinRatio, meta2d_phase, `-log10(HMP)`)) %>% 
    #subset(HMP_LD < 0.05 | HMP_DD<0.05) %>% 
    ggplot(aes(`-log10(HMP)_DD`, `-log10(HMP)_LD`, color=Group_Sample, alpha = phase_offset_score)) +
    theme_minimal() +
    theme(panel.grid.minor = element_blank(), axis.text=element_text(color='black')) +
    geom_point(aes(size=maxNLogHMP)) +
    geom_text_repel(aes(label = Gene), size = 3, max.overlaps = 9) +  # Repelled text labels
    scale_size_continuous(range=c(0,4)) +
    scale_x_continuous(limits=c(0, 12), expand=c(0.05,0)) +
    scale_y_continuous(limits=c(0, 12), expand=c(0.05,0)) +
    scale_color_manual(values=c('#E69F00', '#0072B2', 'red', 'green')) +
    stat_function(fun = function(x) x, color = "black", size = 0.5) +  # y = x line
    #lims(x=c(0, 1), y=c(0, 1)) +
    #geom_smooth(method = "lm", fullrange = TRUE, se=FALSE) +
    facet_grid(~Group, scales='free'), 
  filename='reps_adding_pool/5_plots/scatter.pdf', width=11, height=5)


ARGs_DF_scatter %>%
  group_by(Group) %>% mutate(mean_phase = mean(meta2d_phase)) %>% ungroup %>%
  mutate(phase_offset_score = (cos( (mean_phase - meta2d_phase) / 24 * 2*pi)+1)/2 ) %>% 
  group_by(Gene, Sample, Group) %>% 
  mutate(phase_offset_score = mean(phase_offset_score), 
         meanNLogHMP = mean(`-log10(HMP)`), 
         meanAmp = mean(amp2)) %>% 
  ungroup %>%
  mutate(Group_Sample = paste(Group, Sample)) %>%
  pivot_wider(names_from=Mode, values_from=c(amp2, HMP, MinRatio, meta2d_phase, `-log10(HMP)`)) %>% 
  #subset(HMP_LD < 0.05 | HMP_DD<0.05) %>% 
  ggplot(aes(`-log10(HMP)_DD`, `-log10(HMP)_LD`, color=Group_Sample, alpha = phase_offset_score)) +
  theme_minimal() +
  theme(panel.grid.minor = element_blank(), axis.text=element_text(color='black')) +
  geom_point(aes(size=meanNLogHMP)) +
  geom_text_repel(aes(label = Gene), size = 3, max.overlaps = 9) +  # Repelled text labels
  scale_size_continuous(range=c(0,4)) +
  scale_alpha(range=c(0,1)) +
  geom_abline(slope=1, intercept=0, color='black', size=0.25) +
  scale_x_continuous(limits=c(0, 12), expand=c(0.05,0)) +
  scale_y_continuous(limits=c(0, 12), expand=c(0.05,0)) +
  #lims(x=c(0, 1), y=c(0, 1)) +
  #geom_smooth(method = "lm", fullrange = TRUE, se=FALSE) +
  facet_grid(~Group, scales='free')

plot_grid(ncol=1, rel_heights = c(1, 0.75, 0.75, 1),
  ARGs_DF_scatter %>%
    group_by(Group) %>% mutate(mean_phase = mean(meta2d_phase)) %>% ungroup %>%
    mutate(phase_offset_score = (cos( (mean_phase - meta2d_phase) / 24 * 2*pi)+1)/2 ) %>% 
    group_by(Gene, Sample, Group) %>% 
    mutate(phase_offset_score = mean(phase_offset_score), 
           meanNLogHMP = mean(`-log10(HMP)`), 
           meanAmp = mean(amp2)) %>% 
    ungroup %>%
    mutate(Group_Sample = paste(Group, Sample)) %>%
    pivot_wider(names_from=Mode, values_from=c(amp2, HMP, MinRatio, meta2d_phase, `-log10(HMP)`)) %>% 
    #subset(HMP_LD < 0.05 | HMP_DD<0.05) %>% 
    ggplot(aes(`-log10(HMP)_DD`, `-log10(HMP)_LD`, color=Group_Sample, alpha = phase_offset_score)) +
    theme_minimal() +
    theme(panel.grid.minor = element_blank(), axis.text=element_text(color='black')) +
    geom_point(aes(size=meanNLogHMP)) +
    geom_text_repel(aes(label = Gene), size = 3, max.overlaps = 9) +  # Repelled text labels
    scale_size_continuous(range=c(0,4)) +
    scale_x_continuous(limits=c(0, 12), expand=c(0.05,0)) +
    scale_y_continuous(limits=c(0, 12), expand=c(0.05,0)) +
    geom_abline(slope=1, intercept=0, color='black', size=0.25) +
    #lims(x=c(0, 1), y=c(0, 1)) +
    #geom_smooth(method = "lm", fullrange = TRUE, se=FALSE) +
    facet_grid(~Group, scales='free'), 
  plot_grid(nrow=1, rel_widths=c(1, 0.05, 1), labels = c("ARGs C2", '', "ARGs C3"), label_x = 0.35, label_y = 1,
            
            plot_grid(ncol=1, rel_heights=c(0.1, 1), ggplot()+theme_void(), 
                      subset(ARGs_DF_scatter, Group == 'ARGs' & Sample == 'C2' & Mode == 'DD') %>% 
                        arrange(HMP) %>%
                        {.[1:9, 'Gene']} %>% 
                        lapply(function(G) {
                          repplotter_ofSummary(G, Sample0='C2', 
                                               DF=mutate(avg_data_ALL_adding_pool, Time = as.integer(Time)+72*as.numeric(Mode=='DD')), 
                                               include_SE_bars=TRUE, rib_alpha=0.5, errorWidth=0, errorThick=0.25, hard_ys=NULL)+
                            ggtitle(G) + theme(plot.title = element_text(face = "bold"))}) %>%
                        plot_grid(plotlist=.)),
            
            ggplot()+theme_void(),
            
            plot_grid(ncol=1, rel_heights=c(0.1, 1), ggplot()+theme_void(), 
                      subset(ARGs_DF_scatter, Group == 'ARGs' & Sample == 'C3' & Mode == 'DD') %>% 
                        arrange(HMP) %>%
                        {.[1:9, 'Gene']} %>% 
                        lapply(function(G) {
                          repplotter_ofSummary(G, Sample0='C3', 
                                               DF=mutate(avg_data_ALL_adding_pool, Time = as.integer(Time)+72*as.numeric(Mode=='DD')),
                                               include_SE_bars=TRUE, rib_alpha=0.5, errorWidth=0, errorThick=0.25, hard_ys=NULL) +
                            ggtitle(G) + theme(plot.title = element_text(face = "bold"))}) %>%
                        plot_grid(plotlist=.))),
  
  plot_grid(nrow=1, rel_widths=c(1, 0.05, 1), 
            
            plot_grid(ncol=1, rel_heights=c(0.1, 1), ggplot()+theme_void(), 
                      subset(ARGs_DF_scatter, Group == 'ARGs' & Sample == 'C2' & Mode == 'DD') %>% 
                        arrange(HMP) %>%
                        {.[1:9, 'Gene']} %>% 
                        lapply(function(G) {
                          repplotter_ofSummary(G, Sample0='C2', 
                                               DF=avg_data_ALL_adding_pool_6d_1, 
                                               include_SE_bars=TRUE, rib_alpha=0.5, errorWidth=0, errorThick=0.25, hard_ys=c(0,1))+
                            ggtitle(G) + theme(plot.title = element_text(face = "bold"))
                        }) %>%
                        plot_grid(plotlist=.)),
            
            ggplot()+theme_void(),
            
            plot_grid(ncol=1, rel_heights=c(0.1, 1), ggplot()+theme_void(), 
                      subset(ARGs_DF_scatter, Group == 'ARGs' & Sample == 'C3' & Mode == 'DD') %>% 
                        arrange(HMP) %>%
                        {.[1:9, 'Gene']} %>% 
                        lapply(function(G) {
                          repplotter_ofSummary(G, Sample0='C3', 
                                               DF=avg_data_ALL_adding_pool_6d_1,
                                               include_SE_bars=TRUE, rib_alpha=0.5, errorWidth=0, errorThick=0.25, hard_ys=c(0,1))+
                            ggtitle(G) + theme(plot.title = element_text(face = "bold"))}) %>%
                        plot_grid(plotlist=.))), 
  ARGs_DF_scatter %>%
    group_by(Group) %>% mutate(mean_phase = mean(meta2d_phase)) %>% ungroup %>%
    mutate(phase_offset_score = (cos( (mean_phase - meta2d_phase) / 24 * 2*pi)+1)/2 ) %>% 
    group_by(Gene, Mode, Group) %>% 
    mutate(phase_offset_score = mean(phase_offset_score), 
           meanNLogHMP = mean(`-log10(HMP)`), 
           meanAmp = mean(amp2)) %>% 
    ungroup %>%
    mutate(Group_Mode = paste(Group, Mode)) %>%
    pivot_wider(names_from=Sample, values_from=c(amp2, HMP, MinRatio, meta2d_phase, `-log10(HMP)`)) %>% 
    #subset(HMP_LD < 0.05 | HMP_DD<0.05) %>% 
    ggplot(aes(`-log10(HMP)_C2`, `-log10(HMP)_C3`, color=Group_Mode, alpha = phase_offset_score)) +
    theme_minimal() +
    theme(panel.grid.minor = element_blank(), axis.text=element_text(color='black')) +
    geom_point(aes(size=meanNLogHMP)) +
    geom_text_repel(aes(label = Gene), size = 3, max.overlaps = 9) +  # Repelled text labels
    scale_size_continuous(range=c(0,4)) +
    scale_x_continuous(limits=c(0, 12), expand=c(0.05,0)) +
    scale_y_continuous(limits=c(0, 12), expand=c(0.05,0)) +
    geom_abline(slope=1, intercept=0, color='black', size=0.25) +
    #lims(x=c(0, 1), y=c(0, 1)) +
    #geom_smooth(method = "lm", fullrange = TRUE, se=FALSE) +
    facet_grid(~Group, scales='free')
) %>% 
  ggsave(filename='reps_adding_pool/5_plots/scatter_full.pdf', width=11, height=12.5)

 # inverting
ARGs_DF_scatter %>%
  subset(Mode == 'DD' & Group == 'ARGs') %>%
  group_by(Group) %>% mutate(mean_phase = mean(meta2d_phase)) %>% ungroup %>%
  mutate(phase_offset_score = (cos( (mean_phase - meta2d_phase) / 24 * 2*pi)+1)/2 ) %>% 
  group_by(Gene, Mode, Group) %>% 
  mutate(phase_offset_score = mean(phase_offset_score), 
         meanNLogHMP = mean(`-log10(HMP)`), 
         meanAmp = mean(amp2)) %>% 
  ungroup %>%
  mutate(Group_Mode = paste(Group, Mode)) %>%
  pivot_wider(names_from=Sample, values_from=c(amp2, HMP, MinRatio, meta2d_phase, `-log10(HMP)`)) %>% 
  ggplot(aes(`-log10(HMP)_C3`, `-log10(HMP)_C2`, color=Group_Mode)) +
  theme_minimal() +
  theme(panel.grid.minor = element_blank(), axis.text=element_text(color='black'), 
        legend.position='null') +
  geom_point(aes(size=meanNLogHMP)) +
  geom_text_repel(aes(label = Gene), size = 3, max.overlaps = 9) +  # Repelled text labels
  scale_size_continuous(range=c(0,4)) +
  scale_x_continuous(limits=c(0, 6), expand=c(0.05,0)) +
  scale_y_continuous(limits=c(0, 6), expand=c(0.05,0)) +
  geom_abline(slope=1, intercept=0, color='black', size=0.25) +
  #lims(x=c(0, 1), y=c(0, 1)) +
  #geom_smooth(method = "lm", fullrange = TRUE, se=FALSE) +
  facet_grid(~Group, scales='free')



# Print top genes in each category
lapply(c('C2', 'C3'), function(S) {
  Dir = paste0('reps_adding_pool/5_plots/Top_ARGs_', S, 'DD/')
  Width = 1.05+0.8
  Height = 0.5
    
  subset(ARGs_DF_scatter, Group == 'ARGs' & Sample == S & Mode == 'DD') %>% 
    arrange(HMP) %>%
    {.[1:9, 'Gene']} %>% 
    {lapply(seq_along(.), function(i) {
      G = .[[i]]
      
      ggsave(filename=paste0(Dir, 'scaled/', i, '_', G, '_', 'C3', '.pdf'), width=Width, height=Height, 
             repplotter_ofSummary(G, Sample0='C3', 
                                  DF=mutate(avg_data_ALL_adding_pool, Time = as.integer(Time)+72*as.numeric(Mode=='DD')),
                                  include_SE_bars=TRUE, rib_alpha=0.5, errorWidth=0, errorThick=0.25, hard_ys=c(0,1)))
      ggsave(filename=paste0(Dir, 'scaled/', i, '_', G, '_', 'C2', '.pdf'), width=Width, height=Height, 
             repplotter_ofSummary(G, Sample0='C2',  
                                  DF=mutate(avg_data_ALL_adding_pool, Time = as.integer(Time)+72*as.numeric(Mode=='DD')),
                                  include_SE_bars=TRUE, rib_alpha=0.5, errorWidth=0, errorThick=0.25, hard_ys=c(0,1)))
      
      ggsave(filename=paste0(Dir, 'unscaled/', i, '_', G, '_', 'C3', '.pdf'), width=Width, height=Height, 
             repplotter_ofSummary(G, Sample0='C3', 
                                  DF=mutate(avg_data_ALL_adding_pool, Time = as.integer(Time)+72*as.numeric(Mode=='DD')),
                                  include_SE_bars=TRUE, rib_alpha=0.5, errorWidth=0, errorThick=0.25))
      ggsave(filename=paste0(Dir, 'unscaled/', i, '_', G, '_', 'C2', '.pdf'), width=Width, height=Height, 
             repplotter_ofSummary(G, Sample0='C2',  
                                  DF=mutate(avg_data_ALL_adding_pool, Time = as.integer(Time)+72*as.numeric(Mode=='DD')),
                                  include_SE_bars=TRUE, rib_alpha=0.5, errorWidth=0, errorThick=0.25))
      
      ggsave(filename=paste0(Dir, 'unscaled_noY/', i, '_', G, '_', 'C3', '.pdf'), width=Width, height=Height, 
             repplotter_ofSummary(G, Sample0='C3', 
                                  DF=mutate(avg_data_ALL_adding_pool, Time = as.integer(Time)+72*as.numeric(Mode=='DD')),
                                  include_SE_bars=TRUE, rib_alpha=0.5, errorWidth=0, errorThick=0.25, NO_Y_TEXT = TRUE))
      ggsave(filename=paste0(Dir, 'unscaled_noY/', i, '_', G, '_', 'C2', '.pdf'), width=Width, height=Height, 
             repplotter_ofSummary(G, Sample0='C2',  
                                  DF=mutate(avg_data_ALL_adding_pool, Time = as.integer(Time)+72*as.numeric(Mode=='DD')),
                                  include_SE_bars=TRUE, rib_alpha=0.5, errorWidth=0, errorThick=0.25, NO_Y_TEXT = TRUE))
      
    })}
})


# Strongest ARGs
lapply(c('C2', 'C3'), function(S) {
  Dir = paste0('reps_adding_pool/5_plots/Strongest_ARGs/')
  Width = 1.05+0.8
  Height = 0.5
  
  {lapply(ARGs_strongest_cyc$Gene, function(G) {
      ggsave(filename=paste0(Dir, 'unscaled/', G, '_', S, '.pdf'), width=Width, height=Height, 
             repplotter_ofSummary(G, Sample0=S, 
                                  DF=mutate(avg_data_ALL_adding_pool, Time = as.integer(Time)+72*as.numeric(Mode=='DD')),
                                  include_SE_bars=TRUE, rib_alpha=0.5, errorWidth=0, errorThick=0.25))

      ggsave(filename=paste0(Dir, 'unscaled_noY/', G, '_', S, '.pdf'), width=Width, height=Height, 
             repplotter_ofSummary(G, Sample0=S, 
                                  DF=mutate(avg_data_ALL_adding_pool, Time = as.integer(Time)+72*as.numeric(Mode=='DD')),
                                  include_SE_bars=TRUE, rib_alpha=0.5, errorWidth=0, errorThick=0.25, NO_Y_TEXT = TRUE))
    })}
})


HSRomega <- cbind(read.table('HSRomega/HSRomega_full_CBs_first.txt', fill=TRUE), 
                  read.table('HSRomega/CB_full.txt', fill=TRUE)[1]) %>%
  `colnames<-`(c('Readtag', 'Chr', 'Start', 'Barcode'))


plot_grid(ncol=1,
HSRomega %>%
  merge(data_C2C3[[c('sample', 'mode1', 'time2', 'time1')]] %>% 
          mutate(Barcode = rownames(.))) %>% 
  subset(Start > 21296621-150 & Start < 21310706+150) %>%
  mutate(Sample_Mode = paste(sample, mode1)) %>%
  ggplot(aes(x=Start, fill=time1)) +
  geom_histogram(position='identity', binwidth=10, alpha=0.3),
HSRomega %>%
  merge(data_C2C3[[c('sample', 'mode1', 'time2', 'time1')]] %>% 
          mutate(Barcode = rownames(.))) %>% 
  subset(Start > 21296621-150 & Start < 21310706+150) %>%
  mutate(Sample_Mode = paste(sample, mode1)) %>%
  ggplot(aes(x=Start, color=time1)) +
  geom_density())




# boxplots
subset(ARGs_DF_scatter, select=c(Gene, Sample, Mode, Group, amp2, `-log10(HMP)`, HMP)) %>%
  group_by(Gene, Sample, Group) %>% 
  pivot_wider(names_from=Mode, values_from=c(amp2, `-log10(HMP)`, HMP)) %>% 
  #subset(HMP_LD < 0.05 & HMP_DD<0.05) %>% 
  mutate(Change_amp2 = amp2_LD-amp2_DD, 
            Change_nlog10HMP = `-log10(HMP)_LD` - `-log10(HMP)_DD`) %>%
  mutate(Group_Sample = paste(Group, Sample)) %>% 
  ggplot(aes(y=Change_nlog10HMP, x=Group_Sample, color=Group_Sample)) +
  geom_boxplot()

ARGs_DF_scatter%>%
  group_by(Gene, Sample, Group) %>% 
  mutate(meanNLogHMP = mean(`-log10(HMP)`), 
         meanAmp = mean(amp2)) %>% 
  ungroup %>%
  subset(meanNLogHMP > 1.3) %>%
  mutate(Group_Sample = paste(Group, Sample, Mode)) %>%
  ggplot(aes(y=`-log10(HMP)`, x=Group_Sample, color=Group_Sample)) +
  geom_boxplot()

####    ####    ####    ####    ####    ####    ####    ####    ####    ####   


# 2.8. For supplementary diagram (?)
################################################################################
repplotter_ofSummary(Gene0='sr', DF = subset(avg_data_ALL_adding_pool, Rep=='R0'), Sample0='C3', Mode0='LD')


lapply(c('R0', 'R1', 'R2'), function(R) {
  repplotter_ofSummary(Gene0='Hr38', 
                       DF = subset(avg_data_ALL_adding_pool, Rep==R), 
                       Sample0='C3', Mode0='LD', NO_Y_TEXT = TRUE, mean_point_size = 0.25, mean_linesize = 0.25)
}) %>%
  plot_grid(plotlist=., ncol=1) %>%
  ggsave(filename='oh.pdf', width=0.5, height=1.5)

reps_adding_pool$res$AGG %>% do.call(what=rbind) %>%
  subset(Gene == 'Hr38' & Sample=='C3' & Mode=='LD')


####    ####    ####    ####    ####    ####    ####    ####    ####    ####   

# 2.9. Yerbol Rp's
################################################################################
reps_adding_pool$res$AGG %>% do.call(what=rbind) %>%
  group_by(Gene, Sample, Mode) %>%
  summarise(HMP = hmp.stat(meta2d_pvalue), MinRatio = max(Ratio), meta2d_phase=Arg(mean(complex(argument = (2*pi*meta2d_phase)/24)))*24/(2*pi)) %>%
  mutate(meta2d_phase = meta2d_phase %% 24) %>%
  merge(., subset(to_filter, second_pct>0.1, select=c(Gene, Sample, Mode, second_pct, max, ampl))) %>% 
  mutate(Mean = (2*max-ampl)/2) %>%
  mutate(S_M = paste(Sample, Mode)) %>%
  
  {ggplot(., aes(Mean, -log10(HMP))) +
  scale_x_continuous(trans='log10') +
  geom_smooth(colour = 'black') +
  #geom_point() +
  #geom_point(data=subset(., Gene %in% merge(passQC, GO_ribo_cell)$Gene), 
  #           aes(Mean, -log10(HMP)), color='red', size=0.5) +
  geom_smooth(data=subset(., Gene %in% merge(passQC, GO_ribo_cell)$Gene), 
             aes(Mean, -log10(HMP)), color='red', size=0.5) +
  facet_wrap(~S_M, scales='free')}

Filtereds$DFs$reps_adding_pool %>%
  mutate(S_M = paste(Sample, Mode)) %>%
  mutate(mean_binned = round(Mean*10)/10) %>%
  group_by(S_M, mean_binned) %>%
  summarise(nlog10_HMP = mean(-log10(HMP)), MinRatio = mean(MinRatio)) %>%
  ggplot(aes(mean_binned, nlog10_HMP, color=S_M)) +
  geom_smooth()

####    ####    ####    ####    ####    ####    ####    ####    ####    ####   


# 2.9. GROUP PROMOTER ANALYSIS
################################################################################
set.seed(1)
FIMO_subsets = Subset_Stats$subsets[c('GO_ribo_mito', 'GO_ribo_cell', 'GO_ribo_noPre', 'ETC_mi', 'ARGs', 'ARGs2FC2x')] %>% 
  c(list(mt = passQC[grepl('^mt:', passQC$Gene),'Gene'],
         junk100 = sample(rownames(data_C2C3), 100),
         junk100_2 = sample(rownames(data_C2C3), 100))) %>% 
  lapply(function(G) {
    WD_split = '/Users/Teddy/Desktop/motif_databases/FLY/split_files'
    
    lapply(levels(data_C2C3$sample_mode), function(S_M) {
      S = strsplit(S_M, '_')[[1]][1]
      M = strsplit(S_M, '_')[[1]][2]
      
      
      df = subset(Filtereds$DFs$reps_adding_pool, Sample == S & Mode == M & Gene %in% G)

      if (nrow(df) <= 2) {
        return(list())
      } else {
        make_FIMO_run(df, motifPath = paste0(WD_split, '/onlyExpr_', gsub('_', '\\.', S_M), '.meme'),
                      cyclers = df$Gene, QNOTP = TRUE, thres=0.05)
      }
    }) %>% 
      `names<-`(levels(data_C2C3$sample_mode))
  })

FIMO_subsets_list$DF %>%
  {.[grepl('fru', .$TFs), ]}
  

Subset_Stats$subsets[c('GO_ribo_mito', 'GO_ribo_cell', 'GO_ribo_noPre', 'ETC_mi', 'ARGs', 'ARGs2FC2x')] %>% 
  lapply(length)

reps_adding_pool$res$AGG %>% do.call(what=rbind) %>%
  group_by(Gene, Sample, Mode) %>%
  summarise(HMP = hmp.stat(meta2d_pvalue), MinRatio = max(Ratio), meta2d_phase=Arg(mean(complex(argument = (2*pi*meta2d_phase)/24)))*24/(2*pi)) %>%
  merge(., subset(to_filter, second_pct>0.1, select=c(Gene, Sample, Mode))) %>% 
  subset(MinRatio > 1.5) %>%
  subset(Gene %in% ARG$converted$min2 & HMP<0.05) %>%
  group_by(Sample, Mode) %>%
  summarise(n())

{FIMO_subsets_list <- list(
  DF = lapply(FIMO_subsets, function(L2) {
    lapply(L2, function(L3) {
      list(TFs = paste0(L3$results$TF %>% unique %>% {.[. %in% passQC$Gene]} %>% sort, collapse='/'))
    }) %>% as.list %>%
      rbindlist(idcol='condition') %>%
      as.data.frame
  }) %>%
    as.list %>%
    rbindlist(idcol='Group')
)

FIMO_subsets_list$List <- FIMO_subsets_list$DF %>%
  {split(., list(.$Group))} %>%
  lapply(function(DF) {
    
    subset(DF, nchar(TFs) !=0)$TFs %>%
      lapply(function(TF) {
        TF %>% strsplit('/') %>% unlist %>% sort
      }) %>%
      unlist %>% {.[duplicated(.)]} %>% unique
  })


# in either day and not in either night
FIMO_subsets_list$FIMO_day_either <- unlist(FIMO_subsets_list$List$ARGs2FC2x, 
       FIMO_subsets_list$List$ETC_mi) %>% 
  {.[. %!in% unlist(
    FIMO_subsets_list$List$GO_ribo_cell, 
    FIMO_subsets_list$List$GO_ribo_mito)]} %>% unique
# in either night and not in either day
FIMO_subsets_list$FIMO_night_either <- unlist(FIMO_subsets_list$List$GO_ribo_cell, 
       FIMO_subsets_list$List$GO_ribo_mito) %>% 
{.[. %!in% unlist(FIMO_subsets_list$List$ARGs2FC2x, 
                  FIMO_subsets_list$List$ETC_mi)]} %>% unique

# in both day and not in both night
FIMO_subsets_list$FIMO_day_both <- c(FIMO_subsets_list$List$ARGs2FC2x, 
  FIMO_subsets_list$List$ETC_mi) %>% {.[duplicated(.)]} %>% unique %>%
  {.[. %!in% ( c(FIMO_subsets_list$List$GO_ribo_cell, 
                 FIMO_subsets_list$List$GO_ribo_mito) %>% 
                 {.[duplicated(.)]} %>% unique)
  ]} %>% unique
# in both night and not in both day
FIMO_subsets_list$FIMO_night_both <- c(FIMO_subsets_list$List$GO_ribo_cell, 
  FIMO_subsets_list$List$GO_ribo_mito) %>% {.[duplicated(.)]} %>% unique %>%
  {.[. %!in% ( c(FIMO_subsets_list$List$ARGs2FC2x, 
                 FIMO_subsets_list$List$ETC_mi) %>% 
                 {.[duplicated(.)]} %>% unique)
  ]} %>% unique
}

plot_grid(ncol=1, rel_heights=c(0.4, 
                                length(FIMO_subsets_list$FIMO_day_both),
                                0.4,
                                length(FIMO_subsets_list$FIMO_night_both)),
          plot_grid(labels=c('Night both, C2', 'Night both, C3'), nrow=1, ggplot()+theme_void(), ggplot()+theme_void()),
  plot_grid(nrow=1, 
            lapply(FIMO_subsets_list$FIMO_day_both, function(G) {
              repplotter_ofSummary(G, Sample0='C2', 
                                   DF=mutate(avg_data_ALL_adding_pool, Time = as.integer(Time)+72*as.numeric(Mode=='DD')),
                                   include_SE_bars=TRUE, rib_alpha=0.5, errorWidth=0, errorThick=0.25) + ggtitle(G)
            }) %>%
              plot_grid(plotlist=., ncol=1),
            lapply(FIMO_subsets_list$FIMO_day_both, function(G) {
              repplotter_ofSummary(G, Sample0='C3', 
                                   DF=mutate(avg_data_ALL_adding_pool, Time = as.integer(Time)+72*as.numeric(Mode=='DD')),
                                   include_SE_bars=TRUE, rib_alpha=0.5, errorWidth=0, errorThick=0.25) + ggtitle(G)
            }) %>%
              plot_grid(plotlist=., ncol=1)
  ),
  plot_grid(labels=c('Night both, C2', 'Night both, C3'), nrow=1, ggplot()+theme_void(), ggplot()+theme_void()),
  plot_grid(nrow=1, 
            lapply(FIMO_subsets_list$FIMO_night_both, function(G) {
              repplotter_ofSummary(G, Sample0='C2', 
                                   DF=mutate(avg_data_ALL_adding_pool, Time = as.integer(Time)+72*as.numeric(Mode=='DD')),
                                   include_SE_bars=TRUE, rib_alpha=0.5, errorWidth=0, errorThick=0.25) + ggtitle(G)
            }) %>%
              plot_grid(plotlist=., ncol=1),
            lapply(FIMO_subsets_list$FIMO_night_both, function(G) {
              repplotter_ofSummary(G, Sample0='C3', 
                                   DF=mutate(avg_data_ALL_adding_pool, Time = as.integer(Time)+72*as.numeric(Mode=='DD')),
                                   include_SE_bars=TRUE, rib_alpha=0.5, errorWidth=0, errorThick=0.25) + ggtitle(G)
            }) %>%
              plot_grid(plotlist=., ncol=1)
  )
) %>% ggsave(filename='reps_adding_pool/Tables/FIMO_exclusive.pdf', width=4, height=10)

####    ####    ####    ####    ####    ####    ####    ####    ####    ####   


# 2.10. CWO ChIP, FRU
################################################################################
WD_ChIP <- "/Users/Teddy/Desktop/CCG_distribution/cwo_ChIP.xlsx"

ChIP_Hardin <- list(CWO = rbind(readxl::read_excel(WD_ChIP, 
                                                   sheet = 1, col_names = TRUE) %>% cbind(ZT = 'ZT2'),
                                readxl::read_excel(WD_ChIP, 
                                                   sheet = 2, col_names = TRUE) %>% cbind(ZT = 'ZT14')),
                    CLK = rbind(readxl::read_excel(WD_ChIP, 
                                                   sheet = 3, col_names = TRUE) %>% cbind(ZT = 'ZT2'), 
                                readxl::read_excel(WD_ChIP, 
                                                   sheet = 4, col_names = TRUE) %>% cbind(ZT = 'ZT14'))) %>% # NOTE, RESULTS ARE GENE-MAPPED, AFTER REMOVAL OF INTERGENIC MAPPINGS; adult heads
  {c(., list(CWO_common_genes = .$CWO %>% 
               group_by(`Gene Name`) %>% filter(n()>1) %>% .[['Gene Name']] %>% unique, 
             CLK_common_genes = .$CLK %>%
               group_by(`Gene Name`) %>% filter(n()>1) %>% .[['Gene Name']] %>% unique,
             CWO_either_genes = .$CWO[['Gene Name']] %>% unique,
             CLK_either_genes = .$CLK[['Gene Name']] %>% unique, 
             both_either = intersect(.$CWO[['Gene Name']] %>% unique, .$CLK[['Gene Name']] %>% unique)))}


 %>%
  {.$`Gene Name`}

lapply(levels(data_C2C3$sample_mode), function(S_M) {
  S = gsub('_.*', '', S_M)
  M = gsub('.*_', '', S_M)
  
  U0 = subset(to_filter, second_pct > 0.1 & Sample == S & Mode == M)$Gene
  B0 = ChIP_Hardin$CWO$`Gene Name`
  B0 = B0[B0 %in% U0]
  
  A0 = subset(Filtereds$DFs$reps_adding_pool, Sample == S & Mode == M)$Gene
  
  run_fisher_general(A=A0, B=B0, U=U0)
  
}) %>% `names<-`(levels(data_C2C3$sample_mode))




# WHAT DIFFERENTIATES LD?

TFs_LD <- subset(Filtereds$DFs$reps_adding_pool, Mode == 'LD' & Gene %in% merge(passQC, GO_TRs)$Gene)$Gene

TFs_DD <- subset(Filtereds$DFs$reps_adding_pool, Mode == 'DD' & Gene %in% merge(passQC, GO_TRs)$Gene)$Gene

TFs_LD_inBoth <- TFs_LD[duplicated(TFs_LD)]
TFs_DD_inBoth <- TFs_LD[duplicated(TFs_DD)] # none
TFs_LD_inBoth_notDDEither <- TFs_LD_inBoth[TFs_LD_inBoth %!in% TFs_DD]

TFs_LD_inBoth_notDDEither[ TFs_LD_inBoth_notDDEither %in% ARG$converted$min2]


# FRU:
plot_grid(nrow=4,
          repplotter_ofSummary('fru', Sample0='C2', 
                               DF=mutate(avg_data_ALL_adding_pool, Time = as.integer(Time)+72*as.numeric(Mode=='DD')),
                               include_SE_bars=TRUE, rib_alpha=0.5, errorWidth=0, errorThick=0.25, NO_Y_TEXT = FALSE) + ggtitle('fru C2'),
          
          repplotter_ofSummary('fru', Sample0='C3', 
                               DF=mutate(avg_data_ALL_adding_pool, Time = as.integer(Time)+72*as.numeric(Mode=='DD')),
                               include_SE_bars=TRUE, rib_alpha=0.5, errorWidth=0, errorThick=0.25, NO_Y_TEXT = FALSE) + ggtitle('fru C3'), 
          
          repplotter_ofSummary('Rh7', Sample0='C2', 
                               DF=mutate(avg_data_ALL_adding_pool, Time = as.integer(Time)+72*as.numeric(Mode=='DD')),
                               include_SE_bars=TRUE, rib_alpha=0.5, errorWidth=0, errorThick=0.25, NO_Y_TEXT = FALSE) + ggtitle('Rh7 C2'),
          
          repplotter_ofSummary('Rh7', Sample0='C3', 
                               DF=mutate(avg_data_ALL_adding_pool, Time = as.integer(Time)+72*as.numeric(Mode=='DD')),
                               include_SE_bars=TRUE, rib_alpha=0.5, errorWidth=0, errorThick=0.25, NO_Y_TEXT = FALSE) + ggtitle('Rh7 C3'), 
          
          repplotter_ofSummary('cwo', Sample0='C2', 
                               DF=mutate(avg_data_ALL_adding_pool, Time = as.integer(Time)+72*as.numeric(Mode=='DD')),
                               include_SE_bars=TRUE, rib_alpha=0.5, errorWidth=0, errorThick=0.25, NO_Y_TEXT = FALSE) + ggtitle('cwo C2'),
          
          repplotter_ofSummary('cwo', Sample0='C3', 
                               DF=mutate(avg_data_ALL_adding_pool, Time = as.integer(Time)+72*as.numeric(Mode=='DD')),
                               include_SE_bars=TRUE, rib_alpha=0.5, errorWidth=0, errorThick=0.25, NO_Y_TEXT = FALSE) + ggtitle('cwo C3'),
          
          repplotter_ofSummary('usp', Sample0='C2', 
                               DF=mutate(avg_data_ALL_adding_pool, Time = as.integer(Time)+72*as.numeric(Mode=='DD')),
                               include_SE_bars=TRUE, rib_alpha=0.5, errorWidth=0, errorThick=0.25, NO_Y_TEXT = FALSE) + ggtitle('usp C2'),
          repplotter_ofSummary('usp', Sample0='C3', 
                               DF=mutate(avg_data_ALL_adding_pool, Time = as.integer(Time)+72*as.numeric(Mode=='DD')),
                               include_SE_bars=TRUE, rib_alpha=0.5, errorWidth=0, errorThick=0.25, NO_Y_TEXT = FALSE) + ggtitle('usp C3')
)

subset(Filtereds$DFs$reps_adding_pool, Gene == 'Rh7')$meta2d_phase # 9

Filtereds$DFs$reps_adding_pool %>%
  mutate(phaseBin = round(meta2d_phase*2)/2) %>%
  group_by(Sample, Mode, phaseBin) %>% 
  summarise(binCount = n()) %>% 
  mutate(normBinCount = binCount/sum(binCount), binCount=NULL) %>%
  pivot_wider(names_from = Mode, values_from = normBinCount) %>% 
  mutate(ratioBinCount = LD/DD) %>%
  ggplot(aes(phaseBin, ratioBinCount, fill=Sample)) +
  geom_bar(stat = 'identity')

####    ####    ####    ####    ####    ####    ####    ####    ####    ####   


# 2.11. Printing supplementary tables:
################################################################################

reps_adding_pool$res$AGG %>% do.call(what=rbind) %>%
  group_by(Gene, Sample, Mode) %>%
  summarise(HMP = hmp.stat(meta2d_pvalue), MinRatio = max(Ratio), 
            meta2d_phase=Arg(mean(complex(argument = (2*pi*meta2d_phase)/24)))*24/(2*pi)) %>%
  merge(subset(to_filter, select=c(Gene, Sample, Mode, Cell_proportion, second_pct))) %>% 
  mutate(DETECTED = second_pct>0.1) %>%
  left_join(passQC) %>%
  {left_join(., BH_procedure(subset(., DETECTED==TRUE), 
                             Qcol=NULL, Pcol='HMP', Q=0.05, CUT = FALSE))} %>%
  mutate(HMP = ifelse(DETECTED, HMP, NA), 
         meta2d_phase = ifelse(DETECTED, meta2d_phase, NA), 
         MinRatio = ifelse(DETECTED, MinRatio, NA)) %>%
  `colnames<-`(c('Gene', 'Cell_Type', 'Mode', 'Cycling_HMP_value', 'Ratio', 'Cycling_phase', 
                 'Avg_cell_percent_expressing', 'second_most_cell_percent_expressing', 
                 'DETECTED', 'ENSEMBL_ID', 'Cycling_BH.Q_value', 'PASS_CYCLING')) %>%
    {.[, c(1,2,3, 7,8,9,10, 5,4,6, 11,12)]} %>%
  write.csv(file='reps_adding_pool/Tables/6_GeneStats.csv', quote=FALSE, row.names=FALSE)


list.files('reps_adding_pool/2_GO/Results') %>%
  lapply(function(N) {
    SM = gsub('_.*', '', N)
    GO_type = gsub('.txt', '', N) %>%
      gsub('.*_', '', .)
    read_delim('reps_adding_pool/2_GO/Results/C3DD_CC.txt', skip = 11) %>%
      {cbind(Condition = SM, 
             Category = GO_type, 
             .)}
  }) %>%
  do.call(what='rbind') %>% 
  write.csv(file='reps_adding_pool/Tables/7_GO.csv', quote=FALSE)


Subset_Stats$subsets %>%
  {.[c('ARGs', 'ARGs2FC2x', 'GO_ribo_cell', 'GO_ribo_mito', 'ETC_mi')]} %>%
  stack %>%
  `colnames<-`(c('Gene', 'Group')) %>%
  merge(data.frame(Group = c('ARGs', 'ARGs2FC2x', 'GO_ribo_cell', 'GO_ribo_mito', 'ETC_mi'),
                   Origin = c('Chen_2016', 'Chen_2016', 'GOCC:0022626', 'GOCC:0005761', 'GOBP:0022900 ^mt:'))) %>%
  merge(subset(to_filter, select=c('Gene', 'Sample', 'Mode', 'second_pct'))) %>%
  mutate(DETECTED = second_pct > 0.1, .keep='unused') %>%
  mutate(CYCLE = (paste0(Gene, Sample, Mode) %in% paste0(Filtereds$DFs$reps_adding_pool$Gene,
                                                         Filtereds$DFs$reps_adding_pool$Sample,
                                                         Filtereds$DFs$reps_adding_pool$Mode) )) %>%
  write.csv(file='reps_adding_pool/Tables/8_Groups.csv', quote=FALSE, row.names=FALSE)


lapply(c('R0', 'R1', 'R2'), function(R) {
  cbind(reps_adding_pool[[R]]$avg[[c('sample', 'mode1', 'time3')]], Rep=R)
}) %>%
  do.call(what=rbind) %>% 
  write.csv(file='reps_adding_pool/Tables/9_pRep_assignment.csv', quote=FALSE)


FIMO_subsets_list$DF %>% 
  write.csv('reps_adding_pool/Tables/X_FIMO_subsets.csv', row.names=FALSE, quote=FALSE)


FIMO_subsets_list[lapply(FIMO_subsets_list, class) == 'character'] %>%
  lapply(function(L) {
    list(TFs = paste0(L, collapse='/'))
  }) %>%
  do.call(what=rbind) %>%
  write.csv('reps_adding_pool/Tables/XX_FIMO_exclusive.csv', quote=FALSE)



####    ####    ####    ####    ####    ####    ####    ####    ####    ####   


# Seurat: AGG and AGG norms

# functions RData
# to_filter
# pool3_AGG
# Filtereds
# ALL avg_data...




save.image('AGG_to_analyze/functions.RData')

# %!in% 
# add_cell_subsets 
# add_rectangles_to_plots 
# avg_plotter_concat 
# avg_plotter_concat_mode 
# BH_procedure
# counter
# daynorm
# from_corrplot_meanmaker
# gene_hclusterer
# hist_phase
# hmap_phase
# make_corrplot
# make_corrs
# make_heatmap
# make_panel
# make_pdf_filtered
# make_pdf_targets
# make_phase_from_GOdf
# make_plot
# oneRect
# pattern_to_gene_helper
# pdf_combine_and_delete
# pdf_combine_and_delete_N
# pdf_gene_lists
# plot_annotation_table
# plotter_geneCombine
# plotter_geneCombine_modes
# run_fisher_general
# split_list
# standard_plotter
# standard_plotter_2
# subplotter_concat
# subset_one_cluster_timepoint
# substatter
# take_percentiles
# trim_TFclusters

c('to_filter',
  'avg_data', 
  'avg_data_down_times2_AGG',
  'avg_data_time2',
  'avg_data_times2',
  'avg_data_times2_01',
  'avg_data_times2_01_2d',
  'avg_data_times2_AGG',
  'avg_data_times2_max1',
  
  'added_C2C3', 
  'added_C2C3_AGG',
  #'added_C2C3_pool',
  #'added_C2C3_pool3',
  #'added_C2C3_pool_AGG',
  
  'ARG') %>%
  lapply(function(Name) {saveRDS(get(Name), paste0('AGG_to_analyze/', Name, '.rds'))})

####    ####    ####    ####    ####    ####    ####    ####    ####    ####   



# Ambient

FeaturePlot(data_C2C3, features=c('Gs2', 'e', 'Eaat1', 'per', 'cwo', 'repo', 'ninaE', 'Rh2', 'Rh4', 'Rh5', 'Rh6'))
o <- data_C2C3@assays$RNA@counts[c('per', 'cwo'),] %>% t %>% as.data.frame

run_fisher_general(A = rownames(subset(o, per != 0)), 
                   B = rownames(subset(o, cwo != 0)), 
                   U = rownames(o))


Ines_bulk <- read.table('/Users/Teddy/Desktop/to_GEO/Bulk/GSE233184_PolyAnormtomax_6tp.txt', header=TRUE)





