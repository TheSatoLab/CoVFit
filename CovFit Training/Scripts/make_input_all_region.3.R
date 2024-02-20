#!/usr/bin/env R

library(tidyverse)
library(data.table)
library(rbin)
library(cmdstanr)
library(patchwork)
library(RColorBrewer)
library(nnet)
library(broom)
library(seqinr)
library(Biostrings)
library(tictoc)
library(MASS)
library(foreach)
library(doParallel)
library(dplyr)


#Transmissibility
bin.size <- 1
generation_time <- 2.1
core.num <- 12
date.max <- "2023-10-15"
ref_hap_clade <- "22B"

setwd("/Users/jumpeiito/Desktop/analysis/Sato_analysis/fitness_prediction/dataset/231102")

UK.v <- c("England","Scotland","Wales")

metadata.original.name <- "metadata.tsv"
metadata.original <- fread(metadata.original.name,header=T,sep="\t",quote="",check.names=T)

metadata.original <- metadata.original %>% filter(Host == "Human")

metadata.original <- metadata.original %>% mutate(Virus.name2 = Virus.name)  %>% separate(Virus.name2,c("virus", "country", "name","year"), sep="/")
  
metadata.name <- "20231101.full.nextclade.tsv"
metadata <- fread(metadata.name,header=T,sep="\t",quote="",check.names=T)

data.count_N.name <- "gene_S.translation.count_table.txt"
data.count_N <- fread(data.count_N.name,header=T,sep="\t",quote="",check.names=T)

data.count_N.filtered <- data.count_N %>% dplyr::rename(seqName = name) %>% mutate(prop.N = count_N / length) %>% filter(prop.N<0.01,count_stop==1)

metadata.merged <- metadata %>% inner_join(data.count_N.filtered,by="seqName")
metadata.merged <- metadata.merged %>% distinct(seqName,.keep_all=T)
metadata.merged <- metadata.merged %>% mutate(seqName2 = seqName) %>% separate(seqName2,c("virus", "country", "name","time"), sep="/")
metadata.merged <- metadata.merged %>% mutate(time2 = time) %>% separate(time2,c("year","Collection.date","Submission.date"), sep="\\|")
metadata.merged <- metadata.merged %>% filter(!is.na(Collection.date),str_length(Collection.date) == 10, virus == "hCoV-19")
metadata.merged <- metadata.merged %>% mutate(Collection.date = as.Date(Collection.date)) %>% filter(Collection.date <= as.Date(date.max))
metadata.merged <- metadata.merged %>% inner_join(metadata.original %>% dplyr::select(name),by="name")

metadata.merged <- metadata.merged %>% mutate(country = ifelse(country %in% UK.v,"UK",as.character(country)))

country.interest.df <- metadata.merged %>% group_by(country) %>% summarize(count.country = n()) %>% arrange(desc(count.country)) %>% filter(count.country > 100000)
country.interest.v <- country.interest.df$country

metadata.merged <- metadata.merged %>% filter(country %in% country.interest.df$country)

mut.info.name <- "20231101.full.nextclade.mut_long.tsv"
mut.info <- fread(mut.info.name,header=T,sep="\t",quote="",check.names=T)

mut.info <- mut.info %>% filter(grepl("S:",mut))

#mut.info.old <- mut.info %>% inner_join(metadata.merged %>% dplyr::select(seqName),by="seqName")
mut.info <- mut.info %>% right_join(metadata.merged %>% dplyr::select(seqName),by="seqName")

count.df.mut <- mut.info %>% group_by(mut) %>% summarize(count.mut = n())
count.df.mut <- count.df.mut %>% arrange(desc(count.mut)) %>% filter(count.mut > 100)

mut.info <- mut.info %>% inner_join(count.df.mut %>% dplyr::select(mut),by="mut")
mut.info <- mut.info %>% unique()
mut.info <- mut.info %>% mutate(mut = as.factor(mut))

mut.info.sum <- mut.info %>% group_by(seqName) %>% summarize(haplotype = toString(mut))
mut.info.sum <- mut.info.sum %>% mutate(hap_Id = paste("hap_",as.character(as.numeric(as.factor(haplotype))),sep=""))
mut.info.sum <- mut.info.sum %>% inner_join(metadata.merged %>% dplyr::select(seqName,country),by="seqName")

hap.info.df <- mut.info.sum %>% ungroup() %>% dplyr::select(hap_Id, haplotype) %>% unique()

count.df.hap <- mut.info.sum %>% group_by(hap_Id,country) %>% summarize(count.hap = n()) %>% arrange(desc(count.hap))

count.df.hap <- count.df.hap %>% filter(count.hap>=20)

mut.info.sum <- mut.info.sum %>% filter(hap_Id %in% as.character(count.df.hap$hap_Id))

metadata.merged <- metadata.merged %>% inner_join(mut.info.sum %>% dplyr::select(seqName,hap_Id),by="seqName")
metadata.merged <- metadata.merged %>% distinct(seqName,.keep_all = T)

metadata.merged <- metadata.merged %>% filter(!is.na(Collection.date))

metadata.filtered.interest <- metadata.merged %>%
  mutate(Collection.date = as.Date(Collection.date),
         date.num = as.numeric(Collection.date) - min(as.numeric(Collection.date))  + 1,
         date.bin = cut(date.num,seq(0,max(date.num),bin.size)),
         date.bin.num = as.numeric(date.bin))

metadata.filtered.interest <- metadata.filtered.interest %>% filter(!is.na(date.bin.num))

count.hap.mat <- metadata.filtered.interest %>%
  filter(clade == ref_hap_clade) %>%
  group_by(country,hap_Id) %>%
  summarize(count.hap = n()) %>%
  spread(key=country,value=count.hap) %>% na.omit()

count.df.hap2 <- count.hap.mat %>% gather(key=country,value=count,-hap_Id)
count.df.hap2.min <- count.df.hap2 %>% group_by(hap_Id) %>% summarize(count.min = min(count)) %>% arrange(desc(count.min))

hap_Id.ref <- count.df.hap2.min$hap_Id[1] 

print(count.hap.mat %>% filter(hap_Id == hap_Id.ref))


rm(metadata.original);gc();gc()
rm(metadata);gc();gc()
rm(mut.info);gc();gc()



# Register the parallel backend to use many processors

registerDoParallel(cores=core.num)

data.coef.df <- data.frame()

# Parallelize the for loop using foreach
results <- foreach(country.interest = country.interest.df$country, .combine = rbind) %dopar% {
  print(country.interest)

  metadata.filtered.interest.country <- metadata.filtered.interest %>% filter(country == country.interest)
  
  metadata.filtered.interest.bin <- metadata.filtered.interest.country %>% group_by(date.bin.num,hap_Id) %>% summarize(count = n()) %>% ungroup()
  
  hap.count.df <- metadata.filtered.interest.bin %>% group_by(hap_Id) %>% summarize(total = sum(count)) %>% filter(total >= 20) %>% arrange(desc(total))
  
  metadata.filtered.interest.bin.spread <- metadata.filtered.interest.bin %>% filter(hap_Id %in% as.character(hap.count.df$hap_Id)) %>% spread(key=hap_Id,value = count)
  metadata.filtered.interest.bin.spread[is.na(metadata.filtered.interest.bin.spread)] <- 0
  
  X <- as.matrix(data.frame(X0 = 1, X1 = metadata.filtered.interest.bin.spread$date.bin.num))
  
  Y <- metadata.filtered.interest.bin.spread %>% dplyr::select(- date.bin.num)
  
  hap_Id.order.v <- c(hap_Id.ref,hap.count.df$hap_Id[-which(hap.count.df$hap_Id == hap_Id.ref)])
  
  Y <- Y %>% dplyr::select(all_of(hap_Id.order.v))
  
  max.date <- max(metadata.filtered.interest$date.bin.num)
  
  date.v <- X[,2] / max.date
  Y <- as.matrix(Y)
  
  fit.multinom <- multinom(Y ~ date.v, maxit = 10000, MaxNWts = 25000, model =T)
  
  
  fit.coef <- coef(fit.multinom)
  colnames(fit.coef) <- c("intercept","slope")
  
  fit.coef <- fit.coef %>% as.data.frame() %>% mutate(hap_Id = rownames(fit.coef))
  fit.coef <- fit.coef %>% mutate(relative_Re = exp(((slope / max.date) / bin.size)  * generation_time))
  fit.coef <- fit.coef %>% dplyr::select(hap_Id,relative_Re)
  
  hap_Id.reference <- hap_Id.order.v[1]
  
  fit.coef <- rbind(data.frame(hap_Id = hap_Id.reference, relative_Re = 1),fit.coef)
  fit.coef <- fit.coef %>% mutate(country = country.interest)
  
  return(fit.coef)
}

# Combine all results into a single dataframe
data.coef.df <- do.call(rbind, results)

# Stop the parallel cluster after the job is done
stopImplicitCluster()


g <- ggplot(results,aes(x=country,y=relative_Re))
g <- g + geom_violin()
g


data.coef.df <- results

#mut.info.name <- #was "20230330.nextclade.ProtBert.mut_long.tsv", but this is old and not found in the directory. 
mut.info <- fread(mut.info.name,header=T,sep="\t",quote="",check.names=T)

mut.info <- mut.info %>% filter(grepl("S:",mut))

mut.info <- mut.info %>% right_join(metadata.filtered.interest %>% dplyr::select(seqName,hap_Id), by = "seqName")

mut.info <- mut.info %>% unique()
mut.info <- mut.info %>% mutate(mut = as.factor(mut))

mut.info.sum.all_S <- mut.info %>% group_by(seqName) %>% summarize(haplotype.all = toString(mut))

mut.info.sum.all_S <- mut.info.sum.all_S %>% mutate(hap_Id.all_mut = paste("hap_all_",as.character(as.numeric(as.factor(haplotype.all))),sep=""))

metadata.filtered.interest <- metadata.filtered.interest %>% inner_join(mut.info.sum.all_S %>% dplyr::select(seqName,hap_Id.all_mut),by="seqName")

date.quantile.df <- metadata.filtered.interest %>% mutate(Collection.date = as.Date(Collection.date)) %>% group_by(hap_Id) %>% summarize(date.first = quantile(Collection.date,0.01,type=1))
date.quantile.df <- date.quantile.df %>% arrange(date.first)

all_S.major.hap.df <- metadata.filtered.interest %>% group_by(hap_Id,hap_Id.all_mut) %>% summarize(count = n()) %>% ungroup() %>% group_by(hap_Id) %>% slice_max(count,n=1,with_ties = F)

metadata.filtered.interest.representative <- metadata.filtered.interest %>% inner_join(all_S.major.hap.df %>% dplyr::select(hap_Id,hap_Id.all_mut),by=c("hap_Id","hap_Id.all_mut")) %>% group_by(hap_Id) %>% sample_n(1)

metadata.filtered.interest.representative <- metadata.filtered.interest.representative %>% dplyr::select(seqName,clade,Nextclade_pango,hap_Id)
metadata.filtered.interest.representative <- metadata.filtered.interest.representative %>% inner_join(hap.info.df,by="hap_Id")

data.coef.df.merged <- data.coef.df %>% inner_join(metadata.filtered.interest.representative,by="hap_Id")
data.coef.df.merged <- data.coef.df.merged %>% arrange(desc(relative_Re))

data.coef.df.merged <- data.coef.df.merged %>% left_join(date.quantile.df,by="hap_Id")


out.name <- "metadata.representative.all_countries.with_date.v2.txt"
write.table(data.coef.df.merged,out.name,col.names=T,row.names=F,sep="\t",quote=F)



count.df.hap.all_regions <- metadata.filtered.interest %>% group_by(hap_Id) %>% summarize(count.total = n())

data.coef.df.merged <- data.coef.df.merged %>% inner_join(count.df.hap,by=c("hap_Id","country"))


g <- ggplot(data.coef.df.merged, aes(x=date.first,y=log(relative_Re), size = count.hap, fill = clade))
g <- g + geom_point(shape = 21, alpha = 0.7)
g <- g + theme(panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.background = element_blank(),
               strip.text = element_text(size=8)
)
g <- g + theme_set(theme_classic(base_size = 12, base_family = "Helvetica"))
g <- g + theme(panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.background = element_blank(),
               strip.text = element_text(size=8)
)
g <- g + geom_hline(yintercept=c(-1,-0.5,0,0.25,0.5,1),color="gray70")
g <- g + facet_wrap(~country)
g <- g + xlab("Date") + ylab("Relative Re (log)")
#g <- g + scale_y_continuous(lim=c(0,2))

pdf.name <- "time_vs_fitness_each_country_log.pdf"
pdf(pdf.name, width = 15, height = 10)
print(g)
dev.off()

g <- ggplot(data.coef.df.merged, aes(x=date.first,y=log(relative_Re+0.1), size = count.hap, fill = clade))
g <- g + geom_point(shape = 21, alpha = 0.7)
g <- g + theme(panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.background = element_blank(),
               strip.text = element_text(size=8)
)
g <- g + theme_set(theme_classic(base_size = 12, base_family = "Helvetica"))
g <- g + theme(panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.background = element_blank(),
               strip.text = element_text(size=8)
)
g <- g + geom_hline(yintercept=c(-1,-0.5,0,0.25,0.5,1),color="gray70")
g <- g + facet_wrap(~country)
g <- g + xlab("Date") + ylab("Relative Re (log)")
#g <- g + scale_y_continuous(lim=c(0,2))

pdf.name <- "time_vs_fitness_each_country_log_with_peudo.pdf"
pdf(pdf.name, width = 15, height = 10)
print(g)
dev.off()



g <- ggplot(data.coef.df.merged, aes(x=date.first,y=relative_Re, size = count.hap, fill = clade))
g <- g + geom_point(shape = 21, alpha = 0.7)
g <- g + theme(panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.background = element_blank(),
               strip.text = element_text(size=8)
)
g <- g + theme_set(theme_classic(base_size = 12, base_family = "Helvetica"))
g <- g + theme(panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.background = element_blank(),
               strip.text = element_text(size=8)
)
g <- g + geom_hline(yintercept=c(1,2,3,4),color="gray70")
g <- g + facet_wrap(~country)
g <- g + xlab("Date") + ylab("Relative Re")
#g <- g + scale_y_continuous(lim=c(0,2.5))

pdf.name <- "time_vs_fitness_each_country.pdf"
pdf(pdf.name, width = 15, height = 10)
print(g)
dev.off()


haplotype.mut_profile <- metadata.filtered.interest.representative %>% dplyr::select(hap_Id,haplotype) %>% separate_rows(haplotype, sep = ",\\s*")
haplotype.mut_profile <- haplotype.mut_profile %>% mutate(pos = gsub("[^0-9]","",haplotype), pos = as.numeric(pos))
haplotype.mut_profile <- haplotype.mut_profile %>% inner_join(metadata.filtered.interest.representative %>% dplyr::select(hap_Id,clade),by="hap_Id")
haplotype.mut_profile <- haplotype.mut_profile %>% mutate(hap_Id = factor(hap_Id, levels=date.quantile.df$hap_Id)) %>% arrange(hap_Id)

count.hap_with_mut_each_pos <- haplotype.mut_profile %>% group_by(pos) %>% summarize(count.hap = n()) %>% arrange(pos)

g2 <- ggplot(count.hap_with_mut_each_pos,aes(x=pos,y=log(count.hap,10)))
g2 <- g2 + geom_point(size=0.01)
g2 <- g2 + xlab("") + ylab("log10(count)")

g1 <- ggplot(haplotype.mut_profile,aes(x=pos,y=hap_Id, color = clade))
g1 <- g1 + geom_point(size=0.01)

g1 <- g1 + theme(panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.background = element_blank(),
               strip.text = element_text(size=8)
)
g1 <- g1 + theme_set(theme_classic(base_size = 12, base_family = "Helvetica"))
g1 <- g1 + theme(panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.background = element_blank(),
               strip.text = element_text(size=8)
)
g1 <- g1 + theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())
g1 <- g1 + xlab("Amino acid position") + ylab(paste("haplotype (",as.character(nrow(date.quantile.df)),")",sep=""))

pdf.name <- "mut_heatmap.pdf"
pdf(pdf.name, width = 10, height = 5)
print(g2 / g1 + plot_layout(ncol = 1, heights = c(1, 3)))
dev.off()

