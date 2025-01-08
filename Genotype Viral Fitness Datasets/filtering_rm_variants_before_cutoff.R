#!/usr/bin/env R

library(tidyverse)
library(data.table)
library(patchwork)
library(viridis)
library(scales)

setwd("your directory")

metadata.name <- "dataset/metadata.representative.all_countries.with_date.v2.with_seq_231102.txt"
metadata <- read.table(metadata.name,header=T,sep="\t")
metadata <- metadata %>% mutate(date.first = as.Date(date.first))
metadata <- metadata %>% mutate(clade = ifelse(str_detect(Nextclade_pango,"BA\\.2\\.86") | str_detect(Nextclade_pango,"JN\\."), "BA.2.86",as.character(clade)))


clade_cutoff_df.name <- "dataset/clade_cutoff_for_cutoff_sliding.txt"
clade_cutoff_df <- read.table(clade_cutoff_df.name,header=T,sep="\t")

dim(metadata)
metadata <- metadata %>% left_join(clade_cutoff_df, by = "clade")
dim(metadata)

metadata.filtered <- metadata %>% filter(date.first >= cutoff | is.na(cutoff))

out.name <- "dataset/metadata.representative.all_countries.with_date.v2.with_seq_231102_wo_variants_before_cutof.txt"
write.table(metadata.filtered, out.name, col.names=T,row.names=F, sep = "\t", quote = F)



metadata.filtered.interest <- metadata.filtered %>% select(hap_Id,clade,date.first) %>% unique()

clade_first_date.df <- metadata.filtered.interest %>% group_by(clade) %>% summarize(date.first = min(date.first), month = as.Date(format(date.first, "%Y-%m-01")))

month.v <- c("2022-01-01","2022-02-01","2022-03-01","2022-04-01","2022-05-01","2022-06-01","2022-07-01","2022-08-01","2022-09-01","2022-10-01","2022-11-01","2022-12-01",
             "2023-01-01","2023-02-01","2023-03-01","2023-04-01","2023-05-01","2023-06-01","2023-07-01","2023-08-01","2023-09-01")


clade.v <- rev(c("21K","21L","22C","22A","22B","22D","22E","22F","23C","23A","23F","BA.2.86"))
clade_name.v <- rev(c("21K (BA.1)","21L (BA.2)","22C (BA.2.12)","22A (BA.4)","22B (BA.5)","22D (BA.2.75)","22E (BQ.1)","22F (XBB)","23C (CH.1.1)","23A (XBB.1.5)","23F (EG.5.1)","BA.2.86"))

clade_name.df <- data.frame(clade = clade.v, clade_name = clade_name.v)

presence_df <- data.frame()

for(clade_interest in clade.v){
  for(month_interest in month.v){
    #clade_interest <- clade.v[1]
    #month <- month.v[1]
    cutoff <- clade_first_date.df %>% filter(clade == clade_interest) %>% pull(month)
    value <- 0
    if(as.Date(month_interest) > as.Date(cutoff)){
      value <- 1
    }
    temp_df <- data.frame(clade = clade_interest, month = month_interest, value = value)
    presence_df <- rbind(presence_df,temp_df)
  }
}

presence_df <- presence_df %>% left_join(clade_name.df, by = "clade")
presence_df <- presence_df %>% mutate(clade_name = factor(clade_name,levels=clade_name.v))

presence_df$month <- format(as.Date(presence_df$month, format = "%Y-%m-%d"), "%Y-%m")

g <- ggplot(presence_df, aes(x=month, y=clade_name, fill=as.factor(value))) +
  geom_tile() +
  theme_classic() +
  scale_fill_manual(values = c("0" = "gray80", "1" = "navy")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

g

pdf.name <- "dataset/clade_cutoff_for_cutoff_sliding.pdf"
pdf(pdf.name, width =5.5, height = 2.7)
plot(g)
dev.off()

#


