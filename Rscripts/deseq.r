library("tidyverse")
library ("phyloseq")
library("ggpubr")
library("DESeq2")
library("qiime2R")

set.seed(711)

dir = as.character(read.table("dir.txt"))
setwd(dir)

load("./16s/physeq_16s.RData") # Load 16s data
load("./ITS/physeq_ITS.RData") # Load ITS data

# Obtain ASVs (16s)
seq_16s = read_qza("./16s/dada2/FeatureData[Sequence]_16s.qza")
seq_16s = as.data.frame(seq_16s$data)
seq_16s$ASV = row.names(seq_16s)
colnames(seq_16s) = c("Sequence", "ASV")
seq_16s = as_tibble(seq_16s)
#

# Obtain ASVs (ITS)
seq_ITS = read_qza("./ITS/dada2/FeatureData[Sequence]_ITS.qza")
seq_ITS = as.data.frame(seq_ITS$data)
seq_ITS$ASV = row.names(seq_ITS)
colnames(seq_ITS) = c("Sequence", "ASV")
seq_ITS = as_tibble(seq_ITS)
#

plasmopara = merge_samples(physeq_ITS, "Sample", fun = mean)
plasmopara = transform_sample_counts(plasmopara, function(x) x / sum(x) * 100) # Transform abundance data to relative abundances

# Obtain relative abundances of Plasmopara viticola in samples
plasmopara = psmelt(plasmopara) %>%
  group_by(., Genus, Sample) %>%
  summarise(., Abundance = sum(Abundance))
plasmopara = filter(plasmopara, Genus %in% grep("Plasmopara", Genus, value = TRUE))
plasmopara$Genus = "Plasmopara viticola"
plasmopara = plasmopara %>%
	group_by(., Genus, Sample) %>%
	summarise(., Abundance = sum(Abundance)) %>%
    ungroup()
#

# Obtain Plasmopara viticola presence factor
plasmopara = plasmopara %>%
    mutate(Plasmopara = if_else(Abundance>0, "Yes", "No")) %>%
    select(Sample, Plasmopara) %>%
    data.frame()
#

row.names(plasmopara) = plasmopara$Sample

physeq_16s = merge_samples(physeq_16s, "Sample", fun = mean)
sample_data(physeq_16s) = sample_data(plasmopara)


# DESeq2 process (16s)
deseq2_16s = phyloseq_to_deseq2(physeq_16s, ~ Plasmopara)
deseq2_16s = DESeq(deseq2_16s, test="Wald", fitType="parametric")
#

# Filter significant results of DESeq2 (16s)
res = results(deseq2_16s, cooksCutoff = FALSE)
res = as.data.frame(res)
de_deseq_significant = filter(res, padj < 0.01)
de_deseq_significant$ASV = row.names(de_deseq_significant)

significant_tax_16s = as.data.frame(tax_table(physeq_16s)[de_deseq_significant$ASV,])
significant_tax_16s$ASV = row.names(significant_tax_16s)
de_deseq_significant = left_join(de_deseq_significant, significant_tax_16s)

de_deseq_significant$Genus = factor(de_deseq_significant$Genus, 
    unique(de_deseq_significant$Genus[order(de_deseq_significant$log2FoldChange, decreasing = TRUE)]))

seq_16s = filter(seq_16s, ASV %in% significant_tax_16s$ASV)
significant_seq_16s = left_join(seq_16s, de_deseq_significant)
write.table(significant_seq_16s, file="Supporting Information Table S4.txt")
#

de_deseq_significant=mutate(de_deseq_significant, Plasmopara=if_else(log2FoldChange<0, "No", "Yes"))

# Figure 4. Identified by the DESeq2 tool significantly different abundant (adjusted p < 0.01) bacterial ASVs between grape samples which grouped based on presence of Plasmopara viticola
ggplot(data=de_deseq_significant) +
  geom_col(aes(x=Genus, y=0, fill=Plasmopara),alpha=0.3) +
  scale_fill_manual(values = c("#401CE6","#e50914")) +
  geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=0,ymax=-Inf),alpha=0.005, color="black", fill="#401CE6") +
  geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=0,ymax=Inf),alpha=0.005, color="black", fill="#e50914") +
  geom_point(aes(x=Genus, y=log2FoldChange, color=Class), size=4) + 
  theme_linedraw() +
  scale_color_manual(values = c("#490a3d","#bd1550","#e97f02","#8a9b0f")) +
  theme(axis.text.y = element_text(color="black", face="italic"), legend.text=element_text(face="italic"), 
    legend.title = element_text(face = "italic"), legend.position="bottom", legend.box="vertical"
  ) + 
    guides(color=guide_legend(ncol=2, order = 1)) +
  labs(fill="Presence of Plasmopara viticola", y="log2 fold change") +
  coord_flip() 

ggsave("Figure 4.png", width = 8, height = 6)
#


physeq_ITS = merge_samples(physeq_ITS, "Sample", fun = mean)
sample_data(physeq_ITS) = sample_data(plasmopara)

# DESeq2 process (ITS)
deseq2_ITS = phyloseq_to_deseq2(physeq_ITS, ~ Plasmopara)
deseq2_ITS = DESeq(deseq2_ITS, test="Wald", fitType="parametric")
#

# Filter significant results of DESeq2 (ITS)
res = results(deseq2_ITS, cooksCutoff = FALSE)
res = as.data.frame(res)
de_deseq_significant = filter(res, padj < 0.01)
de_deseq_significant$ASV = row.names(de_deseq_significant)

significant_tax_ITS = as.data.frame(tax_table(physeq_ITS)[de_deseq_significant$ASV,])
significant_tax_ITS$ASV = row.names(significant_tax_ITS)
de_deseq_significant = left_join(de_deseq_significant, significant_tax_ITS)

de_deseq_significant$Genus = factor(de_deseq_significant$Genus, 
    unique(de_deseq_significant$Genus[order(de_deseq_significant$log2FoldChange, decreasing = TRUE)]))


seq_ITS = filter(seq_ITS, ASV %in% significant_tax_ITS$ASV)
significant_seq_ITS = left_join(seq_ITS, de_deseq_significant)
write.table(significant_seq_ITS, file="Supporting Information Table S5.txt")
#

de_deseq_significant=mutate(de_deseq_significant, Plasmopara=if_else(log2FoldChange<0, "No", "Yes"))

# Figure 5. Identified by the DESeq2 tool significantly different abundant (adjusted p < 0.01) fungal ASVs between grape samples which grouped based on presence of Plasmopara viticola
ggplot(data=de_deseq_significant) +
  geom_col(aes(x=Genus, y=0, fill=Plasmopara),alpha=0.3)+
  scale_fill_manual(values = c("#401CE6","#e50914")) +
  geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=0,ymax=-Inf),alpha=0.021, color="black", fill="#401CE6")+
  geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=0,ymax=Inf),alpha=0.021, color="black", fill="#e50914")+
  geom_point(aes(x=Genus, y=log2FoldChange, color=Class), size=4) + 
  theme_linedraw() +
  scale_color_manual(values = c("#655643","#e50914","#e6ac27","#bf4d28")) +
  theme(axis.text.y = element_text(color="black", face="italic"), legend.text=element_text(face="italic"), 
    legend.title = element_text(face = "italic"), legend.position="bottom", legend.box="vertical"
  ) + 
  guides(color=guide_legend(ncol=2, order = 1)) +
  labs(fill="Presence of Plasmopara viticola", y="log2 fold change") +
  coord_flip() 

ggsave("Figure 5.png", width = 7, height = 6)