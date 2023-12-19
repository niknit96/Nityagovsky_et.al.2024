library("tidyverse")
library ("phyloseq")

set.seed(711)

dir = as.character(read.table("dir.txt"))
setwd(dir)

load("./16s/physeq_16s.RData") # Load 16s data
load("./ITS/physeq_ITS.RData") # Load ITS data


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

# Resample an ASV table such that all samples have the same library size
physeq_16s = rarefy_even_depth(physeq_16s)
physeq_ITS = rarefy_even_depth(physeq_ITS)
#

alpha_diversity_16s = merge_samples(physeq_16s, "Sample", fun = mean)
sample_data(alpha_diversity_16s)$Sample = sample_names(alpha_diversity_16s)

# Performs a number of standard alpha diversity estimates for 16s data
alpha_diversity_16s = estimate_richness(alpha_diversity_16s)
#

# Obtain Pielou’s evenness index in 16s data
alpha_diversity_16s$Pielou = alpha_diversity_16s$Shannon/log(alpha_diversity_16s$Observed)
#

alpha_diversity_16s$Sample = rownames(alpha_diversity_16s)
alpha_diversity_16s = left_join(alpha_diversity_16s, plasmopara)

alpha_diversity_ITS = merge_samples(physeq_ITS, "Sample", fun = mean)
sample_data(alpha_diversity_ITS)$Sample = sample_names(alpha_diversity_ITS)

# Performs a number of standard alpha diversity estimates for ITS data
alpha_diversity_ITS = estimate_richness(alpha_diversity_ITS)
#

# Obtain Pielou’s evenness index in ITS data
alpha_diversity_ITS$Pielou = alpha_diversity_ITS$Shannon/log(alpha_diversity_ITS$Observed)
#

alpha_diversity_ITS$Sample = rownames(alpha_diversity_ITS)
alpha_diversity_ITS = left_join(alpha_diversity_ITS, plasmopara)

# Figure 2b. Pielou's evenness index (16s)
ggplot(data = alpha_diversity_16s, aes(x = Plasmopara, y = Pielou)) + 
	geom_boxplot() +
	geom_jitter(width = 0.2) +
	theme_bw() +
	scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1), limits = c(0,1.25)) +
	theme(text = element_text(size=25), axis.title.x=element_text(face = "italic"),
    	axis.text.x = element_text(size=30, color="#000000"),
    	title=element_text(size=23, face = "italic")
	) +
	labs(x = "Presence of Plasmopara viticola", y = "") +
	ggtitle("Pielou's evenness index (16s)")
ggsave("Figure 2b.png", width = 7, height = 7)
#

# Figure 2d. Pielou's evenness index (ITS)
ggplot(data = alpha_diversity_ITS, aes(x = Plasmopara, y = Pielou)) + 
	geom_boxplot() +
	geom_jitter(width = 0.2) +
	theme_bw() +
	scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1), limits = c(0,1.25)) +
	theme(text = element_text(size=25), axis.title.x=element_text(face = "italic"),
    	axis.text.x = element_text(size=30, color="#000000"),
    	title=element_text(size=23, face = "italic")
	) +
	labs(x = "Presence of Plasmopara viticola", y = "") +
	ggtitle("Pielou's evenness index (ITS)")
ggsave("Figure 2d.png", width = 7, height = 7)
#

# Figure 2a. Number of ASVs (16s)
ggplot(data = alpha_diversity_16s, aes(x = Plasmopara, y = Observed)) + 
	geom_boxplot() +
	geom_jitter(width = 0.2) +
	theme_bw() +
	scale_y_continuous(breaks = c(0,100,200), limits = c(0,300)) +
	theme(text = element_text(size=25), axis.title.x=element_text(face = "italic"),
    	axis.text.x = element_text(size=30, color="#000000"),
    	title=element_text(size=23, face = "italic")
	) +
	labs(x = "Presence of Plasmopara viticola", y = "") +
	ggtitle("Number of ASVs (16s)")
ggsave("Figure 2a.png", width = 7, height = 7)
#

# Figure 2c. Number of ASVs (ITS)
ggplot(data = alpha_diversity_ITS, aes(x = Plasmopara, y = Observed)) + 
	geom_boxplot() +
	geom_jitter(width = 0.2) +
	theme_bw() +
	scale_y_continuous(breaks = c(0,100,200), limits = c(0,300)) +
	theme(text = element_text(size=25), axis.title.x=element_text(face = "italic"),
    	axis.text.x = element_text(size=30, color="#000000"),
    	title=element_text(size=23, face = "italic")
	) +
	labs(x = "Presence of Plasmopara viticola", y = "") +
	ggtitle("Number of ASVs (ITS)")
ggsave("Figure 2c.png", width = 7, height = 7)


print("Wilcoxon rank sum test. Pielou's evenness index (16s)")
pairwise.wilcox.test(alpha_diversity_16s$Pielou, alpha_diversity_16s$Plasmopara, p.adjust.method = "none", paired=FALSE)

print("Wilcoxon rank sum test. Number of ASVs (16s)")
pairwise.wilcox.test(alpha_diversity_16s$Observed, alpha_diversity_16s$Plasmopara, p.adjust.method = "none", paired=FALSE)

print("Wilcoxon rank sum test. Pielou's evenness index (ITS)")
pairwise.wilcox.test(alpha_diversity_ITS$Pielou, alpha_diversity_ITS$Plasmopara, p.adjust.method = "none", paired=FALSE)

print("Wilcoxon rank sum test. Number of ASVs (ITS)")
pairwise.wilcox.test(alpha_diversity_ITS$Observed, alpha_diversity_ITS$Plasmopara, p.adjust.method = "none", paired=FALSE)