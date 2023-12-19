library("tidyverse")
library ("phyloseq")
library("vegan")
library("ggpubr")

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

physeq_16s = merge_samples(physeq_16s, "Sample", fun = mean)

# Transform to even sampling depth (16s data)
total = median(sample_sums(physeq_16s))
standf = function(x, t=total) round(t * (x / sum(x)))
physeq_16s = transform_sample_counts(physeq_16s, standf)
#

# Preprocess for Bray–Curtis dissimilarity index computetion (16s data)
physeq_16s = psmelt(physeq_16s) %>%
	group_by(., Genus, Sample) %>%
	summarise(., Abundance = sum(Abundance))

physeq_16s = pivot_wider(physeq_16s, id_cols = Genus, names_from = "Sample", values_from = "Abundance")
physeq_16s = as.data.frame(physeq_16s)
row.names(physeq_16s) = physeq_16s[,1]
physeq_16s = physeq_16s[,-1]
physeq_16s = t(physeq_16s)
#

vdist = vegdist(physeq_16s, method="bray") # The function computes Bray–Curtis dissimilarity index (16s data)

vdist_NMDS = metaMDS(physeq_16s, distance = "bray", k = 2, maxit = 999,  trymax = 500, wascores = TRUE) # Converting to nonmetric multidimensional scaling (NMDS)
NMDS = as.data.frame(scores(vdist_NMDS, display="site"))
NMDS$Sample = rownames(NMDS)

NMDS = left_join(NMDS, plasmopara, by ="Sample")

mycolors = c("#401CE6", "#FF1202")

# Figure 3a. Bray–Curtis beta diversity NMDS plot of grape endophytic bacteria
ggscatter(NMDS, x = "NMDS1", y = "NMDS2",
    color = "Plasmopara",
    mean.point = TRUE, ellipse = TRUE) +
    stat_stars(aes(color = Plasmopara)) +
    scale_colour_manual(values = mycolors) +
    fill_palette(mycolors) +
    theme(legend.text = element_text(size = 15, face = "italic"),
        legend.title = element_text(face = "italic"),
        title = element_text(face = "italic"),
        plot.title = element_text(hjust=0.5),
        legend.position="bottom"
    ) + 
    labs(title="Bray–Curtis beta diversity (16s)", 
        color = "Presence of Plasmopara viticola", 
        fill = "Presence of Plasmopara viticola")

ggsave("Figure 3a.png", width = 6, height = 6)
#



physeq_ITS = merge_samples(physeq_ITS, "Sample", fun = mean)

# Transform to even sampling depth (ITS data)
total = median(sample_sums(physeq_ITS))
standf = function(x, t=total) round(t * (x / sum(x)))
physeq_ITS = transform_sample_counts(physeq_ITS, standf)
#

# Preprocess for Bray–Curtis dissimilarity index computetion (ITS data)
physeq_ITS = psmelt(physeq_ITS) %>%
	group_by(., Genus, Sample) %>%
	summarise(., Abundance = sum(Abundance))

physeq_ITS = pivot_wider(physeq_ITS, id_cols = Genus, names_from = "Sample", values_from = "Abundance")
physeq_ITS = as.data.frame(physeq_ITS)
row.names(physeq_ITS) = physeq_ITS[,1]
physeq_ITS = physeq_ITS[,-1]
physeq_ITS = t(physeq_ITS)
#

vdist = vegdist(physeq_ITS, method="bray") # The function computes Bray–Curtis dissimilarity index (ITS data)

vdist_NMDS = metaMDS(physeq_ITS, distance = "bray", k = 2, maxit = 999,  trymax = 500, wascores = TRUE) # Converting to nonmetric multidimensional scaling (NMDS)
NMDS = as.data.frame(scores(vdist_NMDS, display="site"))
NMDS$Sample = rownames(NMDS)

NMDS = left_join(NMDS, plasmopara, by ="Sample")

mycolors = c("#401CE6", "#FF1202")


# Figure 3b. Bray-Curtis beta diversity NMDS plot of grape endophytic fungi
ggscatter(NMDS, x = "NMDS1", y = "NMDS2",
    color = "Plasmopara",
    mean.point = TRUE, ellipse = TRUE) +
    stat_stars(aes(color = Plasmopara)) +
    scale_colour_manual(values = mycolors) +
    fill_palette(mycolors) +
    theme(legend.text = element_text(size = 15, face = "italic"),
        legend.title = element_text(face = "italic"),
        title = element_text(face = "italic"),
        plot.title = element_text(hjust=0.5),
        legend.position="bottom"
    ) + 
    labs(title="Bray–Curtis beta diversity (ITS)", 
        color = "Presence of Plasmopara viticola", 
        fill = "Presence of Plasmopara viticola")

ggsave("Figure 3b.png", width = 6, height = 6)
#


# PERMANOVA results (16s)
print("PERMANOVA results (16s)")
adonis2(physeq_16s ~ Plasmopara, data = NMDS, permutations = 999, method = "bray")
#

# PERMANOVA results (ITS)
print("PERMANOVA results (ITS)")
adonis2(physeq_ITS ~ Plasmopara, data = NMDS, permutations = 999, method = "bray")
#