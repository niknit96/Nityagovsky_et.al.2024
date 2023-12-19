library("tidyverse")
library ("phyloseq")

dir = as.character(read.table("dir.txt"))
setwd(dir)

load("./ITS/physeq_ITS.RData") # Load ITS data

# Obtain relative abundances of Plasmopara viticola in samples
Barplot_plasmopara = merge_samples(physeq_ITS, "Sample", fun = mean)
Barplot_plasmopara = transform_sample_counts(Barplot_plasmopara, function(x) x / sum(x) * 100) # Transform abundance data to relative abundances

Barplot_plasmopara = psmelt(Barplot_plasmopara) %>%
	group_by(., Genus, Sample) %>%
	summarise(., Abundance = sum(Abundance))
Barplot_plasmopara = filter(Barplot_plasmopara, Genus %in% grep("Plasmopara", Genus, value = TRUE))
Barplot_plasmopara$Genus = "Plasmopara viticola"
Barplot_plasmopara = Barplot_plasmopara %>%
	group_by(., Genus, Sample) %>%
	summarise(., Abundance = sum(Abundance))

Barplot_plasmopara = left_join(Barplot_plasmopara, sampledata_ITS, by="Sample")
#

# For Plasmopara abundance decreasing sorting in barplot
Barplot_plasmopara = Barplot_plasmopara[order(Barplot_plasmopara$Abundance, decreasing = TRUE), ]
Barplot_plasmopara$Sample = factor(x = Barplot_plasmopara$Sample, levels = Barplot_plasmopara$Sample)
#

# For location sorting
Barplot_plasmopara = Barplot_plasmopara %>%
	mutate(Organ = if_else(Organ.material == "Leaf", "L", "S"))

Location_short = c("Rikord Island (P-4)",
"Russky Island (P-3)",
"Vladivostok (Gh)",
"Vladivostok (P-1, P-2)",
"Vineyard Makarevich (M, M-dm, Ad, Muk)",
"Ivanovka village (P-5)",
"Vineyard PRIM ORGANICA (Alfa, Pr-St)",
'The Verkhne-Ussuriysky Research Station (P-6)',
"Near the city of Nevelsk (S-3)",
"Near the city of Kholmsk (S-2)",
"The botanical garden on Sakhalin Island (S-Va, S-1)",
"Litovko village (Kh-1)",
"Silinsky forest (Kh-2)")

Places = c("Rikord Island (P-4)",
"Russky Island (P-3)",
"Vladivostok (Gh)",
"Vladivostok (P-1, P-2)",
"Vineyard Makarevich\n(M, M-dm, Ad, Muk)",
"Ivanovka village (P-5)",
"Vineyard PRIM ORGANICA\n(Alfa, Pr-St)",
'The Verkhne-Ussuriysky\nResearch Station (P-6)',
"Near the city of Nevelsk (S-3)",
"Near the city of Kholmsk (S-2)",
"The botanical garden\non Sakhalin Island (S-Va, S-1)",
"Litovko village (Kh-1)",
"Silinsky forest (Kh-2)")

Places_num = paste(1:13, ". ", Places, sep="")
Places_for_left_join = cbind(Location_short, Places_num)

colnames(Places_for_left_join) = c("Location_short", "Places")
Places_for_left_join = as.data.frame(Places_for_left_join)


Barplot_plasmopara = left_join(Barplot_plasmopara, Places_for_left_join, by="Location_short")

Barplot_plasmopara$Places = factor(Barplot_plasmopara$Places,Places_num)
#

# Figure 1b. Relative abundance of P. viticola in samples
ggplot(Barplot_plasmopara) +
	geom_col(aes(x = Sample, y = Abundance, fill = Genus), colour = "black", linetype = "solid", size = 0.1) +
	scale_fill_manual(values="#e50914") +
	facet_wrap(
		vars(Places),
		ncol =3,
		scales = "free_x",
		shrink = TRUE,
		labeller = "label_value",
		as.table = TRUE,
		switch = NULL,
		drop = TRUE
	) +
	guides(fill = guide_legend(ncol = 1, byrow = TRUE)) +
	theme_bw() +
	theme(
		legend.text = element_text(size=28, face="bold.italic"),
		legend.position="top",
		legend.title = element_text(size=18),
		axis.text.x = element_text(size=1, angle=0, vjust = 0.5, color = "#00000000"),
		axis.text.y = element_text(size=18, color = "black"),
		strip.text = element_text(size=25),
		axis.title.x = element_text(size=18),
		axis.title.y = element_text(size=28)
	) +
	labs(x = "", y = "Relative abundance, %", fill = "") +
	scale_y_continuous(labels = scales::percent_format(scale = 1)) +
	geom_text(aes(x = Sample, y = Abundance, 
		label = ifelse(Abundance > 0, scales::percent(Abundance, scale = 1, accuracy = 0.1), ""), 
		group = Genus), 
		inherit.aes = FALSE, position = position_nudge(x = 0, y = 2), 
		size = 3.2, angle=0, col = "black"
	) + 
	geom_text(aes(x = Sample, y = -2, label = Organ), 
		size = 4, angle=0, col = "black")

ggsave("Figure 1b.png", width = 20, height = 20)
#