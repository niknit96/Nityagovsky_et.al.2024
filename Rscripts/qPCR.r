library("tidyverse")
library("multcompView")
library("ggpubr")
library ("phyloseq")

dir = as.character(read.table("dir.txt"))
setwd(dir)

qPCR = readRDS('./qPCR.RDS') # load qPCR data
plasmopara = readRDS('./plasmopara_NGS.RDS') # load Plasmopara abundance data

# ANOVA for qPCR data ('Plant' factor)
final_tab = NULL
for(amplicon in c("PvITS1_1", "PvITS1_2", "PvCox1_1")) {

    tab_for_anova = filter(qPCR, Amplicon %in% amplicon)

    mod = aov(Relative.expression ~ Plant, data = tab_for_anova)
    tukey = TukeyHSD(exp_aov <- mod)
    anova_res = multcompLetters4(exp_aov, tukey, reversed = TRUE)

    anova_res = as.data.frame(anova_res$Plant$Letters)
    anova_res$Plant = row.names(anova_res)
    colnames(anova_res) = c("Anova", "Plant")

    tab_for_anova = left_join(tab_for_anova, anova_res, by=c("Plant"))
    final_tab = bind_rows(final_tab, tab_for_anova)
}
#

final_tab$Plant = str_replace_all(final_tab$Plant, "__", "-")

# Measurement of mean relative amplification level and standart error in qPCR data
qPCR = final_tab %>%
    group_by(Plant, Amplicon) %>%
    summarise(Relative.expression_mean = mean(Relative.expression), se = sd(Relative.expression)/sqrt(length((Relative.expression))), Anova = Anova) %>%
    group_by(Plant, Amplicon, Relative.expression_mean, se, Anova) %>%
    summarise()
#

plasmopara$Plant = str_replace_all(plasmopara$Plant, "__", "-")

# Measurement of mean Plasmopara viticola abundance and standart error in NGS data
plasmopara = plasmopara %>%
    group_by(Plant) %>%
    summarise(Plasmopara_mean = mean(Plasmopara), Plasmopara_se = sd(Plasmopara)/sqrt(length((Plasmopara)))) %>%
    unique()
#

final_tab$Plant = str_replace_all(final_tab$Plant, "__", "-")

final_tab = left_join(qPCR, plasmopara, by="Plant") # General table

sample_order=c("Gh","P-1","P-2",
    "P-3","P-4","P-5","P-6",
    "Kh-1","Kh-2","S-Va","S-1",
    "S-2","S-3","Ad","Muk","Alfa",
    "Pr-St","M-dm","Mildew","Negative (mix)")

final_tab$Plant = factor(final_tab$Plant, sample_order)
final_tab$Amplicon = factor(final_tab$Amplicon, c("PvITS1_1", "PvITS1_2", "PvCox1_1"))

# Figure 6a. Quantification of the amplification of PvITS1_1, PvITS1_2 and PvCox1_1 in the DNAs of grapevine samples performed by qPCR
ggplot(final_tab) +
	geom_col(aes(x = Plant, y = Relative.expression_mean), colour = "black", 
        linetype = "solid", size = 0.1, position = position_dodge(width=0.9)
    ) +
    facet_grid(vars(Amplicon)) +
    geom_errorbar(aes(x = Plant, ymin = Relative.expression_mean-se, ymax = Relative.expression_mean+se), 
        position = position_dodge(width=0.9), width = 0.25
    ) +
    theme_bw() +
    labs(title="qPCR SybrGreen detecting of Plasmopara viticola", 
        x = "Plants", 
        y="Relative level of amplification of\nP. viticola amplicon, r.u."
    ) +
    theme(title = element_text(face = "italic"), strip.text = element_text(size=14),
        axis.text.y = element_text(face = "italic"),
        axis.text.x = element_text(angle=40, vjust = 0.6, hjust = 0.7)
    ) +
    geom_text(aes(x = Plant, y = Relative.expression_mean+se+0.1, label = Anova), inherit.aes = FALSE, position = position_dodge(width=0.9))


ggsave("Figure 6a.png", width = 7, height = 7, dpi=300)
#

# Figure 6b. The relative representation of Plasmopara viticola in NGS samples
ggplot(final_tab) +
	geom_col(aes(x = Plant, y = Plasmopara_mean), colour = "black", 
        linetype = "solid", size = 0.1, position = position_dodge(width=0.9)
    ) +
    geom_errorbar(aes(x = Plant, ymin = Plasmopara_mean-Plasmopara_se, ymax = Plasmopara_mean+Plasmopara_se), 
        position = position_dodge(width=0.9), width = 0.25
    ) +
    geom_text(data=filter(final_tab, Plant %in% c("Mildew", "Negative (mix)")), aes(x=Plant, y=2, label="n.m."), size=3) +
    theme_bw() +
    labs(title="NGS detecting of Plasmopara viticola", 
        x = "Plants", 
        y="Relative abundance of P. viticola, %"
    ) +
    theme(title = element_text(face = "italic"),
        axis.text.y = element_text(face = "italic"),
        axis.text.x = element_text(angle=40, vjust = 0.6, hjust = 0.7)
    )

ggsave("Figure 6b.png", width = 7, height = 5, dpi=300)
#

# Figure 7. Pearson correlation coefficient
ggscatter(final_tab, y = "Relative.expression_mean", x = "Plasmopara_mean", size = 0.3,
          color = "Amplicon", palette = "jco",
          facet.by = "Amplicon",
          add = "reg.line", conf.int = FALSE
        ) +
        stat_cor(method = "pearson",  p.accuracy = 0.001) +
        theme(title = element_text(face = "italic"),
            axis.text.x = element_text(face = "italic"),
            axis.text.y = element_text(face = "italic")
        ) +
        labs(color = "Amplicon", 
            y = "Mean relative level of amplification of\nP. viticola amplicons in qPCR data, r.u.", 
            x="Mean relative abundance of\nP. viticola in NGS data, %"
        )

ggsave("Figure 7.png", width = 9, height = 5, dpi=300)
#

