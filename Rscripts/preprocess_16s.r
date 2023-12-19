library("tidyverse")
library("qiime2R")
library ("phyloseq")

### Load data to R
## Begin

dir = as.character(read.table("dir.txt"))
setwd(dir)



asv_meta = read_qza("./16s/dada2/FeatureTable[Frequency]_16s.qza") # Table of metagenome data reads for each ASV by samples

asv_meta = cbind(row.names(asv_meta$data), asv_meta$data)
asv_meta = as.data.frame(asv_meta)
colnames(asv_meta)[1] = "Species"
row.names(asv_meta) = asv_meta[,1]
asv_meta = asv_meta[,-1]
asv_meta[] = apply(asv_meta[], 2, as.numeric)

tax_meta = read_qza("./16s/feature-classifier_classify-sklearn/FeatureData[Taxonomy]_16s.qza") # Taxonomy table for metagenome data
tax_meta = parse_taxonomy(tax_meta$data, trim_extra=FALSE)
tax_meta[is.na(tax_meta)] <- "uncultured"

# Swap values ​​from Swap_unidentified in tax_meta to previous taxon level
Swap_unidentified = c("uncultured", "unidentified", "metagenome", "bacteriap25", "Unknown")
for(Swap in Swap_unidentified) {
    for(j in 1:7) {
        for (i in 1:nrow(tax_meta)) {
            if(grepl(Swap, tax_meta[i,j])) {
                tax_meta[i,j] = tax_meta[i,j-1] }
                else { tax_meta[i,j] = tax_meta[i,j]
            }       
        }
    }
}


#

# Filtering metagenome data from non-significant taxa
tax_meta = filter(tax_meta, 
Phylum != "d__Bacteria" & 
Phylum != "Unassigned" &
Kingdom != "d__Eukaryota" & 
Phylum != "k__Fungi" & 
Kingdom != "d__Archaea" &
Kingdom != "k__Viridiplantae" &
Genus != "g__Mitochondria" &
Genus != "g__Chloroplast")
#

tax_meta["Other",] = c("Other","Other","Other","Other","Other","Other","Other")

# Name ASVs to Genus level
tax_meta[grep("f__", tax_meta[,"Genus"]),"Genus"] = make.unique(tax_meta[grep("f__", tax_meta[,"Genus"]), "Genus"])
tax_meta[grep("o__", tax_meta[,"Genus"]),"Genus"] = make.unique(tax_meta[grep("o__", tax_meta[,"Genus"]), "Genus"])
tax_meta[grep("c__", tax_meta[,"Genus"]),"Genus"] = make.unique(tax_meta[grep("c__", tax_meta[,"Genus"]), "Genus"])
tax_meta[grep("p__", tax_meta[,"Genus"]),"Genus"] = make.unique(tax_meta[grep("p__", tax_meta[,"Genus"]), "Genus"])
tax_meta[,"Genus"] = str_replace_all(tax_meta[,"Genus"], "(\\d+)$", function(x) as.numeric(x) + 1)
#

sampledata_16s = read.table(file="./Metagenome metadata/sampledata_16s_meta.txt", header = TRUE, colClasses = "character", sep="\t") # Sample metadata for the metagenome
row.names(sampledata_16s) = sampledata_16s$SRA


# Creating a phyloseq object for metagenome data
sampledata_meta = sample_data(sampledata_16s)


tax_16s = tax_table(as.matrix(tax_meta))
asv_16s = otu_table(asv_meta, taxa_are_rows = TRUE)
physeq_16s = phyloseq(asv_16s, tax_16s, sampledata_meta)

#
## End

filtered_data = as.data.frame(colSums(otu_table(physeq_16s)))
filtered_data = as.data.frame(t(t(filtered_data)))
colnames(filtered_data) = "Sequences after filtration used in analysis"
filtered_data$Sample = row.names(filtered_data)
raw_data = read_qza("./16s/dada2/16s-stats-dada2.qza")$data
raw_data$Sample = row.names(raw_data) 
data_16s = left_join(raw_data,filtered_data, by="Sample")
data_16s = data_16s[,c("Sample", "input", "Sequences after filtration used in analysis")]
colnames(data_16s) = c("SRA", "Raw paired-end reads", "Sequences after filtration used in analysis")
data_16s$Sum_raw = sum(data_16s[,"Raw paired-end reads"])
data_16s$Mean_raw = mean(data_16s[,"Raw paired-end reads"])
data_16s$Median_raw = median(data_16s[,"Raw paired-end reads"])
data_16s$Sum_filtered = sum(data_16s[,"Sequences after filtration used in analysis"])
data_16s$Mean_filtered = mean(data_16s[,"Sequences after filtration used in analysis"])
data_16s$Median_filtered = median(data_16s[,"Sequences after filtration used in analysis"])

sampledata_16s = left_join(sampledata_16s, data_16s, by="SRA")

write.table(sampledata_16s, file="Sampledata_16s_with_summary.txt", quote=FALSE, sep="\t")

save.image(file='./16s/physeq_16s.RData')