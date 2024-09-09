## Rafael Valentin
## InterBIO project 2024
# Code generated to create a series of Manhattan plots for clinical data


################################################################################
# Manhattan plots
################################################################################

## ===================================
# Run through M/B-values for Kruskal-Wallis
# testing for sig CpGs for all grps
## ===================================

# Clear workspace to free up memory
rm(list = setdiff(ls(), c("Kennedy_CpG.M_values", "Kennedy_CpG.B_values", "Output")))
gc()

# Read in additional covariate data
Covariate_data = read.csv("INTERBIO_cord_ewas_covars.csv")
Master_sheet = read.csv("20230726_kennedy_master_sample_sheet.csv")

##** DOESN"T WORK SINCE NO COMMON ID NAMES!!!! **##
#Covariate_data = merge(Covariate_data, Master_sheet[,c(2,10)], by.x = "Sample_ID", by.y = "Sample_Name") 

# Fix "Output" to contain additional data just read in
##** Nothing to add until Master and Data align **##

# Create object (M-values)
test.frame = as.data.frame(t(Kennedy_CpG.M_values))
test.frame = merge(Output[,c(1,12)], test.frame, by.x = "SampleID", by.y = "row.names")

# Scaled M-value version
test.frame2 = caret::preProcess(test.frame, method = c("scale"))
test.frame2 = predict(test.frame2, test.frame)

test.frame3 = as.data.frame(t(Kennedy_CpG.B_values))
test.frame3 = merge(Output[,c(1,12)], test.frame3, by.x = "SampleID", by.y = "row.names")
test.frame3 = test.frame3[,colnames(test.frame3) %in% colnames(test.frame)]
gc()

# Run analysis (no correction - unscaled, scaled, and B-values)
Kruskal.test_M.vals = apply(test.frame[3:ncol(test.frame)], 2, function(i){
  kruskal.test(i ~ test.frame$phenotype)
})

Kruskal.test_B.vals = apply(test.frame3[3:ncol(test.frame3)], 2, function(i){
  kruskal.test(i ~ test.frame3$phenotype)
})
gc()


## -----------------------
# Pull out just p-vals
## -----------------------
library(dplyr)

# Unscaled M-values
Sig.results_M.vals = lapply(Kruskal.test_M.vals, "[[", 3)
Sig.results_M.vals = as.data.frame(t(dplyr::bind_rows(Sig.results_M.vals)))
Sig.results_M.vals = merge(Sig.results_M.vals, IlluminaHumanMethylationEPICanno.ilm10b5.hg38::Locations, by = "row.names")
colnames(Sig.results_M.vals) = c("CpG Site", "p.value", "Chromosome", "Position", "Strand")
Sig.results_M.vals = Sig.results_M.vals[Sig.results_M.vals$Chromosome != "chrY",]
Sig.results_M.vals$Chromosome = factor(Sig.results_M.vals$Chromosome,
                                       levels = c("chr1", "chr2", "chr3",
                                                  "chr4", "chr5", "chr6",
                                                  "chr7", "chr8", "chr9",
                                                  "chr10", "chr11", "chr12",
                                                  "chr13", "chr14", "chr15",
                                                  "chr16", "chr17", "chr18",
                                                  "chr19", "chr20", "chr21",
                                                  "chr22", "chrX"))
Sig.results_M.vals$p.adj = p.adjust(Sig.results_M.vals$p.value, method = "bonferroni", n=nrow(Sig.results_M.vals))
Sig.results_M.vals = Sig.results_M.vals %>% relocate(p.adj, .after=p.value)

Chromosome_Lengths = Sig.results_M.vals %>% group_by(Chromosome) %>%
  summarise(chr_len = max(Position)) %>%
  mutate(tot = cumsum(as.numeric(chr_len))-chr_len) %>% select(-chr_len)

Sig.results_M.vals = merge(Sig.results_M.vals, Chromosome_Lengths, by = "Chromosome")
rm(Chromosome_Lengths)

Sig.results_M.vals$BPcum = Sig.results_M.vals$Position + Sig.results_M.vals$tot
Sig.results_M.vals$is_sig.base = ifelse(Sig.results_M.vals$p.value <= 0.05/nrow(Sig.results_M.vals), "yes", "no")
Sig.results_M.vals$is_sig.adj = ifelse(Sig.results_M.vals$p.adj <= 0.05, "yes", "no")

Sig.results_M.vals = merge(Sig.results_M.vals, IlluminaHumanMethylationEPICanno.ilm10b5.hg38::Other, by.x = "CpG Site", by.y = "row.names")


# B-values
Sig.results_B.vals = lapply(Kruskal.test_B.vals, "[[", 3)
Sig.results_B.vals = as.data.frame(t(dplyr::bind_rows(Sig.results_B.vals)))
Sig.results_B.vals = merge(Sig.results_B.vals, IlluminaHumanMethylationEPICanno.ilm10b5.hg38::Locations, by = "row.names")
colnames(Sig.results_B.vals) = c("CpG Site", "p.value", "Chromosome", "Position", "Strand")
Sig.results_B.vals = Sig.results_B.vals[Sig.results_B.vals$Chromosome != "chrY",]
Sig.results_B.vals$Chromosome = factor(Sig.results_B.vals$Chromosome,
                                       levels = c("chr1", "chr2", "chr3",
                                                  "chr4", "chr5", "chr6",
                                                  "chr7", "chr8", "chr9",
                                                  "chr10", "chr11", "chr12",
                                                  "chr13", "chr14", "chr15",
                                                  "chr16", "chr17", "chr18",
                                                  "chr19", "chr20", "chr21",
                                                  "chr22", "chrX"))
Sig.results_B.vals$p.adj = p.adjust(Sig.results_B.vals$p.value, method = "bonferroni", n=nrow(Sig.results_B.vals))
Sig.results_B.vals = Sig.results_B.vals %>% relocate(p.adj, .after=p.value)

Chromosome_Lengths = Sig.results_B.vals %>% group_by(Chromosome) %>%
  summarise(chr_len = max(Position)) %>%
  mutate(tot = cumsum(as.numeric(chr_len))-chr_len) %>% select(-chr_len)

Sig.results_B.vals = merge(Sig.results_B.vals, Chromosome_Lengths, by = "Chromosome")
rm(Chromosome_Lengths)

Sig.results_B.vals$BPcum = Sig.results_B.vals$Position + Sig.results_B.vals$tot
Sig.results_B.vals$is_sig.base = ifelse(Sig.results_B.vals$p.value <= 0.05/nrow(Sig.results_B.vals), "yes", "no")
Sig.results_B.vals$is_sig.adj = ifelse(Sig.results_B.vals$p.adj <= 0.05, "yes", "no")

Sig.results_B.vals = merge(Sig.results_B.vals, IlluminaHumanMethylationEPICanno.ilm10b5.hg38::Other, by.x = "CpG Site", by.y = "row.names")
gc()


## -----------------------
# Generate Manhattan plot
## -----------------------
library(ggplot2)
library(dplyr)

##** Unscaled M-Value**##
axisdf = Sig.results_M.vals %>% group_by(Chromosome) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
ggplot(Sig.results_M.vals, aes(x=BPcum, y = log10(p.value)*-1)) +
  geom_point(aes(color = Chromosome), alpha =0.8, size=1.3) +
  scale_color_manual(values = rep(c("lightblue", "steelblue"), length(levels(Sig.results_M.vals$Chromosome)))) +
  scale_x_continuous(labels = axisdf$Chromosome, breaks = axisdf$center) +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw() + theme(legend.position = "none",
                     axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10),
                     axis.text.y = element_text(size = 10)) +
  geom_hline(yintercept = -log10(0.05/nrow(Sig.results_M.vals)), linetype = 2, color = "red") +
  geom_hline(yintercept = -log10(1/nrow(Sig.results_M.vals)), linetype = 2, color = "black") +
  ggrepel::geom_label_repel(data = subset(Sig.results_M.vals, is_sig.adj=="yes"), aes(label = `CpG Site`), size = 3, max.overlaps = 15) +
  labs(x = "Chromosome", y = "-Log10(P)",
       title = "Significant CpGs across all Fetal Growth Trajectories",
       subtitle = "Kruskal-Wallis test - Unscaled M-values")
rm(axisdf)


##** B-value **##
axisdf = Sig.results_B.vals %>% group_by(Chromosome) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
ggplot(Sig.results_B.vals, aes(x=BPcum, y = log10(p.value)*-1)) +
  geom_point(aes(color = Chromosome), alpha =0.8, size=1.3) +
  scale_color_manual(values = rep(c("lightblue", "steelblue"), length(levels(Sig.results_B.vals$Chromosome)))) +
  scale_x_continuous(labels = axisdf$Chromosome, breaks = axisdf$center) +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw() + theme(legend.position = "none",
                     axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10),
                     axis.text.y = element_text(size = 10)) +
  geom_hline(yintercept = -log10(0.05/nrow(Sig.results_B.vals)), linetype = 2, color = "red") +
  geom_hline(yintercept = -log10(1/nrow(Sig.results_B.vals)), linetype = 2, color = "black") +
  ggrepel::geom_label_repel(data = subset(Sig.results_B.vals, is_sig.adj=="yes"), aes(label = `CpG Site`), size = 3, max.overlaps = 15) +
  labs(x = "Chromosome", y = "-Log10(P)",
       title = "Significant CpGs across all Fetal Growth Trajectories",
       subtitle = "Kruskal-Wallis test - B-values")
rm(axisdf)
