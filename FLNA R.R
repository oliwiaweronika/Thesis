install.packages('BiocManager')
BiocManager::install("phyloseq")
library("phyloseq")
library("ggplot2")
library("scales")
install.packages("ggpubr")
library("ggpubr")
install.packages("vegan")
library("vegan")
library("dplyr")
BiocManager::install("DESeq2")
library(DESeq2)
install.packages("remotes")
remotes::install_github("vmikk/metagMisc")
library(metagMisc)
library(ggsci)

#define colors for each tretament to be used in most charts
genotypec<-c("#DC9C03","#965499","#00A087FF")
treatmentc<-c("#4DBBD5FF", "#E64B35FF")

otu1 <- otu_table(ASV, taxa_are_rows=TRUE)

sam1 <- sample_data(metadata) 

tax1 <- tax_table(taxonomy)
## tax1 <- data.table::as.data.table(taxonomy)

#make sure rows of taxa match the rows of OTU table
colnames(tax1) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus")
rownames(tax1)<- rownames(otu1)

#build a class data phyloseq object (named FLNA) that contains OTU table, taxonomy info and metadata info for all samples
FLNA <- phyloseq(otu1,tax1,sam1)

#check data to see if all data makes sense/readable
most_abundant_taxa <- sort(taxa_sums(FLNA), TRUE)[1:25]
ex2 <- prune_taxa(names(most_abundant_taxa), FLNA)
topFamilies <- tax_table(ex2)[, "Family"]
as(topFamilies, "vector")

#######################################################################################
#As a first analysis, we will look at the distribution of read counts from our samples
sample_sum_df <- data.frame(sum = sample_sums(FLNA))
ggplot(sample_sum_df, aes(x = sum)) + 
  geom_histogram(color = "black", fill = "indianred", binwidth = 2500) +
  ggtitle("Distribution of sample sequencing depth") + 
  xlab("Read counts") +
  theme(axis.title.y = element_blank())

# output mean, max and min of sample read counts
smin <- min(sample_sums(FLNA))
smin
smean <- mean(sample_sums(FLNA))
smean
ssd<- sd(sample_sums(FLNA))
ssd
smax <- max(sample_sums(FLNA))
smax

#calculate coverage for each sample - note: you can do this before or after rarefying the data
phyloseq_coverage(FLNA,correct_singletons = FALSE, add_attr = T)


#rarefy data (if desired; if not jump this step and normalize the easy way(below))

FLNA_rare<-rarefy_even_depth(FLNA, sample.size = min(sample_sums(FLNA)),
                                 rngseed = 771, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)

# Convert phyloseq to a data frame
df <- psmelt(FLNA_rare)

# Calculate min, max, mean, and sd
min_rare <- min(df$Abundance, na.rm = TRUE)
max_rare <- max(df$Abundance, na.rm = TRUE)
mean_rare <- mean(df$Abundance, na.rm = TRUE)
sd_rare <- sd(df$Abundance, na.rm = TRUE)

# Print the results
print(paste("Min: ", min_rare))
print(paste("Max: ", max_rare))
print(paste("Mean: ", mean_rare))
print(paste("SD: ", sd_rare))

#quick look into the new OTU table
otu_FLNA_rare<-otu_table(FLNA_rare)
otu_FLNA_rare

write.csv(otu_FLNA_rare, file = "otu_FLNA_rare.csv")

#fix the order at which samples appear on charts: first control. etc..
sample_data(FLNA_rare)$genotype<- factor(sample_data(FLNA_rare)$genotype,levels = c("wt", "Q","R"))
sample_data(FLNA_rare)$treatment<- factor(sample_data(FLNA_rare)$treatment, levels=c("antibiotics", "none"))

#subset samples by genotype
FLNA_rare_wt<-subset_samples(FLNA_rare, genotype=="wt")
FLNA_rare_Q<-subset_samples(FLNA_rare, genotype=="Q")
FLNA_rare_R<-subset_samples(FLNA_rare, genotype=="R")

# Subset samples by genotype and treatment
FLNA_rare_wt_none <- subset_samples(FLNA_rare, genotype == "wt" & treatment == "none")
FLNA_rare_Q_none <- subset_samples(FLNA_rare, genotype == "Q" & treatment == "none")
FLNA_rare_R_none <- subset_samples(FLNA_rare, genotype == "R" & treatment == "none")

# Combine the subset samples into a phyloseq object
FLNA_rare_none <- merge_phyloseq(FLNA_rare_wt_none, FLNA_rare_Q_none, FLNA_rare_R_none)

######## Create alpha diversity plot for untreated samples
alpha_none <- plot_richness(FLNA_rare_none, x = "genotype", color = "genotype", measures = c("InvSimpson", "Shannon", "Observed")) + scale_color_manual(values = genotypec)
alpha_none + geom_boxplot(lwd = 1, width = 0.5) + geom_point() + stat_compare_means(method = 'wilcox', label = 'p.signif', bracket.size = 0.1, label.y.npc = "top", tip.length = 0.01, hide.ns = TRUE, label.x.npc = "center") + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9))

# Subset samples by genotype and treatment
FLNA_rare_Q_abx <- subset_samples(FLNA_rare, genotype == "Q" & treatment == "antibiotics")
FLNA_rare_R_abx <- subset_samples(FLNA_rare, genotype == "R" & treatment == "antibiotics")

# Combine the subset samples into a phyloseq object
FLNA_rare_abx <- merge_phyloseq(FLNA_rare_Q_abx, FLNA_rare_R_abx)

####### Create alpha diversity plot for antibiotic treated samples

genotypec<-c("#965499","#00A087FF")

alpha_abx <- plot_richness(FLNA_rare_abx, x = "genotype", color = "genotype", measures = c("InvSimpson", "Shannon", "Observed")) + scale_color_manual(values = genotypec)
alpha_abx + geom_boxplot(lwd = 1, width = 0.5) + geom_point() + stat_compare_means(method = 'wilcox', label = 'p.signif', bracket.size = 0.1, label.y.npc = "top", tip.length = 0.01, hide.ns = TRUE, label.x.npc = "center") + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9))

#extract important data for a new alpha div simple table - genotype (untreated)
richness.genotype.rare<- cbind(estimate_richness(FLNA_rare_none, 
                                                 measures = c('Observed','Shannon',"Simpson","InvSimpson",'genotype')),
                               sample_data(FLNA_rare_none)$genotype)

write.csv(richness.genotype.rare, file = "alphaDivFLNA_geno.csv")


colnames(richness.genotype.rare) <- c('Observed','Shannon',"Simpson","InvSimpson",'genotype')
richness.genotype.rare

#additional stats: note if data is normally distributed, you can use a t-test instead of wilcox - genotype (untreated)
compare_means(Shannon ~ genotype,  data = richness.genotype.rare,
              ref.group = "wt", method = "wilcox")
compare_means(InvSimpson ~ genotype,  data = richness.genotype.rare,
              ref.group = "wt", method = "wilcox")
compare_means(Observed ~ genotype,  data = richness.genotype.rare,
              ref.group = "wt", method = "wilcox")

###### the above is good to compare FLNQ and FLNR to wt, but if we want to compare FLNQ to FLNR use below stats
###### but only look at the ones with R as group 2 as already have data comparing wt
compare_means(Shannon ~ genotype,  data = richness.genotype.rare,
              ref.group = "Q", method = "wilcox")
compare_means(InvSimpson ~ genotype,  data = richness.genotype.rare,
              ref.group = "Q", method = "wilcox")
compare_means(Observed ~ genotype,  data = richness.genotype.rare,
              ref.group = "Q", method = "wilcox")

#extract important data for a new alpha div simple table - treatment (Q/R)
richness.treatedgenotype.rare<- cbind(estimate_richness(FLNA_rare_abx, 
                                        measures = c('Observed','Shannon',"Simpson","InvSimpson",'genotype')),
                      sample_data(FLNA_rare_abx)$genotype)

write.csv(richness.treatedgenotype.rare, file = "alphaDivFLNA_treatment.csv")

colnames(richness.treatedgenotype.rare) <- c('Observed','Shannon',"Simpson","InvSimpson",'genotype')
richness.treatedgenotype.rare

#additional stats: note if data is normally distributed, you can use a t-test instead of wilcox - treat 
compare_means(Shannon ~ genotype,  data = richness.treatedgenotype.rare,
              ref.group = "Q", method = "wilcox")
compare_means(InvSimpson ~ genotype,  data = richness.treatedgenotype.rare,
              ref.group = "Q", method = "wilcox")
compare_means(Observed ~ genotype,  data = richness.treatedgenotype.rare,
              ref.group = "Q", method = "wilcox")

###################################Calculate rel abundance and filter out some data###########################

#calculate relative abundance 

FLNA_rare_Rel  = transform_sample_counts(FLNA_rare, function(x) x / sum(x) )

ASV_Rel = as(otu_table(FLNA_rare_Rel), "matrix")
sample_data(FLNA_rare_Rel)$ID<- factor(sample_data(FLNA_rare_Rel)$Mouse_ID,levels = c("4003","4004","4038","4039","4037","4040","4044",
                               "4046","4022","4020","4021","4023","7659","7703","7643","7647","7653","7657","7656","7654","7658",
                               "7704","7644","7646","7701","7700","7702","7757","7623","7622","7650","7651","7707","7706","7758",
                               "7636","7633","7769","7767","7766","7649","7648"))

######################################STACKED BARPLOTS##########################################################

# Find out what our top 10 genus and family are 

genus_count <- table(taxonomy$Genus)
sorted_genus <- sort(genus_count, decreasing = TRUE)
top_10_genus <- head(sorted_genus, 10)
print(top_10_genus)

family_count <- table(taxonomy$Family)
sorted_family <- sort(family_count, decreasing = TRUE)
top_10_family <- head(sorted_family, 10)
print(top_10_family)

######################################Family level########################################
FLNA_Rel_family <- tax_glom(FLNA_rare_Rel, taxrank = 'Family',NArm=FALSE) # agglomerate taxa
FamilyRel <- psmelt(FLNA_Rel_family ) # create dataframe from phyloseq object
FamilyRel$Family <- as.character(FamilyRel$Family) #convert to character
FamilyRel$Family[FamilyRel$Abundance < 0.05] <- "Other" #rename taxa with < 5% abundance

ggplot(FamilyRel, aes(x =ID, y = Abundance, fill = Family))+geom_bar(stat="identity") +ylab("Relative abundance")+xlab("")+theme_bw()+theme(legend.text = element_text(size=8),legend.title = element_text(size=10),axis.title.y=element_text(size=12),axis.title.x = element_text(size=12),axis.text.y = element_text(size=8),axis.text.x = element_text(angle = 45, hjust=1,size=7),strip.text = element_text(size = 10))+facet_wrap(genotype ~ treatment,scales= "free_x",nrow=1)

#stats to analyse
analyze_family <- function(Family) {
  avg_abundance <- aggregate(Abundance ~ Family + genotype + treatment, FamilyRel, mean)
  avg_abundance[avg_abundance$Family == Family, ]
}

families_present <- unique(FamilyRel$Family)
results_family <- lapply(families_present, analyze_family)
results_family

results_family_df <- do.call(rbind, results_family)
write.csv(results_family_df, file = "results_family.csv")


##################################Phylum Level###################################

FLNA_Rel_phyla <- tax_glom(FLNA_rare_Rel, taxrank = 'Phylum') # agglomerate taxa
PhylaRel<- psmelt(FLNA_Rel_phyla) # create dataframe from phyloseq object
PhylaRel$Phylum<- as.character(PhylaRel$Phylum) #convert to character
PhylaRel$Phylum[PhylaRel$Abundance < 0.01] <- "Other" #rename taxa with < 1% abundance

ggplot(PhylaRel, aes(x =ID, y = Abundance, fill = Phylum))+geom_bar(stat="identity") +ylab("Relative abundance")+xlab("")+theme_bw()+theme(legend.text = element_text(size=8),legend.title = element_text(size=10),axis.title.y=element_text(size=12),axis.title.x = element_text(size=12),axis.text.y = element_text(size=8),axis.text.x = element_text(angle = 45, hjust=1,size=7),strip.text = element_text(size = 10))+facet_wrap(genotype ~ treatment,scales= "free_x",nrow=1)#+scale_fill_manual(values=colours_htot_phyla)

#stats to analyse
analyze_phylum <- function(Phylum) {
  avg_abundance <- aggregate(Abundance ~ Phylum + genotype + treatment, PhylaRel, mean)
  avg_abundance[avg_abundance$Phylum == Phylum, ]
}

phylum_present <- unique(PhylaRel$Phylum)
results_phylum <- lapply(phylum_present, analyze_phylum)
results_phylum

results_phylum_df <- do.call(rbind, results_phylum)
write.csv(results_phylum_df, file = "results_phylum.csv")

######################################Genus level########################################
FLNA_Rel_genus <- tax_glom(FLNA_rare_Rel, taxrank = 'Genus',NArm=FALSE) # agglomerate taxa
GenusRel <- psmelt(FLNA_Rel_genus ) # create dataframe from phyloseq object
GenusRel$Genus <- as.character(GenusRel$Genus) #convert to character
GenusRel$Genus[GenusRel$Abundance < 0.05] <- "Other" #rename taxa with < 5% abundance

ggplot(GenusRel, aes(x =ID, y = Abundance, fill = Genus))+geom_bar(stat="identity") +ylab("Relative abundance")+xlab("")+theme_bw()+theme(legend.text = element_text(size=8),legend.title = element_text(size=10),axis.title.y=element_text(size=12),axis.title.x = element_text(size=12),axis.text.y = element_text(size=8),axis.text.x = element_text(angle = 45, hjust=1,size=7),strip.text = element_text(size = 10))+facet_wrap(genotype ~ treatment,scales= "free_x",nrow=1)#+scale_fill_manual(values=colours_htot)

#stats to analyse
analyze_genus <- function(Genus) {
  avg_abundance <- aggregate(Abundance ~ Genus + genotype + treatment, GenusRel, mean)
  avg_abundance[avg_abundance$Genus == Genus, ]
}

genus_present <- unique(GenusRel$Genus)
results_genus <- lapply(genus_present, analyze_genus)
results_genus

results_genus_df <- do.call(rbind, results_genus)
write.csv(results_genus_df, file = "results_genus.csv")

####################################BETA DIVERSITY ANALYSIS#########################################################


FLNA_rare_none <- subset_samples(FLNA_rare, treatment == "none")
bc_dist_none <- phyloseq::distance(FLNA_rare_none, method="bray", weighted=F)
results_none <- adonis2(bc_dist_none ~ sample_data(FLNA_rare_none)$genotype)
bc_dist_none
results_none

bc.nmds_none <- metaMDS(bc_dist_none)
stressplot(bc.nmds_none)

write.csv(bc_dist_none, file = "permanova_none")

# For NMDS plot on genotype with "none" in treatment
geno.scrs <- as.data.frame(scores(bc.nmds, display = "sites"))
geno.scrs <- cbind(geno.scrs, sites = ifelse(ord.meta$treatment == "none", ord.meta$genotype, NA))
geno.scrs <- geno.scrs[!is.na(geno.scrs$sites),]
head(geno.scrs)


NMDS.mean = aggregate(geno.scrs[,c("NMDS1", "NMDS2")], list(group = geno.scrs$sites), mean)
NMDS.mean


veganCovEllipse <- function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

df_ell.dune.treatment <- data.frame() #sets up a data frame before running the function.
for(g in levels(geno.scrs$sites)){
  df_ell.dune.treatment <- rbind(df_ell.dune.treatment, cbind(as.data.frame(with(site.scrs [site.scrs$site==g,],
                                                                                 veganCovEllipse(cov.wt(cbind(NMDS1,NMDS2),wt=rep(1/length(NMDS1),length(NMDS1)))$cov,center=c(mean(NMDS1),mean(NMDS2))))) ,site=g))
}

geno.scrs$sites <- as.factor(geno.scrs$sites)
sites <- c("#DC9C03","#965499","#00A087FF")

#draw the NMDS ordination plot with elipses for each dataset
NMDS_plot<-ggplot(data = geno.scrs, aes(NMDS1, NMDS2)) + geom_point(aes(color = sites),show.legend=NA,size=1.5)+scale_color_manual(values=sites)
NMDS_plot

NMDS_plot_elipse <- ggplot(data = geno.scrs, aes(NMDS1, NMDS2)) +
  geom_point(aes(color = sites), show.legend=NA, size=1.5) +
  scale_color_manual(values=sites, labels =c("wt","Q", "R")) +
  stat_ellipse(data = geno.scrs, aes(NMDS1, NMDS2, color = sites), type = "t", level = 0.95)
NMDS_plot_elipse



##########################
#########################

# For NMDS plot on genotypes Q and R but antibiotic treated

FLNA_rare_treat <- subset_samples(FLNA_rare, treatment == "antibiotics")
bc_dist_treat <- phyloseq::distance(FLNA_rare_treat, method="bray", weighted=F)
results_treat <- adonis2(bc_dist_treat ~ sample_data(FLNA_rare_treat)$genotype)
bc_dist_treat
results_treat

bc.nmds_treat <- metaMDS(bc_dist_treat)
stressplot(bc.nmds_treat)

write.csv(bc_dist_treat, file = "permanova_treat")

geno.scrs <- as.data.frame(scores(bc.nmds, display = "sites"))
geno.scrs <- cbind(geno.scrs, sites = ifelse(ord.meta$treatment == "antibiotics", ord.meta$genotype, NA))
geno.scrs <- geno.scrs[!is.na(geno.scrs$sites),]
head(geno.scrs)

NMDS.mean = aggregate(geno.scrs[,c("NMDS1", "NMDS2")], list(group = geno.scrs$sites), mean)
NMDS.mean

veganCovEllipse <- function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

df_ell.dune.treatment <- data.frame() #sets up a data frame before running the function.
for(g in levels(geno.scrs$sites)){
  df_ell.dune.treatment <- rbind(df_ell.dune.treatment, cbind(as.data.frame(with(site.scrs [site.scrs$site==g,],
                                                                                 veganCovEllipse(cov.wt(cbind(NMDS1,NMDS2),wt=rep(1/length(NMDS1),length(NMDS1)))$cov,center=c(mean(NMDS1),mean(NMDS2))))) ,site=g))
}

geno.scrs$sites <- as.factor(geno.scrs$sites)
sites <- c("#965499","#00A087FF")

#draw the NMDS ordination plot with elipses for each dataset
NMDS_plot<-ggplot(data = geno.scrs, aes(NMDS1, NMDS2)) + geom_point(aes(color = sites),show.legend=NA,size=1.5)+scale_color_manual(values=sites)
NMDS_plot

NMDS_plot_elipse <- ggplot(data = geno.scrs, aes(NMDS1, NMDS2)) +
  geom_point(aes(color = sites), show.legend=NA, size=1.5) +
  scale_color_manual(values=sites, labels =c("Q", "R")) +
  stat_ellipse(data = geno.scrs, aes(NMDS1, NMDS2, color = sites), type = "t", level = 0.95)
NMDS_plot_elipse



###########DESEQ2 anaysis#################


###################Comparison treatments: Q vs R in None###################

FLNA_prune = prune_taxa(taxa_sums(FLNA) > 10, FLNA)

FLNA_prune_treat<- subset_samples(FLNA_prune, treatment %in% c("none")) 

sample_data(FLNA_prune_treat)$genotype<- factor(sample_data(FLNA_prune_treat)$genotype,levels = c("Q", "wt","R"))

treatdss = phyloseq_to_deseq2(FLNA_prune_treat, ~genotype)

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(treatdss), 1, gm_mean)
treatdss = estimateSizeFactors(treatdss, geoMeans = geoMeans)
treatdss = DESeq(treatdss, fitType="local")

res = results(treatdss)
res = res[order(res$padj, na.last=NA), ]
alpha = 0.05
sigtab = res[(res$padj < alpha), ]
sigtab = res[(res$padj < alpha) & (abs(res$log2FoldChange) > 2), ]
sigtab
sigtab = cbind(as(sigtab, "data.frame"), 
               as(tax_table(FLNA_prune)[rownames(sigtab), ], "matrix"))
sigtabgen = subset(sigtab, !is.na(Genus))

write.csv(sigtab, file = "DeSeq2_Q_R_none.csv")

x = tapply(sigtabgen$log2FoldChange, sigtabgen$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels=names(x))

x = tapply(sigtabgen$log2FoldChange, sigtabgen$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Genus = factor(as.character(sigtabgen$Genus), levels=names(x))

x = tapply(sigtabgen$log2FoldChange, sigtabgen$Order, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Order = factor(as.character(sigtabgen$Order), levels=names(x))



qr1<-ggplot(sigtabgen, aes(y=Genus, x=log2FoldChange, color=Family)) + 
  geom_vline(xintercept = 0.0, color = "gray", linewidth = 0.5) +
  geom_point(size=4) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+theme_bw()+scale_color_d3(palette="category20")
qr1


qr2<-ggplot(sigtabgen, aes(y=Genus, x=log2FoldChange, color=Phylum)) + 
  geom_vline(xintercept = 0.0, color = "gray", linewidth = 0.5) +
  geom_point(size=4) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+theme_bw()+scale_color_d3(palette="category20")
qr2

###################Comparison treatments: wt vs R in None###################

FLNA_prune = prune_taxa(taxa_sums(FLNA) > 10, FLNA)

FLNA_prune_treat<- subset_samples(FLNA_prune, treatment %in% c("none")) 

sample_data(FLNA_prune_treat)$genotype<- factor(sample_data(FLNA_prune_treat)$genotype,levels = c("wt", "Q","R"))

treatdss = phyloseq_to_deseq2(FLNA_prune_treat, ~genotype)

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(treatdss), 1, gm_mean)
treatdss = estimateSizeFactors(treatdss, geoMeans = geoMeans)
treatdss = DESeq(treatdss, fitType="local")

res = results(treatdss)
res = res[order(res$padj, na.last=NA), ]
alpha = 0.05
sigtab = res[(res$padj < alpha), ]
sigtab = res[(res$padj < alpha) & (abs(res$log2FoldChange) > 2), ]
sigtab
sigtab = cbind(as(sigtab, "data.frame"), 
               as(tax_table(FLNA_prune)[rownames(sigtab), ], "matrix"))
sigtabgen = subset(sigtab, !is.na(Genus))

write.csv(sigtab, file = "DeSeq2_WT_R_None.csv")

x = tapply(sigtabgen$log2FoldChange, sigtabgen$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels=names(x))

x = tapply(sigtabgen$log2FoldChange, sigtabgen$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Genus = factor(as.character(sigtabgen$Genus), levels=names(x))

x = tapply(sigtabgen$log2FoldChange, sigtabgen$Order, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Order = factor(as.character(sigtabgen$Order), levels=names(x))



wtr1<-ggplot(sigtabgen, aes(y=Genus, x=log2FoldChange, color=Family)) + 
  geom_vline(xintercept = 0.0, color = "gray", linewidth = 0.5) +
  geom_point(size=4) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+theme_bw()+scale_color_d3(palette="category20")
wtr1

wtr2<-ggplot(sigtabgen, aes(y=Genus, x=log2FoldChange, color=Phylum)) + 
  geom_vline(xintercept = 0.0, color = "gray", linewdith = 0.5) +
  geom_point(size=4) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+theme_bw()+scale_color_d3(palette="category20")
wtr2

###################Comparison treatments: wt vs Q in None###################

FLNA_prune = prune_taxa(taxa_sums(FLNA) > 10, FLNA)

FLNA_prune_treat<- subset_samples(FLNA_prune, treatment %in% c("none")) 

sample_data(FLNA_prune_treat)$genotype<- factor(sample_data(FLNA_prune_treat)$genotype,levels = c("wt", "R","Q"))

treatdss = phyloseq_to_deseq2(FLNA_prune_treat, ~genotype)

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(treatdss), 1, gm_mean)
treatdss = estimateSizeFactors(treatdss, geoMeans = geoMeans)
treatdss = DESeq(treatdss, fitType="local")

res = results(treatdss)
res = res[order(res$padj, na.last=NA), ]
alpha = 0.05
sigtab = res[(res$padj < alpha), ]
sigtab = res[(res$padj < alpha) & (abs(res$log2FoldChange) > 2), ]
sigtab
sigtab = cbind(as(sigtab, "data.frame"), 
               as(tax_table(FLNA_prune)[rownames(sigtab), ], "matrix"))
sigtabgen = subset(sigtab, !is.na(Genus))

write.csv(sigtab, file = "DeSeq2_WT_Q_None.csv")

x = tapply(sigtabgen$log2FoldChange, sigtabgen$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels=names(x))

x = tapply(sigtabgen$log2FoldChange, sigtabgen$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Genus = factor(as.character(sigtabgen$Genus), levels=names(x))

x = tapply(sigtabgen$log2FoldChange, sigtabgen$Order, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Order = factor(as.character(sigtabgen$Order), levels=names(x))


wtq1<-ggplot(sigtabgen, aes(y=Genus, x=log2FoldChange, color=Family)) + 
  geom_vline(xintercept = 0.0, color = "gray", linewidth = 0.5) +
  geom_point(size=4) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+theme_bw()+scale_color_d3(palette="category20")
wtq1

wtq2<-ggplot(sigtabgen, aes(y=Genus, x=log2FoldChange, color=Phylum)) + 
  geom_vline(xintercept = 0.0, color = "gray", linewidth = 0.5) +
  geom_point(size=4) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+theme_bw()+scale_color_d3(palette="category20")
wtq2

###################Comparison treatments: Q vs R in Antibiotics###################
################### don't include figures in diss but include data in main text ###########

FLNA_prune = prune_taxa(taxa_sums(FLNA) > 10, FLNA)

FLNA_prune_treat<- subset_samples(FLNA_prune, treatment %in% c("antibiotics")) 

sample_data(FLNA_prune_treat)$genotype<- factor(sample_data(FLNA_prune_treat)$genotype,levels = c("Q", "wt","R"))

treatdss = phyloseq_to_deseq2(FLNA_prune_treat, ~genotype)

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(treatdss), 1, gm_mean)
treatdss = estimateSizeFactors(treatdss, geoMeans = geoMeans)
treatdss = DESeq(treatdss, fitType="local")

res = results(treatdss)
res = res[order(res$padj, na.last=NA), ]
alpha = 0.05
sigtab = res[(res$padj < alpha), ]
sigtab = res[(res$padj < alpha) & (abs(res$log2FoldChange) > 2), ]
sigtab
sigtab = cbind(as(sigtab, "data.frame"), 
               as(tax_table(FLNA_prune)[rownames(sigtab), ], "matrix"))
sigtabgen = subset(sigtab, !is.na(Genus))

write.csv(sigtab, file = "DeSeq2_Q_R_Treat.csv")

x = tapply(sigtabgen$log2FoldChange, sigtabgen$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels=names(x))

x = tapply(sigtabgen$log2FoldChange, sigtabgen$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Genus = factor(as.character(sigtabgen$Genus), levels=names(x))

x = tapply(sigtabgen$log2FoldChange, sigtabgen$Order, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Order = factor(as.character(sigtabgen$Order), levels=names(x))


color_vector <- rainbow(2)

qra1<-ggplot(sigtabgen, aes(y=Genus, x=log2FoldChange, color=Family)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(size=4) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+theme_bw()+scale_color_d3(palette="category20")
qra1

qra2<-ggplot(sigtabgen, aes(y=Genus, x=log2FoldChange, color=Phylum)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(size=4) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))+theme_bw()+scale_color_d3(palette="category20")
qra2



