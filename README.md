# Soldier-microbiome
This repository contains the R code used for the analysis and visualisation of data generated in the Gurkha-microbiome study.

```R
#### Analyses and Final Plots ####

############################################# START ##################################################

#### QC plots START ####

library(ggplot2)

# This code loads in the required data and assigns them to variables. 

insert <- read.csv("insert_stats.tsv", sep='\t')
reads <- read.csv("read_counts.tsv", sep="\t")
length <- read.csv("read_length_stats.tsv", sep="\t", row.names = 1)

gen <- read.csv("../genus_transposed_3g.csv", stringsAsFactors = T) # Used for group info
group <- gen[,c(1,94)]
reads$Group <- with(group, Group[match(reads$Sample,Genus)]) # Assigns groups to sample IDs.
reads <- reads[!(reads$Sample %in% c("Undetermined")), ] # Ignores reads without barcode.

# This code creates a plot of the total bases

reads <- reads[grep("QC", reads$Step),]
ggplot(reads, aes(x=Sample, y=Total_Bases, fill=Group)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=c("Control" = "#57bda3", "Gurkha" = "#6aa2e9","Non-Gurkha" = "#d99a46")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 5700000000)) +
  theme_minimal() + ylab("Total Bases") +
  theme(panel.grid=element_blank(), axis.line=element_line(colour="black"),
        axis.title = element_text(family="montserrat"),
        axis.text.y = element_text(family="montserrat"),
        axis.text.x = element_text(family="montserrat", angle=90),
        legend.text = element_text(family="montserrat"),
        legend.title = element_blank())
ggsave("total_bases.png", device="png", dpi=300)

#### QC plots END ####

#### Bubble plot of genus abundance START ####

# This code loads the ggplot2, reshape2, and showtext libraries into R, which are used for data visualization and font handling.

library(ggplot2)
library(reshape2)
library(showtext)

# This code adds the Montserrat font from Google Fonts to the showtext font library, and sets showtext to automatically use this font for plotting.

font_add_google(name="montserrat", family="montserrat")
showtext_auto()

# This code reads in a CSV file named "genus_transposed_3g.csv" as a dataframe called "gen", with the row names set to the first column of the CSV file. The "stringsAsFactors" argument is set to TRUE, which will convert any character columns in the CSV file to factor variables in the dataframe.

gen <- read.csv("genus_transposed_3g.csv", stringsAsFactors = T, row.names = 1)

# This code creates a character vector called "row.names.remove", which contains the row names that should be removed from the "gen" dataframe.

row.names.remove <- c("NC1", "OC", "DES", "PC", "PC1")

# This code removes the rows with the names listed in "row.names.remove" from the "gen" dataframe using subsetting notation.

gen <- gen[!row.names(gen) %in% row.names.remove,]

# This code removes columns 93 to 95 from the "gen" dataframe using subsetting notation and stores the result in a new dataframe called "gen2".

gen2 <- gen[,-c(93:95)]

# This code removes any columns in "gen2" where all values are less than 1, using the "keep" function from the dplyr package with a lambda function. The resulting dataframe is stored in "gen3".

gen3 <- gen2 %>% 
  keep(~ any(. >= 1))

# This code rearranges the columns in "gen3" to be in alphabetical order, adds the "Group" and "Sample" columns from the original "gen" dataframe, and replaces underscores in the column names with spaces. The resulting dataframe is stored in "gen4".

gen4 <- gen3[,order(colnames(gen3))]
gen4$Group <- gen$Group
gen4$Sample <- rownames(gen)
colnames(gen4) <- gsub("_"," ", colnames(gen4))

# This code uses the "melt" function from the reshape2 package to convert the "gen4" dataframe from wide format to long format, using the "Sample" and "Group" columns as the IDs.

gen5 <- melt(gen4, id = c("Sample", "Group"))

# This code uses ggplot2 to create a bubble plot from the gen5 dataframe.

bubble_plot <- ggplot(gen5, aes(x = Sample, y = variable)) + 
  geom_point(aes(size = value, fill = Group), alpha = 1, shape = 21) + 
  scale_size_continuous(limits = c(0.000001, 100), range = c(1,10), breaks = c(1,10,50,75)) + 
  labs( x= "", y = "", size = "Relative Abundance (%)", fill = "Group")  + 
  theme(legend.key=element_blank(), 
        axis.text.x = element_text(colour = "black", size = 12, angle = 45, hjust = 1, family="montserrat"), 
        axis.text.y = element_text(colour = "black", size = 13, family="montserrat", face = "italic"), 
        legend.text = element_text(size = 13, colour ="black", family="montserrat"), 
        legend.title = element_text(size = 13, family="montserrat", face = "bold"), panel.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
        legend.position = "right", panel.grid.major.y = element_line(colour = "grey95")) +  
  scale_fill_manual(values = c("#6aa2e9", "#d99a46"), guide = guide_legend(override.aes = list(size=5))) + 
  scale_y_discrete(limits = rev(levels(gen5$variable))) 

# This code saves the bubble plot to a .tiff file at a resolution of 300 dpi. 

tiff(filename = "Genus_abundance_bubble_plot.tiff", width = 14, height = 8, units = "in", res = 300)

bubble_plot

dev.off()

#### Bubble plot of genus abundance END ####

######################################################################################################

#### Alpha diversity START ####

# This code loads three libraries: "ggplot2" for data visualization, "vegan" for ecological analysis, and "ggpubr" for publication-ready plots.

library(ggplot2)
library(vegan)
library(ggpubr)

#Species diversity (shannon)

dft <- read.csv("species_transposed_3g.csv", row.names = 1) # This reads a CSV file named "species_transposed_3g.csv" and assigns its contents to a variable named dft.The argument "species_transposed_3g.csv" specifies the file name, while row.names = 1 indicates that the first column of the CSV file should be used as row names.
row.names.remove <- c("NC1", "OC", "DES") # This code creates a character vector called "row.names.remove", which contains the row names that should be removed from the "dft" dataframe.
dft <- dft[!row.names(dft) %in% row.names.remove,] # This code removes the rows with the names listed in "row.names.remove" from the "dft" dataframe using subsetting notation.
species_df <- dft[,-c(233:235)] # Removes descriptive variables

shannon <- as.data.frame(diversity(species_df)) # Calculate shannon diversity of each sample
shannon$Group <- dft$Group # Append the descriptive variables back on to the dataframe
shannon$Group <- as.factor(shannon$Group) # Convert the Group variable column to factors
names(shannon)[1] <- "shannon_diversity" # Rename the diversity value column
No_control <- shannon[(shannon$Group == "Gurkha" | shannon$Group == "Non-Gurkha"), ] # Removes the two control samples

# Plot dot and box plot of data 

shann_spec_plot <- ggboxplot(No_control, x="Group", y="shannon_diversity", add="jitter", color = "Group", palette=c("#6aa2e9", "#d99a46"))
my_comparisons <- list(c("Gurkha","Non-Gurkha"))
shann_spec_plot_2 <- shann_spec_plot + stat_compare_means(comparisons = my_comparisons) +
  theme_minimal() +
  theme(axis.title = element_text(family="montserrat", size = 12),
        axis.title.x = element_blank(),
        axis.text = element_text(family="montserrat", size = 12),
        legend.position ="none") +
  labs(y="Shannon Diversity")

# Species diversity (simpson)

spec <- read.csv("species_transposed_3g.csv", row.names = 1) # This reads a CSV file named "species_transposed_3g.csv" and assigns its contents to a variable named spec.The argument "species_transposed_3g.csv" specifies the file name, while row.names = 1 indicates that the first column of the CSV file should be used as row names.
row.names.remove <- c("NC1", "OC", "DES") # This code creates a character vector called "row.names.remove", which contains the row names that should be removed from the "spec" dataframe.
spec <- spec[!row.names(spec) %in% row.names.remove,] # This code removes the rows with the names listed in "row.names.remove" from the "spec" dataframe using subsetting notation.
spec_df <- spec[,-c(233:235)] # Remove descriptive variables

simpson <- as.data.frame(diversity(spec_df, index = "simpson")) # Calculate simpson diversity of each sample
simpson$Group <- spec$Group # Append the descriptive variables back on to the dataframe
simpson$Group <- as.factor(simpson$Group) # Convert the Group variable column to factors
names(simpson)[1] <- "simpson_diversity" # Rename the diversity value column
No_control_simpson <- simpson[(simpson$Group == "Gurkha" | simpson$Group == "Non-Gurkha"), ] # Removes the two control samples

# Plot dot and box plot of data 

simp_spec_plot <- ggboxplot(No_control_simpson, x="Group", y="simpson_diversity", add="jitter", color = "Group", palette=c("#6aa2e9", "#d99a46"))
my_comparisons <- list(c("Gurkha","Non-Gurkha"))
simp_spec_plot_2 <- simp_spec_plot + stat_compare_means(comparisons = my_comparisons) +
  theme_minimal() +
  theme(axis.title = element_text(family="montserrat", size = 12),
        axis.title.x = element_blank(),
        axis.text = element_text(family="montserrat", size = 12),
        legend.position ="none") +
  labs(y="Simpson Diversity")

# Genus diversity (Shannon)

dft2 <- read.csv("genus_transposed_3g.csv", row.names = 1) # This reads a CSV file named "genus_transposed_3g.csv" and assigns its contents to a variable named dft2.The argument "genus_transposed_3g.csv" specifies the file name, while row.names = 1 indicates that the first column of the CSV file should be used as row names.
row.names.remove <- c("NC1", "OC", "DES") # This code creates a character vector called "row.names.remove", which contains the row names that should be removed from the "dft2" dataframe.
dft2 <- dft2[!row.names(dft2) %in% row.names.remove,] # This code removes the rows with the names listed in "row.names.remove" from the "dft2" dataframe using subsetting notation.
genus_df <- dft2[,-c(93,94,95)]# Remove descriptive variables

shannon_g <- as.data.frame(diversity(genus_df)) # Calculate shannon diversity of each sample
shannon_g$Group <- dft2$Group # Append the descriptive variables back on to the dataframe
shannon_g$Group <- as.factor(shannon_g$Group) # Convert the Group variable column to factors
names(shannon_g)[1] <- "shannon_diversity" # Rename the diversity value column
No_control_g <- shannon_g[(shannon_g$Group == "Gurkha" | shannon_g$Group == "Non-Gurkha"), ] # Removes the two control samples

# Plot dot and box plot of data

shann_gen_plot <- ggboxplot(No_control_g, x="Group", y="shannon_diversity", add="jitter", color = "Group", palette=c("#6aa2e9", "#d99a46"))
my_comparisons <- list(c("Gurkha","Non-Gurkha"))
shann_gen_plot_2 <- shann_gen_plot + stat_compare_means(comparisons = my_comparisons) +
  theme_minimal() +
  theme(axis.title = element_text(family="montserrat", size = 12),
        axis.title.x = element_blank(),
        axis.text = element_text(family="montserrat", size = 12),
        legend.position ="none") +
  labs(y="Shannon Diversity")

# Genus diversity (simpson)

gen_2 <- read.csv("genus_transposed_3g.csv", row.names = 1) # This reads a CSV file named "genus_transposed_3g.csv" and assigns its contents to a variable named gen_2.The argument "genus_transposed_3g.csv" specifies the file name, while row.names = 1 indicates that the first column of the CSV file should be used as row names.
row.names.remove <- c("NC1", "OC", "DES") # This code creates a character vector called "row.names.remove", which contains the row names that should be removed from the "gen2" dataframe.
gen_2 <- gen_2[!row.names(gen_2) %in% row.names.remove,] # This code removes the rows with the names listed in "row.names.remove" from the "gen_2" dataframe using subsetting notation.
gen_2_df <- gen_2[,-c(93,94,95)]# Remove descriptive variables

simpson_g <- as.data.frame(diversity(gen_2_df, index = "simpson")) # Calculate simpson diversity of each sample
simpson_g$Group <- gen_2$Group # Append the descriptive variables back on to the dataframe
simpson_g$Group <- as.factor(simpson_g$Group) # Convert the Group variable column to factors
names(simpson_g)[1] <- "simpson_diversity" # Rename the diversity value column
No_control_simpson_g <- simpson_g[(simpson_g$Group == "Gurkha" | simpson_g$Group == "Non-Gurkha"), ] # Removes the two control samples

# Plot dot and box plot of data 

simp_gen_plot <- ggboxplot(No_control_simpson_g, x="Group", y="simpson_diversity", add="jitter", color = "Group", palette=c("#6aa2e9", "#d99a46"))
my_comparisons <- list(c("Gurkha","Non-Gurkha"))
simp_gen_plot_2 <- simp_gen_plot + stat_compare_means(comparisons = my_comparisons) +
  theme_minimal() +
  theme(axis.title = element_text(family="montserrat", size = 12),
        axis.title.x = element_blank(),
        axis.text = element_text(family="montserrat", size = 12),
        legend.position ="none") +
  labs(y="Simpson Diversity")

### The following code creates multipanel diversity plots.

### Shannon species & genus

tiff(filename = "Shannon_diversity_SG.tiff", width = 8, height = 7, units = "in", res = 300)

ggarrange(shann_spec_plot_2, shann_gen_plot_2, 
          labels = c("A", "B"),
          ncol = 2, nrow = 1)

dev.off()

### Simpson species & genus

tiff(filename = "Simpson_diversity_SG.tiff", width = 8, height = 7, units = "in", res = 300)

ggarrange(simp_spec_plot_2, simp_gen_plot_2, 
          labels = c("A", "B"),
          ncol = 2, nrow = 1)

dev.off()

### Shannon and Simpson species

tiff(filename = "Shannon_and_Simpson_diversity_S.tiff", width = 8, height = 7, units = "in", res = 300)

ggarrange(shann_spec_plot_2, simp_spec_plot_2, 
          labels = c("A", "B"),
          ncol = 2, nrow = 1)

dev.off()

### Shannon and Simpson genus

tiff(filename = "Shannon_and_Simpson_diversity_G.tiff", width = 8, height = 7, units = "in", res = 300)

ggarrange(shann_gen_plot_2, simp_gen_plot_2, 
          labels = c("A", "B"),
          ncol = 2, nrow = 1)

dev.off()

#### Alpha diversity END ####

######################################################################################################

#### Vegan Final Stats START ####

library(vegan)
library(ggrepel)

# Create Bray-Curtis matrices for genus and species

gen <- read.csv("genus_transposed_3g.csv", stringsAsFactors = T, row.names = 1) # This reads a CSV file named "genus_transposed_3g.csv" and assigns its contents to a variable named gen.The argument "genus_transposed_3g.csv" specifies the file name, while row.names = 1 indicates that the first column of the CSV file should be used as row names.
row.names.remove <- c("NC1", "OC", "DES", "PC", "PC1") # This code creates a character vector called "row.names.remove", which contains the row names that should be removed from the "gen" dataframe.
gen <- gen[!row.names(gen) %in% row.names.remove,] # This code removes the rows with the names listed in "row.names.remove" from the "gen" dataframe using subsetting notation.
gen.dist <- vegdist(gen[,-c(93,94,95)], "bray") # This code creates a bray-curtis matrix from the genus data, excluding the descriptive variable columns.
spe <- read.csv("species_transposed_3g.csv", stringsAsFactors = T, row.names = 1) # This reads a CSV file named "species_transposed_3g.csv" and assigns its contents to a variable named spe.The argument "species_transposed_3g.csv" specifies the file name, while row.names = 1 indicates that the first column of the CSV file should be used as row names.
row.names.remove <- c("NC1", "OC", "DES", "PC", "PC1") # This code creates a character vector called "row.names.remove", which contains the row names that should be removed from the "spe" dataframe.
spe <- spe[!row.names(spe) %in% row.names.remove,] # This code removes the rows with the names listed in "row.names.remove" from the "spe" dataframe using subsetting notation.
spe.dist <-vegdist(spe[,-c(233,234,235)], "bray") # This code creates a bray-curtis matrix from the species data, excluding the descriptive variable columns.

# This code computes PERMANOVAs for group, age and base from the bray-curtis matrices using 10,000 permutations

set.seed(2)
stats.gen.group <- adonis(gen.dist ~ Group, data = gen, method = "bray", permutations=10000)
print(stats.gen.group) # SIG
stats.gen.age <- adonis(gen.dist ~ Age, data = gen, method = "bray", permutations=10000)
print(stats.gen.age)
stats.gen.base <- adonis(gen.dist ~ Base, data = gen, method = "bray", permutations=10000)
print(stats.gen.base)

set.seed(2)
stats.spe.group <- adonis(spe.dist ~ Group, data = spe, method = "bray", permutations=10000)
print(stats.spe.group) # SIG
stats.spe.age <- adonis(spe.dist ~ Age, data = spe, method = "bray", permutations=10000)
print(stats.spe.age)
stats.spe.base <- adonis(spe.dist ~ Base, data = spe, method = "bray", permutations=10000)
print(stats.spe.base)

# The following code creates NMDS plots for genus and species using bray-curtis matrices. 

# Genus

NMDS_gen <- gen <- read.csv("genus_transposed_3g.csv", stringsAsFactors = T, row.names = 1) # This reads a CSV file named "genus_transposed_3g.csv" and assigns its contents to a variable named NMDS_gen.The argument "genus_transposed_3g.csv" specifies the file name, while row.names = 1 indicates that the first column of the CSV file should be used as row names.
row.names.remove <- c("NC1", "OC", "DES") # This code creates a character vector called "row.names.remove", which contains the row names that should be removed from the "NMDS_gen" dataframe.
NMDS_gen <- NMDS_gen <- NMDS_gen[!row.names(NMDS_gen) %in% row.names.remove,] # This code removes the rows with the names listed in "row.names.remove" from the "NMDS_gen" dataframe using subsetting notation.

### This code performs an NMDS analysis on the bray-curtis matrices.

set.seed(2)
all.mds <- metaMDS(NMDS_gen[,-c(93,94,95)], trymax = 1000)
data.scores <- as.data.frame(scores(all.mds))
data.scores$group <- NMDS_gen$Group # Assign groups to samples.
data.scores$sample <- row.names(data.scores)

# The following code uses ggplot to create a genus-level NMDS plot

genus_NMDS <- ggplot(data=data.scores,aes(x=NMDS1,y=NMDS2,colour=group, label=sample, family="montserrat"))+
  stat_ellipse(aes(x=NMDS1,y=NMDS2,colour=group),level = 0.5, show.legend = F) +
  geom_point(size=1) +
  scale_color_manual(values=c("Control" = "#57bda3", "Gurkha" = "#6aa2e9","Non-Gurkha" = "#d99a46")) +
  geom_text_repel(family="montserrat", show.legend = F, size = 4)+
  theme_minimal() +
  theme(axis.title = element_text(family="montserrat", size = 16),
        axis.text = element_text(family="montserrat", size = 14),
        legend.text = element_text(family="montserrat", size = 12),
        legend.title = element_blank(),
        legend.box.spacing = unit(0,"pt"),
        legend.margin = margin(0,10,0,0, unit = "pt"))

# Species

NMDS_spe <- NMDS_spe <- read.csv("species_transposed_3g.csv", stringsAsFactors = T, row.names = 1) # This reads a CSV file named "species_transposed_3g.csv" and assigns its contents to a variable named NMDS_spe.The argument "species_transposed_3g.csv" specifies the file name, while row.names = 1 indicates that the first column of the CSV file should be used as row names.
row.names.remove <- c("NC1", "OC", "DES") # This code creates a character vector called "row.names.remove", which contains the row names that should be removed from the "NMDS_spe" dataframe.
NMDS_spe <- NMDS_spe <- NMDS_spe[!row.names(NMDS_spe) %in% row.names.remove,] # This code removes the rows with the names listed in "row.names.remove" from the "NMDS_spe" dataframe using subsetting notation.

### This code performs an NMDS analysis on the bray-curtis matrices.

set.seed(2)
all.mds <- metaMDS(NMDS_spe[,-c(233,234,235)], trymax = 1000)
data.scores <- as.data.frame(scores(all.mds))
data.scores$group <- NMDS_spe$Group
data.scores$sample <- row.names(data.scores)

# The following code uses ggplot to create a species-level NMDS plot

species_NMDS <- ggplot(data=data.scores,aes(x=NMDS1,y=NMDS2,colour=group, label=sample, family="montserrat"))+
  stat_ellipse(aes(x=NMDS1,y=NMDS2,colour=group),level = 0.50, show.legend = F) +
  geom_point(size=1) +
  scale_color_manual(values=c("Control" = "#57bda3", "Gurkha" = "#6aa2e9","Non-Gurkha" = "#d99a46")) +
  geom_text_repel(family="montserrat", show.legend = F, size = 4)+
  theme_minimal() +
  theme(axis.title = element_text(family="montserrat", size = 16),
        axis.text = element_text(family="montserrat", size = 14),
        legend.text = element_text(family="montserrat", size = 12),
        legend.title = element_blank(),
        legend.box.spacing = unit(0,"pt"),
        legend.margin = margin(0,10,0,0, unit = "pt"))

### This code creates a multipanel NMDS figure containing both the species and genus NMDS plots. 

tiff(filename = "NMDS_SG.tiff", width = 8, height = 12, units = "in", res = 300)

ggarrange(species_NMDS, genus_NMDS, 
          labels = c("A", "B"),
          ncol = 1, nrow = 2)

dev.off()

#### Vegan Final Stats END ####

######################################################################################################

### Differential abundance analysis START ###

library(Maaslin2)

# Species

DA_spec <- read.csv("species_transposed_3g.csv", row.names = 1) # This reads a CSV file named "species_transposed_3g.csv" and assigns its contents to a variable named DA_spec.The argument "species_transposed_3g.csv" specifies the file name, while row.names = 1 indicates that the first column of the CSV file should be used as row names.
row.names.remove <- c("PC", "PC1", "NC1", "OC", "DES") # This code creates a character vector called "row.names.remove", which contains the row names that should be removed from the "DA_spec" dataframe.
DA_spec <- DA_spec[!row.names(DA_spec) %in% row.names.remove,] # This code removes the rows with the names listed in "row.names.remove" from the "DA_spec" dataframe using subsetting notation.
DA_spec_meta <- DA_spec[,c(233:235)] # This code creates a dataframe called DA_spec_meta containing the descriptive variables.
DA_spec <- DA_spec[,-c(233:235)] # Remove descriptive variables

#The following code performs an association analysis using maaslin2

fit_data = Maaslin2(input_data     = DA_spec, 
                     input_metadata = DA_spec_meta,
                     output         = "DA_spec", 
                     fixed_effects  = "Group")

# Genus

DA_gen <- read.csv("genus_transposed_3g.csv", row.names = 1) # This reads a CSV file named "genus_transposed_3g.csv" and assigns its contents to a variable named DA_gen.The argument "genus_transposed_3g.csv" specifies the file name, while row.names = 1 indicates that the first column of the CSV file should be used as row names.
row.names.remove <- c("PC", "PC1", "NC1", "OC", "DES") # This code creates a character vector called "row.names.remove", which contains the row names that should be removed from the "DA_gen" dataframe.
DA_gen <- DA_gen[!row.names(DA_gen) %in% row.names.remove,] # This code removes the rows with the names listed in "row.names.remove" from the "DA_gen" dataframe using subsetting notation.
DA_gen_meta <- DA_gen[,c(93,94,95)] # This code creates a dataframe called DA_gen_meta containing the descriptive variables.
DA_gen <- DA_gen[,-c(93,94,95)] # Remove descriptive variables

##The following code performs an association analysis using maaslin2

fit_data2 = Maaslin2(input_data     = DA_gen, 
                     input_metadata = DA_gen_meta,
                     output         = "DA_gen",
                     fixed_effects = "Group")

#### Differential abundance analysis END ####

############################################## END ###################################################
```
