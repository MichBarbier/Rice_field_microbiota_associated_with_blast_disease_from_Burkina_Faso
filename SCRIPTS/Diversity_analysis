##########################################################################################
##    This script allows to asses the alpha diversity and to plot the beta diversity    ##
##########################################################################################


setwd("C:/[[YOUR PATH]]/Rice_field_microbiota_associated_with_blast_disease_from_Burkina_Faso")


#### PACKAGES ####

library("data.table")
library("ggnewscale")
library("ggpubr")
library("ggtext")
library("patchwork")
library("tidyverse")

library("microeco") # Allows to asses the alpha diversity
library("vegan") # Allows to performe the PERMANOVA and the NMDS



#### ALPHA DIVERSITY ####

Data <- read.csv("DATA/ASV_TABLE.csv", sep = ";", dec = ".", header = TRUE)

Metadata <- read.csv("DATA/METADATA.csv", sep = ";", dec = ".", header = TRUE)

Taxa <- Data[,1:7]
ASVs <- Data[,9:95]

# Calculate different alpha diversity indices 
Alpha <- microtable$new(otu = ASVs, tax = Taxa)
Alpha <- Alpha$cal_alphadiv()

Alpha <- as.data.frame(Alpha$alpha_diversity)
Alpha[,"Disease"] <- Metadata[,"Disease"]

Colors <- c("orangered", "steelblue")

# Specific richness 
p1 <- ggplot(Alpha)+
  geom_boxplot(aes(x = Disease, y = Observed, fill = Disease, color = Disease), alpha = 0.5)+
  scale_fill_manual("Disease", values = Colors, name = "", guide = "none")+
  scale_color_manual("Disease", values = Colors, name = "", guide = "none")+
  geom_point(aes(x = Disease, y = Observed, color = Disease))+
  scale_color_manual("Disease", values = Colors, name = "", guide = "none")+
  geom_signif(aes(x = Disease, y = Observed), test = wilcox.test, comparisons = list(c("Healthy", "Diseased")), map_signif_level = TRUE, y_position = 1100, textsize = 8)+
  stat_summary(aes(x = Disease, y = Observed), fun = mean, colour = "black", geom = "point", shape = "+", size = 6.5)+
  ggtitle("Observed richness")+
  ylab("Number of ASVs")+
  scale_y_continuous(limits = c(0, 1250), breaks = seq(0, 1250, 250))+
  theme_bw()+
  theme(plot.title = element_text(face = "bold", size = 13, hjust = 0.5), axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 12), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.text = element_text(size = 12))
p1

# Shannon index
p2 <- ggplot(Alpha)+
  geom_boxplot(aes(x = Disease, y = Shannon, fill = Disease, color = Disease), alpha = 0.5)+
  scale_fill_manual("Disease", values = Colors, name = "", guide = "none")+
  scale_color_manual("Disease", values = Colors, name = "", guide = "none")+
  geom_point(aes(x = Disease, y = Shannon, color = Disease))+
  scale_color_manual("Disease", values = Colors, name = "", guide = "none")+
  geom_signif(aes(x = Disease, y = Shannon), test = wilcox.test, comparisons = list(c("Healthy", "Diseased")), map_signif_level = TRUE, y_position = 6.56, textsize = 8)+
  stat_summary(aes(x = Disease, y = Shannon), fun = mean, colour = "black", geom = "point", shape = "+", size = 6.5)+
  ggtitle("Shannon index")+
  ylab("Index")+
  scale_y_continuous(limits = c(1.5, 7.5), breaks = seq(1.5, 7.5, 1))+
  theme_bw()+
  theme(plot.title = element_text(face = "bold", size = 13, hjust = 0.5), axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 12), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.text = element_text(size = 12))
p2



#### BETA DIVERSITY ####

Data <- read.csv("DATA/ASV_TABLE.csv", sep = ";", dec = ".", header = TRUE)

Metadata <- read.csv("DATA/METADATA.csv", sep = ";", dec = ".", header = TRUE)

rownames(Data) <- Data[,"ASV"]
Data <- Data[,-1:-8]
Data <- as.data.frame(t(Data))

Data[,"Field"] <- Metadata[,"Field"]
Data[,"Site"] <- Metadata[,"Site"]
Data[,"Disease"] <- Metadata[,"Disease"]

# Performe a NMDS based on the Bray-Curtis distances with 1000 permutations
NMDS <- metaMDS(Data[,-14699:-14701], distance = "bray", try = 1000, trymax = 1000)
NMDS_scores <- scores(NMDS)
NMDS_scores_sites <- as.data.frame(NMDS_scores[["sites"]])

NMDS$stress
stressplot(NMDS)
Stress_value <- paste0("Stress: ", round(NMDS$stress, 4))

Distances <- vegdist(Data[,-14699:-14701], method = "bray")

# Performe a PERMANOVA to asses the impact of the sampling site
PERMANOVA <- adonis2(Distances~Data[,"Site"], data = Data[,-14699:-14701], permutations = 1000, method = "bray")
PERMANOVA_SITE <- paste0("PERMANOVA (Site): p = ", round(PERMANOVA$`Pr(>F)`[1], 4), "; F-ratio = ", round(PERMANOVA$`F`[1], 4), "; R2 = ", round(PERMANOVA$`R2`[1], 4))

# Performe a PERMANOVA to asses the impact of the field health 
PERMANOVA <- adonis2(Distances~Data[,"Disease"], data = Data[,-14699:-14701], permutations = 1000, method = "bray")
PERMANOVA_DISEASE <- paste0("PERMANOVA (Disease): p = ", round(PERMANOVA$`Pr(>F)`[1], 4), "; F-ratio = ", round(PERMANOVA$`F`[1], 4), "; R2 = ", round(PERMANOVA$`R2`[1], 4))

# Performe a PERMANOVA to asses the synergic impact of the sampling site and the field health 
PERMANOVA <- adonis2(Distances~Data[,"Site"]*Data[,"Disease"], data = Data[,-14699:-14701], permutations = 1000, method = "bray")
PERMANOVA_BOTH <- paste0("PERMANOVA (Site x Disease): p = ", round(PERMANOVA$`Pr(>F)`[1], 4), "; F-ratio = ", round(PERMANOVA$`F`[1], 4), "; R2 = ", round(PERMANOVA$`R2`[1], 4))

Colors <- c("orangered", "steelblue")
Shapes <- c(21,22,24,25)

p3 <- ggplot(NMDS_scores_sites, aes(NMDS1,NMDS2))+
  stat_ellipse(data = NMDS_scores_sites, aes(group = Data$Disease, color = Data$Disease), type = "norm", show.legend = FALSE)+
  geom_point(aes(shape = Data$Site, fill = Data$Disease, color = Data$Disease), alpha = 0.75, size = 3)+
  scale_color_manual(Data$Disease, values = Colors, name = "")+
  scale_fill_manual(Data$Disease, values = Colors, name = "")+
  scale_shape_manual(Data$Disease, values = Shapes, name = "Sampling sites")+
  scale_x_continuous(limits = c(-4, 3), breaks = seq(-4, 3, 1))+
  scale_y_continuous(limits = c(-2.5, 3.25), breaks = seq(-2, 3, 1))+
  theme_bw()+
  annotate("text", x = -4, y = 3, label = Stress_value, size = 3.5, hjust= 0)+
  annotate("text", x = -4, y = 2.75, label = PERMANOVA_SITE, size = 3.5, hjust= 0)+
  annotate("text", x = -4, y = 2.5, label = PERMANOVA_DISEASE, size = 3.5, hjust= 0)+
  annotate("text", x = -4, y = 2.25, label = PERMANOVA_BOTH, size = 3.5, hjust= 0)+
  theme(axis.title.y = element_text(size = 12), axis.text.y = element_text(size = 12), axis.title.x = element_text(size = 12), axis.text.x = element_text(size = 12), legend.text = element_text(size = 12), legend.title = element_text(face = "bold", size = 12))
p3



#### SUPPLEMENTARY FIGURE 4 ####

Model <- "ABCCC"

p10 <- p1 + p2 + p3 + plot_layout(design = Model, guides = "collect") 
p10

ggsave(plot = p10, dpi = 1000, device = "jpeg", width = 12, height = 8, filename = "FIGURES/SUPPLEMENTARY_FIGURE_4.jpeg")
