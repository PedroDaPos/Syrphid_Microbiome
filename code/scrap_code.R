install.packages("aPCoA")
library(aPCoA)
apcoa_otu <- t(as(otu_table(Allograpta_NowolbNS %>% subset_samples(MigratoryStatus == "No")), "matrix"))
apcoa_bray <- vegdist(apcoa_otu, method = "bray")
apcoa_meta <- data.frame(sample_data(Allograpta_NowolbNS %>% subset_samples(MigratoryStatus == "No")))
rownames(apcoa_meta)<-rownames(as.matrix(apcoa_bray))
opar<-par(mfrow=c(1,2),
          mar=c(3.1, 4.1, 4.1, 3.1),
          mgp=c(2, 0.5, 0),
          oma=c(0, 0, 0, 4))
result<-aPCoA(apcoa_bray~Year,apcoa_meta,S_W_N, drawCenter = FALSE, drawEllipse = FALSE)
par(opar)
results <- aPCoA(apcoa_bray ~ Year, data = apcoa_meta, maincov = State, drawCenter = FALSE, drawEllipse = FALSE)


, col = c("red", "green", "blue", "black"))

apcoa_otu <- t(as(otu_table(Eupeodes_Nowolb), "matrix"))
apcoa_bray <- vegdist(apcoa_otu, method = "bray")
apcoa_meta <- data.frame(sample_data(Eupeodes_Nowolb))

results <- aPCoA(apcoa_bray ~ Sex, data = apcoa_meta, maincov = Sex, drawCenter = F, drawEllipse = F)

pseq_beta <- Allograpta_NowolbNS %>% phyloseq::subset_samples(!Season == "Spring")
dist = phyloseq::distance(pseq_beta, method="bray")
ordination = ordinate(pseq_beta, method="PCoA", distance=dist)
plot_ordination(pseq_beta, ordination, color="S_W_N", shape = "MigratoryStatus") + 
  theme_classic() +
  geom_point(size=14, alpha = 0.5) + scale_colour_brewer(type="qual", palette="Set1") +
  theme(strip.background = element_blank(), text = element_text(size = 16)) 
#ggsave("./Figures/PCA_euclidean.pdf", width=15, height=10)


# PCoA with density curves on the x and y-axes
## Allograpta with No Wolbachia, wunifrac
unifrac_dist <- ANoWolb_unifrac_dist_d_1
ps_object <- Allograpta_NowolbNS
pcoa <- cmdscale(as.matrix(unifrac_dist), eig = TRUE, k = 2)
scores <- as.data.frame(pcoa$points)
colnames(scores) <- c("MDS1", "MDS2")

# Add metadata
metadata <- data.frame(sample_data(ps_object))
metadata$MigratoryStatus <- factor(metadata$MigratoryStatus, levels = c("No", "Yes", "Unknown"))
plot_data <- cbind(scores, metadata)

# Create variance variable for the plot axis labels
var_explained <- round(100 * pcoa$eig / sum(pcoa$eig), 2)

# Generate figure (PCoA and boxplots along both axis)
p1 <- ggplot(plot_data, aes(x = MDS1, y = MDS2, color = S_W_N, shape = MigratoryStatus)) +
  geom_point(size = 10, alpha = 0.5) +
  scale_color_brewer(palette = "Set2") +
  labs(
    x = paste0("MDS1 (", var_explained[1], "%)"),
    y = paste0("MDS2 (", var_explained[2], "%)"),
    color = "Population",
    shape = "Migratory?",
    title = "<i>Allograpta</i>, Weighted Unifrac"
  ) +
  theme_minimal() + theme(legend.text = element_text(size = 14), legend.title = element_text(size = 14), plot.title = element_markdown()) +
  scale_shape_manual(values= 15:17) +
  ggside::geom_xsideboxplot(aes(fill = S_W_N, y = S_W_N, group = S_W_N), orientation = "y", show.legend = FALSE) +
  ggside::geom_ysideboxplot(aes(fill = S_W_N, x = S_W_N, group = S_W_N), orientation = "x", show.legend = FALSE) +
  ggside::scale_xsidey_discrete(labels = NULL) +
  ggside::scale_ysidex_discrete(labels = NULL) +
  ggside::theme_ggside_void()+
  scale_fill_brewer(palette = "Set2")

#Allograpta with no Wolbachia, unweighted unifrac

unifrac_dist <- ANoWolb_unifrac_distUW
ps_object <- Allograpta_NowolbNS
pcoa <- cmdscale(as.matrix(unifrac_dist), eig = TRUE, k = 2)
scores <- as.data.frame(pcoa$points)
colnames(scores) <- c("MDS1", "MDS2")

# Add metadata
metadata <- data.frame(sample_data(ps_object))
metadata$MigratoryStatus <- factor(metadata$MigratoryStatus, levels = c("No", "Yes", "Unknown"))
plot_data <- cbind(scores, metadata)

# Create variance variable for the plot axis labels
var_explained <- round(100 * pcoa$eig / sum(pcoa$eig), 2)

# Generate figure (PCoA and boxplots along both axis)
p2 <- ggplot(plot_data, aes(x = MDS1, y = MDS2, color = S_W_N, shape = MigratoryStatus)) +
  geom_point(size = 10, alpha = 0.5) +
  scale_color_brewer(palette = "Set2") +
  labs(
    x = paste0("MDS1 (", var_explained[1], "%)"),
    y = paste0("MDS2 (", var_explained[2], "%)"),
    color = "Population",
    shape = "Migratory?",
    title = "<i>Allograpta</i>, Unweighted Unifrac"
  ) +
  theme_minimal() + theme(legend.text = element_text(size = 14), legend.title = element_text(size = 14), plot.title = element_markdown()) +
  scale_shape_manual(values= 15:17) +
  ggside::geom_xsideboxplot(aes(fill = S_W_N, y = S_W_N, group = S_W_N), orientation = "y", show.legend = FALSE) +
  ggside::geom_ysideboxplot(aes(fill = S_W_N, x = S_W_N, group = S_W_N), orientation = "x", show.legend = FALSE) +
  ggside::scale_xsidey_discrete(labels = NULL) +
  ggside::scale_ysidex_discrete(labels = NULL) +
  ggside::theme_ggside_void()+
  scale_fill_brewer(palette = "Set2")

## Eupeodes with No Wolbachia, wunifrac
unifrac_dist <- EaNoWolb_unifrac_dist_d_1
ps_object <- Eupeodes_Nowolb
pcoa <- cmdscale(as.matrix(unifrac_dist), eig = TRUE, k = 2)
scores <- as.data.frame(pcoa$points)
colnames(scores) <- c("MDS1", "MDS2")

# Add metadata
metadata <- data.frame(sample_data(ps_object))
metadata$MigratoryStatus <- factor(metadata$MigratoryStatus, levels = c("No", "Yes", "Unknown"))
plot_data <- cbind(scores, metadata)

# Create variance variable for the plot axis labels
var_explained <- round(100 * pcoa$eig / sum(pcoa$eig), 2)

# Generate figure (PCoA and boxplots along both axis)
p3 <- ggplot(plot_data, aes(x = MDS1, y = MDS2, color = S_W_N, shape = MigratoryStatus)) +
  geom_point(size = 10, alpha = 0.5) +
  scale_color_brewer(palette = "Set2") +
  labs(
    x = paste0("MDS1 (", var_explained[1], "%)"),
    y = paste0("MDS2 (", var_explained[2], "%)"),
    color = "Population",
    shape = "Migratory?",
    title = "<i>Eupeodes</i>, Weighted Unifrac"
  ) +
  theme_minimal() + theme(legend.text = element_text(size = 14), legend.title = element_text(size = 14), plot.title = element_markdown()) +
  scale_shape_manual(values= 15:17) +
  ggside::geom_xsideboxplot(aes(fill = S_W_N, y = S_W_N, group = S_W_N), orientation = "y", show.legend = FALSE) +
  ggside::geom_ysideboxplot(aes(fill = S_W_N, x = S_W_N, group = S_W_N), orientation = "x", show.legend = FALSE) +
  ggside::scale_xsidey_discrete(labels = NULL) +
  ggside::scale_ysidex_discrete(labels = NULL) +
  ggside::theme_ggside_void()+
  scale_fill_brewer(palette = "Set2")

#Euepeodes with no Wolbachia, unweighted unifrac

unifrac_dist <- EaNoWolb_unifrac_distUW
ps_object <- Eupeodes_Nowolb
pcoa <- cmdscale(as.matrix(unifrac_dist), eig = TRUE, k = 2)
scores <- as.data.frame(pcoa$points)
colnames(scores) <- c("MDS1", "MDS2")

# Add metadata
metadata <- data.frame(sample_data(ps_object))
metadata$MigratoryStatus <- factor(metadata$MigratoryStatus, levels = c("No", "Yes", "Unknown"))
plot_data <- cbind(scores, metadata)

# Create variance variable for the plot axis labels
var_explained <- round(100 * pcoa$eig / sum(pcoa$eig), 2)

# Generate figure (PCoA and boxplots along both axis)
p4 <- ggplot(plot_data, aes(x = MDS1, y = MDS2, color = S_W_N, shape = MigratoryStatus)) +
  geom_point(size = 10, alpha = 0.5) +
  scale_color_brewer(palette = "Set2") +
  labs(
    x = paste0("MDS1 (", var_explained[1], "%)"),
    y = paste0("MDS2 (", var_explained[2], "%)"),
    color = "Population",
    shape = "Migratory?",
    title = "<i>Eupeodes</i>, Unweighted Unifrac"
  ) +
  theme_minimal() + theme(legend.text = element_text(size = 14), legend.title = element_text(size = 14), plot.title = element_markdown()) +
  scale_shape_manual(values= 15:17) +
  ggside::geom_xsideboxplot(aes(fill = S_W_N, y = S_W_N, group = S_W_N), orientation = "y", show.legend = FALSE) +
  ggside::geom_ysideboxplot(aes(fill = S_W_N, x = S_W_N, group = S_W_N), orientation = "x", show.legend = FALSE) +
  ggside::scale_xsidey_discrete(labels = NULL) +
  ggside::scale_ysidex_discrete(labels = NULL) +
  ggside::theme_ggside_void()+
  scale_fill_brewer(palette = "Set2")

p5 <- plot_grid(p1 + theme(legend.position="none"),
          p3 + theme(legend.position="none"),
          p2 + theme(legend.position="none"), 
          p4 + theme(legend.position="none"), labels = c("A", "C", "B", "D"), ncol = 2, nrows = 2)

# extract the legend from one of the plots
legend <- get_legend(
  # create some space to the left of the legend
  p1 + theme(legend.box.margin = margin(0, 1, 1, 12))
)

plot_grid(p5, legend, rel_widths = c(6, 1))

library(microbiome)
library(DirichletMultinomial)
library(parallel)
numCores <- detectCores()
data_focus <- Allograpta_NowolbNS %>% subset_samples(MigratoryStatus == "No")
data_genus <- tax_glom(data_focus, "Genus")
data_genuscomp <- microbiome::transform(data_genus, "compositional")
taxa <- core_members(data_genuscomp, detection = 0.1/100, prevalence = 20/100)
data_genuscomp <- prune_taxa(taxa, data_genuscomp)
trim_data <- transform_sample_counts(Wolb_rr_up_filt_genus, function(x) x / sum(x))
trim_data_filt <- filter_taxa(trim_data, function(x) mean(x) < 0.0001, TRUE)
trim_data_filt_lean <- prune_samples(sample_sums(trim_data_filt) > 0, trim_data_filt)

dat <- abundances(data_focus)
count <- as.matrix(t(dat))
#fit <- lapply(1:10, dmn, count = count, verbose=TRUE)
fit <- mclapply(1:5, dmn, count = count, verbose=TRUE, mc.cores = numCores)
laplace <- sapply(fit, laplace)
plot(laplace, type="b", xlab = "Number of Dirichlet Components", ylab="Laplace",lwd = 2)


dbRDA <- capscale(apcoa_otu ~ S_W_N + Sex + Year + Season, data=apcoa_meta, distance = "bray")
plot(dbRDA)
plot(dbRDA, type = "none")
points(dbRDA, display = "sites", 
       col = apcoa_meta$S_W_N, cex = 1.4)
with(dbRDA, legend("topleft", legend=levels(apcoa_meta$S_W_N), cex =1.4, col =c("red", "green"[as.numeric(apcoa_meta$S_W_N)])))



# MDS visualization
nmds <- metaMDS(ANoWolb_unifrac_dist_d_05, k =2, trymax = 100)
nmds_df <- data.frame(NMDS1 = nmds$points[,1], NMDS2 = nmds$points[,2], Population = metadata_ANoWolb$S_W_N, MigratoryStatus = )
ggplot(nmds_df, aes(x=NMDS1, y=NMDS2, color = Population)) +
  geom_point(size = 4) +
  theme_minimal()



# boxplot of data dispersion
# Allograpta No Wolbachia; are non-local individuals pushing data dispersion for this mostly non-migratory species?
disper_df <- data.frame(Population = metadata_ANoWolb$S_W_N, Distance = disper_NS_ANoWolb$distances, MigratoryStatus = metadata_ANoWolb$MigratoryStatus)
ggplot(disper_df, aes(x= Population, y = Distance, fill = Population)) + geom_boxplot() + 
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), strip.text = element_text(face = "italic", size = 14),
        legend.text = element_markdown(size = 14), legend.title = element_markdown(size = 14)) +
        geom_jitter(aes(color = MigratoryStatus), width = 0.2) +
        scale_fill_manual(values = alpha(c("North" = "orange", "South" = "purple"), .3), 
                           labels = c("North", "South")) +
  scale_color_manual(values = c("No" = "brown", "Unknown" = "black", "Yes" = "lightgray"),
                     labels = c("Non-migratory", "Unknown", "Migratory"))
ggsave("./Documents/Syrphids_microbiome/Syrphid_Microbiome/results/Figures/Boxplot_disper.pdf")

# excluding unknown and non-local samples
ggplot(disper_df %>% filter(MigratoryStatus == "No"), aes(x= Population, y = Distance, fill = Population)) + geom_boxplot(outlier.shape = NA) + 
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), strip.text = element_text(face = "italic", size = 14),
        legend.text = element_markdown(size = 14), legend.title = element_markdown(size = 14)) +
  geom_jitter(aes(color = MigratoryStatus), width = 0.2) +
  scale_fill_manual(values = alpha(c("North" = "orange", "South" = "purple"), .3), 
                    labels = c("North", "South")) +
  scale_color_manual(values = c("No" = "brown"), 
                     labels = "Non-migratory")
ggsave("./Documents/Syrphids_microbiome/Syrphid_Microbiome/results/Figures/Boxplot_disper_NonMigratory.pdf")


# boxplot of data dispersion
# Eupeodes No Wolbachia;
disper_df <- data.frame(Population = metadata_EaNoWolb$S_W_N, Distance = disper_NS_EaNoWolb$distances, MigratoryStatus = metadata_EaNoWolb$MigratoryStatus)
ggplot(disper_df, aes(x= Population, y = Distance, fill = Population)) + geom_boxplot() + 
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), strip.text = element_text(face = "italic", size = 14),
        legend.text = element_markdown(size = 14), legend.title = element_markdown(size = 14)) +
  geom_jitter(aes(color = MigratoryStatus), width = 0.2) +
  scale_fill_manual(values = alpha(c("North" = "orange", "South" = "purple"), .3), 
                    labels = c("North", "South")) +
  scale_color_manual(values = c("No" = "brown", "Unknown" = "black", "Yes" = "darkgray"),
                     labels = c("Non-migratory", "Unknown", "Migratory"))

ggsave("./Documents/Syrphids_microbiome/Syrphid_Microbiome/results/Figures/Eupeodes_Boxplot_disper.pdf")

# focusing on only migrating individuals
ggplot(disper_df %>% filter(MigratoryStatus == "Yes"), aes(x= Population, y = Distance, fill = Population)) + geom_boxplot(outlier.shape = NA) + 
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), strip.text = element_text(face = "italic", size = 14),
        legend.text = element_markdown(size = 14), legend.title = element_markdown(size = 14)) +
  geom_jitter(aes(color = MigratoryStatus), width = 0.2) +
  scale_fill_manual(values = alpha(c("North" = "orange", "South" = "purple"), .3), 
                    labels = c("North", "South")) +
  scale_color_manual(values = c("Yes" = "darkgray"), 
                     labels = "Migratory")
ggsave("./Documents/Syrphids_microbiome/Syrphid_Microbiome/results/Figures/Eupeodes_Boxplot_disper_only_Migratory.pdf")
