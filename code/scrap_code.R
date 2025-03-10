install.packages("aPCoA")
library(aPCoA)
apcoa_otu <- t(as(otu_table(Allograpta_NowolbNS %>% subset_samples(MigratoryStatus == "No")), "matrix"))
apcoa_bray <- vegdist(apcoa_otu, method = "bray")
apcoa_meta <- data.frame(sample_data(Allograpta_NowolbNS %>% subset_samples(MigratoryStatus == "No")))
rownames(apcoa_meta)<-rownames(as.matrix(apcoa_bray))
opar<-par(mfrow=c(1,2),
          mar=c(3.1, 3.1, 3.1, 5.1),
          mgp=c(2, 0.5, 0),
          oma=c(0, 0, 0, 4))
result<-aPCoA(apcoa_bray~Year,apcoa_meta,S_W_N, drawCenter = TRUE, drawEllipse = FALSE)
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
plot_ordination(pseq_beta, ordination, color="S_W_N") + 
  theme_classic() +
  geom_point(size=8, alpha = 0.5) + scale_colour_brewer(type="qual", palette="Set1") +
  theme(strip.background = element_blank(), text = element_text(size = 16)) 
#ggsave("./Figures/PCA_euclidean.pdf", width=15, height=10)



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




table(meta(Allograpta_NowolbNS)$S_W_N[1:10], meta(Allograpta_NowolbNS)$ID[1:10])


nmds <- metaMDS()

# boxplot of data dispersion
disper_df <- data.frame(Population = metadata_ANoWolb$S_W_N, Distance = disper_NS_ANoWolb$distances, MigratoryStatus = metadata_ANoWolb$MigratoryStatus)
ggplot(disper_df, aes(x= Population, y = Distance, fill = Population)) + geom_boxplot() + 
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), strip.text = element_text(face = "italic", size = 14),
        legend.text = element_markdown(size = 14), legend.title = element_markdown(size = 14)) +
        geom_jitter(aes(color = MigratoryStatus), width = 0.2) +
        scale_fill_manual(values = alpha(c("North" = "orange", "South" = "purple"), .3), 
                           labels = c("North", "South")) +
  scale_color_manual(values = c("No" = "brown", "Unknown" = "black", "Yes" = "lightgray"))
ggsave("./Documents/Syrphids_microbiome/Syrphid_Microbiome/results/Figures/Boxplot_disper.pdf")



ggplot(disper_df %>% filter(MigratoryStatus == "No"), aes(x= Population, y = Distance, fill = Population, shape = MigratoryStatus)) + geom_boxplot() + theme_minimal() + geom_jitter()

