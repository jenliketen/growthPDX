setwd("~/Desktop/growth/data/results")

library(ggplot2)


# Load the PDXE dataset
pdxun <- readRDS("~/Desktop/growth/data/pdxe_untreated.Rda")


# Compute the Pearson correlation between the growth features:
# doubling time (timeToDouble_published), survival (time.last_published), and slope
## Store as variables if you want, but the results are within expected ranges
corGrowthFeature <- function(feature1, feature2, eSet=pData(pdxun)) {
  return(cor(x=eSet[, feature1], y=eSet[, feature2], method="pearson"))
}

corGrowthFeature("timeToDouble_published", "time.last_published")
corGrowthFeature("timeToDouble_published", "slope")
corGrowthFeature("time.last_published", "slope")


# For inspection: how many samples for each tissue type?
number_of_samples <- c()
tissue <- c()
for (i in unique(pData(pdxun)$tumor.type)) {
  a <- nrow(pData(pdxun)[pData(pdxun)$tumor.type==i, ])
  b <- i
  number_of_samples <- c(number_of_samples, a)
  tissue <- c(tissue, b)
}

countSamples <- data.frame(tissue, number_of_samples)
countSamples <- countSamples[order(countSamples$tissue), ]


# Initial diagnostic plots to visualize the mean values of the growth features
growthByTissue <- pData(pdxun)[, c("tumor.type", "timeToDouble_published", "time.last_published", "slope")]
growthByTissue <- growthByTissue[order(growthByTissue$tumor.type), ]
colnames(growthByTissue) <- c("Tissue", "DoublingTime", "Survival", "Slope")

anovaPrelimPlot <- function(tissue, feature) {
  plot <- ggplot(growthByTissue, aes_string(x=tissue, y=feature)) +
    geom_boxplot(fill="grey80", color="blue") +
    scale_x_discrete() + xlab(tissue) +
    ylab(feature) +
    theme_bw() +
    theme(axis.line=element_line(color="black"),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          panel.border=element_blank(),
          panel.background=element_blank())
  return(plot)
}

pdf("anova_prelim_plots.pdf",
    width=12, height=8)
for (i in unique(colnames(growthByTissue))) {
  if (i != "Tissue") {
    print(anovaPrelimPlot("Tissue", i))
  }
}
dev.off()


# One-way ANOVA
## Focus mainly on doubling time
doublingTime.anova <- aov(DoublingTime ~ Tissue, data=growthByTissue)
doublingTime.anova.summary <- summary(doublingTime.anova)
survival.anova <- aov(Survival ~ Tissue, data=growthByTissue)
survival.anova.summary <- summary(survival.anova)
slope.anova <- aov(Slope ~ Tissue, data=growthByTissue)
slope.anova.summary <- summary(slope.anova)

growth.anova <- list("Doubling Time"=doublingTime.anova.summary,
                     "Survival"=survival.anova.summary,
                     "Slope"=slope.anova.summary)

saveRDS(growth.anova, file="growth_anova.Rda")


# Pariwise comparisons
## Again, focus mainly on doubling time
doublingTime.pairwise <- TukeyHSD(doublingTime.anova)
survival.pairwise <- TukeyHSD(survival.anova)
slope.pairwise <- TukeyHSD(slope.anova)

growth.pairwise <- list("Doubling Time"=doublingTime.pairwise,
                        "Survival"=survival.pairwise,
                        "Slope"=slope.pairwise)

saveRDS(growth.pairwise, file="growth_pairwise.Rda")