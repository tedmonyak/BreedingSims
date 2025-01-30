library("ggpubr")
library("grid")
setwd("~/Documents/CSU/R/BreedingSims")

output_dir <- file.path(getwd(), "AggregateQtl")
if (!dir.exists(output_dir)) dir.create(output_dir)

ne50qtl2h20 <- read.csv("Output/QtlMonteCarlo/Ne_50_qtl_2_selProp_0.1_h2_0.2_gens_100_2025-01-26_14_54/result_dataframe.csv") %>%
  mutate(h2=0.2,qtl=2,ne=50,sel=0.1)
ne50qtl2h50 <- read.csv("Output/QtlMonteCarlo/Ne_50_qtl_2_selProp_0.25_h2_0.5_gens_100_2025-01-24_12_35/result_dataframe.csv") %>%
  mutate(h2=0.5,qtl=2,ne=50,sel=0.25)
ne50qtl20h20 <- read.csv("Output/QtlMonteCarlo/Ne_50_qtl_20_selProp_0.1_h2_0.2_gens_100_2025-01-27_11_27/result_dataframe.csv") %>%
  mutate(h2=0.2,qtl=20,ne=50,sel=0.1)
ne50qtl20h50 <- read.csv("Output/QtlMonteCarlo/Ne_50_qtl_20_selProp_0.25_h2_0.5_gens_100_2025-01-24_15_00/result_dataframe.csv") %>%
  mutate(h2=0.5,qtl=20,ne=50,sel=0.25)
ne500qtl2h20 <- read.csv("Output/QtlMonteCarlo/Ne_500_qtl_2_selProp_0.1_h2_0.2_gens_100_2025-01-26_16_46/result_dataframe.csv") %>%
  mutate(h2=0.2,qtl=2,ne=500,sel=0.1)
ne500qtl2h50 <- read.csv("Output/QtlMonteCarlo/Ne_500_qtl_2_selProp_0.25_h2_0.5_gens_100_2025-01-24_17_42/result_dataframe.csv") %>%
  mutate(h2=0.5,qtl=2,ne=500,sel=0.25)
ne500qtl20h20 <- read.csv("Output/QtlMonteCarlo/Ne_500_qtl_20_selProp_0.1_h2_0.2_gens_100_2025-01-27_13_34/result_dataframe.csv") %>%
  mutate(h2=0.2,qtl=20,ne=500,sel=0.1)
ne500qtl20h50 <- read.csv("Output/QtlMonteCarlo/Ne_500_qtl_20_selProp_0.25_h2_0.5_gens_100_2025-01-25_06_06/result_dataframe.csv") %>%
  mutate(h2=0.5,qtl=20,ne=500,sel=0.25)

olig_small.df <- rbind(ne50qtl2h20,ne50qtl2h50)
olig_large.df <- rbind(ne500qtl2h20,ne500qtl2h50)
poly_small.df <- rbind(ne50qtl20h20,ne50qtl20h50)
poly_large.df <- rbind(ne500qtl20h20,ne500qtl20h50)

ymax = max(c(olig_small.df$nSigQtl,
           olig_large.df$nSigQtl,
           poly_small.df$nSigQtl,
           poly_large.df$n.SigQtl))

theme <- theme(
  axis.title.x = element_text(family="Helvetica", size=10),
  axis.title.y = element_text(family="Helvetica", size=10),
  plot.title = element_text(family="Helvetica", size=12, hjust = 0.5),
  legend.text = element_text(family="Helvetica", size=12),
  legend.title = element_text(family="Helvetica", size=12),
  legend.key = element_rect(linewidth=0.05),
  plot.caption = element_text(family="Helvetica", size=10, hjust = 0.5),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_rect(fill = "white", color = "black"),
  aspect.ratio = 1)

olig_small <- ggplot(data=olig_small.df, aes(x=as.factor(h2), y=nSigQtl, fill=type)) +
  geom_boxplot(outlier.shape = NA) +
  ylim(0,ymax) +
  scale_fill_manual(name = "Parental Cross Type",
                    values = c("Inter" = "orange",
                               "Intra" = "dodgerblue2"),
                    labels = c("Inter-population", "Intra-population")) +
  theme +
  labs(x="h2",
       y="Significant QTL")

olig_large <- ggplot(data=olig_large.df, aes(x=as.factor(h2), y=nSigQtl, fill=type)) +
  geom_boxplot(outlier.shape = NA) +
  ylim(0,ymax) +
  scale_fill_manual(name = "Parental Cross Type",
                    values = c("Inter" = "orange",
                               "Intra" = "dodgerblue2"),
                    labels = c("Inter-population", "Intra-population")) +
  theme +
  labs(x="h2",
       y="Significant QTL")

poly_small <- ggplot(data=poly_small.df, aes(x=as.factor(h2), y=nSigQtl, fill=type)) +
  geom_boxplot(outlier.shape = NA) +
  ylim(0,ymax) +
  scale_fill_manual(name = "Parental Cross Type",
                    values = c("Inter" = "orange",
                               "Intra" = "dodgerblue2"),
                    labels = c("Inter-population", "Intra-population")) +
  theme +
  labs(x="h2",
       y="Significant QTL")

poly_large <- ggplot(data=poly_large.df, aes(x=as.factor(h2), y=nSigQtl, fill=type)) +
  geom_boxplot(outlier.shape = NA) +
  ylim(0,ymax) +
  scale_fill_manual(name = "Parental Cross Type",
                    values = c("Inter" = "orange",
                               "Intra" = "dodgerblue2"),
                    labels = c("Inter-population", "Intra-population")) +
  theme +
  labs(x="h2",
       y="Significant QTL")
label_size = 8
olig_label <- wrap_elements(panel = textGrob('Oligogenic\nTrait Architecture\n(QTL per Trait = 20',
                                             rot=90,
                                             gp=gpar(fontsize=10)))
poly_label <- wrap_elements(panel = textGrob('Polygenic\nTrait Architecture\n(QTL per Trait = 200)',
                                             rot=90,
                                             gp=gpar(fontsize=10)))
small_label <- wrap_elements(panel = textGrob('Small Subopulation\n(Ne = 50)',
                                              gp=gpar(fontsize=10)))
large_label <- wrap_elements(panel = textGrob('Large Subpopulation\n(Ne = 500)',
                                              gp=gpar(fontsize=10)))

p <- plot_spacer()+small_label+large_label+olig_label+olig_small+olig_large+poly_label+poly_small+poly_large + plot_layout(guides='collect', axes='collect', widths=c(1,2,2), heights=c(1,2,2))& theme(legend.position="right")

fname <- file.path(output_dir, "aggregate_qtl.pdf")
ggplot2::ggsave(filename = fname,
                device = "pdf",
                height=7,
                width=8)
