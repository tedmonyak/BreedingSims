library("ggpubr")
library("grid")
setwd("~/Documents/CSU/R/BreedingSims")

subpop50_qtl20 <- read.csv(file.path(getwd(), "Output/AverageEffectSize/subpop50_qtl20/effSize.csv"))
subpop50_qtl2 <- read.csv(file.path(getwd(), "Output/AverageEffectSize/subpop50_qtl2/effSize.csv"))
subpop500_qtl2 <- read.csv(file.path(getwd(), "Output/AverageEffectSize/subpop500_qtl2/effSize.csv"))
dfs <- list(subpop50_qtl20, subpop50_qtl2, subpop500_qtl2)

output_dir1 <- "~/Documents/CSU/R/BreedingSims/Output/AverageEffectSize/subpop50_qtl20"
output_dir2 <- "~/Documents/CSU/R/BreedingSims/Output/AverageEffectSize/subpop50_qtl2"
output_dir3 <- "~/Documents/CSU/R/BreedingSims/Output/AverageEffectSize/subpop500_qtl2"
output_dirs <- c(output_dir1, output_dir2, output_dir3)

removeOutliers <- "false"
for (i in 1:3) {
  df <- dfs[[i]]
  for (N in c(100)) {
    # ORDER
    order.df <- df %>%
      group_by(orderFixed) %>%
      #filter(!(abs(effectSize - median(effectSize)) > 2*sd(effectSize))) %>%
      summarize(meanEffectSize = mean(effectSize), n=n()) %>%
      filter(n > N)
    
    gO <- ggplot(order.df, aes(x=orderFixed, y=meanEffectSize)) +
      geom_point() +
      theme +
      labs(title="Figure 1", x="Order Fixed", y="Mean Effect Size") +
      geom_smooth(formula = y ~ x, method="loess")
    
    ggplot2::ggsave(filename = paste0("effect_size_order.pdf"),
                    path=output_dirs[i],
                    device = "pdf")
    
    # FITNESS
    fit.df <- df %>%
      mutate(smoothedFit=round(fitness,digits=1)/-4) %>%
      group_by(smoothedFit) %>%
      #filter(!(abs(effectSize - median(effectSize)) > 2*sd(effectSize))) %>%
      summarize(meanEffectSize = mean(effectSize), n=n()) %>%
      filter(n > N)
    
    gF <- ggplot(fit.df, aes(x=smoothedFit, y=meanEffectSize)) +
      geom_point() +
      scale_x_reverse() +
      theme +
      labs(title="Figure 1", x="Normalized Distance from Fitness Optimum", y="Mean Effect Size") +
      geom_smooth(formula = y ~ x, method="loess")
    
    ggplot2::ggsave(filename = paste0("effect_size_fitness.pdf"),
                    path=output_dirs[i],
                    device = "pdf")
    
    # GENERATION
    gen.df <- df %>%
      group_by(gen) %>%
      #filter(!(abs(effectSize - median(effectSize)) > 2*sd(effectSize))) %>%
      summarize(meanEffectSize = mean(effectSize), n=n()) %>%
      filter(n > N)
    
    gG <- ggplot(gen.df, aes(x=gen, y=meanEffectSize)) +
      geom_point() +
      theme +
      labs(title="Figure 1", x="Generation", y="Mean Effect Size") +
      geom_smooth(formula = y ~ x, method="loess")
    
    ggplot2::ggsave(filename = paste0("effect_size_generation.pdf"),
                    path=output_dirs[i],
                    device = "pdf")
    
  }
}



