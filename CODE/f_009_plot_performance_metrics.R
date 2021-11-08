################################################################################
################################################################################
# ADNI NIGHTINGALE (NG) LONGITUDINAL ## V001 ## MAFF ###########################
# F009 DISTRIBUTION PERFORMANCE METRICS:
#    + BY PHENOTYPE
# ---------------------------------------------------------------------------- #

library(ggplot2)
library(data.table)
library(tidyverse)


temp = readRDS("./RESULTS/BPDIA_cv10_lmx.rds")[, c(301:307)]
temp %>%
  group_by(metabolites) %>%
  dplyr::summarize(Mean = mean(RMSE, na.rm=TRUE)) %>%
  dplyr::arrange(Mean)


L_HDL_C <- subset(temp, metabolites == "L_HDL_C")

pdf(file = "./RESULTS/plot_performance_metrics_BPDIA_L_HDL_C.pdf")
L_HDL_C %>%
  gather(.,MSE:BIC,key ="Metric",value = "Value")%>%
  ggplot(aes(x=Metric,y=Value,fill=Metric)) +
  geom_boxplot()+coord_flip() +
  facet_wrap(~Metric,ncol=1,scales="free") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(color = "blue")
    ) +
  labs(title = "Pheno::BPDIA",
       subtitle = "L_HDL_C",
       caption = "Data source: Nightingale longitudinal")
dev.off()