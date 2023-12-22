  ### Figures from experiments with the simulated data
library(ggplot2)
library(ggpubr)
library(ggtext)
library(ggh4x)
library(ggsci)
library(gridExtra)
library(patchwork)
library(ggthemes)
source("utils.R")

#######################################################################
### Main Figure 2
load(file = "data/sim/fig_data/fig2.RData")
fig2_df$Hypothesis <- paste0("cophescan ", fig2_df$Hypothesis)
fig2_df$Hypothesis <- factor(fig2_df$Hypothesis,
levels = c("cophescan Hn", "cophescan Ha", "cophescan Hc"))
g2a <- ggplot(fig2_df, aes(y = type, x = value, fill = Hypothesis)) +
  geom_col(width = 0.8) +
  facet_nested(Priors + BF + rg ~  truth, scales = "free",
  space = "free", remove_labels = "y") +
  theme_minimal()  %+replace%
  theme(strip.placement = "outside", text = element_text(size = 15,
    colour = "grey15"), legend.position = "bottom",
    axis.text.y =  element_blank(),
    plot.title = element_text(hjust = 0.5),
    strip.text.x = element_textbox(size = 15, color = "black",
    fill = "#e6e6e6ff", box.color = "#e6e6e6ff",
    halign =  0.5, linetype = 1,
    r = unit(5, "pt"), width = unit(1, "npc"),
    padding = margin(2, 0, 1, 0), margin = margin(3, 3, 3, 3),
      face = "bold"),
    strip.text.y = element_textbox(orientation = "right-rotated", size = 12,
    color = "black", fill = "aliceblue",
    box.color = "cornflowerblue", halign  =  0.5,
    linetype = 1, r = unit(5, "pt"), width = unit(1, "npc"),
    padding = margin(2, 0, 1, 0),
    margin = margin(3, 3, 3, 3), face = "bold"),
    axis.title.x = element_blank(), axis.text.x = element_blank(),
    axis.ticks.x = element_blank()) +
    labs(x = "", y = "") +
  scale_fill_viridis_d(name = "Predicted finding",
    alpha = 0.7, option =  "D") +
  geom_text(aes(x = value, label = value_label, color = Hypothesis),
  hjust = 1.1, fontface = "bold", size = 5) +
  scale_color_manual(values = c("grey90", "grey30", "grey22"), guide = "none")

load(file = "data/sim/fig_data/fig2b.RData")
fig2b_melt$variable <- factor(fig2b_melt$variable,
  levels = c("notAssoc", "Assoc"), labels =
    c("FDR not significant", "FDR significant"))
g2b <- ggplot(fig2b_melt, aes(y = type, x = value, fill = variable)) +
  geom_col(width = 0.8) +
  facet_nested(type ~ truth, scales = "free", space = "free",
  remove_labels = "y") +
  theme_minimal() %+replace%
  theme(strip.placement = "outside", text = element_text(size = 15,
    colour = "grey15"), legend.position = "bottom",
    axis.text.y = element_blank(), plot.title = element_text(hjust = 0.5),
    strip.text.x = element_blank(),
    strip.text.y = element_textbox(orientation = "right-rotated",
    size = 12, color = "black", fill = "aliceblue",
    box.color = "cornflowerblue", halign = 0.5, linetype = 1,
    r = unit(5, "pt"), width = unit(1, "npc"),
    padding = margin(2, 0, 1, 0), margin = margin(3, 3, 3, 3),
    face = "bold")) +
  labs(x = "Percentage of traits analysed", y = "") +
  scale_fill_viridis_d(name = "", alpha = 0.7, option = "A") +
  geom_text(aes(x = value, label = value_label, color = variable),
  hjust = 1.1, fontface = "bold", size = 5) +
  scale_color_manual(values = c("grey90", "grey30", "grey22"), guide = "none")

g2 <- g2a / g2b +
  plot_layout(heights = unit(c(12, 3.5), c("cm", "cm"))) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom", plot.margin = unit(c(1, 1, 0, 1), "cm"))
# g2
ggsave("figures/final_figures/fig2ab.png",
  width = 13, height = 10,
  bg = "white", units = "in", dpi = 300)
ggsave("figures/final_figures/fig2ab.pdf",
  width = 13, height = 10, bg = "white", units = "in")

#######################################################################
### Main Figure 3
#### Figure 3a
load(file = "data/sim/fig_data/fig3a.RData")
fig3a_df$hc <- fig3a_df$x < 0.6 & fig3a_df$y >= 0.6
fig3a_df$nohc <- fig3a_df$x < 0.6
####   (59+21)/(2203+664)*100
####  [1] 2.790373
####   (1/4637+4627) * 100
####  [1] 0.01079447
####   (77/87715) * 100
####  [1] 0.0877843

fig3a_df$truth_grped <- gsub("True |2", "", fig3a_df$truth)
g3a <- ggplot(aes(x = x, y = y), data = fig3a_df) +
  facet_wrap(~truth, scales = "free") +
  xlab(expression(paste("ppHc - Without r")[g])) +
  ylab(expression(paste("ppHc - With r")[g])) +
  labs(color = "r<sub>g</sub>")  +
  geom_point(alpha = 0.8, aes(color = rg), size = 2.5) +
  geom_abline(color = "grey66") + theme_minimal()  %+replace%
  theme(text = element_text(size = 13),
        legend.text = element_text(size = 12),
        legend.title = element_markdown(size = 12),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13, angle = 90, hjust = 0.5),
        strip.text.x = element_textbox(size = 12, color = "black",
        fill = "#e6e6e6ff", box.color = "#e6e6e6ff",
        halign =  0.5, linetype = 1, r = unit(5, "pt"),
        width = unit(1, "npc"),
        padding = margin(2, 0, 1, 0),
        margin = margin(3, 3, 3, 3))) +
  scale_color_gradient(low = "grey", high = "purple", na.value = NA) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
  labels = c("0", "0.25", "0.5", "0.75", "1")) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
    labels = c("0", "0.25", "0.5", "0.75", "1"))
  g3a

#### Figure 3b
load(file = "data/sim/fig_data/fig3b_fdr.RData")
threeb_lab <- c(
  expression(paste("ABF|No.r")[g]), expression(paste("SuSIE.BF|No.r")[g]),
  expression(paste("ABF|With.r")[g]), expression(paste("SuSIE.BF|With.r")[g]))
names(threeb_lab) <- unique(fig3b_df_fdr$`BF|Covariate`)
fig3b_df_fdr$`BF|Covariate` <- factor(fig3b_df_fdr$`BF|Covariate`,
  levels = unique(fig3b_df_fdr$`BF|Covariate`))
ci_fig3b <- Hmisc::binconf(fig3b_df_fdr$FP, fig3b_df_fdr$TP +
  fig3b_df_fdr$FP, method = "wilson")
fig3b_df_fdr <- cbind(fig3b_df_fdr, ci_fig3b)
g3b <- ggplot(fig3b_df_fdr, aes(true.perc.Hc, FDR)) +
  geom_hline(yintercept = 0.05, linetype = "dashed", col = "grey55") +
  geom_vline(xintercept = unique(fig3b_df_fdr$true.num.Hc),
  linetype = "dotted", col = "grey66") +
  geom_pointrange(mapping = aes(ymin = Lower, ymax = Upper,
    shape = Priors, fill = `BF|Covariate`,
    color = `BF|Covariate`), size = 0.8, alpha = 0.85,
    linewidth = 1.1) +
  theme_bw() %+replace%
  theme(text = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    axis.text = element_text(size = 14)) +
  scale_color_pander(labels = threeb_lab) +
  scale_fill_pander(labels = threeb_lab) +
  xlab("True Hc %") +
  ylab("FDR") +
  scale_x_continuous(breaks = seq(0, 6, 1), limits = c(0, 6)) +
  scale_y_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1))
g3b

##### Combine 3a and 3b
g3a / g3b +
  plot_layout(heights = c(1, 1.4)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 12))

ggsave("figures/final_figures/fig3.pdf", width = 9, height = 10, units = "in")
ggsave("figures/final_figures/fig3.png",
  width = 9, height = 10, units = "in",
  bg = "white")

####################################################################### 0
### Supplementary figure 2
### Distribution of rg values in simulated data
load("data/sim/sim_data//sim_data.RData")

suppf2 <- simdata_sing %>%
mutate(truth_comb = case_when(
    grepl("Hn", truth)  ~ "Hn",
    grepl("Ha", truth) ~ "Ha",
    grepl("Hc", truth) ~ "Hc")) %>%
    mutate(truth_comb = factor(truth_comb, levels = c("Hn", "Ha", "Hc"),
     labels = c("True~Hn", "True~Ha+Ha2", "True~Hc+Hc2"))) %>%
     ggplot(aes(x = rg)) +
     facet_wrap(~truth_comb, nrow = 2, labeller = label_parsed,
     ncol = 3, scales = "free_y") +
     geom_histogram(alpha = 0.7, binwidth = 0.05,
     fill = "grey66", col = "grey44") +
         theme_bw() +
         xlab(expression(paste("Simulated r")[g])) +
         ylab("Trait counts") +
         theme(text = element_text(size = 14), strip.placement = "outside",
             axis.text = element_text(size = 16),
             axis.title = element_text(size = 14),
             legend.text = element_text(size = 14),
             legend.title = element_text(size = 14),
             strip.text.x = element_text(size = 14, color = "black"),
             strip.background.x = element_rect(colour = "black",
             fill = "white"))
suppf2

ggsave("figures/final_figures/suppfig2.pdf",
    bg = "white", width = 12, height = 5, units = "in")
ggsave("figures/final_figures/suppfig2.png",
    bg = "white", width = 12, height = 5, units = "in")

#######################################################################
### Supplementary Figure 3
source("sim05_run_diagnostics.R")

#######################################################################
### Supplementary Figure 4
load(file = "data/sim/fig_data/suppfig4df.RData")

suppfig4df  %>% ggplot(aes(Hc.cutoff, FDR, color = label)) +
  geom_point(size = 3, alpha = 0.6) +
  facet_wrap(~type, nrow = 2, scales = "free_x") +
  theme_bw() %+replace%
  theme(strip.placement = "outside", axis.text = element_text(size = 14),
  axis.title = element_text(size = 14), legend.text = element_text(size = 14),
  legend.title = element_text(size = 14), strip.text.x = element_text(size = 14,
  color = "black"), strip.background.x = element_rect(colour = "black",
  fill = "white"), legend.position = "bottom") +
  scale_y_continuous(limits = c(0, 1)) +
    geom_hline(yintercept = 0.05, linetype = "dashed", col = "grey44") +
    ylab("FDR") +
    scale_color_colorblind() + guides(color = guide_legend(title = ""))

ggsave("figures/final_figures/suppfig4.pdf",
       bg = "white", width = 14, height = 9, units = "in")
ggsave("figures/final_figures/suppfig4.png",
       bg = "white", width = 14, height = 9, units = "in")

#######################################################################
