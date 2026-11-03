library(tidyverse)
library(dplyr)
library(ggpointdensity)
library(viridis)
library(multcomp)
library(reshape2)
library(scales)
library(mgcv)


theme_pie <-   theme_classic(base_size = 7, base_family = "Helvetica") +
  theme(strip.background = element_blank(),
        axis.text.y = element_text(color = "black"),
        axis.text.x = element_text(color = "black"),
        text = element_text(family = "Helvetica"),
        legend.text=element_text(size=7),
        legend.key.size = unit(0.3, "cm"),
        plot.background = element_blank())

colours_dip <- c("#F9B233","#BE1622","#1D71B8","#A0D2F3","#006633","#7BA68F","#82368C","#AC8BBF","#F39200","#F9CDAB","#B2B2B2", "#4C5454","#E6007E", "#DEDC00") #Deletion # Homology Arm, # guides, 


path = "/Users/jw38/repos/dipper/"
outdir = paste0(path, "plots_figures/")

# Features (and add some context features)
data <- read_csv("/Users/jw38/repos/dipper/files/51-14_counts_featurized.csv") 
data2 <- read_csv("/Users/jw38/repos/dipper/files/51-14_counts_featurized_core.csv") 


########################################## Figure 2

# Length plot (not genomic position)
long_df <- data2 %>%
  dplyr::select("LFC_S2_HAP_DPr04", "LFC_S2_HAP_DPr08", "LFC_S2_HAP_BGH01","length_dpr04", "length_dpr08", "length_bgh01_hap") %>%
  pivot_longer(
    cols = starts_with("length"),
    names_to = "category",
    values_to = "length") %>%
  pivot_longer(
    cols = starts_with("LFC"),
    names_to = "data_category",
    values_to = "data") %>%
  mutate(
    data_category = tolower(data_category),
    category_suffix = tolower(str_replace(str_extract(category, "(dpr04|dpr08|bgh01_hap)$"), "bgh01_hap", "bgh01")),
    data_category_suffix = tolower(str_extract(data_category, "(dpr04|dpr08|bgh01)$"))
  ) %>%
  filter(category_suffix == data_category_suffix)

long_df2 <- data2 %>%
  dplyr::select("LFC_S2_HAP_DPr04", "LFC_S2_HAP_DPr08", "LFC_S2_HAP_BGH01","length_dpr04") %>%
  pivot_longer(
    cols = starts_with("length"),
    names_to = "category",
    values_to = "length") %>%
  pivot_longer(
    cols = starts_with("LFC"),
    names_to = "data_category",
    values_to = "data") %>%
  mutate(
    data_category = tolower(data_category),
    category_suffix = tolower(str_replace(str_extract(category, "(dpr04|dpr08|bgh01_hap)$"), "bgh01_hap", "bgh01")),
    data_category_suffix = tolower(str_extract(data_category, "(dpr04|dpr08|bgh01)$"))
  )


gam_model <- gam(data ~ exp(-length/1000), data = long_df, family = gaussian())
summary(gam_model)
gam_model2 <- gam(data ~ exp(-length), data = long_df2, family = gaussian())
summary(gam_model2)

linear_model <- lm(data ~ length, data = long_df2)
summary(linear_model)


ggplot(long_df %>% filter(length < 10000), aes(x = length/1000 , y = data)) +
  geom_point(shape = 16, alpha = 0.8, size = 0.5) +
  geom_smooth(method = "gam", formula = y ~ exp(-x), se = TRUE) +  # Adding a GAM fit
  labs(x = "Deletion length [kb]", y = "Log2 fold change") +
  theme_pie

ggplot(long_df2 %>% filter(length < 10000), aes(x = length/1000 , y = data)) +
  geom_point(shape = 16, alpha = 0.8, size = 0.5) +
  geom_smooth(method = "lm", se = TRUE) +  # Adding a GAM fit
  #geom_smooth(method = "gam", formula = y ~ exp(-x), se = TRUE) +  # Adding a GAM fit
  labs(x = "Deletion length [kb]", y = "Log2 fold change") +
  theme_pie







p <- ggplot(data2 %>% filter(length_dpr04 > -100), aes(x = length_dpr04/1000 , y = L2FC_S2_HAP_DPr04)) +
  geom_point(shape = 16, color = "#000000", alpha = 0.4, size = 0.2) +
  labs(x = "Deletion length [kb]", y = "Log2 fold change") +
  scale_x_log10(labels = c("0.1", "1", "10", "100", "1000"), limits = c(0.01, 1500), breaks = c(0.1, 1, 10, 100, 1000), expand = c(0,0)) + 
  coord_cartesian(ylim = c(-4, 7.5), expand = 0) +
  theme_pie
ggsave(paste0(outdir, 'Fig2_length_l2fc_log.pdf'), p, width = 9, height = 4, units = "cm", dpi = 600, useDingbats = F)


p <- ggplot(data2, aes(x = length_dpr04/1000 , y = L2FC)) +
  geom_point(shape = 16, color = "#000000", alpha = 0.4, size = 0.2) +
  labs(x = "Deletion length [kb]", y = "Log2 fold change") +
  coord_cartesian(xlim = c(-0.5, 10), ylim = c(-2.5, 5), expand = 0) +
  theme_pie
ggsave(paste0(outdir, 'Fig2_length_lfc_s1_10kb.pdf'), p, width = 4, height = 4, units = "cm", dpi = 600, useDingbats = F)

## Does it correlate more with length or with genomic position?
p <- ggplot(data2 %>% filter(length_dpr04 > -100), aes(x = start_genomic , y = L2FC_S2_HAP_DPr04)) +
  geom_point(shape = 16, color = "#000000", alpha = 0.4, size = 0.2) +
  labs(x = "Deletion length [kb]", y = "Log2 fold change") +
  scale_x_log10(labels = c("0.1", "1", "10", "100", "1000"), limits = c(0.01, 1500), breaks = c(0.1, 1, 10, 100, 1000), expand = c(0,0)) + 
  coord_cartesian(ylim = c(-4, 7.5), expand = 0) +
  theme_pie
ggsave(paste0(outdir, 'Fig2_length_l2fc_log.pdf'), p, width = 9, height = 4, units = "cm", dpi = 600, useDingbats = F)




# mead editing rates per length bin
data_bins <- data2 %>%
  mutate(length_bin = cut(length_dpr04, 
                          breaks = c(-Inf, 0 ,1000, 4999, 9999, 99999, 499999, 999999, Inf),
                          labels = c("neg ctrl", "<1", "1-5", "5-10", "10-100", 
                                     "100-500", "500-1000", ">1000")))



data_bins %>%
  group_by(length_bin) %>%
  summarise(mean = mean(L2FC,na.rm = TRUE))

test <- data2 %>% filter(length_bin == "neg ctrl")












p <- ggplot(data_bins, aes(x=length_bin, y=LFC_s1)) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.1,  fill = "#F1EEE4") +
  labs( x="Length category [kb]", y="Log2 fold change") +
  coord_cartesian(xlim = c(0.5, 8.5), ylim = c(-2.5, 5), expand = 0) +
  theme_pie + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(paste0(outdir, 'Fig2_length_bins.pdf'), p, width = 4, height = 5, units = "cm", dpi = 600, useDingbats = F)

# PBS length
pbs <-  read_csv("/Users/jw38/repos/dipper/files/DP102_mergeddata2_pbs.csv")
df_wide_pbs <- pbs %>%
  select(length_pbs, LFC_s1, spc_scfld_rtt) %>% 
  pivot_wider(names_from = length_pbs, values_from = LFC_s1) 

p <- ggplot(df_wide_pbs, aes(x = `11`, y = `16`)) +
  geom_smooth(method = "lm", col = "grey", size = 0.1,  se = FALSE) +
  geom_point(shape = 16,alpha = 0.8, size = 0.3, color = "#1D71B8") +
  annotate(geom = "text", x = 4, y = -0.6, label = paste0("R = ", round(sqrt(summary(lm(`16` ~ `11`, df_wide_pbs))$adj.r.squared), 2)), size = 2) +
  annotate(geom = "text", x = 4, y = -1, label = paste0("y ~ ", round((summary(lm(`16` ~ `11`, df_wide_pbs))$coefficients[2]), 2), "x"), size = 2) +
  labs(x = "11nt PBS (Log2FC)", y = "16nt PBS (Log2FC)") +
  coord_cartesian(xlim = c(-1.5, 5), ylim = c(-1.5, 5), expand = 0) +
  theme_pie
ggsave(paste0(outdir, 'Fig2_pbs_scatter_1116.pdf'), p, width = 4, height = 4, units = "cm", dpi = 600, useDingbats = F)

p <- ggplot(df_wide_pbs, aes(x = `11`, y = `13`)) +
  geom_smooth(method = "lm", col = "grey", size = 0.1,  se = FALSE) +
  geom_point(shape = 16,alpha = 0.8, size = 0.3, color = "#1D71B8") +
  annotate(geom = "text", x = 4, y = -0.6, label = paste0("R = ", round(sqrt(summary(lm(`13` ~ `11`, df_wide_pbs))$adj.r.squared), 2)), size = 2) +
  annotate(geom = "text", x = 4, y = -1, label = paste0("y ~ ", round((summary(lm(`13` ~ `11`, df_wide_pbs))$coefficients[2]), 2), "x"), size = 2) +
  labs(x = "11nt PBS (Log2FC)", y = "13nt PBS (Log2FC)") +
  coord_cartesian(xlim = c(-1.5, 5), ylim = c(-1.5, 5), expand = 0) +
  theme_pie
ggsave(paste0(outdir, 'Fig2_pbs_scatter_1113.pdf'), p, width = 4, height = 4, units = "cm", dpi = 600, useDingbats = F)

p <- ggplot(df_wide_pbs, aes(x = `11`, y = `9`)) +
  geom_smooth(method = "lm", col = "grey", size = 0.1,  se = FALSE) +
  geom_point(shape = 16,alpha = 0.8, size = 0.3, color = "#1D71B8") +
  annotate(geom = "text", x = 4, y = -0.6, label = paste0("R = ", round(sqrt(summary(lm(`9` ~ `11`, df_wide_pbs))$adj.r.squared), 2)), size = 2) +
  annotate(geom = "text", x = 4, y = -1, label = paste0("y ~ ", round((summary(lm(`9` ~ `11`, df_wide_pbs))$coefficients[2]), 2), "x"), size = 2) +
  labs(x = "11nt PBS (Log2FC)", y = "9nt PBS (Log2FC)") +
  coord_cartesian(xlim = c(-1.5, 5), ylim = c(-1.5, 5), expand = 0) +
  theme_pie 
ggsave(paste0(outdir, 'Fig2_pbs_scatter_119.pdf'), p, width = 4, height = 4, units = "cm", dpi = 600, useDingbats = F)


p <- ggplot(pbs, aes(x=as.factor(length_pbs), y=LFC_s1, color = length_pbs)) +
  geom_boxplot(outlier.size=0.05, color = "black",  linewidth = 0.2, fill = "#F1EEE4") +
  #scale_color_gradient(low = "#A0D2F3", high = "#144D7B") +
  labs(x="PBS length", y="Log2 fold change") +
  coord_cartesian(xlim = c(0.5, 4.5), ylim = c(-2.5, 5), expand = 0) +
  theme_pie + theme(legend.position="none")
ggsave(paste0(outdir, 'Fig2_pbs_boxplot.pdf'), p, width = 4, height = 4, units = "cm", dpi = 600, useDingbats = F)

# Scaffolds
scaffolds <-  read_csv("/Users/jw38/repos/dipper/files/DP102_mergeddata2_scaffoldsubset.csv")
df_wide_scaffolds <- scaffolds %>%
  dplyr::select(scaffold_name, LFC_s1, spc_ext) %>% 
  pivot_wider(names_from = scaffold_name, values_from = LFC_s1) 

p <- ggplot(df_wide_scaffolds, aes(x = improved, y = paste)) +
  geom_abline(intercept = 0, slope = 1, col = "grey", size = 0.1) +
  geom_point(shape = 16,alpha = 0.8, size = 0.3, color = '#4C5454') +
  #annotate(geom = "text", x = 3.5, y = -1.5, label = paste0("R = ", round(sqrt(summary(lm(paste ~ improved, df_wide_scaffolds))$adj.r.squared), 2)), size = 2) +
  annotate(geom = "text", x = 3.5, y = -1.5, label = paste0("y ~ x") )+
  labs(x = "Improved scaffold", y = "PASTE scaffold") +
  coord_cartesian(xlim = c(-2.5, 5), ylim = c(-2.5, 5), expand = 0) +
  theme_pie
ggsave(paste0(outdir, 'Fig2_scaffold_scatter.pdf'), p, width = 4, height = 4, units = "cm", dpi = 600, useDingbats = F)

p <- ggplot(df_wide_scaffolds, aes(x = paste, y = improved)) +
  geom_smooth(method = "lm", col = "grey", size = 0.1,  se = FALSE) +
  geom_point(shape = 16,alpha = 0.8, size = 0.3, color = '#4C5454') +
  annotate(geom = "text", x = 3.5, y = -1.5, label = paste0("R = ", round(sqrt(summary(lm(improved ~ paste, df_wide_scaffolds))$adj.r.squared), 2)), size = 2) +
  annotate(geom = "text", x = 3.5, y = -2, label = paste0("y ~ ", round(coef(lm(improved ~ paste, df_wide_scaffolds))[2], 2), "x"), size = 2) +
  labs(x = "PASTE scaffold", y = "Improved scaffold") +
  coord_cartesian(xlim = c(-2.5, 5), ylim = c(-2.5, 5), expand = 0) +
  theme_pie

p <- ggplot(scaffolds, aes(x=scaffold_name, y=LFC_s1)) +
  geom_boxplot(outlier.size=0.1,  linewidth = 0.1) +
  labs( x="scaffold", y="LFC_s1") +
  theme_pie 
ggsave(paste0(outdir, 'Fig2_scaffold_boxplot.pdf'), p, width = 4, height = 4, units = "cm", dpi = 600, useDingbats = F)


# Overhangs
overhangs <-  read_csv("/Users/jw38/repos/dipper/files/DP102_mergeddata2_overhangs.csv") %>%
  filter(ha_name != 'mismatch1') %>%
  filter(ha_name != 'mismatch2') %>%
  mutate(ha_name2 = ifelse(ha_name == "twinPE_34", "34nt", ifelse(ha_name == "twinPE_29", "29nt", "38nt")) ) 

df_wide_rtt <- overhangs %>%
  select(length, LFC_s1, spc_scfld_pbs, ha_name) %>% 
  pivot_wider(names_from = ha_name, values_from = LFC_s1) %>%
  drop_na()

keep_overhangs <- df_wide_rtt$spc_scfld_pbs
overhangs <- overhangs %>%
  filter(spc_scfld_pbs %in% keep_overhangs)

p <- ggplot(overhangs, aes(x=ha_name2, y=LFC_s1)) +
  geom_boxplot(outlier.size=0.05, color = "black", linewidth = 0.2, fill = "#F1EEE4") +
  labs( x=" ", y="Log2 fold change") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  coord_cartesian(xlim = c(0.5, 3.5), ylim = c(-2.5, 5),  expand = c(0, 0)) +
  theme_pie 
ggsave(paste0(outdir, 'Fig2_overhangs_boxplot.pdf'), p, width = 4, height = 4, units = "cm", dpi = 600, useDingbats = F)

p <-ggplot(df_wide_rtt, aes(x = twinPE_34, y = twinPE_29)) +
  geom_smooth(method = "lm", col = "grey", size = 0.1,  se = FALSE) +
  geom_point(shape = 16,alpha = 0.8, size = 0.3) +
  annotate(geom = "text", x = 3, y = -0.75, label = paste0("R = ", round(sqrt(summary(lm(twinPE_29 ~ twinPE_34, df_wide_rtt))$adj.r.squared), 2)), size = 2) +
  annotate(geom = "text", x = 3, y = -1, label = paste0("y ~ ", round((summary(lm(twinPE_29 ~ twinPE_34, df_wide_rtt))$coefficients[2]), 2), "x"), size = 2) +
  labs(x = "twinPE_34", y = "twinPE_29") +
  coord_cartesian(xlim = c(-2, 4.5), ylim = c(-2, 4.5)) +
  theme_pie
ggsave(paste0(outdir, 'Fig2_overhangs_scatter_twin.pdf'), p, width = 10, height = 4, units = "cm", dpi = 600, useDingbats = F)

p <-ggplot(df_wide_rtt, aes(x = twinPE_34  , y = full_38nt)) +
  geom_smooth(method = "lm", col = "grey", size = 0.1,  se = FALSE) +
  geom_point(shape = 16,alpha = 0.8, size = 0.3) +
  annotate(geom = "text", x = 3, y = -0.75, label = paste0("R = ", round(sqrt(summary(lm(full_38nt ~ twinPE_34, df_wide_rtt))$adj.r.squared), 2)), size = 2) +
  annotate(geom = "text", x = 3, y = -1, label = paste0("y ~ ", round((summary(lm(full_38nt ~ twinPE_34, df_wide_rtt))$coefficients[2]), 2), "x"), size = 2) +
  labs(x = "twinPE_34", y = "full_38nt") +
  coord_cartesian(xlim = c(-2, 4.5), ylim = c(-2, 4.5)) +
  theme_pie
ggsave(paste0(outdir, 'Fig2_overhangs_scatter_twin34paste.pdf'), p, width = 10, height = 4, units = "cm", dpi = 600, useDingbats = F)

p <- ggplot(df_wide_rtt, aes(x = twinPE_34  , y = mismatch2)) +
  geom_smooth(method = "lm", col = "grey", size = 0.1,  se = FALSE) +
  geom_point(shape = 16,alpha = 0.8, size = 0.3) +
  annotate(geom = "text", x = 3, y = -0.75, label = paste0("R = ", round(sqrt(summary(lm(mismatch2 ~ twinPE_34, df_wide_rtt))$adj.r.squared), 2)), size = 2) +
  annotate(geom = "text", x = 3, y = -1, label = paste0("y ~ ", round((summary(lm(mismatch2 ~ twinPE_34, df_wide_rtt))$coefficients[2]), 2), "x"), size = 2) +
  labs(x = "twinPE_34", y = "mismatch2") +
  coord_cartesian(xlim = c(-2, 4.5), ylim = c(-2, 4.5)) +
  theme_pie
ggsave(paste0(outdir, 'Fig2_overhangs_scatter_twin34mismatch2.pdf'), p, width = 10, height = 4, units = "cm", dpi = 600, useDingbats = F)



# PAM
p <- ggplot(data2%>% filter(length_dpr04 < 10000), aes(x=pam, y=LFC)) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.2, fill = "#F1EEE4") +
  labs( x="PAM", y="Log2 fold change") +
  coord_cartesian(xlim = c(0.5, 4.5), ylim = c(-2.5, 5),  expand = 0) +
  theme_pie 
ggsave(paste0(outdir, 'Fig2_pam_boxplot.pdf'), p, width = 4, height = 4, units = "cm", dpi = 600, useDingbats = F)

p <- ggplot(data2%>% filter(length_dpr04 < 10000), aes(x=pam, y=lenres)) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.2, fill = "#F1EEE4") +
  labs( x="PAM", y="Length residual") +
  coord_cartesian(xlim = c(0.5, 4.5), ylim = c(-2.5, 2.5),  expand = 0) +
  theme_pie 
ggsave(paste0(outdir, 'Fig2_pam_boxplot.pdf'), p, width = 4, height = 4, units = "cm", dpi = 600, useDingbats = F)


######## Spacer properties 
#####(How can we regress out properties that are generally good for prime editing vs properties that are particularly good for twinPE?)

# %GC
p <- ggplot(data2 %>% filter(length_dpr04 < 100000), aes(x = percGC_spc, y = LFC)) +
  geom_point(shape = 16, color = "#144D7B", alpha = 0.6, size = 0.5) +
  geom_smooth(method = "lm", col = "#545454", size = 0.1) +
  scale_color_viridis() +
  annotate(geom = "text", x = 0.8, y = -2, label = paste0("R = ", round(sqrt(summary(lm(percGC_spc ~ LFC, data2 %>% filter(length_dpr04 < 100000)))$adj.r.squared), 2)), size = 2) +
  labs(x = "GC content of spacer [%]", y = "Log2 fold change") +
  coord_cartesian(xlim = c(0, 1), ylim = c(-2.5, 5), expand = 0) +
  theme_pie
ggsave(paste0(outdir, 'Fig2_spacer_GC.pdf'), p, width = 4, height = 4, units = "cm", dpi = 600, useDingbats = F)

ggplot(data2 %>% filter(length_dpr04 < 100000), aes(x = percGC_spc, y = lenres_mean)) +
  geom_point(shape = 16, color = "#144D7B", alpha = 0.6, size = 0.5) +
  geom_smooth(method = "lm", col = "#545454", size = 0.1) +
  scale_color_viridis() +
  annotate(geom = "text", x = 0.8, y = -2, label = paste0("R = ", round(sqrt(summary(lm(percGC_spc ~ lenres_mean, data2 %>% filter(length_dpr04 < 100000)))$adj.r.squared), 2)), size = 2) +
  labs(x = "GC content of spacer [%]", y = "Length residual") +
  coord_cartesian(xlim = c(0, 1), ylim = c(-2.5, 5), expand = 0) +
  theme_pie

p <- ggplot(data2 %>% filter(length_dpr04 < 100000, length_dpr04 > 0), aes(x = percGC_pbs, y = lenres_mean)) +
  geom_pointdensity(shape = 16, size = 0.1, show.legend = FALSE) +
  geom_smooth(method = "lm", col = "#545454", size = 0.1) +
  scale_color_viridis() +
  annotate(geom = "text", x = 0.8, y = -3, label = paste0("R = ", round(sqrt(summary(lm(percGC_pbs ~ lenres_mean, data2 %>% filter(length_dpr04 < 100000)))$adj.r.squared), 2)), size = 2) +
  labs(x = "GC content of PBS [%]", y = "Length residual mean") +
  coord_cartesian(xlim = c(0, 1), ylim = c(-4, 4)) +
  theme_pie
ggsave(paste0(outdir, 'Fig2_pbs_GC.pdf'), p, width = 4, height = 4, units = "cm", dpi = 600, useDingbats = F)



# nucleotides
cols_context <- c("mN17", "mN16", "mN15", "mN14", "mN13", "mN12", "mN11", "mN10", 
                  "mN9", "mN8", "mN7", "mN6", "mN5", "mN4","mN3", "mN2", "mN1", 
                  "N1", "N2", "N3", "N4", "N5", "N6", "N7", "N8", "N9", "N10",
                  "N11", "N12", "N13", "N14", "N15", "N16", "N17", "N18", "N19", "N20"
                  #  "N21", "N22", "N23", "N24", "N25", "N26", "N27", "N28", "N29", "N30",
                  #  "N31", "N32", "N33", "N34", "N35", "N36", "N37", "N38", "N39", "N40",
                  #  "N41", "N42", "N43", "N44", "N45", "N46", "N47", "N48", "N49", "N50"
)

melted_data <- data2 %>% 
  #filter(length_dpr04 < 100000) %>%
  #filter(!is.na(N30))%>%
  pivot_longer(cols = all_of(cols_context), names_to = "position", values_to = "nucleotide") %>% 
  group_by(position, nucleotide) %>% 
  summarise(mean_value = mean(lenres_mean, na.rm = TRUE)) 


p <- ggplot(melted_data, aes(x = factor(position, levels = cols_context), y = nucleotide, fill = mean_value)) +
  geom_tile() +
  scale_fill_gradient2(low = "#1D71B8", mid = "white", high = "#BE1622", midpoint = 0) +
  labs(x = "Position relative to nick site", 
       y = "Nucleotide", 
       fill = "Mean lenresidual") +
  theme_pie
ggsave(paste0(outdir, 'Fig2_nucleotides_heatmap.pdf'), p, width = 16, height = 2.5, units = "cm", dpi = 600, useDingbats = F)


# other lengths
# look up means from main data set for lenres mutate(lenres_mean = LFC_s1 - mean_LFC_s1, lenres = LFC_s1 - median_LFC_s1)
mean_values <- data2 %>%
  group_by(target) %>%
  summarize(mean_LFC_s1 = mean(mean_LFC_s1, na.rm = TRUE))

melted_data2 <- overhangs %>% 
  left_join(mean_values, by = "target") %>%
  mutate(lenres_mean = LFC_s1 - mean_LFC_s1) %>%
  filter(scaffold_name == "improved") %>%
  filter(ha_name == 'full_38nt') %>%
  filter(length < 100000) %>%
  mutate( mN4 = substr(target, start = 14, stop = 14),
          mN3 = substr(target, start = 15, stop = 15),
          mN2 = substr(target, start = 16, stop = 16),
          mN1 = substr(target, start = 17, stop = 17),
          N1 = substr(target, start = 18, stop = 18),
          N2 = substr(target, start = 19, stop = 19),
          N3 = substr(target, start = 20, stop = 20),
          N4 = substr(target, start = 21, stop = 21),
          N5 = substr(target, start = 22, stop = 22),
          N6 = substr(target, start = 23, stop = 23)
  )  %>%
  pivot_longer(cols = all_of(cols_context2), names_to = "position", values_to = "nucleotide") %>% 
  group_by(position, nucleotide) %>% 
  summarise(mean_value = mean(lenres_mean, na.rm = TRUE)) 

cols_context2 <- c("mN4","mN3", "mN2", "mN1", 
                   "N1", "N2", "N3", "N4", "N5", "N6"
)

p <- ggplot(melted_data2, aes(x = factor(position, levels = cols_context2), y = nucleotide, fill = mean_value)) +
  geom_tile() +
  scale_fill_gradient2(low = "#1D71B8", mid = "white", high = "#BE1622", midpoint = 0) +
  labs(x = "Position relative to nick site", 
       y = "Nucleotide", 
       fill = "Mean lenres_mean") +
  theme_pie
ggsave(paste0(outdir, 'Fig2_nucleotides_heatmap_29.pdf'), p, width = 8, height = 2.5, units = "cm", dpi = 600, useDingbats = F)




####Other properties
p <- ggplot(data2 %>% filter(length_abs < 10000), aes(x=pos1compl, y=LFC)) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.2, fill = "#F1EEE4") +
  labs(x="pos1compl", y="Log2 fold change") +
  coord_cartesian(xlim = c(0.5, 2.5), ylim = c(-2.5, 5),  expand = 0) +
  theme_pie 

ggplot(data2 %>% filter(length_abs < 10000), aes(x=pos1compl, y=lenres)) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.2, fill = "#F1EEE4") +
  labs(x="pos1compl", y="Log2 fold change") +
  coord_cartesian(xlim = c(0.5, 2.5), ylim = c(-2.5, 5),  expand = 0) +
  theme_pie 

ggsave(paste0(outdir, 'Fig2_pos1compl_boxplot.pdf'), p, width = 4, height = 4, units = "cm", dpi = 600, useDingbats = F)

p <- ggplot(data %>% filter(length_abs < 10000), aes(x=VF_spc, y=LFC)) +
  geom_point() +
  labs( x="VF spacer", y="LFC_s1") +
  theme_pie 
ggsave(paste0(outdir, 'Fig2_pos1compl_boxplot.pdf'), p, width = 4, height = 4, units = "cm", dpi = 600, useDingbats = F)

######### Is length and guide effect all that counts? How much are we missing?

library(caret)
data2_clean <- data2 %>% filter(!is.na(L2FC)) %>%
  filter(length_dpr04 > 10000 & length_dpr04 < 100000 )
set.seed(123)
train_indices <- sample(seq_len(nrow(data2_clean)), size = 0.8 * nrow(data2_clean))

# Split the data into training and test sets based on the sampled indices
train_set <- data2_clean[train_indices, ]
test_set <- data2_clean[-train_indices, ]

model_ideal <- lm(L2FC ~ length_dpr04 + pridict_score + DeepCas9, train_set)
predictions_ideal <- predict(model_ideal, newdata = test_set)
test_set$pred_ideal <- predictions_ideal


p <- ggplot(test_set, aes(x = L2FC, y = pred_ideal)) +
  geom_point(shape = 16, color = "#35597A", alpha = 0.8, size = 0.3) +
  geom_smooth(method = "lm", col = "black", size = 0.5) +
  annotate(geom = "text", x = 4, y = 1, label = paste0("R = ", round(sqrt(summary(lm(pred_ideal ~ L2FC, test_set))$adj.r.squared), 2), size = 1)) +
  labs(x = "Measured L2FC", y = "Predicted L2FC") +
  theme_pie
ggsave(paste0(outdir, 'prediction_testset.pdf'), p, width = 4, height = 4, units = "cm", dpi = 600, useDingbats = F)


rsqr_10kb <- sqrt(summary(lm(L2FC ~ pred_ideal, test_set))$adj.r.squared) # 694
rsqr_100kb <- sqrt(summary(lm(L2FC ~ pred_ideal, test_set))$adj.r.squared) # 28
rsqr_1000kb <- sqrt(summary(lm(L2FC ~ pred_ideal, test_set))$adj.r.squared) # 17


Models <- data.frame(Rsq = c(rsqr_10kb, rsqr_100kb, rsqr_1000kb),
                     label = c('0 - 10kb', '10 - 100 kb', '100 - 1000 kb'),
                     order = c(1, 2, 3))

p <- ggplot(Models, aes(x = reorder(label, order), y = Rsq)) +
  geom_bar(stat="identity") +
  labs(x = "Deletion length subset", y = "Variance explained in test set [R2]") +
  coord_cartesian(ylim = c(0, 1)) +
  theme_pie

ggsave(paste0(outdir, 'ExplainingVariation.pdf'), p, width = 5, height = 4, units = "cm", dpi = 600, useDingbats = F)




######### Epigenetic marks:
# normalize by spacer efficiency?

ggplot(data2 %>% filter(H3K27me3_ENCFF139KZL_average_value_100 > 0), aes(x = H3K27me3_ENCFF139KZL_average_value_100, y = L2FC)) +
  geom_point(shape = 16, color = "#000000", alpha = 1, size = 0.5) +
  geom_smooth(method = "lm", col = "#FFD700", size = 0.1,  se = FALSE) +
  annotate(geom = "text", x = 0.1, y = -2, label = paste0("R = ", round(sqrt(summary(lm(L2FC ~ H3K27me3_ENCFF139KZL_average_value_100, data2 %>% filter(H3K27me3_ENCFF139KZL_average_value_100 >0)))$adj.r.squared), 3)), size = 2) +
  theme_pie

ggplot(data2 %>% filter(length_dpr04 < 50000), aes(x = `Dnase-seq_ENCFF972GVB_average_value_100`, y = L2FC)) +
  geom_point(shape = 16, color = "#000000", alpha = 1, size = 0.5) +
  geom_smooth(method = "lm", col = "#FFD700", size = 0.1,  se = FALSE) +
  annotate(geom = "text", x = 0.1, y = -2, label = paste0("R = ", round(sqrt(summary(lm(L2FC ~ `Dnase-seq_ENCFF972GVB_average_value_100`, data2 %>% filter(length_dpr04 < 1000)))$adj.r.squared), 3)), size = 2) +
  theme_pie


ggplot(data2, aes(x = HEK293T_H3K9me3_signalC, y = lenres)) +
  geom_point(shape = 16, color = "#000000", alpha = 1, size = 0.5) +
  geom_smooth(method = "lm", col = "#FFD700", size = 0.1,  se = FALSE) +
  annotate(geom = "text", x = 4, y = -2, label = paste0("R = ", round(sqrt(summary(lm(lenres ~ HEK293T_H3K9me3_signalC, data2))$adj.r.squared), 3)), size = 2) +
  theme_pie

ggplot(data2, aes(x = HEK293T_H3K27me3_mean, y = lenres)) +
  geom_point(shape = 16, color = "#000000", alpha = 1, size = 0.5) +
  geom_smooth(method = "lm", col = "#FFD700", size = 0.1,  se = FALSE) +
  annotate(geom = "text", x = 4, y = -0.1, label = paste0("R = ", round(sqrt(summary(lm(lenres ~ HEK293T_H3K27me3_mean, data2))$adj.r.squared), 3)), size = 2) +
  theme_pie

ggplot(data2 %>% filter(HEK293T_CTCF_signalC > 0 & HEK293T_CTCF_signalC < 10), aes(x = HEK293T_CTCF_signalC, y = lenres)) +
  geom_point(shape = 16, color = "#000000", alpha = 1, size = 0.5) +
  geom_smooth(method = "lm", col = "#FFD700", size = 0.1,  se = FALSE) +
  annotate(geom = "text", x = 4, y = -2, label = paste0("R = ", round(sqrt(summary(lm(lenres ~ HEK293T_CTCF_mean, data2))$adj.r.squared), 3)), size = 2) +
  theme_pie

ggplot(data2, aes(x = HEK293T_DNase_signalC, y = lenres)) +
  geom_point(shape = 16, color = "#000000", alpha = 1, size = 0.5) +
  geom_smooth(method = "lm", col = "#FFD700", size = 0.1,  se = FALSE) +
  annotate(geom = "text", x = 4, y = -2, label = paste0("R = ", round(sqrt(summary(lm(lenres ~ HEK293T_DNase_mean, data2))$adj.r.squared), 3)), size = 2) +
  theme_pie

ggplot(data2 %>% filter(HEK293T_DNase_signalC < 10), aes(x = HEK293T_DNase_signalC, y = lenres)) +
  geom_point(shape = 16, color = "#000000", alpha = 1, size = 0.5) +
  geom_smooth(method = "lm", col = "#FFD700", size = 0.1,  se = FALSE) +
  annotate(geom = "text", x = 4, y = -2, label = paste0("R = ", round(sqrt(summary(lm(lenres ~ HEK293T_DNase_mean, data2 %>% filter(HEK293T_DNase_signalC < 10)))$adj.r.squared), 3)), size = 2) +
  theme_pie

ggplot(data2 %>% filter(HEK293T_ATAC_signalC < 10), aes(x = HEK293T_DNase_signalC, y = lenres)) +
  geom_point(shape = 16, color = "#000000", alpha = 1, size = 0.5) +
  geom_smooth(method = "lm", col = "#FFD700", size = 0.1,  se = FALSE) +
  annotate(geom = "text", x = 4, y = -2, label = paste0("R = ", round(sqrt(summary(lm(lenres ~ HEK293T_DNase_mean, data2 %>% filter(HEK293T_DNase_signalC < 10)))$adj.r.squared), 3)), size = 2) +
  theme_pie

ggplot(data2 %>% filter(HEK293T_H3K27ac_mean < 20), aes(x = HEK293T_H3K27ac_mean, y = lenres)) +
  geom_point(shape = 16, color = "#000000", alpha = 1, size = 0.5) +
  geom_smooth(method = "lm", col = "#FFD700", size = 0.1,  se = FALSE) +
  annotate(geom = "text", x = 6, y = -1, label = paste0("R = ", round(sqrt(summary(lm(lenres ~ HEK293T_H3K27ac_mean, data2 %>% filter(HEK293T_H3K27ac_mean < 20)))$adj.r.squared), 3)), size = 2) +
  theme_pie

ggplot(data2, aes(x = HEK293T_H3K4me1_mean, y = lenres)) +
  geom_point(shape = 16, color = "#000000", alpha = 1, size = 0.5) +
  geom_smooth(method = "lm", col = "#FFD700", size = 0.1,  se = FALSE) +
  annotate(geom = "text", x = 4, y = -2, label = paste0("R = ", round(sqrt(summary(lm(lenres ~ HEK293T_H3K4me1_mean, data2))$adj.r.squared), 3)), size = 2) +
  theme_pie

########## ChromHMM: "emission_name_HEK"   "abbreviation_HEK" "emission_HEK"   
# Change some categories as the number per category where low
data2 <- data2 %>%
  mutate(emission_name_juliane = case_when(
    emission_name_HEK %in% c("Active TSS", "Flanking TSS", "Poised TSS", "Transcription strong", "Transcription weak","Gene body") ~ "Transcription",
    emission_name_HEK %in% c("Enhancer weak", "Genic Enhancer", "Strong enhancer") ~ "Enhancer",
    emission_name_HEK %in% c("Quiescent", "Quiescent Gene") ~ "Quiescent",
    TRUE ~ as.character(emission_name_HEK)
  ))

# Plot the observations per group
data2 %>%
  filter(length_dpr04 < 100000) %>%
  group_by(emission_name_juliane) %>%
  summarise(n = n()) %>%
  ggplot(aes(x = emission_name_juliane, y = n)) +
  geom_bar(stat = "identity", fill = "#56B4E9") +
  geom_text(aes(label = n), vjust = -0.5, size = 3) +
  labs(x = "ChromHMM label", y = "Number of Observations") +
  labs(x = "ChromHMM label", y = "Number of Observations") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Boxplot of ChromHMM

# Perform ANOVA and TukeyHSD
anova_res <- aov(lenres ~ emission_name_juliane, data = data2 %>% filter(length_dpr04 < 50000))
tukey_res <- TukeyHSD(anova_res)
p_values <- as.data.frame(tukey_res$emission_name_juliane)
p_values$comparison <- rownames(p_values)
significant_p_values <- p_values %>% filter(`p adj` < 0.05)

ggplot(data2 %>% filter(emission_name_juliane != "Repressed polycomb" & length_dpr04 < 50000), aes(x=emission_name_juliane, y=lenres)) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.1,  fill = "#F1EEE4") +
  geom_point() +
  labs(x="ChromHMM label", y="Lenres") +
  #geom_text(data = significant_p_values, aes(x = comparison, y = 20, label = sprintf("p = %.3f", `p adj`)), angle = 90, vjust = 1.5) +
  coord_cartesian(xlim = c(0.5, 5.5), ylim = c(-1.5, 1.5), expand = 0) +
  theme_pie + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot(data2 %>% filter(emission_name_juliane != "Repressed polycomb" & length_dpr04 < 50000), aes(x=emission_name_juliane, y=lenres)) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.1,  fill = "#F1EEE4") +
  labs(x="ChromHMM label", y="Lenres") +
  geom_text(data = significant_p_values, aes(x = comparison, y = 20, label = sprintf("p = %.3f", `p adj`)), angle = 90, vjust = 1.5) +
  coord_cartesian(xlim = c(0.5, 4.5), ylim = c(-1.5, 1.5), expand = 0) +
  theme_pie + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


########## Predictions

# Cas9 score
p <- ggplot(data2 %>% filter(., length_abs < 10000), aes(x = L2FC, y = DeepCas9)) +
  geom_point(shape = 16, color = "#000000", alpha = 0.7, size = 0.2) +
  geom_smooth(method = "lm", col = "#1D71B8", size = 0.2) +
  annotate(geom = "text", x = 4, y = 10, label = paste0("R = ", round(sqrt(summary(lm(L2FC ~ DeepCas9, data2%>% filter(., length_abs < 10000)))$adj.r.squared), 3)), size = 2) +
  labs(x = "L2FC", y = "DeepCas9", title = "DeepCas9 for deletions < 10kb") +
  #coord_cartesian(xlim = c(0.1, 100), ylim = c(0.1, 100)) +
  #geom_text(hjust = 0, nudge_x = 0.1, size = 2) +
  theme_pie
ggsave(paste0(outdir, 'Fig4_DeepCas9.pdf'), p, width = 4, height = 4, units = "cm", dpi = 600, useDingbats = F)
ggsave(paste0(outdir, 'Fig4_DeepCas9b.pdf'), p, width = 6, height = 6, units = "cm", dpi = 600, useDingbats = F)


p <- ggplot(data %>% filter(., length_abs < 2000), aes(x = LFC, y = DeepCas9)) +
  geom_point(shape = 16, color = "#000000", alpha = 1, size = 0.5) +
  geom_smooth(method = "lm", col = "#FFD700", size = 0.1,  se = FALSE) +
  annotate(geom = "text", x = 4, y = 10, label = paste0("R = ", round(sqrt(summary(lm(LFC ~ DeepCas9, data%>% filter(., length_abs < 2000)))$adj.r.squared), 3)), size = 2) +
  labs(x = "LFC Screen 1", y = "DeepCas9", title = "DeepCas9 for deletions < 2kb") +
  #coord_cartesian(xlim = c(0.1, 100), ylim = c(0.1, 100)) +
  #geom_text(hjust = 0, nudge_x = 0.1, size = 2) +
  theme_pie
ggsave(paste0(outdir, 'Fig4_DeepCas9_2kb.pdf'), p, width = 4, height = 4, units = "cm", dpi = 600, useDingbats = F)

# PRIDICT score
p <- ggplot(data %>% filter(., length_abs < 10000), aes(x = LFC_s1, y = pridict1_score)) +
  geom_point(shape = 16, color = "#000000", alpha = 1, size = 0.5) +
  geom_smooth(method = "lm", col = "#FFD700", size = 0.1,  se = FALSE) +
  annotate(geom = "text", x = 4, y = 1, label = paste0("R = ", round(sqrt(summary(lm(LFC_s1 ~ pridict1_score, data%>% filter(., length_abs < 10000)))$adj.r.squared), 3)), size = 2) +
  labs(x = "LFC Screen 1", y = "PRIDICT score", title = "PRIDICT score for deletions < 10kb using attB insertion") +
  theme_pie
ggsave(paste0(outdir, 'Fig4_scores_pridict.pdf'), p, width = 4, height = 4, units = "cm", dpi = 600, useDingbats = F)

p <- ggplot(data2 %>% filter(., length_abs < 10000), aes(x = L2FC, y = pridict_score)) +
  geom_point(shape = 16, color = "#000000", alpha = 0.7, size = 0.2) +
  geom_smooth(method = "lm", col = "#1D71B8", size = 0.2) +
  annotate(geom = "text", x = 4, y = 1, label = paste0("R = ", round(sqrt(summary(lm(L2FC ~ pridict_score, data2%>% filter(., length_abs < 10000)))$adj.r.squared), 3)), size = 2) +
  labs(x = "L2FC", y = "PRIDICT2 score", title = "PRIDICT2 score for deletions < 10kb using attB insertion") +
  theme_pie
ggsave(paste0(outdir, 'Fig4_scores_pridict.pdf'), p, width = 4, height = 4, units = "cm", dpi = 600, useDingbats = F)
ggsave(paste0(outdir, 'Fig4_scores_pridictb.pdf'), p, width = 6, height = 6, units = "cm", dpi = 600, useDingbats = F)

ggplot(data2 %>% filter(length_dpr04 < 10000), aes(x = lenres, y = pridict2_score_HEK)) +
  geom_point(shape = 16, color = "#000000", alpha = 1, size = 0.5) +
  geom_smooth(method = "lm", col = "#FFD700", size = 0.1,  se = FALSE) +
  annotate(geom = "text", x = 4, y = 1, label = paste0("R = ", round(sqrt(summary(lm(lenres ~ pridict2_score_HEK, data2%>% filter(., length_dpr04 < 10000)))$adj.r.squared), 3)), size = 2) +
  labs(x = "lenres", y = "PRIDICT2 score", title = "PRIDICT score for deletions < 10kb using attB insertion") +
  theme_pie

ggplot(data2 %>% filter(length_dpr04 < 10000), aes(x = LFC_s1, y = pridict2_score_HEK)) +
  geom_point(shape = 16, color = "#000000", alpha = 1, size = 0.5) +
  geom_smooth(method = "lm", col = "#FFD700", size = 0.1,  se = FALSE) +
  annotate(geom = "text", x = 4, y = 1, label = paste0("R = ", round(sqrt(summary(lm(LFC_s1 ~ pridict2_score_HEK, data2%>% filter(., length_dpr04 < 10000)))$adj.r.squared), 3)), size = 2) +
  labs(x = "LFC Screen 1", y = "PRIDICT2 score", title = "PRIDICT score for deletions < 10kb using attB insertion") +
  theme_pie

p <- ggplot(data %>% filter(., length_abs < 2000), aes(x = LFC_s1, y = pridict_score)) +
  geom_point(shape = 16, color = "#000000", alpha = 1, size = 0.5) +
  geom_smooth(method = "lm", col = "#FFD700", size = 0.1,  se = FALSE) +
  annotate(geom = "text", x = 4, y = 1, label = paste0("R = ", round(sqrt(summary(lm(LFC_s1 ~ pridict_score, data%>% filter(., length_abs < 2000)))$adj.r.squared), 3)), size = 2) +
  labs(x = "LFC Screen 1", y = "PRIDICT score", title = "PRIDICT score for deletions < 2kb") +
  theme_pie
ggsave(paste0(outdir, 'Fig4_scores_pridict_2kb.pdf'), p, width = 4, height = 4, units = "cm", dpi = 600, useDingbats = F)

# ePRIDICT score
p <- ggplot(data %>% filter(., length_abs < 100000), aes(x = LFC, y = `ePRIDICT score (light model)`)) +
  geom_point(shape = 16, color = "#000000", alpha = 1, size = 0.5) +
  geom_smooth(method = "lm", col = "#FFD700", size = 0.1,  se = FALSE) +
  annotate(geom = "text", x = 2, y = 50, label = paste0("R = ", round(sqrt(summary(lm(LFC ~ `ePRIDICT score (light model)`, data%>% filter(., length_abs < 100000)))$adj.r.squared), 3)), size = 2) +
  labs(x = "LFC Screen 1", y = "ePRIDICT score", title = "ePRIDICT score") +
  theme_pie
ggsave(paste0(outdir, 'Fig4_scores_epridict.pdf'), p, width = 4, height = 4, units = "cm", dpi = 600, useDingbats = F)


# DeepPrime score
p <- ggplot(data %>% filter(., length_abs < 10000), aes(x = LFC, y = score_deepprime)) +
  geom_point(shape = 16, color = "#000000", alpha = 1, size = 0.5) +
  geom_smooth(method = "lm", col = "#FFD700", size = 0.1,  se = FALSE) +
  annotate(geom = "text", x = 4, y = 1, label = paste0("R = ", round(sqrt(summary(lm(LFC ~ score_deepprime, data%>% filter(., length_abs < 10000)))$adj.r.squared), 3)), size = 2) +
  labs(x = "LFC Screen 1 (HEK293T)", y = "deepprime score", title = "deepprime score for deletions < 10kb using attB insertion") +
  theme_pie
ggsave(paste0(outdir, 'SuppFig3_scores_deepprime.pdf'), p, width = 4, height = 4, units = "cm", dpi = 600, useDingbats = F)

p <- ggplot(data %>% filter(., length_abs < 10000), aes(x = LFC_s1, y = Doench2014OnTarget)) +
  geom_point(shape = 16, color = "#000000", alpha = 1, size = 0.5) +
  geom_smooth(method = "lm", col = "#FFD700", size = 0.1,  se = FALSE) +
  annotate(geom = "text", x = 4, y = 1, label = paste0("R = ", round(sqrt(summary(lm(LFC_s1 ~ DeepPExCFD, data%>% filter(., length_abs < 10000)))$adj.r.squared), 3)), size = 2) +
  labs(x = "LFC Screen 1 (HEK293T)", y = "Doench 2014 On Target score", title = "Doench 2014 On Target score for deletions < 10kb") +
  #coord_cartesian(xlim = c(0.1, 100), ylim = c(0.1, 100)) +
  #geom_text(hjust = 0, nudge_x = 0.1, size = 2) +
  theme_pie
ggsave(paste0(outdir, 'SuppFig3_scores_Doench2014On.pdf'), p, width = 4, height = 4, units = "cm", dpi = 600, useDingbats = F)

p <- ggplot(data %>% filter(., length_abs < 10000), aes(x = LFC, y = DeepSpCas9)) +
  geom_point(shape = 16, color = "#000000", alpha = 1, size = 0.5) +
  geom_smooth(method = "lm", col = "#FFD700", size = 0.1,  se = FALSE) +
  annotate(geom = "text", x = 4, y = 10, label = paste0("R = ", round(sqrt(summary(lm(LFC ~ DeepSpCas9, data%>% filter(., length_abs < 10000)))$adj.r.squared), 3)), size = 2) +
  labs(x = "LFC Screen 1 (HEK293T)", y = "DeepSpCas9", title = "DeepSpCas9 for deletions < 10kb") +
  #coord_cartesian(xlim = c(0.1, 100), ylim = c(0.1, 100)) +
  #geom_text(hjust = 0, nudge_x = 0.1, size = 2) +
  theme_pie
ggsave(paste0(outdir, 'SuppFig3_scores_DeepSpCas9.pdf'), p, width = 3, height = 4, units = "cm", dpi = 600, useDingbats = F)

