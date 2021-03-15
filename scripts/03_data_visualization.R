

# Load packages

library(tidyverse)
library(cowplot)
library(ggforce)
library(ggtext)
library(dotwhisker)
library(DiagrammeR)
library(viridis)

#==============================================================================


# Import the raw differential expression data for the quantitative summary and
# generate derived data columns

d <- read_csv("data/quantitative_summary/wild_expression_data.csv") %>%
  mutate(
    total_tested = as.integer(total_tested),
    de_total_tested = as.integer(de_total_tested),
    total_annotated = as.integer(total_annotated),
    de_total_annotated = as.integer(de_total_annotated),
    prop_de = round(de_total_tested/total_tested, 2),
    prop_de_unrounded = de_total_tested/total_tested,
    log10_prop_de = log10(prop_de_unrounded),
    log10_prop_de = ifelse(log10_prop_de < -2, -2.5, log10_prop_de),
    log2_prop_de = log2(prop_de_unrounded),
    log2_prop_de = ifelse(log2_prop_de <= -8, -8, log2_prop_de),
    study_mod = paste0(
      study, "\n",
      "Assay: ", assay, "\n",
      "Host class: ", host_class, "\n",
      "Pathogen/parasite: ", pathogen
    ),
    study_mod_alt = paste0(
      study, "\n",
      "Assay: ", assay, ", ",
      "Host class: ", host_class, ", ",
      "Pathogen/parasite: ", pathogen
    ),
    study_mod_simple = paste0(
      study, "<br>",
      assay, "<br>",
      host_class, "<br>",
      "*", pathogen, "*"
    ),
    study = as.factor(study),
    tissue_mod = paste0("Tissue assayed: ", tissue),
    time_point_mod = paste0("Post-exposure time point: ", time_point),
    second_line = ifelse(
      is.na(pathogen_strain), 
      tissue, 
      paste(pathogen_strain, tissue, sep = " - ")
    ),
    second_line = paste(second_line, time_point, sep = " - "),
    group_de = paste(comparison, second_line, sep = "\n")
    %>%
      as.factor(),
    susceptibility_binary = ifelse(susceptibility == "susceptible", 1, 0)
  )

#==============================================================================


# Figure 1a code

# Generate labels that will appear along the bottom of the figure and add to
# the larger data frame

group.labels <- d %>%
  arrange(study_mod_simple, group_de) %>%
  filter(susceptibility == "susceptible") %>%
  group_by(study) %>%
  mutate(
    first_line = paste0(substr(study, 1, 2), as.character(row_number()))
  ) %>%
  ungroup() %>%
  mutate(
    time_point_short = str_replace_all(time_point, " ", ""),
    time_point_short = substr(time_point_short, 1, 3),
    time_point_short = 
      ifelse(str_detect(time_point_short, "opp"), "op", time_point_short),
    time_point_short =
      ifelse(str_detect(time_point_short, "atc"), "cs", time_point_short),
    tissue_short = substr(tissue, 1, 2),
    bottom_label = 
      paste0(first_line, "\n", time_point_short)
  ) %>%
  select(group_de, bottom_label)

d <- left_join(d, group.labels, by = "group_de")

# Generate the segment lengths that will connect disease-resistant and 
# disease-susceptible species comparisons

segment.limits <- d %>%
  group_by(bottom_label) %>%
  summarize(
    min_measure = min(prop_de_unrounded),
    max_measure = max(prop_de_unrounded)
  )

# Plot and save Figure 1a

fig1a <- d %>%
  left_join(., segment.limits, by = "bottom_label") %>%
  ggplot(
    aes(x = bottom_label, y = prop_de_unrounded, 
        color = susceptibility, 
        fill = susceptibility,
        shape = tissue)
  ) +
  geom_segment(
    aes(x = bottom_label, xend = bottom_label, 
        y = min_measure, yend = max_measure),
    color = "darkgrey", size = 2
  ) +
  geom_point(size = 8) +
  ylim(0, 0.6) +
  scale_color_manual(values = alpha(c("dodgerblue", "firebrick2"), 0.9)) +
  scale_fill_manual(values = alpha(c("dodgerblue", "firebrick2"), 0.9)) +
  scale_shape_manual(values = c(17, 15, 16, 18, 25)) +
  ylab("Proportion of genes/contigs/probes differentially expressed") +
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.text.x = element_text(angle = 0, size = 22),
    axis.text = element_text(size = 24),
    axis.title.y = element_text(size = 35),
    legend.text = element_text(size = 35),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.title = element_blank()
  ) +
  guides(color = guide_legend(order = 1)) +
  facet_wrap(~study_mod_simple, scales = "free_x", nrow = 1) +
  theme(strip.text.x = element_markdown()) +
  guides(
    fill = FALSE,
    shape = guide_legend(override.aes = list(fill = "black")) 
  )

fig1a

ggsave("outputs/fig1a.pdf", width = 30, height = 15)

# Output reference table

d %>%
  group_by(bottom_label) %>%
  distinct(
    study, assay, pathogen, host_class, comparison, 
    tissue, time_point, pathogen_strain) %>%
  select(bottom_label, everything()) %>%
  ungroup() %>%
  mutate(bottom_label = substr(bottom_label, 1, 3)) %>%
  arrange(bottom_label) %>%
  rename(
    `Comparison Label` = bottom_label,
    `Study` = study,
    `Assay` = assay,
    `Pathogen` = pathogen,
    `Host Class` = host_class,
    `Host Comparison` = comparison,
    `Tissue Sampled` = tissue,
    `Sampling Time Point` = time_point,
    `Pathogen Strain` = pathogen_strain
  ) %>%
  write_csv(., "outputs/fig1_supplementary_table.csv")

# Alternative version of Figure 1 (horizontal layout)

d %>%
  left_join(., segment.limits, by = "bottom_label") %>%
  ggplot(aes(
    x = prop_de,
    y = reorder(bottom_label, desc(bottom_label)), 
    color = susceptibility)) +
  geom_segment(
    aes(y = bottom_label, yend = bottom_label, 
        x = min_measure, xend = max_measure),
    color = "darkgrey", size = 2
  ) +
  geom_point(size = 8) +
  scale_color_manual(values = alpha(c("dodgerblue", "firebrick2"), 0.9)) +
  xlab("Proportion of genes/contigs differentially expressed") +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.title.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    strip.text.y = element_text(angle = 0, size = 18),
    axis.text = element_text(size = 20),
    axis.title.x = element_text(size = 24)
  ) +
  facet_wrap(~study_mod_simple, scales = "free_y", ncol = 1, strip.position = "right") +
  theme(strip.text.y = element_markdown())

# ggsave("outputs/fig1alt.jpeg", width = 12, height = 14)

#==============================================================================


# Visualize modeling results

# Load m.all model

m.all <- readRDS("data/saved_models/m.all.RDS")
m.all.stanfit <- m.all@stanfit
m.all.out <- rstan::extract(m.all@stanfit)

# Generate parameter traceplot

traceplot <- rstan::traceplot(m.all.stanfit, ncol = 4)

levels(traceplot$data$parameter) <- c(
  "Intercept", "Host Species Type",
  "Bracamonte et al. 2019", "Davy et al. 2017", "Ellison et al. 2015",
  "Eskew et al. 2018", "Fuess et al. 2017", "Poorten and Rosenblum 2016",
  "Sutherland et al. 2014",
  "σ (Study)"
)

traceplot + 
  scale_color_viridis(discrete = TRUE, option = "plasma", end = 0.85) +
  ggtitle("Parameter trace plots for hierarchical Bayesian model of differential gene expression (all data)") +
  theme(
    text = element_text(size = 14, color = "black"),
    plot.title = element_text(size = 18),
    strip.text.x = element_text(size = 14, face = "bold"),
    legend.position = "none"
  )

ggsave("outputs/figS2.jpeg", width = 15, height = 10)

# Generate parameter dotchart

m.all.summary <- as.data.frame(m.all.out) %>%
  select(-11) %>%
  pivot_longer(cols = 1:10, names_to = "parameter") %>%
  group_by(parameter) %>%
  summarize(
    mean = mean(value),
    median = median(value),
    lower = HPDI(value, prob = 0.95)[1],
    upper = HPDI(value, prob = 0.95)[2]
  )

m.all.summary
precis(m.all, prob = 0.95, depth = 2)

m.all.summary %>%
  rename(
    term = parameter,
    estimate = mean, 
    conf.low = lower, 
    conf.high = upper
  ) %>%
  relabel_predictors(
    a = "Intercept",
    bS = "Host Species Type",
    sigma_study = "σ (Study)",
    a_study.1 = "Bracamonte et al. 2019", 
    a_study.2 = "Davy et al. 2017", 
    a_study.3 = "Ellison et al. 2015",
    a_study.4 = "Eskew et al. 2018", 
    a_study.5 = "Fuess et al. 2017", 
    a_study.6 = "Poorten and Rosenblum 2016",
    a_study.7 = "Sutherland et al. 2014"
  ) %>%
  dw_plot(
    dot_args = list(size = 3, color = "darkgreen"),
    whisker_args = list(color = "forestgreen")
  ) +
  xlab("Parameter Estimate") +
  xlim(-7.5, 7.5) +
  theme_minimal() + 
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    text = element_text(size = 14, color = "black"),
    legend.position = "none"
  ) +
  geom_vline(xintercept = 0, colour = "black", linetype = 2)

ggsave("outputs/figS3.jpeg", width = 6, height = 4)

# Load m.anno model

m.anno <- readRDS("data/saved_models/m.anno.RDS")
m.anno.stanfit <- m.anno@stanfit
m.anno.out <- rstan::extract(m.anno@stanfit)

# Generate parameter traceplot

traceplot <- rstan::traceplot(m.anno.stanfit, ncol = 4)

levels(traceplot$data$parameter) <- c(
  "Global Intercept", "Host Species Type",
  "Bracamonte et al. 2019", "Davy et al. 2017", "Ellison et al. 2015",
  "Eskew et al. 2018", "Fuess et al. 2017", "Poorten and Rosenblum 2016",
  "Sutherland et al. 2014",
  "σ (Study)"
)

traceplot + 
  scale_color_viridis(discrete = TRUE, option = "plasma", end = 0.85) +
  ggtitle("Parameter trace plots for hierarchical Bayesian model of differential gene expression (only annotated data)") +
  theme(
    text = element_text(size = 14, color = "black"),
    plot.title = element_text(size = 18),
    strip.text.x = element_text(size = 14, face = "bold"),
    legend.position = "none"
  )

ggsave("outputs/figS4.jpeg", width = 15, height = 10)

# Generate parameter dotchart

m.anno.summary <- as.data.frame(m.anno.out) %>%
  select(-11) %>%
  pivot_longer(cols = 1:10, names_to = "parameter") %>%
  group_by(parameter) %>%
  summarize(
    mean = mean(value),
    median = median(value),
    lower = HPDI(value, prob = 0.95)[1],
    upper = HPDI(value, prob = 0.95)[2]
  )

m.anno.summary
precis(m.anno, prob = 0.95, depth = 2)

m.anno.summary %>%
  rename(
    term = parameter,
    estimate = mean, 
    conf.low = lower, 
    conf.high = upper
  ) %>%
  relabel_predictors(
    a = "Intercept",
    bS = "Host Species Type",
    sigma_study = "σ (Study)",
    a_study.1 = "Bracamonte et al. 2019", 
    a_study.2 = "Davy et al. 2017", 
    a_study.3 = "Ellison et al. 2015",
    a_study.4 = "Eskew et al. 2018", 
    a_study.5 = "Fuess et al. 2017", 
    a_study.6 = "Poorten and Rosenblum 2016",
    a_study.7 = "Sutherland et al. 2014"
  ) %>%
  dw_plot(
    dot_args = list(size = 3, color = "darkgreen"),
    whisker_args = list(color = "forestgreen")
  )+
  xlab("Parameter Estimate") +
  xlim(-7.5, 7.5) +
  theme_minimal() + 
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    text = element_text(size = 14, color = "black"),
    legend.position = "none"
  ) +
  geom_vline(xintercept = 0, colour = "black", linetype = 2)

ggsave("outputs/figS5.jpeg", width = 6, height = 4)

#==============================================================================


# Figure 1b code

# Generate implied probability of differential expression values for both
# species types using posterior samples from the full fit model

preds <- data.frame(
  species_type = c(
    rep("non-susceptible", length(m.all.out$a)), 
    rep("susceptible", length(m.all.out$a))
  ),
  probs = c(logistic(m.all.out$a), logistic(m.all.out$a + m.all.out$bS))
)

# Summarize these implied probability values

preds.summary <- preds %>%
  group_by(species_type) %>%
  summarize(
    median = median(probs),
    mean = mean(probs),
    lower = HPDI(probs, prob = 0.95)[1],
    upper = HPDI(probs, prob = 0.95)[2]
  )

# Plot and save Figure 1b

# Sina plot version

fig1b <- preds %>%
  ggplot(aes(x = species_type, y = probs, group = species_type)) +
  ylab("Implied probability of differential gene\nexpression given pathogen exposure") +
  xlab("Host species type") +
  ylim(0, 0.4) +
  geom_sina(aes(color = species_type), alpha = 0.05) +
  geom_pointrange(
    data = preds.summary,
    aes(y = median, ymin = lower, ymax = upper, group = species_type),
    size = 1
  ) +
  scale_color_manual(values = alpha(c("dodgerblue", "firebrick2"), 0.7)) +
  theme_minimal() +
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position = "none",
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 18),
  )

fig1b

# Density plot version

fig1b <- preds %>%
  ggplot(aes(x = probs, group = species_type)) +
  xlab("\nImplied probability of differential gene expression given pathogen exposure") +
  ylab("Posterior density") +
  xlim(0, 0.15) +
  ylim(0, 125) +
  geom_density(aes(fill = species_type, color = species_type), size = 2) +
  scale_color_manual(
    values = alpha(c("dodgerblue", "firebrick2"), 0.9),
    name = "Host species type"
  ) +
  scale_fill_manual(
    values = alpha(c("dodgerblue", "firebrick2"), 0.2),
    name = "Host species type"
  ) +
  theme_minimal() +
  theme(
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text = element_text(size = 32),
    axis.title = element_text(size = 35),
    legend.position = "none",
  ) +
  annotate(
    geom = "text",
    label = "non-susceptible host species",
    x = 0.04, y = 80,
    size = 14,
    color = alpha("dodgerblue", 0.9)
  ) +
  annotate(
    geom = "text",
    label = "susceptible host species",
    x = 0.05, y = 27,
    size = 14,
    color = alpha("firebrick2", 0.9)
  )

fig1b

ggsave("outputs/fig1b.pdf", width = 20, height = 10)

# Save Figure 1a and 1b together as Figure 1

plot_grid(
  fig1a, fig1b, nrow = 2, hjust = -1,
  labels = c("(a)", "(b)"), label_size = 40
)

ggsave("outputs/fig1.pdf", width = 30, height = 30)

#==============================================================================


# Generate Figure 2 (have to manually export from RStudio to get a jpeg)

grViz("data/figures/flowchart.dot")

# Save as a pdf automatically

grViz("data/figures/flowchart.dot") %>%
  DiagrammeRsvg::export_svg() %>%
  charToRaw() %>%
  rsvg::rsvg_pdf(file = "outputs/fig2.pdf", width = 1400, height = 900)

#==============================================================================


# Old code for generating heatmaps (not used)

d %>%
  ggplot(aes(x = susceptibility, y = forcats::fct_rev(group_de))) +
  geom_tile(aes(fill = log2_prop_de)) +
  geom_text(aes(label = prop_de, fontface = "bold")) +
  xlab("Disease Susceptibility of Host") +
  scale_x_discrete(position = "top") +
  scale_fill_gradient(
    low = "white", high = "firebrick3", na.value = "gainsboro",
    limits = c(-8, 7)
  ) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 22),
    axis.text.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_blank(),
    strip.text.x = element_text(size = 17),
    legend.position = "none"
  ) +
  facet_wrap(~study_mod, scales = "free_y", ncol = 3)

# ggsave("outputs/fig1.jpeg", width = 22, height = 12)

d %>%
  ggplot(aes(x = susceptibility, y = forcats::fct_rev(group_de))) +
  geom_tile(aes(fill = log2_prop_de)) +
  geom_text(aes(label = prop_de, fontface = "bold")) +
  xlab("Disease Susceptibility of Host") +
  scale_x_discrete(position = "top") +
  scale_fill_gradient(
    low = "white", high = "firebrick3", na.value = "gainsboro",
    limits = c(-8, 7)
  ) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 22),
    axis.text.x = element_text(size = 20, face = "bold"),
    axis.title.y = element_blank(),
    strip.text.x = element_text(size = 17),
    legend.position = "none",
    axis.text.y = element_blank()
  ) +
  facet_wrap(~study_mod, scales = "free_y", ncol = 3)

# ggsave("outputs/fig1alt.jpeg", width = 14, height = 12)