

# Load packages

library(tidyverse)
library(rethinking)

#==============================================================================


# Import quantitative summary data

d <- read_csv("data/quantitative_summary/wild_expression_data.csv") %>%
  mutate(
    total_tested = as.integer(total_tested),
    de_total_tested = as.integer(de_total_tested),
    total_annotated = as.integer(total_annotated),
    de_total_annotated = as.integer(de_total_annotated),
    prop_de = round(de_total_tested/total_tested, 2),
    prop_de_unrounded = de_total_tested/total_tested,
    study = as.factor(study),
    second_line = ifelse(
      is.na(pathogen_strain), 
      tissue, 
      paste(pathogen_strain, tissue, sep = " - ")
    ),
    second_line = paste(second_line, time_point, sep = " - "),
    group_de = paste(comparison, second_line, sep = "\n")
    %>%
      as.factor(),
    susceptibility_binary = ifelse(susceptibility == "susceptible", 1, 0),
    # Binary variable indicating whether or not the pathogen was Bd
    bd_binary = ifelse(pathogen == "B. dendrobatidis", 1, 0),
    # Binary variable indicating whether or not the assay used was a microarray
    microarray_binary = ifelse(assay == "microarray", 1, 0)
  )

# Subset down to create a table of disease-resistant vs. disease-susceptible
# species comparisons

comps <- d %>%
  select(
    study, susceptibility, group_de, prop_de_unrounded,
    bd_binary, microarray_binary
  ) %>%
  pivot_wider(names_from = susceptibility, values_from = prop_de_unrounded) %>%
  mutate(
    resistant_less = `non-susceptible` < susceptible,
    resistant_less_or_equal = `non-susceptible` <= susceptible
  ) %>%
  filter(!is.na(`non-susceptible`))

# How many disease-resistant vs. disease-susceptible comparisons are 
# there total?

n_total <- nrow(comps)
n_total

# In how many of the comparisons do the disease-resistant species have less
# gene expression change than the disease-susceptible species?

n_less <- sum(comps$resistant_less)
n_less

# Conduct a Chi-squared test

chisq.test(c(n_less, n_total - n_less), p = c(0.5, 0.5))

# Generate summary tables stratifying by the "bd_binary" and 
# "microarray_binary" variables

comps %>%
  group_by(bd_binary) %>%
  summarize(
    n_less = sum(resistant_less),
    prop_n_less = n_less/n()*100,
    n_equal_or_more = n() - sum(resistant_less),
    prop_n_equal_or_more = (n() - sum(resistant_less))/n()*100
  )

comps %>%
  group_by(microarray_binary) %>%
  summarize(
    n_less = sum(resistant_less),
    prop_n_less = n_less/n()*100,
    n_equal_or_more = n() - sum(resistant_less),
    prop_n_equal_or_more = (n() - sum(resistant_less))/n()*100
  )

#==============================================================================


# Subset to create a dataset appropriate for statistical modeling (i.e., each
# row is a single observation)

d.all <- d %>%
  filter(!is.na(de_total_tested)) %>%
  select(-comparison, -group_de) %>%
  distinct() %>%
  mutate(
    time_point_numeric = str_remove(time_point, " days"),
    time_point_numeric = as.numeric(time_point_numeric)
  )

nrow(d.all)

#==============================================================================


# Fit and save a hierarchical Bayesian model with a main effect of species type
# and a varying effect of study

m.all <- map2stan(
  data = d.all,
  alist(
    de_total_tested ~ dbinom(size = total_tested, prob = p),
    logit(p) <-   
      a + bS*susceptibility_binary + a_study[study],
    a ~ dnorm(0, 5),
    bS ~ dnorm(0, 5),
    a_study[study] ~ dnorm(0, sigma_study),
    sigma_study ~ dexp(1)
  ),
  chains = 4,
  iter = 7500,
  warmup = 2500,
  rng_seed = 8,
  WAIC = FALSE,
  control = list(adapt_delta = 0.99)
)

# Save the fit model

saveRDS(m.all, "data/saved_models/m.all.RDS")

#==============================================================================


# Fit and save a hierarchical Bayesian model with a main effect of species type
# and a varying effect of study, but only using gene expression data from
# annotated genes/contigs

# First generate a data subset that only contains annotated genes/contigs

d.anno <- filter(d.all, !is.na(de_total_annotated)) %>%
  mutate(
    study = as.factor(as.character(study))
  )

nrow(d.anno)

# Visualize the relationship between proportion of differential expression
# in the full gene set versus the annotated gene set

d.anno %>%
  ggplot(
    aes(
      x = de_total_tested/total_tested, 
      y = de_total_annotated/total_annotated, 
      color = susceptibility)
  ) +
  geom_point() +
  theme_minimal()

# Fit the model with this data

m.anno <- map2stan(
  data = d.anno,
  alist(
    de_total_annotated ~ dbinom(size = total_annotated, prob = p),
    logit(p) <-   
      a + bS*susceptibility_binary + a_study[study],
    a ~ dnorm(0, 5),
    bS ~ dnorm(0, 5),
    a_study[study] ~ dnorm(0, sigma_study),
    sigma_study ~ dexp(1)
  ),
  chains = 4,
  iter = 7500,
  warmup = 2500,
  rng_seed = 8,
  WAIC = FALSE,
  control = list(adapt_delta = 0.99)
)

# Save the fit model

saveRDS(m.anno, "data/saved_models/m.anno.RDS")

#==============================================================================


# Fit and save a hierarchical Bayesian model with a main effect of species type
# and a varying effect of study, adding on the main effect of pathogen type 
# (Bd or not) which interacts with species type

m.pathogen.type <- map2stan(
  data = d.all,
  alist(
    de_total_tested ~ dbinom(size = total_tested, prob = p),
    logit(p) <-   
      a + 
      bS*susceptibility_binary + bB*bd_binary + 
      bI*susceptibility_binary*bd_binary + 
      a_study[study],
    a ~ dnorm(0, 5),
    bS ~ dnorm(0, 5),
    bB ~ dnorm(0, 5),
    bI ~ dnorm(0, 5),
    a_study[study] ~ dnorm(0, sigma_study),
    sigma_study ~ dexp(1)
  ),
  chains = 4,
  iter = 7500,
  warmup = 2500,
  rng_seed = 8,
  WAIC = FALSE,
  control = list(adapt_delta = 0.99)
)

# Save the fit model

saveRDS(m.pathogen.type, "data/saved_models/m.pathogen.type.RDS")

#==============================================================================


# Fit and save a hierarchical Bayesian model with a main effect of species type
# and a varying effect of study, adding on the main effect of assay type 
# (microarray or not) which interacts with species type

m.assay.type <- map2stan(
  data = d.all,
  alist(
    de_total_tested ~ dbinom(size = total_tested, prob = p),
    logit(p) <-   
      a + 
      bS*susceptibility_binary + bM*microarray_binary + 
      bI*susceptibility_binary*microarray_binary + 
      a_study[study],
    a ~ dnorm(0, 5),
    bS ~ dnorm(0, 5),
    bM ~ dnorm(0, 5),
    bI ~ dnorm(0, 5),
    a_study[study] ~ dnorm(0, sigma_study),
    sigma_study ~ dexp(1)
  ),
  chains = 4,
  iter = 7500,
  warmup = 2500,
  rng_seed = 8,
  WAIC = FALSE,
  control = list(adapt_delta = 0.99)
)

# Save the fit model

saveRDS(m.assay.type, "data/saved_models/m.assay.type.RDS")
