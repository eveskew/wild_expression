# This script contains summary analyses for cases where simple inspection 
# of relevant literature was insufficient to derive the metrics needed 
# for the meta-analysis (such as the number of genes/contigs/probes in an
# analysis or the number of such annotated entities)


# Load packages

library(tidyverse)

#==============================================================================


# Supplemental analyses for Bracamonte et al. 2019
# https://doi.org/10.1002/ece3.5728

# Calculate the number of DE contigs and number of annotated DE contigs for
# A. japonica
readxl::read_xlsx(
  "data/original_source_supp_info/Bracamonte_etal_2019/Bracamonte_et_al_2019_DEG_A-japonica.xlsx",
  sheet = 1,
  skip = 1) %>%
  slice(1:(n()-1)) %>%
  group_by(Time) %>%
  summarize(
    de_total_contigs_analyzed = n(),
    de_contigs_annotated = sum(!is.na(Name))
  )

# Calculate the number of DE contigs and number of annotated DE contigs for
# A. anguilla
readxl::read_xlsx(
  "data/original_source_supp_info/Bracamonte_etal_2019/Bracamonte_et_al_2019_DEG_A-anguilla.xlsx",
  sheet = 1,
  skip = 1) %>%
  slice(1:(n()-1)) %>%
  group_by(Time) %>%
  summarize(
    de_total_contigs_analyzed = n(),
    de_contigs_annotated = sum(!is.na(Name))
  )

#==============================================================================


# Supplemental analyses for Eskew et al. 2018
# https://doi.org/10.1098/rsos.170910

# Part 1 - Reshaping and saving data

# For this analysis, I didn't originally derive a complete list of annotated
# contigs in the reference transcriptome (unfortunately) or in fact record
# the differentially expressed contigs for each of the relevant experimental 
# conditions (very unfortunately). The first problem isn't a death knell, but 
# the second is frustrating because it appears the internals of edgeR DE 
# calculations may have changed in the interim. Thus, I had to re-run 
# "bioinformatics_analyses.R" from the paper repository 
# (https://github.com/eveskew/frog_chytrid_transcriptomics)
# locally in order to regenerate the "goframeData" and "data" objects. The
# number of DE contigs for each condition are very slightly different than 
# what was reported in the paper as a result. In "wild_expression_data.csv"
# I've opted to keep the overall DE results as reported in the paper, 
# supplemented with the contig annotation information newly calculated here.

# Derive a table listing the contigs used in the paper that had associated
# GO terms (i.e., those that were annotated)
annotated.contigs <- goframeData %>%
  select(Contig) %>%
  distinct() %>%
  mutate(Annotated_w_GO_term = rep(1, nrow(.)))

# How many contigs in the total set of 50,249 contigs were annotated?
nrow(annotated.contigs) # should be 22,940
# How many of these annotated contigs came from Bd?
sum(grepl("BDET", annotated.contigs$Contig))
# How many of these annotated contigs came from frogs?
sum(grepl("Lithobates|Lcla", annotated.contigs$Contig))

# I didn't originally append "data$genes" with the Lithobates clamitans AMP
# contigs, but this is easy enough to fix. Note that "data" is of length
# 41,646, representing the host-specific contig set that was analyzed in the
# paper
assertthat::assert_that(
  length(data$genes$Contig_Name) == length(rownames(data$counts)))
data$genes$Contig_Name <- rownames(data$counts)

# Join the contig GO term annotation information into the existing 
# contig "annotation" data. Add on columns indicating whether a given contig
# was differentially expressed in a given condition
table <- data$genes %>%
  left_join(., annotated.contigs, by = c("Contig_Name" = "Contig")) %>%
  mutate(
    Annotated_w_GO_term = ifelse(
      is.na(Annotated_w_GO_term), 
      0, 
      Annotated_w_GO_term
    ),
    Bull_CarVCon_3 = ifelse(Contig_Name %in% genes_Bull_CarVCon_3, 1, 0),
    Bull_SecVCon_3 = ifelse(Contig_Name %in% genes_Bull_SecVCon_3, 1, 0),
    Bull_SecVCon_7 = ifelse(Contig_Name %in% genes_Bull_SecVCon_7, 1, 0),
    Bull_SecVCon_10 = ifelse(Contig_Name %in% genes_Bull_SecVCon_10, 1, 0),
    Wood_CarVCon_3 = ifelse(Contig_Name %in% genes_Wood_CarVCon_3, 1, 0),
    Wood_CarVCon_7 = ifelse(Contig_Name %in% genes_Wood_CarVCon_7, 1, 0),
    Wood_CarVCon_10 = ifelse(Contig_Name %in% genes_Wood_CarVCon_10, 1, 0),
    Wood_SecVCon_3 = ifelse(Contig_Name %in% genes_Wood_SecVCon_3, 1, 0),
    Wood_SecVCon_7 = ifelse(Contig_Name %in% genes_Wood_SecVCon_7, 1, 0),
    Wood_SecVCon_10 = ifelse(Contig_Name %in% genes_Wood_SecVCon_10, 1, 0)
  )

# Save this table for later reference, if needed
write_csv(table, "data/original_source_supp_info/Eskew_etal_2018/Eskew_etal_2018_annotation_summary.csv")

# Part 2 - Output relevant calculations

# Read back in data
table <- read_csv("data/original_source_supp_info/Eskew_etal_2018/Eskew_etal_2018_annotation_summary.csv")

# How many host-specific contigs had GO term annotations?
sum(table$Annotated_w_GO_term == 1)

# Output summary calculations showing the number of DE contigs and number of
# annotated DE contigs for all relevant conditions
sum(table$Bull_CarVCon_3 == 1)
sum(table$Annotated_w_GO_term == 1 & table$Bull_CarVCon_3 == 1)
sum(table$Bull_SecVCon_3 == 1)
sum(table$Annotated_w_GO_term == 1 & table$Bull_SecVCon_3 == 1)
sum(table$Bull_SecVCon_7 == 1)
sum(table$Annotated_w_GO_term == 1 & table$Bull_SecVCon_7 == 1)
sum(table$Bull_SecVCon_10 == 1)
sum(table$Annotated_w_GO_term == 1 & table$Bull_SecVCon_10 == 1)

sum(table$Wood_CarVCon_3 == 1)
sum(table$Annotated_w_GO_term == 1 & table$Wood_CarVCon_3 == 1)
sum(table$Wood_CarVCon_7 == 1)
sum(table$Annotated_w_GO_term == 1 & table$Wood_CarVCon_7 == 1)
sum(table$Wood_CarVCon_10 == 1)
sum(table$Annotated_w_GO_term == 1 & table$Wood_CarVCon_10 == 1)
sum(table$Wood_SecVCon_3 == 1)
sum(table$Annotated_w_GO_term == 1 & table$Wood_SecVCon_3 == 1)
sum(table$Wood_SecVCon_7 == 1)
sum(table$Annotated_w_GO_term == 1 & table$Wood_SecVCon_7 == 1)
sum(table$Wood_SecVCon_10 == 1)
sum(table$Annotated_w_GO_term == 1 & table$Wood_SecVCon_10 == 1)

#==============================================================================


# Supplemental analysis for Poorten and Rosenblum 2016
# https://doi.org/10.1111/mec.13871

# This paper's Table S2 contains information regarding the number of
# probesets analyzed and the number of those that were annotated (with GO
# term information). I calculate relevant summary metrics, using the authors'
# original criteria for DE significance (< 0.1 for cane toads, < 0.05 for 
# boreal toads)

# cane toad - skin
row1 <- readxl::read_xlsx(
  "data/original_source_supp_info/Poorten_and_Rosenblum_2016/mec13871-sup-0002-tables2.xlsx",
  sheet = 1
) %>%
  summarize(
    n_total_probesets_analyzed = n(),
    de_total_probesets_analyzed = sum(p.value.adj.BH.Cane < 0.1),
    n_probesets_annotated = sum(Anno_Xtrop_Ensembl_GOterm != "NA"),
    de_probesets_annotated = 
      sum(p.value.adj.BH.Cane < 0.1 & Anno_Xtrop_Ensembl_GOterm != "NA")
  )

# cane toad - liver
row2 <- readxl::read_xlsx(
  "data/original_source_supp_info/Poorten_and_Rosenblum_2016/mec13871-sup-0002-tables2.xlsx",
  sheet = 2
) %>%
  summarize(
    n_total_probesets_analyzed = n(),
    de_total_probesets_analyzed = sum(p.value.adj.BH.Cane < 0.1),
    n_probesets_annotated = sum(Anno_Xtrop_Ensembl_GOterm != "NA"),
    de_probesets_annotated = 
      sum(p.value.adj.BH.Cane < 0.1 & Anno_Xtrop_Ensembl_GOterm != "NA")
  )

# cane toad - spleen
row3 <- readxl::read_xlsx(
  "data/original_source_supp_info/Poorten_and_Rosenblum_2016/mec13871-sup-0002-tables2.xlsx",
  sheet = 3
) %>%
  summarize(
    n_total_probesets_analyzed = n(),
    de_total_probesets_analyzed = sum(p.value.adj.BH.Cane < 0.1),
    n_probesets_annotated = sum(Anno_Xtrop_Ensembl_GOterm != "NA"),
    de_probesets_annotated = 
      sum(p.value.adj.BH.Cane < 0.1 & Anno_Xtrop_Ensembl_GOterm != "NA")
  )

# wood frog - skin
row4 <- readxl::read_xlsx(
  "data/original_source_supp_info/Poorten_and_Rosenblum_2016/mec13871-sup-0002-tables2.xlsx",
  sheet = 4
) %>%
  summarize(
    n_total_probesets_analyzed = n(),
    de_total_probesets_analyzed = sum(p.value.adj.BH.boreal < 0.05),
    n_probesets_annotated = sum(Anno_Xtrop_Ensembl_GOterm != "NA"),
    de_probesets_annotated = 
      sum(p.value.adj.BH.boreal < 0.05 & Anno_Xtrop_Ensembl_GOterm != "NA")
  )

# wood frog - liver
row5 <- readxl::read_xlsx(
  "data/original_source_supp_info/Poorten_and_Rosenblum_2016/mec13871-sup-0002-tables2.xlsx",
  sheet = 5
) %>%
  summarize(
    n_total_probesets_analyzed = n(),
    de_total_probesets_analyzed = sum(p.value.adj.BH.boreal < 0.05),
    n_probesets_annotated = sum(Anno_Xtrop_Ensembl_GOterm != "NA"),
    de_probesets_annotated = 
      sum(p.value.adj.BH.boreal < 0.05 & Anno_Xtrop_Ensembl_GOterm != "NA")
  )

# wood frog - spleen
row6 <- readxl::read_xlsx(
  "data/original_source_supp_info/Poorten_and_Rosenblum_2016/mec13871-sup-0002-tables2.xlsx",
  sheet = 6
) %>%
  summarize(
    n_total_probesets_analyzed = n(),
    de_total_probesets_analyzed = sum(p.value.adj.BH.boreal < 0.05),
    n_probesets_annotated = sum(Anno_Xtrop_Ensembl_GOterm != "NA"),
    de_probesets_annotated = 
      sum(p.value.adj.BH.boreal < 0.05 & Anno_Xtrop_Ensembl_GOterm != "NA")
  )

# Bind into a summary table
bind_rows(row1, row2, row3, row4, row5, row6)
