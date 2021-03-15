# Host gene expression in wildlife disease: making sense of species-level responses

This repository contains code, data, and figures that support:

Eskew, E.A., D. Fraser, M.J. Vonhof, M.L. Pinsky, and B. Maslo. In review. Host gene expression in wildlife disease: making sense of species-level responses.

--- 

### Repository Structure

- [`/data`](/data) contains all raw data files for the analysis
	- [`/figures`](/data/figures) contains `.dot` file needed to generate Figure 2
	- [`/original_source_supp_info`](/data/original_source_supp_info) contains, for ease of access, all supplementary material from the original papers analyzed in this work 
	- [`/quantitative_summary`](/data/quantitative_summary) contains the curated gene expression data that was harvested from original research papers
	- [`/saved_models`](/data/saved_models) contains fit Bayesian model objects
	- [`/web_of_science_query_results`](/data/web_of_science_query_results) contains documentation of all papers returned from the Web of Science Core Collection searches conducted during this work
- [`/outputs`](/outputs) contains all figures and tables output from the [`03_data_visualization.R`](/scripts/03_data_visualization.R) script
- [`/scripts`](/scripts) contains the primary analysis scripts, numbered in order of execution
