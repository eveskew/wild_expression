This folder contains read and alignment statistics for each sample and the lists of differentially expressed genes (DEG) and overrepresented Gene Ontology (GO) terms for the article:

    Bracamonte SE, Johnston PR, Monaghan MT, Knopf K. Gene expression response to a nematode parasite in novel and native eel hosts. Ecology and Evolution

For questions and requests please contact Seraina E. Bracamonte: sebracamonte@gmx.net

European eels and Japanese eels were experimentally infected with Anguillicola crassus larvae or sham-infected with PBS and sampled at 3 days post-infection or 23 days post-infection (n = 5 for each species, time point, and treatment).
mRNA was paired-end sequenced on an Illumina HiSeq2500 or HiSeq4000 and reads were de novo assembled with Trinity. Raw reads and assemblies used to create the the count data are available from the BioProjects PRJNA419718 (European eel 3 dpi), PRJNA546508 (European eel 23 dpi & assembly), and PRJNA546510 (Japanese eel 3 & 23 dpi & assembly).
Differentially expressed genes between infected and control eels were identified separately for species and time points using DESeq2 v.1.14.0 on read counts obtained with RSEM v.1.3.0 and a mean coverage cut-off >= 10.
GO term enrichment analysis was done with GOstats v2.48.0 using custom background annotations obtained by blasting against UniProt/Swiss-Prot and searching the Pfam database integrated in Trinotate v3.2.0.

Files:

Bracamonte_et_al_2019_sample_sequencing_stats.xlsx: Number of raw reads, aligned reads, and alignment rate for each sample
    Species = A. anguilla (European eel) or A. japonica (Japanese eel)
    Sample = Sample ID
    Time = sampled at 3 or 23 day post-infection
    Treatment = control or infected
    # raw reads = number of raw sequencing reads of the specific sample
    # 1x aligned = number or reads that aligned concordantly exactly once to the reference transcriptome
    # >1x aligned = number or reads that aligned concordantly more than once to the reference transcriptome
    rate (%) = percent of reads that were aligned to the reference transcriptome

Bracamonte_et_al_2019_DEG_A-anguilla.xlsx: Differentially expressed genes of the European eel (Anguilla anguilla) with UniProt, RefSeq, or nr annotation at 3 and 23 days post-infection
    Gene = Gene identifier assigned by Trinity de novo assembly
    Base mean = mean expression across all samples
    Log2 FC = Log2 fold expression change between infected and control samples
    Wald stat = Wald statistic contrasting infected vs control
    Adj. p-value = Benjamini-Hochberg-adjusted p-value
    Name = annotation retrieved from UniProt/Swiss-Prot, RefSeq, Pfam, or NCBI nr databases
    Dir. = up-regulation (↑) or down-regulation (↓) in infected samples
    Time = 3 or 23 dpi

Bracamonte_et_al_2019_DEG_A-japonica.xlsx: Differentially expressed genes of the Japanese eel (Anguilla japonica) with UniProt, RefSeq, or nr annotation at 3 and 23 days post-infection
    Gene = Gene identifier assigned by Trinity de novo assembly
    Base mean = mean expression across all samples
    Log2 FC = Log2 fold expression change between infected and control samples
    Wald stat = Wald statistic contrasting infected vs control
    Adj. p-value = Benjamini-Hochberg-adjusted p-value
    Name = annotation retrieved from UniProt/Swiss-Prot, RefSeq, Pfam, or NCBI nr databases
    Dir. = up-regulation (↑) or down-regulation (↓) in infected samples
    Time = 3 or 23 dpi

Bracamonte_et_al_2019_GO_A-anguilla.xlsx: Overrepresented Gene Ontology (GO) terms of the European eel (Anguilla anguilla) at 3 and 23 days post-infection
    GO = Gene Ontology ID
    Term = Name of the GO term
    Expected = expected numbers of differentially expressed genes of the GO term
    Count = observed numbers of differentially expressed genes of the GO term
    Size = number of genes of the GO term in the reference transcriptome
    p-value
    Dir. = overrepresentation in up-regulated (↑) or down-regulated (↓) genes
    Time = 3 or 23 dpi

Bracamonte_et_al_2019_GO_A-japonica.xlsx: Overrepresented Gene Ontology (GO) terms of the Japanese eel (Anguilla japonica) at 3 and 23 days post-infection
    GO = Gene Ontology ID
    Term = Name of the GO term
    Expected = expected numbers of differentially expressed genes of the GO term
    Count = observed numbers of differentially expressed genes of the GO term
    Size = number of genes of the GO term in the reference transcriptome
    p-value
    Dir. = overrepresentation in up-regulated (↑) or down-regulated (↓) genes
    Time = 3 or 23 dpi