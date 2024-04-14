# Metagenomic Data Analysis
## Context

Microorganisms represent the largest reservoir of genetic diversity on Earth, with estimates ranging from 4 to 6 × 10^30 cells. They are divided into two distinct domains: bacteria (or Eubacteria) and Archaea. Their diverse metabolic functions attest to their involvement in multiple processes such as organic matter mineralization, wastewater treatment, and methane emission. However, the fact that over 99% of environmental prokaryotes cannot be cultured in the laboratory significantly limits our knowledge of these communities. It is crucial to better characterize the composition of these communities and their role in the ecological balance of our planet. In this context, the "metagenomic" approach constitutes an important new methodology for deepening our understanding of prokaryotic diversity and its functional significance within continental aquatic ecosystems.

## Data

Sampling was conducted in the euphotic zone of Lake Bourget at a depth of 1.5 meters. These water samples were pre-filtered on different membranes to remove eukaryotic organisms (zooplankton, picoeukaryotes, etc.) and then concentrated by ultrafiltration. DNA was then extracted and quantified by fluorescence. A DNA library consisting of 7772 clones was created from 30 to 40 kb fragments. For each of these clones, terminal fragments were sequenced, resulting in 15,744 fragments. Here, only 500 clones, hence 1000 terminal fragments, will be studied. These sequences are available in the file: fragments1000.fasta.

## Comparison of Sequences to Annotated Genes

To study the phylogenetic diversity and metabolic composition of this ecosystem (euphotic zone of Lake Bourget), we need to compare the sequenced fragments to the known and annotated genes in the KEGG database. After downloading the file containing all genes from fully sequenced genomes, the file was formatted as follows:
```bash
formatdb -i genes.pep -n KEGG
```
The BLAST software allows us to compare our genomic DNA fragments to a protein database. This comparison between our fragments and the proteins of fully sequenced genomes in KEGG was done using BLASTX:
```bash
blastall -p blastx -d KEGG -i fragments1000.fasta -o resultat_Blast1000.tab -e 10.0 -m 8 -b 5 -a 4
```
The similarities found are available in the file resultats_Blast1000.tab. This file contains one line per significant alignment (up to 5 KEGG proteins per fragment).

## Known Information about KEGG Genes

  * taxonomy: Contains information about organisms present in KEGG.
  * annotation_genes: Contains information about KEGG proteins. Each line includes a protein name, its function, the orthologous group to which the gene encoding this protein belongs (known as KO), and finally, the list of metabolic pathways in which this protein is involved.
  * pathways: Contains, for each metabolic pathway, the pathway number (ko...), as well as its description in three levels.

## Analysis of Results

The objective is to identify the taxonomy and function of each sequenced DNA fragment to better understand the diversity and metabolic potential present in this ecosystem.

## Implementation

This project offers two different approaches for resolving the analysis:

   * Python: Utilizing Python scripts for data processing and analysis.
   * Bash: Using Bash scripts for data manipulation and computation.

