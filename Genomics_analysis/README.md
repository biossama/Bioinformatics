## Genome Analysis 

#### This script performs various analyses on the genome data of Fusarium graminearum. Below is a summary of the operations carried out:

   * Decompression:
        Decompresses the genome file while retaining the compressed version.
   * Identifying Sequences:
        Identifies sequences in the genome file along with their identifiers.
   * Genome Size and Nucleotide Counts:
        Calculates the size of the genome and the counts of A/T/G/C bases.
   * Potential Restriction Sites:
        Identifies potential restriction sites for the BamH1 enzyme.
   * Downloading Annotation Data:
        Downloads the structural annotation data of Fusarium graminearum.
   * Analyzing Annotation Data:
        Determines the types of features described and the number of annotated genes and mRNA.
   * Converting to BED Format:
        Converts annotation data to BED format, sorts it, and adds a prefix to chromosome names.
   * Calculating Gene Sizes:
        Calculates the size of each gene and the average gene size.
    Analyzing Gene Sizes:
        Counts the number of genes larger and smaller than 5kb.
    Extracting Promoter Coordinates:
        Extracts promoter coordinates for genes and creates a file containing these coordinates.
