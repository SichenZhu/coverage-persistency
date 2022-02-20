# coverage-persistency
This is a pipeline to analyze the abundance of low persistent index gene in microbial samples from different parts of human body. It is intended to verify, quantify, and visualize a general pattern (see sample [figure](SRR5723843-meta.pdf)) existing universtally in whole genome sequencing samples: the genes that were less persistent in the evolution process normally had a higher coverage (abundance) that is comparable to the highly persistent genes.

This pipeline applies bioinformatics computational tools such as [DIAMOND](https://github.com/bbuchfink/diamond), [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml), [SPAdes](https://github.com/ablab/spades), [Samtools](https://github.com/samtools/samtools), and [BedTools2](https://github.com/arq5x/bedtools2) to the whole genome sequencing datasets from different parts of the human body like blood, gut, nasopharyngeal swabs, middle ear, etc.

The Persistent Index (PI), which represents whether a large amount of species inherit this protein sequence from their common ancestor, is calculated using the [EggNOG database](http://eggnog5.embl.de/#/app/home). The Transcripts Per Million (TPM) is calculated using the information of coverage and protein alignment for each nucleotide.

This project may be helpful to understand the transferability of genes and mechanisms of some diseases.


Reference:

1. Buchfink, Benjamin, Chao Xie, and Daniel H. Huson. [Fast and sensitive protein alignment using DIAMOND](https://www.nature.com/articles/nmeth.3176). *Nature methods* 12.1 (2015): 59-60.

2. Langmead B, Salzberg S. [Fast gapped-read alignment with Bowtie 2](https://www.nature.com/articles/nmeth.1923). *Nature Methods*. 2012, 9:357-359.

3. Nurk, Sergey, et al. [metaSPAdes: a new versatile metagenomic assembler](https://pubmed.ncbi.nlm.nih.gov/28298430/). *Genome research* 27.5 (2017): 824-834.

4. Quinlan AR and Hall IM, 2010. [BEDTools: a flexible suite of utilities for comparing genomic features](https://academic.oup.com/bioinformatics/article/26/6/841/244688). *Bioinformatics*. 26, 6, pp. 841–842.

5. Danecek, Petr, et al. [Twelve years of SAMtools and BCFtools](https://academic.oup.com/gigascience/article/10/2/giab008/6137722?login=true). *Gigascience* 10.2 (2021): giab008.

6. Huerta-Cepas, Jaime, et al. [eggNOG 5.0: a hierarchical, functionally and phylogenetically annotated orthology resource based on 5090 organisms and 2502 viruses](https://academic.oup.com/nar/article/47/D1/D309/5173662?login=true). Nucleic acids research 47.D1 (2019): D309-D314.

7. The sample data for this study have been deposited in the European Nucleotide Archive (ENA) at EMBL-EBI under accession number: https://www.ebi.ac.uk/ena/browser/view/PRJNA384716
