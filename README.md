# Non-canonical DNA in human and other ape telomere-to-telomere genomes
<https://www.biorxiv.org/content/10.1101/2024.09.02.610891v2>  
Linnéa Smeds, Kaivan Kamali, Iva Kejnovská, Eduard Kejnovský, Francesca Chiaromonte, Kateryna D. Makova  
Department of Biology, Penn State University, University Park, PA 16802  
Correspondence to Kateryna Makova (kdm16@psu.edu)  

***

## This repository includes the following files and directories:

### UNIX CODE 
- **general_commands.sh:** Includes a list of software, commands for creating bed files of non-B DNA motifs, generate statistics and find overlaps between the motif types. 

- **density_commands.sh:** Commands for calculating non-B DNA motif density and GC content in 100kb windows along the genome.

- **new_sequence_commands.sh:**
Analysis of new sequence in the T2T assemblies vs. old sequence (non T2T reference genomes), including code used for aligning new vs old. 

- **repeat_commands.sh:** Commands for enrichment analysis in repeats (including paths to publicliy available repeat annotations). Also contains detailed analysis of the satellites Walusat, LSAU and SST1. 

- **functional_commands.sh:** Commands for analysis of functional regions.

- **methylation_commands.sh:** Commands for methylation analysis of certain repeats and satellites.

- **centromere_commands.sh:** Analysis of centromeres, including centromeric satellites, CENP-B and SF analysis. Also includes paths to the publicly available annotations of these elements.

- **circos_commands.sh:** Scripts for generating input files for circos plots.

### DIRECTORIES
- **circos/:** Examples of circos scripts for human (including all configuration files and non-B DNA motifs data files).

- **python/:** In-house python scripts used in the study. 

- **R/:** R scripts to generate figures and statistics.

- **helpfiles/:** Textfiles with lists of species, files, repeats, and other information that are necessary to run most of the code in the different command files.
