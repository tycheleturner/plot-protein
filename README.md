plot-protein
============

Plot Protein: Visualization of Mutations

Author: Tychele N. Turner, Ph.D.

Licenses: GNU General Public License version 3.0 (GPLv3), MIT License

Short Description: Protein Plotting Script to Visualize Amino Acid Changes

Programming Language: R

Current version: 3.0.0

Readme Date: 06/11/2016

Description: This script takes mutation information at the protein level and plots out the mutation above the schematic of the protein. It also plots the domains. It now has additional features for specifying the tick size of the x-axis, capability to show labels, and also the option to zoom to a particular region of the protein. 

NOTE: All files should be referring to the same isoform of the protein. This is imperative for drawing the plot correctly.

Required files:

*Mutation file: tab-delimited file containing 5 columns (ProteinId, GeneName, ProteinPositionOfMutation, ReferenceAminoAcid, AlternateAminoAcid) NO HEADER FOR NEEDED FOR THIS FILE

*Protein architecture file: tab-delimited file containing 3 columns (architecture_name, start_site, end_site). This file NEEDS the header and it is the same as what was previously written. This information can be downloaded from the HPRD (http://hprd.org/). Although the most recent files are quite old so looking in the web browser you can get much more up to date information.

*Post-translational modification file: This is a tab-delimited file with only one column and that is the site. This file NEEDS a header and is as previously written.


Basic Usage:
==================================================
```
Rscript plotProtein.R -m psen1_mutation_file.txt -a psen1_architecture_file.txt -p psen1_post_translation_file.txt -l 463
```

Advanced usage:
==================================================
```
Rscript plotProtein.R -m psen1_mutation_file.txt -a psen1_architecture_file.txt -p psen1_post_translation_file.txt -l 464 -n Disease -t 25 -s yes -z yes -b 50 -c 100
```

Running high throughput Plot Protein using snakemake:
==================================================
The set of instructions and input files for this section are different than that shown above. 

Need to have the mutation file formatted as follows:

```
    Column 1: GENE_HUGO_ID (Can use NA if unavailable)
    Column 2: PROTEIN_ID Required: Must match Protein Ids provided in the protein length file
    Column 3: STUDY_NAME (Can use NA if unavailable)
    Column 4: AMINO_ACID_POSITION Required: Amino Acid position of the variant
    Column 5: CHROM (Can use NA if unavailable)
    Column 6: POSITION (Can use NA if unavailable)
    Column 7: REF Allele (Can use NA if unavailable)
    Column 8: ALT Allele (Can use NA if unavailable)
    Column 9: ALLELE_FREQUENCY Optional column. (Can use NA if unavailable) 
```

 * Generate the list of proteins to plot: 

to plot all you could use the command below or you could just make your own list:
```
cut -f2 <mutation file> | sort | uniq > proteins_to_plot.txt
```

 * Get the github repository

```

git clone https://github.com/tycheleturner/plot-protein.git
cd plot-protein/high_throughput/
```

 * Fill out the config file. You'll need a post-translational modification file and a domain file. These can be downloaded from HPRD or you could make your own. Required information is shown below.

	* Post translational modification file has a column 4 with the protein id matching that of the mutation file and column 5 is the site.
	* Domain file has a column 3 with the protein id matching that of the mutation file, column 5 with the domain name, column 7 is the starting amino acid of the domain, and column 8 is the ending amino acid of the domain

 * Run the snakemake either locally OR by submitting to the cluster. Examples of both are shown below. ===

Local
```
snakemake
```

Submitting to cluster 
```
snakemake --cluster 'qsub {params.sge_opts}' -j 100 -w 30 -k
```

