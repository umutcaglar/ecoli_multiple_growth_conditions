# Data and code for: The E. coli molecular phenotype under different growth conditions

*Mehmet U. Caglar, John R. Houser, Craig S. Barnhart, Daniel R. Boutz, Sean M. Carroll, Aurko Dasgupta, Walter F. Lenoir, Bartram L. Smith, Viswanadham Sridhara, Dariya K. Sydykova, Drew Vander Wood, Christopher J. Marx, Edward M. Marcotte, Jeffrey E. Barrick, Claus O. Wilke*


The repository consists primarily of five different groups of folders, labeled with `a_`, `b_`, `c_`, `d_`, and `e_`, representing five different stages of the analysis. For each group, there is a `code` folder that holds all the code, a `results` folder that holds all generated output datafiles, and a `figures` folder that holds all generated figures.

As an example, the folder `c_figures` includes the figures generated by the code in the folder `c_code_change_wrt_variables_RNA&Protein`.

In addition to these five groups of folders, we have data folders named:
 
- `RNA_reads_06_13_2016`: contains RNA abundances
- `Protein_reads_06_16_2016`: contains protein abundances
- `DoublingTimeRawData`: contains doubling times
- `Flux_reads_06_21_2016`: contains flux data

We also have a folder `text` that contains the manuscript and a folder `dictionary` that includes files that map "YP" names for proteins to "ECB" names for RNAs.

Protein and RNA raw files are txt files that includes counts and gene / protein names.

The analysis pipeline follows alphabetical order in the code folders:

a) `dataPreperation_RNA&Protein`: The code in this folder combines distinct count files into single unique ones for mRNA and proteins; then it generates the meta data file associated with them.

	`a_results` file contains:
	- "rnaMatrix_mRNA.csv" contains raw mRNA data combined
	- "resDf_protein_trT_set00_StcYtcNasAgrNgrMgh_SYAN_baseMgAllMg_baseNaAllNa_ExpAllPhase_noMatchFilter_p1Sf_noNorm.csv" contains mRNA data after DeSeq2 size factor normalization
	- "resDf_protein_trT_set00_StcYtcNasAgrNgrMgh_SYAN_baseMgAllMg_baseNaAllNa_ExpAllPhase_noMatchFilter_p1Sf_vst.csv" contains mRNA data after DeSeq2 size factor normalization and vst trasformation (similar to log transformation)
	- "rnaMatrix_meta.csv" is the meta mRNA file
	
	- "proteinMatrix.csv" contains raw protein data combined
	- "resDf_protein_trT_set00_StcYtcNasAgrNgrMgh_SYAN_baseMgAllMg_baseNaAllNa_ExpAllPhase_noMatchFilter_p1Sf_noNorm.csv" contains protein data after DeSeq2 size factor normalization
	- "resDf_protein_trT_set00_StcYtcNasAgrNgrMgh_SYAN_baseMgAllMg_baseNaAllNa_ExpAllPhase_noMatchFilter_p1Sf_vst.csv" contains protein data after DeSeq2 size factor normalization and vst trasformation (similar to log transformation)
	- "metaProtein.csv" is the meta protein file
	

b) `histogram_RNA&Protein`: This folder generates heat map figures for proteins and for RNA. Also calculates table 1.

c) `change_wrt_variables_RNA&Protein`: The code in this folder finds differentiall expressed genes using DeSeq2.

d) `Mf_Pathway_Analyze`:  The code in this folder generates figures and edited tables from results of DAVID runs.

e) `Flux_Analyze`: The code in this folder performs flux data analysis.  `FluxSaltStressAnalyze.R generates` concentration vs flux figures for RNA and proteins.

More details about the analysis pipeline can be found in the file  `instructions_for_pipeline.docx`.






