# Data and code for: The E. coli molecular phenotype under different growth conditions

*Mehmet U. Caglar, John R. Houser, Craig S. Barnhart, Daniel R. Boutz, Sean M. Carroll, Aurko Dasgupta, Walter F. Lenoir, Bartram L. Smith, Viswanadham Sridhara, Dariya K. Sydykova, Drew Vander Wood, Christopher J. Marx, Edward M. Marcotte, Jeffrey E. Barrick, Claus O. Wilke*

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
	
