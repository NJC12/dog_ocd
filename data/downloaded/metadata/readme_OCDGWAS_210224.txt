10.05.2025 nora strom
md5checksum:
9ab61cb7412c42a0ba207856af6da0f8  daner_OCD_full_wo23andMe_190522.gz
ebf837d5556450243a13ce4d2ca90cf6  daner_OCD_full_incl23andMe_clump_W-3000_R-0p2_P-0p05_hrcALL.clumped_10000.gz
30066ad398c2bc5949647c82a4ec8ef8  daner_OCDmeta_wo23andMe_LOOUKBB_080425.gz


These summary statistics relate to the following publication:
Strom, N.I., Gerring, Z.F., Galimberti, M. et al. Genome-wide analyses identify 30 loci associated with obsessive–compulsive disorder. Nat Genet (2025). https://doi.org/10.1038/s41588-025-02189-z

Data without 23andMe:
daner_OCD_full_wo23andMe_190522.gz are the summary statistics of the OCD-3 GWAS, excluding 23andMe data. It contains 23493 European ancestry OCD cases and 1114613 European ancestry controls. Details of individual included studies and cohorts can be found in basic.OCD_full_wo23andMe_190522.num.xls
A list of independent GWAS results are in daner_OCD_full_wo23andMe_190522.gz.p4.clump.areator.sorted.1mhc.xls Manhattanplot can be found in manhattan.nog2.OCD_full_wo23andMe_190522.pdf and a QQ-plot in qq.OCD_full_wo23andMe_190522.pdf

Data without 23andMe and without UKBB:
daner_OCDmeta_wo23andMe_LOOUKBB_080425.gz are the summarystatistics excluding 23andMe as well as UKBB data. It contains 22717 European ancestry OCD cases and 988884 European ancestry controls. Details of individual included studies and cohorts can be found in basic.OCDmeta_wo23andMe_LOOUKBB_080425.num.xls
A list of independent GWAS results are in daner_OCDmeta_wo23andMe_LOOUKBB_080425.gz.p4.clump.areator.sorted.1mhc.xls Manhattanplot can be found in manhattan.nog2.OCDmeta_wo23andMe_LOOUKBB_080425.pdf and a QQ-plot in qq.OCDmeta_wo23andMe_LOOUKBB_080425.pdf

Data including 23andMe (10,000 SNPs)
daner_OCD_full_incl23andMe_clump_W-3000_R-0p2_P-0p05_hrcALL.clumped_10000.gz contains data of the full dataset, including 23andMe data (53660 European ancestry OCD cases and 2044417 European ancestry controls). The data was filtered to include 10,000 SNPs. The main result file was clumped with plink and then filtered to the top 10,000 SNPs. Clumping was done with the following command:
PLINK v1.9.0-b.7.7 64-bit (22 Oct 2024)
Options in effect:
  --bfile hrc_chrALL
  --clump daner_OCD_fullv6_240822_shortenednames_Beta_Neff.gz
  --clump-kb 3000
  --clump-p1 0.05
  --clump-r2 0.2
  --out clump_W-3000_R-0p2_P-0p05_hrcALL
  

The meta-analysis was done using metal within ricopili.
The files are in ricopili-daner-format (see here for details: https://docs.google.com/document/d/1TWIhr8-qpCXB13WCXcU1_HDio8lC_MeWoAg2jlggrtU/):

CHR: chromosome	
SNP: rs-numbers	
BP: base-pair position on genome-build GRCh37
A1: Effect allele	
A2: Allele 2	
FRQ_A_Ncases: A1 allele frequency of affected individuals, the max number of cases is appended to the header *
FRQ_U_Ncontrols: A1 allele frequency of unaffected individuals, the max number of controls is appended to the header *	
INFO: Imputation quality metric**	
OR: Effect size estimate with respect to A1	
SE: Standard error of the effect size estimate (referring to log(OR))	
P: p-value	
ngt: number of datasets where the SNP was genotyped directly (i.e. not imputed)***	
Direction: String of characters +,-,? Each character represents the direction of the estimated effect (+ is one direction, - the other, and ? means that this variant was missing in the original study).	
HetISqt: The I² statsitic for estimating heterogeneity among study effects. 
HetDf: Degrees of freedom for the heterogeneity test	
HetPVa: P-value of the test of heteogeneity among estimated effects in the studies	
Nca: The number of cases used in the test. 	
Nco: The number of controls in the overall meta-analysis. 	
Neff_half: The (half) effective sample size. For single cohorts this is defined as 4*Nca*Nco/(2*(Nca+Nco)). In meta analyses these numbers are summed over cohorts. 

* Note that not all individual studies provided separate allele frequencies for cases and controls.
** Also, not all of the individual datasets provided INFO scores. If they did not, the SNP infoscore was set to 1 for that study.
*** Not all individual datasets provided this information. If they did not, ngt was set to 0 for that study.


