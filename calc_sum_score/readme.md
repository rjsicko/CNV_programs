calc_sum_score V1.0
Author: Robert Sicko - R.J.Sicko@gmail.com

****DESCRIPTION****
This program compares a list of CNVs to a list of regions of interest
and outputs the original list with two additional columns appended.
The two columns appended are 'Overlapped_Genes' and 'Score_Overlapped'
if a CNV overlaps more than one region, the regions will be seperated
by a ; as will the percent overlaps for each
****END DESCRIPTION****

****REQUIRED FORMAT****
	NOTE: No commas in numbers (ie. use 1000000 not 1,000,000)
	NOTE2: Keep headers in your files. The program ignores the first line (don't have spaces in header id's though).
	
	CNV List (tab-delimited txt):
	OurID	Chr		Start	Stop	Case/Control	Copy_Number
	
	Region List(tab-delimited txt):
	Database_From	Unique_ID_1(gene)	Score	Chr		Start	Stop
****END REQUIRED FORMAT****

****CHANGE LOG****
v1.0- Modified from region of interest v1.4
	- Making this program to sum the haploinsufficiency and exac cnv scores for genes our CNVs overlap
		note - 161125_exac_final_cnv_gene_scores.txt created from:	 ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3.1/cnv/exac-final-cnv.gene.scores071316
			- 161125_haploinsuff_v3.txt created from: HI_Predictions_Version3.bed from https://decipher.sanger.ac.uk/about#downloads/data
****END CHANGE LOG****

****TO DO****
	- Add other checks for config file parameters outside expected range.
	- Comment out code!
	- Implement option to not write transcript.
	- Add check for output being too long for a single excel cell. In that case, use secondary column for output.
	- Fix case/control count so it doesn't just count the number of CNVs that overlap a gene, instead counts individuals
****END TO DO****
*********************************************************/