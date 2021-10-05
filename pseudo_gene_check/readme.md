/*****************************************************
pseudo_gene_check_v1.2
Author: Robert Sicko - R.J.Sicko@gmail.com

****DESCRIPTION****
This program compares a list of CNVs to a list of genes(pseudo in this case).
This program compares each CNV in the list to each gene and will output which
CNVs are overlapped(above a threshold) by a gene or multiple genes.

The four columns appended are 'gene_ids' 'gene_percents' and 'total_overlap'
if a CNV overlaps more than one gene, the gene_ids will be seperated
by a ; as will the percent overlaps for each
****END DESCRIPTION****

****REQUIRED FORMAT****
	NOTE: No commas in numbers (ie. use 1000000 not 1,000,000)
	NOTE2: Keep headers in your files. The program ignores the first line (don't have spaces in header id's though).
	
	CNV List (tab-delimited txt):
	OurID	Chr		Start	Stop	Case/Control	Copy_Number
	
	Gene List(tab-delimited txt):
	database_from(or_type)	gene_name	chr		start	stop	
	
****END REQUIRED FORMAT****

****CHANGE LOG****
v1.0- Bugs fixed, programs appears to be functioning as expected
v1.1- Fixed a bug in the 'total_overlap' calculation
	- v1.0 actually addressed the issue that 'total_overlap' could be greater than 100% if you simply sum
	  each OVERLAP in the vector of overlaps. However, the way v1.0 corrected it was also faulty. 
	  It calculated the correct overlap the majority of the time, but over estimated overlap when multiple 
	  seperate overlaps overlapped with each other.
	- To see how this was corrected, see function CoalesceOverlaps... briefly, after all overlaps
	  to a particular CNV were added to the vector of OVERLAPS, these overlaps were coalesced,
	  leaving non overlapping OVERLAPS that were stored in a list. This list was then passed
	  to CalcTotalOverlap for calculating total overlap with our CNV.
v1.2- Added runlog output

****END CHANGE LOG****

****TO DO****
	
****END TO DO****
*********************************************************/