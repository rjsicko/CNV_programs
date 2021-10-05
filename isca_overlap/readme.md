/*****************************************************
isca_overlap_v1.0
Author: Robert Sicko - R.J.Sicko@gmail.com

****DESCRIPTION****
a program to output our CNVs overlap with ISCA entries. Outputs our CNV and the following additional columns:
ISCA Entry Overlapped
ISCA Percent Overlap
ISCA Phenotype(s)
ISCA Type(s)

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


****END CHANGE LOG****

****TO DO****
	
****END TO DO****
*********************************************************/
