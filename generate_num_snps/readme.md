/*****************************************************
generate_num_snps v1.0
Author: Robert Sicko - R.J.Sicko@gmail.com

****DESCRIPTION****
This program takes a list of CNVpartition CNVs (in tab format, with headers) 
and using the PFB file specified, generates the approximate number of SNPs for each CNV.
Only approximate because we use all markers in the PFB file, including those omitted from analysis.
Could be corrected by modifying PFB file to include only those used.

important note: this PFB file is specific to the omni2.5, you'll have to create your own if another array was used.
important note2: unzip the PFB file before use.
****END DESCRIPTION****

****REQUIRED FORMAT****
	NOTE: No commas in numbers (ie. use 1000000 not 1,000,000)
	NOTE2: Keep headers in your files. The program ignores the first line (don't have spaces in header id's though).
	
	CNV List (tab-delimited txt):
	chr  start  end  copy_num  sampleID  startSNP  endSNP  conf  num_SNP
	
	NOTE: this program functions correctly only if the PFB file is sorted by Chr & Position
	PFB file (tab-delimited txt):
	Name	Chr	Position	PFB
****END REQUIRED FORMAT****

****CHANGE LOG****
v1.0- Bugs fixed, programs appears to be functioning as expected

	
****END CHANGE LOG****

****TO DO****
	- Comment out code!
	- See if we can get exact number of SNPs
****END TO DO****
*********************************************************/
