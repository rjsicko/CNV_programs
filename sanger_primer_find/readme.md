/*****************************************************
sanger_primer_find v1.0
Author: Robert Sicko - R.J.Sicko@gmail.com

****DESCRIPTION****
This program takes a list of variants and spits out the primer pair(s) to use to Sanger validate
****END DESCRIPTION****

****REQUIRED FORMAT****
	NOTE: No commas in numbers (ie. use 1000000 not 1,000,000)
	NOTE2: Keep headers in your files. The program ignores the first line (don't have spaces in header id's though).
	
	variant_list (tab-delimited txt):
	VarID	Chr	Start	Stop
	
	primer list(tab-delimited txt):
	Primer_Name		Start_hg19	End_hg19
	//note: file should be structured so pairs of primers are together (line 2/3, 4/5, 6/7 etc)
****END REQUIRED FORMAT****

****CHANGE LOG****
v1.0- Bugs fixed, programs appears to be functioning as expected

****END CHANGE LOG****

****TO DO****
- add a check if the variant falls in a primer binding region
- fix double output of last line - may be issue with reading in?
****END TO DO****
*********************************************************/
