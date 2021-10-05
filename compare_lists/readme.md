/*****************************************************
compare_lists_v1.2
Author: Robert Sicko - R.J.Sicko@gmail.com

****DESCRIPTION****
This program compares two CNV lists and outputs the first file with two additional
fields appended 'num_cases' & 'num_controls' that are the counts of the numbers of cases
and controls in file2 that (CNV) overlap with those in file1 by more than the threshold
specified in the config file.
****END DESCRIPTION****

****REQUIRED FORMAT****
	NOTE: No commas in numbers (ie. use 1000000 not 1,000,000)
	NOTE2: Keep headers in your files. The program ignores the first line (don't have spaces in header id's though).
	
	Both files should be tab-delimited in the following format:
	OurID	Chr		Start	Stop	Case/Control	Copy_Number
****END REQUIRED FORMAT****

****CHANGE LOG****
v1.0- Bugs fixed, programs appears to be functioning as expected
v1.1- Added runlog output
v1.2- Changed overlap type to CNV overlap instead of union
****END CHANGE LOG****

****TO DO****
	- add overlap type selection
****END TO DO****

*********************************************************/