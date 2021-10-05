Algorithm_Compare_v1.3
Author: Robert Sicko - R.J.Sicko@gmail.com

****DESCRIPTION****
This program compares output from seperate CNV calling algorithms
and writes a file with two appended columns. The two additional
columns are 'percents_overlapped' and 'total_overlap' for the percent
of our CNV overlapped by each other call and the total percent of
our CNV overlapped by all other calls. NOTE: these other calls
are in the same Subject_id!

****END DESCRIPTION****

****REQUIRED FORMAT****
	NOTE: No commas in numbers (ie. use 1000000 not 1,000,000).
	NOTE2: Keep headers in your files. The program ignores the first line (don't have spaces in header id's though).
	
	Both files should be tab-delimited in the following format(tab-delimited txt):
	OurID	Subject_id		Chr		Start	Stop	Case/Control	Copy_Number
****END REQUIRED FORMAT****


****CHANGE LOG****
v1.0- Bugs fixed, programs appears to be functioning as expected
v1.1- Output changed to include 0,1,2 or 3
		0=not in other algorithm
		1=exact match in other algorithm
		2=loose match in other algorithm (>=90% union overlap)
		3=exact match in other algorithm, but type differs (dup vs del)
		4=loose match in other algorithm, but type differs (dup vs del)
v1.2- Added runlog output	
v1.3- Overlap is now calculated as in the Pseudogene program.
	- Calculate total percent of our CNV overlapped by other calls (pileup other calls on our CNV)
	- Changed chr variable to string from int. Allows X, Y, XY to be read instead of converting to 23, 24, 25
****END CHANGE LOG****		
		

*********************************************************/