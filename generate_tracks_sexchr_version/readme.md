/*****************************************************
Generate_Tracks_v1.0_sexchr
Author: Robert Sicko - R.J.Sicko@gmail.com

****DESCRIPTION****
This program generates UCSC and DGV browser tracks from our data
****END DESCRIPTION****

****REQUIRED FORMAT****
	NOTE: No commas in numbers (ie. use 1000000 not 1,000,000).
	NOTE2: Keep headers in your files. The program ignores the first line (don't have spaces in header id's though).
	
	File should be tab-delimited in the following format(tab-delimited txt):
	OurID	Subject_id		Chr		Start	Stop	Case/Control	Copy_Number	Sex
****END REQUIRED FORMAT****


****CHANGE LOG****
v1.0- Bugs fixed, programs appears to be functioning as expected
v1.0_sexchr- Modified for sex chr cnvs
		   - Losses colored red, gains colored blue, LOH/normal colored green
****END CHANGE LOG****		
		

*********************************************************/
