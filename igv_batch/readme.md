/*****************************************************
igv_batch_v1.0
Author: Robert Sicko - R.J.Sicko@gmail.com

****DESCRIPTION****
	takes a list of variants, creates a batch file that you can use in IGV.
	batch file will load each sample's BAM, go to variant location and save the snapshot to a png file named: "SAMPLE_CHR_LOCATION_REF_ALT_ZYGOSITY.png"
	note: currently need to find/replace "load BAM_FOR_" with "load BAM_location" in the output batch file.
****END DESCRIPTION****

****REQUIRED FORMAT****
	NOTE: No commas in numbers (ie. use 1000000 not 1,000,000)
	NOTE2: Keep headers in your files. The program ignores the first line (don't have spaces in header id's though).
	
	File should be tab-delimited in the following format:
	CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	genotype	match_chr	match_start	match_end	sample
****END REQUIRED FORMAT****

****CHANGE LOG****
v1.0- Bugs fixed, programs appears to be functioning as expected

****END CHANGE LOG****

****TO DO****
	add a SAMPLE	BAM map to the config file
****END TO DO****

*********************************************************/
