****DESCRIPTION****
Genomic Variant Annotation Tool (GVAT) - an attempt to roll all the CNV annotation tools into a single program. I believe this was mostly working, but please check before relying on the results.
****END DESCRIPTION****

****REQUIRED FORMAT****
	NOTE: No commas in numbers (ie. use 1000000 not 1,000,000)
	NOTE2: Keep headers in your files. The program ignores the first line (don't have spaces in header id's though).
	
	File should be tab-delimited in the following format:
	CNV_ID	sample	phenotype	case_control	chr	start	stop	cn
****END REQUIRED FORMAT****

****CHANGE LOG****
0.0.1 	- Working version
0.0.2	- Lost main.cpp, have to recreate it
		- Now reading from a config file, no longer hardcoded
		- Fix size category to have  # and text instead of just #
		- Fixed 'Cumulative Database CNVs Overlapped With(X)' field, output NA if nothing overlapped
		- Added UCSC, IGV and DGV coordinate fields to the output file
		- Fixed CheckCumulativeDatabaseOverlap function to account for LOH in the file or the cumulativeDB
****END CHANGE LOG****

****TO DO****
	- Handle duplicates! Only count as single subject...
	- Ability to only run portions of the analysis
	- 
****END TO DO****