/*****************************************************
our_database_check_v1.3_sexchr
Author: Robert Sicko - R.J.Sicko@gmail.com

****DESCRIPTION****
program to check overlap with our CNV database, specific for sex chromosomes. For each CNV input we output the CNV with these additional columns:
'Exact_Unique_Times_Seen'
the following two columns:
	'Total % Our CNV With Cumulative Database Overlap(same phenotype, same CNV type)'
	'Cumulative Database CNVs Overlapped With (same phenotype, same CNV type)'
- these two columns repeat with:
	'similar phenotype, same CNV type'
	'different phenotype, same CNV type'
	'same phenotype, different CNV type'
	'similar phenotype, different CNV type'
	'different phenotype, different CNV type'

****END DESCRIPTION****

****REQUIRED FORMAT****
	NOTE: No commas in numbers (ie. use 1000000 not 1,000,000)
	NOTE2: Keep headers in your files. The program ignores the first line (don't have spaces in header id's though).
	
	File should be tab-delimited in the following format:
	ourID	sample	phenotype	case/control	chr	start	stop	copy_num
****END REQUIRED FORMAT****

****CHANGE LOG****
v1.0- Bugs fixed, programs appears to be functioning as expected
v1.1- Added runlog output
	- Added options 5-8 to output, see above.
	- Added output for Number_100, Number_90, Number_50, Number_1BP
	  instead of a single Number_overlapped column.
v1.2- New output, more informative.
	- Similar overlap calculation as in pseudogene check and dgv2 check.
	- Total % of our CNV covered by 1 or more cumulative database calls is now the output.
	- Previous versions didn't take into account that more than one cumulative DB
	  CNV could overlap our CNV. Now, we sum all overlap.
	- Also, ignore self in cumulative DB.
	  That is, don't count 100% overlap just because we matched ourself in the database.
	- Finally, now have three fields for overlap: exact, similar, different
	  These types will take into account the type of CNV.
v1.3- source completely changed
	- source a merger of case_control_v1.2 and CheckCumulativeDatabaseOverlap function from GVAT-0.0.2 and new additions and cleanup
	- output to changed again. We now output 13 columns:
		- first column still 'Exact_Unique_Times_Seen'
		- remaining 12 columns are the following two columns:
			'Total % Our CNV With Cumulative Database Overlap(same phenotype, same CNV type)'
			'Cumulative Database CNVs Overlapped With (same phenotype, same CNV type)'
		- these two columns repeat with:
			'similar phenotype, same CNV type'
			'different phenotype, same CNV type'
			'same phenotype, different CNV type'
			'similar phenotype, different CNV type'
			'different phenotype, different CNV type'
v1.3_sexchr- modified source for sex chr cnvs

****END CHANGE LOG****

****TO DO****

****END TO DO****

*********************************************************/
