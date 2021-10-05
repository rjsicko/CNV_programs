/*****************************************************
case_control_v1.2
Author: Robert Sicko - R.J.Sicko@gmail.com

****DESCRIPTION****
	takes a list of CNVs and outputs the list + appends the following extra columns:
			% Overlap Control(same CNV type)
			% Overlap Control(diff CNV type)
			Num Cases Overlapped (same phenotype, same CNV type)[Threshold: olap_threshold]
			Num Cases Overlapped (similar phenotype, same CNV type)[Threshold: olap_threshold]
			Num Cases Overlapped (different phenotype, same CNV type)[Threshold: olap_threshold]
			Num Cases Overlapped (same phenotype, diff CNV type)[Threshold: olap_threshold]
			Num Cases Overlapped (similar phenotype, diff CNV type)[Threshold: olap_threshold]
			Num Cases Overlapped (different phenotype, diff CNV type)[Threshold: olap_threshold]
****END DESCRIPTION****

****REQUIRED FORMAT****
	NOTE: No commas in numbers (ie. use 1000000 not 1,000,000)
	NOTE2: Keep headers in your files. The program ignores the first line (don't have spaces in header id's though).
	
	File should be tab-delimited in the following format:
	ourID	sample	phenotype	case/control	chr	start	stop	copy_num
****END REQUIRED FORMAT****

****CHANGE LOG****
v1.0- Bugs fixed, programs appears to be functioning as expected
v1.1- Added ability to have multiple groups of similar phenotypes
v1.2- Added a threshold for overlap to count cases that overlap

****END CHANGE LOG****

****TO DO****

****END TO DO****

*********************************************************/