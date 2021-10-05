#ifndef CONFIG_H
#define CONFIG_H

#include <string>
 
		
std::string master_cnv_file;	//this is the master cnv list
												//this list will be used in our output file with the
												//"number_overlap_cases" and "number_overlap_controls" columns added
												
												
std::string comparison_cnv_file;	//this is our test file
													//all entries in master_cnv_file will be tested
													//for overlap with all entries in this file		


		
std::string out_file; //file name to write output to
													//this file will be a tab-delimited file with
													//all entries from master_cnv_file and
													//"number_overlap_cases" and "number_overlap_controls" columns added

 
#endif