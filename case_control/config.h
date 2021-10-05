#ifndef CONFIG_H
#define CONFIG_H

#include <string>
#include <vector>
			
std::string master_cnv_file;	//this is the master cnv list
												//this list will be used in our output file with the
												//"number_overlap_cases" and "number_overlap_controls" columns added
												
												
std::vector<std::string> similar_phenotypes;	
			
double olap_threshold;
		
std::string out_file; //file name to write output to
													//this file will be a tab-delimited file with
													//all entries from master_cnv_file and
													//"number_overlap_cases" and "number_overlap_controls" columns added

 
#endif