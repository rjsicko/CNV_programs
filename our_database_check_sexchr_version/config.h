#ifndef CONFIG_H
#define CONFIG_H

#include <string>
#include <vector>
			
std::string master_cnv_file;	//this is the master cnv list

std::string cumulative_database_file; //cumulative db file						
												
std::vector<std::string> similar_phenotypes;	
			
std::string out_file; //file name to write output to
						//this file will be a tab-delimited file with
						//all entries from master_cnv_file

 
#endif