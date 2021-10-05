#ifndef CONFIG_H
#define CONFIG_H

#include <string>
 
float percent_overlap_cutoff;	//Used to specify cutoff for overlap (0-1)
							
int overlap_type; //Used to specify the type of overlap calculation used. 1=union, 2=region, 3=cnv

int write_transcript_id; //Write transcript_id?
							//1=write transcript_id in output file, 0=do not write transcript_id
							
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