#ifndef CONFIG_H
#define CONFIG_H

#include <string>
 
std::string regions_file;	//this is the regions list
std::string out_file; //file name to write igv batch file to
std::string snapshot_dir;//this is the snapshot directory to use in the batch file
std::string genome;//this is the genome to use in the batch file
#endif