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



#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include "config.h"
#include <algorithm>
#include <sstream>
#include <ctime>
#include <iomanip>


using namespace std;
const string VERSION = "igv_batch_v1.0";

//A struct for our regions
struct REGION{
	//CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	genotype	match_chr	match_start	match_end	sample
	string chr, pos, ref, alt, genotype, sample;
};

bool compareBySample(const REGION &a, const REGION &b)
{
	return a.sample < b.sample;
}
bool ReadFile(vector<REGION>& our_regions)
{
	ifstream infile(regions_file.c_str());
	if(!infile.is_open())
		return false;
	REGION dummy_region;
	bool first_line=true;
	string line;
	while (getline(infile, line))
	{
		istringstream ss(line);
		string throwaway;
		//CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	genotype	match_chr	match_start	match_end	sample
		ss >> dummy_region.chr >> dummy_region.pos
			>> throwaway //ID
			>> dummy_region.ref >> dummy_region.alt
			>> throwaway //QUAL
			>> throwaway //FILTER
			>> throwaway //INFO
			>> throwaway //FORMAT
			>> dummy_region.genotype
			>> throwaway //match_chr
			>> throwaway //match_start
			>> throwaway //match_end
			>> dummy_region.sample;
			
		if(first_line)
			first_line=false;
		else
			our_regions.push_back(dummy_region);
		
	}
	infile.close();
	return true;
}

bool WriteBatch(vector<REGION>& our_regions)
{
	ofstream outfile(out_file.c_str());
	if(!outfile){
		 cerr << "Can not open the results file for writing!\nIt may be in use. Close it and try again\n";
		 return false;
	}
	sort(our_regions.begin(), our_regions.end(), compareBySample);
	outfile << "new" << endl; 
			<< "genome " << genome << endl
			<< "snapshotDirectory " << snapshot_dir << endl;
			
	string sample = "dummy_start";
	for (vector<REGION>::iterator me = our_regions.begin();me != our_regions.end();me++)
	{
		if(sample != me->sample)//outputing a new sample, load a new bam
		{
			sample = me->sample;
			outfile << "new" << endl << "load BAM_FOR_" << me->sample << endl;
		}
		outfile << "goto " << me->chr << ":" << me->pos << "-" << me->pos << endl;
		string filename = me->sample+"_"+me->chr+"_"+me->pos+"_"+me->ref+"_"+me->alt+"_";
		if(me->genotype == "1/1") filename+="hom.png";
		else if(me->genotype == "0/1") filename+="het.png";
		outfile << "snapshot " << filename << endl;
	}
	
	outfile.close();

	return true;
}

bool WriteLog(int num_regions,bool write_batch)
{
	string out = out_file + ".Log";
	ofstream outfile(out.c_str());
	if(!outfile) return false;
	outfile << "Program: " <<  VERSION << endl;
	outfile << "Number regions read: " << num_regions << endl;
	outfile << "Regions file used: " << regions_file << endl;
	time_t t = time(0);   // get time now
    struct tm * now = localtime( & t );
	if(write_batch)
		outfile << "Analysis successfully ran: ";
	else
		outfile << "Analysis failed (unable to write batch file): ";
	outfile << (now->tm_year + 1900) << '-' 
			<< (now->tm_mon + 1) << '-'
			<<  now->tm_mday
			<< endl;
}

bool DefineConfig()
{
	string configname = VERSION + ".config";
	ifstream infile( configname.c_str() );	//open file
	if(!infile.is_open())
		return false;
	
	string line;
	int line_num=1;
	while (getline(infile, line))
	{
		istringstream ss(line);
		if(line_num ==1)
			ss >> regions_file;
		if(line_num ==3)
			ss >> snapshot_dir;
		if(line_num ==5)
			ss >> genome;
		if(line_num ==7)
			ss >> out_file;
		line_num++;
	}
	infile.close();
	return true;
}



int main()
{
	vector<REGION> our_regions;	//main vector of regions 

	string filename;
	cout << "***********************\n***********************\n"
	     << VERSION << endl
		 << "\n***********************\n***********************\n";
	cout << "\nLists to compare should be tab-delimited files in the following format\n";
	cout << "//CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tgenotype\tmatch_chr\tmatch_start\tmatch_end\tsample\n";
	cout << "\nKeep headers in your input files! WARNING: For now remove commas from numbers. ie dont use 1,000,000, use 1000000 in input files\n";
	cout << "Press enter to continue\n";
	cin.get();
	
	if(!DefineConfig()){
		cerr << "Unable to read the config file!\n";
		return 1;
	}
	if(!ReadFile(our_regions)){
		cerr << "Looks like there was a problem reading " << regions_file;
		cerr << "Make sure you updated the config file \n";
		return 1;
	}
	
	cout << "Running now\n";
	int num_regions = our_regions.size();
	bool write_batch = false;
	if(WriteBatch(our_regions)){
		write_batch = true;
		cout << "All set. Check you're output file: " << out_file << endl;
	}
	if(!WriteLog(num_regions, write_batch))
		cerr << "Couldn't open log file for writing\n";
	
	return 0;
}
	
	