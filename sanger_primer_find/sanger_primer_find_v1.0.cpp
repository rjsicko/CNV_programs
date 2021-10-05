/*****************************************************
sanger_primer_find v1.0
Author: Robert Sicko - R.J.Sicko@gmail.com

****DESCRIPTION****
This program takes a list of variants and spits out the primer pair(s) to use to Sanger validate
****END DESCRIPTION****

****REQUIRED FORMAT****
	NOTE: No commas in numbers (ie. use 1000000 not 1,000,000)
	NOTE2: Keep headers in your files. The program ignores the first line (don't have spaces in header id's though).
	
	variant_list (tab-delimited txt):
	VarID	Chr	Start	Stop
	
	primer list(tab-delimited txt):
	Primer_Name		Start_hg19	End_hg19
	//note: file should be structured so pairs of primers are together (line 2/3, 4/5, 6/7 etc)
****END REQUIRED FORMAT****

****CHANGE LOG****
v1.0- Bugs fixed, programs appears to be functioning as expected

****END CHANGE LOG****

****TO DO****
- add a check if the variant falls in a primer binding region
- fix double output of last line - may be issue with reading in?
****END TO DO****
*********************************************************/




#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include "config.h"
#include <algorithm>
#include <sstream>
#include <iomanip>
#include <map>
#include <ctime>




using namespace std;

const string VERSION = "sanger_primer_find_v1.0";

//A struct for our CNVs (or regions)


struct PRIMER{
	string primer_id, chr;	//ints for chromosome, copy numbe and case or control (1, 2)
	long start, stop;	//start and end points for each cnv
	
};

struct VARIANT{
	string var_id, chr;
	long start, stop;	//start and end points for each cnv
	vector<pair<PRIMER,PRIMER> > primers_olap;
	
};



bool ReadFileVar(string filename,vector<VARIANT>& our_vars)
{
	ifstream infile(filename.c_str());
	if(!infile.is_open())
		return false;
	VARIANT dummy_var;
	bool first_line=true;
	string line;
	while (getline(infile, line))
	{
		istringstream ss(line);
		ss >> dummy_var.var_id >> dummy_var.chr
			>> dummy_var.start >> dummy_var.stop;
		
		if(first_line)
			first_line=false;
		else
			our_vars.push_back (dummy_var);
		
	}
	infile.close();
	return true;
}

bool ReadFilePrimers(string filename,vector<pair<PRIMER,PRIMER> >& primers)
{
	//cout << "im in readfileprimers " << endl;
	ifstream infile(filename.c_str());
	if(!infile.is_open())
		return false;
	PRIMER dummy_primer;
	int line_num=1;
	string line;
	pair<PRIMER,PRIMER> dummy_pair;
	while (getline(infile, line))
	{
		//cout << line_num << endl;
		//cout << "line_num % 2 = " << line_num % 2 << endl;	
		istringstream ss(line);
		//primer_id	chr	Start_hg19	End_hg19
		ss >> dummy_primer.primer_id >> dummy_primer.chr
		>> dummy_primer.start >> dummy_primer.stop;
		
		
		if(line_num==1)
		{
			line_num++;
			continue;
		}
		else if(line_num % 2 == 0)
		{//start of a pair
			dummy_pair.first = dummy_primer;
			//cout << dummy_primer.primer_id << endl;
		}
		else
		{//second member of pair
			dummy_pair.second = dummy_primer;
			primers.push_back (dummy_pair);
			//cout << dummy_pair.first.primer_id << " " << dummy_pair.second.primer_id << endl;
		}
		//cout << line_num << ":" << line_num % 2 << endl;
		line_num++;
	}
	if(line_num % 2 != 0) //we didn't end on an odd line. note: we check that's it's even because we ++'ed before end of loop
		cerr << "WARNING: your primer file didn't have an even number of primers. Unpaired primer left over\nCheck your primer file\n";
	infile.close();
	return true;
}

double CalcOverlap(double start_cnv, double start_region, double end_cnv, double end_region)
{
	return ((min(end_cnv,end_region) - max(start_cnv,start_region))
							/
					(end_cnv-start_cnv));
}

	

//function to compare the master list with a test or comparison list
//modifies master_list_cnvs 'num_cases_overlap_with' variable to keep track of the number of overlaps each entry in master list has with the comparison list 
void CheckRegions(vector<VARIANT>& our_vars, vector<pair<PRIMER,PRIMER> >& our_primers)
{
	for (vector<VARIANT>::iterator me = our_vars.begin();me != our_vars.end();me++)
	{
		for(vector<pair<PRIMER,PRIMER> >::iterator them = our_primers.begin();them != our_primers.end();them++)
		{
			if(me->chr==them->first.chr)
			{
				//cout << "their chromosomes were the same\n";
				double me_start = me->start, me_stop = me->stop;
				double them_start = min(them->first.start,them->second.start);
				double them_stop = max(them->first.stop, them->second.stop);
				double overlap;
				overlap = CalcOverlap(me_start,them_start,me_stop,them_stop);
				
				if(overlap > 0)
				{
					//cout << "overlap: " << overlap << endl;
					me->primers_olap.push_back(*them);
					cout << me->primers_olap.size() << endl;
				}
			}
		}
	}
}

bool WriteVars(vector<VARIANT>& our_vars)
{
	//cout << "I got to writeCNVs" << endl;
	ofstream outfile(out_file.c_str());
	if(!outfile){
		 cerr << "Can not open the results file for writing!\nIt may be in use. Close it and try again\n";
		 return false;
	}
	//Chr	Start	Stop
	outfile << "VarID\tChr\tStart\tStop\tprimers_to_use\n";
	for (vector<VARIANT>::iterator me = our_vars.begin();me != our_vars.end();me++)
	{
		outfile << me->var_id << "\t" << me->chr << "\t" 
		<< me->start << "\t" << me->stop << "\t";
		string primer_string;
		for(vector<pair<PRIMER,PRIMER> >::iterator primers = me->primers_olap.begin();primers != me->primers_olap.end();primers++)
		{
			if(primer_string!="")	//if our gene list is not empty
				primer_string = primer_string + ";";	//add ; at end
			primer_string = primer_string + primers->first.primer_id +"&" + primers->second.primer_id;
		}
		outfile << primer_string << endl;
	}
	outfile.close();	//close the file
	return true;
}

bool WriteLog(int num_vars, int num_primers)
{
	string out_trans = out_file + ".Log";
	ofstream outfile(out_trans.c_str());
	if(!outfile){
		 cerr << "Can not open the log file for writing!\nIt may be in use. Close it and try again\n";
		 return false;
	}
	outfile << "Program: " <<  VERSION << endl;
	outfile << "Number variants read: " << num_vars << endl;
	outfile << "Number primers read: " << num_primers << endl;
	outfile << "Variants file used: " << variant_file << endl;
	outfile << "Primer file used: " << primer_file << endl;
	time_t t = time(0);   // get time now
    struct tm * now = localtime( & t );
	outfile	<< "Analysis ran: " 
			<< (now->tm_year + 1900) << '-' 
			<< (now->tm_mon + 1) << '-'
			<<  now->tm_mday
			<< endl;
}
	

//The following function reads our configuration settings from a file and stores the data in variables defined in config.h
bool DefineConfig()
{
	string configname = VERSION + ".config";
	ifstream infile( configname.c_str() );	//open file
	if(!infile.is_open())	//open failed
		return false;
	
	string line;
	int line_num=1;
	while (getline(infile, line))	//lets read one line at a time then parse each line
	{
		istringstream ss(line);	//istringstream for parsing fields
		//string temp;
		
		//I know... BAD practice to hardcode line numbers here
		//but it's all I could come up with to make it work in a resonable amount of time
		if(line_num ==1)
			ss >> variant_file;
		if(line_num ==6)
			ss >> primer_file;
		if(line_num ==9)
			ss >> out_file;
		line_num++;
	}
	if (line_num != 13){
		cout << line_num << endl;
		cerr << "Looks like the config file is corrupt! Did you change what lines the data was on?\n";
		return false;
	}
	infile.close();
	return true;
}



int main()
{
	vector<VARIANT> our_vars;	//main vector of variants
	vector<pair<PRIMER,PRIMER> > our_primers;

	cout << "***********************\n***********************\n"
	     << VERSION << endl
		 << "\n***********************\n***********************\n";
	cout << "Press enter to continue\n";
	cin.get();
	
	if(!DefineConfig()){
		cerr << "Unable to read the config file!\n";
		return 1;
	}
	if(ReadFileVar(variant_file,our_vars))
	{
		//cout << "Thanks! Now specify second file to read (Including path)\n";
		//cin >> filename;
		if(!ReadFilePrimers(primer_file,our_primers))
		{
			//throw runtime_error("Cannot find config item.");
			
			cerr << "Looks like there was a problem reading " << primer_file;
			cerr << " Make sure you updated the compare_lists.config file with the two files you want to compare\n";
			cerr << "Also, the cnv files need to be in the same directory as the program\n";
			//cerr << "\nPerhaps you didn't specify the path?\n";
			return 1;
		}
	}
	else
	{
		cerr << "Looks like there was a problem reading " << variant_file;
		cerr << " Make sure you updated the compare_lists.config file with the two files you want to compare\n";
		cerr << "Also, the cnv files need to be in the same directory as the program\n";		
		//cerr << "\nPerhaps you didn't specify the path?\n";
		return 1;
	}
	
	cout << "Running now\n";
	CheckRegions(our_vars, our_primers);	
	if(WriteVars(our_vars))
		cout << "All set. Check you're output file: " << out_file << endl;
	int num_vars = our_vars.size(),
		num_primers = our_primers.size();
	WriteLog(num_vars, num_primers);			
	//ConsolePrintCnvs(master_list_cnvs);
	return 0;
}
	
	