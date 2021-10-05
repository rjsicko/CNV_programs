/*****************************************************
isca_overlap_v1.0
Author: Robert Sicko - R.J.Sicko@gmail.com

****DESCRIPTION****
a program to output our CNVs overlap with ISCA entries. Outputs our CNV and the following additional columns:
ISCA Entry Overlapped
ISCA Percent Overlap
ISCA Phenotype(s)
ISCA Type(s)

****END DESCRIPTION****

****REQUIRED FORMAT****
	NOTE: No commas in numbers (ie. use 1000000 not 1,000,000)
	NOTE2: Keep headers in your files. The program ignores the first line (don't have spaces in header id's though).
	
	CNV List (tab-delimited txt):
	OurID	Chr		Start	Stop	Case/Control	Copy_Number
	
	Gene List(tab-delimited txt):
	database_from(or_type)	gene_name	chr		start	stop	
	
****END REQUIRED FORMAT****

****CHANGE LOG****
v1.0- Bugs fixed, programs appears to be functioning as expected


****END CHANGE LOG****

****TO DO****
	
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
#include <list>
#include <ctime>

using namespace std;

struct Overlap{
	string isca_id, phenotype, olap_percentage, type;
};


//A struct for our CNVs (or regions)
struct CNV{
	//ourID	chr	start	stop	case/control	copy_num
	string our_id, chr;
	int copy_num, case_control;
	long start, stop;	//start and end points for each cnv
	vector<Overlap> isca_entries_overlapped;
};


struct IscaEntry{
	//database_from	name	phenotype	chr	chromStart	chromEnd	CN
	string isca_id, phenotype, chr;	//block identifier
	long start, stop;	//start and end points for each cnv
	int copy_num;

};

double CalcOverlapUnion(double start_cnv, double start_region, double end_cnv, double end_region)
{
	return ((min(end_cnv,end_region) - max(start_cnv,start_region))
							/
					(max(end_cnv,end_region) - min(start_cnv,start_region)));
}

double CalcOverlapRegion(double start_cnv, double start_region, double end_cnv, double end_region)
{
	return ((min(end_cnv,end_region) - max(start_cnv,start_region))
							/
					(end_region-start_region));
}

double CalcOverlapCNV(double start_cnv, double start_region, double end_cnv, double end_region)
{
	return ((min(end_cnv,end_region) - max(start_cnv,start_region))
							/
					(end_cnv-start_cnv));
}



bool ReadFileCnv(string filename,vector<CNV>& our_cnvs)
{
	ifstream infile(filename.c_str());
	if(!infile.is_open())
		return false;

	CNV dummy_cnv;
	bool first_line=true;
	string line;
	while (getline(infile, line))
	{
		istringstream ss(line);
		//ourID	chr	start	stop	case/control	copy_num
		ss >> dummy_cnv.our_id >> dummy_cnv.chr
			>> dummy_cnv.start >> dummy_cnv.stop
			>> dummy_cnv.case_control >> dummy_cnv.copy_num;

		//cout << dummy_cnv.our_id << "\t" << dummy_cnv.chr << "\t" << dummy_cnv.start << endl;
		if(first_line)
			first_line=false;
		else
			our_cnvs.push_back (dummy_cnv);
		
	}
		
		//cout << "I got to the while loop for file reading\n";
	infile.close();
	//ConsolePrintCnvs(our_cnvs);
	return true;
}

bool ReadFileRegion(string filename,vector<IscaEntry>& isca_regions)
{
	ifstream infile(filename.c_str());
	if(!infile.is_open())
		return false;
	IscaEntry dummy_isca;
	bool first_line=true;
	string line;
	while (getline(infile, line))
	{
		istringstream ss(line);
		//string filler;
		//ss >> filler;	//ignore database_from field
		string token;
		int token_num = 0;
		while(getline(ss, token, '\t'))
		{
			if(token_num == 1)
				dummy_isca.isca_id = token;
			else if(token_num == 2)
				dummy_isca.phenotype = token;
			else if(token_num == 3)
				dummy_isca.chr = token;
			else if(token_num == 4)
				dummy_isca.start = atoi( token.c_str() );
			else if(token_num == 5)
				dummy_isca.stop = atoi( token.c_str() );
			else if(token_num == 6)
				dummy_isca.copy_num = atoi( token.c_str() );
			token_num++;	
		}
		//database_from	name	phenotype	chr	chromStart	chromEnd	CN
		//ss >> dummy_isca.isca_id >> dummy_isca.phenotype >> dummy_isca.chr
			//>> dummy_isca.start >> dummy_isca.stop >> dummy_isca.copy_num;
		if(first_line)
			first_line=false;
		else
			isca_regions.push_back (dummy_isca);
		
	}
		//cout << "I got to the while loop for file reading\n";
	infile.close();
	return true;
}
	

//function to compare the master list with a test or comparison list
//modifies master_list_cnvs 'num_cases_overlap_with' variable to keep track of the number of overlaps each entry in master list has with the comparison list 
void CheckRegions(vector<CNV>& master_list_cnvs, vector<IscaEntry>& isca_regions)
{
	for (vector<CNV>::iterator me = master_list_cnvs.begin();me != master_list_cnvs.end();me++)
	{
		for(vector<IscaEntry>::iterator them = isca_regions.begin();them != isca_regions.end();them++)
		{
			//cout << them->isca_id << "\t" << them->phenotype << "\t" << them->chr << "\t" << them->start << "\t" << them->stop << endl;
			if(me->chr==them->chr)
			{
				double me_start = me->start, me_stop = me->stop, them_start = them->start, them_stop = them->stop;	
					//double overlap_cnv = CalcOverlapCNV(me_start,them_start,me_stop,them_stop);
				double overlap_union = CalcOverlapUnion(me_start,them_start,me_stop,them_stop);
					//double overlap_region = CalcOverlapRegion(me_start,them_start,me_stop,them_stop);
			
				if( overlap_union > percent_overlap_cutoff )
				{	//If there is union overlap > out cutoff then lets keep track of it
				
					//cout << "there was an overlap" << endl;
					Overlap my_overlap;
					
					my_overlap.isca_id = them->isca_id;
					my_overlap.phenotype = them->phenotype;
					if(them->copy_num == 0 || them->copy_num == 1)
						my_overlap.type = "Loss";
					else if(them->copy_num == 3 || them->copy_num == 4)
						my_overlap.type = "Gain";
					
					ostringstream strs;
					overlap_union = overlap_union * 100;
					strs << fixed << setprecision (2) << overlap_union;
					string olap = strs.str();	
					if(olap == "100.00") olap = "100";
					
					my_overlap.olap_percentage = olap;
					
					me->isca_entries_overlapped.push_back(my_overlap);
				}
			}
		}
	}
}



bool WriteCnvs(vector<CNV> master_list_cnvs)
{
	//cout << "I got to writeCNVs" << endl;
	ofstream outfile(out_file.c_str());
	if(!outfile){
		 cerr << "Can not open the results file for writing!\nIt may be in use. Close it and try again\n";
		 return false;
	}
	outfile << "OurID\tChr\tStart\tStop\tCase/Control\tCopy_Number\tISCA Entry Overlapped\tISCA Percent Overlap\tISCA Phenotype(s)\tISCA Type(s)\n";
	
	for (vector<CNV>::iterator me = master_list_cnvs.begin();me != master_list_cnvs.end();me++)
	{
		int region_count=0;
		outfile << me->our_id << "\t" << me->chr << "\t" 
				<< me->start << "\t" << me->stop << "\t" 
				<< me->case_control << "\t" << me->copy_num << "\t";
		
		bool first_entry = true;
		string entries, percentages, phenotypes, types;
		for(vector<Overlap>::iterator isca_entries_itr = me->isca_entries_overlapped.begin();isca_entries_itr != me->isca_entries_overlapped.end(); isca_entries_itr++)
		{
			if(!first_entry)
			{
				entries = entries + ";";
				percentages = percentages + ";";
				phenotypes = phenotypes + ";";
				types = types + ";";
			}
			else
				first_entry = false;
			entries = entries + isca_entries_itr->isca_id;
			percentages = percentages + isca_entries_itr->olap_percentage;
			phenotypes = phenotypes + isca_entries_itr->phenotype;
			types = types + isca_entries_itr->type;		
		}
		if(entries == "")
			outfile << "NA\tNA\tNA\tNA" << endl;
		else
			outfile << entries << "\t" << percentages << "\t" << phenotypes << "\t" << types << endl;

	}
	outfile.close();	//close the file
	return true;
}

bool WriteLog(int num_cnvs, int num_regions)
{
	string out_trans = out_file + ".Log";
	ofstream outfile(out_trans.c_str());
	if(!outfile){
		 cerr << "Can not open the log file for writing!\nIt may be in use. Close it and try again\n";
		 return false;
	}
	outfile << "Program: isca_overlap_v1.0\n";
	outfile << "Number CNVs read: " << num_cnvs << endl;
	outfile << "Number regions(blocks) read: " << num_regions << endl;
	outfile << "CNV file used: " << master_cnv_file << endl;
	outfile << "Regions (blocks) file used: " << comparison_cnv_file << endl;
	outfile << "Percent overlap cutoff: " << percent_overlap_cutoff << endl;
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
	ifstream infile("isca_overlap_v1.0.config");	//open file
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
			ss >> percent_overlap_cutoff;
		if(line_num ==5)
			ss >> master_cnv_file;
		if(line_num ==10)
			ss >> comparison_cnv_file;
		if(line_num ==15)
			ss >> out_file;
		line_num++;
	}
		
		//cout << "I got to the while loop for file reading\n";
	infile.close();
	//ConsolePrintCnvs(our_cnvs);
	return true;
}



int main()
{
	vector<CNV> master_list_cnvs;	//main vectors of CNVs for comparing
	vector<IscaEntry> isca_regions;

	string filename;
	cout << "***********************\n***********************\nisca_overlap_v1.0\n***********************\n***********************\n";
	cout << "\nCNV list file should be a tab-delimited file in the following format:\n";
	cout << "OurID\tChr\tStart\tStop\tCase/Control\tCopy_Number\n";
	cout << "\nISCA file should be a tab-delimited file in the following format:\n";
	cout << "database_from\tname\tphenotype\tchr\tchromStart\tchromEnd\tCN\n";
	cout << "\nKeep headers in your input files! WARNING: For now remove commas from numbers. ie dont use 1,000,000, use 1000000 in input files\n";
	cout << "I will try to deal with commas in numbers in later versions, but for now leave them out\n";
	cout << "Press enter to continue\n";
	cin.get();

	//cout << "Specify first file to read (Including path)\n";
	//cin >> filename;
	
	if(!DefineConfig())
	{
		cerr << "Unable to read the config file!\n";
		return 1;
	}
	if(ReadFileCnv(master_cnv_file,master_list_cnvs))
	{

		if(!ReadFileRegion(comparison_cnv_file,isca_regions))
		{
			cerr << "Looks like there was a problem reading " << comparison_cnv_file;
			cerr << " Make sure you updated the config file with the two files you want to compare\n";
			return 1;
		}
	}
	else{
		cerr << "Looks like there was a problem reading " << master_cnv_file;
		cerr << " Make sure you updated the config file with the two files you want to compare\n";
		return 1;
	}
	
	//ConsolePrint(regions);
	cout << "Running now\n";
	CheckRegions(master_list_cnvs, isca_regions);
	if(WriteCnvs(master_list_cnvs))
		cout << "All set. Check you're output file: " << out_file << endl;
		
	int num_cnvs = master_list_cnvs.size(),
		num_regions = isca_regions.size();
	WriteLog(num_cnvs, num_regions);		
	//ConsolePrintCnvs(master_list_cnvs);
	return 0;
}
	
	