/*****************************************************
compare_lists_exact_v1.1
Author: Robert Sicko - R.J.Sicko@gmail.com

****DESCRIPTION****
This program compares two CNV lists and outputs the first file with two additional
fields appended 'num_cases' & 'num_controls' that are the counts of the numbers of cases
and controls in file2 that (CNV) overlap with those in file1 by more than the threshold
specified in the config file. Only counts exact type overlaps.
****END DESCRIPTION****

****REQUIRED FORMAT****
	NOTE: No commas in numbers (ie. use 1000000 not 1,000,000)
	NOTE2: Keep headers in your files. The program ignores the first line (don't have spaces in header id's though).
	
	Both files should be tab-delimited in the following format:
	OurID	Chr		Start	Stop	Case/Control	Copy_Number
****END REQUIRED FORMAT****

****CHANGE LOG****
v1.0- Bugs fixed, programs appears to be functioning as expected
v1.1- Believe  1.0 was not comparing type of CNV correctly
	- Now we accurately compare CNV type
	- Also made the counts not include self
	

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
#include <ctime>


using namespace std;

//A struct for our CNVs (or regions)
struct CNV{
	string our_id;	//our identifier
	int chr, copy_num, case_control;	//ints for chromosome, copy numbe and case or control (1, 2)
	long start, stop;	//start and end points for each cnv
	int num_cases_overlap_with, num_controls_overlap_with; // count variables, +1 each time we overlap with a case or control
	
};

//Mostly for debugging
void ConsolePrintCnvs(vector<CNV> any_list_of_cnvs)
{
	//cout << "I got to ConsolePrintCnvs function\nAnd my cnvs list contains: ";
	//cout << any_list_of_cnvs.size();
	//cout << flush;
	cout << "OurID\tChr\tStart\tStop\tCase/Control\tCopy_Number\tOverlap_Count_Cases\tOverlap_Count_Controls\n";
	for (vector<CNV>::iterator me = any_list_of_cnvs.begin();me != any_list_of_cnvs.end();me++)
			cout << me->our_id << "\t" << me->chr << "\t" 
				<< me->start << "\t" << me->stop << "\t" 
				<< me->case_control << "\t" << me->copy_num 
				<< "\t" << me->num_cases_overlap_with
				<< "\t" << me->num_controls_overlap_with << endl;
		
}

bool ReadFile(string filename,vector<CNV>& our_cnvs)
{
	ifstream infile(filename.c_str());
	if(!infile.is_open())
		return false;
		//int marker=0;	//marker for where we are in the file //not needed?
	
									/*while(infile >> dummy_cnv.our_id){
									//start reading from our file
									//file should be formatted:
									//our_id		chr		start		stop		case_control		copy_number
									//char place_holder;
									int chr;
									infile >> place_holder >> dummy_cnv.chr >> place_holder 
											>> dummy_cnv.start >> place_holder >> dummy_cnv.stop
											>> place_holder >> dummy_cnv.case_control
											>> place_holder >> dummy_cnv.copy_num;
									infile >> chr;*/
	CNV dummy_cnv;
	bool first_line=true;
	string line;
	while (getline(infile, line))
	{
		istringstream ss(line);
		//string name;
		//int var1, var2, var3;
		ss >> dummy_cnv.our_id >> dummy_cnv.chr
			>> dummy_cnv.start >> dummy_cnv.stop
			>> dummy_cnv.case_control >> dummy_cnv.copy_num;
		//ss >> name >> var1 >> var2 >> var3;
		//first_line = false;
		dummy_cnv.num_cases_overlap_with = 0;
		dummy_cnv.num_controls_overlap_with = 0;
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


double CalcOverlapCNV(double start_cnv, double start_region, double end_cnv, double end_region)
{
	return ((min(end_cnv,end_region) - max(start_cnv,start_region))
							/
					(end_cnv-start_cnv));
}
	

//function to compare the master list with a test or comparison list
//modifies master_list_cnvs 'num_cases_overlap_with' variable to keep track of the number of overlaps each entry in master list has with the comparison list 
void CompareCnvs(vector<CNV>& master_list_cnvs, vector<CNV> comparison_list_cnvs)
{	
	//cout << type_specific << endl;
	for (vector<CNV>::iterator me = master_list_cnvs.begin();me != master_list_cnvs.end();me++)
		for(vector<CNV>::iterator them = comparison_list_cnvs.begin();them != comparison_list_cnvs.end();them++)
		{
			if(me->our_id==them->our_id)
				continue;	//skip to next CNV, don't want to count ourself as an overlap
			if(me->chr==them->chr)
			{
				double me_start = me->start, me_stop = me->stop, them_start = them->start, them_stop = them->stop;	
				double overlap = CalcOverlapCNV(me_start, them_start, me_stop, them_stop);
				/*if(me->our_id=="PUV_0561" && them->our_id=="PUV_0564")
					cout << overlap << endl;*/
				bool type_specific_pass = true;
				if(me->copy_num==them->copy_num)
					type_specific_pass = true;
				else
					type_specific_pass = false;
				if(type_specific_pass)
				{
					//cout << overlap << "&" << percent_overlap_cutoff << endl;
					if(overlap >= percent_overlap_cutoff)
					{
						if(them->case_control==1)
						{
							me->num_cases_overlap_with++;
							/*if(me->our_id=="PUV_0561")
							{
								cout << "case overlap:" << them->our_id << endl;	
							}*/
						}
						else
						{
							me->num_controls_overlap_with++;
							/*if(me->our_id=="PUV_0561")
							{
								cout << "control overlap:"  << them->our_id << endl;
							}*/
						}
					}
				}
			}
		}
}

bool WriteCnvs(vector<CNV> master_list_cnvs)
{
	ofstream outfile(out_file.c_str());
	if(!outfile){
		 cerr << "Can not open the results file for writing!\nIt may be in use. Close it and try again\n";
		 return false;
	}
	outfile << "OurID\tChr\tStart\tStop\tCase/Control\tCopy_Number\tOverlap_Count_Cases\tOverlap_Count_Controls\n";
	for (vector<CNV>::iterator me = master_list_cnvs.begin();me != master_list_cnvs.end();me++)
		outfile << me->our_id << "\t" << me->chr << "\t" 
				<< me->start << "\t" << me->stop << "\t" 
				<< me->case_control << "\t" << me->copy_num 
				<< "\t" << me->num_cases_overlap_with
				<< "\t" << me->num_controls_overlap_with << endl;
		
	outfile.close();
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
	outfile << "Program: compare_lists_exact_v1.1\n";
	outfile << "Number CNVs read: " << num_cnvs << endl;
	outfile << "Number comparison CNVs read: " << num_regions << endl;
	outfile << "CNV file used: " << master_cnv_file << endl;
	outfile << "Comparision file used: " << comparison_cnv_file << endl;
	outfile << "Percent overlap cutoff: " << percent_overlap_cutoff << endl;
	outfile << "Type specific on always in compare_lists_exact" << endl;
	time_t t = time(0);   // get time now
    struct tm * now = localtime( & t );
	outfile	<< "Analysis ran: " 
			<< (now->tm_year + 1900) << '-' 
			<< (now->tm_mon + 1) << '-'
			<<  now->tm_mday
			<< endl;
}

bool DefineConfig()
{
	ifstream infile("compare_lists_exact_v1.1.config");
	if(!infile.is_open())
		return false;
	
	string line;
	int line_num=1;
	while (getline(infile, line))
	{
		istringstream ss(line);
		//string temp;
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
	vector<CNV> master_list_cnvs, comparison_list_cnvs;	//main vectors of CNVs for comparing 

	string filename;
	cout << "***********************\n***********************\nCompare_Lists_Exact v1.1\n***********************\n***********************\n";
	cout << "\nLists to compare should be tab-delimited files in the following format\n";
	cout << "OurID\tChr\tStart\tStop\tCase/Control\tCopy_Number\n";
	cout << "\nKeep headers in your input files! WARNING: For now remove commas from numbers. ie dont use 1,000,000, use 1000000 in input files\n";
	cout << "I will try to deal with commas in numbers in later versions, but for now leave them out\n";
	cout << "Press enter to continue\n";
	cin.get();

	//cout << "Specify first file to read (Including path)\n";
	//cin >> filename;
	
	if(!DefineConfig())	cerr << "Unable to read the config file!\n";
	if(ReadFile(master_cnv_file,master_list_cnvs)){
		//cout << "Thanks! Now specify second file to read (Including path)\n";
		//cin >> filename;
		if(ReadFile(comparison_cnv_file,comparison_list_cnvs)){}
		else{
			//throw runtime_error("Cannot find config item.");
			
			cerr << "Looks like there was a problem reading " << comparison_cnv_file;
			cerr << "Make sure you updated the compare_lists.config file with the two files you want to compare\n";
			cerr << "Also, the cnv files need to be in the same directory as the program\n";
			//cerr << "\nPerhaps you didn't specify the path?\n";
			return 1;
		}
	}
	else{
		cerr << "Looks like there was a problem reading " << master_cnv_file;
		cerr << "Make sure you updated the compare_lists.config file with the two files you want to compare\n";
		cerr << "Also, the cnv files need to be in the same directory as the program\n";		
		//cerr << "\nPerhaps you didn't specify the path?\n";
		return 1;
	}
	//cout << "your overlap threshold is " << percent_overlap_cutoff << endl;
	
	
	cout << "Running compare_cnvs now\n";
	CompareCnvs(master_list_cnvs, comparison_list_cnvs);
	if(WriteCnvs(master_list_cnvs))
		cout << "All set. Check you're output file: " << out_file << endl;
		
	int num_cnvs = master_list_cnvs.size(),
		num_regions = comparison_list_cnvs.size();
	WriteLog(num_cnvs, num_regions);	
	//ConsolePrintCnvs(master_list_cnvs);
	return 0;
}
	
	