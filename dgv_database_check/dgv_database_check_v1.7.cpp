/*****************************************************
DGV_database_check_v1.7
Author: Robert Sicko - R.J.Sicko@gmail.com

****DESCRIPTION****
This program compares a list of CNVs to DGV2(or other) CNV database
and outputs the original list with an additional column appended.
'overlap_info'


Overlap_info:

		'% Our CNV'	= percentage of our CNV overlapped 
		'% DGV CNV' = percentage of DGV CNV overlapped
		'%Union' 	= percentage shared overlap in a union overlap calculation
						     ------|-----|---------------------   Our CNV
							 ----------+++++++++++++++---------   DGV2 CNV
							 ---++++++-------------------------   DGV2 CNV(2)							 
							 ----------XXX---------------------   Intersection
							 ------XXXXXXXXXXXXXXXXXXX---------   Denominator for union overlap
The DGV2 database can be found here: http://dgv.tcag.ca/dgv/app/home?ref=
****END DESCRIPTION****

****REQUIRED FORMAT****
	NOTE: No commas in numbers (ie. use 1000000 not 1,000,000)
	NOTE2: Keep headers in your files. The program ignores the first line (don't have spaces in header id's though).
	
	CNV List (tab-delimited txt):
	OurID	Chr		Start	Stop	Case/Control	Copy_Number
	
	DGV Database File(tab-delimited txt):
	DGV_date	chr	start	end	variantsubtype	frequency
****END REQUIRED FORMAT****

****CHANGE LOG****
v1.0- Bugs fixed, programs appears to be functioning as expected
v1.2- Input file for DGV database modified, now includes copy_number in last column
	- Now our check compares copy number to compare only dups with dups and dels with dels
v1.2- Added runlog output
v1.3- Output changed. Now six columns output:
		'% Our CNV'	= percentage of our CNV overlapped 
		'% DGV CNV' = percentage of DGV CNV overlapped
		'%Union' 	= percentage shared overlap in a union overlap calculation
		these three fields for both, same and different type of cnvs
v1.4- New output, more informative.
	- Similar overlap calculation as in pseudogene check.
	- Total % of our CNV covered by 1 or more DGV2 calls is now the output.
	- See overlap_info description above, previous versions didn't take into account that more than one DGV2
	  CNV could overlap our CNV. Now, we sum all overlap.
v1.5- Now use frequency data available for DGV2 entries
	- Mix the output form v1.3 and v1.4
v1.6- Fixed freq so the highest freq overlapped is output
v1.7- Fixed a bug - the highest frequency of the greatest overlap is now output instead of the overall highest frequency
****END CHANGE LOG****

****TO DO****
-fix freq so that it is the highest freq of a DGV2 overlap. 
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
#include <list>


using namespace std;

string VERSION = "v1.7";

struct OVERLAP{
	string id;
	long start, stop;
	//double freq;
	string from_block;
	bool is_overlapped;
};

bool OVERLAPStartSortPredicate(const OVERLAP& d1, const OVERLAP& d2)
{
	return d1.start < d2.start;
}

		/*//A struct for our CNVs (or regions)
		struct CNV{
			string our_id;	//our identifier
			int chr, copy_num, case_control;	//ints for chromosome, copy numbe and case or control (1, 2)
			long start, stop;	//start and end points for each cnv
			//int num_cases_overlap_with, num_controls_overlap_with; // count variables, +1 each time we overlap with a case or control
			double overlap_cnv_sametype, overlap_dgv_sametype, overlap_union_sametype,
				overlap_cnv_difftype, overlap_dgv_difftype, overlap_union_difftype;
			
			
		};*/


struct BLOCK{

	//	Block_ID	chr		start	stop	freq
	string block_id;	//block identifier
	string chr;	//ints for chromosome, copy numbe and case or control (1, 2)
	int copyNum;
	long start, stop;	//start and end points for each cnv
	double freq;
};


//A struct for our CNVs (or regions)
struct CNV{
	string our_id;	//our identifier
	string chr;
	int copy_num, case_control;	//ints for chromosome, copy number and case or control (1, 2)
	long start, stop;	//start and end points for each cnv
	//int num_cases_overlap_with, num_controls_overlap_with; // count variables, +1 each time we overlap with a case or control
	//string regions_overlapped, percent_overlapped;
	//int in_other;
	vector<OVERLAP> overlap_same;
	list<OVERLAP> coalesced_overlaps_same;
	
	vector<OVERLAP> overlap_diff;
	list<OVERLAP> coalesced_overlaps_diff;
	
	double percent_overlapped_same, percent_overlapped_diff;
	double overlap_single_same, overlap_single_diff, champ_freq_same, champ_freq_diff;
	
};


//Mostly for debugging
void ConsolePrint(vector<CNV> any_list_of_cnvs)
{
	//cout << "I got to ConsolePrintCnvs function\nAnd my cnvs list contains: ";
	//cout << any_list_of_cnvs.size();
	//cout << flush;
	cout << "OurID\tChr\tStart\tStop\tCase/Control\tCopy_Number\tOverlap_Regions\tPercent_Overlapped\n";
	for (vector<CNV>::iterator me = any_list_of_cnvs.begin();me != any_list_of_cnvs.end();me++)
			cout << me->our_id << "\t" << me->chr << "\t" 
				<< me->start << "\t" << me->stop << "\t" 
				<< me->case_control << "\t" << me->copy_num << endl;
								/*
								<< "\t" << me->regions_overlapped
								<< "\t" << me->percent_overlapped << endl;
								*/
		
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
		//string name;
		//int var1, var2, var3;
		ss >> dummy_cnv.our_id >> dummy_cnv.chr
			>> dummy_cnv.start >> dummy_cnv.stop
			>> dummy_cnv.case_control >> dummy_cnv.copy_num;

		dummy_cnv.percent_overlapped_same = 0;
		dummy_cnv.percent_overlapped_diff = 0;
		dummy_cnv.overlap_single_same = 0;
		dummy_cnv.overlap_single_diff = 0;
		dummy_cnv.champ_freq_same = 0;
		dummy_cnv.champ_freq_diff = 0;
		
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

bool ReadFileVariation(string filename,vector<BLOCK>& blocks)
{
	ifstream infile(filename.c_str());
	if(!infile.is_open())
		return false;
	BLOCK dummy_cnv;
	bool first_line=true;
	string line;
	while (getline(infile, line))	//get a line at a time
	{
			
		istringstream ss(line);
	
		//	DB file as of 131127 = DGV_date	chr	start	end	variantsubtype	frequency
		ss >> dummy_cnv.block_id >> dummy_cnv.chr	>> dummy_cnv.start
			>> dummy_cnv.stop >> dummy_cnv.copyNum >> dummy_cnv.freq;
		
	

		if(first_line)	//skip header
			first_line=false;
		else
			blocks.push_back (dummy_cnv);
		
	}
		
		//cout << "I got to the while loop for file reading\n";
	infile.close();
	//ConsolePrintCnvs(our_cnvs);
	return true;
}

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

void CoalesceOverlaps(list<OVERLAP>& overlap_regions)
{
	overlap_regions.sort(OVERLAPStartSortPredicate);
	
	//iterators for keeping track of the two nodes we are comparing
	list<OVERLAP>::iterator me = overlap_regions.begin(),
							thee = overlap_regions.begin(),
                               end = overlap_regions.end();
							   
	
	if ( me != end ) // Treat empty list
		thee++;	//sets it to the second element
		if(thee!=end)	//Treat list with one element
			while(thee != end)	//lets keep comparing until we are at the end
			{
				if(thee->start <= me->stop)	//hit to coalesce them
				{
					long temp_start = min(thee->start,me->start),temp_stop = max(thee->stop,me->stop);
					OVERLAP temp_region;
					temp_region.start = temp_start;
					temp_region.stop = temp_stop;

					overlap_regions.erase(me);
			
					list<OVERLAP>::iterator temp_itr = overlap_regions.erase(thee);
							
					me = overlap_regions.insert(temp_itr,temp_region);
					thee = temp_itr;
				}
				else{
					me++;
					thee++;
				}
			}
}

double CalcTotalOverlap(CNV our_cnv,list<OVERLAP>& overlap_regions)
{
	//cout << "Begining of CalcTotalOverlap function\n";
	double total_olap=0;
	for(list<OVERLAP>::iterator me = overlap_regions.begin();me != overlap_regions.end();me++)
	{
		double me_start = me->start,
				me_stop = me->stop,
				cnv_start = our_cnv.start,
				cnv_stop = our_cnv.stop;
		double temp_olap = CalcOverlapCNV(cnv_start,me_start,cnv_stop,me_stop);
		//cout << temp_olap << endl;
		total_olap = total_olap + temp_olap;
		//cout << total_olap << endl;
	}
	//cout << "End of CalcTotalOverlap function\n";
	return total_olap;
}

	

		


void CheckCNV(vector<CNV>& master_list_cnvs, vector<BLOCK> otherCNVs)
{
	//cout << type_specific << endl;
	for (vector<CNV>::iterator me = master_list_cnvs.begin();me != master_list_cnvs.end();me++)
	{
		for(vector<BLOCK>::iterator them = otherCNVs.begin();them != otherCNVs.end();them++)
		{
			if(me->chr == them->chr)
			{
				
				double me_start = me->start, me_stop = me->stop, them_start = them->start, them_stop = them->stop;	
				double overlap_cnv = CalcOverlapCNV(me_start,them_start,me_stop,them_stop);
			
				if(overlap_cnv>0)
				{	//If there is some overlap with our CNV, then lets keep track of it
				
					//cout << "there was an overlap" << endl;
					OVERLAP my_overlap;
					long dummy_start, dummy_stop;
					dummy_start = max (me->start,them->start);
					dummy_stop = min (me->stop,them->stop);
					my_overlap.start = dummy_start;
					my_overlap.stop = dummy_stop;
					my_overlap.is_overlapped = false;
					if(me->copy_num == them->copyNum || me->copy_num == them->copyNum - 1 || me->copy_num == them-> copyNum + 1){
						//if we match del to del or dup to dup
						me->overlap_same.push_back(my_overlap);
						me->percent_overlapped_same = me->percent_overlapped_same + overlap_cnv;
					}
					else{
						me->overlap_diff.push_back(my_overlap);
						me->percent_overlapped_diff = me->percent_overlapped_diff + overlap_cnv;
					}
				}
			}
		}
		//for same type overlaps
		//copy our overlap_bases vector to our coalesced_overlaps list
		copy(me->overlap_same.begin(), me->overlap_same.end(), back_inserter(me->coalesced_overlaps_same));
		
		//pass our coalesced_overlaps vector to the coalesce function
		CoalesceOverlaps(me->coalesced_overlaps_same);	
		
		//calculate our percent overlapped with the coalesced overlaps
		me->percent_overlapped_same = CalcTotalOverlap(*me,me->coalesced_overlaps_same);
		
		
		//for diff type overlaps
		//copy our overlap_bases vector to our coalesced_overlaps list
		copy(me->overlap_diff.begin(), me->overlap_diff.end(), back_inserter(me->coalesced_overlaps_diff));
		
		//pass our coalesced_overlaps vector to the coalesce function
		CoalesceOverlaps(me->coalesced_overlaps_diff);	
		
		//calculate our percent overlapped with the coalesced overlaps
		me->percent_overlapped_diff = CalcTotalOverlap(*me,me->coalesced_overlaps_diff);
		
		
	}
}

void CheckSingleOverlap(vector<CNV>& master_list_cnvs, vector<BLOCK> otherCNVs)
{
		
	//cout << type_specific << endl;
	for (vector<CNV>::iterator me = master_list_cnvs.begin();me != master_list_cnvs.end();me++)
		for(vector<BLOCK>::iterator them = otherCNVs.begin();them != otherCNVs.end();them++)
		{
			bool same_type = false;
			if(me->copy_num==(them->copyNum+1)
				||me->copy_num==(them->copyNum-1)
				||me->copy_num==them->copyNum)
						same_type = true;
			
			if(me->chr == them->chr 
				&& me->start == them->start
				&& me->stop == them->stop)
			{//exact match
				//cout << "them with 100%? " << them->our_id << endl;
				if(same_type)
				{
					if(me->overlap_single_same == 1)	//already have a 100% overlap, just check freq for a new champ
					{
						if(me->champ_freq_same < them->freq)
							me->champ_freq_same = them->freq;
					}
					else	//new champ for overlap
					{
						me->overlap_single_same = 1;
						me->champ_freq_same = them->freq;
					}
				}
				else	//different type
				{
					if(me->overlap_single_diff == 1)	//already have a 100% overlap, just check freq for a new champ
					{
						if(me->champ_freq_diff < them->freq)
							me->champ_freq_diff = them->freq;
					}
					else
					{
						me->overlap_single_diff = 1;
						me->champ_freq_diff = them->freq;
					}
				}
				continue;
			
			}
			else if(me->chr == them->chr)
			{					
				double overlap_CNV;
				double me_start = me->start, me_stop = me->stop, them_start = them->start, them_stop = them->stop;
				
				overlap_CNV = CalcOverlapCNV(me_start,them_start,me_stop,them_stop);
				if(overlap_CNV > 0)
				{//not exact match, but some overlap
					//cout << "them with overlap? " << "overlap = " << overlap_CNV << endl;
									
					if(same_type)
					{
						//cout << "and same type\n";
						//cout << "current champ overlap =" << me->overlap_single_same << endl;
						//cout << "current overlap here =" << overlap_CNV << endl;
						//if the overlap we are checking is the same as our current overlap, just check for greater frequency
						if(overlap_CNV == me->overlap_single_same)
						{
							//cout << "surely you're not counting these as equal..." << endl;
							if(me->champ_freq_same < them->freq)	//lets check their freq and replace ours with theirs if theirs is higher
								me->champ_freq_same = them->freq;
						}
						else if(overlap_CNV > me->overlap_single_same)
						{//we have a new champ
							//cout << "new champ\n";
							me->overlap_single_same = overlap_CNV;
							me->champ_freq_same = them->freq;	//new champ so set the freq
						}
					}
					else	//diff type
					{
						//cout << "and different type\n";
						if(overlap_CNV == me->overlap_single_diff)
						{
							if(me->champ_freq_diff < them->freq)	//lets check their freq and replace ours with theirs if theirs is higher
								me->champ_freq_diff = them->freq;
						}
						else if(overlap_CNV > me->overlap_single_diff)
						{//we have a new champ
							me->overlap_single_diff = overlap_CNV;
							me->champ_freq_diff = them->freq;	//new champ so set the freq
						}
						
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
	//added fields 	'% Our CNV'	'% DGV CNV' '%Union' 	
	
	outfile << "OurID\tChr\tStart\tStop\tCase/Control\tCopy_Number\t"
			<< "Total % Our CNV with all DGV2 overlap(same type)\t" << "Total % Our CNV with all DGV2 overlap(diff. type)\t"
			<< "% Our CNV with SINGLE DGV2(same type)\tSINGLE DGV2 Frequency(same type)\t% Our CNV with SINGLE DGV2(diff type)\t(SINGLE DGV2 Frequency(diff type)\n";

	/*	int overlap_cnv_sametype, overlap_dgv_sametype, overlap_union_sametype,
		overlap_cnv_difftype, overlap_dgv_difftype, overlap_union_difftype;*/
		
	for (vector<CNV>::iterator me = master_list_cnvs.begin();me != master_list_cnvs.end();me++)
	{
		outfile << me->our_id << "\t" << me->chr << "\t" 
				<< me->start << "\t" << me->stop << "\t"
				<< me->case_control << "\t" << me->copy_num << "\t";
		
		
		//double percent_overlapped_same, percent_overlapped_diff;
		
		if(me->percent_overlapped_same != 0){
			double olap_percent = me->percent_overlapped_same*100;
			ostringstream strs;
			strs << fixed << setprecision (2) << olap_percent;
			string olap = strs.str();	
			if(olap == "100.00") olap = "100";
			outfile << olap << "\t";
		}
		else
			outfile << "NA\t";
			
		if(me->percent_overlapped_diff != 0){
			double olap_percent = me->percent_overlapped_diff*100;
			ostringstream strs;
			strs << fixed << setprecision (2) << olap_percent;
			string olap = strs.str();	
			if(olap == "100.00") olap = "100";
			outfile << olap << "\t";
		}
		else	
			outfile << "NA\t";
		
		if(me->overlap_single_same != 0){
			double olap_percent = me->overlap_single_same*100;
			ostringstream strs;
			strs << fixed << setprecision (2) << olap_percent;
			string olap = strs.str();	
			if(olap == "100.00") olap = "100";
			outfile << olap << "\t";
			
			outfile << me->champ_freq_same << "\t";
			
		}
		else
			outfile << "NA\tNA\t";
		if(me->overlap_single_diff != 0){
			double olap_percent = me->overlap_single_diff*100;
			ostringstream strs;
			strs << fixed << setprecision (2) << olap_percent;
			string olap = strs.str();	
			if(olap == "100.00") olap = "100";
			outfile << olap << "\t";
			
			outfile << me->champ_freq_diff << "\t";
			
		}
		else
			outfile << "NA\tNA\t";

		outfile << "\n";
		
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
	outfile << "Program: dgv_database_check_" << VERSION << endl;
	outfile << "Number CNVs read: " << num_cnvs << endl;
	outfile << "Number DGV entries read: " << num_regions << endl;
	outfile << "CNV file used: " << master_cnv_file << endl;
	outfile << "DGV file used: " << comparison_cnv_file << endl;
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
	string in_name = "dgv_database_check_" + VERSION + ".config";
	ifstream infile(in_name.c_str());	//open file
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
			ss >> master_cnv_file;
		if(line_num ==6)
			ss >> comparison_cnv_file;
		if(line_num ==11)
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
	vector<BLOCK> cnv_database_list;

	string filename;
	cout << "***********************\n***********************\nDGV_database_check_" << VERSION << "\n***********************\n***********************\n";
	cout << "\nCNV list file should be a tab-delimited file in the following format:\n";
	cout << "OurID\tChr\tStart\tStop\tCase/Control\tCopy_Number\n";
	cout << "\nDGV (or any other) database file should be a tab-delimited file in the following format:";
	cout << "ID\tChr\tStart\tEnd\tCopy_num\tFrequency\n";
	cout << "\nKeep headers in your input files! WARNING: For now remove commas from numbers. ie dont use 1,000,000, use 1000000 in input files\n";
	cout << "I will try to deal with commas in numbers in later versions, but for now leave them out\n";
	cout << "Press enter to continue\n";
	cin.get();

	//cout << "Specify first file to read (Including path)\n";
	//cin >> filename;
	
	if(!DefineConfig()){
		cerr << "Unable to read the config file!\n";
		return 1;
	}
	if(ReadFileCnv(master_cnv_file,master_list_cnvs)){
		//cout << "Thanks! Now specify second file to read (Including path)\n";
		//cin >> filename;
		if(!ReadFileVariation(comparison_cnv_file,cnv_database_list)){
			//throw runtime_error("Cannot find config item.");
			
			cerr << "Looks like there was a problem reading " << comparison_cnv_file;
			cerr << " Make sure you updated the compare_lists.config file with the two files you want to compare\n";
			cerr << "Also, the cnv files need to be in the same directory as the program\n";
			//cerr << "\nPerhaps you didn't specify the path?\n";
			return 1;
		}
	}
	else{
		cerr << "Looks like there was a problem reading " << master_cnv_file;
		cerr << " Make sure you updated the compare_lists.config file with the two files you want to compare\n";
		cerr << "Also, the cnv files need to be in the same directory as the program\n";		
		//cerr << "\nPerhaps you didn't specify the path?\n";
		return 1;
	}
	//cout << "your overlap threshold is " << percent_overlap_cutoff << endl;
	
	//ConsolePrint(regions);
	cout << "Running compare_cnvs now\n";
	//ConsolePrint(cnv_database_list);
	CheckCNV(master_list_cnvs, cnv_database_list);
	CheckSingleOverlap(master_list_cnvs, cnv_database_list);
	if(WriteCnvs(master_list_cnvs))
		cout << "All set. Check you're output file: " << out_file << endl;
		
	int num_cnvs = master_list_cnvs.size(),
		num_regions = cnv_database_list.size();
	WriteLog(num_cnvs, num_regions);		
	//ConsolePrintCnvs(master_list_cnvs);
	return 0;
}
	
	