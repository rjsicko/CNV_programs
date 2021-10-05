/*****************************************************
Algorithm_Compare_v1.3
Author: Robert Sicko - R.J.Sicko@gmail.com

****DESCRIPTION****
This program compares output from seperate CNV calling algorithms
and writes a file with two appended columns. The two additional
columns are 'percents_overlapped' and 'total_overlap' for the percent
of our CNV overlapped by each other call and the total percent of
our CNV overlapped by all other calls. NOTE: these other calls
are in the same Subject_id!

****END DESCRIPTION****

****REQUIRED FORMAT****
	NOTE: No commas in numbers (ie. use 1000000 not 1,000,000).
	NOTE2: Keep headers in your files. The program ignores the first line (don't have spaces in header id's though).
	
	Both files should be tab-delimited in the following format(tab-delimited txt):
	OurID	Subject_id		Chr		Start	Stop	Case/Control	Copy_Number
****END REQUIRED FORMAT****


****CHANGE LOG****
v1.0- Bugs fixed, programs appears to be functioning as expected
v1.1- Output changed to include 0,1,2 or 3
		0=not in other algorithm
		1=exact match in other algorithm
		2=loose match in other algorithm (>=90% union overlap)
		3=exact match in other algorithm, but type differs (dup vs del)
		4=loose match in other algorithm, but type differs (dup vs del)
v1.2- Added runlog output	
v1.3- Overlap is now calculated as in the Pseudogene program.
	- Calculate total percent of our CNV overlapped by other calls (pileup other calls on our CNV)
	- Changed chr variable to string from int. Allows X, Y, XY to be read instead of converting to 23, 24, 25
****END CHANGE LOG****		
		

*********************************************************/




#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include "config.h"
#include <algorithm>
#include <sstream>
#include <iomanip>
#include <ctime>
#include <list>


using namespace std;

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

//A struct for our CNVs (or regions)
struct CNV{
	string our_id, subject_id;	//our identifier
	string chr;
	int copy_num, case_control;	//ints for chromosome, copy number and case or control (1, 2)
	long start, stop;	//start and end points for each cnv
	//int num_cases_overlap_with, num_controls_overlap_with; // count variables, +1 each time we overlap with a case or control
	//string regions_overlapped, percent_overlapped;
	//int in_other;
	vector<OVERLAP> overlap_bases;
	list<OVERLAP> coalesced_overlaps;
	
	double percent_overlapped;
	
};



struct BLOCK{

	//	Block_ID	chr		start	stop	freq
	string block_id;	//block identifier
	string chr;	//ints for chromosome, copy numbe and case or control (1, 2)
	int copyNum;
	long start, stop;	//start and end points for each cnv
	//double freq;
};


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

double CalcOverlapCNV(double start_cnv, double start_region, double end_cnv, double end_region)
{
	return ((min(end_cnv,end_region) - max(start_cnv,start_region))
							/
					(end_cnv-start_cnv));
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

bool ReadFileCnv(string filename,vector<CNV>& our_cnvs)
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
		ss >> dummy_cnv.our_id >> dummy_cnv.subject_id >> dummy_cnv.chr
			>> dummy_cnv.start >> dummy_cnv.stop
			>> dummy_cnv.case_control >> dummy_cnv.copy_num;
		//ss >> name >> var1 >> var2 >> var3;
		//first_line = false;
		//dummy_cnv.in_other=0;
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
bool ReadFileOtherCNVs(string filename,vector<BLOCK>& blocks)
{
	//OurID	Subject_id		Chr		Start	Stop	Case/Control	Copy_Number
	ifstream infile(filename.c_str());
	if(!infile.is_open())
		return false;
	BLOCK dummy_block;
	bool first_line=true;
	string line;
	while (getline(infile, line))
	{
			
		istringstream ss(line);
		//string name;
		//int var1, var2, var3;
		//	Block_ID	chr		start	stop
		string filler;
		ss >> filler;	//ignore ourID field, doesnt matter
		
		ss >> dummy_block.block_id >> dummy_block.chr	//subjectID then Chr
			>> dummy_block.start >> dummy_block.stop;	//start then stop
		
		ss >> filler;	//read past cast/control
		ss >> dummy_block.copyNum;	//copynum
		
		if(first_line)
			first_line=false;
		else
			blocks.push_back (dummy_block);
		
	}
		//cout << "I got to the while loop for file reading\n";
	infile.close();
	//ConsolePrintCnvs(our_cnvs);
	return true;
}

		//other overlap types not used
		/*double CalcOverlapUnion(double start_cnv, double start_region, double end_cnv, double end_region)
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
		}*/





				
			/* 		//old compare function (v1.2). replaced below
			//function to compare the master list with a test or comparison list
			//modifies master_list_cnvs 'num_cases_overlap_with' variable to keep track of the number of overlaps each entry in master list has with the comparison list 
			void CheckCNV(vector<CNV>& master_list_cnvs, vector<CNV> compare_list_cnvs)
			{
					
				//cout << type_specific << endl;
				for (vector<CNV>::iterator me = master_list_cnvs.begin();me != master_list_cnvs.end();me++)
					for(vector<CNV>::iterator them = compare_list_cnvs.begin();them != compare_list_cnvs.end();them++){
						//cout << compare_type << endl;
						if(me->chr == them->chr 
						&& me->subject_id == them->subject_id
						&& me->start == them->start
						&& me->stop == them->stop)
						{
								if(me->copy_num==them->copy_num||me->copy_num==them->copy_num+1||me->copy_num==them->copy_num-1)
										me->in_other=1;	//exact match 
								else 
									me->in_other==3;	//exact match but type differs
						}
						else if(me->chr == them->chr 
						&& me->subject_id == them->subject_id)
						{
								double overlap;
								double me_start = me->start, me_stop = me->stop, them_start = them->start, them_stop = them->stop;
								overlap = CalcOverlapUnion(me_start,them_start,me_stop,them_stop);
								//cout << overlap << endl;
								if(overlap >= 0.90)
									if(me->copy_num==them->copy_num||me->copy_num==them->copy_num+1||me->copy_num==them->copy_num-1)
										me->in_other=2;	//loose match, type same
									else
										me->in_other=4;	//loose match and type differs
						}
							/*
											//cout << me->our_id << "\t" << me->copy_num << endl;
											//cin.get();
											//double me_start = me->start, me_stop = me->stop, them_start = them->start, them_stop = them->stop;	
											//double overlap;
											if(overlap_type == 1)	overlap = CalcOverlapUnion(me_start,them_start,me_stop,them_stop);
											if(overlap_type == 2)	overlap = CalcOverlapRegion(me_start,them_start,me_stop,them_stop);
											if(overlap_type == 3)	overlap = CalcOverlapCNV(me_start,them_start,me_stop,them_stop);
											
								
											//cout << overlap << "&" << percent_overlap_cutoff << endl;
											//checks that the start and end of two CNVS are equal
										
											//cout << me->our_id << "\t" << me->copy_num << endl;
											//cin.get();
											if(overlap>percent_overlap_cutoff){
											me->regions_overlapped = me->regions_overlapped + them->unique_id;
											me->percent_overlapped = me->percent_overlapped + olap;
							*/

//function to compare the master list with a test or comparison list
//modifies master_list_cnvs 'num_cases_overlap_with' variable to keep track of the number of overlaps each entry in master list has with the comparison list 
void CheckCNV(vector<CNV>& master_list_cnvs, vector<BLOCK> otherCNVs)
{
	//cout << type_specific << endl;
	for (vector<CNV>::iterator me = master_list_cnvs.begin();me != master_list_cnvs.end();me++)
	{
		for(vector<BLOCK>::iterator them = otherCNVs.begin();them != otherCNVs.end();them++)
		{
			if(me->chr == them->chr && me->subject_id == them->block_id)
			{
				if(me->copy_num == them->copyNum) //only checking exact calls for now
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
						
						//my_overlap.from_block = them->block_id;
						my_overlap.start = dummy_start;
						my_overlap.stop = dummy_stop;
						//my_overlap.id = them->block_id;
						//my_overlap.freq = them->freq;
						my_overlap.is_overlapped = false;
						me->overlap_bases.push_back(my_overlap);
												
						me->percent_overlapped = me->percent_overlapped + overlap_cnv;
					}
				}
			}
		}
		//copy our overlap_bases vector to our coalesced_overlaps list
		copy(me->overlap_bases.begin(), me->overlap_bases.end(), back_inserter(me->coalesced_overlaps));
		
		//pass our coalesced_overlaps vector to the coalesce function
		CoalesceOverlaps(me->coalesced_overlaps);	
		
		//calculate our percent overlapped with the coalesced overlaps
		me->percent_overlapped = CalcTotalOverlap(*me,me->coalesced_overlaps);
		
		
		
	}
}



		//replaced by WriteCnvs below!
		/*bool WriteCnvs(vector<CNV> master_list_cnvs)
		{
			ofstream outfile(out_file.c_str());
			if(!outfile){
				 cerr << "Can not open the results file for writing!\nIt may be in use. Close it and try again\n";
				 return false;
			}
			outfile << "OurID\tSubject_id\tChr\tStart\tStop\tCase/Control\tCopy_Number\tCalled_in_Other_Algorithm\n";
			for (vector<CNV>::iterator me = master_list_cnvs.begin();me != master_list_cnvs.end();me++){
				outfile << me->our_id << "\t" << me->subject_id << "\t" << me->chr << "\t" 
						<< me->start << "\t" << me->stop << "\t" 
						<< me->case_control << "\t" << me->copy_num
						<< "\t" << me->in_other << endl;
			}
			outfile.close();	//close the file
			return true;
		}*/

bool WriteCnvs(vector<CNV> master_list_cnvs)
{
	//cout << "I got to writeCNVs" << endl;
	ofstream outfile(out_file.c_str());
	if(!outfile){
		 cerr << "Can not open the results file for writing!\nIt may be in use. Close it and try again\n";
		 return false;
	}
	outfile << "OurID\tSubjectID\tChr\tStart\tStop\tCase/Control\tCopy_Number\toverlap_percents\ttotal_overlap\n";
	
	for (vector<CNV>::iterator me = master_list_cnvs.begin();me != master_list_cnvs.end();me++)
	{
		int region_count=0;
		outfile << me->our_id << "\t" << me->subject_id << "\t" << me->chr << "\t" 
				<< me->start << "\t" << me->stop << "\t" 
				<< me->case_control << "\t" << me->copy_num << "\t";
				
		//'percents' and 'total_overlap'
		string overlap_freq, overlap_percent;
		//cout << me->percent_overlapped << endl;
		if(me->percent_overlapped>0)
		{
			for(vector<OVERLAP>::iterator my_overlap = me->overlap_bases.begin();
				my_overlap != me->overlap_bases.end();
				my_overlap++)	
			{	//for all overlapped bases in our CNV
				/*if(overlapped_ids != "")
					overlapped_ids = overlapped_ids + ";";
				overlapped_ids = overlapped_ids + my_overlap->from_block;
				*/
				
				if(overlap_percent!= "")
					overlap_percent = overlap_percent + ";";
				
					
				//convert overlap to a string
				double temp_overlap = CalcOverlapCNV(me->start, my_overlap->start, me->stop, my_overlap->stop);
				double olap_percent = temp_overlap*100;
				ostringstream strs;
				strs << fixed << setprecision (2) << olap_percent;
				string olap = strs.str();	
				if(olap == "100.00") olap = "100";
				
				overlap_percent = overlap_percent + olap;
			}
				
			ostringstream strs_nother_temp;
			strs_nother_temp << fixed << setprecision (2) << me->percent_overlapped*100;
			string temp_total_overlap = strs_nother_temp.str();
			outfile << overlap_percent << "\t" << temp_total_overlap << endl;
		}
		else
			outfile << "NA\tNA" << endl;

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
	outfile << "Program: algorithm_compare_v1.3\n";
	outfile << "Number CNVs read: " << num_cnvs << endl;
	outfile << "Number comparison CNVs read: " << num_regions << endl;
	outfile << "CNV file used: " << master_cnv_file << endl;
	outfile << "Comparison CNVs file used: " << comparison_cnv_file << endl;
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
	ifstream infile("algorithm_compare_v1.3.config");	//open file
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
	vector<BLOCK> compare_list_cnvs;

	string filename;
	cout << "***********************\n***********************\ncompare_algorithms v1.3\n***********************\n***********************\n";
	cout << "\nCNV list files should be a tab-delimited file in the following format(NOTE ADDITIONAL 'subject_id' FIELD):\n";
	cout << "OurID\tSubject_id\tChr\tStart\tStop\tCase/Control\tCopy_Number\n";
	cout << "\nKeep headers in your input files! WARNING: For now remove commas from numbers. ie dont use 1,000,000, use 1000000 in input files\n";
	cout << "I will try to deal with commas in numbers in later versions, but for now leave them out\n";
	cout << "Press enter to continue\n";
	cin.get();

	//cout << "Specify first file to read (Including path)\n";
	//cin >> filename;
	
	if(!DefineConfig())	cerr << "Unable to read the config file!\n";
	if(ReadFileCnv(master_cnv_file,master_list_cnvs)){
		//cout << "Thanks! Now specify second file to read (Including path)\n";
		//cin >> filename;
		if(!ReadFileOtherCNVs(comparison_cnv_file,compare_list_cnvs)){
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
	CheckCNV(master_list_cnvs, compare_list_cnvs);
	if(WriteCnvs(master_list_cnvs))
		cout << "All set. Check you're output file: " << out_file << endl;
		
	int num_cnvs = master_list_cnvs.size(),
		num_regions = compare_list_cnvs.size();
	WriteLog(num_cnvs, num_regions);		
	//ConsolePrintCnvs(master_list_cnvs);
	return 0;
}
	
	