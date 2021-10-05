/*****************************************************
pseudo_gene_check_v1.2
Author: Robert Sicko - R.J.Sicko@gmail.com

****DESCRIPTION****
This program compares a list of CNVs to a list of genes(pseudo in this case).
This program compares each CNV in the list to each gene and will output which
CNVs are overlapped(above a threshold) by a gene or multiple genes.

The four columns appended are 'gene_ids' 'gene_percents' and 'total_overlap'
if a CNV overlaps more than one gene, the gene_ids will be seperated
by a ; as will the percent overlaps for each
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
v1.1- Fixed a bug in the 'total_overlap' calculation
	- v1.0 actually addressed the issue that 'total_overlap' could be greater than 100% if you simply sum
	  each OVERLAP in the vector of overlaps. However, the way v1.0 corrected it was also faulty. 
	  It calculated the correct overlap the majority of the time, but over estimated overlap when multiple 
	  seperate overlaps overlapped with each other.
	- To see how this was corrected, see function CoalesceOverlaps... briefly, after all overlaps
	  to a particular CNV were added to the vector of OVERLAPS, these overlaps were coalesced,
	  leaving non overlapping OVERLAPS that were stored in a list. This list was then passed
	  to CalcTotalOverlap for calculating total overlap with our CNV.
v1.2- Added runlog output

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
	string our_id;	//our identifier
	int chr, copy_num, case_control;	//ints for chromosome, copy numbe and case or control (1, 2)
	long start, stop;	//start and end points for each cnv
	//int num_cases_overlap_with, num_controls_overlap_with; // count variables, +1 each time we overlap with a case or control
	//string regions_overlapped, percent_overlapped;
	//multimap<string,string> genes_overlapped;	//this multimap is <gene,transcript>
	//map<string,string> percent_overlapped;	//this map is <transcript,percent_overlap>
	vector<OVERLAP> overlap_bases;
	list<OVERLAP> coalesced_overlaps;
	
	double percent_overlapped;
	
};


struct BLOCK{

	//	Block_ID	chr		start	stop	freq
	string block_id;	//block identifier
	int chr;	//ints for chromosome, copy numbe and case or control (1, 2)
	long start, stop;	//start and end points for each cnv
	//double freq;
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

		//cout << dummy_cnv.our_id << "\t" << dummy_cnv.chr << "\t" << dummy_cnv.start << endl;
		dummy_cnv.percent_overlapped = 0;
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

bool ReadFileRegion(string filename,vector<BLOCK>& blocks)
{
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
		//	database_from(or_type)	chr		start	stop	
		string filler;
		ss >> filler;	//ignore database_from field
		
		ss >> dummy_block.block_id >> dummy_block.chr
			>> dummy_block.start >> dummy_block.stop;

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



		/*double CheckOverlapsForOverlap(vector<OVERLAP> overlap_bases, long cnv_length)
		{
				//cout << "got to begining of checkoverlaps" << endl;
			double total_overlap=0;
			for(vector<OVERLAP>::iterator me = overlap_bases.begin();me != overlap_bases.end();me++)
			{	
				if(me->is_overlapped)
					continue;
				for(vector<OVERLAP>::iterator thee = overlap_bases.begin();thee != overlap_bases.end();thee++)
				{
					if(me->id == thee->id)
						continue;
						//cout << thee->start << endl;
					double first_start = me->start,
							second_start = thee->start,
							first_stop = me->stop,
							second_stop = thee->stop;
					//cout << "first start " << first_start << endl;
					double temp_olap = (min(first_stop,second_stop) - max(first_start,second_start))
										/
										(cnv_length);
					//cout << "min " << min(first_stop,second_stop) << endl;
					//cout << "Temp_olap " << temp_olap << endl;
					if(temp_olap >0)
					{
						total_overlap = total_overlap + temp_olap;
						thee->is_overlapped = true;
					}
				}
				
			}
				//cout << "got to end of checkoverlaps" << endl;
			//cout << "total_overlap " << total_overlap << endl;
			return total_overlap;
		}*/
	

//function to compare the master list with a test or comparison list
//modifies master_list_cnvs 'num_cases_overlap_with' variable to keep track of the number of overlaps each entry in master list has with the comparison list 
void CheckRegions(vector<CNV>& master_list_cnvs, vector<BLOCK> blocks)
{
	//cout << type_specific << endl;
	for (vector<CNV>::iterator me = master_list_cnvs.begin();me != master_list_cnvs.end();me++)
	{
		for(vector<BLOCK>::iterator them = blocks.begin();them != blocks.end();them++)
		{
			if(me->chr==them->chr)
			{
				double me_start = me->start, me_stop = me->stop, them_start = them->start, them_stop = them->stop;	
				double overlap_cnv = CalcOverlapCNV(me_start,them_start,me_stop,them_stop);
					//if(overlap_type == 1)	overlap = CalcOverlapUnion(me_start,them_start,me_stop,them_stop);
					//if(overlap_type == 2)	overlap = CalcOverlapRegion(me_start,them_start,me_stop,them_stop);
			
				if(overlap_cnv>0)
				{	//If there is some overlap with our CNV, then lets keep track of it
				
					//cout << "there was an overlap" << endl;
					OVERLAP my_overlap;
					long dummy_start, dummy_stop;
					dummy_start = max (me->start,them->start);
					dummy_stop = min (me->stop,them->stop);
					
					my_overlap.from_block = them->block_id;
					my_overlap.start = dummy_start;
					my_overlap.stop = dummy_stop;
					my_overlap.id = them->block_id;
					//my_overlap.freq = them->freq;
					my_overlap.is_overlapped = false;
					me->overlap_bases.push_back(my_overlap);
											
					me->percent_overlapped = me->percent_overlapped + overlap_cnv;

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



bool WriteCnvs(vector<CNV> master_list_cnvs)
{
	//cout << "I got to writeCNVs" << endl;
	ofstream outfile(out_file.c_str());
	if(!outfile){
		 cerr << "Can not open the results file for writing!\nIt may be in use. Close it and try again\n";
		 return false;
	}
	outfile << "OurID\tChr\tStart\tStop\tCase/Control\tCopy_Number\tgene_ids\toverlap_percents\ttotal_overlap\n";
	
	for (vector<CNV>::iterator me = master_list_cnvs.begin();me != master_list_cnvs.end();me++)
	{
		int region_count=0;
		outfile << me->our_id << "\t" << me->chr << "\t" 
				<< me->start << "\t" << me->stop << "\t" 
				<< me->case_control << "\t" << me->copy_num << "\t";
				
		//'gene_ids' 'gene_percents' and 'total_overlap'
		string overlapped_ids, overlap_freq, overlap_percent;
		//cout << me->percent_overlapped << endl;
		if(me->percent_overlapped>percent_overlap_cutoff)
		{	//the total overlap is greater than our cutoff
			for(vector<OVERLAP>::iterator my_overlap = me->overlap_bases.begin();
				my_overlap != me->overlap_bases.end();
				my_overlap++)	
			{	//for all overlapped bases in our CNV
				if(overlapped_ids != "")
					overlapped_ids = overlapped_ids + ";";
				overlapped_ids = overlapped_ids + my_overlap->from_block;
				
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
			outfile << overlapped_ids << "\t" << overlap_percent << "\t" << temp_total_overlap << endl;
		}
		else
			outfile << "NA\tNA\tNA" << endl;

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
	outfile << "Program: pseudo_gene_check_v1.2\n";
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
	ifstream infile("pseudo_gene_check_v1.2.config");	//open file
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
	vector<BLOCK> blocks;

	string filename;
	cout << "***********************\n***********************\nPseudo Gene Check_v1.2\n***********************\n***********************\n";
	cout << "\nCNV list file should be a tab-delimited file in the following format:\n";
	cout << "OurID\tChr\tStart\tStop\tCase/Control\tCopy_Number\n";
	cout << "\nCHOP block list file should be a tab-delimited file in the following format:";
	cout << "BlockID\tChr\tStart\tStop\tFreq\n";
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
		if(!ReadFileRegion(comparison_cnv_file,blocks)){
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
	CheckRegions(master_list_cnvs, blocks);
	if(WriteCnvs(master_list_cnvs))
		cout << "All set. Check you're output file: " << out_file << endl;
		
	int num_cnvs = master_list_cnvs.size(),
		num_regions = blocks.size();
	WriteLog(num_cnvs, num_regions);		
	//ConsolePrintCnvs(master_list_cnvs);
	return 0;
}
	
	