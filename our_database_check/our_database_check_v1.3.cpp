/*****************************************************
our_database_check_v1.3
Author: Robert Sicko - R.J.Sicko@gmail.com

****DESCRIPTION****
program to check overlap with our CNV database. For each CNV input we output the CNV with these additional columns:
'Exact_Unique_Times_Seen'
the following two columns:
	'Total % Our CNV With Cumulative Database Overlap(same phenotype, same CNV type)'
	'Cumulative Database CNVs Overlapped With (same phenotype, same CNV type)'
- these two columns repeat with:
	'similar phenotype, same CNV type'
	'different phenotype, same CNV type'
	'same phenotype, different CNV type'
	'similar phenotype, different CNV type'
	'different phenotype, different CNV type'

****END DESCRIPTION****

****REQUIRED FORMAT****
	NOTE: No commas in numbers (ie. use 1000000 not 1,000,000)
	NOTE2: Keep headers in your files. The program ignores the first line (don't have spaces in header id's though).
	
	File should be tab-delimited in the following format:
	ourID	sample	phenotype	case/control	chr	start	stop	copy_num
****END REQUIRED FORMAT****

****CHANGE LOG****
v1.0- Bugs fixed, programs appears to be functioning as expected
v1.1- Added runlog output
	- Added options 5-8 to output, see above.
	- Added output for Number_100, Number_90, Number_50, Number_1BP
	  instead of a single Number_overlapped column.
v1.2- New output, more informative.
	- Similar overlap calculation as in pseudogene check and dgv2 check.
	- Total % of our CNV covered by 1 or more cumulative database calls is now the output.
	- Previous versions didn't take into account that more than one cumulative DB
	  CNV could overlap our CNV. Now, we sum all overlap.
	- Also, ignore self in cumulative DB.
	  That is, don't count 100% overlap just because we matched ourself in the database.
	- Finally, now have three fields for overlap: exact, similar, different
	  These types will take into account the type of CNV.
v1.3- source completely changed
	- source a merger of case_control_v1.2 and CheckCumulativeDatabaseOverlap function from GVAT-0.0.2 and new additions and cleanup
	- output changed again. We now output 13 columns:
		- first column still 'Exact_Unique_Times_Seen'
		- remaining 12 columns are the following two columns:
			'Total % Our CNV With Cumulative Database Overlap(same phenotype, same CNV type)'
			'Cumulative Database CNVs Overlapped With (same phenotype, same CNV type)'
		- these two columns repeat with:
			'similar phenotype, same CNV type'
			'different phenotype, same CNV type'
			'same phenotype, different CNV type'
			'similar phenotype, different CNV type'
			'different phenotype, different CNV type'

****END CHANGE LOG****

****TO DO****

****END TO DO****

*********************************************************/



#include <iostream>
#include <vector>
#include <list>
#include <fstream>
#include <string>
#include "config.h"
#include <algorithm>
#include <sstream>
#include <ctime>
#include <iomanip>

using namespace std;
const string VERSION = "our_database_check_v1.3";



//A struct for our CNVs (used for both the CNV file and the cumulative db file)
struct CNV{
	//cumulative DB format: Our_CNV_ID	Study_ID	Disorder	Case_or_Control	Chromosome	Start	End	Copy_Number
	//GVAT format: ourID	subjectID	phenotype	case/control	chr	start	stop	copy_num
	string our_id, study_id, phenotype, chr;
	int copy_num, case_or_control;	
	double start, stop;	//start and end points for each cnv

	int exact_unique_times_seen;
	//variables for overlap percentages	
	double percent_olap_same_disorder_same_type, percent_olap_similar_disorder_same_type, percent_olap_diff_disorder_same_type;
	double percent_olap_same_disorder_diff_type, percent_olap_similar_disorder_diff_type, percent_olap_diff_disorder_diff_type;
	
	//variables for the id's we overlap
	string ids_olap_same_disorder_same_type, ids_olap_similar_disorder_same_type, ids_olap_diff_disorder_same_type;
	string ids_olap_same_disorder_diff_type, ids_olap_similar_disorder_diff_type, ids_olap_diff_disorder_diff_type;
	
};

//struct for overlaps
struct OVERLAP{
	double start, stop;
};

bool OVERLAPStartSortPredicate(const OVERLAP& d1, const OVERLAP& d2)
{
	return d1.start < d2.start;
}


bool ReadFile(string filename,vector<CNV>& some_cnvs)
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
		//ourID	subjectID	phenotype	case/control	chr	start	stop	copy_num
		ss >> dummy_cnv.our_id >> dummy_cnv.study_id
			>> dummy_cnv.phenotype >> dummy_cnv.case_or_control
			>> dummy_cnv.chr >> dummy_cnv.start
			>> dummy_cnv.stop >> dummy_cnv.copy_num;
		
		dummy_cnv.exact_unique_times_seen = 0;
		//set our percent overlap to 0
		dummy_cnv.percent_olap_same_disorder_same_type = 0;
		dummy_cnv.percent_olap_similar_disorder_same_type = 0;
		dummy_cnv.percent_olap_diff_disorder_same_type = 0;
		dummy_cnv.percent_olap_same_disorder_diff_type = 0;
		dummy_cnv.percent_olap_similar_disorder_diff_type = 0;
		dummy_cnv.percent_olap_diff_disorder_diff_type = 0;
		
		//no need to set overlapped id's to "", by default a string is empty
		
		if(first_line)
			//if we are reading the header line, change our flag for first_line so we know next run through we aren't on the header line
			first_line=false;
		else
			//wasn't a header so push it back
			some_cnvs.push_back(dummy_cnv);
	}
		//cout << "I got to the while loop for file reading\n";
	infile.close();
	return true;
}

//returns a formatted string for outputting overlap from a double containing the calculated overlap
string ConvertOlapToPercent(double olap)
{	
	double percent_olap = olap * 100;
	ostringstream strs;
	strs << fixed << setprecision (2) << percent_olap;
	string final_olap = strs.str();	
	if(final_olap == "100.00") final_olap = "100";
	return final_olap;
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
							   	
	if ( me != end ) // Treat empty list, if not empty it's safe to make thee the second element
		thee++;	//sets it to the second element
	if(thee!=end)	//Treat list with one element, if its not a list with one element its safe to do the while loop
		while(thee != end)	//lets keep comparing until we are at the end
		{
			if(thee->start <= me->stop)	//hit to coalesce them
			{
				double temp_start = min(thee->start,me->start), temp_stop = max(thee->stop,me->stop);
				OVERLAP temp_region;
				temp_region.start = temp_start;
				temp_region.stop = temp_stop;

				overlap_regions.erase(me);
		
				list<OVERLAP>::iterator temp_itr = overlap_regions.erase(thee);
						
				me = overlap_regions.insert(temp_itr,temp_region);
				thee = temp_itr;
			}
			else
			{
				me++;
				thee++;
			}
		}
}
	
//Here we calculate total overlap percentage from our coalesced overlaps
double CalcTotalOverlap(CNV our_cnv,list<OVERLAP>& coalesced_regions)
{
	//cout << "Begining of CalcTotalOverlap function\n";
	double total_olap=0;
	for(list<OVERLAP>::iterator me = coalesced_regions.begin();me != coalesced_regions.end();me++)
	{
		double temp_olap = CalcOverlapCNV(our_cnv.start,me->start,our_cnv.stop,me->stop);
		total_olap = total_olap + temp_olap;
	}
	return total_olap;
}	

//function to compare the master list with a test or comparison list
void CompareCnvs(vector<CNV>& master_list_cnvs, vector<CNV>& comparison_list_cnvs)
{		
	for (vector<CNV>::iterator me = master_list_cnvs.begin();me != master_list_cnvs.end();me++)
	{
	//each CNV in our list
		//make our temporary overlap lists
		list<OVERLAP> my_overlaps_same_pheno_same_cnv, my_overlaps_similar_pheno_same_cnv, my_overlaps_diff_pheno_same_cnv;
		list<OVERLAP> my_overlaps_same_pheno_diff_cnv, my_overlaps_similar_pheno_diff_cnv, my_overlaps_diff_pheno_diff_cnv;
	
		for(vector<CNV>::iterator them = comparison_list_cnvs.begin();them != comparison_list_cnvs.end();them++)
		{
		//for each cumulative db entry
			if(me->study_id == them->study_id)	//lets not count ourself and skip if we are looking at ourself
				continue;
			if(me->chr != them->chr)	//limit what we look at based on chromosome
				continue;
			double overlap = CalcOverlapCNV(me->start, them->start, me->stop, them->stop);								
			if(overlap > 0)	//again limit what we look at, this time only look at cnvs that overlap
			{
			//lets not compare LOH to CNV
				if( (me->copy_num == 2 && them->copy_num != 2) || (them->copy_num == 2 && me->copy_num !=2) )
					continue;
			//lets set our similar phenotype flag
				bool similar_phenotype_flag = false;
				for(vector<string>::iterator itr = similar_phenotypes.begin();
					itr != similar_phenotypes.end();itr++)
				{
					size_t found = itr->find(me->phenotype.c_str());
					if(found != string::npos)	//if we are in one of the strings in the similar phenotypes vector
					{
						found = itr->find(them->phenotype.c_str());
						if(found != string::npos)	//if them phenotype was in the same similar phenotype string as "me" was in
							similar_phenotype_flag = true;
					}
				}					
			//lets set our same cnv type flag
				bool cnv_same_flag = true;
				if(me->copy_num==them->copy_num	//ok we don't check loh vs cnv here because we already continued on that earlier
					||me->copy_num==(them->copy_num+1)
					||me->copy_num==(them->copy_num-1))
				{
					cnv_same_flag = true;

				}
				else
					cnv_same_flag = false;
			
			//now that we have our flags, we know which list to add our overlap to
				OVERLAP temp_olap;
				temp_olap.start = them->start;
				temp_olap.stop = them->stop;
	/*
	double percent_olap_with_same_disorder_same_type, percent_olap_other_similar_disorder_same_type, percent_olap_other_diff_disorder_same_type;
	double percent_olap_with_same_disorder_diff_type, percent_olap_other_similar_disorder_diff_type, percent_olap_other_diff_disorder_diff_type;
	
	//variables for the id's we overlap
	string ids_olap_with_same_disorder_same_type, ids_olap_other_similar_disorder_same_type, ids_olap_other_diff_disorder_same_type;
	string ids_olap_with_same_disorder_diff_type, ids_olap_other_similar_disorder_diff_type, ids_olap_other_diff_disorder_diff_type;
	*/				
				if(me->phenotype == them->phenotype)	//same phenotype
				{
					if(cnv_same_flag)
					{
						my_overlaps_same_pheno_same_cnv.push_back(temp_olap);
						if(me->ids_olap_same_disorder_same_type != "")
							me->ids_olap_same_disorder_same_type = me->ids_olap_same_disorder_same_type + ";";
						me->ids_olap_same_disorder_same_type = me->ids_olap_same_disorder_same_type + them->our_id;
					}
					else
					{
						my_overlaps_same_pheno_diff_cnv.push_back(temp_olap);
						if(me->ids_olap_same_disorder_diff_type != "")
							me->ids_olap_same_disorder_diff_type = me->ids_olap_same_disorder_diff_type + ";";
						me->ids_olap_same_disorder_diff_type = me->ids_olap_same_disorder_diff_type + them->our_id;
					}
				}
				else if(similar_phenotype_flag)	//similar phenotype
				{
					if(cnv_same_flag)
					{
						my_overlaps_similar_pheno_same_cnv.push_back(temp_olap);
						if(me->ids_olap_similar_disorder_same_type != "")
							me->ids_olap_similar_disorder_same_type = me->ids_olap_similar_disorder_same_type + ";";
						me->ids_olap_similar_disorder_same_type = me->ids_olap_similar_disorder_same_type + them->our_id;
					}
					else
					{
						my_overlaps_similar_pheno_diff_cnv.push_back(temp_olap);
						if(me->ids_olap_similar_disorder_diff_type != "")
							me->ids_olap_similar_disorder_diff_type = me->ids_olap_similar_disorder_diff_type + ";";
						me->ids_olap_similar_disorder_diff_type = me->ids_olap_similar_disorder_diff_type + them->our_id;						
					}
				}
				else	//different phenotype
				{
					if(cnv_same_flag)
					{
						my_overlaps_diff_pheno_same_cnv.push_back(temp_olap);
						if(me->ids_olap_diff_disorder_same_type != "")
							me->ids_olap_diff_disorder_same_type = me->ids_olap_diff_disorder_same_type + ";";
						me->ids_olap_diff_disorder_same_type = me->ids_olap_diff_disorder_same_type + them->our_id;							
					}
					else
					{
						my_overlaps_diff_pheno_diff_cnv.push_back(temp_olap);
						if(me->ids_olap_diff_disorder_diff_type != "")
							me->ids_olap_diff_disorder_diff_type = me->ids_olap_diff_disorder_diff_type + ";";
						me->ids_olap_diff_disorder_diff_type = me->ids_olap_diff_disorder_diff_type + them->our_id;							
					}
				}
			//now lets check for exact unique
				//exact breakpoints check
				if(me->start == them->start && me->stop == them->stop && me->copy_num==them->copy_num)
					//now don't count if a same phenotype 
					//don't count if similar phenotype either
					if(me->phenotype != them->phenotype && !similar_phenotype_flag)
					{
						//cout << "we've got an exact unique\n";
						me->exact_unique_times_seen++;
					}
			}
		}
	//now lets deal with collapsing the overlaps
		//pass our overlaps to the coalesce function
		CoalesceOverlaps(my_overlaps_same_pheno_same_cnv);
		CoalesceOverlaps(my_overlaps_similar_pheno_same_cnv);
		CoalesceOverlaps(my_overlaps_diff_pheno_same_cnv);
		//calculate our percent overlapped with the coalesced overlaps
		me->percent_olap_same_disorder_same_type = CalcTotalOverlap(*me,my_overlaps_same_pheno_same_cnv);
		me->percent_olap_similar_disorder_same_type = CalcTotalOverlap(*me,my_overlaps_similar_pheno_same_cnv);
		me->percent_olap_diff_disorder_same_type = CalcTotalOverlap(*me,my_overlaps_diff_pheno_same_cnv);
		
		//pass our overlaps to the coalesce function
		CoalesceOverlaps(my_overlaps_same_pheno_diff_cnv);
		CoalesceOverlaps(my_overlaps_similar_pheno_diff_cnv);
		CoalesceOverlaps(my_overlaps_diff_pheno_diff_cnv);
		//calculate our percent overlapped with the coalesced overlaps
		me->percent_olap_same_disorder_diff_type = CalcTotalOverlap(*me,my_overlaps_same_pheno_diff_cnv);
		me->percent_olap_similar_disorder_diff_type = CalcTotalOverlap(*me,my_overlaps_similar_pheno_diff_cnv);
		me->percent_olap_diff_disorder_diff_type = CalcTotalOverlap(*me,my_overlaps_diff_pheno_diff_cnv);
		
		//check for those with no overlap
		if(me->ids_olap_same_disorder_same_type == "")
			me->ids_olap_same_disorder_same_type = "NA";
		if(me->ids_olap_similar_disorder_same_type == "")
			me->ids_olap_similar_disorder_same_type = "NA";
		if(me->ids_olap_diff_disorder_same_type == "")
			me->ids_olap_diff_disorder_same_type = "NA";
			
		if(me->ids_olap_same_disorder_diff_type == "")
			me->ids_olap_same_disorder_diff_type = "NA";
		if(me->ids_olap_similar_disorder_diff_type == "")
			me->ids_olap_similar_disorder_diff_type = "NA";
		if(me->ids_olap_diff_disorder_diff_type == "")
			me->ids_olap_diff_disorder_diff_type = "NA";
		
	}				
}

bool WriteCnvs(vector<CNV> master_list_cnvs)
{
	ofstream outfile(out_file.c_str());
	if(!outfile)
		 return false;
		 
	//ourID\tstudyID\tphenotype\tcase/control(1,2)\tchr\tstart\tstop\tcopy_num\n";
	//output headers
	outfile << "ourID\tstudyID\tphenotype\tcase/control(1,2)\tchr\tstart\tstop\tcopy_num\t"	 //original fields
			
			<< "Number exact unique times seen (other disorders)\t"
			
			<< "Total % Our CNV With Cumulative Database Overlap(same phenotype, same CNV type)\t"
			<< "Cumulative Database CNVs Overlapped With(same phenotype, same CNV type)\t"
			
			<< "Total % Our CNV With Cumulative Database Overlap(similar phenotype, same CNV type)\t"
			<< "Cumulative Database CNVs Overlapped With(similar phenotype, same CNV type)\t"
			
			<< "Total % Our CNV With Cumulative Database Overlap(different phenotype, same CNV type)\t"
			<< "Cumulative Database CNVs Overlapped With(different phenotype, same CNV type)\t"
			
			
			<< "Total % Our CNV With Cumulative Database Overlap(same phenotype, different CNV type)\t"
			<< "Cumulative Database CNVs Overlapped With(same phenotype, different CNV type)\t"
			
			<< "Total % Our CNV With Cumulative Database Overlap(similar phenotype, different CNV type)\t"
			<< "Cumulative Database CNVs Overlapped With(similar phenotype, different CNV type)\t"
			
			<< "Total % Our CNV With Cumulative Database Overlap(different phenotype, different CNV type)\t"
			<< "Cumulative Database CNVs Overlapped With(different phenotype, different CNV type)\n";
			
	for (vector<CNV>::iterator me = master_list_cnvs.begin();me != master_list_cnvs.end();me++)
	{
	
		//"ourID\tstudyID\tphenotype\tcase/control(1,2)\tchr\tstart\tstop\tcopy_num\t
		outfile << me->our_id << "\t" << me->study_id << "\t" 
				<< me->phenotype << "\t" << me->case_or_control << "\t" 
				<< me->chr << "\t" << me->start << "\t" << me->stop
				<< "\t" << me->copy_num << "\t" 
				
				<< me->exact_unique_times_seen << "\t"
			
				<< ConvertOlapToPercent(me->percent_olap_same_disorder_same_type) << "\t"
				<< me->ids_olap_same_disorder_same_type << "\t"
				<< ConvertOlapToPercent(me->percent_olap_similar_disorder_same_type) << "\t"
				<< me->ids_olap_similar_disorder_same_type << "\t"
				<< ConvertOlapToPercent(me->percent_olap_diff_disorder_same_type) << "\t"
				<< me->ids_olap_diff_disorder_same_type << "\t"
				
				<< ConvertOlapToPercent(me->percent_olap_same_disorder_diff_type) << "\t"
				<< me->ids_olap_same_disorder_diff_type << "\t"
				<< ConvertOlapToPercent(me->percent_olap_similar_disorder_diff_type) << "\t"
				<< me->ids_olap_similar_disorder_diff_type << "\t"
				<< ConvertOlapToPercent(me->percent_olap_diff_disorder_diff_type) << "\t"
				<< me->ids_olap_diff_disorder_diff_type << "\n";			
	}
	
	
	outfile.close();
	return true;
}

bool WriteLog(int num_cnvs, int num_cumulative)
{
	string out_trans = out_file + ".Log";
	ofstream outfile(out_trans.c_str());
	if(!outfile){
		 cerr << "Can not open the log file for writing!\nIt may be in use. Close it and try again\n";
		 return false;
	}
	outfile << "Program: " << VERSION << endl;
	outfile << "Number CNVs read: " << num_cnvs << endl;
	outfile << "CNV file used: " << master_cnv_file << endl;
	outfile << "Number cumulative database entries read: "<< num_cumulative << endl;
	outfile << "Cumulative database file used: " << cumulative_database_file << endl;
	outfile << "Phenotypes called similar: ";
	for(vector<string>::iterator itr = similar_phenotypes.begin();
		itr != similar_phenotypes.end();itr++)
	{
		outfile << *itr;
		outfile << "\n";
	}
	

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
	string configname = VERSION + ".config";
	ifstream infile( configname.c_str() );
	if(!infile.is_open())
		return false;
	
	string line;
	int line_num=1;
	while (getline(infile, line))
	{
		istringstream ss(line);
		
		//string temp;
		if(line_num ==1)
			ss >> master_cnv_file;
		if(line_num ==4)
			ss >> cumulative_database_file;	
		if(line_num ==7)
		{
			string s = line;
			string delimiter = "|";
			size_t pos = 0;
			string token;
			while ((pos = s.find(delimiter)) != string::npos)
			{
				token = s.substr(0, pos);
				similar_phenotypes.push_back(token);
				//cout << token << endl;
				s.erase(0, pos + delimiter.length());
			}
			similar_phenotypes.push_back(s);	//push back the final group of phenotypes
		}
	
		if(line_num ==11)
			ss >> out_file;
		line_num++;
	}
		
		//cout << "I got to the while loop for file reading\n";
	infile.close();
	return true;
}



int main()
{
	vector<CNV> master_list_cnvs, comparison_list_cnvs;	//main vectors of CNVs for comparing 

	string filename;
	cout << "***********************\n***********************\n"
	     << VERSION << endl
		 << "\n***********************\n***********************\n";
	cout << "\nLists to compare should be tab-delimited files in the following format\n";
	cout << "ourID\tstudyID\tphenotype\tcase/control(1,2)\tchr\tstart\tstop\tcopy_num\n";
	cout << "\nKeep headers in your input files! WARNING: For now remove commas from numbers. ie dont use 1,000,000, use 1000000 in input files\n";
	cout << "I will try to deal with commas in numbers in later versions, but for now leave them out\n";
	cout << "Press enter to continue\n";
	cin.get();
	
	cout << "Reading config file...\n";
	if(!DefineConfig())	cerr << "Unable to read the config file!\n";
	
	cout << "Reading CNV file...\n";
	if(!ReadFile(master_cnv_file,master_list_cnvs))
	{
		cerr << "Looks like there was a problem reading " << master_cnv_file;
		cerr << "Make sure you updated the compare_lists.config file with the two files you want to compare\n";	
		cerr << "\nPerhaps you didn't specify the correct path?\n";
		return 1;
	}

	cout << "Reading Cumulative Database file...\n";
	if(!ReadFile(cumulative_database_file,comparison_list_cnvs))
	{
		cerr << "Looks like there was a problem reading " << cumulative_database_file;
		cerr << "Make sure you updated the compare_lists.config file with the two files you want to compare\n";	
		cerr << "\nPerhaps you didn't specify the correct path?\n";
		return 1;
	}
		
	int num_cnvs = master_list_cnvs.size(), num_cumulative = comparison_list_cnvs.size();	
	
	cout << "Running now...\n";
	CompareCnvs(master_list_cnvs,comparison_list_cnvs);
	
	cout << "Writing log file...\n";
	if(!WriteLog(num_cnvs, num_cumulative))
		cerr << "Couldn't open log file for writing\n";

	if(!WriteCnvs(master_list_cnvs))
		cerr << "Couldn't open output file for writing\n";
	else
		cout << "All set. Check you're output file: " << out_file << endl;
	
	return 0;
}
	
	