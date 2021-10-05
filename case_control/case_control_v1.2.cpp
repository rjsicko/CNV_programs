/*****************************************************
case_control_v1.2
Author: Robert Sicko - R.J.Sicko@gmail.com

****DESCRIPTION****
	takes a list of CNVs and outputs the list + appends the following extra columns:
			% Overlap Control(same CNV type)
			% Overlap Control(diff CNV type)
			Num Cases Overlapped (same phenotype, same CNV type)[Threshold: olap_threshold]
			Num Cases Overlapped (similar phenotype, same CNV type)[Threshold: olap_threshold]
			Num Cases Overlapped (different phenotype, same CNV type)[Threshold: olap_threshold]
			Num Cases Overlapped (same phenotype, diff CNV type)[Threshold: olap_threshold]
			Num Cases Overlapped (similar phenotype, diff CNV type)[Threshold: olap_threshold]
			Num Cases Overlapped (different phenotype, diff CNV type)[Threshold: olap_threshold]
****END DESCRIPTION****

****REQUIRED FORMAT****
	NOTE: No commas in numbers (ie. use 1000000 not 1,000,000)
	NOTE2: Keep headers in your files. The program ignores the first line (don't have spaces in header id's though).
	
	File should be tab-delimited in the following format:
	ourID	sample	phenotype	case/control	chr	start	stop	copy_num
****END REQUIRED FORMAT****

****CHANGE LOG****
v1.0- Bugs fixed, programs appears to be functioning as expected
v1.1- Added ability to have multiple groups of similar phenotypes
v1.2- Added a threshold for overlap to count cases that overlap

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
#include <iomanip>


using namespace std;

//A struct for our CNVs (or regions)
struct CNV{
	//ourID	sample	phenotype	case/control	chr	start	stop	copy_num
	string our_id, sample_id, phenotype, chr;	//our identifier
	int copy_num, case_or_control;	
	long start, stop;	//start and end points for each cnv
	int num_cases_overlap_with_same, num_other_similar_disorders_same, num_other_diff_disorders_same;
	int num_cases_overlap_with_diff, num_other_similar_disorders_diff, num_other_diff_disorders_diff;
	double percent_control_overlap_same, percent_control_overlap_diff;
	
};

//Mostly for debugging
//outdated!
/*void ConsolePrintCnvs(vector<CNV> any_list_of_cnvs)
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
		
}*/

bool ReadFile(string filename,vector<CNV>& our_cnvs)
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
		//ourID	sample	phenotype	case/control	chr	start	stop	copy_num
		ss >> dummy_cnv.our_id >> dummy_cnv.sample_id
			>> dummy_cnv.phenotype >> dummy_cnv.case_or_control
			>> dummy_cnv.chr >> dummy_cnv.start
			>> dummy_cnv.stop >> dummy_cnv.copy_num;
		
		//first_line = false;
		dummy_cnv.num_cases_overlap_with_same = 0;
		dummy_cnv.num_other_similar_disorders_same = 0;
		dummy_cnv.num_other_diff_disorders_same = 0;
		dummy_cnv.percent_control_overlap_same = 0;
				
		dummy_cnv.num_cases_overlap_with_diff = 0;
		dummy_cnv.num_other_similar_disorders_diff = 0;
		dummy_cnv.num_other_diff_disorders_diff = 0;
		dummy_cnv.percent_control_overlap_diff = 0;
		//cout << dummy_cnv.our_id << "\t" << dummy_cnv.chr << "\t" << dummy_cnv.start << endl;
		//cout << dummy_cnv.phenotype << "\t" << dummy_cnv.case_or_control;
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
void CompareCnvs(vector<CNV>& master_list_cnvs)
{
			/*string ls = line;
			string delimiter = "\t";
			string token;
			int pos = 0;
			while ((pos = ls.find(delimiter)) != std::string::npos)
			{
				token = ls.substr(0, pos);
				ls.erase(0, pos + delimiter.length());
			}
			
			
			string token = ls.substr(0, s.find(delimiter));*/
		
	//cout << type_specific << endl;
	for (vector<CNV>::iterator me = master_list_cnvs.begin();me != master_list_cnvs.end();me++)
		for(vector<CNV>::iterator them = master_list_cnvs.begin();them != master_list_cnvs.end();them++)
		{
			if(me->our_id == them->our_id)	//lets not count ourself and skip if we are looking at ourself
				continue;
			if(me->chr==them->chr)
			{
				double me_start = me->start, me_stop = me->stop, them_start = them->start, them_stop = them->stop;	
				double overlap = CalcOverlapCNV(me_start, them_start, me_stop, them_stop);								
				if(overlap > 0)
				{
					//cout << overlap << endl;
					bool type_specific_pass = true;
					if(me->copy_num==them->copy_num
						||me->copy_num==(them->copy_num+1)
						||me->copy_num==(them->copy_num-1))
					{
						type_specific_pass = true;

					}
					else
						type_specific_pass = false;

					if(type_specific_pass)
					{
						if(them->case_or_control == 2)	//they are a control
						{
							//percent_control_overlap_same, percent_control_overlap_diff
							if(me->percent_control_overlap_same < overlap)
							{
								//cout << overlap << endl;
								me->percent_control_overlap_same = overlap;	//new champ
							}
						}
						else if(me->phenotype == them->phenotype)
						{//same phenotype as us
							//	int num_cases_overlap_with_same, num_other_similar_disorders_same;
							//	int num_cases_overlap_with_diff, num_other_similar_disorders_diff;
							if(overlap > olap_threshold)
								me->num_cases_overlap_with_same++;
							
						}
						else if(overlap > olap_threshold)
						//either a similar disorder or a different disorder
						//and the threshold is met
						{
							bool similar_pheno = false;	//flag to figure out if we compared two CNVs that are of a similar CNV type
							for(vector<string>::iterator itr = similar_phenotypes.begin();
								itr != similar_phenotypes.end();itr++)
							{
								size_t found = itr->find(me->phenotype.c_str());
								if(found != string::npos)	//if we are in one of the strings in the similar phenotypes vector
								{
									found = itr->find(them->phenotype.c_str());
									if(found != string::npos)	//if them phenotype was in the same similar phenotype string as "me" was in
										similar_pheno = true;
								}
							}
							if(similar_pheno == true)
								me->num_other_similar_disorders_same++;
							else
								me->num_other_diff_disorders_same++;
						}
					}
					else	// diff type
					{
						if(them->case_or_control == 2)
						{
							//percent_control_overlap_same, percent_control_overlap_diff
							if(me->percent_control_overlap_diff < overlap)
								me->percent_control_overlap_diff = overlap;
						}
						else if(me->phenotype == them->phenotype)
						{
							//	int num_cases_overlap_with_same, num_other_similar_disorders_same;
							//	int num_cases_overlap_with_diff, num_other_similar_disorders_diff;
							if(overlap > olap_threshold)
								me->num_cases_overlap_with_diff++;
						}
						else if(overlap > olap_threshold)
						//either a similar disorder or a different disorder
						//and the threshold is met
						{
							bool similar_pheno = false;	//flag to figure out if we compared two CNVs that are of a similar CNV type
							for(vector<string>::iterator itr = similar_phenotypes.begin();
								itr != similar_phenotypes.end();itr++)
							{
								size_t found = itr->find(me->phenotype.c_str());
								if(found != string::npos)	//if we are in one of the strings in the similar phenotypes vector
								{
									found = itr->find(them->phenotype.c_str());
									if(found != string::npos)	//if them phenotype was in the same similar phenotype string as "me" was in
										similar_pheno = true;
								}
							}
							if(similar_pheno == true)
								me->num_other_similar_disorders_diff++;
							else
								me->num_other_diff_disorders_diff++;
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
	////ourID	sample	phenotype	case/control	chr	start	stop	copy_num
	outfile << "OurID\tSampleID\tPhenotype\tCase/Control\tChr\tStart\tStop\tCopy_Number\t"
			<< "% Overlap Control(same CNV type)\t% Overlap Control(diff CNV type)\t"
			<< "Num Cases Overlapped (same phenotype, same CNV type)[Threshold:" << olap_threshold << "]"
			<< "\tNum Cases Overlapped (similar phenotype, same CNV type)[Threshold:" << olap_threshold << "]"
			<< "\tNum Cases Overlapped (different phenotype, same CNV type)[Threshold:" << olap_threshold << "]"
			<< "\tNum Cases Overlapped (same phenotype, diff CNV type)[Threshold:" << olap_threshold << "]"
			<< "\tNum Cases Overlapped (similar phenotype, diff CNV type)[Threshold:" << olap_threshold << "]"
			<< "\tNum Cases Overlapped (different phenotype, diff CNV type)[Threshold:" << olap_threshold << "]"
			<< endl;
	for (vector<CNV>::iterator me = master_list_cnvs.begin();me != master_list_cnvs.end();me++)
	{
		double olap_percent_same = me->percent_control_overlap_same*100;
		ostringstream strs;
		strs << fixed << setprecision (2) << olap_percent_same;
		string olap_same = strs.str();	
		if(olap_same == "100.00") olap_same = "100";
		
		double olap_percent_diff = me->percent_control_overlap_diff*100;
		ostringstream strs_diff;
		strs_diff << fixed << setprecision (2) << olap_percent_diff;
		string olap_diff = strs_diff.str();	
		if(olap_diff == "100.00") olap_diff = "100";
		
		outfile << me->our_id << "\t" << me->sample_id << "\t" 
				<< me->phenotype << "\t" << me->case_or_control << "\t" 
				<< me->chr << "\t" << me->start << "\t" << me->stop
				<< "\t" << me->copy_num << "\t" 
				<< olap_same << "\t"
				<< olap_diff << "\t"
				<< me->num_cases_overlap_with_same << "\t"
				<< me->num_other_similar_disorders_same << "\t"
				<< me->num_other_diff_disorders_same << "\t"
				<< me->num_cases_overlap_with_diff << "\t"
				<< me->num_other_similar_disorders_diff << "\t"
				<< me->num_other_diff_disorders_diff << "\n";
	}
	
	/*	int num_cases_overlap_with_same, num_other_similar_disorders_same, num_other_diff_disorders_same;
	int num_cases_overlap_with_diff, num_other_similar_disorders_diff, num_other_diff_disorders_diff;
	double percent_control_overlap_same, percent_control_overlap_diff;*/
	
	outfile.close();
	return true;
}

bool WriteLog(int num_cnvs)
{
	string out_trans = out_file + ".Log";
	ofstream outfile(out_trans.c_str());
	if(!outfile){
		 cerr << "Can not open the log file for writing!\nIt may be in use. Close it and try again\n";
		 return false;
	}
	outfile << "Program: case_control_v1.2\n";
	outfile << "Number CNVs read: " << num_cnvs << endl;
	outfile << "CNV file used: " << master_cnv_file << endl;
	outfile << "Phenotypes called similar: ";
	for(vector<string>::iterator itr = similar_phenotypes.begin();
		itr != similar_phenotypes.end();itr++)
	{
		outfile << *itr;
		outfile << "\n";
	}
	outfile << "\nThreshold: " << olap_threshold << "Cases with overlap > than this counted" << endl;

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
	ifstream infile("case_control_v1.2.config");
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
		if(line_num ==6)
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
		if(line_num ==10)
			ss >> olap_threshold;
		if(line_num ==14)
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
	cout << "***********************\n***********************\nCase_Control v1.2\n***********************\n***********************\n";
	cout << "\nLists to compare should be tab-delimited files in the following format\n";
	cout << "OurID\tChr\tStart\tStop\tCase/Control\tCopy_Number\n";
	cout << "\nKeep headers in your input files! WARNING: For now remove commas from numbers. ie dont use 1,000,000, use 1000000 in input files\n";
	cout << "I will try to deal with commas in numbers in later versions, but for now leave them out\n";
	cout << "Press enter to continue\n";
	cin.get();

	//cout << "Specify first file to read (Including path)\n";
	//cin >> filename;
	
	if(!DefineConfig())	cerr << "Unable to read the config file!\n";
	if(ReadFile(master_cnv_file,master_list_cnvs)){	}
	else{
		cerr << "Looks like there was a problem reading " << master_cnv_file;
		cerr << "Make sure you updated the compare_lists.config file with the two files you want to compare\n";
		cerr << "Also, the cnv files need to be in the same directory as the program\n";		
		//cerr << "\nPerhaps you didn't specify the path?\n";
		return 1;
	}
	//cout << "your overlap threshold is " << percent_overlap_cutoff << endl;
	
	
	cout << "Running compare_cnvs now\n";
	CompareCnvs(master_list_cnvs);
	int num_cnvs = master_list_cnvs.size();
	if(!WriteLog(num_cnvs))
		cerr << "Couldn't open log file for writing\n";
	if(WriteCnvs(master_list_cnvs))
		cout << "All set. Check you're output file: " << out_file << endl;
	//ConsolePrintCnvs(master_list_cnvs);
	return 0;
}
	
	