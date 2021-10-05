/*****************************************************
case_control_v1.2_sexchr
Author: Robert Sicko - R.J.Sicko@gmail.com

****DESCRIPTION****
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
	
	used case_control program_v1.2, but modified it to handle sex chromosome calls 
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
****END DESCRIPTION****

****REQUIRED FORMAT****
	NOTE: No commas in numbers (ie. use 1000000 not 1,000,000)
	NOTE2: Keep headers in your files. The program ignores the first line (don't have spaces in header id's though).
	
	File should be tab-delimited in the following format:
	ourID	sample	phenotype	case/control	chr	start	stop	copy_num	sex
****END REQUIRED FORMAT****

****CHANGE LOG****
v1.0- Bugs fixed, programs appears to be functioning as expected
v1.1- Added ability to have multiple groups of similar phenotypes
v1.2- Added a threshold for overlap to count cases that overlap
v1.2_sexchr- modified program for sex chromosomes...

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
	//ourID	sample	phenotype	case/control	chr	start	stop	copy_num	sex
	string our_id, sample_id, phenotype, chr, sex;	//our identifier
	int copy_num, case_or_control;	
	long start, stop;	//start and end points for each cnv
	int num_cases_overlap_with_exact, num_other_similar_disorders_exact, num_other_diff_disorders_exact;
	int num_cases_overlap_with_similar, num_other_similar_disorders_similar, num_other_diff_disorders_similar;
	double percent_control_overlap_exact, percent_control_overlap_similar;
	string description, internal_type;
};

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
		//handle blank line
		if(line == "") continue;
		istringstream ss(line);
		//ourID	sample	phenotype	case/control	chr	start	stop	copy_num
		ss >> dummy_cnv.our_id >> dummy_cnv.sample_id
			>> dummy_cnv.phenotype >> dummy_cnv.case_or_control
			>> dummy_cnv.chr >> dummy_cnv.start
			>> dummy_cnv.stop >> dummy_cnv.copy_num >> dummy_cnv.sex;
		
		//first_line = false;
		dummy_cnv.num_cases_overlap_with_exact = 0;
		dummy_cnv.num_other_similar_disorders_exact = 0;
		dummy_cnv.num_other_diff_disorders_exact = 0;
		dummy_cnv.percent_control_overlap_exact = 0;

		dummy_cnv.num_cases_overlap_with_similar = 0;
		dummy_cnv.num_other_similar_disorders_similar = 0;
		dummy_cnv.num_other_diff_disorders_similar = 0;
		dummy_cnv.percent_control_overlap_similar = 0;
		
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
	/*
	int num_cases_overlap_with_exact, num_other_similar_disorders_exact, num_other_diff_disorders_exact;
	int num_cases_overlap_with_similar, num_other_similar_disorders_similar, num_other_diff_disorders_similar;
	double percent_control_overlap_exact, percent_control_overlap_similar;
	*/
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
					if(me->sex == them->sex && me->copy_num == them->copy_num)	//exact type
					{
						if(them->case_or_control == 2)	//they are a control
						{
							if(me->percent_control_overlap_exact < overlap)
							{
								//cout << overlap << endl;
								me->percent_control_overlap_exact = overlap;	//new champ
							}
						}
						else if(me->phenotype == them->phenotype)
						{//same phenotype as us
							//	int num_cases_overlap_with_exact, num_other_similar_disorders_exact;
							//	int num_cases_overlap_with_diff, num_other_similar_disorders_diff;
							if(overlap > olap_threshold)
								me->num_cases_overlap_with_exact++;
							
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
								me->num_other_similar_disorders_exact++;
							else
								me->num_other_diff_disorders_exact++;
						}
					}
					else if(me->internal_type == them->internal_type)	//similar type
					{
						if(them->case_or_control == 2)	//they are a control
						{
							if(me->percent_control_overlap_similar < overlap)
							{
								//cout << overlap << endl;
								me->percent_control_overlap_similar = overlap;	//new champ
							}
						}
						else if(me->phenotype == them->phenotype)
						{//same phenotype as us
							//	int num_cases_overlap_with_exact, num_other_similar_disorders_exact;
							//	int num_cases_overlap_with_diff, num_other_similar_disorders_diff;
							if(overlap > olap_threshold)
								me->num_cases_overlap_with_similar++;
							
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
								me->num_other_similar_disorders_similar++;
							else
								me->num_other_diff_disorders_similar++;
						}
					}
				}
			}
		}
}

void GenerateType(vector<CNV>& master_list_cnvs)
{
	//CIDR notes with X project list PAR and XTR coordinates that are mis-annotated in some illumina products
	//These coordinates are used here
	
	for (vector<CNV>::iterator me = master_list_cnvs.begin();me != master_list_cnvs.end();me++)
	{
		//illumina PAR regions
		if(me->chr == "XY")
		{
			if(me->copy_num == 0)
			{
				me->description = "Loss of both copies (XY region)";
				me->internal_type = "loss";
			}
			if(me->copy_num == 1)
			{
				me->description = "Loss of one copy (XY region)";
				me->internal_type = "loss";
			}
			if(me->copy_num == 2)
			{
				me->description = "LOH (XY region)";
				me->internal_type = "normal";
			}
			if(me->copy_num == 3)
			{
				me->description = "Gain of one copy (XY region)";
				me->internal_type = "gain";
			}
			if(me->copy_num == 4)
			{
				me->description = "Gain of two copies (XY region)";
				me->internal_type = "gain";
			}
		}
		else if(me->sex == "M")
		{
			if(me->chr == "X")
			{
				if(me->copy_num == 0)
				{
					me->description = "Loss of one(only) copy";
					me->internal_type = "loss";
				}
				if(me->copy_num == 1)
				{
					me->description = "Normal";
					me->internal_type = "normal";
				}
				if(me->copy_num == 2)
				{
					me->description = "Gain of one copy";
					me->internal_type = "gain";
				}
				if(me->copy_num == 3)
				{
					me->description = "Gain of two copies";
					me->internal_type = "gain";
				}
				if(me->copy_num == 4)
				{
					me->description = "Gain of three copies";
					me->internal_type = "gain";
				}
			}
			if(me->chr == "Y")
			{
				if(me->copy_num == 0)
				{
					me->description = "Loss of one(only) copy";
					me->internal_type = "loss";
				}
				if(me->copy_num == 1)
				{
					me->description = "Normal";
					me->internal_type = "normal";
				}
				if(me->copy_num == 2)
				{
					me->description = "Gain of one copy";
					me->internal_type = "gain";
				}
				if(me->copy_num == 3)
				{
					me->description = "Gain of two copies";
					me->internal_type = "gain";
				}
				if(me->copy_num == 4)
				{
					me->description = "Gain of three copies";
					me->internal_type = "gain";
				}	
			}
		}	
		else if(me->sex == "F")
		{
			if(me->chr == "X")
			{
				if(me->copy_num == 0)
				{
					me->description = "Loss of both copies";
					me->internal_type = "loss";
				}
				if(me->copy_num == 1)
				{
					me->description = "Loss of one copy";
					me->internal_type = "loss";
				}
				if(me->copy_num == 2)
				{
					me->description = "LOH";
					me->internal_type = "normal";
				}
				if(me->copy_num == 3)
				{
					me->description = "Gain of one copy";
					me->internal_type = "gain";
				}
				if(me->copy_num == 4)
				{
					me->description = "Gain of two copies";	
					me->internal_type = "gain";
				}
			}
			if(me->chr == "Y")
			{
				if(me->copy_num == 0)
				{	
					me->description = "Normal";
					me->internal_type = "normal";
				}
				if(me->copy_num == 1)
				{
					me->description = "Gain of one copy";
					me->internal_type = "gain";
				}
				if(me->copy_num == 2)
				{	
					me->description = "Gain of two copies";
					me->internal_type = "gain";
				}
				if(me->copy_num == 3)
				{
					me->description = "Gain of three copies";
					me->internal_type = "gain";
				}
				if(me->copy_num == 4)
				{
					me->description = "Gain of four copies";
					me->internal_type = "gain";
				}				
			}
		}
		//now annotate the CNVs not annotated PAR by illumina
		//now PAR/XTR check for illumina mis-annotations
		if(me->chr == "X")
		{
			//PAR1
			if( CalcOverlapCNV(me->start, 60001, me->stop, 2699520) > 0 )
			{
				if(me->copy_num == 0)
				{
					me->description = "Loss of both copies (PAR1)";
					me->internal_type = "loss";
				}
				if(me->copy_num == 1)
				{
					me->description = "Loss of one copy (PAR1)";
					me->internal_type = "loss";
				}
				if(me->copy_num == 2)
				{
					me->description = "Normal (PAR1)";
					me->internal_type = "normal";
				}
				if(me->copy_num == 3)
				{
					me->description = "Gain of one copy (PAR1)";
					me->internal_type = "gain";
				}
				if(me->copy_num == 4)
				{
					me->description = "Gain of two copies (PAR1)";
					me->internal_type = "gain";
				}
			}
			//PAR2
			else if( CalcOverlapCNV(me->start, 154931044, me->stop, 155260560) > 0  )
			{
				if(me->copy_num == 0)
				{
					me->description = "Loss of both copies (PAR2)";
					me->internal_type = "loss";
				}
				if(me->copy_num == 1)
				{
					me->description = "Loss of one copy (PAR2)";
					me->internal_type = "loss";
				}
				if(me->copy_num == 2)
				{
					me->description = "Normal (PAR2)";
					me->internal_type = "normal";
				}
				if(me->copy_num == 3)
				{
					me->description = "Gain of one copy (PAR2)";
					me->internal_type = "gain";
				}
				if(me->copy_num == 4)
				{
					me->description = "Gain of two copies (PAR2)";
					me->internal_type = "gain";
				}
			}
			//XTR
			else if( CalcOverlapCNV(me->start, 88395830, me->stop, 92583067) > 0  )
			{
				if(me->copy_num == 0)
				{
					me->description = "Loss of both copies (XTR)";
					me->internal_type = "loss";
				}
				if(me->copy_num == 1)
				{
					me->description = "Loss of one copy (XTR)";
					me->internal_type = "loss";
				}
				if(me->copy_num == 2)
				{
					me->description = "Normal (XTR)";
					me->internal_type = "normal";
				}
				if(me->copy_num == 3)
				{
					me->description = "Gain of one copy (XTR)";
					me->internal_type = "gain";
				}
				if(me->copy_num == 4)
				{
					me->description = "Gain of two copies (XTR)";
					me->internal_type = "gain";
				}
			}
		}
		if(me->chr == "Y")
		{
			//PAR1
			if( CalcOverlapCNV(me->start, 10001, me->stop, 2649520) > 0  )
			{
				if(me->copy_num == 0)
				{
					me->description = "Loss of both copies (PAR1)";
					me->internal_type = "loss";
				}
				if(me->copy_num == 1)
				{
					me->description = "Loss of one copy (PAR1)";
					me->internal_type = "loss";
				}
				if(me->copy_num == 2)
				{
					me->description = "Normal (PAR1)";
					me->internal_type = "normal";
				}
				if(me->copy_num == 3)
				{
					me->description = "Gain of one copy (PAR1)";
					me->internal_type = "gain";
				}
				if(me->copy_num == 4)
				{
					me->description = "Gain of two copies (PAR1)";
					me->internal_type = "gain";
				}
			}
			//PAR2
			else if( CalcOverlapCNV(me->start, 59034050, me->stop, 59363566) > 0  )
			{
				if(me->copy_num == 0)
				{
					me->description = "Loss of both copies (PAR2)";
					me->internal_type = "loss";
				}
				if(me->copy_num == 1)
				{
					me->description = "Loss of one copy (PAR2)";
					me->internal_type = "loss";
				}
				if(me->copy_num == 2)
				{
					me->description = "Normal (PAR2)";
					me->internal_type = "normal";
				}
				if(me->copy_num == 3)
				{
					me->description = "Gain of one copy (PAR2)";
					me->internal_type = "gain";
				}
				if(me->copy_num == 4)
				{
					me->description = "Gain of two copies (PAR2)";
					me->internal_type = "gain";
				}
			}
			//XTR
			else if( CalcOverlapCNV(me->start, 2926383, me->stop, 6543100) > 0 )
			{
				if(me->copy_num == 0)
				{
					me->description = "Loss of both copies (XTR)";
					me->internal_type = "loss";
				}
				if(me->copy_num == 1)
				{
					me->description = "Loss of one copy (XTR)";
					me->internal_type = "loss";
				}
				if(me->copy_num == 2)
				{
					me->description = "Normal (XTR)";
					me->internal_type = "normal";
				}
				if(me->copy_num == 3)
				{
					me->description = "Gain of one copy (XTR)";
					me->internal_type = "gain";
				}
				if(me->copy_num == 4)
				{
					me->description = "Gain of two copies (XTR)";
					me->internal_type = "gain";
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
	////ourID	sample	phenotype	case/control	chr	start	stop	copy_num	sex
	outfile << "OurID\tSampleID\tPhenotype\tCase/Control\tChr\tStart\tStop\tCopy_Number\tSex\tDescription\t"
			<< "% Overlap Control(exact CNV type)\t% Overlap Control(similar CNV type)\t"
			<< "Num Cases Overlapped (same phenotype, exact CNV type)[Threshold:" << olap_threshold << "]"
			<< "\tNum Cases Overlapped (similar phenotype, exact CNV type)[Threshold:" << olap_threshold << "]"
			<< "\tNum Cases Overlapped (different phenotype, exact CNV type)[Threshold:" << olap_threshold << "]"
			<< "\tNum Cases Overlapped (same phenotype, similar CNV type)[Threshold:" << olap_threshold << "]"
			<< "\tNum Cases Overlapped (similar phenotype, similar CNV type)[Threshold:" << olap_threshold << "]"
			<< "\tNum Cases Overlapped (different phenotype, similar CNV type)[Threshold:" << olap_threshold << "]"
			<< endl;
	for (vector<CNV>::iterator me = master_list_cnvs.begin();me != master_list_cnvs.end();me++)
	{
		double olap_percent_exact = me->percent_control_overlap_exact*100;
		ostringstream strs;
		strs << fixed << setprecision (2) << olap_percent_exact;
		string olap_exact = strs.str();	
		if(olap_exact == "100.00") olap_exact = "100";
		
		double olap_percent_sim = me->percent_control_overlap_similar*100;
		ostringstream strs_sim;
		strs_sim << fixed << setprecision (2) << olap_percent_sim;
		string olap_sim = strs_sim.str();	
		if(olap_sim == "100.00") olap_sim = "100";
		
		outfile << me->our_id << "\t" << me->sample_id << "\t" 
				<< me->phenotype << "\t" << me->case_or_control << "\t" 
				<< me->chr << "\t" << me->start << "\t" << me->stop
				<< "\t" << me->copy_num << "\t" 
				<< me->sex << "\t" << me->description << "\t"
				<< olap_exact << "\t"
				<< olap_sim << "\t"
				<< me->num_cases_overlap_with_exact << "\t"
				<< me->num_other_similar_disorders_exact << "\t"
				<< me->num_other_diff_disorders_exact << "\t"
				<< me->num_cases_overlap_with_similar << "\t"
				<< me->num_other_similar_disorders_similar << "\t"
				<< me->num_other_diff_disorders_similar << "\n";
	}
	
	/*	int num_cases_overlap_with_exact, num_other_similar_disorders_exact, num_other_diff_disorders_exact;
	int num_cases_overlap_with_similar, num_other_similar_disorders_similar, num_other_diff_disorders_similar;
	double percent_control_overlap_exact, percent_control_overlap_similar;
	string description;*/
	
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
	ifstream infile("case_control_v1.2_sexchr.config");
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
	cout << "***********************\n***********************\nCase_Control v1.2_sexchr\n***********************\n***********************\n";
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
	GenerateType(master_list_cnvs);
	if(WriteCnvs(master_list_cnvs))
		cout << "All set. Check you're output file: " << out_file << endl;
	//ConsolePrintCnvs(master_list_cnvs);
	return 0;
}
	
	