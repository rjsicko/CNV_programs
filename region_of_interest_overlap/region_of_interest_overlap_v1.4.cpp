/*****************************************************
region_of_interest_overlap v1.4
Author: Robert Sicko - R.J.Sicko@gmail.com

****DESCRIPTION****
This program compares a list of CNVs to a list of regions of interest
and outputs the original list with two additional columns appended.
The two columns appended are 'Overlapped_Regions' and 'Percent_Overlap'
if a CNV overlaps more than one region, the regions will be seperated
by a ; as will the percent overlaps for each
****END DESCRIPTION****

****REQUIRED FORMAT****
	NOTE: No commas in numbers (ie. use 1000000 not 1,000,000)
	NOTE2: Keep headers in your files. The program ignores the first line (don't have spaces in header id's though).
	
	CNV List (tab-delimited txt):
	OurID	Chr		Start	Stop	Case/Control	Copy_Number
	
	Region List(tab-delimited txt):
	Database_From	Unique_ID_1(gene)	Unique_ID_2(transcript)		Chr		Start	Stop
****END REQUIRED FORMAT****

****CHANGE LOG****
v1.0- Bugs fixed, programs appears to be functioning as expected
v1.1- Added "NA" to 'Overlapped_Regions' and 'Percent_Overlap' fields if no regions were detected for a particular CNV
	- Added option in 'region_of_interest_overlap.config' to choose the type of overlap used
		Note: see 'Explanation_Overlap_Options.txt' for more information on the types of overlap
	- Commented a lot more of the code (I'll try to comment it all eventually)
	- Fixed a the 'Percent_Overlap' output so now allows 2 decimels
v1.2- Added option for writing 'transcript_id' field in output file
		As a consequence, the input file format now contains an extra 'transcript_id' field
		If you do not need this field, simply use a dummy code in this column and in config file specify 0 for write_transcript_id
		NOTE: Did not implement the do not write_transcript_id option yet. Therefore, write_transcript_id has to be used.
			  Will make option to not write transcript_id in v1.3!
	- Since I was already modifying the config file for the above option, I cleaned it up and removed entries not used
	- As a consequence of adding the write transcript_id option, the storage format has changed for gene_list and overlap_percent.
	  Each CNV now has a multimap<string,string> & a map<string,string> where the multimap is <gene_name,transcript_ids>
	  and the map is <transcript_id,overlap_percentage>.
	  Since the storage was changed to maps, the output is now alphabetically sorted.
	- "Bug" fixed when importing program output into Excel. In previous versions output was surrounded by (). For single transcript genes,
	  when the data was imported into Excel, the percent_overlap would be a negative number. This is because Excel automatically converts
	  (#) into -#. This can be changed when importing the data by telling Excel to treat the field as text. However, the easiest solution 
	  was to change the surrounding ()'s into []'s.
	- Added column for number_regions_overlapped
v1.3- Added 3 output files. A gene list file , a transcript output file and a log file.
	- The gene list file is a list of genes overlapped and the number cases and number controls that overlap it
	- The transcript list file is a list of genes overlapped and the number cases and number controls that overlap it
	- The log file is a file describing the particular run of the program
	- Added check for overlap_cutoff > 1; just incase user specifies 90 instead of 0.90 in config file!
	- Added check for config file manipulation.
		- count total line numbers and make sure user didnt remove any lines
v1.4- Changed chr variable from int to string
	- Added const string VERSION variable for tracking what version we are
****END CHANGE LOG****

****TO DO****
	- Add other checks for config file parameters outside expected range.
	- Comment out code!
	- Implement option to not write transcript.
	- Add check for output being too long for a single excel cell. In that case, use secondary column for output.
	- Fix case/control count so it doesn't just count the number of CNVs that overlap a gene, instead counts individuals
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

const string VERSION = "region_of_interest_overlap_v1.4";

//A struct for our CNVs (or regions)
struct CNV{
	string our_id, chr;	//our identifier
	int copy_num, case_control;	//ints for chromosome, copy numbe and case or control (1, 2)
	long start, stop;	//start and end points for each cnv
	//int num_cases_overlap_with, num_controls_overlap_with; // count variables, +1 each time we overlap with a case or control
	//string regions_overlapped, percent_overlapped;
	multimap<string,string> genes_overlapped;	//this multimap is <gene,transcript>
	map<string,string> percent_overlapped;	//this map is <transcript,percent_overlap>
	
};

struct REGION{

	//"Database_From\tUnique_ID\tChr\tStart\tStop\n";
	string database_id, unique_id, unique_id_2;	//our identifier
	string chr;	//ints for chromosome, copy numbe and case or control (1, 2)
	long start, stop;	//start and end points for each cnv
	//int num_cases_overlap_with, num_controls_overlap_with; // count variables, +1 each time we overlap with a case or control
	
};

struct REGION_COUNT{
	int num_cases, num_controls;
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

void ConsolePrint(vector<REGION> any_regions)
{
	//cout << "I got to ConsolePrintCnvs function\nAnd my cnvs list contains: ";
	//cout << any_list_of_cnvs.size();
	//cout << flush;
	cout << "Database_ID\tUnique_ID\tChr\tStart\tStop\n";
	for (vector<REGION>::iterator me = any_regions.begin();me != any_regions.end();me++)
			cout << me->database_id<< "\t" << me->unique_id << "\t" 
				<< me->chr << "\t" << me->start << "\t" 
				<< me->stop << "\t" << endl;
		
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
		ss >> dummy_cnv.our_id >> dummy_cnv.chr
			>> dummy_cnv.start >> dummy_cnv.stop
			>> dummy_cnv.case_control >> dummy_cnv.copy_num;
		//ss >> name >> var1 >> var2 >> var3;
		//first_line = false;
		//dummy_cnv.regions_overlapped = "";
		//dummy_cnv.percent_overlapped = "";
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

bool ReadFileRegion(string filename,vector<REGION>& regions)
{
	ifstream infile(filename.c_str());
	if(!infile.is_open())
		return false;
	REGION dummy_region;
	bool first_line=true;
	string line;
	while (getline(infile, line))
	{
			
		istringstream ss(line);
		//string name;
		//int var1, var2, var3;
		ss >> dummy_region.database_id >> dummy_region.unique_id >> dummy_region.unique_id_2 
			>> dummy_region.chr	>> dummy_region.start >> dummy_region.stop;

		if(first_line)
			first_line=false;
		else
			regions.push_back (dummy_region);
		
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

	

//function to compare the master list with a test or comparison list
//modifies master_list_cnvs 'num_cases_overlap_with' variable to keep track of the number of overlaps each entry in master list has with the comparison list 
void CheckRegions(vector<CNV>& master_list_cnvs, vector<REGION> regions)
{
	//cout << type_specific << endl;
	for (vector<CNV>::iterator me = master_list_cnvs.begin();me != master_list_cnvs.end();me++)
	{
		//vector<REGION>::iterator previous_region;
		//bool first_run=true;
		for(vector<REGION>::iterator them = regions.begin();them != regions.end();them++)
		{
			if(me->chr==them->chr)
			{
				double me_start = me->start, me_stop = me->stop, them_start = them->start, them_stop = them->stop;	
				double overlap;
				if(overlap_type == 1)	overlap = CalcOverlapUnion(me_start,them_start,me_stop,them_stop);
				if(overlap_type == 2)	overlap = CalcOverlapRegion(me_start,them_start,me_stop,them_stop);
				if(overlap_type == 3)	overlap = CalcOverlapCNV(me_start,them_start,me_stop,them_stop);
		
				if(overlap>percent_overlap_cutoff)
				{
											/*if(me->regions_overlapped != ""){
												me->regions_overlapped = me->regions_overlapped + ";";
												me->percent_overlapped = me->percent_overlapped + ";";
											}
											double olap_percent = overlap*100;
											ostringstream strs;
											strs << setprecision (4) << olap_percent;
											string olap = strs.str();
											if(write_transcript_id==0)
												me->regions_overlapped = me->regions_overlapped + them->unique_id;
											if(write_transcript_id==1){
												if(first_run){
													me->regions_overlapped = me->regions_overlapped + them->unique_id + "(" + them->unique_id_2;
													me->percent_overlapped = me->percent_overlapped + "(" + olap;
												}
												else{
													first_run=false;
													if(them->unique_id==previous_region->unique_id){	//if we match previous gene
														me->regions_overlapped = me->regions_overlapped + "," + them->unique_id_2;
														me->percent_overlapped = me->percent_overlapped + "," + olap;
													}
													else{	//we don't match previous gene
														me->regions_overlapped = me->regions_overlapped
														+ "),"	//close previous ')'
														+ them->unique_id + "(" + them->unique_id_2;	//new 'gene(transcript' added
														me->percent_overlapped = me->percent_overlapped
														+ ")," 
														+ "(" + olap;
													}
												}
												previous_region=them;//for the next cycle through, previous_region points to previous region
												*/
						//get down to 4 digits for olap
						double olap_percent = overlap*100;
						ostringstream strs;
						strs << fixed << setprecision (2) << olap_percent;
						string olap = strs.str();	
						if(olap == "100.00") olap = "100";
						//insert gene and transcript. if gene's there it will just insert another transcript
						me->genes_overlapped.insert(pair<string, string>(them->unique_id,them->unique_id_2));
						//now insert overlap percentage for that transcript
						me->percent_overlapped.insert(pair<string, string>(them->unique_id_2,olap));
						//cout << "I just finished inserting pairs in the maps" << endl;
				
				}
			}
		}
	}
}

int writeGeneList(vector<CNV> master_list_cnvs, vector<REGION> regions)
{
	map<string,REGION_COUNT> overall_region_list;
	
	//loop through regions list make a map entry for each gene in region list
	for (vector<REGION>::iterator reg_itr = regions.begin();reg_itr != regions.end();reg_itr++)
	{
		//REGION_COUNT temp_reg;
		overall_region_list[reg_itr->unique_id];
	}	
	//now loop through all CNVs and their associated regions and count number of occurences of each gene	
	for (vector<CNV>::iterator me = master_list_cnvs.begin();me != master_list_cnvs.end();me++)
	{
		for(multimap<string,string>::iterator gene_list = me->genes_overlapped.begin();
			gene_list != me->genes_overlapped.end();
			gene_list=me->genes_overlapped.upper_bound(gene_list->first))	//for all genes
			{
				map<string,REGION_COUNT>::iterator found_itr;
				found_itr = overall_region_list.find(gene_list->first);
				if(found_itr == overall_region_list.end())
					return -1;
				else{
					if(me->case_control == 1)
					{
						found_itr->second.num_cases = found_itr->second.num_cases + 1;
					}
					if(me->case_control == 2)
						found_itr->second.num_controls = found_itr->second.num_controls + 1;
				}
			}
	}
	
	string out_gene  = out_file + ".GeneList.txt";
	
	ofstream outfile(out_gene.c_str());
	if(!outfile){
		 cerr << "Can not open the results file for writing!\nIt may be in use. Close it and try again\n";
		 return -2;
	}
	
	outfile << "region(gene)\tcase_count\tcontrol_count\n";
	
	for(map<string,REGION_COUNT>::iterator gene_list_itr = overall_region_list.begin();
		gene_list_itr != overall_region_list.end();
		gene_list_itr++)
	{
		if(gene_list_itr->second.num_cases ==0 && gene_list_itr->second.num_controls ==0)
			continue;
		else
			outfile << gene_list_itr->first << "\t"
					<< gene_list_itr->second.num_cases
					<< "\t" << gene_list_itr->second.num_controls << endl;
	}
		
	
	outfile.close();
	return 0;
}

int writeTransList(vector<CNV> master_list_cnvs, vector<REGION> regions)
{
	map<string,REGION_COUNT> overall_region_list;
	
	//loop through regions list make a map entry for each transcript in region list
	for (vector<REGION>::iterator reg_itr = regions.begin();reg_itr != regions.end();reg_itr++)
	{
		//REGION_COUNT temp_reg;
		overall_region_list[reg_itr->unique_id_2];
	}	
	//now loop through all CNVs and their associated regions and count number of occurences of each transcripts	
	for (vector<CNV>::iterator me = master_list_cnvs.begin();me != master_list_cnvs.end();me++)
	{
		for(map<string,string>::iterator trans_list = me->percent_overlapped.begin();
			trans_list != me->percent_overlapped.end();
			trans_list++)	//for all transcripts
			{
				map<string,REGION_COUNT>::iterator found_itr;
				found_itr = overall_region_list.find(trans_list->first);
				if(found_itr == overall_region_list.end())
					return -1;
				else{
					if(me->case_control == 1)
					{
						found_itr->second.num_cases = found_itr->second.num_cases + 1;
					}
					if(me->case_control == 2)
						found_itr->second.num_controls = found_itr->second.num_controls + 1;
				}
			}
	}
	
	string out_trans = out_file + ".TransList.txt";
	
	ofstream outfile(out_trans.c_str());
	if(!outfile){
		 cerr << "Can not open the results file for writing!\nIt may be in use. Close it and try again\n";
		 return -2;
	}
	
	outfile << "TranscriptID\tcase_count\tcontrol_count\n";
	
	for(map<string,REGION_COUNT>::iterator trans_list_itr = overall_region_list.begin();
		trans_list_itr != overall_region_list.end();
		trans_list_itr++)
	{
		if(trans_list_itr->second.num_cases ==0 && trans_list_itr->second.num_controls ==0)
			continue;
		else
			outfile << trans_list_itr->first << "\t"
					<< trans_list_itr->second.num_cases
					<< "\t" << trans_list_itr->second.num_controls << endl;
	}
		
	
	outfile.close();
	return 0;
}


bool WriteCnvs(vector<CNV> master_list_cnvs)
{
	//cout << "I got to writeCNVs" << endl;
	ofstream outfile(out_file.c_str());
	if(!outfile){
		 cerr << "Can not open the results file for writing!\nIt may be in use. Close it and try again\n";
		 return false;
	}
	outfile << "OurID\tChr\tStart\tStop\tCase/Control\tCopy_Number\tOverlapped_Regions\tPercent_Overlap\tNumber_Regions_Overlapped\n";
	if(write_transcript_id==0)
	{
		cout << "no write transcript id" << endl;
		cout << "Sorry haven't implemented this yet" << endl;
		cout << "Try running the program again with write_transcript_id turned on and then"
				<< "search and replace in excel for [*] (replacing with nothing)... Not perfect, but I'll get to output with only gene name soon\n";
	}
	if(write_transcript_id==1)
	{	
		
		//cout << "I got to the write_transcript_id if statement in the writeCNVs function" << endl;
		for (vector<CNV>::iterator me = master_list_cnvs.begin();me != master_list_cnvs.end();me++){
			int region_count=0;
			outfile << me->our_id << "\t" << me->chr << "\t" 
					<< me->start << "\t" << me->stop << "\t" 
					<< me->case_control << "\t" << me->copy_num << "\t";
			string gene_string="", overlap_string="";
			multimap<string,string>::iterator gene_list = me->genes_overlapped.begin();
			map<string,string>::iterator percent_overlap_list;
			//cout << "I defined my iterators\n";
			for(;gene_list != me->genes_overlapped.end();gene_list=me->genes_overlapped.upper_bound(gene_list->first))	//for all genes
			{
				region_count++;
				//cout << "I got to the for loop inside write\n";
				if(gene_string!="")	//if our gene list is not empty
					gene_string = gene_string + ";";	//add ; at end
				gene_string = gene_string + gene_list->first//add gene name to our gene_string
							+ "[";							//and get ready for transcript ids	
				if(overlap_string!="")
					overlap_string = overlap_string + ";";
								
				//temp iterator for the next for
				multimap<string,string>::iterator temp_gene_list=gene_list;
				
				//now i=0, keep going until i gets to the number of transcripts for this gene
				//cout << me->genes_overlapped.count(gene_list->first) << endl;
				for(int i = 0;i < me->genes_overlapped.count(gene_list->first);i++,temp_gene_list++)
				{
					//cout << temp_gene_list->first << "(" << temp_gene_list->second << ")" << endl;
					//cout << (me->percent_overlapped.find(temp_gene_list->second))->second << endl;
					//cout << ")" << me->percent_overlapped.find(temp_gene_list->second) << endl;
					if(i==0)	//first transcript
						overlap_string = overlap_string + "[";
					else	//not our first transcript
					{									
						gene_string = gene_string + ",";	//append comma
						overlap_string = overlap_string + ",";
					}
					gene_string = gene_string + temp_gene_list->second;	//append transcript name
					percent_overlap_list = me->percent_overlapped.find(temp_gene_list->second);	//get iterator for this transcripts map
					
						/*if(percent_overlap_list->first == "ENST00000240189.2")
							cout << percent_overlap_list->second << endl;*/
						
					overlap_string = overlap_string + percent_overlap_list->second;
				}
				//must be done with transcripts for this gene, now lets close the )
				gene_string = gene_string + "]";
				overlap_string = overlap_string + "]";
	
			}
			if(gene_string == "")
				gene_string = "NA";
			if(overlap_string == "")
				overlap_string = "NA";
				
			outfile << gene_string << "\t" << overlap_string << "\t" << region_count << endl;	
		}
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
	outfile << "Program: " <<  VERSION << endl;
	outfile << "Number CNVs read: " << num_cnvs << endl;
	outfile << "Number regions(genes) read: " << num_regions << endl;
	outfile << "CNV file used: " << master_cnv_file << endl;
	outfile << "Region of interest file used: " << comparison_cnv_file << endl;
	outfile << "Overlap type specified: " << overlap_type << endl;
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
		{
			float temp_overlap_cutoff;
			ss >> temp_overlap_cutoff;
			//cout << temp_overlap_cutoff << endl;
			if(temp_overlap_cutoff > 1){
				cerr << "WHOA! check your config file. I think you might have used 0-100 for overlap cutoff instead of 0-1... I'll convert it\n";
				cerr << "If something other than you using 0-100 for overlap cutoff happened, fix the config file and rerun the program\n";
				temp_overlap_cutoff = temp_overlap_cutoff/100;
			}
			 percent_overlap_cutoff = temp_overlap_cutoff;
		}
		if(line_num ==5)
			ss >> overlap_type;
		if(line_num ==9)
			ss >> write_transcript_id;
		if(line_num ==13)
			ss >> master_cnv_file;
		if(line_num ==18)
			ss >> comparison_cnv_file;
		if(line_num ==23)
			ss >> out_file;
		line_num++;
	}
	if (line_num != 31){
		//cout << line_num << endl;
		cerr << "Looks like the config file is corrupt! Did you change what lines the data was on?\n";
		return false;
	}
		//cout << "I got to the while loop for file reading\n";
	infile.close();
	//ConsolePrintCnvs(our_cnvs);
	return true;
}



int main()
{
	vector<CNV> master_list_cnvs;	//main vectors of CNVs for comparing
	vector<REGION> regions;

	string filename;
	cout << "***********************\n***********************\n"
	     << VERSION << endl
		 << "\n***********************\n***********************\n";
	cout << "\nCNV list file should be a tab-delimited file in the following format:\n";
	cout << "OurID\tChr\tStart\tStop\tCase/Control\tCopy_Number\n";
	cout << "\nRegion of interest file should be a tab-delimited file in the following format:"
			<< "\t\tNOTE:ADDED A FIELD FOR UNIQUE_ID_2, THIS CAN BE USED FOR TRANSCRIPTS\n";
	cout << "Database_From\tUnique_ID_1(gene)\tUnique_ID_2(transcript)\tChr\tStart\tStop\n";
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
		if(!ReadFileRegion(comparison_cnv_file,regions)){
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
	CheckRegions(master_list_cnvs, regions);	
	if(WriteCnvs(master_list_cnvs))
		cout << "All set. Check you're output file: " << out_file << endl;
	if(writeGeneList(master_list_cnvs,regions) == -1)
		cerr << "ERROR writing Aux file... Looks like we have gene name in our list that wasn't in our regions list!"
				<< "\nThis should never happen, so something must be going on.\n";
	if(writeTransList(master_list_cnvs,regions) == -1)
		cerr << "ERROR writing Aux file... Looks like we have gene name in our list that wasn't in our regions list!"
				<< "\nThis should never happen, so something must be going on.\n";
	int num_cnvs = master_list_cnvs.size(),
		num_regions = regions.size();
	WriteLog(num_cnvs, num_regions);			
	//ConsolePrintCnvs(master_list_cnvs);
	return 0;
}
	
	