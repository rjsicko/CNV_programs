/*****************************************************
Generate_Tracks_v1.0
Author: Robert Sicko - R.J.Sicko@gmail.com

****DESCRIPTION****
This program generates UCSC and DGV browser tracks from our data
****END DESCRIPTION****

****REQUIRED FORMAT****
	NOTE: No commas in numbers (ie. use 1000000 not 1,000,000).
	NOTE2: Keep headers in your files. The program ignores the first line (don't have spaces in header id's though).
	
	File should be tab-delimited in the following format(tab-delimited txt):
	OurID	Subject_id		Chr		Start	Stop	Case/Control	Copy_Number
****END REQUIRED FORMAT****


****CHANGE LOG****
v1.0- Bugs fixed, programs appears to be functioning as expected
v1.1- Changed chr to strign to handle sex chr
	- Added const string VERSION variable for tracking what version we are
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


using namespace std;

const string VERSION = "generate_tracks_v1.1";

//A struct for our CNVs (or regions)
struct CNV{
	string our_id, subject_id, chr;	//our identifier
	int copy_num, case_control;	//ints for chromosome, copy numbe and case or control (1, 2)
	long start, stop;	//start and end points for each cnv
};

template <typename T>
  string NumberToString ( T Number )
  {
     ostringstream ss;
     ss << Number;
     return ss.str();
  }


bool WriteLog(int num_cnvs)
{
	string out_trans = out_file + ".Log";
	ofstream outfile(out_trans.c_str());
	if(!outfile){
		 cerr << "Can not open the log file for writing!\nIt may be in use. Close it and try again\n";
		 return false;
	}
	outfile << "Program: generate_tracks_v1.0\n";
	outfile << "Number CNVs read: " << num_cnvs << endl;
	outfile << "CNV file used: " << master_cnv_file << endl;
	time_t t = time(0);   // get time now
    struct tm * now = localtime( & t );
	outfile	<< "Analysis ran: " 
			<< (now->tm_year + 1900) << '-' 
			<< (now->tm_mon + 1) << '-'
			<<  now->tm_mday
			<< endl;
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
		ss >> dummy_cnv.our_id >> dummy_cnv.subject_id >> dummy_cnv.chr
			>> dummy_cnv.start >> dummy_cnv.stop
			>> dummy_cnv.case_control >> dummy_cnv.copy_num;


		if(first_line)
			first_line=false;
		else
			our_cnvs.push_back (dummy_cnv);
		
	}
	infile.close();
	return true;
}

string ThrowString(CNV the_cnv)
{
	if(the_cnv.copy_num == 0)	//homodel
		return "HomoDel";    	
	if(the_cnv.copy_num == 1)	//hetdel
		return "HetDel";     	
	if(the_cnv.copy_num == 2)	//loh
		return "LOH";        	
	if(the_cnv.copy_num == 3)	//dupl
		return "Dupl";       	
	if(the_cnv.copy_num == 4)	//tripl
		return "Tripl";
}


string ThrowColor(CNV the_cnv)	
{
	//returns a string for r,g,b numbers for each cnv
	if(the_cnv.copy_num == 0)	//homodel
		return "255,154,154";
	if(the_cnv.copy_num == 1)	//hetdel
		return "255,128,0";
	if(the_cnv.copy_num == 2)	//loh
		return "0,255,0";
	if(the_cnv.copy_num == 3)	//dupl
		return "0,0,255";
	if(the_cnv.copy_num == 4)	//tripl
		return "102,0,204";
}

bool WriteUCSCTrack(vector<CNV> master_list_cnvs)
{

	string out_ucsc = out_file + ".UCSC_Track_" + NumberToString(master_list_cnvs.size()) + ".txt";
	ofstream outfile(out_ucsc.c_str());
	if(!outfile){
		 cerr << "Can not open the UCSC track file for writing!\nIt may be in use. Close it and try again\n";
		 return false;
	}
	
	
	outfile << "track name=\"" << out_ucsc << "\" description=\"" << out_ucsc << "\" visibility=2 itemRgb=\"On\"" << endl;
	for (vector<CNV>::iterator me = master_list_cnvs.begin();me != master_list_cnvs.end();me++){
		string combined_name = me->our_id + "_" + me->subject_id + "_";
		if(me->case_control == 1)
			combined_name = combined_name + "Case";
		else
			combined_name = combined_name + "Control";
		combined_name = combined_name + "_" + ThrowString(*me);
		//UCSC uses BED format which is 0 based coordinates so subtract 1 from our coordinates
		outfile << "chr" << me->chr << "\t" << (me->start)-1
				<< "\t" << (me->stop)-1 << "\t" << combined_name
				<< "\t0\t+\t" << (me->start)-1 << "\t" << (me->stop)-1 
				<< "\t" << ThrowColor(*me) << endl;

	}

	outfile.close();	//close the file
	return true;
}

bool WriteDGVTrack(vector<CNV> master_list_cnvs)
{

	string out_dgv = out_file + ".DGV_Track_" + NumberToString(master_list_cnvs.size()) + ".txt";
	ofstream outfile(out_dgv.c_str());
	if(!outfile){
		 cerr << "Can not open the DGV track file for writing!\nIt may be in use. Close it and try again\n";
		 return false;
	}
	
	outfile << "[HomoDel]\nbgcolor= pink\n\n[HetDel]\nbgcolor= orange\n\n[Dupl]\nbgcolor= blue\n\n[Tripl]\nbgcolor= purple\n\n[LOH]\nbgcolor= green\n\n";
	
	for (vector<CNV>::iterator me = master_list_cnvs.begin();me != master_list_cnvs.end();me++){
		string combined_name = me->our_id + "_" + me->subject_id + "_";
		if(me->case_control == 1)
			combined_name = combined_name + "Case";
		else
			combined_name = combined_name + "Control";
		combined_name = combined_name + "_" + ThrowString(*me);
		string combined_location = "chr" +  me->chr + ":" + NumberToString(me->start) + ".." + NumberToString(me->stop); 
		outfile << ThrowString(*me) << "\t" << combined_name
				<< "\t" << combined_location << endl;
	}	
	
	outfile.close();	//close the file
	return true;
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
		//I know... BAD practice to hardcode line numbers here
		//but it's all I could come up with to make it work in a resonable amount of time
		if(line_num ==1)
			ss >> master_cnv_file;
		if(line_num ==4)
			ss >> out_file;
		line_num++;
	}
	infile.close();
	//ConsolePrintCnvs(our_cnvs);
	return true;
}



int main()
{
	vector<CNV> master_list_cnvs;	//main vectors of CNVs for comparing

	string filename;
	cout << "***********************\n***********************\n"
	     << VERSION << endl
		 << "\n***********************\n***********************\n";
	cout << "\nCNV list files should be a tab-delimited file in the following format(NOTE ADDITIONAL 'subject_id' FIELD):\n";
	cout << "OurID\tSubject_id\tChr\tStart\tStop\tCase/Control\tCopy_Number\n";
	cout << "\nKeep headers in your input files! WARNING: For now remove commas from numbers. ie dont use 1,000,000, use 1000000 in input files\n";
	cout << "I will try to deal with commas in numbers in later versions, but for now leave them out\n";
	cout << "Press enter to continue\n";
	cin.get();

	if(!DefineConfig())	cerr << "Unable to read the config file!\n";
	if(ReadFileCnv(master_cnv_file,master_list_cnvs))
	{
		cout << "Running now\n";
		WriteUCSCTrack(master_list_cnvs);
		WriteDGVTrack(master_list_cnvs);
	}
	else{
		cerr << "Looks like there was a problem reading " << master_cnv_file;
		cerr << " Make sure you updated the compare_lists.config file with the two files you want to compare\n";
		cerr << "Also, the cnv files need to be in the same directory as the program\n";		
		return 1;
	}
	

	
	cout << "All set. Check you're output files: " << out_file << endl;
		
	int num_cnvs = master_list_cnvs.size();
	WriteLog(num_cnvs);		

	return 0;
}