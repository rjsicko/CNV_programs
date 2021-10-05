/*****************************************************
generate_num_snps v1.0
Author: Robert Sicko - R.J.Sicko@gmail.com

****DESCRIPTION****
This program takes a list of CNVpartition CNVs (in tab format, with headers) 
and using the PFB file specified, generates the approximate number of SNPs for each CNV.
Only approximate because we use all markers in the PFB file, including those omitted from analysis.
Could be corrected by modifying PFB file to include only those used.

****END DESCRIPTION****

****REQUIRED FORMAT****
	NOTE: No commas in numbers (ie. use 1000000 not 1,000,000)
	NOTE2: Keep headers in your files. The program ignores the first line (don't have spaces in header id's though).
	
	CNV List (tab-delimited txt):
	chr  start  end  copy_num  sampleID  startSNP  endSNP  conf  num_SNP
	
	NOTE: this program functions correctly only if the PFB file is sorted by Chr & Position
	PFB file (tab-delimited txt):
	Name	Chr	Position	PFB
****END REQUIRED FORMAT****

****CHANGE LOG****
v1.0- Bugs fixed, programs appears to be functioning as expected

	
****END CHANGE LOG****

****TO DO****
	- Comment out code!
	- See if we can get exact number of SNPs
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



//A struct for our CNVs (or regions)
//"chr \t start \t end \t copy_num \t sampleID \t startSNP \t endSNP \t conf \t num_SNP\n";
struct CNV{
	string startSNP, endSNP, sampleID;	//our identifier
	int chr, copy_num, num_SNP;	//ints for chromosome, copy number and case or control (1, 2)
	long start, stop;
	double conf;
};

struct MARKER{
	string snp_id;
	bool included;
	MARKER(string temp_id) : snp_id(temp_id), included(false) { }
	 bool operator==(const MARKER& r) const
	{
       return this->snp_id == r.snp_id;
    }
  
};


bool ReadFileCnv(string filename,vector<CNV>& our_cnvs)
{
		//"chr \t start \t end \t copy_num \t sampleID \t startSNP \t endSNP \t conf \t num_SNP\n";
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
		ss >> dummy_cnv.chr >> dummy_cnv.start
			>> dummy_cnv.stop >> dummy_cnv.copy_num
			>> dummy_cnv.sampleID >> dummy_cnv.startSNP
			>> dummy_cnv.endSNP >> dummy_cnv.conf
			>> dummy_cnv.num_SNP;

		
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

bool ReadFileSNPS(string filename,vector< vector<MARKER> > &SNPS )
{
	ifstream infile(filename.c_str());
	if(!infile.is_open())
		return false;
	bool first_line=true;
	
	string line;
	while (getline(infile, line))
	{
		string snp_id;
		int chr;
		istringstream ss(line);
		//string name;
		//int var1, var2, var3;
		ss >> snp_id >> chr;

		if(first_line)
			first_line=false;
		else{
			MARKER temp_snp(snp_id);
			SNPS[chr].push_back(temp_snp);
		}
	}
		
		//cout << "I got to the while loop for file reading\n";
	infile.close();
	//ConsolePrintCnvs(our_cnvs);
	return true;
}

			//NOT FUNCTIONING
			/*bool TrimSNPList(string filename,vector< vector<MARKER> > &SNPS )
			{
				//Name(marker)	C06075.Log R Ratio	C06075.B Allele Freq
				ifstream infile(filename.c_str());
				if(!infile.is_open())
					return false;
				bool first_line=true;
				
				string line;
				int line_num=0;
				while (getline(infile, line))
				{
					if(line_num == 500,000)
						cout << "Almost half way there... 500,000 lines read!\n";
					line_num++;
					string snp_id;
					istringstream ss(line);
					//string name;
					//int var1, var2, var3;
					ss >> snp_id;

					if(first_line)
						first_line=false;
					else{
							//vector<MARKER>::iterator found = SNPS[0].end()
						//look through all chr positions
						int i = 0;
						//cout << "I'm just before the for loop\n";
						for (vector< vector<MARKER> >::iterator itr = SNPS.begin();itr != SNPS.end();itr++, i++)
							//and all markers in the vector at each positions
							for (vector<MARKER>::iterator j = SNPS[i].begin();j != SNPS[i].end();j++)
							{
								if(snp_id == j->snp_id)
									j->included = true;
							}
										
					}
				}	
					//cout << "I got to the while loop for file reading\n";
				infile.close();
				//ConsolePrintCnvs(our_cnvs);
				return true;
			}*/

	


void CountSNPS(vector<CNV>& master_list_cnvs, vector< vector<MARKER> > SNPS)
{
	//cout << type_specific << endl;
	for (vector<CNV>::iterator me = master_list_cnvs.begin();me != master_list_cnvs.end();me++)
	{
		vector<MARKER>::iterator foundSNP = find(SNPS[me->chr].begin(), SNPS[me->chr].end(), me->startSNP);
		if(foundSNP == SNPS[me->chr].end())
			me->num_SNP = -1;	//startSNP wasn't found so some error
		else{
			int snp_count=0;
			while(foundSNP != SNPS[me->chr].end() && foundSNP->snp_id != me->endSNP)
			{
				//if(foundSNP->included)
				snp_count++;
				foundSNP++;
			}
			snp_count++;
			me->num_SNP = snp_count;
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
		//"chr \t start \t end \t copy_num \t sampleID \t startSNP \t endSNP \t conf \t num_SNP\n";
	
	outfile << "Chr\tStart\tStop\tCopy_Number\tsampleID\tstartSNP\tendSNP\tconf\tapprox_num_SNPs\n";
		
		//cnv struct
		/*string startSNP, endSNP, sampleID;	//our identifier
		int chr, copy_num, num_SNP;	//ints for chromosome, copy number and case or control (1, 2)
		long start, stop, conf;	//start and end points for each cnv*/		
	
	for (vector<CNV>::iterator me = master_list_cnvs.begin();me != master_list_cnvs.end();me++){
		outfile << me->chr << "\t" << me->start << "\t" 
				<< me->stop << "\t" << me->copy_num << "\t" 
				<< me->sampleID << "\t" << me->startSNP << "\t"
				<< me->endSNP << "\t" << me->conf << "\t"
				<< me->num_SNP << endl;
	}
	outfile.close();	//close the file
	return true;
}

bool WriteLog(int num_cnvs)
{
	//cout << "in WriteLog\n";
	string out_trans = out_file + ".Log";
	ofstream outfile(out_trans.c_str());
	if(!outfile){
		 cerr << "Can not open the log file for writing!\nIt may be in use. Close it and try again\n";
		 return false;
	}
	outfile << "Program: generate_num_snps_v1.0\n";
	outfile << "Number CNVs read: " << num_cnvs << endl;
	outfile << "CNV file used: " << master_cnv_file << endl;
	outfile << "PFB file used: " << pfb_file << endl;
	outfile << "Sample file used (determines which markers from PFB to use): " << sample_file << endl;
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
	ifstream infile("generate_num_snps_v1.0.config");	//open file
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
		if(line_num ==4)
			ss >> pfb_file;
		//if(line_num ==7)
			//ss >> sample_file;
		if(line_num ==10)
			ss >> out_file;
		line_num++;
	}
	if (line_num != 15){
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
	vector< vector<MARKER> > SNPS (27);

	string filename;
	cout << "***********************\n***********************\nGenerate Number SNPS v1.0\n***********************\n***********************\n";
	cout << "WARNING: Different file format from our other programs\n";
	cout << "\nCNV list file should be a tab-delimited file in the following format (PennCNV tab format with headers):\n";
	cout << "chr \t start \t end \t copy_num \t sampleID \t startSNP \t endSNP \t conf \t num_SNP\n";
	cout << "PFB file should be a tab-delimited file in the following format:";
	cout << "Name \t Chr \t Position \t PFB\n";
	cout << "NOTE: this program functions correctly only if the PFB file is sorted by Chr & Position\n";
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
		if(!ReadFileSNPS(pfb_file,SNPS ))
		{
			//throw runtime_error("Cannot find config item.");
			
			cerr << "Looks like there was a problem reading " << pfb_file;
			cerr << " Make sure you updated the generate_num_snps_v1.0.config file\n";
			cerr << "Also, the cnv files need to be in the same directory as the program\n";
			//cerr << "\nPerhaps you didn't specify the path?\n";
			return 1;
		}
	}
	else{
		cerr << "Looks like there was a problem reading " << master_cnv_file;
		cerr << " Make sure you updated the generate_num_snps_v1.0.config file\n";
		cerr << "Also, the cnv files need to be in the same directory as the program\n";		
		//cerr << "\nPerhaps you didn't specify the path?\n";
		return 1;
	}
	
	/*if (!TrimSNPList(sample_file,SNPS))
	{
		cerr << "Looks like there was a problem reading " << sample_file;
		cerr << " Make sure you updated the generate_num_snps_v1.0.config file\n";
	}*/	
	
	
	cout << "Running now\n";
	
	CountSNPS(master_list_cnvs, SNPS);
	if(WriteCnvs(master_list_cnvs))
		cout << "All set. Check you're output file: " << out_file << endl;
	
	int num_cnvs = master_list_cnvs.size();
	WriteLog(num_cnvs);			
	//ConsolePrintCnvs(master_list_cnvs);
	return 0;
}
	
	