#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include "cnv.h"
#include "config.h"
#include "genomicregion.h"

using namespace std;

//NOTE LOST MAIN.CPP for GVAT0.0.1 have to recreate it! this is a partially working version

//make config global since we use it all over
Config ourconfig("src/GVAT-0.0.2.config");

string cnv_filename = ourconfig.GetFilename("cnv_filename");
string DGV2_filename = ourconfig.GetFilename("DGV2_filename");
string other_algorithm_filename = ourconfig.GetFilename("other_algorithm_filename");
string gencode_filename = ourconfig.GetFilename("gencode_filename");
string ccds_filename = ourconfig.GetFilename("ccds_filename");
string omim_filename = ourconfig.GetFilename("omim_filename");
string cumulative_database_filename = ourconfig.GetFilename("cumulative_database_filename");
string chop_cnvs_filename = ourconfig.GetFilename("chop_cnvs_filename");
string hapmap_cnvs_filename = ourconfig.GetFilename("hapmap_cnvs_filename");
string pathway_genes_filename = ourconfig.GetFilename("pathway_genes_filename");
string pseudo_gene_filename = ourconfig.GetFilename("pseudo_gene_filename");
string pseudo_parent_filename = ourconfig.GetFilename("pseudo_parent_filename");
string out_filename = ourconfig.GetFilename("out_filename");
double case_control_overlap_threshold = 0.50;
	
bool ReadDGV2File(vector<GenomicRegion>& dgv2_vector)
{
	ifstream infile(DGV2_filename.c_str());
	if(!infile.is_open())
		return false;	//file failed to open

	bool first_line = true;
	string line;
	while (getline(infile, line))	//while there are still lines left to read
	{
		//DGV_date	chr	start	end	variantsubtype	frequency
		string chromosome, junk;
		int start, stop, copy_num;
		double freq;
		
		istringstream ss(line);
		ss >> junk >> chromosome
		    >> start >> stop >> copy_num >> freq;
		
		if(first_line)	//ignore header line
			first_line=false;
		
		else{	//we already skipped header, start adding the entries
			//GenomicRegion(string const& description, string const& unique_id, string const& chromosome, int start, int stop, int copy_num)
			GenomicRegion dummy_dgv2(chromosome, start, stop, copy_num, freq);
			dgv2_vector.push_back(dummy_dgv2);
		}	
	}
	infile.close();
	return true;
}


//generic read region function
bool ReadRegionFile(vector<GenomicRegion>& region_of_interest_vector, const string& region_of_interest_filename)
{
	ifstream infile(region_of_interest_filename.c_str());
	if(!infile.is_open())
		return false;	//file failed to open

	bool first_line = true;
	string line;
	while (getline(infile, line))	//while there are still lines left to read
	{
		string database_id, unique_id, unique_id_2, chromosome;
		int start, stop;
		
				//datbase_from	type	gene_name	chr	start	end
					//>> dummy_region.chr	>> dummy_region.start >> dummy_region.stop;
		//database_from	gene	gene_id	chr	txStart	txEnd	
		istringstream ss(line);
		ss >> database_id >> unique_id >> unique_id_2 
		   >> chromosome >> start >> stop;
		
		if(first_line)	//ignore header line
			first_line=false;
		
		else{	//we already skipped header, start adding the entries
			//GenomicRegion(string const& description, string const& unique_id, string const& unique_id_2, string const& chromosome, int start, int stop)
			GenomicRegion dummy_region(database_id, unique_id, unique_id_2, chromosome, start, stop);
			region_of_interest_vector.push_back(dummy_region);
		}	
	}
	infile.close();
	return true;
}

/*
bool ReadRegionPseudo(vector<GenomicRegion>& region_of_interest_vector, const string& region_of_interest_filename)
{
	ifstream infile(region_of_interest_filename.c_str());
	if(!infile.is_open())
		return false;	//file failed to open

	bool first_line = true;
	string line;
	while (getline(infile, line))	//while there are still lines left to read
	{
		string database_id, unique_id, unique_id_2, chromosome;
		int start, stop;
		//Database_from	ID	Chromosome	Start_Coordinate	Stop_Coordinate
			
		istringstream ss(line);
		ss >> database_id >> unique_id  
		   >> chromosome >> start >> stop;
		
		unique_id_2 = "";
		if(first_line)	//ignore header line
			first_line=false;
		
		else{	//we already skipped header, start adding the entries
			//GenomicRegion(string const& description, string const& unique_id, string const& unique_id_2, string const& chromosome, int start, int stop)
			GenomicRegion dummy_region(database_id, unique_id, unique_id_2, chromosome, start, stop);
			region_of_interest_vector.push_back(dummy_region);
		}	
	}
	infile.close();
	return true;
}*/


//read function for cumulative DB
bool CumulativeDBFile(vector<GenomicRegion>& cumulative_db)
{
	ifstream infile(cumulative_database_filename.c_str());
	if(!infile.is_open())
		return false;	//file failed to open

	bool first_line = true;
	string line;
	while (getline(infile, line))	//while there are still lines left to read
	{
		//GenomicRegion(string const& our_id, string const& chromosome, int start, int stop, int copy_num, int case_or_control);
		string our_id, chromosome, dummy;
		int start, stop, copy_num, case_or_control;
		
		istringstream ss(line);
		//Our_CNV_ID	Study_ID	Chromosome	Start	End	Case_or_Control	Copy_Number
		ss >> our_id >> dummy >> chromosome
			>> start >> stop >> case_or_control >> copy_num;
		
		if(first_line)	//ignore header line
			first_line=false;
		
		else{	//we already skipped header, start adding the entries
			//constructor for CNVs from cumulative DB
			//GenomicRegion(string const& our_id, string const& chromosome, int start, int stop, int copy_num, int case_or_control);
			GenomicRegion dummy_region(our_id, chromosome, start, stop, copy_num, case_or_control);
			cumulative_db.push_back(dummy_region);
		}	
	}
	infile.close();
	return true;
}

//read function for hapmap/chop
bool ReadBlock(vector<GenomicRegion>& region_of_interest_vector, const string& region_of_interest_filename)
{
	ifstream infile(region_of_interest_filename.c_str());
	if(!infile.is_open())
		return false;	//file failed to open

	bool first_line = true;
	string line;
	while (getline(infile, line))	//while there are still lines left to read
	{
		string unique_id, chromosome, freq;
		int start, stop;
		
		////Block_ID	chr	start	stop	freq
			//Database_from	ID	Chromosome	Start_Coordinate	Stop_Coordinate
			
		istringstream ss(line);
		ss >> unique_id >> chromosome >> start >> stop >> freq;
		
		if(first_line)	//ignore header line
			first_line=false;
		
		else{	//we already skipped header, start adding the entries
			//(string const& block_id, string const& chromosome, int start, int stop, string const& freq)
			GenomicRegion dummy_region(unique_id, chromosome, start, stop, freq);
			region_of_interest_vector.push_back(dummy_region);
		}	
	}
	infile.close();
	return true;
}

bool ReadOtherAlgorithm(vector<GenomicRegion>& region_of_interest_vector, const string& region_of_interest_filename)
{
	ifstream infile(region_of_interest_filename.c_str());
	if(!infile.is_open())
		return false;	//file failed to open

	bool first_line = true;
	string line;
	while (getline(infile, line))	//while there are still lines left to read
	{
		string algorithm, our_id, study_id, chromosome;
		int start, stop, case_or_control, copy_num;

		//algorithm	our_id	sampleID	case/control	chr	start	end	copy_num	
		istringstream ss(line);
		ss >> algorithm >> our_id
			>> study_id >> case_or_control >> chromosome 
			>> start >> stop >> copy_num;

		if(first_line)	//ignore header line
			first_line=false;
		
		else{	//we already skipped header, start adding the entries
			//algorithm	our_id	sampleID	case/control	chr	start	end	copy_num
			GenomicRegion dummy_region(algorithm, our_id, study_id, case_or_control, chromosome, start, stop, copy_num);	
			region_of_interest_vector.push_back(dummy_region);
		}	
	}
	infile.close();
	return true;
	//algorithm, string const& our_id, string const& study_id, int case_or_control, string const& chromosome, int start, int stop, int copy_num)
}

int main()
{
	cout << "******************************************\nGVAT-0.0.2\n******************************************\n";
	cout << "Press enter to continue\n";
	cin.get();
	CnvVector main_cnv_list;
	vector<GenomicRegion> other_algorithm_cnvs;
	vector<GenomicRegion> dgv2_vector;
	vector<GenomicRegion> gencode_genes;
	vector<GenomicRegion> ccds_genes;
	vector<GenomicRegion> omim_genes;
	vector<GenomicRegion> cumulative_db;
	vector<GenomicRegion> chop_regions;
	vector<GenomicRegion> hapmap_regions;
	vector<GenomicRegion> pathway_genes;
	vector<GenomicRegion> pseudo_genes;
	vector<GenomicRegion> pseudo_parent_genes;
	
	
	if(!main_cnv_list.ReadCnvFile(cnv_filename))
	{
		cerr << "Looks like there was a problem reading " << cnv_filename;
		return 1;
	}
	//other algorithm list
	else if(!ReadOtherAlgorithm(other_algorithm_cnvs, other_algorithm_filename))
	{
		cerr << "Looks like there was a problem reading " <<  other_algorithm_filename;
		return 1;
	}
	//dgv2
	else if(!ReadDGV2File(dgv2_vector))	//dgv2 file reading failed
	{
		cerr << "Looks like there was a problem reading " << DGV2_filename;
		return 1;
	}
	//gencode
	else if(!ReadRegionFile(gencode_genes, gencode_filename))	//failed to open gencode gene file
	{
		cerr << "Looks like there was a problem reading " << gencode_filename;
		return 1;
	}
	//ccds
	else if(!ReadRegionFile(ccds_genes, ccds_filename))		//failed to open ccds gene file
	{
		cerr << "Looks like there was a problem reading " << ccds_filename;
		return 1;
	}
	//omim
	else if(!ReadRegionFile(omim_genes, omim_filename))		//failed to open omim gene file
	{
		cerr << "Looks like there was a problem reading " << omim_filename;
		return 1;
	}
	//cumulativeDB
	else if(!CumulativeDBFile(cumulative_db))	//failed to open cumulative DB file
	{
		cerr << "Looks like there was a problem reading " << cumulative_database_filename;
		return 1;
	}
	//chop
	else if(!ReadBlock(chop_regions, chop_cnvs_filename))	
	{
		cerr << "Looks like there was a problem reading " << chop_cnvs_filename;
		return 1;
	}
	//hapmap
	else if(!ReadBlock(hapmap_regions, hapmap_cnvs_filename))	
	{
		cerr << "Looks like there was a problem reading " << hapmap_cnvs_filename;
		return 1;
	}
	//pathway genes
	else if(!ReadRegionFile(pathway_genes, pathway_genes_filename))	
	{
		cerr << "Looks like there was a problem reading " << pathway_genes_filename;
		return 1;
	}
	//pseudogenes
	else if(!ReadBlock(pseudo_genes, pseudo_gene_filename))
	{
		cerr << "Looks like there was a problem reading " << pseudo_gene_filename;
		return 1;
	}
	//parent genes
	else if(!ReadBlock(pseudo_parent_genes, pseudo_parent_filename))	
	{
		cerr << "Looks like there was a problem reading " << pseudo_parent_filename;
		return 1;
	}
	
	
	main_cnv_list.CalculateCaseControl(case_control_overlap_threshold);
	main_cnv_list.CalculateDGV2Overlap(dgv2_vector);
	main_cnv_list.CompareToAnotherAlgorithm(other_algorithm_cnvs);
	//main_cnv_list.ReadCnvFile(cnv_filename);
	main_cnv_list.CheckRegionOfInterest(ccds_genes, "ccds_genes", 2);
	main_cnv_list.CheckRegionOfInterest(omim_genes, "omim_genes", 2);
	main_cnv_list.CheckRegionOfInterest(gencode_genes, "gencode_genes", 2);
	main_cnv_list.CheckRegionOfInterest(pathway_genes, "pathway_genes", 2);
	main_cnv_list.CheckCumulativeDatabaseOverlap(cumulative_db);
	main_cnv_list.CheckCnvr(chop_regions, "chop");
	main_cnv_list.CheckCnvr(hapmap_regions, "hapmap");
	main_cnv_list.CheckCnvr(pseudo_genes, "pseudogene");
	main_cnv_list.CheckCnvr(pseudo_parent_genes, "parentgene");
	main_cnv_list.WriteSelfToFile(out_filename);
	
	
	if(!main_cnv_list.WriteSelfToFile(out_filename))
	{
		cerr << "Looks like there was a problem reading " << out_filename;
		return 1;
	}		
	
}
		