#ifndef GENOMICREGION_H
#define GENOMICREGION_H

#include <string>

using namespace std;
 
class GenomicRegion{

public:

	//DGV2 entry constructor
		//GenomicRegion(string const& description, string const& unique_id, string const& chromosome, int start, int stop, int copy_num)
	GenomicRegion(string const&, int, int, int, double);
	
	
	//constructor for CNVs from cumulative DB
	GenomicRegion(string const& our_id, string const& chromosome, int start, int stop, int copy_num, int case_or_control);
	
	//constructor for one of our CNVs for compare algorithms
	//				our_id >>     study_id >> case_or_control >> chromosome >> start >> stop >> copy_num;
	GenomicRegion(string const&, string const&, string const&, int, string const&, int, int, int);
	
	//Gene or transcript constructor
		//GenomicRegion::GenomicRegion(string const& description, string const& unique_id, string const& unique_id_2, string const& chromosome, int start, int stop)
	GenomicRegion(string const&, string const&, string const&, string const&, int, int);
		
	//constructor for cnvr block
	////Block_ID	chr	start	stop	freq
	GenomicRegion(string const& block_id, string const& chromosome, int start, int stop, string const& freq);
	
	//constructor for most basic type: chr, start, stop
	GenomicRegion(string const&, int, int);
	
	
	//operators >, <, >=, <= based on the start coordinate for each
	friend bool operator<(const GenomicRegion& lhs, const GenomicRegion& rhs);
	friend bool operator>(const GenomicRegion& lhs, const GenomicRegion& rhs);
	friend bool operator<=(const GenomicRegion& lhs, const GenomicRegion& rhs);
	friend bool operator>=(const GenomicRegion& lhs, const GenomicRegion& rhs);
	
	unsigned int GetStart() const;
	unsigned int GetStop() const;
	string GetChromosome() const;
	int GetCopyNum() const;
	int GetCaseControlStatus() const;
	string GetUniqueId() const;
	string GetUniqueId2() const;
	string GetDescription() const;
		
private:
	string description_;	//identifier of type
	string unique_id_;
	string unique_id_2_;
	string chromosome_;
	double freq_;
	int copy_num_, case_or_control_;
	unsigned int start_, stop_;	//start and end points for each region
	unsigned int length_;	//length of region

};



#endif