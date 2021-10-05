#ifndef CNV_H
#define CNV_H

#include <string>
#include <vector>
#include <list>
#include <map>
#include "genomicregion.h"

using namespace std;
 
class Cnv{

public:
	Cnv(string const&, string const&, string const&, int, string const&, int, int, int);
	unsigned int GetStart() const;
	unsigned int GetStop() const;
	string GetOurId() const;
	bool CentromereSpanCheck() const;
	
private:
	string our_id_;	//our identifier "disorder_0001"
	string study_id_;	//study identifier "A#####"
	string chromosome_;
	string cnv_type_;	//string for Het Del, Homo Del, etc.
		//string case_or_control_;	//case or control
	string phenotype_;	//disorder
	int copy_num_;
	int case_or_control_;
	unsigned int start_, stop_;	//start and end points for each cnv
	unsigned int length_;	//length of CNV
	string size_category_;	//for annotating size
	bool major_chromosomal_abnormality_;
	bool spans_centromere_;
	
	//case control count variables
	int num_cases_overlap_with_same_, num_controls_overlap_with_same_; // count variables, +1 each time we overlap with a case or control
	int num_cases_overlap_with_diff_, num_controls_overlap_with_diff_;
	int num_cases_exact_, num_controls_exact_;
	
	//dgv2 overlap variables
	vector<GenomicRegion> dgv2_overlaps_same_;
	list<GenomicRegion> dgv2_coalesced_overlaps_same_;
	vector<GenomicRegion> dgv2_overlaps_diff_;
	list<GenomicRegion> dgv2_coalesced_overlaps_diff_;
	double dgv2_percent_overlapped_same_, dgv2_percent_overlapped_diff_;
	
	//other algorithm variables
	vector<GenomicRegion> other_algorithm_overlaps_same_;
	list<GenomicRegion> other_algorithm_coalesced_overlaps_same_;
	double other_algorithm_percent_overlapped_same_;
	
	//cumulative database variables
	vector<GenomicRegion> cumulativeDB_overlaps_exacttype_;
	list<GenomicRegion> cumulativeDB_coalesced_overlaps_exacttype_;
	vector<GenomicRegion> cumulativeDB_overlaps_simtype_;
	list<GenomicRegion> cumulativeDB_coalesced_overlaps_simtype_;
	vector<GenomicRegion> cumulativeDB_overlaps_difftype_;
	list<GenomicRegion> cumulativeDB_coalesced_overlaps_difftype_;
	double cumulativeDB_percent_overlapped_exacttype_, cumulativeDB_percent_overlapped_simtype_, cumulativeDB_percent_overlapped_difftype_;
	string which_cumulativeDB_overlaps_exacttype_, which_cumulativeDB_overlaps_simtype_, which_cumulativeDB_overlaps_difftype_;
	int num_exact_unique_;
	
	//cnvr variables
	vector<GenomicRegion> chop_overlaps_;
	list<GenomicRegion> chop_coalesced_overlaps_;
	vector<GenomicRegion> hapmap_overlaps_;
	list<GenomicRegion> hapmap_coalesced_overlaps_;
	string chop_overlaps_ids_, chop_percentages_, chop_frequencies_;
	string chop_total_percent_;
	string hapmap_overlaps_ids_, hapmap_percentages_, hapmap_frequencies_;
	string hapmap_total_percent_;
	
	//transcripts/genes maps
	multimap<string,string> regions_overlapped_;	//this multimap is <gene,transcript>
	map<string,string> percent_regions_overlapped_;	//this map is <transcript,percent_overlap>
	//ccds
	string ccds_gene_string_, ccds_overlap_string_;
	int ccds_region_count_;
	//omim
	string omim_gene_string_, omim_overlap_string_;
	int omim_region_count_;
	//gencode
	string gencode_gene_string_, gencode_overlap_string_;
	int gencode_region_count_;
	//pathway genes
	string pathway_gene_string_, pathway_gene_overlap_string_;
	
	vector<GenomicRegion> pseudo_overlaps_;
	list<GenomicRegion> pseudo_coalesced_overlaps_;
	vector<GenomicRegion> parent_overlaps_;
	list<GenomicRegion> parent_coalesced_overlaps_;
	string pseudo_overlaps_ids_, pseudo_percentages_;
	string pseudo_total_percent_;
	string parent_overlaps_ids_, parent_percentages_;
	string parent_total_percent_;
	
	
	
	friend class CnvVector;
};

class CnvVector{

public:
	//void AddCnv(Cnv const& a_cnv);
	void CalculateCaseControl(double case_control_overlap_threshold);
	void CalculateDGV2Overlap(vector<GenomicRegion>& dgv2_vector);
	void CompareToAnotherAlgorithm(vector<GenomicRegion>& other_algorithm_cnvs);
	
	bool ReadCnvFile(string const& cnv_filename);
	bool WriteSelfToFile(string const& out_filename);
	
	void CheckRegionOfInterest(vector<GenomicRegion>& region_of_interest_vector, string const& what_region, int overlap_type);
	void SetWriteRegions(string const& what_regions);
	void CheckCumulativeDatabaseOverlap(vector<GenomicRegion>& cumulative_cnvs);
	void CheckCnvr(vector<GenomicRegion>& some_blocks, const string& which_cnvr);
	
	
private:
	vector<Cnv> cnv_vector_;
	string FormatOverlap(const double percent_overlapped) const;
	string GetIndividualCompareAlgorithmPercents(Cnv our_cnv) const;
	int OutputLongString(string const& my_long_string, ofstream & outfile);
	string algorithm_compared_to_;
	int num_columns_for_ccds_;
	int num_columns_for_gencode_;
	int num_columns_for_omim_;
};

#endif