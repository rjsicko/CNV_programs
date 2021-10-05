#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <map>
#include <algorithm>
#include "overlap.h"
#include "genomicregion.h"
#include "cnv.h"

using namespace std;


				/*private:
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
					int size_category_;	//for annotating size
					int num_cases_overlap_with_same_, num_controls_overlap_with_same_; // count variables, +1 each time we overlap with a case or control
					int num_cases_overlap_with_diff_, num_controls_overlap_with_diff_;
					int num_cases_exact_, num_controls_exact_;
					
					friend class CnvVector;
				};*/

			//one cnv
			
	//Cnv(string const&,       string const&,          string const&,           int,                 string const&,            int,       int,      int);
Cnv::Cnv(string const& our_id, string const& study_id, string const& phenotype, int case_or_control, string const& chromosome, int start, int stop, int copy_num)
{
	//assign all of our private variables based on input from constructor
	our_id_ = our_id;
	study_id_ = study_id;
	chromosome_ = chromosome;
	copy_num_ = copy_num;
	start_ = start;
	stop_ = stop;
	phenotype_ = phenotype;
	case_or_control_ = case_or_control;
	//define cnv_type_ variable based on copy_num
	if(copy_num == 0)
		cnv_type_ = "Homo Del";
	else if(copy_num == 1)
		cnv_type_ = "Het Del";
	else if(copy_num == 3)
		cnv_type_ = "Dupl";
	else if(copy_num == 4)
		cnv_type_ = "Tripl";
	
	//calculate our length
	length_ = stop - start + 1;
	//define size_category_ variable based on size			
	if(length_ < 1000)
		size_category_ = "1_<1000-bp";
	else if(length_ >= 1000 && length_ < 5000)
		size_category_ = "2_[1000-5000)-bp";
	else if(length_ >= 5000 && length_ < 25000)
		size_category_ = "3_[5-25)-Kb";
	else if(length_ >= 25000 && length_ < 100000)
		size_category_ = "4_[25-100)-Kb";
	else if(length_ >= 100000 && length_ < 500000)
		size_category_ = "5_[100-500)-Kb";
	else if(length_ >= 500000 && length_ < 1000000)
		size_category_ = "6_[0.5-1)-Mb";
	else if(length_ >= 100000)
		size_category_ = "7_=1-Mb";	
	
	num_exact_unique_ = 0;	
	//set all overlap numbers to 0 
	//these are our overlap > threshold counts
	num_cases_overlap_with_same_ = 0;
	num_controls_overlap_with_same_ = 0;
	num_cases_overlap_with_diff_ = 0;
	num_controls_overlap_with_diff_ = 0;
	
	//these are our overlap exact counts
	num_cases_exact_ = 0;
	num_controls_exact_ = 0;
	
	//some basic annotation
	if(length_ >= 5000000)	//major chr abn
		major_chromosomal_abnormality_ = true;
	else 
		major_chromosomal_abnormality_ = false;
	
	if(CentromereSpanCheck())
		spans_centromere_ = true;
	else
		spans_centromere_ = false;
}


unsigned int Cnv::GetStart() const
{
	return start_;
}

unsigned int Cnv::GetStop() const
{
	return stop_;
}
string Cnv::GetOurId() const
{
	return our_id_;
}

bool Cnv::CentromereSpanCheck() const
{
	unsigned int cent_start, cent_stop;
	//assign the cent_start and cent_stop variables depending on which chromosome our cnv is on
	if(chromosome_ == "1") cent_start = 125600000;
	if(chromosome_ == "2") cent_start = 93300000;
	if(chromosome_ == "3") cent_start = 91300000;
	if(chromosome_ == "4") cent_start = 50000000;
	if(chromosome_ == "5") cent_start = 48000000;
	if(chromosome_ == "6") cent_start = 61000000;
	if(chromosome_ == "7") cent_start = 60200000;
	if(chromosome_ == "8") cent_start = 45400000;
	if(chromosome_ == "9") cent_start = 49000000;
	if(chromosome_ == "10") cent_start = 40000000;
	if(chromosome_ == "11") cent_start = 54000000;
	if(chromosome_ == "12") cent_start = 35500000;
	if(chromosome_ == "13") cent_start = 17500000;
	if(chromosome_ == "14") cent_start = 17500000;
	if(chromosome_ == "15") cent_start = 18800000;
	if(chromosome_ == "16") cent_start = 36500000;
	if(chromosome_ == "17") cent_start = 23800000;
	if(chromosome_ == "18") cent_start = 17000000;
	if(chromosome_ == "19") cent_start = 26500000;
	if(chromosome_ == "20") cent_start = 27500000;
	if(chromosome_ == "21") cent_start = 13000000;
	if(chromosome_ == "22") cent_start = 14600000;
	if(chromosome_ == "X") cent_start = 60500000;
	if(chromosome_ == "Y") cent_start = 12300000;
	
	cent_stop = cent_start + 1;
	
	if( CalcOverlapCNV(start_, cent_start, stop_, cent_stop) > 0)	//if we overlapped the centromere at all
		return true;
	else
		return false;

}

	/*void CnvVector::AddCnv(Cnv const& a_cnv)
	{
		cnv_vector_.push_back(a_cnv);
	}*/

void CnvVector::CalculateCaseControl(double case_control_overlap_threshold)
{		
	//for each cnv in our list
	for (vector<Cnv>::iterator me = cnv_vector_.begin();me != cnv_vector_.end();me++)
		//go through each OTHER cnv and compare ourself to them
		for(vector<Cnv>::iterator them = cnv_vector_.begin();them != cnv_vector_.end();them++)
		{
			if(me->our_id_ == them->our_id_)	//comparing the same CNV so skip, don't count ourself
				continue;
			else if(me->chromosome_ == them->chromosome_)
			{	//the same chromosome?
				bool phenotype_match_flag = (me->phenotype_ == them->phenotype_);
				bool copy_number_match_flag = (me->copy_num_ == them->copy_num_);
				bool real_control_flag = (me->case_or_control_ == 2);
				
				if(me->start_ == them->start_ && me->stop_ == them->stop_
				   && copy_number_match_flag)	
				//exact match, so only worry about if the types of CNVs match
				{
					//first deal with if we are a real control
					//if so, we only have two options overlap a case or overlap another real control
					
					//real control matched to another real control
					if(real_control_flag && them->case_or_control_ == 2){
						me->num_controls_exact_++;
						me->num_controls_overlap_with_same_++;	//increment the above-threshold count to since 100% overlap will be above the threshold
					}
					
					//otherwise, im a real control and i overlap a case
					else if(real_control_flag){
						me->num_cases_exact_++;
						me->num_cases_overlap_with_same_++;	//increment the above-threshold count to since 100% overlap will be above the threshold
					}
					
					//on to the case situations
					//im a case, overlap a real control
					else if( (!real_control_flag) && (them->case_or_control_ == 2) ){
						me->num_controls_exact_++;
						me->num_controls_overlap_with_same_++;	//increment the above-threshold count to since 100% overlap will be above the threshold
					}
						
					//im a case, overlap another disorder (pseudo-control)
					else if( (!real_control_flag) && (!phenotype_match_flag) ){
						me->num_controls_exact_++;
						me->num_controls_overlap_with_same_++;	//increment the above-threshold count to since 100% overlap will be above the threshold
					}
					
					//im a case, overlap another case
					else if( (!real_control_flag) && phenotype_match_flag){
						me->num_cases_exact_++;
						me->num_cases_overlap_with_same_++;	//increment the above-threshold count too since 100% overlap will be above the threshold
					}
					
					//otherwise we forgot to account for all possibilities, lets spit out a negative -999 for debugging
					else{
						me->num_cases_exact_ = -999;
						me->num_controls_exact_ = -999;
					}
				}	//finished with the exact matches
				
					/*if(me->our_id_ == "PUV_1480")
					{
							std::cout << "Me: " << me->our_id_ << endl;
							std::cout << "Them: " << them->our_id_ << endl;
							std::cout << "Calculated overlap: " << std::setprecision(3)
									 << double(CalcOverlapCNV(me->start_, them->start_, me->stop_, them->stop_)) << endl;
					 }*/
					
				//now take case of overlap > threshold situations	
				else if(CalcOverlapCNV(me->start_, them->start_, me->stop_, them->stop_) > case_control_overlap_threshold)
				{	
					//specific case for debugging
						/*if(me->our_id_ == "PUV_1480")
						{
							std::cout << "Me: " << me->our_id_ << endl;
							std::cout << "Them: " << them->our_id_ << endl;
							std::cout << "Calculated overlap: "
								 << double(CalcOverlapCNV(me->start_, them->start_, me->stop_, them->stop_)) << endl;
						}*/
					
					//lets take care of the same type of cnv first
					if(copy_number_match_flag
					   || (me->copy_num_ == them->copy_num_ + 1)
					   || (me->copy_num_ == them->copy_num_ - 1) )
					{
						
							//if(me->chromosome_ == "7" && me->start_ == 70421644 && me->stop_ == 70423234)
							//std::cout << "copy_number_match_flag: " << copy_number_match_flag << endl;
						//first deal with if we are a real control
						//if so, we only have two options overlap a case or overlap another real control
						
						//real control matched to another real control
						if(real_control_flag && them->case_or_control_ == 2)
							me->num_controls_overlap_with_same_++;
						
						//otherwise, im a real control and i overlap a case
						else if(real_control_flag)
							me->num_cases_overlap_with_same_++;
						
						//on to the case situations
						//im a case, overlap a real control
						else if( (!real_control_flag) && (them->case_or_control_ == 2) )
							me->num_controls_overlap_with_same_++;
							
						//im a case, overlap another disorder (pseudo-control)
						else if( (!real_control_flag) && (!phenotype_match_flag) )
							me->num_controls_overlap_with_same_++;
						
						//im a case, overlap another case
						else if( (!real_control_flag) && phenotype_match_flag)
						{
								/*if(me->our_id_ == "PUV_1480")
								{
									std::cout << "got here\n";
								}*/
							me->num_cases_overlap_with_same_++;
						}
						
						//otherwise we forgot to account for all possibilities, lets spit out a negative -5000 for debugging
						else{
							me->num_controls_overlap_with_same_ = -5000;
							me->num_cases_overlap_with_same_ = -5000;
						}
					}
					
					//then, lets take care of different type of cnv
					else
					{		
						//first deal with if we are a real control
						//if so, we only have two options overlap a case or overlap another real control
						
						//real control matched to another real control
						if(real_control_flag && them->case_or_control_ == 2)
							me->num_controls_overlap_with_diff_++;
						
						//otherwise, im a real control and i overlap a case
						else if(real_control_flag)
							me->num_cases_overlap_with_diff_++;
						
						//on to the case situations
						//im a case, overlap a real control
						else if( (!real_control_flag) && (them->case_or_control_ == 2) )
							me->num_controls_overlap_with_diff_++;
							
						//im a case, overlap another disorder (pseudo-control)
						else if( (!real_control_flag) && (!phenotype_match_flag) )
							me->num_controls_overlap_with_diff_++;
						
						//im a case, overlap another case
						else if( (!real_control_flag) && phenotype_match_flag)
							me->num_cases_overlap_with_diff_++;
						
						//otherwise we forgot to account for all possibilities, lets spit out a negative -2222 for debugging
						else{
							me->num_cases_overlap_with_diff_ = -2222;
							me->num_controls_overlap_with_diff_ = -2222;
						}
					}
				}
			}
		}
}

void CnvVector::CalculateDGV2Overlap(vector<GenomicRegion>& dgv2_vector)
{
	//for each cnv in our list
	for (vector<Cnv>::iterator me = cnv_vector_.begin();me != cnv_vector_.end();me++)
	{
		//go through each OTHER cnv and compare ourself to them
		for(vector<GenomicRegion>::iterator them = dgv2_vector.begin();them != dgv2_vector.end();them++)
		{
			//only bother calculating overlap if we are on the same chromosome
			if(me->chromosome_ == them->GetChromosome())
			{
				double overlap_cnv = CalcOverlapCNV(me->start_, them->GetStart(), me->stop_, them->GetStop());
				if(overlap_cnv>0)
				{	//If there is some overlap with our CNV, then lets keep track of it
					
					unsigned int dummy_start, dummy_stop;
					dummy_start = max (me->start_,them->GetStart());
					dummy_stop = min (me->stop_,them->GetStop());
					
					GenomicRegion my_overlap(me->chromosome_,dummy_start,dummy_stop);
					
					//if we match del to del or dup to dup
					if(me->copy_num_ == them->GetCopyNum() || me->copy_num_ == them->GetCopyNum() - 1 || me->copy_num_ == them->GetCopyNum() + 1)
					{
						me->dgv2_overlaps_same_.push_back(my_overlap);
						me->dgv2_percent_overlapped_same_ = me->dgv2_percent_overlapped_same_ + overlap_cnv;
					}
					else
					{
						me->dgv2_overlaps_diff_.push_back(my_overlap);
						me->dgv2_percent_overlapped_diff_ = me->dgv2_percent_overlapped_diff_ + overlap_cnv;
					}
				}
			}
		}
		//for same type overlaps
		//copy our overlap_bases vector to our coalesced_overlaps list
		copy(me->dgv2_overlaps_same_.begin(), me->dgv2_overlaps_same_.end(), back_inserter(me->dgv2_coalesced_overlaps_same_));
		
			//used for debugging
			/*if(me->our_id_ == "HPPE_Control_0003")
			{
				cout << me->dgv2_coalesced_overlaps_same_.size() << endl
					 << CalcTotalOverlap(*me,me->dgv2_coalesced_overlaps_same_) << endl;
			}*/
		
		//pass our coalesced_overlaps vector to the coalesce function
		CoalesceOverlaps(me->dgv2_coalesced_overlaps_same_);	
		
			//used for debugging
			/*if(me->our_id_ == "HPPE_Control_0003")
			{
				cout << me->dgv2_coalesced_overlaps_same_.size() << endl
					 << "overlap in cnv.cpp: " << CalcTotalOverlap(*me,me->dgv2_coalesced_overlaps_same_) << endl;
			}*/
		
		//calculate our percent overlapped with the coalesced overlaps
		me->dgv2_percent_overlapped_same_ = CalcTotalOverlap(*me,me->dgv2_coalesced_overlaps_same_);
		
		
		//for diff type overlaps
		//copy our overlap_bases vector to our coalesced_overlaps list
		copy(me->dgv2_overlaps_diff_.begin(), me->dgv2_overlaps_diff_.end(), back_inserter(me->dgv2_coalesced_overlaps_diff_));
		
		//pass our coalesced_overlaps vector to the coalesce function
		CoalesceOverlaps(me->dgv2_coalesced_overlaps_diff_);	
		
		//calculate our percent overlapped with the coalesced overlaps
		me->dgv2_percent_overlapped_diff_ = CalcTotalOverlap(*me,me->dgv2_coalesced_overlaps_diff_);
	}
		
}


				//vector<GenomicRegion> other_algorithm_overlaps_same_;
				//list<GenomicRegion> other_algorithm_coalesced_overlaps_same_;
void CnvVector::CompareToAnotherAlgorithm(vector<GenomicRegion>& other_algorithm_cnvs)
{
	bool first_run = true;
	//for each cnv in our list
	for (vector<Cnv>::iterator me = cnv_vector_.begin();me != cnv_vector_.end();me++)
	{
		//go through each OTHER cnv and compare ourself to them
		for(vector<GenomicRegion>::iterator them = other_algorithm_cnvs.begin();them != other_algorithm_cnvs.end();them++)
		{
			if(first_run)
			{
				algorithm_compared_to_ = them->GetDescription();
				first_run = false;
			}
			
			//only bother calculating overlap if we are on the same chromosome
			//and we are the same subject (since we're comparing algorithms)
			if(me->chromosome_ == them->GetChromosome() && me->study_id_ == them->GetUniqueId2())
			{
				if(me->copy_num_ == them->GetCopyNum()) //only checking exact calls for now
				{
					double overlap_cnv = CalcOverlapCNV(me->start_, them->GetStart(), me->stop_, them->GetStop());
					if(overlap_cnv>0)
					{	//If there is some overlap with our CNV, then lets keep track of it
					
						unsigned int dummy_start, dummy_stop;
						dummy_start = max (me->start_,them->GetStart());
						dummy_stop = min (me->stop_,them->GetStop());
					
						GenomicRegion my_overlap(me->chromosome_,dummy_start,dummy_stop);
						
						me->other_algorithm_overlaps_same_.push_back(my_overlap);
						me->other_algorithm_percent_overlapped_same_ = me->other_algorithm_percent_overlapped_same_ + overlap_cnv;

					}
				}
			}
		}
		//copy our overlap_bases vector to our coalesced_overlaps list
		copy(me->other_algorithm_overlaps_same_.begin(), me->other_algorithm_overlaps_same_.end(), back_inserter(me->other_algorithm_coalesced_overlaps_same_));

		//pass our coalesced_overlaps vector to the coalesce function
		CoalesceOverlaps(me->other_algorithm_coalesced_overlaps_same_);	
		
		//calculate our percent overlapped with the coalesced overlaps
		me->other_algorithm_percent_overlapped_same_ = CalcTotalOverlap(*me,me->other_algorithm_coalesced_overlaps_same_);	
	}
}
	

//our_id	sampleID	disorder	case/control	chr	start	end	copy_num	
bool CnvVector::ReadCnvFile(string const& cnv_filename)
{
	ifstream infile(cnv_filename.c_str());
	if(!infile.is_open())
		return false;	//file failed to open

	bool first_line = true;
	string line;
	while (getline(infile, line))	//while there are still lines left to read
	{
		string our_id, study_id, phenotype, chromosome;
		int start, stop, copy_num, case_or_control;
		
		istringstream ss(line);
		ss >> our_id >> study_id >> phenotype >> case_or_control
		   >> chromosome >> start >> stop >> copy_num;
		
		if(first_line)	//ignore header line
			first_line=false;
		
		else{	//we already skipped header, start adding the CNVs
			//Cnv(string const& our_id, string const& study_id, string const& phenotype, int case_or_control, string const& chromosome, int start, int stop, int copy_num)
			//Cnv::Cnv(string const& our_id, string const& study_id, string const& phenotype, int case_or_control, string const& chromosome, int start, int stop, int copy_num)
			Cnv dummy_cnv(our_id, study_id, phenotype, case_or_control, chromosome, start, stop, copy_num);
			cnv_vector_.push_back(dummy_cnv);
		}	
	}
	infile.close();
	
	num_columns_for_gencode_ = 1;
	num_columns_for_ccds_ = 1;
	num_columns_for_omim_ = 1;
	return true;
}

string CnvVector::FormatOverlap(const double percent_overlapped) const
{
	if(percent_overlapped != 0)
	{
		double olap_percent = percent_overlapped * 100;
		ostringstream strs;
		strs << fixed << setprecision (2) << olap_percent;
		string olap = strs.str();	
		if(olap == "100.00") olap = "100";
		return olap;
	}
	else	
		return "NA";
}

string CnvVector::GetIndividualCompareAlgorithmPercents(Cnv our_cnv) const
{
	string overlap_percent;
	for(vector<GenomicRegion>::iterator my_overlap = our_cnv.other_algorithm_overlaps_same_.begin();
				my_overlap != our_cnv.other_algorithm_overlaps_same_.end();
				my_overlap++)	
	{	
				
		if(overlap_percent!= "")
			overlap_percent = overlap_percent + ";";
		
			
		//convert overlap to a string
		double temp_overlap = CalcOverlapCNV(our_cnv.start_, my_overlap->GetStart(), our_cnv.stop_, my_overlap->GetStop());
		
		overlap_percent = overlap_percent + FormatOverlap(temp_overlap);
	}
	if(overlap_percent == "") overlap_percent = "NA";
	return overlap_percent;
}


int CnvVector::OutputLongString(string const& my_long_string, ofstream & outfile)
{
	int column_count = 1;
	if(my_long_string.size() > 32000)
	{
		string first_part = my_long_string.substr (0, 32000);
		outfile << first_part << "\t";
		
		if(my_long_string.size() >  32000)		//going to need to split output into 2 columns
		{
			column_count++;
			string first_sub = my_long_string.substr (32000, 32000);
			outfile << first_sub << "\t";
		}
		if(my_long_string.size() > 64000)		//going to need to split output into 3 columns
		{
			column_count++;
			string second_sub = my_long_string.substr (64000, 32000);
			outfile << second_sub << "\t";
		}
		if(my_long_string.size() > 96000)		//going to need to split output into 4 columns
		{
			column_count++;
			string third_sub = my_long_string.substr (96000, 32000);
			outfile << third_sub << "\t";
		}
		if(my_long_string.size() > 128000) 	//going to need to split output into 5 columns
		{
			column_count++;
			string fourth_sub = my_long_string.substr (128000, 32000);
			outfile << fourth_sub << "\t";		
		}
	}
	else
		outfile << my_long_string << "\t";
	return column_count;
}
	

bool CnvVector::WriteSelfToFile(string const& out_filename)
{
	ofstream outfile(out_filename.c_str());
	if(!outfile)
		 return false;
	
	//print header line
	outfile << "OurID\t"
			<< "StudyID\t"
			<< "Phenotype\t"
			<< "Case/Control\t"
			<< "Chromosome\t"
			<< "Start\t"
			<< "Stop\t"
			<< "Copy_Number\t"
			<< "Type\t"
			<< "Size\t"
			<< "Size Category\t"
			<< "Overlap Count Cases(Same CNV Type) > 0.50 Threshold\t"
			<< "Overlap Count Controls(Same CNV Type) > 0.50 Threshold\t"
			<< "Overlap Count Cases(Different CNV Type) > 0.50 Threshold\t"
			<< "Overlap Count Controls(Different CNV Type) > 0.50 Threshold\t"
			<< "Overlap Count Cases (Exact)\t"
			<< "Overlap Count Controls (Exact)\t"
			<< "Total % Our CNV With DGV2 Overlap(Same Type)\t"
			<< "Total % Our CNV With DGV2 Overlap(Different Type)\t"
			<< "% Our CNV Overlapped by Call(s) in " << algorithm_compared_to_ << "\t"
			<< "Total % Our CNV Overlapped by " << algorithm_compared_to_ << " Call(s)\t"
			<< "Major Chromosomal Abnormality\t"
			<< "Spans Centromere\t"
			<< "Number Exact Unique Times Seen\t"
			<< "Total % Our CNV With Cumulative Database Overlap(Exact Type)\t"
			<< "Cumulative Database CNVs Overlapped With(Exact Type)\t"			
			<< "Total % Our CNV With Cumulative Database Overlap(Similar Type)\t"
			<< "Cumulative Database CNVs Overlapped With(Similar Type)\t"
			<< "Total % Our CNV With Cumulative Database Overlap(Different Type)\t"
			<< "Cumulative Database CNVs Overlapped With(Different Type)\t"
			<< "HapMap Ids\t"
			<< "HapMap Est. Freq\t"
			<< "Percent Each HapMap\t"
			<< "Total % Our CNV Overlapped by HapMap\t"
			<< "CHOP Ids\t"
			<< "CHOP Freq\t"
			<< "Percent Each CHOP\t"
			<< "Total % Our CNV Overlapped by CHOP\t"
			<< "Gencode Transcript(s) Overlapped\t";
	if(num_columns_for_gencode_ > 1)
		for(int i = 1; i < num_columns_for_gencode_; i++)
			outfile << "Gencode Transcripts Overflow (" << i+1 << ")\t";
	outfile << "Transcript Coverage\t"
			<< "Number Transcripts Covered\t"
			<< "CCDS Genes Overlapped\t";
	if(num_columns_for_ccds_ > 1)
		for(int i = 1; i < num_columns_for_ccds_; i++)
			outfile << "CCDS Overflow (" << i+1 << ")\t";
	outfile << "Number CCDS Genes Covered\t"
			<< "Number CCDS Genes Unique\t"
			<< "OMIM Genes Overlapped\t";
	if(num_columns_for_omim_ > 1)
		for(int i = 1; i < num_columns_for_omim_; i++)
			outfile << "OMIM Genes Overflow (" << i+1 << ")\t";
	outfile << "Percent OMIM Coverage\t"
			<< "Number OMIM Overlapped\t"
			<< "Pathway Gene Name(s)\t"
			<< "Pathway Percent Overlap\t"
			<< "Pseudogene ID(s)\t"
			<< "Percent Each Pseudogene\t"
			<< "Total % Our CNV Overlapped by Pseudogene(s)\t"
			<< "Parentgene ID(s)\t"
			<< "Percent Each Parentgene\t"
			<< "Total % Our CNV Overlapped by Parentgene(s)\t"
			<< "DGV Coordinates\t"
			<< "Illumina IGV Coordinates\t"
			<< "UCSC Coordinates\n";
	

	//for (vector<Cnv>::iterator me = cnv_vector_.begin();me != cnv_vector_.end();me++)		
	for (vector<Cnv>::iterator me = cnv_vector_.begin();me != cnv_vector_.end();me++)
	{
		outfile << me->our_id_ << "\t"
				<< me->study_id_ << "\t" 
				<< me->phenotype_ << "\t"
				<< me->case_or_control_ << "\t" 
				<< me->chromosome_ << "\t" 
				<< me->start_ << "\t" 
				<< me->stop_ << "\t" 
				<< me->copy_num_ << "\t"
				<< me->cnv_type_ << "\t"
				<< me->length_ << "\t"
				<< me->size_category_ << "\t"
				<< me->num_cases_overlap_with_same_ << "\t"
				<< me->num_controls_overlap_with_same_ << "\t"
				<< me->num_cases_overlap_with_diff_ << "\t"
				<< me->num_controls_overlap_with_diff_ << "\t"
				<< me->num_cases_exact_ << "\t"
				<< me->num_controls_exact_ << "\t"
				<< FormatOverlap(me->dgv2_percent_overlapped_same_) << "\t"
				<< FormatOverlap(me->dgv2_percent_overlapped_diff_) << "\t"
				<< GetIndividualCompareAlgorithmPercents(*me) << "\t"
				<< FormatOverlap(me->other_algorithm_percent_overlapped_same_) << "\t"
				<< me->major_chromosomal_abnormality_ << "\t"
				<< me->spans_centromere_ << "\t"
				<< me->num_exact_unique_ << "\t"
				<< FormatOverlap(me->cumulativeDB_percent_overlapped_exacttype_) << "\t";
		if(me->which_cumulativeDB_overlaps_exacttype_ == "")
			outfile << "NA" << "\t";
		else
			outfile << me->which_cumulativeDB_overlaps_exacttype_ << "\t";
		outfile << FormatOverlap(me->cumulativeDB_percent_overlapped_simtype_) << "\t";
		if(me->which_cumulativeDB_overlaps_simtype_ == "")
			outfile << "NA" << "\t";
		else
			outfile	<< me->which_cumulativeDB_overlaps_simtype_ << "\t";
		outfile	<< FormatOverlap(me->cumulativeDB_percent_overlapped_difftype_) << "\t";
		if(me->which_cumulativeDB_overlaps_difftype_ == "")
			outfile << "NA" << "\t";
		else
			outfile << me->which_cumulativeDB_overlaps_difftype_ << "\t";
		outfile << me->hapmap_overlaps_ids_ << "\t"
				<< me->hapmap_frequencies_ << "\t"
				<< me->hapmap_percentages_ << "\t"
				<< me->hapmap_total_percent_ << "\t"
				<< me->chop_overlaps_ids_ << "\t"
				<< me->chop_frequencies_ << "\t"
				<< me->chop_percentages_ << "\t"
				<< me->chop_total_percent_ << "\t";				
				
		//output gencode stuff
		if(num_columns_for_gencode_ > 1)
		{
					/*if(me->our_id_ == "HTXY_B_3562")
						cout << "HTXY_B_3562: and got to num_columns > 1.\n";*/
			int column_count = OutputLongString(me->gencode_gene_string_, outfile);
					/*if(me->our_id_ == "HTXY_B_3562")
						cout << "got back from outputlongstring. column_count =" << column_count << " and num columns for gencode = " << num_columns_for_gencode_ << endl;*/
			while(column_count < num_columns_for_gencode_)
			{
				outfile << "NA\t";
				column_count++;
			}
		
		}
		else
			outfile << me->gencode_gene_string_ << "\t";
		outfile << me->gencode_overlap_string_ << "\t"
				<< me->gencode_region_count_ << "\t";

		//output ccds stuff
		if(num_columns_for_ccds_ > 1)
		{
			int column_count = OutputLongString(me->ccds_gene_string_, outfile);
			while(column_count < num_columns_for_ccds_)
			{
				outfile << "NA\t";
				column_count++;
			}
		}
		else
			outfile << me->ccds_gene_string_ << "\t";	
		outfile << "CalculateInExcel!\t"
				<< me->ccds_region_count_ << "\t";
		
		//output omim
		if(num_columns_for_omim_ > 1)
		{
			int column_count = OutputLongString(me->omim_gene_string_, outfile);
			while(column_count < num_columns_for_omim_)
			{
				outfile << "NA\t";
				column_count++;
			}
		}
		else
			outfile << me->omim_gene_string_ << "\t";
		outfile << me->omim_overlap_string_ << "\t"
				<< me->omim_region_count_ << "\t"
				<< me->pathway_gene_string_ << "\t"
				<< me->pathway_gene_overlap_string_ << "\t"
				<< me->pseudo_overlaps_ids_ << "\t"
				<< me->pseudo_percentages_ << "\t"
				<< me->pseudo_total_percent_ << "\t"
				<< me->parent_overlaps_ids_ << "\t"
				<< me->parent_percentages_ << "\t"
				<< me->parent_total_percent_ << "\t";		
		outfile << "chr" << me->chromosome_ << ":" << (me->start_-10000) << ".." << (me->stop_+10000) << "\t"
				<< (me->start_-200000)<< " - " << (me->stop_+200000) << "\t"
				<< "chr" << me->chromosome_ << ":" << (me->start_-10000) << "-" << (me->stop_+10000) << "\n";
	}		
		
	outfile.close();
	return true;
}


void CnvVector::SetWriteRegions(string const& what_regions)
{
	for (vector<Cnv>::iterator me = cnv_vector_.begin();me != cnv_vector_.end();me++)
	{
		int region_count=0;
		string gene_string="", overlap_string="";
		multimap<string,string>::iterator regions_list = me->regions_overlapped_.begin();
		map<string,string>::iterator percent_regions_overlapped;
			
		for(;regions_list != me->regions_overlapped_.end();regions_list = me->regions_overlapped_.upper_bound(regions_list->first))	//for all genes
		{
			region_count++;
			if(gene_string!="")	//if our gene list is not empty
				gene_string = gene_string + ";";	//add ; at end
			gene_string = gene_string + regions_list->first//add gene name to our gene_string
						+ "[";							//and get ready for transcript ids	
			if(overlap_string!="")
				overlap_string = overlap_string + ";";
								
			//temp iterator for the next for
			multimap<string,string>::iterator temp_region_list = regions_list;
			//now i=0, keep going until i gets to the number of transcripts for this gene
			for(int i = 0;i < me->regions_overlapped_.count(regions_list->first);i++,temp_region_list++)
			{
				if(i==0)	//first transcript
					overlap_string = overlap_string + "[";
				else	//not our first transcript
				{									
					gene_string = gene_string + ",";	//append comma
					overlap_string = overlap_string + ",";
				}
				gene_string = gene_string + temp_region_list->second;	//append transcript name
				percent_regions_overlapped = me->percent_regions_overlapped_.find(temp_region_list->second);	//get iterator for this transcripts map	
				overlap_string = overlap_string + percent_regions_overlapped->second;
			}
			//must be done with transcripts for this gene, now lets close the )
			gene_string = gene_string + "]";
			overlap_string = overlap_string + "]";
	
		}
		if(gene_string == "")
			gene_string = "NA";
		if(overlap_string == "")
			overlap_string = "NA";
		
		if(what_regions == "ccds_genes")
		{
			me->ccds_gene_string_ = gene_string;
			me->ccds_overlap_string_ = overlap_string;
			me->ccds_region_count_ = region_count;
			//checks below for exceeding max string length in excel
			if(gene_string.size() > 32000 && num_columns_for_ccds_ < 2)		//going to need to split output into 2 columns
				num_columns_for_ccds_ = 2;
			if(gene_string.size() > 64000 && num_columns_for_ccds_ < 3)		//going to need to split output into 3 columns
				num_columns_for_ccds_ = 3;
			if(gene_string.size() > 96000 && num_columns_for_ccds_ < 4)		//going to need to split output into 4 columns
				num_columns_for_ccds_ = 4;
			if(gene_string.size() > 128000 && num_columns_for_ccds_ < 5) 	//going to need to split output into 5 columns
				num_columns_for_ccds_ = 5;			
		}		
		else if(what_regions == "omim_genes")
		{
			me->omim_gene_string_ = gene_string;
			me->omim_overlap_string_ = overlap_string;
			me->omim_region_count_ = region_count;
			//checks below for exceeding max string length in excel
			if(gene_string.size() > 32000 && num_columns_for_omim_ < 2)		//going to need to split output into 2 columns
				num_columns_for_omim_ = 2;
			if(gene_string.size() > 64000 && num_columns_for_omim_ < 3)		//going to need to split output into 3 columns
				num_columns_for_omim_ = 3;
			if(gene_string.size() > 96000 && num_columns_for_omim_ < 4)		//going to need to split output into 4 columns
				num_columns_for_omim_ = 4;
			if(gene_string.size() > 128000 && num_columns_for_omim_ < 5) 	//going to need to split output into 5 columns
				num_columns_for_omim_ = 5;						
		}
		else if(what_regions == "gencode_genes")
		{
			/*if(me->our_id_ == "HTXY_B_3562")
				cout << gene_string.size() << endl;*/
				
			me->gencode_gene_string_ = gene_string;
			me->gencode_overlap_string_ = overlap_string;
			me->gencode_region_count_ = region_count;
			//checks below for exceeding max string length in excel
			if(gene_string.size() > 32000 && num_columns_for_gencode_ < 2)		//going to need to split output into 2 columns
				num_columns_for_gencode_ = 2;
			if(gene_string.size() > 64000 && num_columns_for_gencode_ < 3)		//going to need to split output into 3 columns
				num_columns_for_gencode_ = 3;
			if(gene_string.size() > 96000 && num_columns_for_gencode_ < 4)		//going to need to split output into 4 columns
				num_columns_for_gencode_ = 4;
			if(gene_string.size() > 128000 && num_columns_for_gencode_ < 5) 	//going to need to split output into 5 columns
				num_columns_for_gencode_ = 5;
		}
		else if(what_regions == "pathway_genes")
		{
			me->pathway_gene_string_ = gene_string;
			me->pathway_gene_overlap_string_ = overlap_string;
		}
	}
}

void CnvVector::CheckRegionOfInterest(vector<GenomicRegion>& region_of_interest_vector, string const& what_region, int overlap_type)
{
	
	for (vector<Cnv>::iterator me = cnv_vector_.begin();me != cnv_vector_.end();me++)
	{
		//for all Cnv clear their maps incase they were just used by a different comparison
		me->regions_overlapped_.clear();
		me->percent_regions_overlapped_.clear();
		for(vector<GenomicRegion>::iterator them = region_of_interest_vector.begin();them != region_of_interest_vector.end();them++)
		{
			//only bother calculating overlap if we are on the same chromosome
			if(me->chromosome_ == them->GetChromosome())
			{	
				double overlap;
				if(overlap_type == 1)	overlap = CalcOverlapUnion(me->start_, them->GetStart(), me->stop_, them->GetStop());
				if(overlap_type == 2)	overlap = CalcOverlapRegion(me->start_, them->GetStart(), me->stop_, them->GetStop());
				if(overlap_type == 3)	overlap = CalcOverlapCNV(me->start_, them->GetStart(), me->stop_, them->GetStop());					
									
				if(overlap > 0)
				{
						string olap = FormatOverlap(overlap);
						
						//insert gene and transcript. if gene's there it will just insert another transcript
						me->regions_overlapped_.insert( pair<string, string>( them->GetUniqueId(),them->GetUniqueId2() ) );
						//now insert overlap percentage for that transcript
						me->percent_regions_overlapped_.insert(pair<string, string>(them->GetUniqueId2(),olap));
						//cout << "I just finished inserting pairs in the maps" << endl;
				
				}
			}
		}
	}
	SetWriteRegions(what_region);
}


void CnvVector::CheckCumulativeDatabaseOverlap(vector<GenomicRegion>& cumulative_cnvs)
{
	//for each cnv in our list
	for (vector<Cnv>::iterator me = cnv_vector_.begin();me != cnv_vector_.end();me++)
	{
		//go through each OTHER cnv and compare ourself to them
		for(vector<GenomicRegion>::iterator them = cumulative_cnvs.begin();them != cumulative_cnvs.end();them++)
		{
			//only bother calculating overlap if we are on the same chromosome and we're not counting ourself
			if(me->chromosome_ == them->GetChromosome() && me->GetOurId() != them->GetUniqueId())
			{
				double overlap_cnv = CalcOverlapCNV(me->start_, them->GetStart(), me->stop_, them->GetStop());
				
				//# Exact Unique Times Seen tracking
				//LOH is already ok here...
				if(me->start_ == them->GetStart() && me->stop_ == them->GetStop()
				   && (me->copy_num_ == them->GetCopyNum()) && me->our_id_ != them->GetUniqueId() )
				{
					//first lets isolate only the disorder part of the CNV id.
					//to do this, search for first _ then take substring from 0- first _
					
					string temp_our_id = me->our_id_;
					string temp_their_id = them->GetUniqueId();
					
					unsigned int position_found = me->our_id_.find("_");
					if (position_found!=std::string::npos)
						temp_our_id.erase(position_found);
					position_found = them->GetUniqueId().find("_");
					if (position_found!=std::string::npos)
						temp_their_id.erase(position_found);						
						
					if(temp_our_id == temp_their_id)	//same disorder
					{
						if(them->GetCaseControlStatus() == 2)	//it is this disorder, but it's a control
						{
							me->num_exact_unique_++;
								//if(me->our_id_== "HTXY_B_3668") cout << them->GetUniqueId() << endl;
						}
					}
					else{
						me->num_exact_unique_++;
								//if(me->our_id_== "HTXY_B_3668") cout << them->GetUniqueId() << endl;
					}
				}
				if(overlap_cnv>0)
				{	//If there is some overlap with our CNV, then lets keep track of it
					
					unsigned int dummy_start, dummy_stop;
					dummy_start = max (me->start_,them->GetStart());
					dummy_stop = min (me->stop_,them->GetStop());
					
					GenomicRegion my_overlap(me->chromosome_,dummy_start,dummy_stop);
					
								/*
								vector<GenomicRegion> cumulativeDB_overlaps_exacttype_;
								list<GenomicRegion> cumulativeDB_coalesced_overlaps_exacttype_;
								vector<GenomicRegion> cumulativeDB_overlaps_simtype_;
								list<GenomicRegion> cumulativeDB_coalesced_overlaps_simtype_;
								vector<GenomicRegion> cumulativeDB_overlaps_difftype_;
								list<GenomicRegion> cumulativeDB_coalesced_overlaps_difftype_;
								double cumulativeDB_percent_overlapped_exacttype_, cumulativeDB_percent_overlapped_simtype_, cumulativeDB_percent_overlapped_difftype_;
								string which_cumulativeDB_overlaps_exacttype_, which_cumulativeDB_overlaps_simtype_, which_cumulativeDB_overlaps_difftype_;
								*/
					
					//exact type
					//This is already LOH ok
					if(me->copy_num_ == them->GetCopyNum())
					{
						me->cumulativeDB_overlaps_exacttype_.push_back(my_overlap);
						me->cumulativeDB_percent_overlapped_exacttype_ = me->cumulativeDB_percent_overlapped_exacttype_ + overlap_cnv;
						if(me->which_cumulativeDB_overlaps_exacttype_ != "")
							me->which_cumulativeDB_overlaps_exacttype_ = me->which_cumulativeDB_overlaps_exacttype_ + ";";
						me->which_cumulativeDB_overlaps_exacttype_ = me->which_cumulativeDB_overlaps_exacttype_ + them->GetUniqueId();
						
					}
					//similar type
					else if (me->copy_num_ == (them->GetCopyNum() - 1) || me->copy_num_ == (them->GetCopyNum() + 1) )
					{
						//LOH check
						if(me->copy_num_ != 2 && them->GetCopyNum() != 2)
						{
							me->cumulativeDB_overlaps_simtype_.push_back(my_overlap);
							me->cumulativeDB_percent_overlapped_simtype_ = me->cumulativeDB_percent_overlapped_simtype_ + overlap_cnv;
							if(me->which_cumulativeDB_overlaps_simtype_ != "")
								me->which_cumulativeDB_overlaps_simtype_ = me->which_cumulativeDB_overlaps_simtype_ + ";";
							me->which_cumulativeDB_overlaps_simtype_ = me->which_cumulativeDB_overlaps_simtype_ + them->GetUniqueId();		
						}
					}
					//different type
					else
					{
						//LOH check
						if(me->copy_num_  != 2 && them->GetCopyNum() != 2)
						{
							me->cumulativeDB_overlaps_difftype_.push_back(my_overlap);
							me->cumulativeDB_percent_overlapped_difftype_ = me->cumulativeDB_percent_overlapped_difftype_ + overlap_cnv;
							if(me->which_cumulativeDB_overlaps_difftype_ != "")
								me->which_cumulativeDB_overlaps_difftype_ = me->which_cumulativeDB_overlaps_difftype_ + ";";
							me->which_cumulativeDB_overlaps_difftype_ = me->which_cumulativeDB_overlaps_difftype_ + them->GetUniqueId();						
						}
					}
				}
			}
		}
		//for exact type overlaps
		//copy our overlap_bases vector to our coalesced_overlaps list
		copy(me->cumulativeDB_overlaps_exacttype_.begin(), me->cumulativeDB_overlaps_exacttype_.end(), back_inserter(me->cumulativeDB_coalesced_overlaps_exacttype_));
		//pass our coalesced_overlaps vector to the coalesce function
		CoalesceOverlaps(me->cumulativeDB_coalesced_overlaps_exacttype_);	
		//calculate our percent overlapped with the coalesced overlaps
		me->cumulativeDB_percent_overlapped_exacttype_ = CalcTotalOverlap(*me,me->cumulativeDB_coalesced_overlaps_exacttype_);
		
		//for similar type overlaps
		//copy our overlap_bases vector to our coalesced_overlaps list
		copy(me->cumulativeDB_overlaps_simtype_.begin(), me->cumulativeDB_overlaps_simtype_.end(), back_inserter(me->cumulativeDB_coalesced_overlaps_simtype_));
		//pass our coalesced_overlaps vector to the coalesce function
		CoalesceOverlaps(me->cumulativeDB_coalesced_overlaps_simtype_);	
		//calculate our percent overlapped with the coalesced overlaps
		me->cumulativeDB_percent_overlapped_simtype_ = CalcTotalOverlap(*me,me->cumulativeDB_coalesced_overlaps_simtype_);
		
		//for different type overlaps
		//copy our overlap_bases vector to our coalesced_overlaps list
		copy(me->cumulativeDB_overlaps_difftype_.begin(), me->cumulativeDB_overlaps_difftype_.end(), back_inserter(me->cumulativeDB_coalesced_overlaps_difftype_));
		//pass our coalesced_overlaps vector to the coalesce function
		CoalesceOverlaps(me->cumulativeDB_coalesced_overlaps_difftype_);	
		//calculate our percent overlapped with the coalesced overlaps
		me->cumulativeDB_percent_overlapped_difftype_ = CalcTotalOverlap(*me,me->cumulativeDB_coalesced_overlaps_difftype_);
	}
	
}

void CnvVector::CheckCnvr(vector<GenomicRegion>& some_blocks, const string& which_cnvr)
{
	for (vector<Cnv>::iterator me = cnv_vector_.begin();me != cnv_vector_.end();me++)
	{
		for(vector<GenomicRegion>::iterator them = some_blocks.begin();them != some_blocks.end();them++)
		{
			if(me->chromosome_ == them->GetChromosome())
			{
				double overlap_cnv = CalcOverlapCNV(me->start_, them->GetStart(), me->stop_, them->GetStop());
			
				if(overlap_cnv>0)
				{	//If there is some overlap with our CNV, then lets keep track of it		
					unsigned int dummy_start, dummy_stop;
					dummy_start = max (me->start_,them->GetStart());
					dummy_stop = min (me->stop_,them->GetStop());
					
					GenomicRegion my_overlap(me->chromosome_,dummy_start,dummy_stop);
							/*
							string chop_overlaps_ids_, chop_percentages_, chop_frequencies_;
							double chop_total_percent_;
							string hapmap_overlaps_ids_, hapmap_percentages_, hapmap_frequencies_;
							double hapmap_total_percent_;
							*/					
					if(which_cnvr == "chop")
					{
						if(me->chop_overlaps_ids_ != "")	me->chop_overlaps_ids_ = me->chop_overlaps_ids_ + ";";	//append ;
						me->chop_overlaps_ids_ = me->chop_overlaps_ids_ + them->GetUniqueId();
						if(me->chop_percentages_ != "")	me->chop_percentages_ = me->chop_percentages_ + ";";	//append ;
						me->chop_percentages_ = me->chop_percentages_ + FormatOverlap(overlap_cnv);
						if(me->chop_frequencies_ != "")	me->chop_frequencies_ = me->chop_frequencies_ + ";";	//append ;
						me->chop_frequencies_ = me->chop_frequencies_ + them->GetUniqueId2();
						me->chop_overlaps_.push_back(my_overlap);
						
					}
					else if(which_cnvr == "hapmap")
					{
						if(me->hapmap_overlaps_ids_ != "")	me->hapmap_overlaps_ids_ = me->hapmap_overlaps_ids_ + ";";	//append ;
						me->hapmap_overlaps_ids_ = me->hapmap_overlaps_ids_ + them->GetUniqueId();
						if(me->hapmap_percentages_ != "")	me->hapmap_percentages_ = me->hapmap_percentages_ + ";";	//append ;
						me->hapmap_percentages_ = me->hapmap_percentages_ + FormatOverlap(overlap_cnv);
						if(me->hapmap_frequencies_ != "")	me->hapmap_frequencies_ = me->hapmap_frequencies_ + ";";	//append ;
						me->hapmap_frequencies_ = me->hapmap_frequencies_ + them->GetUniqueId2();
						me->hapmap_overlaps_.push_back(my_overlap);
					
					}
					else if(which_cnvr == "pseudogene")
					{
						if(me->pseudo_overlaps_ids_ != "")	me->pseudo_overlaps_ids_ = me->pseudo_overlaps_ids_ + ";";	//append ;
						me->pseudo_overlaps_ids_ = me->pseudo_overlaps_ids_ + them->GetUniqueId();
						if(me->pseudo_percentages_ != "")	me->pseudo_percentages_ = me->pseudo_percentages_ + ";";	//append ;
						me->pseudo_percentages_ = me->pseudo_percentages_ + FormatOverlap(overlap_cnv);
						me->pseudo_overlaps_.push_back(my_overlap);
					}
					else if(which_cnvr == "parentgene")
					{
						if(me->parent_overlaps_ids_ != "")	me->parent_overlaps_ids_ = me->parent_overlaps_ids_ + ";";	//append ;
						me->parent_overlaps_ids_ = me->parent_overlaps_ids_ + them->GetUniqueId();
						if(me->parent_percentages_ != "")	me->parent_percentages_ = me->parent_percentages_ + ";";	//append ;
						me->parent_percentages_ = me->parent_percentages_ + FormatOverlap(overlap_cnv);
						me->parent_overlaps_.push_back(my_overlap);
					}						
				}
			}
		}
		if(which_cnvr == "chop")
		{
			//copy our overlap_bases vector to our coalesced_overlaps list
			copy(me->chop_overlaps_.begin(), me->chop_overlaps_.end(), back_inserter(me->chop_coalesced_overlaps_));
			
			//pass our coalesced_overlaps vector to the coalesce function
			CoalesceOverlaps(me->chop_coalesced_overlaps_);	
			
			//calculate our percent overlapped with the coalesced overlaps
			me->chop_total_percent_ = FormatOverlap( CalcTotalOverlap(*me,me->chop_coalesced_overlaps_) );
			
			if(me->chop_overlaps_ids_ == "")	me->chop_overlaps_ids_ = "NA";
			if(me->chop_percentages_ == "")	me->chop_percentages_ = "NA";
			if(me->chop_frequencies_ == "")	me->chop_frequencies_ = "NA";
		}
		else if(which_cnvr == "hapmap")
		{
			//copy our overlap_bases vector to our coalesced_overlaps list
			copy(me->hapmap_overlaps_.begin(), me->hapmap_overlaps_.end(), back_inserter(me->hapmap_coalesced_overlaps_));
			
			//pass our coalesced_overlaps vector to the coalesce function
			CoalesceOverlaps(me->hapmap_coalesced_overlaps_);	
			
			//calculate our percent overlapped with the coalesced overlaps
			me->hapmap_total_percent_ =  FormatOverlap( CalcTotalOverlap(*me,me->hapmap_coalesced_overlaps_) );		
			
			if(me->hapmap_overlaps_ids_ == "")	me->hapmap_overlaps_ids_ = "NA";
			if(me->hapmap_percentages_ == "")	me->hapmap_percentages_ = "NA";
			if(me->hapmap_frequencies_ == "")	me->hapmap_frequencies_ = "NA";			
		}
		else if(which_cnvr == "pseudogene")
		{
			//copy our overlap_bases vector to our coalesced_overlaps list
			copy(me->pseudo_overlaps_.begin(), me->pseudo_overlaps_.end(), back_inserter(me->pseudo_coalesced_overlaps_));
			
			//pass our coalesced_overlaps vector to the coalesce function
			CoalesceOverlaps(me->pseudo_coalesced_overlaps_);	
			
			//calculate our percent overlapped with the coalesced overlaps
			me->pseudo_total_percent_ =  FormatOverlap( CalcTotalOverlap(*me,me->pseudo_coalesced_overlaps_) );		
			
			if(me->pseudo_overlaps_ids_ == "")	me->pseudo_overlaps_ids_ = "NA";
			if(me->pseudo_percentages_ == "")	me->pseudo_percentages_ = "NA";
		}
		else if(which_cnvr == "parentgene")
		{
			//copy our overlap_bases vector to our coalesced_overlaps list
			copy(me->parent_overlaps_.begin(), me->parent_overlaps_.end(), back_inserter(me->parent_coalesced_overlaps_));
			
			//pass our coalesced_overlaps vector to the coalesce function
			CoalesceOverlaps(me->parent_coalesced_overlaps_);	
			
			//calculate our percent overlapped with the coalesced overlaps
			me->parent_total_percent_ =  FormatOverlap( CalcTotalOverlap(*me,me->parent_coalesced_overlaps_) );		
			
			if(me->parent_overlaps_ids_ == "")	me->parent_overlaps_ids_ = "NA";
			if(me->parent_percentages_ == "")	me->parent_percentages_ = "NA";
		}
	
	}
}








