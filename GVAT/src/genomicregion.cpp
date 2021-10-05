#include "genomicregion.h"

using namespace std;

		/*class GenomicRegion{

		public:

			//DGV2 entry constructor
			GenomicRegion(string const&, string const&, string const&, string const&, int, int, int);
			
			//Gene or transcript constructor
			GenomicRegion(string const&, string const&, string const&, int, int, int);
			
			
			
		private:
			string description_;	//identifier of type
			string unique_id_;
			string unique_id_2;
			string chromosome_;
			int copy_num_;
			unsigned int start_, stop_;	//start and end points for each region
			unsigned int length_;	//length of region
		};*/
		

//constructor for DGV2 entry		
GenomicRegion::GenomicRegion(string const& chromosome, int start, int stop, int copy_num, double freq)
{
	//DGV_date	chr	start	end	variantsubtype	frequency
	//assign all of our private variables based on input from constructor
	//description_ = description;
	//unique_id_ = unique_id;
	chromosome_ = chromosome;
	copy_num_ = copy_num;
	start_ = start;
	stop_ = stop;
	freq_ = freq;
	
}	

//OurID		Chr		Start	Stop	Copy_num	case/control
	//constructor for CNVs from cumulative DB
GenomicRegion::GenomicRegion(string const& our_id, string const& chromosome, int start, int stop, int copy_num, int case_or_control)
{
	unique_id_ = our_id;
	chromosome_ = chromosome;
	start_ = start;
	stop_ = stop;
	copy_num_ = copy_num;
	case_or_control_ = case_or_control;
}

//constructor for one of our CNVs for compare algorithms
	//			algorithm >> our_id >>     study_id >> case_or_control >> chromosome >> start >> stop >> copy_num;
GenomicRegion::GenomicRegion(string const& algorithm, string const& our_id, string const& study_id, int case_or_control, string const& chromosome, int start, int stop, int copy_num)
{
	description_ = algorithm;
	unique_id_ = our_id;
	unique_id_2_ = study_id;
	case_or_control_ = case_or_control;
	chromosome_ = chromosome;
	copy_num_ = copy_num;
	start_ = start;
	stop_ = stop;
}	

//constructor for gene
GenomicRegion::GenomicRegion(string const& description, string const& unique_id, string const& unique_id_2, string const& chromosome, int start, int stop)
{
	//assign all of our private variables based on input from constructor
	description_ = description;
	unique_id_ = unique_id;
	unique_id_2_ = unique_id_2;
	chromosome_ = chromosome;
	start_ = start;
	stop_ = stop;
	
}


//constructor for cnvr block
////Block_ID	chr	start	stop	freq
GenomicRegion::GenomicRegion(string const& block_id, string const& chromosome, int start, int stop, string const& freq)
{
	unique_id_ = block_id;
	chromosome_ = chromosome;
	start_ = start;
	stop_ = stop;
	unique_id_2_ = freq;
}

//constructor for most basic type: chr, start, stop
GenomicRegion::GenomicRegion(string const& chromosome, int start, int stop)
{
	//assign all of our private variables based on input from constructor
	chromosome_ = chromosome;
	start_ = start;
	stop_ = stop;
}
	
bool operator<(const GenomicRegion& lhs, const GenomicRegion& rhs)
{
	return lhs.start_ < rhs.start_;
}

bool operator>(const GenomicRegion& lhs, const GenomicRegion& rhs)
{
	return rhs < lhs;
}

bool operator<=(const GenomicRegion& lhs, const GenomicRegion& rhs)
{
	return !(rhs < lhs);
}

bool operator>=(const GenomicRegion& lhs, const GenomicRegion& rhs)
{
	return !(lhs < rhs);
}

unsigned int GenomicRegion::GetStart() const
{
	return start_;
}

unsigned int GenomicRegion::GetStop() const
{
	return stop_;
}

int GenomicRegion::GetCaseControlStatus() const
{
	return case_or_control_;
}

string GenomicRegion::GetChromosome() const
{
	return chromosome_;
}

int GenomicRegion::GetCopyNum() const
{
	return copy_num_;
}

string GenomicRegion::GetUniqueId() const
{
	return unique_id_;
}

string GenomicRegion::GetUniqueId2() const
{
	return unique_id_2_;
}

string GenomicRegion::GetDescription() const
{
	return description_;
}

