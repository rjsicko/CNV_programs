#ifndef OVERLAP_H
#define OVERLAP_H

#include <iostream>
#include <algorithm>
#include "cnv.h"

using namespace std;
 
double CalcOverlapUnion(int start_cnv, int start_region, int end_cnv, int end_region)
{
	return double( ( double( min(end_cnv,end_region) ) - max(start_cnv,start_region) ) )
							/
					( double( max(end_cnv,end_region) ) - min(start_cnv,start_region) );
}

double CalcOverlapRegion(int start_cnv, int start_region, int end_cnv, int end_region)
{
	return double ( ( double( min(end_cnv,end_region) ) - max(start_cnv,start_region) ) )
							/
					( double(end_region) - start_region );
}

double CalcOverlapCNV(int start_cnv, int start_region, int end_cnv, int end_region)
{
	return double( ( double( min(end_cnv,end_region) ) - max(start_cnv,start_region) ) )
							/
					( double(end_cnv) - start_cnv);
}

void CoalesceOverlaps(list<GenomicRegion>& overlap_regions)
{
	overlap_regions.sort();	//std sort since we defined >,< operators for GenomicRegion class
	
	//iterators for keeping track of the two nodes we are comparing
	list<GenomicRegion>::iterator me = overlap_regions.begin(),
							thee = overlap_regions.begin(),
                               end = overlap_regions.end();
							   	
	if ( me != end ) // Treat empty list, if not empty it's safe to make thee the second element
		thee++;	//sets it to the second element
	if(thee!=end)	//Treat list with one element, if its not a list with one element its safe to do the while loop
		while(thee != end)	//lets keep comparing until we are at the end
		{
			if(thee->GetStart() <= me->GetStop())	//hit to coalesce them
			{
				unsigned int temp_start = min(thee->GetStart(),me->GetStart()), temp_stop = max(thee->GetStop(),me->GetStop());
				GenomicRegion temp_region("1",temp_start,temp_stop);	//basic region constructor, used 1 for a dummy chromosome value

				overlap_regions.erase(me);
		
				list<GenomicRegion>::iterator temp_itr = overlap_regions.erase(thee);
						
				me = overlap_regions.insert(temp_itr,temp_region);
				thee = temp_itr;
			}
			else
			{
				me++;
				thee++;
			}
		}
}

//Here we calculate total overlap percentage from our coalesced overlaps
double CalcTotalOverlap(Cnv our_cnv,list<GenomicRegion>& overlap_regions)
{
	//cout << "Begining of CalcTotalOverlap function\n";
	double total_olap=0;
	for(list<GenomicRegion>::iterator me = overlap_regions.begin();me != overlap_regions.end();me++)
	{
		double temp_olap = CalcOverlapCNV(our_cnv.GetStart(),me->GetStart(),our_cnv.GetStop(),me->GetStop());
			//used for debugging
			/*if(our_cnv.GetOurId() == "HPPE_Control_0003")
			{
				cout << "our_cnv start: " << our_cnv.GetStart() << endl
					 << "genomic region start: " << me->GetStart() << endl
					 << "our_cnv start: " << our_cnv.GetStop() << endl
					 << "genomic region start: " << me->GetStop() << endl;
			}*/
		//cout << temp_olap << endl;
		total_olap = total_olap + temp_olap;
		//cout << total_olap << endl;
	}
	//cout << "End of CalcTotalOverlap function\n";
	return total_olap;
}	


 
#endif