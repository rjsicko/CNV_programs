/*****************************************************
cnvr_check_v1.2
Author: Robert Sicko - R.J.Sicko@gmail.com

****DESCRIPTION****
This program compares a list of CNVs to a list of "blocks" from the CHOP
database and outputs the original list with four additional columns appended.
The CHOP blocks each have an associated frequency. This program compares each 
CNV in the list to each block and will output which CNVs are overlapped
(above a threshold) by a block or multiple blocks.

The four columns appended are 'block_ids' 'block_freq' 'block_percents' and 'total_overlap'
if a CNV overlaps more than one block, the block_freq will be seperated
by a ; as will the percent overlaps for each
****END DESCRIPTION****

****REQUIRED FORMAT****
	NOTE: No commas in numbers (ie. use 1000000 not 1,000,000)
	NOTE2: Keep headers in your files. The program ignores the first line (don't have spaces in header id's though).
	
	CNV List (tab-delimited txt):
	OurID	Chr		Start	Stop	Case/Control	Copy_Number
	
	CHOP Block List(tab-delimited txt):
	Block_ID	chr		start	stop	freq
	
****END REQUIRED FORMAT****

****CHANGE LOG****
v1.0- Bugs fixed, programs appears to be functioning as expected
v1.1- Fixed a semi-bug that allowed 'total_overlap' to be greater than 100%.
	- This was caused by simply summing the percent of our CNV overlapped by each "block"
	- To see how this was corrected, see function CoalesceOverlaps... briefly, after all overlaps
	  to a particular CNV were added to the vector of OVERLAPS, these overlaps were coalesced,
	  leaving non overlapping OVERLAPS that were stored in a list. This list was then passed
	  to CalcTotalOverlap for calculating total overlap with our CNV.
v1.2- Added runlog output		  

****END CHANGE LOG****

****TO DO****
	
****END TO DO****
*********************************************************/

