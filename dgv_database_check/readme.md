/*****************************************************
DGV_database_check_v1.7
Author: Robert Sicko - R.J.Sicko@gmail.com

****DESCRIPTION****
This program compares a list of CNVs to DGV2(or other) CNV database
and outputs the original list with an additional column appended.
'overlap_info'


Overlap_info:

		'% Our CNV'	= percentage of our CNV overlapped 
		'% DGV CNV' = percentage of DGV CNV overlapped
		'%Union' 	= percentage shared overlap in a union overlap calculation
						         ------|-----|---------------------   Our CNV
							 ----------+++++++++++++++---------   DGV2 CNV
							 ---++++++-------------------------   DGV2 CNV(2)							 
							 ----------XXX---------------------   Intersection
							 ------XXXXXXXXXXXXXXXXXXX---------   Denominator for union overlap

The DGV2 database can be found here: http://dgv.tcag.ca/dgv/app/home?ref=
****END DESCRIPTION****

****REQUIRED FORMAT****
	NOTE: No commas in numbers (ie. use 1000000 not 1,000,000)
	NOTE2: Keep headers in your files. The program ignores the first line (don't have spaces in header id's though).
	
	CNV List (tab-delimited txt):
	OurID	Chr		Start	Stop	Case/Control	Copy_Number
	
	DGV Database File(tab-delimited txt):
	DGV_date	chr	start	end	variantsubtype	frequency
****END REQUIRED FORMAT****

****CHANGE LOG****
v1.0- Bugs fixed, programs appears to be functioning as expected
v1.2- Input file for DGV database modified, now includes copy_number in last column
	- Now our check compares copy number to compare only dups with dups and dels with dels
v1.2- Added runlog output
v1.3- Output changed. Now six columns output:
		'% Our CNV'	= percentage of our CNV overlapped 
		'% DGV CNV' = percentage of DGV CNV overlapped
		'%Union' 	= percentage shared overlap in a union overlap calculation
		these three fields for both, same and different type of cnvs
v1.4- New output, more informative.
	- Similar overlap calculation as in pseudogene check.
	- Total % of our CNV covered by 1 or more DGV2 calls is now the output.
	- See overlap_info description above, previous versions didn't take into account that more than one DGV2
	  CNV could overlap our CNV. Now, we sum all overlap.
v1.5- Now use frequency data available for DGV2 entries
	- Mix the output form v1.3 and v1.4
v1.6- Fixed freq so the highest freq overlapped is output
v1.7- Fixed a bug - the highest frequency of the greatest overlap is now output instead of the overall highest frequency
****END CHANGE LOG****

****TO DO****
-fix freq so that it is the highest freq of a DGV2 overlap. 
****END TO DO****


*********************************************************/
