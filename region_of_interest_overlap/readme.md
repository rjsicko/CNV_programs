/*****************************************************
region_of_interest_overlap v1.4
Author: Robert Sicko - R.J.Sicko@gmail.com

****DESCRIPTION****
This program compares a list of CNVs to a list of regions of interest
and outputs the original list with two additional columns appended.
The two columns appended are 'Overlapped_Regions' and 'Percent_Overlap'
if a CNV overlaps more than one region, the regions will be seperated
by a ; as will the percent overlaps for each
Note: the majority of the annotation databases were downloaded from the UCSC table browser (https://genome.ucsc.edu/cgi-bin/hgTables) and formatted for use. 

****END DESCRIPTION****

****REQUIRED FORMAT****
	NOTE: No commas in numbers (ie. use 1000000 not 1,000,000)
	NOTE2: Keep headers in your files. The program ignores the first line (don't have spaces in header id's though).
	
	CNV List (tab-delimited txt):
	OurID	Chr		Start	Stop	Case/Control	Copy_Number
	
	Region List(tab-delimited txt):
	Database_From	Unique_ID_1(gene)	Unique_ID_2(transcript)		Chr		Start	Stop
****END REQUIRED FORMAT****

****CHANGE LOG****
v1.0- Bugs fixed, programs appears to be functioning as expected
v1.1- Added "NA" to 'Overlapped_Regions' and 'Percent_Overlap' fields if no regions were detected for a particular CNV
	- Added option in 'region_of_interest_overlap.config' to choose the type of overlap used
		Note: see 'Explanation_Overlap_Options.txt' for more information on the types of overlap
	- Commented a lot more of the code (I'll try to comment it all eventually)
	- Fixed a the 'Percent_Overlap' output so now allows 2 decimels
v1.2- Added option for writing 'transcript_id' field in output file
		As a consequence, the input file format now contains an extra 'transcript_id' field
		If you do not need this field, simply use a dummy code in this column and in config file specify 0 for write_transcript_id
		NOTE: Did not implement the do not write_transcript_id option yet. Therefore, write_transcript_id has to be used.
			  Will make option to not write transcript_id in v1.3!
	- Since I was already modifying the config file for the above option, I cleaned it up and removed entries not used
	- As a consequence of adding the write transcript_id option, the storage format has changed for gene_list and overlap_percent.
	  Each CNV now has a multimap<string,string> & a map<string,string> where the multimap is <gene_name,transcript_ids>
	  and the map is <transcript_id,overlap_percentage>.
	  Since the storage was changed to maps, the output is now alphabetically sorted.
	- "Bug" fixed when importing program output into Excel. In previous versions output was surrounded by (). For single transcript genes,
	  when the data was imported into Excel, the percent_overlap would be a negative number. This is because Excel automatically converts
	  (#) into -#. This can be changed when importing the data by telling Excel to treat the field as text. However, the easiest solution 
	  was to change the surrounding ()'s into []'s.
	- Added column for number_regions_overlapped
v1.3- Added 3 output files. A gene list file , a transcript output file and a log file.
	- The gene list file is a list of genes overlapped and the number cases and number controls that overlap it
	- The transcript list file is a list of genes overlapped and the number cases and number controls that overlap it
	- The log file is a file describing the particular run of the program
	- Added check for overlap_cutoff > 1; just incase user specifies 90 instead of 0.90 in config file!
	- Added check for config file manipulation.
		- count total line numbers and make sure user didnt remove any lines
v1.4- Changed chr variable from int to string
	- Added const string VERSION variable for tracking what version we are
****END CHANGE LOG****

****TO DO****
	- Add other checks for config file parameters outside expected range.
	- Comment out code!
	- Implement option to not write transcript.
	- Add check for output being too long for a single excel cell. In that case, use secondary column for output.
	- Fix case/control count so it doesn't just count the number of CNVs that overlap a gene, instead counts individuals
****END TO DO****
*********************************************************/
