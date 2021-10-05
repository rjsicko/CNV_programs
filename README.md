# CNV_programs
 
Various programs mostly related to copy number variant (CNV) studies. These programs do things like output overlap with various databases/regions/genes etc.
General use of these = edit program.config to specify input, output and reference files. save config (do not change config file name; programs are hardcoded to check for the specified config). Then run the program.   

Notes:
- These were written before I knew about file format standards - they use custom file formats instead of parsing standard format (i.e. vcf or bed files).
- Sorry, the input is not standard across all programs. Some programs include extra columns (phenotype, sex etc.) I created each for specifically what I needed at the time.
- If you use the code in your project, it would probably be best to adapt it to use standard file formats and then grab what you need.
