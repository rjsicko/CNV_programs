opened OMIM_mim2gene in csved, removed duplicates and deleted the center columns (leaving omim# and gene).
saved the omim_gene_ids.txt file that contains only these OMIM ids and genes

I'll use this file in GVAT. Instead of searching for overlap with these genes again... I'll simply search the genecode genes overlapped by a particular CNV for the gene name in this file.