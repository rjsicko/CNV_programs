#include<iostream>
#include<string>
#include "config.h"

using namespace std;


int main()
{
		/*
		cnv_filename	
		DGV2_filename	data_store/annotation_dbs/131121-DGV2_hg19.txt
		other_algorithm_filename	
		gencode_filename	data_store/annotation_dbs/131121-GENCODE_V17_HG19.txt
		ccds_filename	data_store/annotation_dbs/131121-CCDS_HG19.txt
		omim_filename	data_store/annotation_dbs/131121-OMIM_gene_ids.txt
		cumulative_database_filename	data_store/annotation_dbs/131121-Cumulative_CLEANED_DB.txt
		chop_cnvs_filename	data_store/annotation_dbs/CNVBlocks_hg19_freq1_1048.txt
		hapmap_cnvs_filename	data_store/annotation_dbs/HapMap3_CNPs_hg19_freq1_850.txt
		pathway_genes_filename	data_store/annotation_dbs/131121-SigTrans&DevBio_from_reactome.txt
		pseudo_gene_filename	data_store/annotation_dbs/PseudoGene_build69_modified_for_pseudo_gene_check_v1.0.txt
		pseudo_parent_filename	
		*/
	Config ourconfig("src/GVAT-0.0.2.config");	
	cout << ourconfig.GetFilename("DGV2_filename") << endl;
	
}