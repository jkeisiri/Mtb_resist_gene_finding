# Mtb_resist_gene_finding
Using for finding resistant gene by searching SNP in variant calling files (VCF)

## Installation:
 -       pip install pandas

## Installation:
You have to add list of genes or list of mutations that you need to find in gene_list.csv file. The important parameters (Start, Stop and Strand) need to be filled. You can see an example from gene_list.csv file. For input files, you can use VCCF file with out heading information (example files attached). 

Usage:
 -       mtbtranslation.py -i [path of input] -o [path of output] -g [path of gene_list.csv] -a [path of Amino.csv] -r [path of Ref_Mtb.csv]

Example:
 -       mtbtranslation.py -i ./input/*.vcf -o ./output -g ./files/gene_list.csv -a ./files/Amino.csv -r ./files/Ref_Mtb.csv
  
