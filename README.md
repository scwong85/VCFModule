VCFModule
==========
This module process the vcf file from samtools and output a column delimited file with genotype information for each locus in the vcff file.


Files
-----
* parseVCF.py
* README.md


Using the module
-----------------
There are 2 modes for this module.
Mode 1: parsing vcf file
Mode 2: parsing vcf file with filtering criteria

Example:

Mode 1

    python processVCF.py <vcf_file> <output_file> mode
    python processVCF.py myvcf.vcf myoutput.txt 1

Mode 2

    python processVCF.py <vcf_file> <output_file> mode qualtity_score callrate maf
    #python processVCF.py myvcf.vcf myoutput.txt 2 900 0.8 0.2
    

Output
-------
Output column: 
Chromosome	Position	Ref	Alt	Quality	GenotypeCount	SampleSize	CallRate	MAF	HWEChiSq	inHWE
GenotypeCount is a dictionary data structure with genotype as key and count as item
