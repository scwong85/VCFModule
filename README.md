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

  python processVCF.py myvcf.vcf myoutput.txt 1
