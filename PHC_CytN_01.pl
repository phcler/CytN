#!/usr/bin/perl 

#This script takes vcf files as input.
#Allelic selection is encoded by a string where "0" stands for the reference allele, "1" for the alternative allele and "\d" for a digit (either 0 or 1). 
#The order of the characters in this string follows the order of the 30 homokaryons used for parallel SNP calling (Clergeot, Rode et al. 2019). 

@vcf_lines = <STDIN>;
$i = 59;
foreach (@vcf_lines){
	chomp;
	@snp_data = split (/\s+/,$vcf_lines[$i]);
	$snp_string = ();
	for ($j = 9; $j <= 38; $j++) {
		@snp_info = split (/\:/,$snp_data[$j]);
		$snp_string .= $snp_info[0];
	}
	if ($snp_string =~ /00\d0\d\d\d\d\d\d\d0\d\d10\d0000\d\d\d\d\d\d\d10/) {
		print "$vcf_lines[$i]";
	}	
	$i++;
}