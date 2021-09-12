#!/usr/bin/perl 

#This script retrieves positions and unitig names from a list of SNPs in simplified Hapmap format.

@lines = <STDIN>;
$i=1;
foreach (@lines){
	@snp_data = split (/\s+/,$_);
	print "$snp_data[3]\t","$snp_data[2]\n";
	$i++;
}
