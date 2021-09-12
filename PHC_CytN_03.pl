#!/usr/bin/perl 

#This script filters a list of SNPs in simplified Hapmap format for their QUAL phred-scale quality score Q.

@SHap_lines = <STDIN>;
$i = 0;
foreach (@SHap_lines){
	@snp_data = split (/\s+/,$SHap_lines[$i]);
	if ($snp_data[0] > 100) {
		print "$SHap_lines[$i]";
	}
	$i++;
}