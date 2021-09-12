#!/usr/bin/perl 

#This script retrieves the position and unitig name of SNPs identified as causing a non-synonymous mutation in an exon.
#The script requires three arguments to run: 1) the rank in the file of filtered SNP; 2) the SNP position; 3) the SNP unitig name 

open IN, "snp_9cytn32_11.txt";
@lines = <IN>;
$arg1 = shift@ARGV;
$arg2 = shift@ARGV;
$arg3 = shift@ARGV;
$l = 0;
foreach (@lines){
	@syn_data = split(/\s+/,$lines[$l]);
	if ($syn_data[0] == $arg1) {
		print "$arg2 $arg3\n";
		last;
	} else {
		$l++;
	}
}