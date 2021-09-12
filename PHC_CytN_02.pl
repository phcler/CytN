#!/usr/bin/perl 

#This script converts a vcf file to a list of SNPs in simplified Hapmap format: Q score; reference allele/alternative allele; unitig number; position; allele borne by each isolate (all separated by tabs; the reference isolate is the last in the list).

@vcf_lines = <STDIN>;
$lines = scalar(@vcf_lines);

for ($j = 0; $j < $lines; $j++) {
	@snp_data = split (/\s+/,$vcf_lines[$j]);
	print "$snp_data[5]\t","$snp_data[3]","/","$snp_data[4]\t","$snp_data[0]\t","$snp_data[1]\t";
	for ($i = 9; $i < 38; $i++) {
		@sub_data = split (/\:/,$snp_data[$i]);
		if ($sub_data[0] == 0) {
			print "$snp_data[3]\t";
		} else {
			print "$snp_data[4]\t";
		}
	}
	print "$snp_data[3]\n";
}