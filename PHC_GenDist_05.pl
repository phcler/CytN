#!/usr/bin/perl 

#This script takes vcf files as input. It removes SNPs found in the mitochondrial genome sequence and borne by unitig_36 and unitig_37.

@vcf_lines = <STDIN>;
$lines = scalar(@vcf_lines);

for ($j = 0; $j < $lines; $j++) {
	@snp_data = split(/\s+/,$vcf_lines[$j]);
	unless ($snp_data[0] eq "unitig_36" || $snp_data[0] eq "unitig_37") {
		print "$vcf_lines[$j]";
	}
}

