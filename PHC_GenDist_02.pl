#!/usr/bin/perl 

#This script takes vcf files as input. It removes SNPs whose number of reads corresponding to the reference allele (RO) and to the alternative allele (AO) are both equal to 0. 

@vcf_lines = <STDIN>;

for ($j = 0; $j < scalar(@vcf_lines); $j++) {
    chomp;
    @snp_data = split (/\s+/,$vcf_lines[$j]);
    $false = 0;
	for ($i = 9; $i<= 38; $i++) {
        @isd = split (/\:/,$snp_data[$i]);
        if ($isd[2] == 0 && $isd[4] == 0){
          	$false++;
        }
    }
    unless ($false != 0){
        print "$vcf_lines[$j]";
    }
}
