#!/usr/bin/perl 

#This script takes vcf files as input. It removes SNPs whose genotype is supported, for any isolate, by less than the fraction of the reads (between 0 and 1) given as argument on the command line.

@vcf_lines = <STDIN>;
$arg1 = shift@ARGV;

for ($j = 0; $j < scalar(@vcf_lines); $j++) {
    chomp;
    @snp_data = split (/\s+/,$vcf_lines[$j]);
    $count = 0;
    for ($i = 9; $i<= 38; $i++) {
        @isd = split (/\:/,$snp_data[$i]);
        $fRO = $isd[2]/$isd[1];
        $fAO = $isd[4]/$isd[1];
        if ($isd[0] == 0){
                if ($fRO >= $arg1) {
                        $count++;
                }
        } elsif ($isd[0] == 1){
                if ($fAO >= $arg1) {
                        $count++;
                }
        }
    }
    if ($count == 30){
       print "$vcf_lines[$j]";
       }
}

