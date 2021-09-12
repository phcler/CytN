#!/usr/bin/perl 

#This script retrieves a piece of genomic sequence surrounding a SNP. 
#The reference genome sequence in fasta format is used as input. 
#The script run with three arguments: 1) the position of the SNP; 2) the full name of the unitig bearing it; 3) the length of the fragment desired (with the SNP in the middle).

use List::Util 'min';
open IN,"pb_hgap_illumina_corrected_dp50_ao0.9_paired0.8.fasta";

@fasta_lines = <IN>;
$arg1 = shift@ARGV;
$arg2 = shift@ARGV;
$arg3 = shift@ARGV;
$i = 0;

if ($arg2 eq "unitig_0") {
    foreach (@fasta_lines) {
        if ($fasta_lines[$i] =~ /$arg2/ ) {
            $j = $i + 1;
                foreach (@fasta_lines) {
                        if (substr($fasta_lines[$j],0,1) eq ">") {
                                last;
                        } else {
                                chomp $fasta_lines[$j];
                                $unitig .= $fasta_lines[$j];    
                        }
                        $j++
                }
        }
        $i++;
    }
} else {
    $unitig_number = substr($arg2,7,3);
    foreach (@fasta_lines) {
       if (substr($fasta_lines[$i],8,3) == $unitig_number) {
        $j = $i + 1;
                foreach (@fasta_lines) {
                        if (substr($fasta_lines[$j],0,1) eq ">") {
                                last;
                        } else {
                                chomp $fasta_lines[$j];
                                $unitig .= $fasta_lines[$j];    
                        }
                        $j++
                }
        }
        $i++;
    }
}

print "".(substr($unitig,$arg1 - min(($arg3/2),$arg1), min($arg3,length($unitig)-$arg1+$arg3/2)))."\n";
