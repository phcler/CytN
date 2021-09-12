#!/usr/bin/perl 

#This script identifies SNPs causing non-synonymous mutations using the information provided by script PHC_CytN_08.pl.
#The script requires three arguments for each SNP to run: 1) the rank in the file of filtered SNP; 2) the reference allele; 3) the alternative allele.

open IN,"snp_9cytn32_09.txt";
@lines = <IN>;

#snp rank:
$arg1 = shift@ARGV;
#reference allele:
$arg2 = shift@ARGV;
#laternative allele:
$arg3 = shift@ARGV;

$gen_code{"F"} = ("TTT TTC");
$gen_code{"L"} = ("TTA TTG CTT CTC CTA CTG");
$gen_code{"I"} = ("ATT ATC ATA");
$gen_code{"M"} = ("ATG");
$gen_code{"V"} = ("GTT GTC GTA GTG");
$gen_code{"S"} = ("TCT TCC TCA TCG AGT AGC");
$gen_code{"P"} = ("CCT CCC CCA CCG");
$gen_code{"T"} = ("ACT ACC ACA ACG");
$gen_code{"A"} = ("GCT GCC GCA GCG");
$gen_code{"Y"} = ("TAT TAC");
$gen_code{"stop"} = ("TAG TAA TGA");
$gen_code{"H"} = ("CAT CAC");
$gen_code{"Q"} = ("CAA CAG");
$gen_code{"N"} = ("AAT AAC");
$gen_code{"K"} = ("AAA AAG");
$gen_code{"D"} = ("GAT GAC");
$gen_code{"E"} = ("GAA GAG");
$gen_code{"C"} = ("TGT TGC");
$gen_code{"W"} = ("TGG");
$gen_code{"R"} = ("CGT CGC CGA CGG AGA AGG");
$gen_code{"G"} = ("GGT GGC GGA GGG");

foreach (@lines){
	@ref_codons = ();
	@alt_codons = ();
	$aa_ref = ();
	$aa_alt = ();
	@snp_data = split (/\s+/,$_);
	@snp_number = split (/\./,$snp_data[0]);
	if ($snp_number[0] == $arg1 && $snp_data[2] eq "cds") {
		print "$arg1\t$snp_data[1]\t";
		if ($snp_data[4] =~ /frameshift/) {
			print "frameshift detected\n";
		} elsif ($snp_data[5] =~ /manual/) {
			print "manual check required\n";
		} else {
			$r = int((length($arg2)+1)/3)+1;
			$a = int((length($arg3)+1)/3)+1;
			#snp on positive strand (+)
			if ($snp_data[5] eq "plus") {
				print "+\t";
				#reference DNA sequence ref_seq (+)
				$ref_seq = substr($snp_data[9],0,3*$r);
				$ref_seq =~ tr/atgc/ATGC/;
				#reference amino-acids aa_ref (+)
				for ($n=0; $n<3*$r; $n+=3) {
					push(@ref_codons,substr($ref_seq,$n,3));
				}
				$i=0;
				foreach (@ref_codons) {
					foreach (keys %gen_code) {
						if ($gen_code{$_} =~ /$ref_codons[$i]/) {
							$aa_ref .= $_;
						}
					}
					$i++;
				}
				#alternative DNA sequence alt_seq (+)
				if ($snp_data[7] == 1) {	
					$alt_seq = $arg3.substr($snp_data[9],length($arg2),3*$a-length($arg3));
				} elsif ($snp_data[7] == 2) {
					$alt_seq = substr($snp_data[9],0,1).$arg3;
					$alt_seq .= substr($snp_data[9],length($arg2)+1,3*$a-length($arg3)-1);
				} elsif ($snp_data[7] == 3) {
					$alt_seq = substr($snp_data[9],0,2).$arg3;
					$alt_seq .= substr($snp_data[9],length($arg2)+2,3*$a-length($arg3)-2);
				}
				#alternative amino-acids aa_alt (+)
				for ($n=0; $n<3*$a; $n+=3) {
					push(@alt_codons,substr($alt_seq,$n,3));
				}
				$i=0;
				foreach (@alt_codons) {
					foreach (keys %gen_code) {
						if ($gen_code{$_} =~ /$alt_codons[$i]/) {
							$aa_alt .= $_;
						}
					}
					$i++;
				}
			#snp on negative strand (-)
			} elsif ($snp_data[5] eq "minus") {
				print "-\t";
				$anti_snp_data[9] = ();
				for ($n=0; $n<length($snp_data[9]);$n++) {
					$anti_snp_data[9] .= substr($snp_data[9],(length($snp_data[9])-1-$n),1);
				}
				$anti_snp_data[9] =~ tr/atgcATGC/tacgTACG/;
				#reference DNA sequence ref_seq (-)
				$anti_ref_seq = substr($anti_snp_data[9],0,3*$r);
				$ref_seq = ();
				for ($n=0; $n<length($anti_ref_seq);$n++) {
					$ref_seq .= substr($anti_ref_seq,(length($anti_ref_seq)-1-$n),1);
				}
				$ref_seq =~ tr/atgcATGC/tacgTACG/;
				$ref_seq =~ tr/atgc/ATGC/;
				#reference amino-acids aa_ref (-)
				for ($n=0; $n<3*$r; $n+=3) {
					push(@ref_codons,substr($ref_seq,$n,3));
				}
				$i=0;
				foreach (@ref_codons) {
					foreach (keys %gen_code) {
						if ($gen_code{$_} =~ /$ref_codons[$i]/) {
							$aa_ref .= $_;
						}
					}
					$i++;
				}
				#alternative DNA sequence alt_seq (-)
				$anti_alt_seq = ();
				if ($snp_data[7] == 3) {	
					$anti_alt_seq = $arg3.substr($anti_snp_data[9],length($arg2),3*$a-length($arg3));
				} elsif ($snp_data[7] == 2) {
					$anti_alt_seq = substr($anti_snp_data[9],0,1).$arg3;
					$anti_alt_seq .= substr($anti_snp_data[9],length($arg2)+1,3*$a-length($arg3)-1);
				} elsif ($snp_data[7] == 1) {
					$anti_alt_seq = substr($anti_snp_data[9],0,2).$arg3;
					$anti_alt_seq .= substr($anti_snp_data[9],length($arg2)+2,3*$a-length($arg3)-2);
				}
				$alt_seq = ();
				for ($n=0; $n<length($anti_alt_seq);$n++) {
					$alt_seq .= substr($anti_alt_seq,(length($anti_alt_seq)-1-$n),1);
				}
				$alt_seq =~ tr/atgcATGC/tacgTACG/;
				$alt_seq =~ tr/atgc/ATGC/;
				#alternative amino-acids aa_alt (-)
				for ($n=0; $n<3*$a; $n+=3) {
					push(@alt_codons,substr($alt_seq,$n,3));
				}
				$i=0;
				foreach (@alt_codons) {
					foreach (keys %gen_code) {
						if ($gen_code{$_} =~ /$alt_codons[$i]/) {
							$aa_alt .= $_;
						}
					}
					$i++;
				}
			}
			print "$ref_seq\t$snp_data[7]\t$aa_ref\t$alt_seq\t$aa_alt\t";
			if ($aa_alt eq $aa_ref) {
				print "synonymous\t";
			} else {
				print "non-synonymous\t";
			}
			print "$r $a\n";		
		}	
	}
}
