#!/usr/bin/perl 

#This script analyzes nucleotide blast results files in order to identify SNPs localized in an exon of a gene. 

#The script requires two arguments to run: 
#1) SNP position within the Query sequence; 
#2) maximum number of codons potentially modified by the alternative allele
#The script has seven possible outputs:

#1) "no hit found";
#2) "alignment score too low": cutoff = 180;
#3) "snp in 3' (or 5') UTR": UTR = untranslated regions in transcript sequence;
#4) "snp into intron";
#5) "cds frameshift": SNP possibly in an exon, but there is a frame shift between Subject and Query sequence;
#6) "cds manual check required": case requires manual sequence edition to conlude;
#7) "cds	738	T	plus	2219	2	740	CTACACAATCGC" (for example): SNP in an exon; alignment score between Query and Subject sequences; nucleotide at the position expected for the SNP in the Query sequence; strand on which the Query sequence aligns with the Subject sequence; position in the Subject sequence corresponding to the SNP in the Query sequence; frame used to translate the Query sequence in order to align with the translation of the Subject sequence; codon in the Subject sequence corresponding to the codon bearing the SNP in the Query sequence; codons in the Query sequence potentially modified by SNP allele substitution (here four).


@lines = <STDIN>;
$arg1 = shift@ARGV;
$arg2 = shift@ARGV;
$i = 0;
foreach (@lines) {
	if (/\*/) {
		print "no hit found";
		last;
	} elsif (/^>/) {
		/>(.*)/;
		chomp;
		print "$_\t";
		$j = $i + 7;
		$k = 0;
		$intS = 0;
		@match_endQ = ();
		@match_endS = ();
		@seqQ = ();
		@seqS = ();
		@strand = ();
		@lengthS = ();
		@score = split(/\s+/,$lines[$i+3]);
		if ($score[3] < 180) {
			print "alignment score too low\t";
			last;
		}
		until ($lines[$j] =~ /Lambda/ || $lines[$j] =~ /^>/) {
			if ($lines[$j] =~ /Query/) {
				@alignQ = split (/\s+/,$lines[$j]);
				@alignS = split (/\s+/,$lines[$j+2]);
				$match_endQ[2*$k] = $alignQ[1];
				$seqQ[$k] = $alignQ[2];
				$match_endQ[2*$k+1] = $alignQ[3];
				$match_endS[2*$k] = $alignS[1];
				$seqS[$k] = $alignS[2];
				$match_endS[2*$k+1] = $alignS[3];
				$k++;
			}
			$j++;
		}
		for ($n=0; $n<$k; $n++) {
			if (($seqS[$n] =~ /\-{20}$/ && $seqS[$n+1] =~ /^\-/) || ($seqS[$n] =~ /\-$/ && $seqS[$n+1] =~ /\-{20}/)) {
				$seqS[$n] =~ s/\-/\i/g;
				$seqS[$n+1] =~ s/\-/\i/g;
			} elsif ($seqS[$n] =~ /\-{20}/) {
				$seqS[$n] =~ s/\-/\i/g;	
			}
		}
		$startQ = 0;
		$startS = 0;
		$Q = ();
		$S = ();
		$endQ = 0;
		$endS = 0;
		for ($n=0; $n<$k; $n++) {
			if ($match_endQ[2*$n+1] >= $arg1 && $match_endQ[2*$n+1] < $arg1 + 3 && $n<$k-1 && $match_endQ[2*$n+1] + 1 == $match_endQ[2*$n+2]) {
				&double_sequence;
			} elsif ($match_endQ[2*$n+1] <= $arg1 && $match_endQ[2*$n+1] > $arg1 - 3 && $n<$k-1 && $match_endQ[2*$n+1] + 1 == $match_endQ[2*$n+2]) {
				&double_sequence;
			} elsif ($match_endQ[2*$n] <= $arg1 - 2 && $match_endQ[2*$n+1] >= $arg1 + 2) {
				&sequence;
			} elsif ($k>0 && $match_endQ[2*$k-2] <= $arg1 && $match_endQ[2*$k-1] > $arg1 + 3) {
				&sequence;
			} elsif ($match_endQ[0] >= $arg1 && $match_endQ[1] < $arg1 + 3) {
				&sequence;
			}
		}
		if ($Q ne ()) {
			$intS = ($S =~ tr/\i/\i/);
			if ($intS != 0 && substr($S,$arg1-$startQ,1) eq "i") {
				print "snp into intron (i)\t";
				last;
			} else {
				print "cds\t$score[3]\t";
			}
			$gapsQ = ($Q =~ tr/\-/\-/);
			$gapsS = ($S =~ tr/\-/\-/);
			if ($gapsS != $gapsQ) {
				$Dgap = ($gapsS - $gapsQ)/3;
				unless ($Dgap =~ /^-?\d+\z/) {
					print "frameshift detected";
					last;
				}
			}
			$Q =~ s/\-//g;
			@strand = split (/\//,$lines[$i+5]);
			if ($strand[1] eq "Plus\n") {
				$positionSpl3 = ($startS + $gapsS + $arg1 - $startQ + $gapsQ)/3;
				$positionSpl2 = ($startS + $gapsS + $arg1 - $startQ + $gapsQ +1)/3;
				$positionSpl1 = ($startS + $gapsS + $arg1 - $startQ + $gapsQ +2)/3;
				if ($positionSpl3 =~ /^\d+\z/) {
					$snp_codon = substr($Q,($arg1-2-$startQ),3*$arg2);
					print "".substr($Q,($arg1-$startQ),1)." ","plus ".3*$positionSpl3."\t3\t","$positionSpl3\t","$snp_codon\t";
				} elsif ($positionSpl2 =~ /^\d+\z/) {
					$snp_codon = substr($Q,($arg1-1-$startQ),3*$arg2);
					print "".substr($Q,($arg1-$startQ),1)." ","plus ".3*$positionSpl3."\t2\t","$positionSpl2\t","$snp_codon\t"; 
				} elsif ($positionSpl1 =~ /^\d+\z/) {
					$snp_codon = substr($Q,($arg1-$startQ),3*$arg2);
					print "".substr($Q,($arg1-$startQ),1)." ","plus ".3*$positionSpl3."\t1\t","$positionSpl1\t","$snp_codon\t"; 
				}
			} else {
				$positionSmn3 = ($endS + $endQ + $gapsS + $gapsQ - $arg1)/3;
				$positionSmn2 = ($endS + $endQ + $gapsS + $gapsQ - $arg1 + 1)/3;
				$positionSmn1 = ($endS + $endQ + $gapsS + $gapsQ - $arg1 + 2)/3;
				if ($positionSmn1 =~ /^\d+\z/) {
					$snp_anticodon = substr($Q,($arg1 - $startQ - 2),3*$arg2+2);
					&antiparallel;
					print "".substr($Q,($arg1-$startQ),1)." ","minus ".3*$positionSmn3."\t1\t","$positionSmn1\t","$snp_codon\t";
				} elsif ($positionSmn2 =~ /^\d+\z/) {
					$snp_anticodon = substr($Q,($arg1 - $startQ - 1),3*$arg2+1);
					&antiparallel;
					print "".substr($Q,($arg1-$startQ),1)." ","minus ".3*$positionSmn3."\t2\t","$positionSmn2\t","$snp_codon\t"; 
				} elsif ($positionSmn3 =~ /^\d+\z/) {
					$snp_anticodon = substr($Q,($arg1 - $startQ),3*$arg2);
					&antiparallel;
					print "".substr($Q,($arg1-$startQ),1)." ","minus ".3*$positionSmn3."\t3\t","$positionSmn3\t","$snp_codon\t"; 
				}	
			}
		} else {
			chomp $lines[$i+1];
			@lengthS = split(/\=/,$lines[$i+1]);
			@strand = split (/\//,$lines[$i+5]);
			@sortedQ = sort {$a <=> $b} @match_endQ;
			@sortedS = sort {$a <=> $b} @match_endS;
			$startQ = $sortedQ[0];
			$endQ = $sortedQ[2*$k-1];
			if ($strand[1] eq "Plus\n") {
				$startS = $sortedS[0];
				$endS = $sortedS[2*$k-1];
				if ($endQ < $arg1 && $arg1 - $endQ > $lengthS[1] - $endS) {
					print "snp in 3'UTR\t";
				} elsif ($startQ > $arg1 && $startQ - $arg1 > $startS - 1) {
					print "snp in 5'UTR\t";
				} elsif (($startS > 1 && $startQ > $arg1 && $startS - 1 > $startQ - $arg1) || ($endS < $lengthS[1] && $endQ < $arg1 && $arg1 - $endQ < $lengthS[1] - $endS)) {
					print "cds\t$score[3]\t+\tmanual check required\t";
				} elsif ($startQ < $arg1 && $endQ > $arg1) {
					for ($n=0; $n<$k; $n++) {
						if ($sortedQ[2*$n+2] > $arg1 && $sortedQ[2*$n+1] < $arg1) {
							if ($sortedQ[2*$n+2] != $sortedQ[2*$n] + 60) {
								if ($sortedS[2*$n+2] == $sortedS[2*$n+1] + 1 || $sortedS[2*$n+2] - $sortedS[2*$n+1] <= 20) {
									print "snp into intron (prob)\t";
								} elsif ($sortedS[2*$n+2] - $sortedS[2*$n+1] > 20) {
									print "cds\t$score[3]\t+\tmanual check required\t";
								}
							}
						} elsif ($sortedQ[2*$n+1] >= $arg1 && $sortedQ[2*$n+1] < $arg1 + 2) {
							print "cds\t$score[3]\t+\tmanual check required\t";
						}
					}
				}
			} elsif ($strand[1] eq "Minus\n") {
				@revsortedS = reverse @sortedS;
				$startS = $revsortedS[0];
				$endS = $revsortedS[2*$k-1];
				if ($endQ < $arg1 && $arg1 - $endQ > $endS - 1) {
					print "snp in 5'UTR\t";
				} elsif ($startQ > $arg1 && $startQ - $arg1 > $lengthS[1] - $startS) {
					print "snp in 3'UTR\t";
				} elsif (($endS > 1 && $endQ < $arg1 && $arg1 - $endQ < $endS - 1) || ($startS < $lengthS[1] && $startQ > $arg1 && $startQ - $arg1 < $lengthS[1] - $startS)) {
					print "cds\t$score[3]\t-\tmanual check required\t";
				} elsif ($startQ < $arg1 && $endQ > $arg1) {
					for ($n=0; $n<$k; $n++) {
						if ($sortedQ[2*$n+2] > $arg1 && $sortedQ[2*$n+1] < $arg1) {
							if ($sortedQ[2*$n+2] != $sortedQ[2*$n] + 60) {
								if ($revsortedS[2*$n+2] + 1 == $revsortedS[2*$n+1] || $revsortedS[2*$n+1] - $revsortedS[2*$n+2] <= 20) {
									print "snp into intron (prob)\t";
								} elsif ($revsortedS[2*$n+1] - $revsortedS[2*$n+2] > 20) {
									print "cds\t$score[3]\t-\tmanual check required\t";
								}
							}
						} elsif ($sortedQ[2*$n+1] >= $arg1 && $sortedQ[2*$n+1] < $arg1 + 2) {
							print "cds\t$score[3]\t-\tmanual check required\t";
						}
					}
				}
			} 
		}
		last;
	}
	$i++;
}
print "\n";

sub sequence {
	$startQ = $match_endQ[2*$n];
	$Q = $seqQ[$n];
	$endQ = $match_endQ[2*$n+1];
	$startS = $match_endS[2*$n];
	$S = $seqS[$n];
	$endS = $match_endS[2*$n+1];
}

sub double_sequence {
	$startQ = $match_endQ[2*$n];
	$Q = $seqQ[$n].$seqQ[$n+1];
	$endQ = $match_endQ[2*$n+3];
	$startS = $match_endS[2*$n];
	$S = $seqS[$n].$seqS[$n+1];
	$endS = $match_endS[2*$n+3];
}

sub antiparallel {
	$snp_codon = ();
	for ($n=0; $n<length($snp_anticodon);$n++) {
		$snp_codon .= substr($snp_anticodon,(length($snp_anticodon)-1-$n),1);
	}
	$snp_codon =~ tr/atgcATGC/tacgTACG/;
}