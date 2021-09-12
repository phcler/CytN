#!/usr/bin/perl

#This script detects ORFs in a list of transcript sequences in fasta format.

open IN, "Hannosum_v2.FilteredModels1.transcripts.fsa";

%transcripts = ();

while (<IN>) {
 	chomp;
	if (/^>/) {
		/>(.*)/;
		$transcript_name = $_;
	} else {
		$transcripts{$transcript_name} .= $_;
	}
}

foreach (keys %transcripts) {
	print "$_\n";

	%length_ORF_strand = ();
	%sequence_ORF_strand = ();

	for ($s=0; $s<=2; $s++) {
		$j = 0;
		@strand = ();
		for ($i=0+$s; $i<length($transcripts{$_})-2; $i+=3) {
			$strand[$j] = substr($transcripts{$_},$i,3);
			$j++;
		}

		$start_codon = 0;
		for ($k = 0; $k<scalar(@strand); $k++) {
			if ($strand[$k] eq "ATG") {
				$start_codon = $k + 1;
				last;
			}
		}

		$stop_codon = 0;
		if ($start_codon == 0) {
			for ($l = 0; $l<scalar(@strand); $l++) {
				if ($strand[$l] eq "TAA" || $strand[$l] eq "TAG" || $strand[$l] eq "TGA") {
					$stop_codon = $l + 1;
					last;
				}
			}
		} else {
			for ($l = $start_codon; $l<scalar(@strand); $l++) {
				if ($strand[$l] eq "TAA" || $strand[$l] eq "TAG" || $strand[$l] eq "TGA") {
					$stop_codon = $l + 1;
					last;
				}
			}
		}
		if ($stop_codon == 0) {
		}

		if ($stop_codon != 0) {
			$length_ORF_strand{$s+1} = ($stop_codon-$start_codon)*3;
			if ($start_codon != 0) {
				for ($m = $k; $m<=$l; $m++) {
					$sequence_ORF_strand{$s+1} .= $strand[$m]; 
				}
			} else {
				for ($m = 0; $m<=$l; $m++) {
					$sequence_ORF_strand{$s+1} .= $strand[$m]; 
				}
			}
		} else {
			$length_ORF_strand{$s+1} = length($transcripts{$_})-($start_codon*3);
			if ($start_codon != 0) {
				for ($m = $k; $m<scalar(@strand); $m++) {
					$sequence_ORF_strand{$s+1} .= $strand[$m]; 
				}
			} else {
				for ($m = 0; $m<scalar(@strand); $m++) {
					$sequence_ORF_strand{$s+1} .= $strand[$m]; 
				}
			}
		}
	}

	if ($length_ORF_strand{1} > $length_ORF_strand{2}) {
		if ($length_ORF_strand{1} > $length_ORF_strand{3}) {
			print "$sequence_ORF_strand{1}\n";
		} else {
			print "$sequence_ORF_strand{3}\n";
		}
	} else {
		if ($length_ORF_strand{2} > $length_ORF_strand{3}) {
			print "$sequence_ORF_strand{2}\n";
		} else {
			print "$sequence_ORF_strand{3}\n";
		}
	}
}