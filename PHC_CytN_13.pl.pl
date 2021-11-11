#!/usr/bin/perl 

@lines = <STDIN>;
$i = 0;
$j = 1;
$count = 0;
@first_line_data = split (/\s+/,$lines[0]);
foreach (@lines){
	if ($lines[$i] =~ /intron/) {
		if ($count == 0) {
			@line_data = split (/\s+/,$lines[$i]);
			$ORF_ID = $line_data[1];
			@ORF_num = split(/\>/,$line_data[1]);
			print "$ORF_num[1]\t";
			$count++;
		} else {
			@line_data = split (/\s+/,$lines[$i]);
			if ($line_data[1] ne $ORF_ID) {
				$count++;
				print "$j\n";
				$ORF_ID = $line_data[1];
				$j = 1;
				@ORF_num = split(/\>/,$line_data[1]);
				print "$ORF_num[1]\t";
			} else {
				$j++;
			}	
		}
	}
	$i++;
}
print "$j\n";
print "$count\n";