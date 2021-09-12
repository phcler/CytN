#!/usr/bin/perl 
open IN, "snp_9cytn32_04.txt";
@lines = <IN>;
$arg1 = shift@ARGV;
$arg2 = shift@ARGV;
$l = 0;
foreach (@lines){
	@snp_data = split(/\s+/,$lines[$l]);
	if ($snp_data[2] eq $arg2 && $snp_data[3] == $arg1) {
		print "$lines[$l]";
		last;
	} else {
		$l++;
	}
}