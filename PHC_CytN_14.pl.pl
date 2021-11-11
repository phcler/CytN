#!/usr/bin/perl 

@lines = <STDIN>;
$arg1 = shift@ARGV;
$i = 0;
for ($i=0;$i<scalar(@lines);$i+=2){
	@Numb_ORF = split(/\>/,$lines[$i]);
	if ($Numb_ORF[1] == $arg1) {
		print "$lines[$i]","$lines[$i+1]";
		last;
	} 
}

