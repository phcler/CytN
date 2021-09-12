#!/usr/bin/perl 

#This script shorthens the names of a list of jgi sequences in fasta format to their identification numbers only.

@lines = <STDIN>;
foreach (@lines){
	if(/^>/){
		@est_data = split (/\|/,$_);
		print ">","$est_data[2]\n";
	} else {
		print "$_";
	}
}
print "$i\n";