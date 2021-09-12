#!/usr/bin/perl

@lines = <STDIN>;
$l = 0;

foreach (@lines){
	if ($lines[$l] =~ /non-synonymous/) {
		print "$lines[$l]";
	}
	$l++;
}