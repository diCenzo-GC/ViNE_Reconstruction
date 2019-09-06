#!usr/bin/perl
use 5.010;

while(<>) {
	if(/>/) {
		@line = split(' ', $_);
		@line[0] =~ s/lcl|//;
		say(@line[0]);
	}
	else {
		print($_);
	}
}
