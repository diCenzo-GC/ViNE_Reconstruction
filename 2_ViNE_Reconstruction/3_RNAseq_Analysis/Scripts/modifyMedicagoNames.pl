#!usr/bin/perl
use 5.010;

while(<>) {
	if(/>/) {
		chomp;
		@line = split(' \| ', $_);
		@positionA = split(':', @line[3]);
		@positionB = split('-', @positionA[1]);
		say("@line[0]__@positionB[0]__@positionB[1]");
	}
	else {
		print("$_");
	}
}
