#!usr/bin/perl
use 5.010;

while(<>) {
    if(/>/) {
        @line = split('locus_tag=', $_);
	@line2 = split(']', @line[1]);
        say(">@line2[0]");
    }
    else {
        print("$_");
    }
}
