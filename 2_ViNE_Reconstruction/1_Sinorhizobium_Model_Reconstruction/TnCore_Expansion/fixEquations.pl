#!usr/bin/perl
use 5.010;

while(<>) {
    chomp;
    @line = split(' ', $_);
    $n = 1;
    foreach $i (@line) {
        if(index($i, 'cpd') != -1) {
            print("$i ($n)");
            $n = 1;
        }
        elsif($i eq '+') {
            print(' + ');
        }
        elsif($i eq '<==>') {
            print(' <==> ');
        }
        elsif($i eq '==>') {
            print(' ==> ');
        }
        else {
            $n = $i;
        }
    }
    print("\n");
}
