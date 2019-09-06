#!usr/bin/perl
use 5.010;

$length = 0;
$end = 0;

while(<>) {
    @line = split("\t", $_);
    if(@line[2] eq 'gene') {
        say("$name\t$length");
        @info = split(';', @line[8]);
        $name = @info[0];
        $name =~ s/gene_id=//;
        $length = 0;
        $end = 0;
    }
    elsif(@line[2] eq 'five_prime_UTR') {
        if(@line[3] > $end) {
            if(@line[4] - @line[3] + 1 > 0) {
                $length = $length + (@line[4] - @line[3] + 1);
            }
        }
        else {
            if(@line[4] - $end > 0) {
                $length = $length + (@line[4] - $end);
            }
        }
        $end = @line[4];
    }
    elsif(@line[2] eq 'exon') {
        if(@line[3] > $end) {
            if(@line[4] - @line[3] + 1 > 0) {
                $length = $length + (@line[4] - @line[3] + 1);
            }
        }
        else {
            if(@line[4] - $end > 0) {
                $length = $length + (@line[4] - $end);
            }
        }
        $end = @line[4];
    }
    elsif(@line[2] eq 'three_prime_UTR') {
        if(@line[3] > $end) {
            if(@line[4] - @line[3] + 1 > 0) {
                $length = $length + (@line[4] - @line[3] + 1);
            }
        }
        else {
            if(@line[4] - $end > 0) {
                $length = $length + (@line[4] - $end);
            }
        }
        $end = @line[4];
    }
}
say("$name\t$length");

