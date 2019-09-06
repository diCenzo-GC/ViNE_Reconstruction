#!usr/bin/perl
use 5.010;

# File variables
$sizeFile = @ARGV[1];
$countFile = @ARGV[0];

# Store gene lengths
$geneNumber = -1;
open($sizes, '<', $sizeFile);
while(<$sizes>) {
    @line = split("\t", $_);
    push(@lengths, @line[1]);
    push(@geneNames, @line[0]);
    $geneNumber++;
    push(@geneNumber, $geneNumber);
}
close($sizes);

# Store counts
open($htcount, '<', $countFile);
while(<$htcount>) {
    @line = split("\t", $_);
    push(@counts, @line[1]);
}
close($htcount);

# Calculate RPK value and add all RPK values together
$totalrpk = 0;
foreach $i (@geneNumber) {
    $rpk = @counts[$i] / (@lengths[$i] / 1000);
    push(@rpk, $rpk);
    $totalrpk = $totalrpk + $rpk;
}
$totalrpk = $totalrpk / 1000000;

# Calculate TPM values
foreach $i (@geneNumber) {
    $tpm = @rpk[$i] / $totalrpk;
    say("@geneNames[$i]\t$tpm");
}
