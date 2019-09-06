#!usr/bin/perl
use 5.010;

$lastLocus = 'XXX';
$start = 0;
$end = 0;
$startTemp = 0;
$endTemp = 0;
$identity = 0;
$totalLength = 0;
$geneLength = 0;
$chromosome = 0;
$strand = 0;

while(<>) {
	chomp;
	@line = split("\t", $_);
	if(@line[0] eq $lastLocus) {
		if(@line[14] eq $strand) {
			if(@line[14] eq 'plus') {
				if(@line[10] < $start) {
					$startTemp = @line[10];
				}
				if(@line[11] > $end) {
					$endTemp = @line[11];
				}
				if(abs($startTemp - $endTemp) < ($origLength + 2500)) {
					$start = $startTemp;
					$end = $endTemp;
					$identity = $identity + (@line[2] * @line[3]);
					$totalLength = $totalLength + @line[3];
				}
			}
			elsif(@line[14] eq 'minus') {
				if(@line[10] > $start) {
					$startTemp = @line[10];
				}
				if(@line[11] < $end) {
					$endTemp = @line[11];
				}
				if(abs($startTemp - $endTemp) < ($origLength + 2500)) {
					$start = $startTemp;
					$end = $endTemp;
					$identity = $identity + (@line[2] * @line[3]);
					$totalLength = $totalLength + @line[3];
				}
			}
			$geneLength = @line[6];
			$chromosome = @line[1];
			$strand = @line[14];
			$lastLocus = @line[0];
		}
	}
	elsif(@line[0] ne @lastLocus) {
		unless($lastLocus eq 'XXX') {
			$overalIdentity = $identity / $totalLength;
			$percentLength = 100 * $totalLength / $geneLength;
			$newLength = abs($start - $end);
			say("$lastLocus\t$chromosome\t$start\t$end\t$strand\t$overalIdentity\t$percentLength\t$totalLength\t$origLength\t$newLength");
			$start = 0;
			$end = 0;
			$startTemp = 0;
			$endTemp = 0;
			$identity = 0;
			$totalLength = 0;
			$geneLength = 0;
			$chromosome = 0;
			$strand = 0;
		}
		$lastLocus = @line[0];
		$start = @line[10];
		$end = @line[11];
		$startTemp = @line[10];
		$endTemp = @line[11];
		$identity = $identity + (@line[2] * @line[3]);
		$totalLength = $totalLength + @line[3];
		$geneLength = @line[6];
		$chromosome = @line[1];
		$strand = @line[14];
		@name = split('__', @line[0]);
		$origLength = abs(@name[1] - $name[2]);
	}
}		
$overalIdentity = 100 * $identity / $totalLength;
$percentLength = $totalLength / $geneLength;
$newLength = abs($start - $end);
say("$lastLocus\t$chromosome\t$start\t$end\t$strand\t$overalIdentity\t$percentLength\t$totalLength\t$origLength\t$newLength");
