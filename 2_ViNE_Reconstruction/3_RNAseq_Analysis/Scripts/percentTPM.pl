#!usr/bin/perl
use 5.010;

while(<>) {
	chomp;
	if(/zoneI/) {
		say("$_");
	}
	else {
		@line = split("\t", $_);
		if(@line[1] + @line[2] + @line[3] + @line[4] + @line[5] > 0) {
			$zoneI = 100 * (@line[1] / (@line[1] + @line[2] + @line[3] + @line[4] + @line[5]));
			$zoneIId = 100 * (@line[2] / (@line[1] + @line[2] + @line[3] + @line[4] + @line[5]));
			$zoneIIp = 100 * (@line[3] / (@line[1] + @line[2] + @line[3] + @line[4] + @line[5]));
			$zoneIZ = 100 * (@line[4] / (@line[1] + @line[2] + @line[3] + @line[4] + @line[5]));
			$zoneIII = 100 * (@line[5] / (@line[1] + @line[2] + @line[3] + @line[4] + @line[5]));
			say("@line[0]\t$zoneI\t$zoneIId\t$zoneIIp\t$zoneIZ\t$zoneIII");
		}
		else {
			say("@line[0]\t0\t0\t0\t0\t0");
		}
	}
}
