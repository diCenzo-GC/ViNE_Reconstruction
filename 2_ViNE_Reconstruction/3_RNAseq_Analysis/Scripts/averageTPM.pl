#!usr/bin/perl
use 5.010;

while(<>) {
	if(/zoneI_1/) {
		print("Gene\tzoneI\tzoneIId\tzoneIIp\tzoneIZ\tzoneIII\n");
	}
	else {
		chomp;
		@line = split("\t", $_);
		$zoneI = (@line[1] + @line[2] + @line[3]) / 3;
		$zoneIId = (@line[4] + @line[5] + @line[6]) / 3;
		$zoneIIp = (@line[7] + @line[8] + @line[9]) / 3;
		$zoneIZ = (@line[10] + @line[11] + @line[12]) / 3;
		$zoneIII = (@line[13] + @line[14] + @line[15]) / 3;
		say("@line[0]\t$zoneI\t$zoneIId\t$zoneIIp\t$zoneIZ\t$zoneIII");
	}
}
