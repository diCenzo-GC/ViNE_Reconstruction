#!usr/bin/perl
use 5.010;

# Input files
$medicagoFile = 'combinedGenome/blastAnalysis/Medicago_truncatula_blast_output_parsed.txt';
$sinoFile = 'combinedGenome/blastAnalysis/Sinorhizobium_meliloti_blast_output_parsed.txt';

# Print the header
say('##gff-version 3');

# Make the Medicago part of the gff file
open($medic, '<', $medicagoFile);
while(<$medic>) {
	chomp;
	@line = split("\t", $_);
	@name = split('\.', @line[0]);
	if(@line[4] eq 'plus') {
		say("@line[1]\tRefSeq\tgene\t@line[2]\t@line[3]\t.\t+\t.\tgene_id=@name[0];length=@line[7]");
	}
	else {
		say("@line[1]\tRefSeq\tgene\t@line[3]\t@line[2]\t.\t-\t.\tgene_id=@name[0];length=@line[7]");
        }
}
close($medic);

# Make the Sinorhizobium part of the gff file
open($sino, '<', $sinoFile);
while(<$sino>) {
	chomp;
	@line = split("\t", $_);
	if(@line[14] eq 'plus') {
		say("@line[1]\tRefseq\tgene\t@line[10]\t@line[11]\t.\t+\t.\tgene_id=@line[0]");
	}
	else {
		say("@line[1]\tRefseq\tgene\t@line[11]\t@line[10]\t.\t-\t.\tgene_id=@line[0]");
	}
}
close(@sino);
