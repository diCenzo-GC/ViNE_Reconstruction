
#!/usr/bin/perl

my $file = $ARGV[0];

print "\n working on file $file \n";


open(OPEN, $file);
@array = <OPEN>;

close OPEN;

#print @array;


foreach $line(@array){

	@SplitLine = split('\t', $line);
	#print "\n$SplitLine[0]";
	
	push(@SharedSinoMNXcode, "$SplitLine[0]\[e0\]\n");	

	@SplitLine2 = split(';', $SplitLine[1]);

	$cpd4EX = substr($SplitLine2[0], 0, -2);
	
	push(@EXListofMapped, "EX_$cpd4EX\_e0\n");

		foreach $oldline(@SplitLine2){
		chomp $oldline;
	
#identify the compartment and adjust the model name accordingly. First identify the right position of the compartment specification (different for the sino model (-1) and the truncatula model (-1))
	
		if ($file =~ /Sino/ ){

		$compartment_position = -2;	

		}

		else {	$compartment_position = -1}

		$compartment = substr $oldline, $compartment_position;
		#print "\n\nthis is the compartment -_$compartment- ";
		push(@FullArray, "$SplitLine[0]\t$oldline_$compartment\n");	

#if we are working with the sino metabolic model and if the compound is 'b' than replace it with the cytosol compartment of the plant. Also, in truncatula model, the compartment is appended at the end of mets as "[C]" and not as "_C". This is why I have introduced the following if statement.

		if ($file =~ /Sino/ ){
		
		if ($compartment =~ 'b') {

		$compartment = 'C';	

		}
	
		push(@NewNames, "$SplitLine[0]\[$compartment\]\n");


		push(@OldNames, "$oldline\n");

		}

		else {
	
		push(@NewNames, "$SplitLine[0]\[$compartment\]\n");

		$oldline_edit = substr($oldline, 0, -2);

		push(@OldNames, "$oldline_edit\[$compartment\]\n");


		}

		}

}


open(WRITE, ">$file.NewMetNames");
print WRITE @NewNames;
close WRITE;

open(WRITE, ">$file.OldMetNames");
print WRITE @OldNames;
close WRITE;

open(WRITE, ">$file.ExListofShared");
print WRITE @EXListofMapped;
close WRITE;

open(WRITE, ">$file.ExListofShared_MNXcode");
print WRITE @SharedSinoMNXcode;
close WRITE;




