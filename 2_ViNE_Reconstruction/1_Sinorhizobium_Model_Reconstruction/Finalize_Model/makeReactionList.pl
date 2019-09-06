#!usr/bin/perl
use 5.010;

#start the reaction list
say('<listOfReactions>');
print("\n");

#prepare the reaction list
while(<>) {
	if(/TRUE/ or /FALSE/) {
		s/TRUE/true/g;
		s/FALSE/false/g;
		chomp;
		@reaction = split("\t",$_);
		#prepare header line of the reaction
		print('<reaction id="');
		print("@reaction[0]");
		print('" name="');
		print("@reaction[1]");
		print('" reversible="');
		print("@reaction[4]");
		say('">');
		#start notes section
		say("\t<notes>");
		#add genes
		print("\t\t");
		print('<html:p>GENE_ASSOCIATION:');
		print("@reaction[7]");
		say('</html:p>');
		#add proteins
		print("\t\t");
		print('<html:p>PROTEIN_ASSOCIATION:');
		print("@reaction[8]");
		say('</html:p>');
		#add KEGG_RID
		print("\t\t");
		print('<html:p>KEGG_RID:');
		print("@reaction[5]");
		say('</html:p>');
		#add PROTEIN_CLASS
		print("\t\t");
		print('<html:p>PROTEIN_CLASS:');
		print("@reaction[6]");
		say('</html:p>');
		#add support level
		print("\t\t");
		print('<html:p>SUPPORT_LEVEL:');
		print("@reaction[9]");
		say('</html:p>');
		#add references
		print("\t\t");
		print('<html:p>REFERENCES:');
		print("@reaction[10]");
		say('</html:p>');
		#add notes
		print("\t\t");
		print('<html:p>NOTES:');
		print("@reaction[11]");
		say('</html:p>');
		#end notes section
		say("\t</notes>");
		#prepare reaction
		@reaction[2] =~ s/\ <//g;
		@reaction[2] =~ s/>\ //g;
		@reaction[2] =~ s/\+/XXX/g;
		@compounds = split('==', @reaction[2]);
		@reactants = split(' XXX ', @compounds[0]);
		@products = split(' XXX ', @compounds[1]);
		#add reactants
		say("\t<listOfReactants>");
		foreach $reactant (@reactants) {
			print("\t\t");
			print('<speciesReference species="');
			@reactant2 = split("\ ", $reactant);
			print("@reactant2[0]");
			print('" stoichiometry="');
			@reactant2[1] =~ s/\(//g;
			@reactant2[1] =~ s/\)//g;
			print("@reactant2[1]");
			say('"/>');	
		}
		say("\t</listOfReactants>");
		#add products
		say("\t<listOfProducts>");
		foreach $product (@products) {
			print("\t\t");
			print('<speciesReference species="');
			@product2 = split("\ ", $product);
			print("@product2[0]");
			print('" stoichiometry="');
			@product2[1] =~ s/\(//g;
			@product2[1] =~ s/\)//g;
			print("@product2[1]");
			say('"/>');	
		}
		say("\t</listOfProducts>");
		#start kinetic laws
		say("\t<kineticLaw>");
		#add math
		print("\t\t");
		say('<math xmlns="http://www.w3.org/1998/Math/MathML">');
		print("\t\t\t");
		say('<ci> FLUX_VALUE </ci>');
		print("\t\t");
		say('</math>');
		#add parameters
		print("\t\t");
		say('<listOfParameters>');
		print("\t\t\t");
		if(/true/) {
			say('<parameter id="LOWER_BOUND" value="-1000" units="mmol_per_gDW_per_hr"/>');
		}
		if(/false/) {
			say('<parameter id="LOWER_BOUND" value="0" units="mmol_per_gDW_per_hr"/>');
		}
		print("\t\t\t");
		say('<parameter id="UPPER_BOUND" value="1000" units="mmol_per_gDW_per_hr"/>');
		print("\t\t\t");
		say('<parameter id="OBJECTIVE_COEFFICIENT" value="0"/>');
		print("\t\t\t");
		say('<parameter id="FLUX_VALUE" value="0.000000" units="mmol_per_gDW_per_hr"/>');
		print("\t\t");
		say('</listOfParameters>');
		print("\t");
		say('</kineticLaw>');
		#finish reaction
		say('</reaction>');
		print("\n");
	}
	else {
		chomp;
		@subheading = split("\t",$_);
		print('<!--');
		print("$subheading[0]");
		say("-->\n");
	}
}

#finish reaction list
say('</listOfReactions>');
print("\n");

#finish file
say('</model>');
say('</sbml>');
