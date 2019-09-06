#!usr/bin/perl
use 5.010;

#prepare the header section of the model
say('<?xml version="1.0" encoding="UTF-8"?>');
say('<sbml xmlns="http://www.sbml.org/sbml/level2" level="2" version="1" xmlns:html="http://www.w3.org/1999/xhtml">');
say('<model id="iGD1575" name="Sinorhizobium_meliloti_1021">');
print("\n");
say('<listOfUnitDefinitions>');
print("\t");
say('<unitDefinition id="mmol_per_gDW_per_hr">');
print("\t\t");
say('<listOfUnits>');
print("\t\t\t");
say('<unit kind="mole" scale="-3"/>');
print("\t\t\t");
say('<unit kind="gram" exponent="-1"/>');
print("\t\t\t");
say('<unit kind="second" multiplier=".00027777" exponent="-1"/>');
print("\t\t");
say('</listOfUnits>');
print("\t");
say('</unitDefinition>');
say('</listOfUnitDefinitions>');
print("\n");
say('<listOfCompartments>');
print("\t");
say('<compartment id="c0" name="Cytosol_0" />');
print("\t");
say('<compartment id="e0" name="Extracellular_0" />');
say('</listOfCompartments>');
print("\n");
say('<listOfSpecies>');
print("\n");

#prepare the metabolite list
while(<>) {
	if(/TRUE/ or /FALSE/) {
		s/TRUE/true/g;
		s/FALSE/false/g;
		chomp;
		@metabolites = split("\t", $_);
		print('<species id="');
		print("@metabolites[0]");
		print('" name="');
		print("@metabolites[1]");
		print('" compartment="');
		print("@metabolites[2]");
        print('" hasOnlySubstanceUnits="false" boundaryCondition="');
        print("@metabolites[5]");
        say('" constant="false">');
        say("\t<notes>");
        print("\t\t");
        say('<body xmlns="http://www.w3.org/1999/xhtml">');
        say("\t\t\t<p>FORMULA: @metabolites[4]</p>");
        say("\t\t\t<p>CHARGE: @metabolites[3]</p>");
        say("\t\t</body>\n\t</notes>\n</species>\n");
	}
}

#finish metabolite list
print("\n");
say('</listOfSpecies>');
print("\n");
