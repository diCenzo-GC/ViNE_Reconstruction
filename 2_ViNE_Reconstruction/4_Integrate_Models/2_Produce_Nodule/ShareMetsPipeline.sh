
#/bin/bash
#  a bash pipeline to get shared compounds between sino boundary mets and truncatula cellular (cytosolic mets)
# initial compounds lists are in Sino_boundary.mets and Truncatula_Cytosol.mets
# MNX reference file is chem_xref.tsv
# final output file is in Mapped_in_shared

echo
echo "*** Starting bash pipeline for getting shared mets between sino and truncatula ***"
echo

# Get MNX notation for boundary sino mets and cellular truncatula mets
echo ">Get MNX notation for boundary sino mets and cellular truncatula mets"
./GetMappedMets_sino.sh
./GetMappedMets_trunca.sh
echo done
echo

echo ">Identify shared mets between the two"
# Identify shared mets between the two
./GetShared_MNX.sh
echo done
echo

echo ">Get functional annotation of the shared mets"
# Get functional annotation of the shared boundary mets
./GetMappedMets_shared.sh
echo done
echo


echo ">Get sino unshared mets"
# Get sinorhizobium boundary mets that are not shared with the plant
./GetUnsharedSino.sh
echo done
echo

echo ">Get sino met names of shared mets"
# Get the  sino names for the metabolites that are shared with the plant
./GetMappingOfSharedSino.sh
echo done
echo

echo ">Get truncatula met names of shared mets"
# do the same for the trunca model
./GetMappingOfSharedTrunca.sh
echo done
echo

echo ">Prepare new namings for importing in matlab -sino model"
# create two arrays with new and old names for use in Matlab
perl PrepareMetsSubstituion_Sino.pl SinoNamesMapping_of_shared
echo done
echo

echo ">Prepare new namings for importing in Matlab -truncatula model"
# create two arrays with new and old names for use in Matlab 
perl PrepareMetsSubstituion_Trunca.pl TruncaNamesMapping_of_shared
echo
echo
echo "*** Done, end of bash pipeline, going back to matlab code ***"
echo 
echo
