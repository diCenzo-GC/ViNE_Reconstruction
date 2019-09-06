#Convert the metabolite list and the reaction list into a properly formatted metabolic reconstruction in SBML format.
#In the first script, change the name of the model accordlingly.

perl makeMetaboliteList.pl metaboliteList.txt > testModel.xml
perl makeReactionList.pl reactionList.txt >> testModel.xml

