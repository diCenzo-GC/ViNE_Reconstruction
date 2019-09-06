
#/bin/bash

rm  Mapped_in_Sino
cut -d'|' -f 1 Sino_boundary.mets | while read pat; do grep "$pat"  MappedCompoundsSino >> Mapped_in_Sino; if [[ `grep -wc "$pat" MappedCompoundsSino` > 1 ]]; then echo $pat; fi; done
