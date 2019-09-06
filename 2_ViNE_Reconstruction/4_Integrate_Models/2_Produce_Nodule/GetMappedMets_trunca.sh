
#/bin/bash

rm Mapped_in_Trunca
cut -d'|' -f 1 Truncatula_Cytosol.mets | while read pat; do grep -w "$pat" MappedCompoundsTrunca >> Mapped_in_Trunca; if [[ `grep -cw "$pat" MappedCompoundsTrunca` > 1 ]]; then echo $pat; fi; done
