#/bin/bash

grep -vxF -f shared_MNXM Sino_MNXM  > unshared_bounday_sino

rm Mapped_in_unshared;
cut -d'|' -f 1 unshared_bounday_sino | while read pat; do grep -w "$pat" chem_prop.tsv >> Mapped_in_unshared; if [[ `grep -cw "$pat" shared_MNXM` > 1 ]]; then echo $pat; fi; done


