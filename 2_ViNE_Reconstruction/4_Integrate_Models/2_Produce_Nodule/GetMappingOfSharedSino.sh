#/bin/bash
rm SinoNamesMapping_of_shared
cut -d'|' -f 1 shared_MNXM | while read pat; do grep -w "$pat" Mapped_in_Sino >> SinoNamesMapping_of_shared ; if [[ `grep -cw "$pat" shared_MNXM` > 1 ]]; then echo $pat; fi; done

