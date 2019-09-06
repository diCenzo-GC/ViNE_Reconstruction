#/bin/bash
rm TruncaNamesMapping_of_shared
cut -d'|' -f 1 shared_MNXM | while read pat; do grep -w "$pat" Mapped_in_Trunca >> TruncaNamesMapping_of_shared ; if [[ `grep -cw "$pat" shared_MNXM` > 1 ]]; then echo $pat; fi; done

