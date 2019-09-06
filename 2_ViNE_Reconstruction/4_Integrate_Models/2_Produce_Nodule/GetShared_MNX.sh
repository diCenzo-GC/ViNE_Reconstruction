#!/bin/bash

awk '{print $1}' Mapped_in_Sino > Sino_MNXM
awk '{print $1}' Mapped_in_Trunca > Truncatula_MNXM

cat Sino_MNXM | sort > one_list.tmp
cat Truncatula_MNXM | sort > two_list.tmp


comm -12 one_list.tmp two_list.tmp > shared_MNXM

rm *.tmp
