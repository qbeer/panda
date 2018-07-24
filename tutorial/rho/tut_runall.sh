#!/bin/bash

nev=100
prefix="signal"
input="pp_jpsi2pi_jpsi_mumu.dec"
pbeam=6.231552

if test "$1" != ""; then
  nev=$1
fi

if test "$2" != ""; then
  prefix=$2
fi

if test "$3" != ""; then
  input=$3
fi

if test "$4" != ""; then
  pbeam=$4
fi

root -l -b -q  tut_sim.C\($nev,\"$prefix\",\"$input\",$pbeam\)
root -l -b -q  tut_aod.C\($nev,\"$prefix\"\)
