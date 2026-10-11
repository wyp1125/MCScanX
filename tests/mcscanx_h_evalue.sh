#!/bin/sh
# MCScanX and MCScanX_h must report the same e-value for the same block.
set -e
cd "$(dirname "$0")/.."
T=$(mktemp -d)
trap 'rm -rf "$T"' EXIT
${CXX:-g++} struct.cc mcscan.cc read_data.cc out_utils.cc dagchainer.cc msa.cc permutation.cc -o "$T/MCScanX"
${CXX:-g++} struct.cc mcscan_h.cc read_homology.cc out_homology.cc dagchainer.cc msa.cc permutation.cc -o "$T/MCScanX_h"
for i in 1 2 3 4 5 6 7; do
  printf 'aa1\ta%d\t%d\t%d\n' $i ${i}000 ${i}500
  printf 'bb1\tb%d\t%d\t%d\n' $i ${i}000 ${i}500
done > "$T/t.gff"
printf 'a%d\tb%d\n' 1 1 2 2 3 3 4 4 5 5 6 6 > "$T/t.homology"
awk -v OFS='\t' '{print $1,$2,100,100,0,0,1,100,1,100,"1e-50",200}' "$T/t.homology" > "$T/t.blast"
"$T/MCScanX" "$T/t" -a -e 1 > /dev/null
x=$(grep -o 'e_value=[^ ]*' "$T/t.collinearity")
"$T/MCScanX_h" "$T/t" -a -e 1 > /dev/null
h=$(grep -o 'e_value=[^ ]*' "$T/t.collinearity")
if [ "$x" != "$h" ]; then
  echo "mcscanx_h_evalue: MCScanX $x, MCScanX_h $h"
  exit 1
fi
echo "mcscanx_h_evalue: ok ($x)"
