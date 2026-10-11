#!/bin/sh
# dag_main must not index score past its end after print_chains shrinks it.
# _GLIBCXX_ASSERTIONS turns the out-of-bounds operator[] into an abort.
set -e
cd "$(dirname "$0")/.."
T=$(mktemp -d)
trap 'rm -rf "$T"' EXIT
${CXX:-g++} -D_GLIBCXX_ASSERTIONS struct.cc mcscan_h.cc read_homology.cc out_homology.cc dagchainer.cc msa.cc permutation.cc -o "$T/MCScanX_h"
for i in 1 2 3 4 5 6 7; do
  printf 'aa1\ta%d\t%d\t%d\n' $i ${i}000 ${i}500
  printf 'bb1\tb%d\t%d\t%d\n' $i ${i}000 ${i}500
done > "$T/t.gff"
printf 'a%d\tb%d\n' 1 1 2 2 3 3 4 4 5 5 6 6 6 7 > "$T/t.homology"
"$T/MCScanX_h" "$T/t" -a > /dev/null
echo "dag_main_bounds: ok"
