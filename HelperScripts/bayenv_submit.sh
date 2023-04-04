#/!/bin/bash

for d in ./run*
do
 (cd "$d" && ./calc_bf.sh 241forbayenv.txt 241stan4envirofile.txt 241covarmatrix_hwe.txt 5 500000 4)
done
