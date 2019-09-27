#!/bin/sh

for file in *.pdf; do
noext=`echo $file | sed 's#.pdf##'`
./pdf2eps.sh $noext
done
