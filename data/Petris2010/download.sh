#!/bin/bash
# Download data files
URL=https://www.jstatsoft.org/index.php/jss/article/downloadSuppFile/v036i12
FILES="Folland.dat HL.dat IPCONGD.txt"

for f in $FILES
do
  wget $URL/$f
done
