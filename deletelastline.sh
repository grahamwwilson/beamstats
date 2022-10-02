#!/bin/sh

FN=$1

echo 'Deleting last line of file named '${FN}
echo 'Current file length: '
wc -l ${FN}
echo 'Last 3 lines of current file'
tail -3 ${FN}

head -n -1 ${FN} >temp.txt
mv temp.txt ${FN}

echo 'New file length'
wc -l ${FN}
echo 'Last 2 lines of current file'
tail -2 ${FN}

exit
