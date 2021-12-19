#!/bin/sh
#
# launchlumistats.sh
#

lumifile=$1
echo 'Using lumi file '${lumifile}
wc -l ${lumifile} >LumiFile_LineCount.dat

# Make symbolic link to the lumi file
if [ -L "lumifile.ini" ]
then
   rm lumifile.ini 
   ln -s ${lumifile} lumifile.ini
else
   ln -s ${lumifile} lumifile.ini
fi

# Run the already compiled executable

./checklumicorr

exit
