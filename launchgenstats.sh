#!/bin/sh
#
# launchgenstats.sh
#

genfile=$1
echo 'Using genfile '${genfile}
wc -l ${genfile} >GenFile_LineCount.dat

# Make symbolic link to the beam file
if [ -L "genfile.ini" ]
then
   rm genfile.ini 
   ln -s ${genfile} genfile.ini
else
   ln -s ${genfile} genfile.ini
fi

# Run the already compiled executable

./checkgencorr

exit
