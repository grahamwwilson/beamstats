#!/bin/sh
#
# launchmcstats.sh
#

mcfile=$1
echo 'Using mcfile '${mcfile}

# Make symbolic link to the MC file
if [ -L "mcfile.ini" ]
then
   rm mcfile.ini 
   ln -s ${mcfile} mcfile.ini
else
   ln -s ${mcfile} mcfile.ini
fi

# Run the already compiled executable

./checkmccorr

exit
