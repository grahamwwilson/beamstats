#!/bin/sh
#
# launchstats.sh
#

beamfile=$1
echo 'Using beamfile '${beamfile}
wc -l ${beamfile} >BeamFile_LineCount.dat

# Make symbolic link to the beam file
if [ -L "beamfile.ini" ]
then
   rm beamfile.ini 
   ln -s ${beamfile} beamfile.ini
else
   ln -s ${beamfile} beamfile.ini
fi

# Run the already compiled executable

./checkcorr

exit
