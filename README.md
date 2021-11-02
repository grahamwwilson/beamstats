# beamstats
Beam file statistics

Read in beam file prepared for Guinea-PIG run and 
report statistics

# Compilation instructions
make   #Compile the code (needs gfortran)

# Running the code
./launchstats.sh input-beam-file.ini

# How this works
This runs the launchstats.sh unix shell script. 
The specified input file (eg. electron.ini) 
is symbolically linked to beamfile.ini. 
The information on (E, x, y, z, x', y') for each beam particle 
is read into an array and statistics on the mean and rms of each 
quantity and correlations are computed.

## Limitations
Currently the array size used is set 
to a maximum of 1,000,000 particles per file
