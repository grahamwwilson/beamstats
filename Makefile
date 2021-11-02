FC:=gfortran

all: checkcorr.f arraysize.inc
	$(FC) -o checkcorr checkcorr.f
