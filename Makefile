FC:=gfortran

all: checkcorr.f checklumicorr.f checkgencorr.f checkmccorr.f arraysize.inc
	$(FC) -o checkcorr checkcorr.f
	$(FC) -o checklumicorr checklumicorr.f
	$(FC) -o checkgencorr checkgencorr.f
	$(FC) -o checkmccorr checkmccorr.f	
