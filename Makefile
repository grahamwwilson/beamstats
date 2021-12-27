FC:=gfortran

all: checkcorr.f checklumicorr.f checkgencorr.f arraysize.inc
	$(FC) -o checkcorr checkcorr.f
	$(FC) -o checklumicorr checklumicorr.f
	$(FC) -o checkgencorr checkgencorr.f		
