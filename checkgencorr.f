      program checkgencorr
*
* Check statistics and correlations amongst 
* DumpBeamEnergies.dat quantities on the ILD DSTs.
*         gplikefile << mcpElectronEnergy[2] << " " << mcpElectronEnergy[3] << " " << *mcpAfterBSM << " " << *mcpAfterISRM << " " << *mcpFinalStateM << " " << *mcpMuonicM << endl;
*
*                Graham W. Wilson,    20-DEC-2021
*
      implicit none
      
      integer i,j,i1,i2
      integer N
      include 'arraysize.inc'
      double precision X(6,N)
      double precision xm(6),xr(6)
      double precision xmean,xrms
      integer nevents
      
      print *,'Array size pre-dimensioned to N = ',N
      
* Check how many particles are available to read
      open(unit=11,file='GenFile_LineCount.dat',status='OLD')      
      read(11,*)nevents
      close(11)
      
      print *,'Trying to read ',nevents,'events'
      
      if(nevents.le.N)then
         print *,'Array size is sufficient for requested particle count'
         print *,'This will WORK!'
      else
         print *,'Array size is too small '
         print *,'Update arraysize.inc appropriately and recompile'
         print *,'Will STOP!'
         stop
      endif
      
* Read initial beam file

      open(unit=21,file='genfile.ini',status='OLD')      
      
      do i=1,nevents
         read(21,*)x(1,i),x(2,i),x(3,i),x(4,i),x(5,i),x(6,i)
* check a couple of lines to guard against potential file reading errors        
         if(i.le.2.or.i.eq.nevents)then
            print *,i,x(1,i),x(2,i),x(3,i),x(4,i),x(5,i),x(6,i)
         endif
      enddo
      
      close(21)
      
* Calculate statistics for the various "columns"      
      
      do j=1,6 
        call meanandrms(nevents,x,j,xmean,xrms)
        xm(j) = xmean
        xr(j) = xrms
      enddo

      do i1=1,5
         do i2=i1+1,6
            call correlations(nevents,x,i1,i2,xm,xr)
         enddo
      enddo
      
      end
      
      subroutine meanandrms(nevents,x,j,xmean,xrms)
* Calculate mean and rms for quantity with index j 
*         gplikefile << mcpElectronEnergy[2] << " " << mcpElectronEnergy[3] << " " << *mcpAfterBSM << " " << *mcpAfterISRM << " " << *mcpFinalStateM << " " << *mcpMuonicM << endl;     
      implicit none
      integer N
      include 'arraysize.inc'      
      double precision X(6,N)
      
      integer i,j,nevents
      double precision xsum,xxsum,varx,sx
      double precision xmin,xmax
      double precision xmean,xrms
* Nominal beam energy
      double precision ebeam
      parameter (ebeam = 125.0d0)
* Nominal half crossing-angle
      double precision ha
      double precision enorm(3)
      
      parameter (ha=0.007d0)
      character*13 cvalues(6)
      data cvalues/'E1      (GeV)',
     +             'E2      (GeV)',
     +             'ECM     (GeV)',
     +             'ECMI    (GeV)',
     +             'M-FS    (GeV)',
     +             'M-mumu  (GeV)'/  
      print *,' '    
      print *,'Statistics for ',cvalues(j)
      
      xsum = 0.0d0
      xxsum = 0.0d0
      xmin =  1.0d10
      xmax = -1.0d10
      
      enorm(1) = ebeam/cos(ha)
      enorm(2) = ebeam/cos(ha)
      enorm(3) = 2.0*ebeam
      
      do i=1,nevents
         xsum = xsum + x(j,i)
         xxsum = xxsum + x(j,i)*x(j,i)
         if(x(j,i).gt.xmax)xmax = x(j,i)
         if(x(j,i).lt.xmin)xmin = x(j,i)
      enddo
      
      xsum = xsum/dble(nevents)
      xxsum = xxsum/dble(nevents)
      
      varx = xxsum - xsum*xsum
      sx = sqrt(varx)
      print *,'Mean     ',xsum
      print *,'Rms      ',sx
* Extrema values (both absolute and in normalized deviations)
      print *,'xmin     ',xmin,'    ',(xmin-xsum)/sx
      print *,'xmax     ',xmax,'    ',(xmax-xsum)/sx
      if(j.le.3)then
         print *,'Rms/Mean ',sx/xsum,'Nom: ',enorm(j),
     +           ' Mean/Nom ',xsum/enorm(j),
     +           ' (1-(mean/nom))',1.0d0-(xsum/enorm(j)),
     +           'Rms/Norm ',sx/enorm(j)
      endif
      
      xmean = xsum
      xrms  = sx
      
      end
      
      subroutine correlations(nevents,x,i,j,xm,xr)
* Calculate correlation (if any) between i-th quantity and j-th quantity
      implicit none
      integer N
      include 'arraysize.inc'
      double precision X(6,N)
      double precision xm(6),xr(6)
      
      integer i,j,iev,nevents
      double precision xysum,covxy,rho
      character*13 cvalues(6)
      data cvalues/'E1      (GeV)',
     +             'E2      (GeV)',
     +             'ECM     (GeV)',
     +             'ECMI    (GeV)',
     +             'M-FS    (GeV)',
     +             'M-mumu  (GeV)'/
     
      print *,' '
      
      xysum = 0.0d0
      
      do iev=1,nevents
         xysum = xysum + x(i,iev)*x(j,iev)
      enddo
      xysum = xysum/dble(nevents)
      
* Calculate covariance and correlation coefficient 
* making use of precomputed means and rms of each variable
      covxy = xysum - xm(i)*xm(j)
      rho = covxy/(xr(i)*xr(j))
      
*      print *,'Covariance ',covxy
*      print *,'Vari       ',xr(i)*xr(i)
*      print *,'Varj       ',xr(j)*xr(j)

* Correlation Coefficient
*      print *,'Rho        ',rho
      print *,'Correlation Coefficient for ',
     + cvalues(i),' with ',cvalues(j),' of ',rho
      
      end
