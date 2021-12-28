      program checkmccorr
*
* Check statistics and correlations amongst 
* (x1, x2) quantities in the reweighting MC toy data-sets
*
* Use similar logic to ~/BeamSpectrumGenerator/rwbingp.f
*
*                Graham W. Wilson,    27-DEC-2021
*
      implicit none
      
      integer i,j,i1,i2
      integer N
      include 'arraysize.inc'
      integer npar
      parameter (npar=2)
      double precision X(npar,N)
      double precision xm(npar),xr(npar)
      double precision xmean,xrms
      
      integer ntoread
      integer nevents
      integer nread
      
      double precision ebeam
      parameter (ebeam = 125.0d0)
      
      integer iev,rtype,ibin1,ibin2,icomb
      double precision x1,x2,y1,y2,z1,z2,weight,weightp,rwt
      
      integer lun
      parameter (lun = 22)
      
      print *,'Array size pre-dimensioned to N = ',N
      
* Initialize input file
      open(unit=lun,file='mcfile.ini',status='OLD')

* Read the header
      call readheader(lun,ntoread)
      
      print *,'Trying to read ',ntoread,' events'
      
      if(ntoread.le.N)then
         print *,'Array size is sufficient for requested event count'
         print *,'This will WORK!'
      else
         print *,'Array size is too small '
         print *,'Update arraysize.inc appropriately and recompile'
         print *,'Will STOP!'
         stop
      endif
      
      nevents = 0
      
      do i=1,ntoread
         read(lun,*)iev,rtype,x1,x2,y1,y2,z1,z2,
     +              weight,weightp,rwt,ibin1,ibin2,icomb      
         nevents = nevents + 1           
         x(1,i) = x1
         x(2,i) = x2
      enddo
      
      close(lun)
      
      print *,'nevents analyzed ',nevents
      
* Calculate statistics for the various "columns"      
      
      do j=1,npar
        call meanandrms(nevents,x,j,xmean,xrms)
        xm(j) = xmean
        xr(j) = xrms
        if(j.le.2)then
*           print *,'Mean energy loss relative to ',Ebeam,' of ',
*     +              (Ebeam-xm(j))/Ebeam
           print *,'Mean energy loss relative to ',Ebeam,' of ',
     +              (1.0d0-xm(j))
        endif
      enddo

      do i1=1,1
         do i2=i1+1,npar
            call correlations(nevents,x,i1,i2,xm,xr)
         enddo
      enddo
      
      end
            
      subroutine meanandrms(nevents,x,j,xmean,xrms)
* Calculate mean and rms for quantity with index j      
      implicit none
      integer N
      include 'arraysize.inc'

      integer npar
      parameter (npar=2)            
      double precision X(npar,N)
      
      integer i,j,nevents
      double precision xsum,xxsum,varx,sx
      double precision xmin,xmax
      double precision xmean,xrms
      character*9 cvalues(npar)
      data cvalues/'E1/Eb    ',
     +             'E2/Eb    '/
       
      print *,' '    
      print *,'Statistics for ',cvalues(j)
      
      xsum = 0.0d0
      xxsum = 0.0d0
      xmin =  1.0d10
      xmax = -1.0d10
      
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
      if(j.le.2.or.j.eq.6)then
         print *,'Rms/Mean ',sx/xsum
      endif
      
      xmean = xsum
      xrms  = sx
      
      end
      
      subroutine correlations(nevents,x,i,j,xm,xr)
* Calculate correlation (if any) between i-th quantity and j-th quantity
      implicit none
      integer N
      include 'arraysize.inc'
      
      integer npar
      parameter (npar=2)            
      double precision X(npar,N)      
      
      double precision xm(2),xr(2)
      integer i,j,iev,nevents
      double precision xysum,covxy,rho
      
      character*9 cvalues(npar)
      data cvalues/'E1/Eb    ',
     +             'E2/Eb    '/     
      print *,' '
      
      xysum = 0.0d0
      
      do iev=1,nevents
         xysum = xysum + x(i,iev)*x(j,iev)
      enddo
      xysum = xysum/dble(nevents)
      
* Calculate covariance and correlation coefficient 
* making use of precomputed means and rms of each variable
      covxy = xysum - xm(i)*xm(j)
      rho = 0.0d0
      if(xr(i).gt.0.0d0.and.xr(j).gt.0.0d0)then     
         rho = covxy/(xr(i)*xr(j))
      endif
      
*      print *,'Covariance ',covxy
*      print *,'Vari       ',xr(i)*xr(i)
*      print *,'Varj       ',xr(j)*xr(j)

* Correlation Coefficient
*      print *,'Rho        ',rho
      print *,'Correlation Coefficient for ',
     + cvalues(i),' with ',cvalues(j),' of ',rho
      
      end
      
      include 'readheader.f'
