      program checklumicorr
*
* Check statistics and correlations amongst 
* quantities in the lumi.ee.out files produced by Guinea-PIG.
*
* Update to work better with different size files
*
*                Graham W. Wilson,    18-DEC-2021
*
      implicit none
      
      integer i,j,i1,i2
      integer N
      include 'arraysize.inc'
      integer npar
      parameter (npar=17)
      double precision X(npar,N)
      double precision xm(npar),xr(npar)
      double precision xmean,xrms
      
      integer ntoread
      integer nevents
      integer nread
      integer iev
      
      double precision ebeam
      parameter (ebeam = 125.0d0)
      
      double precision x1,x2,xpv,ypv,zpv,x1p,y1p,x2p,y2p
      integer itime,sx1,sy1,sz1,sx2,sy2,sz2,iorder
      
      print *,'Array size pre-dimensioned to N = ',N
      
* Check how many particles are available to read
      open(unit=11,file='LumiFile_LineCount.dat',status='OLD')      
      read(11,*)ntoread
      close(11)
      
      print *,'Trying to read ',ntoread,' events'
      
      if(ntoread.le.N)then
         print *,'Array size is sufficient for requested particle count'
         print *,'This will WORK!'
      else
         print *,'Array size is too small '
         print *,'Update arraysize.inc appropriately and recompile'
         print *,'Will STOP!'
         stop
      endif
      
* Read initial beam file

      open(unit=21,file='lumifile.ini',status='OLD')      
      
      nevents = 0
      iev = 0
      
      do i=1,ntoread
         read(21,*)x1,x2,xpv,ypv,zpv,itime,x1p,y1p,x2p,y2p,
     +             sx1,sy1,sz1,sx2,sy2,sz2,iorder
* Rescale to Ebeam
         x1 = x1/ebeam
         x2 = x2/ebeam
* Various potential cuts - select one and only one condition
* to store the event. FIXME - write this more elegantly.
*         if(zpv.ge.0.0d0)then
*         if(zpv.lt.0.0d0)then
         if(abs(zpv).ge.0.0d0)then
         iev = iev + 1
         nevents = nevents + 1           
         call myfill(iev,x,x1,x2,xpv,ypv,zpv,itime,x1p,y1p,x2p,y2p,
     +             sx1,sy1,sz1,sx2,sy2,sz2,iorder)
         endif
* check a couple of lines to guard against potential file reading errors        
*         if(i.le.2.or.i.eq.nevents)then
*            print *,i,x(1,i),x(2,i),x(3,i),x(4,i),x(5,i),x(6,i)
*         endif
      enddo
      
      close(21)
      
      print *,'nevents analyzed ',nevents
      
* Calculate statistics for the various "columns"      
      
      do j=1,npar
        if(j.le.10.or.j.eq.npar)then 
           call meanandrms(nevents,x,j,xmean,xrms)
        endif
        xm(j) = xmean
        xr(j) = xrms
        if(j.le.2)then
*           print *,'Mean energy loss relative to ',Ebeam,' of ',
*     +              (Ebeam-xm(j))/Ebeam
           print *,'Mean energy loss relative to ',Ebeam,' of ',
     +              (1.0d0-xm(j))
        endif
      enddo

      do i1=1,10
         do i2=i1+1,npar
            if(i2.lt.11.or.i2.eq.17)then
               call correlations(nevents,x,i1,i2,xm,xr)
            endif
         enddo
      enddo
      
      end
      
      subroutine myfill(i,x,x1,x2,xpv,ypv,zpv,itime,x1p,y1p,x2p,y2p,
     +             sx1,sy1,sz1,sx2,sy2,sz2,iorder)
      implicit none
      integer i
      
      integer N
      include 'arraysize.inc'
      integer npar
      parameter (npar=17)
      double precision X(npar,N)      
      
      double precision x1,x2,xpv,ypv,zpv,x1p,y1p,x2p,y2p
      integer itime,sx1,sy1,sz1,sx2,sy2,sz2,iorder      
      
      X(1,i) = x1
      X(2,i) = x2
      X(3,i) = xpv
      X(4,i) = ypv
      X(5,i) = zpv
      X(6,i) = dble(itime)
      X(7,i) = x1p
      X(8,i) = y1p
      X(9,i) = x2p
      X(10,i) = y2p
      X(11,i) = dble(sx1)
      X(12,i) = dble(sy1)      
      X(13,i) = dble(sz1)      
      X(14,i) = dble(sx2)
      X(15,i) = dble(sy2)      
      X(16,i) = dble(sz2)
      X(17,i) = dble(iorder)
      
      end
      
      subroutine meanandrms(nevents,x,j,xmean,xrms)
* Calculate mean and rms for quantity with index j      
      implicit none
      integer N
      include 'arraysize.inc'

      integer npar
      parameter (npar=17)            
      double precision X(npar,N)
      
      integer i,j,nevents
      double precision xsum,xxsum,varx,sx
      double precision xmin,xmax
      double precision xmean,xrms
      character*9 cvalues(npar)
      data cvalues/'E1/Eb    ',
     +             'E2/Eb    ',
     +             'xpv  (nm)',
     +             'ypv  (nm)',
     +             'zpv  (um)',
     +             'time  (?)',
     +             'x1p (rad)',
     +             'y1p (rad)',
     +             'x2p (rad)',
     +             'y2p (rad)',
     +             'sx1      ',
     +             'sy1      ',
     +             'sz1      ',
     +             'sx2      ',
     +             'sy2      ',
     +             'sz2      ',
     +             'order    '/
       
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
      parameter (npar=17)            
      double precision X(npar,N)      
      
      double precision xm(6),xr(6)
      integer i,j,iev,nevents
      double precision xysum,covxy,rho
      
      character*9 cvalues(npar)
      data cvalues/'E1/Eb    ',
     +             'E2/Eb    ',
     +             'xpv  (nm)',
     +             'ypv  (nm)',
     +             'zpv  (um)',
     +             'time  (?)',
     +             'x1p (rad)',
     +             'y1p (rad)',
     +             'x2p (rad)',
     +             'y2p (rad)',
     +             'sx1      ',
     +             'sy1      ',
     +             'sz1      ',
     +             'sx2      ',
     +             'sy2      ',
     +             'sz2      ',
     +             'order    '/      
     
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
