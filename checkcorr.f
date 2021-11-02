      program checkcorr
*
* Check statistics and correlations amongst initial beam file quantities
*
* Update to work better with different size files
*
*                Graham W. Wilson,    02-NOV-2021
*
      implicit none
      
      integer i,j,i1,i2
      integer N
      include 'arraysize.inc'
      double precision X(6,N)
      double precision xm(6),xr(6)
      double precision xmean,xrms
      integer nparticles
      
      print *,'Array size pre-dimensioned to N = ',N
      
* Check how many particles are available to read
      open(unit=11,file='BeamFile_LineCount.dat',status='OLD')      
      read(11,*)nparticles
      close(11)
      
      print *,'Trying to read ',nparticles,'particles'
      
      if(nparticles.le.N)then
         print *,'Array size is sufficient for requested particle count'
         print *,'This will WORK!'
      else
         print *,'Array size is too small '
         print *,'Update arraysize.inc appropriately and recompile'
         print *,'Will STOP!'
         stop
      endif
      
* Read initial beam file

      open(unit=21,file='beamfile.ini',status='OLD')      
      
      do i=1,nparticles
         read(21,*)x(1,i),x(2,i),x(3,i),x(4,i),x(5,i),x(6,i)
* check a couple of lines to guard against potential file reading errors        
         if(i.le.2.or.i.eq.nparticles)then
            print *,i,x(1,i),x(2,i),x(3,i),x(4,i),x(5,i),x(6,i)
         endif
      enddo
      
      close(21)
      
* Calculate statistics for the various "columns"      
      
      do j=1,6 
        call meanandrms(nparticles,x,j,xmean,xrms)
        xm(j) = xmean
        xr(j) = xrms
      enddo

      do i1=1,5
         do i2=i1+1,6
            call correlations(nparticles,x,i1,i2,xm,xr)
         enddo
      enddo
      
      end
      
      subroutine meanandrms(nparticles,x,j,xmean,xrms)
* Calculate mean and rms for quantity with index j      
      implicit none
      integer N
      include 'arraysize.inc'      
      double precision X(6,N)
      
      integer i,j,nparticles
      double precision xsum,xxsum,varx,sx
      double precision xmin,xmax
      double precision xmean,xrms
      character*9 cvalues(6)
      data cvalues/'E   (GeV)',
     +             'x    (um)',
     +             'y    (um)',
     +             'z    (um)',
     +             'xp (urad)',
     +             'yp (urad)'/  
      print *,' '    
      print *,'Statistics for ',cvalues(j)
      
      xsum = 0.0d0
      xxsum = 0.0d0
      xmin =  1.0d10
      xmax = -1.0d10
      
      do i=1,nparticles
         xsum = xsum + x(j,i)
         xxsum = xxsum + x(j,i)*x(j,i)
         if(x(j,i).gt.xmax)xmax = x(j,i)
         if(x(j,i).lt.xmin)xmin = x(j,i)
      enddo
      
      xsum = xsum/dble(nparticles)
      xxsum = xxsum/dble(nparticles)
      
      varx = xxsum - xsum*xsum
      sx = sqrt(varx)
      print *,'Mean     ',xsum
      print *,'Rms      ',sx
* Extrema values (both absolute and in normalized deviations)
      print *,'xmin     ',xmin,'    ',(xmin-xsum)/sx
      print *,'xmax     ',xmax,'    ',(xmax-xsum)/sx
      if(j.eq.1)then
         print *,'Rms/Mean ',sx/xsum
      endif
      
      xmean = xsum
      xrms  = sx
      
      end
      
      subroutine correlations(nparticles,x,i,j,xm,xr)
* Calculate correlation (if any) between i-th quantity and j-th quantity
      implicit none
      integer N
      include 'arraysize.inc'
      double precision X(6,N)
      double precision xm(6),xr(6)
      
      integer i,j,iev,nparticles
      double precision xysum,covxy,rho
      character*9 cvalues(6)
      data cvalues/'E   (GeV)',
     +             'x    (um)',
     +             'y    (um)',
     +             'z    (um)',
     +             'xp (urad)',
     +             'yp (urad)'/ 
     
      print *,' '
      
      xysum = 0.0d0
      
      do iev=1,nparticles
         xysum = xysum + x(i,iev)*x(j,iev)
      enddo
      xysum = xysum/dble(nparticles)
      
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
