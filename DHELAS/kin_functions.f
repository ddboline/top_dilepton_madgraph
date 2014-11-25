c************************************************************************
c  THIS FILE CONTAINS THE DEFINITIONS OF USEFUL FUNCTIONS OF MOMENTA:
c  
c  DOT(p1,p2)         : 4-Vector Dot product
c  R2(p1,p2)          : distance in eta,phi between two particles
c  SumDot(P1,P2,dsign): invariant mass of 2 particles
c  eta(p)             : rapidity of particle in the lab frame
c  DELTA_PHI(P1, P2)  : separation in phi of two particles 
c  ET(p)              : transverse energy of particle
c  PT(p)              : transverse momentum of particle
c  DJ(p1,p2)          : y (Durham) value for two partons
c  DJ2(p1,p2)         : scalar product squared
c************************************************************************

      DOUBLE PRECISION FUNCTION R2(P1,P2)
c************************************************************************
c     Distance in eta,phi between two particles.
c************************************************************************
      IMPLICIT NONE
c
c     Arguments
c
      double precision p1(0:3),p2(0:3)
c
c     External
c
      double precision ETA,DELTA_PHI
      external eta,delta_phi
c-----
c  Begin Code
c-----
      R2 = (DELTA_PHI(P1,P2))**2+(ETA(p1)-ETA(p2))**2
      RETURN
      END

      DOUBLE PRECISION FUNCTION SumDot(P1,P2,dsign)
c************************************************************************
c     Invarient mass of 2 particles
c************************************************************************
      IMPLICIT NONE
c
c     Arguments
c
      double precision p1(0:3),p2(0:3),dsign
c
c     Local
c      
      integer i
      double precision ptot(0:3)
c
c     External
c
      double precision dot
      external dot
c-----
c  Begin Code
c-----

      do i=0,3
         ptot(i)=p1(i)+dsign*p2(i)
      enddo
      SumDot = dot(ptot,ptot)
      RETURN
      END

      DOUBLE PRECISION  FUNCTION eta(p)
c************************************************************************
c     Returns rapidity of particle in the lab frame
c************************************************************************
      IMPLICIT NONE
c
c     Arguments
c
      double precision  p(0:3)
c
c     Local
c
      double precision pm
c
c     Global
c
      integer                                        lpp(2)
      double precision    ebeam(2), xbk(2),q2fact(2)
      common/to_collider/ ebeam   , xbk   ,q2fact,   lpp
c-----
c  Begin Code
c-----
c      pm=dsqrt(p(1)**2+p(2)**2+p(3)**2)
      pm = p(0)
      eta = .5d0*dlog((pm+p(3))/(pm-p(3)))+
     $     .5d0*dlog(xbk(1)*ebeam(1)/(xbk(2)*ebeam(2)))
      end

      DOUBLE PRECISION FUNCTION DELTA_PHI(P1, P2)
c************************************************************************
c     Returns separation in phi of two particles p1,p2
c************************************************************************
      IMPLICIT NONE
c
c     Arguments
c
      double precision p1(0:3),p2(0:3)
c
c     Local
c
      REAL*8 DENOM, TEMP
c-----
c  Begin Code
c-----
      DENOM = SQRT(P1(1)**2 + P1(2)**2) * SQRT(P2(1)**2 + P2(2)**2)
      TEMP = MAX(-0.99999999D0, (P1(1)*P2(1) + P1(2)*P2(2)) / DENOM)
      TEMP = MIN( 0.99999999D0, TEMP)
      DELTA_PHI = ACOS(TEMP)
      END



      double precision function et(p)
c************************************************************************
c     Returns transverse energy of particle
c************************************************************************
      IMPLICIT NONE
c
c     Arguments
c
      double precision p(0:3)
c
c     Local
c
      double precision pt
c-----
c  Begin Code
c-----
      pt = dsqrt(p(1)**2+p(2)**2)
      if (pt .gt. 0d0) then
         et = p(0)*pt/dsqrt(pt**2+p(3)**2)
      else
         et = 0d0
      endif
      end

      double precision function pt(p)
c************************************************************************
c     Returns transverse momentum of particle
c************************************************************************
      IMPLICIT NONE
c
c     Arguments
c
      double precision p(0:3)
c-----
c  Begin Code
c-----

      pt = dsqrt(p(1)**2+p(2)**2)

      return
      end

      double precision function DJ(p1,p2)
c***************************************************************************
c     Uses Durham algorythm to calculate the y value for two partons
c***************************************************************************
      implicit none
c
c     Arguments
c
      double precision p1(0:3),p2(0:3)
c
c     Local
c
      double precision pm1,pm2,costh
      integer j
c-----
c  Begin Code
c-----
      pm1 = dsqrt(p1(1)**2+p1(2)**2+p1(3)**2)
      pm2 = dsqrt(p2(1)**2+p2(2)**2+p2(3)**2)
      if (pm1*pm2 .ne. 0d0) then
         costh = (p1(1)*p2(1)+p1(2)*p2(2)+p1(3)*p2(3))/(pm1*pm2)
         dj = 2d0*min(p1(0)**2,p2(0)**2)*(1d0-costh)   !Durham
c         dj = 2d0*p1(0)*p2(0)*(1d0-costh)    !JADE
      else
         print*,'Warning 0 momentum in Durham algorythm'
         write(*,'(4e15.5)') (p1(j),j=0,3)
         write(*,'(4e15.5)') (p2(j),j=0,3)
         dj = 0d0
      endif
      end

      double precision function DJ2(p1,p2)
c***************************************************************************
c     Uses Lorentz
c***************************************************************************
      implicit none
c
c     Arguments
c
      double precision p1(0:3),p2(0:3)
c
c     Local
c
      integer j
c
c     External
c
      double precision dot
c-----
c  Begin Code
c-----
      dj2 = dot(p1,p1)+2d0*dot(p1,p2)+dot(p2,p2)
      return
      end

      subroutine switchmom(p1,p,ic,jc,nexternal)
c**************************************************************************
c     Changes stuff for crossings
c**************************************************************************
      implicit none
      integer nexternal
      integer jc(nexternal),ic(nexternal)
      real*8 p1(0:3,nexternal),p(0:3,nexternal)
      integer i,j
c-----
c Begin Code
c-----
      do i=1,nexternal
         do j=0,3
            p(j,ic(i))=p1(j,i)
         enddo
      enddo
      do i=1,nexternal
         jc(i)=1
      enddo
      jc(ic(1))=-1
      jc(ic(2))=-1
      end

      double precision function dot(p1,p2)
C****************************************************************************
C     4-Vector Dot product
C****************************************************************************
      implicit none
      double precision p1(0:3),p2(0:3)
      dot=p1(0)*p2(0)-p1(1)*p2(1)-p1(2)*p2(2)-p1(3)*p2(3)
      end

