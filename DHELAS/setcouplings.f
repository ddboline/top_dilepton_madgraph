c
      SUBROUTINE SETCOUPLINGS()
c
c******************************************************************
c     sets up masses and coupling constants for Helas
c******************************************************************
      implicit none
c
c Constants
c
      double precision  sw2
      parameter        (sw2 = .2220d0)
c
c   masses and widths of fermions
c
      double precision tmass,      bmass,     cmass
      parameter       (tmass=174.3d0,bmass=4.7d0, cmass=1.420d0)
      double precision smass,      umass,     dmass
      parameter       (smass=0d0,  umass=0d0, dmass=0d0)
      double precision twidth,     bwidth,    cwidth
      parameter       (twidth=1.508d0, bwidth=0d0,cwidth=0d0)
      double precision swidth,     uwidth,    dwidth
      parameter       (swidth=0d0, uwidth=0d0,dwidth=0d0)
      double precision emass,      mumass,    taumass
      parameter       (emass=0d0,  mumass=0d0,taumass=1.777d0)
      double precision ewidth,     muwidth,    tauwidth
      parameter       (ewidth=0d0, muwidth=0d0,tauwidth=0d0)
c
c   masses and widths of bosons
c
c
c Local Variables
c
      integer i
c
c Global Variables
c
c
C GLOBAL VARIABLES
C  
      REAL*8         alpha,fourpi,ee,sw,cw,sc2,v,ee2,als
      COMMON /COUPBASIC/ alpha,ee, sw, cw, als   
      REAL*8         GW, GWWA, GWWZ
      COMMON /COUP1/ GW, GWWA, GWWZ
      REAL*8         GAL(2),GAU(2),GAD(2),GWF(2)
      COMMON /COUP2A/GAL,   GAU,   GAD,   GWF
      REAL*8         GZN(2),GZL(2),GZU(2),GZD(2),G1(2)
      COMMON /COUP2B/GZN,   GZL,   GZU,   GZD,   G1
      REAL*8         GWWH,GZZH,GHHH,GWWHH,GZZHH,GHHHH
      COMMON /COUP3/ GWWH,GZZH,GHHH,GWWHH,GZZHH,GHHHH
      COMPLEX*16     GH(2,12)
      COMMON /COUP4/ GH
      REAL*8         WMASS,WWIDTH,ZMASS,ZWIDTH
      COMMON /VMASS1/WMASS,WWIDTH,ZMASS,ZWIDTH
      REAL*8         AMASS,AWIDTH,HMASS,HWIDTH
      COMMON /VMASS2/AMASS,AWIDTH,HMASS,HWIDTH
      REAL*8           GG(2), G
      COMMON /COUPQCD/ GG,    G
      REAL*8            FMASS(12), FWIDTH(12)
      COMMON /FERMIONS/ FMASS,     FWIDTH
c
c
      REAL*8 ez, ey
      REAL*8 r_zero, r_half, r_one, r_two, r_three, r_four, r_ote
      REAL*8 r_pi, r_ialph, r_alpha
      PARAMETER( r_zero=0.0d0, r_half=0.5d0, r_one=1.0d0, r_two=2.0d0 )
      PARAMETER( r_three=3.0d0 )
      PARAMETER( r_four=4.0d0, r_ote=128.0d0 )
      PARAMETER( r_alpha=132.50698 )
      PARAMETER( r_pi=3.14159265358979323846d0, r_ialph=137.0359895d0 )
c
c-----------
c Begin Code
c-----------
c      Include 'constants.inc'
c          
      GW=0
      GWWA=0
      GWWZ=0
      GAL(1)=0
      GAL(2)=0
      GAU(1)=0
      GAU(2)=0
      GAD(1)=0
      GAD(2)=0
      GWF(1)=0
      GWF(2)=0
      GZN(1)=0
      GZN(2)=0
      GZL(1)=0
      GZL(2)=0
      GZU(1)=0
      GZU(2)=0
      GZD(1)=0
      GZD(2)=0
      G1(1)=0
      G1(2)=0
      WMASS=0
      WWIDTH=0
      ZMASS=0
      ZWIDTH=0
      HMASS=0
      HWIDTH=0
      DO I=0,12
         FMASS(I)=0
         FWIDTH(I)=0
      ENDDO
c

      fmass(1) = emass
      fmass(2) = 0d0
      fmass(3) = umass
      fmass(4) = dmass
      fmass(5) = mumass
      fmass(6) = 0d0
      fmass(7) = cmass
      fmass(8) = smass
      fmass(9) = taumass
      fmass(10)= 0d0
      fmass(11)= tmass
      fmass(12)= bmass
c
      fwidth(1) = ewidth
      fwidth(2) = 0d0
      fwidth(3) = uwidth
      fwidth(4) = dwidth
      fwidth(5) = muwidth
      fwidth(6) = 0d0
      fwidth(7) = cwidth
      fwidth(8) = swidth
      fwidth(9) = tauwidth
      fwidth(10)= 0d0
      fwidth(11)= twidth
      fwidth(12)= bwidth
c
c      wmass1=wm
c      zmass1=zm
c      amass1=amass
c      hmass1=hm
c      wwidth1=ww
c      zwidth1=zw
c      awidth1=awidth
c      hwidth1=hw
c
c  Calls to Helas routines to set couplings
c
C  
c
c
c

c
      WMASS=80.419
      WWIDTH=2.087
      ZMASS=91.1888
      ZWIDTH=2.4974
      HMASS=120.0
      HWIDTH=0.001
c
c
c
C
c
      alpha = r_alpha
c     alpha = r_one / r_ialph
      fourpi = r_four * r_pi
      ee=sqrt( alpha * fourpi )
      ee2=alpha*fourpi
      sw=sqrt( sw2 )
      cw=sqrt( r_one - sw2 )
      sc2=sw2*( r_one - sw2 )
      v = r_two * zmass*sqrt(sc2)/sqrt(ee2)
c
cc
      GW    =  ee/sw
      GWWA  =  ee
      GWWZ  =  ee*cw/sw
c     WRITE(*,*) GW, GWWA, GWWZcccc
c
      ez=ee/(sw*cw)
      ey=ee*(sw/cw)
c
      GAL(1) =  ee
      GAL(2) =  ee
      GAU(1) = -ee*r_two/r_three
      GAU(2) = -ee*r_two/r_three
      GAD(1) =  ee   /r_three
      GAD(2) =  ee   /r_three
      GWF(1) = -ee/sqrt(r_two*sw2)
      GWF(2) =  r_zero
      GZN(1) = -ez*  r_half
      GZN(2) =  r_zero

      GZL(1) = -ez*(-r_half+sw2)
      GZL(2) = -ey

c      write(*,*) 'ez: ' , ez
c      write(*,*) 'r_half: ', r_half
c      write(*,*) 'sin2w: ', sw2
c      write(*,*) 'ey: ', ey
c      write(*,*) 'sw2*2/3: ', sw2*r_two/r_three
      GZU(1) = -ez*( r_half-sw2*r_two/r_three)
      GZU(2) =  ey*          r_two/r_three


c      write(*,*) 'gzu: ', gzu
      GZD(1) = -ez*(-r_half+sw2   /r_three)
      GZD(2) = -ey             /r_three
c
      GWWH  =   ee2/sw2*r_half*v
      GZZH  =   ee2/sc2*r_half*v
      GHHH  =  -hmass**2/v*r_three
      GWWHH =   ee2/sw2*r_half
      GZZHH =   ee2/sc2*r_half
      GHHHH = -(hmass/v)**2*r_three
c
C
      EZ=sqrt(alpha*fourpi)/sqrt(sw2*(1.d0-sw2))
c      write(*,*) 'ez', ez
      DO I=0,12
         G=EZ*FMASS(I)*0.5d0/ZMASS
C
         GH(1,I) = DCMPLX( -G )
         GH(2,I) = DCMPLX( -G )
      ENDDO

c      write(*,*) 'gh', gh

      G1(1)  =  r_one
      G1(2)  =  r_one
c
c  QCD couplings
c
      g = 1.21712451d0
      als = g*g/fourpi
      gg(1)=-g
      gg(2)=-g

      write(*,*) 'Why am I here'
      stop
      write(*,*) 'G = ', g

c      write(*,*)  'sw2 ', 'gw ', 'gwwa ', 'gwwz '
c      write(*,*) sw2,gw,gwwa,gwwz
c      write(*,*)  'sw2 ', 'gal ', 'gau ', 'gad ', 'gwf '
c      write(*,*) sw2,gal,gau,gad,gwf
c      write(*,*)  'gzn ', 'gzl ', 'gzu ', 'gzd ', 'g1 '
c      write(*,*) gzn,gzl,gzu,gzd,g1
c      write(*,*)  'sw2 ', 'zmass1 ', 'hmass1 ', 'gwwh '
c      write(*,*) sw2,zmass,hmass,gwwh
c      write(*,*)  'gzzh ', 'ghhh ', 'gwwhh ', 'gzzhh ', 'ghhhh '
c      write(*,*) gzzh,ghhh,gwwhh,gzzhh,ghhhh
c      write(*,*)  'sw2 ', 'zmass1 ', 'fmass(i) ', 'gchf(1,i) '
c      do i=1,12
c         write(*,*) sw2,zmass,fmass(i),gh(1,i)
c      enddo

c      stop

      END

       
