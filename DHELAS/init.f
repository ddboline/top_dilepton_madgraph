      SUBROUTINE INIT(var)

      implicit none

c====================================================================
c
c  Define common block containing all coupling constants and masses.
c
c====================================================================
c
      character*200        var

      double complex       gg(2)
      double precision     alpha,ee, sin2w, gfermi, alfas,g
      common /COUPL_BASIC/ alpha,ee, sin2w, gfermi, alfas,g,gg
c
      double precision     hmass, wmass, zmass, amass,
     &                     tmass, bmass, lmass, cmass
      common /COUPL_MASS/  hmass, wmass, zmass, amass,
     &                     tmass, bmass, lmass, cmass

      double precision     fmass(12), fwidth(12)      
      common /COUPL_MASS/  fmass, fwidth

      double precision     hwidth, wwidth, zwidth, 
     &                     twidth, lwidth, awidth
      common /COUPL_WIDTH/ hwidth, wwidth, zwidth, 
     &                     twidth, lwidth, awidth

      double complex       gh(2,12)
      common /COUP4/       gh

      double complex       gal(2), gad(2), gau(2), gwf(2),
     &                     gzn(2), gzl(2), gzd(2), gzu(2), g1(2)
      double precision     gw, gwwa, gwwz
      common /COUPL_GAUGE/ gal   , gad   , gau   , gwf   ,
     &                     gzn   , gzl   , gzd   , gzu   ,
     &                     gw, gwwa, gwwz, g1

      double complex       gwfc(2),  gwfs(2)
      common /coupl_ckm/   gwfc,     gwfs

      double complex       gwwh, gzzh, gwwhh, gzzhh, ghhh, ghhhh
      common /COUPL_SCAL/  gwwh, gzzh, gwwhh, gzzhh, ghhh, ghhhh

      double complex       ghtop(2), ghbot(2), ghtau(2), ghcha(2)
      common /COUPL_YUK/   ghtop   , ghbot   , ghtau   , ghcha

      double precision     xzmass, xwmass
      common /COUPL_XMASS/ xzmass, xwmass

      double complex       xzl(2) , xzb(2) , xzt(2) ,
     &     xwpq(2), xwmq(2), xwpl(2), xwml(2)
      common /COUPL_XFFV/  xzl    , xzb    , xzt    ,
     &                     xwpq   , xwmq   , xwpl   , xwml

      double complex       xzhz, xwhwp, xwhwm
      common /COUPL_XVSS/  xzhz, xwhwp, xwhwm

      double complex       xwzwp, xwzwm, xwawp, xwawm
      common /COUPL_XVVS/  xwzwp, xwzwm, xwawp, xwawm

      double complex       xwzhwp, xwzhwm, xwahwp, xwahwm
      common /COUPL_XVVSS/ xwzhwp, xwzhwm, xwahwp, xwahwm

      REAL*8 ez, ey
      REAL*8 r_zero, r_half, r_one, r_two, r_three, r_four, r_ote
      REAL*8 r_pi, r_ialph, fourpi
      PARAMETER( r_zero=0.0d0, r_half=0.5d0, r_one=1.0d0, r_two=2.0d0 )
      PARAMETER( r_three=3.0d0 )
      PARAMETER( r_four=4.0d0, r_ote=128.9d0 )
      PARAMETER( r_pi=3.14159265358979323846d0, r_ialph=137.0359895d0 )

      logical lhdecay,ltdecay,lwdecay,lzdecay
      common/to_set_par/lhdecay,ltdecay,lwdecay,lzdecay
c
      data amass/0d0/  ! set the photon mass to zero

      integer icollider
      integer                                        lpp(2)
      double precision    ebeam(2), xbk(2),q2fact(2)
      common/to_collider/ ebeam   , xbk   ,q2fact,   lpp

      real*8          scale
      logical               fixed_ren_scale,fixed_fac_scale
      common/to_scale/scale,fixed_ren_scale,fixed_fac_scale

c
c     Input values 
c
      double precision  hgf,hale, hss,wm,zm !EW input parameters
      double precision  mtMS,mbMS,mcMS      !MSbar masses
      double precision  Vud,Vus             !CKM matrix elements
      integer  idef                         !EW couplings option
      common/input_values/hgf,hale,hss,wm,zm,
     &                    mtMS,mbMS,mcMS,
     &                    Vud,Vus, 
     &                    idef


c     Calculated intermediate values that are printed out
c
      double precision  decw, w_w_nl, w_w_tau, w_w_ud, w_w_cs
      double precision  decz, w_z_nn, w_z_ll, w_z_tau
      double precision  w_z_uu, w_z_dd, w_z_cc, w_z_bb
      double precision  mt_h,mb_h,mc_h
      double precision  w_h_tt,w_h_bb,w_h_tau,w_h_cc
      double precision  w_h_ww,w_h_zz,w_h_WLWL,w_h_ZLZL
      common/calc_values/
     &                   decw, w_w_nl, w_w_tau, w_w_ud, w_w_cs,
     &                   decz, w_z_nn, w_z_ll, w_z_tau,
     &                   w_z_uu, w_z_dd, w_z_cc, w_z_bb,
     &                   mt_h, mb_h, mc_h,
     &                   w_h_tt,w_h_bb,w_h_tau,w_h_cc,
     &                   w_h_ww,w_h_zz,w_h_WLWL,w_h_ZLZL



      double precision GetG
      External GetG

      double precision  sw2, sw, cw, sc2, v, ee2
      double precision  g_temp
      integer i

c***********************************************************************
c     this files contains the common blocks for the 
c     pdf and the alpha_s settings
c
c     pdlabel= string identifying the pdf    
c     asmz   = alpha_s(Mz) is set based on the pdf chosen in setcuts.f
c     nloop  = order of the running of alpha_s based on the pdf chosen  
c***********************************************************************
      double precision asmz
      common/couple/asmz
      integer nloop
      common/running/nloop
      character*7 pdlabel
      common/to_pdf/pdlabel
      
      zmass = 91.188D0
c      scale=zmass   
c      asmz =0.130d0  
      pdlabel='cteq6_l'
      asmz =0.118d0   
      nloop=2

c************************************************************************
c     EW couplings: **INPUT** values. Output values can be different!   *
c************************************************************************

      hgf    = 1.16639d-5         !G_F  
      hale   = 1d0/128.9d0        !alpha_EM(m_Z)  
      hss    = 0.2312d0           !sin^2(theta_W) 
      wm     = 80.419d0           !W  mass 
      zm     = 91.188d0           !Z  mass 

c************************************************************************
c     Z and W widths                                                    *
c************************************************************************

      lzdecay =.true.   ! if true calculates the Z width @ LO 
      zwidth  = 2.495d0  ! PDG best value. Overidden if lzdecay=.true.
      lwdecay =.true.   ! if true calculates the W width @ LO 
      wwidth  = 2.12d0   ! PDG best value. Overidden if lwdecay=.true.

c************************************************************************
c     Higgs                                                             *
c************************************************************************

      hmass   = 120d0    !Higgs  mass
      lhdecay =.true.    !if true calls hdecay to calculate width and brs
      hwidth  = 3.7d-3   !Higgs width. Overidden if lhdecay=.true.

c************************************************************************     
c     Fermions                                                          *
c************************************************************************     

      tmass   = 174.3d0  !top pole mass
      ltdecay =.true.   !if true calculates the top width @ LO 
      twidth  = 1.6d0    !top width . Overridden if ltdecay=.true. 
      bmass   = 4.7d0    !bottom pole mass
      cmass   = 1.42d0   !charm mass
      lmass   = 1.777d0  !tau mass
      lwidth  = 2.36d-12 !tau width, PDG value

c************************************************************************     
c     MSbar masses for Yukawa couplings                                 *
c************************************************************************     

      mtMS = 164.5d0     !top    MSbar mass
      mbMS = 4.2d0       !botton MSbar mass
      mcMS = 1.42d0      !charm  MSbar mass

c************************************************************************     
c     CKM matrix elements                                               *
c     symmetric 3x3 matrix, Vud=Vcs, Vus=Vcd, Vcb=Vub=0                 *
c************************************************************************     
c     >>>>>>>>>>>>>>>****NOTE****<<<<<<<<<<<<<<<<<<<<<<<<<
c     these couplings matter ONLY when interaction_CKM.dat
c     has been used to generate all the diagrams with off-diagonal
c     couplings. The default of MadEvent is a diagonal 
c     CKM matrix AND any value for Vud and Vus is IGNORED!

      Vud=0.975d0
      Vus=dsqrt(1d0-Vud**2)

c************************************************************************     
c     Calculation scheme for EW couplings                               *
c************************************************************************     
c
c     idef=0   : Old MadEvent default (= AlpGen with iewopt=2) 
c                input values = sin^2(theta_W),alpha(m_Z),m_Z
c
c     idef=1   : New Madevent default, "G_mu scheme"
c                = LUSIFER and AlpGen (iewopt=3) defaults 
c                input values = G_F,m_Z,m_W
c
c     idef=2   : input  values = G_F,sin^2(theta_W),alpha_(m_Z) 
c                output values = m_W,m_Z.
c
c     idef=3   : User choice. Edit the subroutine coupsm in couplings.f
c                You have to know what you're doing.
c
      idef=1
      call ewcoup   !calculate the EW couplings

c************************************************************************     
c     Renormalization and factorization scales                          *
c************************************************************************     
c
c     Set the following flags to .false. if an event-by-event
c     scale choice is requested. In this case edit the 
c     user subroutines set_ren_scale and set_fac_scale below. 
c
      fixed_ren_scale=.true.  
      fixed_fac_scale=.true.
c     
c     Fixed RENORMALIZATION scale.   
c     For a dynamic choice edit the subroutine set_ren_scale below.      
c
      scale=wmass/2    ! ren scale. Overidden if fixed_ren_scale=.false.
c
c     Fixed FACTORIZATION scale 
c     For a dynamic choice edit the subroutine set_fac_scale below.       
c
c     If fixed_fac_scale=.false. the following values will be overridden.
c
      q2fact(1) = zmass**2      ! fact scale**2 for pdf1
      q2fact(2) = zmass**2      ! fact scale**2 for pdf2     


c************************************************************************     
c    Collider energy and type                                           *
c************************************************************************     
c
c     lpp  = -1 (antiproton), 0 (no pdf), 1 (proton)
c     ebeam= energy of each beam in GeV
c
c     predefined choices
c
c     icollider=0    User choice
c     icollider=1    e+ e- @   500 GeV  
c     icollider=2    P  Pb @  2000 GeV
c     icollider=3    P  P  @ 14000 GeV
c     icollider=4    Pb Pb @ 14000 GeV
c
      icollider = 2
c
      if (icollider .eq. 0) then  !User choice
         lpp(1) = 0
         lpp(2) = 0 
         ebeam(1) = 250d0
         ebeam(2) = ebeam(1)
      elseif (icollider .eq. 1) then  !e+e-         
         lpp(1) = 0
         lpp(2) = 0 
         ebeam(1) = 250d0
         ebeam(2) = ebeam(1)
      elseif (icollider .eq. 2) then  !Tevatron
         lpp(1) = 1       !Proton
         lpp(2) = -1      !Anti-proton
         ebeam(1) = 980d0
         ebeam(2) = ebeam(1)
      elseif (icollider .eq. 3) then  !LHC
         lpp(1) = 1       !Proton
         lpp(2) = 1       !Proton
         ebeam(1) = 7000d0
         ebeam(2) = ebeam(1)
      elseif (icollider .eq. 4) then  !LHC-BAR
         lpp(1) =-1       !Antiproton
         lpp(2) =-1       !Antiproton
         ebeam(1) = 7000d0
         ebeam(2) = ebeam(1)
      endif

c************************************************************************     
c    Sets up the couplings for HELAS                                    *
c************************************************************************     

      call coupsm(0)    ! to set HELAS couplings 

      write(*,*) '** Done setting couplings'
C      call writecoup(var)
      return

      END

      SUBROUTINE set_tmass(new_mtop)

      implicit none
      double precision new_mtop

      double precision     hmass, wmass, zmass, amass,
     &                     tmass, bmass, lmass, cmass
      common /COUPL_MASS/  hmass, wmass, zmass, amass,
     &                     tmass, bmass, lmass, cmass

      double precision     fmass(12), fwidth(12)      
      common /COUPL_MASS/  fmass, fwidth

      tmass = new_mtop

      return
      end
