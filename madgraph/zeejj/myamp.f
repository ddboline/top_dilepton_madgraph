      double precision function testamp(p)
c*****************************************************************************
c     Approximates matrix element by propagators
c*****************************************************************************
      implicit none
c
c     Constants
c     
      include 'genps.inc'
      double precision   zero
      parameter (zero = 0d0)
c
c     Arguments
c
      double precision p(0:3,nexternal)
c      integer iconfig
c
c     Local
c
      double precision xp(0:3,-nexternal:nexternal)
      double precision mpole(-nexternal:0),shat,tsgn
      integer i,j,iconfig

      double precision pmass(-nexternal:0,lmaxconfigs)
      double precision pwidth(-nexternal:0,lmaxconfigs)
      integer pow(-nexternal:0,lmaxconfigs)
      logical first_time
c
c     Global
c
      integer iforest(2,-max_branch:-1,lmaxconfigs)
      common/to_forest/ iforest
      integer            mapconfig(0:lmaxconfigs), this_config
      common/to_mconfigs/mapconfig, this_config
      
      include 'coupl.inc'
c
c     External
c
      double precision dot

      save pmass,pwidth,pow
      data first_time /.true./
c-----
c  Begin Code
c-----      
      iconfig = this_config
      if (first_time) then
c         include 'props.inc'
         first_time=.false.
      endif

      do i=1,nexternal
         mpole(-i)=0d0
         do j=0,3
            xp(j,i)=p(j,i)
         enddo
      enddo
c      mpole(-3) = 174**2
c      shat = dot(p(0,1),p(0,2))/(1800)**2
      shat = dot(p(0,1),p(0,2))/(10)**2
c      shat = 1d0
      testamp = 1d0
      tsgn    = +1d0
      do i=-1,-(nexternal-3),-1              !Find all the propagotors
         if (iforest(1,i,iconfig) .eq. 1) tsgn=-1d0
         do j=0,3
            xp(j,i) = xp(j,iforest(1,i,iconfig))
     $           +tsgn*xp(j,iforest(2,i,iconfig))
         enddo
         if (pwidth(i,iconfig) .ne. 0d0 .and. .false.) then
            testamp=testamp/((dot(xp(0,i),xp(0,i))
     $                        -pmass(i,iconfig)**2)**2
     $         -(pmass(i,iconfig)*pwidth(i,iconfig))**2)
         else
            testamp = testamp/((dot(xp(0,i),xp(0,i)) -
     $                          pmass(i,iconfig)**2)
     $                          **(pow(i,iconfig)))
         endif
        testamp=testamp*shat**(pow(i,iconfig))
c        write(*,*) i,iconfig,pow(i,iconfig),pmass(i,iconfig)
      enddo
c      testamp = 1d0/dot(xp(0,-1),xp(0,-1))
      testamp=abs(testamp)
c      testamp = 1d0
      end

      logical function cut_bw(p)
c*****************************************************************************
c     Approximates matrix element by propagators
c*****************************************************************************
      implicit none
c
c     Constants
c     
      include 'genps.inc'
      double precision   zero
      parameter (zero = 0d0)
c
c     Arguments
c
      double precision p(0:3,nexternal)
c
c     Local
c
      double precision xp(0:3,-nexternal:nexternal)
      double precision mpole(-nexternal:0),shat,tsgn
      integer i,j,iconfig

      double precision pmass(-nexternal:0,lmaxconfigs)
      double precision pwidth(-nexternal:0,lmaxconfigs)
      integer pow(-nexternal:0,lmaxconfigs)
      logical first_time, onshell
      double precision xmass
      integer nbw
c
c     Global
c
      integer iforest(2,-max_branch:-1,lmaxconfigs)
      common/to_forest/ iforest
      integer            mapconfig(0:lmaxconfigs), this_config
      common/to_mconfigs/mapconfig, this_config

      integer        lbw(0:nexternal)  !Use of B.W.
      common /to_BW/ lbw
      
      include 'coupl.inc'
c
c     External
c
      double precision dot

      save pmass,pwidth,pow
      data first_time /.true./
c-----
c  Begin Code
c-----      
      cut_bw = .false.    !Default is we passed the cut
      iconfig = this_config
      if (first_time) then
         include 'props.inc'
         nbw = 0
         tsgn = 1d0
         do i=-1,-(nexternal-3),-1
            if (iforest(1,i,iconfig) .eq. 1) tsgn=-1d0
            if (pwidth(i,iconfig) .gt. 0 .and. tsgn .eq. 1d0) then
               nbw=nbw+1
               if (lbw(nbw) .eq. 1) then
                  write(*,*) 'Requiring BW ',i,nbw
               elseif(lbw(nbw) .eq. 2) then
                  write(*,*) 'Excluding BW ',i,nbw
               else
                  write(*,*) 'No cut BW ',i,nbw
               endif
            endif
         enddo
         first_time=.false.
      endif

      do i=1,nexternal
         mpole(-i)=0d0
         do j=0,3
            xp(j,i)=p(j,i)
         enddo
      enddo
      nbw = 0
      tsgn    = +1d0
      do i=-1,-(nexternal-3),-1              !Loop over propagators
         if (iforest(1,i,iconfig) .eq. 1) tsgn=-1d0
         do j=0,3
            xp(j,i) = xp(j,iforest(1,i,iconfig))
     $           +tsgn*xp(j,iforest(2,i,iconfig))
         enddo
         if (tsgn .gt. 0d0 .and. pwidth(i,iconfig) .gt. 0d0 ) then !This is B.W.
            nbw = nbw+1
c            write(*,*) 'Checking BW',nbw
            xmass = sqrt(dot(xp(0,i),xp(0,i)))
c            write(*,*) 'xmass',xmass,pmass(i,iconfig)
            onshell = (abs(xmass - pmass(i,iconfig)) .lt.
     $           5.*pwidth(i,iconfig))
            if (onshell .and. (lbw(nbw).eq. 2) ) cut_bw=.true.
            if (.not. onshell .and. (lbw(nbw).eq. 1)) cut_bw=.true.
c            write(*,*) nbw,xmass,onshell,lbw(nbw),cut_bw
         endif
      enddo
      end


      subroutine set_peaks
c*****************************************************************************
c     Attempts to determine peaks for this configuration
c*****************************************************************************
      implicit none
c
c     Constants
c     
      include 'genps.inc'
      double precision   zero
      parameter (zero = 0d0)
c
c     Arguments
c
c
c     Local
c
      double precision  xm(-nexternal:nexternal)
      double precision  xe(-nexternal:nexternal)
      double precision tsgn, xo, a
      double precision x1,x2,xk(nexternal)
      double precision dr,mtot,etot, stot
      integer i, iconfig, l1, l2, j, nt, nbw

      double precision pmass(-nexternal:0,lmaxconfigs)
      double precision pwidth(-nexternal:0,lmaxconfigs)
      integer pow(-nexternal:0,lmaxconfigs)

c
c     Global
c
      integer iforest(2,-max_branch:-1,lmaxconfigs)
      common/to_forest/ iforest

      integer            mapconfig(0:lmaxconfigs), this_config
      common/to_mconfigs/mapconfig, this_config

      real*8         emass(nexternal)
      common/to_mass/emass

      integer                                        lpp(2)
      double precision    ebeam(2), xbk(2),q2fact(2)
      common/to_collider/ ebeam   , xbk   ,q2fact,   lpp

      double precision etmin(3:nexternal),etamax(3:nexternal)
      double precision                    r2min(3:nexternal,3:nexternal)
      double precision s_min(nexternal,nexternal)
      common/to_cuts/  etmin     ,etamax     , r2min, s_min

      double precision      spole(maxinvar),swidth(maxinvar),bwjac
      common/to_brietwigner/spole          ,swidth          ,bwjac

      integer        lbw(0:nexternal)  !Use of B.W.
      common /to_BW/ lbw
      
      include 'coupl.inc'

c
c     External
c

c-----
c  Begin Code
c-----      
      include 'props.inc'
      stot = 4d0*ebeam(1)*ebeam(2)
c      etmin = 10
      nt = 0
      iconfig = this_config
      mtot = 0d0
      etot = 0d0   !Total energy needed
      do i=3,nexternal  !assumes 2 incoming
         xm(i)=emass(i)
         xe(i)=max(emass(i),etmin(i))
         xk(i)= 0d0
         etot = etot+xe(i)
         mtot=mtot+xm(i)         
      enddo
      tsgn    = +1d0
      nbw = 0
      do i=-1,-(nexternal-3),-1              !Find all the propagotors
         if (iforest(1,i,iconfig) .eq. 1) tsgn=-1d0
         if (tsgn .eq. 1d0) then                         !s channel
            xm(i) = xm(iforest(1,i,iconfig))+xm(iforest(2,i,iconfig))
            xe(i) = xe(iforest(1,i,iconfig))+xe(iforest(2,i,iconfig))
            mtot = mtot - xm(i)
            etot = etot - xe(i)
            if (iforest(1,i,iconfig) .gt. 0
     &           .and. iforest(2,i,iconfig) .gt. 0) then
               xm(i)=max(xm(i),
     &           sqrt(s_min(iforest(1,i,iconfig),iforest(2,i,iconfig))))
               xe(i)=max(xe(i),xm(i))
            endif
c            write(*,*) 'iconfig,i',iconfig,i
c            write(*,*) pwidth(i,iconfig),pmass(i,iconfig)
            if (pwidth(i,iconfig) .gt. 0) nbw=nbw+1
            if (pwidth(i,iconfig) .gt. 0 .and. (lbw(nbw) .eq. 1 .or.
     $           lbw(0) .eq. 0)) then         !B.W.

c               nbw = nbw +1

               if (i .eq. -(nexternal-3)) then  !This is s-hat
                  j = 3*(nexternal-2)-4+1    !set i to ndim+1
                  write(*,*) 'Setting PDF BW',j,nbw,pmass(i,iconfig)
                  spole(j)=pmass(i,iconfig)*pmass(i,iconfig)/stot
                  swidth(j) = pwidth(i,iconfig)*pmass(i,iconfig)/stot
                  xm(i) = pmass(i,iconfig)
               else
                  write(*,*) 'Setting BW',i,nbw,pmass(i,iconfig)
                  spole(-i)=pmass(i,iconfig)*pmass(i,iconfig)/stot
                  swidth(-i) = pwidth(i,iconfig)*pmass(i,iconfig)/stot
                  xm(i) = pmass(i,iconfig)
               endif
            else                                  !1/x^pow
               if (xm(i) - pmass(i,iconfig) .le. 0d0) then !Can hit pole
c                  write(*,*) 'Setting new min',i,xm(i),pmass(i,iconfig)
                  l1 = iforest(1,i,iconfig)                  !need dr cut
                  l2 = iforest(2,i,iconfig)
                  if (l2 .lt. l1) then
                     j = l1
                     l1 = l2
                     l2 = j
                  endif
                  dr = 0
                  if (l1 .gt. 0) dr = r2min(l2,l1)**2 !dr only for external
c                  write(*,*) 'using r2min',l2,l1,sqrt(dr)
                  dr = dr*.8d0                        !0.8 to hit peak hard
                  xo = 0.5d0*pmass(i,iconfig)**2      !0.5 to hit peak hard
                  if (dr .gt. 0d0) xo = max(xo,etmin(l2)*etmin(l1)*dr)
c                  write(*,*) 'Got dr',i,l1,l2,dr
                  xo = xo+pmass(i,iconfig)**2
c                  write(*,*) 'Setting xm',i,xm(i),sqrt(xo)
                  xm(i) = sqrt(xo)    !Reset xm to realistic minimum
                  xo = xo/stot
c                  xo = sqrt(pmass(i,iconfig)**2+ pmass(i,iconfig)**2)
c                  xo = pmass(i,iconfig)+etmin
               else
c                  write(*,*) 'Using xm',i,xm(i)
                  xo = xm(i)**2/stot
               endif
               a=pmass(i,iconfig)**2/stot
c               call setgrid(-i,xo,a,pow(i,iconfig))

c               write(*,*) 'Enter minimum for ',-i, xo
c               read(*,*) xo

               if (pwidth(i,iconfig) .eq. 0) call setgrid(-i,xo,a,1)
               if (pwidth(i,iconfig) .gt. 0) then
                  write(*,*) 'Using flat grid for BW',i,nbw,
     $                 pmass(i,iconfig)
               endif
            endif
            xe(i) = max(xm(i),xe(i))               
            etot = etot+xe(i)
            mtot=mtot+max(xm(i),xm(i))
c            write(*,*) 'New mtot',i,mtot,xm(i)
         else                                        !t channel
c
c     Check closest to p1
c
            nt = nt+1
            l2 = iforest(2,i,iconfig) !need dr cut
            x1 = 0            
            if (l2 .gt. 0) x1 = etmin(l2)
            x1 = max(x1, xm(l2)/1d0)
            if (nt .gt. 1) x1 = max(x1,xk(nt-1))
            xk(nt)=x1
c            write(*,*) 'Using 1',l2,x1

c
c     Check closest to p2
c
            j = i-1
            l2 = iforest(2,j,iconfig)
            x2 = 0
            if (l2 .gt. 0) x2 = etmin(l2)
            x2 = max(x2, xm(l2)/1d0)
c            if (nt .gt. 1) x2 = max(x2,xk(nt-1))
            
c            write(*,*) 'Using 2',l2,x2

            xo = min(x1,x2)

            xo = xo*xo/stot
            a=-pmass(i,iconfig)**2/stot
c            call setgrid(-i,xo,a,pow(i,iconfig))

c               write(*,*) 'Enter minimum for ',-i, xo
c               read(*,*) xo

             call setgrid(-i,xo,a,1)
         endif
      enddo
      if (abs(lpp(1)) .eq. 1 .or. abs(lpp(2)) .eq. 1) then
         write(*,*) 'etot',etot,nexternal
         xo = etot**2/stot
         i = 3*(nexternal-2) - 4 + 1
         if (swidth(i) .eq. 0d0) then
            call setgrid(i,xo,0d0,1)
         endif
      endif

      i=-8
c      write(*,*) 'Enter minimum for ',-i, xo
c      read(*,*) xo      
c      if (xo .gt. 0)      call setgrid(-i,xo,a,1)

      i=-10
c      write(*,*) 'Enter minimum for ',-i, xo
c      read(*,*) xo      
c      if (xo .gt. 0) call setgrid(-i,xo,a,1)

      end



      subroutine find_matches(iconfig,isym)
c*****************************************************************************
c     Finds all of the matches to see what gains may be possible
c*****************************************************************************
      implicit none
c
c     Constants
c     
      include 'genps.inc'
      double precision   zero
      parameter (zero = 0d0)
c
c     Arguments
c
      integer iconfig,isym(0:10)
c
c     Local
c
      integer i,j, nc, imatch
      double precision tprop(3,nexternal),xprop(3,nexternal)
c
c     Global
c
      integer iforest(2,-max_branch:-1,lmaxconfigs)
      common/to_forest/ iforest

      integer            mapconfig(0:lmaxconfigs), this_config
      common/to_mconfigs/mapconfig, this_config

c      include 'coupl.inc'
c      include 'configs.inc'
c
c     External
c

c-----
c  Begin Code
c-----      
c
c     First get peaks for configuration interested in
c
      if (.false.) then
         isym(0)=1
         isym(1)=iconfig
      else
         call get_peaks(iconfig,tprop)
         isym(0)=0
         do i=1,mapconfig(0)
            call get_peaks(i,xprop)
c        call sort_prop(xprop)
            call match_peak(tprop,xprop,imatch)
            if (imatch .eq. 1) then
               write(*,*) 'Found Match',iconfig,i
               isym(0)=isym(0)+1
               isym(isym(0)) = i
            endif
         enddo
      endif
      write(*,'(20i4)') (isym(i),i=0,isym(0))
      end

      subroutine find_all_matches()
c*****************************************************************************
c     Finds all of the matches to see what gains may be possible
c*****************************************************************************
      implicit none
c
c     Constants
c     
      include 'genps.inc'
      double precision   zero
      parameter (zero = 0d0)
c
c     Arguments
c
      double precision tprop(3,nexternal),xprop(3,nexternal)
      integer iconfig
c
c     Local
c
      integer i,j, nc, imatch
      logical gm(lmaxconfigs)
c
c     Global
c
      integer iforest(2,-max_branch:-1,lmaxconfigs)
      common/to_forest/ iforest

      integer            mapconfig(0:lmaxconfigs), this_config
      common/to_mconfigs/mapconfig, this_config

c      include 'coupl.inc'
c      include 'configs.inc'
c
c     External
c

c-----
c  Begin Code
c-----      
      nc = 0
      do i=1,mapconfig(0)
         gm(i)=.false.
      enddo
      do j=1,mapconfig(0)
         if (.not. gm(j)) then
            nc=nc+1
c            write(*,*) 'Need config ',j
            call get_peaks(j,tprop)
c            call sort_prop(tprop)
            write(*,'(i4,4e12.4)') j,(tprop(1,i), i=1,4)
            do i=j+1,mapconfig(0)
               call get_peaks(i,xprop)
c               call sort_prop(xprop)
               call match_peak(tprop,xprop,imatch)
               if (imatch .eq. 1) then
                  write(*,*) 'Found Match',j,i
                  gm(i)=.true.
               endif
            enddo
         endif
      enddo
      write(*,*) 'Found matches',mapconfig(0),nc
      stop
      end

      subroutine sort_prop(xprop)
c*****************************************************************************
c     Sort props in order from min to max based on 1st component only
c*****************************************************************************
      implicit none
c
c     Constants
c     
      include 'genps.inc'
      double precision   zero
      parameter (zero = 0d0)
c
c     Arguments
c
      double precision xprop(3,nexternal)
c
c     Local
c
      integer i,j,imin
      double precision temp(3,nexternal),xmin
      logical used(nexternal)
c
c     Global
c
c
c     External
c
c-----
c  Begin Code
c-----      
      do i=1,nexternal-3
         used(i)=.false.
      enddo
      do j=1,nexternal-3
         do i=1,nexternal-3
            xmin = 2d0
            if (.not. used(i) .and. xprop(1,i) .lt. xmin) then
               xmin = xprop(1,i)
               imin = i
            endif
         enddo
         do i=1,3
            temp(i,j)=xprop(i,imin)
         enddo
         used(i)=.true.
      enddo
      do i=1,nexternal-3
         do j=1,3
            xprop(j,i)=temp(j,i)
         enddo
      enddo
      end


      subroutine match_peak(tprop,xprop,imatch)
c*****************************************************************************
c     Determines if two sets of peaks are equivalent
c*****************************************************************************
      implicit none
c
c     Constants
c     
      include 'genps.inc'
      double precision   zero
      parameter (zero = 0d0)
c
c     Arguments
c
      double precision xprop(3,nexternal),tprop(3,nexternal)
      integer imatch
c
c     Local
c
      integer i,j
c
c     Global
c
c
c     External
c
c-----
c  Begin Code
c-----      
      imatch = 1                     !By default assume match
      do i=1,nexternal-3
         do j=1,3
            if (tprop(j,i) .ne. xprop(j,i)) imatch=0
         enddo
      enddo
      end

      subroutine get_peaks(iconfig,xt)
c*****************************************************************************
c     Attempts to determine peaks for this configuration
c*****************************************************************************
      implicit none
c
c     Constants
c     
      include 'genps.inc'
      double precision   zero
      parameter (zero = 0d0)
c
c     Arguments
c
      double precision xt(3,nexternal)
      integer iconfig

c
c     Local
c
      double precision  xm(-nexternal:nexternal)
      double precision tsgn, xo, a
      double precision x1,x2
      double precision dr,mtot, stot
      integer i, l1, l2, j

      double precision pmass(-nexternal:0,lmaxconfigs)
      double precision pwidth(-nexternal:0,lmaxconfigs)
      integer pow(-nexternal:0,lmaxconfigs)

      integer imatch(0:lmaxconfigs)

c
c     Global
c
      integer iforest(2,-max_branch:-1,lmaxconfigs)
      common/to_forest/ iforest

      integer            mapconfig(0:lmaxconfigs), this_config
      common/to_mconfigs/mapconfig, this_config

      real*8         emass(nexternal)
      common/to_mass/emass

      integer                                        lpp(2)
      double precision    ebeam(2), xbk(2),q2fact(2)
      common/to_collider/ ebeam   , xbk   ,q2fact,   lpp

      double precision etmin(3:nexternal),etamax(3:nexternal)
      double precision                    r2min(3:nexternal,3:nexternal)
      double precision s_min(nexternal,nexternal)
      common/to_cuts/  etmin     ,etamax     , r2min, s_min

      double precision      spole(maxinvar),swidth(maxinvar),bwjac
      common/to_brietwigner/spole          ,swidth          ,bwjac
      
      include 'coupl.inc'
c      include 'props.inc'

c
c     External
c

c-----
c  Begin Code
c-----      
      stot = 4d0*ebeam(1)*ebeam(2)
c      iconfig = this_config
      mtot = 0d0
      do i=1,nexternal
         xm(i)=emass(i)
         mtot=mtot+xm(i)
      enddo
      tsgn    = +1d0
      do i=-1,-(nexternal-3),-1              !Find all the propagotors
         if (iforest(1,i,iconfig) .eq. 1) tsgn=-1d0
         if (tsgn .eq. 1d0) then                         !s channel
            xm(i) = xm(iforest(1,i,iconfig))+xm(iforest(2,i,iconfig))
            mtot = mtot - xm(i)
            if (pwidth(i,iconfig) .gt. 0) then    !B.W.
               write(*,*) 'Setting BW',i,pmass(i,iconfig)
               spole(-i)=pmass(i,iconfig)*pmass(i,iconfig)/stot
               swidth(-i) = pwidth(i,iconfig)*pmass(i,iconfig)/stot
               xm(i) = pmass(i,iconfig)
               xt(1,-i) = spole(-i)
               xt(2,-i) = swidth(-i)
               xt(3,-i) = 2
            else                                  !1/x^pow
               if (xm(i) - pmass(i,iconfig) .le. 0d0) then !Can hit pole
c                  write(*,*) 'Setting new min',i,xm(i),pmass(i,iconfig)
                  l1 = iforest(1,i,iconfig)                  !need dr cut
                  l2 = iforest(2,i,iconfig)
                  if (l2 .lt. l1) then
                     j = l1
                     l1 = l2
                     l2 = j
                  endif
                  dr = 0
                  if (l1 .gt. 0) dr = r2min(l2,l1)**2 !dr only for external
                  dr = dr*.8d0                        !0.8 to hit peak hard
                  xo = 0.5d0*pmass(i,iconfig)**2      !0.5 to hit peak hard
                  if (dr .gt. 0d0) xo = max(xo,etmin(l2)*etmin(l1)*dr)
c                  write(*,*) 'Got dr',i,l1,l2,dr
                  xo = xo+pmass(i,iconfig)**2
c                  write(*,*) 'Setting xm',i,xm(i),sqrt(xo)
                  xm(i) = sqrt(xo)    !Reset xm to realistic minimum
                  xo = xo/stot
c                  xo = sqrt(pmass(i,iconfig)**2+ pmass(i,iconfig)**2)
c                  xo = pmass(i,iconfig)+etmin
               else
c                  write(*,*) 'Using xm',i,xm(i)
                  xo = xm(i)**2/stot
               endif
               a=pmass(i,iconfig)**2/stot
c               call setgrid(-i,xo,a,pow(i,iconfig))
c               call setgrid(-i,xo,a,1)
               xt(1,-i) = xo
               xt(2,-i) = a
               xt(3,-i) = 1
            endif
            mtot=mtot+xm(i)
c            write(*,*) 'New mtot',i,mtot,xm(i)
         else                                        !t channel
c
c     Check closest to p1
c
            l2 = iforest(2,i,iconfig) !need dr cut
            x1 = 0
            if (l2 .gt. 0) x1 = etmin(l2)
            x1 = max(x1, xm(l2)/2d0)
c            write(*,*) 'Using 1',l2,x1

c
c     Check closest to p2
c
            j = i-1
            l2 = iforest(2,j,iconfig)
            x2 = 0
            if (l2 .gt. 0) x2 = etmin(l2)
            x2 = max(x2, xm(l2)/2d0)
            
c            write(*,*) 'Using 2',l2,x2

            xo = min(x1,x2)

            xo = xo*xo/stot
            a=-pmass(i,iconfig)**2/stot
c            call setgrid(-i,xo,a,pow(i,iconfig))
c            call setgrid(-i,xo,a,1)
               xt(1,-i) = xo
               xt(2,-i) = a
               xt(3,-i) = 1
         endif
      enddo
      if (abs(lpp(1)) .eq. 1 .or. abs(lpp(2)) .eq. 1) then
         xo = mtot**2/stot
         i = 3*(nexternal-2) - 4 + 1
c         call setgrid(i,xo,0d0,1)
      endif
      end


      subroutine writeamp(p)
c
c     Constants
c     
      include 'genps.inc'
c
c     Arguments
c
      double precision p(0:3,nexternal)
c
c     Local
c
      double precision xp(0:3,-nexternal:nexternal)
      integer i,j,iconfig
c
c     Global
c
      integer iforest(2,-max_branch:-1,lmaxconfigs)
      common/to_forest/ iforest

      integer            mapconfig(0:lmaxconfigs), this_config
      common/to_mconfigs/mapconfig, this_config

c
c     External
c
      double precision dot
c-----
c  Begin Code
c-----
      iconfig = this_config
      do i=1,nexternal
         do j=0,3
            xp(j,i)=p(j,i)
         enddo
      enddo
      shat = dot(p(0,1),p(0,2))
      testamp = 1d0
      tsgn    = +1d0
      do i=-1,-(nexternal-3),-1              !Find all the propagotors
         if (iforest(1,i,iconfig) .eq. 1) tsgn=-1d0
         do j=0,3
            xp(j,i) = xp(j,iforest(1,i,iconfig))
     $           +tsgn*xp(j,iforest(2,i,iconfig))
         enddo
c         testamp=testamp*shat/(dot(xp(0,i),xp(0,i)))**2
c         testamp=testamp**2
      enddo
      write(*,'(A,4e15.5)') 'V',(dot(xp(0,-i),xp(0,-i)),i=1,nexternal-3)
c      testamp=abs(testamp)
      end


      subroutine histamp(p,dwgt)
c
c     Constants
c     
      include 'genps.inc'
c
c     Arguments
c
      double precision p(0:3,nexternal),dwgt
c
c     Local
c
      double precision xp(0:3,-nexternal:nexternal)
      integer i,j,iconfig
      real wgt
c
c     Global
c
      integer iforest(2,-max_branch:-1,lmaxconfigs)
      common/to_forest/ iforest

      integer            mapconfig(0:lmaxconfigs), this_config
      common/to_mconfigs/mapconfig, this_config

c
c     External
c
      double precision dot
c-----
c  Begin Code
c-----
      wgt = dwgt
      iconfig = this_config
      do i=1,nexternal
         do j=0,3
            xp(j,i)=p(j,i)
         enddo
      enddo
      shat = dot(p(0,1),p(0,2))
      testamp = 1d0
      tsgn    = +1d0
      do i=-1,-(nexternal-3),-1              !Find all the propagotors
         if (iforest(1,i,iconfig) .eq. 1) tsgn=-1d0
         do j=0,3
            xp(j,i) = xp(j,iforest(1,i,iconfig))
     $           +tsgn*xp(j,iforest(2,i,iconfig))
         enddo
      enddo
      do i=1,nexternal-3
c         write(*,*) sqrt(abs(dot(xp(0,-i),xp(0,-i)))),wgt
         call hfill(10+i,real(sqrt(abs(dot(xp(0,-i),xp(0,-i))))),0.,wgt)
      enddo
      end


      subroutine check_limits(p,xlimit, iconfig)
c*************************************************************************
c     Checks limits on all of the functions being integrated
c*************************************************************************
      implicit none
c
c     Constants
c
      include 'genps.inc'
c
c     Arguments
c
      double precision p(0:3,nexternal), xlimit(2,nexternal)
      integer iconfig
c
c     Local
c
      double precision xp(0:3,-nexternal:nexternal), tsgn
      double precision sm
      integer ik(4)
      integer i,j, k1,k2
c
c     Global
c
      integer iforest(2,-max_branch:-1,lmaxconfigs)
      common/to_forest/ iforest
      integer            mapconfig(0:lmaxconfigs), this_config
      common/to_mconfigs/mapconfig, this_config
c
c     External
c
      double precision dot
      external dot

c      data ik/1,2,3,0/
c-----
c  Begin Code
c-----
c
c     Transform from rambo(1:4) format to helas (0:3)
c     
      do i=1,nexternal
         do j=0,3
            xp(j,i)=p(j,i)
         enddo
      enddo
c
c     Now build propagators
c      
      tsgn=+1d0
      do i=-1,-(nexternal-3),-1
         if (iforest(1,i,iconfig) .eq. 1) tsgn=-1d0
         k1=iforest(1,i,iconfig)
         k2=iforest(2,i,iconfig)
         do j=0,3
            xp(j,i)=xp(j,k1)+tsgn*xp(j,k2)
         enddo
         sm = tsgn*dot(xp(0,i),xp(0,i))
         if (sm .lt. xlimit(1,-i)) then
            xlimit(1,-i)=sm
c            write(*,*) 'New limit',-i,sm
         endif
         if (sm .gt. xlimit(2,-i)) then
            xlimit(2,-i)=sm
c            write(*,*) 'New limit',-i,sm
         endif
      enddo
      end
