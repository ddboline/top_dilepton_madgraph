      subroutine SetSilent(isilent)
      implicit none

      integer lhasilent,isilent
      common/lhasilent/lhasilent

      lhasilent=isilent
      return

      end
