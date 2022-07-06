************************************************************************
*
*     getters.f
*
************************************************************************
      subroutine getperturbativeorder(iptout)
      implicit none
      include "../commons/ipt.h"
      integer iptout
      iptout = ipt
      return
      end
************************************************************************      
      subroutine getperturbativeorderalpha(iptalphaout)
      implicit none
      include "../commons/ipt.h"
      integer iptalphaout
      iptalphaout = iptalpha
      return
      end
************************************************************************ 
      subroutine getflavourschemeint(nsintout)
      implicit none
      include "../commons/ns.h"
      integer nsintout
      if (ns.eq."FFNS") then
         nsintout = 0
      elseif (ns.eq."VFNS") then
         nsintout = 1
      endif
      return
      end
************************************************************************       
      subroutine getfactorisationschemeint(fsintout)
      implicit none
      include "../commons/facscheme.h"      
      integer fsintout
      if (facscheme.eq."MSBAR") then
         fsintout = 0
      elseif (facscheme.eq."DELTA") then
         fsintout = 1
      endif
      return
      end
************************************************************************ 
      subroutine getrenormalisationschemeint(rsintout)
      implicit none
      include "../commons/renscheme.h"      
      integer rsintout
      if (renscheme.eq."MSBAR") then
         rsintout = 0
      elseif (renscheme.eq."FIXED") then
         rsintout = 1
      elseif (renscheme.eq."ALPMZ") then
         rsintout = 2
      elseif (renscheme.eq."ALGMU") then
         rsintout = 3
      endif
      return
      end
************************************************************************      
      subroutine getactiveflav(nlmaxout,numaxout,ndmaxout)
      implicit none
      include "../commons/activeflavours.h"
      integer nlmaxout,numaxout,ndmaxout
      nlmaxout = nlmax
      numaxout = numax
      ndmaxout = ndmax
      return
      end
************************************************************************
      subroutine getactiveflavaem(nlmaxaemout,numaxaemout,ndmaxaemout)
      implicit none
      include "../commons/activeflavours.h"
      integer nlmaxaemout,numaxaemout,ndmaxaemout
      nlmaxaemout = nlmaxaem
      numaxaemout = numaxaem
      ndmaxaemout = ndmaxaem
      return
      end
************************************************************************
      subroutine getwaem(waemout)
      implicit none
      include "../commons/activeflavours.h"
      integer waemout
      waemout = waem
      return
      end
************************************************************************ 
      subroutine getalpharef(arefout)
      implicit none
      include "../commons/alpha.h"      
      double precision arefout
      arefout = aref
      return
      end
************************************************************************      
      subroutine getalphaqref(qrefout)
      implicit none
      include "../commons/alpha.h"      
      double precision qrefout
      qrefout = dsqrt(q2ref)
      return
      end
************************************************************************            
      subroutine getthresholds2(q2thrs)
      implicit none
      include "../commons/massthrs.h"      
      double precision q2thrs(11)
      q2thrs=q2th
      return
      end
************************************************************************
      subroutine getnphotot(varout)
      implicit none
      include "../commons/tecparam.h"
      integer varout
      varout = nphotot
      return
      end
************************************************************************
      subroutine getnphomin(varout)
      implicit none
      include "../commons/tecparam.h"
      integer varout
      varout = nphomin
      return
      end      
************************************************************************
      subroutine getnmatexp(varout)
      implicit none
      include "../commons/tecparam.h"
      integer varout
      varout = nmatexp
      return
      end
************************************************************************
      subroutine getnstpaem(varout)
      implicit none
      include "../commons/tecparam.h"
      integer varout
      varout = nstpaem
      return
      end
************************************************************************
      subroutine getminvmel(varout)
      implicit none
      include "../commons/tecparam.h"
      integer varout
      varout = minvmel
      return
      end
************************************************************************
      subroutine getrinvmel(varout)
      implicit none
      include "../commons/tecparam.h"
      integer varout
      varout = rinvmel
      return
      end
************************************************************************            
      subroutine getdeltagmu(varout)
      implicit none
      include "../commons/tecparam.h"
      double precision varout
      varout = deltagmu
      return
      end
************************************************************************            
      subroutine getafixsolint(varout)
      implicit none
      include "../commons/tecparam.h"
      integer varout
      if (afixsol.eq."PATHOR") then      
         varout = 0
      elseif (afixsol.eq."MAGNUS") then
         varout = 1
      endif         
      return
      end      
************************************************************************      
      subroutine getb0(b0)
      implicit none
      include "../commons/beta.h"
      include "../commons/consts.h"      
      double precision b0(1:11)
      b0(:) = -beta0(:)/4d0/pi
      return
      end      
************************************************************************
      subroutine getb1(b1)
      implicit none
      include "../commons/beta.h"
      include "../commons/consts.h"      
      double precision b1(1:11)
      b1(:) = -beta1(:)/4d0/4d0/pi/pi
      return
      end      
************************************************************************
      subroutine getc2(c2)
      implicit none
      include "../commons/nfsum.h"
      double precision c2(1:11)
      c2(:) = nfsum2(:)
      return
      end
************************************************************************
      subroutine getc4(c4)
      implicit none
      include "../commons/nfsum.h"
      double precision c4(1:11)
      c4(:) = nfsum4(:)
      return
      end
************************************************************************
      subroutine getmz2(mz2out)      
      implicit none
      include "../commons/massthrs.h"      
      double precision mz2out
      mz2out=q2th(10)
      return
      end
************************************************************************
      subroutine getmw2(mw2out)      
      implicit none
      include "../commons/massthrs.h"      
      double precision mw2out
      mw2out=q2th(9)
      return
      end
************************************************************************
      subroutine geta0(a0)
      implicit none
      include "../commons/alpha.h"
      double precision a0
      a0 = ath(1)
      return
      end

      
      
