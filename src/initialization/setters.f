************************************************************************
*
*     setters.f
*
************************************************************************
      subroutine setperturbativeorder(iptin)
      implicit none
      include "../commons/ipt.h"
      integer iptin
      if ( (iptin.lt.0) .or. (iptin.gt.1) ) then
         write(6,*) "in setperturbativeorder:"
         write(6,*) "invalid value"
         call exit(-10)         
      else
         ipt = iptin
      endif
      return
      end
************************************************************************            
      subroutine setperturbativeorderalpha(iptalphain)
      implicit none
      include "../commons/ipt.h"
      integer iptalphain
      if ( (iptalphain.lt.0) .or. (iptalphain.gt.1) ) then
         write(6,*) "in setperturbativeorderalpha:"
         write(6,*) "invalid value"
         call exit(-10)         
      else
         iptalpha = iptalphain
      endif
      return
      end
************************************************************************            
      subroutine setflavourscheme(nsin)
      implicit none
      include "../commons/ns.h"
      character*4 nsin
      if ( (ns.ne."VFNS") .and. (ns.ne."FFNS") ) then
         write(6,*) "in setflavourscheme:"
         write(6,*) "invalid value"
         call exit(-10)         
      else 
         ns = nsin
      endif
      return
      end
************************************************************************
      subroutine setflavourschemeint(nsinint)
      implicit none
      include "../commons/ns.h"
      integer nsinint
      if (nsinint.eq.0) then
         ns = "FFNS"
      elseif (nsinint.eq.1) then
         ns = "VFNS"
      else
         write(6,*) "in setflavourschemeint:"
         write(6,*) "invalid value"
         call exit(-10)
      endif 
      return
      end
************************************************************************
      subroutine setfactorisationscheme(fsin)
      implicit none
      include "../commons/facscheme.h"
      character*5 fsin
      if ( (fsin.ne."MSBAR") .and. (fsin.ne."DELTA") ) then
         write(6,*) "in setfactorisationscheme:"
         write(6,*) "invalid value"
         call exit(-10)         
      else
         facscheme = fsin
      endif
      return
      end
************************************************************************
      subroutine setfactorisationschemeint(fsinint)
      implicit none
      include "../commons/facscheme.h"
      integer fsinint
      if (fsinint.eq.0) then
         facscheme = "MSBAR"
      elseif (fsinint.eq.1) then
         facscheme = "DELTA"
      else
         write(6,*) "in setfactorisationschemeint:"
         write(6,*) "invalid value"
         call exit(-10)
      endif 
      return
      end
************************************************************************      
      subroutine setrenormalisationscheme(rsin)
      implicit none
      include "../commons/renscheme.h"
      character*5 rsin
      if ( (rsin.ne."MSBAR") .and. (rsin.ne."FIXED") 
     .     .and. (rsin.ne."ALPMZ") .and. (rsin.ne."ALGMU") ) then
         write(6,*) "in setrenormalisationscheme:"
         write(6,*) "invalid value"
         call exit(-10)         
      else      
         renscheme = rsin
      endif
      return
      end
************************************************************************
      subroutine setrenormalisationschemeint(rsinint)
      implicit none
      include "../commons/renscheme.h"
      integer rsinint
      if (rsinint.eq.0) then
         renscheme = "MSBAR"
      elseif (rsinint.eq.1) then
         renscheme = "FIXED"
      elseif (rsinint.eq.2) then
         renscheme = "ALPMZ"
      elseif (rsinint.eq.3) then
         renscheme = "ALGMU"
      else
         write(6,*) "in setrenormalisationscheme:"
         write(6,*) "invalid value"
         call exit(-10)
      endif 
      return
      end
************************************************************************      
      subroutine setactiveflav(nlmaxin,numaxin,ndmaxin)
      implicit none
      include "../commons/activeflavours.h"
      integer nlmaxin,numaxin,ndmaxin
      if ( (nlmaxin.lt.1) .or. (nlmaxin.gt.3) ) then
         write(6,*) "in setactiveflav:"
         write(6,*) "invalid value"
         call exit(-10)         
      else
         nlmax = nlmaxin
      endif
      if ( (numaxin.lt.0) .or. (numaxin.gt.2) ) then      
         write(6,*) "in setactiveflav:"
         write(6,*) "invalid value"
         call exit(-10)         
      else
         numax = numaxin         
      endif
      if ( (ndmaxin.lt.0) .or. (ndmaxin.gt.3) ) then      
         write(6,*) "in setactiveflav:"
         write(6,*) "invalid value"
         call exit(-10)         
      else
         ndmax = ndmaxin
      endif      
      return
      end
************************************************************************
      subroutine setactiveflavaem(nlmaxaemin,numaxaemin,ndmaxaemin)
      implicit none
      include "../commons/activeflavours.h"
      integer nlmaxaemin,numaxaemin,ndmaxaemin
      if ( (nlmaxaemin.lt.1) .or. (nlmaxaemin.gt.3) ) then
         write(6,*) "in setactiveflav:"
         write(6,*) "invalid value"
         call exit(-10)         
      else
         nlmaxaem = nlmaxaemin
      endif
      if ( (numaxaemin.lt.0) .or. (numaxaemin.gt.2) ) then      
         write(6,*) "in setactiveflav:"
         write(6,*) "invalid value"
         call exit(-10)         
      else
         numaxaem = numaxaemin         
      endif
      if ( (ndmaxaemin.lt.0) .or. (ndmaxaemin.gt.3) ) then      
         write(6,*) "in setactiveflav:"
         write(6,*) "invalid value"
         call exit(-10)         
      else
         ndmaxaem = ndmaxaemin
      endif      
      return
      end
************************************************************************      
      subroutine setwaem(waemin)
      implicit none
      include "../commons/activeflavours.h"
      integer waemin
      if ( (waemin.lt.0) .or. (waemin.gt.1) ) then
         write(6,*) "in setwaem:"
         write(6,*) "invalid value"
         call exit(-10)         
      else
         waem = waemin
      endif         
      return
      end
************************************************************************            
      subroutine setalpha(ain,qin)
      implicit none
      include "../commons/alpha.h"
      double precision ain, qin
      aref = ain
      q2ref = qin**2
      return
      end
************************************************************************
      subroutine setthresholds(me,mu,md,ms,mm,mc,mt,mb,mw,mz)
      implicit none
      include "../commons/massthrs.h"
      double precision me,mu,md,ms,mm,mc,mt,mb,mw,mz
      q2th(1) = me**2
      q2th(2) = mu**2
      q2th(3) = md**2
      q2th(4) = ms**2
      q2th(5) = mm**2
      q2th(6) = mc**2
      q2th(7) = mt**2
      q2th(8) = mb**2
      q2th(9) = mw**2
      q2th(10)= mz**2
      return
      end
************************************************************************
      subroutine setnphotot(varin)
      implicit none
      include "../commons/tecparam.h"
      integer varin
      nphotot = varin
      return
      end
************************************************************************     
      subroutine setnphomin(varin)
      implicit none
      include "../commons/tecparam.h"
      integer varin
      nphomin = varin
      return
      end
************************************************************************
      subroutine setnmatexp(varin)
      implicit none
      include "../commons/tecparam.h"
      integer varin
      nmatexp = varin
      return
      end
************************************************************************
      subroutine setnstpaem(varin)
      implicit none
      include "../commons/tecparam.h"
      integer varin
      nstpaem = varin
      return
      end
************************************************************************
      subroutine setminvmel(varin)
      implicit none
      include "../commons/tecparam.h"
      integer varin
      minvmel = varin
      return
      end
************************************************************************
      subroutine setrinvmel(varin)
      implicit none
      include "../commons/tecparam.h"
      integer varin
      rinvmel = varin
      return
      end
************************************************************************      
      subroutine setafixsolint(intin)
      implicit none
      include "../commons/tecparam.h"
      integer intin
      if (intin.eq.0) then
         afixsol = "PATHOR"
      elseif (intin.eq.1) then
         afixsol = "MAGNUS"
      else
         write(6,*) "in setafixsolint:"
         write(6,*) "invalid value"
         call exit(-10)
      endif 
      return
      end
************************************************************************
      subroutine setdeltagmu(deltagmuin)
      implicit none
      include "../commons/tecparam.h"
      double precision deltagmuin
      deltagmu = deltagmuin
      return
      end
      
