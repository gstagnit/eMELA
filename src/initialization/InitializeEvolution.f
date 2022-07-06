************************************************************************
*
*     InitializeEvolution.f:
*
*     Initialization for the production of the evolution tables.
*
************************************************************************
      subroutine InitializeEvolution
*
      implicit none
*
      include "../commons/ipt.h"
      include "../commons/ns.h"
      include "../commons/facscheme.h"
      include "../commons/renscheme.h"
      include "../commons/activeflavours.h"
      include "../commons/alpha.h"      
      include "../commons/massthrs.h"
      include "../commons/tecparam.h"
**
*     Internal variables
*
      integer nf
*
*     Initialize the coefficient of the beta functions in QED
*
      call initBeta
*
*     Alpha at Q20 and at the heavy lepton mass thresholds
*
      call initCouplings
*
*     Evolution mode
*
      write(6,*) "  "
      write(6,*) "======= Evolution parameters ======="
      write(6,*) "  "
*
*     Perturbative order
*
      if(ipt.eq.0)then
         write(6,*) "Perturbative order: LO"
      elseif(ipt.eq.1)then
         write(6,*) "Perturbative order: NLO"
      else
         write(6,*) "In InitializeEvolution.f:"
         write(6,*) "Invalid perturbative order, ipt = ",ipt
         call exit(-10)
      endif
*
*     Factorisation scheme
*
      if(FACSCHEME.eq."MSBAR")then
         write(6,*) "Factorisation scheme: MSbar"
      elseif(FACSCHEME.eq."DELTA")then
         write(6,*) "Factorisation scheme: Delta"
      else
         write(6,*) "In InitializeEvolution.f:"
         write(6,*) "Unknown factorisation scheme = ",FACSCHEME
         call exit(-10)
      endif
*      
*     Renormalisation scheme
*
      if(RENSCHEME.eq."MSBAR")then
         write(6,*) "Renormalisation scheme: MSbar"
         if(iptalpha.eq.0)then
            write(6,*) " Perturbative order alpha: LO"
         elseif(iptalpha.eq.1)then
            write(6,*) " Perturbative order alpha: NLO"
         else
            write(6,*) "In InitializeEvolution.f:"
            write(6,*) "Invalid perturbative order, iptalpha=",iptalpha
            call exit(-10)
         endif
      elseif(RENSCHEME.eq."FIXED")then
         write(6,*) "Renormalisation scheme: alpha fixed"
      elseif((RENSCHEME.eq."ALPMZ").or.(RENSCHEME.eq."ALGMU")
     .        .or.(RENSCHEME.eq."AFAKE"))then
         if(RENSCHEME.eq."ALPMZ")then
            write(6,*) "Renormalisation scheme: alpha(MZ)"
         endif
         if(RENSCHEME.eq."ALGMU")then
            write(6,*) "Renormalisation scheme: alphaGmu"
         endif
         if(RENSCHEME.eq."AFAKE")then
            write(6,*) "Renormalisation scheme: alphaGmu FAKE"
         endif         
         if(AFIXSOL.eq."PATHOR")then
            write(6,*) "  method: path_ordering"
         elseif(AFIXSOL.eq."MAGNUS")then
            write(6,*) "  method: magnus"
         else
            write(6,*) "In InitializeEvolution.f:"
            write(6,*) "Unknown alpha(MZ)/alphaGmu method = ",AFIXSOL
            call exit(-10)
         endif
      else
         write(6,*) "In InitializeEvolution.f:"
         write(6,*) "Unknown renormalisation scheme = ",RENSCHEME
         call exit(-10)
      endif
*
*     Evolution scheme
*
      if((NLMAX.lt.0).or.(NUMAX.lt.0).or.(NDMAX.lt.0))then
         write(6,*) "In InitializeEvolution.f:"
         write(6,*) "Wrong (NLMAX,NUMAX,NDMAX) = ",
     .        NLMAX,NUMAX,NDMAX
      endif
*      
      if(NS.eq."VFNS")then
         write(6,*) "Evolution scheme: VFNS"
         write(6,"(a,i1,a,i1,a,i1)") " NLMAX = ",NLMAX,
     .        ", NUMAX = ",NUMAX,", NDMAX = ",NDMAX
      elseif(NS.eq."FFNS")then
         write(6,*) "Evolution scheme: FFNS"
         write(6,"(a,i1,a,i1,a,i1)") " NL = ",NLMAX,
     .        ", NU = ",NUMAX,", ND = ",NDMAX         
      else
         write(6,*) "In InitializeEvolution.f:"
         write(6,*) "Unknown mass scheme = ",NS
         call exit(-10)
      endif
**            
*     Thresholds
*
      write(6,*) "Mass thresholds:"
      write(6,"(a,f14.9,a)") "   me  =",dsqrt(q2th(1 ))," GeV"
      write(6,"(a,f14.9,a)") "   mu  =",dsqrt(q2th(2 ))," GeV"
      write(6,"(a,f14.9,a)") "   md  =",dsqrt(q2th(3 ))," GeV"
      write(6,"(a,f14.9,a)") "   ms  =",dsqrt(q2th(4 ))," GeV"
      write(6,"(a,f14.9,a)") "   mm  =",dsqrt(q2th(5 ))," GeV"
      write(6,"(a,f14.9,a)") "   mc  =",dsqrt(q2th(6 ))," GeV"
      write(6,"(a,f14.9,a)") "   mt  =",dsqrt(q2th(7 ))," GeV"
      write(6,"(a,f14.9,a)") "   mb  =",dsqrt(q2th(8 ))," GeV"
      write(6,"(a,f14.9,a)") "   MW  =",dsqrt(q2th(9 ))," GeV"
      write(6,"(a,f14.9,a)") "   MZ  =",dsqrt(q2th(10))," GeV"
**      write(6,"(a,f14.9,a)") "   mtp =",dsqrt(q2th(11))," GeV"      
*      
*     Check that thresholds are ordered
*
      do nf = 2, 11
         if (dsqrt(q2th(nf - 1)).gt.dsqrt(q2th(nf))) then
            write(6,*) "In InitializeEvolution.f:"
            write(6,*) "Fermion masses are not ordered"
            call exit(-10)
         endif
      enddo
**
*     Alpha reference values
*
      if(RENSCHEME.ne."MSBAR")then
         write(6,*) "Alpha FIXED with reference value:"
         write(6,"(a,f14.9)")   "   Alpha(Qref) = ",AREF         
      else
         write(6,*) "Alpha reference value:"
         write(6,"(a,f14.9,a)") "   Qref = ",sqrt(Q2REF)," GeV"
         write(6,"(a,f14.9)")   "   Alpha(Qref) = ",AREF
         if((NLMAXAEM.lt.0).or.(NUMAXAEM.lt.0).or.(NDMAXAEM.lt.0))then
            write(6,*) "In InitializeEvolution.f:"
            write(6,*) "Wrong (NLMAXAEM,NUMAXAEM,NDMAXAEM) = ",
     .           NLMAXAEM,NUMAXAEM,NDMAXAEM
         endif         
         if(NS.eq."VFNS")then
            write(6,*) "Alpha evolution scheme: VFNS"
            write(6,"(a,i1,a,i1,a,i1)") " NLMAXAEM = ",NLMAXAEM,
     .           ", NUMAXAEM = ",NUMAXAEM,", NDMAXAEM = ",NDMAXAEM
         elseif(NS.eq."FFNS")then
            write(6,*) "Alpha evolution scheme: FFNS"
            write(6,"(a,i1,a,i1,a,i1)") " NLAEM = ",NLMAXAEM,
     .           ", NUAEM = ",NUMAXAEM,", NDAEM = ",NDMAXAEM         
         else
            write(6,*) "In InitializeEvolution.f:"
            write(6,*) "Unknown mass scheme = ",NS
            call exit(-10)
         endif
      endif
*
      if(waem.eq.0)then
         write(6,*) "W effects NOT included in evolution"
         if(RENSCHEME.eq."ALGMU")then
            write(6,*) "WARNING: GMU SCHEME WITH WOFF MAY BE UNPHYSICAL"
         endif
      elseif(waem.eq.1)then
         write(6,*) "W effects included in evolution"
      else
         write(6,*) "In InitializeEvolution.f:"
         write(6,*) "Unknown Waem value, Waem = ",waem
         call exit(-10)
      endif
*
      return
      end
