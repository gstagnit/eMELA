************************************************************************
*
*     initBeta.f:
*
*     Initialize the coefficients of the QED beta function.
*
************************************************************************
      subroutine initBeta
*
      implicit none
*
      include "../commons/charges.h"      
      include "../commons/activeflavours.h"
      include "../commons/beta.h"
      include "../commons/nfsum.h"
      include "../commons/ns.h"         
**
*     Internal Variables
*
      integer i
      integer nlaem(1:11), nuaem(1:11), ndaem(1:11)
*     
*     fill b0 and b1 according to nlmaxaem,numaxaem,ndmaxaem,waem
*
      call setnlnund(nlmaxaem,numaxaem,ndmaxaem,nlaem,nuaem,ndaem)
      do i = 1,11
         beta0(i) = -4d0/3d0 * ( el2 * nlaem(i)
     .        + nc * ( eu2 * nuaem(i) + ed2 * ndaem(i) ) )
         beta1(i) = -4d0 * ( el4 * nlaem(i)
     .        + nc * ( eu4 * nuaem(i) + ed4 * ndaem(i) ) )
      enddo
*
*     Both in FFNS/VFNS W boson effects only above the W mass
      do i = 9,11
         beta0(i) = -4d0/3d0 * (-3d0/4d0*beta0(i) - 21d0/4d0*WAEM )
      enddo
*     
*     fill nfsum2 and nfsum4 according to nlmax,numax,ndmax
*     fill nl,nu,nd according to nlmax,numax,ndmax
*
      call setnlnund(nlmax,numax,ndmax,nl,nu,nd)
      do i = 1,11
         nfsum2(i) = ( el2 * nl(i)
     .        + nc * ( eu2 * nu(i) + ed2 * nd(i) ) )
         nfsum4(i) = ( el4 * nl(i)
     .        + nc * ( eu4 * nu(i) + ed4 * nd(i) ) )
      enddo
*
      return
      end
************************************************************************      
*     Utility subroutine to fill the nl,nu,nd vectors according to
*     the maximum number of allowed leptons, down- and up-quarks.
      subroutine setnlnund(nlmax,numax,ndmax,nl,nu,nd)
      implicit none
      include "../commons/ns.h"      
*     Input
      integer nlmax, numax, ndmax
*     Output
      integer nl(1:11), nu(1:11), nd(1:11)
*     Internal
      integer i
*
*     if FFNS, than all the elements of the vectors are set to
*     the same values nlmax,numax,ndmax      
      if(ns.eq."FFNS")then
         do i=1,11
            nl(i) = nlmax
            nu(i) = numax
            nd(i) = ndmax
         enddo
         return
      endif
*      
      if(ns.eq."VFNS")then
         do i=1,11
            nl(i) = 0
            nu(i) = 0
            nd(i) = 0
         enddo
*     
         if (nlmax.ge.1) then
            do i=1,11
               nl(i) = 1
            enddo
         endif
         if (nlmax.ge.2) then
            do i=5,11
               nl(i) = 2
            enddo
         endif
         if (nlmax.ge.3) then
            do i=7,11
               nl(i) = 3
            enddo
         endif
*
         if (numax.ge.1) then
            do i=2,11
               nu(i) = 1
            enddo
         endif
         if (numax.ge.2) then
            do i=6,11
               nu(i) = 2
            enddo
         endif
         if (numax.ge.3) then
            do i=11,11
               nu(i) = 3
            enddo
         endif
*     
         if (ndmax.ge.1) then
            do i=3,11
               nd(i) = 1
            enddo
         endif
         if (ndmax.ge.2) then
            do i=4,11
               nd(i) = 2
            enddo
         endif
         if (ndmax.ge.3) then
            do i=8,11
               nd(i) = 3
            enddo
         endif
*         
      endif
*
      return
      end
