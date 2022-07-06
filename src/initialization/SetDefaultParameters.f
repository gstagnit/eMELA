************************************************************************
*
*     setdefaultparameters.f:
*
*     it sets the default parameters of the evolution.
*
************************************************************************
      subroutine setdefaultparameters
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
*
*     default values
*
      ipt         = 1
      iptalpha    = 1
      ns          = "VFNS"
      facscheme   = "MSBAR"
      renscheme   = "MSBAR"
      nlmax       = 3
      numax       = 2
      ndmax       = 3
      nlmaxaem    = 3
      numaxaem    = 2
      ndmaxaem    = 3
      waem        = 1
      q2ref       = 0.000510998928d0**2
      aref        = 0.0072973525693d0
      q2th(1)     = 0.000510998928d0**2
      q2th(2)     = 0.00216d0**2
      q2th(3)     = 0.00467d0**2
      q2th(4)     = 0.093d0**2
      q2th(5)     = 0.10566d0**2
      q2th(6)     = 1.27d0**2
      q2th(7)     = 1.77686d0**2
      q2th(8)     = 4.18d0**2
      q2th(9)     = 80.379d0**2
      q2th(10)    = 91.1876d0**2
**********************************************************************     
*     SETTING THE TOP THRESHOLD TO VERY LARGE VALUE
      q2th(11)    = 1000000000d0**2
**********************************************************************
      nphotot     = 300
      nphomin     = 50
      nmatexp     = 30
      nstpaem     = 10
      minvmel     = 51
      rinvmel     = 20
      afixsol     = "MAGNUS"
c     WARNING: ASSOCIATED TO VALUE OF ALPHA = 1/132.18289853516
      deltagmu    = 4.09293586102687773d0
*
      return
      end
