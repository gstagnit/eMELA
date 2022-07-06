************************************************************************
*     Convert from the LHA convention ordering:
*
*      -9  -8  -7  -6  -5  -4  -3  -2  -1   0   1   2   3   4   5   6   7   8   9
*      bb  sb  db  tb  cb  ub  ta+ mu+ e+   gm  e-  mu- tu- u   c   t   d   s   b
*
*     to the EVOLUTION convention:
*
*      1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19
*     gm Sgl T1l T2l  Vl V1l V2l Sgu T1u T2u  Vu V1u V2u Sgd T1d T2d  Vd V1d V2d
*
************************************************************************
      subroutine lha2evln(fin,fout)
*
      implicit none
*
      include "../commons/rotations.h"
**
*     Input Variables
*
      double precision fin(-9:9)
**
*     Internal Variables
*
      integer i,j
**
*     Output Variables
*
      double precision fout(19)
*
*     Intialise to zero
*
      do i=1,19
         fout(i) = 0d0
      enddo

*     Gamma
      fout(1) = fin(0)
*
      do i=1,3
         do j=1,3
*     Leptons
            fout(i+1) = fout(i+1) + p2e(j,i) * (fin(j) + fin(-j))
            fout(i+4) = fout(i+4) + p2e(j,i) * (fin(j) - fin(-j))
*     Up-quarks
            fout(i+7)  = fout(i+7)  + p2e(j,i) * (fin(j+3) + fin(-j-3))
            fout(i+10) = fout(i+10) + p2e(j,i) * (fin(j+3) - fin(-j-3))
*     Down-quarks
            fout(i+13) = fout(i+13) + p2e(j,i) * (fin(j+6) + fin(-j-6))
            fout(i+16) = fout(i+16) + p2e(j,i) * (fin(j+6) - fin(-j-6))
         enddo
      enddo
*
      return
      end
*
************************************************************************
*     Convert from the EVOLUTION convention ordering:
*
*      1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19
*     gm Sgl T1l T2l  Vl V1l V2l Sgu T1u T2u  Vu V1u V2u Sgd T1d T2d  Vd V1d V2d
*
*     to the LHA convention:
*
*      -9  -8  -7  -6  -5  -4  -3  -2  -1   0   1   2   3   4   5   6   7   8   9
*      bb  sb  db  tb  cb  ub  ta+ mu+ e+   gm  e-  mu- tu- u   c   t   d   s   b
*
************************************************************************
      subroutine evln2lha(fin,fout)
*
      implicit none
*
      include "../commons/rotations.h"
**
*     Input Variables
*
      double precision fin(19)
**
*     Internal Variables
*
      integer i,j
**
*     Output Variables
*
      double precision fout(-9:9)
*
*     Intialise to zero
*
      do i=-9,9
         fout(i) = 0d0
      enddo

*     Gamma
      fout(0) = fin(1)
*
      do i=1,3
         do j=1,3
*     Leptons
            fout(i)  = fout(i)  + e2p(j,i) * (fin(j+1) + fin(j+4))
            fout(-i) = fout(-i) + e2p(j,i) * (fin(j+1) - fin(j+4))
*     Up-quarks
            fout(i+3)  = fout(i+3)  + e2p(j,i) * (fin(j+7) + fin(j+10))
            fout(-i-3) = fout(-i-3) + e2p(j,i) * (fin(j+7) - fin(j+10))
*     Down-quarks
            fout(i+6)  = fout(i+6)  + e2p(j,i) * (fin(j+13) + fin(j+16))
            fout(-i-6) = fout(-i-6) + e2p(j,i) * (fin(j+13) - fin(j+16))
         enddo
      enddo
*
      return
      end
*
************************************************************************
*     Convert from the EVOLUTION convention ordering (complex version):
*
*      1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19
*     gm Sgl T1l T2l  Vl V1l V2l Sgu T1u T2u  Vu V1u V2u Sgd T1d T2d  Vd V1d V2d
*
*     to the LHA convention:
*
*      -9  -8  -7  -6  -5  -4  -3  -2  -1   0   1   2   3   4   5   6   7   8   9
*      bb  sb  db  tb  cb  ub  ta+ mu+ e+   gm  e-  mu- tu- u   c   t   d   s   b
*
************************************************************************
      subroutine evln2lhac(fin,fout)
*
      implicit none
*
      include "../commons/rotations.h"
**
*     Input Variables
*
      double complex fin(19)
**
*     Internal Variables
*
      integer i,j
**
*     Output Variables
*
      double complex fout(-9:9)
*
*     Intialise to zero
*
      do i=-9,9
         fout(i) = 0d0
      enddo

*     Gamma
      fout(0) = fin(1)
*
      do i=1,3
         do j=1,3
*     Leptons
            fout(i)  = fout(i)  + e2p(j,i) * (fin(j+1) + fin(j+4))
            fout(-i) = fout(-i) + e2p(j,i) * (fin(j+1) - fin(j+4))
*     Up-quarks
            fout(i+3)  = fout(i+3)  + e2p(j,i) * (fin(j+7) + fin(j+10))
            fout(-i-3) = fout(-i-3) + e2p(j,i) * (fin(j+7) - fin(j+10))
*     Down-quarks
            fout(i+6)  = fout(i+6)  + e2p(j,i) * (fin(j+13) + fin(j+16))
            fout(-i-6) = fout(-i-6) + e2p(j,i) * (fin(j+13) - fin(j+16))
         enddo
      enddo
*
      return
      end
