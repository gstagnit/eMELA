*     -*-fortran-*-
*
*     active flavours of leptons, up-type and down-type quarks
*
      integer nlmax, numax, ndmax
      integer nlmaxaem, numaxaem, ndmaxaem
      integer nl(1:11), nu(1:11), nd(1:11)
      integer waem
*
      common / activeflavoursmela / nl,nu,nd
      common / nmaxmela / nlmax,numax,ndmax,nlmaxaem,numaxaem,ndmaxaem
      common / waemmela / waem
