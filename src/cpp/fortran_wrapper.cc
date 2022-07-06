#include <string>
#include <cstring>
#include <memory>
#include <cmath>

#include "eMELA/eMELA.hh"
#include "eMELA/bspdf.hh"

namespace eMELA {

  /// a namespace for the fortran-wrapper which contains commonly-used
  /// structures and means to transfer fortran <-> C++
  namespace fwrapper {

#ifdef USE_LHAPDF
    std::shared_ptr<bspdf> ptr_bs;
#endif
    
    /// copied from src/LHAGlue.cc of LHAPDF
    /// Fortran-string -> C++-string converter
    std::string fstr_to_ccstr(const char* fstring, const std::size_t fstring_len,
			      bool spcpad=false) {
      // Allocate space for an equivalent C-string (with an extra terminating null byte)
      char* s = new char[fstring_len+1];
      // Copy all characters and add the terminating null byte
      strncpy(s, fstring, fstring_len);
      s[fstring_len] = '\0';
      // Replace all trailing spaces with null bytes unless explicitly stopped
      if (!spcpad) {
	for (int i = fstring_len-1; i >= 0; --i) {
	  if (s[i] != ' ') break;
	  s[i] = '\0';
	}
      }
      std::string rtn(s); //< copy the result to a C++ string
      delete[] s; //< clean up the dynamic array
      return rtn;
    }

  }
}
  

using namespace eMELA::fwrapper;

extern "C" {   

#ifdef USE_LHAPDF
  /// Initialize the PDF object from grid
  ///
  /// Corresponds to the following Fortran subroutine interface structure:
  ///
  ///    subroutine initfromgrid_name(pdfname)
  ///    character*N      pdfname
  ///
  /// where on input
  ///
  ///     pdfname   name of the pdf 
  ///               (to be searched in $LHAPDF_DATA_PATH)
  ///
  ///  The second argument is not present in the Fortran call
  ///  (Fortran passes an additional argument with the length).
  ///
  void initfromgrid_name_(const char* pdfname, int pdfname_length)
  {
    std::string p = fstr_to_ccstr(pdfname, pdfname_length);
    std::cout << "Reading grid with name: " << p << std::endl;
    eMELA::InitializeFromGrid(p);
  }

  /// Initialize the PDF object from grid
  ///
  /// Corresponds to the following Fortran subroutine interface structure:
  ///
  ///    subroutine initfromgrid_lhaid(lhaid)
  ///    integer      lhaid
  ///
  /// where on input
  ///
  ///     lhaid   index of the pdf in the pdfsets.index file
  ///             (to be searched in $LHAPDF_DATA_PATH)
  ///  
  void initfromgrid_lhaid_(int const& lhaid)
  {
    eMELA::InitializeFromGrid(lhaid);
  }
#endif
  
  /// Retrieve the value of alpha(Q2)
  ///
  /// Corresponds to the following Fortran subroutine interface structure:
  ///
  ///    subroutine alphaQ2(Q2,aem)
  ///    double precision Q2,aem
  ///
  void alphaq2_(double const& Q2, double & aem)
  {
    aem = eMELA::aQED(Q2);
  }		 
  
  /// Retrieve the electron PDF matched with analytics in the z->1 region
  ///
  /// Corresponds to the following Fortran subroutine interface structure:
  ///
  ///    subroutine elpdfq2(usegrids,idx,x,omx,q2,gamma,f)
  ///    integer usegrids,idx
  ///    double precision x,omx,q2,gamma,f
  ///
  /// where on input
  ///
  ///    usegrids   0: use grids, 1: call directly the numerical/analytical code  
  ///    idx        11: electron, 22: photon, -11: positron
  ///    x          momentum fraction of the PDF
  ///    omx        1-x (used in the asymptotic region)
  ///    q2         final scale squared
  ///    gamma      the PDF returned is f(x)*(1-x)^(1-gamma)
  ///
  /// and on output
  ///
  ///    f          the value of the PDF f(x,q2)
  ///
  void elpdfq2_(int    const& usegrids,
		int    const& idx,
		double const& x,
		double const& omx,
		double const& q2,
		double const& gamma,
		double      & f)
  {
    const bool bgrids = (usegrids == 0) ? true : false;
    if (bgrids) {
      f = eMELA::GridPdf(idx, x, omx, sqrt(q2), gamma);
    } else {
      f = eMELA::CodePdf(idx, x, omx, sqrt(q2), gamma);
    }
  }

#ifdef USE_LHAPDF
  ////////////////////////////////////////////////////////////////////////////////
  /// BREMSSTRAHLUNG
  ////////////////////////////////////////////////////////////////////////////////
  void bs_initfromgrid_name_(const char* pdfname, int pdfname_length)
  {
    std::string p = fstr_to_ccstr(pdfname, pdfname_length);
    std::cout << "Reading grid with bremsstrahlung with name: " << p << std::endl;
    ptr_bs = make_bspdf(p);    
  }

  void bs_initfromgrid_lhaid_(int const& lhaid)
  {
    ptr_bs = make_bspdf(lhaid);    
  }
    
  void bs_elpdfq2_(int    const& icom,
		   int    const& pid,
		   int    const& ibeam,		   
		   double const& x,
		   double const& omx,
		   double const& q2,
		   double const& gamma,
		   double      & f)
  {
    f = ptr_bs->get_pdf(icom, pid, ibeam, x, omx, sqrt(q2), gamma);
  }
#endif
  
}




