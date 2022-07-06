#include <cmath>
#include <memory>
#include <iostream>

#include "eMELA/eMELA.hh"
#include "eMELA/version.h"
#include "eMELA/analytics.hh"
#include "eMELA/numerics.hh"

#ifdef USE_LHAPDF
#include "eMELA/grid.hh"
#endif

using namespace std;

extern "C"
{
  void initializeevolution_(void);
}  


namespace eMELA {

  // Show banner just the first time
  // TODO: When do we show the banner?
  static bool banner_first_time = false;
  
  /// Struct to store parameters for switching
  struct switchparams
  {
    double omxcut = 1e-8;    
    int kswitch = 2;
  };
  
  /// Pointer to the local objects
#ifdef USE_LHAPDF  
  unique_ptr<GridMELA> grid = nullptr;
#endif
  
  unique_ptr<vector<double>> athrs = nullptr;
  unique_ptr<switchparams> swpar = nullptr;
  
  //_________________________________________________________________________________
  /// local function
  /// When building the Grid object, save the values of alpha at the thresholds
  /// instead of recomputing them at any time.  
  void save_athrs()
  {
    const std::vector<double> q2thrs = GetThresholds2();
    athrs.reset(new vector<double>);
    for (unsigned i=0; i < q2thrs.size(); i++) {
      athrs->push_back(aQED(q2thrs[i]));
    }
  }

  //_________________________________________________________________________________
  void InitializeEvolution()
  {
    initializeevolution_();    
    save_athrs();
    swpar.reset(new switchparams);
  }

  //_________________________________________________________________________________
  void QuickInitialize(std::string const& pert_order,
		       std::string const& fac_scheme,
		       std::string const& ren_scheme,
		       double      const& alpha_mz)
  {
    SetDefaultParameters();
    if (pert_order == "NLL")     SetPerturbativeOrder(1);
    else if (pert_order == "LL") SetPerturbativeOrder(0);
    else throw std::runtime_error("[QuickInitialize]: wrong perturbative order");
    SetFactorisationScheme(fac_scheme);
    SetRenormalisationScheme(ren_scheme);
    SetAlpha(alpha_mz, 91.1876);
    InitializeEvolution();    
  }

#ifdef USE_LHAPDF
  //_________________________________________________________________________________
  void InitializeFromGrid(std::string const& pdfname)
  {
    SetDefaultParameters();
    grid.reset(new GridMELA(pdfname));
    InitializeEvolution();    
  }
  
  //_________________________________________________________________________________
  void InitializeFromGrid(int const& lhaid)
  {
    SetDefaultParameters();       
    grid.reset(new GridMELA(lhaid));
    InitializeEvolution();    
  }
#endif
  
  //_________________________________________________________________________________
  double AsyPdf(const double & omx, const double & Q, const double & gamma)
  {
    return Asymptotic(omx, Q, gamma, *athrs);
  }
  
  //_________________________________________________________________________________
  double AsyInt(const double & omx, const double & Q, const double & gamma, const double & kappa)
  {
    return Asymptotic(omx, Q, gamma, *athrs, true, kappa);
  }

  //_________________________________________________________________________________
  double GridNumPdf(const int & idx, const double & x, const double & Q,
		    const double & gamma)
  {
#ifndef USE_LHAPDF
    throw std::runtime_error("[GridNumPdf]: usage of grids requires LHAPDF");
#else
    if (grid == nullptr)
      throw std::runtime_error("[GridNumPdf] No grid has been initialised");
    
    return grid->NumGrid(idx, x, Q, gamma);
#endif
  }
  
  //_________________________________________________________________________________
  int pdg2mela(int ipdg) {
    int imela;
    switch (ipdg) {
    case  -5:  imela = -9; break; // bbar     
    case  -3:  imela = -8; break; // sbar     
    case  -1:  imela = -7; break; // dbar     
    case  -6:  imela = -6; break; // tbar       
    case  -4:  imela = -5; break; // cbar     
    case  -2:  imela = -4; break; // ubar     	
    case -15:  imela = -3; break; // tau+
    case -13:  imela = -2; break; // mu+
    case -11:  imela = -1; break; // e+
    case  22:  imela =  0; break; // photon
    case  11:  imela =  1; break; // e- 
    case  13:  imela =  2; break; // mu-
    case  15:  imela =  3; break; // tau-
    case   2:  imela =  4; break; // u
    case   4:  imela =  5; break; // c
    case   6:  imela =  6; break; // t
    case   1:  imela =  7; break; // d
    case   3:  imela =  8; break; // s     
    case   5:  imela =  9; break; // b
    default: throw std::runtime_error("[pdg2mela] unknown PDG code");
    }
    return imela;
  }
    
  //_________________________________________________________________________________  
  double CodeNumPdf(const int & idx, const double & x, const double & Q,
		    const double & gamma)
  {
    std::map<int, double> xpdf = xDistributions(x,Q);
    const double power = pow(1.0-x, 1.0-gamma);
    int imela = pdg2mela(idx);
    return xpdf.at(imela)/x*power;
  }

  //_________________________________________________________________________________
  void SetSwitchParameters(const double & omxcut, const int & kswitch)
  {
    swpar->omxcut = omxcut;
    swpar->kswitch = kswitch;
  }
  
  //_________________________________________________________________________________
  double _SwitchPdf(const bool & usegrid,
		    const int & idx, const double & x, const double & omx, const double & Q,
		    const double & gamma)
  {
    double res = 0.0;
    if (idx == 11 && omx < swpar->omxcut) {
      /// Analytic with rescaling if electron and x close to 1
      double xcut = 1.0-swpar->omxcut;
      double num0 = (usegrid) ? GridNumPdf(idx, xcut, Q, gamma) : CodeNumPdf(idx, xcut, Q, gamma);
      double asy0 = AsyPdf(swpar->omxcut, Q, gamma);
      double asy  = AsyPdf(omx, Q, gamma);
      // Try new switch with smooth function
      if (swpar->kswitch == 0) {
	res = asy*(num0/asy0);
      } else {
	double sfun = pow(omx/swpar->omxcut,swpar->kswitch);
	res = asy*( sfun*num0/asy0 + 1 - sfun );
      }
      // std::cout << "ASYMTC x: " << x << " omx: " << omx << " elpdf: " << res << std::endl;      
    } else {
      /// Numerical
      res = (usegrid) ? GridNumPdf(idx, x, Q, gamma) : CodeNumPdf(idx, x, Q, gamma);
      // std::cout << "NUM x: " << x << " omx: " << omx << " elpdf: " << res << std::endl;  
    }    
    return res;   
  }
  
  //_________________________________________________________________________________    
  double GridPdf(const int & idx, const double & x, const double & omx, const double & Q,
		 const double & gamma)
  {
#ifndef USE_LHAPDF
    throw std::runtime_error("[GridPdf]: usage of grids requires LHAPDF");
#else
    if (grid == nullptr)
      throw std::runtime_error("[GridPdf] No grid has been initialised");
    
    return _SwitchPdf(true, idx, x, omx, Q, gamma);    
#endif
  }

  //_________________________________________________________________________________      
  double CodePdf(const int & idx, const double & x, const double & omx, const double & Q,
		 const double & gamma)
  {
    return _SwitchPdf(false, idx, x, omx, Q, gamma);
  }
  
  //_________________________________________________________________________________      
  void Banner()
  {
    if (banner_first_time)  return;
    banner_first_time = true;
    
    //std::cout << "\033[1;33m\n";
    std::cout << "# ----------------------------------------------------------------------------\n";
    std::cout << "#                               ***************                               \n";    
    std::cout << "#                               **** eMELA ****                               \n";
    std::cout << "#                                     v" << VERSION <<                       "\n";
    std::cout << "#                               ***************                               \n";
    std::cout << "# A library for the evolution of the electron PDFs in QED at NLL              \n";
    std::cout << "# Webpage: https://github.com/gstagnit/eMELA                                  \n";
    std::cout << "# If you use this code for a scientific publication, please cite:             \n";
    std::cout << "# arXiv:1909.03886, arXiv:1911.12040, arXiv:2105.06688, arXiv:2206.XXXXX      \n";
    std::cout << "# ----------------------------------------------------------------------------\n";
    //std::cout << "\033[39m\n";
  }
  
  
}
