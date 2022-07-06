#pragma once

#include <string>
#include <map>
#include <vector>

#include "eMELA/bspdf.hh"

namespace eMELA {

  //////////////////////////////////////////////////////////////////////
  /// Main methods
  //////////////////////////////////////////////////////////////////////  
  /// Set default parameters
  void SetDefaultParameters(void);
  /// Initialize evolution
  void InitializeEvolution(void);

  // Set parameters from grid
  void InitializeFromGrid(std::string const& pdfname);    
  void InitializeFromGrid(int const& lhaid);
  
  // one-call function to initialize by using default parameters or arguments
  void QuickInitialize(std::string const& pert_order = "NLL",
		       std::string const& fac_scheme = "DELTA",
		       std::string const& ren_scheme = "MSBAR",
		       double      const& alpha_mz   = 1.0/127.95471413988483012);
  
  /// Alpha at Q2
  double aQED(double q2);
  /// Alpha at the electron mass
  double aQEDme();

  // Full PDF (switching between numerical from code and analytic z->1)
  double CodePdf(const int & idx, const double & x, const double & omx, const double & Q,
		 const double & gamma = 1.0);
  // Numerical PDF (directly from code)
  double CodeNumPdf(const int & idx, const double & x, const double & Q, const double & gamma = 1.0);

  // Full PDF (switching between numerical from grid and analytic z->1)
  double GridPdf(const int & idx, const double & x, const double & omx, const double & Q,
		 const double & gamma = 1.0);
  // Numerical PDF (reading from grid, if not initialised it throws error)
  double GridNumPdf(const int & idx, const double & x, const double & Q, const double & gamma = 1.0);
  
  // Analytic asymptotic PDF
  double AsyPdf(const double & omx, const double & Q, const double & gamma = 1.0);  

  // Legacy LL PDFs used in the literature see e.g.
  // https://arxiv.org/pdf/hep-ph/0006307.pdf
  // or
  // https://iopscience.iop.org/article/10.1209/0295-5075/17/2/007/pdf
  //
  // betaE is the variable appearing in the asymptotic solution,
  // except betaS which is the coefficient of 3/4 in the exponent.
  // betaH is the variable appearing in the recursive solution.
  //
  // index labelling LL PDFs:
  // - 0: collinear
  // - 1: beta
  // - 2: eta
  // - 3: mixed
  //
  // the normalization is such that in "beta" scheme:
  // - betaE = betaS = betaH = alpha/pi * (L - 1)
  // whereas in the "collinear" scheme i.e. the LL PDF of https://arxiv.org/pdf/1911.12040.pdf
  // - betaE = betaS = betaH = alpha/pi * L
  // where L = log(mu^2/mu0^2)
  //
  // the PDF is multiplied by (1-z)^(1-gamma)
  //  
  // Asymptotic + recursive (without double counting)  
  double LLPDF(int const& index, double const& x, double const& omx, double const& Q,
	       double const& gamma = 1.0);
  // Asymptotic only
  double AsyLLPDF(int const& index, double const& omx, double const& Q,
		  double const& gamma = 1.0);

  /// Write grids of the current PDF
  // name: name of the grid
  // nx: number of x points
  // xmin: lower extreme
  // xmed: separation between lin and log grid
  // ymax: upper extreme 1-10^(-ymax)
  // frac: # points in lin grid = nx * frac, # points in log grid = nx * (1-frac)
  // nQ: number of Q points
  // Qmin: lower extreme
  // Qmax: upper extreme
  // Qlin: use a linear grid in Q
  // bspdf: pointer to an object bspdf for BREMSSTRAHLUNG.
  void WriteGrid(std::string const& name,
		 bspdf*             brem = NULL,		 
		 int         const& nx   = 500,
		 double      const& xmin = 0.001,
		 double      const& xmed = 0.9,
		 double      const& ymax = 9,
		 double      const& frac = 0.4,
		 int         const& nQ   = 100,
		 double      const& Qmin = 0.1,
		 double      const& Qmax = 1000.0,
		 bool        const& Qlin = false);
  
  /// Banner utilities
  void Banner();

  //////////////////////////////////////////////////////////////////////
  /// BREMSSTRAHLUNG
  //////////////////////////////////////////////////////////////////////  
  // Utility function 
  double AsyInt(const double & omx, const double & Q, const double & beta, const double & kappa);
  
  //////////////////////////////////////////////////////////////////////
  /// Setters
  //////////////////////////////////////////////////////////////////////
  /// Perturbative order of the solution
  /// 0: LL, 1:NLL
  void SetPerturbativeOrder(int);
  /// Perturbative order of alpha (when running)
  /// 0: LL, 1:NLL
  void SetPerturbativeOrderAlpha(int);
  /// Flavour scheme (VFNS, FFNS)
  void SetFlavourScheme(std::string const&);
  /// 0: FFNS, 1: VFNS
  void SetFlavourSchemeInt(int);
  /// Factorisation scheme (MSBAR, DELTA)    
  void SetFactorisationScheme(std::string const&);
  /// 0: MSBAR, 1: DELTA
  void SetFactorisationSchemeInt(int);
  /// Renormalisation scheme (MSBAR, FIXED, ALPMZ, ALGMU)  
  void SetRenormalisationScheme(std::string const&);
  /// 0: MSBAR, 1: FIXED, 2: ALPMZ, 3: ALGMU
  void SetRenormalisationSchemeInt(int);
  /// Number of flavours in the evolution NL, NU, ND
  /// (turned on at thresholds in the VFNS, or kept fixed in the FFNS)
  void SetActiveFlavours(int nlmax, int numax, int ndmax);
  /// Number of flavours in the evolution of alpha NL, NU, ND
  /// (turned on at thresholds in the VFNS, or kept fixed in the FFNS)
  void SetActiveFlavoursAlpha(int nlmaxaem, int numaxaem, int ndmaxaem);
  /// Same as above, but with a vector
  void SetActiveFlavours(std::vector<int> active);
  /// Same as above, but with a vector
  void SetActiveFlavoursAlpha(std::vector<int> active);
  /// 0: OFF, 1: ON
  void SetWalpha(int);
  /// Initial value of alpha ain at reference scale Qin
  void SetAlpha(double ain, double Qin);
  /// Set the value of the thresholds
  void SetThresholds(double me, double mu, double md, double ms, double mm,
  		     double mc, double mt, double mb, double MW, double MZ);
  /// As above, but with vector
  void SetThresholds(std::vector<double> thresholds);
  /// TECHNICAL PARAMETERS
  /// Number of total steps in the path ordering solution  
  void SetNPHOTOT(int);
  /// Number of minimum steps in each mass range in the path ordering solution  
  void SetNPHOMIN(int);
  /// Number of terms in the expasion of 4x4 matrices
  void SetNMATEXP(int);
  /// Number of steps in the numerical solution of alpha (at NLO)
  void SetNSTPAEM(int);
  /// Technical parameter for Mellin inverse  
  void SetMINVMEL(int);
  /// Technical parameter for Mellin inverse  
  void SetRINVMEL(int);
  /// Solution to be used in the alpha(MZ) or GMU scheme (MAGNUS, PATHOR)
  void SetAFIXSOLint(int);
  /// Value of deltaGMU, to be used in the GMU scheme, associated to alpha value
  /// WARNING: DEFAULT VALUE ASSOCIATED TO VALUE OF ALPHA = 1/132.18289853516
  /// IF YOU USE A DIFFERENT ALPHA, YOU HAVE TO CHANGE THIS VALUE
  void SetDELTAGMU(double);

  //////////////////////////////////////////////////////////////////////
  /// Getters
  //////////////////////////////////////////////////////////////////////
  /// See the setters above for explanation
  int GetPerturbativeOrder();
  int GetPerturbativeOrderAlpha();
  int GetFlavourSchemeInt();
  int GetFactorisationSchemeInt();
  int GetRenormalisationSchemeInt();
  std::vector<int> GetActiveFlavours();
  std::vector<int> GetActiveFlavoursAlpha();
  int GetWalpha();  
  double GetAlphaRef();
  double GetAlphaQref();
  std::vector<double> GetThresholds();
  std::vector<double> GetThresholds2();
  int GetNPHOTOT();
  int GetNPHOMIN();
  int GetNMATEXP();
  int GetNSTPAEM();
  int GetMINVMEL();
  int GetRINVMEL();
  int GetAFIXSOLint();
  double GetDELTAGMU();

}
