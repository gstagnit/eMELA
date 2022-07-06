#include <vector>
#include <math.h>

#include "eMELA/eMELA.hh"
#include "eMELA/analytics.hh"
#include "eMELA/numerics.hh"

extern"C"
{
  // initialise to default
  void setdefaultparameters_(void);

  // get alpha and PDFs
  double aqed_(double*);
  void xdistributions_(double*, double*, double*);
  void xdistributionsev_(double*, double*, double*);    
  
  // setters
  void setperturbativeorder_(int*);
  void setperturbativeorderalpha_(int*);
  void setflavourscheme_(char*);
  void setflavourschemeint_(int*);
  void setfactorisationscheme_(char*);
  void setfactorisationschemeint_(int*);
  void setrenormalisationscheme_(char*);
  void setrenormalisationschemeint_(int*);
  void setactiveflav_(int*, int*, int*);
  void setactiveflavaem_(int*, int*, int*);    
  void setwaem_(int*);
  void setalpha_(double*, double*);
  void setthresholds_(double*, double*, double*, double*, double*,
		      double*, double*, double*, double*, double*);
  void setnphotot_(int*);
  void setnmatexp_(int*);
  void setnphomin_(int*);
  void setnstpaem_(int*);
  void setminvmel_(int*);
  void setrinvmel_(int*);
  void setafixsolint_(int*);
  void setdeltagmu_(double*);
  
  // getters
  void getperturbativeorder_(int*);    
  void getperturbativeorderalpha_(int*);
  void getflavourschemeint_(int*);
  void getfactorisationschemeint_(int*);
  void getrenormalisationschemeint_(int*);
  void getactiveflav_(int*, int*, int*);
  void getactiveflavaem_(int*, int*, int*);
  void getwaem_(int*);
  void getalpharef_(double*);
  void getalphaqref_(double*);
  void getthresholds2_(double*);
  void getnphotot_(int*);
  void getnmatexp_(int*);
  void getnphomin_(int*);
  void getnstpaem_(int*);
  void getminvmel_(int*);
  void getrinvmel_(int*);
  void getafixsolint_(int*);
  void getdeltagmu_(double*);

  // used in analytical expressions
  void getc2_(double*);
  void getc4_(double*);
  void getb0_(double*);
  void getb1_(double*);
  void getdk_(int*, double*);  
  void getmz2_(double*);
  void getmw2_(double*);

  void geta0_(double*);
}

namespace eMELA
{
  void SetDefaultParameters(void)
  {
    setdefaultparameters_();
  }

  double aQED(double q2)
  {
    return aqed_(&q2)*4.0*M_PI;
  }  

  double aQEDme()
  {
    double a0;
    geta0_(&a0);
    return a0*4.0*M_PI;
  }
  
  /// PDFs (from -9 to 9)
  std::map<int, double> xDistributions(double x, double Q)
  {
    double* xfph = new double[19];
    xdistributions_(&x, &Q, xfph);
    std::map<int, double> xfout;
    for (int i = -9; i < 10; i++)
      xfout.insert({i, xfph[i+9]});
    delete[] xfph;
    return xfout;
  }

  /// PDFs (from 1 to 19)
  std::map<int, double> xDistributionsEv(double x, double Q)
  {
    double* xfev = new double[19];
    xdistributionsev_(&x, &Q, xfev);
    std::map<int, double> xfout;
    for (int i = 1; i < 20; i++)
      xfout.insert({i, xfev[i-1]});
    delete[] xfev;
    return xfout;
  }

  
  void SetPerturbativeOrder(int iptin)
  {
    setperturbativeorder_(&iptin);
  }

  void SetPerturbativeOrderAlpha(int iptalphain)
  {
    setperturbativeorderalpha_(&iptalphain);
  }
  
  int GetPerturbativeOrder()
  {
    int pto;
    getperturbativeorder_(&pto);
    return pto;
  }

  int GetPerturbativeOrderAlpha()
  {
    int iptalpha;
    getperturbativeorderalpha_(&iptalpha);
    return iptalpha;
  }

  void SetFlavourScheme(std::string const& fnsin)
  {
    std::vector<char> cstr(fnsin.c_str(), fnsin.c_str() + fnsin.size() + 1);
    setflavourscheme_(cstr.data());
  }

  void SetFlavourSchemeInt(int fnsin)
  {
    setflavourschemeint_(&fnsin);
  }
  
  int GetFlavourSchemeInt()
  {
    int fscheme;
    getflavourschemeint_(&fscheme);
    return fscheme;
  }

  void SetFactorisationScheme(std::string const& fsin)
  {
    std::vector<char> cstr(fsin.c_str(), fsin.c_str() + fsin.size() + 1);
    setfactorisationscheme_(cstr.data());
  }

  void SetFactorisationSchemeInt(int fnsinint)
  {
    setfactorisationschemeint_(&fnsinint);
  }

  int GetFactorisationSchemeInt()
  {
    int fscheme;
    getfactorisationschemeint_(&fscheme);
    return fscheme;
  }

  void SetRenormalisationScheme(std::string const& rdin)
  {
    std::vector<char> cstr(rdin.c_str(), rdin.c_str() + rdin.size() + 1);
    setrenormalisationscheme_(cstr.data());
  }

  void SetRenormalisationSchemeInt(int rnsinint)
  {
    setrenormalisationschemeint_(&rnsinint);
  }

  int GetRenormalisationSchemeInt()
  {
    int rdcheme;
    getrenormalisationschemeint_(&rdcheme);
    return rdcheme;
  }

  void SetActiveFlavours(int nlmax, int numax, int ndmax)
  {
    setactiveflav_(&nlmax, &numax, &ndmax);    
  }

  void SetActiveFlavoursAlpha(int nlmaxaem, int numaxaem, int ndmaxaem)
  {
    setactiveflavaem_(&nlmaxaem, &numaxaem, &ndmaxaem);    
  }
  
  void SetActiveFlavours(std::vector<int> active)
  {
    int & nlmax = active[0];
    int & numax = active[1];
    int & ndmax = active[2];    
    setactiveflav_(&nlmax, &numax, &ndmax);
  }
  
  void SetActiveFlavoursAlpha(std::vector<int> active)
  {
    int & nlmaxaem = active[0];
    int & numaxaem = active[1];
    int & ndmaxaem = active[2];    
    setactiveflavaem_(&nlmaxaem, &numaxaem, &ndmaxaem);
  }
  
  std::vector<int> GetActiveFlavours()
  {
    int nlmax, numax, ndmax;
    getactiveflav_(&nlmax, &numax, &ndmax);
    std::vector<int> res;
    res.push_back(nlmax);
    res.push_back(numax);
    res.push_back(ndmax);
    return res;
  }

  std::vector<int> GetActiveFlavoursAlpha()
  {
    int nlmax, numax, ndmax;
    getactiveflavaem_(&nlmax, &numax, &ndmax);
    std::vector<int> res;
    res.push_back(nlmax);
    res.push_back(numax);
    res.push_back(ndmax);
    return res;
  }    

  void SetWalpha(int walpha)
  {
    setwaem_(&walpha);
  }
  
  int GetWalpha()
  {
    int walpha;
    getwaem_(&walpha);
    return walpha;
  }

  void SetAlpha(double ain, double Qin)
  {
    setalpha_(&ain, &Qin);
  }

  double GetAlphaRef()
  {
    double res;
    getalpharef_(&res);
    return res;
  }

  double GetAlphaQref()
  {
    double res;
    getalphaqref_(&res);
    return res;
  }
  
  void SetThresholds(double me, double mu, double md, double ms, double mm,
		     double mc, double mt, double mb,
		     double MW, double MZ)
  {
    setthresholds_(&me, &mu, &md, &ms, &mm, &mc, &mt, &mb, &MW, &MZ);
  }

  void SetThresholds(std::vector<double> thrs)
  {
    double & me = thrs[0];
    double & mu = thrs[1];
    double & md = thrs[2];
    double & ms = thrs[3];
    double & mm = thrs[4];
    double & mc = thrs[5];
    double & mt = thrs[6];
    double & mb = thrs[7];
    double & MW = thrs[8];
    double & MZ = thrs[9];
    
    setthresholds_(&me, &mu, &md, &ms, &mm, &mc, &mt, &mb, &MW, &MZ);
  }
  
  std::vector<double> GetThresholds()
  {
    double* q2thrsf = new double[11];
    getthresholds2_(q2thrsf);
    std::vector<double> thrs;
    for (int i = 0; i < 11; i++)
      thrs.push_back(sqrt(q2thrsf[i]));
    delete[] q2thrsf;
    return thrs;
  }  

  std::vector<double> GetThresholds2()
  {
    double* q2thrsf = new double[11];
    getthresholds2_(q2thrsf);
    std::vector<double> q2thrs;
    for (int i = 0; i < 11; i++)
      q2thrs.push_back(q2thrsf[i]);
    delete[] q2thrsf;
    return q2thrs;
  }  

  void SetNPHOTOT(int res)
  {
    setnphotot_(&res);
  }

  int GetNPHOTOT()
  {
    int res;
    getnphotot_(&res);
    return res;
  }

  void SetNPHOMIN(int res)
  {
    setnphomin_(&res);
  }

  int GetNPHOMIN()
  {
    int res;
    getnphomin_(&res);
    return res;
  }

  void SetNMATEXP(int res)
  {
    setnmatexp_(&res);
  }

  int GetNMATEXP()
  {
    int res;
    getnmatexp_(&res);
    return res;
  }

  void SetNSTPAEM(int res)
  {
    setnstpaem_(&res);
  }

  int GetNSTPAEM()
  {
    int res;
    getnstpaem_(&res);
    return res;
  }

  void SetMINVMEL(int res)
  {
    setminvmel_(&res);
  }

  int GetMINVMEL()
  {
    int res;
    getminvmel_(&res);
    return res;
  }

  void SetRINVMEL(int res)
  {
    setrinvmel_(&res);
  }

  int GetRINVMEL()
  {
    int res;
    getrinvmel_(&res);
    return res;
  }

  void SetAFIXSOLint(int intin)
  {
    setafixsolint_(&intin);
  }

  int GetAFIXSOLint()
  {
    int res;
    getafixsolint_(&res);
    return res;
  }

  void SetDELTAGMU(double din)
  {
    setdeltagmu_(&din);
  }

  double GetDELTAGMU()
  {
    double res;
    getdeltagmu_(&res);
    return res;
  }
  
  int GetRegionMU2(double mu2)
  {
    std::vector<double> q2thrs = GetThresholds2();
    int region = 10;
    region = region - 1;
    while (region >= 0) {
      if (mu2 > q2thrs[region]) {
    	break;
      } else {
    	region = region - 1;
      }
    }
    return region;
  }
  
  double GetC2(int region)
  {
    double* C2f = new double[10];
    getc2_(C2f);
    double res =  C2f[region];
    delete[] C2f;
    return res;
  }

  double GetC4(int region)
  {
    double* C4f = new double[10];
    getc4_(C4f);
    double res = C4f[region];
    delete[] C4f;
    return res;
  }

  double Getb0(int region)
  {
    double* b0f = new double[10];
    getb0_(b0f);
    double res = b0f[region];
    delete[] b0f;
    return res;
  }

  double Getb1(int region)
  {
    double* b1f = new double[10];
    getb1_(b1f);
    double res = b1f[region];
    delete[] b1f;
    return res;
  }

  double GetDk(int region)
  {
    // C++ and Fortran indices are shifted by one
    region += 1;
    double dk;
    getdk_(&region, &dk);
    return dk;
  }
  
  double GetMZ2()
  {
    double mz2;
    getmz2_(&mz2);
    return mz2;
  }

  double GetMW2()
  {
    double mw2;
    getmw2_(&mw2);
    return mw2;
  }
  
}
