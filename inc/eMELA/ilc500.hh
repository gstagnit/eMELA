#pragma once

#include "eMELA/bspdf.hh"

class ilc500:public bspdf
{
public:
  ilc500();
  ilc500(const std::string& name);
  ilc500(int const& lhaid);  

  virtual int components()
  {
    return 4;
  }

  virtual double calc_pdf(int icom, int particleid, int beamid, double x, double omx, double Q, double beta);
  virtual double asy_pdf(int icom, int particleid, int beamid, double x, double omx, double Q, double beta);

  virtual double bs_integral(int icom, int beamid);

  virtual std::string beamspectrum_info();

  virtual std::string bsname() { return "ILC500"; }
   
protected:

  void init();
  
  //constants for ilc500 fit
  double kappa = 1./3.;
  double f11 = 0.5012;
  double f01 = 0.1613;
  double p01 = -8.514;
  double q01 = -5.808;
  double f10 = 0.1613;
  double p10 = -8.505;
  double q10 = -5.823;
  double f001 = 0.2528;
  double p001 = -7.535;
  double q001 = -6.790;
  double f002 = 0.2524;
  double p002 = -7.481;
  double q002 = -6.849;

  double integrate_isr(int particleid, double p, double q, double x, double omx, double Q, double beta);
  double asy_integrate_isr(int particleid, double p, double q, double x, double omx, double Q, double beta);
};
