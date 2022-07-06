#pragma once

namespace eMELA
{
  double Asymptotic( double const& omz, double const& Q, double const& betab,
		     std::vector<double> const& athrs = std::vector<double>(),
		     bool const& bson = false, double const& kappa = 0.0);
  
  //////////////////////////////////////////////////////////////////////
  /// Used in the analytic expressions
  //////////////////////////////////////////////////////////////////////
  // Return the number of the region of mu2
  // i.e. 1 if me2 <= mu2 < mu2
  //      2 if mu2 <= mu2 < md2
  //      ...
  //      8 if mb2 <= mu2 < mw2
  //      9 if mw2 <= mu2 < mz2
  //     10 if mz2 <= mu2 < mtp2
  //     11 if mtp2 <= mu2
  int GetRegionMU2(double mu2);

  // Get each of the expressions in the relevant region
  double GetC2(int region);
  double GetC4(int region);
  double Getb0(int region);
  double Getb1(int region);
  double GetDk(int region);
  double GetMZ2();
  double GetMW2();
  
}



