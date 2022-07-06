#include <vector>
#include <cmath>
#include <iostream>

#include "eMELA/eMELA.hh"
#include "eMELA/constants.hh"
#include "eMELA/specialfunctions.hh"
#include "eMELA/analytics.hh"

namespace eMELA
{

  // low-level function, containing analytic expressions in the asymptotic region
  // multiplied by (1-z)^(1-betab)
  // kappa is used with bremsstrahlung
  double Asymptotic(double const& omz, double const& Q, double const& betab, std::vector<double> const& athrs, bool const& bson, double const& kappa)
  {
    // Alpha is fixed if we are in the FIXED(1) or ALPMZ(2) or ALGMU(3) renormalization scheme
    const bool aemfix = (GetRenormalisationSchemeInt() > 0) ? true : false;

    const double lambda0 = 3.0/4.0;   

    const double MZ2 = GetMZ2();
    // const double MW2 = GetMW2();
    
    double Q2 = Q*Q;
    double alphaMU = aQED(Q2);
    
    double log2mumz = pow(log(Q2/MZ2),2);
    
    const std::vector<double> q2thrs = GetThresholds2();
    const double me2 = q2thrs[0];

    std::vector<double> myathrs;
    if (athrs.empty()) {
      for (unsigned i=0; i < q2thrs.size(); i++) {
	myathrs.push_back(aQED(q2thrs[i]));
      }    
    } else if (athrs.size() != q2thrs.size()) {
      throw std::runtime_error("[Asymptotic]: athrs vector has wrong dimension!");
    } else {
      myathrs = athrs;
    }

    double alphaM0 = myathrs[0]; // aQEDinit();
        
    const double aref = GetAlphaRef();
    const double ln2 = log(Q2/me2);
    const double eta0 = aref/M_PI * ln2;

    // TODO: we are always assuming that mu0 is equal to the electron mass
    const double log2mu0mz = pow(log(me2/MZ2),2);
    const double logmu0 = 0.0;      
    
    //----------------------------------------------------------------------
    if (GetPerturbativeOrder() == 0) { // LL

      double tsum(0.0);
      
      // at LL there is no dependence on thresholds if alpha is fixed
      if (!aemfix) {
	
	if (GetFlavourSchemeInt() == 0 && GetWalpha() == 0) { // FFNS
	
	  int nfa = GetRegionMU2(Q2);
	  double b0 = Getb0(nfa);
	  tsum = 1./(2.*b0*M_PI) * log(alphaMU/alphaM0);
	
	} else if (GetFlavourSchemeInt() == 1 || (GetFlavourSchemeInt() == 0 && GetWalpha() == 1) ) { // VFNS
	
	  unsigned k = GetRegionMU2(Q2);
	  if (k >= q2thrs.size())
	    throw std::runtime_error("[Asymptotic]: ERROR");
	  double alphaMK = myathrs[k]; // aQED(q2thrs[k]);
	  tsum = 1./(2.*Getb0(k)*M_PI) * log(alphaMU/alphaMK);
	  for (unsigned i=0; i<k; i++) {
	    double alphaMI = myathrs[i]; // aQED(q2thrs[i]);
	    double alphaMIP1 = myathrs[i+1]; // aQED(q2thrs[i+1]);
	    double barti = 1./(2.*Getb0(i)*M_PI) * log(alphaMIP1/alphaMI);
	    tsum += barti;
	  }
	  
	} else {
	  throw std::runtime_error("[Asymptotic]: ERROR");
	} // END OF FLAV SCHEME
	
      }
	
      double csi0 = (aemfix) ? eta0 : 2.0*tsum;
      double csihat0 = (aemfix) ? eta0*lambda0 : 3.0/2.0*tsum;
      double power = pow(omz, csi0 - betab);
      ////////////////////////////////////////////////////////////////////////////////
      /// BREMSSTRAHLUNG
      double bsint = (bson) ? pow(omz, kappa) * tgamma(kappa) * tgamma(csi0) / tgamma(kappa+csi0) : 1.0;
      ////////////////////////////////////////////////////////////////////////////////      
      return exp( csihat0 - emc * csi0 ) / tgamma( 1 + csi0 ) * csi0 * power * bsint;

    //----------------------------------------------------------------------      
    } else if (GetPerturbativeOrder() == 1) { // NLL

      double csi1(0.0),csihat1(0.0);
      
      if (GetFlavourSchemeInt() == 0 && GetWalpha() == 0) { // FFNS
	
	unsigned nff = GetRegionMU2(Q2);
	double C2  = GetC2(nff);
	double lambda1 = 3./8. - Pi2/2. + 6.*zeta3 - C2/18. * (3. + 4.*Pi2);

	if (GetRenormalisationSchemeInt() == 1) { // REN FIXED
	  csi1 = eta0 * (1.0 - 5./9.*aref/M_PI * C2);
	  csihat1 = eta0 * (lambda0 + 1.0/4.0*aref/M_PI * lambda1);
	  
	} else if (GetRenormalisationSchemeInt() >= 2) { // REN ALPMZ or ALGMU or AFAKE
	  double b0 = Getb0(nff);
	  double arefsq = aref*aref;

	  ////////////////////////////////////////////////////////////////////////////////
	  // FRAGILE! THE VALUE IS ASSOCIATED TO THE VALUE OF ALPHAGMU
	  double DELTAGMU = GetDELTAGMU();
	  double DMUfac = (GetRenormalisationSchemeInt() == 2) ? 5.0/9.0/M_PI*C2 : DELTAGMU;
	  double C2fac = (GetRenormalisationSchemeInt() == 2) ? 5.0/3.0*C2 : 3.0*M_PI*DELTAGMU;
	  ////////////////////////////////////////////////////////////////////////////////	  

	  // THE PRESENCE OF THE W IS DELICATE, AND EVEN IN THE FFNS IT IS ACTIVE ONLY ABOVE THE W MASS
	  // HENCE, IN CASE OF WON, WE USE THE EQUATION WITH THRESHOLDS
	  
	  // double wboson = (1.0/6.0 + 7.0/4.0*log(MZ2/MW2))*GetWalpha();	  
	  // // If wboson is null, then in csi1 there is a cancellation
	  // csi1 = ( eta0 + arefsq/2.0/M_PI * b0 * (log2mumz - log2mu0mz)
	  // 	   - aref/4.0/M_PI * 4.0*wboson * eta0 );
	  
	  csi1 = ( eta0 + arefsq/2.0/M_PI * b0 * (log2mumz - log2mu0mz)
		   + aref * ( DMUfac - 5.0/9.0/M_PI*C2 ) * eta0 );
	  
	  // csihat1 = ( lambda0*eta0 + 3.0/8.0/M_PI * arefsq * b0 * (log2mumz - log2mu0mz)
	  // 	      + aref/4.0/M_PI * ( lambda1 + C2fac - 3.0*wboson ) * eta0 );
	  csihat1 = ( lambda0*eta0 + 3.0/8.0/M_PI * arefsq * b0 * (log2mumz - log2mu0mz)
		      + aref/4.0/M_PI * ( lambda1 + C2fac ) * eta0 );
	  
	  // std::cout << "C++" << 2.0/3.0 * (5.0/3.0 * C2 + wboson) << std::endl;
	  // std::cout << "log " << eta0/aref*M_PI << std::endl;
	  // std::cout << "A "
	  // 	    << b0 * (log2mumz - log2mu07mz) << std::endl;
	  // std::cout << "csi1 " << csi1 << " csihat1 " << csihat1 << std::endl;	  
	    
 	} else if (GetRenormalisationSchemeInt() == 0) { // REN MSBAR

	  unsigned nfa = GetRegionMU2(Q2);
	  double b0 = Getb0(nfa);
	  double b1 = Getb1(nfa);	  
	  double tk = 1./(2.*b0*M_PI) * log(alphaMU/alphaM0);     
	  double alphafact = alphaMU / (4.*Pi2*b0);
	  double exptfact = 1 - exp(-2*M_PI*b0*tk);
	  csi1 = 2. * tk - alphafact * exptfact * (20./9. * C2 + 4.*M_PI*b1/b0);
	  csihat1 = 3./2. * tk + alphafact * exptfact * (lambda1 - 3.*M_PI*b1/b0);
	  
	} else {
	  throw std::runtime_error("[Asymptotic]: ERROR");
	} // END OF REN SCHEME

      } else if (GetFlavourSchemeInt() == 1 || (GetFlavourSchemeInt() == 0 && GetWalpha() == 1) ) { // VFNS
	
	unsigned k = GetRegionMU2(Q2);
	if (k >= q2thrs.size())
	  throw std::runtime_error("[Asymptotic]: ERROR");

	if (GetRenormalisationSchemeInt() == 1) { // REN FIXED
	  const double logq2mk = log(Q2/q2thrs[k]);
	  const double a2pisq = pow(aref/(2.0*M_PI),2);
	  csi1 = eta0 - 20./9.* a2pisq * GetC2(k) * logq2mk;
	  double lambda1 = 3./8. - Pi2/2. + 6.*zeta3 - GetC2(k)/18. * (3. + 4.*Pi2);    	  
	  csihat1 = eta0 * lambda0 + a2pisq * lambda1 * logq2mk;
	  for (unsigned i=0; i<k; i++) {
	    const double logmip1mi = log(q2thrs[i+1]/q2thrs[i]);
	    csi1 += - 20./9. * a2pisq * GetC2(i) * logmip1mi;
	    double lambda1 = 3./8. - Pi2/2. + 6.*zeta3 - GetC2(i)/18. * (3. + 4.*Pi2);    	 	    
	    csihat1 += a2pisq * lambda1 * logmip1mi;
	  }
	  
	} else if (GetRenormalisationSchemeInt() >= 2) { // REN ALPMZ or ALGMU
	  
	  const double m2k = q2thrs[k];
	  const double m2kp1 = (k >= 9) ? MZ2 : q2thrs[k+1];
	  // const double m2kp1 = q2thrs[k+1];
	  // std::cout << "k " << k << " m2kp1 " << m2kp1 << std::endl;
	  const double logq2m2k = log(Q2/m2k);
	  const double log2q2m2kp1 = pow(log(Q2/m2kp1),2);
	  const double log2m2km2kpi = pow(log(m2k/m2kp1),2);
	  const double difflog2 = log2q2m2kp1 - log2m2km2kpi;
	  double sumlog(0.0), sumb0log2(0.0), sumCDlog(0.0), sumlamlog(0.0), sumDlog(0.0);
	  for (unsigned i=0; i<k; i++) {
	    // std::cout << "i " << i << " mip1 " << q2thrs[i+1] << std::endl;
	    const double lambda1i = 3./8. - Pi2/2. + 6.*zeta3 - GetC2(i)/18. * (3. + 4.*Pi2);
	    const double logmip1mi = log(q2thrs[i+1]/q2thrs[i]);
	    const double log2mimip1 = pow(log(q2thrs[i]/q2thrs[i+1]),2);
	    sumlog += logmip1mi;
	    sumb0log2 += Getb0(i) * log2mimip1;
	    sumCDlog += (10.0/9.0*GetC2(i) - GetDk(i)) * logmip1mi;
	    sumlamlog += lambda1i * logmip1mi;
	    sumDlog += GetDk(i) * logmip1mi;
	    // std::cout << "Di(" << i << ") " << GetDk(i) << std::endl;
	    // std::cout << "Di expr "
	    // 	      << ( 2.0*M_PI*Getb0(i) * log(q2thrs[i+1]/MZ2)
	    // 		   + 10.0/9.0*GetC2(i) ) << std::endl;	    
	    //  sumClog += GetC2(i) * logmip1mi;	    
	  }
	  csi1 = ( aref/M_PI * ( logq2m2k + sumlog )
		   + aref*aref/2.0/M_PI * 
   		   ( Getb0(k) * difflog2 - sumb0log2
		     - 1.0/M_PI * (10.0/9.0*GetC2(k) - GetDk(k))* logq2m2k
		     - 1.0/M_PI * sumCDlog ) );
	  const double lambda1k = 3./8. - Pi2/2. + 6.*zeta3 - GetC2(k)/18. * (3. + 4.*Pi2);
	  csihat1 = ( lambda0*aref/M_PI * ( logq2m2k + sumlog )
		      + 3.0/8.0*aref*aref/M_PI *
		      ( Getb0(k) * difflog2 - sumb0log2 )
		      + aref*aref/4.0/M_PI/M_PI *
		      ( lambda1k * logq2m2k + sumlamlog
	  		+ 3.0/2.0 * ( GetDk(k) * logq2m2k + sumDlog ) ) );
	  
	  // std::cout << "logq2m2k + sumlog " << logq2m2k + sumlog << std::endl;
	  // std::cout << "A " 
	  // 	    << ( Getb0(k) * difflog2 - sumb0log2
	  // 	     - 1.0/M_PI * (10.0/9.0*GetC2(k) - GetDk(k))* logq2m2k
	  // 		 - 1.0/M_PI * sumCDlog ) << std::endl;
	  // std::cout << "csi1 " << csi1 << " csihat1 " << csihat1 << std::endl;
	  
	} else if (GetRenormalisationSchemeInt() == 0) { // REN MSBAR
	  double alphaMK = myathrs[k]; // aQED(q2thrs[k]);
	  double tk = 1./(2.*Getb0(k)*M_PI) * log(alphaMU/alphaMK);        
	  double alphafact = alphaMU / (4.*Pi2*Getb0(k));
	  double exptfact = 1 - exp(-2*M_PI*Getb0(k)*tk);
	  csi1 = 2. * tk - alphafact * exptfact * (20./9. * GetC2(k) + 4.*M_PI*Getb1(k)/Getb0(k));
	  double lambda1 = 3./8. - Pi2/2. + 6.*zeta3 - GetC2(k)/18. * (3. + 4.*Pi2);    
	  csihat1 = 3./2. * tk + alphafact * exptfact * (lambda1 - 3.*M_PI*Getb1(k)/Getb0(k));    
	  for (unsigned i=0; i<k; i++) {
	    double alphaMIP1 = myathrs[i+1]; // aQED(q2thrs[i+1]);	  
	    double alphaMI = myathrs[i]; // aQED(q2thrs[i]);
	    double barti = 1./(2.*Getb0(i)*M_PI) * log(alphaMIP1/alphaMI);
	    double alphafact = alphaMIP1 / (4.*Pi2*Getb0(i));    
	    double exptfact = 1 - exp(-2*M_PI*Getb0(i)*barti);
	    double factcsi1 = 20./9. * GetC2(i) + 4.*M_PI*Getb1(i)/Getb0(i);
	    csi1 += 2. * barti - alphafact * exptfact * factcsi1;
	    double lambda1 = 3./8. - Pi2/2. + 6.*zeta3 - GetC2(i)/18. * (3. + 4.*Pi2);    
	    csihat1 += 3./2. * barti + alphafact * exptfact * (lambda1 - 3.*M_PI*Getb1(i)/Getb0(i));
	  }

	} else {
	  throw std::runtime_error("[Asymptotic]: ERROR");
	} // END OF REN SCHEME
	
      } else {
	throw std::runtime_error("[Asymptotic]: ERROR");
      } // END OF FLAV SCHEME
      
      const double factor = exp( csihat1 - emc * csi1 ) / tgamma( 1 + csi1 ) * csi1;
      const double power = pow(omz, csi1 - betab);

      const double a0opi = (aemfix) ? aref/M_PI : alphaM0/M_PI;
      const double aopi  = alphaMU/M_PI;
      
      const double lz = log(omz);
      
      if (GetFactorisationSchemeInt() == 0) { // MSBAR
	
	////////////////////////////////////////////////////////////////////////////////
	/// BREMSSTRAHLUNG 		
	const double psicsi1 = (bson) ? psi(csi1+kappa) : psi(csi1);
	const double dpsicsi1 = (bson) ? dpsi(csi1+kappa,1) : dpsi(csi1,1);
	////////////////////////////////////////////////////////////////////////////////	
	const double A = - emc - psicsi1;
	const double B = ( 0.5 * emc*emc + Pi2/12. + emc*psicsi1
			   + 0.5 * psicsi1*psicsi1 - 0.5 * dpsicsi1 );
	const double clog0 = 1. + a0opi * ( (logmu0 - 1)*( A + 3./4. ) - 2*B + 7./4. );
	const double clog1 = a0opi * ( (logmu0 - 1 - 2*A) );
	const double clog2 = a0opi * ( - 1.0 );
	
	////////////////////////////////////////////////////////////////////////////////
	/// BREMSSTRAHLUNG 		
	double normbs = (bson) ? tgamma(kappa)*tgamma(csi1)/tgamma(kappa+csi1) : 1.0;
	////////////////////////////////////////////////////////////////////////////////
	return factor * power * normbs * ( clog0 + clog1 * lz + clog2 * lz * lz );
	
      } else if (GetFactorisationSchemeInt() == 1) { // DELTA

	if (aemfix) {
	  
	  ////////////////////////////////////////////////////////////////////////////////
	  /// BREMSSTRAHLUNG 		
	  const double psicsi1 = (bson) ? psi(csi1+kappa) : psi(csi1);
	  ////////////////////////////////////////////////////////////////////////////////	
	  const double A = - emc - psicsi1;
	  ////////////////////////////////////////////////////////////////////////////////
	  /// BREMSSTRAHLUNG	
	  double normBS = (bson) ? pow(omz, kappa) * tgamma(kappa)*tgamma(csi1)/tgamma(kappa+csi1) : 1.0;
	  ////////////////////////////////////////////////////////////////////////////////       
	  return factor * power * normBS * ( 1.0 + aref/M_PI * logmu0 * ( A + lz + 3.0/4.0 ) );

	} else {
	
	  const double denz = 1.0+a0opi*lz*lz;
	  const double lz2 = lz*lz;
	  const double lz4 = lz2*lz2;
	  const double lz6 = lz4*lz2;
	  const double a0opisq = a0opi*a0opi;
	  const double a0opicb = a0opisq*a0opi;
	  const double cc = M_PI*M_PI/6.0 - 1.0;

	  ////////////////////////////////////////////////////////////////////////////////
	  /// BREMSSTRAHLUNG 		    	  
	  const double psi0k = (bson) ?  psi(csi1+kappa)    :  psi(csi1)   ;
	  const double psi1k = (bson) ? dpsi(csi1+kappa, 1) : dpsi(csi1, 1);
	  const double psi2k = (bson) ? dpsi(csi1+kappa, 2) : dpsi(csi1, 2);
	  const double psi3k = (bson) ? dpsi(csi1+kappa, 3) : dpsi(csi1, 3);
	  const double psi4k = (bson) ? dpsi(csi1+kappa, 4) : dpsi(csi1, 4);
	  ////////////////////////////////////////////////////////////////////////////////	  

	  const double Ak = -emc-psi0k;
	  const double Bk = 0.5*emc*emc + Pi2/12.0 + emc*psi0k + 0.5*psi0k*psi0k - 0.5*psi1k;
	  const double Ck = ( -emc*emc*emc/6.0-emc*Pi2/12.0-1.0/6.0*psi0k*(3.0*emc*emc+3.0*emc*psi0k+psi0k*psi0k)
			      +0.5*psi0k*(psi1k-Pi2/6.0)+0.5*emc*psi1k-1.0/6.0*psi2k-zeta3/3.0 );
	  const double Dk = ( pow(emc,4)/24.0+emc*emc*Pi2/24.0+Pi2*Pi2/160.0
			      +1.0/3.0*emc*zeta3+1.0/24.0*(8.0*zeta3+
							   (2.0*emc+psi0k)*(pow(emc+psi0k,2)+emc*emc+Pi2))*psi0k
			      -1.0/4.0*(emc*emc+Pi2/6.0+(2*emc+psi0k)*psi0k-0.5*psi1k)*psi1k+
			      1.0/6.0*(emc+psi0k)*psi2k-1.0/24.0*psi3k );
	
	  const double d1 = - Ak;
	  const double d2 = 2.0*Bk - Pi2/6.0;
	  const double d3 = -6.0*Ck + 0.5*Pi2*Ak - 2.0*zeta3;
	  const double d4 = 24.0*Dk - Pi2*Ak*Ak + 8.0*zeta3*Ak + Pi2*psi1k - 3.0*Pi2*Pi2/20.0;
	  const double d5 = ( pow(emc+psi0k,5) + ( -10*pow(emc+psi0k,3) + 15*(emc+psi0k)*psi1k - 10*psi2k ) * psi1k
			      + 10*pow(emc+psi0k,2)*psi2k - 5*(emc+psi0k)*psi3k + psi4k );
	
	  const double s10h = -1.0/a0opi;
	  const double s11h = (1.0-2.0*d1)*lz;
	  const double s12h =  a0opi*(cc-1.0+3.0*d1-3.0*d2)*lz*lz + cc-d1+d2;
	  const double s13h =  ( pow(a0opi,2)*(1.0-4.0*d1+6.0*d2-4.0*d3-2*cc*(1-2*d1))*pow(lz,3)
				 +2*(a0opi)*(d1-3*d2+2*d3-cc*(1-2*d1))*lz );
	  const double s14h =  ( a0opi*(cc*cc-a0opi*cc*(3-2*cc)*lz2+a0opisq*(1-3*cc+cc*cc)*lz4)
				 - a0opi*(-2*cc+a0opi*(3+8*cc)*lz2-5*a0opisq*(1-2*cc)*lz4) * d1
				 - a0opi*(1+2*cc-a0opi*(13+8*cc)*lz2+10*a0opisq*(1-cc)*lz4) * d2
				 + 2*a0opi*(1-10*a0opi*lz2+5*a0opisq*lz4) * d3
				 - a0opi*(1-10*a0opi*lz2+5*a0opisq*lz4) * d4 );
	  const double s15h =  ( a0opisq*(3*cc*cc-2*a0opi*cc*(2-3*cc)*lz2+a0opisq*(1-4*cc+3*cc*cc)*lz4)*lz
				 - 2*a0opisq*(3*cc*(1+cc)-2*a0opi*(1+3*cc-3*cc*cc)*lz2
					      + 3*a0opisq*(1-3*cc+cc*cc)*lz4)*lz * d1
				 - a0opisq*(-3*(1+6*cc)+2*a0opi*(11+6*cc)*lz2-15*a0opisq*(1-2*cc)*lz4)*lz * d2
				 + 4*a0opisq*(-3*(1+cc)+2*a0opi*(6+cc)*lz2-5*a0opisq*(1-cc)*lz4)*lz * d3
				 + 5*a0opisq*(3-10*a0opi*lz2+3*a0opisq*lz4)*lz * d4
				 - 2*a0opisq*(3-10*a0opi*lz2+3*a0opisq*lz4)*lz * d5 );

	  const double s20h =  1.0/a0opi * lz;
	  const double s21h =  -(1-d1)*lz*lz - 1.0/a0opi*d1;
	  const double s22h =  -(a0opi)*(cc-1+2*d1-d2)*pow(lz,3)-(cc-2*d1+3*d2)*lz;
	  const double s23h =  ( -pow(a0opi,2)*(1-3*d1+3*d2-d3-cc*(2-3*d1))*pow(lz,4)
				 -(a0opi)*(3*d1-8*d2+6*d3-2*cc*(1-d1))*lz*lz+cc*d1-d2+d3 );
	  const double s24h =  ( a0opi*(cc*cc-a0opi*cc*(3-2*cc)*lz2+a0opisq*(1-3*cc+cc*cc)*lz4)*lz
				 + 4*a0opi*(-cc+a0opi*(1+cc)*lz2-a0opisq*(1-2*cc)*lz4)*lz * d1
				 - 3*a0opi*(-1-2*cc+5*a0opi*lz2-2*a0opisq*(1-cc)*lz4)*lz* d2
				 - 4*a0opi*(2-5*a0opi*lz2+a0opisq*lz4)*lz * d3
				 + a0opi*(5-10*a0opi*lz2+a0opisq*lz4)*lz * d4 );
	  const double s25h =  ( a0opisq*(3*cc*cc-2*a0opi*cc*(2-3*cc)*lz2+a0opisq*(1-4*cc+3*cc*cc)*lz4)*lz2
				 - a0opi*(cc*cc-3*a0opi*cc*(3+cc)*lz2+a0opisq*(5+6*cc-9*cc*cc)*lz4
					  - 5*a0opicb*(1-3*cc+cc*cc)*lz6) * d1
				 + 2*a0opi*(cc-3*a0opi*(1+4*cc)*lz2+3*a0opisq*(4-cc)*lz4-5*a0opicb*(1-2*cc)*lz6) * d2
				 + a0opi*(-(1+2*cc)+6*a0opi*(4+3*cc)*lz2
					  -5*a0opisq*(9-2*cc)*lz4+10*a0opicb*(1-cc)*lz6) * d3
				 + a0opi*(2-33*a0opi*lz2+40*a0opisq*lz4-5*a0opicb*lz6) * d4
				 - a0opi*(1-15*a0opi*lz2+15*a0opisq*lz4 - a0opicb*lz6) * d5 );
	
	  const double denz2 = pow(denz,2);
	  const double denz3 = pow(denz,3);
	  const double denz4 = pow(denz,4);
	  const double denz5 = pow(denz,5);
	  const double denz6 = pow(denz,6);
	
	  const double s10 = aopi/a0opi + (aopi-a0opi) * s10h/denz;
	  const double s11 = (aopi-a0opi) * s11h/denz2;
	  const double s12 = (aopi-a0opi) * s12h/denz3;
	  const double s13 = (aopi-a0opi) * s13h/denz4;
	  const double s14 = (aopi-a0opi) * s14h/denz5;
	  const double s15 = (aopi-a0opi) * s15h/denz6;

	  const double s20 = aopi/a0opi * (-lz+d1) + (aopi-a0opi) * s20h/denz;
	  const double s21 = (aopi-a0opi) * s21h/denz2;
	  const double s22 = (aopi-a0opi) * s22h/denz3;
	  const double s23 = (aopi-a0opi) * s23h/denz4;
	  const double s24 = (aopi-a0opi) * s24h/denz5;
	  const double s25 = (aopi-a0opi) * s25h/denz6;
	
	  const double s1sum = s10 + s11 + s12 + s13 + s14 + s15;
	  const double s2sum = s20 + s21 + s22 + s23 + s24 + s25;

	  const double fac_s1 = (1.0 + 3.0/4.0*a0opi * logmu0);
	  const double fac_s2 = - a0opi * logmu0;

	  ////////////////////////////////////////////////////////////////////////////////
	  /// BREMSSTRAHLUNG 		    
	  double bspow = (bson) ? pow(omz, kappa) * tgamma(kappa) * tgamma(csi1) / tgamma(kappa+csi1) : 1.0;
	  ////////////////////////////////////////////////////////////////////////////////    
	  return factor * power * (fac_s1 * s1sum + fac_s2 * s2sum) * bspow;
	  
	}
	  
      } else {
	throw std::runtime_error("[Asymptotic]: ERROR");	
      } // END OF FAC SCHEME
    } else {	
      throw std::runtime_error("[Asymptotic]: ERROR");
    } // END OF ORDER
    
    throw std::runtime_error("[Asymptotic]: ERROR");	    
  }

  //_________________________________________________________________________________
  double local_asy(double const& omz, double const& betaE, double const& betaS, double const& gamma)
  {
    return exp(3.0/4.0*betaS - emc*betaE) / tgamma(1+betaE) * betaE * pow(omz, betaE-gamma);
  }

  double local_rec(double const& z, double const& betaH)
  {
    double omz = 1.-z;
    double lnz = log(z);
    double lnomz = log(omz);    
    
    double rec1z = 1+z;
    double rec2z = 4*(1+z)*lnomz + 5.0 + z + 4.0*lnz/omz - 3.0*(1+z)*lnz;
    double rec3z = ( (1+z) * (6.0*dilog(z) + 12.0*pow(lnomz,2) - 3.0*Pi2)
		     + 1.0/omz * (3.0/2.0*(1.0+8*z+3*z*z)*lnz + 6*(z+5)*omz*lnomz
				  + 12*(1+z*z)*lnz*lnomz - 1.0/2.0*(1+7*z*z)*lnz*lnz
				  + 1.0/4.0*(39-24*z-15*z*z)) );
    
    return - betaH/2.0 * rec1z - pow(betaH,2)/8.0 * rec2z - pow(betaH,3)/48.0 * rec3z;
  }

  void warmup(int const& index, double const& Q, double & betaE, double & betaS, double & betaH)
  {
    double Q2 = Q*Q;
    double aem = aQED(Q2);
    const std::vector<double> q2thrs = GetThresholds2();    
    const double me2 = q2thrs[0];
    const double ln2 = log(Q2/me2);
    const double eta = aem/M_PI * ln2;
    const double beta = aem/M_PI * (ln2-1.0);
    
    if (index == 0) { // nobeta
      betaE = eta;
      betaS = eta;
      betaH = eta;
    } else if (index == 1) { // mixed
      betaE = beta;
      betaS = eta;
      betaH = eta;
    } else if (index == 2) { // eta
      betaE = beta;
      betaS = beta;
      betaH = eta;
    } else if (index == 3) { // beta
      betaE = beta;
      betaS = beta;
      betaH = beta;
    } else {
      throw std::runtime_error("[ElPDFLegacy]: invalid index");      
    }    
  }
  
  double LLPDF(int const& index, double const& z, double const& omz,
	       double const& Q, double const& gamma)
  {
    double betaE, betaS, betaH;
    warmup(index, Q, betaE, betaS, betaH);
    double asy = local_asy(omz, betaE, betaS, gamma);
    double rec = local_rec(z, betaH);
    return asy + rec * pow(omz,1.-gamma);
  }
  
  double AsyLLPDF(int const& index, double const& omz,
		  double const& Q, double const& gamma)
  {
    double betaE, betaS, betaH;
    warmup(index, Q, betaE, betaS, betaH);    
    double asy = local_asy(omz, betaE, betaS, gamma);
    return asy;
  }  

}
