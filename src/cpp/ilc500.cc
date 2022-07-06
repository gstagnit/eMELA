#ifdef USE_LHAPDF

#include "eMELA/ilc500.hh"
#include "eMELA/eMELA.hh"
#include <cmath>
#include <boost/math/quadrature/gauss_kronrod.hpp>

using std::exp;
using std::pow;
using std::sqrt;

ilc500::ilc500()
{
}

void ilc500::init()
{
  auto info = _vptr_lha[0]->info();  
  f11 = boost::lexical_cast<double>(info.get_entry("ILC500_f11"));
  f01 = boost::lexical_cast<double>(info.get_entry("ILC500_f01"));
  p01 = boost::lexical_cast<double>(info.get_entry("ILC500_p01"));
  q01 = boost::lexical_cast<double>(info.get_entry("ILC500_q01"));
  f10 = boost::lexical_cast<double>(info.get_entry("ILC500_f10"));
  p10 = boost::lexical_cast<double>(info.get_entry("ILC500_p10"));
  q10 = boost::lexical_cast<double>(info.get_entry("ILC500_q10"));
  f001 = boost::lexical_cast<double>(info.get_entry("ILC500_f001"));
  p001 = boost::lexical_cast<double>(info.get_entry("ILC500_p001"));
  q001 = boost::lexical_cast<double>(info.get_entry("ILC500_q001"));
  f002 = boost::lexical_cast<double>(info.get_entry("ILC500_f002"));
  p002 = boost::lexical_cast<double>(info.get_entry("ILC500_p002"));
  q002 = boost::lexical_cast<double>(info.get_entry("ILC500_q002"));
}

ilc500::ilc500(const std::string& name) : bspdf(name)
{ 
}

ilc500::ilc500(int const& lhaid) : bspdf(lhaid)
{
}

std::string ilc500::beamspectrum_info()
{
  std::ostringstream oss;
  oss.precision(17);
  oss << "beamspectrum_type: " << "ILC500" << "\n"
      << "ILC500_f11: " << f11 << "\n"
      << "ILC500_f01: " << f01 << "\n"
      << "ILC500_p01: " << p01 << "\n"
      << "ILC500_q01: " << q01 << "\n"
      << "ILC500_f10: " << f10 << "\n"
      << "ILC500_p10: " << p10 << "\n"
      << "ILC500_q10: " << q10 << "\n"
      << "ILC500_f001: " << f001 << "\n"
      << "ILC500_p001: " << p001 << "\n"
      << "ILC500_q001: " << q001 << "\n"
      << "ILC500_f002: " << f002 << "\n"
      << "ILC500_p002: " << p002 << "\n"
      << "ILC500_q002: " << q002 << "\n";
  return oss.str();
}

double ilc500::calc_pdf(int icom, int particleid, int beamid, double x, double omx, double Q, double beta)
{
  int prod;
  if      ((particleid == 11 && beamid == +1) || (particleid == -11 && beamid == -1)) prod = +11;
  else if ((particleid == 11 && beamid == -1) || (particleid == -11 && beamid == +1)) prod = -11;
  else if ((particleid == 22 && beamid == +1) || (particleid ==  22 && beamid == -1)) prod =  22;  
  else throw std::runtime_error("Unknown particle/beam combination, particle = " + std::to_string(particleid)+", beam = " + std::to_string(beamid));

  //NOTE: only the product of both beams are meaningful
  //and hence the normalisation of each beam is customary
  switch (icom)
    {
    case 1:
      return sqrt(f11) * eMELA::GridPdf(prod, x, omx, Q, beta);
      break;
    case 2:
      if (beamid == 1)
	{
	  return f01 / sqrt(f11) * integrate_isr(prod, p01, q01, x, omx, Q, beta);
	}
      else
	{
	  return sqrt(f11) * eMELA::GridPdf(prod, x, omx, Q, beta);
	}
      break;
    case 3:
      if (beamid == -1)
	{
	  return f10 / sqrt(f11) * integrate_isr(prod, p10, q10, x, omx, Q, beta);
	}
      else
	{
	  return sqrt(f11) * eMELA::GridPdf(prod, x, omx, Q, beta);
	}
      break;
    case 4:
      if (beamid == 1)
	{
	  return f001 * integrate_isr(prod, p001, q001, x, omx, Q, beta);
	}
      else
	{
	  return f002 * integrate_isr(prod, p002, q002, x, omx, Q, beta);
	}
      break;
    default:
      throw std::runtime_error("Unimplemented cases: icom=" + std::to_string(icom));
      break;
    }
}

double ilc500::asy_pdf(int icom, int particleid, int beamid, double x, double omx, double Q, double beta)
{
  int prod;
  if      ((particleid == 11 && beamid == +1) || (particleid == -11 && beamid == -1)) prod = +11;
  else if ((particleid == 11 && beamid == -1) || (particleid == -11 && beamid == +1)) prod = -11;
  else if ((particleid == 22 && beamid == +1) || (particleid ==  22 && beamid == -1)) prod =  22;  
  else throw std::runtime_error("Unknown particle/beam combination, particle = " + std::to_string(particleid)+", beam = " + std::to_string(beamid));
  
  //NOTE: only the product of both beams are meaningful
  //and hence the normalisation of each beam is customary
  switch (icom)
    {
    case 1:
      return sqrt(f11) * eMELA::GridPdf(prod, x, omx, Q, beta);
      break;
    case 2:
      if (beamid == 1)
	{
	  return f01 / sqrt(f11) * asy_integrate_isr(prod, p01, q01, x, omx, Q, beta);
	}
      else
	{
	  return sqrt(f11) * eMELA::GridPdf(prod, x, omx, Q, beta);
	}
      break;
    case 3:
      if (beamid == -1)
	{
	  return f10 / sqrt(f11) * asy_integrate_isr(prod, p10, q10, x, omx, Q, beta);
	}
      else
	{
	  return sqrt(f11) * eMELA::GridPdf(prod, x, omx, Q, beta);
	}
      break;
    case 4:
      if (beamid == 1)
	{
	  return f001 * asy_integrate_isr(prod, p001, q001, x, omx, Q, beta);
	}
      else
	{
	  return f002 * asy_integrate_isr(prod, p002, q002, x, omx, Q, beta);
	}
      break;
    default:
      throw std::runtime_error("Unimplemented cases: icom=" + std::to_string(icom));
      break;
    }
}

double ilc500::integrate_isr(int pid, double p, double q, double x, double omx, double Q, double beta)
{
  //auxiallary paramters that controlls the variable transformation(gamma)
  //and domain splitting(mid)
  //The result should be independent on them
  double gamma;
  double mid;
  if (pid == 11)
    {
      gamma = 0.03;
      mid = 0.5;
    }
  else
    {
      gamma = 0.0;
      mid = 1.0;
    }
  //Currently set e+ in e- beam to be 0
  //due to some instability when x->1
  if (pid == -11)
    {
      return 0;
    }
  auto term1 = [&](double r)
    {
      double ra = mid;
      double t = ra * pow(r, 1.0 / kappa);
      double z = t * omx + x;
      double omz = (1.0 - t) * omx;
      double g = exp(p * t * omx / z) * exp(q * pow((t * omx / z), 1.5));
      g *= pow(z, -kappa);
      g *= pow(1.0 - t, -1. + gamma) * pow(ra, kappa) / kappa;
      g *= pow(omx, kappa + gamma - beta);
      g *= eMELA::GridPdf(pid, z, omz, Q, gamma);
      //g *= isrll(z, omz, Q, gamma);
      //std::cout<<std::setprecision(17)<<pid<<" "<<r<<" "<<t<<" "<<x<<" "<<z<<" "<<omz<<" "<<g<<std::endl;
      return g;
    };
  auto term2 = [&](double r)
    {
      double rb = 1. - mid;
      double t = 1. - rb * pow(r, 1. / gamma);
      double z = t * omx + x;
      double omz = rb * pow(r, 1. / gamma) * omx;
      double g = exp(p * t * omx / z) * exp(q * pow((t * omx / z), 1.5));
      g *= pow(z, -kappa);
      g *= pow(t, -1. + kappa) * pow(rb, gamma) / gamma;
      g *= pow(omx, kappa + gamma - beta);
      //g *= isrll(z, omz, Q, gamma);
      //std::cerr<<r<<" "<<t<<" "<<z<<" "<<omz<<" "<<g<<std::endl;
      g *= eMELA::GridPdf(pid, z, omz, Q, gamma);
      //std::cerr<<g<<std::endl;
      return g;
    };
  boost::math::quadrature::gauss_kronrod<double, 31> integrator;
  double err1 = 0;
  double err2 = 0;
  double int1 = 0;
  double int2 = 0;
  int1 =  integrator.integrate(term1, 0, 1., 15, 1e-6, &err1);
  if (pid == 11)
    { int2 = integrator.integrate(term2, 0, 1., 15, 1e-6, &err2); }
  //std::cerr << "x:" << x << " omx:" << omx << " int:" << int1 << "+-" << err1 << " " << int2 << "+-" << err2 << " => " << int1 + int2 << "+-" << err1 + err2 << std::endl;
  return int1 + int2;
}

double ilc500::asy_integrate_isr(int pid, double p, double q, double x, double omx, double Q, double beta)
{
  //in asymptotic region, only electron is nonzero
  if (pid != 11) { return 0; }
  return eMELA::AsyInt(omx, Q, beta, kappa);
}

double ilc500::bs_integral(int icom, int beamid)
{
  auto pqint = [](double p, double q)
    {
      boost::math::quadrature::gauss_kronrod<double, 31> integrator;
      double res = integrator.integrate([&](double t)
    {
      double b = 1. / 3.;
      double omx = pow(t, 1. / b);
      double r = exp(p * omx) * exp(q * pow(omx, 3. / 2.));
      return r / b;
    },
					0., 1., 25, 1e-8);
      return res;
    };
  switch (icom)
    {
    case 1:
      return sqrt(f11);
      break;
    case 2:
      if (beamid == 1)
	{
	  return f01 / sqrt(f11) * pqint(p01, q01);
	}
      else
	{
	  return sqrt(f11);
	}
      break;
    case 3:
      if (beamid == -1)
	{
	  return f10 / sqrt(f11) * pqint(p10, q10);
	}
      else
	{
	  return sqrt(f11);
	}
      break;
    case 4:
      if (beamid == 1)
	{
	  return f001 * pqint(p001, q001);
	}
      else
	{
	  return f002 * pqint(p002, q002);
	}
      break;
    default:
      throw std::runtime_error("Unimplemented cases: icom=" + std::to_string(icom));
      break;
    }
}

#endif
