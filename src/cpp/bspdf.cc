#ifdef USE_LHAPDF

#include <sys/stat.h>
#include "eMELA/bspdf.hh"
#include "eMELA/eMELA.hh"
#include "eMELA/ilc500.hh"
#include <fstream>
#include <sstream>
#include <set>

namespace
{
  std::string get_str(double val, int prec)
  {
    std::stringstream stream;
    stream << std::fixed << std::setprecision(prec) << val;
    return stream.str();
  }
}

bspdf::bspdf()
{
}

bspdf::bspdf(const std::string& name)
{
  LHAPDF::mkPDFs(name, _vptr_lha);  
}

bspdf::bspdf(int const& lhaid)
{
  // the mkPDFs function accepts only the set name as argument
  // ---> finding the name of the PDF set via lookupPDF
  std::pair<std::string, int> set_id = LHAPDF::lookupPDF(lhaid);
  std::string name = set_id.first;
  LHAPDF::mkPDFs(name, _vptr_lha);
}

bspdf::~bspdf()
{
}

double bspdf::get_pdf(int icom, int pid, int beamid, double x, double omx, double Q, double beta)
{
  if (omx < _omxcut)//using asymptotic
    {
      return asy_pdf(icom, pid, beamid, x, omx, Q, beta);
    }
  else if (_use_grid)//using precomputed LHAGrid
    {
      return grid_pdf(icom, pid, beamid, x, omx, Q, beta);
    }
  else
    {
      return calc_pdf(icom, pid, beamid, x, omx, Q, beta);
    }
}

double bspdf::grid_pdf(int icom, int pid, int beamid, double x, double omx, double Q, double beta)
{
  int index = (beamid + 1) / 2 * this->components() + icom;
  double res;
  switch (pid)
    {
    case 11:
      res = _vptr_lha[index]->xfxQ(x, Q).at(11);
      break;
    case 22:
      res = _vptr_lha[index]->xfxQ(x, Q).at(22);
      break;
    case -11:
      res = _vptr_lha[index]->xfxQ(x, Q).at(-11);
      break;      
    default:
      throw std::runtime_error("Unimplemented pid=" + std::to_string(pid));
      break;
    }
  if ((pid == 11 && beamid == 1) || (pid == -11 && beamid == -1))
    { res *= pow(omx, -beta); }
  else
    {
      res *= pow(omx, 1. - beta);
    }
  return res;
}

double bspdf::calc_pdf(int icom, int pid, int beamid, double x, double omx, double Q, double beta)
{
  throw std::runtime_error("Error: calc_pdf of base class shouldn't be called");
  return 0.;
}

double bspdf::asy_pdf(int icom, int pid, int beamid, double x, double omx, double Q, double beta)
{
  throw std::runtime_error("Error: calc_pdf of base class shouldn't be called");
  return 0.;
}

double bspdf::bs_integral(int icom, int beamid)
{
  throw std::runtime_error("Error: bs_integral of base class shouldn't be called");
  return 0.;
}

template <typename T>
std::shared_ptr<bspdf> internal_make_bspdf(T var, LHAPDF::PDFInfo info)
{
  std::shared_ptr<bspdf> ptr;

  if (!(info.has_key("beamspectrum_type"))) {
    throw std::runtime_error("Grid with ISR, but no bremsstrahlung!");
  }
  
  std::string bstype = info.get_entry("beamspectrum_type");
  if (bstype == "ILC500") 
    ptr = std::make_shared<ilc500>(var);
  else
    throw std::runtime_error("Unknown beamspectrum type:" + bstype);

  eMELA::InitializeFromGrid(var);
  
  return ptr;
}

std::shared_ptr<bspdf> make_bspdf(const std::string& name)
{
  std::unique_ptr<LHAPDF::GridPDF> lha(new LHAPDF::GridPDF(name, 0));
  auto info = lha->info();    
  return internal_make_bspdf(name, info);
}

 std::shared_ptr<bspdf> make_bspdf(int const& lhaid)
{
  std::unique_ptr<LHAPDF::GridPDF> lha(new LHAPDF::GridPDF(lhaid));
  auto info = lha->info();      
  return internal_make_bspdf(lhaid, info);
}

#endif
