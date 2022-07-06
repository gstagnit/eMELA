#pragma once

#include <memory>
#ifdef USE_LHAPDF
#include <LHAPDF/GridPDF.h>
#endif

//The main class for the PDF with beamstrahlung
class bspdf
{
public:
  bspdf();
  bspdf(const std::string& name);
  bspdf(int const& lhaid);  
  
  virtual ~bspdf();

#ifdef USE_LHAPDF  
  //return number of components
  virtual int components()
  {
    return _vptr_lha.size() / 2;
  }
#endif 

  //get the PDF
  //icom: the index of component
  //particleid: particle id(11:electron,22:photon,-11:positron)
  //beamid: beam id(1:electron,-1:positron)
  //x: x
  //omx: 1-x
  //Q: the factorisation scale for ISR
  //beta: an auxialiary parameter
  //return value: the PDF result times (1-x)^{1-beta}
  double get_pdf(int icom, int particleid, int beamid, double x, double omx, double Q, double beta = 1);

  //get the PDF from precomputed LHAgrid
  //icom: the index of component
  //particleid: particle id(11:electron,22:photon,-11:positron)
  //beamid: beam id(1:electron,-1:positron)
  //x: x
  //omx: 1-x
  //Q: the factorisation scale for ISR
  //beta: an auxialiary parameter
  //return value: the PDF result times (1-x)^{1-beta}
  double grid_pdf(int icom, int particleid, int beamid, double x, double omx, double Q, double beta = 1);

  //perform the co-evolution of beamstrahlung function and ISR
  //icom: the index of component
  //particleid: particle id(11:electron,22:photon,-11:positron)
  //beamid: beam id(1:electron,-1:positron)
  //x: x
  //omx: 1-x
  //Q: the factorisation scale for ISR
  //beta: an auxialiary parameter
  //return value: the PDF result times (1-x)^{1-beta}
  virtual double calc_pdf(int icom, int particleid, int beamid, double x, double omx, double Q, double beta = 1);

  //compute the co-evolution of beamstrahlung function and ISR when x->1
  //asymptotic analytic expressions are adopted
  //icom: the index of component
  //particleid: particle id(11:electron,22:photon,-11:positron)
  //beamid: beam id(1:electron,-1:positron)
  //x: x
  //omx: 1-x
  //Q: the factorisation scale for ISR
  //beta: an auxialiary parameter
  //return value: the PDF result times (1-x)^{1-beta}
  virtual double asy_pdf(int icom, int particleid, int beamid, double x, double omx, double Q, double beta = 1);

  //Debug only:
  //return the integral of one compoenent without coevolution with ISR
  virtual double bs_integral(int icom, int beamid);

  virtual std::string beamspectrum_info() = 0;

  virtual std::string bsname() = 0;

  void use_grid(bool use) { _use_grid = use; }
  
protected:

#ifdef USE_LHAPDF
  LHAPDF::PDFInfo _lhainfo;
  std::vector<std::shared_ptr<LHAPDF::PDF>> _vptr_lha;
#endif
  
  double _omxcut = 1e-8;
  double _use_grid = true;
};

std::shared_ptr<bspdf> make_bspdf(const std::string& pdfname);
std::shared_ptr<bspdf> make_bspdf(int const& lhaid);
