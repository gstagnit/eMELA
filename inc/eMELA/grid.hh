#pragma once

#ifdef USE_LHAPDF
#include <LHAPDF/GridPDF.h>
#endif

namespace eMELA {
  
  class GridMELA {
    
  public:
    
    GridMELA();    
    GridMELA(std::string const& pdfname);    
    GridMELA(int const& lhaid);    
    GridMELA(const GridMELA& copy);    
    ~GridMELA();
#ifdef USE_LHAPDF
    LHAPDF::PDFInfo get_info() { return _lha->info(); }
#endif
    double NumGrid(const int & idx, const double & x, const double & Q, const double & beta);

  private:
    // Private constructor
    void _GridMELA(LHAPDF::PDFInfo const& info);
#ifdef USE_LHAPDF
    LHAPDF::GridPDF* _lha;
#endif
    double _omxmax;    
  };

}

