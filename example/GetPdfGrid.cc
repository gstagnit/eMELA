///
/// Minimal example to get the PDF from the grid
///
/// Remember to export the path to the grid folder:
///    export LHAPDF_DATA_PATH=/path/to/grid/folder:$LHAPDF_DATA_PATH
///
#include <iostream>

#include "eMELA/eMELA.hh"

using namespace std;
using namespace eMELA;

int main()
{
  InitializeFromGrid("NLL_DELTA_MSBAR");
  
  int pdgid = 11;
  double x = 0.5;
  double omx = 1.0-x;
  double Q = 100.0;
  double gamma = 1.0;
  
  double pdf = GridPdf(pdgid, x, omx, Q, gamma);
  cout << "Electron PDF at (x=" << x << ",Q=" << Q << ") = " << pdf << endl;

  return 0;
}
