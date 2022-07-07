///
/// Minimal example to get the PDF from the code
///
#include <iostream>

#include "eMELA/eMELA.hh"

using namespace std;
using namespace eMELA;

int main()
{
  double aMZ = 1.0/127.95471413988483012;
  QuickInitialize("NLL", "DELTA", "MSBAR", aMZ);

  int pdgid = 11;
  double x = 0.5;
  double omx = 1.0-x;
  double Q = 100.0;
  double gamma = 1.0;
  
  double pdf = CodePdf(pdgid, x, omx, Q, gamma);
  cout << "Electron PDF at (x=" << x << ",Q=" << Q << ") = " << pdf << endl;

  return 0;
}
