///
/// Minimal example to show how to write a grid
///
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>

#include "eMELA/eMELA.hh"

using namespace std;

int main (int argc, char ** argv)
{
  if (argc != 4) {
    cout << "usage: ./CreateGrid PTO FAC REN" << endl;    
    return -1;
  }

  string pto(argv[1]);
  string fac(argv[2]);
  string ren(argv[3]);
  
  double ainit;
  if (ren == "MSBAR")      ainit = 1.0/127.95471413988483012;
  else if (ren == "ALPMZ") ainit = 1.0/128.94002842670017122;
  else if (ren == "ALGMU") ainit = 1.0/132.18289853516;
  else throw runtime_error("wrong ren scheme");
    
  eMELA::QuickInitialize(pto, fac, ren, ainit);

  std::ostringstream sgridname;
  sgridname << pto;
  if (pto == "NLL") sgridname << "_" << fac;
  sgridname << "_" << ren;
  std::string gridname = sgridname.str(); 

  eMELA::WriteGrid(gridname);
  
  return 0;
  
}
