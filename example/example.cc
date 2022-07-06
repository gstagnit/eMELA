///
/// Minimal example
///
/// Compile with:
///    g++ `path/to/eMELA/install/bin/eMELA-config --cppflags --ldflags` -o example example.cc
/// You may need to export:
///    export DYLD_LIBRARY_PATH=path/to/eMELA/install/lib:$DYLD_LIBRARY_PATH
///
#include <iostream>

#include "eMELA/eMELA.hh"

using namespace std;
using namespace eMELA;

int main()
{
  SetDefaultParameters();
  InitializeEvolution();

  double pdf = CodePdf(1, 0.5, 0.5, 100.0, 1.0);
  cout << pdf << endl;

  return 0;
}
