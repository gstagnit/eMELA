#pragma once

namespace eMELA {

  /// PDFs multiplied by x
  std::map<int, double> xDistributions(double x, double Q);
  /// PDFs multiplied by x (evolution basis)
  std::map<int, double> xDistributionsEv(double x, double Q);
  
}
