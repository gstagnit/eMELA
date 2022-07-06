#pragma once

#include <complex>

namespace eMELA
{
  /**
   * @brief function for the computation of the complex digamma function
   * @param z: complex argument
   * @return The digamma function in N
   */
  std::complex<double> psi(std::complex<double> const& z);

  /**
   * @brief function for the computation of the real digamma function
   * @param z: real argument
   * @return The digamma function in z
   */
  inline double psi(double const& z)
  {
    std::complex<double> zc(z, 0.0);
    return psi(zc).real();
  }

  /**
   * @brief function for the computation of the complex polygamma function
   * @param z: complex argument
   * @param m: positive integer degree of the polygamma
   * @return The polygamma of degree "m" function in "z"
   */
  std::complex<double> dpsi(std::complex<double> const& z, int const& m);

  /**
   * @brief function for the computation of the real digamma function
   * @param z: real argument
   * @return The digamma function in z
   */
  inline double dpsi(double const& z, int const& m)
  {
    std::complex<double> zc(z, 0.0);
    return dpsi(zc, m).real();
  }

  /**
   * @brief function for the computation of the Nielsen's generalized dilogs.
   * @param n: integer argument
   * @param p: integer argument
   * @param x: real argument
   * @return \f$\mathrm{S}_{n,p}(x)\f$
   * @note Implementation translated from CERNLIB WGPLG.
   */
  double wgplg(int const& n, int const& p, double const& x);

  /**
   * @brief function for the computation of the polylogarithms
   * @param n: integer argument
   * @param x: real argument
   * @return \f$\mathrm{Li}_{n}(x)\f$
   * @note Proxy for \f$\mathrm{S}_{n-1,1}(x)\f$
   */
  inline double polylog(int const& n, double const& x) { return wgplg(n-1,1,x); }

  /**
   * @brief function for the computation of the dilogarithm
   * @param x: real argument
   * @return \f$\mathrm{Li}_{2}(x)\f$
   */
  inline double dilog(double const& x) { return polylog(2,x); }

  /**
   * @brief function for the computation of the trilogarithm
   * @param x: real argument
   * @return \f$\mathrm{Li}_{3}(x)\f$
   */
  inline double trilog(double const& x) { return polylog(3,x); }
}
