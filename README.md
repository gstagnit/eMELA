# eMELA

``eMELA`` is a library that implements the evolution in pure QED of the unpolarised electron parton distribution functions (PDFs) up to next-to-leading logarithmic (NLL) approximation, according to the papers listed at the end of this file. ``eMELA`` includes the evolution with multiple fermion families and their mass thresholds, and it gives one the possibility of choosing among three different UV-renormalisation schemes (`MSBAR`: $\overline{\rm MS}$, `ALPMZ`: $\alpha(m_Z)$, `ALGMU`: $G_\mu$) and two different factorisation schemes (`MSBAR`: $\overline{\rm MS}$, `DELTA`: $\Delta$).

More in detail, ``eMELA`` is an improved QED-version of [MELA](https://github.com/vbertone/MELA). It consists of a Fortran code responsible for the numerical evolution of the PDFs, and a C++ wrapper that provides one with the analytical solutions in the asymptotic $z \to 1$ region. Moreover, the possibility is given to the user to output the PDFs as grids compliant with the [LHAPDF](https://lhapdf.hepforge.org/index.html) format, that can be employed at a later stage.  Regardless of whether the numerical solution is computed at runtime or read from the grids, eMELA offers the possibility to switch to the analytical solution in the asymptotic region. This must be considered as the default option when using the PDFs for physics simulations.

**Note**: eMELA supersedes [ePDF](https://github.com/gstagnit/ePDF), that was limited to the evolution with a single lepton in the $\overline{\rm MS}$ renormalisation and factorisation schemes.

## Download

``eMELA`` can be directly obtained from the github repository:

https://github.com/gstagnit/eMELA/releases

For the last development branch, one can clone the master code:

```Shell
git clone https://github.com/gstagnit/eMELA.git
```

## Dependencies

In order to install the code, [CMake](https://cmake.org/) is required (usually available by default on Unix systems). To use the grids, [LHAPDF](https://lhapdf.hepforge.org/index.html) is required (compatibility with version 6.3.0 has been tested).

## Installation

The code can be compiled using the following procedure:
```Shell
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=/eMELA/install/ -DWITH_LHAPDF=[ON|OFF] ..
make && make install
cd ..
```
with `/eMELA/install` denoting the folder where `eMELA` will be installed. By the default, if no prefix specification is given, the program will be installed in the `/usr/local` folder on Linux systems.

More configuration options (static vs. dynamic library, path to LHAPDF installation, ...) can be accessed through ``ccmake`` e.g.
```Shell
cd build
ccmake ..
```
Pre-computed grids with different choices of renormalisation and factorisation schemes will be installed in the  ``/eMELA/install/grids/`` folder.

## Examples

In case of successful installation, some test codes will be generated in the ``examples/`` folder.

## Usage

All the relevant functions are defined in a single header file, which has to be included as:
```C++
#include "eMELA/eMELA.hh"
using namespace eMELA;
```
### Initialization
One first needs to initialize the PDFs. This can be done in different ways.

The simplest way is with a call to the function `QuickInitialize`, which accepts four arguments:
```C++
void QuickInitialize(std::string const& pert_order = "NLL",
                     std::string const& fac_scheme = "DELTA",
                     std::string const& ren_scheme = "MSBAR",   
                     double const& alpha_mz = 1.0/127.95471413988483012);
```
where:
- Perturbative order of the evolution: `NLL`, `LL`
- Factorisation scheme: `DELTA` , `MSBAR`
- Renormalisation scheme: `MSBAR`, `ALPMZ`, `ALGMU`
- Value of $\alpha$ at the Z mass $m_Z$ = 91.1876

(The default `alpha_mz` is the value of $\alpha(m_Z)$ in the `MSBAR` renormalisation scheme.)

For instance, one could write:
```C++
QuickInitialize("NLL", "MSBAR", "ALPMZ", 1.0/128.940);
```

If one needs to modify other options, several setters function `SetXXX(...)` are provided, see the inline comments next to each function in the header file. Note that in this case, one should first call the function `SetDefaultParameters()`, modify what needed and at the end call the `InitializeEvolution()` function e.g.:
```C++
SetDefaultParameters();
SetXXX(...);
...
SetXXX(...);
InitializeEvolution();
```

Instead, if one prefers to use a grid, one should use either one of the two following functions:
```C++
void InitializeFromGrid(std::string const& pdfname);
void InitializeFromGrid(int const& lhaid);
```
In the first call, the name of the PDF should be provided; `eMELA` will look for a folder `pdfname` among the paths contained in the `LHAPDF_DATA_PATH` environment variable (see above). In the second call, the index of the PDF as declared in the `eMELA/install/grids/pdfsets.index` file should be provided; again, the `eMELA/install/grids` path has to be in the  `LHAPDF_DATA_PATH` variable. By loading the grid, all the parameters which have been used to build the PDF (which are stored in the info file of the PDF) are also loaded in `eMELA`.

### Get the electron PDFs

At this point, one can get the electron PDFs multiplied by a damping factor $(1-z)^{1-\gamma}$ (see below) by means of the function:
```C++
double CodePdf(const int & idx,
               const double & x, const double & omx,
               const double & Q, const double & gamma = 1.0);
```
or the function:
```C++
double GridPdf(const int & idx,
               const double & x, const double & omx,
               const double & Q, const double & gamma = 1.0);
```
The arguments are the same in the two functions, but in the former the PDFs are calculated at runtime, whereas in the latter the PDFs are read from the grid. Note that, in both cases, the PDF of the electron in the asymptotic region (by default defined as the region $1 - x < 10^{-8}$) is built by means of a switch between the numerical solution (calculated at runtime or read from the grid) and the asymptotic analytical solution.

If one is interested in the numerical solution or the asymptotic analytical solution separately, the functions
```C++
double CodeNumPdf(const int & idx, const double & x, const double & Q, const double & gamma = 1.0);
double GridNumPdf(const int & idx, const double & x, const double & Q, const double & gamma = 1.0);
```
and the function
```C++
double AsyPdf(const double & omx, const double & Q, const double & gamma = 1.0);
```
are provided, respectively.

The arguments are:
- `idx`: PDG index of the parton (e.g. 11 = electron, 22 = photon, -11 = positron, 1 = d quark, etc.)
- `x`: momentum fraction
- `omx`: simply $1-x$ (it is important to use directly this variable in the asymptotic region when the PDFs appear under integration, since values of $x$ very close to 1 are going to be probed).
- `Q`: factorisation scale
- `gamma`: the PDF returned is multiplied by $(1-z)^{1-\gamma}$.

Of course, the positron PDFs are simply given by charge conjugation e.g. the PDF of the positron from the positron beam is equal to the PDF of the electron from the electron beam etc.

### Write your own grid

In order to write a grid, it is first needed to initialize the PDF, as indicated above. Then the following function is provided:
```C++
WriteGrid(std::string const& name, bspdf* brem = NULL, ...)  
```
where the dots denote grid parameters related to the spacing of $x$ and $Q$ points, with sensible default values. See the header file for details.
It will create a folder `name` with the grid inside. The second argument allows one to write a grid with beamstrahlung effects. Its default value is `NULL` i.e. the grid will not have beamstrahlung. If a pointer to a beamstrahlung object is created (see above), then it can be passed as second argument to write the associated grid.

### Legacy LL PDFs

We provide also functions for the LL "legacy" PDFs adopted in the literature, see Appendix A of [arXiv:2206.XXXX](https://arxiv.org/abs/2105.XXX).
They can be called as
```C++
double LLPDF(int const& index, double const& x, double const& omx, double const& Q, double const& gamma = 1.0);
```
where `index` denotes the "scheme" of the LL PDFs:
- 0: collinear
- 1: beta
- 2: eta
- 3: mixed

Note that this function adopts the value of $\alpha$ (and its evolutions) as stored inside eMELA; hence, eMELA PDFs need to be initilised in any case.
All the parameters not related to the evolution of alpha will be ignored when calling ``LLPDF``. 

### Beamstrahlung

`eMELA` provides also PDFs with beamstrahlung effects, according to the procedure presented in [arXiv:2108.10261](https://arxiv.org/pdf/2108.10261.pdf). (All the equations refer to this paper.) Beamstrahlung effects are collider dependent, hence for any collider, one should implement a new class, derived from the `bspdf` class. As an example, the class `ilc500` is provided in the `eMELA` source code, and grids including beamstrahlung effects are included in the `eMELA/install/grids` folder.

Suppose to be interested in using the `NLL_DELTA_MSBAR_ILC500` grid. In order to use beamstrahlung, one needs to include the following headers:
```C++
#include "eMELA/bspdf.hh"
#include "eMELA/ilc500.hh"
```
where the second one has to be replaced with the header related to the specific collider.

Then, one should build a pointer to the `bspdf` object via:
```C++
auto ptr_bs = make_bspdf("NLL_DELTA_MSBAR_ILC500");
```
Finally, a single component can be accessed with:
```C++
double fcom = ptr_bs->get_pdf(icom, pid, bid, x, omx, Q, gamma);
```
where:
- `icom`: index of the component
- `pid`: PDG code of the particle
- `bid`: identifier of the beam (1 = electron, -1 = positron)
- `x` and `omx`: $x$ and $1-x$ respectively
- `Q`: factorisation scale
- `gamma`: as for the PDFs, the beamstrahlung component is multiplied by $(1-x)^{1-\gamma}$.

### Implementing new beamstrahlung type

In practice, other beamstrahlung type may be adopted, and here we document the necessary steps to implement it. 

The first step, clearly, is to prepare the beamstrahlung functions that want to be adopted. Note that such function should be written in a factorisable form, as Eq.(27). In the following we suppose that the function form is similar to Eq. (28), where all those f_01,f_10,f_00+,f_00- functions esstential the same f function with different parameters, as in Eq. (31-34). Other forms of the beamstrahlung functions can be adopted, but require more modification. (For more detail about that, please contact one of the author, see contacts below).

To implement a new beamstrahlung, a new class derived from the `bspdf` base class should be implemented. The `ilc500` class can be recognized as an example. Suppose the new class is called `newbs`, below we list several changes that should be made:
1. This function computes the coevolution between ISR and beamstrahlung functions:
 ```C++
void newbs::integrate_isr(int pid, double p, double q, double x, double omx, double Q, double beta);
```
Such coevolution is done via numerical integration. For numerical accuracy and performance, the domain of integration is separated into two intervals, and for each of them a variable transformation is performed to achieve a flatter integrand. The corresponding integrand are defined as `term1` and `term2` in the `ilc500` code. In particular, the line starting with `double g = ` is the computation of the single f function, and it should be changed into the new f function accordingly.

2. in the header file, corresponding values for parameters should be set accordingly.

3. in `bspdf.cc`, the function `internal_make_bspdf` should be changed to include the new beamspectrum_type.


## Linking `eMELA` to external software

The `eMELA-config` configuration script (installed in `/eMELA/install/bin`) can be used at the compilation stage to link `eMELA` to external software. It can be included in the `PATH` environment variable:
```Shell
export PATH=/eMELA/install/bin:$PATH
```
Then, supposing a `C++` program called `mymain.cc`, one could compile with:
```Shell
g++ `eMELA-config --cppflags` mymain.cc -o mymain `eMELA-config --ldflags`
```
or supposing a Fortran program called `mymain.f`:
```Shell
gfortran mymain.f -o mymain `eMELA-config --ldflags` -lstdc++
```
On some systems, it may be required to add the ``/eMELA/install/lib`` folder into the `LD_LIBRARY_PATH` environment variable (or `DYLD_LIBRARY_PATH` on Mac systems) e.g.
```Shell
export LD_LIBRARY_PATH=/eMELA/install/lib:$LD_LIBRARY_PATH
```

To use grids, one would need to:
- pass the LHAPDF compilear and linker flags (this can be done with the options `--cflags` and `--libs` of the script `lhapdf-config`);
- update the `LD_LIBRARY_PATH` with the path where the LHAPDF library is installed;
- modify the `LHAPDF` environment variable `LHAPDF_DATA_PATH` to include the path to the grids e.g. if one wishes to use the default grids shipped with `eMELA` (installed in `eMELA/install/grids`), hence one would need to execute
```Shell
export LHAPDF_DATA_PATH=/eMELA/install/grids:$LHAPDF_DATA_PATH
```

## References

Please cite the following references when using ``eMELA`` in scientific works:

- S. Frixione, *Initial conditions for electron and photon structure and fragmentation functions*, [arXiv:1909.03886](https://arxiv.org/abs/1909.03886)
- V. Bertone, M. Cacciari, S. Frixione, G. Stagnitto, *The partonic structure of the electron at the next-to-leading logarithmic accuracy in QED*, [arXiv:1911.12040](https://arxiv.org/abs/1911.12040)
- S. Frixione, *On factorisation schemes for the electron parton distribution functions in QED*, [arXiv:2105.06688](https://arxiv.org/abs/2105.06688)
- V. Bertone, M. Cacciari, S. Frixione, G. Stagnitto, X. Zhao, M. Zaro, *Improving methods and predictions at high-energy e+e- colliders within collinear factorisation*, [arXiv:2206.XXXX](https://arxiv.org/abs/2105.XXX)

## Contacts

- Valerio Bertone: valerio.bertone@cern.ch
- Giovanni Stagnitto: giovanni.stagnitto@physik.uzh.ch
- Xiaoran Zhao: xiaoran.zhao@uniroma3.it

## TODO

- Create an Authors file and a LICENSE file
- make sure that with `make install` the grids are copied

## FOLDER TO BE REMOVED FROM THE DISTRIBUTED VERSION

- bspdf: original Xiaoran's beamstrahlung code
- tmp: old unused files
- private: code forh checks, writing down the grids, etc.
- toymodel: calculate toy model cross section and plot it
