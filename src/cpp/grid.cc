#ifdef USE_LHAPDF

#include <fstream>
#include <sys/stat.h>
#include <functional>

#include "eMELA/eMELA.hh"
#include "eMELA/grid.hh"
#include "eMELA/version.h"
#include "eMELA/numerics.hh"
#include "eMELA/bspdf.hh"

namespace eMELA {

  //_________________________________________________________________________________
  /// Real constructor
  void GridMELA::_GridMELA(LHAPDF::PDFInfo const& info)    
  {
    std::cout << "Reading INFO from grid" << std::endl;
    
    if ( info.get_entry("eMELA_version") != VERSION )
      throw std::runtime_error("[GridMELA] Grid written with different version of MELA!");

    SetDefaultParameters();
    
    SetPerturbativeOrder(info.get_entry_as<int>("eMELA_PerturbativeOrder"));
    SetFlavourSchemeInt(info.get_entry_as<int>("eMELA_FlavourSchemeInt"));
    SetRenormalisationSchemeInt(info.get_entry_as<int>("eMELA_RenormalisationSchemeInt"));
    SetFactorisationSchemeInt(info.get_entry_as<int>("eMELA_FactorisationSchemeInt"));
    SetWalpha(info.get_entry_as<int>("eMELA_Walpha"));
    SetPerturbativeOrderAlpha(info.get_entry_as<int>("eMELA_PerturbativeOrderAlpha"));    
    SetNPHOTOT(info.get_entry_as<int>("eMELA_NPHOTOT"));                         
    SetNMATEXP(info.get_entry_as<int>("eMELA_NMATEXP"));                         
    SetNSTPAEM(info.get_entry_as<int>("eMELA_NSTPAEM"));
    SetMINVMEL(info.get_entry_as<int>("eMELA_MINVMEL"));                  
    SetRINVMEL(info.get_entry_as<int>("eMELA_RINVMEL"));
    SetAFIXSOLint(info.get_entry_as<int>("eMELA_AFIXSOLint"));
    SetDELTAGMU(info.get_entry_as<double>("eMELA_DELTAGMU"));
    SetAlpha(info.get_entry_as<double>("eMELA_AlphaRef"),
    	   info.get_entry_as<double>("eMELA_AlphaQref"));
    SetThresholds(info.get_entry_as<std::vector<double>>("eMELA_Thresholds"));
    SetActiveFlavours(info.get_entry_as<std::vector<int>>("eMELA_ActiveFlavours"));
    SetActiveFlavoursAlpha(info.get_entry_as<std::vector<int>>("eMELA_ActiveFlavoursAlpha"));
    
    double ymax = info.get_entry_as<double>("YMax");
    _omxmax = pow(10,-ymax);
  }
  
  //_________________________________________________________________________________  
  GridMELA::GridMELA()
  {
    // set to null the pointer, no grid is used
    _lha = NULL;    
  }

  //_________________________________________________________________________________
  GridMELA::GridMELA(std::string const& pdfname) 
  {
    // Adding path to LHAPDF searching paths    
    // LHAPDF::pathsAppend(LHAPDF::dirname(LHAPDF::dirname(infopath)));
    // Better to use the $LHAPDF_DATA_PATH option
    
    _lha = new LHAPDF::GridPDF(pdfname, 0);
    LHAPDF::PDFInfo info = _lha->info();
    _GridMELA(info);
  }

  //_________________________________________________________________________________
  GridMELA::GridMELA(int const& lhaid) 
  {
    _lha = new LHAPDF::GridPDF(lhaid);
    LHAPDF::PDFInfo info = _lha->info();
    _GridMELA(info);
  }    
  
  //_________________________________________________________________________________
  GridMELA::~GridMELA()
  {
    // base class distructor already called by default
    delete _lha; _lha = NULL;
  }

   //_________________________________________________________________________________
  double GridMELA::NumGrid(const int & idx, const double & x, const double & Q,
			   const double & gamma)
  {
    if (_lha == NULL)
      throw std::runtime_error("[NumGrid] Using NumGrid, but no grid has been initialised");
    
    // We turn off the PDF for x greater than the last node in the grid to avoid errors like
    // Error in LHAPDF::ContinuationExtrapolator, x > xMax (last x knot)
    const double omx = 1.0-x;
    if (omx < _omxmax) {
      return 0.0;
    } else {
      if (idx == 11)  return _lha->xfxQ(x, Q).at(11)  * pow(omx, -gamma);
      if (idx == 22)  return _lha->xfxQ(x, Q).at(22)  * pow(omx, 1.0-gamma);
      if (idx == -11) return _lha->xfxQ(x, Q).at(-11) * pow(omx, 1.0-gamma);
      throw std::runtime_error("[NumGrid] idx should be 11, 22 or -11. idx = "+std::to_string(idx));
    }
  }  

  //_________________________________________________________________________________    
  // Utility function to get a double as a string with arbitrary precision
  std::string get_str(const double & val, const int & prec)
  {
    std::stringstream stream;
    stream << std::fixed << std::setprecision(prec) << val;
    return stream.str();
  }
  
  //_________________________________________________________________________________  
  void WriteGrid(std::string const& name,
		 bspdf*      brem,		 
		 int         const& nx   ,
		 double      const& xmin ,
		 double      const& xmed ,
		 double      const& ymax ,
		 double      const& frac ,
		 int         const& nQ   ,
		 double      const& Qmin ,
		 double      const& Qmax ,
		 bool        const& Qlin )
  {
    // Create a folder with the name of the setup and write the info file.
    if (mkdir(name.c_str(), 0777) != 0)
      throw std::runtime_error("[WriteGrid] cannot create folder for the set "+name+". If existing, not attempting overwriting.");
    
    // Information string.
    std::string info;

    // LHAPDF common names, trying to fill with useful info, but also to avoid errors
    info += "SetDesc: 'Set generated with eMELA. Name: " + name + "'\n";
    info += "Authors: V. Bertone, M. Cacciari, S. Frixione, G. Stagnitto, X. Zhao, M. Zaro\n";
    info += "References: arXiv:1909.03886, arXiv:1911.12040, arXiv:2105.06688, arXiv:XXXX\n";

    info += "Format: lhagrid1\n";
    info += "DataVersion: 1\n";
    if (brem != NULL)
      info += "NumMembers: " + std::to_string(1 + 2 * brem->components()) + "\n";
    else 
      info += "NumMembers: 1\n";
    info += "Particle: 11\n";
    info += "Flavors: [11, 22, -11]\n";
    info += "NumFlavors: 1\n";
    info += "FlavorScheme: fixed\n";    
    info += "XMin: " + get_str(xmin, 20) + "\n";
    info += "XMax: " + get_str(1-pow(10,-ymax), 20) + "\n";
    info += "YMax: " + get_str(ymax, 20) + "\n";
    info += "QMin: " + get_str(Qmin, 20) + "\n";
    info += "QMax: " + get_str(Qmax, 20) + "\n";

    if (brem != NULL)
      info += brem->beamspectrum_info();
    
    // Write eMELA info required to read the grid and set the parameters    
    std::map<std::string, std::string> eMELA_config;
    
    eMELA_config["eMELA_version"]                  = VERSION;
    eMELA_config["eMELA_PerturbativeOrder"]        = get_str(GetPerturbativeOrder(), 0);
    eMELA_config["eMELA_PerturbativeOrderAlpha"]   = get_str(GetPerturbativeOrderAlpha(), 0);
    eMELA_config["eMELA_FlavourSchemeInt"]         = get_str(GetFlavourSchemeInt(), 0);
    eMELA_config["eMELA_FactorisationSchemeInt"]   = get_str(GetFactorisationSchemeInt(), 0);
    eMELA_config["eMELA_RenormalisationSchemeInt"] = get_str(GetRenormalisationSchemeInt(), 0);    
    eMELA_config["eMELA_Walpha"]                   = get_str(GetWalpha(), 0);
    eMELA_config["eMELA_NPHOTOT"]                  = get_str(GetNPHOTOT(), 0);                  
    eMELA_config["eMELA_NMATEXP"]                  = get_str(GetNMATEXP(), 0);                  
    eMELA_config["eMELA_NSTPAEM"]                  = get_str(GetNSTPAEM(), 0);
    eMELA_config["eMELA_MINVMEL"]                  = get_str(GetMINVMEL(), 0);               
    eMELA_config["eMELA_RINVMEL"]                  = get_str(GetRINVMEL(), 0);
    eMELA_config["eMELA_AFIXSOLint"]               = get_str(GetAFIXSOLint(), 0);
    eMELA_config["eMELA_DELTAGMU"]                 = get_str(GetDELTAGMU(), 15);
    eMELA_config["eMELA_AlphaRef"]                 = get_str(GetAlphaRef(), 15);              
    eMELA_config["eMELA_AlphaQref"]                = get_str(GetAlphaQref(), 15);             

    for (auto const& x : eMELA_config) {
      info +=  x.first + ": " + x.second + "\n";
    }

    info += "eMELA_ActiveFlavours: [";
    std::vector<int> activeflavours = GetActiveFlavours();
    for (unsigned i = 0; i < activeflavours.size(); i++)
      info += get_str(activeflavours[i], 15) + ", ";
    //info = info.substr(0, info.size() - 2);
    info += "]\n";

    info += "eMELA_ActiveFlavoursAlpha: [";
    std::vector<int> activeflavoursalpha = GetActiveFlavoursAlpha();
    for (unsigned i = 0; i < activeflavoursalpha.size(); i++)
      info += get_str(activeflavoursalpha[i], 15) + ", ";
    //info = info.substr(0, info.size() - 2);
    info += "]\n";
    
    info += "eMELA_Thresholds: [";
    std::vector<double> thrs = GetThresholds();
    for (unsigned i = 0; i < thrs.size(); i++)
      info += get_str(thrs[i], 15) + ", ";
    //info = info.substr(0, info.size() - 2);
    info += "]\n";

    // Generate grid in Q.
    std::vector<double> Qgrid{Qmin};    
    if (Qlin == true && (nQ > 1)) {
      // We use a linearly spaced grid.
      const double Qstep = (Qmax - Qmin) / ( nQ - 1);
      for (int i = 1; i < nQ; i++)
	Qgrid.push_back(Qgrid.back() + Qstep);
    } else {
      // We use a logarithmically spaced grid.
      const double Qstep = exp( log(Qmax / Qmin) / ( nQ - 1 ) );
      for (int i = 1; i < nQ; i++)
	Qgrid.push_back(Qgrid.back() * Qstep);
    }
    
    info += "Qs_grid: [";
    for (unsigned i = 0; i < Qgrid.size(); i++)
      info += get_str(Qgrid[i], 15) + ", ";
    //info = info.substr(0, info.size() - 2);
    info += "]\n";
    
    //**********************************************************************
    // The stuf written below are not used, need to check what is actually
    // needed in order to prevent any LHAPDF error    
    // Quark masses to zero, to avoid 'LHAPDF::MetadataError'
    info += "MUp: 0\n";
    info += "MDown: 0\n";
    info += "MStrange: 0\n";
    info += "MCharm: 0\n";
    info += "MBottom: 0\n";
    info += "MTop: 0\n";    
    info += "MZ: 0.00000\n";
    info += "AlphaS_MZ: 0.0000\n";
    info += "AlphaS_OrderQCD: 0.00000\n";
    info += "AlphaS_Type: ipol\n";
    
    info += "AlphaS_Qs: [";
    for (unsigned i=0; i < Qgrid.size(); i++) info += "0.0, ";
    info = info.substr(0, info.size() - 2);
    info += "]\n";
    
    info += "AlphaS_Vals: [";
    for (unsigned i=0; i < Qgrid.size(); i++) info += "0.0, ";   
    info = info.substr(0, info.size() - 2);
    info += "]\n";
    
    // Open info file and print the information in it.
    std::ofstream out(name + "/" + name + ".info");
    out << info;
    out.close();

    // Generate grid in x.
    std::vector<double> xgrid{xmin};

    // Compute number of point in the linear and logarithmic grids
    const int nx1 = nx * frac;
    const int nx2 = nx - nx1;

    // Linear grid
    const double xstp1 = ( xmed - xmin ) / nx1;
    for (int ix = 1; ix < nx1; ix++)
      xgrid.push_back(xgrid.back() + xstp1);

    // Logarithmic grid
    const double howclosetoone = ymax;    
    const std::function<double(double const&)> fy = [] (double const& x) -> double{ return - log10( 1 - x ); };
    const std::function<double(double const&)> fx = [] (double const& y) -> double{ return 1 - pow(10, - y); };
    const double xstp2 = ( howclosetoone - fy(xmed) ) / ( nx2 - 1 );
    std::vector<double> txgrid{fy(xmed)};
    for (int ix = 1; ix < nx2; ix++)
      txgrid.push_back(txgrid.back() + xstp2);

    // Now compute grid in x.
    for (double tx : txgrid)
      xgrid.push_back(fx(tx));

    // If BS, then we have several grid files
    int upindex = (brem != NULL) ? 2 * brem->components() : 0;

    // Debug
    //int lowindex = (brem != NULL) ? 1 : 0;
    int lowindex = 0;
    
    for (int index = lowindex; index <= upindex; ++index) {
      
      std::ostringstream oss;
      oss << name << "_" << std::setfill('0') << std::setw(4) << index << ".dat";
      std::ofstream ofs(name + "/" + oss.str());
      ofs << std::scientific;
      ofs << std::setprecision(17);
      ofs << "PdfType: replica\n"
	  << "Format: lhagrid1\n"
	  << "---\n";

      // Write x grid.
      for (double x : xgrid)
	ofs << x << " ";
      ofs << "\n";
      
      // Write Q grid.
      for (double q : Qgrid)
	ofs << q << " ";
      ofs << "\n";
      
      // Write flavour indices.
      ofs << "11 22 -11\n";
      
      // Keep track of progress
      unsigned int totpoints = xgrid.size() * Qgrid.size();
      unsigned int counter = 0;      
      std::cout << "Start writing a total of " << totpoints
		<< " points in the grid " << oss.str() << std::endl;
      
      for (double x : xgrid) {
	double omx = 1.0-x;
	
	for (double Q : Qgrid) {

	  if (index == 0) { //original PDF without BS	    
	    //// Compute PDFs
	    std::map<int, double> xpdf = xDistributions(x, Q);
	    //// Multiply the electron PDF by (1.0-x), then photon
	    ofs << xpdf.at(1)/x*(1-x) << "\t" << xpdf.at(0)/x << "\t" << xpdf.at(-1)/x << std::endl;
	  } else {
	    int beamid = 2 * ((index - 1) / brem->components()) - 1;
	    int icom = (index - 1) % brem->components() + 1;	    
	    double r1 = brem->calc_pdf(icom, 11, beamid, x, omx, Q, ((beamid == 1) ? 0 : 1));
	    double r0 = brem->calc_pdf(icom, 22, beamid, x, omx, Q, 1);
	    double rm1 = brem->calc_pdf(icom, -11, beamid, x, omx, Q, ((beamid == -1) ? 0 : 1));
	    ofs << r1 << " " << r0 << " " << rm1 << std::endl;
	  }
	  
	  counter += 1;
	  if (counter % ((int) totpoints/10) == 0)
	    std::cout << "Written " << ((double) counter)/totpoints*100.0
		      << "% of points in the grid..." << std::endl;
	}
      }
    
      ofs << "---\n";
      ofs.close();
      
    } // end loop over index
    
    std::cout << "Grid written successfully!" << std::endl;
    
  }
  
}

#endif
