#pragma once
#include <molecule.hpp>
#include <gauxc/basisset.hpp>
#include <gauxc/molecule.hpp>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <map>
#include <iterator>
#include <cerr.hpp>
#include <molecule.hpp>
#include <basisset.hpp>
#include <Eigen/Dense>
#include <quantum.hpp>
#include <gauxc/xc_integrator.hpp>
#include <gauxc/xc_integrator/impl.hpp>

namespace ChronusQ {


  class GauXCUtils {

    public:

      GauXC::BasisSet<double> gbasis;
      GauXC::BasisSet<double> gpbasis;
      GauXC::Molecule gmol;
      std::vector<GauXC::XCTask> * lb_tasks;
      std::shared_ptr<GauXC::RuntimeEnvironment> grt;
      std::shared_ptr<GauXC::XCIntegrator<Eigen::MatrixXd>> integrator_pointer=nullptr;

      // Default ctors
      GauXCUtils() = default;
      GauXCUtils( GauXCUtils && )      = default;
      GauXCUtils( const GauXCUtils & ) = default;
      GauXCUtils& operator=(const GauXCUtils&) = default;

      static ExchCXX::Functional get_functional(std::string fname);
      static ExchCXX::Functional get_epcfunctional(std::string fname);
      static GauXC::Molecule make_gmol(const Molecule& molecule);
      static GauXC::BasisSet<double> make_gbasis(const BasisSet& basis);

  };


  /**
   * A struct to hold information pertaining to controlling GauXC
   * This may be more appropriate to put into singleslater/base.hpp???
   */
  struct GauXCOptions {
    bool useGPU = false;                
    float gpuMemFrac = 0.95;           
    double basisTol  = 1e-10;
    double pbasisTol = 1e-10;
    size_t batchSize = 4096;
    std::string funcName;
    std::string prot_funcName;
    ExchCXX::Spin xcSpin;
    ExchCXX::Backend xcBackend         = ExchCXX::Backend::libxc;
    GauXC::AtomicGridSizeDefault grid  = GauXC::AtomicGridSizeDefault::UltraFineGrid;  
    GauXC::PruningScheme pruningScheme = GauXC::PruningScheme::Unpruned; 
    GauXC::XCWeightAlg xcWeightAlg     = GauXC::XCWeightAlg::SSF;
    GauXC::RadialQuad radialQuad       = GauXC::RadialQuad::MurrayHandyLaming;
    std::string intKernel              = "Default";

    static const std::map<std::string, GauXC::AtomicGridSizeDefault> mg_map;
    static const std::map<std::string, GauXC::PruningScheme> prune_map;
    static const std::map<std::string, GauXC::XCWeightAlg> xcweight_map;
    static const std::map<std::string, GauXC::RadialQuad> radialquad_map;

    // Build a GauXCUtils class given using the parse GAUXCOptions
    std::shared_ptr<GauXCUtils> buildGauXCUtils( const std::shared_ptr<const BasisSet>& basis, const std::shared_ptr<const BasisSet>& basis2,
        const Molecule& mol, MPI_Comm comm);

    void printGauXCSettings(std::ostream&out);
  };

}
