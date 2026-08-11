#pragma once
#include <molecule.hpp>
#include <gauxc/basisset.hpp>
#include <gauxc/molecule.hpp>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <map>
#include <vector>
#include <iterator>
#include <cerr.hpp>
#include <molecule.hpp>
#include <basisset.hpp>
#include <Eigen/Dense>
#include <quantum.hpp>
#include <gauxc/xc_integrator.hpp>
#include <gauxc/xc_integrator/impl.hpp>

namespace ChronusQ {

  struct QuantumSubsystem;

  class GauXCUtils {

    public:

      GauXC::BasisSet<double> gbasis;
      std::vector<GauXC::BasisSet<double>> gbases;
      GauXC::Molecule gmol;
      std::vector<GauXC::XCTask> * lb_tasks;
      std::shared_ptr<GauXC::RuntimeEnvironment> grt;
      std::shared_ptr<GauXC::LoadBalancer> load_balancer;
      std::shared_ptr<GauXC::XCIntegrator<Eigen::MatrixXd>> integrator_pointer=nullptr;
      GauXC::MultiParticleFunctionalSpec multiparticle_functional_spec;

      // Hybrid functional coefficients through ExchCXX
      // ExchCXX HybCoeffs definition for reference:
      // double alpha = 0.0; // the coefficient of HF part for global hybrid functionals or the coefficient of long-range HF part for range-separated functionals
      // double beta  = 0.0; // the deduction of the short-range part. Following the notation of libxc. So the real coefficient of short-range HF is alpha + beta
      // double omega = 0.0; // the range-separation parameter
      ExchCXX::HybCoeffs hybridCoefficients;

      // Using the same definition as ExchCXX
      bool isRangeSeparatedHybrid() const {
        return hybridCoefficients.beta != 0.0 or hybridCoefficients.omega != 0.0;
      }

      // Get xHFX from GauXC functional definition or input.
      double xHFX = 0.;

      // Default ctors
      GauXCUtils() = default;
      GauXCUtils( GauXCUtils && )      = default;
      GauXCUtils( const GauXCUtils & ) = default;
      GauXCUtils& operator=(const GauXCUtils&) = default;

      static ExchCXX::Functional get_functional(std::string fname);
      static ExchCXX::Functional get_epcfunctional(std::string fname);
      static bool is_range_separated(std::string fname);
      static GauXC::Molecule make_gmol(const Molecule& molecule);
      static GauXC::BasisSet<double> make_gbasis(const BasisSet& basis);
      static ExchCXX::XCKernel get_xckernel(std::string kernel, ExchCXX::Spin xcSpin, ExchCXX::Backend xcBackend );

  };


  /**
   * A struct to hold information pertaining to controlling GauXC
   * This may be more appropriate to put into singleslater/base.hpp???
   */
  struct GauXCOptions {
    bool useGPU = false;                
    float gpuMemFrac = 0.95;           
    double basisTol  = 1e-10;
    double otherBasisTol = 1e-10;
    size_t batchSize = 4096;
    std::string funcName;
    ExchCXX::Spin xcSpin;
    ExchCXX::Backend xcBackend         = ExchCXX::Backend::libxc;
    GauXC::AtomicGridSizeDefault grid  = GauXC::AtomicGridSizeDefault::UltraFineGrid;  
    GauXC::PruningScheme pruningScheme = GauXC::PruningScheme::Unpruned; 
    GauXC::XCWeightAlg xcWeightAlg     = GauXC::XCWeightAlg::SSF;
    GauXC::RadialQuad radialQuad       = GauXC::RadialQuad::MurrayHandyLaming;
    std::string intKernel              = "Default";
    GauXC::RadialSize nrad             = GauXC::RadialSize(99);
    GauXC::AngularSize nang            = GauXC::AngularSize(590);
    bool custom_grid                   = false;

    static const std::map<std::string, GauXC::AtomicGridSizeDefault> mg_map;
    static const std::map<std::string, GauXC::PruningScheme> prune_map;
    static const std::map<std::string, GauXC::XCWeightAlg> xcweight_map;
    static const std::map<std::string, GauXC::RadialQuad> radialquad_map;

    // BYO Functional
    ExchCXX::HybCoeffs hyb_coeffs = {0.0, 0.0, 0.0};
    std::string kernels;

    // Build a GauXCUtils class given using the parse GAUXCOptions
    std::shared_ptr<GauXCUtils> buildGauXCUtils( const std::shared_ptr<const BasisSet>& basis,
        const Molecule& mol, MPI_Comm comm);
    std::shared_ptr<GauXCUtils> buildGauXCUtils( const std::vector<QuantumSubsystem>& quantumSubsystems,
        const Molecule& mol, MPI_Comm comm);

    void printGauXCSettings(std::ostream&out,
        const std::vector<QuantumSubsystem>* quantumSubsystems = nullptr);
  };

}
