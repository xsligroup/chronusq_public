/*
 *  This file is part of the Chronus Quantum (ChronusQ) software package
 *  
 *  Copyright (C) 2014-2024 Li Research Group (University of Washington)
 *  
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *  
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License along
 *  with this program; if not, write to the Free Software Foundation, Inc.,
 *  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *  
 *  Contact the Developers:
 *    E-Mail: xsli@uw.edu
 *  
 */
#include <cxxapi/options.hpp>
#include <cerr.hpp>
#include <gauxcutils.hpp>
#include <gauxc/xc_integrator/integrator_factory.hpp>
#include <gauxc/molgrid/defaults.hpp>
#include <gauxc/molecular_weights.hpp>
#include <exchcxx/enums/spin.hpp>


namespace ChronusQ {


  /**
   *
   *  Check valid keywords in the section.
   *
  */
  std::set<std::string> CQGAUXC_VALID(const std::map<std::string, std::string>& inputSection) {

    // Allowed keywords
    std::set<std::string> allowedKeywords = {
        "GPU",                      // True or False
        "GPUMEMFRAC",               // Float between 0 and 1
        "BATCHSIZE",                // size_t
        "BASISTOL",                 // double
        "OTHERBASISTOL",            // double
        "PBASISTOL",                // double, legacy alias for OTHERBASISTOL
        "GRID",                     // string: fine, ultrafine, superfine, GM3, GM5
        "PRUNINGSCHEME",            // string: unpruned, robust, treutler
        "XCWEIGHTALG",              // string: Becke, SSF, LKO
        "RADIALQUAD",               // string: MuraKnowles,MurrayHandyLaming,TreutlerAhlrichs
        "XCBACKEND",                // string: libxc, builtin
        "INTKERNEL",                 // string: default, shellbatched (gpu only), incore (gpu only), reference (cpu only)
        "FUNCTIONAL",               // Build functional from list of allowed x and c kernels. 
        "XHFX",                     // Global hybrid scaling factor
    };

    return CQInvalidKeywords(allowedKeywords, inputSection);
  } //CQGAUXC_VALID

  /**
   *
   * Construct GauXCOptions using the input file
   *
   */
  GauXCOptions CQGauXCOptions(std::ostream&out , CQInputFile &input, SingleSlaterOptions &ssOptions,
    const std::vector<QuantumSubsystem>* quantumSubsystems) {
    
    out << "\nGauXC Settings:\n" << BannerTop << "\n\n" ;
    out << "  Parsing Input Options:" << std::endl;
    out << bannerMid << std::endl;
    GauXCOptions gauxcOpts;

    OPTOPT( gauxcOpts.useGPU     = input.getData<bool>("GAUXC/GPU"); )
    OPTOPT( gauxcOpts.gpuMemFrac = input.getData<double>("GAUXC/GPUMEMFRAC"); )
    // pin to [0,1]
    if (gauxcOpts.gpuMemFrac > 1.0)       gauxcOpts.gpuMemFrac = 1.0; 
    else if (gauxcOpts.gpuMemFrac < 0.0)  gauxcOpts.gpuMemFrac = 0.0; // TODO: test what happens at extremes
    
    // Parse basisset tolerance 
    OPTOPT( gauxcOpts.basisTol      = input.getData<double>("GAUXC/BASISTOL"); )
    OPTOPT( gauxcOpts.otherBasisTol = input.getData<double>("GAUXC/PBASISTOL"); )
    OPTOPT( gauxcOpts.otherBasisTol = input.getData<double>("GAUXC/OTHERBASISTOL"); )
    
    // Parse batch size
    OPTOPT( gauxcOpts.batchSize = input.getData<size_t>("GAUXC/BATCHSIZE"); )
    gauxcOpts.batchSize = gauxcOpts.batchSize > 1 ? gauxcOpts.batchSize : 1; // pin to [1,inf)

    // Initialize strings into which we will read users input
    std::string inputGrid, inputPruningScheme, inputXCWeightAlg, inputRadialQuad, inputIntKernel;

    // Parse atomic grid size, default to ultrafine
    OPTOPT( inputGrid = input.getData<std::string>("GAUXC/GRID"); )
    if (GauXCOptions::mg_map.find(inputGrid) == GauXCOptions::mg_map.end()){
      // Try to convert from grid specified in DFTInts section to GauXC grid
      if (inputGrid.empty()){
        std::pair<size_t, size_t> dftGrid(ssOptions.intParam.nRad, ssOptions.intParam.nAng);
        if (dftGrid == std::pair<size_t, size_t>(35, 110)){
          inputGrid = "GM3";
        } else if (dftGrid == std::pair<size_t, size_t>(50, 302)){
          inputGrid = "GM5";
        } else if (dftGrid == std::pair<size_t, size_t>(75, 302)){
          inputGrid = "FINE";
        } else if (dftGrid == std::pair<size_t, size_t>(99, 590)){
          inputGrid = "ULTRAFINE";
        } else if (dftGrid == std::pair<size_t, size_t>(250, 974)){
          inputGrid = "SUPERFINE";
        } else { 
          out << "  Warning: Grid not set in GauXC section: Using Custom Grid. " << std::endl;
          gauxcOpts.custom_grid = true;
        }
      } else{
        out << "  " << std::setw(39) << "Invalid GAUXC/GRID Keyword!"; 
        out << "  Set to default (Ultrafine)" << std::endl;
        inputGrid = "ULTRAFINE"; 
      }
    } 
    
    if(gauxcOpts.custom_grid){
      
      // Check if nAng input is a valid Lebedev number.
      size_t N = ssOptions.intParam.nAng;
      if( N!=6 && N!=14 && N!=26 && N!=38 && N!=50 && N!=74 && N!=86 && N!=110 && N!=146 &&
        N!=170 && N!=194 && N!=230 && N!=266 && N!=302 && N!=590 && N!=974 )
          CErr("Number of Angular Points NYI.");

      GauXC::RadialSize nRad = GauXC::RadialSize(ssOptions.intParam.nRad);
      GauXC::AngularSize nAng = GauXC::AngularSize(ssOptions.intParam.nAng);
      gauxcOpts.nrad = nRad;
      gauxcOpts.nang = nAng;
    } else{
      gauxcOpts.grid = GauXCOptions::mg_map.at(inputGrid);
    }

    // Parse pruning scheme, default to unpruned
    OPTOPT( inputPruningScheme = input.getData<std::string>("GAUXC/PRUNINGSCHEME"); )
    if (GauXCOptions::prune_map.find(inputPruningScheme) == GauXCOptions::prune_map.end()) {
      out << "  " << std::setw(39) << "PruningScheme not set or unrecognized;";
      out << "Set to default (Unpruned)" << std::endl;
      inputPruningScheme = "UNPRUNED"; } 
    gauxcOpts.pruningScheme = GauXCOptions::prune_map.at(inputPruningScheme);

    // Parse XC weight algorithm, default to SSF
    OPTOPT( inputXCWeightAlg = input.getData<std::string>("GAUXC/XCWEIGHTALG"); )
    if (GauXCOptions::xcweight_map.find(inputXCWeightAlg) == GauXCOptions::xcweight_map.end()) {
      out << "  " << std::setw(39) << "XCWeightAlg not set or unrecognized;";
      out << "Set to default (SSF)" << std::endl;
      inputXCWeightAlg = "SSF"; }
    gauxcOpts.xcWeightAlg = GauXCOptions::xcweight_map.at(inputXCWeightAlg);

    // Parse radial quadruture, default to MurrayHandyLaming
    OPTOPT( inputRadialQuad = input.getData<std::string>("GAUXC/RADIALQUAD"); )
    if (GauXCOptions::radialquad_map.find(inputRadialQuad) == GauXCOptions::radialquad_map.end()) {
      out << "  " << std::setw(39) << "RadialQuad not set or unrecognized;";
      out << "Set to default (MurrayHandyLaming)" << std::endl;
      inputRadialQuad = "MURRAYHANDYLAMING"; }
    gauxcOpts.radialQuad = GauXCOptions::radialquad_map.at(inputRadialQuad);

    // Parse integrator kernel. Input is sanitized on GauXC side
    OPTOPT( inputIntKernel = input.getData<std::string>("GAUXC/INTKERNEL"); )
    if( !inputIntKernel.empty() )
      gauxcOpts.intKernel = inputIntKernel;


    // Get XC functional info parsed in ssOptions
    gauxcOpts.funcName        =  ssOptions.refOptions.funcName;
    // Parse spin for ExchCXX
    gauxcOpts.xcSpin = ExchCXX::Spin::Unpolarized;
    if (ssOptions.refOptions.refType != isRRef) gauxcOpts.xcSpin = ExchCXX::Spin::Polarized; // UKS, GKS, 2C/X2C KS
    
    // Parse xc evalution backend for ExchCXX
    std::string inputXCBackend = "LIBXC";
    OPTOPT( inputXCBackend = input.getData<std::string>("GAUXC/XCBACKEND"); )
    if (not inputXCBackend.compare("LIBXC")) gauxcOpts.xcBackend = ExchCXX::Backend::libxc;
    else if (not inputXCBackend.compare("BUILTIN")) gauxcOpts.xcBackend = ExchCXX::Backend::builtin;
    else CErr(inputXCBackend + " not a valid GAUXC/XCBACKEND Keyword",out);
    
    out << std::endl;

    
    // Parse CUSTOM Functional from kernel list and hybridization
    OPTOPT( gauxcOpts.hyb_coeffs.alpha      = input.getData<double>("GAUXC/XHFX"); )
    OPTOPT( gauxcOpts.kernels = input.getData<std::string>("GAUXC/FUNCTIONAL"); )

    // Print full GauXC settings 
    gauxcOpts.printGauXCSettings(out, quantumSubsystems);

    out << std::endl << BannerEnd << std::endl;

    return gauxcOpts;
  } //CQGauXCOptions

  /**
   *
   * Construct GauXCUtils object from parsed options
   *
   */
  std::shared_ptr<GauXCUtils> GauXCOptions::buildGauXCUtils( const std::shared_ptr<const BasisSet>& basis,
    const Molecule& mol, MPI_Comm comm )
    {
      
      std::shared_ptr<GauXCUtils> gauxcUtils = std::make_shared<GauXCUtils>();

      // Generate molecule and basis 
      gauxcUtils->gmol   = gauxcUtils->make_gmol(mol);
      gauxcUtils->gbasis = gauxcUtils->make_gbasis(*basis);

      // Generate GauXC Runtime and select execution space / Kernel
      GauXC::ExecutionSpace exec_space = useGPU ? GauXC::ExecutionSpace::Device : GauXC::ExecutionSpace::Host;
      #ifdef CQ_ENABLE_CUDA
        #ifdef CQ_ENABLE_MPI
        gauxcUtils->grt = useGPU ? std::make_shared<GauXC::DeviceRuntimeEnvironment>(comm, gpuMemFrac) : std::make_shared<GauXC::RuntimeEnvironment>(comm);
        #else
        gauxcUtils->grt = useGPU ? std::make_shared<GauXC::DeviceRuntimeEnvironment>(gpuMemFrac) : std::make_shared<GauXC::RuntimeEnvironment>();
        #endif
      #else
        if(useGPU) CErr("useGPU enabled but CQ not compiled with CUDA! Set CQ_ENABLE_CUDA=ON");
        #ifdef CQ_ENABLE_MPI
        gauxcUtils->grt = std::make_shared<GauXC::RuntimeEnvironment>(comm);
        #else
        gauxcUtils->grt = std::make_shared<GauXC::RuntimeEnvironment>();
        #endif
      #endif

      // Set up molecular grid
      GauXC::MolGrid mg = custom_grid ? GauXC::MolGridFactory::create_default_molgrid(
            gauxcUtils->gmol, pruningScheme, GauXC::BatchSize(batchSize), radialQuad, nrad, nang) :
            GauXC::MolGridFactory::create_default_molgrid(
            gauxcUtils->gmol, pruningScheme, GauXC::BatchSize(batchSize), radialQuad, grid);

      // Screen shells 
      for( auto& sh : gauxcUtils->gbasis  ){ sh.set_shell_tolerance( basisTol ); }

      // Setup Load Balancer
      GauXC::LoadBalancerFactory lb_factory(exec_space, "Default");
      auto lb = lb_factory.get_shared_instance(*(gauxcUtils->grt), gauxcUtils->gmol, mg, gauxcUtils->gbasis);
      gauxcUtils->load_balancer = lb;
      auto& tasks = lb->get_tasks(); // Pregenerate tasks
      gauxcUtils->lb_tasks = &tasks; 

      // Transform weights
      GauXC::MolecularWeightsSettings mw_settings;
      mw_settings.weight_alg = xcWeightAlg;
      GauXC::MolecularWeightsFactory mw_factory( exec_space, "Default", mw_settings );    
      auto mw = mw_factory.get_instance();
      mw.modify_weights(*lb);

      // Build GauXC Integrator
      GauXC::XCIntegratorFactory<Eigen::MatrixXd> integrator_factory(exec_space, "Replicated", intKernel, "Default", "Default");  
      // Setup XC functional and build integrator
      GauXC::functional_type func;

      if (!funcName.compare("CUSTOM")) {
          std::vector<std::pair<double,ExchCXX::XCKernel>> kernel_list;
          std::stringstream ss(kernels);
          std::string kernel;
          double scale;

          while ( ss >> kernel >> scale){
            kernel_list.push_back({scale, gauxcUtils->get_xckernel(kernel, xcSpin, xcBackend)});
          }
          func = GauXC::functional_type( kernel_list, hyb_coeffs);

      } else {
          func = GauXC::functional_type( xcBackend, gauxcUtils->get_functional(funcName), xcSpin);
          hyb_coeffs = func.hyb_exx();
      }
      

      gauxcUtils->xHFX = hyb_coeffs.alpha;

      gauxcUtils->integrator_pointer = integrator_factory.get_shared_instance(func, lb);
      
      return gauxcUtils;

    } // End GauXCUtils Builder

  std::shared_ptr<GauXCUtils> GauXCOptions::buildGauXCUtils(
    const std::vector<QuantumSubsystem>& quantumSubsystems,
    const Molecule& mol, MPI_Comm comm )
    {

      std::vector<std::shared_ptr<const BasisSet>> bases;
      bases.reserve(quantumSubsystems.size());

      GauXC::MultiParticleFunctionalSpec functional_spec;
      functional_spec.intra_functionals.resize(quantumSubsystems.size());

      auto make_intra_functional = [&](const std::string& name, ExchCXX::Spin spin) {
        if(!name.compare("CUSTOM")) {
          std::vector<std::pair<double,ExchCXX::XCKernel>> kernel_list;
          std::stringstream ss(kernels);
          std::string kernel;
          double scale;

          while(ss >> kernel >> scale)
            kernel_list.push_back({scale, GauXCUtils::get_xckernel(kernel, spin, xcBackend)});

          return std::make_shared<GauXC::functional_type>(kernel_list, hyb_coeffs);
        } else {
          auto func = std::make_shared<GauXC::functional_type>(
            xcBackend, GauXCUtils::get_functional(name), spin);
          hyb_coeffs = func->hyb_exx();
          return func;
        }
      };

      // Procedural builds the electronic subsystem first. Its KS functional is
      // the intra-particle term; EPC refs on later subsystems become inter terms.
      for(size_t i = 0; i < quantumSubsystems.size(); ++i) {
        const auto& sys = quantumSubsystems[i];
        const auto& ref = sys.ssOptions.refOptions;
        bases.push_back(sys.basis);

        if(i == 0 and ref.isKSRef and not ref.isEPCRef) {
          ExchCXX::Spin spin = ref.refType != isRRef ?
            ExchCXX::Spin::Polarized : ExchCXX::Spin::Unpolarized;
          functional_spec.intra_functionals[i].push_back(make_intra_functional(ref.funcName, spin));
        }

        if(ref.isEPCRef) {
          GauXC::MultiParticlePairFunctional pair_spec;
          pair_spec.electron = 0;
          pair_spec.particle = i;
          pair_spec.functionals.push_back(std::make_shared<GauXC::functional_type>(
            ExchCXX::Backend::builtin, GauXCUtils::get_epcfunctional(ref.funcName),
            ExchCXX::Spin::Polarized));
          functional_spec.inter_functionals.push_back(std::move(pair_spec));
        }
      }

      if(bases.empty()) CErr("Cannot build MultiParticle GauXCUtils without basis sets");

      std::shared_ptr<GauXCUtils> gauxcUtils = std::make_shared<GauXCUtils>();
      gauxcUtils->multiparticle_functional_spec = functional_spec;
      gauxcUtils->xHFX = hyb_coeffs.alpha;

      gauxcUtils->gmol = gauxcUtils->make_gmol(mol);
      gauxcUtils->gbases.reserve(bases.size());
      for(const auto& basis : bases) {
        if(!basis) CErr("Cannot build MultiParticle GauXCUtils with null basis");
        gauxcUtils->gbases.emplace_back(gauxcUtils->make_gbasis(*basis));
      }
      if(!gauxcUtils->gbases.empty())
        gauxcUtils->gbasis = gauxcUtils->gbases.front();

      GauXC::ExecutionSpace exec_space = useGPU ? GauXC::ExecutionSpace::Device : GauXC::ExecutionSpace::Host;
      #ifdef CQ_ENABLE_CUDA
        #ifdef CQ_ENABLE_MPI
        gauxcUtils->grt = useGPU ? std::make_shared<GauXC::DeviceRuntimeEnvironment>(comm, gpuMemFrac) : std::make_shared<GauXC::RuntimeEnvironment>(comm);
        #else
        gauxcUtils->grt = useGPU ? std::make_shared<GauXC::DeviceRuntimeEnvironment>(gpuMemFrac) : std::make_shared<GauXC::RuntimeEnvironment>();
        #endif
      #else
        if(useGPU) CErr("useGPU enabled but CQ not compiled with CUDA! Set CQ_ENABLE_CUDA=ON");
        #ifdef CQ_ENABLE_MPI
        gauxcUtils->grt = std::make_shared<GauXC::RuntimeEnvironment>(comm);
        #else
        gauxcUtils->grt = std::make_shared<GauXC::RuntimeEnvironment>();
        #endif
      #endif

      GauXC::MolGrid mg = custom_grid ? GauXC::MolGridFactory::create_default_molgrid(
            gauxcUtils->gmol, pruningScheme, GauXC::BatchSize(batchSize), radialQuad, nrad, nang) :
            GauXC::MolGridFactory::create_default_molgrid(
            gauxcUtils->gmol, pruningScheme, GauXC::BatchSize(batchSize), radialQuad, grid);

      for(size_t i = 0; i < gauxcUtils->gbases.size(); ++i) {
        const double tol = i == 0 ? basisTol : otherBasisTol;
        for(auto& sh : gauxcUtils->gbases[i]) sh.set_shell_tolerance(tol);
      }

      GauXC::LoadBalancerFactory lb_factory(exec_space, "Default");
      auto lb = lb_factory.get_shared_instance(*(gauxcUtils->grt), gauxcUtils->gmol, mg, gauxcUtils->gbases);
      gauxcUtils->load_balancer = lb;
      auto& tasks = lb->get_tasks();
      gauxcUtils->lb_tasks = &tasks;

      GauXC::MolecularWeightsSettings mw_settings;
      mw_settings.weight_alg = xcWeightAlg;
      GauXC::MolecularWeightsFactory mw_factory( exec_space, "Default", mw_settings );
      auto mw = mw_factory.get_instance();
      mw.modify_weights(*lb);

      // GauXC needs a representative functional to construct the integrator.
      // The full intra/inter-particle spec is passed at each XC evaluation.
      std::shared_ptr<GauXC::functional_type> seed_func;
      for(const auto& funcs : functional_spec.intra_functionals) {
        if(!funcs.empty()) {
          seed_func = funcs.front();
          break;
        }
      }
      if(!seed_func) {
        for(const auto& pair : functional_spec.inter_functionals) {
          if(!pair.functionals.empty()) {
            seed_func = pair.functionals.front();
            break;
          }
        }
      }
      if(!seed_func) CErr("Cannot build MultiParticle GauXC integrator without functionals");

      GauXC::XCIntegratorFactory<Eigen::MatrixXd> integrator_factory(exec_space, "Replicated", intKernel, "Default", "Default");
      gauxcUtils->integrator_pointer = integrator_factory.get_shared_instance(seed_func, lb);

      return gauxcUtils;

    } // End MultiParticle GauXCUtils Builder

  void GauXCOptions::printGauXCSettings(std::ostream &out,
    const std::vector<QuantumSubsystem>* quantumSubsystems){

    size_t width = 28;
    bool doMultiParticle = quantumSubsystems and not quantumSubsystems->empty();

    out << "  Full GauXC Settings:" << std::endl;
    out << bannerMid << std::endl;

    out << "  " << std::setw(width) << "Use GPU:";
    out << (useGPU ? "True" : "False") << std::endl;
    
    out << "  " << std::setw(width) << "GPU Memory Fraction:";
    out << gpuMemFrac << std::endl;
    
    out << "  " << std::setw(width) << "Batch Size:";
    out << batchSize << std::endl;

    if(doMultiParticle) {
      auto labelOf = [](const QuantumSubsystem& sys) -> const std::string& {
        return sys.label.empty() ? sys.inputLabel : sys.label;
      };
      const std::string electronLabel = labelOf(quantumSubsystems->front());

      for(size_t i = 0; i < quantumSubsystems->size(); ++i) {
        const auto& sys = quantumSubsystems->at(i);
        const auto& ref = sys.ssOptions.refOptions;

        if(i == 0 and ref.isKSRef and not ref.isEPCRef) {
          out << "  " << std::setw(width) << "Intra Functional:";
          out << labelOf(sys) << " = " << ref.funcName << std::endl;
        }

        if(ref.isEPCRef) {
          out << "  " << std::setw(width) << "Inter Functional:";
          out << electronLabel << "-" << labelOf(sys) << " = " << ref.funcName << std::endl;
        }
      }
    } else {
      out << "  " << std::setw(width) << "XC Functional:";
      out << funcName << std::endl;
    }
    
    out << "  " << std::setw(width) << "XC Spin:";
    out << (xcSpin==ExchCXX::Spin::Unpolarized ? "Unpolarized" : "Polarized")  << std::endl;

    out << "  " << std::setw(width) << "XC Backend:";
    out << (xcBackend==ExchCXX::Backend::libxc ? "Libxc" : "Builtin")  << std::endl;

    out << "  " << std::setw(width) << "Integrator Kernel:";
    out << intKernel << std::endl;
    
    out << "  " << std::setw(width) << "Basis Tolerance:";
    out << basisTol << std::endl;

    if(doMultiParticle) {
      out << "  " << std::setw(width) << "Other Basis Tolerance:";
      out << otherBasisTol << std::endl;
    }

    out << "  " << std::setw(width) << "Grid:";
    out << ( custom_grid ? "Custom ("+std::to_string(nrad.get())+","+std::to_string(nang.get())+")"
        :   grid==GauXC::AtomicGridSizeDefault::UltraFineGrid ? "UltraFineGrid (99,590)" 
        :   grid==GauXC::AtomicGridSizeDefault::SuperFineGrid ? "SuperFineGrid (175(Z<2) or 250(Z>=2), 974)" 
        :   grid==GauXC::AtomicGridSizeDefault::FineGrid ?      "FineGrid (75,302)" 
        :   grid==GauXC::AtomicGridSizeDefault::GM3 ?           "GM3 (35,110)" 
        :                                                       "GM5 (50,302)" )  << std::endl;

    out << "  " << std::setw(width) << "Pruning Scheme:";
    out << (pruningScheme==GauXC::PruningScheme::Unpruned ? "Unpruned" 
        :   pruningScheme==GauXC::PruningScheme::Robust ?   "Robust" 
        :   "Treutler")  << std::endl;

    out << "  " << std::setw(width) << "XC Weight Algorithm:";
    out << (xcWeightAlg==GauXC::XCWeightAlg::SSF ? "SSF" 
        :   xcWeightAlg==GauXC::XCWeightAlg::Becke ? "Becke" 
        :   "LKO")  << std::endl;

    out << "  " << std::setw(width) << "Radial Quadrature:";
    out << (radialQuad==GauXC::RadialQuad::MurrayHandyLaming ? "MurrayHandyLaming" 
        :   radialQuad==GauXC::RadialQuad::MuraKnowles ?       "MuraKnowles" 
        :   "TreutlerAhlrichs")  << std::endl;

  }

}; // namespace ChronusQ
