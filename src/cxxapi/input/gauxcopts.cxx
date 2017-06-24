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
  void CQGAUXC_VALID( std::ostream &out, CQInputFile &input ) {

    // Allowed keywords
    std::vector<std::string> allowedKeywords = {
        "GPU",                      // True or False
        "GPUMEMFRAC",               // Float between 0 and 1
        "BATCHSIZE",                // size_t
        "BASISTOL",                 // double
        "PBASISTOL",                // double
        "GRID",                     // string: fine, ultrafine, superfine, GM3, GM5
        "PRUNINGSCHEME",            // string: unpruned, robust, treutler
        "XCWEIGHTALG",              // string: Becke, SSF, LKO
        "RADIALQUAD",               // string: MuraKnowles,MurrayHandyLaming,TreutlerAldrichs
        "XCBACKEND",                // string: libxc, builtin
        "INTKERNEL"                 // string: default, shellbatched (gpu only), incore (gpu only), reference (cpu only)
    };

    // Specified keywords
    std::vector<std::string> intsKeywords = input.getDataInSection("GAUXC");

    // Make sure all of basisKeywords in allowedKeywords
    for( auto &keyword : intsKeywords ) {
      auto ipos = std::find(allowedKeywords.begin(),allowedKeywords.end(),keyword);
      if( ipos == allowedKeywords.end() ) 
        CErr("Keyword GAUXC." + keyword + " is not recognized",std::cout);// Error
    }

  } //CQGAUXC_VALID

  /**
   *
   * Construct GauXCOptions using the input file
   *
   */
  GauXCOptions CQGauXCOptions(std::ostream&out , CQInputFile &input, SingleSlaterOptions &ssOptions,
    SingleSlaterOptions &prot_ssOptions) {
    
    out << "\nGauXC Settings:\n" << BannerTop << "\n\n" ;
    out << "  Parsing Input Options:" << std::endl;
    out << bannerMid << std::endl;
    GauXCOptions gauxcOpts;

    OPTOPT( gauxcOpts.useGPU     = input.getData<bool>("GAUXC.GPU"); )
    OPTOPT( gauxcOpts.gpuMemFrac = input.getData<double>("GAUXC.GPUMEMFRAC"); )
    // pin to [0,1]
    if (gauxcOpts.gpuMemFrac > 1.0)       gauxcOpts.gpuMemFrac = 1.0; 
    else if (gauxcOpts.gpuMemFrac < 0.0)  gauxcOpts.gpuMemFrac = 0.0; // TODO: test what happens at extremes
    
    // Parse basisset tolerance 
    OPTOPT( gauxcOpts.basisTol      = input.getData<double>("GAUXC.BASISTOL"); )
    OPTOPT( gauxcOpts.pbasisTol     = input.getData<double>("GAUXC.PBASISTOL"); )
    
    // Parse batch size
    OPTOPT( gauxcOpts.batchSize = input.getData<size_t>("GAUXC.BATCHSIZE"); )
    gauxcOpts.batchSize = gauxcOpts.batchSize > 1 ? gauxcOpts.batchSize : 1; // pin to [1,inf)

    // Initialize strings into which we will read users input
    std::string inputGrid, inputPruningScheme, inputXCWeightAlg, inputRadialQuad, inputIntKernel;

    // Parse atomic grid size, default to ultrafine
    OPTOPT( inputGrid = input.getData<std::string>("GAUXC.GRID"); )
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
          out << "  Warning: Grid not set in GauXC section " << std::endl;
          out << "           and Grid specified in 'DFTINTS' section is not available for GauXC" << std::endl;
          out << "  Set to default (Ultrafine)" << std::endl;
          inputGrid = "ULTRAFINE"; 
        }
      } else{
        out << "  " << std::setw(39) << "Invalid GAUXC.GRID Keyword!"; 
        out << "  Set to default (Ultrafine)" << std::endl;
        inputGrid = "ULTRAFINE"; 
      }
    } 
    gauxcOpts.grid = GauXCOptions::mg_map.at(inputGrid);

    // Parse pruning scheme, default to unpruned
    OPTOPT( inputPruningScheme = input.getData<std::string>("GAUXC.PRUNINGSCHEME"); )
    if (GauXCOptions::prune_map.find(inputPruningScheme) == GauXCOptions::prune_map.end()) {
      out << "  " << std::setw(39) << "PruningScheme not set or unrecognized;";
      out << "Set to default (Unpruned)" << std::endl;
      inputPruningScheme = "UNPRUNED"; } 
    gauxcOpts.pruningScheme = GauXCOptions::prune_map.at(inputPruningScheme);

    // Parse XC weight algorithm, default to SSF
    OPTOPT( inputXCWeightAlg = input.getData<std::string>("GAUXC.XCWEIGHTALG"); )
    if (GauXCOptions::xcweight_map.find(inputXCWeightAlg) == GauXCOptions::xcweight_map.end()) {
      out << "  " << std::setw(39) << "XCWeightAlg not set or unrecognized;";
      out << "Set to default (SSF)" << std::endl;
      inputXCWeightAlg = "SSF"; }
    gauxcOpts.xcWeightAlg = GauXCOptions::xcweight_map.at(inputXCWeightAlg);

    // Parse radial quadruture, default to MurrayHandyLaming
    OPTOPT( inputRadialQuad = input.getData<std::string>("GAUXC.RADIALQUAD"); )
    if (GauXCOptions::radialquad_map.find(inputRadialQuad) == GauXCOptions::radialquad_map.end()) {
      out << "  " << std::setw(39) << "RadialQuad not set or unrecognized;";
      out << "Set to default (MurrayHandyLaming)" << std::endl;
      inputRadialQuad = "MURRAYHANDYLAMING"; }
    gauxcOpts.radialQuad = GauXCOptions::radialquad_map.at(inputRadialQuad);

    // Parse integrator kernel. Input is sanitized on GauXC side
    OPTOPT( inputIntKernel = input.getData<std::string>("GAUXC.INTKERNEL"); )
    if( !inputIntKernel.empty() )
      gauxcOpts.intKernel = inputIntKernel;


    // Get XC functional info parsed in ssOptions
    gauxcOpts.funcName        =  ssOptions.refOptions.funcName;
    gauxcOpts.prot_funcName   =  prot_ssOptions.refOptions.funcName; 
    // Parse spin for ExchCXX
    gauxcOpts.xcSpin = ExchCXX::Spin::Unpolarized;
    if (ssOptions.refOptions.refType != isRRef) gauxcOpts.xcSpin = ExchCXX::Spin::Polarized; // UKS, GKS, 2C/X2C KS
    
    // Parse xc evalution backend for ExchCXX
    std::string inputXCBackend = "LIBXC";
    OPTOPT( inputXCBackend = input.getData<std::string>("GAUXC.XCBACKEND"); )
    if (not inputXCBackend.compare("LIBXC")) gauxcOpts.xcBackend = ExchCXX::Backend::libxc;
    else if (not inputXCBackend.compare("BUILTIN")) gauxcOpts.xcBackend = ExchCXX::Backend::builtin;
    else CErr(inputXCBackend + " not a valid GAUXC.XCBACKEND Keyword",out);
    
    out << std::endl;

    // Print full GauXC settings 
    gauxcOpts.printGauXCSettings(out);

    out << std::endl << BannerEnd << std::endl;

    return gauxcOpts;
  } //CQGauXCOptions

  /**
   *
   * Construct GauXCUtils object from parsed options
   *
   */
  std::shared_ptr<GauXCUtils> GauXCOptions::buildGauXCUtils( const std::shared_ptr<const BasisSet>& basis, const std::shared_ptr<const BasisSet>& pbasis,
    const Molecule& mol, MPI_Comm comm )
    {
      
      std::shared_ptr<GauXCUtils> gauxcUtils = std::make_shared<GauXCUtils>();
      
      // Decide whether to do NEO  
      bool doNEO = pbasis ? true : false; 

      // Generate molecule and basis 
      gauxcUtils->gmol   = gauxcUtils->make_gmol(mol);
      gauxcUtils->gbasis = gauxcUtils->make_gbasis(*basis);
      if(doNEO) gauxcUtils->gpbasis = gauxcUtils->make_gbasis(*pbasis);

      // Generate GauXC Runtime and select execution space / Kernel
      GauXC::ExecutionSpace exec_space = useGPU ? GauXC::ExecutionSpace::Device : GauXC::ExecutionSpace::Host;
      if(useGPU) CErr("GPU GauXC Not Yet Implemented!");
      #ifdef CQ_ENABLE_MPI
        gauxcUtils->grt = std::make_shared<GauXC::RuntimeEnvironment>(comm);
      #else
        gauxcUtils->grt = std::make_shared<GauXC::RuntimeEnvironment>();
      #endif

      // Set up molecular grid
      auto mg = GauXC::MolGridFactory::create_default_molgrid(
          gauxcUtils->gmol, pruningScheme, GauXC::BatchSize(batchSize), radialQuad, grid);

      // Screen shells 
      for( auto& sh : gauxcUtils->gbasis  ){ sh.set_shell_tolerance( basisTol ); }
      if(doNEO) for( auto& sh : gauxcUtils->gpbasis ){ sh.set_shell_tolerance( pbasisTol ); }

      // Setup Load Balancer
      GauXC::LoadBalancerFactory lb_factory(exec_space, "Default");
      std::shared_ptr<GauXC::LoadBalancer> lb;
      lb = doNEO ? lb_factory.get_shared_instance(*(gauxcUtils->grt), gauxcUtils->gmol, mg, gauxcUtils->gbasis, gauxcUtils->gpbasis)
          : lb_factory.get_shared_instance(*(gauxcUtils->grt), gauxcUtils->gmol, mg, gauxcUtils->gbasis);
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
      GauXC::functional_type func( xcBackend, gauxcUtils->get_functional(funcName), xcSpin);
      if (doNEO){
        // EPC functionals only has builtin implementations, and always be polarized
        GauXC::functional_type epcfunc( ExchCXX::Backend::builtin, gauxcUtils->get_epcfunctional(prot_funcName), ExchCXX::Spin::Polarized );
        gauxcUtils->integrator_pointer = integrator_factory.get_shared_instance(func, epcfunc, lb);
      } else {
        gauxcUtils->integrator_pointer = integrator_factory.get_shared_instance(func, lb);
      }
      
      return gauxcUtils;

    } // End GauXCUtils Builder

  void GauXCOptions::printGauXCSettings(std::ostream &out){

    size_t width = 28;
    bool doNEO = !prot_funcName.empty();

    out << "  Full GauXC Settings:" << std::endl;
    out << bannerMid << std::endl;

    out << "  " << std::setw(width) << "Use GPU:";
    out << (useGPU ? "True" : "False") << std::endl;
    
    out << "  " << std::setw(width) << "GPU Memory Fraction:";
    out << gpuMemFrac << std::endl;
    
    out << "  " << std::setw(width) << "Batch Size:";
    out << batchSize << std::endl;
    
    out << "  " << std::setw(width) << "XC Functional:";
    out << funcName << std::endl;
    
    out << "  " << std::setw(width) << "XC Spin:";
    out << (xcSpin==ExchCXX::Spin::Unpolarized ? "Unpolarized" : "Polarized")  << std::endl;

    out << "  " << std::setw(width) << "XC Backend:";
    out << (xcBackend==ExchCXX::Backend::libxc ? "Libxc" : "Builtin")  << std::endl;

    out << "  " << std::setw(width) << "Integrator Kernel:";
    out << intKernel << std::endl;
    
    out << "  " << std::setw(width) << "Basis Tolerance:";
    out << basisTol << std::endl;

    if(doNEO){
      out << "  " << std::setw(width) << "EPC Functional:";
      out << prot_funcName << std::endl;

      out << "  " << std::setw(width) << "Proton Basis Tolerance:";
      out << pbasisTol << std::endl;
    }

    out << "  " << std::setw(width) << "Grid:";
    out << (grid==GauXC::AtomicGridSizeDefault::UltraFineGrid ? "UltraFineGrid (99,590)" 
        :   grid==GauXC::AtomicGridSizeDefault::SuperFineGrid ? "SuperFineGrid (175(Z<2) or 250(Z>=2), 974)" 
        :   grid==GauXC::AtomicGridSizeDefault::FineGrid ?      "FineGrid (75,302)" 
        :   grid==GauXC::AtomicGridSizeDefault::GM3 ?           "GM3 (35,110)" 
        :   "GM5 (50,302)")  << std::endl;

    out << "  " << std::setw(width) << "Pruning Scheme:";
    out << (pruningScheme==GauXC::PruningScheme::Unpruned ? "Unpruned" 
        :   pruningScheme==GauXC::PruningScheme::Robust ?   "Robust" 
        :   "Treutler")  << std::endl;

    out << "  " << std::setw(width) << "XC Weight Algorithm:";
    out << (xcWeightAlg==GauXC::XCWeightAlg::SSF ? "SSF" 
        :   xcWeightAlg==GauXC::XCWeightAlg::Becke ? "Becke" 
        :   "LKO")  << std::endl;

    out << "  " << std::setw(width) << "Radial Quadruture:";
    out << (radialQuad==GauXC::RadialQuad::MurrayHandyLaming ? "MurrayHandyLaming" 
        :   radialQuad==GauXC::RadialQuad::MuraKnowles ?       "MuraKnowles" 
        :   "TreutlerAldrichs")  << std::endl;

  }

}; // namespace ChronusQ
