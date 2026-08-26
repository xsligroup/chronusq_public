#include <gauxcutils.hpp>
#include <gauxc/xc_integrator/integrator_factory.hpp>
#include <gauxc/molgrid/defaults.hpp>
#include <gauxc/molecular_weights.hpp>
#include <exchcxx/enums/spin.hpp>
namespace ChronusQ {

const std::map< std::string, GauXC::AtomicGridSizeDefault > GauXCOptions::mg_map = {
    {"FINE",          GauXC::AtomicGridSizeDefault::FineGrid},
    {"ULTRAFINE",     GauXC::AtomicGridSizeDefault::UltraFineGrid},
    {"SUPERFINE",     GauXC::AtomicGridSizeDefault::SuperFineGrid},
    {"GM3",           GauXC::AtomicGridSizeDefault::GM3},
    {"GM5",           GauXC::AtomicGridSizeDefault::GM5}
};
const std::map< std::string, GauXC::PruningScheme > GauXCOptions::prune_map = {
    {"UNPRUNED", GauXC::PruningScheme::Unpruned},
    {"ROBUST",   GauXC::PruningScheme::Robust},
    {"TREUTLER", GauXC::PruningScheme::Treutler}
};
const std::map< std::string, GauXC::XCWeightAlg > GauXCOptions::xcweight_map = {
    {"BECKE",  GauXC::XCWeightAlg::Becke},
    {"SSF",    GauXC::XCWeightAlg::SSF},
    {"LKO",    GauXC::XCWeightAlg::LKO}
};
const std::map< std::string, GauXC::RadialQuad > GauXCOptions::radialquad_map = {
    {"MURAKNOWLES",         GauXC::RadialQuad::MuraKnowles},
    {"MURRAYHANDYLAMING",   GauXC::RadialQuad::MurrayHandyLaming},
    {"TREUTLERAHLRICHS",    GauXC::RadialQuad::TreutlerAhlrichs}
};


ExchCXX::Functional GauXCUtils::get_functional(std::string fname) {

  //std::transform(fname.begin(), fname.end(), fname.begin(),std::toupper);

  if (!fname.compare("BLYP")) {
    return ExchCXX::Functional::BLYP;
  } else if (!fname.compare("B3LYP")) {
    return ExchCXX::Functional::B3LYP;
  } else if (!fname.compare("PBEXPBEC")) {
    return ExchCXX::Functional::PBE;
  } else if (!fname.compare("PBE")) {
    return ExchCXX::Functional::PBE;
  } else if (!fname.compare("revPBE")) {
    return ExchCXX::Functional::revPBE;
  } else if (!fname.compare("PBE0")) {
    return ExchCXX::Functional::PBE0;
  } else if (!fname.compare("SVWN5")) {
    return ExchCXX::Functional::SVWN5;
  } else if (!fname.compare("SVWN3")) {
    return ExchCXX::Functional::SVWN3;
  } else if (!fname.compare("LDA")) {
    return ExchCXX::Functional::LDA; 
  } else if (!fname.compare("BHANDH")) {
    return ExchCXX::Functional::BHANDH; 
  } else if (!fname.compare("CAMB3LYP")) {
    return ExchCXX::Functional::CAMB3LYP;
  } else if (!fname.compare("HSE06")) {
    return ExchCXX::Functional::HSE06;
  } else if (!fname.compare("LRCWPBE")) {
    return ExchCXX::Functional::LRCwPBE;
  } else if (!fname.compare("LCWPBE")) {
    return ExchCXX::Functional::LCwPBE;
  } else if (!fname.compare("WB97")) {
    return ExchCXX::Functional::wB97;
  } else if (!fname.compare("WB97X")) {
    return ExchCXX::Functional::wB97X;
  } else if (!fname.compare("PW91")) {
    return ExchCXX::Functional::PW91;
  } else if (!fname.compare("CUSTOM")) {
    CErr("Custom Requires Functional Definition under GAUXC header.");
  
  // MGGA functionals in GauXC temporarily not enabled for CQ
  // } else if (!fname.compare("SCAN")) {
  //   return ExchCXX::Functional::SCAN;
  // } else if (!fname.compare("R2SCAN")) {
  //   return ExchCXX::Functional::R2SCAN;
  // } else if (!fname.compare("R2SCANL")) {
  //   return ExchCXX::Functional::R2SCANL;
  } else {
    CErr("Invalid Functional for Gauxc");
  }

}

ExchCXX::Functional GauXCUtils::get_epcfunctional(std::string fname) {

  //std::transform(fname.begin(), fname.end(), fname.begin(),std::toupper);

  if (!fname.compare("EPC17_1")) {
    return ExchCXX::Functional::EPC17_1;
  } else if (!fname.compare("EPC17_2")) {
    return ExchCXX::Functional::EPC17_2;
  } else if (!fname.compare("EPC18_1")) {
    return ExchCXX::Functional::EPC18_1;
  } else if (!fname.compare("EPC18_2")) {
    return ExchCXX::Functional::EPC18_2;
  } else {
    CErr("Invalid EPCFunctional for Gauxc");
  }

}

bool GauXCUtils::is_range_separated(std::string fname) {
  // Derive from ExchCXX's own definition
  if (not ExchCXX::functional_map.key_exists(fname)) return false;
  try {
    // is_range_separated() only reads hyb_coefs_, which is spin-independent.
    // Build with Polarized so EPC kernels (which require polarized spin, e.g. EPC17_2) can be constructed
    // treat any functional that cannot be built as not range-separated.
    return ExchCXX::XCFunctional(ExchCXX::Backend::builtin, ExchCXX::functional_map.value(fname), ExchCXX::Spin::Polarized).is_range_separated();
  } catch (const std::exception&) {
    return false;
  }
}

ExchCXX::XCKernel GauXCUtils::get_xckernel(std::string kernel, ExchCXX::Spin xcSpin, ExchCXX::Backend xcBackend) {

  //std::transform(fname.begin(), fname.end(), fname.begin(),std::toupper);

  if (!kernel.compare("PBE_X")) {
    return ExchCXX::XCKernel( xcBackend, ExchCXX::Kernel::PBE_X, xcSpin);
  } else if (!kernel.compare("PBE_C")) {
    return ExchCXX::XCKernel( xcBackend, ExchCXX::Kernel::PBE_C, xcSpin);
  } else if (!kernel.compare("B88")) {
    return ExchCXX::XCKernel( xcBackend, ExchCXX::Kernel::B88, xcSpin);
  } else if (!kernel.compare("LYP")) {
    return ExchCXX::XCKernel( xcBackend, ExchCXX::Kernel::LYP, xcSpin);
  } else if (!kernel.compare("SLATEREXCHANGE")) {
    return ExchCXX::XCKernel( xcBackend, ExchCXX::Kernel::SlaterExchange, xcSpin);
  } else if (!kernel.compare("VWN3")) {
    return ExchCXX::XCKernel( xcBackend, ExchCXX::Kernel::VWN3, xcSpin);
  } else if (!kernel.compare("VWN5")) {
    return ExchCXX::XCKernel( xcBackend, ExchCXX::Kernel::VWN5, xcSpin);
  } else if (!kernel.compare("VWN")) {
    return ExchCXX::XCKernel( xcBackend, ExchCXX::Kernel::VWN, xcSpin);
  } else if (!kernel.compare("REVPBE_X")) {
    return ExchCXX::XCKernel( xcBackend, ExchCXX::Kernel::revPBE_X, xcSpin);
  } else if (!kernel.compare("OPTX_X")) {
    return ExchCXX::XCKernel( xcBackend, ExchCXX::Kernel::OPTX_X, xcSpin);
  } else if (!kernel.compare("PW91_X")) {
    return ExchCXX::XCKernel( xcBackend, ExchCXX::Kernel::PW91_X, xcSpin);
  } else if (!kernel.compare("PW91_C")) {
    return ExchCXX::XCKernel( xcBackend, ExchCXX::Kernel::PW91_C, xcSpin);
  } else if (!kernel.compare("PW91_LDA")) {
    return ExchCXX::XCKernel( xcBackend, ExchCXX::Kernel::PW91_LDA, xcSpin);
  } else if (!kernel.compare("PZ81")) {
    return ExchCXX::XCKernel( xcBackend, ExchCXX::Kernel::PZ81, xcSpin);
  } else if (!kernel.compare("ITYH_X")) {
    return ExchCXX::XCKernel( xcBackend, ExchCXX::Kernel::ITYH_X, xcSpin);
  } else if (!kernel.compare("P86_C")) {
    return ExchCXX::XCKernel( xcBackend, ExchCXX::Kernel::P86_C, xcSpin);
  } else if (!kernel.compare("B97_D")) {
    return ExchCXX::XCKernel( xcBackend, ExchCXX::Kernel::B97_D, xcSpin);
  } else {
    CErr("Invalid Kernel for Gauxc");
  }

}



GauXC::Molecule GauXCUtils::make_gmol(const Molecule& molecule) {

  GauXC::Molecule gmol;
  for (size_t i =0; i< molecule.nAtoms; i++) {

    gmol.emplace_back(GauXC::AtomicNumber(molecule.atoms[i].nucCharge),
        molecule.atoms[i].coord[0],
        molecule.atoms[i].coord[1],
        molecule.atoms[i].coord[2] );


  }

  return gmol;
}


GauXC::BasisSet<double> GauXCUtils::make_gbasis(const BasisSet& basis) {

  std::vector<GauXC::Shell<double>> gauxc_shell_vec;
  gauxc_shell_vec.reserve(basis.nShell);

  GauXC::Shell<double>::prim_array gauxc_alpha_arr;
  GauXC::Shell<double>::prim_array gauxc_coeff_arr;
  GauXC::Shell<double>::cart_array gauxc_cart_arr;

  for(const auto& shell : basis.shells) {
    for (auto a1 = 0; a1 < shell.nprim(); a1++) {
      gauxc_alpha_arr[a1] = shell.alpha[a1];
      gauxc_coeff_arr[a1] = shell.contr[0].coeff[a1];
    }
    for (auto a1 = 0; a1 < 3; a1++) {
      gauxc_cart_arr[a1] = shell.O[a1];
    }

    GauXC::Shell<double> add_gauxc_shell(
        GauXC::PrimSize(shell.nprim()),
        GauXC::AngularMomentum(shell.contr[0].l),
        GauXC::SphericalType(shell.contr[0].pure), gauxc_alpha_arr,
        gauxc_coeff_arr, gauxc_cart_arr, false);
    gauxc_shell_vec.emplace_back(std::move(add_gauxc_shell));
  }

  return GauXC::BasisSet<double>(gauxc_shell_vec);
} 

}
