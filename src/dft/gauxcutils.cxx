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
    {"TREUTLERALDRICHS",    GauXC::RadialQuad::TreutlerAldrichs}
};


ExchCXX::Functional GauXCUtils::get_functional(std::string fname) {

  //std::transform(fname.begin(), fname.end(), fname.begin(),std::toupper);

  if (!fname.compare("BLYP")) {
    return ExchCXX::Functional::BLYP;
  } else if (!fname.compare("B3LYP")) {
    return ExchCXX::Functional::B3LYP;
  } else if (!fname.compare("PBEXPBEC")) {
    return ExchCXX::Functional::PBE;
  } else if (!fname.compare("revPBE")) {
    return ExchCXX::Functional::revPBE;
  } else if (!fname.compare("PBE0")) {
    return ExchCXX::Functional::PBE0;
  } else if (!fname.compare("SVWN5")) {
    return ExchCXX::Functional::SVWN5;
  } else if (!fname.compare("SVWN3")) {
    return ExchCXX::Functional::SVWN3;
  // MGGA functionals in GauXC temporarily not enabled for CQ
  // } else if (!fname.compare("SCAN")) {
  //   return ExchCXX::Functional::SCAN;
  // } else if (!fname.compare("R2SCAN")) {
  //   return ExchCXX::Functional::R2SCAN;
  // } else if (!fname.compare("R2SCANL")) {
  //   return ExchCXX::Functional::R2SCANL;
  }else {
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
