#include "d3utils.hpp"
#include "dftd3.h" 
#include <molecule.hpp>
#include <algorithm>
#include <cctype>
#include <cerr.hpp>

namespace ChronusQ {
  
  // helper to turn dftd3_error into std::string
  static std::string d3_error_string(dftd3_error err) {
    constexpr int kBuf = 1024;
    char buf[kBuf] = {};
    int n = kBuf - 1;
    dftd3_get_error(err, buf, &n);
    return std::string(buf);
  }

  // helper for param loading
  dftd3_param load_param(dftd3_error err, ChronusQ::D3Utils::Scheme scheme, const std::string& method, bool atm) {
    auto m = const_cast<char*>(method.c_str());
    switch(scheme) {
      case ChronusQ::D3Utils::Scheme::D3BJ:    return dftd3_load_rational_damping(err, m, atm);
      case ChronusQ::D3Utils::Scheme::D3Zero:  return dftd3_load_zero_damping(err, m, atm);
      case ChronusQ::D3Utils::Scheme::D3MBJ:   return dftd3_load_mrational_damping(err, m, atm);
      case ChronusQ::D3Utils::Scheme::D3MZero: return dftd3_load_mzero_damping(err, m, atm);
      case ChronusQ::D3Utils::Scheme::D3OP:    return dftd3_load_optimizedpower_damping(err, m, atm);
    }
    return nullptr;
  } 

  void D3Utils::setSchemeFromString(const std::string& name_in) {
    std::string s(name_in);
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c){ return std::tolower(c); });
    
    if(s == "d3bj")    { scheme_ = Scheme::D3BJ;    return; }
    if(s == "d3zero")  { scheme_ = Scheme::D3Zero;  return; }
    if(s == "d3bjm")   { scheme_ = Scheme::D3MBJ;   return; }
    if(s == "d3zerom") { scheme_ = Scheme::D3MZero; return; }
    if(s == "d3op")    { scheme_ = Scheme::D3OP;    return; }
  
    CErr("Unknown D3 scheme: " + name_in, std::cout);
  }
  
  void D3Utils::evaluate(const Molecule& mol) {
    const int N = static_cast<int>(mol.nAtoms);
    std::vector<int> Z(mol.nAtoms);
    for (size_t i = 0; i < mol.nAtoms; ++i) Z[i] = mol.atoms[i].atomicNumber;
    std::vector<double> xyz = mol.getTotalCoordinates();
    evaluate(N, Z.data(), xyz.data());
  }

  // ---------------- Evaluate ----------------
  void D3Utils::evaluate(int natoms, const int* Z, const double* xyz) {

    dftd3_error     err   = dftd3_new_error();
    dftd3_structure mol   = dftd3_new_structure(err, natoms, Z, xyz, nullptr, nullptr);
    if(dftd3_check_error(err)) CErr("DFTD3 structure error: " + d3_error_string(err), std::cout);

    dftd3_model     model = dftd3_new_d3_model(err, mol);
    if(dftd3_check_error(err)) CErr("DFTD3 model error: " + d3_error_string(err), std::cout);

    std::string parameterRef = ref_;
    dftd3_param param = load_param(err, scheme_, parameterRef, useATM_);

    if (dftd3_check_error(err)) {
      parameterRef = "pbe";

      std::cout << "DFTD3 parameters for '" << ref_
                << "' were not found; defaulting to PBE parameters.\n";

      param = load_param(err, scheme_, parameterRef, useATM_);

      if (dftd3_check_error(err)) {
        CErr("DFTD3 PBE parameter error: " + d3_error_string(err),
            std::cout);
      }
    }
    
    // Pointers for derivative outputs
    double* gradPtr  = nullptr;
    double* sigmaPtr = nullptr;
    std::array<double,9> sigmaBuf{};  // sigma (strain/virial), must provide when gradient is requested

    if (computeGrad_) {
      result_.gradient.assign(3 * natoms, 0.0);
      gradPtr  = result_.gradient.data();
      sigmaPtr = sigmaBuf.data();
    }

    dftd3_get_dispersion(err, mol, model, param, &result_.energy, gradPtr, sigmaPtr);
    if(dftd3_check_error(err)) CErr("DFTD3 dispersion error: " + d3_error_string(err), std::cout);

    if (computePairwise_) {
      result_.pair2.assign(static_cast<size_t>(natoms) * static_cast<size_t>(natoms), 0.0);
      dftd3_get_pairwise_dispersion(err, mol, model, param, result_.pair2.data(), nullptr);
      if(dftd3_check_error(err)) CErr("DFTD3 pairwise error: " + d3_error_string(err), std::cout);
    }

    dftd3_delete_param(&param);
    dftd3_delete_model(&model);
    dftd3_delete_structure(&mol);
    dftd3_delete_error(&err);

}



} // namespace ChronusQ
