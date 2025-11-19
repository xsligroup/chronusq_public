#pragma once
#include <memory>
#include <vector>
#include <array>
#include <string>
#include <stdexcept>

namespace ChronusQ {
  
  // Forward declaration
  struct Molecule;

  class D3Utils {
    public:
      enum class Scheme {
        D3BJ,     // "d3bj"    -> dftd3_load_rational_damping
        D3Zero,   // "d3zero"  -> dftd3_load_zero_damping
        D3MBJ,    // "d3bjm"   -> dftd3_load_mrational_damping
        D3MZero,  // "d3zerom" -> dftd3_load_mzero_damping
        D3OP      // "d3op"    -> dftd3_load_optimizedpower_damping
      };
    
      struct Result {
        double energy = 0.0;             
        std::vector<double> pair2;    // N×N row-major pairwise 2-body, compute if computePairwise_=true
        std::vector<double> gradient;
      };
    
      D3Utils() = default;
      D3Utils(std::string ref, Scheme scheme = Scheme::D3BJ,
              bool useATM = true, bool computePairwise = false, bool computeGrad = true)
      : ref_(std::move(ref)), scheme_(scheme),
        useATM_(useATM), computePairwise_(computePairwise), computeGrad_(computeGrad) {}
      
      D3Utils(std::string ref, std::string schemeString,
              bool useATM = true, bool computePairwise = false, bool computeGrad = true)
      : ref_(std::move(ref)), model_(schemeString),
        useATM_(useATM), computePairwise_(computePairwise), computeGrad_(computeGrad) {
          setSchemeFromString(schemeString);
        }
    
      // getters and setters
      const std::string& ref()   const { return ref_;   }
      const std::string& model() const { return model_; }
      Scheme          scheme()   const { return scheme_; }
      bool            useATM()   const { return useATM_; }
      bool   computePairwise()   const { return computePairwise_; }
      bool       computeGrad()   const { return computeGrad_; }
      const Result& result()     const { return result_; }
      void setRef(std::string r)       { ref_ = std::move(r); }
      void setModel(std::string m)     { model_ = m; }
      void setScheme(Scheme s)         { scheme_ = s; }
      void setSchemeFromString(const std::string& name);
      void setUseATM(bool b)           { useATM_ = b; }
      void setComputePairwise(bool b)  { computePairwise_ = b; }
      void setComputeGrad(bool b)      { computeGrad_ = b; }
    
      // Evaluate dispersion correction (non-PBC).
      // natoms: number of atoms
      // Z     : pointer to atomic numbers (length N)  
      // xyz   : pointer to interleaved positions (length 3N), in Bohr
      void evaluate(int natoms, const int* Z, const double* xyz);
      void evaluate(const Molecule& mol);
    
    private:
      std::string ref_;
      std::string model_;
      Scheme      scheme_;
      bool        useATM_ = true;
      bool        computePairwise_ = false;
      bool        computeGrad_ = true;
      Result      result_;
    };

}