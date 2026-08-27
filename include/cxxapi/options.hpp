/* 
 *  This file is part of the Chronus Quantum (ChronusQ) software package
 *  
 *  Copyright (C) 2014-2022 Li Research Group (University of Washington)
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
#pragma once

#include <cxxapi/input.hpp>
#include <cxxapi/procedural.hpp>
#include <quantumsubsystems.hpp>
#include <molecule.hpp>
#include <basisset.hpp>
#include <integrals.hpp>
#include <singleslater.hpp>
#include <realtime.hpp>
#include <response.hpp>
#include <coupledcluster.hpp>
#include <posthartreefock.hpp>
#include <mcscf.hpp>
#include <regex>
#include <cubegen.hpp> 
#include <perturb.hpp>
#include <newperturb.hpp>
#include <mp.hpp>
#include <memory>
#include <physcon.hpp>

// Preprocessor directive to aid the digestion of optional 
// input arguments
#define OPTOPT(x) try{ x; } catch(...) { ; }

namespace ChronusQ {

  // TODO: Remove duplicated JobType enum
  enum class CQJobType {
    SCF,
    RT,
    LR,
    CC,
    EOMCC,
    MR,
    BOMD,
    EHRENFEST
  };

  /*************/
  /* OLD CODES */
  /*************/

  // Type of Job
  enum class JobType {
    SCF,
    RT,
    RESP,
    CC,
    EOMCC,
    CI,
    PT,
    MP2,
    BOMD,
    EHRENFEST,
    UNKNOWN
  };



  // Tedious, but there isn't an easier way to do this
  inline JobType parseJob(std::string jobStr) {
    JobType job;
    if( jobStr == "SCF" ) {
      job = JobType::SCF;
    }
    else if( jobStr == "RT" ) {
      job = JobType::RT;
    }
    else if( jobStr == "RESP" ) {
      job = JobType::RESP;
    }
    else if( jobStr == "CC" ) {
      job = JobType::CC;
    }
    else if( jobStr == "EOMCC" ) {
      job = JobType::EOMCC;
    }
    else if( jobStr == "BOMD" ) {
      job = JobType::BOMD;
    }
    else if( jobStr == "EHRENFEST" ) {
      job = JobType::EHRENFEST;
    }
    else if( jobStr == "CI" ) {
      job = JobType::CI;
    }
    else if( jobStr == "PERTURB" ) {
      job = JobType::PT;
    }
    else if( jobStr == "MP2" ) {
      job = JobType::MP2;
    }
    else {
      job = JobType::UNKNOWN;
    }
    return job;
  };

  // Function definitions ofr option parsing. 
  // See src/cxxapi/input/*opts.cxx for documentation

  // Parse the options relating to the Molecule object
  Molecule CQMoleculeOptions(std::ostream &, CQInputFile &, std::string &);

  std::set<std::string> CQMOLECULE_VALID(const std::map<std::string, std::string>& inputSection);

  void parseGeomInp(Molecule &, std::string &, std::ostream &, bool, bool);

  void parseGeomFchk(Molecule &, std::string &, std::ostream &, bool);

  RefOptions parseRef(std::ostream &, Molecule &, std::vector<std::string> &);

  void buildFunclist(std::vector<std::shared_ptr<DFTFunctional>> &,
    std::string);

  void parseIntParam(std::ostream &, CQInputFile &, IntegrationParam &);

  void parseHamiltonianOptions(std::ostream &, CQInputFile &, 
    BasisSet &basis, RefOptions &refOptions, HamiltonianOptions &hamiltonianOptions, std::string);

  bool parseAtomicType(std::ostream &, CQInputFile &, ATOMIC_X2C_TYPE &, std::string);

  // Parse the options relating to the BasisSet
  std::shared_ptr<BasisSet> CQBasisSetOptions(std::ostream &, CQInputFile &,
    Molecule &, std::string, const std::vector<size_t>& atomIndices = {});

  std::set<std::string> CQBASIS_VALID(const std::map<std::string, std::string>& inputSection);


  // Parse the options relating to the SingleSlaterOptions
  SingleSlaterOptions CQSingleSlaterOptions(
      std::ostream &, CQInputFile &, Molecule &, BasisSet &);

  template <typename MatsT, typename IntsT>
  std::shared_ptr<SingleSlaterBase> buildMultiParticleSS(
    std::ostream& out,
    Molecule& mol,
    std::vector<QuantumSubsystem>& quantumSubsystems,
    std::vector<QuantumPairInteraction>& quantumPairInteractions,
    const std::vector<std::shared_ptr<SingleSlaterBase>>& subSS,
    SCFControls scfControls);

  std::shared_ptr<SingleSlaterBase> CQMultiParticleSSOptions(
    std::ostream& out,
    CQInputFile& input,
    Molecule& mol,
    std::vector<QuantumSubsystem>& quantumSubsystems,
    std::vector<QuantumPairInteraction>& quantumPairInteractions,
    SCFControls scfControls);

  /// Resolve SCF guess options parsed from the input for each subsystem in a multiparticle calculation
  void resolveSubsystemGuessOptions(std::vector<QuantumSubsystem>& quantumSubsystems, const SCFControls& scfControls);

  std::set<std::string> CQQM_VALID(const std::map<std::string, std::string>& inputSection);
  std::set<std::string> CQQUANTUMSUBSYSTEMQM_VALID(const std::map<std::string, std::string>& inputSection);
  std::set<std::string> CQDFTINT_VALID(const std::map<std::string, std::string>& inputSection);

  // Parse RT options
  std::shared_ptr<TDEMFieldBase> parseRTField(std::string&, std::ostream& );

  void HandleRTInitState(std::ostream&, CQInputFile&, std::shared_ptr<RealTimeMultiSlaterVectorManagerBase>&);

  std::shared_ptr<RealTimeBase> CQRealTimeOptions(
    std::ostream &, CQInputFile &, std::shared_ptr<SingleSlaterBase> &,
    std::shared_ptr<MCWaveFunctionBase> &,
    std::shared_ptr<TDEMPerturbation>& ,
    EMPerturbation &
  );
  std::shared_ptr<RealTimeBase> CQRealTimeMultiSlaterOptions(
    std::ostream &, CQInputFile &, std::shared_ptr<SingleSlaterBase> &,
    std::shared_ptr<MCWaveFunctionBase> &,
    EMPerturbation &
  );

  std::set<std::string> CQRT_VALID(const std::map<std::string, std::string>& inputSection);
  void CQRTExpandSubsystemAlgorithms(TDSCFOptions&, const std::vector<QuantumSubsystem>&);

  // Parse Response options
  std::shared_ptr<ResponseBase> CQResponseOptions(
    std::ostream &, CQInputFile &, std::shared_ptr<SingleSlaterBase> &,
    EMPerturbation &
  );

  std::set<std::string> CQRESPONSE_VALID(const std::map<std::string, std::string>& inputSection);
  std::set<std::string> CQMOR_VALID(const std::map<std::string, std::string>& inputSection);

  // Parse integral options
  std::shared_ptr<IntegralsBase> CQIntsOptions(std::ostream &, 
    CQInputFile &, Molecule &,
    std::shared_ptr<BasisSet>, std::shared_ptr<BasisSet>,
    std::shared_ptr<BasisSet>, std::string int_sec = "INTS");

  std::set<std::string> CQINTS_VALID(const std::map<std::string, std::string>& inputSection);

  // Parse the external field options
  inline void handleField(const std::string& fieldInputStr, EMPerturbation& parsedField, const EMPerturbation& otherField = EMPerturbation()) {
      auto const regexOFF = std::regex("false|off",std::regex_constants::icase);
      if( std::regex_search(fieldInputStr, regexOFF) ) {
          return;
      } else if( fieldInputStr.empty() ) {
        parsedField.addField(otherField);
        return;
      }

      std::vector<std::string> tokens;
      split(tokens,fieldInputStr);

      if( tokens.size() < 4 )
        CErr(fieldInputStr + "is not a valid Field specification");

      std::string fieldTypeStr = tokens[0];

      EMFieldTyp fieldType;
      if( not fieldTypeStr.compare("ELECTRIC") )
        fieldType = Electric;
      else if( not fieldTypeStr.compare("MAGNETIC") )
        fieldType = Magnetic;
      else
        CErr(fieldTypeStr + "not a valid Field type");

      if( tokens.size() == 4 ) {
        cart_t field = {std::stod(tokens[1]), std::stod(tokens[2]),
                        std::stod(tokens[3])};
        parsedField.addField(fieldType,field);
      } else
        CErr("Non Dipole fields NYI");
  };

  // Parse the SCF options
  SCFControls CQSCFOptions(std::ostream&, CQInputFile&, EMPerturbation &);

  void HandleOrbitalSwaps(std::ostream&, CQInputFile&, SingleSlaterBase&, std::string);

  std::set<std::string> CQSCF_VALID(const std::map<std::string, std::string>& inputSection);

  // Parse Davidson energy specific settings
  size_t HandleNRootsInput(std::string,
                  std::vector<std::pair<double, size_t>> &);


  // Parse CC options
#ifdef CQ_HAS_TA
  CoupledClusterSettings CQCCOptions(std::ostream &, CQInputFile &);
  EOMSettings CQEOMCCOptions(std::ostream &, CQInputFile &);
#endif
  std::set<std::string> CQCC_VALID(const std::map<std::string, std::string>& inputSection);
  std::set<std::string> CQEOMCC_VALID(const std::map<std::string, std::string>& inputSection);

  // Parse geometry modifier options
  JobType CQGeometryOptions(std::ostream& out, CQInputFile& input, SafeFile& rstFile,
    JobType job, Molecule& mol, std::shared_ptr<SingleSlaterBase> ss,
    std::shared_ptr<MCWaveFunctionBase> mcscf,
    std::shared_ptr<RealTimeBase>& rt,
    std::shared_ptr<TDEMPerturbation>& tdPert,
    std::vector<QuantumSubsystem>& quantumSubsystems,
    std::vector<QuantumPairInteraction>& quantumPairInteractions,
    EMPerturbation& emPert, TDSCFOptions& tdSCFOptions,
    const std::shared_ptr<GauXCOptions>& gauxcOptions = nullptr);

  JobType CQDynamicsOptions(std::ostream& out, CQInputFile& input, SafeFile& rstFile,
    JobType job, Molecule& mol, std::shared_ptr<SingleSlaterBase> ss, std::shared_ptr<MCWaveFunctionBase> mcscf,
    std::shared_ptr<RealTimeBase>& rt,
    std::shared_ptr<TDEMPerturbation>& tdPert,
    std::vector<QuantumSubsystem>& quantumSubsystems,
    std::vector<QuantumPairInteraction>& quantumPairInteractions,
    EMPerturbation& emPert, TDSCFOptions& tdSCFOptions,
    const std::shared_ptr<GauXCOptions>& gauxcOptions = nullptr);

  std::set<std::string> CQDYNAMICS_VALID(const std::map<std::string, std::string>& inputSection);

  // Parse MCSCF options
  struct MCSCFJobType;
  std::shared_ptr<MCWaveFunctionBase> CQBuildMCWaveFunction(std::ostream &,
     CQInputFile &, std::shared_ptr<SingleSlaterBase> &, EMPerturbation &, std::shared_ptr<CubeGen>, const std::string &,
     std::shared_ptr<MCSCFJobType>&,std::shared_ptr<MCSCFSettings>);
  std::shared_ptr<MCSCFBase> CQMCSCFOptions(std::ostream &,
     CQInputFile &, std::shared_ptr<SingleSlaterBase> &, std::shared_ptr<MCWaveFunctionBase> &, EMPerturbation &, std::shared_ptr<CubeGen>, bool isNEO, std::vector<QuantumSubsystem>, std::vector<QuantumPairInteraction>);
  std::shared_ptr<MCSCFSettings> CQGetMCSCFSettings(std::ostream &out, CQInputFile &input, std::shared_ptr<MCSCFJobType> &mcscfjobtype, std::string prefix = "");
  std::shared_ptr<MCSCFBase> CQBuildMCSCFOptions(std::ostream &,
     CQInputFile &, std::shared_ptr<MCWaveFunctionBase> &, EMPerturbation &, std::shared_ptr<CubeGen>, std::string &,
     std::shared_ptr<MCSCFJobType>&mcscfjob,std::shared_ptr<MCSCFSettings>&);

  std::set<std::string> CQMCSCF_VALID(const std::map<std::string, std::string>& inputSection);
  
  void HandlePostHFProperties(std::ostream &, CQInputFile &,
    std::shared_ptr<PostHartreeFockBase> & postHF,
    std::string postHFSection);

   void HandleSavePDMSPostHF(std::ostream &, CQInputFile &,
    std::shared_ptr<PostHartreeFockBase> & postHF,
    std::string postHFSection);
  
  void HandlePostHFRDMPrinting(std::ostream &, CQInputFile &,
    std::shared_ptr<PostHartreeFockBase> & postHF,
    std::string postHFSection);

  void HandlePostHFOrbitalSwaps(std::ostream &out, CQInputFile &input,
    std::shared_ptr<SingleSlaterBase> &ss, 
    std::shared_ptr<PostHartreeFockBase> & postHF,
    std::string postHFSection);
  
  void ConstructActiveSpaces(std::ostream & out, CQInputFile & input,
                             const std::vector<size_t> & nActOs,
                             size_t nActE, size_t MOOffset,
                             int maxInterspaceEX,
                             std::vector<ActiveSpaceParameters> & actS,
                             std::vector<std::vector<size_t>> & refOcc,
                             std::string postHFSection);
  
  void ReadReferenceOcc(std::ostream & out, CQInputFile & input,
    std::vector<std::vector<size_t>> & refOcc, std::string postHFSection);

  std::shared_ptr<PostHartreeFockBase> CQCIOptions(std::ostream &,
    CQInputFile &, std::shared_ptr<SingleSlaterBase> &, EMPerturbation &, std::shared_ptr<CubeGen> cu);
  
  // Parse GauXC options                                                           
  GauXCOptions CQGauXCOptions(std::ostream&, CQInputFile &input, SingleSlaterOptions &ssOptions,
    const std::vector<QuantumSubsystem>* quantumSubsystems = nullptr);

  std::set<std::string> CQCI_VALID(const std::map<std::string, std::string>& inputSection);
  
  // Save reference info
  void saveRefs(SingleSlaterOptions &, std::shared_ptr<SingleSlaterBase> &);
  // NewPerturb Options
  std::shared_ptr<PostHartreeFockBase> CQMRPTSettings(std::ostream &,
            CQInputFile &, std::shared_ptr<PostHartreeFockBase> &);
  std::set<std::string> CQMRPT_VALID(const std::map<std::string, std::string>& inputSection);
  // Parse Perturb options
  std::shared_ptr<MCWaveFunctionBase> CQPerturbOptions(std::ostream &,
            CQInputFile &, std::shared_ptr<MCWaveFunctionBase> &);
  // Parse MP2 options
  std::shared_ptr<MP2Base> CQMP2Options(std::ostream &, CQInputFile &, 
              std::shared_ptr<SingleSlaterBase> &);


  std::set<std::string> CQMP2_VALID(const std::map<std::string, std::string>& inputSection);

  std::set<std::string> CQPERTURB_VALID(const std::map<std::string, std::string>& inputSection);

  void CQMiscOptions(std::ostream &, CQInputFile &);

  std::set<std::string> CQMISC_VALID(const std::map<std::string, std::string>& inputSection);

  std::shared_ptr<CubeGen> CQCUBEOptions(std::ostream&, CQInputFile&,
    std::shared_ptr<Molecule> mol, std::shared_ptr<BasisSet> &, EMPerturbation &, double );

  void CQCUBEOptionalKeywords(std::ostream&, CQInputFile&,
    CubeGenOptions&, std::string);

  void ParseCubeSubsection(std::ostream&, CQInputFile&, std::string,
    CubeGenOptions &, std::shared_ptr<CubeGen> cu);

  void ParseOrbitalPropSubsection(std::ostream&, CQInputFile&,
    std::shared_ptr<SingleSlaterBase> ss);

  std::set<std::string> CQCUBE_VALID(const std::map<std::string, std::string>& inputSection);

  std::set<std::string> CQORBPROP_VALID(const std::map<std::string, std::string>& inputSection);

  std::set<std::string> CQGAUXC_VALID(const std::map<std::string, std::string>& inputSection);

  std::set<std::string> CQPHYSCON_VALID(const std::map<std::string, std::string>& inputSection);

  void CQPhysConSetOptions(std::ostream &out, CQInputFile &input );

  void printInvalidKeys(const std::set<std::string> &invalidKeywords,
                        const std::string &prefix);

  std::set<std::string> CQInvalidKeywords(
      const std::set<std::string> &allowedKeywords,
      const std::map<std::string, std::string>& inputSection);

  void CQINPUT_VALID(std::ostream &out, CQInputFile &input);

  void copyIntermediates();

};




