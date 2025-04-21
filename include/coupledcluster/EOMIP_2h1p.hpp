#pragma once 
#include <chronusq_sys.hpp>
#include <coupledcluster.hpp>

namespace ChronusQ{
  template <typename MatsT>
  EOMIP_2h1p<MatsT>::EOMIP_2h1p(const SafeFile &savFile,
    CCIntermediates<MatsT> &intermediates,
    const EOMSettings &eomSettings,
    const CoupledClusterSettings &ccSettings):
    EOMCCBase<MatsT>(savFile, intermediates,eomSettings,ccSettings),
    vLabel_(intermediates.vLabel), oLabel_(intermediates.oLabel),
    T1_(intermediates.T->get_tensor("OneBody")), T2_(intermediates.T->get_tensor("TwoBody")),
    tau(intermediates.tau),
    F_ae(intermediates.F_ae),
    F_mi(intermediates.F_mi),
    F_me(intermediates.F_me),
    W_mnie(intermediates.W_mnie),
    W_mbij(intermediates.W_mbij),
    W_mnij(intermediates.W_mnij),
    W_mbej(intermediates.W_mbej),
    W_amef(intermediates.W_amef) {
    TAManager &TAmanager = TAManager::get();

    // without L, we do not need D
    if (eomSettings.oscillator_strength == false) {
      if(intermediates.D_ai)   TAmanager.free("vo", std::move(intermediates.D_ai), true);
      if(intermediates.D_abij) TAmanager.free("vvoo", std::move(intermediates.D_abij), true);
    }

    nV_ = TAmanager.getRange(vLabel_).extent();
    nO_ = TAmanager.getRange(oLabel_).extent();
    nO2shift_ = nO_ * (nO_ - 1) / 2;
    this->Hbar_dimension_offsets.emplace("OneBody", 0);
    this->Hbar_dimension_offsets.emplace("TwoBody", nO_);
    this->Hbar_dim = nO_ + nV_ * nO2shift_;
    this->outOfBound_ = (this->Hbar_dim + 1) * (this->Hbar_dim + 1);
    ijIndices_.clear();
    ijIndices_.resize(nO_, std::vector<size_t>(nO_, outOfBound_));
    size_t idx = 0;
    for (size_t j = 0; j < nO_; j++) {
      for (size_t i = 0; i < std::min(j, nO_); i++) {
        ijIndices_[i][j] = idx;
        ijIndices_[j][i] = idx++;
      }
    }

    this->tensor_builder_.push_back(std::string(""));
    this->tensor_builder_.push_back(std::string({intermediates.oLabel}));
    this->tensor_builder_.push_back(std::string("OneBody"));
    this->tensor_builder_.push_back(std::string({intermediates.vLabel}));
    this->tensor_builder_.push_back(std::string({intermediates.oLabel,intermediates.oLabel}));
    this->tensor_builder_.push_back(std::string("TwoBody"));

  }

  template <typename MatsT>
  inline size_t EOMIP_2h1p<MatsT>::toCompoundS(size_t a) const {
    return a;
  }

  template <typename MatsT>
  inline size_t EOMIP_2h1p<MatsT>::toCompoundD(size_t a, size_t i, size_t j) const {
    size_t ij = ijIndices_[i][j];
    if (ij == outOfBound_)
      return outOfBound_;
    return a + ij * nV_;
  }

  template <typename MatsT>
  void EOMIP_2h1p<MatsT>::initializeEOMCC() {

    TAManager &TAmanager = TAManager::get();

    if (not W_mnie.is_initialized()){
      W_mnie = TAmanager.malloc<MatsT>("ooov");
    }

    if (not W_amef.is_initialized()){
      W_amef = TAmanager.malloc<MatsT>("vovv");
    }

    if (not W_mbij.is_initialized()){
      W_mbij = TAmanager.malloc<MatsT>("ovoo");
    }

  }

  template <typename MatsT>
  void EOMIP_2h1p<MatsT>::formF_ae() {
    F_ae("a,e") -= 0.5 * this->T1_("a,m") * F_me("m,e");
  }

  template <typename MatsT>
  void EOMIP_2h1p<MatsT>::formF_mi() {
    F_mi("m,i") += 0.5 * this->T1_("e,i") * F_me("m,e");
  }

  template <typename MatsT>
  void EOMIP_2h1p<MatsT>::formW_mnij() {
    W_mnij("m,n,i,j") += 0.25 * tau("e,f,i,j") * conj(this->antiSymMoints["vvoo"]("e,f,m,n"));
  }

  template <typename MatsT>
  void EOMIP_2h1p<MatsT>::formW_mbej() {
    W_mbej("m,b,e,j") -= 0.5 * this->T2_("f,b,j,n") * conj(this->antiSymMoints["vvoo"]("e,f,m,n"));
  }

  template <typename MatsT>
  void EOMIP_2h1p<MatsT>::formW_mnie() {
    W_mnie("m,n,i,e") = - conj(this->antiSymMoints["vooo"]("e,i,m,n")) + this->T1_("f,i") * conj(this->antiSymMoints["vvoo"]("f,e,m,n"));
  }

  template <typename MatsT>
  void EOMIP_2h1p<MatsT>::formW_amef() {
    W_amef("a,m,e,f") = conj(this->antiSymMoints["vvvo"]("e,f,a,m")) - this->T1_("a,n") * conj(this->antiSymMoints["vvoo"]("e,f,n,m"));
  }
  template <typename MatsT>
  void EOMIP_2h1p<MatsT>::formW_mbij() {
    W_mbij("m,b,i,j") = - this->antiSymMoints["vooo"]("b,m,i,j") - F_me("m,e") * this->T2_("b,e,i,j");
    W_mbij("m,b,i,j") += - this->T1_("b,n") * W_mnij("m,n,i,j");
    W_mbij("m,b,i,j") += - 0.5 * conj(this->antiSymMoints["vvvo"]("e,f,b,m")) * tau("e,f,i,j");
    W_mbij("m,b,i,j") += - conj(this->antiSymMoints["vooo"]("e,i,m,n")) * this->T2_("b,e,j,n");
    W_mbij("m,b,i,j") += conj(this->antiSymMoints["vooo"]("e,j,m,n")) * this->T2_("b,e,i,n");
    TArray tmp = TAManager::get().malloc<MatsT>("ovvo");
    tmp("m,b,e,j") = - this->antiSymMoints["vovo"]("b,m,e,j") - this->T2_("b,f,n,j") * conj(this->antiSymMoints["vvoo"]("e,f,m,n"));
    W_mbij("m,b,i,j") += this->T1_("e,i") * tmp("m,b,e,j");
    W_mbij("m,b,i,j") += - this->T1_("e,j") * tmp("m,b,e,i");
    TAManager::get().free("ovvo", std::move(tmp));
  }

  template <typename MatsT>
  void EOMIP_2h1p<MatsT>::formEOMIntermediates() {

    TAManager &TAmanager = TAManager::get();

    // If DFCCSD was used for ground state, initiate necessary slices
    if (this->ccSettings_.cctype == CC_TYPE::DFCCSD) {
      // Build slices of ERI needed to build intermediates
      this->antiSymMoints["vooo"] = TAmanager.malloc<MatsT>("vooo");
      this->antiSymMoints["vooo"]("a,n,i,j")  = this->riMoints["bvo"]("Q,a,i") * this->riMoints["boo"]("Q,n,j");
      this->antiSymMoints["vooo"]("a,n,i,j") -= this->antiSymMoints["vooo"]("a,n,j,i");

      this->antiSymMoints["vvoo"] = TAmanager.malloc<MatsT>("vvoo");
      this->antiSymMoints["vvoo"]("a,b,i,j")  = this->riMoints["bvo"]("Q,a,i") * this->riMoints["bvo"]("Q,b,j");
      this->antiSymMoints["vvoo"]("a,b,i,j") -= this->antiSymMoints["vvoo"]("a,b,j,i");

      this->antiSymMoints["vovo"] = TAmanager.malloc<MatsT>("vovo");
      this->antiSymMoints["vovo"]("a,m,e,i")  = this->riMoints["bvv"]("Q,a,e") * this->riMoints["boo"]("Q,m,i");
      this->antiSymMoints["vovo"]("a,m,e,i") -= this->riMoints["bvo"]("Q,a,i") * this->riMoints["bov"]("Q,m,e");

      this->antiSymMoints["vvvo"] = TAmanager.malloc<MatsT>("vvvo");
      this->antiSymMoints["vvvo"]("a,b,e,j")  = this->riMoints["bvv"]("Q,a,e") * this->riMoints["bvo"]("Q,b,j");
      this->antiSymMoints["vvvo"]("a,b,e,j") -= this->antiSymMoints["vvvo"]("b,a,e,j");
    }

    formF_ae();
    formF_mi();
    formW_mnij();
//    formW_abef();
    formW_mbej();
    formW_mnie();
    formW_amef();
    formW_mbij();
//    formW_abei();

    if (this->ccSettings_.cctype == CC_TYPE::DFCCSD) {
      // "vvoo" slice still needed in formR2_tilde and buildRightZeroBody, don't clear here
      TAmanager.free("vooo",std::move(this->antiSymMoints["vooo"]));
      TAmanager.free("vovo",std::move(this->antiSymMoints["vovo"]));
      TAmanager.free("vvvo",std::move(this->antiSymMoints["vvvo"]));
      TA::get_default_world().gop.fence();
    }
  }

  template <typename MatsT>
  void EOMIP_2h1p<MatsT>::buildSigma(const MBExpansion<MatsT> &V, MBExpansion<MatsT> &HV, EOMCCEigenVecType vecType) const {

    const TArray &V2 = V.get_tensor("TwoBody");
    TArray &HV2 = HV.get_tensor("TwoBody");
    const TArray &V1 = V.get_tensor("OneBody");
    TArray &HV1 = HV.get_tensor("OneBody");

    switch (vecType) {
      case EOMCCEigenVecType::RIGHT:
        formR_tilde(V1, V2, HV1, HV2);
        break;
      case EOMCCEigenVecType::LEFT:
        formL_tilde(V1, V2, HV1, HV2);
        break;
    }
  }

  template <typename MatsT>
  void EOMIP_2h1p<MatsT>::formR_tilde(const TArray &R1, const TArray &R2, TArray &sigmaR1, TArray &sigmaR2) const {
    TAManager &TAmanager = TAManager::get();// sigmaR1 = +1.00 f(i,i) R1(m)  // flops: o1v0L1 = o1v0L1 | mem: o1v0L1 = o1v0L1

    // Eq. (A7) without R_{bc}^{ijk}
    sigmaR1("i")  = -F_mi("m,i") * R1("m");
    sigmaR1("i") += F_me("m,e") * R2("e,i,m");
    sigmaR1("i") += -0.5 * W_mnie("m,n,i,e") * R2("e,m,n");

    // Eq. (A8) without R_{bc}^{ijk}
    TArray tmp_R2 = TAmanager.malloc<MatsT>("voo");
    tmp_R2("b,i,j")  = -0.5 * W_mbij("m,b,i,j") * R1("m");
    tmp_R2("b,i,j") += -F_mi("m,i") * R2("b,m,j");
    tmp_R2("b,i,j") += 0.5 * F_ae("b,e") * R2("e,i,j");
    tmp_R2("b,i,j") += 0.25 * W_mnij("m,n,i,j") * R2("b,m,n");
    tmp_R2("b,i,j") += -W_mbej("m,b,e,j") * R2("e,m,i");

    TArray I_e = TAmanager.malloc<MatsT>("v");
    I_e("e") = -0.5 * conj(this->antiSymMoints["vvoo"]("e,f,m,n")) * R2("f,m,n");
    tmp_R2("b,i,j") += 0.5 * I_e("e") * T2_("e,b,i,j");
    TAmanager.free("v", std::move(I_e));

    // A(ij)
    sigmaR2("b,i,j")  = tmp_R2("b,i,j");
    sigmaR2("b,i,j") -= tmp_R2("b,j,i");
    TAmanager.free("voo", std::move(tmp_R2));

  }

  template <typename MatsT>
  void EOMIP_2h1p<MatsT>::buildDiag(MatsT * diag, std::vector<double> eps) const {

    TAManager &TAmanager = TAManager::get();
    size_t n_v = TAmanager.getRange(vLabel_).extent();
    size_t n_o = TAmanager.getRange(oLabel_).extent();
    std::fill_n(diag, this->Hbar_dim, MatsT(0.0));
    MatsT * diag_i = diag;
    for (auto i = 0; i < n_o; i++){
      diag_i[toCompoundS(i)] = -eps[i];
    }
    MatsT * diag_aij = diag + this->Hbar_dimension_offsets.at("TwoBody");
    for (auto a = 0; a < n_v; a++){
      for (auto i = 0; i < n_o; i++){
        for (auto j = 0; j < i; j++){
          diag_aij[toCompoundD(a,i,j)] = eps[a+n_o] - eps[i] - eps[j];
        }
      }
    }
  }



  template <typename MatsT>
  void EOMIP_2h1p<MatsT>::runLambda(){} 


  template <typename MatsT>
  typename Davidson<dcomplex>::VecsGen_t EOMIP_2h1p<MatsT>::EmptyDavidsonVectorBuilder(){
      // Algorithm with implicit Hbar matrix
      typename Davidson<dcomplex>::VecsGen_t vecsGenEOM;
      if (this->eomSettings.hbar_type == EOM_HBAR_TYPE::IMPLICIT) {
        vecsGenEOM = [this](size_t nVec)->std::shared_ptr<SolverVectors<dcomplex>> {
          return std::make_shared<MBExpansionSet<dcomplex>>(this->tensor_builder_, nVec, this->savFile_);
        }; // implicit vecsGenerator

        return vecsGenEOM;
      }
      return vecsGenEOM;
  }

  template <typename MatsT>
  typename Davidson<dcomplex>::LinearTrans_t EOMIP_2h1p<MatsT>::DavidsonResidualBuilder(EOMCCEigenVecType &eigenVecType){
      if (this->eomSettings.hbar_type == EOM_HBAR_TYPE::IMPLICIT) {
        this->funcEOM = [this, &eigenVecType]( size_t nVec, SolverVectors<dcomplex> &V,
            SolverVectors<dcomplex> &AV) {

          MBExpansionSet<dcomplex> *V_ptr = nullptr, *AV_ptr = nullptr;
          size_t Vshift = 0, AVshift = 0;
          try {
            V_ptr = &dynamic_cast<MBExpansionSet<dcomplex>&>(V);
          } catch(const std::bad_cast& e) {
            SolverVectorsView<dcomplex>& V_view = dynamic_cast<SolverVectorsView<dcomplex>&>(V);
            V_ptr = &dynamic_cast<MBExpansionSet<dcomplex>&>(V_view.getVecs());
            Vshift = V_view.shift();
          }

          try {
            AV_ptr = &dynamic_cast<MBExpansionSet<dcomplex>&>(AV);
          } catch(const std::bad_cast& e) {
            SolverVectorsView<dcomplex>& AV_view = dynamic_cast<SolverVectorsView<dcomplex>&>(AV);
            AV_ptr = &dynamic_cast<MBExpansionSet<dcomplex>&>(AV_view.getVecs());
            AVshift = AV_view.shift();
          }

          for (size_t i = 0; i < nVec; i++) {
            const MBExpansion<dcomplex> &Vi = V_ptr->get(i + Vshift);
            MBExpansion<dcomplex> &AVi = AV_ptr->get(i + AVshift);
            buildSigma(Vi, AVi, eigenVecType);
            TA::get_default_world().gop.fence();
            AVi.enforceSymmetry();
          }

        }; // implicit sigmaBuilder
        return this->funcEOM;
      }
      return this->funcEOM;

  }
  template <typename MatsT>
  typename Davidson<dcomplex>::LinearTrans_t EOMIP_2h1p<MatsT>::DavidsonPreconditionerBuilder(dcomplex * curEig, dcomplex * eomDiag){

      double PCsmall = this->eomSettings.davidson_preCond_small;

      this->PCEOM = [this, eomDiag, curEig, PCsmall]( size_t nVec, SolverVectors<dcomplex> &V,
          SolverVectors<dcomplex> &AV) {

        AV.set_data(0, nVec, V, 0);

        MBExpansionSet<dcomplex> *AV_ptr = nullptr;
        size_t AVshift = 0;

        try {
          AV_ptr = &dynamic_cast<MBExpansionSet<dcomplex>&>(AV);
        } catch(const std::bad_cast& e) {
          SolverVectorsView<dcomplex>& AV_view = dynamic_cast<SolverVectorsView<dcomplex>&>(AV);
          AV_ptr = &dynamic_cast<MBExpansionSet<dcomplex>&>(AV_view.getVecs());
          AVshift = AV_view.shift();
        }

        for (size_t iVec = 0; iVec < nVec; iVec++) {

          MBExpansion<dcomplex> &curB = AV_ptr->get(iVec + AVshift);

          MatsT * Diag_A = eomDiag;
          TA::foreach_inplace(curB.get_tensor("OneBody"), [iVec, curEig, Diag_A, this, PCsmall](TA::TensorZ &tile){
            const auto& lobound = tile.range().lobound();
            const auto& upbound = tile.range().upbound();

            dcomplex denom = 0.0;
            std::vector<std::size_t> x{0};
            for(x[0] = lobound[0]; x[0] < upbound[0]; ++x[0]) {
              denom = curEig[iVec] - Diag_A[toCompoundS(x[0])];
              if (std::abs(denom) >= PCsmall) tile[x] /= denom;
            }
          });
          TA::get_default_world().gop.fence();

          MatsT * Diag_AIJ = eomDiag + this->Hbar_dimension_offsets.at("TwoBody");
          TA::foreach_inplace(curB.get_tensor("TwoBody"), [iVec, curEig, Diag_AIJ, this, PCsmall](TA::TensorZ &tile){
            const auto& lobound = tile.range().lobound();
            const auto& upbound = tile.range().upbound();

            dcomplex denom = 0.0;
            std::vector<std::size_t> x{0, 0, 0};
            for(x[0] = lobound[0]; x[0] < upbound[0]; ++x[0]){
              size_t a = x[0];
              for(x[1] = lobound[1]; x[1] < upbound[1]; ++x[1]){
                size_t i = x[1];
                for(x[2] = lobound[2]; x[2] < upbound[2]; ++x[2]){
                  if (x[1] == x[2])
                    continue;
                  size_t j = x[2];
                  denom = curEig[iVec] - Diag_AIJ[toCompoundD(a, i, j)];
                  if (std::abs(denom) >= PCsmall) tile[x] /= denom;
                }
              }
            }
          });
          TA::get_default_world().gop.fence();

          curB.enforceSymmetry();
        }
      }; // implicit preConditioner
      return this->PCEOM;
  }

  

  template <typename MatsT>
  EOMIP_2h1p<MatsT>::~EOMIP_2h1p() {

  }


}
