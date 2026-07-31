#pragma once 
#include <chronusq_sys.hpp>
#include <coupledcluster.hpp>

namespace ChronusQ{
  template <typename MatsT>
  EOMEA<MatsT>::EOMEA(const SafeFile &savFile,
    CCIntermediates<MatsT> &intermediates,
    const EOMSettings &eomSettings,
    const CoupledClusterSettings &ccSettings):
    EOMCCBase<MatsT>(savFile, intermediates,eomSettings,ccSettings),
    vLabel_(intermediates.vLabel), oLabel_(intermediates.oLabel),
    reuse_tmps_(intermediates.sigmaOps),
    tmps_(intermediates.tempOps),
    T1_(intermediates.T->get_tensor("OneBody")), T2_(intermediates.T->get_tensor("TwoBody")) {
    TAManager &TAmanager = TAManager::get();

    // without L, we do not need D
    if (eomSettings.oscillator_strength == false) {
      if(intermediates.D_ai)   TAmanager.free("vo", std::move(intermediates.D_ai), true);
      if(intermediates.D_abij) TAmanager.free("vvoo", std::move(intermediates.D_abij), true);
    }

    nV_ = TAmanager.getRange(vLabel_).extent();
    nO_ = TAmanager.getRange(oLabel_).extent();
    nV2shift_ = nV_ * (nV_ - 1) / 2;
    this->Hbar_dimension_offsets.emplace("OneBody", 0);
    this->Hbar_dimension_offsets.emplace("TwoBody", nV_);
    this->Hbar_dim = nV_ + nO_ * nV2shift_;
    this->outOfBound_ = (this->Hbar_dim + 1) * (this->Hbar_dim + 1);
    abIndices_.clear();
    abIndices_.resize(nV_, std::vector<size_t>(nV_, outOfBound_));
    size_t idx = 0;
    for (size_t b = 0; b < nV_; b++) {
      for (size_t a = 0; a < std::min(b, nV_); a++) {
        abIndices_[a][b] = idx;
        abIndices_[b][a] = idx++;
      }
    }

    this->tensor_builder_.push_back(std::string({intermediates.vLabel}));
    this->tensor_builder_.push_back(std::string(""));
    this->tensor_builder_.push_back(std::string("OneBody"));
    this->tensor_builder_.push_back(std::string({intermediates.vLabel,intermediates.vLabel}));
    this->tensor_builder_.push_back(std::string({intermediates.oLabel}));
    this->tensor_builder_.push_back(std::string("TwoBody"));

  }

  template <typename MatsT>
  inline size_t EOMEA<MatsT>::toCompoundS(size_t a) const {
    return a;
  }

  template <typename MatsT>
  inline size_t EOMEA<MatsT>::toCompoundD(size_t a, size_t b, size_t i) const {
    size_t ab = abIndices_[a][b];
    if (ab == outOfBound_)
      return outOfBound_;
    return ab + i * nV2shift_;
  }

  template <typename MatsT>
  void EOMEA<MatsT>::initializeEOMCC() {

    TAManager &TAmanager = TAManager::get();

reuse_tmps_.emplace(std::make_pair("vvoo_1", TAmanager.malloc<MatsT>("vvoo")));
    reuse_tmps_.emplace(std::make_pair("vvoo_2", TAmanager.malloc<MatsT>("vvoo")));
    reuse_tmps_.emplace(std::make_pair("vo_30", TAmanager.malloc<MatsT>("vo")));
    reuse_tmps_.emplace(std::make_pair("vv_3", TAmanager.malloc<MatsT>("vv")));
    reuse_tmps_.emplace(std::make_pair("vvvv_4", TAmanager.malloc<MatsT>("vvvv")));
    reuse_tmps_.emplace(std::make_pair("vvvo_26", TAmanager.malloc<MatsT>("vvvo")));
    reuse_tmps_.emplace(std::make_pair("vvoo_5", TAmanager.malloc<MatsT>("vvoo")));
    reuse_tmps_.emplace(std::make_pair("vvoo_6", TAmanager.malloc<MatsT>("vvoo")));
    reuse_tmps_.emplace(std::make_pair("vvoo_7", TAmanager.malloc<MatsT>("vvoo")));
    reuse_tmps_.emplace(std::make_pair("vo_27", TAmanager.malloc<MatsT>("vo")));
    reuse_tmps_.emplace(std::make_pair("vv_8", TAmanager.malloc<MatsT>("vv")));
    reuse_tmps_.emplace(std::make_pair("vo_28", TAmanager.malloc<MatsT>("vo")));
    reuse_tmps_.emplace(std::make_pair("vv_9", TAmanager.malloc<MatsT>("vv")));
    reuse_tmps_.emplace(std::make_pair("oo_10", TAmanager.malloc<MatsT>("oo")));
    reuse_tmps_.emplace(std::make_pair("vo_11", TAmanager.malloc<MatsT>("vo")));
    reuse_tmps_.emplace(std::make_pair("vooo_12", TAmanager.malloc<MatsT>("vooo")));
    reuse_tmps_.emplace(std::make_pair("vooo_13", TAmanager.malloc<MatsT>("vooo")));
    reuse_tmps_.emplace(std::make_pair("oo_25", TAmanager.malloc<MatsT>("oo")));
    reuse_tmps_.emplace(std::make_pair("vo_33", TAmanager.malloc<MatsT>("vo")));
    reuse_tmps_.emplace(std::make_pair("vo_14", TAmanager.malloc<MatsT>("vo")));
    reuse_tmps_.emplace(std::make_pair("oo_15", TAmanager.malloc<MatsT>("oo")));
    reuse_tmps_.emplace(std::make_pair("vo_32", TAmanager.malloc<MatsT>("vo")));
    reuse_tmps_.emplace(std::make_pair("oo_16", TAmanager.malloc<MatsT>("oo")));
    reuse_tmps_.emplace(std::make_pair("vo_29", TAmanager.malloc<MatsT>("vo")));
    reuse_tmps_.emplace(std::make_pair("vvvo_17", TAmanager.malloc<MatsT>("vvvo")));
    reuse_tmps_.emplace(std::make_pair("vo_18", TAmanager.malloc<MatsT>("vo")));
    reuse_tmps_.emplace(std::make_pair("vvvo_19", TAmanager.malloc<MatsT>("vvvo")));
    reuse_tmps_.emplace(std::make_pair("oo_20", TAmanager.malloc<MatsT>("oo")));
    reuse_tmps_.emplace(std::make_pair("vo_31", TAmanager.malloc<MatsT>("vo")));
    reuse_tmps_.emplace(std::make_pair("vo_21", TAmanager.malloc<MatsT>("vo")));
    reuse_tmps_.emplace(std::make_pair("vo_22", TAmanager.malloc<MatsT>("vo")));
    reuse_tmps_.emplace(std::make_pair("vo_23", TAmanager.malloc<MatsT>("vo")));
    reuse_tmps_.emplace(std::make_pair("vo_24", TAmanager.malloc<MatsT>("vo")));

  }




  template <typename MatsT>
  void EOMEA<MatsT>::formEOMIntermediates() {

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

      this->antiSymMoints["vvvv"] = TAmanager.malloc<MatsT>("vvvv");
      this->antiSymMoints["vvvv"]("a,b,e,f")  = this->riMoints["bvv"]("Q,a,e") * this->riMoints["bvv"]("Q,b,f");
      this->antiSymMoints["vvvv"]("a,b,e,f") -= this->antiSymMoints["vvvv"]("a,b,f,e");
    }

    reuse_tmps_["vvoo_1"]("b,f,m,i")  = conj(this->antiSymMoints["vvvo"]("a,b,f,i")) * this->T1_("a,m");
    reuse_tmps_["vvoo_2"]("a,e,m,i")  = conj(this->antiSymMoints["vvvo"]("a,b,e,i")) * this->T1_("b,m");
    reuse_tmps_["vo_30"]("e,m")  = this->T1_("a,i") * reuse_tmps_["vvoo_2"]("a,e,m,i");
    reuse_tmps_["vv_3"]("b,e")  = conj(this->antiSymMoints["vvvo"]("a,b,e,i")) * this->T1_("a,i");
    reuse_tmps_["vvvv_4"]("f,e,b,a")  = conj(this->antiSymMoints["vvoo"]("a,b,j,i")) * this->T2_("e,f,j,i");
    reuse_tmps_["vvvo_26"]("b,e,f,m")  = this->T1_("a,m") * reuse_tmps_["vvvv_4"]("f,e,b,a");
    reuse_tmps_["vvoo_5"]("f,b,m,j")  = conj(this->antiSymMoints["vvoo"]("a,b,j,i")) * this->T2_("a,f,m,i");
    reuse_tmps_["vvoo_6"]("e,b,m,i")  = conj(this->antiSymMoints["vvoo"]("a,b,j,i")) * this->T2_("a,e,m,j");
    reuse_tmps_["vvoo_7"]("f,a,m,i")  = conj(this->antiSymMoints["vvoo"]("a,b,j,i")) * this->T2_("b,f,m,j");
    reuse_tmps_["vo_27"]("f,m")  = this->T1_("a,i") * reuse_tmps_["vvoo_7"]("f,a,m,i");
    reuse_tmps_["vv_8"]("f,a")  = 0.50 * conj(this->antiSymMoints["vvoo"]("a,b,j,i")) * this->T2_("b,f,j,i");
    reuse_tmps_["vo_28"]("e,m")  = this->T1_("a,m") * reuse_tmps_["vv_8"]("e,a");
    reuse_tmps_["vv_9"]("e,b")  = 0.50 * conj(this->antiSymMoints["vvoo"]("a,b,j,i")) * this->T2_("a,e,j,i");
    reuse_tmps_["oo_10"]("m,j")  = conj(this->antiSymMoints["vvoo"]("a,b,j,i")) * this->T2_("a,b,m,i");
    reuse_tmps_["vo_11"]("f,m")  = 0.50 * conj(this->antiSymMoints["vvvo"]("a,b,f,i")) * this->T2_("a,b,m,i");
    reuse_tmps_["vooo_12"]("b,m,i,j")  = conj(this->antiSymMoints["vvoo"]("a,b,j,i")) * this->T1_("a,m");
    reuse_tmps_["vooo_13"]("a,m,i,j")  = conj(this->antiSymMoints["vvoo"]("a,b,j,i")) * this->T1_("b,m");
    reuse_tmps_["oo_25"]("j,m")  = this->T1_("a,i") * reuse_tmps_["vooo_13"]("a,m,i,j");
    reuse_tmps_["vo_33"]("e,m")  = this->T1_("e,j") * reuse_tmps_["oo_25"]("j,m");
    reuse_tmps_["vo_14"]("b,j")  = conj(this->antiSymMoints["vvoo"]("a,b,j,i")) * this->T1_("a,i");
    reuse_tmps_["oo_15"]("m,j")  = conj(this->antiSymMoints["vooo"]("a,m,j,i")) * this->T1_("a,i");
    reuse_tmps_["vo_32"]("f,m")  = this->T1_("f,j") * reuse_tmps_["oo_15"]("m,j");
    reuse_tmps_["oo_16"]("m,i")  = 0.50 * conj(this->antiSymMoints["vvoo"]("a,b,j,i")) * this->T2_("a,b,m,j");
    reuse_tmps_["vo_29"]("f,m")  = this->T1_("f,i") * reuse_tmps_["oo_16"]("m,i");
    reuse_tmps_["vvvo_17"]("f,e,a,m")  = conj(this->antiSymMoints["vooo"]("a,m,j,i")) * this->T2_("e,f,j,i");
    reuse_tmps_["vo_18"]("f,m")  = 0.50 * conj(this->antiSymMoints["vooo"]("a,m,j,i")) * this->T2_("a,f,j,i");
    reuse_tmps_["vvvo_19"]("b,f,e,m")  = this->antiSymMoints["vvvv"]("e,f,a,b") * this->T1_("a,m");
    reuse_tmps_["oo_20"]("m,i")  = this->fockMatrix_ta["ov"]("i,a") * this->T1_("a,m");
    reuse_tmps_["vo_31"]("e,m")  = this->T1_("e,i") * reuse_tmps_["oo_20"]("m,i");
    reuse_tmps_["vo_21"]("e,m")  = this->antiSymMoints["vovo"]("e,i,a,m") * this->T1_("a,i");
    reuse_tmps_["vo_22"]("f,m")  = this->fockMatrix_ta["ov"]("i,a") * this->T2_("a,f,m,i");
    reuse_tmps_["vo_23"]("e,m")  = this->fockMatrix_ta["vv"]("e,a") * this->T1_("a,m");
    reuse_tmps_["vo_24"]("e,m")  = this->fockMatrix_ta["oo"]("i,m") * this->T1_("e,i");
    
  }

  template <typename MatsT>
  void EOMEA<MatsT>::buildSigma(const MBExpansion<MatsT> &V, MBExpansion<MatsT> &HV, EOMCCEigenVecType vecType) const {

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
  void EOMEA<MatsT>::formR_tilde(const TArray &R1, const TArray &R2, TArray &sigmaR1, TArray &sigmaR2) const {
    TAManager &TAmanager = TAManager::get();// sigmaR1 = +1.00 f(i,i) R1(e)  // flops: o0v1L1 = o0v1L1 | mem: o0v1L1 = o0v1L1
    
    // sigmaR1 += +1.00 f(i,a) this->T1_(a,i) R1(e)  // flops: o0v1L1 += o0v1L1 | mem: o0v1L1 += o0v1L1
    
    // sigmaR1 += -0.50 <j,i||j,i> R1(e)  // flops: o0v1L1 += o0v1L1 | mem: o0v1L1 += o0v1L1
    
    // sigmaR1 += -0.50 <j,i||a,b> this->T1_(a,i) this->T1_(b,j) R1(e)  // flops: o0v1L1 += o0v1L1 | mem: o0v1L1 += o0v1L1
    
    // sigmaR1 += +0.25 <j,i||a,b> this->T2_(a,b,j,i) R1(e)  // flops: o0v1L1 += o0v1L1 | mem: o0v1L1 += o0v1L1
    
    // sigmaR1 += +1.00 <i,e||a,b> this->T1_(a,i) R1(b)  // flops: o0v1L1 += o0v2L1 | mem: o0v1L1 += o0v1L1
    sigmaR1("e") = -reuse_tmps_["vv_3"]("b,e") * R1("b");
    
    // sigmaR1 += -0.50 <j,i||a,b> this->T2_(a,e,j,i) R1(b)  // flops: o0v1L1 += o0v2L1 | mem: o0v1L1 += o0v1L1
    sigmaR1("e") -= reuse_tmps_["vv_9"]("e,b") * R1("b");
    
    // sigmaR1 += +1.00 f(e,a) R1(a)  // flops: o0v1L1 += o0v2L1 | mem: o0v1L1 += o0v1L1
    sigmaR1("e") += this->fockMatrix_ta["vv"]("e,a") * R1("a");
    
    // sigmaR1 += +1.00 <j,i||a,b> this->T1_(a,i) R2(b,e,j)  // flops: o0v1L1 += o1v2L1 | mem: o0v1L1 += o0v1L1
    sigmaR1("e") += reuse_tmps_["vo_14"]("b,j") * R2("b,e,j");
    
    // sigmaR1 += -1.00 f(i,a) R2(a,e,i)  // flops: o0v1L1 += o1v2L1 | mem: o0v1L1 += o0v1L1
    sigmaR1("e") -= this->fockMatrix_ta["ov"]("i,a") * R2("a,e,i");
    
    // sigmaR2 = -1.00 P(e,f) f(e,m) R1(f)  // flops: o1v2L1 = o1v2L1 | mem: o1v2L1 = o1v2L1
    sigmaR2("e,f,m")  = -1.00 * this->fockMatrix_ta["vo"]("e,m") * R1("f");
    sigmaR2("e,f,m") += reuse_tmps_["vo_28"]("f,m") * R1("e");
    sigmaR2("e,f,m") += reuse_tmps_["vo_29"]("f,m") * R1("e");
    sigmaR2("e,f,m") += reuse_tmps_["vo_27"]("f,m") * R1("e");
    
    // sigmaR2 += +0.25 <j,i||a,b> this->T2_(a,b,j,i) R2(e,f,m)  // flops: o1v2L1 += o1v2L1 | mem: o1v2L1 += o1v2L1
    
    // sigmaR2 += +0.50 P(e,f) <i,e||a,b> this->T2_(a,b,m,i) R1(f)  // flops: o1v2L1 += o1v2L1 | mem: o1v2L1 += o1v2L1
    sigmaR2("e,f,m") -= reuse_tmps_["vo_11"]("e,m") * R1("f");
    
    // sigmaR2 += +1.00 P(e,f) f(i,m) this->T1_(e,i) R1(f)  // flops: o1v2L1 += o1v2L1 | mem: o1v2L1 += o1v2L1
    sigmaR2("e,f,m") += reuse_tmps_["vo_24"]("e,m") * R1("f");
    
    // sigmaR2 += +1.00 f(i,i) R2(e,f,m)  // flops: o1v2L1 += o1v2L1 | mem: o1v2L1 += o1v2L1
    
    // sigmaR2 += -1.00 P(e,f) <i,e||a,b> this->T1_(a,i) this->T1_(b,m) R1(f)  // flops: o1v2L1 += o1v2L1 | mem: o1v2L1 += o1v2L1
    sigmaR2("e,f,m") += reuse_tmps_["vo_30"]("e,m") * R1("f");
    sigmaR2("e,f,m") -= reuse_tmps_["vo_18"]("f,m") * R1("e");
    sigmaR2("e,f,m") += reuse_tmps_["vo_33"]("f,m") * R1("e");
    sigmaR2("e,f,m") -= reuse_tmps_["vo_22"]("f,m") * R1("e");
    
    // sigmaR2 += +0.50 P(e,f) <j,i||a,m> this->T2_(a,e,j,i) R1(f)  // flops: o1v2L1 += o1v2L1 | mem: o1v2L1 += o1v2L1
    sigmaR2("e,f,m") += reuse_tmps_["vo_18"]("e,m") * R1("f");
    
    // sigmaR2 += -0.50 <j,i||j,i> R2(e,f,m)  // flops: o1v2L1 += o1v2L1 | mem: o1v2L1 += o1v2L1
    sigmaR2("e,f,m") += this->fockMatrix_ta["vo"]("f,m") * R1("e");
    sigmaR2("e,f,m") += reuse_tmps_["vo_11"]("f,m") * R1("e");
    
    // sigmaR2 += -1.00 P(e,f) <i,e||a,m> this->T1_(a,i) R1(f)  // flops: o1v2L1 += o1v2L1 | mem: o1v2L1 += o1v2L1
    sigmaR2("e,f,m") += reuse_tmps_["vo_21"]("e,m") * R1("f");
    
    // sigmaR2 += +1.00 f(i,a) this->T1_(a,i) R2(e,f,m)  // flops: o1v2L1 += o1v2L1 | mem: o1v2L1 += o1v2L1
    sigmaR2("e,f,m") -= reuse_tmps_["vo_24"]("f,m") * R1("e");
    sigmaR2("e,f,m") -= reuse_tmps_["vo_21"]("f,m") * R1("e");
    sigmaR2("e,f,m") -= reuse_tmps_["vo_30"]("f,m") * R1("e");
    sigmaR2("e,f,m") += reuse_tmps_["vo_23"]("f,m") * R1("e");
    sigmaR2("e,f,m") -= reuse_tmps_["vo_31"]("f,m") * R1("e");
    
    // sigmaR2 += +1.00 P(e,f) f(i,a) this->T1_(a,m) this->T1_(e,i) R1(f)  // flops: o1v2L1 += o1v2L1 | mem: o1v2L1 += o1v2L1
    sigmaR2("e,f,m") += reuse_tmps_["vo_31"]("e,m") * R1("f");
    
    // sigmaR2 += -1.00 P(e,f) <j,i||a,m> this->T1_(a,i) this->T1_(e,j) R1(f)  // flops: o1v2L1 += o1v2L1 | mem: o1v2L1 += o1v2L1
    sigmaR2("e,f,m") -= reuse_tmps_["vo_32"]("e,m") * R1("f");
    
    // sigmaR2 += +1.00 P(e,f) f(i,a) this->T2_(a,e,m,i) R1(f)  // flops: o1v2L1 += o1v2L1 | mem: o1v2L1 += o1v2L1
    sigmaR2("e,f,m") += reuse_tmps_["vo_22"]("e,m") * R1("f");
    
    // sigmaR2 += -0.50 P(e,f) <j,i||a,b> this->T1_(e,i) this->T2_(a,b,m,j) R1(f)  // flops: o1v2L1 += o1v2L1 | mem: o1v2L1 += o1v2L1
    sigmaR2("e,f,m") -= reuse_tmps_["vo_29"]("e,m") * R1("f");
    
    // sigmaR2 += -0.50 P(e,f) <j,i||a,b> this->T1_(a,m) this->T2_(b,e,j,i) R1(f)  // flops: o1v2L1 += o1v2L1 | mem: o1v2L1 += o1v2L1
    sigmaR2("e,f,m") -= reuse_tmps_["vo_28"]("e,m") * R1("f");
    
    // sigmaR2 += -1.00 P(e,f) <j,i||a,b> this->T1_(a,i) this->T2_(b,e,m,j) R1(f)  // flops: o1v2L1 += o1v2L1 | mem: o1v2L1 += o1v2L1
    sigmaR2("e,f,m") -= reuse_tmps_["vo_27"]("e,m") * R1("f");
    
    // sigmaR2 += -0.50 <j,i||a,b> this->T1_(a,i) this->T1_(b,j) R2(e,f,m)  // flops: o1v2L1 += o1v2L1 | mem: o1v2L1 += o1v2L1
    
    // sigmaR2 += -1.00 P(e,f) <j,i||a,b> this->T1_(a,i) this->T1_(b,m) this->T1_(e,j) R1(f)  // flops: o1v2L1 += o1v2L1 | mem: o1v2L1 += o1v2L1
    sigmaR2("e,f,m") -= reuse_tmps_["vo_33"]("e,m") * R1("f");
    sigmaR2("e,f,m") += reuse_tmps_["vo_32"]("f,m") * R1("e");
    
    // sigmaR2 += -1.00 P(e,f) f(e,a) this->T1_(a,m) R1(f)  // flops: o1v2L1 += o1v2L1 | mem: o1v2L1 += o1v2L1
    sigmaR2("e,f,m") -= reuse_tmps_["vo_23"]("e,m") * R1("f");
    
    // sigmaR2 += +1.00 <j,i||a,m> this->T1_(a,i) R2(e,f,j)  // flops: o1v2L1 += o2v2L1 | mem: o1v2L1 += o1v2L1
    sigmaR2("e,f,m") += reuse_tmps_["oo_15"]("m,j") * R2("e,f,j");
    
    // sigmaR2 += -1.00 f(i,a) this->T1_(a,m) R2(e,f,i)  // flops: o1v2L1 += o2v2L1 | mem: o1v2L1 += o1v2L1
    sigmaR2("e,f,m") -= reuse_tmps_["oo_20"]("m,i") * R2("e,f,i");
    
    // sigmaR2 += +1.00 <j,i||a,b> this->T1_(a,i) this->T1_(b,m) R2(e,f,j)  // flops: o1v2L1 += o2v2L1 | mem: o1v2L1 += o1v2L1
    sigmaR2("e,f,m") += reuse_tmps_["oo_25"]("j,m") * R2("e,f,j");
    
    // sigmaR2 += -1.00 f(i,m) R2(e,f,i)  // flops: o1v2L1 += o2v2L1 | mem: o1v2L1 += o1v2L1
    sigmaR2("e,f,m") -= this->fockMatrix_ta["oo"]("i,m") * R2("e,f,i");
    
    // sigmaR2 += -0.50 <j,i||a,b> this->T2_(a,b,m,i) R2(e,f,j)  // flops: o1v2L1 += o2v2L1 | mem: o1v2L1 += o1v2L1
    sigmaR2("e,f,m") -= 0.50 * reuse_tmps_["oo_10"]("m,j") * R2("e,f,j");
    
    // sigmaR1 += -0.50 <i,e||a,b> R2(a,b,i)  // flops: o0v1L1 += o1v3L1 | mem: o0v1L1 += o0v1L1
    sigmaR1("e") += 0.50 * conj(this->antiSymMoints["vvvo"]("a,b,e,i")) * R2("a,b,i");
    
    // sigmaR2 += +0.50 <j,i||a,m> this->T2_(e,f,j,i) R1(a)  // flops: o1v2L1 += o1v3L1 | mem: o1v2L1 += o1v2L1
    sigmaR2("e,f,m") += 0.50 * reuse_tmps_["vvvo_17"]("f,e,a,m") * R1("a");
    
    // sigmaR2 += -1.00 <e,f||a,b> this->T1_(a,m) R1(b)  // flops: o1v2L1 += o1v3L1 | mem: o1v2L1 += o1v2L1
    sigmaR2("e,f,m") -= reuse_tmps_["vvvo_19"]("b,f,e,m") * R1("b");
    
    // sigmaR2 += -0.50 <j,i||a,b> this->T1_(a,m) this->T2_(e,f,j,i) R1(b)  // flops: o1v2L1 += o1v3L1 | mem: o1v2L1 += o1v2L1
    sigmaR2("e,f,m") -= 0.50 * reuse_tmps_["vvvo_26"]("b,e,f,m") * R1("b");
    
    // sigmaR2 += +1.00 <e,f||a,m> R1(a)  // flops: o1v2L1 += o1v3L1 | mem: o1v2L1 += o1v2L1
    sigmaR2("e,f,m") += this->antiSymMoints["vvvo"]("e,f,a,m") * R1("a");
    
    // sigmaR2 += +0.50 <e,f||a,b> R2(a,b,m)  // flops: o1v2L1 += o1v4L1 | mem: o1v2L1 += o1v2L1
    sigmaR2("e,f,m") += 0.50 * this->antiSymMoints["vvvv"]("e,f,a,b") * R2("a,b,m");
    
    // sigmaR2 += +1.00 <j,i||a,b> this->T1_(a,m) this->T1_(e,i) this->T1_(f,j) R1(b)  // flops: o1v2L1 += o3v1L1 o3v1L1 o2v2L1 | mem: o1v2L1 += o3v0L1 o2v1L1 o1v2L1
    sigmaR2("e,f,m") += reuse_tmps_["vooo_12"]("b,m,i,j") * R1("b") * this->T1_("e,i") * this->T1_("f,j");
    
    // sigmaR2 += -1.00 <j,i||a,m> this->T1_(e,i) this->T1_(f,j) R1(a)  // flops: o1v2L1 += o3v1L1 o3v1L1 o2v2L1 | mem: o1v2L1 += o3v0L1 o2v1L1 o1v2L1
    sigmaR2("e,f,m") -= conj(this->antiSymMoints["vooo"]("a,m,j,i")) * R1("a") * this->T1_("e,i") * this->T1_("f,j");
    
    // tmps_[1_Lvoo](e,m,i) = 0.50 eri[vovv](e,i,a,b) * R2(a,b,m) // flops: o2v1L1 = o2v3L1 | mem: o2v1L1 = o2v1L1
    tmps_.emplace(std::make_pair("voo_1", TAmanager.malloc<MatsT>("voo")));
    tmps_["voo_1"]("e,m,i")  = 0.50 * conj(this->antiSymMoints["vvvo"]("a,b,e,i")) * R2("a,b,m");
    
    // tmps_[28_Lvvo](e,f,m) = 1.00 this->T1_(f,i) * eri[vovv](e,i,a,b) * R2(a,b,m) // flops: o1v2L1 = o2v2L1 | mem: o1v2L1 = o1v2L1
    tmps_.emplace(std::make_pair("vvo_28", TAmanager.malloc<MatsT>("vvo")));
    tmps_["vvo_28"]("e,f,m")  = this->T1_("f,i") * tmps_["voo_1"]("e,m,i");
    TAmanager.free("voo", std::move(tmps_["voo_1"]));
    sigmaR2("e,f,m") += tmps_["vvo_28"]("f,e,m");
    
    // sigmaR2 += +0.50 P(e,f) <i,e||a,b> this->T1_(f,i) R2(a,b,m)  // flops: o1v2L1 += o1v2L1 | mem: o1v2L1 += o1v2L1
    sigmaR2("e,f,m") -= tmps_["vvo_28"]("e,f,m");
    TAmanager.free("vvo", std::move(tmps_["vvo_28"]));
    
    // tmps_[2_Lvvo](f,e,m) = 1.00 eri[vovo](e,i,a,m) * R2(a,f,i) // flops: o1v2L1 = o2v3L1 | mem: o1v2L1 = o1v2L1
    tmps_.emplace(std::make_pair("vvo_2", TAmanager.malloc<MatsT>("vvo")));
    tmps_["vvo_2"]("f,e,m")  = this->antiSymMoints["vovo"]("e,i,a,m") * R2("a,f,i");
    sigmaR2("e,f,m") += tmps_["vvo_2"]("e,f,m");
    
    // sigmaR2 += +1.00 P(e,f) <i,e||a,m> R2(a,f,i)  // flops: o1v2L1 += o1v2L1 | mem: o1v2L1 += o1v2L1
    sigmaR2("e,f,m") -= tmps_["vvo_2"]("f,e,m");
    TAmanager.free("vvo", std::move(tmps_["vvo_2"]));
    
    // tmps_[3_Lvvo](e,f,m) = 1.00 eri[vovv](f,i,a,b) * this->T1_(a,m) * R2(b,e,i) // flops: o1v2L1 = o2v3L1 | mem: o1v2L1 = o1v2L1
    tmps_.emplace(std::make_pair("vvo_3", TAmanager.malloc<MatsT>("vvo")));
    tmps_["vvo_3"]("e,f,m")  = reuse_tmps_["vvoo_1"]("b,f,m,i") * R2("b,e,i");
    sigmaR2("e,f,m") -= tmps_["vvo_3"]("e,f,m");
    
    // sigmaR2 += -1.00 P(e,f) <i,e||a,b> this->T1_(a,m) R2(b,f,i)  // flops: o1v2L1 += o1v2L1 | mem: o1v2L1 += o1v2L1
    sigmaR2("e,f,m") += tmps_["vvo_3"]("f,e,m");
    TAmanager.free("vvo", std::move(tmps_["vvo_3"]));
    
    // tmps_[4_Lvvo](e,f,m) = 1.00 R2(b,f,j) * eri[oovv](j,i,a,b) * this->T2_(a,e,m,i) // flops: o1v2L1 = o2v3L1 | mem: o1v2L1 = o1v2L1
    tmps_.emplace(std::make_pair("vvo_4", TAmanager.malloc<MatsT>("vvo")));
    tmps_["vvo_4"]("e,f,m")  = R2("b,f,j") * reuse_tmps_["vvoo_5"]("e,b,m,j");
    
    // sigmaR2 += +1.00 P(e,f) <j,i||a,b> this->T2_(a,e,m,i) R2(b,f,j)  // flops: o1v2L1 += o1v2L1 | mem: o1v2L1 += o1v2L1
    sigmaR2("e,f,m") += tmps_["vvo_4"]("e,f,m");
    sigmaR2("e,f,m") -= tmps_["vvo_4"]("f,e,m");
    TAmanager.free("vvo", std::move(tmps_["vvo_4"]));
    
    // tmps_[5_vvvo](f,b,e,m) = 1.00 eri[vovv](e,i,a,b) * this->T2_(a,f,m,i) // flops: o1v3 = o2v4 | mem: o1v3 = o1v3
    tmps_.emplace(std::make_pair("vvvo_5", TAmanager.malloc<MatsT>("vvvo")));
    tmps_["vvvo_5"]("f,b,e,m")  = conj(this->antiSymMoints["vvvo"]("a,b,e,i")) * this->T2_("a,f,m,i");
    
    // tmps_[20_Lvvo](f,e,m) = 1.00 R1(b) * eri[vovv](f,i,a,b) * this->T2_(a,e,m,i) // flops: o1v2L1 = o1v3L1 | mem: o1v2L1 = o1v2L1
    tmps_.emplace(std::make_pair("vvo_20", TAmanager.malloc<MatsT>("vvo")));
    tmps_["vvo_20"]("f,e,m")  = R1("b") * tmps_["vvvo_5"]("e,b,f,m");
    TAmanager.free("vvvo", std::move(tmps_["vvvo_5"]));
    sigmaR2("e,f,m") -= tmps_["vvo_20"]("f,e,m");
    
    // sigmaR2 += -1.00 P(e,f) <i,e||a,b> this->T2_(a,f,m,i) R1(b)  // flops: o1v2L1 += o1v2L1 | mem: o1v2L1 += o1v2L1
    sigmaR2("e,f,m") += tmps_["vvo_20"]("e,f,m");
    TAmanager.free("vvo", std::move(tmps_["vvo_20"]));
    
    // tmps_[6_Looo](m,i,j) = 1.00 eri[oovv](j,i,a,b) * R2(a,b,m) // flops: o3v0L1 = o3v2L1 | mem: o3v0L1 = o3v0L1
    tmps_.emplace(std::make_pair("ooo_6", TAmanager.malloc<MatsT>("ooo")));
    tmps_["ooo_6"]("m,i,j")  = conj(this->antiSymMoints["vvoo"]("a,b,j,i")) * R2("a,b,m");
    
    // sigmaR2 += +0.25 <j,i||a,b> this->T2_(e,f,j,i) R2(a,b,m)  // flops: o1v2L1 += o3v2L1 | mem: o1v2L1 += o1v2L1
    sigmaR2("e,f,m") += 0.25 * tmps_["ooo_6"]("m,i,j") * this->T2_("e,f,j,i");
    
    // sigmaR2 += -0.50 <j,i||a,b> this->T1_(e,i) this->T1_(f,j) R2(a,b,m)  // flops: o1v2L1 += o3v1L1 o2v2L1 | mem: o1v2L1 += o2v1L1 o1v2L1
    sigmaR2("e,f,m") -= 0.50 * tmps_["ooo_6"]("m,i,j") * this->T1_("e,i") * this->T1_("f,j");
    TAmanager.free("ooo", std::move(tmps_["ooo_6"]));
    
    // tmps_[7_Lvoo](f,m,i) = 1.00 eri[oovo](j,i,a,m) * R2(a,f,j) // flops: o2v1L1 = o3v2L1 | mem: o2v1L1 = o2v1L1
    tmps_.emplace(std::make_pair("voo_7", TAmanager.malloc<MatsT>("voo")));
    tmps_["voo_7"]("f,m,i")  = conj(this->antiSymMoints["vooo"]("a,m,j,i")) * R2("a,f,j");
    
    // tmps_[25_Lvvo](e,f,m) = 1.00 this->T1_(f,i) * eri[oovo](j,i,a,m) * R2(a,e,j) // flops: o1v2L1 = o2v2L1 | mem: o1v2L1 = o1v2L1
    tmps_.emplace(std::make_pair("vvo_25", TAmanager.malloc<MatsT>("vvo")));
    tmps_["vvo_25"]("e,f,m")  = this->T1_("f,i") * tmps_["voo_7"]("e,m,i");
    TAmanager.free("voo", std::move(tmps_["voo_7"]));
    
    // sigmaR2 += -1.00 P(e,f) <j,i||a,m> this->T1_(e,i) R2(a,f,j)  // flops: o1v2L1 += o1v2L1 | mem: o1v2L1 += o1v2L1
    sigmaR2("e,f,m") -= tmps_["vvo_25"]("f,e,m");
    sigmaR2("e,f,m") += tmps_["vvo_25"]("e,f,m");
    TAmanager.free("vvo", std::move(tmps_["vvo_25"]));
    
    // tmps_[8_Lvoo](f,i,m) = 1.00 R2(b,f,j) * eri[oovv](j,i,a,b) * this->T1_(a,m) // flops: o2v1L1 = o3v2L1 | mem: o2v1L1 = o2v1L1
    tmps_.emplace(std::make_pair("voo_8", TAmanager.malloc<MatsT>("voo")));
    tmps_["voo_8"]("f,i,m")  = R2("b,f,j") * reuse_tmps_["vooo_12"]("b,m,i,j");
    
    // tmps_[22_Lvvo](e,f,m) = 1.00 this->T1_(f,i) * R2(b,e,j) * eri[oovv](j,i,a,b) * this->T1_(a,m) // flops: o1v2L1 = o2v2L1 | mem: o1v2L1 = o1v2L1
    tmps_.emplace(std::make_pair("vvo_22", TAmanager.malloc<MatsT>("vvo")));
    tmps_["vvo_22"]("e,f,m")  = this->T1_("f,i") * tmps_["voo_8"]("e,i,m");
    TAmanager.free("voo", std::move(tmps_["voo_8"]));
    sigmaR2("e,f,m") -= tmps_["vvo_22"]("e,f,m");
    
    // sigmaR2 += +1.00 P(e,f) <j,i||a,b> this->T1_(a,m) this->T1_(e,i) R2(b,f,j)  // flops: o1v2L1 += o1v2L1 | mem: o1v2L1 += o1v2L1
    sigmaR2("e,f,m") += tmps_["vvo_22"]("f,e,m");
    TAmanager.free("vvo", std::move(tmps_["vvo_22"]));
    
    // tmps_[9_Lvvo](e,f,m) = 1.00 f[vv](f,a) * R2(a,e,m) // flops: o1v2L1 = o1v3L1 | mem: o1v2L1 = o1v2L1
    tmps_.emplace(std::make_pair("vvo_9", TAmanager.malloc<MatsT>("vvo")));
    tmps_["vvo_9"]("e,f,m")  = this->fockMatrix_ta["vv"]("f,a") * R2("a,e,m");
    sigmaR2("e,f,m") -= tmps_["vvo_9"]("e,f,m");
    
    // sigmaR2 += +1.00 P(e,f) f(e,a) R2(a,f,m)  // flops: o1v2L1 += o1v2L1 | mem: o1v2L1 += o1v2L1
    sigmaR2("e,f,m") += tmps_["vvo_9"]("f,e,m");
    TAmanager.free("vvo", std::move(tmps_["vvo_9"]));
    
    // tmps_[10_Lvvo](e,f,m) = 1.00 R2(b,f,m) * eri[vovv](e,i,a,b) * this->T1_(a,i) // flops: o1v2L1 = o1v3L1 | mem: o1v2L1 = o1v2L1
    tmps_.emplace(std::make_pair("vvo_10", TAmanager.malloc<MatsT>("vvo")));
    tmps_["vvo_10"]("e,f,m")  = R2("b,f,m") * reuse_tmps_["vv_3"]("b,e");
    sigmaR2("e,f,m") += tmps_["vvo_10"]("f,e,m");
    
    // sigmaR2 += +1.00 P(e,f) <i,e||a,b> this->T1_(a,i) R2(b,f,m)  // flops: o1v2L1 += o1v2L1 | mem: o1v2L1 += o1v2L1
    sigmaR2("e,f,m") -= tmps_["vvo_10"]("e,f,m");
    TAmanager.free("vvo", std::move(tmps_["vvo_10"]));
    
    // tmps_[11_Lvvo](e,f,m) = 1.00 R2(b,f,m) * eri[oovv](j,i,a,b) * this->T2_(a,e,j,i) // flops: o1v2L1 = o1v3L1 | mem: o1v2L1 = o1v2L1
    tmps_.emplace(std::make_pair("vvo_11", TAmanager.malloc<MatsT>("vvo")));
    tmps_["vvo_11"]("e,f,m")  = R2("b,f,m") * reuse_tmps_["vv_9"]("e,b");
    
    // sigmaR2 += -0.50 P(e,f) <j,i||a,b> this->T2_(a,e,j,i) R2(b,f,m)  // flops: o1v2L1 += o1v2L1 | mem: o1v2L1 += o1v2L1
    sigmaR2("e,f,m") -= tmps_["vvo_11"]("e,f,m");
    sigmaR2("e,f,m") += tmps_["vvo_11"]("f,e,m");
    TAmanager.free("vvo", std::move(tmps_["vvo_11"]));
    
    // tmps_[12_Lo](i) = 0.50 eri[oovv](j,i,a,b) * R2(a,b,j) // flops: o1v0L1 = o2v2L1 | mem: o1v0L1 = o1v0L1
    tmps_.emplace(std::make_pair("o_12", TAmanager.malloc<MatsT>("o")));
    tmps_["o_12"]("i")  = 0.50 * conj(this->antiSymMoints["vvoo"]("a,b,j,i")) * R2("a,b,j");
    
    // sigmaR1 += +0.50 <j,i||a,b> this->T1_(e,i) R2(a,b,j)  // flops: o0v1L1 += o1v1L1 | mem: o0v1L1 += o0v1L1
    sigmaR1("e") += tmps_["o_12"]("i") * this->T1_("e,i");
    
    // sigmaR2 += -0.50 <j,i||a,b> this->T2_(e,f,m,i) R2(a,b,j)  // flops: o1v2L1 += o2v2L1 | mem: o1v2L1 += o1v2L1
    sigmaR2("e,f,m") -= tmps_["o_12"]("i") * this->T2_("e,f,m,i");
    TAmanager.free("o", std::move(tmps_["o_12"]));
    
    // tmps_[13_Lvoo](e,i,m) = 1.00 R1(b) * eri[oovv](j,i,a,b) * this->T2_(a,e,m,j) // flops: o2v1L1 = o2v2L1 | mem: o2v1L1 = o2v1L1
    tmps_.emplace(std::make_pair("voo_13", TAmanager.malloc<MatsT>("voo")));
    tmps_["voo_13"]("e,i,m")  = R1("b") * reuse_tmps_["vvoo_6"]("e,b,m,i");
    
    // tmps_[27_Lvvo](e,f,m) = 1.00 this->T1_(f,i) * R1(b) * eri[oovv](j,i,a,b) * this->T2_(a,e,m,j) // flops: o1v2L1 = o2v2L1 | mem: o1v2L1 = o1v2L1
    tmps_.emplace(std::make_pair("vvo_27", TAmanager.malloc<MatsT>("vvo")));
    tmps_["vvo_27"]("e,f,m")  = this->T1_("f,i") * tmps_["voo_13"]("e,i,m");
    TAmanager.free("voo", std::move(tmps_["voo_13"]));
    
    // sigmaR2 += +1.00 P(e,f) <j,i||a,b> this->T1_(e,i) this->T2_(a,f,m,j) R1(b)  // flops: o1v2L1 += o1v2L1 | mem: o1v2L1 += o1v2L1
    sigmaR2("e,f,m") += tmps_["vvo_27"]("f,e,m");
    sigmaR2("e,f,m") -= tmps_["vvo_27"]("e,f,m");
    TAmanager.free("vvo", std::move(tmps_["vvo_27"]));
    
    // tmps_[14_Lvoo](e,m,i) = 1.00 eri[vovo](e,i,a,m) * R1(a) // flops: o2v1L1 = o2v2L1 | mem: o2v1L1 = o2v1L1
    tmps_.emplace(std::make_pair("voo_14", TAmanager.malloc<MatsT>("voo")));
    tmps_["voo_14"]("e,m,i")  = this->antiSymMoints["vovo"]("e,i,a,m") * R1("a");
    
    // tmps_[26_Lvvo](e,f,m) = 1.00 this->T1_(f,i) * eri[vovo](e,i,a,m) * R1(a) // flops: o1v2L1 = o2v2L1 | mem: o1v2L1 = o1v2L1
    tmps_.emplace(std::make_pair("vvo_26", TAmanager.malloc<MatsT>("vvo")));
    tmps_["vvo_26"]("e,f,m")  = this->T1_("f,i") * tmps_["voo_14"]("e,m,i");
    TAmanager.free("voo", std::move(tmps_["voo_14"]));
    sigmaR2("e,f,m") += tmps_["vvo_26"]("f,e,m");
    
    // sigmaR2 += +1.00 P(e,f) <i,e||a,m> this->T1_(f,i) R1(a)  // flops: o1v2L1 += o1v2L1 | mem: o1v2L1 += o1v2L1
    sigmaR2("e,f,m") -= tmps_["vvo_26"]("e,f,m");
    TAmanager.free("vvo", std::move(tmps_["vvo_26"]));
    
    // tmps_[15_Lvoo](f,m,i) = 1.00 f[ov](i,a) * R2(a,f,m) // flops: o2v1L1 = o2v2L1 | mem: o2v1L1 = o2v1L1
    tmps_.emplace(std::make_pair("voo_15", TAmanager.malloc<MatsT>("voo")));
    tmps_["voo_15"]("f,m,i")  = this->fockMatrix_ta["ov"]("i,a") * R2("a,f,m");
    
    // tmps_[24_Lvvo](e,f,m) = 1.00 this->T1_(f,i) * f[ov](i,a) * R2(a,e,m) // flops: o1v2L1 = o2v2L1 | mem: o1v2L1 = o1v2L1
    tmps_.emplace(std::make_pair("vvo_24", TAmanager.malloc<MatsT>("vvo")));
    tmps_["vvo_24"]("e,f,m")  = this->T1_("f,i") * tmps_["voo_15"]("e,m,i");
    TAmanager.free("voo", std::move(tmps_["voo_15"]));
    sigmaR2("e,f,m") += tmps_["vvo_24"]("e,f,m");
    
    // sigmaR2 += -1.00 P(e,f) f(i,a) this->T1_(e,i) R2(a,f,m)  // flops: o1v2L1 += o1v2L1 | mem: o1v2L1 += o1v2L1
    sigmaR2("e,f,m") -= tmps_["vvo_24"]("f,e,m");
    TAmanager.free("vvo", std::move(tmps_["vvo_24"]));
    
    // tmps_[16_Lvoo](f,j,m) = 1.00 R2(b,f,m) * eri[oovv](j,i,a,b) * this->T1_(a,i) // flops: o2v1L1 = o2v2L1 | mem: o2v1L1 = o2v1L1
    tmps_.emplace(std::make_pair("voo_16", TAmanager.malloc<MatsT>("voo")));
    tmps_["voo_16"]("f,j,m")  = R2("b,f,m") * reuse_tmps_["vo_14"]("b,j");
    
    // tmps_[21_Lvvo](e,f,m) = 1.00 this->T1_(f,j) * R2(b,e,m) * eri[oovv](j,i,a,b) * this->T1_(a,i) // flops: o1v2L1 = o2v2L1 | mem: o1v2L1 = o1v2L1
    tmps_.emplace(std::make_pair("vvo_21", TAmanager.malloc<MatsT>("vvo")));
    tmps_["vvo_21"]("e,f,m")  = this->T1_("f,j") * tmps_["voo_16"]("e,j,m");
    TAmanager.free("voo", std::move(tmps_["voo_16"]));
    sigmaR2("e,f,m") -= tmps_["vvo_21"]("e,f,m");
    
    // sigmaR2 += +1.00 P(e,f) <j,i||a,b> this->T1_(a,i) this->T1_(e,j) R2(b,f,m)  // flops: o1v2L1 += o1v2L1 | mem: o1v2L1 += o1v2L1
    sigmaR2("e,f,m") += tmps_["vvo_21"]("f,e,m");
    TAmanager.free("vvo", std::move(tmps_["vvo_21"]));
    
    // tmps_[17_Lvoo](f,i,m) = 1.00 R1(b) * eri[vovv](f,i,a,b) * this->T1_(a,m) // flops: o2v1L1 = o2v2L1 | mem: o2v1L1 = o2v1L1
    tmps_.emplace(std::make_pair("voo_17", TAmanager.malloc<MatsT>("voo")));
    tmps_["voo_17"]("f,i,m")  = R1("b") * reuse_tmps_["vvoo_1"]("b,f,m,i");
    
    // tmps_[23_Lvvo](f,e,m) = 1.00 R1(b) * eri[vovv](e,i,a,b) * this->T1_(a,m) * this->T1_(f,i) // flops: o1v2L1 = o2v2L1 | mem: o1v2L1 = o1v2L1
    tmps_.emplace(std::make_pair("vvo_23", TAmanager.malloc<MatsT>("vvo")));
    tmps_["vvo_23"]("f,e,m")  = tmps_["voo_17"]("e,i,m") * this->T1_("f,i");
    TAmanager.free("voo", std::move(tmps_["voo_17"]));
    sigmaR2("e,f,m") -= tmps_["vvo_23"]("e,f,m");
    
    // sigmaR2 += -1.00 P(e,f) <i,e||a,b> this->T1_(a,m) this->T1_(f,i) R1(b)  // flops: o1v2L1 += o1v2L1 | mem: o1v2L1 += o1v2L1
    sigmaR2("e,f,m") += tmps_["vvo_23"]("f,e,m");
    TAmanager.free("vvo", std::move(tmps_["vvo_23"]));
    
    // tmps_[18_Lo](i) = 1.00 f[ov](i,a) * R1(a) // flops: o1v0L1 = o1v1L1 | mem: o1v0L1 = o1v0L1
    tmps_.emplace(std::make_pair("o_18", TAmanager.malloc<MatsT>("o")));
    tmps_["o_18"]("i")  = this->fockMatrix_ta["ov"]("i,a") * R1("a");
    
    // sigmaR1 += -1.00 f(i,a) this->T1_(e,i) R1(a)  // flops: o0v1L1 += o1v1L1 | mem: o0v1L1 += o0v1L1
    sigmaR1("e") -= tmps_["o_18"]("i") * this->T1_("e,i");
    
    // sigmaR2 += +1.00 f(i,a) this->T2_(e,f,m,i) R1(a)  // flops: o1v2L1 += o2v2L1 | mem: o1v2L1 += o1v2L1
    sigmaR2("e,f,m") += tmps_["o_18"]("i") * this->T2_("e,f,m,i");
    TAmanager.free("o", std::move(tmps_["o_18"]));
    
    // tmps_[19_Lo](j) = 1.00 R1(b) * eri[oovv](j,i,a,b) * this->T1_(a,i) // flops: o1v0L1 = o1v1L1 | mem: o1v0L1 = o1v0L1
    tmps_.emplace(std::make_pair("o_19", TAmanager.malloc<MatsT>("o")));
    tmps_["o_19"]("j")  = R1("b") * reuse_tmps_["vo_14"]("b,j");
    
    // sigmaR1 += +1.00 <j,i||a,b> this->T1_(a,i) this->T1_(e,j) R1(b)  // flops: o0v1L1 += o1v1L1 | mem: o0v1L1 += o0v1L1
    sigmaR1("e") += tmps_["o_19"]("j") * this->T1_("e,j");
    
    // sigmaR2 += -1.00 <j,i||a,b> this->T1_(a,i) this->T2_(e,f,m,j) R1(b)  // flops: o1v2L1 += o2v2L1 | mem: o1v2L1 += o1v2L1
    sigmaR2("e,f,m") -= tmps_["o_19"]("j") * this->T2_("e,f,m,j");
    TAmanager.free("o", std::move(tmps_["o_19"]));
    

  }

  template <typename MatsT>
  void EOMEA<MatsT>::buildDiag(MatsT * diag, const std::vector<double> &eps) const {

    TAManager &TAmanager = TAManager::get();
    size_t n_v = TAmanager.getRange(vLabel_).extent();
    size_t n_o = TAmanager.getRange(oLabel_).extent();
    std::fill_n(diag, this->Hbar_dim, MatsT(0.0));

    MatsT * diag_a = diag;
    for (auto a = 0; a < n_v; a++){
      diag_a[toCompoundS(a)] = eps[a+n_o];
    }
    MatsT * diag_abi = diag + this->Hbar_dimension_offsets.at("TwoBody");
    for (auto a = 0; a < n_v; a++){
      for (auto b = 0; b < a; b++){
        for (auto i = 0; i < n_o; i++){
          diag_abi[toCompoundD(a,b,i)] = eps[a+n_o] + eps[b+n_o] - eps[i];
        }
      }
    }
  }



  template <typename MatsT>
  void EOMEA<MatsT>::runLambda(){} 


  template <typename MatsT>
  typename Davidson<MatsT>::VecsGen_t EOMEA<MatsT>::EmptyDavidsonVectorBuilder(){
      // Algorithm with implicit Hbar matrix
      typename Davidson<MatsT>::VecsGen_t vecsGenEOM;
      if (this->eomSettings.hbar_type == EOM_HBAR_TYPE::IMPLICIT) {
        vecsGenEOM = [this](size_t nVec)->std::shared_ptr<SolverVectors<MatsT>> {
          return std::make_shared<MBExpansionSet<MatsT>>(this->tensor_builder_, nVec, this->savFile_);
        }; // implicit vecsGenerator

        return vecsGenEOM;
      }
      return vecsGenEOM;
  }

  template <typename MatsT>
  typename Davidson<MatsT>::LinearTrans_t EOMEA<MatsT>::DavidsonResidualBuilder(EOMCCEigenVecType &eigenVecType){
      if (this->eomSettings.hbar_type == EOM_HBAR_TYPE::IMPLICIT) {
        this->funcEOM = [this, &eigenVecType]( size_t nVec, SolverVectors<MatsT> &V,
            SolverVectors<MatsT> &AV) {

          MBExpansionSet<MatsT> *V_ptr = nullptr, *AV_ptr = nullptr;
          size_t Vshift = 0, AVshift = 0;
          try {
            V_ptr = &dynamic_cast<MBExpansionSet<MatsT>&>(V);
          } catch(const std::bad_cast& e) {
            SolverVectorsView<MatsT>& V_view = dynamic_cast<SolverVectorsView<MatsT>&>(V);
            V_ptr = &dynamic_cast<MBExpansionSet<MatsT>&>(V_view.getVecs());
            Vshift = V_view.shift();
          }

          try {
            AV_ptr = &dynamic_cast<MBExpansionSet<MatsT>&>(AV);
          } catch(const std::bad_cast& e) {
            SolverVectorsView<MatsT>& AV_view = dynamic_cast<SolverVectorsView<MatsT>&>(AV);
            AV_ptr = &dynamic_cast<MBExpansionSet<MatsT>&>(AV_view.getVecs());
            AVshift = AV_view.shift();
          }

          for (size_t i = 0; i < nVec; i++) {
            const MBExpansion<MatsT> &Vi = V_ptr->get(i + Vshift);
            MBExpansion<MatsT> &AVi = AV_ptr->get(i + AVshift);
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
  typename Davidson<MatsT>::LinearTrans_t EOMEA<MatsT>::DavidsonPreconditionerBuilder(dcomplex * curEig, MatsT * eomDiag){

      double PCsmall = this->eomSettings.davidson_preCond_small;

      this->PCEOM = [this, eomDiag, curEig, PCsmall]( size_t nVec, SolverVectors<MatsT> &V,
          SolverVectors<MatsT> &AV) {

        AV.set_data(0, nVec, V, 0);

        MBExpansionSet<MatsT> *AV_ptr = nullptr;
        size_t AVshift = 0;

        try {
          AV_ptr = &dynamic_cast<MBExpansionSet<MatsT>&>(AV);
        } catch(const std::bad_cast& e) {
          SolverVectorsView<MatsT>& AV_view = dynamic_cast<SolverVectorsView<MatsT>&>(AV);
          AV_ptr = &dynamic_cast<MBExpansionSet<MatsT>&>(AV_view.getVecs());
          AVshift = AV_view.shift();
        }

        for (size_t iVec = 0; iVec < nVec; iVec++) {

          MBExpansion<MatsT> &curB = AV_ptr->get(iVec + AVshift);
          MatsT curEigI = 0.0;
          if constexpr (std::is_same_v<MatsT, double>) {
            curEigI = curEig[iVec].real();
          } else {
            curEigI = curEig[iVec];
          }

          MatsT * Diag_A = eomDiag;
          TA::foreach_inplace(curB.get_tensor("OneBody"), [iVec, curEigI, Diag_A, this, PCsmall](TA::Tensor<MatsT> &tile){
            const auto& lobound = tile.range().lobound();
            const auto& upbound = tile.range().upbound();

            MatsT denom = 0.0;
            std::vector<std::size_t> x{0};
            for(x[0] = lobound[0]; x[0] < upbound[0]; ++x[0]) {
              denom = curEigI - Diag_A[toCompoundS(x[0])];
              if (std::abs(denom) >= PCsmall) tile[x] /= denom;
            }
          });
          TA::get_default_world().gop.fence();

          MatsT * Diag_ABI = eomDiag + this->Hbar_dimension_offsets.at("TwoBody");
          TA::foreach_inplace(curB.get_tensor("TwoBody"), [iVec, curEigI, Diag_ABI, this, PCsmall](TA::Tensor<MatsT> &tile){
            const auto& lobound = tile.range().lobound();
            const auto& upbound = tile.range().upbound();

            MatsT denom = 0.0;
            std::vector<std::size_t> x{0, 0, 0};
            for(x[0] = lobound[0]; x[0] < upbound[0]; ++x[0]){
              size_t a = x[0];
              for(x[1] = lobound[1]; x[1] < upbound[1]; ++x[1]){
                if (x[0] == x[1])
                  continue;
                size_t b = x[1];
                for(x[2] = lobound[2]; x[2] < upbound[2]; ++x[2]){
                  size_t i = x[2];
                  denom = curEigI - Diag_ABI[toCompoundD(a, b, i)];
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
  EOMEA<MatsT>::~EOMEA() {

    TAManager &TAmanager = TAManager::get();

    TAmanager.free("vvoo",std::move(reuse_tmps_["vvoo_1"]), true); 
    TAmanager.free("vvoo",std::move(reuse_tmps_["vvoo_2"]), true); 
    TAmanager.free("vo",std::move(reuse_tmps_["vo_30"]), true); 
    TAmanager.free("vv",std::move(reuse_tmps_["vv_3"]), true); 
    TAmanager.free("vvvv",std::move(reuse_tmps_["vvvv_4"]), true); 
    TAmanager.free("vvvo",std::move(reuse_tmps_["vvvo_26"]), true); 
    TAmanager.free("vvoo",std::move(reuse_tmps_["vvoo_5"]), true); 
    TAmanager.free("vvoo",std::move(reuse_tmps_["vvoo_6"]), true); 
    TAmanager.free("vvoo",std::move(reuse_tmps_["vvoo_7"]), true); 
    TAmanager.free("vo",std::move(reuse_tmps_["vo_27"]), true); 
    TAmanager.free("vv",std::move(reuse_tmps_["vv_8"]), true); 
    TAmanager.free("vo",std::move(reuse_tmps_["vo_28"]), true); 
    TAmanager.free("vv",std::move(reuse_tmps_["vv_9"]), true); 
    TAmanager.free("oo",std::move(reuse_tmps_["oo_10"]), true); 
    TAmanager.free("vo",std::move(reuse_tmps_["vo_11"]), true); 
    TAmanager.free("vooo",std::move(reuse_tmps_["vooo_12"]), true); 
    TAmanager.free("vooo",std::move(reuse_tmps_["vooo_13"]), true); 
    TAmanager.free("oo",std::move(reuse_tmps_["oo_25"]), true); 
    TAmanager.free("vo",std::move(reuse_tmps_["vo_33"]), true); 
    TAmanager.free("vo",std::move(reuse_tmps_["vo_14"]), true); 
    TAmanager.free("oo",std::move(reuse_tmps_["oo_15"]), true); 
    TAmanager.free("vo",std::move(reuse_tmps_["vo_32"]), true); 
    TAmanager.free("oo",std::move(reuse_tmps_["oo_16"]), true); 
    TAmanager.free("vo",std::move(reuse_tmps_["vo_29"]), true); 
    TAmanager.free("vvvo",std::move(reuse_tmps_["vvvo_17"]), true); 
    TAmanager.free("vo",std::move(reuse_tmps_["vo_18"]), true); 
    TAmanager.free("vvvo",std::move(reuse_tmps_["vvvo_19"]), true); 
    TAmanager.free("oo",std::move(reuse_tmps_["oo_20"]), true); 
    TAmanager.free("vo",std::move(reuse_tmps_["vo_31"]), true); 
    TAmanager.free("vo",std::move(reuse_tmps_["vo_21"]), true); 
    TAmanager.free("vo",std::move(reuse_tmps_["vo_22"]), true); 
    TAmanager.free("vo",std::move(reuse_tmps_["vo_23"]), true); 
    TAmanager.free("vo",std::move(reuse_tmps_["vo_24"]), true); 
  }


}
