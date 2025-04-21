#pragma once 
#include <chronusq_sys.hpp>
#include <coupledcluster.hpp>

namespace ChronusQ{
  template <typename MatsT>
  EOMDIP_4h2pCCSDT<MatsT>::EOMDIP_4h2pCCSDT(const SafeFile &savFile,
    CCIntermediates<MatsT> &intermediates,
    const EOMSettings &eomSettings,
    const CoupledClusterSettings &ccSettings):
    EOMCCBase<MatsT>(savFile, intermediates,eomSettings,ccSettings),
    vLabel_(intermediates.vLabel), oLabel_(intermediates.oLabel),
    reused_(intermediates.sigmaOps),
    tmps_(intermediates.tempOps),
    T1_(intermediates.T->get_tensor("OneBody")), T2_(intermediates.T->get_tensor("TwoBody")), T3_(intermediates.T->get_tensor("ThreeBody")) {
    TAManager &TAmanager = TAManager::get();

    // without L, we do not need D
    if (eomSettings.oscillator_strength == false) {
      if(intermediates.D_ai)   TAmanager.free("vo", std::move(intermediates.D_ai), true);
      if(intermediates.D_abij) TAmanager.free("vvoo", std::move(intermediates.D_abij), true);
    }

    nV_ = TAmanager.getRange(vLabel_).extent();
    nO_ = TAmanager.getRange(oLabel_).extent();
    nV2shift_ = nV_ * (nV_ - 1) / 2;
    nO2shift_ = nO_ * (nO_ - 1) / 2;
    nO3shift_ = nO_ * (nO_ - 1) * (nO_ - 2) / 6;
    nO4shift_ = nO_ * (nO_ - 1) * (nO_ - 2) * (nO_ - 3) / 24;
    this->Hbar_dimension_offsets.emplace("TwoBody", 0);
    this->Hbar_dimension_offsets.emplace("ThreeBody", nO2shift_);
    this->Hbar_dimension_offsets.emplace("FourBody", nO2shift_ + nO3shift_*nV_);
    this->Hbar_dim = nO2shift_ + nO3shift_*nV_ + nO4shift_*nV2shift_;
    this->outOfBound_ = (this->Hbar_dim + 1) * (this->Hbar_dim + 1);
    abIndices_.clear();
    abIndices_.resize(nV_, std::vector<size_t>(nV_, outOfBound_));
    ijIndices_.clear();
    ijIndices_.resize(nO_, std::vector<size_t>(nO_, outOfBound_));
    ijkIndices_.clear();
    ijkIndices_.resize(nO_, std::vector<std::vector<size_t>>(nO_, std::vector<size_t>(nO_, outOfBound_)));
    ijklIndices_.clear();
    ijklIndices_.resize(nO_, std::vector<std::vector<std::vector<size_t>>>(nO_, std::vector<std::vector<size_t>>(nO_, std::vector<size_t>(nO_, outOfBound_))));
    size_t idx = 0;
    for (size_t b = 0; b < nV_; b++) {
      for (size_t a = 0; a < b; a++) {
        abIndices_[a][b] = idx;
        abIndices_[b][a] = idx++;
      }
    }

    idx = 0;
    for (size_t j = 0; j < nO_; j++) {
      for (size_t i = 0; i < j; i++) {
        ijIndices_[i][j] = idx;
        ijIndices_[j][i] = idx++;
      }
    }

    idx = 0;
    for (size_t k = 0; k < nO_; k++) {
      for (size_t j = 0; j < k; j++) {
        for (size_t i = 0; i < j; i++) {
          ijkIndices_[i][j][k] = idx;
          ijkIndices_[j][i][k] = idx;
          ijkIndices_[i][k][j] = idx;
          ijkIndices_[j][k][i] = idx;
          ijkIndices_[k][i][j] = idx;
          ijkIndices_[k][j][i] = idx++;
        }
      }
    }

    idx = 0;
    for (size_t l = 0; l < nO_; l++) {
      for (size_t k = 0; k < l; k++) {
        for (size_t j = 0; j < k; j++) {
          for (size_t i = 0; i < j; i++) {
            ijklIndices_[i][j][k][l] = idx;
            ijklIndices_[j][i][k][l] = idx;
            ijklIndices_[i][k][j][l] = idx;
            ijklIndices_[j][k][i][l] = idx;
            ijklIndices_[k][i][j][l] = idx;
            ijklIndices_[k][j][i][l] = idx;
            ijklIndices_[i][j][l][k] = idx;
            ijklIndices_[j][i][l][k] = idx;
            ijklIndices_[i][k][l][j] = idx;
            ijklIndices_[j][k][l][i] = idx;
            ijklIndices_[k][i][l][j] = idx;
            ijklIndices_[k][j][l][i] = idx;
            ijklIndices_[i][l][j][k] = idx;
            ijklIndices_[j][l][i][k] = idx;
            ijklIndices_[i][l][k][j] = idx;
            ijklIndices_[j][l][k][i] = idx;
            ijklIndices_[k][l][i][j] = idx;
            ijklIndices_[k][l][j][i] = idx;
            ijklIndices_[l][i][j][k] = idx;
            ijklIndices_[l][j][i][k] = idx;
            ijklIndices_[l][i][k][j] = idx;
            ijklIndices_[l][j][k][i] = idx;
            ijklIndices_[l][k][i][j] = idx;
            ijklIndices_[l][k][j][i] = idx++;
          }
        }
      }
    }

    this->tensor_builder_.push_back(std::string(""));
    this->tensor_builder_.push_back(std::string({intermediates.oLabel,intermediates.oLabel}));
    this->tensor_builder_.push_back(std::string("TwoBody"));
    this->tensor_builder_.push_back(std::string({intermediates.vLabel}));
    this->tensor_builder_.push_back(std::string({intermediates.oLabel,intermediates.oLabel,intermediates.oLabel}));
    this->tensor_builder_.push_back(std::string("ThreeBody"));

    this->tensor_builder_.push_back(std::string({intermediates.vLabel,intermediates.vLabel}));
    this->tensor_builder_.push_back(std::string({intermediates.oLabel,intermediates.oLabel,intermediates.oLabel,intermediates.oLabel}));
    this->tensor_builder_.push_back(std::string("FourBody"));

  }

  template <typename MatsT>
  inline size_t EOMDIP_4h2pCCSDT<MatsT>::toCompoundS(size_t i, size_t j) const {
    if (i >= nO_ or j >= nO_)
      return outOfBound_;
    return ijIndices_[i][j];
  }

  template <typename MatsT>
  inline size_t EOMDIP_4h2pCCSDT<MatsT>::toCompoundD(size_t a, size_t i, size_t j, size_t k) const {
  size_t ijk = ijkIndices_[i][j][k];
  if (ijk == outOfBound_)
    return outOfBound_;
  return a + ijk * nV_;
  }

  template <typename MatsT>
  inline size_t EOMDIP_4h2pCCSDT<MatsT>::toCompoundT(size_t a, size_t b, size_t i, size_t j, size_t k, size_t l) const {
  size_t ijkl = ijklIndices_[i][j][k][l];
  size_t ab = abIndices_[a][b];
  if (ijkl == outOfBound_ or ab == outOfBound_)
      return outOfBound_;
    return ab + ijkl * nV2shift_;
  }

  template <typename MatsT>
  void EOMDIP_4h2pCCSDT<MatsT>::initializeEOMCC() {

    TAManager &TAmanager = TAManager::get();



    reused_.emplace(std::make_pair("1_vovo", TAmanager.malloc<MatsT>("vovo")));
    reused_.emplace(std::make_pair("2_oovv", TAmanager.malloc<MatsT>("oovv")));
    reused_.emplace(std::make_pair("3_ovov", TAmanager.malloc<MatsT>("ovov")));
    reused_.emplace(std::make_pair("4_ovov", TAmanager.malloc<MatsT>("ovov")));
    reused_.emplace(std::make_pair("5_oovo", TAmanager.malloc<MatsT>("oovo")));
    reused_.emplace(std::make_pair("6_oovo", TAmanager.malloc<MatsT>("oovo")));
    reused_.emplace(std::make_pair("7_oooo", TAmanager.malloc<MatsT>("oooo")));
    reused_.emplace(std::make_pair("8_voov", TAmanager.malloc<MatsT>("voov")));
    reused_.emplace(std::make_pair("9_vovv", TAmanager.malloc<MatsT>("vovv")));
    reused_.emplace(std::make_pair("10_vovv", TAmanager.malloc<MatsT>("vovv")));
    reused_.emplace(std::make_pair("11_vvoo", TAmanager.malloc<MatsT>("vvoo")));
    reused_.emplace(std::make_pair("12_vvoo", TAmanager.malloc<MatsT>("vvoo")));
    reused_.emplace(std::make_pair("13_vo", TAmanager.malloc<MatsT>("vo")));
    reused_.emplace(std::make_pair("14_voov", TAmanager.malloc<MatsT>("voov")));
    reused_.emplace(std::make_pair("15_vovo", TAmanager.malloc<MatsT>("vovo")));
    reused_.emplace(std::make_pair("16_vvvo", TAmanager.malloc<MatsT>("vvvo")));
    reused_.emplace(std::make_pair("17_vooo", TAmanager.malloc<MatsT>("vooo")));
    reused_.emplace(std::make_pair("18_voov", TAmanager.malloc<MatsT>("voov")));
    reused_.emplace(std::make_pair("19_vooo", TAmanager.malloc<MatsT>("vooo")));
    reused_.emplace(std::make_pair("20_voov", TAmanager.malloc<MatsT>("voov")));
    reused_.emplace(std::make_pair("21_vovo", TAmanager.malloc<MatsT>("vovo")));
    reused_.emplace(std::make_pair("22_vooo", TAmanager.malloc<MatsT>("vooo")));
    reused_.emplace(std::make_pair("23_vo", TAmanager.malloc<MatsT>("vo")));
    reused_.emplace(std::make_pair("24_ooov", TAmanager.malloc<MatsT>("ooov")));
    reused_.emplace(std::make_pair("25_ov", TAmanager.malloc<MatsT>("ov")));
    reused_.emplace(std::make_pair("26_ov", TAmanager.malloc<MatsT>("ov")));
    reused_.emplace(std::make_pair("27_ov", TAmanager.malloc<MatsT>("ov")));
    reused_.emplace(std::make_pair("28_vo", TAmanager.malloc<MatsT>("vo")));
    reused_.emplace(std::make_pair("29_vo", TAmanager.malloc<MatsT>("vo")));
    reused_.emplace(std::make_pair("30_oo", TAmanager.malloc<MatsT>("oo")));
    reused_.emplace(std::make_pair("31_vo", TAmanager.malloc<MatsT>("vo")));
    reused_.emplace(std::make_pair("32_oo", TAmanager.malloc<MatsT>("oo")));
    reused_.emplace(std::make_pair("33_ov", TAmanager.malloc<MatsT>("ov")));
    reused_.emplace(std::make_pair("34_vo", TAmanager.malloc<MatsT>("vo")));
    reused_.emplace(std::make_pair("35_vo", TAmanager.malloc<MatsT>("vo")));
    reused_.emplace(std::make_pair("36_oo", TAmanager.malloc<MatsT>("oo")));
    reused_.emplace(std::make_pair("37_vo", TAmanager.malloc<MatsT>("vo")));
    reused_.emplace(std::make_pair("38_vo", TAmanager.malloc<MatsT>("vo")));
    reused_.emplace(std::make_pair("39_vo", TAmanager.malloc<MatsT>("vo")));
    reused_.emplace(std::make_pair("40_vo", TAmanager.malloc<MatsT>("vo")));
    reused_.emplace(std::make_pair("41_vo", TAmanager.malloc<MatsT>("vo")));
    reused_.emplace(std::make_pair("42_oo", TAmanager.malloc<MatsT>("oo")));
    reused_.emplace(std::make_pair("43_oo", TAmanager.malloc<MatsT>("oo")));
    reused_.emplace(std::make_pair("44_oo", TAmanager.malloc<MatsT>("oo")));
    reused_.emplace(std::make_pair("45_oooo", TAmanager.malloc<MatsT>("oooo")));
    reused_.emplace(std::make_pair("46_vooo", TAmanager.malloc<MatsT>("vooo")));
    reused_.emplace(std::make_pair("47_vooo", TAmanager.malloc<MatsT>("vooo")));
    reused_.emplace(std::make_pair("48_ovoo", TAmanager.malloc<MatsT>("ovoo")));
    reused_.emplace(std::make_pair("49_ooov", TAmanager.malloc<MatsT>("ooov")));
    reused_.emplace(std::make_pair("50_vooo", TAmanager.malloc<MatsT>("vooo")));
    reused_.emplace(std::make_pair("51_vooo", TAmanager.malloc<MatsT>("vooo")));
    reused_.emplace(std::make_pair("52_vooo", TAmanager.malloc<MatsT>("vooo")));
    reused_.emplace(std::make_pair("53_oooo", TAmanager.malloc<MatsT>("oooo")));
    reused_.emplace(std::make_pair("54_vooo", TAmanager.malloc<MatsT>("vooo")));
    reused_.emplace(std::make_pair("55_oovo", TAmanager.malloc<MatsT>("oovo")));
    reused_.emplace(std::make_pair("56_vvoo", TAmanager.malloc<MatsT>("vvoo")));
    reused_.emplace(std::make_pair("57_voov", TAmanager.malloc<MatsT>("voov")));
    reused_.emplace(std::make_pair("58_vvoo", TAmanager.malloc<MatsT>("vvoo")));
    reused_.emplace(std::make_pair("59_vvoo", TAmanager.malloc<MatsT>("vvoo")));
    reused_.emplace(std::make_pair("60_vooo", TAmanager.malloc<MatsT>("vooo")));
    reused_.emplace(std::make_pair("61_voov", TAmanager.malloc<MatsT>("voov")));
    reused_.emplace(std::make_pair("62_vvoo", TAmanager.malloc<MatsT>("vvoo")));
    reused_.emplace(std::make_pair("63_vvoo", TAmanager.malloc<MatsT>("vvoo")));
    reused_.emplace(std::make_pair("64_vvoo", TAmanager.malloc<MatsT>("vvoo")));
    reused_.emplace(std::make_pair("65_vvoo", TAmanager.malloc<MatsT>("vvoo")));
    reused_.emplace(std::make_pair("66_vovv", TAmanager.malloc<MatsT>("vovv")));
    reused_.emplace(std::make_pair("67_oovv", TAmanager.malloc<MatsT>("oovv")));
    reused_.emplace(std::make_pair("68_oovv", TAmanager.malloc<MatsT>("oovv")));
    reused_.emplace(std::make_pair("69_oovv", TAmanager.malloc<MatsT>("oovv")));
    reused_.emplace(std::make_pair("70_ovvo", TAmanager.malloc<MatsT>("ovvo")));
    reused_.emplace(std::make_pair("71_vovo", TAmanager.malloc<MatsT>("vovo")));
    reused_.emplace(std::make_pair("72_ovvo", TAmanager.malloc<MatsT>("ovvo")));
    reused_.emplace(std::make_pair("73_ovvo", TAmanager.malloc<MatsT>("ovvo")));
    reused_.emplace(std::make_pair("74_ovvo", TAmanager.malloc<MatsT>("ovvo")));
    reused_.emplace(std::make_pair("75_ovvo", TAmanager.malloc<MatsT>("ovvo")));
    reused_.emplace(std::make_pair("76_vovo", TAmanager.malloc<MatsT>("vovo")));
    reused_.emplace(std::make_pair("77_ovvo", TAmanager.malloc<MatsT>("ovvo")));
    reused_.emplace(std::make_pair("78_vovo", TAmanager.malloc<MatsT>("vovo")));
    reused_.emplace(std::make_pair("79_vvoo", TAmanager.malloc<MatsT>("vvoo")));
    reused_.emplace(std::make_pair("80_vvoo", TAmanager.malloc<MatsT>("vvoo")));
    reused_.emplace(std::make_pair("81_vvoo", TAmanager.malloc<MatsT>("vvoo")));
    reused_.emplace(std::make_pair("82_voov", TAmanager.malloc<MatsT>("voov")));
    reused_.emplace(std::make_pair("83_vvoo", TAmanager.malloc<MatsT>("vvoo")));
    reused_.emplace(std::make_pair("84_vvoo", TAmanager.malloc<MatsT>("vvoo")));
    reused_.emplace(std::make_pair("85_voov", TAmanager.malloc<MatsT>("voov")));
    reused_.emplace(std::make_pair("86_vvoo", TAmanager.malloc<MatsT>("vvoo")));
    reused_.emplace(std::make_pair("87_vv", TAmanager.malloc<MatsT>("vv")));
    reused_.emplace(std::make_pair("88_vooo", TAmanager.malloc<MatsT>("vooo")));
    reused_.emplace(std::make_pair("89_vv", TAmanager.malloc<MatsT>("vv")));
    reused_.emplace(std::make_pair("90_vo", TAmanager.malloc<MatsT>("vo")));
    reused_.emplace(std::make_pair("91_vo", TAmanager.malloc<MatsT>("vo")));
    reused_.emplace(std::make_pair("92_vo", TAmanager.malloc<MatsT>("vo")));
    reused_.emplace(std::make_pair("93_vv", TAmanager.malloc<MatsT>("vv")));
    reused_.emplace(std::make_pair("94_vv", TAmanager.malloc<MatsT>("vv")));
    reused_.emplace(std::make_pair("95_oooo", TAmanager.malloc<MatsT>("oooo")));
    reused_.emplace(std::make_pair("96_vvov", TAmanager.malloc<MatsT>("vvov")));
    reused_.emplace(std::make_pair("97_vvov", TAmanager.malloc<MatsT>("vvov")));
    reused_.emplace(std::make_pair("98_oovv", TAmanager.malloc<MatsT>("oovv")));
    reused_.emplace(std::make_pair("99_ovvo", TAmanager.malloc<MatsT>("ovvo")));
    reused_.emplace(std::make_pair("100_vvoo", TAmanager.malloc<MatsT>("vvoo")));
    reused_.emplace(std::make_pair("101_oovvoo", TAmanager.malloc<MatsT>("oovvoo")));
    reused_.emplace(std::make_pair("102_voovoo", TAmanager.malloc<MatsT>("voovoo")));
    reused_.emplace(std::make_pair("103_oovvoo", TAmanager.malloc<MatsT>("oovvoo")));
    reused_.emplace(std::make_pair("104_ovvo", TAmanager.malloc<MatsT>("ovvo")));
    reused_.emplace(std::make_pair("105_vvoo", TAmanager.malloc<MatsT>("vvoo")));
    reused_.emplace(std::make_pair("106_ovov", TAmanager.malloc<MatsT>("ovov")));
    reused_.emplace(std::make_pair("107_ovov", TAmanager.malloc<MatsT>("ovov")));
    reused_.emplace(std::make_pair("108_vooooo", TAmanager.malloc<MatsT>("vooooo")));
    reused_.emplace(std::make_pair("109_oovvoo", TAmanager.malloc<MatsT>("oovvoo")));
    reused_.emplace(std::make_pair("110_ovvooo", TAmanager.malloc<MatsT>("ovvooo")));
    reused_.emplace(std::make_pair("111_vvvo", TAmanager.malloc<MatsT>("vvvo")));
    reused_.emplace(std::make_pair("112_oovoov", TAmanager.malloc<MatsT>("oovoov")));
    reused_.emplace(std::make_pair("113_oovoov", TAmanager.malloc<MatsT>("oovoov")));
    reused_.emplace(std::make_pair("114_vvoooo", TAmanager.malloc<MatsT>("vvoooo")));
    reused_.emplace(std::make_pair("115_vovooo", TAmanager.malloc<MatsT>("vovooo")));
    reused_.emplace(std::make_pair("116_voovoo", TAmanager.malloc<MatsT>("voovoo")));
    reused_.emplace(std::make_pair("117_voovoo", TAmanager.malloc<MatsT>("voovoo")));
    reused_.emplace(std::make_pair("118_voovoo", TAmanager.malloc<MatsT>("voovoo")));
    reused_.emplace(std::make_pair("119_vvoooo", TAmanager.malloc<MatsT>("vvoooo")));

  }


  template <typename MatsT>
  void EOMDIP_4h2pCCSDT<MatsT>::formEOMIntermediates() {

    TAManager &TAmanager = TAManager::get();

    reused_["1_vovo"]("a,l,b,k")  = conj(this->antiSymMoints["vvvo"]("c,d,a,m")) * this->T1_("c,l") * this->T2_("d,b,k,m");
    reused_["2_oovv"]("i,l,a,b")  = this->T2_("a,b,n,m") * this->antiSymMoints["oooo"]("m,n,i,l");
    reused_["3_ovov"]("j,b,l,a")  = this->T2_("c,b,l,n") * conj(this->antiSymMoints["vooo"]("c,j,m,n")) * this->T1_("a,m");
    reused_["4_ovov"]("l,b,i,a")  = conj(this->antiSymMoints["vvoo"]("c,d,m,n")) * this->T1_("c,l") * this->T2_("d,b,i,n") * this->T1_("a,m");
    reused_["5_oovo"]("m,j,a,k")  = this->T2_("b,a,k,l") * conj(this->antiSymMoints["vooo"]("b,j,l,m"));
    reused_["6_oovo"]("i,m,a,k")  = conj(this->antiSymMoints["vvoo"]("b,c,l,m")) * this->T1_("b,i") * this->T2_("c,a,k,l");
    reused_["7_oooo"]("i,j,l,m")  = conj(this->antiSymMoints["vvoo"]("b,c,l,m")) * this->T2_("b,c,i,j");
    reused_["8_voov"]("b,j,m,d")  = conj(this->antiSymMoints["vvoo"]("c,d,m,n")) * this->T2_("c,b,j,n");
    reused_["9_vovv"]("d,j,a,b")  = conj(this->antiSymMoints["vvoo"]("c,d,m,n")) * this->T1_("c,j") * this->T2_("a,b,n,m");
    reused_["10_vovv"]("c,l,a,b")  = this->T2_("a,b,n,m") * conj(this->antiSymMoints["vooo"]("c,l,m,n"));
    reused_["11_vvoo"]("a,b,i,l")  = this->fockMatrix_ta["ov"]("m,c") * this->T3_("c,a,b,i,l,m");
    reused_["12_vvoo"]("a,b,i,l")  = conj(this->antiSymMoints["vvoo"]("c,d,m,n")) * this->T1_("c,m") * this->T3_("d,a,b,i,l,n");
    reused_["13_vo"]("a,k")  = conj(this->antiSymMoints["vvoo"]("b,c,l,m")) * this->T3_("b,c,a,k,m,l");
    reused_["14_voov"]("a,k,m,c")  = conj(this->antiSymMoints["vvoo"]("b,c,l,m")) * this->T2_("b,a,k,l");
    reused_["15_vovo"]("a,i,b,l")  = this->antiSymMoints["vovo"]("a,m,c,i") * this->T2_("c,b,l,m");
    reused_["16_vvvo"]("a,b,d,j")  = this->antiSymMoints["vvvv"]("a,b,c,d") * this->T1_("c,j");
    reused_["17_vooo"]("a,l,j,k")  = conj(this->antiSymMoints["vvvo"]("b,c,a,l")) * this->T2_("b,c,j,k");
    reused_["18_voov"]("a,k,l,b")  = reused_["17_vooo"]("a,m,k,l") * this->T1_("b,m");
    reused_["19_vooo"]("a,l,j,k")  = this->antiSymMoints["vovo"]("a,l,b,j") * this->T1_("b,k");
    reused_["20_voov"]("a,l,k,b")  = reused_["19_vooo"]("a,m,l,k") * this->T1_("b,m");
    reused_["21_vovo"]("a,l,c,k")  = conj(this->antiSymMoints["vvvo"]("b,c,a,l")) * this->T1_("b,k");
    reused_["22_vooo"]("a,l,k,j")  = reused_["21_vovo"]("a,l,c,k") * this->T1_("c,j");
    reused_["23_vo"]("a,k")  = conj(this->antiSymMoints["vvvo"]("b,c,a,l")) * this->T2_("b,c,k,l");
    reused_["24_ooov"]("j,k,l,b")  = conj(this->antiSymMoints["vvoo"]("a,b,k,l")) * this->T1_("a,j");
    reused_["25_ov"]("i,a")  = conj(this->antiSymMoints["vvoo"]("b,c,l,m")) * this->T2_("b,c,i,m") * this->T1_("a,l");
    reused_["26_ov"]("i,a")  = this->T2_("b,a,m,l") * conj(this->antiSymMoints["vooo"]("b,i,l,m"));
    reused_["27_ov"]("l,b")  = conj(this->antiSymMoints["vvoo"]("a,b,k,l")) * this->T1_("a,k");
    reused_["28_vo"]("a,k")  = this->fockMatrix_ta["ov"]("l,b") * this->T2_("b,a,k,l");
    reused_["29_vo"]("a,i")  = this->antiSymMoints["vovo"]("a,l,b,i") * this->T1_("b,l");
    reused_["30_oo"]("l,j")  = conj(this->antiSymMoints["vooo"]("a,j,k,l")) * this->T1_("a,k");
    reused_["31_vo"]("a,k")  = this->fockMatrix_ta["vv"]("a,b") * this->T1_("b,k");
    reused_["32_oo"]("j,k")  = this->fockMatrix_ta["ov"]("k,a") * this->T1_("a,j");
    reused_["33_ov"]("k,a")  = this->T1_("a,l") * this->fockMatrix_ta["oo"]("l,k");
    reused_["34_vo"]("a,i")  = reused_["29_vo"]("a,i");
    reused_["34_vo"]("a,i") -= reused_["31_vo"]("a,i");
    reused_["34_vo"]("a,i") += reused_["28_vo"]("a,i");
    reused_["34_vo"]("a,i") += reused_["33_ov"]("i,a");
    reused_["34_vo"]("a,i") += 0.25 * reused_["13_vo"]("a,i");
    reused_["34_vo"]("a,i") -= 0.50 * reused_["23_vo"]("a,i");
    reused_["34_vo"]("a,i") -= 0.50 * reused_["26_ov"]("i,a");
    reused_["34_vo"]("a,i") += 0.50 * reused_["25_ov"]("i,a");
    reused_["35_vo"]("a,i")  = this->T2_("c,a,m,l") * reused_["24_ooov"]("i,l,m,c");
    reused_["36_oo"]("j,l")  = this->T1_("b,j") * reused_["27_ov"]("l,b");
    reused_["37_vo"]("a,i")  = this->T1_("a,l") * reused_["32_oo"]("i,l");
    reused_["38_vo"]("a,i")  = this->T1_("a,m") * reused_["30_oo"]("m,i");
    reused_["39_vo"]("a,i")  = reused_["35_vo"]("a,i");
    reused_["39_vo"]("a,i") += 2.00 * reused_["34_vo"]("a,i");
    reused_["39_vo"]("a,i") += 2.00 * reused_["38_vo"]("a,i");
    reused_["39_vo"]("a,i") += 2.00 * reused_["37_vo"]("a,i");
    reused_["40_vo"]("a,i")  = reused_["36_oo"]("i,m") * this->T1_("a,m");
    reused_["41_vo"]("a,i")  = reused_["40_vo"]("a,i");
    reused_["41_vo"]("a,i") += 0.50 * reused_["39_vo"]("a,i");
    reused_["42_oo"]("j,l")  = conj(this->antiSymMoints["vvoo"]("a,b,k,l")) * this->T2_("a,b,j,k");
    reused_["43_oo"]("m,i")  = reused_["30_oo"]("m,i");
    reused_["43_oo"]("m,i") -= 0.50 * reused_["42_oo"]("i,m");
    reused_["44_oo"]("m,i")  = reused_["43_oo"]("m,i");
    reused_["44_oo"]("m,i") += reused_["36_oo"]("i,m");
    reused_["45_oooo"]("i,j,k,l")  = this->T1_("b,i") * reused_["24_ooov"]("j,k,l,b");
    reused_["46_vooo"]("a,j,k,m")  = reused_["45_oooo"]("j,k,l,m") * this->T1_("a,l");
    reused_["47_vooo"]("a,j,k,m")  = conj(this->antiSymMoints["vvoo"]("b,c,l,m")) * this->T3_("b,c,a,j,k,l");
    reused_["48_ovoo"]("m,a,i,j")  = conj(this->antiSymMoints["vvoo"]("b,c,l,m")) * this->T1_("b,l") * this->T2_("c,a,i,j");
    reused_["49_ooov"]("m,i,j,a")  = this->T1_("a,l") * this->antiSymMoints["oooo"]("l,m,i,j");
    reused_["50_vooo"]("a,i,j,m")  = reused_["47_vooo"]("a,i,j,m");
    reused_["50_vooo"]("a,i,j,m") += 2.00 * reused_["49_ooov"]("m,i,j,a");
    reused_["51_vooo"]("a,i,j,m")  = this->T1_("a,l") * reused_["7_oooo"]("i,j,l,m");
    reused_["52_vooo"]("a,i,j,m")  = reused_["51_vooo"]("a,i,j,m");
    reused_["52_vooo"]("a,i,j,m") += 2.00 * reused_["48_ovoo"]("m,a,i,j");
    reused_["52_vooo"]("a,i,j,m") += reused_["50_vooo"]("a,i,j,m");
    reused_["53_oooo"]("k,l,m,j")  = conj(this->antiSymMoints["vooo"]("b,j,l,m")) * this->T1_("b,k");
    reused_["54_vooo"]("a,j,m,k")  = this->T1_("a,l") * reused_["53_oooo"]("j,l,m,k");
    reused_["55_oovo"]("m,j,a,k")  = reused_["5_oovo"]("m,j,a,k");
    reused_["55_oovo"]("m,j,a,k") += reused_["54_vooo"]("a,k,m,j");
    reused_["56_vvoo"]("a,b,i,l")  = conj(this->antiSymMoints["vvvo"]("c,d,a,m")) * this->T3_("c,d,b,i,l,m");
    reused_["57_voov"]("b,i,l,a")  = conj(this->antiSymMoints["vvoo"]("c,d,m,n")) * this->T3_("c,d,b,i,l,n") * this->T1_("a,m");
    reused_["58_vvoo"]("a,b,k,l")  = conj(this->antiSymMoints["vvvo"]("c,d,a,m")) * this->T1_("c,m") * this->T2_("d,b,k,l");
    reused_["59_vvoo"]("a,b,k,l")  = this->fockMatrix_ta["vv"]("a,c") * this->T2_("c,b,k,l");
    reused_["60_vooo"]("a,i,j,l")  = this->fockMatrix_ta["ov"]("l,b") * this->T2_("b,a,i,j");
    reused_["61_voov"]("a,i,l,b")  = this->antiSymMoints["vooo"]("a,m,i,l") * this->T1_("b,m");
    reused_["62_vvoo"]("a,b,k,l")  = reused_["56_vvoo"]("a,b,k,l");
    reused_["62_vvoo"]("a,b,k,l") += 2.00 * reused_["58_vvoo"]("a,b,k,l");
    reused_["62_vvoo"]("a,b,k,l") -= reused_["57_voov"]("b,k,l,a");
    reused_["62_vvoo"]("a,b,k,l") += 2.00 * reused_["61_voov"]("a,k,l,b");
    reused_["62_vvoo"]("a,b,k,l") -= 2.00 * reused_["59_vvoo"]("a,b,k,l");
    reused_["63_vvoo"]("a,b,k,l")  = this->T1_("a,n") * reused_["48_ovoo"]("n,b,k,l");
    reused_["64_vvoo"]("a,b,j,k")  = this->T1_("a,m") * reused_["60_vooo"]("b,j,k,m");
    reused_["65_vvoo"]("a,b,k,l")  = reused_["63_vvoo"]("a,b,k,l");
    reused_["65_vvoo"]("a,b,k,l") += 0.50 * reused_["62_vvoo"]("a,b,k,l");
    reused_["65_vvoo"]("a,b,k,l") += reused_["64_vvoo"]("a,b,k,l");
    reused_["66_vovv"]("d,k,a,b")  = reused_["9_vovv"]("d,k,a,b");
    reused_["66_vovv"]("d,k,a,b") -= 2.00 * reused_["16_vvvo"]("a,b,d,k");
    reused_["67_oovv"]("j,k,a,b")  = this->T1_("d,j") * reused_["66_vovv"]("d,k,a,b");
    reused_["68_oovv"]("k,l,a,b")  = this->T1_("d,k") * reused_["24_ooov"]("l,m,n,d") * this->T1_("a,m") * this->T1_("b,n");
    reused_["69_oovv"]("j,k,a,b")  = reused_["67_oovv"]("j,k,a,b");
    reused_["69_oovv"]("j,k,a,b") -= 2.00 * reused_["68_oovv"]("j,k,a,b");
    reused_["70_ovvo"]("j,a,b,l")  = this->T3_("c,a,b,l,n,m") * conj(this->antiSymMoints["vooo"]("c,j,m,n"));
    reused_["71_vovo"]("a,j,b,l")  = conj(this->antiSymMoints["vvoo"]("c,d,m,n")) * this->T2_("c,a,j,m") * this->T2_("d,b,l,n");
    reused_["72_ovvo"]("l,a,b,i")  = conj(this->antiSymMoints["vvoo"]("c,d,m,n")) * this->T1_("c,m") * this->T1_("d,l") * this->T2_("a,b,i,n");
    reused_["73_ovvo"]("l,a,b,k")  = conj(this->antiSymMoints["vooo"]("c,l,m,n")) * this->T1_("c,m") * this->T2_("a,b,k,n");
    reused_["74_ovvo"]("l,a,b,i")  = this->fockMatrix_ta["ov"]("m,c") * this->T1_("c,l") * this->T2_("a,b,i,m");
    reused_["75_ovvo"]("k,a,b,l")  = this->T2_("a,b,l,m") * this->fockMatrix_ta["oo"]("m,k");
    reused_["76_vovo"]("a,k,b,l")  = reused_["71_vovo"]("a,k,b,l");
    reused_["76_vovo"]("a,k,b,l") += reused_["75_ovvo"]("k,a,b,l");
    reused_["76_vovo"]("a,k,b,l") += reused_["74_ovvo"]("k,a,b,l");
    reused_["76_vovo"]("a,k,b,l") += 0.50 * reused_["70_ovvo"]("k,a,b,l");
    reused_["76_vovo"]("a,k,b,l") += reused_["72_ovvo"]("k,a,b,l");
    reused_["76_vovo"]("a,k,b,l") += reused_["73_ovvo"]("k,a,b,l");
    reused_["77_ovvo"]("i,a,b,l")  = reused_["10_vovv"]("c,i,a,b") * this->T1_("c,l");
    reused_["78_vovo"]("a,k,b,l")  = reused_["76_vovo"]("a,k,b,l");
    reused_["78_vovo"]("a,k,b,l") += 0.50 * reused_["77_ovvo"]("k,a,b,l");
    reused_["79_vvoo"]("a,b,j,k")  = this->antiSymMoints["vvvv"]("a,b,c,d") * this->T2_("c,d,j,k");
    reused_["80_vvoo"]("b,a,k,l")  = conj(this->antiSymMoints["vvoo"]("c,d,m,n")) * this->T2_("d,b,n,m") * this->T2_("c,a,k,l");
    reused_["81_vvoo"]("a,b,i,l")  = conj(this->antiSymMoints["vvoo"]("c,d,m,n")) * this->T2_("c,a,n,m") * this->T2_("d,b,i,l");
    reused_["82_voov"]("b,j,k,a")  = this->T1_("a,m") * this->antiSymMoints["oooo"]("m,n,j,k") * this->T1_("b,n");
    reused_["83_vvoo"]("a,b,k,l")  = reused_["12_vvoo"]("a,b,k,l");
    reused_["83_vvoo"]("a,b,k,l") += 0.50 * reused_["81_vvoo"]("a,b,k,l");
    reused_["83_vvoo"]("a,b,k,l") -= 0.50 * reused_["2_oovv"]("k,l,a,b");
    reused_["83_vvoo"]("a,b,k,l") += 0.50 * reused_["80_vvoo"]("b,a,k,l");
    reused_["83_vvoo"]("a,b,k,l") += 0.50 * reused_["79_vvoo"]("a,b,k,l");
    reused_["83_vvoo"]("a,b,k,l") += reused_["82_voov"]("b,k,l,a");
    reused_["83_vvoo"]("a,b,k,l") += reused_["11_vvoo"]("a,b,k,l");
    reused_["84_vvoo"]("a,b,k,l")  = this->T2_("a,b,n,m") * reused_["7_oooo"]("k,l,m,n");
    reused_["85_voov"]("a,k,l,b")  = this->T1_("a,m") * reused_["7_oooo"]("k,l,m,n") * this->T1_("b,n");
    reused_["86_vvoo"]("a,b,k,l")  = reused_["83_vvoo"]("a,b,k,l");
    reused_["86_vvoo"]("a,b,k,l") -= 0.25 * reused_["84_vvoo"]("a,b,k,l");
    reused_["86_vvoo"]("a,b,k,l") += 0.50 * reused_["85_voov"]("a,k,l,b");
    reused_["87_vv"]("c,a")  = conj(this->antiSymMoints["vvoo"]("b,c,l,m")) * this->T2_("b,a,m,l");
    reused_["88_vooo"]("a,i,j,m")  = -1.00 * reused_["52_vooo"]("a,i,j,m");
    reused_["89_vv"]("a,c")  = conj(this->antiSymMoints["vvvo"]("b,c,a,l")) * this->T1_("b,l");
    reused_["90_vo"]("a,k")  = reused_["27_ov"]("m,c") * this->T2_("c,a,k,m");
    reused_["91_vo"]("a,k")  = reused_["89_vv"]("a,c") * this->T1_("c,k");
    reused_["92_vo"]("a,i")  = reused_["90_vo"]("a,i");
    reused_["92_vo"]("a,i") += reused_["91_vo"]("a,i");
    reused_["93_vv"]("a,c")  = conj(this->antiSymMoints["vvoo"]("b,c,l,m")) * this->T2_("b,a,m,l");
    reused_["94_vv"]("a,c")  = reused_["93_vv"]("a,c");
    reused_["94_vv"]("a,c") -= 2.00 * reused_["89_vv"]("a,c");
    reused_["95_oooo"]("i,j,k,l")  = reused_["7_oooo"]("i,j,k,l");
    reused_["95_oooo"]("i,j,k,l") -= 2.00 * reused_["45_oooo"]("i,j,k,l");
    reused_["96_vvov"]("a,b,l,d")  = conj(this->antiSymMoints["vvoo"]("c,d,m,n")) * this->T3_("c,a,b,l,n,m");
    reused_["97_vvov"]("a,b,l,d")  = reused_["96_vvov"]("a,b,l,d");
    reused_["97_vvov"]("a,b,l,d") += reused_["66_vovv"]("d,l,a,b");
    reused_["98_oovv"]("l,k,a,b")  = this->T1_("c,k") * conj(this->antiSymMoints["vooo"]("c,l,m,n")) * this->T1_("a,m") * this->T1_("b,n");
    reused_["99_ovvo"]("j,a,b,i")  = conj(this->antiSymMoints["vvoo"]("c,d,m,n")) * this->T2_("c,d,j,m") * this->T2_("a,b,i,n");
    reused_["100_vvoo"]("a,b,i,l")  = this->antiSymMoints["vvvo"]("a,b,c,i") * this->T1_("c,l");
    reused_["101_oovvoo"]("n,k,a,b,j,l")  = this->T3_("c,a,b,j,l,m") * conj(this->antiSymMoints["vooo"]("c,k,m,n"));
    reused_["102_voovoo"]("a,j,n,b,i,k")  = conj(this->antiSymMoints["vvoo"]("c,d,m,n")) * this->T2_("c,a,j,m") * this->T2_("d,b,i,k");
    reused_["103_oovvoo"]("n,k,a,b,j,l")  = reused_["101_oovvoo"]("n,k,a,b,j,l");
    reused_["103_oovvoo"]("n,k,a,b,j,l") += reused_["102_voovoo"]("a,k,n,b,j,l");
    reused_["104_ovvo"]("k,a,b,l")  = conj(this->antiSymMoints["vvoo"]("c,d,m,n")) * this->T1_("c,k") * this->T3_("d,a,b,l,n,m");
    reused_["105_vvoo"]("a,b,l,k")  = reused_["100_vvoo"]("a,b,l,k");
    reused_["105_vvoo"]("a,b,l,k") += 0.50 * reused_["104_ovvo"]("l,a,b,k");
    reused_["105_vvoo"]("a,b,l,k") += reused_["98_oovv"]("l,k,a,b");
    reused_["105_vvoo"]("a,b,l,k") += 0.50 * reused_["99_ovvo"]("l,a,b,k");
    reused_["106_ovov"]("l,b,k,a")  = reused_["3_ovov"]("l,b,k,a");
    reused_["106_ovov"]("l,b,k,a") += reused_["1_vovo"]("a,l,b,k");
    reused_["107_ovov"]("k,b,l,a")  = reused_["4_ovov"]("k,b,l,a");
    reused_["107_ovov"]("k,b,l,a") += reused_["15_vovo"]("a,k,b,l");
    reused_["108_vooooo"]("a,i,j,l,m,n")  = conj(this->antiSymMoints["vvoo"]("c,d,m,n")) * this->T3_("c,d,a,i,j,l");
    reused_["109_oovvoo"]("l,n,a,b,j,k")  = conj(this->antiSymMoints["vvoo"]("c,d,m,n")) * this->T1_("c,l") * this->T3_("d,a,b,j,k,m");
    reused_["110_ovvooo"]("n,a,b,i,j,l")  = conj(this->antiSymMoints["vvoo"]("c,d,m,n")) * this->T1_("c,m") * this->T3_("d,a,b,i,j,l");
    reused_["111_vvvo"]("a,d,b,j")  = conj(this->antiSymMoints["vvvo"]("c,d,a,m")) * this->T2_("c,b,j,m");
    reused_["112_oovoov"]("n,l,b,i,j,a")  = this->T2_("c,b,i,j") * conj(this->antiSymMoints["vooo"]("c,l,m,n")) * this->T1_("a,m");
    reused_["113_oovoov"]("i,n,b,k,l,a")  = conj(this->antiSymMoints["vvoo"]("c,d,m,n")) * this->T1_("c,i") * this->T2_("d,b,k,l") * this->T1_("a,m");
    reused_["114_vvoooo"]("a,b,j,k,l,m")  = this->fockMatrix_ta["ov"]("m,c") * this->T3_("c,a,b,j,k,l");
    reused_["115_vovooo"]("a,m,b,i,j,l")  = conj(this->antiSymMoints["vvvo"]("c,d,a,m")) * this->T3_("c,d,b,i,j,l");
    reused_["116_voovoo"]("a,m,j,b,k,l")  = this->antiSymMoints["vovo"]("a,m,c,j") * this->T2_("c,b,k,l");
    reused_["117_voovoo"]("a,m,k,b,j,l")  = conj(this->antiSymMoints["vvvo"]("c,d,a,m")) * this->T1_("c,k") * this->T2_("d,b,j,l");
    reused_["118_voovoo"]("b,l,n,a,i,j")  = conj(this->antiSymMoints["vvoo"]("c,d,m,n")) * this->T2_("d,b,l,m") * this->T2_("c,a,i,j");
    reused_["119_vvoooo"]("a,b,j,k,l,n")  = this->T1_("a,m") * reused_["108_vooooo"]("b,j,k,l,m,n");

}

  template <typename MatsT>
  void EOMDIP_4h2pCCSDT<MatsT>::buildSigma(const MBExpansion<MatsT> &V, MBExpansion<MatsT> &HV, EOMCCEigenVecType vecType) const {
    const TArray &V4 = V.get_tensor("FourBody");
    TArray &HV4 = HV.get_tensor("FourBody");
    const TArray &V3 = V.get_tensor("ThreeBody");
    TArray &HV3 = HV.get_tensor("ThreeBody");
    const TArray &V2 = V.get_tensor("TwoBody");
    TArray &HV2 = HV.get_tensor("TwoBody");

    switch (vecType) {
      case EOMCCEigenVecType::RIGHT:
        formR_tilde(V2, V3, V4, HV2, HV3, HV4);
        break;
      case EOMCCEigenVecType::LEFT:
        formL_tilde(V2, V3, V4, HV2, HV3, HV4);
        break;
    }
  }

  template <typename MatsT>
  void EOMDIP_4h2pCCSDT<MatsT>::formR_tilde(const TArray &R2, const TArray &R3, const TArray &R4, TArray &sigmaR2, TArray &sigmaR3, TArray &sigmaR4) const {
    TAManager &TAmanager = TAManager::get();    // sigmaR2  = -1.00 P(i,j) f(k,j) R2(i,k) 
    // flops: o2v0L1  = o3v0L1
    //  mems: o2v0L1  = o2v0L1
    tmps_["perm_Loo"] = TAmanager.malloc<MatsT>("oo");
    tmps_["perm_Loo"]("i,j")  = R2("i,k") * this->fockMatrix_ta["oo"]("k,j");
    sigmaR2("i,j")  = -1.00 * tmps_["perm_Loo"]("i,j");
    sigmaR2("i,j") += tmps_["perm_Loo"]("j,i");
    TAmanager.free("oo", std::move(tmps_["perm_Loo"]));
    
    // sigmaR3  = +1.00 f(a,b) R3(b,i,j,k) 
    // flops: o3v1L1  = o3v2L1
    //  mems: o3v1L1  = o3v1L1
    sigmaR3("a,i,j,k")  = this->fockMatrix_ta["vv"]("a,b") * R3("b,i,j,k");
    
    // sigmaR4  = +1.00 P(a,b) f(a,c) R4(c,b,i,j,k,l) 
    // flops: o4v2L1  = o4v3L1
    //  mems: o4v2L1  = o4v2L1
    tmps_["perm_Lvvoooo"] = TAmanager.malloc<MatsT>("vvoooo");
    tmps_["perm_Lvvoooo"]("a,b,i,j,k,l")  = this->fockMatrix_ta["vv"]("a,c") * R4("c,b,i,j,k,l");
    sigmaR4("a,b,i,j,k,l")  = tmps_["perm_Lvvoooo"]("a,b,i,j,k,l");
    sigmaR4("a,b,i,j,k,l") -= tmps_["perm_Lvvoooo"]("b,a,i,j,k,l");
    TAmanager.free("vvoooo", std::move(tmps_["perm_Lvvoooo"]));
    
    // sigmaR2 += -1.00 P(i,j) f(k,a) R2(i,k) this->T1_(a,j) 
    // flops: o2v0L1 += o3v0L1
    //  mems: o2v0L1 += o2v0L1
    tmps_["perm_Loo"] = TAmanager.malloc<MatsT>("oo");
    tmps_["perm_Loo"]("i,j")  = R2("i,k") * reused_["32_oo"]("j,k");
    sigmaR2("i,j") -= tmps_["perm_Loo"]("i,j");
    sigmaR2("i,j") += tmps_["perm_Loo"]("j,i");
    TAmanager.free("oo", std::move(tmps_["perm_Loo"]));
    
    // sigmaR2 += +1.00 P(i,j) <l,k||a,j> R2(i,l) this->T1_(a,k) 
    //            += -0.50 P(i,j) <l,k||a,b> R2(i,l) this->T2_(a,b,j,k) 
    //            += +1.00 P(i,j) <l,k||a,b> R2(i,l) this->T1_(a,k) this->T1_(b,j) 
    // flops: o2v0L1 += o3v0L1
    //  mems: o2v0L1 += o2v0L1
    tmps_["perm_Loo"] = TAmanager.malloc<MatsT>("oo");
    tmps_["perm_Loo"]("i,j")  = R2("i,l") * reused_["44_oo"]("l,j");
    sigmaR2("i,j") -= tmps_["perm_Loo"]("i,j");
    sigmaR2("i,j") += tmps_["perm_Loo"]("j,i");
    TAmanager.free("oo", std::move(tmps_["perm_Loo"]));
    
    // sigmaR2 += +0.50 <l,k||i,j> R2(l,k) 
    // flops: o2v0L1 += o4v0L1
    //  mems: o2v0L1 += o2v0L1
    sigmaR2("i,j") -= 0.50 * this->antiSymMoints["oooo"]("k,l,i,j") * R2("l,k");
    
    // sigmaR2 += +0.50 P(i,j) <l,k||a,j> R2(l,k) this->T1_(a,i) 
    // flops: o2v0L1 += o4v0L1
    //  mems: o2v0L1 += o2v0L1
    tmps_["perm_Loo"] = TAmanager.malloc<MatsT>("oo");
    tmps_["perm_Loo"]("i,j")  = 0.50 * R2("l,k") * reused_["53_oooo"]("i,k,l,j");
    sigmaR2("i,j") -= tmps_["perm_Loo"]("i,j");
    sigmaR2("i,j") += tmps_["perm_Loo"]("j,i");
    TAmanager.free("oo", std::move(tmps_["perm_Loo"]));
    
    // sigmaR2 += +0.25 <l,k||a,b> R2(l,k) this->T2_(a,b,i,j) 
    //            += -0.50 <l,k||a,b> R2(l,k) this->T1_(a,j) this->T1_(b,i) 
    // flops: o2v0L1 += o4v0L1
    //  mems: o2v0L1 += o2v0L1
    sigmaR2("i,j") -= 0.25 * R2("l,k") * reused_["95_oooo"]("i,j,k,l");
    
    // sigmaR2 += +1.00 f(k,a) R3(a,i,j,k) 
    // flops: o2v0L1 += o3v1L1
    //  mems: o2v0L1 += o2v0L1
    sigmaR2("i,j") += this->fockMatrix_ta["ov"]("k,a") * R3("a,i,j,k");
    
    // sigmaR2 += -1.00 <l,k||a,b> R3(b,i,j,l) this->T1_(a,k) 
    // flops: o2v0L1 += o3v1L1
    //  mems: o2v0L1 += o2v0L1
    sigmaR2("i,j") += R3("b,i,j,l") * reused_["27_ov"]("l,b");
    
    // sigmaR2 += +0.50 P(i,j) <l,k||a,j> R3(a,i,l,k) 
    // flops: o2v0L1 += o4v1L1
    //  mems: o2v0L1 += o2v0L1
    tmps_["perm_Loo"] = TAmanager.malloc<MatsT>("oo");
    tmps_["perm_Loo"]("i,j")  = 0.50 * R3("a,i,l,k") * conj(this->antiSymMoints["vooo"]("a,j,k,l"));
    sigmaR2("i,j") -= tmps_["perm_Loo"]("i,j");
    sigmaR2("i,j") += tmps_["perm_Loo"]("j,i");
    TAmanager.free("oo", std::move(tmps_["perm_Loo"]));
    
    // sigmaR2 += -0.50 P(i,j) <l,k||a,b> R3(b,i,l,k) this->T1_(a,j) 
    // flops: o2v0L1 += o4v1L1
    //  mems: o2v0L1 += o2v0L1
    tmps_["perm_Loo"] = TAmanager.malloc<MatsT>("oo");
    tmps_["perm_Loo"]("i,j")  = 0.50 * R3("b,i,l,k") * reused_["24_ooov"]("j,k,l,b");
    sigmaR2("i,j") += tmps_["perm_Loo"]("i,j");
    sigmaR2("i,j") -= tmps_["perm_Loo"]("j,i");
    TAmanager.free("oo", std::move(tmps_["perm_Loo"]));
    
    // sigmaR2 += +0.25 <l,k||a,b> R4(a,b,i,j,l,k) 
    // flops: o2v0L1 += o4v2L1
    //  mems: o2v0L1 += o2v0L1
    sigmaR2("i,j") -= 0.25 * conj(this->antiSymMoints["vvoo"]("a,b,k,l")) * R4("a,b,i,j,l,k");
    
    // sigmaR3 += -0.50 <m,l||b,c> R3(c,i,j,k) this->T2_(b,a,m,l) 
    //              += +1.00 <l,a||b,c> R3(c,i,j,k) this->T1_(b,l) 
    // flops: o3v1L1 += o3v2L1
    //  mems: o3v1L1 += o3v1L1
    sigmaR3("a,i,j,k") += 0.50 * R3("c,i,j,k") * reused_["94_vv"]("a,c");
    
    // sigmaR3 += +0.25 <m,l||b,c> R2(m,l) this->T3_(b,c,a,i,j,k) 
    // flops: o3v1L1 += o5v1L1
    //  mems: o3v1L1 += o3v1L1
    sigmaR3("a,i,j,k") -= 0.25 * R2("m,l") * reused_["108_vooooo"]("a,i,j,k,l,m");
    
    // sigmaR3 += -1.00 f(l,b) R4(b,a,i,j,k,l) 
    // flops: o3v1L1 += o4v2L1
    //  mems: o3v1L1 += o3v1L1
    sigmaR3("a,i,j,k") -= this->fockMatrix_ta["ov"]("l,b") * R4("b,a,i,j,k,l");
    
    // sigmaR3 += +1.00 <m,l||b,c> R4(c,a,i,j,k,m) this->T1_(b,l) 
    // flops: o3v1L1 += o4v2L1
    //  mems: o3v1L1 += o3v1L1
    sigmaR3("a,i,j,k") -= R4("c,a,i,j,k,m") * reused_["27_ov"]("m,c");
    
    // sigmaR3 += -0.50 <l,a||b,c> R4(b,c,i,j,k,l) 
    // flops: o3v1L1 += o4v3L1
    //  mems: o3v1L1 += o3v1L1
    sigmaR3("a,i,j,k") += 0.50 * conj(this->antiSymMoints["vvvo"]("b,c,a,l")) * R4("b,c,i,j,k,l");
    
    // sigmaR3 += -1.00 f(l,b) R3(b,i,j,k) this->T1_(a,l) 
    // flops: o3v1L1 += o4v1L1 o4v1L1
    //  mems: o3v1L1 += o4v0L1 o3v1L1
    sigmaR3("a,i,j,k") -= this->fockMatrix_ta["ov"]("l,b") * R3("b,i,j,k") * this->T1_("a,l");
    
    // sigmaR3 += +1.00 <m,l||b,c> R3(c,i,j,k) this->T1_(a,m) this->T1_(b,l) 
    // flops: o3v1L1 += o4v1L1 o4v1L1
    //  mems: o3v1L1 += o4v0L1 o3v1L1
    sigmaR3("a,i,j,k") -= R3("c,i,j,k") * reused_["27_ov"]("m,c") * this->T1_("a,m");
    
    // sigmaR4 += -0.50 P(a,b) <n,m||c,d> R4(d,b,i,j,k,l) this->T2_(c,a,n,m) 
    //                += +1.00 P(a,b) <m,a||c,d> R4(d,b,i,j,k,l) this->T1_(c,m) 
    // flops: o4v2L1 += o4v3L1
    //  mems: o4v2L1 += o4v2L1
    tmps_["perm_Lvvoooo"] = TAmanager.malloc<MatsT>("vvoooo");
    tmps_["perm_Lvvoooo"]("a,b,i,j,k,l")  = 0.50 * R4("d,b,i,j,k,l") * reused_["94_vv"]("a,d");
    sigmaR4("a,b,i,j,k,l") += tmps_["perm_Lvvoooo"]("a,b,i,j,k,l");
    sigmaR4("a,b,i,j,k,l") -= tmps_["perm_Lvvoooo"]("b,a,i,j,k,l");
    TAmanager.free("vvoooo", std::move(tmps_["perm_Lvvoooo"]));
    
    // sigmaR4 += +0.50 <a,b||c,d> R4(c,d,i,j,k,l) 
    // flops: o4v2L1 += o4v4L1
    //  mems: o4v2L1 += o4v2L1
    sigmaR4("a,b,i,j,k,l") += 0.50 * this->antiSymMoints["vvvv"]("a,b,c,d") * R4("c,d,i,j,k,l");
    
    // sigmaR4 += -1.00 P(a,b) f(m,c) R4(c,b,i,j,k,l) this->T1_(a,m) 
    // flops: o4v2L1 += o5v2L1 o5v2L1
    //  mems: o4v2L1 += o5v1L1 o4v2L1
    tmps_["perm_Lvvoooo"] = TAmanager.malloc<MatsT>("vvoooo");
    tmps_["perm_Lvvoooo"]("a,b,i,j,k,l")  = this->fockMatrix_ta["ov"]("m,c") * R4("c,b,i,j,k,l") * this->T1_("a,m");
    sigmaR4("a,b,i,j,k,l") -= tmps_["perm_Lvvoooo"]("a,b,i,j,k,l");
    sigmaR4("a,b,i,j,k,l") += tmps_["perm_Lvvoooo"]("b,a,i,j,k,l");
    TAmanager.free("vvoooo", std::move(tmps_["perm_Lvvoooo"]));
    
    // sigmaR4 += +1.00 P(a,b) <n,m||c,d> R4(d,b,i,j,k,l) this->T1_(a,n) this->T1_(c,m) 
    // flops: o4v2L1 += o5v2L1 o5v2L1
    //  mems: o4v2L1 += o5v1L1 o4v2L1
    tmps_["perm_Lvvoooo"] = TAmanager.malloc<MatsT>("vvoooo");
    tmps_["perm_Lvvoooo"]("a,b,i,j,k,l")  = R4("d,b,i,j,k,l") * reused_["27_ov"]("n,d") * this->T1_("a,n");
    sigmaR4("a,b,i,j,k,l") -= tmps_["perm_Lvvoooo"]("a,b,i,j,k,l");
    sigmaR4("a,b,i,j,k,l") += tmps_["perm_Lvvoooo"]("b,a,i,j,k,l");
    TAmanager.free("vvoooo", std::move(tmps_["perm_Lvvoooo"]));
    
    // sigmaR4 += +0.50 P(a,b) <m,a||c,d> R4(c,d,i,j,k,l) this->T1_(b,m) 
    // flops: o4v2L1 += o5v3L1 o5v2L1
    //  mems: o4v2L1 += o5v1L1 o4v2L1
    tmps_["perm_Lvvoooo"] = TAmanager.malloc<MatsT>("vvoooo");
    tmps_["perm_Lvvoooo"]("a,b,i,j,k,l")  = 0.50 * conj(this->antiSymMoints["vvvo"]("c,d,a,m")) * R4("c,d,i,j,k,l") * this->T1_("b,m");
    sigmaR4("a,b,i,j,k,l") -= tmps_["perm_Lvvoooo"]("a,b,i,j,k,l");
    sigmaR4("a,b,i,j,k,l") += tmps_["perm_Lvvoooo"]("b,a,i,j,k,l");
    TAmanager.free("vvoooo", std::move(tmps_["perm_Lvvoooo"]));
    
    // flops: o3v1L1  = o3v1L1 o3v1L1 o4v2L1 o3v1L1 o4v2L1 o3v1L1 o4v2L1 o3v1L1 o4v1L1 o3v1L1 o3v1L1 o3v2L1 o3v1L1 o4v1L1 o3v1L1 o5v1L1 o4v1L1 o3v1L1 o3v1L1 o3v1L1 o5v2L1 o3v1L1 o5v1L1 o4v1L1 o3v1L1 o3v1L1 o3v2L1 o3v1L1 o4v1L1 o3v1L1 o5v2L1 o3v1L1 o3v1L1
    //  mems: o3v1L1  = o3v1L1 o3v1L1 o3v1L1 o3v1L1 o3v1L1 o3v1L1 o3v1L1 o3v1L1 o3v1L1 o3v1L1 o1v1L1 o3v1L1 o3v1L1 o3v1L1 o3v1L1 o4v0L1 o3v1L1 o3v1L1 o3v1L1 o3v1L1 o3v1L1 o3v1L1 o4v0L1 o3v1L1 o3v1L1 o1v1L1 o3v1L1 o3v1L1 o3v1L1 o3v1L1 o3v1L1 o3v1L1 o3v1L1
    tmps_.emplace(std::make_pair("1_vooo", TAmanager.malloc<MatsT>("vooo")));
    tmps_["1_vooo"]("a,i,j,k")  = R2("j,k") * reused_["92_vo"]("a,i");
    tmps_["1_vooo"]("a,i,j,k") -= this->fockMatrix_ta["vo"]("a,i") * R2("j,k");
    tmps_["1_vooo"]("a,i,j,k") -= reused_["21_vovo"]("a,l,c,i") * R3("c,j,k,l");
    tmps_["1_vooo"]("a,i,j,k") += reused_["14_voov"]("a,i,m,c") * R3("c,j,k,m");
    tmps_["1_vooo"]("a,i,j,k") += this->antiSymMoints["vovo"]("a,l,b,i") * R3("b,j,k,l");
    tmps_["1_vooo"]("a,i,j,k") += R3("a,j,k,m") * reused_["44_oo"]("m,i");
    tmps_["1_vooo"]("a,i,j,k") -= 0.50 * conj(this->antiSymMoints["vooo"]("b,i,l,m")) * R2("m,l") * this->T2_("b,a,j,k");
    tmps_["1_vooo"]("a,i,j,k") += this->fockMatrix_ta["oo"]("l,i") * R3("a,j,k,l");
    tmps_["1_vooo"]("a,i,j,k") += reused_["24_ooov"]("i,l,m,c") * R3("c,j,k,m") * this->T1_("a,l");
    tmps_["1_vooo"]("a,i,j,k") += reused_["41_vo"]("a,i") * R2("j,k");
    tmps_["1_vooo"]("a,i,j,k") -= 0.50 * conj(this->antiSymMoints["vooo"]("b,i,l,m")) * R4("b,a,j,k,m,l");
    tmps_["1_vooo"]("a,i,j,k") -= conj(this->antiSymMoints["vooo"]("b,i,l,m")) * R3("b,j,k,m") * this->T1_("a,l");
    tmps_["1_vooo"]("a,i,j,k") += 0.50 * R2("m,l") * reused_["24_ooov"]("i,l,m,c") * this->T2_("c,a,j,k");
    tmps_["1_vooo"]("a,i,j,k") += R3("a,j,k,l") * reused_["32_oo"]("i,l");
    tmps_["1_vooo"]("a,i,j,k") += 0.50 * R4("c,a,j,k,m,l") * reused_["24_ooov"]("i,l,m,c");
    
    // sigmaR3 += +1.00 <m,l||b,c> R2(j,k) this->T2_(c,a,i,m) this->T1_(b,l) 
    //              += +1.00 <l,a||b,c> R2(j,k) this->T1_(b,l) this->T1_(c,i) 
    //              += +1.00 f(a,i) R2(j,k) 
    //              += -1.00 <l,a||b,c> R3(c,j,k,l) this->T1_(b,i) 
    //              += +1.00 <m,l||b,c> R3(c,j,k,m) this->T2_(b,a,i,l) 
    //              += +1.00 <l,a||b,i> R3(b,j,k,l) 
    //              += +1.00 <m,l||b,i> R3(a,j,k,m) this->T1_(b,l) 
    //              += -0.50 <m,l||b,c> R3(a,j,k,m) this->T2_(b,c,i,l) 
    //              += +1.00 <m,l||b,c> R3(a,j,k,m) this->T1_(b,l) this->T1_(c,i) 
    //              += -0.50 <m,l||b,i> R2(m,l) this->T2_(b,a,j,k) 
    //              += -1.00 f(l,i) R3(a,j,k,l) 
    //              += +1.00 <m,l||b,c> R3(c,j,k,m) this->T1_(a,l) this->T1_(b,i) 
    //              += +1.00 <m,l||b,c> R2(j,k) this->T1_(a,m) this->T1_(b,l) this->T1_(c,i) 
    //              += +0.50 <m,l||b,c> R2(j,k) this->T2_(c,a,m,l) this->T1_(b,i) 
    //              += +1.00 <l,a||b,i> R2(j,k) this->T1_(b,l) 
    //              += +1.00 f(a,b) R2(j,k) this->T1_(b,i) 
    //              += -1.00 f(l,b) R2(j,k) this->T2_(b,a,i,l) 
    //              += -1.00 f(l,i) R2(j,k) this->T1_(a,l) 
    //              += +0.25 <m,l||b,c> R2(j,k) this->T3_(b,c,a,i,m,l) 
    //              += -0.50 <l,a||b,c> R2(j,k) this->T2_(b,c,i,l) 
    //              += -0.50 <m,l||b,i> R2(j,k) this->T2_(b,a,m,l) 
    //              += +0.50 <m,l||b,c> R2(j,k) this->T1_(a,l) this->T2_(b,c,i,m) 
    //              += +1.00 <m,l||b,i> R2(j,k) this->T1_(a,m) this->T1_(b,l) 
    //              += -1.00 f(l,b) R2(j,k) this->T1_(a,l) this->T1_(b,i) 
    //              += -0.50 <m,l||b,i> R4(b,a,j,k,m,l) 
    //              += -1.00 <m,l||b,i> R3(b,j,k,m) this->T1_(a,l) 
    //              += +0.50 <m,l||b,c> R2(m,l) this->T2_(c,a,j,k) this->T1_(b,i) 
    //              += -1.00 f(l,b) R3(a,j,k,l) this->T1_(b,i) 
    //              += +0.50 <m,l||b,c> R4(c,a,j,k,m,l) this->T1_(b,i) 
    sigmaR3("a,i,j,k") -= tmps_["1_vooo"]("a,i,j,k");
    
    // sigmaR3 += +1.00 P(j,k) <m,l||b,c> R2(i,j) this->T2_(c,a,k,m) this->T1_(b,l) 
    //              += +1.00 P(j,k) <l,a||b,c> R2(i,j) this->T1_(b,l) this->T1_(c,k) 
    //              += +1.00 P(j,k) f(a,k) R2(i,j) 
    //              += -1.00 P(j,k) <l,a||b,c> R3(c,i,j,l) this->T1_(b,k) 
    //              += +1.00 P(j,k) <m,l||b,c> R3(c,i,j,m) this->T2_(b,a,k,l) 
    //              += +1.00 P(j,k) <l,a||b,k> R3(b,i,j,l) 
    //              += +1.00 P(j,k) <m,l||b,k> R3(a,i,j,m) this->T1_(b,l) 
    //              += -0.50 P(j,k) <m,l||b,c> R3(a,i,j,m) this->T2_(b,c,k,l) 
    //              += +1.00 P(j,k) <m,l||b,c> R3(a,i,j,m) this->T1_(b,l) this->T1_(c,k) 
    //              += -0.50 P(j,k) <m,l||b,k> R2(m,l) this->T2_(b,a,i,j) 
    //              += -1.00 P(j,k) f(l,k) R3(a,i,j,l) 
    //              += +1.00 P(j,k) <m,l||b,c> R3(c,i,j,m) this->T1_(a,l) this->T1_(b,k) 
    //              += +1.00 P(j,k) <m,l||b,c> R2(i,j) this->T1_(a,m) this->T1_(b,l) this->T1_(c,k) 
    //              += +0.50 P(j,k) <m,l||b,c> R2(i,j) this->T2_(c,a,m,l) this->T1_(b,k) 
    //              += +1.00 P(j,k) <l,a||b,k> R2(i,j) this->T1_(b,l) 
    //              += +1.00 P(j,k) f(a,b) R2(i,j) this->T1_(b,k) 
    //              += -1.00 P(j,k) f(l,b) R2(i,j) this->T2_(b,a,k,l) 
    //              += -1.00 P(j,k) f(l,k) R2(i,j) this->T1_(a,l) 
    //              += +0.25 P(j,k) <m,l||b,c> R2(i,j) this->T3_(b,c,a,k,m,l) 
    //              += -0.50 P(j,k) <l,a||b,c> R2(i,j) this->T2_(b,c,k,l) 
    //              += -0.50 P(j,k) <m,l||b,k> R2(i,j) this->T2_(b,a,m,l) 
    //              += +0.50 P(j,k) <m,l||b,c> R2(i,j) this->T1_(a,l) this->T2_(b,c,k,m) 
    //              += +1.00 P(j,k) <m,l||b,k> R2(i,j) this->T1_(a,m) this->T1_(b,l) 
    //              += -1.00 P(j,k) f(l,b) R2(i,j) this->T1_(a,l) this->T1_(b,k) 
    //              += -0.50 P(j,k) <m,l||b,k> R4(b,a,i,j,m,l) 
    //              += -1.00 P(j,k) <m,l||b,k> R3(b,i,j,m) this->T1_(a,l) 
    //              += +0.50 P(j,k) <m,l||b,c> R2(m,l) this->T2_(c,a,i,j) this->T1_(b,k) 
    //              += -1.00 P(j,k) f(l,b) R3(a,i,j,l) this->T1_(b,k) 
    //              += +0.50 P(j,k) <m,l||b,c> R4(c,a,i,j,m,l) this->T1_(b,k) 
    sigmaR3("a,i,j,k") -= tmps_["1_vooo"]("a,k,i,j");
    sigmaR3("a,i,j,k") += tmps_["1_vooo"]("a,j,i,k");
    TAmanager.free("vooo", std::move(tmps_["1_vooo"]));
    
    // flops: o3v1L1  = o4v1L1 o4v1L1 o3v1L1 o4v1L1 o3v1L1 o3v2L1 o3v2L1 o3v1L1 o5v1L1 o3v1L1 o5v1L1 o3v1L1 o5v1L1 o3v1L1 o4v1L1 o3v1L1 o4v1L1 o3v1L1 o4v1L1 o3v1L1
    //  mems: o3v1L1  = o3v1L1 o3v1L1 o3v1L1 o3v1L1 o3v1L1 o1v1L1 o3v1L1 o3v1L1 o3v1L1 o3v1L1 o3v1L1 o3v1L1 o3v1L1 o3v1L1 o3v1L1 o3v1L1 o3v1L1 o3v1L1 o3v1L1 o3v1L1
    tmps_.emplace(std::make_pair("2_vooo", TAmanager.malloc<MatsT>("vooo")));
    tmps_["2_vooo"]("a,i,j,k")  = reused_["17_vooo"]("a,l,i,j") * R2("k,l");
    tmps_["2_vooo"]("a,i,j,k") += 2.00 * reused_["46_vooo"]("a,i,j,m") * R2("k,m");
    tmps_["2_vooo"]("a,i,j,k") -= 2.00 * reused_["22_vooo"]("a,l,j,i") * R2("k,l");
    tmps_["2_vooo"]("a,i,j,k") += conj(this->antiSymMoints["vvoo"]("b,c,l,m")) * R3("c,k,m,l") * this->T2_("b,a,i,j");
    tmps_["2_vooo"]("a,i,j,k") -= 0.50 * R3("a,k,m,l") * reused_["7_oooo"]("i,j,l,m");
    tmps_["2_vooo"]("a,i,j,k") -= this->antiSymMoints["oooo"]("l,m,i,j") * R3("a,k,m,l");
    tmps_["2_vooo"]("a,i,j,k") += R3("a,k,m,l") * reused_["45_oooo"]("i,j,l,m");
    tmps_["2_vooo"]("a,i,j,k") -= 2.00 * reused_["60_vooo"]("a,i,j,l") * R2("k,l");
    tmps_["2_vooo"]("a,i,j,k") += 2.00 * this->antiSymMoints["vooo"]("a,l,i,j") * R2("k,l");
    tmps_["2_vooo"]("a,i,j,k") += R2("k,m") * reused_["88_vooo"]("a,i,j,m");
    
    // sigmaR3 += -0.50 P(i,j) <l,a||b,c> R2(i,l) this->T2_(b,c,j,k) 
    //              += -1.00 P(i,j) <m,l||b,c> R2(i,m) this->T1_(a,l) this->T1_(b,k) this->T1_(c,j) 
    //              += +1.00 P(i,j) <l,a||b,c> R2(i,l) this->T1_(b,k) this->T1_(c,j) 
    //              += -0.50 P(i,j) <m,l||b,c> R3(c,i,m,l) this->T2_(b,a,j,k) 
    //              += +0.25 P(i,j) <m,l||b,c> R3(a,i,m,l) this->T2_(b,c,j,k) 
    //              += +0.50 P(i,j) <m,l||j,k> R3(a,i,m,l) 
    //              += -0.50 P(i,j) <m,l||b,c> R3(a,i,m,l) this->T1_(b,k) this->T1_(c,j) 
    //              += -1.00 P(i,j) f(l,b) R2(i,l) this->T2_(b,a,j,k) 
    //              += -1.00 P(i,j) <l,a||j,k> R2(i,l) 
    //              += +0.50 P(i,j) <m,l||b,c> R2(i,m) this->T1_(a,l) this->T2_(b,c,j,k) 
    //              += +1.00 P(i,j) <m,l||b,c> R2(i,m) this->T2_(c,a,j,k) this->T1_(b,l) 
    //              += +0.50 P(i,j) <m,l||b,c> R2(i,m) this->T3_(b,c,a,j,k,l) 
    //              += +1.00 P(i,j) <m,l||j,k> R2(i,m) this->T1_(a,l) 
    sigmaR3("a,i,j,k") += 0.50 * tmps_["2_vooo"]("a,j,k,i");
    sigmaR3("a,i,j,k") -= 0.50 * tmps_["2_vooo"]("a,i,k,j");
    
    // sigmaR3 += -0.50 <l,a||b,c> R2(k,l) this->T2_(b,c,i,j) 
    //              += -1.00 <m,l||b,c> R2(k,m) this->T1_(a,l) this->T1_(b,j) this->T1_(c,i) 
    //              += +1.00 <l,a||b,c> R2(k,l) this->T1_(b,j) this->T1_(c,i) 
    //              += -0.50 <m,l||b,c> R3(c,k,m,l) this->T2_(b,a,i,j) 
    //              += +0.25 <m,l||b,c> R3(a,k,m,l) this->T2_(b,c,i,j) 
    //              += +0.50 <m,l||i,j> R3(a,k,m,l) 
    //              += -0.50 <m,l||b,c> R3(a,k,m,l) this->T1_(b,j) this->T1_(c,i) 
    //              += -1.00 f(l,b) R2(k,l) this->T2_(b,a,i,j) 
    //              += -1.00 <l,a||i,j> R2(k,l) 
    //              += +0.50 <m,l||b,c> R2(k,m) this->T1_(a,l) this->T2_(b,c,i,j) 
    //              += +1.00 <m,l||b,c> R2(k,m) this->T2_(c,a,i,j) this->T1_(b,l) 
    //              += +0.50 <m,l||b,c> R2(k,m) this->T3_(b,c,a,i,j,l) 
    //              += +1.00 <m,l||i,j> R2(k,m) this->T1_(a,l) 
    sigmaR3("a,i,j,k") += 0.50 * tmps_["2_vooo"]("a,i,j,k");
    TAmanager.free("vooo", std::move(tmps_["2_vooo"]));
    
    // flops: o3v1L1  = o4v1L1 o5v1L1 o3v1L1
    //  mems: o3v1L1  = o3v1L1 o3v1L1 o3v1L1
    tmps_.emplace(std::make_pair("3_oovo", TAmanager.malloc<MatsT>("oovo")));
    tmps_["3_oovo"]("i,j,a,k")  = R2("i,m") * reused_["55_oovo"]("m,j,a,k");
    tmps_["3_oovo"]("i,j,a,k") += 0.50 * R3("a,i,m,l") * reused_["53_oooo"]("k,l,m,j");
    
    // sigmaR3 += -1.00 P(i,k) <m,l||b,j> R2(i,m) this->T2_(b,a,k,l) 
    //              += -1.00 P(i,k) <m,l||b,j> R2(i,m) this->T1_(a,l) this->T1_(b,k) 
    //              += -0.50 P(i,k) <m,l||b,j> R3(a,i,m,l) this->T1_(b,k) 
    sigmaR3("a,i,j,k") += tmps_["3_oovo"]("i,j,a,k");
    sigmaR3("a,i,j,k") -= tmps_["3_oovo"]("k,j,a,i");
    
    // sigmaR3 += +1.00 P(i,j) <m,l||b,k> R2(i,m) this->T2_(b,a,j,l) 
    //              += +1.00 P(i,j) <m,l||b,k> R2(i,m) this->T1_(a,l) this->T1_(b,j) 
    //              += +0.50 P(i,j) <m,l||b,k> R3(a,i,m,l) this->T1_(b,j) 
    sigmaR3("a,i,j,k") -= tmps_["3_oovo"]("i,k,a,j");
    sigmaR3("a,i,j,k") += tmps_["3_oovo"]("j,k,a,i");
    
    // sigmaR3 += +1.00 P(j,k) <m,l||b,i> R2(j,m) this->T2_(b,a,k,l) 
    //              += +1.00 P(j,k) <m,l||b,i> R2(j,m) this->T1_(a,l) this->T1_(b,k) 
    //              += +0.50 P(j,k) <m,l||b,i> R3(a,j,m,l) this->T1_(b,k) 
    sigmaR3("a,i,j,k") -= tmps_["3_oovo"]("j,i,a,k");
    sigmaR3("a,i,j,k") += tmps_["3_oovo"]("k,i,a,j");
    TAmanager.free("oovo", std::move(tmps_["3_oovo"]));
    
    // flops: o3v1L1  = o4v1L1 o4v1L1 o3v1L1
    //  mems: o3v1L1  = o3v1L1 o3v1L1 o3v1L1
    tmps_.emplace(std::make_pair("4_vooo", TAmanager.malloc<MatsT>("vooo")));
    tmps_["4_vooo"]("a,i,k,j")  = R2("j,l") * reused_["19_vooo"]("a,l,i,k");
    tmps_["4_vooo"]("a,i,k,j") += reused_["6_oovo"]("i,m,a,k") * R2("j,m");
    
    // sigmaR3 += -1.00 P(i,j) <l,a||b,k> R2(i,l) this->T1_(b,j) 
    //              += -1.00 P(i,j) <m,l||b,c> R2(i,m) this->T2_(c,a,j,l) this->T1_(b,k) 
    sigmaR3("a,i,j,k") += tmps_["4_vooo"]("a,k,j,i");
    sigmaR3("a,i,j,k") -= tmps_["4_vooo"]("a,k,i,j");
    
    // sigmaR3 += -1.00 P(j,k) <l,a||b,i> R2(j,l) this->T1_(b,k) 
    //              += -1.00 P(j,k) <m,l||b,c> R2(j,m) this->T2_(c,a,k,l) this->T1_(b,i) 
    sigmaR3("a,i,j,k") += tmps_["4_vooo"]("a,i,k,j");
    sigmaR3("a,i,j,k") -= tmps_["4_vooo"]("a,i,j,k");
    
    // sigmaR3 += +1.00 P(i,k) <l,a||b,j> R2(i,l) this->T1_(b,k) 
    //              += +1.00 P(i,k) <m,l||b,c> R2(i,m) this->T2_(c,a,k,l) this->T1_(b,j) 
    sigmaR3("a,i,j,k") -= tmps_["4_vooo"]("a,j,k,i");
    sigmaR3("a,i,j,k") += tmps_["4_vooo"]("a,j,i,k");
    TAmanager.free("vooo", std::move(tmps_["4_vooo"]));
    
    // flops: o4v0L1  = o5v2L1
    //  mems: o4v0L1  = o4v0L1
    tmps_.emplace(std::make_pair("5_oooo", TAmanager.malloc<MatsT>("oooo")));
    tmps_["5_oooo"]("i,j,k,l")  = conj(this->antiSymMoints["vvoo"]("b,c,l,m")) * R4("b,c,i,j,k,m");
    
    // sigmaR3 += +0.50 <m,l||b,c> R4(b,c,i,j,k,m) this->T1_(a,l) 
    // flops: o3v1L1 += o4v1L1
    //  mems: o3v1L1 += o3v1L1
    sigmaR3("a,i,j,k") -= 0.50 * this->T1_("a,l") * tmps_["5_oooo"]("i,j,k,l");
    
    // flops: o4v2L1  = o2v2 o4v2L1 o5v2L1 o5v2L1 o4v2L1 o4v2L1
    //  mems: o4v2L1  = o2v2 o4v2L1 o4v2L1 o4v2L1 o4v2L1 o4v2L1
    tmps_.emplace(std::make_pair("6_voovoo", TAmanager.malloc<MatsT>("voovoo")));
    tmps_["6_voovoo"]("a,i,l,b,j,k")  = (reused_["20_voov"]("a,i,l,b") + reused_["107_ovov"]("i,b,l,a")) * R2("j,k");
    tmps_["6_voovoo"]("a,i,l,b,j,k") += R3("b,j,k,m") * reused_["19_vooo"]("a,m,i,l");
    tmps_["6_voovoo"]("a,i,l,b,j,k") += R3("b,j,k,n") * reused_["6_oovo"]("i,n,a,l");
    
    // sigmaR4 += -1.00 P(k,l) P(a,b) <m,a||c,i> R2(j,k) this->T1_(b,m) this->T1_(c,l) 
    //                += -1.00 P(k,l) P(a,b) <n,m||c,d> R2(j,k) this->T1_(a,m) this->T2_(d,b,l,n) this->T1_(c,i) 
    //                += -1.00 P(k,l) P(a,b) <m,a||c,i> R2(j,k) this->T2_(c,b,l,m) 
    //                += -1.00 P(k,l) P(a,b) <m,a||c,i> R3(b,j,k,m) this->T1_(c,l) 
    //                += -1.00 P(k,l) P(a,b) <n,m||c,d> R3(b,j,k,n) this->T2_(d,a,l,m) this->T1_(c,i) 
    sigmaR4("a,b,i,j,k,l") += tmps_["6_voovoo"]("a,i,l,b,j,k");
    sigmaR4("a,b,i,j,k,l") -= tmps_["6_voovoo"]("a,i,k,b,j,l");
    sigmaR4("a,b,i,j,k,l") -= tmps_["6_voovoo"]("b,i,l,a,j,k");
    sigmaR4("a,b,i,j,k,l") += tmps_["6_voovoo"]("b,i,k,a,j,l");
    
    // sigmaR4 += -1.00 P(j,l) P(a,b) <m,a||c,k> R2(i,j) this->T1_(b,m) this->T1_(c,l) 
    //                += -1.00 P(j,l) P(a,b) <n,m||c,d> R2(i,j) this->T1_(a,m) this->T2_(d,b,l,n) this->T1_(c,k) 
    //                += -1.00 P(j,l) P(a,b) <m,a||c,k> R2(i,j) this->T2_(c,b,l,m) 
    //                += -1.00 P(j,l) P(a,b) <m,a||c,k> R3(b,i,j,m) this->T1_(c,l) 
    //                += -1.00 P(j,l) P(a,b) <n,m||c,d> R3(b,i,j,n) this->T2_(d,a,l,m) this->T1_(c,k) 
    sigmaR4("a,b,i,j,k,l") += tmps_["6_voovoo"]("a,k,l,b,i,j");
    sigmaR4("a,b,i,j,k,l") -= tmps_["6_voovoo"]("a,k,j,b,i,l");
    sigmaR4("a,b,i,j,k,l") -= tmps_["6_voovoo"]("b,k,l,a,i,j");
    sigmaR4("a,b,i,j,k,l") += tmps_["6_voovoo"]("b,k,j,a,i,l");
    
    // sigmaR4 += +1.00 P(k,l) P(a,b) <m,a||c,j> R2(i,k) this->T1_(b,m) this->T1_(c,l) 
    //                += +1.00 P(k,l) P(a,b) <n,m||c,d> R2(i,k) this->T1_(a,m) this->T2_(d,b,l,n) this->T1_(c,j) 
    //                += +1.00 P(k,l) P(a,b) <m,a||c,j> R2(i,k) this->T2_(c,b,l,m) 
    //                += +1.00 P(k,l) P(a,b) <m,a||c,j> R3(b,i,k,m) this->T1_(c,l) 
    //                += +1.00 P(k,l) P(a,b) <n,m||c,d> R3(b,i,k,n) this->T2_(d,a,l,m) this->T1_(c,j) 
    sigmaR4("a,b,i,j,k,l") -= tmps_["6_voovoo"]("a,j,l,b,i,k");
    sigmaR4("a,b,i,j,k,l") += tmps_["6_voovoo"]("a,j,k,b,i,l");
    sigmaR4("a,b,i,j,k,l") += tmps_["6_voovoo"]("b,j,l,a,i,k");
    sigmaR4("a,b,i,j,k,l") -= tmps_["6_voovoo"]("b,j,k,a,i,l");
    
    // sigmaR4 += +1.00 P(k,l) P(a,b) <m,a||c,l> R2(j,k) this->T1_(b,m) this->T1_(c,i) 
    //                += +1.00 P(k,l) P(a,b) <n,m||c,d> R2(j,k) this->T1_(a,m) this->T2_(d,b,i,n) this->T1_(c,l) 
    //                += +1.00 P(k,l) P(a,b) <m,a||c,l> R2(j,k) this->T2_(c,b,i,m) 
    //                += +1.00 P(k,l) P(a,b) <m,a||c,l> R3(b,j,k,m) this->T1_(c,i) 
    //                += +1.00 P(k,l) P(a,b) <n,m||c,d> R3(b,j,k,n) this->T2_(d,a,i,m) this->T1_(c,l) 
    sigmaR4("a,b,i,j,k,l") -= tmps_["6_voovoo"]("a,l,i,b,j,k");
    sigmaR4("a,b,i,j,k,l") += tmps_["6_voovoo"]("a,k,i,b,j,l");
    sigmaR4("a,b,i,j,k,l") += tmps_["6_voovoo"]("b,l,i,a,j,k");
    sigmaR4("a,b,i,j,k,l") -= tmps_["6_voovoo"]("b,k,i,a,j,l");
    
    // sigmaR4 += +1.00 P(i,j) P(a,b) <m,a||c,j> R2(k,l) this->T1_(b,m) this->T1_(c,i) 
    //                += +1.00 P(i,j) P(a,b) <n,m||c,d> R2(k,l) this->T1_(a,m) this->T2_(d,b,i,n) this->T1_(c,j) 
    //                += +1.00 P(i,j) P(a,b) <m,a||c,j> R2(k,l) this->T2_(c,b,i,m) 
    //                += +1.00 P(i,j) P(a,b) <m,a||c,j> R3(b,k,l,m) this->T1_(c,i) 
    //                += +1.00 P(i,j) P(a,b) <n,m||c,d> R3(b,k,l,n) this->T2_(d,a,i,m) this->T1_(c,j) 
    sigmaR4("a,b,i,j,k,l") -= tmps_["6_voovoo"]("a,j,i,b,k,l");
    sigmaR4("a,b,i,j,k,l") += tmps_["6_voovoo"]("a,i,j,b,k,l");
    sigmaR4("a,b,i,j,k,l") += tmps_["6_voovoo"]("b,j,i,a,k,l");
    sigmaR4("a,b,i,j,k,l") -= tmps_["6_voovoo"]("b,i,j,a,k,l");
    
    // sigmaR4 += +1.00 P(j,k) P(a,b) <m,a||c,l> R2(i,j) this->T1_(b,m) this->T1_(c,k) 
    //                += +1.00 P(j,k) P(a,b) <n,m||c,d> R2(i,j) this->T1_(a,m) this->T2_(d,b,k,n) this->T1_(c,l) 
    //                += +1.00 P(j,k) P(a,b) <m,a||c,l> R2(i,j) this->T2_(c,b,k,m) 
    //                += +1.00 P(j,k) P(a,b) <m,a||c,l> R3(b,i,j,m) this->T1_(c,k) 
    //                += +1.00 P(j,k) P(a,b) <n,m||c,d> R3(b,i,j,n) this->T2_(d,a,k,m) this->T1_(c,l) 
    sigmaR4("a,b,i,j,k,l") -= tmps_["6_voovoo"]("a,l,k,b,i,j");
    sigmaR4("a,b,i,j,k,l") += tmps_["6_voovoo"]("a,l,j,b,i,k");
    sigmaR4("a,b,i,j,k,l") += tmps_["6_voovoo"]("b,l,k,a,i,j");
    sigmaR4("a,b,i,j,k,l") -= tmps_["6_voovoo"]("b,l,j,a,i,k");
    TAmanager.free("voovoo", std::move(tmps_["6_voovoo"]));
    
    // flops: o4v2L1  = o5v2L1 o4v2L1 o4v3L1 o4v2L1 o5v2L1 o4v2L1
    //  mems: o4v2L1  = o4v2L1 o2v2L1 o4v2L1 o4v2L1 o4v2L1 o4v2L1
    tmps_.emplace(std::make_pair("7_vovooo", TAmanager.malloc<MatsT>("vovooo")));
    tmps_["7_vovooo"]("a,l,b,i,j,k")  = reused_["117_voovoo"]("a,m,l,b,i,j") * R2("k,m");
    tmps_["7_vovooo"]("a,l,b,i,j,k") += 0.50 * R3("b,k,n,m") * reused_["24_ooov"]("l,m,n,d") * this->T2_("d,a,i,j");
    tmps_["7_vovooo"]("a,l,b,i,j,k") += R2("k,n") * reused_["112_oovoov"]("n,l,b,i,j,a");
    
    // sigmaR4 += -1.00 P(i,k) P(a,b) <m,a||c,d> R2(i,m) this->T2_(d,b,k,l) this->T1_(c,j) 
    //                += -0.50 P(i,k) P(a,b) <n,m||c,d> R3(b,i,n,m) this->T2_(d,a,k,l) this->T1_(c,j) 
    //                += -1.00 P(i,k) P(a,b) <n,m||c,j> R2(i,n) this->T1_(a,m) this->T2_(c,b,k,l) 
    sigmaR4("a,b,i,j,k,l") += tmps_["7_vovooo"]("a,j,b,k,l,i");
    sigmaR4("a,b,i,j,k,l") -= tmps_["7_vovooo"]("a,j,b,i,l,k");
    sigmaR4("a,b,i,j,k,l") -= tmps_["7_vovooo"]("b,j,a,k,l,i");
    sigmaR4("a,b,i,j,k,l") += tmps_["7_vovooo"]("b,j,a,i,l,k");
    
    // sigmaR4 += +1.00 P(i,j) P(a,b) <m,a||c,d> R2(i,m) this->T2_(d,b,j,l) this->T1_(c,k) 
    //                += +0.50 P(i,j) P(a,b) <n,m||c,d> R3(b,i,n,m) this->T2_(d,a,j,l) this->T1_(c,k) 
    //                += +1.00 P(i,j) P(a,b) <n,m||c,k> R2(i,n) this->T1_(a,m) this->T2_(c,b,j,l) 
    sigmaR4("a,b,i,j,k,l") -= tmps_["7_vovooo"]("a,k,b,j,l,i");
    sigmaR4("a,b,i,j,k,l") += tmps_["7_vovooo"]("a,k,b,i,l,j");
    sigmaR4("a,b,i,j,k,l") += tmps_["7_vovooo"]("b,k,a,j,l,i");
    sigmaR4("a,b,i,j,k,l") -= tmps_["7_vovooo"]("b,k,a,i,l,j");
    
    // sigmaR4 += -1.00 P(i,j) P(a,b) <m,a||c,d> R2(l,m) this->T2_(d,b,i,k) this->T1_(c,j) 
    //                += -0.50 P(i,j) P(a,b) <n,m||c,d> R3(b,l,n,m) this->T2_(d,a,i,k) this->T1_(c,j) 
    //                += -1.00 P(i,j) P(a,b) <n,m||c,j> R2(l,n) this->T1_(a,m) this->T2_(c,b,i,k) 
    sigmaR4("a,b,i,j,k,l") += tmps_["7_vovooo"]("a,j,b,i,k,l");
    sigmaR4("a,b,i,j,k,l") -= tmps_["7_vovooo"]("a,i,b,j,k,l");
    sigmaR4("a,b,i,j,k,l") -= tmps_["7_vovooo"]("b,j,a,i,k,l");
    sigmaR4("a,b,i,j,k,l") += tmps_["7_vovooo"]("b,i,a,j,k,l");
    
    // sigmaR4 += +1.00 P(j,k) P(a,b) <m,a||c,d> R2(j,m) this->T2_(d,b,k,l) this->T1_(c,i) 
    //                += +0.50 P(j,k) P(a,b) <n,m||c,d> R3(b,j,n,m) this->T2_(d,a,k,l) this->T1_(c,i) 
    //                += +1.00 P(j,k) P(a,b) <n,m||c,i> R2(j,n) this->T1_(a,m) this->T2_(c,b,k,l) 
    sigmaR4("a,b,i,j,k,l") -= tmps_["7_vovooo"]("a,i,b,k,l,j");
    sigmaR4("a,b,i,j,k,l") += tmps_["7_vovooo"]("a,i,b,j,l,k");
    sigmaR4("a,b,i,j,k,l") += tmps_["7_vovooo"]("b,i,a,k,l,j");
    sigmaR4("a,b,i,j,k,l") -= tmps_["7_vovooo"]("b,i,a,j,l,k");
    
    // sigmaR4 += -1.00 P(k,l) P(a,b) <m,a||c,d> R2(k,m) this->T2_(d,b,i,j) this->T1_(c,l) 
    //                += -0.50 P(k,l) P(a,b) <n,m||c,d> R3(b,k,n,m) this->T2_(d,a,i,j) this->T1_(c,l) 
    //                += -1.00 P(k,l) P(a,b) <n,m||c,l> R2(k,n) this->T1_(a,m) this->T2_(c,b,i,j) 
    sigmaR4("a,b,i,j,k,l") += tmps_["7_vovooo"]("a,l,b,i,j,k");
    sigmaR4("a,b,i,j,k,l") -= tmps_["7_vovooo"]("a,k,b,i,j,l");
    sigmaR4("a,b,i,j,k,l") -= tmps_["7_vovooo"]("b,l,a,i,j,k");
    sigmaR4("a,b,i,j,k,l") += tmps_["7_vovooo"]("b,k,a,i,j,l");
    
    // sigmaR4 += -1.00 P(i,j) P(a,b) <m,a||c,d> R2(i,m) this->T2_(d,b,j,k) this->T1_(c,l) 
    //                += -0.50 P(i,j) P(a,b) <n,m||c,d> R3(b,i,n,m) this->T2_(d,a,j,k) this->T1_(c,l) 
    //                += -1.00 P(i,j) P(a,b) <n,m||c,l> R2(i,n) this->T1_(a,m) this->T2_(c,b,j,k) 
    sigmaR4("a,b,i,j,k,l") += tmps_["7_vovooo"]("a,l,b,j,k,i");
    sigmaR4("a,b,i,j,k,l") -= tmps_["7_vovooo"]("a,l,b,i,k,j");
    sigmaR4("a,b,i,j,k,l") -= tmps_["7_vovooo"]("b,l,a,j,k,i");
    sigmaR4("a,b,i,j,k,l") += tmps_["7_vovooo"]("b,l,a,i,k,j");
    TAmanager.free("vovooo", std::move(tmps_["7_vovooo"]));
    
    // flops: o4v2L1  = o4v2L1 o4v3L1 o5v2L1 o4v2L1 o5v2L1 o4v2L1
    //  mems: o4v2L1  = o2v2L1 o4v2L1 o4v2L1 o4v2L1 o4v2L1 o4v2L1
    tmps_.emplace(std::make_pair("8_ovovoo", TAmanager.malloc<MatsT>("ovovoo")));
    tmps_["8_ovovoo"]("i,b,j,a,k,l")  = R3("b,j,n,m") * conj(this->antiSymMoints["vooo"]("c,i,m,n")) * this->T2_("c,a,k,l");
    tmps_["8_ovovoo"]("i,b,j,a,k,l") += 2.00 * reused_["113_oovoov"]("i,n,b,k,l,a") * R2("j,n");
    tmps_["8_ovovoo"]("i,b,j,a,k,l") += 2.00 * R2("j,m") * reused_["116_voovoo"]("a,m,i,b,k,l");
    
    // sigmaR4 += +0.50 P(i,j) P(a,b) <n,m||c,l> R3(b,i,n,m) this->T2_(c,a,j,k) 
    //                += +1.00 P(i,j) P(a,b) <n,m||c,d> R2(i,n) this->T1_(a,m) this->T2_(d,b,j,k) this->T1_(c,l) 
    //                += +1.00 P(i,j) P(a,b) <m,a||c,l> R2(i,m) this->T2_(c,b,j,k) 
    sigmaR4("a,b,i,j,k,l") -= 0.50 * tmps_["8_ovovoo"]("l,b,i,a,j,k");
    sigmaR4("a,b,i,j,k,l") += 0.50 * tmps_["8_ovovoo"]("l,b,j,a,i,k");
    sigmaR4("a,b,i,j,k,l") += 0.50 * tmps_["8_ovovoo"]("l,a,i,b,j,k");
    sigmaR4("a,b,i,j,k,l") -= 0.50 * tmps_["8_ovovoo"]("l,a,j,b,i,k");
    
    // sigmaR4 += +0.50 P(i,j) P(a,b) <n,m||c,j> R3(b,l,n,m) this->T2_(c,a,i,k) 
    //                += +1.00 P(i,j) P(a,b) <n,m||c,d> R2(l,n) this->T1_(a,m) this->T2_(d,b,i,k) this->T1_(c,j) 
    //                += +1.00 P(i,j) P(a,b) <m,a||c,j> R2(l,m) this->T2_(c,b,i,k) 
    sigmaR4("a,b,i,j,k,l") -= 0.50 * tmps_["8_ovovoo"]("j,b,l,a,i,k");
    sigmaR4("a,b,i,j,k,l") += 0.50 * tmps_["8_ovovoo"]("i,b,l,a,j,k");
    sigmaR4("a,b,i,j,k,l") += 0.50 * tmps_["8_ovovoo"]("j,a,l,b,i,k");
    sigmaR4("a,b,i,j,k,l") -= 0.50 * tmps_["8_ovovoo"]("i,a,l,b,j,k");
    
    // sigmaR4 += +0.50 P(i,k) P(a,b) <n,m||c,j> R3(b,i,n,m) this->T2_(c,a,k,l) 
    //                += +1.00 P(i,k) P(a,b) <n,m||c,d> R2(i,n) this->T1_(a,m) this->T2_(d,b,k,l) this->T1_(c,j) 
    //                += +1.00 P(i,k) P(a,b) <m,a||c,j> R2(i,m) this->T2_(c,b,k,l) 
    sigmaR4("a,b,i,j,k,l") -= 0.50 * tmps_["8_ovovoo"]("j,b,i,a,k,l");
    sigmaR4("a,b,i,j,k,l") += 0.50 * tmps_["8_ovovoo"]("j,b,k,a,i,l");
    sigmaR4("a,b,i,j,k,l") += 0.50 * tmps_["8_ovovoo"]("j,a,i,b,k,l");
    sigmaR4("a,b,i,j,k,l") -= 0.50 * tmps_["8_ovovoo"]("j,a,k,b,i,l");
    
    // sigmaR4 += +0.50 P(k,l) P(a,b) <n,m||c,l> R3(b,k,n,m) this->T2_(c,a,i,j) 
    //                += +1.00 P(k,l) P(a,b) <n,m||c,d> R2(k,n) this->T1_(a,m) this->T2_(d,b,i,j) this->T1_(c,l) 
    //                += +1.00 P(k,l) P(a,b) <m,a||c,l> R2(k,m) this->T2_(c,b,i,j) 
    sigmaR4("a,b,i,j,k,l") -= 0.50 * tmps_["8_ovovoo"]("l,b,k,a,i,j");
    sigmaR4("a,b,i,j,k,l") += 0.50 * tmps_["8_ovovoo"]("k,b,l,a,i,j");
    sigmaR4("a,b,i,j,k,l") += 0.50 * tmps_["8_ovovoo"]("l,a,k,b,i,j");
    sigmaR4("a,b,i,j,k,l") -= 0.50 * tmps_["8_ovovoo"]("k,a,l,b,i,j");
    
    // sigmaR4 += -0.50 P(i,j) P(a,b) <n,m||c,k> R3(b,i,n,m) this->T2_(c,a,j,l) 
    //                += -1.00 P(i,j) P(a,b) <n,m||c,d> R2(i,n) this->T1_(a,m) this->T2_(d,b,j,l) this->T1_(c,k) 
    //                += -1.00 P(i,j) P(a,b) <m,a||c,k> R2(i,m) this->T2_(c,b,j,l) 
    sigmaR4("a,b,i,j,k,l") += 0.50 * tmps_["8_ovovoo"]("k,b,i,a,j,l");
    sigmaR4("a,b,i,j,k,l") -= 0.50 * tmps_["8_ovovoo"]("k,b,j,a,i,l");
    sigmaR4("a,b,i,j,k,l") -= 0.50 * tmps_["8_ovovoo"]("k,a,i,b,j,l");
    sigmaR4("a,b,i,j,k,l") += 0.50 * tmps_["8_ovovoo"]("k,a,j,b,i,l");
    
    // sigmaR4 += -0.50 P(j,k) P(a,b) <n,m||c,i> R3(b,j,n,m) this->T2_(c,a,k,l) 
    //                += -1.00 P(j,k) P(a,b) <n,m||c,d> R2(j,n) this->T1_(a,m) this->T2_(d,b,k,l) this->T1_(c,i) 
    //                += -1.00 P(j,k) P(a,b) <m,a||c,i> R2(j,m) this->T2_(c,b,k,l) 
    sigmaR4("a,b,i,j,k,l") += 0.50 * tmps_["8_ovovoo"]("i,b,j,a,k,l");
    sigmaR4("a,b,i,j,k,l") -= 0.50 * tmps_["8_ovovoo"]("i,b,k,a,j,l");
    sigmaR4("a,b,i,j,k,l") -= 0.50 * tmps_["8_ovovoo"]("i,a,j,b,k,l");
    sigmaR4("a,b,i,j,k,l") += 0.50 * tmps_["8_ovovoo"]("i,a,k,b,j,l");
    TAmanager.free("ovovoo", std::move(tmps_["8_ovovoo"]));
    
    // flops: o4v2L1  = o4v2L1 o6v2L1 o4v2L1 o5v1L1 o5v2L1 o4v2L1
    //  mems: o4v2L1  = o4v2L1 o4v2L1 o4v2L1 o4v0L1 o4v2L1 o4v2L1
    tmps_.emplace(std::make_pair("9_oovovo", TAmanager.malloc<MatsT>("oovovo")));
    tmps_["9_oovovo"]("j,k,a,i,b,l")  = reused_["78_vovo"]("a,i,b,l") * R2("j,k");
    tmps_["9_oovovo"]("j,k,a,i,b,l") += 0.50 * R4("a,b,j,k,n,m") * reused_["53_oooo"]("l,m,n,i");
    tmps_["9_oovovo"]("j,k,a,i,b,l") += reused_["24_ooov"]("i,m,n,d") * R3("d,j,k,n") * this->T2_("a,b,l,m");
    
    // sigmaR4 += +1.00 P(j,k) <n,m||c,d> R2(i,j) this->T2_(c,a,l,m) this->T2_(d,b,k,n) 
    //                += -1.00 P(j,k) f(m,l) R2(i,j) this->T2_(a,b,k,m) 
    //                += -1.00 P(j,k) f(m,c) R2(i,j) this->T2_(a,b,k,m) this->T1_(c,l) 
    //                += +0.50 P(j,k) <n,m||c,l> R2(i,j) this->T3_(c,a,b,k,n,m) 
    //                += +1.00 P(j,k) <n,m||c,d> R2(i,j) this->T2_(a,b,k,n) this->T1_(c,m) this->T1_(d,l) 
    //                += +1.00 P(j,k) <n,m||c,l> R2(i,j) this->T2_(a,b,k,n) this->T1_(c,m) 
    //                += +0.50 P(j,k) <n,m||c,l> R2(i,j) this->T2_(a,b,n,m) this->T1_(c,k) 
    //                += +0.50 P(j,k) <n,m||c,l> R4(a,b,i,j,n,m) this->T1_(c,k) 
    //                += +1.00 P(j,k) <n,m||c,d> R3(d,i,j,n) this->T2_(a,b,k,m) this->T1_(c,l) 
    sigmaR4("a,b,i,j,k,l") -= tmps_["9_oovovo"]("i,j,a,l,b,k");
    sigmaR4("a,b,i,j,k,l") += tmps_["9_oovovo"]("i,k,a,l,b,j");
    
    // sigmaR4 += -1.00 P(k,l) <n,m||c,d> R2(j,k) this->T2_(c,a,i,m) this->T2_(d,b,l,n) 
    //                += +1.00 P(k,l) f(m,i) R2(j,k) this->T2_(a,b,l,m) 
    //                += +1.00 P(k,l) f(m,c) R2(j,k) this->T2_(a,b,l,m) this->T1_(c,i) 
    //                += -0.50 P(k,l) <n,m||c,i> R2(j,k) this->T3_(c,a,b,l,n,m) 
    //                += -1.00 P(k,l) <n,m||c,d> R2(j,k) this->T2_(a,b,l,n) this->T1_(c,m) this->T1_(d,i) 
    //                += -1.00 P(k,l) <n,m||c,i> R2(j,k) this->T2_(a,b,l,n) this->T1_(c,m) 
    //                += -0.50 P(k,l) <n,m||c,i> R2(j,k) this->T2_(a,b,n,m) this->T1_(c,l) 
    //                += -0.50 P(k,l) <n,m||c,i> R4(a,b,j,k,n,m) this->T1_(c,l) 
    //                += -1.00 P(k,l) <n,m||c,d> R3(d,j,k,n) this->T2_(a,b,l,m) this->T1_(c,i) 
    sigmaR4("a,b,i,j,k,l") += tmps_["9_oovovo"]("j,k,a,i,b,l");
    sigmaR4("a,b,i,j,k,l") -= tmps_["9_oovovo"]("j,l,a,i,b,k");
    
    // sigmaR4 += +1.00 P(i,j) <n,m||c,d> R2(k,l) this->T2_(c,a,j,m) this->T2_(d,b,i,n) 
    //                += -1.00 P(i,j) f(m,j) R2(k,l) this->T2_(a,b,i,m) 
    //                += -1.00 P(i,j) f(m,c) R2(k,l) this->T2_(a,b,i,m) this->T1_(c,j) 
    //                += +0.50 P(i,j) <n,m||c,j> R2(k,l) this->T3_(c,a,b,i,n,m) 
    //                += +1.00 P(i,j) <n,m||c,d> R2(k,l) this->T2_(a,b,i,n) this->T1_(c,m) this->T1_(d,j) 
    //                += +1.00 P(i,j) <n,m||c,j> R2(k,l) this->T2_(a,b,i,n) this->T1_(c,m) 
    //                += +0.50 P(i,j) <n,m||c,j> R2(k,l) this->T2_(a,b,n,m) this->T1_(c,i) 
    //                += +0.50 P(i,j) <n,m||c,j> R4(a,b,k,l,n,m) this->T1_(c,i) 
    //                += +1.00 P(i,j) <n,m||c,d> R3(d,k,l,n) this->T2_(a,b,i,m) this->T1_(c,j) 
    sigmaR4("a,b,i,j,k,l") -= tmps_["9_oovovo"]("k,l,a,j,b,i");
    sigmaR4("a,b,i,j,k,l") += tmps_["9_oovovo"]("k,l,a,i,b,j");
    
    // sigmaR4 += -1.00 P(j,l) <n,m||c,d> R2(i,j) this->T2_(c,a,k,m) this->T2_(d,b,l,n) 
    //                += +1.00 P(j,l) f(m,k) R2(i,j) this->T2_(a,b,l,m) 
    //                += +1.00 P(j,l) f(m,c) R2(i,j) this->T2_(a,b,l,m) this->T1_(c,k) 
    //                += -0.50 P(j,l) <n,m||c,k> R2(i,j) this->T3_(c,a,b,l,n,m) 
    //                += -1.00 P(j,l) <n,m||c,d> R2(i,j) this->T2_(a,b,l,n) this->T1_(c,m) this->T1_(d,k) 
    //                += -1.00 P(j,l) <n,m||c,k> R2(i,j) this->T2_(a,b,l,n) this->T1_(c,m) 
    //                += -0.50 P(j,l) <n,m||c,k> R2(i,j) this->T2_(a,b,n,m) this->T1_(c,l) 
    //                += -0.50 P(j,l) <n,m||c,k> R4(a,b,i,j,n,m) this->T1_(c,l) 
    //                += -1.00 P(j,l) <n,m||c,d> R3(d,i,j,n) this->T2_(a,b,l,m) this->T1_(c,k) 
    sigmaR4("a,b,i,j,k,l") += tmps_["9_oovovo"]("i,j,a,k,b,l");
    sigmaR4("a,b,i,j,k,l") -= tmps_["9_oovovo"]("i,l,a,k,b,j");
    
    // sigmaR4 += +1.00 P(k,l) <n,m||c,d> R2(i,k) this->T2_(c,a,j,m) this->T2_(d,b,l,n) 
    //                += -1.00 P(k,l) f(m,j) R2(i,k) this->T2_(a,b,l,m) 
    //                += -1.00 P(k,l) f(m,c) R2(i,k) this->T2_(a,b,l,m) this->T1_(c,j) 
    //                += +0.50 P(k,l) <n,m||c,j> R2(i,k) this->T3_(c,a,b,l,n,m) 
    //                += +1.00 P(k,l) <n,m||c,d> R2(i,k) this->T2_(a,b,l,n) this->T1_(c,m) this->T1_(d,j) 
    //                += +1.00 P(k,l) <n,m||c,j> R2(i,k) this->T2_(a,b,l,n) this->T1_(c,m) 
    //                += +0.50 P(k,l) <n,m||c,j> R2(i,k) this->T2_(a,b,n,m) this->T1_(c,l) 
    //                += +0.50 P(k,l) <n,m||c,j> R4(a,b,i,k,n,m) this->T1_(c,l) 
    //                += +1.00 P(k,l) <n,m||c,d> R3(d,i,k,n) this->T2_(a,b,l,m) this->T1_(c,j) 
    sigmaR4("a,b,i,j,k,l") -= tmps_["9_oovovo"]("i,k,a,j,b,l");
    sigmaR4("a,b,i,j,k,l") += tmps_["9_oovovo"]("i,l,a,j,b,k");
    
    // sigmaR4 += +1.00 P(k,l) <n,m||c,d> R2(j,k) this->T2_(c,a,l,m) this->T2_(d,b,i,n) 
    //                += -1.00 P(k,l) f(m,l) R2(j,k) this->T2_(a,b,i,m) 
    //                += -1.00 P(k,l) f(m,c) R2(j,k) this->T2_(a,b,i,m) this->T1_(c,l) 
    //                += +0.50 P(k,l) <n,m||c,l> R2(j,k) this->T3_(c,a,b,i,n,m) 
    //                += +1.00 P(k,l) <n,m||c,d> R2(j,k) this->T2_(a,b,i,n) this->T1_(c,m) this->T1_(d,l) 
    //                += +1.00 P(k,l) <n,m||c,l> R2(j,k) this->T2_(a,b,i,n) this->T1_(c,m) 
    //                += +0.50 P(k,l) <n,m||c,l> R2(j,k) this->T2_(a,b,n,m) this->T1_(c,i) 
    //                += +0.50 P(k,l) <n,m||c,l> R4(a,b,j,k,n,m) this->T1_(c,i) 
    //                += +1.00 P(k,l) <n,m||c,d> R3(d,j,k,n) this->T2_(a,b,i,m) this->T1_(c,l) 
    sigmaR4("a,b,i,j,k,l") -= tmps_["9_oovovo"]("j,k,a,l,b,i");
    sigmaR4("a,b,i,j,k,l") += tmps_["9_oovovo"]("j,l,a,k,b,i");
    TAmanager.free("oovovo", std::move(tmps_["9_oovovo"]));
    
    // flops: o4v2L1  = o5v2L1 o3v2 o4v2L1 o5v2L1 o4v2L1 o4v2L1
    //  mems: o4v2L1  = o4v2L1 o2v2 o4v2L1 o4v2L1 o4v2L1 o4v2L1
    tmps_.emplace(std::make_pair("10_voovoo", TAmanager.malloc<MatsT>("voovoo")));
    tmps_["10_voovoo"]("a,l,i,b,j,k")  = R3("b,j,k,m") * reused_["22_vooo"]("a,m,l,i");
    tmps_["10_voovoo"]("a,l,i,b,j,k") += reused_["22_vooo"]("a,m,l,i") * this->T1_("b,m") * R2("j,k");
    tmps_["10_voovoo"]("a,l,i,b,j,k") -= R3("b,j,k,n") * reused_["46_vooo"]("a,i,l,n");
    
    // sigmaR4 += -1.00 P(k,l) P(a,b) <m,a||c,d> R3(b,j,k,m) this->T1_(c,l) this->T1_(d,i) 
    //                += -1.00 P(k,l) P(a,b) <m,a||c,d> R2(j,k) this->T1_(b,m) this->T1_(c,l) this->T1_(d,i) 
    //                += +1.00 P(k,l) P(a,b) <n,m||c,d> R3(b,j,k,n) this->T1_(a,m) this->T1_(c,l) this->T1_(d,i) 
    sigmaR4("a,b,i,j,k,l") += tmps_["10_voovoo"]("a,l,i,b,j,k");
    sigmaR4("a,b,i,j,k,l") -= tmps_["10_voovoo"]("a,k,i,b,j,l");
    sigmaR4("a,b,i,j,k,l") -= tmps_["10_voovoo"]("b,l,i,a,j,k");
    sigmaR4("a,b,i,j,k,l") += tmps_["10_voovoo"]("b,k,i,a,j,l");
    
    // sigmaR4 += -1.00 P(a,b) <m,a||c,d> R3(b,k,l,m) this->T1_(c,j) this->T1_(d,i) 
    //                += -1.00 P(a,b) <m,a||c,d> R2(k,l) this->T1_(b,m) this->T1_(c,j) this->T1_(d,i) 
    //                += +1.00 P(a,b) <n,m||c,d> R3(b,k,l,n) this->T1_(a,m) this->T1_(c,j) this->T1_(d,i) 
    sigmaR4("a,b,i,j,k,l") += tmps_["10_voovoo"]("a,j,i,b,k,l");
    sigmaR4("a,b,i,j,k,l") -= tmps_["10_voovoo"]("b,j,i,a,k,l");
    
    // sigmaR4 += -1.00 P(a,b) <m,a||c,d> R3(b,i,l,m) this->T1_(c,k) this->T1_(d,j) 
    //                += -1.00 P(a,b) <m,a||c,d> R2(i,l) this->T1_(b,m) this->T1_(c,k) this->T1_(d,j) 
    //                += +1.00 P(a,b) <n,m||c,d> R3(b,i,l,n) this->T1_(a,m) this->T1_(c,k) this->T1_(d,j) 
    sigmaR4("a,b,i,j,k,l") += tmps_["10_voovoo"]("a,k,j,b,i,l");
    sigmaR4("a,b,i,j,k,l") -= tmps_["10_voovoo"]("b,k,j,a,i,l");
    
    // sigmaR4 += -1.00 P(j,k) P(a,b) <m,a||c,d> R3(b,i,j,m) this->T1_(c,l) this->T1_(d,k) 
    //                += -1.00 P(j,k) P(a,b) <m,a||c,d> R2(i,j) this->T1_(b,m) this->T1_(c,l) this->T1_(d,k) 
    //                += +1.00 P(j,k) P(a,b) <n,m||c,d> R3(b,i,j,n) this->T1_(a,m) this->T1_(c,l) this->T1_(d,k) 
    sigmaR4("a,b,i,j,k,l") += tmps_["10_voovoo"]("a,l,k,b,i,j");
    sigmaR4("a,b,i,j,k,l") -= tmps_["10_voovoo"]("a,l,j,b,i,k");
    sigmaR4("a,b,i,j,k,l") -= tmps_["10_voovoo"]("b,l,k,a,i,j");
    sigmaR4("a,b,i,j,k,l") += tmps_["10_voovoo"]("b,l,j,a,i,k");
    TAmanager.free("voovoo", std::move(tmps_["10_voovoo"]));
    
    // flops: o4v2L1  = o1v1 o4v2L1 o6v2L1 o5v2L1 o4v2L1 o5v2L1 o5v2L1 o4v2L1 o5v3L1 o4v2L1 o5v2L1 o5v2L1 o4v2L1 o5v3L1 o4v2L1 o5v2L1 o5v2L1 o4v2L1 o5v3L1 o4v2L1 o4v3L1 o4v2L1 o6v2L1 o5v2L1 o4v2L1 o4v2L1 o4v2L1
    //  mems: o4v2L1  = o1v1 o4v2L1 o5v1L1 o4v2L1 o4v2L1 o5v1L1 o4v2L1 o4v2L1 o4v2L1 o4v2L1 o5v1L1 o4v2L1 o4v2L1 o4v2L1 o4v2L1 o5v1L1 o4v2L1 o4v2L1 o4v2L1 o4v2L1 o4v2L1 o4v2L1 o5v1L1 o4v2L1 o4v2L1 o4v2L1 o4v2L1
    tmps_.emplace(std::make_pair("11_vovooo", TAmanager.malloc<MatsT>("vovooo")));
    tmps_["11_vovooo"]("a,l,b,i,j,k")  = (reused_["41_vo"]("a,l") + -1.00 * this->fockMatrix_ta["vo"]("a,l")) * R3("b,i,j,k");
    tmps_["11_vovooo"]("a,l,b,i,j,k") += conj(this->antiSymMoints["vooo"]("c,l,m,n")) * R4("c,b,i,j,k,n") * this->T1_("a,m");
    tmps_["11_vovooo"]("a,l,b,i,j,k") += R3("d,i,j,k") * reused_["21_vovo"]("a,m,d,l") * this->T1_("b,m");
    tmps_["11_vovooo"]("a,l,b,i,j,k") -= R4("d,b,i,j,k,n") * reused_["14_voov"]("a,l,n,d");
    tmps_["11_vovooo"]("a,l,b,i,j,k") -= R3("d,i,j,k") * reused_["8_voov"]("b,l,m,d") * this->T1_("a,m");
    tmps_["11_vovooo"]("a,l,b,i,j,k") -= this->antiSymMoints["vovo"]("a,m,c,l") * R4("c,b,i,j,k,m");
    tmps_["11_vovooo"]("a,l,b,i,j,k") -= this->antiSymMoints["vovo"]("a,m,c,l") * R3("c,i,j,k") * this->T1_("b,m");
    tmps_["11_vovooo"]("a,l,b,i,j,k") += R4("d,b,i,j,k,m") * reused_["21_vovo"]("a,m,d,l");
    tmps_["11_vovooo"]("a,l,b,i,j,k") += R3("d,i,j,k") * reused_["111_vvvo"]("a,d,b,l");
    tmps_["11_vovooo"]("a,l,b,i,j,k") -= R4("d,b,i,j,k,n") * reused_["24_ooov"]("l,m,n,d") * this->T1_("a,m");
    tmps_["11_vovooo"]("a,l,b,i,j,k") += R3("b,i,j,k") * reused_["92_vo"]("a,l");
    
    // sigmaR4 += -1.00 P(k,l) P(a,b) <n,m||c,d> R3(b,i,j,k) this->T1_(a,n) this->T1_(c,m) this->T1_(d,l) 
    //                += -0.50 P(k,l) P(a,b) <n,m||c,d> R3(b,i,j,k) this->T2_(d,a,n,m) this->T1_(c,l) 
    //                += -1.00 P(k,l) P(a,b) <m,a||c,l> R3(b,i,j,k) this->T1_(c,m) 
    //                += -1.00 P(k,l) P(a,b) f(a,c) R3(b,i,j,k) this->T1_(c,l) 
    //                += +1.00 P(k,l) P(a,b) f(m,c) R3(b,i,j,k) this->T2_(c,a,l,m) 
    //                += +1.00 P(k,l) P(a,b) f(m,l) R3(b,i,j,k) this->T1_(a,m) 
    //                += -0.25 P(k,l) P(a,b) <n,m||c,d> R3(b,i,j,k) this->T3_(c,d,a,l,n,m) 
    //                += +0.50 P(k,l) P(a,b) <m,a||c,d> R3(b,i,j,k) this->T2_(c,d,l,m) 
    //                += +0.50 P(k,l) P(a,b) <n,m||c,l> R3(b,i,j,k) this->T2_(c,a,n,m) 
    //                += -0.50 P(k,l) P(a,b) <n,m||c,d> R3(b,i,j,k) this->T1_(a,m) this->T2_(c,d,l,n) 
    //                += -1.00 P(k,l) P(a,b) <n,m||c,l> R3(b,i,j,k) this->T1_(a,n) this->T1_(c,m) 
    //                += +1.00 P(k,l) P(a,b) f(m,c) R3(b,i,j,k) this->T1_(a,m) this->T1_(c,l) 
    //                += -1.00 P(k,l) P(a,b) f(a,l) R3(b,i,j,k) 
    //                += -1.00 P(k,l) P(a,b) <n,m||c,l> R4(c,b,i,j,k,n) this->T1_(a,m) 
    //                += -1.00 P(k,l) P(a,b) <m,a||c,d> R3(d,i,j,k) this->T1_(b,m) this->T1_(c,l) 
    //                += +1.00 P(k,l) P(a,b) <n,m||c,d> R4(d,b,i,j,k,n) this->T2_(c,a,l,m) 
    //                += +1.00 P(k,l) P(a,b) <n,m||c,d> R3(d,i,j,k) this->T1_(a,m) this->T2_(c,b,l,n) 
    //                += +1.00 P(k,l) P(a,b) <m,a||c,l> R4(c,b,i,j,k,m) 
    //                += +1.00 P(k,l) P(a,b) <m,a||c,l> R3(c,i,j,k) this->T1_(b,m) 
    //                += -1.00 P(k,l) P(a,b) <m,a||c,d> R4(d,b,i,j,k,m) this->T1_(c,l) 
    //                += -1.00 P(k,l) P(a,b) <m,a||c,d> R3(d,i,j,k) this->T2_(c,b,l,m) 
    //                += +1.00 P(k,l) P(a,b) <n,m||c,d> R4(d,b,i,j,k,n) this->T1_(a,m) this->T1_(c,l) 
    //                += -1.00 P(k,l) P(a,b) <n,m||c,d> R3(b,i,j,k) this->T2_(d,a,l,n) this->T1_(c,m) 
    //                += -1.00 P(k,l) P(a,b) <m,a||c,d> R3(b,i,j,k) this->T1_(c,m) this->T1_(d,l) 
    sigmaR4("a,b,i,j,k,l") += tmps_["11_vovooo"]("a,l,b,i,j,k");
    sigmaR4("a,b,i,j,k,l") -= tmps_["11_vovooo"]("a,k,b,i,j,l");
    sigmaR4("a,b,i,j,k,l") -= tmps_["11_vovooo"]("b,l,a,i,j,k");
    sigmaR4("a,b,i,j,k,l") += tmps_["11_vovooo"]("b,k,a,i,j,l");
    
    // sigmaR4 += -1.00 P(i,j) P(a,b) <n,m||c,d> R3(b,i,k,l) this->T1_(a,n) this->T1_(c,m) this->T1_(d,j) 
    //                += -0.50 P(i,j) P(a,b) <n,m||c,d> R3(b,i,k,l) this->T2_(d,a,n,m) this->T1_(c,j) 
    //                += -1.00 P(i,j) P(a,b) <m,a||c,j> R3(b,i,k,l) this->T1_(c,m) 
    //                += -1.00 P(i,j) P(a,b) f(a,c) R3(b,i,k,l) this->T1_(c,j) 
    //                += +1.00 P(i,j) P(a,b) f(m,c) R3(b,i,k,l) this->T2_(c,a,j,m) 
    //                += +1.00 P(i,j) P(a,b) f(m,j) R3(b,i,k,l) this->T1_(a,m) 
    //                += -0.25 P(i,j) P(a,b) <n,m||c,d> R3(b,i,k,l) this->T3_(c,d,a,j,n,m) 
    //                += +0.50 P(i,j) P(a,b) <m,a||c,d> R3(b,i,k,l) this->T2_(c,d,j,m) 
    //                += +0.50 P(i,j) P(a,b) <n,m||c,j> R3(b,i,k,l) this->T2_(c,a,n,m) 
    //                += -0.50 P(i,j) P(a,b) <n,m||c,d> R3(b,i,k,l) this->T1_(a,m) this->T2_(c,d,j,n) 
    //                += -1.00 P(i,j) P(a,b) <n,m||c,j> R3(b,i,k,l) this->T1_(a,n) this->T1_(c,m) 
    //                += +1.00 P(i,j) P(a,b) f(m,c) R3(b,i,k,l) this->T1_(a,m) this->T1_(c,j) 
    //                += -1.00 P(i,j) P(a,b) f(a,j) R3(b,i,k,l) 
    //                += -1.00 P(i,j) P(a,b) <n,m||c,j> R4(c,b,i,k,l,n) this->T1_(a,m) 
    //                += -1.00 P(i,j) P(a,b) <m,a||c,d> R3(d,i,k,l) this->T1_(b,m) this->T1_(c,j) 
    //                += +1.00 P(i,j) P(a,b) <n,m||c,d> R4(d,b,i,k,l,n) this->T2_(c,a,j,m) 
    //                += +1.00 P(i,j) P(a,b) <n,m||c,d> R3(d,i,k,l) this->T1_(a,m) this->T2_(c,b,j,n) 
    //                += +1.00 P(i,j) P(a,b) <m,a||c,j> R4(c,b,i,k,l,m) 
    //                += +1.00 P(i,j) P(a,b) <m,a||c,j> R3(c,i,k,l) this->T1_(b,m) 
    //                += -1.00 P(i,j) P(a,b) <m,a||c,d> R4(d,b,i,k,l,m) this->T1_(c,j) 
    //                += -1.00 P(i,j) P(a,b) <m,a||c,d> R3(d,i,k,l) this->T2_(c,b,j,m) 
    //                += +1.00 P(i,j) P(a,b) <n,m||c,d> R4(d,b,i,k,l,n) this->T1_(a,m) this->T1_(c,j) 
    //                += -1.00 P(i,j) P(a,b) <n,m||c,d> R3(b,i,k,l) this->T2_(d,a,j,n) this->T1_(c,m) 
    //                += -1.00 P(i,j) P(a,b) <m,a||c,d> R3(b,i,k,l) this->T1_(c,m) this->T1_(d,j) 
    sigmaR4("a,b,i,j,k,l") += tmps_["11_vovooo"]("a,j,b,i,k,l");
    sigmaR4("a,b,i,j,k,l") -= tmps_["11_vovooo"]("a,i,b,j,k,l");
    sigmaR4("a,b,i,j,k,l") -= tmps_["11_vovooo"]("b,j,a,i,k,l");
    sigmaR4("a,b,i,j,k,l") += tmps_["11_vovooo"]("b,i,a,j,k,l");
    TAmanager.free("vovooo", std::move(tmps_["11_vovooo"]));
    
    // flops: o4v2L1  = o5v2L1 o6v1L1 o6v1L1 o5v2L1 o4v2L1 o5v2L1 o4v2L1 o3v1L1 o4v3L1 o4v2L1 o4v1L1 o5v2L1 o4v2L1 o6v1L1 o6v1L1 o5v2L1 o4v2L1 o4v1L1 o5v2L1 o4v2L1 o5v2L1 o4v2L1 o3v1L1 o4v3L1 o4v2L1 o4v3L1 o4v2L1 o4v3L1 o4v2L1 o4v3L1 o4v2L1 o5v2L1 o4v2L1
    //  mems: o4v2L1  = o4v2L1 o6v0L1 o5v1L1 o4v2L1 o4v2L1 o4v2L1 o4v2L1 o1v1L1 o4v2L1 o4v2L1 o4v0L1 o4v2L1 o4v2L1 o6v0L1 o5v1L1 o4v2L1 o4v2L1 o4v0L1 o4v2L1 o4v2L1 o4v2L1 o4v2L1 o1v1L1 o4v2L1 o4v2L1 o4v2L1 o4v2L1 o4v2L1 o4v2L1 o4v2L1 o4v2L1 o4v2L1 o4v2L1
    tmps_.emplace(std::make_pair("12_ovvooo", TAmanager.malloc<MatsT>("ovvooo")));
    tmps_["12_ovvooo"]("l,a,b,i,j,k")  = R4("a,b,i,j,k,m") * this->fockMatrix_ta["oo"]("m,l");
    tmps_["12_ovvooo"]("l,a,b,i,j,k") += R3("d,i,j,k") * reused_["24_ooov"]("l,m,n,d") * this->T1_("a,m") * this->T1_("b,n");
    tmps_["12_ovvooo"]("l,a,b,i,j,k") += R4("a,b,i,j,k,m") * reused_["32_oo"]("l,m");
    tmps_["12_ovvooo"]("l,a,b,i,j,k") += 0.50 * conj(this->antiSymMoints["vooo"]("c,l,m,n")) * R2("n,m") * this->T3_("c,a,b,i,j,k");
    tmps_["12_ovvooo"]("l,a,b,i,j,k") -= R3("d,i,j,k") * reused_["27_ov"]("n,d") * this->T2_("a,b,l,n");
    tmps_["12_ovvooo"]("l,a,b,i,j,k") -= R3("c,i,j,k") * conj(this->antiSymMoints["vooo"]("c,l,m,n")) * this->T1_("a,m") * this->T1_("b,n");
    tmps_["12_ovvooo"]("l,a,b,i,j,k") -= this->fockMatrix_ta["ov"]("m,c") * R3("c,i,j,k") * this->T2_("a,b,l,m");
    tmps_["12_ovvooo"]("l,a,b,i,j,k") += R4("a,b,i,j,k,n") * reused_["44_oo"]("n,l");
    tmps_["12_ovvooo"]("l,a,b,i,j,k") -= 0.50 * R2("n,m") * reused_["24_ooov"]("l,m,n,d") * this->T3_("d,a,b,i,j,k");
    tmps_["12_ovvooo"]("l,a,b,i,j,k") -= 0.50 * reused_["97_vvov"]("a,b,l,d") * R3("d,i,j,k");
    tmps_["12_ovvooo"]("l,a,b,i,j,k") -= this->antiSymMoints["vvvo"]("a,b,c,l") * R3("c,i,j,k");
    tmps_["12_ovvooo"]("l,a,b,i,j,k") += 0.50 * R3("c,i,j,k") * reused_["10_vovv"]("c,l,a,b");
    tmps_["12_ovvooo"]("l,a,b,i,j,k") -= 0.50 * this->T2_("a,b,l,m") * tmps_["5_oooo"]("i,j,k,m");
    TAmanager.free("oooo", std::move(tmps_["5_oooo"]));
    
    // sigmaR4 += -1.00 P(i,j) f(m,j) R4(a,b,i,k,l,m) 
    //                += +1.00 P(i,j) <n,m||c,d> R3(d,i,k,l) this->T1_(a,m) this->T1_(b,n) this->T1_(c,j) 
    //                += -1.00 P(i,j) f(m,c) R4(a,b,i,k,l,m) this->T1_(c,j) 
    //                += +0.50 P(i,j) <n,m||c,j> R2(n,m) this->T3_(c,a,b,i,k,l) 
    //                += -1.00 P(i,j) <n,m||c,d> R3(d,i,k,l) this->T2_(a,b,j,n) this->T1_(c,m) 
    //                += -1.00 P(i,j) <n,m||c,j> R3(c,i,k,l) this->T1_(a,m) this->T1_(b,n) 
    //                += +1.00 P(i,j) f(m,c) R3(c,i,k,l) this->T2_(a,b,j,m) 
    //                += +1.00 P(i,j) <n,m||c,j> R4(a,b,i,k,l,n) this->T1_(c,m) 
    //                += -0.50 P(i,j) <n,m||c,d> R4(a,b,i,k,l,n) this->T2_(c,d,j,m) 
    //                += +1.00 P(i,j) <n,m||c,d> R4(a,b,i,k,l,n) this->T1_(c,m) this->T1_(d,j) 
    //                += -0.50 P(i,j) <n,m||c,d> R2(n,m) this->T3_(d,a,b,i,k,l) this->T1_(c,j) 
    //                += -0.50 P(i,j) <n,m||c,d> R3(d,i,k,l) this->T3_(c,a,b,j,n,m) 
    //                += -0.50 P(i,j) <n,m||c,d> R3(d,i,k,l) this->T2_(a,b,n,m) this->T1_(c,j) 
    //                += -1.00 P(i,j) <a,b||c,d> R3(d,i,k,l) this->T1_(c,j) 
    //                += +1.00 P(i,j) <a,b||c,j> R3(c,i,k,l) 
    //                += +0.50 P(i,j) <n,m||c,j> R3(c,i,k,l) this->T2_(a,b,n,m) 
    //                += -0.50 P(i,j) <n,m||c,d> R4(c,d,i,k,l,n) this->T2_(a,b,j,m) 
    sigmaR4("a,b,i,j,k,l") -= tmps_["12_ovvooo"]("j,a,b,i,k,l");
    sigmaR4("a,b,i,j,k,l") += tmps_["12_ovvooo"]("i,a,b,j,k,l");
    
    // sigmaR4 += -1.00 P(k,l) f(m,l) R4(a,b,i,j,k,m) 
    //                += +1.00 P(k,l) <n,m||c,d> R3(d,i,j,k) this->T1_(a,m) this->T1_(b,n) this->T1_(c,l) 
    //                += -1.00 P(k,l) f(m,c) R4(a,b,i,j,k,m) this->T1_(c,l) 
    //                += +0.50 P(k,l) <n,m||c,l> R2(n,m) this->T3_(c,a,b,i,j,k) 
    //                += -1.00 P(k,l) <n,m||c,d> R3(d,i,j,k) this->T2_(a,b,l,n) this->T1_(c,m) 
    //                += -1.00 P(k,l) <n,m||c,l> R3(c,i,j,k) this->T1_(a,m) this->T1_(b,n) 
    //                += +1.00 P(k,l) f(m,c) R3(c,i,j,k) this->T2_(a,b,l,m) 
    //                += +1.00 P(k,l) <n,m||c,l> R4(a,b,i,j,k,n) this->T1_(c,m) 
    //                += -0.50 P(k,l) <n,m||c,d> R4(a,b,i,j,k,n) this->T2_(c,d,l,m) 
    //                += +1.00 P(k,l) <n,m||c,d> R4(a,b,i,j,k,n) this->T1_(c,m) this->T1_(d,l) 
    //                += -0.50 P(k,l) <n,m||c,d> R2(n,m) this->T3_(d,a,b,i,j,k) this->T1_(c,l) 
    //                += -0.50 P(k,l) <n,m||c,d> R3(d,i,j,k) this->T3_(c,a,b,l,n,m) 
    //                += -0.50 P(k,l) <n,m||c,d> R3(d,i,j,k) this->T2_(a,b,n,m) this->T1_(c,l) 
    //                += -1.00 P(k,l) <a,b||c,d> R3(d,i,j,k) this->T1_(c,l) 
    //                += +1.00 P(k,l) <a,b||c,l> R3(c,i,j,k) 
    //                += +0.50 P(k,l) <n,m||c,l> R3(c,i,j,k) this->T2_(a,b,n,m) 
    //                += -0.50 P(k,l) <n,m||c,d> R4(c,d,i,j,k,n) this->T2_(a,b,l,m) 
    sigmaR4("a,b,i,j,k,l") -= tmps_["12_ovvooo"]("l,a,b,i,j,k");
    sigmaR4("a,b,i,j,k,l") += tmps_["12_ovvooo"]("k,a,b,i,j,l");
    TAmanager.free("ovvooo", std::move(tmps_["12_ovvooo"]));
    
    // flops: o4v2L1  = o2v2 o4v2L1 o6v2L1 o4v2L1 o4v2L1 o5v3L1 o4v2L1 o6v2L1 o4v2L1 o2v2L1 o2v3L1 o4v3L1 o4v2L1
    //  mems: o4v2L1  = o2v2 o4v2L1 o4v2L1 o4v2L1 o3v1L1 o4v2L1 o4v2L1 o4v2L1 o4v2L1 o0v2L1 o2v2L1 o4v2L1 o4v2L1
    tmps_.emplace(std::make_pair("13_vvoooo", TAmanager.malloc<MatsT>("vvoooo")));
    tmps_["13_vvoooo"]("a,b,i,l,j,k")  = (reused_["86_vvoo"]("a,b,i,l") + this->antiSymMoints["vvoo"]("a,b,i,l")) * R2("j,k");
    tmps_["13_vvoooo"]("a,b,i,l,j,k") -= 0.50 * this->antiSymMoints["oooo"]("m,n,i,l") * R4("a,b,j,k,n,m");
    tmps_["13_vvoooo"]("a,b,i,l,j,k") += conj(this->antiSymMoints["vvoo"]("c,d,m,n")) * R3("d,j,k,n") * this->T3_("c,a,b,i,l,m");
    tmps_["13_vvoooo"]("a,b,i,l,j,k") -= 0.25 * R4("a,b,j,k,n,m") * reused_["7_oooo"]("i,l,m,n");
    tmps_["13_vvoooo"]("a,b,i,l,j,k") += 0.50 * conj(this->antiSymMoints["vvoo"]("c,d,m,n")) * R2("n,m") * this->T2_("c,a,i,l") * this->T2_("d,b,j,k");
    
    // sigmaR4 += -1.00 P(j,k) <n,m||c,d> R2(i,j) this->T3_(d,a,b,k,l,n) this->T1_(c,m) 
    //                += -0.50 P(j,k) <n,m||c,d> R2(i,j) this->T2_(c,a,n,m) this->T2_(d,b,k,l) 
    //                += +0.50 P(j,k) <n,m||k,l> R2(i,j) this->T2_(a,b,n,m) 
    //                += -0.50 P(j,k) <n,m||c,d> R2(i,j) this->T2_(c,a,k,l) this->T2_(d,b,n,m) 
    //                += +0.50 P(j,k) <a,b||c,d> R2(i,j) this->T2_(c,d,k,l) 
    //                += -1.00 P(j,k) <n,m||k,l> R2(i,j) this->T1_(a,m) this->T1_(b,n) 
    //                += +1.00 P(j,k) f(m,c) R2(i,j) this->T3_(c,a,b,k,l,m) 
    //                += +0.25 P(j,k) <n,m||c,d> R2(i,j) this->T2_(a,b,n,m) this->T2_(c,d,k,l) 
    //                += -0.50 P(j,k) <n,m||c,d> R2(i,j) this->T1_(a,m) this->T1_(b,n) this->T2_(c,d,k,l) 
    //                += +1.00 P(j,k) <a,b||k,l> R2(i,j) 
    //                += +0.50 P(j,k) <n,m||k,l> R4(a,b,i,j,n,m) 
    //                += -1.00 P(j,k) <n,m||c,d> R3(d,i,j,n) this->T3_(c,a,b,k,l,m) 
    //                += +0.25 P(j,k) <n,m||c,d> R4(a,b,i,j,n,m) this->T2_(c,d,k,l) 
    //                += -0.50 P(j,k) <n,m||c,d> R2(n,m) this->T2_(c,a,k,l) this->T2_(d,b,i,j) 
    sigmaR4("a,b,i,j,k,l") += tmps_["13_vvoooo"]("a,b,k,l,i,j");
    sigmaR4("a,b,i,j,k,l") -= tmps_["13_vvoooo"]("a,b,j,l,i,k");
    
    // sigmaR4 += -1.00 P(k,l) <n,m||c,d> R2(j,k) this->T3_(d,a,b,i,l,n) this->T1_(c,m) 
    //                += -0.50 P(k,l) <n,m||c,d> R2(j,k) this->T2_(c,a,n,m) this->T2_(d,b,i,l) 
    //                += +0.50 P(k,l) <n,m||i,l> R2(j,k) this->T2_(a,b,n,m) 
    //                += -0.50 P(k,l) <n,m||c,d> R2(j,k) this->T2_(c,a,i,l) this->T2_(d,b,n,m) 
    //                += +0.50 P(k,l) <a,b||c,d> R2(j,k) this->T2_(c,d,i,l) 
    //                += -1.00 P(k,l) <n,m||i,l> R2(j,k) this->T1_(a,m) this->T1_(b,n) 
    //                += +1.00 P(k,l) f(m,c) R2(j,k) this->T3_(c,a,b,i,l,m) 
    //                += +0.25 P(k,l) <n,m||c,d> R2(j,k) this->T2_(a,b,n,m) this->T2_(c,d,i,l) 
    //                += -0.50 P(k,l) <n,m||c,d> R2(j,k) this->T1_(a,m) this->T1_(b,n) this->T2_(c,d,i,l) 
    //                += +1.00 P(k,l) <a,b||i,l> R2(j,k) 
    //                += +0.50 P(k,l) <n,m||i,l> R4(a,b,j,k,n,m) 
    //                += -1.00 P(k,l) <n,m||c,d> R3(d,j,k,n) this->T3_(c,a,b,i,l,m) 
    //                += +0.25 P(k,l) <n,m||c,d> R4(a,b,j,k,n,m) this->T2_(c,d,i,l) 
    //                += -0.50 P(k,l) <n,m||c,d> R2(n,m) this->T2_(c,a,i,l) this->T2_(d,b,j,k) 
    sigmaR4("a,b,i,j,k,l") += tmps_["13_vvoooo"]("a,b,i,l,j,k");
    sigmaR4("a,b,i,j,k,l") -= tmps_["13_vvoooo"]("a,b,i,k,j,l");
    
    // sigmaR4 += -1.00 P(i,k) <n,m||c,d> R2(i,l) this->T3_(d,a,b,j,k,n) this->T1_(c,m) 
    //                += -0.50 P(i,k) <n,m||c,d> R2(i,l) this->T2_(c,a,n,m) this->T2_(d,b,j,k) 
    //                += +0.50 P(i,k) <n,m||j,k> R2(i,l) this->T2_(a,b,n,m) 
    //                += -0.50 P(i,k) <n,m||c,d> R2(i,l) this->T2_(c,a,j,k) this->T2_(d,b,n,m) 
    //                += +0.50 P(i,k) <a,b||c,d> R2(i,l) this->T2_(c,d,j,k) 
    //                += -1.00 P(i,k) <n,m||j,k> R2(i,l) this->T1_(a,m) this->T1_(b,n) 
    //                += +1.00 P(i,k) f(m,c) R2(i,l) this->T3_(c,a,b,j,k,m) 
    //                += +0.25 P(i,k) <n,m||c,d> R2(i,l) this->T2_(a,b,n,m) this->T2_(c,d,j,k) 
    //                += -0.50 P(i,k) <n,m||c,d> R2(i,l) this->T1_(a,m) this->T1_(b,n) this->T2_(c,d,j,k) 
    //                += +1.00 P(i,k) <a,b||j,k> R2(i,l) 
    //                += +0.50 P(i,k) <n,m||j,k> R4(a,b,i,l,n,m) 
    //                += -1.00 P(i,k) <n,m||c,d> R3(d,i,l,n) this->T3_(c,a,b,j,k,m) 
    //                += +0.25 P(i,k) <n,m||c,d> R4(a,b,i,l,n,m) this->T2_(c,d,j,k) 
    //                += -0.50 P(i,k) <n,m||c,d> R2(n,m) this->T2_(c,a,j,k) this->T2_(d,b,i,l) 
    sigmaR4("a,b,i,j,k,l") += tmps_["13_vvoooo"]("a,b,j,k,i,l");
    sigmaR4("a,b,i,j,k,l") -= tmps_["13_vvoooo"]("a,b,j,i,k,l");
    TAmanager.free("vvoooo", std::move(tmps_["13_vvoooo"]));
    
    // flops: o4v2L1  = o4v2L1 o5v1L1 o5v2L1 o4v2L1
    //  mems: o4v2L1  = o4v2L1 o4v0L1 o4v2L1 o4v2L1
    tmps_.emplace(std::make_pair("14_vvoooo", TAmanager.malloc<MatsT>("vvoooo")));
    tmps_["14_vvoooo"]("a,b,l,i,j,k")  = R2("j,k") * reused_["105_vvoo"]("a,b,l,i");
    tmps_["14_vvoooo"]("a,b,l,i,j,k") += R3("c,j,k,n") * conj(this->antiSymMoints["vooo"]("c,l,m,n")) * this->T2_("a,b,i,m");
    
    // sigmaR4 += +1.00 P(k,l) <a,b||c,j> R2(i,k) this->T1_(c,l) 
    //                += -0.50 P(k,l) <n,m||c,d> R2(i,k) this->T3_(d,a,b,l,n,m) this->T1_(c,j) 
    //                += -1.00 P(k,l) <n,m||c,j> R2(i,k) this->T1_(a,m) this->T1_(b,n) this->T1_(c,l) 
    //                += -0.50 P(k,l) <n,m||c,d> R2(i,k) this->T2_(a,b,l,n) this->T2_(c,d,j,m) 
    //                += -1.00 P(k,l) <n,m||c,j> R3(c,i,k,n) this->T2_(a,b,l,m) 
    sigmaR4("a,b,i,j,k,l") += tmps_["14_vvoooo"]("a,b,j,l,i,k");
    sigmaR4("a,b,i,j,k,l") -= tmps_["14_vvoooo"]("a,b,j,k,i,l");
    
    // sigmaR4 += +1.00 P(k,l) <a,b||c,l> R2(j,k) this->T1_(c,i) 
    //                += -0.50 P(k,l) <n,m||c,d> R2(j,k) this->T3_(d,a,b,i,n,m) this->T1_(c,l) 
    //                += -1.00 P(k,l) <n,m||c,l> R2(j,k) this->T1_(a,m) this->T1_(b,n) this->T1_(c,i) 
    //                += -0.50 P(k,l) <n,m||c,d> R2(j,k) this->T2_(a,b,i,n) this->T2_(c,d,l,m) 
    //                += -1.00 P(k,l) <n,m||c,l> R3(c,j,k,n) this->T2_(a,b,i,m) 
    sigmaR4("a,b,i,j,k,l") += tmps_["14_vvoooo"]("a,b,l,i,j,k");
    sigmaR4("a,b,i,j,k,l") -= tmps_["14_vvoooo"]("a,b,k,i,j,l");
    
    // sigmaR4 += -1.00 P(j,l) <a,b||c,k> R2(i,j) this->T1_(c,l) 
    //                += +0.50 P(j,l) <n,m||c,d> R2(i,j) this->T3_(d,a,b,l,n,m) this->T1_(c,k) 
    //                += +1.00 P(j,l) <n,m||c,k> R2(i,j) this->T1_(a,m) this->T1_(b,n) this->T1_(c,l) 
    //                += +0.50 P(j,l) <n,m||c,d> R2(i,j) this->T2_(a,b,l,n) this->T2_(c,d,k,m) 
    //                += +1.00 P(j,l) <n,m||c,k> R3(c,i,j,n) this->T2_(a,b,l,m) 
    sigmaR4("a,b,i,j,k,l") -= tmps_["14_vvoooo"]("a,b,k,l,i,j");
    sigmaR4("a,b,i,j,k,l") += tmps_["14_vvoooo"]("a,b,k,j,i,l");
    
    // sigmaR4 += +1.00 P(j,k) <a,b||c,l> R2(i,j) this->T1_(c,k) 
    //                += -0.50 P(j,k) <n,m||c,d> R2(i,j) this->T3_(d,a,b,k,n,m) this->T1_(c,l) 
    //                += -1.00 P(j,k) <n,m||c,l> R2(i,j) this->T1_(a,m) this->T1_(b,n) this->T1_(c,k) 
    //                += -0.50 P(j,k) <n,m||c,d> R2(i,j) this->T2_(a,b,k,n) this->T2_(c,d,l,m) 
    //                += -1.00 P(j,k) <n,m||c,l> R3(c,i,j,n) this->T2_(a,b,k,m) 
    sigmaR4("a,b,i,j,k,l") += tmps_["14_vvoooo"]("a,b,l,k,i,j");
    sigmaR4("a,b,i,j,k,l") -= tmps_["14_vvoooo"]("a,b,l,j,i,k");
    
    // sigmaR4 += +1.00 P(i,j) <a,b||c,j> R2(k,l) this->T1_(c,i) 
    //                += -0.50 P(i,j) <n,m||c,d> R2(k,l) this->T3_(d,a,b,i,n,m) this->T1_(c,j) 
    //                += -1.00 P(i,j) <n,m||c,j> R2(k,l) this->T1_(a,m) this->T1_(b,n) this->T1_(c,i) 
    //                += -0.50 P(i,j) <n,m||c,d> R2(k,l) this->T2_(a,b,i,n) this->T2_(c,d,j,m) 
    //                += -1.00 P(i,j) <n,m||c,j> R3(c,k,l,n) this->T2_(a,b,i,m) 
    sigmaR4("a,b,i,j,k,l") += tmps_["14_vvoooo"]("a,b,j,i,k,l");
    sigmaR4("a,b,i,j,k,l") -= tmps_["14_vvoooo"]("a,b,i,j,k,l");
    
    // sigmaR4 += -1.00 P(k,l) <a,b||c,i> R2(j,k) this->T1_(c,l) 
    //                += +0.50 P(k,l) <n,m||c,d> R2(j,k) this->T3_(d,a,b,l,n,m) this->T1_(c,i) 
    //                += +1.00 P(k,l) <n,m||c,i> R2(j,k) this->T1_(a,m) this->T1_(b,n) this->T1_(c,l) 
    //                += +0.50 P(k,l) <n,m||c,d> R2(j,k) this->T2_(a,b,l,n) this->T2_(c,d,i,m) 
    //                += +1.00 P(k,l) <n,m||c,i> R3(c,j,k,n) this->T2_(a,b,l,m) 
    sigmaR4("a,b,i,j,k,l") -= tmps_["14_vvoooo"]("a,b,i,l,j,k");
    sigmaR4("a,b,i,j,k,l") += tmps_["14_vvoooo"]("a,b,i,k,j,l");
    TAmanager.free("vvoooo", std::move(tmps_["14_vvoooo"]));
    
    // flops: o4v2L1  = o4v2L1 o5v2L1 o4v2L1
    //  mems: o4v2L1  = o4v2L1 o4v2L1 o4v2L1
    tmps_.emplace(std::make_pair("15_ooovov", TAmanager.malloc<MatsT>("ooovov")));
    tmps_["15_ooovov"]("j,k,l,b,i,a")  = R2("j,k") * reused_["106_ovov"]("l,b,i,a");
    tmps_["15_ooovov"]("j,k,l,b,i,a") += R3("b,j,k,n") * reused_["55_oovo"]("n,l,a,i");
    
    // sigmaR4 += -1.00 P(i,j) P(a,b) <n,m||c,j> R2(k,l) this->T1_(a,m) this->T2_(c,b,i,n) 
    //                += -1.00 P(i,j) P(a,b) <m,a||c,d> R2(k,l) this->T2_(d,b,i,m) this->T1_(c,j) 
    //                += -1.00 P(i,j) P(a,b) <n,m||c,j> R3(b,k,l,n) this->T2_(c,a,i,m) 
    //                += -1.00 P(i,j) P(a,b) <n,m||c,j> R3(b,k,l,n) this->T1_(a,m) this->T1_(c,i) 
    sigmaR4("a,b,i,j,k,l") += tmps_["15_ooovov"]("k,l,j,b,i,a");
    sigmaR4("a,b,i,j,k,l") -= tmps_["15_ooovov"]("k,l,i,b,j,a");
    sigmaR4("a,b,i,j,k,l") -= tmps_["15_ooovov"]("k,l,j,a,i,b");
    sigmaR4("a,b,i,j,k,l") += tmps_["15_ooovov"]("k,l,i,a,j,b");
    
    // sigmaR4 += -1.00 P(k,l) P(a,b) <n,m||c,j> R2(i,k) this->T1_(a,m) this->T2_(c,b,l,n) 
    //                += -1.00 P(k,l) P(a,b) <m,a||c,d> R2(i,k) this->T2_(d,b,l,m) this->T1_(c,j) 
    //                += -1.00 P(k,l) P(a,b) <n,m||c,j> R3(b,i,k,n) this->T2_(c,a,l,m) 
    //                += -1.00 P(k,l) P(a,b) <n,m||c,j> R3(b,i,k,n) this->T1_(a,m) this->T1_(c,l) 
    sigmaR4("a,b,i,j,k,l") += tmps_["15_ooovov"]("i,k,j,b,l,a");
    sigmaR4("a,b,i,j,k,l") -= tmps_["15_ooovov"]("i,l,j,b,k,a");
    sigmaR4("a,b,i,j,k,l") -= tmps_["15_ooovov"]("i,k,j,a,l,b");
    sigmaR4("a,b,i,j,k,l") += tmps_["15_ooovov"]("i,l,j,a,k,b");
    
    // sigmaR4 += -1.00 P(j,k) P(a,b) <n,m||c,l> R2(i,j) this->T1_(a,m) this->T2_(c,b,k,n) 
    //                += -1.00 P(j,k) P(a,b) <m,a||c,d> R2(i,j) this->T2_(d,b,k,m) this->T1_(c,l) 
    //                += -1.00 P(j,k) P(a,b) <n,m||c,l> R3(b,i,j,n) this->T2_(c,a,k,m) 
    //                += -1.00 P(j,k) P(a,b) <n,m||c,l> R3(b,i,j,n) this->T1_(a,m) this->T1_(c,k) 
    sigmaR4("a,b,i,j,k,l") += tmps_["15_ooovov"]("i,j,l,b,k,a");
    sigmaR4("a,b,i,j,k,l") -= tmps_["15_ooovov"]("i,k,l,b,j,a");
    sigmaR4("a,b,i,j,k,l") -= tmps_["15_ooovov"]("i,j,l,a,k,b");
    sigmaR4("a,b,i,j,k,l") += tmps_["15_ooovov"]("i,k,l,a,j,b");
    
    // sigmaR4 += -1.00 P(k,l) P(a,b) <n,m||c,l> R2(j,k) this->T1_(a,m) this->T2_(c,b,i,n) 
    //                += -1.00 P(k,l) P(a,b) <m,a||c,d> R2(j,k) this->T2_(d,b,i,m) this->T1_(c,l) 
    //                += -1.00 P(k,l) P(a,b) <n,m||c,l> R3(b,j,k,n) this->T2_(c,a,i,m) 
    //                += -1.00 P(k,l) P(a,b) <n,m||c,l> R3(b,j,k,n) this->T1_(a,m) this->T1_(c,i) 
    sigmaR4("a,b,i,j,k,l") += tmps_["15_ooovov"]("j,k,l,b,i,a");
    sigmaR4("a,b,i,j,k,l") -= tmps_["15_ooovov"]("j,l,k,b,i,a");
    sigmaR4("a,b,i,j,k,l") -= tmps_["15_ooovov"]("j,k,l,a,i,b");
    sigmaR4("a,b,i,j,k,l") += tmps_["15_ooovov"]("j,l,k,a,i,b");
    
    // sigmaR4 += +1.00 P(k,l) P(a,b) <n,m||c,i> R2(j,k) this->T1_(a,m) this->T2_(c,b,l,n) 
    //                += +1.00 P(k,l) P(a,b) <m,a||c,d> R2(j,k) this->T2_(d,b,l,m) this->T1_(c,i) 
    //                += +1.00 P(k,l) P(a,b) <n,m||c,i> R3(b,j,k,n) this->T2_(c,a,l,m) 
    //                += +1.00 P(k,l) P(a,b) <n,m||c,i> R3(b,j,k,n) this->T1_(a,m) this->T1_(c,l) 
    sigmaR4("a,b,i,j,k,l") -= tmps_["15_ooovov"]("j,k,i,b,l,a");
    sigmaR4("a,b,i,j,k,l") += tmps_["15_ooovov"]("j,l,i,b,k,a");
    sigmaR4("a,b,i,j,k,l") += tmps_["15_ooovov"]("j,k,i,a,l,b");
    sigmaR4("a,b,i,j,k,l") -= tmps_["15_ooovov"]("j,l,i,a,k,b");
    
    // sigmaR4 += +1.00 P(j,l) P(a,b) <n,m||c,k> R2(i,j) this->T1_(a,m) this->T2_(c,b,l,n) 
    //                += +1.00 P(j,l) P(a,b) <m,a||c,d> R2(i,j) this->T2_(d,b,l,m) this->T1_(c,k) 
    //                += +1.00 P(j,l) P(a,b) <n,m||c,k> R3(b,i,j,n) this->T2_(c,a,l,m) 
    //                += +1.00 P(j,l) P(a,b) <n,m||c,k> R3(b,i,j,n) this->T1_(a,m) this->T1_(c,l) 
    sigmaR4("a,b,i,j,k,l") -= tmps_["15_ooovov"]("i,j,k,b,l,a");
    sigmaR4("a,b,i,j,k,l") += tmps_["15_ooovov"]("i,l,k,b,j,a");
    sigmaR4("a,b,i,j,k,l") += tmps_["15_ooovov"]("i,j,k,a,l,b");
    sigmaR4("a,b,i,j,k,l") -= tmps_["15_ooovov"]("i,l,k,a,j,b");
    TAmanager.free("ooovov", std::move(tmps_["15_ooovov"]));
    
    // flops: o4v2L1  = o5v2L1
    //  mems: o4v2L1  = o4v2L1
    tmps_.emplace(std::make_pair("16_ovvooo", TAmanager.malloc<MatsT>("ovvooo")));
    tmps_["16_ovvooo"]("j,a,b,i,k,l")  = R2("l,n") * reused_["109_oovvoo"]("j,n,a,b,i,k");
    
    // sigmaR4 += -1.00 P(i,j) <n,m||c,d> R2(l,n) this->T3_(d,a,b,i,k,m) this->T1_(c,j) 
    sigmaR4("a,b,i,j,k,l") += tmps_["16_ovvooo"]("j,a,b,i,k,l");
    sigmaR4("a,b,i,j,k,l") -= tmps_["16_ovvooo"]("i,a,b,j,k,l");
    
    // sigmaR4 += -1.00 P(i,k) <n,m||c,d> R2(i,n) this->T3_(d,a,b,k,l,m) this->T1_(c,j) 
    sigmaR4("a,b,i,j,k,l") += tmps_["16_ovvooo"]("j,a,b,k,l,i");
    sigmaR4("a,b,i,j,k,l") -= tmps_["16_ovvooo"]("j,a,b,i,l,k");
    
    // sigmaR4 += +1.00 P(i,j) <n,m||c,d> R2(i,n) this->T3_(d,a,b,j,l,m) this->T1_(c,k) 
    sigmaR4("a,b,i,j,k,l") -= tmps_["16_ovvooo"]("k,a,b,j,l,i");
    sigmaR4("a,b,i,j,k,l") += tmps_["16_ovvooo"]("k,a,b,i,l,j");
    
    // sigmaR4 += +1.00 P(j,k) <n,m||c,d> R2(j,n) this->T3_(d,a,b,k,l,m) this->T1_(c,i) 
    sigmaR4("a,b,i,j,k,l") -= tmps_["16_ovvooo"]("i,a,b,k,l,j");
    sigmaR4("a,b,i,j,k,l") += tmps_["16_ovvooo"]("i,a,b,j,l,k");
    
    // sigmaR4 += -1.00 P(i,j) <n,m||c,d> R2(i,n) this->T3_(d,a,b,j,k,m) this->T1_(c,l) 
    sigmaR4("a,b,i,j,k,l") += tmps_["16_ovvooo"]("l,a,b,j,k,i");
    sigmaR4("a,b,i,j,k,l") -= tmps_["16_ovvooo"]("l,a,b,i,k,j");
    
    // sigmaR4 += -1.00 P(k,l) <n,m||c,d> R2(k,n) this->T3_(d,a,b,i,j,m) this->T1_(c,l) 
    sigmaR4("a,b,i,j,k,l") += tmps_["16_ovvooo"]("l,a,b,i,j,k");
    sigmaR4("a,b,i,j,k,l") -= tmps_["16_ovvooo"]("k,a,b,i,j,l");
    TAmanager.free("ovvooo", std::move(tmps_["16_ovvooo"]));
    
    // flops: o4v2L1  = o5v2L1
    //  mems: o4v2L1  = o4v2L1
    tmps_.emplace(std::make_pair("17_ovvooo", TAmanager.malloc<MatsT>("ovvooo")));
    tmps_["17_ovvooo"]("j,a,b,i,k,l")  = R2("l,n") * reused_["103_oovvoo"]("n,j,a,b,i,k");
    
    // sigmaR4 += +1.00 P(i,j) <n,m||c,j> R2(l,n) this->T3_(c,a,b,i,k,m) 
    //                += +1.00 P(i,j) <n,m||c,d> R2(l,n) this->T2_(c,a,j,m) this->T2_(d,b,i,k) 
    sigmaR4("a,b,i,j,k,l") -= tmps_["17_ovvooo"]("j,a,b,i,k,l");
    sigmaR4("a,b,i,j,k,l") += tmps_["17_ovvooo"]("i,a,b,j,k,l");
    
    // sigmaR4 += -1.00 P(i,j) <n,m||c,k> R2(i,n) this->T3_(c,a,b,j,l,m) 
    //                += -1.00 P(i,j) <n,m||c,d> R2(i,n) this->T2_(c,a,k,m) this->T2_(d,b,j,l) 
    sigmaR4("a,b,i,j,k,l") += tmps_["17_ovvooo"]("k,a,b,j,l,i");
    sigmaR4("a,b,i,j,k,l") -= tmps_["17_ovvooo"]("k,a,b,i,l,j");
    
    // sigmaR4 += -1.00 P(j,k) <n,m||c,i> R2(j,n) this->T3_(c,a,b,k,l,m) 
    //                += -1.00 P(j,k) <n,m||c,d> R2(j,n) this->T2_(c,a,i,m) this->T2_(d,b,k,l) 
    sigmaR4("a,b,i,j,k,l") += tmps_["17_ovvooo"]("i,a,b,k,l,j");
    sigmaR4("a,b,i,j,k,l") -= tmps_["17_ovvooo"]("i,a,b,j,l,k");
    
    // sigmaR4 += +1.00 P(k,l) <n,m||c,l> R2(k,n) this->T3_(c,a,b,i,j,m) 
    //                += +1.00 P(k,l) <n,m||c,d> R2(k,n) this->T2_(c,a,l,m) this->T2_(d,b,i,j) 
    sigmaR4("a,b,i,j,k,l") -= tmps_["17_ovvooo"]("l,a,b,i,j,k");
    sigmaR4("a,b,i,j,k,l") += tmps_["17_ovvooo"]("k,a,b,i,j,l");
    
    // sigmaR4 += +1.00 P(i,j) <n,m||c,l> R2(i,n) this->T3_(c,a,b,j,k,m) 
    //                += +1.00 P(i,j) <n,m||c,d> R2(i,n) this->T2_(c,a,l,m) this->T2_(d,b,j,k) 
    sigmaR4("a,b,i,j,k,l") -= tmps_["17_ovvooo"]("l,a,b,j,k,i");
    sigmaR4("a,b,i,j,k,l") += tmps_["17_ovvooo"]("l,a,b,i,k,j");
    
    // sigmaR4 += +1.00 P(i,k) <n,m||c,j> R2(i,n) this->T3_(c,a,b,k,l,m) 
    //                += +1.00 P(i,k) <n,m||c,d> R2(i,n) this->T2_(c,a,j,m) this->T2_(d,b,k,l) 
    sigmaR4("a,b,i,j,k,l") -= tmps_["17_ovvooo"]("j,a,b,k,l,i");
    sigmaR4("a,b,i,j,k,l") += tmps_["17_ovvooo"]("j,a,b,i,l,k");
    TAmanager.free("ovvooo", std::move(tmps_["17_ovvooo"]));
    
    // flops: o6v0L1  = o6v2L1
    //  mems: o6v0L1  = o6v0L1
    tmps_.emplace(std::make_pair("18_oooooo", TAmanager.malloc<MatsT>("oooooo")));
    tmps_["18_oooooo"]("i,j,k,l,m,n")  = conj(this->antiSymMoints["vvoo"]("c,d,m,n")) * R4("c,d,i,j,k,l");
    
    // sigmaR4 += +0.25 <n,m||c,d> R4(c,d,i,j,k,l) this->T2_(a,b,n,m) 
    // flops: o4v2L1 += o6v2L1
    //  mems: o4v2L1 += o4v2L1
    sigmaR4("a,b,i,j,k,l") -= 0.25 * this->T2_("a,b,n,m") * tmps_["18_oooooo"]("i,j,k,l,m,n");
    
    // sigmaR4 += -0.50 <n,m||c,d> R4(c,d,i,j,k,l) this->T1_(a,m) this->T1_(b,n) 
    // flops: o4v2L1 += o6v1L1 o5v2L1
    //  mems: o4v2L1 += o5v1L1 o4v2L1
    sigmaR4("a,b,i,j,k,l") += 0.50 * this->T1_("a,m") * tmps_["18_oooooo"]("i,j,k,l,m,n") * this->T1_("b,n");
    TAmanager.free("oooooo", std::move(tmps_["18_oooooo"]));
    
    // flops: o4v2L1  = o5v0L1 o5v2L1 o5v0L1 o5v2L1 o4v2L1
    //  mems: o4v2L1  = o4v0L1 o4v2L1 o4v0L1 o4v2L1 o4v2L1
    tmps_.emplace(std::make_pair("19_ooovvo", TAmanager.malloc<MatsT>("ooovvo")));
    tmps_["19_ooovvo"]("i,j,k,a,b,l")  = R2("k,n") * reused_["7_oooo"]("i,j,m,n") * this->T2_("a,b,l,m");
    tmps_["19_ooovvo"]("i,j,k,a,b,l") += 2.00 * this->antiSymMoints["oooo"]("m,n,i,j") * R2("k,n") * this->T2_("a,b,l,m");
    
    // sigmaR4 += -0.50 P(k,l) <n,m||c,d> R2(k,n) this->T2_(a,b,l,m) this->T2_(c,d,i,j) 
    //                += -1.00 P(k,l) <n,m||i,j> R2(k,n) this->T2_(a,b,l,m) 
    sigmaR4("a,b,i,j,k,l") += 0.50 * tmps_["19_ooovvo"]("i,j,k,a,b,l");
    sigmaR4("a,b,i,j,k,l") -= 0.50 * tmps_["19_ooovvo"]("i,j,l,a,b,k");
    
    // sigmaR4 += +0.50 P(j,l) <n,m||c,d> R2(j,n) this->T2_(a,b,l,m) this->T2_(c,d,i,k) 
    //                += +1.00 P(j,l) <n,m||i,k> R2(j,n) this->T2_(a,b,l,m) 
    sigmaR4("a,b,i,j,k,l") -= 0.50 * tmps_["19_ooovvo"]("i,k,j,a,b,l");
    sigmaR4("a,b,i,j,k,l") += 0.50 * tmps_["19_ooovvo"]("i,k,l,a,b,j");
    
    // sigmaR4 += -0.50 P(j,k) <n,m||c,d> R2(j,n) this->T2_(a,b,k,m) this->T2_(c,d,i,l) 
    //                += -1.00 P(j,k) <n,m||i,l> R2(j,n) this->T2_(a,b,k,m) 
    sigmaR4("a,b,i,j,k,l") += 0.50 * tmps_["19_ooovvo"]("i,l,j,a,b,k");
    sigmaR4("a,b,i,j,k,l") -= 0.50 * tmps_["19_ooovvo"]("i,l,k,a,b,j");
    
    // sigmaR4 += -0.50 P(i,l) <n,m||c,d> R2(i,n) this->T2_(a,b,l,m) this->T2_(c,d,j,k) 
    //                += -1.00 P(i,l) <n,m||j,k> R2(i,n) this->T2_(a,b,l,m) 
    sigmaR4("a,b,i,j,k,l") += 0.50 * tmps_["19_ooovvo"]("j,k,i,a,b,l");
    sigmaR4("a,b,i,j,k,l") -= 0.50 * tmps_["19_ooovvo"]("j,k,l,a,b,i");
    
    // sigmaR4 += -0.50 P(i,j) <n,m||c,d> R2(i,n) this->T2_(a,b,j,m) this->T2_(c,d,k,l) 
    //                += -1.00 P(i,j) <n,m||k,l> R2(i,n) this->T2_(a,b,j,m) 
    sigmaR4("a,b,i,j,k,l") += 0.50 * tmps_["19_ooovvo"]("k,l,i,a,b,j");
    sigmaR4("a,b,i,j,k,l") -= 0.50 * tmps_["19_ooovvo"]("k,l,j,a,b,i");
    
    // sigmaR4 += +0.50 P(i,k) <n,m||c,d> R2(i,n) this->T2_(a,b,k,m) this->T2_(c,d,j,l) 
    //                += +1.00 P(i,k) <n,m||j,l> R2(i,n) this->T2_(a,b,k,m) 
    sigmaR4("a,b,i,j,k,l") -= 0.50 * tmps_["19_ooovvo"]("j,l,i,a,b,k");
    sigmaR4("a,b,i,j,k,l") += 0.50 * tmps_["19_ooovvo"]("j,l,k,a,b,i");
    TAmanager.free("ooovvo", std::move(tmps_["19_ooovvo"]));
    
    // flops: o4v2L1  = o5v2L1 o5v2L1 o4v2L1 o3v2L1 o4v3L1 o4v2L1
    //  mems: o4v2L1  = o4v2L1 o4v2L1 o4v2L1 o1v1L1 o4v2L1 o4v2L1
    tmps_.emplace(std::make_pair("20_vvoooo", TAmanager.malloc<MatsT>("vvoooo")));
    tmps_["20_vvoooo"]("a,b,i,j,l,k")  = R2("k,n") * reused_["110_ovvooo"]("n,a,b,i,j,l");
    tmps_["20_vvoooo"]("a,b,i,j,l,k") += R2("k,m") * reused_["114_vvoooo"]("a,b,i,j,l,m");
    tmps_["20_vvoooo"]("a,b,i,j,l,k") -= 0.50 * conj(this->antiSymMoints["vvoo"]("c,d,m,n")) * R3("d,k,n,m") * this->T3_("c,a,b,i,j,l");
    
    // sigmaR4 += +1.00 P(k,l) <n,m||c,d> R2(k,n) this->T3_(d,a,b,i,j,l) this->T1_(c,m) 
    //                += -1.00 P(k,l) f(m,c) R2(k,m) this->T3_(c,a,b,i,j,l) 
    //                += -0.50 P(k,l) <n,m||c,d> R3(d,k,n,m) this->T3_(c,a,b,i,j,l) 
    sigmaR4("a,b,i,j,k,l") -= tmps_["20_vvoooo"]("a,b,i,j,l,k");
    sigmaR4("a,b,i,j,k,l") += tmps_["20_vvoooo"]("a,b,i,j,k,l");
    
    // sigmaR4 += +1.00 P(i,j) <n,m||c,d> R2(i,n) this->T3_(d,a,b,j,k,l) this->T1_(c,m) 
    //                += -1.00 P(i,j) f(m,c) R2(i,m) this->T3_(c,a,b,j,k,l) 
    //                += -0.50 P(i,j) <n,m||c,d> R3(d,i,n,m) this->T3_(c,a,b,j,k,l) 
    sigmaR4("a,b,i,j,k,l") -= tmps_["20_vvoooo"]("a,b,j,k,l,i");
    sigmaR4("a,b,i,j,k,l") += tmps_["20_vvoooo"]("a,b,i,k,l,j");
    TAmanager.free("vvoooo", std::move(tmps_["20_vvoooo"]));
    
    // flops: o4v2L1  = o5v0L1 o5v2L1
    //  mems: o4v2L1  = o4v0L1 o4v2L1
    tmps_.emplace(std::make_pair("21_ooovvo", TAmanager.malloc<MatsT>("ooovvo")));
    tmps_["21_ooovvo"]("i,j,k,a,b,l")  = R2("i,n") * reused_["53_oooo"]("j,m,n,k") * this->T2_("a,b,l,m");
    
    // sigmaR4 += -1.00 P(i,k) <n,m||c,j> R2(i,n) this->T2_(a,b,k,m) this->T1_(c,l) 
    sigmaR4("a,b,i,j,k,l") += tmps_["21_ooovvo"]("i,l,j,a,b,k");
    sigmaR4("a,b,i,j,k,l") -= tmps_["21_ooovvo"]("k,l,j,a,b,i");
    
    // sigmaR4 += +1.00 P(j,l) <n,m||c,k> R2(j,n) this->T2_(a,b,l,m) this->T1_(c,i) 
    sigmaR4("a,b,i,j,k,l") -= tmps_["21_ooovvo"]("j,i,k,a,b,l");
    sigmaR4("a,b,i,j,k,l") += tmps_["21_ooovvo"]("l,i,k,a,b,j");
    
    // sigmaR4 += +1.00 P(i,l) <n,m||c,j> R2(i,n) this->T2_(a,b,l,m) this->T1_(c,k) 
    sigmaR4("a,b,i,j,k,l") -= tmps_["21_ooovvo"]("i,k,j,a,b,l");
    sigmaR4("a,b,i,j,k,l") += tmps_["21_ooovvo"]("l,k,j,a,b,i");
    
    // sigmaR4 += +1.00 P(k,l) <n,m||c,i> R2(k,n) this->T2_(a,b,l,m) this->T1_(c,j) 
    sigmaR4("a,b,i,j,k,l") -= tmps_["21_ooovvo"]("k,j,i,a,b,l");
    sigmaR4("a,b,i,j,k,l") += tmps_["21_ooovvo"]("l,j,i,a,b,k");
    
    // sigmaR4 += +1.00 P(i,j) <n,m||c,k> R2(i,n) this->T2_(a,b,j,m) this->T1_(c,l) 
    sigmaR4("a,b,i,j,k,l") -= tmps_["21_ooovvo"]("i,l,k,a,b,j");
    sigmaR4("a,b,i,j,k,l") += tmps_["21_ooovvo"]("j,l,k,a,b,i");
    
    // sigmaR4 += -1.00 P(i,l) <n,m||c,k> R2(i,n) this->T2_(a,b,l,m) this->T1_(c,j) 
    sigmaR4("a,b,i,j,k,l") += tmps_["21_ooovvo"]("i,j,k,a,b,l");
    sigmaR4("a,b,i,j,k,l") -= tmps_["21_ooovvo"]("l,j,k,a,b,i");
    
    // sigmaR4 += -1.00 P(i,j) <n,m||c,l> R2(i,n) this->T2_(a,b,j,m) this->T1_(c,k) 
    sigmaR4("a,b,i,j,k,l") += tmps_["21_ooovvo"]("i,k,l,a,b,j");
    sigmaR4("a,b,i,j,k,l") -= tmps_["21_ooovvo"]("j,k,l,a,b,i");
    
    // sigmaR4 += +1.00 P(i,k) <n,m||c,l> R2(i,n) this->T2_(a,b,k,m) this->T1_(c,j) 
    sigmaR4("a,b,i,j,k,l") -= tmps_["21_ooovvo"]("i,j,l,a,b,k");
    sigmaR4("a,b,i,j,k,l") += tmps_["21_ooovvo"]("k,j,l,a,b,i");
    
    // sigmaR4 += -1.00 P(k,l) <n,m||c,j> R2(k,n) this->T2_(a,b,l,m) this->T1_(c,i) 
    sigmaR4("a,b,i,j,k,l") += tmps_["21_ooovvo"]("k,i,j,a,b,l");
    sigmaR4("a,b,i,j,k,l") -= tmps_["21_ooovvo"]("l,i,j,a,b,k");
    
    // sigmaR4 += -1.00 P(j,l) <n,m||c,i> R2(j,n) this->T2_(a,b,l,m) this->T1_(c,k) 
    sigmaR4("a,b,i,j,k,l") += tmps_["21_ooovvo"]("j,k,i,a,b,l");
    sigmaR4("a,b,i,j,k,l") -= tmps_["21_ooovvo"]("l,k,i,a,b,j");
    
    // sigmaR4 += -1.00 P(j,k) <n,m||c,l> R2(j,n) this->T2_(a,b,k,m) this->T1_(c,i) 
    sigmaR4("a,b,i,j,k,l") += tmps_["21_ooovvo"]("j,i,l,a,b,k");
    sigmaR4("a,b,i,j,k,l") -= tmps_["21_ooovvo"]("k,i,l,a,b,j");
    
    // sigmaR4 += +1.00 P(j,k) <n,m||c,i> R2(j,n) this->T2_(a,b,k,m) this->T1_(c,l) 
    sigmaR4("a,b,i,j,k,l") -= tmps_["21_ooovvo"]("j,l,i,a,b,k");
    sigmaR4("a,b,i,j,k,l") += tmps_["21_ooovvo"]("k,l,i,a,b,j");
    TAmanager.free("ooovvo", std::move(tmps_["21_ooovvo"]));
    
    // flops: o4v2L1  = o6v2L1 o5v2L1 o4v2L1 o5v2L1 o4v2L1
    //  mems: o4v2L1  = o4v2L1 o4v2L1 o4v2L1 o4v2L1 o4v2L1
    tmps_.emplace(std::make_pair("22_vooovo", TAmanager.malloc<MatsT>("vooovo")));
    tmps_["22_vooovo"]("a,i,j,l,b,k")  = R3("b,k,n,m") * reused_["108_vooooo"]("a,i,j,l,m,n");
    tmps_["22_vooovo"]("a,i,j,l,b,k") += 2.00 * R2("k,m") * reused_["115_vovooo"]("a,m,b,i,j,l");
    tmps_["22_vooovo"]("a,i,j,l,b,k") -= 2.00 * reused_["119_vvoooo"]("a,b,i,j,l,n") * R2("k,n");
    
    // sigmaR4 += -0.25 P(k,l) P(a,b) <n,m||c,d> R3(b,k,n,m) this->T3_(c,d,a,i,j,l) 
    //                += -0.50 P(k,l) P(a,b) <m,a||c,d> R2(k,m) this->T3_(c,d,b,i,j,l) 
    //                += +0.50 P(k,l) P(a,b) <n,m||c,d> R2(k,n) this->T1_(a,m) this->T3_(c,d,b,i,j,l) 
    sigmaR4("a,b,i,j,k,l") += 0.25 * tmps_["22_vooovo"]("a,i,j,l,b,k");
    sigmaR4("a,b,i,j,k,l") -= 0.25 * tmps_["22_vooovo"]("a,i,j,k,b,l");
    sigmaR4("a,b,i,j,k,l") -= 0.25 * tmps_["22_vooovo"]("b,i,j,l,a,k");
    sigmaR4("a,b,i,j,k,l") += 0.25 * tmps_["22_vooovo"]("b,i,j,k,a,l");
    
    // sigmaR4 += -0.25 P(i,j) P(a,b) <n,m||c,d> R3(b,i,n,m) this->T3_(c,d,a,j,k,l) 
    //                += -0.50 P(i,j) P(a,b) <m,a||c,d> R2(i,m) this->T3_(c,d,b,j,k,l) 
    //                += +0.50 P(i,j) P(a,b) <n,m||c,d> R2(i,n) this->T1_(a,m) this->T3_(c,d,b,j,k,l) 
    sigmaR4("a,b,i,j,k,l") += 0.25 * tmps_["22_vooovo"]("a,j,k,l,b,i");
    sigmaR4("a,b,i,j,k,l") -= 0.25 * tmps_["22_vooovo"]("a,i,k,l,b,j");
    sigmaR4("a,b,i,j,k,l") -= 0.25 * tmps_["22_vooovo"]("b,j,k,l,a,i");
    sigmaR4("a,b,i,j,k,l") += 0.25 * tmps_["22_vooovo"]("b,i,k,l,a,j");
    TAmanager.free("vooovo", std::move(tmps_["22_vooovo"]));
    
    // flops: o4v2L1  = o2v2 o4v2L1 o5v2L1 o4v2L1 o5v2L1 o4v2L1 o5v2L1 o4v2L1 o4v3L1 o4v3L1 o4v2L1 o3v3L1 o4v3L1 o4v2L1 o4v2L1 o5v2L1 o5v2L1 o4v2L1 o5v2L1 o4v2L1
    //  mems: o4v2L1  = o2v2 o4v2L1 o4v2L1 o4v2L1 o4v2L1 o4v2L1 o4v2L1 o4v2L1 o2v2L1 o4v2L1 o4v2L1 o2v2L1 o4v2L1 o4v2L1 o3v1L1 o5v1L1 o4v2L1 o4v2L1 o4v2L1 o4v2L1
    tmps_.emplace(std::make_pair("23_voovoo", TAmanager.malloc<MatsT>("voovoo")));
    tmps_["23_voovoo"]("a,i,l,b,j,k")  = (reused_["18_voov"]("a,i,l,b") + 2.00 * reused_["65_vvoo"]("a,b,i,l")) * R2("j,k");
    tmps_["23_voovoo"]("a,i,l,b,j,k") += R3("b,j,k,m") * reused_["17_vooo"]("a,m,i,l");
    tmps_["23_voovoo"]("a,i,l,b,j,k") += 2.00 * this->antiSymMoints["vooo"]("a,m,i,l") * R3("b,j,k,m");
    tmps_["23_voovoo"]("a,i,l,b,j,k") -= 2.00 * R3("b,j,k,m") * reused_["60_vooo"]("a,i,l,m");
    tmps_["23_voovoo"]("a,i,l,b,j,k") -= conj(this->antiSymMoints["vvoo"]("c,d,m,n")) * R4("d,b,j,k,n,m") * this->T2_("c,a,i,l");
    tmps_["23_voovoo"]("a,i,l,b,j,k") -= 2.00 * conj(this->antiSymMoints["vvvo"]("c,d,a,m")) * R3("d,j,k,m") * this->T2_("c,b,i,l");
    tmps_["23_voovoo"]("a,i,l,b,j,k") += 2.00 * conj(this->antiSymMoints["vvoo"]("c,d,m,n")) * R3("d,j,k,n") * this->T2_("c,b,i,l") * this->T1_("a,m");
    tmps_["23_voovoo"]("a,i,l,b,j,k") += R3("b,j,k,n") * reused_["88_vooo"]("a,i,l,n");
    
    // sigmaR4 += +0.50 P(i,k) P(a,b) <m,a||c,d> R2(i,l) this->T1_(b,m) this->T2_(c,d,j,k) 
    //                += +1.00 P(i,k) P(a,b) <n,m||c,d> R2(i,l) this->T1_(a,n) this->T2_(d,b,j,k) this->T1_(c,m) 
    //                += +0.50 P(i,k) P(a,b) <m,a||c,d> R2(i,l) this->T3_(c,d,b,j,k,m) 
    //                += +1.00 P(i,k) P(a,b) <m,a||c,d> R2(i,l) this->T2_(d,b,j,k) this->T1_(c,m) 
    //                += -0.50 P(i,k) P(a,b) <n,m||c,d> R2(i,l) this->T1_(a,m) this->T3_(c,d,b,j,k,n) 
    //                += +1.00 P(i,k) P(a,b) <m,a||j,k> R2(i,l) this->T1_(b,m) 
    //                += +1.00 P(i,k) P(a,b) f(a,c) R2(i,l) this->T2_(c,b,j,k) 
    //                += -1.00 P(i,k) P(a,b) f(m,c) R2(i,l) this->T1_(a,m) this->T2_(c,b,j,k) 
    //                += +0.50 P(i,k) P(a,b) <m,a||c,d> R3(b,i,l,m) this->T2_(c,d,j,k) 
    //                += +1.00 P(i,k) P(a,b) <m,a||j,k> R3(b,i,l,m) 
    //                += +1.00 P(i,k) P(a,b) f(m,c) R3(b,i,l,m) this->T2_(c,a,j,k) 
    //                += -0.50 P(i,k) P(a,b) <n,m||c,d> R4(d,b,i,l,n,m) this->T2_(c,a,j,k) 
    //                += -1.00 P(i,k) P(a,b) <m,a||c,d> R3(d,i,l,m) this->T2_(c,b,j,k) 
    //                += +1.00 P(i,k) P(a,b) <n,m||c,d> R3(d,i,l,n) this->T1_(a,m) this->T2_(c,b,j,k) 
    //                += -0.50 P(i,k) P(a,b) <n,m||c,d> R3(b,i,l,n) this->T1_(a,m) this->T2_(c,d,j,k) 
    //                += -1.00 P(i,k) P(a,b) <n,m||c,d> R3(b,i,l,n) this->T2_(d,a,j,k) this->T1_(c,m) 
    //                += -0.50 P(i,k) P(a,b) <n,m||c,d> R3(b,i,l,n) this->T3_(c,d,a,j,k,m) 
    //                += -1.00 P(i,k) P(a,b) <n,m||j,k> R3(b,i,l,n) this->T1_(a,m) 
    sigmaR4("a,b,i,j,k,l") -= 0.50 * tmps_["23_voovoo"]("a,j,k,b,i,l");
    sigmaR4("a,b,i,j,k,l") += 0.50 * tmps_["23_voovoo"]("a,j,i,b,k,l");
    sigmaR4("a,b,i,j,k,l") += 0.50 * tmps_["23_voovoo"]("b,j,k,a,i,l");
    sigmaR4("a,b,i,j,k,l") -= 0.50 * tmps_["23_voovoo"]("b,j,i,a,k,l");
    
    // sigmaR4 += +0.50 P(k,l) P(a,b) <m,a||c,d> R2(j,k) this->T1_(b,m) this->T2_(c,d,i,l) 
    //                += +1.00 P(k,l) P(a,b) <n,m||c,d> R2(j,k) this->T1_(a,n) this->T2_(d,b,i,l) this->T1_(c,m) 
    //                += +0.50 P(k,l) P(a,b) <m,a||c,d> R2(j,k) this->T3_(c,d,b,i,l,m) 
    //                += +1.00 P(k,l) P(a,b) <m,a||c,d> R2(j,k) this->T2_(d,b,i,l) this->T1_(c,m) 
    //                += -0.50 P(k,l) P(a,b) <n,m||c,d> R2(j,k) this->T1_(a,m) this->T3_(c,d,b,i,l,n) 
    //                += +1.00 P(k,l) P(a,b) <m,a||i,l> R2(j,k) this->T1_(b,m) 
    //                += +1.00 P(k,l) P(a,b) f(a,c) R2(j,k) this->T2_(c,b,i,l) 
    //                += -1.00 P(k,l) P(a,b) f(m,c) R2(j,k) this->T1_(a,m) this->T2_(c,b,i,l) 
    //                += +0.50 P(k,l) P(a,b) <m,a||c,d> R3(b,j,k,m) this->T2_(c,d,i,l) 
    //                += +1.00 P(k,l) P(a,b) <m,a||i,l> R3(b,j,k,m) 
    //                += +1.00 P(k,l) P(a,b) f(m,c) R3(b,j,k,m) this->T2_(c,a,i,l) 
    //                += -0.50 P(k,l) P(a,b) <n,m||c,d> R4(d,b,j,k,n,m) this->T2_(c,a,i,l) 
    //                += -1.00 P(k,l) P(a,b) <m,a||c,d> R3(d,j,k,m) this->T2_(c,b,i,l) 
    //                += +1.00 P(k,l) P(a,b) <n,m||c,d> R3(d,j,k,n) this->T1_(a,m) this->T2_(c,b,i,l) 
    //                += -0.50 P(k,l) P(a,b) <n,m||c,d> R3(b,j,k,n) this->T1_(a,m) this->T2_(c,d,i,l) 
    //                += -1.00 P(k,l) P(a,b) <n,m||c,d> R3(b,j,k,n) this->T2_(d,a,i,l) this->T1_(c,m) 
    //                += -0.50 P(k,l) P(a,b) <n,m||c,d> R3(b,j,k,n) this->T3_(c,d,a,i,l,m) 
    //                += -1.00 P(k,l) P(a,b) <n,m||i,l> R3(b,j,k,n) this->T1_(a,m) 
    sigmaR4("a,b,i,j,k,l") -= 0.50 * tmps_["23_voovoo"]("a,i,l,b,j,k");
    sigmaR4("a,b,i,j,k,l") += 0.50 * tmps_["23_voovoo"]("a,i,k,b,j,l");
    sigmaR4("a,b,i,j,k,l") += 0.50 * tmps_["23_voovoo"]("b,i,l,a,j,k");
    sigmaR4("a,b,i,j,k,l") -= 0.50 * tmps_["23_voovoo"]("b,i,k,a,j,l");
    
    // sigmaR4 += +0.50 P(j,k) P(a,b) <m,a||c,d> R2(i,j) this->T1_(b,m) this->T2_(c,d,k,l) 
    //                += +1.00 P(j,k) P(a,b) <n,m||c,d> R2(i,j) this->T1_(a,n) this->T2_(d,b,k,l) this->T1_(c,m) 
    //                += +0.50 P(j,k) P(a,b) <m,a||c,d> R2(i,j) this->T3_(c,d,b,k,l,m) 
    //                += +1.00 P(j,k) P(a,b) <m,a||c,d> R2(i,j) this->T2_(d,b,k,l) this->T1_(c,m) 
    //                += -0.50 P(j,k) P(a,b) <n,m||c,d> R2(i,j) this->T1_(a,m) this->T3_(c,d,b,k,l,n) 
    //                += +1.00 P(j,k) P(a,b) <m,a||k,l> R2(i,j) this->T1_(b,m) 
    //                += +1.00 P(j,k) P(a,b) f(a,c) R2(i,j) this->T2_(c,b,k,l) 
    //                += -1.00 P(j,k) P(a,b) f(m,c) R2(i,j) this->T1_(a,m) this->T2_(c,b,k,l) 
    //                += +0.50 P(j,k) P(a,b) <m,a||c,d> R3(b,i,j,m) this->T2_(c,d,k,l) 
    //                += +1.00 P(j,k) P(a,b) <m,a||k,l> R3(b,i,j,m) 
    //                += +1.00 P(j,k) P(a,b) f(m,c) R3(b,i,j,m) this->T2_(c,a,k,l) 
    //                += -0.50 P(j,k) P(a,b) <n,m||c,d> R4(d,b,i,j,n,m) this->T2_(c,a,k,l) 
    //                += -1.00 P(j,k) P(a,b) <m,a||c,d> R3(d,i,j,m) this->T2_(c,b,k,l) 
    //                += +1.00 P(j,k) P(a,b) <n,m||c,d> R3(d,i,j,n) this->T1_(a,m) this->T2_(c,b,k,l) 
    //                += -0.50 P(j,k) P(a,b) <n,m||c,d> R3(b,i,j,n) this->T1_(a,m) this->T2_(c,d,k,l) 
    //                += -1.00 P(j,k) P(a,b) <n,m||c,d> R3(b,i,j,n) this->T2_(d,a,k,l) this->T1_(c,m) 
    //                += -0.50 P(j,k) P(a,b) <n,m||c,d> R3(b,i,j,n) this->T3_(c,d,a,k,l,m) 
    //                += -1.00 P(j,k) P(a,b) <n,m||k,l> R3(b,i,j,n) this->T1_(a,m) 
    sigmaR4("a,b,i,j,k,l") -= 0.50 * tmps_["23_voovoo"]("a,k,l,b,i,j");
    sigmaR4("a,b,i,j,k,l") += 0.50 * tmps_["23_voovoo"]("a,j,l,b,i,k");
    sigmaR4("a,b,i,j,k,l") += 0.50 * tmps_["23_voovoo"]("b,k,l,a,i,j");
    sigmaR4("a,b,i,j,k,l") -= 0.50 * tmps_["23_voovoo"]("b,j,l,a,i,k");
    TAmanager.free("voovoo", std::move(tmps_["23_voovoo"]));
    
    // flops: o4v2L1  = o5v2 o4v2 o5v2L1
    //  mems: o4v2L1  = o4v2 o4v2 o4v2L1
    tmps_.emplace(std::make_pair("24_vovooo", TAmanager.malloc<MatsT>("vovooo")));
    tmps_["24_vovooo"]("b,l,a,i,k,j")  = (reused_["118_voovoo"]("b,l,n,a,i,k") + this->T2_("a,b,l,m") * reused_["45_oooo"]("i,k,m,n")) * R2("j,n");
    
    // sigmaR4 += +1.00 P(k,l) <n,m||c,d> R2(k,n) this->T2_(c,a,i,j) this->T2_(d,b,l,m) 
    //                += +1.00 P(k,l) <n,m||c,d> R2(k,n) this->T2_(a,b,l,m) this->T1_(c,j) this->T1_(d,i) 
    sigmaR4("a,b,i,j,k,l") -= tmps_["24_vovooo"]("b,l,a,i,j,k");
    sigmaR4("a,b,i,j,k,l") += tmps_["24_vovooo"]("b,k,a,i,j,l");
    
    // sigmaR4 += +1.00 P(i,j) <n,m||c,d> R2(i,n) this->T2_(c,a,k,l) this->T2_(d,b,j,m) 
    //                += +1.00 P(i,j) <n,m||c,d> R2(i,n) this->T2_(a,b,j,m) this->T1_(c,l) this->T1_(d,k) 
    sigmaR4("a,b,i,j,k,l") -= tmps_["24_vovooo"]("b,j,a,k,l,i");
    sigmaR4("a,b,i,j,k,l") += tmps_["24_vovooo"]("b,i,a,k,l,j");
    
    // sigmaR4 += +1.00 P(i,l) <n,m||c,d> R2(i,n) this->T2_(c,a,j,k) this->T2_(d,b,l,m) 
    //                += +1.00 P(i,l) <n,m||c,d> R2(i,n) this->T2_(a,b,l,m) this->T1_(c,k) this->T1_(d,j) 
    sigmaR4("a,b,i,j,k,l") -= tmps_["24_vovooo"]("b,l,a,j,k,i");
    sigmaR4("a,b,i,j,k,l") += tmps_["24_vovooo"]("b,i,a,j,k,l");
    
    // sigmaR4 += -1.00 P(i,k) <n,m||c,d> R2(i,n) this->T2_(c,a,j,l) this->T2_(d,b,k,m) 
    //                += -1.00 P(i,k) <n,m||c,d> R2(i,n) this->T2_(a,b,k,m) this->T1_(c,l) this->T1_(d,j) 
    sigmaR4("a,b,i,j,k,l") += tmps_["24_vovooo"]("b,k,a,j,l,i");
    sigmaR4("a,b,i,j,k,l") -= tmps_["24_vovooo"]("b,i,a,j,l,k");
    
    // sigmaR4 += +1.00 P(j,k) <n,m||c,d> R2(j,n) this->T2_(c,a,i,l) this->T2_(d,b,k,m) 
    //                += +1.00 P(j,k) <n,m||c,d> R2(j,n) this->T2_(a,b,k,m) this->T1_(c,l) this->T1_(d,i) 
    sigmaR4("a,b,i,j,k,l") -= tmps_["24_vovooo"]("b,k,a,i,l,j");
    sigmaR4("a,b,i,j,k,l") += tmps_["24_vovooo"]("b,j,a,i,l,k");
    
    // sigmaR4 += -1.00 P(j,l) <n,m||c,d> R2(j,n) this->T2_(c,a,i,k) this->T2_(d,b,l,m) 
    //                += -1.00 P(j,l) <n,m||c,d> R2(j,n) this->T2_(a,b,l,m) this->T1_(c,k) this->T1_(d,i) 
    sigmaR4("a,b,i,j,k,l") += tmps_["24_vovooo"]("b,l,a,i,k,j");
    sigmaR4("a,b,i,j,k,l") -= tmps_["24_vovooo"]("b,j,a,i,k,l");
    TAmanager.free("vovooo", std::move(tmps_["24_vovooo"]));
    
    // flops: o4v2L1  = o4v2L1 o6v2L1 o4v2L1
    //  mems: o4v2L1  = o4v2L1 o4v2L1 o4v2L1
    tmps_.emplace(std::make_pair("25_oooovv", TAmanager.malloc<MatsT>("oooovv")));
    tmps_["25_oooovv"]("k,l,i,j,a,b")  = reused_["69_oovv"]("i,j,a,b") * R2("k,l");
    tmps_["25_oooovv"]("k,l,i,j,a,b") += R4("a,b,k,l,n,m") * reused_["45_oooo"]("i,j,m,n");
    
    // sigmaR4 += -0.50 <n,m||c,d> R2(k,l) this->T2_(a,b,n,m) this->T1_(c,j) this->T1_(d,i) 
    //                += -1.00 <a,b||c,d> R2(k,l) this->T1_(c,j) this->T1_(d,i) 
    //                += +1.00 <n,m||c,d> R2(k,l) this->T1_(a,m) this->T1_(b,n) this->T1_(c,j) this->T1_(d,i) 
    //                += -0.50 <n,m||c,d> R4(a,b,k,l,n,m) this->T1_(c,j) this->T1_(d,i) 
    sigmaR4("a,b,i,j,k,l") += 0.50 * tmps_["25_oooovv"]("k,l,i,j,a,b");
    
    // sigmaR4 += -0.50 <n,m||c,d> R2(i,l) this->T2_(a,b,n,m) this->T1_(c,k) this->T1_(d,j) 
    //                += -1.00 <a,b||c,d> R2(i,l) this->T1_(c,k) this->T1_(d,j) 
    //                += +1.00 <n,m||c,d> R2(i,l) this->T1_(a,m) this->T1_(b,n) this->T1_(c,k) this->T1_(d,j) 
    //                += -0.50 <n,m||c,d> R4(a,b,i,l,n,m) this->T1_(c,k) this->T1_(d,j) 
    sigmaR4("a,b,i,j,k,l") += 0.50 * tmps_["25_oooovv"]("i,l,j,k,a,b");
    
    // sigmaR4 += -0.50 P(j,k) <n,m||c,d> R2(i,j) this->T2_(a,b,n,m) this->T1_(c,l) this->T1_(d,k) 
    //                += -1.00 P(j,k) <a,b||c,d> R2(i,j) this->T1_(c,l) this->T1_(d,k) 
    //                += +1.00 P(j,k) <n,m||c,d> R2(i,j) this->T1_(a,m) this->T1_(b,n) this->T1_(c,l) this->T1_(d,k) 
    //                += -0.50 P(j,k) <n,m||c,d> R4(a,b,i,j,n,m) this->T1_(c,l) this->T1_(d,k) 
    sigmaR4("a,b,i,j,k,l") += 0.50 * tmps_["25_oooovv"]("i,j,k,l,a,b");
    sigmaR4("a,b,i,j,k,l") -= 0.50 * tmps_["25_oooovv"]("i,k,j,l,a,b");
    
    // sigmaR4 += -0.50 P(k,l) <n,m||c,d> R2(j,k) this->T2_(a,b,n,m) this->T1_(c,l) this->T1_(d,i) 
    //                += -1.00 P(k,l) <a,b||c,d> R2(j,k) this->T1_(c,l) this->T1_(d,i) 
    //                += +1.00 P(k,l) <n,m||c,d> R2(j,k) this->T1_(a,m) this->T1_(b,n) this->T1_(c,l) this->T1_(d,i) 
    //                += -0.50 P(k,l) <n,m||c,d> R4(a,b,j,k,n,m) this->T1_(c,l) this->T1_(d,i) 
    sigmaR4("a,b,i,j,k,l") += 0.50 * tmps_["25_oooovv"]("j,k,i,l,a,b");
    sigmaR4("a,b,i,j,k,l") -= 0.50 * tmps_["25_oooovv"]("j,l,i,k,a,b");
    TAmanager.free("oooovv", std::move(tmps_["25_oooovv"]));

/////////////////////////////////////////////////////////


  }

  template <typename MatsT>
  void EOMDIP_4h2pCCSDT<MatsT>::buildDiag(MatsT * diag, std::vector<double> eps) const {

    TAManager &TAmanager = TAManager::get();
    size_t n_v = TAmanager.getRange(vLabel_).extent();
    size_t n_o = TAmanager.getRange(oLabel_).extent();
    std::fill_n(diag, this->Hbar_dim, MatsT(0.0));

    for (auto j = 0; j < n_o; j++){
      for (auto i = 0; i < j; i++){
        diag[toCompoundS(i,j)] = - eps[i] - eps[j];
      }
    }
    MatsT * diag3 = diag + this->Hbar_dimension_offsets.at("ThreeBody");
    for (auto a = 0; a < n_v; a++){
      for (auto k = 0; k < n_o; k++){
        for (auto j = 0; j < k; j++){
          for (auto i = 0; i < j; i++){
            diag3[toCompoundD(a,i,j,k)] = eps[a+n_o] - eps[i] - eps[j] - eps[k];
          }
        }
      }
    }
    MatsT * diag4 = diag + this->Hbar_dimension_offsets.at("FourBody");
    for (auto b = 0; b < n_v; b++){
      for (auto a = 0; a < b; a++){
        for (auto l = 0; l < n_o; l++){
          for (auto k = 0; k < l; k++){
            for (auto j = 0; j < k; j++){
              for (auto i = 0; i < j; i++){
                diag4[toCompoundT(a,b,i,j,k,l)] = eps[a+n_o] + eps[b+n_o] - eps[i] - eps[j] - eps[k] - eps[l];
              }
            }
          }
        }
      }
    }
  }



  template <typename MatsT>
  void EOMDIP_4h2pCCSDT<MatsT>::runLambda(){} 


  template <typename MatsT>
  typename Davidson<dcomplex>::VecsGen_t EOMDIP_4h2pCCSDT<MatsT>::EmptyDavidsonVectorBuilder(){
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
  typename Davidson<dcomplex>::LinearTrans_t EOMDIP_4h2pCCSDT<MatsT>::DavidsonResidualBuilder(EOMCCEigenVecType &eigenVecType){
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
  typename Davidson<dcomplex>::LinearTrans_t EOMDIP_4h2pCCSDT<MatsT>::DavidsonPreconditionerBuilder(dcomplex * curEig, dcomplex * eomDiag){

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

          MatsT * Diag_IJ = eomDiag;
          TA::foreach_inplace(curB.get_tensor("TwoBody"), [iVec, curEig, Diag_IJ, this, PCsmall](TA::TensorZ &tile){
            const auto& lobound = tile.range().lobound();
            const auto& upbound = tile.range().upbound();

            dcomplex denom = 0.0;
            std::vector<std::size_t> x{0, 0};
            for(x[0] = lobound[0]; x[0] < upbound[0]; ++x[0]) {
              for(x[1] = lobound[1]; x[1] < upbound[1]; ++x[1]) {
                if (x[0] == x[1]) continue;
                denom = curEig[iVec] - Diag_IJ[toCompoundS(x[0], x[1])];
                if (std::abs(denom) >= PCsmall) tile[x] /= denom;
              }
            }
          });
          TA::get_default_world().gop.fence();

          dcomplex *diagD = eomDiag + nV_;

          MatsT * Diag_AIJK = eomDiag + this->Hbar_dimension_offsets.at("ThreeBody");
          TA::foreach_inplace(curB.get_tensor("ThreeBody"), [iVec, curEig, Diag_AIJK, this, PCsmall](TA::TensorZ &tile){
            const auto& lobound = tile.range().lobound();
            const auto& upbound = tile.range().upbound();

            dcomplex denom = 0.0;
            std::vector<std::size_t> x{0, 0, 0, 0};
            for(x[0] = lobound[0]; x[0] < upbound[0]; ++x[0]){
              size_t a = x[0];
              for(x[1] = lobound[1]; x[1] < upbound[1]; ++x[1]){
                for(x[2] = lobound[2]; x[2] < upbound[2]; ++x[2]){
                  for(x[3] = lobound[3]; x[3] < upbound[3]; ++x[3]){
                    if (x[1] == x[2] or x[2] == x[3] or x[3] == x[1])
                      continue;
                    size_t i = x[1];
                    size_t j = x[2];
                    size_t k = x[3];
                    denom = curEig[iVec] - Diag_AIJK[toCompoundD(a, i, j, k)];
                    if (std::abs(denom) >= PCsmall) tile[x] /= denom;
                  }
                }
              }
            }
          });
          TA::get_default_world().gop.fence();

          MatsT * Diag_ABIJKL = eomDiag + this->Hbar_dimension_offsets.at("FourBody");
          TA::foreach_inplace(curB.get_tensor("ThreeBody"), [iVec, curEig, Diag_ABIJKL, this, PCsmall](TA::TensorZ &tile){
            const auto& lobound = tile.range().lobound();
            const auto& upbound = tile.range().upbound();

            dcomplex denom = 0.0;
            std::vector<std::size_t> x{0, 0, 0, 0, 0, 0};
            for(x[0] = lobound[0]; x[0] < upbound[0]; ++x[0]){
              size_t a = x[0];
              for(x[1] = lobound[1]; x[1] < upbound[1]; ++x[1]){
                if (x[0] == x[1])
                  continue;
                size_t b = x[1];
                for(x[2] = lobound[2]; x[2] < upbound[2]; ++x[2]){
                  size_t i = x[2];
                  for(x[3] = lobound[3]; x[3] < upbound[3]; ++x[3]){
                    if (x[2] == x[3]) continue;
                    size_t j = x[3];
                    for(x[4] = lobound[4]; x[4] < upbound[4]; ++x[4]){
                       if (x[4] == x[2] or x[4] == x[2]) continue; 
                      size_t k = x[4];
                      for(x[5] = lobound[5]; x[5] < upbound[5]; ++x[5]){
                        if (x[2] == x[5] or x[3] == x[5] or x[4] == x[5]) continue;
                        size_t l = x[5];
                        denom = curEig[iVec] - Diag_ABIJKL[toCompoundT(a, b, i, j, k, l)];
                        if (std::abs(denom) >= PCsmall) tile[x] /= denom;
                      }
                    }
                  }
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
  EOMDIP_4h2pCCSDT<MatsT>::~EOMDIP_4h2pCCSDT() {

    TAManager &TAmanager = TAManager::get();

    TAmanager.free("vovo",std::move(reused_["1_vovo"]), true); 
    TAmanager.free("oovv",std::move(reused_["2_oovv"]), true); 
    TAmanager.free("ovov",std::move(reused_["3_ovov"]), true); 
    TAmanager.free("ovov",std::move(reused_["4_ovov"]), true); 
    TAmanager.free("oovo",std::move(reused_["5_oovo"]), true); 
    TAmanager.free("oovo",std::move(reused_["6_oovo"]), true); 
    TAmanager.free("oooo",std::move(reused_["7_oooo"]), true); 
    TAmanager.free("voov",std::move(reused_["8_voov"]), true); 
    TAmanager.free("vovv",std::move(reused_["9_vovv"]), true); 
    TAmanager.free("vovv",std::move(reused_["10_vovv"]), true); 
    TAmanager.free("vvoo",std::move(reused_["11_vvoo"]), true); 
    TAmanager.free("vvoo",std::move(reused_["12_vvoo"]), true); 
    TAmanager.free("vo",std::move(reused_["13_vo"]), true); 
    TAmanager.free("voov",std::move(reused_["14_voov"]), true); 
    TAmanager.free("vovo",std::move(reused_["15_vovo"]), true); 
    TAmanager.free("vvvo",std::move(reused_["16_vvvo"]), true); 
    TAmanager.free("vooo",std::move(reused_["17_vooo"]), true); 
    TAmanager.free("voov",std::move(reused_["18_voov"]), true); 
    TAmanager.free("vooo",std::move(reused_["19_vooo"]), true); 
    TAmanager.free("voov",std::move(reused_["20_voov"]), true); 
    TAmanager.free("vovo",std::move(reused_["21_vovo"]), true); 
    TAmanager.free("vooo",std::move(reused_["22_vooo"]), true); 
    TAmanager.free("vo",std::move(reused_["23_vo"]), true); 
    TAmanager.free("ooov",std::move(reused_["24_ooov"]), true); 
    TAmanager.free("ov",std::move(reused_["25_ov"]), true); 
    TAmanager.free("ov",std::move(reused_["26_ov"]), true); 
    TAmanager.free("ov",std::move(reused_["27_ov"]), true); 
    TAmanager.free("vo",std::move(reused_["28_vo"]), true); 
    TAmanager.free("vo",std::move(reused_["29_vo"]), true); 
    TAmanager.free("oo",std::move(reused_["30_oo"]), true); 
    TAmanager.free("vo",std::move(reused_["31_vo"]), true); 
    TAmanager.free("oo",std::move(reused_["32_oo"]), true); 
    TAmanager.free("ov",std::move(reused_["33_ov"]), true); 
    TAmanager.free("vo",std::move(reused_["34_vo"]), true); 
    TAmanager.free("vo",std::move(reused_["35_vo"]), true); 
    TAmanager.free("oo",std::move(reused_["36_oo"]), true); 
    TAmanager.free("vo",std::move(reused_["37_vo"]), true); 
    TAmanager.free("vo",std::move(reused_["38_vo"]), true); 
    TAmanager.free("vo",std::move(reused_["39_vo"]), true); 
    TAmanager.free("vo",std::move(reused_["40_vo"]), true); 
    TAmanager.free("vo",std::move(reused_["41_vo"]), true); 
    TAmanager.free("oo",std::move(reused_["42_oo"]), true); 
    TAmanager.free("oo",std::move(reused_["43_oo"]), true); 
    TAmanager.free("oo",std::move(reused_["44_oo"]), true); 
    TAmanager.free("oooo",std::move(reused_["45_oooo"]), true); 
    TAmanager.free("vooo",std::move(reused_["46_vooo"]), true); 
    TAmanager.free("vooo",std::move(reused_["47_vooo"]), true); 
    TAmanager.free("ovoo",std::move(reused_["48_ovoo"]), true); 
    TAmanager.free("ooov",std::move(reused_["49_ooov"]), true); 
    TAmanager.free("vooo",std::move(reused_["50_vooo"]), true); 
    TAmanager.free("vooo",std::move(reused_["51_vooo"]), true); 
    TAmanager.free("vooo",std::move(reused_["52_vooo"]), true); 
    TAmanager.free("oooo",std::move(reused_["53_oooo"]), true); 
    TAmanager.free("vooo",std::move(reused_["54_vooo"]), true); 
    TAmanager.free("oovo",std::move(reused_["55_oovo"]), true); 
    TAmanager.free("vvoo",std::move(reused_["56_vvoo"]), true); 
    TAmanager.free("voov",std::move(reused_["57_voov"]), true); 
    TAmanager.free("vvoo",std::move(reused_["58_vvoo"]), true); 
    TAmanager.free("vvoo",std::move(reused_["59_vvoo"]), true); 
    TAmanager.free("vooo",std::move(reused_["60_vooo"]), true); 
    TAmanager.free("voov",std::move(reused_["61_voov"]), true); 
    TAmanager.free("vvoo",std::move(reused_["62_vvoo"]), true); 
    TAmanager.free("vvoo",std::move(reused_["63_vvoo"]), true); 
    TAmanager.free("vvoo",std::move(reused_["64_vvoo"]), true); 
    TAmanager.free("vvoo",std::move(reused_["65_vvoo"]), true); 
    TAmanager.free("vovv",std::move(reused_["66_vovv"]), true); 
    TAmanager.free("oovv",std::move(reused_["67_oovv"]), true); 
    TAmanager.free("oovv",std::move(reused_["68_oovv"]), true); 
    TAmanager.free("oovv",std::move(reused_["69_oovv"]), true); 
    TAmanager.free("ovvo",std::move(reused_["70_ovvo"]), true); 
    TAmanager.free("vovo",std::move(reused_["71_vovo"]), true); 
    TAmanager.free("ovvo",std::move(reused_["72_ovvo"]), true); 
    TAmanager.free("ovvo",std::move(reused_["73_ovvo"]), true); 
    TAmanager.free("ovvo",std::move(reused_["74_ovvo"]), true); 
    TAmanager.free("ovvo",std::move(reused_["75_ovvo"]), true); 
    TAmanager.free("vovo",std::move(reused_["76_vovo"]), true); 
    TAmanager.free("ovvo",std::move(reused_["77_ovvo"]), true); 
    TAmanager.free("vovo",std::move(reused_["78_vovo"]), true); 
    TAmanager.free("vvoo",std::move(reused_["79_vvoo"]), true); 
    TAmanager.free("vvoo",std::move(reused_["80_vvoo"]), true); 
    TAmanager.free("vvoo",std::move(reused_["81_vvoo"]), true); 
    TAmanager.free("voov",std::move(reused_["82_voov"]), true); 
    TAmanager.free("vvoo",std::move(reused_["83_vvoo"]), true); 
    TAmanager.free("vvoo",std::move(reused_["84_vvoo"]), true); 
    TAmanager.free("voov",std::move(reused_["85_voov"]), true); 
    TAmanager.free("vvoo",std::move(reused_["86_vvoo"]), true); 
    TAmanager.free("vv",std::move(reused_["87_vv"]), true); 
    TAmanager.free("vooo",std::move(reused_["88_vooo"]), true); 
    TAmanager.free("vv",std::move(reused_["89_vv"]), true); 
    TAmanager.free("vo",std::move(reused_["90_vo"]), true); 
    TAmanager.free("vo",std::move(reused_["91_vo"]), true); 
    TAmanager.free("vo",std::move(reused_["92_vo"]), true); 
    TAmanager.free("vv",std::move(reused_["93_vv"]), true); 
    TAmanager.free("vv",std::move(reused_["94_vv"]), true); 
    TAmanager.free("oooo",std::move(reused_["95_oooo"]), true); 
    TAmanager.free("vvov",std::move(reused_["96_vvov"]), true); 
    TAmanager.free("vvov",std::move(reused_["97_vvov"]), true); 
    TAmanager.free("oovv",std::move(reused_["98_oovv"]), true); 
    TAmanager.free("ovvo",std::move(reused_["99_ovvo"]), true); 
    TAmanager.free("vvoo",std::move(reused_["100_vvoo"]), true); 
    TAmanager.free("oovvoo",std::move(reused_["101_oovvoo"]), true); 
    TAmanager.free("voovoo",std::move(reused_["102_voovoo"]), true); 
    TAmanager.free("oovvoo",std::move(reused_["103_oovvoo"]), true); 
    TAmanager.free("ovvo",std::move(reused_["104_ovvo"]), true); 
    TAmanager.free("vvoo",std::move(reused_["105_vvoo"]), true); 
    TAmanager.free("ovov",std::move(reused_["106_ovov"]), true); 
    TAmanager.free("ovov",std::move(reused_["107_ovov"]), true); 
    TAmanager.free("vooooo",std::move(reused_["108_vooooo"]), true); 
    TAmanager.free("oovvoo",std::move(reused_["109_oovvoo"]), true); 
    TAmanager.free("ovvooo",std::move(reused_["110_ovvooo"]), true); 
    TAmanager.free("vvvo",std::move(reused_["111_vvvo"]), true); 
    TAmanager.free("oovoov",std::move(reused_["112_oovoov"]), true); 
    TAmanager.free("oovoov",std::move(reused_["113_oovoov"]), true); 
    TAmanager.free("vvoooo",std::move(reused_["114_vvoooo"]), true); 
    TAmanager.free("vovooo",std::move(reused_["115_vovooo"]), true); 
    TAmanager.free("voovoo",std::move(reused_["116_voovoo"]), true); 
    TAmanager.free("voovoo",std::move(reused_["117_voovoo"]), true); 
    TAmanager.free("voovoo",std::move(reused_["118_voovoo"]), true); 
    TAmanager.free("vvoooo",std::move(reused_["119_vvoooo"]), true); 
  }


}
