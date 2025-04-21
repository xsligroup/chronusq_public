#include <coupledcluster/MBExpansion.hpp>
#include <coupledcluster.hpp>
namespace ChronusQ{
  template <typename MatsT>
  inline std::pair<size_t, double> EOMEA<MatsT>::toCompoundSS(size_t a, size_t b, size_t ldH) const {
    return std::make_pair(a + b * ldH, 1.0);
  }

  template <typename MatsT>
  inline std::pair<size_t, double> EOMEA<MatsT>::toCompoundSD(size_t a, size_t b,
                                                                         size_t c, size_t i, size_t ldH) const {
    size_t bci = toCompoundD(b,c,i);
    if (bci == outOfBound_)
      return std::make_pair(outOfBound_, 0.0);
    double sign = signD(b,c,i);
    return std::make_pair(a + bci * ldH, sign);
  }


  template <typename MatsT>
  inline std::pair<size_t, double> EOMEA<MatsT>::toCompoundDS(size_t a, size_t b, size_t i, size_t c,
                                                                         size_t ldH) const {
    size_t abi = toCompoundD(a,b,i);
    if (abi == outOfBound_)
      return std::make_pair(outOfBound_, 0.0);
    double sign = signD(a,b,i);
    return std::make_pair(abi + c * ldH, sign);
  }

  template <typename MatsT>
  inline std::pair<size_t, double> EOMEA<MatsT>::toCompoundDD(size_t a, size_t b, size_t i,
                                                                         size_t c, size_t d, size_t j, size_t ldH) const {
    size_t abi = toCompoundD(a,b,i), cdj = toCompoundD(c,d,j);
    if (abi == outOfBound_ or cdj == outOfBound_)
      return std::make_pair(outOfBound_, 0.0);
    double sign = signD(a,b,i);
    sign *= signD(c,d,j);
    return std::make_pair(abi + cdj * ldH, sign);
  }

  //for full_diagonalization
   template <typename MatsT>
  cqmatrix::Matrix<MatsT> EOMEA<MatsT>::buildHbar(bool includeGroundState) const{
    if (includeGroundState) CErr("EA should not include Ground State in the Hamiltonian.");

    TAManager &TAmanager = TAManager::get();
    size_t nV = TAmanager.getRange(vLabel_).extent();
    size_t nO = TAmanager.getRange(oLabel_).extent();

    cqmatrix::Matrix<MatsT> fullMat(this->Hbar_dim);
    fullMat.clear();

    MatsT * Hbar = fullMat.pointer();
    size_t ldH = fullMat.dimension();
    size_t nCol = fullMat.dimension();

    MatsT * HbarSS = Hbar;
    MatsT * HbarSD = HbarSS + nV_ * ldH;
    MatsT * HbarDS = HbarSS + nV_;
    MatsT * HbarDD = HbarSD + nV_;

    TArray H_vv = TAmanager.malloc<MatsT>("vv");  
    TArray H_vvvo = TAmanager.malloc<MatsT>("vvvo");  
    TArray H_vvov = TAmanager.malloc<MatsT>("vvvo");  
    TArray H_vvovvo = TAmanager.malloc<MatsT>("vvvvoo");  
 
    buildHbarTA(H_vv, H_vvvo, H_vvov, H_vvovvo);

    TA::foreach_inplace( H_vv, [this, &ldH, HbarSS ](TA::Tensor<MatsT>& tile) {

      const auto& lobound = tile.range().lobound();
      const auto& upbound = tile.range().upbound();

      std::size_t x[] = {0, 0};
      for(x[0] = lobound[0]; x[0] != upbound[0]; ++x[0])
        for(x[1] = lobound[1]; x[1] != upbound[1]; ++x[1]) {
              MatsT v = tile[x];

              size_t a = x[0], e = x[1];
                auto idx_sgn = toCompoundSS(a,e,ldH);
//std::cout<<"SS "<<i<<" "<<j<<" "<<k<<" "<<l<<" "<<idx_sgn.first<<" "<<idx_sgn.second<<" "<<v<<std::endl<<std::flush;
                if (isInBound(idx_sgn.first)) {
                  HbarSS[idx_sgn.first] = idx_sgn.second * v;
                }
        }
    });
    TA::get_default_world().gop.fence();
    TA::foreach_inplace( H_vvvo, [this, &ldH, HbarSD ](TA::Tensor<MatsT>& tile) {

      const auto& lobound = tile.range().lobound();
      const auto& upbound = tile.range().upbound();

      std::size_t x[] = {0, 0, 0, 0};
      for(x[0] = lobound[0]; x[0] != upbound[0]; ++x[0])
        for(x[1] = lobound[1]; x[1] != upbound[1]; ++x[1])
          for(x[2] = lobound[2]; x[2] != upbound[2]; ++x[2]) {
            if (x[1] == x[2])
              continue;
            for(x[3] = lobound[3]; x[3] != upbound[3]; ++x[3]) {
                  MatsT v = tile[x];
    
                  size_t a = x[0], b = x[1], c = x[2], i = x[3];
                  if (b<c){
                    auto idx_sgn = toCompoundSD(a,b,c,i,ldH);
//std::cout<<"SD "<<i<<" "<<j<<" "<<b<<" "<<k<<" "<<l<<" "<<m<<" "<<idx_sgn.first<<" "<<idx_sgn.second<<" "<<v<<std::endl<<std::flush;
                    if (isInBound(idx_sgn.first)) {
                      HbarSD[idx_sgn.first] = idx_sgn.second * v;
                    }
                  }
            }
        }
    });
    TA::get_default_world().gop.fence();
    TA::foreach_inplace( H_vvov, [this, &ldH, HbarDS ](TA::Tensor<MatsT>& tile) {

      const auto& lobound = tile.range().lobound();
      const auto& upbound = tile.range().upbound();

      std::size_t x[] = {0, 0, 0, 0};
      for(x[0] = lobound[0]; x[0] != upbound[0]; ++x[0])
        for(x[1] = lobound[1]; x[1] != upbound[1]; ++x[1]) {
          if (x[0] == x[1])
            continue;
          for(x[2] = lobound[2]; x[2] != upbound[2]; ++x[2]) 
            for(x[3] = lobound[3]; x[3] != upbound[3]; ++x[3]) {
                  MatsT v = tile[x];

                  //size_t i = x[0], j = x[1], k = x[2], a = x[3], l = x[4], m = x[5];
                  size_t a = x[0], b = x[1], c = x[2], i = x[3];
                  if (a<b){
                    auto idx_sgn = toCompoundDS(a,b,i,c,ldH);
//std::cout<<"DS "<<a<<" "<<i<<" "<<j<<" "<<k<<" "<<l<<" "<<m<<" "<<idx_sgn.first<<" "<<idx_sgn.second<<" "<<v<<std::endl<<std::flush;
                    if (isInBound(idx_sgn.first)) {
                      HbarDS[idx_sgn.first] = idx_sgn.second * v;
                    }
                  }
            }
        }
    });
    TA::get_default_world().gop.fence();
    TA::foreach_inplace( H_vvovvo, [this, &ldH, HbarDD ](TA::Tensor<MatsT>& tile) {

      const auto& lobound = tile.range().lobound();
      const auto& upbound = tile.range().upbound();

      std::size_t x[] = {0, 0, 0, 0, 0, 0};
      for(x[0] = lobound[0]; x[0] != upbound[0]; ++x[0])
        for(x[1] = lobound[1]; x[1] != upbound[1]; ++x[1]){ 
          if (x[0] == x[1])
            continue;
          for(x[2] = lobound[2]; x[2] != upbound[2]; ++x[2]) 
            for(x[3] = lobound[3]; x[3] != upbound[3]; ++x[3]) {
              if (x[2] == x[3])
                continue;
              for(x[4] = lobound[4]; x[4] != upbound[4]; ++x[4])
                for(x[5] = lobound[5]; x[5] != upbound[5]; ++x[5]) { 
                      MatsT v = tile[x];

                      //size_t i = x[0], j = x[1], k = x[2], a = x[3], l = x[4], m = x[5], n = x[6], b = x[7];
                      size_t a = x[0], b = x[1], c = x[2], d = x[3], i = x[4], j = x[5];
                      if (a<b && c<d){
                        auto idx_sgn = toCompoundDD(a,b,i,c,d,j,ldH);
//std::cout<<"DD "<<a<<" "<<i<<" "<<j<<" "<<k<<" "<<b<<" "<<l<<" "<<m<<" "<<n<<" "<<idx_sgn.first<<" "<<idx_sgn.second<<" "<<v<<std::endl<<std::flush;
                        if (isInBound(idx_sgn.first)) {
                          HbarDD[idx_sgn.first] = idx_sgn.second * v;
                        }
                      }
                    }
            }
        }
    });

    TA::get_default_world().gop.fence();
    TAmanager.free("vv", std::move(H_vv));
    TAmanager.free("vvvo", std::move(H_vvov));
    TAmanager.free("vvvo", std::move(H_vvvo));
    TAmanager.free("vvvvoo", std::move(H_vvovvo));
    return fullMat;
  }


  //for full_diagonalization
  template <typename MatsT>
  void EOMEA<MatsT>::buildHbarTA(TArray & H_vv, TArray & H_vvvo, TArray & H_vvov, TArray & H_vvovvo) const{

    TAManager &TAmanager = TAManager::get();

    TArray Id_oo = TAmanager.malloc<MatsT>("oo");
    TArray Id_vv = TAmanager.malloc<MatsT>("vv");
    TA::foreach_inplace(Id_oo, [&](TA::Tensor<MatsT> &tile){
      const auto& lobound = tile.range().lobound();
      const auto& upbound = tile.range().upbound();

      std::vector<std::size_t> x{0, 0};
      for(x[0] = lobound[0]; x[0] < upbound[0]; ++x[0])
        for(x[1] = lobound[1]; x[1] < upbound[1]; ++x[1])
          if(x[0]==x[1]) 
            tile[x] = 1.0;
          else 
            tile[x] = 0.0;
    });
    TA::foreach_inplace(Id_vv, [&](TA::Tensor<MatsT> &tile){
      const auto& lobound = tile.range().lobound();
      const auto& upbound = tile.range().upbound();

      std::vector<std::size_t> x{0, 0};
      for(x[0] = lobound[0]; x[0] < upbound[0]; ++x[0])
        for(x[1] = lobound[1]; x[1] < upbound[1]; ++x[1])
          if(x[0]==x[1]) 
            tile[x] = 1.0;
          else 
            tile[x] = 0.0;
    });
    TA::get_default_world().gop.fence();    


 // H_vv = +1.00 f(a,e)  // flops: o0v2 = o0v2 | mem: o0v2 = o0v2
H_vv("a,e")  = this->fockMatrix_ta["vv"]("a,e");

// H_vvov = +1.00 <a,b||e,i>  // flops: o1v3 = o1v3 | mem: o1v3 = o1v3
H_vvov("a,b,e,i")  = this->antiSymMoints["vvvo"]("a,b,e,i");

// H_vvvo = -1.00 <j,a||e,f>  // flops: o1v3 = o1v3 | mem: o1v3 = o1v3
H_vvvo("a,e,f,j")  = conj(this->antiSymMoints["vvvo"]("e,f,a,j"));

// H_vv += -0.50 d(a,e) <j,i||b,c> this->T1_(b,i) this->T1_(c,j)  // flops: o0v2 += o0v2 | mem: o0v2 += o0v2

// H_vv += -0.50 d(a,e) <j,i||j,i>  // flops: o0v2 += o0v2 | mem: o0v2 += o0v2

// H_vv += +1.00 d(a,e) f(i,b) this->T1_(b,i)  // flops: o0v2 += o0v2 | mem: o0v2 += o0v2

// H_vv += +0.25 d(a,e) <j,i||b,c> this->T2_(b,c,j,i)  // flops: o0v2 += o0v2 | mem: o0v2 += o0v2

// H_vv += +1.00 d(a,e) f(i,i)  // flops: o0v2 += o0v2 | mem: o0v2 += o0v2

// H_vvov += +1.00 d(a,e) f(b,i)  // flops: o1v3 += o1v3 | mem: o1v3 += o1v3
H_vvov("a,b,e,i") += Id_vv("a,e") * this->fockMatrix_ta["vo"]("b,i");

// H_vvov += -1.00 d(b,e) f(a,i)  // flops: o1v3 += o1v3 | mem: o1v3 += o1v3
H_vvov("a,b,e,i") -= Id_vv("b,e") * this->fockMatrix_ta["vo"]("a,i");

// H_vvvo += +1.00 d(a,e) f(j,f)  // flops: o1v3 += o1v3 | mem: o1v3 += o1v3
H_vvvo("a,e,f,j") += Id_vv("a,e") * this->fockMatrix_ta["ov"]("j,f");

// H_vvvo += -1.00 d(a,f) f(j,e)  // flops: o1v3 += o1v3 | mem: o1v3 += o1v3
H_vvvo("a,e,f,j") -= Id_vv("a,f") * this->fockMatrix_ta["ov"]("j,e");

// H_vvov += +1.00 f(j,e) this->T2_(a,b,i,j)  // flops: o1v3 += o2v3 | mem: o1v3 += o1v3
H_vvov("a,b,e,i") += this->fockMatrix_ta["ov"]("j,e") * this->T2_("a,b,i,j");

// H_vvov += -1.00 <a,b||c,e> this->T1_(c,i)  // flops: o1v3 += o1v4 | mem: o1v3 += o1v3
H_vvov("a,b,e,i") -= this->antiSymMoints["vvvv"]("a,b,c,e") * this->T1_("c,i");

// H_vvov += +0.50 <k,j||e,i> this->T2_(a,b,k,j)  // flops: o1v3 += o3v3 | mem: o1v3 += o1v3
H_vvov("a,b,e,i") += 0.50 * conj(this->antiSymMoints["vooo"]("e,i,k,j")) * this->T2_("a,b,k,j");

// H_vvovvo = +1.00 d(i,j) <a,b||e,f>  // flops: o2v4 = o2v4 | mem: o2v4 = o2v4
H_vvovvo("a,b,e,f,i,j")  = Id_oo("i,j") * this->antiSymMoints["vvvv"]("a,b,e,f");

// H_vvovvo += -1.00 d(a,f) <j,b||e,i>  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
H_vvovvo("a,b,e,f,i,j") += Id_vv("a,f") * this->antiSymMoints["vovo"]("b,j,e,i");

// H_vvovvo += +1.00 d(b,f) <j,a||e,i>  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
H_vvovvo("a,b,e,f,i,j") -= Id_vv("b,f") * this->antiSymMoints["vovo"]("a,j,e,i");

// H_vvovvo += +1.00 d(a,e) <j,b||f,i>  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
H_vvovvo("a,b,e,f,i,j") -= Id_vv("a,e") * this->antiSymMoints["vovo"]("b,j,f,i");

// H_vvovvo += -1.00 d(b,e) <j,a||f,i>  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
H_vvovvo("a,b,e,f,i,j") += Id_vv("b,e") * this->antiSymMoints["vovo"]("a,j,f,i");

// H_vvovvo += -1.00 <j,k||e,f> this->T2_(a,b,i,k)  // flops: o2v4 += o3v4 | mem: o2v4 += o2v4
H_vvovvo("a,b,e,f,i,j") -= conj(this->antiSymMoints["vvoo"]("e,f,j,k")) * this->T2_("a,b,i,k");

// tmps_[1_vvvo](a,e,b,i) = 1.00 eri[vovv](b,j,c,e) * this->T2_(c,a,i,j) // flops: o1v3 = o2v4 | mem: o1v3 = o1v3
tmps_.emplace(std::make_pair("vvvo_1", TAmanager.malloc<MatsT>("vvvo")));
tmps_["1_vvvo"]("a,e,b,i")  = conj(this->antiSymMoints["vvvo"]("c,e,b,j")) * this->T2_("c,a,i,j");

// H_vvov += -1.00 P(a,b) <j,a||c,e> this->T2_(c,b,i,j)  // flops: o1v3 += o1v3 | mem: o1v3 += o1v3
H_vvov("a,b,e,i") += tmps_["1_vvvo"]("b,e,a,i");
H_vvov("a,b,e,i") -= tmps_["1_vvvo"]("a,e,b,i");
TAmanager.free("vvvo", std::move(tmps_["1_vvvo"]));

// tmps_[2_vvoo](a,e,i,j) = 1.00 eri[oovv](j,k,c,e) * this->T2_(c,a,i,k) // flops: o2v2 = o3v3 | mem: o2v2 = o2v2
tmps_.emplace(std::make_pair("vvoo_2", TAmanager.malloc<MatsT>("vvoo")));
tmps_["2_vvoo"]("a,e,i,j")  = conj(this->antiSymMoints["vvoo"]("c,e,j,k")) * this->T2_("c,a,i,k");

// H_vvovvo += -1.00 d(a,f) <j,k||c,e> this->T2_(c,b,i,k)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
H_vvovvo("a,b,e,f,i,j") -= tmps_["2_vvoo"]("b,e,i,j") * Id_vv("a,f");

// H_vvovvo += +1.00 d(b,f) <j,k||c,e> this->T2_(c,a,i,k)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
H_vvovvo("a,b,e,f,i,j") += tmps_["2_vvoo"]("a,e,i,j") * Id_vv("b,f");

// H_vvovvo += +1.00 d(a,e) <j,k||c,f> this->T2_(c,b,i,k)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
H_vvovvo("a,b,e,f,i,j") += tmps_["2_vvoo"]("b,f,i,j") * Id_vv("a,e");

// H_vvovvo += -1.00 d(b,e) <j,k||c,f> this->T2_(c,a,i,k)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
H_vvovvo("a,b,e,f,i,j") -= tmps_["2_vvoo"]("a,f,i,j") * Id_vv("b,e");
TAmanager.free("vvoo", std::move(tmps_["2_vvoo"]));

// tmps_[3_vvoo](a,e,i,j) = 1.00 eri[oovv](k,j,c,e) * this->T2_(c,a,i,k) // flops: o2v2 = o3v3 | mem: o2v2 = o2v2
tmps_.emplace(std::make_pair("vvoo_3", TAmanager.malloc<MatsT>("vvoo")));
tmps_["3_vvoo"]("a,e,i,j")  = conj(this->antiSymMoints["vvoo"]("c,e,k,j")) * this->T2_("c,a,i,k");

// tmps_[27_vvvo](e,b,a,i) = 1.00 this->T1_(a,j) * eri[oovv](k,j,c,e) * this->T2_(c,b,i,k) // flops: o1v3 = o2v3 | mem: o1v3 = o1v3
tmps_.emplace(std::make_pair("vvvo_27", TAmanager.malloc<MatsT>("vvvo")));
tmps_["27_vvvo"]("e,b,a,i")  = this->T1_("a,j") * tmps_["3_vvoo"]("b,e,i,j");
TAmanager.free("vvoo", std::move(tmps_["3_vvoo"]));
H_vvov("a,b,e,i") -= tmps_["27_vvvo"]("e,a,b,i");

// H_vvov += +1.00 P(a,b) <k,j||c,e> this->T1_(a,j) this->T2_(c,b,i,k)  // flops: o1v3 += o1v3 | mem: o1v3 += o1v3
H_vvov("a,b,e,i") += tmps_["27_vvvo"]("e,b,a,i");
TAmanager.free("vvvo", std::move(tmps_["27_vvvo"]));

// tmps_[4_vvvv](b,a,e,c) = 1.00 eri[oovv](k,j,c,e) * this->T2_(a,b,k,j) // flops: o0v4 = o2v4 | mem: o0v4 = o0v4
tmps_.emplace(std::make_pair("vvvv_4", TAmanager.malloc<MatsT>("vvvv")));
tmps_["4_vvvv"]("b,a,e,c")  = conj(this->antiSymMoints["vvoo"]("c,e,k,j")) * this->T2_("a,b,k,j");

// H_vvovvo += +0.50 d(i,j) <l,k||e,f> this->T2_(a,b,l,k)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
H_vvovvo("a,b,e,f,i,j") += 0.50 * tmps_["4_vvvv"]("b,a,f,e") * Id_oo("i,j");
TAmanager.free("vvvv", std::move(tmps_["4_vvvv"]));

// tmps_[5_vvvv](a,e,c,b) = 1.00 eri[vovv](b,j,c,e) * this->T1_(a,j) // flops: o0v4 = o1v4 | mem: o0v4 = o0v4
tmps_.emplace(std::make_pair("vvvv_5", TAmanager.malloc<MatsT>("vvvv")));
tmps_["5_vvvv"]("a,e,c,b")  = conj(this->antiSymMoints["vvvo"]("c,e,b,j")) * this->T1_("a,j");

// H_vvovvo += +1.00 P(a,b) d(i,j) <k,a||e,f> this->T1_(b,k)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
H_vvovvo("a,b,e,f,i,j") -= tmps_["5_vvvv"]("b,f,e,a") * Id_oo("i,j");
H_vvovvo("a,b,e,f,i,j") += tmps_["5_vvvv"]("a,f,e,b") * Id_oo("i,j");
TAmanager.free("vvvv", std::move(tmps_["5_vvvv"]));

// tmps_[6_vvoo](e,b,i,j) = 1.00 eri[vovv](b,j,c,e) * this->T1_(c,i) // flops: o2v2 = o2v3 | mem: o2v2 = o2v2
tmps_.emplace(std::make_pair("vvoo_6", TAmanager.malloc<MatsT>("vvoo")));
tmps_["6_vvoo"]("e,b,i,j")  = conj(this->antiSymMoints["vvvo"]("c,e,b,j")) * this->T1_("c,i");

// H_vvovvo += +1.00 d(b,e) <j,a||c,f> this->T1_(c,i)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
H_vvovvo("a,b,e,f,i,j") -= tmps_["6_vvoo"]("f,a,i,j") * Id_vv("b,e");

// H_vvovvo += +1.00 d(a,f) <j,b||c,e> this->T1_(c,i)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
H_vvovvo("a,b,e,f,i,j") -= tmps_["6_vvoo"]("e,b,i,j") * Id_vv("a,f");

// H_vvovvo += -1.00 d(a,e) <j,b||c,f> this->T1_(c,i)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
H_vvovvo("a,b,e,f,i,j") += tmps_["6_vvoo"]("f,b,i,j") * Id_vv("a,e");

// H_vvovvo += -1.00 d(b,f) <j,a||c,e> this->T1_(c,i)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
H_vvovvo("a,b,e,f,i,j") += tmps_["6_vvoo"]("e,a,i,j") * Id_vv("b,f");

// tmps_[28_vvvo](a,e,b,i) = 1.00 this->T1_(b,j) * eri[vovv](a,j,c,e) * this->T1_(c,i) // flops: o1v3 = o2v3 | mem: o1v3 = o1v3
tmps_.emplace(std::make_pair("vvvo_28", TAmanager.malloc<MatsT>("vvvo")));
tmps_["28_vvvo"]("a,e,b,i")  = this->T1_("b,j") * tmps_["6_vvoo"]("e,a,i,j");
TAmanager.free("vvoo", std::move(tmps_["6_vvoo"]));

// H_vvov += -1.00 P(a,b) <j,a||c,e> this->T1_(c,i) this->T1_(b,j)  // flops: o1v3 += o1v3 | mem: o1v3 += o1v3
H_vvov("a,b,e,i") += tmps_["28_vvvo"]("a,e,b,i");
H_vvov("a,b,e,i") -= tmps_["28_vvvo"]("b,e,a,i");
TAmanager.free("vvvo", std::move(tmps_["28_vvvo"]));

// tmps_[7_vv](a,e) = 0.50 eri[oovv](l,k,c,e) * this->T2_(c,a,l,k) // flops: o0v2 = o2v3 | mem: o0v2 = o0v2
tmps_.emplace(std::make_pair("vv_7", TAmanager.malloc<MatsT>("vv")));
tmps_["7_vv"]("a,e")  = 0.50 * conj(this->antiSymMoints["vvoo"]("c,e,l,k")) * this->T2_("c,a,l,k");

// H_vv += -0.50 <j,i||b,e> this->T2_(b,a,j,i)  // flops: o0v2 += o0v2 | mem: o0v2 += o0v2
H_vv("a,e") -= tmps_["7_vv"]("a,e");

// tmps_[17_vvoo](f,a,j,i) = 1.00 Id[oo](i,j) * Id[vv](a,f) // flops: o2v2 = o2v2 | mem: o2v2 = o2v2
tmps_.emplace(std::make_pair("vvoo_17", TAmanager.malloc<MatsT>("vvoo")));
tmps_["17_vvoo"]("f,a,j,i")  = Id_oo("i,j") * Id_vv("a,f");

// H_vvovvo += -0.50 d(a,e) d(i,j) <l,k||c,f> this->T2_(c,b,l,k)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
H_vvovvo("a,b,e,f,i,j") -= tmps_["7_vv"]("b,f") * tmps_["17_vvoo"]("e,a,j,i");

// H_vvovvo += +0.50 d(a,f) d(i,j) <l,k||c,e> this->T2_(c,b,l,k)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
H_vvovvo("a,b,e,f,i,j") += tmps_["7_vv"]("b,e") * tmps_["17_vvoo"]("f,a,j,i");

// H_vvovvo += -0.50 d(b,f) d(i,j) <l,k||c,e> this->T2_(c,a,l,k)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
H_vvovvo("a,b,e,f,i,j") -= tmps_["7_vv"]("a,e") * tmps_["17_vvoo"]("f,b,j,i");

// H_vvovvo += +0.50 d(b,e) d(i,j) <l,k||c,f> this->T2_(c,a,l,k)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
H_vvovvo("a,b,e,f,i,j") += tmps_["7_vv"]("a,f") * tmps_["17_vvoo"]("e,b,j,i");
TAmanager.free("vv", std::move(tmps_["7_vv"]));

// tmps_[8_vo](a,i) = 0.50 eri[vovv](a,j,c,d) * this->T2_(c,d,i,j) // flops: o1v1 = o2v3 | mem: o1v1 = o1v1
tmps_.emplace(std::make_pair("vo_8", TAmanager.malloc<MatsT>("vo")));
tmps_["8_vo"]("a,i")  = 0.50 * conj(this->antiSymMoints["vvvo"]("c,d,a,j")) * this->T2_("c,d,i,j");

// H_vvov += +0.50 d(b,e) <j,a||c,d> this->T2_(c,d,i,j)  // flops: o1v3 += o1v3 | mem: o1v3 += o1v3
H_vvov("a,b,e,i") -= tmps_["8_vo"]("a,i") * Id_vv("b,e");

// H_vvov += -0.50 d(a,e) <j,b||c,d> this->T2_(c,d,i,j)  // flops: o1v3 += o1v3 | mem: o1v3 += o1v3
H_vvov("a,b,e,i") += tmps_["8_vo"]("b,i") * Id_vv("a,e");
TAmanager.free("vo", std::move(tmps_["8_vo"]));

// tmps_[9_vvvo](a,f,c,j) = 1.00 eri[oovv](j,k,c,f) * this->T1_(a,k) // flops: o1v3 = o2v3 | mem: o1v3 = o1v3
tmps_.emplace(std::make_pair("vvvo_9", TAmanager.malloc<MatsT>("vvvo")));
tmps_["9_vvvo"]("a,f,c,j")  = conj(this->antiSymMoints["vvoo"]("c,f,j,k")) * this->T1_("a,k");

// H_vvvo += +1.00 <j,i||e,f> this->T1_(a,i)  // flops: o1v3 += o1v3 | mem: o1v3 += o1v3
H_vvvo("a,e,f,j") += tmps_["9_vvvo"]("a,f,e,j");

// H_vvovvo += -1.00 d(i,j) <l,k||e,f> this->T1_(a,k) this->T1_(b,l)  // flops: o2v4 += o1v4 o2v4 | mem: o2v4 += o0v4 o2v4
H_vvovvo("a,b,e,f,i,j") -= tmps_["9_vvvo"]("a,f,e,l") * this->T1_("b,l") * Id_oo("i,j");
TAmanager.free("vvvo", std::move(tmps_["9_vvvo"]));

// tmps_[10_vvvo](b,e,a,i) = 1.00 eri[vovo](a,j,e,i) * this->T1_(b,j) // flops: o1v3 = o2v3 | mem: o1v3 = o1v3
tmps_.emplace(std::make_pair("vvvo_10", TAmanager.malloc<MatsT>("vvvo")));
tmps_["10_vvvo"]("b,e,a,i")  = this->antiSymMoints["vovo"]("a,j,e,i") * this->T1_("b,j");
H_vvov("a,b,e,i") += tmps_["10_vvvo"]("a,e,b,i");

// H_vvov += +1.00 P(a,b) <j,a||e,i> this->T1_(b,j)  // flops: o1v3 += o1v3 | mem: o1v3 += o1v3
H_vvov("a,b,e,i") -= tmps_["10_vvvo"]("b,e,a,i");
TAmanager.free("vvvo", std::move(tmps_["10_vvvo"]));

// tmps_[11_vooo](f,i,k,j) = 1.00 eri[oovv](j,k,c,f) * this->T1_(c,i) // flops: o3v1 = o3v2 | mem: o3v1 = o3v1
tmps_.emplace(std::make_pair("vooo_11", TAmanager.malloc<MatsT>("vooo")));
tmps_["11_vooo"]("f,i,k,j")  = conj(this->antiSymMoints["vvoo"]("c,f,j,k")) * this->T1_("c,i");

// H_vvov += -0.50 <k,j||c,e> this->T1_(c,i) this->T2_(a,b,k,j)  // flops: o1v3 += o3v3 | mem: o1v3 += o1v3
H_vvov("a,b,e,i") -= 0.50 * tmps_["11_vooo"]("e,i,j,k") * this->T2_("a,b,k,j");

// tmps_[29_vvoo](f,a,j,i) = 1.00 this->T1_(a,k) * eri[oovv](j,k,c,f) * this->T1_(c,i) // flops: o2v2 = o3v2 | mem: o2v2 = o2v2
tmps_.emplace(std::make_pair("vvoo_29", TAmanager.malloc<MatsT>("vvoo")));
tmps_["29_vvoo"]("f,a,j,i")  = this->T1_("a,k") * tmps_["11_vooo"]("f,i,k,j");

// H_vvov += +1.00 <k,j||c,e> this->T1_(c,i) this->T1_(a,j) this->T1_(b,k)  // flops: o1v3 += o2v3 | mem: o1v3 += o1v3
H_vvov("a,b,e,i") += tmps_["29_vvoo"]("e,a,k,i") * this->T1_("b,k");

// H_vvovvo += +1.00 d(b,f) <j,k||c,e> this->T1_(c,i) this->T1_(a,k)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
H_vvovvo("a,b,e,f,i,j") += tmps_["29_vvoo"]("e,a,j,i") * Id_vv("b,f");

// H_vvovvo += -1.00 d(a,f) <j,k||c,e> this->T1_(c,i) this->T1_(b,k)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
H_vvovvo("a,b,e,f,i,j") -= tmps_["29_vvoo"]("e,b,j,i") * Id_vv("a,f");

// H_vvovvo += -1.00 d(b,e) <j,k||c,f> this->T1_(c,i) this->T1_(a,k)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
H_vvovvo("a,b,e,f,i,j") -= tmps_["29_vvoo"]("f,a,j,i") * Id_vv("b,e");

// H_vvovvo += +1.00 d(a,e) <j,k||c,f> this->T1_(c,i) this->T1_(b,k)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
H_vvovvo("a,b,e,f,i,j") += tmps_["29_vvoo"]("f,b,j,i") * Id_vv("a,e");
TAmanager.free("vvoo", std::move(tmps_["29_vvoo"]));

// tmps_[30_vo](a,i) = 0.50 this->T2_(d,a,k,j) * eri[oovv](k,j,c,d) * this->T1_(c,i) // flops: o1v1 = o3v2 | mem: o1v1 = o1v1
tmps_.emplace(std::make_pair("vo_30", TAmanager.malloc<MatsT>("vo")));
tmps_["30_vo"]("a,i")  = 0.50 * this->T2_("d,a,k,j") * tmps_["11_vooo"]("d,i,j,k");
TAmanager.free("vooo", std::move(tmps_["11_vooo"]));

// H_vvov += +0.50 d(a,e) <k,j||c,d> this->T1_(c,i) this->T2_(d,b,k,j)  // flops: o1v3 += o1v3 | mem: o1v3 += o1v3
H_vvov("a,b,e,i") += tmps_["30_vo"]("b,i") * Id_vv("a,e");

// H_vvov += -0.50 d(b,e) <k,j||c,d> this->T1_(c,i) this->T2_(d,a,k,j)  // flops: o1v3 += o1v3 | mem: o1v3 += o1v3
H_vvov("a,b,e,i") -= tmps_["30_vo"]("a,i") * Id_vv("b,e");
TAmanager.free("vo", std::move(tmps_["30_vo"]));

// tmps_[12_vvoo](a,f,i,j) = 1.00 eri[oovo](j,k,f,i) * this->T1_(a,k) // flops: o2v2 = o3v2 | mem: o2v2 = o2v2
tmps_.emplace(std::make_pair("vvoo_12", TAmanager.malloc<MatsT>("vvoo")));
tmps_["12_vvoo"]("a,f,i,j")  = conj(this->antiSymMoints["vooo"]("f,i,j,k")) * this->T1_("a,k");

// H_vvov += -1.00 <k,j||e,i> this->T1_(a,j) this->T1_(b,k)  // flops: o1v3 += o2v3 | mem: o1v3 += o1v3
H_vvov("a,b,e,i") -= tmps_["12_vvoo"]("a,e,i,k") * this->T1_("b,k");

// H_vvovvo += -1.00 d(a,e) <j,k||f,i> this->T1_(b,k)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
H_vvovvo("a,b,e,f,i,j") -= tmps_["12_vvoo"]("b,f,i,j") * Id_vv("a,e");

// H_vvovvo += -1.00 d(b,f) <j,k||e,i> this->T1_(a,k)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
H_vvovvo("a,b,e,f,i,j") -= tmps_["12_vvoo"]("a,e,i,j") * Id_vv("b,f");

// H_vvovvo += +1.00 d(a,f) <j,k||e,i> this->T1_(b,k)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
H_vvovvo("a,b,e,f,i,j") += tmps_["12_vvoo"]("b,e,i,j") * Id_vv("a,f");

// H_vvovvo += +1.00 d(b,e) <j,k||f,i> this->T1_(a,k)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
H_vvovvo("a,b,e,f,i,j") += tmps_["12_vvoo"]("a,f,i,j") * Id_vv("b,e");
TAmanager.free("vvoo", std::move(tmps_["12_vvoo"]));

// tmps_[13_oo](i,j) = 0.50 eri[oovv](j,k,c,d) * this->T2_(c,d,i,k) // flops: o2v0 = o3v2 | mem: o2v0 = o2v0
tmps_.emplace(std::make_pair("oo_13", TAmanager.malloc<MatsT>("oo")));
tmps_["13_oo"]("i,j")  = 0.50 * conj(this->antiSymMoints["vvoo"]("c,d,j,k")) * this->T2_("c,d,i,k");

// tmps_[40_vvoo](e,a,j,i) = 1.00 Id[vv](a,e) * eri[oovv](j,k,c,d) * this->T2_(c,d,i,k) // flops: o2v2 = o2v2 | mem: o2v2 = o2v2
tmps_.emplace(std::make_pair("vvoo_40", TAmanager.malloc<MatsT>("vvoo")));
tmps_["40_vvoo"]("e,a,j,i")  = Id_vv("a,e") * tmps_["13_oo"]("i,j");
TAmanager.free("oo", std::move(tmps_["13_oo"]));

// H_vvovvo += -0.50 d(a,e) d(b,f) <j,k||c,d> this->T2_(c,d,i,k)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
H_vvovvo("a,b,e,f,i,j") -= tmps_["40_vvoo"]("e,a,j,i") * Id_vv("b,f");

// H_vvovvo += +0.50 d(b,e) d(a,f) <j,k||c,d> this->T2_(c,d,i,k)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
H_vvovvo("a,b,e,f,i,j") += tmps_["40_vvoo"]("e,b,j,i") * Id_vv("a,f");
TAmanager.free("vvoo", std::move(tmps_["40_vvoo"]));

// tmps_[14_oo](i,j) = 0.50 eri[oovv](k,j,c,d) * this->T2_(c,d,i,k) // flops: o2v0 = o3v2 | mem: o2v0 = o2v0
tmps_.emplace(std::make_pair("oo_14", TAmanager.malloc<MatsT>("oo")));
tmps_["14_oo"]("i,j")  = 0.50 * conj(this->antiSymMoints["vvoo"]("c,d,k,j")) * this->T2_("c,d,i,k");

// tmps_[37_vo](b,i) = 1.00 this->T1_(b,j) * eri[oovv](k,j,c,d) * this->T2_(c,d,i,k) // flops: o1v1 = o2v1 | mem: o1v1 = o1v1
tmps_.emplace(std::make_pair("vo_37", TAmanager.malloc<MatsT>("vo")));
tmps_["37_vo"]("b,i")  = this->T1_("b,j") * tmps_["14_oo"]("i,j");
TAmanager.free("oo", std::move(tmps_["14_oo"]));

// H_vvov += -0.50 d(b,e) <k,j||c,d> this->T1_(a,j) this->T2_(c,d,i,k)  // flops: o1v3 += o1v3 | mem: o1v3 += o1v3
H_vvov("a,b,e,i") -= tmps_["37_vo"]("a,i") * Id_vv("b,e");

// H_vvov += +0.50 d(a,e) <k,j||c,d> this->T1_(b,j) this->T2_(c,d,i,k)  // flops: o1v3 += o1v3 | mem: o1v3 += o1v3
H_vvov("a,b,e,i") += tmps_["37_vo"]("b,i") * Id_vv("a,e");
TAmanager.free("vo", std::move(tmps_["37_vo"]));

// tmps_[15_vo](a,i) = 0.50 eri[oovo](k,j,c,i) * this->T2_(c,a,k,j) // flops: o1v1 = o3v2 | mem: o1v1 = o1v1
tmps_.emplace(std::make_pair("vo_15", TAmanager.malloc<MatsT>("vo")));
tmps_["15_vo"]("a,i")  = 0.50 * conj(this->antiSymMoints["vooo"]("c,i,k,j")) * this->T2_("c,a,k,j");

// H_vvov += +0.50 d(b,e) <k,j||c,i> this->T2_(c,a,k,j)  // flops: o1v3 += o1v3 | mem: o1v3 += o1v3
H_vvov("a,b,e,i") += tmps_["15_vo"]("a,i") * Id_vv("b,e");

// H_vvov += -0.50 d(a,e) <k,j||c,i> this->T2_(c,b,k,j)  // flops: o1v3 += o1v3 | mem: o1v3 += o1v3
H_vvov("a,b,e,i") -= tmps_["15_vo"]("b,i") * Id_vv("a,e");
TAmanager.free("vo", std::move(tmps_["15_vo"]));

// tmps_[16_vv](d,a) = 1.00 eri[vovv](a,j,c,d) * this->T1_(c,j) // flops: o0v2 = o1v3 | mem: o0v2 = o0v2
tmps_.emplace(std::make_pair("vv_16", TAmanager.malloc<MatsT>("vv")));
tmps_["16_vv"]("d,a")  = conj(this->antiSymMoints["vvvo"]("c,d,a,j")) * this->T1_("c,j");

// H_vv += +1.00 <i,a||b,e> this->T1_(b,i)  // flops: o0v2 += o0v2 | mem: o0v2 += o0v2
H_vv("a,e") -= tmps_["16_vv"]("e,a");

// H_vvovvo += -1.00 d(b,e) d(i,j) <k,a||c,f> this->T1_(c,k)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
H_vvovvo("a,b,e,f,i,j") += tmps_["16_vv"]("f,a") * tmps_["17_vvoo"]("e,b,j,i");

// H_vvovvo += +1.00 d(a,e) d(i,j) <k,b||c,f> this->T1_(c,k)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
H_vvovvo("a,b,e,f,i,j") -= tmps_["16_vv"]("f,b") * tmps_["17_vvoo"]("e,a,j,i");

// H_vvovvo += +1.00 d(b,f) d(i,j) <k,a||c,e> this->T1_(c,k)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
H_vvovvo("a,b,e,f,i,j") -= tmps_["16_vv"]("e,a") * tmps_["17_vvoo"]("f,b,j,i");

// H_vvovvo += -1.00 d(a,f) d(i,j) <k,b||c,e> this->T1_(c,k)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
H_vvovvo("a,b,e,f,i,j") += tmps_["16_vv"]("e,b") * tmps_["17_vvoo"]("f,a,j,i");

// tmps_[33_vo](b,i) = 1.00 this->T1_(d,i) * eri[vovv](b,j,c,d) * this->T1_(c,j) // flops: o1v1 = o1v2 | mem: o1v1 = o1v1
tmps_.emplace(std::make_pair("vo_33", TAmanager.malloc<MatsT>("vo")));
tmps_["33_vo"]("b,i")  = this->T1_("d,i") * tmps_["16_vv"]("d,b");
TAmanager.free("vv", std::move(tmps_["16_vv"]));

// H_vvov += -1.00 d(b,e) <j,a||c,d> this->T1_(c,j) this->T1_(d,i)  // flops: o1v3 += o1v3 | mem: o1v3 += o1v3
H_vvov("a,b,e,i") += tmps_["33_vo"]("a,i") * Id_vv("b,e");

// H_vvov += +1.00 d(a,e) <j,b||c,d> this->T1_(c,j) this->T1_(d,i)  // flops: o1v3 += o1v3 | mem: o1v3 += o1v3
H_vvov("a,b,e,i") -= tmps_["33_vo"]("b,i") * Id_vv("a,e");
TAmanager.free("vo", std::move(tmps_["33_vo"]));

// H_vvovvo += +1.00 d(b,f) d(i,j) f(a,e)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
H_vvovvo("a,b,e,f,i,j") += tmps_["17_vvoo"]("f,b,j,i") * this->fockMatrix_ta["vv"]("a,e");

// H_vvovvo += -1.00 d(b,e) d(i,j) f(a,f)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
H_vvovvo("a,b,e,f,i,j") -= tmps_["17_vvoo"]("e,b,j,i") * this->fockMatrix_ta["vv"]("a,f");

// H_vvovvo += -1.00 d(a,f) d(i,j) f(b,e)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
H_vvovvo("a,b,e,f,i,j") -= tmps_["17_vvoo"]("f,a,j,i") * this->fockMatrix_ta["vv"]("b,e");

// H_vvovvo += +1.00 d(a,e) d(i,j) f(b,f)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
H_vvovvo("a,b,e,f,i,j") += tmps_["17_vvoo"]("e,a,j,i") * this->fockMatrix_ta["vv"]("b,f");

// H_vvovvo += +0.50 d(b,e) d(a,f) d(i,j) <l,k||l,k>  // flops: o2v4 += o0v2 o2v4 | mem: o2v4 += o0v2 o2v4

// H_vvovvo += +0.25 d(a,e) d(b,f) d(i,j) <l,k||c,d> this->T2_(c,d,l,k)  // flops: o2v4 += o0v2 o2v4 | mem: o2v4 += o0v2 o2v4

// H_vvovvo += +0.50 d(b,e) d(a,f) d(i,j) <l,k||c,d> this->T1_(c,k) this->T1_(d,l)  // flops: o2v4 += o0v2 o2v4 | mem: o2v4 += o0v2 o2v4

// H_vvovvo += -0.25 d(b,e) d(a,f) d(i,j) <l,k||c,d> this->T2_(c,d,l,k)  // flops: o2v4 += o0v2 o2v4 | mem: o2v4 += o0v2 o2v4

// H_vvovvo += +1.00 d(a,e) d(b,f) d(i,j) f(k,c) this->T1_(c,k)  // flops: o2v4 += o0v2 o2v4 | mem: o2v4 += o0v2 o2v4

// H_vvovvo += -1.00 d(b,e) d(a,f) d(i,j) f(k,k)  // flops: o2v4 += o0v2 o2v4 | mem: o2v4 += o0v2 o2v4

// H_vvovvo += -0.50 d(a,e) d(b,f) d(i,j) <l,k||l,k>  // flops: o2v4 += o0v2 o2v4 | mem: o2v4 += o0v2 o2v4

// H_vvovvo += -0.50 d(a,e) d(b,f) d(i,j) <l,k||c,d> this->T1_(c,k) this->T1_(d,l)  // flops: o2v4 += o0v2 o2v4 | mem: o2v4 += o0v2 o2v4

// H_vvovvo += +1.00 d(a,e) d(b,f) d(i,j) f(k,k)  // flops: o2v4 += o0v2 o2v4 | mem: o2v4 += o0v2 o2v4

// H_vvovvo += -1.00 d(b,e) d(a,f) d(i,j) f(k,c) this->T1_(c,k)  // flops: o2v4 += o0v2 o2v4 | mem: o2v4 += o0v2 o2v4

// tmps_[18_vo](d,k) = 1.00 eri[oovv](k,j,c,d) * this->T1_(c,j) // flops: o1v1 = o2v2 | mem: o1v1 = o1v1
tmps_.emplace(std::make_pair("vo_18", TAmanager.malloc<MatsT>("vo")));
tmps_["18_vo"]("d,k")  = conj(this->antiSymMoints["vvoo"]("c,d,k,j")) * this->T1_("c,j");

// tmps_[32_vv](e,a) = 1.00 this->T1_(a,j) * eri[oovv](j,i,b,e) * this->T1_(b,i) // flops: o0v2 = o1v2 | mem: o0v2 = o0v2
tmps_.emplace(std::make_pair("vv_32", TAmanager.malloc<MatsT>("vv")));
tmps_["32_vv"]("e,a")  = this->T1_("a,j") * tmps_["18_vo"]("e,j");

// H_vvovvo += -1.00 d(b,e) d(i,j) <l,k||c,f> this->T1_(c,k) this->T1_(a,l)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
H_vvovvo("a,b,e,f,i,j") -= tmps_["32_vv"]("f,a") * tmps_["17_vvoo"]("e,b,j,i");

// H_vvovvo += +1.00 d(b,f) d(i,j) <l,k||c,e> this->T1_(c,k) this->T1_(a,l)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
H_vvovvo("a,b,e,f,i,j") += tmps_["32_vv"]("e,a") * tmps_["17_vvoo"]("f,b,j,i");

// H_vvovvo += -1.00 d(a,f) d(i,j) <l,k||c,e> this->T1_(c,k) this->T1_(b,l)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
H_vvovvo("a,b,e,f,i,j") -= tmps_["32_vv"]("e,b") * tmps_["17_vvoo"]("f,a,j,i");

// H_vvovvo += +1.00 d(a,e) d(i,j) <l,k||c,f> this->T1_(c,k) this->T1_(b,l)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
H_vvovvo("a,b,e,f,i,j") += tmps_["32_vv"]("f,b") * tmps_["17_vvoo"]("e,a,j,i");

// tmps_[22_vv](a,e) = 1.00 f[ov](i,e) * this->T1_(a,i) // flops: o0v2 = o1v2 | mem: o0v2 = o0v2
tmps_.emplace(std::make_pair("vv_22", TAmanager.malloc<MatsT>("vv")));
tmps_["22_vv"]("a,e")  = this->fockMatrix_ta["ov"]("i,e") * this->T1_("a,i");

// H_vvovvo += -1.00 d(b,f) d(i,j) f(k,e) this->T1_(a,k)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
H_vvovvo("a,b,e,f,i,j") -= tmps_["22_vv"]("a,e") * tmps_["17_vvoo"]("f,b,j,i");

// H_vvovvo += -1.00 d(a,e) d(i,j) f(k,f) this->T1_(b,k)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
H_vvovvo("a,b,e,f,i,j") -= tmps_["22_vv"]("b,f") * tmps_["17_vvoo"]("e,a,j,i");

// H_vvovvo += +1.00 d(a,f) d(i,j) f(k,e) this->T1_(b,k)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
H_vvovvo("a,b,e,f,i,j") += tmps_["22_vv"]("b,e") * tmps_["17_vvoo"]("f,a,j,i");

// H_vvovvo += +1.00 d(b,e) d(i,j) f(k,f) this->T1_(a,k)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
H_vvovvo("a,b,e,f,i,j") += tmps_["22_vv"]("a,f") * tmps_["17_vvoo"]("e,b,j,i");
TAmanager.free("vvoo", std::move(tmps_["17_vvoo"]));

// H_vvvo += -1.00 d(a,e) <j,i||b,f> this->T1_(b,i)  // flops: o1v3 += o1v3 | mem: o1v3 += o1v3
H_vvvo("a,e,f,j") -= tmps_["18_vo"]("f,j") * Id_vv("a,e");

// H_vvvo += +1.00 d(a,f) <j,i||b,e> this->T1_(b,i)  // flops: o1v3 += o1v3 | mem: o1v3 += o1v3
H_vvvo("a,e,f,j") += tmps_["18_vo"]("e,j") * Id_vv("a,f");

// H_vvov += -1.00 <k,j||c,e> this->T1_(c,j) this->T2_(a,b,i,k)  // flops: o1v3 += o2v3 | mem: o1v3 += o1v3
H_vvov("a,b,e,i") -= tmps_["18_vo"]("e,k") * this->T2_("a,b,i,k");

// tmps_[31_vo](b,i) = 1.00 eri[oovv](k,j,c,d) * this->T1_(c,j) * this->T2_(d,b,i,k) // flops: o1v1 = o2v2 | mem: o1v1 = o1v1
tmps_.emplace(std::make_pair("vo_31", TAmanager.malloc<MatsT>("vo")));
tmps_["31_vo"]("b,i")  = tmps_["18_vo"]("d,k") * this->T2_("d,b,i,k");

// H_vvov += -1.00 d(b,e) <k,j||c,d> this->T1_(c,j) this->T2_(d,a,i,k)  // flops: o1v3 += o1v3 | mem: o1v3 += o1v3
H_vvov("a,b,e,i") -= tmps_["31_vo"]("a,i") * Id_vv("b,e");

// H_vvov += +1.00 d(a,e) <k,j||c,d> this->T1_(c,j) this->T2_(d,b,i,k)  // flops: o1v3 += o1v3 | mem: o1v3 += o1v3
H_vvov("a,b,e,i") += tmps_["31_vo"]("b,i") * Id_vv("a,e");
TAmanager.free("vo", std::move(tmps_["31_vo"]));

// H_vv += +1.00 <j,i||b,e> this->T1_(b,i) this->T1_(a,j)  // flops: o0v2 += o0v2 | mem: o0v2 += o0v2
H_vv("a,e") += tmps_["32_vv"]("e,a");
TAmanager.free("vv", std::move(tmps_["32_vv"]));

// tmps_[34_oo](k,i) = 1.00 this->T1_(d,i) * eri[oovv](k,j,c,d) * this->T1_(c,j) // flops: o2v0 = o2v1 | mem: o2v0 = o2v0
tmps_.emplace(std::make_pair("oo_34", TAmanager.malloc<MatsT>("oo")));
tmps_["34_oo"]("k,i")  = this->T1_("d,i") * tmps_["18_vo"]("d,k");
TAmanager.free("vo", std::move(tmps_["18_vo"]));

// tmps_[41_vo](b,i) = 1.00 this->T1_(b,k) * this->T1_(d,i) * eri[oovv](k,j,c,d) * this->T1_(c,j) // flops: o1v1 = o2v1 | mem: o1v1 = o1v1
tmps_.emplace(std::make_pair("vo_41", TAmanager.malloc<MatsT>("vo")));
tmps_["41_vo"]("b,i")  = this->T1_("b,k") * tmps_["34_oo"]("k,i");

// H_vvov += -1.00 d(b,e) <k,j||c,d> this->T1_(c,j) this->T1_(d,i) this->T1_(a,k)  // flops: o1v3 += o1v3 | mem: o1v3 += o1v3
H_vvov("a,b,e,i") -= tmps_["41_vo"]("a,i") * Id_vv("b,e");

// H_vvov += +1.00 d(a,e) <k,j||c,d> this->T1_(c,j) this->T1_(d,i) this->T1_(b,k)  // flops: o1v3 += o1v3 | mem: o1v3 += o1v3
H_vvov("a,b,e,i") += tmps_["41_vo"]("b,i") * Id_vv("a,e");
TAmanager.free("vo", std::move(tmps_["41_vo"]));

// tmps_[42_vvoo](e,a,i,k) = 1.00 Id[vv](a,e) * this->T1_(d,i) * eri[oovv](k,j,c,d) * this->T1_(c,j) // flops: o2v2 = o2v2 | mem: o2v2 = o2v2
tmps_.emplace(std::make_pair("vvoo_42", TAmanager.malloc<MatsT>("vvoo")));
tmps_["42_vvoo"]("e,a,i,k")  = Id_vv("a,e") * tmps_["34_oo"]("k,i");
TAmanager.free("oo", std::move(tmps_["34_oo"]));

// H_vvovvo += -1.00 d(b,e) d(a,f) <j,k||c,d> this->T1_(c,k) this->T1_(d,i)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
H_vvovvo("a,b,e,f,i,j") -= tmps_["42_vvoo"]("e,b,i,j") * Id_vv("a,f");

// H_vvovvo += +1.00 d(a,e) d(b,f) <j,k||c,d> this->T1_(c,k) this->T1_(d,i)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
H_vvovvo("a,b,e,f,i,j") += tmps_["42_vvoo"]("e,a,i,j") * Id_vv("b,f");
TAmanager.free("vvoo", std::move(tmps_["42_vvoo"]));

// tmps_[19_vo](b,i) = 1.00 f[ov](j,c) * this->T2_(c,b,i,j) // flops: o1v1 = o2v2 | mem: o1v1 = o1v1
tmps_.emplace(std::make_pair("vo_19", TAmanager.malloc<MatsT>("vo")));
tmps_["19_vo"]("b,i")  = this->fockMatrix_ta["ov"]("j,c") * this->T2_("c,b,i,j");

// H_vvov += +1.00 d(b,e) f(j,c) this->T2_(c,a,i,j)  // flops: o1v3 += o1v3 | mem: o1v3 += o1v3
H_vvov("a,b,e,i") += tmps_["19_vo"]("a,i") * Id_vv("b,e");

// H_vvov += -1.00 d(a,e) f(j,c) this->T2_(c,b,i,j)  // flops: o1v3 += o1v3 | mem: o1v3 += o1v3
H_vvov("a,b,e,i") -= tmps_["19_vo"]("b,i") * Id_vv("a,e");
TAmanager.free("vo", std::move(tmps_["19_vo"]));

// tmps_[20_vo](a,i) = 1.00 eri[vovo](a,j,c,i) * this->T1_(c,j) // flops: o1v1 = o2v2 | mem: o1v1 = o1v1
tmps_.emplace(std::make_pair("vo_20", TAmanager.malloc<MatsT>("vo")));
tmps_["20_vo"]("a,i")  = this->antiSymMoints["vovo"]("a,j,c,i") * this->T1_("c,j");

// H_vvov += +1.00 d(a,e) <j,b||c,i> this->T1_(c,j)  // flops: o1v3 += o1v3 | mem: o1v3 += o1v3
H_vvov("a,b,e,i") -= tmps_["20_vo"]("b,i") * Id_vv("a,e");

// H_vvov += -1.00 d(b,e) <j,a||c,i> this->T1_(c,j)  // flops: o1v3 += o1v3 | mem: o1v3 += o1v3
H_vvov("a,b,e,i") += tmps_["20_vo"]("a,i") * Id_vv("b,e");
TAmanager.free("vo", std::move(tmps_["20_vo"]));

// tmps_[21_oo](i,k) = 1.00 eri[oovo](k,j,c,i) * this->T1_(c,j) // flops: o2v0 = o3v1 | mem: o2v0 = o2v0
tmps_.emplace(std::make_pair("oo_21", TAmanager.malloc<MatsT>("oo")));
tmps_["21_oo"]("i,k")  = conj(this->antiSymMoints["vooo"]("c,i,k,j")) * this->T1_("c,j");

// tmps_[35_vo](b,i) = 1.00 this->T1_(b,k) * eri[oovo](k,j,c,i) * this->T1_(c,j) // flops: o1v1 = o2v1 | mem: o1v1 = o1v1
tmps_.emplace(std::make_pair("vo_35", TAmanager.malloc<MatsT>("vo")));
tmps_["35_vo"]("b,i")  = this->T1_("b,k") * tmps_["21_oo"]("i,k");

// H_vvov += +1.00 d(a,e) <k,j||c,i> this->T1_(c,j) this->T1_(b,k)  // flops: o1v3 += o1v3 | mem: o1v3 += o1v3
H_vvov("a,b,e,i") += tmps_["35_vo"]("b,i") * Id_vv("a,e");

// H_vvov += -1.00 d(b,e) <k,j||c,i> this->T1_(c,j) this->T1_(a,k)  // flops: o1v3 += o1v3 | mem: o1v3 += o1v3
H_vvov("a,b,e,i") -= tmps_["35_vo"]("a,i") * Id_vv("b,e");
TAmanager.free("vo", std::move(tmps_["35_vo"]));

// tmps_[38_vvoo](e,a,j,i) = 1.00 Id[vv](a,e) * eri[oovo](j,k,c,i) * this->T1_(c,k) // flops: o2v2 = o2v2 | mem: o2v2 = o2v2
tmps_.emplace(std::make_pair("vvoo_38", TAmanager.malloc<MatsT>("vvoo")));
tmps_["38_vvoo"]("e,a,j,i")  = Id_vv("a,e") * tmps_["21_oo"]("i,j");
TAmanager.free("oo", std::move(tmps_["21_oo"]));

// H_vvovvo += -1.00 d(b,e) d(a,f) <j,k||c,i> this->T1_(c,k)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
H_vvovvo("a,b,e,f,i,j") -= tmps_["38_vvoo"]("e,b,j,i") * Id_vv("a,f");

// H_vvovvo += +1.00 d(a,e) d(b,f) <j,k||c,i> this->T1_(c,k)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
H_vvovvo("a,b,e,f,i,j") += tmps_["38_vvoo"]("e,a,j,i") * Id_vv("b,f");
TAmanager.free("vvoo", std::move(tmps_["38_vvoo"]));

// H_vv += -1.00 f(i,e) this->T1_(a,i)  // flops: o0v2 += o0v2 | mem: o0v2 += o0v2
H_vv("a,e") -= tmps_["22_vv"]("a,e");
TAmanager.free("vv", std::move(tmps_["22_vv"]));

// tmps_[23_vo](a,i) = 1.00 f[vv](a,c) * this->T1_(c,i) // flops: o1v1 = o1v2 | mem: o1v1 = o1v1
tmps_.emplace(std::make_pair("vo_23", TAmanager.malloc<MatsT>("vo")));
tmps_["23_vo"]("a,i")  = this->fockMatrix_ta["vv"]("a,c") * this->T1_("c,i");

// H_vvov += -1.00 d(b,e) f(a,c) this->T1_(c,i)  // flops: o1v3 += o1v3 | mem: o1v3 += o1v3
H_vvov("a,b,e,i") -= tmps_["23_vo"]("a,i") * Id_vv("b,e");

// H_vvov += +1.00 d(a,e) f(b,c) this->T1_(c,i)  // flops: o1v3 += o1v3 | mem: o1v3 += o1v3
H_vvov("a,b,e,i") += tmps_["23_vo"]("b,i") * Id_vv("a,e");
TAmanager.free("vo", std::move(tmps_["23_vo"]));

// tmps_[24_oo](i,j) = 1.00 f[ov](j,c) * this->T1_(c,i) // flops: o2v0 = o2v1 | mem: o2v0 = o2v0
tmps_.emplace(std::make_pair("oo_24", TAmanager.malloc<MatsT>("oo")));
tmps_["24_oo"]("i,j")  = this->fockMatrix_ta["ov"]("j,c") * this->T1_("c,i");

// tmps_[36_vo](b,i) = 1.00 this->T1_(b,j) * f[ov](j,c) * this->T1_(c,i) // flops: o1v1 = o2v1 | mem: o1v1 = o1v1
tmps_.emplace(std::make_pair("vo_36", TAmanager.malloc<MatsT>("vo")));
tmps_["36_vo"]("b,i")  = this->T1_("b,j") * tmps_["24_oo"]("i,j");

// H_vvov += -1.00 d(a,e) f(j,c) this->T1_(c,i) this->T1_(b,j)  // flops: o1v3 += o1v3 | mem: o1v3 += o1v3
H_vvov("a,b,e,i") -= tmps_["36_vo"]("b,i") * Id_vv("a,e");

// H_vvov += +1.00 d(b,e) f(j,c) this->T1_(c,i) this->T1_(a,j)  // flops: o1v3 += o1v3 | mem: o1v3 += o1v3
H_vvov("a,b,e,i") += tmps_["36_vo"]("a,i") * Id_vv("b,e");
TAmanager.free("vo", std::move(tmps_["36_vo"]));

// tmps_[39_vvoo](e,a,j,i) = 1.00 Id[vv](a,e) * f[ov](j,c) * this->T1_(c,i) // flops: o2v2 = o2v2 | mem: o2v2 = o2v2
tmps_.emplace(std::make_pair("vvoo_39", TAmanager.malloc<MatsT>("vvoo")));
tmps_["39_vvoo"]("e,a,j,i")  = Id_vv("a,e") * tmps_["24_oo"]("i,j");
TAmanager.free("oo", std::move(tmps_["24_oo"]));

// H_vvovvo += -1.00 d(a,e) d(b,f) f(j,c) this->T1_(c,i)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
H_vvovvo("a,b,e,f,i,j") -= tmps_["39_vvoo"]("e,a,j,i") * Id_vv("b,f");

// H_vvovvo += +1.00 d(b,e) d(a,f) f(j,c) this->T1_(c,i)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
H_vvovvo("a,b,e,f,i,j") += tmps_["39_vvoo"]("e,b,j,i") * Id_vv("a,f");
TAmanager.free("vvoo", std::move(tmps_["39_vvoo"]));

// tmps_[25_vo](b,i) = 1.00 f[oo](j,i) * this->T1_(b,j) // flops: o1v1 = o2v1 | mem: o1v1 = o1v1
tmps_.emplace(std::make_pair("vo_25", TAmanager.malloc<MatsT>("vo")));
tmps_["25_vo"]("b,i")  = this->fockMatrix_ta["oo"]("j,i") * this->T1_("b,j");

// H_vvov += -1.00 d(a,e) f(j,i) this->T1_(b,j)  // flops: o1v3 += o1v3 | mem: o1v3 += o1v3
H_vvov("a,b,e,i") -= tmps_["25_vo"]("b,i") * Id_vv("a,e");

// H_vvov += +1.00 d(b,e) f(j,i) this->T1_(a,j)  // flops: o1v3 += o1v3 | mem: o1v3 += o1v3
H_vvov("a,b,e,i") += tmps_["25_vo"]("a,i") * Id_vv("b,e");
TAmanager.free("vo", std::move(tmps_["25_vo"]));

// tmps_[26_vvoo](e,a,i,j) = 1.00 Id[vv](a,e) * f[oo](j,i) // flops: o2v2 = o2v2 | mem: o2v2 = o2v2
tmps_.emplace(std::make_pair("vvoo_26", TAmanager.malloc<MatsT>("vvoo")));
tmps_["26_vvoo"]("e,a,i,j")  = Id_vv("a,e") * this->fockMatrix_ta["oo"]("j,i");

// H_vvovvo += +1.00 d(b,e) d(a,f) f(j,i)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
H_vvovvo("a,b,e,f,i,j") += tmps_["26_vvoo"]("e,b,i,j") * Id_vv("a,f");

// H_vvovvo += -1.00 d(a,e) d(b,f) f(j,i)  // flops: o2v4 += o2v4 | mem: o2v4 += o2v4
H_vvovvo("a,b,e,f,i,j") -= tmps_["26_vvoo"]("e,a,i,j") * Id_vv("b,f");
TAmanager.free("vvoo", std::move(tmps_["26_vvoo"]));

    TA::get_default_world().gop.fence();







  }
}
