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

#include <particleintegrals/twopints/incore4indexreleri.hpp>
#include <particleintegrals/twopints/incore4indextpi.hpp>
#include <particleintegrals/twopints/incoreritpi.hpp>
#include <particleintegrals/onepints/relativisticints.hpp>
#include <matrix.hpp>
#include <cqlinalg.hpp>
#include <physcon.hpp>

namespace ChronusQ {

  // InCore4indexTPI::spatialToSpinBlock helper macros
  #define PAULITOSPINOR_I(SPIN_M, SPIN_N, SCALE_) \
    auto SPIN_N = SPIN_M; \
    double SCALE_ = 1.;

  #define PAULITOSPINOR_X(SPIN_M, SPIN_N, SCALE_) \
    auto SPIN_N = 1 - SPIN_M; \
    double SCALE_ = 1.;
     
  #define PAULITOSPINOR_Y(SPIN_M, SPIN_N, SCALE_) \
    auto SPIN_N = 1 - SPIN_M; \
    dcomplex SCALE_; \
    if (SPIN_M == 0) SCALE_ = dcomplex(0.,-1.); \
    else SCALE_ = dcomplex(0.,1.);

  #define PAULITOSPINOR_Z(SPIN_M, SPIN_N, SCALE_) \
    auto SPIN_N = SPIN_M; \
    double SCALE_; \
    if (SPIN_M == 0) SCALE_ = 1.; \
    else SCALE_ = -1.;

  #define PAULITOSPINOR_4INDEXTPI(OPR1, OPR2) \
    for (auto splm = 0ul; splm < 2; splm++) \
    for (auto spmu = 0ul; spmu < 2; spmu++) { \
      PAULITOSPINOR_ ## OPR1(spmu, spnu, scale1) \
      PAULITOSPINOR_ ## OPR2(splm, spsg, scale2) \
      auto scale = scale1 * scale2; \
      for (auto sg = 0ul; sg < NB; sg++) \
      for (auto lm = 0ul; lm < NB; lm++) \
      for (auto nu = 0ul; nu < NB; nu++) \
      for (auto mu = 0ul; mu < NB; mu++)  \
        spinor(spmu*NB+mu, spnu*NB+nu, splm*NB+lm, spsg*NB+sg) = scale * (*this)(mu, nu, lm, sg); \
    } 
  
  #define TRY_PAULITOSPINOR_4INDEXTPI(OP1, OP2) \
    if( TRANS1 == #OP1[0] and TRANS2 == #OP2[0]) \
      PAULITOSPINOR_4INDEXTPI(OP1,OP2)   
  
  template <typename IntsT>
  template <typename IntsU>
  InCore4indexTPI<IntsU> InCore4indexTPI<IntsT>::spatialToSpinBlock(
      char TRANS1, char TRANS2) const {
    
    size_t NB = this->nBasis();
    InCore4indexTPI<IntsU> spinor(2*NB);
    spinor.clear();

    //Try all the combinations
    TRY_PAULITOSPINOR_4INDEXTPI(I,I)
    else TRY_PAULITOSPINOR_4INDEXTPI(I,X)
    else TRY_PAULITOSPINOR_4INDEXTPI(I,Y)
    else TRY_PAULITOSPINOR_4INDEXTPI(I,Z)
    else TRY_PAULITOSPINOR_4INDEXTPI(X,I)
    else TRY_PAULITOSPINOR_4INDEXTPI(X,X)
    else TRY_PAULITOSPINOR_4INDEXTPI(X,Y)
    else TRY_PAULITOSPINOR_4INDEXTPI(X,Z)
    else TRY_PAULITOSPINOR_4INDEXTPI(Y,I)
    else TRY_PAULITOSPINOR_4INDEXTPI(Y,X)
    else TRY_PAULITOSPINOR_4INDEXTPI(Y,Y)
    else TRY_PAULITOSPINOR_4INDEXTPI(Y,Z)
    else TRY_PAULITOSPINOR_4INDEXTPI(Z,I)
    else TRY_PAULITOSPINOR_4INDEXTPI(Z,X)
    else TRY_PAULITOSPINOR_4INDEXTPI(Z,Y)
    else TRY_PAULITOSPINOR_4INDEXTPI(Z,Z)
    else CErr("Not a valid pair of TRANS1 and TRANS2");

    return spinor;
  } // InCore4indexTPI::spatialToSpinBlock

  template<>
  template<>
  InCore4indexTPI<double> InCore4indexTPI<double>::spatialToSpinBlock(
      char TRANS1, char TRANS2) const {
    
    size_t NB = this->nBasis();
    InCore4indexTPI<double> spinor(2*NB);
    spinor.clear();
    
    //Try all the combinations with out Y
    TRY_PAULITOSPINOR_4INDEXTPI(I,I)
    else TRY_PAULITOSPINOR_4INDEXTPI(I,X)
    else TRY_PAULITOSPINOR_4INDEXTPI(I,Z)
    else TRY_PAULITOSPINOR_4INDEXTPI(X,I)
    else TRY_PAULITOSPINOR_4INDEXTPI(X,X)
    else TRY_PAULITOSPINOR_4INDEXTPI(X,Z)
    else TRY_PAULITOSPINOR_4INDEXTPI(Z,I)
    else TRY_PAULITOSPINOR_4INDEXTPI(Z,X)
    else TRY_PAULITOSPINOR_4INDEXTPI(Z,Z)
    else CErr("Not a valid pair of TRANS1 and TRANS2");

    return spinor;
  
  } // InCore4indexTPI::spatialToSpinBlock(MatsT=double, IntsT=double)
  
  // IntsU can only be complex here
  template <typename IntsT>
  template <typename IntsU>
  InCore4indexRelERI<IntsU> InCore4indexRelERI<IntsT>::spatialToSpinBlock() const {
    
    size_t nSpinorRelComp;
    size_t NB    = this->nBasis();
    size_t twoNB  = NB*2;
    size_t twoNB2 = twoNB*twoNB;
    dcomplex scale = 1./(4*SpeedOfLight*SpeedOfLight);
    dcomplex iscale = dcomplex(0.0, 1./(4*SpeedOfLight*SpeedOfLight));

    if (this->nRelComp() == 0)       nSpinorRelComp = 0; // direct Coulomb  (LL|LL)
    else if (this->nRelComp() == 4)  nSpinorRelComp = 1; // + Dirac Coulomb (LL|SS) 
    else if (this->nRelComp() == 23) nSpinorRelComp = 2; // + Gaunt (LS |dot SL) ??
    else CErr("Unrecognizable nRelComponent");
  
    InCore4indexRelERI<IntsU> spinor(twoNB, nSpinorRelComp);    
    
    // LLLL part
    {
      auto tmp = InCore4indexTPI<IntsT>::template spatialToSpinBlock<IntsU>();
      SetMat('N', twoNB2, twoNB2, dcomplex(1.), tmp.pointer(), twoNB2, 
        spinor.pointer(), twoNB2);
    } 

    #define ADD_COMPONENTS_TO_SPINOR_RELERI(TRANS1, TRANS2, SCALE, COMP, SPINOR) \
    { auto tmp = this->components_[COMP].template spatialToSpinBlock<IntsU>(TRANS1, TRANS2);\
      MatAdd('N', 'N', twoNB2, twoNB2, dcomplex(1.), SPINOR.pointer(), twoNB2, \
        SCALE, tmp.pointer(), twoNB2, SPINOR.pointer(), twoNB2); }
    
    // SSLL part
    if (nSpinorRelComp > 0) {
      
      auto & SSLL = spinor.components_[0];
    
      SSLL.clear();
      // SSLL += I(1)I(2) ∇A∙∇B(ij|kl)
      ADD_COMPONENTS_TO_SPINOR_RELERI('I', 'I', scale, 0, SSLL);
      
      // SSLL += I(1)x(2) [∇Ax∇B(ij|kl)]_x 
      ADD_COMPONENTS_TO_SPINOR_RELERI('X', 'I', iscale, 1, SSLL);
      // SSLL += I(1)y(2) [∇Ax∇B(ij|kl)]_y 
      ADD_COMPONENTS_TO_SPINOR_RELERI('Y', 'I', iscale, 2, SSLL);
      // SSLL += I(1)z(2) [∇Ax∇B(ij|kl)]_z 
      ADD_COMPONENTS_TO_SPINOR_RELERI('Z', 'I', iscale, 3, SSLL);
      //SSLL.output(std::cout, "DC SSLL AO Integrals", true);
    } 
  
    if (nSpinorRelComp > 1) { 
      CErr("Guant term spatial to spinor transformation is not implemented");
    }
  
    return spinor;
  
  } // InCore4indexRelERI::spatialToSpinBlock
  
  template <typename IntsT>
  template <typename IntsU>
  InCoreRITPI<IntsU> InCoreRITPI<IntsT>::spatialToSpinBlock() const {
    size_t NB = this->nBasis();
    InCoreRITPI<IntsU> spinor(2*NB, NBRI);
/*
    for ( auto sp = 0ul; sp < 2; sp++)
    for ( auto nu = 0ul; nu < NB; nu++)
    for ( auto mu = 0ul; mu < NB; mu++)
    for ( auto L = 0ul; L < nRIBasis(); L++) {
      spinor(L, sp*NB + mu, sp*NB + nu) = (*this)(L, mu, nu);
    }
*/
    SetMatDiag(NB*NBRI, NB, pointer(), NB*NBRI, spinor.pointer(), 2*NB*NBRI);
    return spinor;
  }

  template <typename IntsT>
  template <typename MatsT>
  cqmatrix::Matrix<MatsT> OnePRelInts<IntsT>::formW() const {
    if (hasSpinOrbit()) {
      size_t NB = this->nBasis();
      cqmatrix::Matrix<MatsT> W(2*NB);
      // W = [ W1  W2 ]
      //     [ W3  W4 ]
      dcomplex *W1 = W.pointer();
      dcomplex *W2 = W1 + 2*NB*NB;
      dcomplex *W3 = W1 + NB;
      dcomplex *W4 = W2 + NB;
      // W1 = pV.p + i (pVxp)(Z)
      MatAdd('N','N',NB,NB,dcomplex(1.),scalar().pointer(),NB,
             dcomplex(0.,1.),SOZ().pointer(),NB,W1,2*NB);
      // W4 = pV.p - i (pVxp)(Z)
      MatAdd('N','N',NB,NB,dcomplex(1.),scalar().pointer(),NB,
             dcomplex(0.,-1.),SOZ().pointer(),NB,W4,2*NB);
      // W2 = (pVxp)(Y) + i (pVxp)(X)
      MatAdd('N','N',NB,NB,dcomplex(1.),SOY().pointer(),NB,
             dcomplex(0.,1.),SOX().pointer(),NB,W2,2*NB);
      // W3 = -(pVxp)(Y) + i (pVxp)(X)
      MatAdd('N','N',NB,NB,dcomplex(-1.),SOY().pointer(),NB,
             dcomplex(0.,1.),SOX().pointer(),NB,W3,2*NB);
      return W;
    } else
      return scalar().matrix().template spatialToSpinBlock<MatsT>();
  }

  template <>
  template <>
  cqmatrix::Matrix<double> OnePRelInts<double>::formW() const {
    if (hasSpinOrbit()) {
      cqmatrix::Matrix<double> dummy(1);
      CErr("W matrix with spin-orbit cannot be real.");
      return dummy;
    } else {
      return scalar().matrix().template spatialToSpinBlock<double>();
    }
  }
 
  template InCore4indexTPI<dcomplex> InCore4indexTPI<double>::spatialToSpinBlock(char, char) const;
  template InCore4indexTPI<dcomplex> InCore4indexTPI<dcomplex>::spatialToSpinBlock(char, char) const;


  template <>
  template <>
  InCore4indexRelERI<double>  InCore4indexRelERI<double>::spatialToSpinBlock() const { 
    CErr("It's not valid to have double for relativitic 2C spinor");
    abort();
  };

  template InCore4indexRelERI<dcomplex> InCore4indexRelERI<double>::spatialToSpinBlock() const;
  template InCore4indexRelERI<dcomplex> InCore4indexRelERI<dcomplex>::spatialToSpinBlock() const;
  
  template InCoreRITPI<double> InCoreRITPI<double>::spatialToSpinBlock() const;
  template InCoreRITPI<dcomplex> InCoreRITPI<double>::spatialToSpinBlock() const;
  template InCoreRITPI<dcomplex> InCoreRITPI<dcomplex>::spatialToSpinBlock() const;

  template cqmatrix::Matrix<dcomplex> OnePRelInts<double>::formW() const;
  template cqmatrix::Matrix<dcomplex> OnePRelInts<dcomplex>::formW() const;

}; // namespace ChronusQ
