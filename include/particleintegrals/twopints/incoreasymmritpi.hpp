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

#include <particleintegrals/twopints.hpp>
#include <particleintegrals/twopints/incore4indextpi.hpp>
#include <particleintegrals/twopints/incoreritpi.hpp>
#include <cqlinalg/blas1.hpp>
#include <cqlinalg/blas3.hpp>
#include <cqlinalg/blasutil.hpp>
#include <cxxapi/output.hpp>

namespace ChronusQ {

  enum class NEO_CD_ALG {
    //Note: For the following 3 _AUX options, if there's no corresponding aux basis available, then will build them on the fly
    ELEC_AUX,              // Will use aux basis from (ee|ee) to approximate (ee|pp)
    PROT_AUX,              // Will use aux basis from (pp|pp) to approximate (ee|pp)
    ELEC_AND_PROT_AUX,     // Will use both aux basis from (ee|ee) and (pp|pp) to approximate (ee|pp)

    AUTO,                  // Default option. Will try to use aux basis by detecting what's available. If can't find one, then default to 4-index 
  };

/*
class InCoreAsymmRITPI :
 This class allows approximation of an asymmetric ERI matrix, by using one or two existing
 auxliary basis.

 For example, for NEO, (ee|pp) integral can be approxiamated by substituting in either 
 electronic or protonic auxiliary basis. Both auxiliary basis can also be used together.   
*/

  template <typename IntsT>
  class InCoreAsymmRITPI : public TwoPInts<IntsT> {

    template <typename IntsU>
    friend class InCoreAsymmRITPI;

  protected:
    std::shared_ptr<InCoreRITPI<IntsT>> aux1_ = nullptr;
    std::shared_ptr<InCoreRITPI<IntsT>> aux2_ = nullptr;
    IntsT *partialTPI_ = nullptr; /// The missing part that needs to be built, to be used together with existing 3-index tensor 

    bool build4I_ = false; // Explicitly build four-index TPI
    bool printError_ = false; // Build full precision 4-index TPI, and build approximate 4-index TPI, then report error
    std::shared_ptr<InCore4indexTPI<IntsT>> eri4I_ = nullptr; // Four-index TPI

  public:

    // CONSTRUCTORS:
    // Disable defualt constructor
    // InCoreAsymmRITPI() = delete;
    // Constructor for when one aux basis is given
    InCoreAsymmRITPI(CQMemManager &mem, std::shared_ptr<InCoreRITPI<IntsT>> aux1, size_t sNB, bool build4I = false):
        TwoPInts<IntsT>(mem, aux1->nBasis(), sNB), aux1_(aux1), build4I_(build4I){}
    InCoreAsymmRITPI(CQMemManager &mem, size_t NB, std::shared_ptr<InCoreRITPI<IntsT>> aux2, bool build4I = false):
        TwoPInts<IntsT>(mem, NB, aux2->nBasis()), aux2_(aux2), build4I_(build4I){}
    // Constructor for when two aux basis are given
    InCoreAsymmRITPI(CQMemManager &mem, std::shared_ptr<InCoreRITPI<IntsT>> aux1, std::shared_ptr<InCoreRITPI<IntsT>> aux2, bool build4I = false):
        TwoPInts<IntsT>(mem, aux1->nBasis(), aux2->nBasis()), aux1_(aux1), aux2_(aux2), build4I_(build4I){}


    // COPY CONSTRUCTOR:
    InCoreAsymmRITPI(const InCoreAsymmRITPI &other): TwoPInts<IntsT>(other),
      aux1_(other.getAux1()), aux2_(other.getAux2()), build4I_(other.build4I_), eri4I_(other.eri4I_) {
      malloc();
      std::copy_n(other.partialTPI_, getPartialTPISize(), this->partialTPI_);
    }
    template <typename IntsU>
    InCoreAsymmRITPI( const InCoreAsymmRITPI<IntsU> &other, int = 0 ): TwoPInts<IntsT>(other)
    {
      if (std::is_same<IntsU, dcomplex>::value
          and std::is_same<IntsT, double>::value)
        CErr("Cannot create a Real InCoreRITPI from a Complex one.");
      CErr("Type conversion NYI in InCoreAsymmRITPI.");
    }

    // DESTRUCTOR:
    virtual ~InCoreAsymmRITPI() {
      if(partialTPI_) this->memManager().free(partialTPI_);
    }

    // MOVE CONSTRUCTOR:
    InCoreAsymmRITPI( InCoreAsymmRITPI &&other ): TwoPInts<IntsT>(std::move(other)),
        aux1_(other.getAux1()), aux2_(other.getAux2()), partialTPI_(other.partialTPI_),build4I_(other.build4I_),
        eri4I_(other.eri4I_) {
      other.partialTPI_ = nullptr;
    }

    // ASSIGNMENT OPERATOR
    InCoreAsymmRITPI& operator=( const InCoreAsymmRITPI &other ) {
      if (this != &other) { // self-assignment check expected
        TwoPInts<IntsT>::operator=(other);// if can't pass, manually 
        this->aux1_ = other.getAux1();
        this->aux2_ = other.getAux2();
        build4I_ = other.build4I_;
        this->eri4I_ = other.eri4I_;
        malloc(); // reallocate memory
        std::copy_n(other.partialTPI_, getPartialTPISize(), this->partialTPI_);
      }
      return *this;
    }

    // MOVE ASSIGNMENT OPERATOR
    InCoreAsymmRITPI& operator=( InCoreAsymmRITPI &&other ) {
      if (this != &other) { // self-assignment check expected
        TwoPInts<IntsT>::operator=(std::move(other));
        this->memManager().free(partialTPI_);
        this->aux1_ = other.getAux1();
        this->aux2_ = other.getAux2();
        build4I_ = other.build4I_;
        this->eri4I_ = other.eri4I_;
        this->partialTPI_ = other.partialTPI_;
        other.partialTPI_ = nullptr;
      }
      return *this;
    }

    // Computation interfaces
    virtual void computeAOInts(BasisSet&, Molecule&, EMPerturbation&,
        OPERATOR, const HamiltonianOptions&) override{
      CErr("AO integral evaluation is NOT implemented in class InCoreAsymmRITPI. This class only handles integral with existing aux basis");
    }

    virtual void computeAOInts(BasisSet&, BasisSet&, Molecule&, EMPerturbation&,
        OPERATOR, const HamiltonianOptions&) override;


    void computeOneCholeskyPartialTPILibint(BasisSet&, BasisSet&);

    void computeTwoCholeskyPartialTPILibint(BasisSet&, BasisSet&);

    void computeOneCholeskyPartialTPIPrebuilt4Index(BasisSet&, BasisSet&);

    void computeTwoCholeskyPartialTPIPrebuilt4Index(BasisSet&, BasisSet&);

    InCore4indexTPI<IntsT> to4indexERI() {

      if(not partialTPI_) CErr ("Missing PartialTPI. Can't convert to 4-index");

      size_t NB1 = this->nBasis();
      size_t NB2 = this->snBasis();
      size_t NB1_Squared = NB1 * NB1;
      size_t NB2_Squared = NB2 * NB2;

      InCore4indexTPI<IntsT> eri4i(this->memManager(), NB1, NB2);
      
      if(aux1_){
        size_t NBRI1 = aux1_->nRIBasis();
        if(aux2_){
          size_t NBRI2 = aux2_->nRIBasis();
          double *SCR = this->memManager().template malloc<IntsT>(NBRI1 * NB2_Squared);
          blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans, NBRI1, NB2_Squared, NBRI2, IntsT(1.), partialTPI_, NBRI1, aux2_->pointer(), NBRI2, IntsT(0.), SCR, NBRI1);
          blas::gemm(blas::Layout::ColMajor,blas::Op::Trans,blas::Op::NoTrans, NB1_Squared, NB2_Squared, NBRI1, IntsT(1.), aux1_->pointer(), NBRI1, SCR, NBRI1, IntsT(0.), eri4i.pointer(), NB1_Squared);
          this->memManager().free(SCR);
        } else{
          blas::gemm(blas::Layout::ColMajor,blas::Op::Trans,blas::Op::NoTrans,NB1_Squared, NB2_Squared, NBRI1, IntsT(1.), aux1_->pointer(), NBRI1, partialTPI_, NBRI1, IntsT(0.), eri4i.pointer(), NB1_Squared);
        }
      } else {
        if(aux2_){
          size_t NBRI2 = aux2_->nRIBasis();
          blas::gemm(blas::Layout::ColMajor,blas::Op::Trans,blas::Op::NoTrans,NB1_Squared, NB2_Squared, NBRI2, IntsT(1.), partialTPI_, NBRI2, aux2_->pointer(), NBRI2, IntsT(0.), eri4i.pointer(), NB1_Squared);
        } else{
          CErr ("No aux basis found. Can't convert to 4-index");
        } 
      }
      return eri4i;
    }
    
    // Determine the size of partialERI
    size_t getPartialTPISize(){
      size_t sizePartialTPI = 0, NB1 = this->nBasis(), NB2 = this->snBasis();
      if (aux1_){
        sizePartialTPI = aux2_ ? aux1_->nRIBasis()*aux2_->nRIBasis() : aux1_->nRIBasis()*NB2*NB2;
      } else{
        sizePartialTPI = aux2_ ? aux2_->nRIBasis()*NB1*NB1 : 0;
      }
      return sizePartialTPI;
    }

    void malloc() {

      size_t sizePartialTPI = getPartialTPISize();
      if (sizePartialTPI == 0) CErr("No Auxiliary basis available to form InCoreAsymmRITPI object");
      
      // delete old partialTPI_ if size is wrong
      if (partialTPI_) {
        if (this->memManager().getSize(partialTPI_) == sizePartialTPI)
          return;
        this->memManager().free(partialTPI_);
      }
      
      try { partialTPI_ = this->memManager().template malloc<IntsT>(sizePartialTPI); }
      catch(...) {
        std::cout << std::fixed;
        std::cout << "Insufficient memory for the full RI-ERI tensor ("
                  << (sizePartialTPI/1e9) * sizeof(double) << " GB)" << std::endl;
        std::cout << std::endl << this->memManager() << std::endl;
        CErr();
      }
    }

    // Single element interfaces
    virtual IntsT operator()(size_t, size_t, size_t, size_t) const override{
      CErr("Single Element Indexing NYI");
    };
    virtual IntsT operator()(size_t, size_t) const override{
      CErr("Single Element Indexing NYI");
    }

    virtual void output(std::ostream& out, const std::string& s= "", 
        bool printFull = false) const override{
      if (s == "")
        out << "  Asymmetric Two Particle Integral:" << std::endl;
      else
        out << "  ERI[" << s << "]:" << std::endl;
      out << "    * Contraction Algorithm: ";
      out << "INCORE Cholesky decomposition (Using Aux Basis)";
      out << std::endl;
      if (printFull) {
        CErr("Print Full for AsymmRITPI NYI");
        out << bannerEnd << std::endl;
      }
    }


    virtual void clear() override{
      std::fill_n(partialTPI_, getPartialTPISize(), IntsT(0.));
    }

    // pointers direct access
    std::shared_ptr<InCoreRITPI<IntsT>> getAux1() const {return aux1_;}
    std::shared_ptr<InCoreRITPI<IntsT>> getAux2() const {return aux2_;}
    
    void setPartialTPI(const IntsT* partialTPI){
      std::copy_n(partialTPI, getPartialTPISize() ,partialTPI_);
    }

    IntsT* pointer() { return partialTPI_; }
    const IntsT* pointer() const { return partialTPI_; }

    void setPrintError(bool printError) { printError_ = printError;}

    void printError(BasisSet &basisSet1, BasisSet &basisSet2, Molecule& mol, EMPerturbation& emPert){

      std::cout << bannerTop <<  std::endl;
      std::cout << "\nCalculating Error Associated with approximate (ee|pp) Integrals: \n" << std::endl;

      size_t NB1 = basisSet1.nBasis;
      size_t NB2 = basisSet2.nBasis;
      
      if(aux1_){
        std::cout<< "     * Error from the approximation of (ee|ee) Integrals" << std::endl;
        size_t lenTPI = NB1*NB1*NB1*NB1;
        calculateDifferece(basisSet1, basisSet1, mol, emPert, {-1., 1.}, ELECTRON_REPULSION, 
            aux1_->to4indexERI().pointer(), lenTPI);
      }
      
      if(aux2_){
        std::cout<< "     * Error from the approximation of (pp|pp) Integrals" << std::endl;
        size_t lenTPI = NB2*NB2*NB2*NB2;
        calculateDifferece(basisSet2, basisSet2, mol, emPert, {1., ProtMassPerE}, ELECTRON_REPULSION, 
            aux2_->to4indexERI().pointer(), lenTPI);
      }
      
      std::cout<< "     * Error from the approximation of (ee|pp) Integrals" << std::endl;
      size_t lenTPI = NB1*NB1*NB2*NB2;
      calculateDifferece(basisSet1, basisSet2, mol, emPert, {1., ProtMassPerE}, EP_ATTRACTION, 
          to4indexERI().pointer(), lenTPI);
 
      std::cout << bannerEnd <<  std::endl;

    }

    void calculateDifferece(BasisSet &basisSet1, BasisSet &basisSet2,  Molecule& mol, EMPerturbation& emPert, 
        Particle P, OPERATOR op, IntsT* approxTPI, size_t lenTPI){
      // Build Exact TPI:
      HamiltonianOptions temp_opt;
      temp_opt.particle = P;
      InCore4indexTPI<IntsT> exactTPI(this->memManager(),basisSet1.nBasis,basisSet2.nBasis);
      exactTPI.computeAOInts(basisSet1, basisSet2, mol, emPert, op, temp_opt);

      // Get the difference between approximate TPI and exact TPI
      IntsT* diff = this->memManager().template malloc<IntsT>(lenTPI);
      std::copy_n(approxTPI,lenTPI,diff);
      blas::axpy(lenTPI, -1.0, exactTPI.pointer(), 1, diff, 1);

      auto minmaxApprox = std::minmax_element(diff, diff+lenTPI);
      IntsT maxElement = std::max(-*minmaxApprox.first, *minmaxApprox.second);
      IntsT norm = blas::nrm2(lenTPI, diff, 1);

      std::cout << "       Max Element: " << maxElement << std::endl; 
      std::cout << "       Norm: " << norm << std::endl; 
      std::cout << std::endl; 

      this->memManager().free(diff);
    }

    
  }; // class InCoreAsymmRITPI

  template <typename MatsT, typename IntsT>
  class InCoreAsymmRITPIContraction : public InCore4indexTPIContraction<MatsT,IntsT> {

    template <typename MatsU, typename IntsU>
    friend class InCoreAsymmRITPIContraction;

  public:

    // Constructors

    InCoreAsymmRITPIContraction() = delete;
    
    InCoreAsymmRITPIContraction(std::shared_ptr<TwoPInts<IntsT>> tpi):
      InCore4indexTPIContraction<MatsT,IntsT>(tpi) {}

    template <typename MatsU>
    InCoreAsymmRITPIContraction(
        const InCoreAsymmRITPIContraction<MatsU,IntsT> &other, int dummy = 0 ):
      InCoreAsymmRITPIContraction(other.ints_) {}
    template <typename MatsU>
    InCoreAsymmRITPIContraction(
        InCoreAsymmRITPIContraction<MatsU,IntsT> &&other, int dummy = 0 ):
      InCoreAsymmRITPIContraction(other.ints_) {}

    InCoreAsymmRITPIContraction( const InCoreAsymmRITPIContraction &other ):
      InCoreAsymmRITPIContraction(other, 0) {}
    InCoreAsymmRITPIContraction( InCoreAsymmRITPIContraction &&other ):
      InCoreAsymmRITPIContraction(std::move(other), 0) {}

    // Computation interfaces
    virtual void JContract(
        MPI_Comm,
        TwoBodyContraction<MatsT>&) const;

    virtual void KContract(
        MPI_Comm,
        TwoBodyContraction<MatsT>&) const
      {CErr("K Contraction for (ee|pp) integrals are not valid");}

    virtual ~InCoreAsymmRITPIContraction() {}

  }; // class InCoreAsymmRITPIContraction

}; // namespace ChronusQ
