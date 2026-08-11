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

#include <chronusq_sys.hpp>
#include <singleslater.hpp>
#include <basisset/basisset_util.hpp>
#include <cqlinalg/blasext.hpp>
#include <util/timer.hpp>
#include <dft.hpp>
#include <gauxcutils.hpp>
#include <physcon.hpp>
#include <fockbuilder/kcoef.hpp>
#include <particleintegrals/twopints/gtodirecttpi.hpp>
#include <particleintegrals/twopints/impl.hpp>

// KS_DEBUG_LEVEL == 1 - Timing
#ifndef KS_DEBUG_LEVEL
#  define KS_DEBUG_LEVEL 0
#endif

namespace ChronusQ {


  /**
   *  \breif The Kohn--Sham class.
   *
   *  Specializes the SingleSlater class for a Kohn--Sham description of the
   *  many-body wave function
   */ 
  template <typename MatsT, typename IntsT>
  class KohnSham : virtual public SingleSlater<MatsT,IntsT>,
    public std::enable_shared_from_this<KohnSham<MatsT,IntsT>> {

  protected:

    // Useful typedefs
    typedef MatsT*                    oper_t;
    typedef std::vector<oper_t>       oper_t_coll;
    typedef std::vector<oper_t_coll>  oper_t_coll2;

  public:

    std::shared_ptr<KohnSham<MatsT,IntsT>> getPtr(){ return this->shared_from_this(); }


    std::vector<std::shared_ptr<DFTFunctional>> functionals; ///< XC kernels
    IntegrationParam intParam; ///< Numerical integration controls

    bool doVXC_ = true; ///< If this object is responsible for forming VXC
    bool isGGA_; ///< Whether or not the XC kernel is within the GGA
    double XCEnergy = 0.; ///< Exchange-correlation energy (intra-particle only)
    std::vector<std::vector<double> > XCGradient; ///< Exchange-correlation energy gradient

    std::shared_ptr<cqmatrix::PauliSpinorMatrices<IntsT>> VXC; ///< VXC terms

    // Allocated during GauXC setup only for range-separated functionals
    std::shared_ptr<cqmatrix::PauliSpinorMatrices<MatsT>> shortRangeExchangeMatrix;
    std::shared_ptr<TPIContractions<MatsT,IntsT>>         shortRangeExchangeContraction;
    std::shared_ptr<GradInts<TwoPInts,IntsT>>             shortRangeExchangeGradIntegrals;

    // Current Timings
    double VXCDur;

    // Inherit ctors from SingleSlater<T>

    template <typename... Args>
    KohnSham(std::string funcName,
      std::vector<std::shared_ptr<DFTFunctional>> funclist,
      MPI_Comm c, IntegrationParam ip,
      Molecule &mol, BasisSet &basis,
      std::shared_ptr<Integrals<IntsT>> aoi, Args... args) : 
      SingleSlater<MatsT,IntsT>(c,mol,basis,aoi,args...),
      WaveFunctionBase(c,mol,basis,args...),
      QuantumBase(c,args...), isGGA_(false),
      functionals(std::move(funclist)),intParam(ip){ 

      // Append HF tags to reference names
      if(this->nC == 1) {
        if(this->iCS) {
          this->refLongName_  += "Restricted " + funcName;
          this->refShortName_ += "R" + funcName;
        } else {
          this->refLongName_  += "Unrestricted " + funcName;
          this->refShortName_ += "U" + funcName;
        }
      } else {
        this->refLongName_  += "Generalized " + funcName;
        this->refShortName_ += "G" + funcName;
      }

      size_t NB = this->basisSet().nBasis;
      if(this->nC > 1)
        VXC = std::make_shared<cqmatrix::PauliSpinorMatrices<IntsT>>(NB, true);
      else if (not this->iCS)
        VXC = std::make_shared<cqmatrix::PauliSpinorMatrices<IntsT>>(NB, false);
      else
        VXC = std::make_shared<cqmatrix::PauliSpinorMatrices<IntsT>>(NB, false, false);
      VXC->clear();

      // initialize the gradients to be zero
      XCGradient.resize(this->molecule_.atoms.size());
      for(size_t ic = 0; ic < this->molecule_.atoms.size(); ic++) {
        for(size_t xyz = 0; xyz < 3; xyz++) {
          XCGradient[ic].push_back(0.);
        }
      }


    }; // KohnSham constructor


    template <typename... Args>
    KohnSham(std::string rL, std::string rS, std::string funcName,
      std::vector<std::shared_ptr<DFTFunctional>> funclist,
      MPI_Comm c, IntegrationParam ip, 
      Molecule &mol, BasisSet &basis,
      std::shared_ptr<Integrals<IntsT>> aoi, Args... args) : 
      SingleSlater<MatsT,IntsT>(c,mol,basis,aoi,args...),
      WaveFunctionBase(c,mol,basis,args...),
      QuantumBase(c,args...), isGGA_(false),
      functionals(std::move(funclist)),intParam(ip) { 

      this->refLongName_  += rL + " " + funcName;
      this->refShortName_ += rS + funcName;

      size_t NB = this->basisSet().nBasis;
      if(this->nC > 1)
        VXC = std::make_shared<cqmatrix::PauliSpinorMatrices<IntsT>>(NB, true);
      else if (not this->iCS)
        VXC = std::make_shared<cqmatrix::PauliSpinorMatrices<IntsT>>(NB, false);
      else
        VXC = std::make_shared<cqmatrix::PauliSpinorMatrices<IntsT>>(NB, false, false);
      VXC->clear();

      // initialize the gradients to be zero
      XCGradient.resize(this->molecule_.atoms.size());
      for(size_t ic = 0; ic < this->molecule_.atoms.size(); ic++) {
        for(size_t xyz = 0; xyz < 3; xyz++) {
          XCGradient[ic].push_back(0.);
        }
      }

    }; // KohnSham constructor


    // Copy and Move ctors
      
    template <typename MatsU> 
      KohnSham(const KohnSham<MatsU,IntsT> &other, int dummy = 0); 
    template <typename MatsU> 
      KohnSham(KohnSham<MatsU,IntsT> &&other, int dummy = 0);
    KohnSham(const KohnSham<MatsT,IntsT> &other);
    KohnSham(KohnSham<MatsT,IntsT> &&other);


    void setupRangeSeparatedHybridExchange(const ExchCXX::HybCoeffs &rangeSeparatedHybridCoefficients) {
      if (this->nC > 2)
        CErr("Range-separated hybrid exchange is currently implemented only for 1C and 2C references.");
      if (this->basisSet().basisType != REAL_GTO)
        CErr("Range-separated hybrid exchange is currently implemented only for real GTOs.");

      // Build the erfc TPI using the same settings as the regular electronic TPI.
      if (not this->aoints_->shortRangeTPI) {
        this->aoints_->shortRangeTPI = this->aoints_->TPI->createWithKernel(TPI_KERNEL::ShortRangeErfc, rangeSeparatedHybridCoefficients.omega);
        EMPerturbation pert;
        auto begin = tick();
        this->aoints_->shortRangeTPI->computeAOInts(this->basisSet(), this->molecule(), pert,
                                                    ELECTRON_REPULSION, this->fockBuilder->getHamiltonianOptions());
        if (MPIRank(this->comm) == 0)
          std::cout << "    ShortRange-Erfc-K integral setup duration = " << tock(begin) << " s" << std::endl;

        if (auto srtpi = std::dynamic_pointer_cast<InCoreRITPI<IntsT>>(this->aoints_->shortRangeTPI))
          if (auto eri3j = std::dynamic_pointer_cast<DistributedERI3J<IntsT>>(srtpi->eri3j()))
            if (srtpi->redistribute()) eri3j->redistributeToSplitNBRI();
      }

      // Set up contraction and output matrix
      shortRangeExchangeContraction = makeTPIContraction<MatsT,IntsT>(this->aoints_->shortRangeTPI);
      shortRangeExchangeContraction->printContractionTiming = this->scfControls.printContractionTiming;
      shortRangeExchangeMatrix = std::make_shared<cqmatrix::PauliSpinorMatrices<MatsT>>(*this->exchangeMatrix);

      // The erfc-K gradient contribution stays direct for now (TODOAL: support Incore integrals)
      std::shared_ptr<DirectTPI<IntsT>> shortRangeExchangeDirectTPI;
      if (auto directTPI = std::dynamic_pointer_cast<DirectTPI<IntsT>>(this->aoints_->shortRangeTPI)) {
        shortRangeExchangeDirectTPI = directTPI;
      } else {
        shortRangeExchangeDirectTPI =
          std::make_shared<DirectTPI<IntsT>>(this->basisSet(), this->basisSet(), this->molecule(), 1e-12,
                                             DirectTPI<IntsT>::Kernel::ShortRangeErfc, rangeSeparatedHybridCoefficients.omega);
        shortRangeExchangeDirectTPI->computeSchwarz();
      }

      std::vector<std::shared_ptr<DirectTPI<IntsT>>> gradIntegrals(3*this->molecule().nAtoms,shortRangeExchangeDirectTPI);
      shortRangeExchangeGradIntegrals = std::make_shared<GradInts<TwoPInts,IntsT>>(this->basisSet().nBasis,this->molecule().nAtoms,
                                                                                   std::move(gradIntegrals));
    }

    void setupRangeSeparatedHybridExchange() override {
      if (not this->gauxcUtils or not this->gauxcUtils->isRangeSeparatedHybrid()) return;
      // Clear the existing erfc tensor so it is rebuilt for the current geometry.
      if (this->aoints_) this->aoints_->shortRangeTPI = nullptr;
      shortRangeExchangeContraction.reset();
      shortRangeExchangeMatrix.reset();
      shortRangeExchangeGradIntegrals.reset();
      setupRangeSeparatedHybridExchange(this->gauxcUtils->hybridCoefficients);
    }

    /**
     *  \brief Form the short-range erfc exchange correction.
     *
     *  K_RSH = alpha K_full + beta K_erfc. 
     */
    void formRangeSeparatedHybridExchange(EMPerturbation &pert, const ExchCXX::HybCoeffs &rangeSeparatedHybridCoefficients) {

      if (not shortRangeExchangeContraction)
        CErr("Range-separated hybrid exchange was not initialized during setup.");

      shortRangeExchangeMatrix->clear();
      auto shortRangeExchangeBegin = tick();

      auto ritpi_incore = std::dynamic_pointer_cast<InCoreRITPIContraction<MatsT,IntsT>>(shortRangeExchangeContraction);
      auto ritpi_dist   = std::dynamic_pointer_cast<DistributedRITPIContraction<MatsT,IntsT>>(shortRangeExchangeContraction);
      const bool useKCoef = this->nC == 1 and (ritpi_incore or (ritpi_dist and ritpi_dist->canUseKCoef()));

      if (useKCoef) {
        contractExchangeKCoef(this->comm, this->iCS, not this->denEqCoeff_, *this,
                              shortRangeExchangeContraction, *shortRangeExchangeMatrix);
      } else {
        std::vector<TwoBodyContraction<MatsT>> contractions;
        contractions.push_back({this->onePDM->S().pointer(),
                                shortRangeExchangeMatrix->S().pointer(),
                                true, EXCHANGE});
        if (shortRangeExchangeMatrix->hasZ())
          contractions.push_back({this->onePDM->Z().pointer(),
                                  shortRangeExchangeMatrix->Z().pointer(),
                                  true, EXCHANGE});
        if (shortRangeExchangeMatrix->hasXY()) {
          contractions.push_back({this->onePDM->Y().pointer(),
                                  shortRangeExchangeMatrix->Y().pointer(),
                                  true, EXCHANGE});
          contractions.push_back({this->onePDM->X().pointer(),
                                  shortRangeExchangeMatrix->X().pointer(),
                                  true, EXCHANGE});
        }
        shortRangeExchangeContraction->twoBodyContract(this->comm, true, contractions, pert);
      }

      if (MPIRank(this->comm) == 0 and this->scfControls.printContractionTiming)
        std::cout << "        ShortRange-Erfc-K (RSH) " << (useKCoef ? "K-coeff" : "density")
                  << " contraction duration = " << tock(shortRangeExchangeBegin) << " s" << std::endl;

      if (MPIRank(this->comm) == 0) {
        const double shortRangeExchangeCoefficient = rangeSeparatedHybridCoefficients.beta;
        *this->twoeH -= shortRangeExchangeCoefficient * *shortRangeExchangeMatrix;
        *this->fockMatrix -= shortRangeExchangeCoefficient * *shortRangeExchangeMatrix;
      }
    }

    std::vector<double> getShortRangeExchangeGrad(
      EMPerturbation &pert, double shortRangeExchangeCoefficient) {

      if(not shortRangeExchangeGradIntegrals)
        CErr("Range-separated exchange gradient integrals were not initialized.");
      shortRangeExchangeGradIntegrals->computeAOInts(this->basisSet(),this->basisSet(),this->molecule(),pert,
                                                     ELECTRON_REPULSION,this->fockBuilder->getHamiltonianOptions());
      return this->fockBuilder->getShortRangeExchangeGrad(*this,pert,*shortRangeExchangeGradIntegrals,
                                                          shortRangeExchangeCoefficient);
    }


    /**
     *  \brief Kohn-Sham specialization of formFock
     *
     *  Compute VXC and increment the fock matrix
     */  
    virtual void formFock(EMPerturbation &pert, bool increment = false, double HFX = 0.) {
      double xHFX;
      const ExchCXX::HybCoeffs *rangeSeparatedHybridCoefficients = nullptr;
      if (not this->intParam.useGauXC) {
        xHFX = functionals.size() != 0 ? functionals.back()->xHFX : 1.;
      } else {
        if (this->gauxcUtils->isRangeSeparatedHybrid()) {
          rangeSeparatedHybridCoefficients = &this->gauxcUtils->hybridCoefficients;
          xHFX = rangeSeparatedHybridCoefficients->alpha;
        } else {
          xHFX = this->gauxcUtils->xHFX;
        }
      }

      SingleSlater<MatsT,IntsT>::formFock(pert,increment,xHFX);
      if (rangeSeparatedHybridCoefficients)
        formRangeSeparatedHybridExchange(pert, *rangeSeparatedHybridCoefficients);

      if( doVXC_ ) {
        ProgramTimer::tick("Form VXC");
        if (not this->intParam.useGauXC) {
          // Using in-house DFT code to calculate VXC
          formVXC(pert);
          ROOT_ONLY(this->comm);
          // Add VXC in Fock matrix
          *this->fockMatrix += *VXC;
        } else {
          // Using GauXC to calculate VXC 
          // Get system info 
          bool is_gks = this->onePDM->hasZ() and this->onePDM->hasXY();
          bool is_uks = this->onePDM->hasZ() and not this->onePDM->hasXY();
          bool is_dks = this->nC == 4;
          bool is_rks = not is_uks and not is_gks and not is_dks; 

          size_t NB = this->basisSet().nBasis; 

          // Convert CQ matrices to be Eigen matrices to feed into GauXC
          Eigen::Matrix<double, -1, -1> Ps, Pz, Py, Px;
          
          // Initialize return values
          double EXC = 0.0;
          Eigen::MatrixXd VXCs, VXCz, VXCx, VXCy;

          // 4 component LL Approximation
          if (is_dks) {
            // Parse 4C-DFT Options
            bool vll  = false;
            bool full = false;
            auto dkstype = this->fockBuilder->hamiltonianOptions_.dksType;
            if(dkstype == DKS_TYPE::VLL){
            vll = true;
          } else if (dkstype == DKS_TYPE::FULL) {
            full = true;
          } 
            Eigen::Matrix<double, -1, -1> Ps_SS, Pz_SS, Py_SS, Px_SS, Ps_imag, Pz_imag, Py_imag, Px_imag,Ps_SS_imag, Pz_SS_imag, Py_SS_imag, Px_SS_imag;
            
            if (full) {
            Eigen::MatrixXd VXCs_raw, VXCz_raw, VXCx_raw, VXCy_raw, VXC_zero;
            Eigen::MatrixXd VXCs_raw_s, VXCz_raw_s, VXCx_raw_s, VXCy_raw_s;
            Eigen::MatrixXd VXCs(2*NB,2*NB), VXCz(2*NB,2*NB), VXCx(2*NB,2*NB), VXCy(2*NB,2*NB);
            Eigen::MatrixXd VXCs_raw_s_im, VXCz_raw_s_im, VXCx_raw_s_im, VXCy_raw_s_im;
            Eigen::MatrixXcd VXCs_im(2*NB,2*NB), VXCz_im(2*NB,2*NB), VXCx_im(2*NB,2*NB), VXCy_im(2*NB,2*NB);
            Eigen::Matrix<MatsT,Eigen::Dynamic,Eigen::Dynamic> VXCs_SS(NB,NB), VXCz_SS(NB,NB), VXCx_SS(NB,NB), VXCy_SS(NB,NB);

            VXC_zero = Eigen::MatrixXd::Zero(NB,NB);
            VXCs_raw = VXCz_raw = VXCx_raw = VXCy_raw = Eigen::MatrixXd::Zero(NB,NB);
            VXCs_raw_s = VXCz_raw_s = VXCx_raw_s = VXCy_raw_s = Eigen::MatrixXd::Zero(NB,NB);
            VXCs_raw_s_im = VXCz_raw_s_im = VXCx_raw_s_im = VXCy_raw_s_im = Eigen::MatrixXd::Zero(NB,NB);
            VXCs_im = VXCz_im = VXCx_im = VXCy_im = Eigen::MatrixXcd::Zero(2*NB,2*NB);

            VXCs_SS = VXCz_SS = VXCx_SS = VXCy_SS = Eigen::Matrix<MatsT,Eigen::Dynamic,Eigen::Dynamic>::Zero(NB,NB);

            // Initialize needed matrices for Component and Spin Scatter
            bool allocateLLMS = true; 
            bool allocateSS   = true;
            bool allocateLSSL = false;

            auto dummy_pauli = std::make_shared<cqmatrix::PauliSpinorMatrices<MatsT>>(0, false, false);

            auto onePDMLLSCR = std::make_shared<cqmatrix::PauliSpinorMatrices<MatsT>>(NB, true, true);
            auto onePDMSSSCR = std::make_shared<cqmatrix::PauliSpinorMatrices<MatsT>>(NB, true, true);
            auto onePDMLSSCR = std::make_shared<cqmatrix::PauliSpinorMatrices<MatsT>>(NB, true, true);
            auto onePDMSLSCR = std::make_shared<cqmatrix::PauliSpinorMatrices<MatsT>>(NB, true, true);

            auto onePDMLL = allocateLLMS ? onePDMLLSCR: dummy_pauli;
            auto onePDMSS = allocateSS   ? onePDMSSSCR: dummy_pauli;
            auto onePDMLS = allocateLSSL ? onePDMLSSCR: dummy_pauli;
            auto onePDMSL = allocateLSSL ? onePDMSLSCR: dummy_pauli;

            // Scatter 1 Particle Density Matrix into component blocks.
            this->onePDM->componentScatter(*onePDMLL, *onePDMLS, *onePDMSL, *onePDMSS);
        
            // Scatter LL block of density into Pauli spin components.
            Ps = Eigen::Map<Eigen::Matrix<double, -1, -1>>(onePDMLL->real_part().S().pointer(), NB, NB);
            Px = Eigen::Map<Eigen::Matrix<double, -1, -1>>(onePDMLL->real_part().X().pointer(), NB, NB);
            Py = Eigen::Map<Eigen::Matrix<double, -1, -1>>(onePDMLL->real_part().Y().pointer(), NB, NB);
            Pz = Eigen::Map<Eigen::Matrix<double, -1, -1>>(onePDMLL->real_part().Z().pointer(), NB, NB);

            // Ps_imag = Eigen::Map<Eigen::Matrix<double, -1, -1>>(onePDMLL->imag_part().S().pointer(), NB, NB);
            // Px_imag = Eigen::Map<Eigen::Matrix<double, -1, -1>>(onePDMLL->imag_part().X().pointer(), NB, NB);
            // Py_imag = Eigen::Map<Eigen::Matrix<double, -1, -1>>(onePDMLL->imag_part().Y().pointer(), NB, NB);
            // Pz_imag = Eigen::Map<Eigen::Matrix<double, -1, -1>>(onePDMLL->imag_part().Z().pointer(), NB, NB);

            // Scatter SS block of density into Pauli spin components.
            Ps_SS = Eigen::Map<Eigen::Matrix<double, -1, -1>>(onePDMSS->real_part().S().pointer(), NB, NB);
            Px_SS = Eigen::Map<Eigen::Matrix<double, -1, -1>>(onePDMSS->real_part().X().pointer(), NB, NB);
            Py_SS = Eigen::Map<Eigen::Matrix<double, -1, -1>>(onePDMSS->real_part().Y().pointer(), NB, NB);
            Pz_SS = Eigen::Map<Eigen::Matrix<double, -1, -1>>(onePDMSS->real_part().Z().pointer(), NB, NB);

            Ps_SS_imag = Eigen::Map<Eigen::Matrix<double, -1, -1>>(onePDMSS->imag_part().S().pointer(), NB, NB);
            Px_SS_imag = Eigen::Map<Eigen::Matrix<double, -1, -1>>(onePDMSS->imag_part().X().pointer(), NB, NB);
            Py_SS_imag = Eigen::Map<Eigen::Matrix<double, -1, -1>>(onePDMSS->imag_part().Y().pointer(), NB, NB);
            Pz_SS_imag = Eigen::Map<Eigen::Matrix<double, -1, -1>>(onePDMSS->imag_part().Z().pointer(), NB, NB);
            
            // // Prints the particle count for LL and SS, use to compare against the dft integrated values.
            // Eigen::Matrix<double, -1, -1> T, S;
            // T = Eigen::Map<Eigen::Matrix<double, -1, -1>>(reinterpret_cast<double*>(this->aoints_->kinetic->pointer()), NB, NB);
            // S = Eigen::Map<Eigen::Matrix<double, -1, -1>>(reinterpret_cast<double*>(this->aoints_->overlap->pointer()), NB, NB);

            // std::cout<<"trace(Ps*S) 1/2c2*trace(Ps_SS * T)"<<std::endl;
            // std::cout<<(Ps*S).trace()<<" "<<(1./(2*SpeedOfLight*SpeedOfLight))*(Ps_SS*T).trace()<<std::endl;
            // //

            std::tie(EXC, VXCs_raw, VXCz_raw, VXCy_raw, VXCx_raw, VXCs_raw_s, VXCz_raw_s, VXCy_raw_s, VXCx_raw_s , VXCs_raw_s_im, VXCz_raw_s_im, VXCy_raw_s_im, VXCx_raw_s_im) 
                  = this->gauxcUtils->integrator_pointer->eval_exc_vxc( Ps, Pz, Py, Px, Ps_SS, Pz_SS, Py_SS, Px_SS, Ps_SS_imag, Pz_SS_imag, Py_SS_imag, Px_SS_imag );
            
            auto VXCLL = std::make_shared<cqmatrix::PauliSpinorMatrices<MatsT>>(NB, true, true);
            auto VXCSS = std::make_shared<cqmatrix::PauliSpinorMatrices<MatsT>>(NB, true, true);

            VXCLL->clear();
            VXCSS->clear();

            VXCLL->S() += VXCs_raw;
            VXCLL->Z() += VXCz_raw;
            VXCLL->Y() += VXCy_raw;
            VXCLL->X() += VXCx_raw;

            VXCs_SS.real() << VXCs_raw_s;
            VXCz_SS.real() << VXCz_raw_s;
            VXCy_SS.real() << VXCy_raw_s;
            VXCx_SS.real() << VXCx_raw_s;
            
            if constexpr (std::is_same_v<typename decltype(VXCs_SS)::Scalar, dcomplex >){
              VXCs_SS.imag() << VXC_zero;
              VXCz_SS.imag() << VXC_zero;
              VXCy_SS.imag() << VXC_zero;
              VXCx_SS.imag() << VXC_zero;

              // VXCs_SS.imag() << VXCs_raw_s_im;
              // VXCz_SS.imag() << VXCz_raw_s_im;
              // VXCy_SS.imag() << VXCy_raw_s_im;
              // VXCx_SS.imag() << VXCx_raw_s_im;
            }

            VXCSS->S() += VXCs_SS;
            VXCSS->Z() += VXCz_SS;
            VXCSS->Y() += VXCy_SS;
            VXCSS->X() += VXCx_SS;


            this->fockMatrix->componentAdd('N',MatsT(2.),"LL",*VXCLL);
            this->fockMatrix->componentAdd('N',MatsT(2.),"SS",*VXCSS);

            // EXC Energy
            this->XCEnergy = EXC;

          } //End full
          else if (vll){
            Eigen::MatrixXd VXCs_raw, VXCz_raw, VXCx_raw, VXCy_raw, VXC_zero;
            Eigen::MatrixXd VXCs(2*NB,2*NB), VXCz(2*NB,2*NB), VXCx(2*NB,2*NB), VXCy(2*NB,2*NB);
            VXC_zero = Eigen::MatrixXd::Zero(NB,NB);

            Ps = Pz = Py = Px = Eigen::MatrixXd::Zero(NB,NB);
            
            // Initialize needed matrices for Component and Spin Scatter
            bool allocateLLMS = true; 
            bool allocateSS   = false;
            bool allocateLSSL = false;

            auto dummy_pauli = std::make_shared<cqmatrix::PauliSpinorMatrices<MatsT>>(0, false, false);

            auto onePDMLLSCR = std::make_shared<cqmatrix::PauliSpinorMatrices<MatsT>>(NB, true, true);
            auto onePDMSSSCR = std::make_shared<cqmatrix::PauliSpinorMatrices<MatsT>>(NB, true, true);
            auto onePDMLSSCR = std::make_shared<cqmatrix::PauliSpinorMatrices<MatsT>>(NB, true, true);
            auto onePDMSLSCR = std::make_shared<cqmatrix::PauliSpinorMatrices<MatsT>>(NB, true, true);

            auto onePDMLL = allocateLLMS ? onePDMLLSCR: dummy_pauli;
            auto onePDMSS = allocateSS   ? onePDMSSSCR: dummy_pauli;
            auto onePDMLS = allocateLSSL ? onePDMLSSCR: dummy_pauli;
            auto onePDMSL = allocateLSSL ? onePDMSLSCR: dummy_pauli;

            // Scatter 1 Particle Density Matrix into component blocks.
            this->onePDM->componentScatter(*onePDMLL, *onePDMLS, *onePDMSL, *onePDMSS);

            // Scatter LL block of density into Pauli spin components.
            Ps = Eigen::Map<Eigen::Matrix<double, -1, -1>>(onePDMLL->real_part().S().pointer(), NB, NB);
            Px = Eigen::Map<Eigen::Matrix<double, -1, -1>>(onePDMLL->real_part().X().pointer(), NB, NB);
            Py = Eigen::Map<Eigen::Matrix<double, -1, -1>>(onePDMLL->real_part().Y().pointer(), NB, NB);
            Pz = Eigen::Map<Eigen::Matrix<double, -1, -1>>(onePDMLL->real_part().Z().pointer(), NB, NB);


            std::tie(EXC, VXCs_raw, VXCz_raw, VXCy_raw, VXCx_raw) = this->gauxcUtils->integrator_pointer->eval_exc_vxc( Ps, Pz, Py, Px );

            // Form VXCLL mat in Pauli form by padding
            VXCs << VXCs_raw, VXC_zero, VXC_zero, VXC_zero;
            VXCz << VXCz_raw, VXC_zero, VXC_zero, VXC_zero;
            VXCy << VXCy_raw, VXC_zero, VXC_zero, VXC_zero;
            VXCx << VXCx_raw, VXC_zero, VXC_zero, VXC_zero;


            // EXC Energy
            this->XCEnergy = EXC;
 
            // Add VXCLL terms to Fock Matrix
            VXCs *= 2.0;
            this->fockMatrix->S() +=  VXCs;
           
            if(!is_rks){
              VXCz *= 2.0;
              this->fockMatrix->Z() +=  VXCz;
              if(is_gks or is_dks){
                VXCy *= 2.0;
                VXCx *= 2.0;
                this->fockMatrix->Y() +=  VXCy;
                this->fockMatrix->X() +=  VXCx;
            
              }
            }
          } // End VLL
        } // End DKS

        else {
            auto VXC_zero = Eigen::MatrixXd::Zero(NB,NB);

          // Call corresonding epc evaluation functions 
          Ps = Eigen::Map<Eigen::Matrix<double, -1, -1>>(this->onePDM->real_part().S().pointer(), NB, NB); 
          if (is_rks) {                                        
            Ps /= 2.0; // Need to scale by 0.5 due to GauXC's RKS logic
            std::tie(EXC, VXCs) = this->gauxcUtils->integrator_pointer->eval_exc_vxc( Ps );
          } else {
            Pz = Eigen::Map<Eigen::Matrix<double, -1, -1>>(this->onePDM->real_part().Z().pointer(), NB, NB); 
            if (is_uks) {              
              std::tie(EXC, VXCs, VXCz) = this->gauxcUtils->integrator_pointer->eval_exc_vxc( Ps, Pz);
            } else {              
              Py = Eigen::Map<Eigen::Matrix<double, -1, -1>>(this->onePDM->real_part().Y().pointer(), NB, NB); 
              Px = Eigen::Map<Eigen::Matrix<double, -1, -1>>(this->onePDM->real_part().X().pointer(), NB, NB); 

              std::tie(EXC, VXCs, VXCz, VXCy, VXCx) = this->gauxcUtils->integrator_pointer->eval_exc_vxc( Ps, Pz, Py, Px);
            }
          }
          // Assign computed EXC and VXC (with a scaling factor of 2)
          this->XCEnergy = EXC;
          VXCs *= 2.0;
          this->fockMatrix->S() +=  VXCs;
          if(!is_rks){
            VXCz *= 2.0;
            this->fockMatrix->Z() +=  VXCz;
            if(is_gks or is_dks){
              VXCy *= 2.0;
              VXCx *= 2.0;
              this->fockMatrix->Y() +=  VXCy;
              this->fockMatrix->X() +=  VXCx;
            }
          } 
        }
        }  // end GauXC
        ProgramTimer::tock("Form VXC");
      } // end VXC

    }; // formFock

    /**
     *  \brief Kohn-Sham specialization of getGrad
     *
     *  Compute EXC gradient and increment the HF gradient
     */
    std::vector<double> getGrad(EMPerturbation& pert, bool equil, bool saveInts, double xHFX = 1.) {

      const ExchCXX::HybCoeffs *rangeSeparatedHybridCoefficients = nullptr;
      if(not this->intParam.useGauXC) {
        xHFX = functionals.size() != 0 ? functionals.back()->xHFX : 1.;
      } else if(this->gauxcUtils->isRangeSeparatedHybrid()) {
        rangeSeparatedHybridCoefficients = &this->gauxcUtils->hybridCoefficients;
        xHFX = rangeSeparatedHybridCoefficients->alpha;
      } else {
        xHFX = this->gauxcUtils->xHFX;
      }

      size_t nAtoms = this->molecule().nAtoms;
      size_t nGrad = 3*nAtoms;

      // Obtain HF gradient
      std::vector<double> gradient(nGrad, 0.);
      gradient = SingleSlater<MatsT,IntsT>::getGrad(pert,equil,saveInts,xHFX);

      if(rangeSeparatedHybridCoefficients) {
        auto shortRangeGradient = getShortRangeExchangeGrad(pert,rangeSeparatedHybridCoefficients->beta);
        std::transform(gradient.begin(),gradient.end(),shortRangeGradient.begin(),gradient.begin(),std::plus<double>());
      }

      for(size_t ic = 0; ic < nAtoms; ic++) 
        for(size_t XYZ = 0; XYZ < 3; XYZ++) 
          this->XCGradient[ic][XYZ] = 0.0;


      // TangDD Remove before merge
      //auto magAmp = pert.getDipoleAmp(Magnetic);
      //std::cout<<"magAmp in ks 0: "<<magAmp[0]<<" 1: "<<magAmp[1]<<" 2: "<<magAmp[2]<<std::endl;

      formEXCGradient(pert);

      //std::cout << "Main XC Gradient:" << std::endl; 
      for(size_t ic = 0; ic < nAtoms; ic++) {
        for(size_t XYZ = 0; XYZ < 3; XYZ++) {
          //std::cout << std::setprecision(8);
          //std::cout << std::setw(16) << this->XCGradient[ic][XYZ] << " ";
          gradient[ic*3+XYZ] += this->XCGradient[ic][XYZ];
        }
        std::cout << std::endl;
      }
      //std::cout << std::setprecision(16) << "XC Energy in computeGradients(): " << this->XCEnergy << std::endl;

      return gradient;
    }


    /**
     *  \brief Kohn-Sham specialization of computeEnergy
     *
     *  Compute EXC and add it to the HF energy 
     */
    using QuantumBase::computeEnergy;
    virtual void computeEnergy() {
      SingleSlater<MatsT,IntsT>::computeEnergy();
      // Add EXC in the total energy
      //std::cout<<"XCEnergy in Hartree "<< XCEnergy <<std::endl;
      this->totalEnergy += XCEnergy;
    }; // computeEnergy



    virtual void printFockTimings(std::ostream &out) {
  
      out << "    Fock Timings:\n";
      out << "      Wall time G[D] = " << std::setw(8)
          << std::setprecision(5)  << std::scientific
          << this->GDDur << " s\n";
      out << "      Wall time VXC  = " << std::setw(8)
          << std::setprecision(5)  << std::scientific
          << VXCDur << " s\n\n";
  
  
    }; // SingleSlater<T>::printFockTimings





    // KS specific functions
    // See include/singleslater/kohnsham/vxc.hpp for docs.

    // VXC
    void formVXC(EMPerturbation&); 

    void formEXCGradient(EMPerturbation&);
    void formEXCGradientInHouse(EMPerturbation&);
    void formEXCGradientGauXC(EMPerturbation&);

    // FXC Terms
    template <typename U>
    void evalTransDen(SHELL_EVAL_TYPE typ, size_t NPts,size_t NBE, size_t NB, 
      std::vector<std::pair<size_t,size_t>> &subMatCut, U *SCR1,
      U *SCR2, U *DENMAT, U *Den, U *GDenX, U *GDenY, U *GDenZ,
      U *BasisScr);

    void loadFXCder(size_t NPts, double *Den, double *sigma, double *EpsEval, double *VRhoEval, 
      double *V2RhoEval, double *VsigmaEval, double *V2sigmaEval, double *V2RhosigmaEval, 
      double *EpsSCR, double *VRhoSCR, double *VsigmaSCR, double *V2RhoEvalSCR, double *V2sigmaEvalSCR,
      double *V2RhosigmaEvalSCR); 

    template <typename U>
    void constructZVarsFXC(DENSITY_TYPE denTyp, bool isGGA, size_t NPts, 
      double* GDenS, double* GDenZ, double* GDenY, double* GDenX,
      U* TS, U* TZ, U* TY, U* TX,
      U* GTS, U* GTZ, U* GTY, U* GTX,
      double *VrhoEval, double *VsigmaEval, 
      double *V2rhoEval, double *V2sigmaEval, double *V2RhosigmaEval, 
      U *ZrhoVar1, U *ZgammaVar1, U *ZgammaVar2, U *ZgammaVar3, U *ZgammaVar4);

    template <typename U>
    void constructZVarsFXC(DENSITY_TYPE denTyp, bool isGGA, size_t NPts, 
      double* GDenS, double* GDenZ, double* GDenY, double* GDenX,
      bool * Msmall, double *Mnorm, 
      double *Kx, double *Ky, double *Kz, 
      double *Hx, double *Hy, double *Hz,
      U* TS, U* TZ, U* TY, U* TX,
      U* GTS, U* GTZ, U* GTY, U* GTX,
      U* gPTss, U* gPTsz, U* gPTsy, U* gPTsx, U* gPTzz, 
      U* gPTyy, U* gPTxx,  
      double *VrhoEval, double *VsigmaEval, 
      double *V2rhoEval, double *V2sigmaEval, double *V2RhosigmaEval, 
      U *ZrhoVar1, U *ZgammaVar1, U *ZgammaVar2, U *ZgammaVar3, U *ZgammaVar4);

    template <typename U>
    void formZ_fxc(DENSITY_TYPE denType, bool isGGA, size_t NPts, size_t NBE, size_t IOff,
      double epsScreen, std::vector<double> &weights,
      U *ZrhoVar1, U *ZgammaVar1, U *ZgammaVar2, U *ZgammaVar3, U *ZgammaVar4,
      double* GDenS, double* GDenZ, double* GDenY, double* GDenX, U* GTS, U* GTZ, U* GTY, U* GTX,
      double *BasisScr, U* ZMAT);

    void formZ_fxc(DENSITY_TYPE denType, bool isGGA, size_t NPts, size_t NBE, size_t IOff,
      double epsScreen, std::vector<double> &weights,
      dcomplex *ZrhoVar1, dcomplex *ZgammaVar1, dcomplex *ZgammaVar2, dcomplex *ZgammaVar3, dcomplex *ZgammaVar4,
      double* GDenS, double* GDenZ, double* GDenY, double* GDenX, dcomplex* GTS, dcomplex* GTZ, dcomplex* GTY, dcomplex* GTX,
      dcomplex *BasisScr, dcomplex* ZMAT);

    // GTO-based TDDFT
    template <typename U>
    void formZ_fxc(DENSITY_TYPE denType, bool isGGA, size_t NPts, size_t NBE, size_t IOff,
      double epsScreen, std::vector<double> &weights,
      U *ZrhoVar1, U *ZgammaVar1, U *ZgammaVar2, U *ZgammaVar3, U *ZgammaVar4,
      bool * Msmall, double *Mnorm, 
      double* DSDMnorm, double* signMD, 
      double* GDenS, double* GDenZ, double* GDenY, double* GDenX, 
      double *Kx, double *Ky, double *Kz, 
      double *Hx, double *Hy, double *Hz,
      U* GTS, U* GTZ, U* GTY, U* GTX,
      U* gPTss, U* gPTsz, U* gPTsy, U* gPTsx, U* gPTzz, 
      U* gPTyy, U* gPTxx,  
      double *BasisScr, U* ZMAT);

    void formZ_fxc(DENSITY_TYPE denType, bool isGGA, size_t NPts, size_t NBE, size_t IOff,
      double epsScreen, std::vector<double> &weights,
      dcomplex *ZrhoVar1, dcomplex *ZgammaVar1, dcomplex *ZgammaVar2, dcomplex *ZgammaVar3, dcomplex *ZgammaVar4,
      bool * Msmall, double *Mnorm, 
//      double* DenS, double* DenZ, double* DenY, double* DenX, 
      double* DSDMnorm, double* signMD, 
      double* GDenS, double* GDenZ, double* GDenY, double* GDenX, 
      double *Kx, double *Ky, double *Kz, 
      double *Hx, double *Hy, double *Hz,
      dcomplex* GTS, dcomplex* GTZ, dcomplex* GTY, dcomplex* GTX,
      dcomplex* gPTss, dcomplex* gPTsz, dcomplex* gPTsy, dcomplex* gPTsx, dcomplex* gPTzz, 
      dcomplex* gPTyy, dcomplex* gPTxx,  
      dcomplex *BasisScr, dcomplex* ZMAT);

    // Calculate gPTss,sx,sy,sz 
    template <typename U>
    void mkgPTVar( 
      size_t NPts, 
      double* GDenS, double* GDenZ, double* GDenY, double* GDenX, 
      U* GTS, U* GTZ, U* GTY, U* GTX,
      U* gPTss, U* gPTsz, U* gPTsy, U* gPTsx, U* gPTzz, 
      U* gPTyy, U* gPTxx  
      );


    template <typename U>
    void formFXC(MPI_Comm c,  std::vector<TwoBodyContraction<U>> &cList, EMPerturbation& );



    // SCF Functions
    void buildOrbitalModifierOptions();
    void computeFullNRStep(MatsT*);
    std::pair<double,MatsT*> getStab();

  }; // class KohnSham


}; // namespace ChronusQ

