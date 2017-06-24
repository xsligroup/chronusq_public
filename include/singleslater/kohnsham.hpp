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
    double XCEnergy; ///< Exchange-correlation energy
    std::vector<std::vector<double> > XCGradient; ///< Exchange-correlation energy gradient

    std::shared_ptr<cqmatrix::PauliSpinorMatrices<IntsT>> VXC; ///< VXC terms

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


    /**
     *  \brief Kohn-Sham specialization of formFock
     *
     *  Compute VXC and increment the fock matrix
     */  
    virtual void formFock(EMPerturbation &pert, bool increment = false, double HFX = 0.) {

      double xHFX = functionals.size() != 0 ? functionals.back()->xHFX : 1.;

      SingleSlater<MatsT,IntsT>::formFock(pert,increment,xHFX);

      ProgramTimer::tick("Form VXC");
      if( doVXC_ ) {
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
          bool is_rks = not is_uks and not is_gks; 
          size_t NB = this->basisSet().nBasis; 

          // Convert CQ matrices to be Eigen matrices to feed into GauXC
          Eigen::Matrix<double, -1, -1> Ps, Pz, Py, Px;
          Ps = Eigen::Map<Eigen::Matrix<double, -1, -1>>(this->onePDM->real_part().S().pointer(), NB, NB); 

          // Initialize return values
          double EXC = 0.0;
          Eigen::MatrixXd VXCs, VXCz, VXCx, VXCy;

          // Call corresonding epc evaluation functions 
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
            if(is_gks){
              VXCy *= 2.0;
              VXCx *= 2.0;
              this->fockMatrix->Y() +=  VXCy;
              this->fockMatrix->X() +=  VXCx;
            }
          } 
        }  // end GauXC
      } // end VXC
      ProgramTimer::tock("Form VXC"); 

    }; // formFock

    /**
     *  \brief Kohn-Sham specialization of getGrad
     *
     *  Compute EXC gradient and increment the HF gradient
     */
    std::vector<double> getGrad(EMPerturbation& pert, bool equil, bool saveInts, double xHFX = 1.) {
      
      xHFX = functionals.size() != 0 ? functionals.back()->xHFX : 1.;

      size_t nAtoms = this->molecule().nAtoms;
      size_t nGrad = 3*nAtoms;

      // Obtain HF gradient
      std::vector<double> gradient(nGrad, 0.);
      gradient = SingleSlater<MatsT,IntsT>::getGrad(pert,equil,saveInts,xHFX);

      for(size_t ic = 0; ic < nAtoms; ic++) 
        for(size_t XYZ = 0; XYZ < 3; XYZ++) 
          this->XCGradient[ic][XYZ] = 0.0;

      formEXCGradient();

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
      this->totalEnergy += XCEnergy;

      //std::cout << std::setprecision(16) << "XC Energy in computeEnergy(): " << this->XCEnergy << std::endl;
        
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

    void formEXCGradient();

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

