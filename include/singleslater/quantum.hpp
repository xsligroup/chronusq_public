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
//#define _DEBUGGIAO 
//#define seperatemag
//#define _DEBUG_EWDM
//#define _DEBUG_4C_ANGULAR

#include <singleslater.hpp>
#include <singleslater/multiparticless.hpp>
#include <cqlinalg/blasext.hpp>
#include <cqlinalg/blasutil.hpp>
#include <cqlinalg/blas3.hpp>
#include <quantum/properties.hpp>
#include <mointstransformer.hpp>
#include <util/timer.hpp>
#include <algorithm>
#include <array>
#include <complex>
#include <type_traits>

namespace ChronusQ {


  template <typename MatsT, typename IntsT>
  void SingleSlater<MatsT,IntsT>::constructRDMBuilder()
  {
    if(this->RDMBuilder != nullptr) return;
    RDM_BUILDER_TYPE rdmType;
    if(this->particle.charge > 0) // Handle protons
      rdmType = this->scfControls.protrdmBuilderType;      
    else
      rdmType = this->scfControls.rdmBuilderType;

    auto* fb = this->fockBuilder.get();
    if (auto* neofb = dynamic_cast<NEOFockBuilder<MatsT, IntsT>*>(fb)) {
      fb = neofb->getNonNEOUpstream();              // still a raw pointer
    }
    bool iRO = (dynamic_cast<ROFock<MatsT, IntsT>*>(fb) != nullptr);
    if(rdmType == RDM_BUILDER_TYPE::MOM && iRO) CErr("MOM with ROHF NYI");

    if(rdmType == RDM_BUILDER_TYPE::MOM)
      this->RDMBuilder = std::make_shared<MOMRDMBuilder<MatsT,IntsT>>(*this);
    else if(rdmType == RDM_BUILDER_TYPE::NEOSTATEAVERAGE)
      this->RDMBuilder = std::make_shared<NEOStateAveragedRDMBuilder<MatsT,IntsT>>(this->scfControls.NEOStateAverageNStates);
    else if(rdmType == RDM_BUILDER_TYPE::FINITETEMP)
      CErr("NEO Finite Temp NYI!");
    else
      this->RDMBuilder = std::make_shared<AufbauRDMBuilder<MatsT,IntsT>>();

  }


  /**
   *  \brief Forms the 1PDM using a set of orbitals 
   *
   *  specialization of Quantum::formDensity. Populates / overwrites
   *  onePDM storage
   * 
   *  NOTE: This function assumes the MO's are in the AO basis
   */ 
  template <typename MatsT, typename IntsT>
  void SingleSlater<MatsT,IntsT>::formDensity() {

    if(MPIRank(comm)==0)
      this->constructRDMBuilder();

    size_t NB  = this->nAlphaOrbital() * nC;
    size_t NB2 = NB * NB;
    auto* fb = this->fockBuilder.get();
    if (auto* neofb = dynamic_cast<NEOFockBuilder<MatsT, IntsT>*>(fb)) {
      fb = neofb->getNonNEOUpstream();              // still a raw pointer
    }
    bool iRO = (dynamic_cast<ROFock<MatsT, IntsT>*>(fb) != nullptr);
    cqmatrix::Matrix<MatsT> temp(NB);

    // ROHF copy modified orbitals to redundant set
    if( iRO ){
      std::copy_n(this->mo[0].pointer(),NB*NB,this->mo[1].pointer());
      std::copy_n(this->eps1,NB,this->eps2);
    }

    // Form the 1RDM on the root MPI process for similar reasons to below
    if(MPIRank(comm)==0)
      this->RDMBuilder->buildRDM(*this);

    // Form the 1PDM on the root MPI process as slave processes
    // do not posses the up-to-date MO coefficients
    if( MPIRank(comm) == 0 ) {

      if(nC == 1) {

        cqmatrix::Matrix<MatsT> DA(NB);

        //this->mo[0].output(std::cout, "mo1", true);

        // DA = CA * RDM * CA**H
        blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans, blas::Op::ConjTrans, NB, NB, NB, MatsT(1.), this->oneRDM->pointer(), NB,
            this->mo[0].pointer(), NB, MatsT(0.), temp.pointer(), NB);
        blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans, blas::Op::NoTrans, NB, NB, NB, MatsT(1.), this->mo[0].pointer(), NB,
            temp.pointer(), NB, MatsT(0.), DA.pointer(), NB);

        if(iCS) {

          // DS = 2 * DA
          *this->onePDM = cqmatrix::PauliSpinorMatrices<MatsT>::spinBlockScatterBuild(DA);

        } else {

          cqmatrix::Matrix<MatsT> DB(NB);
          temp.clear();

          // DB = CB * RDM * CB**H
          blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans, blas::Op::ConjTrans, NB, NB, NB, MatsT(1.), this->oneRDMB->pointer(), NB,
              this->mo[1].pointer(), NB, MatsT(0.), temp.pointer(), NB);
          blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans, blas::Op::NoTrans, NB, NB, NB, MatsT(1.), this->mo[1].pointer(), NB,
              temp.pointer(), NB, MatsT(0.), DB.pointer(), NB);

          // DS = DA + DB
          // DZ = DA - DB
          *this->onePDM = cqmatrix::PauliSpinorMatrices<MatsT>::spinBlockScatterBuild(DA,DB);

        }
      } else {

        // 2C or 4C cases
        cqmatrix::Matrix<MatsT> spinBlockForm(NB);

        if( nC == 2 ) {
          blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans, blas::Op::ConjTrans, NB, NB, NB, MatsT(1.), this->oneRDM->pointer(), NB,
              this->mo[0].pointer(), NB, MatsT(0.), temp.pointer(), NB);
          blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans, blas::Op::NoTrans, NB, NB, NB, MatsT(1.), this->mo[0].pointer(), NB,
              temp.pointer(), NB, MatsT(0.), spinBlockForm.pointer(), NB);
        } else if( nC == 4 ) {
          blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans, blas::Op::ConjTrans, NB, NB, this->nO, MatsT(1.), this->mo[0].pointer()+(NB/2)*NB, NB,
              this->mo[0].pointer()+(NB/2)*NB, NB, MatsT(0.), spinBlockForm.pointer(), NB);
        }

        *this->onePDM = spinBlockForm.template spinScatter<MatsT>();

      }
      ao2orthoDen();
    }


#ifdef CQ_ENABLE_MPI

    // Broadcast the 1PDM to all MPI processes
    if( MPISize(comm) > 1 ) {
      std::cerr  << "  *** Scattering the 1PDM ***\n";
      for(auto p : this->onePDM->SZYXPointers())
        MPIBCast(p,NB*NB/(std::min(2,nC))/(std::min(2,nC)),0,comm);
    }

#endif


#if 0
      print1PDM(std::cerr);
#endif

  }; // SingleSlater<T>::formDensity


  /**
   *  \brief Computes the total field free energy of a single slater determinent
   *
   *  Given a 1PDM and a Fock matrix (specifically the core Hamiltonian and
   *  G[D]), compute the energy.
   *
   *  \warning Assumes that density and Fock matrix have the appropriate form
   *
   *  F(S) = F(A) + F(B)
   *  F(Z) = F(A) - F(B)
   *  ...
   *
   *  \f[
   *     E = \frac{1}{2} \mathrm{Tr}[\mathbf{P}(\mathbf{H} + \mathbf{F})] +
   *     V_{NN}
   *  \f]
   *
   *  Specialization of Quantum<T>::computeEnergy, populates / overwrites 
   *  OBEnergy and MBEnergy
   */ 
  template <typename MatsT, typename IntsT>
  void SingleSlater<MatsT,IntsT>::computeEnergy() {

    ROOT_ONLY(comm);

    // Scalar core hamiltonian contribution to the energy
    this->OBEnergy = 
      this->template computeOBProperty<DENSITY_TYPE::SCALAR>(
         coreH->S().pointer());

#ifdef _DEBUGGIAO
    std::cout<<"One-body energy in Hartree "<< 0.5*std::real(this->OBEnergy) <<std::endl;
#endif

    // One body Spin Orbit
    double SOEnergy = 0.;
    if (coreH->hasZ())
      SOEnergy =
        this->template computeOBProperty<DENSITY_TYPE::MZ>(
           coreH->Z().pointer());
    if (coreH->hasXY()) {
      SOEnergy +=
        this->template computeOBProperty<DENSITY_TYPE::MY>(
           coreH->Y().pointer());
      SOEnergy +=
        this->template computeOBProperty<DENSITY_TYPE::MX>(
           coreH->X().pointer());
    }

 
    this->OBEnergy += SOEnergy;
#ifdef _DEBUGGIAO
    std::cout<<"Spin Zeeman energy in Hartree "<< 0.5*std::real(SOEnergy) <<std::endl;
#endif

    this->OBEnergy *= 0.5;

    // Compute many-body contribution to energy
    // *** These calls are safe as proper zeros are returned by
    // property engine ***
    this->MBEnergy = 
      this->template computeOBProperty<DENSITY_TYPE::SCALAR>(
          twoeH->S().pointer());

    if (twoeH->hasZ())
      this->MBEnergy +=
        this->template computeOBProperty<DENSITY_TYPE::MZ>(
            twoeH->Z().pointer());
    if (twoeH->hasXY()) {
      this->MBEnergy +=
        this->template computeOBProperty<DENSITY_TYPE::MY>(
            twoeH->Y().pointer());
      this->MBEnergy +=
        this->template computeOBProperty<DENSITY_TYPE::MX>(
            twoeH->X().pointer());
    }

    this->MBEnergy *= 0.25;
#ifdef _DEBUGGIAO
    std::cout<<"Many-body energy in Hartree "<< std::real(this->MBEnergy) <<std::endl;
    std::cout<<"Nuc energy in Hartree "<< std::real(this->molecule().nucRepEnergy) <<std::endl;
    std::cout<<"E-field energy in Hartree "<< std::real(this->extraEnergy) <<std::endl;
#endif

    // Assemble total energy
    this->totalEnergy = 
      this->OBEnergy + this->MBEnergy + this->molecule().nucRepEnergy
       + this->extraEnergy;

    // Sanity checks
    assert( not std::isnan(this->OBEnergy) );
    assert( not std::isnan(this->MBEnergy) );
    assert( not std::isinf(this->OBEnergy) );
    assert( not std::isinf(this->MBEnergy) );

  }; // SingleSlater<T>::computeEnergy

  template <typename MatsT, typename IntsT>
  std::vector<double> SingleSlater<MatsT,IntsT>::getEnergySummary() {
    std::vector<double> result = QuantumBase::getEnergySummary();
    result.push_back(this->molecule().nucRepEnergy);
    return result;
  };

  template <typename MatsT, typename IntsT>
  void SingleSlater<MatsT,IntsT>::computeMultipole(EMPerturbation &pert, const std::vector<PROPERTY> &properties) {
    ROOT_ONLY(comm);

    // By default, compute all multipole properties
    std::vector<PROPERTY> propsToCompute = properties;
    if(propsToCompute.empty()) propsToCompute = {ELECTRIC_DIPOLE, ELECTRIC_QUADRUPOLE, ELECTRIC_OCTUPOLE};

    auto hasProperty = [&](PROPERTY p){ return std::find(propsToCompute.begin(), propsToCompute.end(), p) != propsToCompute.end(); };

    if(this->nC == 4 and hasProperty(ELECTRIC_DIPOLE)){
      compute4CDipole(pert);
      return;
    }

    if(this->nC == 2 and pchgDipole_[0] and pchgDipole_[1] and pchgDipole_[2] and hasProperty(ELECTRIC_DIPOLE)){
      computeFockX2CDipole(pert);
      return;
    }

    // Compute elecric contribution to the dipoles
    if(hasProperty(ELECTRIC_DIPOLE)){
    for(auto iXYZ = 0; iXYZ < 3; iXYZ++) 
      this->elecDipole[iXYZ] = -this->template computeOBProperty<SCALAR>((*this->aoints_->lenElectric)[iXYZ]->pointer());


    // Nuclear contributions to the dipoles
    for(auto &atom : this->molecule().atoms){
      if (atom.quantum) continue;  
      MatAdd('N','N',3,1,1.,&this->elecDipole[0],3,atom.nucCharge,
          &atom.coord[0],3,&this->elecDipole[0],3);
    }
    } // if(hasProperty(ELECTRIC_DIPOLE))


#ifdef seperatemag
  // Here output orbital Zeeman energy

    dcomplex onei = dcomplex(0,1);
    auto magAmp = pert.getDipoleAmp(Magnetic);
    // std::cout<<"magnetic field Bx= "<<magAmp[0] <<" By= "<<magAmp[1]
    //   <<" Bz= "<<magAmp[2]<<std::endl;
    //size_t DSize = memManager.template getSize(this->onePDM[SCALAR]); 
    // std::cout<<"DSize="<<DSize<<std::endl;

    // this part add the orbital Zeeman energy contribution  
    double OrbitZeeman = 0.0;  

    const std::array<std::string,3> dipoleList =
      { "X","Y","Z" };

    int NB = 0 + this->basisSet().nBasis;
    //prettyPrintSmart(std::cout,"L",(*this->aoints_->magnetic)[dipoleList[2]]->pointer(),NB,NB,NB);

    for ( auto index = 0 ; index < 3 ; index++ ) {
      OrbitZeeman += std::real( 0.5 * magAmp[index] * onei * this->particle.charge * (1.0/this->particle.mass) * 
        blas::dot(NB*NB, this->onePDM->S().pointer(), 1, (*this->aoints_->magnetic)[dipoleList[index]]->pointer(), 1 ));
    } // for ( auto inde = 0 ; inde < 3 ; inde++ ) 

    std::cout<<"Orbit Zeeman contribution "<<OrbitZeeman<<std::endl; 

// Here we calculate diamgnetic contribution
    double diamagcontrib = 0.0;

    const std::array<std::string,3> diagindex =
      { "XX","YY","ZZ" };

    double diagcoeff[3];
    diagcoeff[0] = 1.0/8.0*(magAmp[1]*magAmp[1]+magAmp[2]*magAmp[2]); 
    diagcoeff[1] = 1.0/8.0*(magAmp[0]*magAmp[0]+magAmp[2]*magAmp[2]);    
    diagcoeff[2] = 1.0/8.0*(magAmp[0]*magAmp[0]+magAmp[1]*magAmp[1]);    

    // add diagonal part

    for ( auto index = 0 ; index < 3 ; index++ ) { 
      diamagcontrib += std::real( diagcoeff[index] * ((-1.0 * this->particle.charge) * 1.0/this->particle.mass) *
        this->template computeOBProperty<DENSITY_TYPE::SCALAR>((*this->aoints_->lenElectric)[diagindex[index]]->pointer()));
    }  

    const std::array<std::string,3> offindex =
      { "XY","XZ","YZ" };

    double offcoeff[3];
    offcoeff[0] = -1.0/4.0*magAmp[0]*magAmp[1];
    offcoeff[1] = -1.0/4.0*magAmp[0]*magAmp[2];
    offcoeff[2] = -1.0/4.0*magAmp[1]*magAmp[2];

    for ( auto index = 0 ; index < 3 ; index++ ) { 
      diamagcontrib += std::real( offcoeff[index] * (-1.0 * this->particle.charge) * (1.0/this->particle.mass) * 
        this->template computeOBProperty<DENSITY_TYPE::SCALAR>((*this->aoints_->lenElectric)[diagindex[index]]->pointer()));
    } 

    std::cout<<"Diamagnetic contribution = "<<diamagcontrib<<std::endl;

#endif

    // Electric contribution to the quadrupoles
    if(hasProperty(ELECTRIC_QUADRUPOLE)){
    for(size_t iXYZ = 0, iX = 0; iXYZ < 3; iXYZ++)
    for(size_t jXYZ = iXYZ     ; jXYZ < 3; jXYZ++, iX++){

      this->elecQuadrupole[iXYZ][jXYZ] = -
        this->template computeOBProperty<SCALAR>((*this->aoints_->lenElectric)[iX+3]->pointer());
      
      this->elecQuadrupole[jXYZ][iXYZ] = this->elecQuadrupole[iXYZ][jXYZ]; 
    }
    
    // Nuclear contributions to the quadrupoles
    for(auto &atom : this->molecule().atoms){
    
    if (atom.quantum) continue;

    for(size_t iXYZ = 0; iXYZ < 3; iXYZ++)
    for(size_t jXYZ = 0; jXYZ < 3; jXYZ++) 
      this->elecQuadrupole[iXYZ][jXYZ] +=
        atom.nucCharge * atom.coord[iXYZ] * atom.coord[jXYZ];
    }
    } // if(hasProperty(ELECTRIC_QUADRUPOLE))

    // Electric contribution to the octupoles
    if(hasProperty(ELECTRIC_OCTUPOLE)){
    for(size_t iXYZ = 0, iX = 0; iXYZ < 3; iXYZ++)
    for(size_t jXYZ = iXYZ     ; jXYZ < 3; jXYZ++)
    for(size_t kXYZ = jXYZ     ; kXYZ < 3; kXYZ++, iX++){

      this->elecOctupole[iXYZ][jXYZ][kXYZ] = -
        this->template computeOBProperty<SCALAR>(
          (*this->aoints_->lenElectric)[iX+9]->pointer());

      this->elecOctupole[iXYZ][kXYZ][jXYZ] = this->elecOctupole[iXYZ][jXYZ][kXYZ]; 

      this->elecOctupole[jXYZ][iXYZ][kXYZ] = this->elecOctupole[iXYZ][jXYZ][kXYZ]; 

      this->elecOctupole[jXYZ][kXYZ][iXYZ] = this->elecOctupole[iXYZ][jXYZ][kXYZ]; 

      this->elecOctupole[kXYZ][iXYZ][jXYZ] = this->elecOctupole[iXYZ][jXYZ][kXYZ]; 

      this->elecOctupole[kXYZ][jXYZ][iXYZ] = this->elecOctupole[iXYZ][jXYZ][kXYZ]; 
    }

    // Nuclear contributions to the octupoles
    for(auto &atom : this->molecule().atoms){

    if (atom.quantum) continue;

    for(size_t iXYZ = 0; iXYZ < 3; iXYZ++)
    for(size_t jXYZ = 0; jXYZ < 3; jXYZ++)
    for(size_t kXYZ = 0; kXYZ < 3; kXYZ++)
      this->elecOctupole[iXYZ][jXYZ][kXYZ] +=
        atom.nucCharge * atom.coord[iXYZ] * atom.coord[jXYZ] *
        atom.coord[kXYZ];
    }
    } // if(hasProperty(ELECTRIC_OCTUPOLE))
    //std::cout << std::string(this->particle.charge>0 ? "Protonic" : "Electronic") << " Subsystem Dipole: " 
    //    << this->elecDipole[0] << " " << this->elecDipole[1] << " " << this->elecDipole[2] << std::endl;
  };


  template <typename MatsT, typename IntsT>
  void SingleSlater<MatsT,IntsT>::computeSpin() {

    ROOT_ONLY(comm);

    auto clampTiny = [](auto &value) {
      using ValueT = std::decay_t<decltype(value)>;
      if (std::abs(value) < 1e-12) value = ValueT(0);
    };

    this->SExpect[0] = 0.5 * this->template computeOBProperty<MX>(
      this->aoints_->overlap->pointer());
    this->SExpect[1] = 0.5 * this->template computeOBProperty<MY>(
      this->aoints_->overlap->pointer());
    this->SExpect[2] = 0.5 * this->template computeOBProperty<MZ>(
      this->aoints_->overlap->pointer());

    if( not this->onePDM->hasZ() and not this->onePDM->hasXY() ) this->SSq = 0;
    else {
      size_t NB = this->basisSet().nBasis;
      MatsT * SCR  = CQMemManager::get().malloc<MatsT>(NB*NB);
      std::fill_n(SCR, NB*NB, MatsT(0.));
      MatsT * SCR2 = CQMemManager::get().malloc<MatsT>(NB*NB);
      std::fill_n(SCR2, NB*NB, MatsT(0.));


      // SCR2 = S * D(S) * S
    /*      
      blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::ConjTrans,NB,NB,NB,T(1.),aoints_->overlap,NB,this->onePDM[SCALAR],NB,
        T(0.),SCR,NB);
      blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::ConjTrans,NB,NB,NB,T(1.),SCR,NB,aoints_->overlap,NB,T(0.),SCR2,NB);
    */
      blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,NB,NB,NB,MatsT(1.),this->aoints_->overlap->pointer(),NB,
           this->onePDM->S().pointer(),NB,MatsT(0.),SCR,NB);
      blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::ConjTrans,NB,NB,NB,MatsT(1.),this->aoints_->overlap->pointer(),NB,
           SCR,NB,MatsT(0.),SCR2,NB);

      
      this->SSq = 3 * this->nO - (3./2.) * 
        this->template computeOBProperty<SCALAR>(SCR2);
  

       if (this->onePDM->hasZ()) {
         // SCR2 = D(Z) * S * D(Z)
         blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,NB,NB,NB,MatsT(1.),this->aoints_->overlap->pointer(),NB,
           this->onePDM->Z().pointer(),NB,MatsT(0.),SCR,NB);
         blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::ConjTrans,NB,NB,NB,MatsT(1.),this->aoints_->overlap->pointer(),NB,
           SCR,NB,MatsT(0.),SCR2,NB);
        
         this->SSq += 0.5 * this->template computeOBProperty<MZ>(SCR2);
       }

      if( this->onePDM->hasXY() ) {
  
        // SCR2 = D(Y) * S * D(Y)
        blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,NB,NB,NB,MatsT(1.),this->aoints_->overlap->pointer(),NB,
             this->onePDM->Y().pointer(),NB,MatsT(0.),SCR,NB);
        blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::ConjTrans,NB,NB,NB,MatsT(1.),this->aoints_->overlap->pointer(),NB,
             SCR,NB,MatsT(0.),SCR2,NB);
        
        this->SSq += 0.5 * this->template computeOBProperty<MY>(SCR2);


        // SCR2 = D(X) * S * D(X)
        blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,NB,NB,NB,MatsT(1.),this->aoints_->overlap->pointer(),NB,
             this->onePDM->X().pointer(),NB,MatsT(0.),SCR,NB);
        blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::ConjTrans,NB,NB,NB,MatsT(1.),this->aoints_->overlap->pointer(),NB,
             SCR,NB,MatsT(0.),SCR2,NB);
        
        this->SSq += 0.5 * this->template computeOBProperty<MX>(SCR2);

      }

      for(auto i = 0; i < 3; i++) 
        this->SSq += 4 * this->SExpect[i] * this->SExpect[i];

      this->SSq *= 0.25;

      for (auto &component : this->SExpect) clampTiny(component);
      clampTiny(this->SSq);

      // std::cout << "Spin Expectation Values: Sx = " << this->SExpect[0] << ", Sy = " << this->SExpect[1] << ", Sz = " << this->SExpect[2] << std::endl;
      // std::cout << "Spin Expectation Value: S^2 = " << this->SSq << std::endl;

      CQMemManager::get().free(SCR,SCR2);
    }

  };
  
template <typename MatsT, typename IntsT>
  void SingleSlater<MatsT,IntsT>::computeSpinAndAngularProperties() {

    computeSpinAndAngularProperties(nullptr, 0ul);

  };

  /*
  This function computes the expectation value of L^2, Lz, Lx, Ly, J^2, Jz, Jx, and Jy.
  1RDM is assumed to be D(p, q) = < a_p^dagger a_q > in the AO basis taken from singleslater::onePDM.
  2RDM is assumed to be 2RDM(p, q, r, s) = < a_p^dagger a_q^dagger a_s a_r > in the MO basis. 
  Note: If 2RDM not provided, the two-body contribution will be approximated using the 1RDM (which is exact for a single determinant but approximate for correlated wavefunctions).
  */
  template <typename MatsT, typename IntsT>
  void SingleSlater<MatsT,IntsT>::computeSpinAndAngularProperties(const InCore4indexTPI<MatsT>* twoRDM) {

    computeSpinAndAngularProperties(twoRDM, 0ul);

  }

  template <typename MatsT, typename IntsT>
  void SingleSlater<MatsT,IntsT>::computeSpinAndAngularProperties(const InCore4indexTPI<MatsT>* twoRDM, size_t activeEndOff) {

    ROOT_ONLY(comm);

    bool has2RDM = (twoRDM != nullptr);
    
    auto* fb = this->fockBuilder.get();
    if (auto* neofb = dynamic_cast<NEOFockBuilder<MatsT, IntsT>*>(fb)) {
      fb = neofb->getNonNEOUpstream();
    }
    const bool isROHF = (dynamic_cast<ROFock<MatsT, IntsT>*>(fb) != nullptr);

    // Make tiny values to zero to avoid numerical noise in the computed properties like making things unphysically negative. 
    auto clampTiny = [](auto &value) {
      using ValueT = std::decay_t<decltype(value)>;
      if (std::abs(value) < 1e-12) value = ValueT(0);
    };

    if (this->aoints_->angmom == nullptr)
      return;

    if (this->nC == 4)
      CErr("4C spin and angular momentum properties have not been implemented yet.", std::cout);

    
    auto trace_matrix = [](const auto &mat) -> dcomplex {
      // Fast real trace helper used by the expectation-value contractions below.
      dcomplex value = 0.0;
      for (size_t i = 0; i < mat.nRows(); ++i)
        value += mat(i, i);
      return value;
    };

    auto matrix_product = [](const cqmatrix::Matrix<dcomplex> &left,
                             const cqmatrix::Matrix<dcomplex> &right) {
      if (left.nColumns() != right.nRows())
        CErr("computeSpinAndAngularProperties: matrix_product dimension mismatch");
      // Dense matrix multiply in column-major storage.
      cqmatrix::Matrix<dcomplex> result(left.nRows(), right.nColumns());
      blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans,
          (int)left.nRows(), (int)right.nColumns(), (int)left.nColumns(),
          dcomplex(1.0, 0.0), left.pointer(), (int)left.nRows(), right.pointer(),
          (int)right.nRows(), dcomplex(0.0, 0.0), result.pointer(), (int)left.nRows());
      return result;
    };

    auto matrix_product_trans_first = [](const cqmatrix::Matrix<dcomplex> &left,
                             const cqmatrix::Matrix<dcomplex> &right) {
      if (left.nRows() != right.nRows())
        CErr("computeSpinAndAngularProperties: matrix_product_trans_first dimension mismatch");
      // Dense matrix multiply in column-major storage.
      cqmatrix::Matrix<dcomplex> result(left.nColumns(), right.nColumns());
      blas::gemm(blas::Layout::ColMajor, blas::Op::ConjTrans, blas::Op::NoTrans,
          (int)left.nColumns(), (int)right.nColumns(), (int)left.nRows(),
          dcomplex(1.0, 0.0), left.pointer(), (int)left.nRows(), right.pointer(),
          (int)right.nRows(), dcomplex(0.0, 0.0), result.pointer(), (int)result.nRows());
      return result;
    };


    auto pair_expectation = [&](const cqmatrix::Matrix<dcomplex> &density,
                                const cqmatrix::Matrix<dcomplex> &op1,
                                const cqmatrix::Matrix<dcomplex> &op2) -> double {

      // General product expectation for two one-body operators.
      auto op12 = matrix_product(op1, op2);
      dcomplex one_body = trace_matrix(matrix_product_trans_first(density, op12));

      dcomplex two_body = 0.0;
      const size_t NB = density.nRows();

      if (has2RDM) {
        const size_t rdmDim = twoRDM->nBasis();
        const size_t requestedEnd = (activeEndOff > 0ul) ? activeEndOff : rdmDim;
        const size_t contractionEnd = std::min(requestedEnd, std::min(NB, rdmDim));
        // If the 2RDM is available, use it to compute the two-body contribution.
        for (size_t p = 0; p < contractionEnd; ++p)
          for (size_t q = 0; q < contractionEnd; ++q)
            for (size_t r = 0; r < contractionEnd; ++r)
              for (size_t s = 0; s < contractionEnd; ++s)
                two_body += op1(p, q) * op2(r, s) * (*twoRDM)(p, q, r, s);
        // std::cout << std::scientific << std::setprecision(10)
        //           << "  twoRDM version\n"   
        //           << "    One-body contribution = " << std::real(one_body) << '\n'             
        //           << "    Two-body contribution = " << std::real(two_body) << '\n'
        //           << "    Component total       = " << std::real(one_body - two_body)
        //           << std::endl;
      }else{      
        // Faster case for a single determinant:
        // Do1(s,q) = sum_p D(p, s) * op1(p, q) = matrix_product(D, op1)
        auto Do1 =  matrix_product_trans_first(density, op1); //(D)^T * op1
        // Do2(q,s) = sum_r D(r, q) * op2(r, s) = matrix_product(D, op2)
        auto Do2 = matrix_product_trans_first(density, op2); // (D)^T * op2
        // one-body term: sum_pq Do1(s,q) * op2(q, s) = trace(Do1 * op2)
        auto one_body_fast = trace_matrix(matrix_product(Do1, op2));
        
        // two-body term: sum_pqrs D(p, s) * op1(p, q) * D(r, q) * op2(r, s) - D(p, s) * op1(p, q) * D(r, s) * op2(r, q)
        // = trace(Do1 * Do2) - trace(Do1) * trace(Do2)
        auto two_body_term1 = trace_matrix(matrix_product(Do1,Do2));
        auto two_body_term2 = trace_matrix(Do1) * trace_matrix(Do2);
        auto two_body_fast = two_body_term1 - two_body_term2;

        // std::cout << std::scientific << std::setprecision(10)
        //           << "  O^2 (eq.143) component fast debug\n"
        //           << "    One-body contribution = " << std::real(one_body_fast) << '\n'                
        //           << "    Two-body contribution = " << std::real(two_body_fast) << '\n'
        //           << "    Component total       = " << std::real(one_body_fast - two_body_fast)
        //           << std::endl;
        two_body = two_body_fast;
        two_body *= -1.0;
      }    

      return std::real(one_body + two_body);
    };

    auto AO2MO = [&](const cqmatrix::Matrix<dcomplex> &AO,
                    const cqmatrix::Matrix<dcomplex> &coefficients) {
      
      // one and two component
      if (this->nC == 1 or this->nC == 2) 
      {
      const size_t NB = AO.nRows();
      cqmatrix::Matrix<dcomplex> tmp(NB);
      cqmatrix::Matrix<dcomplex> MO(NB);

      // MO = C^dagger * AO * C
      blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans,
                 NB, NB, NB, dcomplex(1.0, 0.0), AO.pointer(), NB,
                 coefficients.pointer(), NB, dcomplex(0.0, 0.0), tmp.pointer(), NB);
      blas::gemm(blas::Layout::ColMajor, blas::Op::ConjTrans, blas::Op::NoTrans,
                 NB, NB, NB, dcomplex(1.0, 0.0), coefficients.pointer(), NB,
                 tmp.pointer(), NB, dcomplex(0.0, 0.0), MO.pointer(), NB);
        return MO;
      }
      else 
      {

        const size_t nAO = coefficients.nRows();
        cqmatrix::Matrix<dcomplex> AO_work = [&]() {
          if (AO.nRows() == nAO && AO.nColumns() == nAO)
            return cqmatrix::Matrix<dcomplex>(AO);

          if (AO.nRows() * 2 == nAO && AO.nColumns() * 2 == nAO) {
            cqmatrix::Matrix<dcomplex> LS(AO.nRows());
            cqmatrix::Matrix<dcomplex> SL(AO.nRows());
            cqmatrix::Matrix<dcomplex> SS(AO.nRows());
            LS.clear();
            SL.clear();
            SS.clear();

            cqmatrix::Matrix<dcomplex> promoted(nAO);
            promoted.componentGather(AO, LS, SL, SS, false);
            return promoted;
          }

          CErr("computeSpinAndAngularProperties: AO2MO incompatible AO/operator and MO-coefficient dimensions");
          return cqmatrix::Matrix<dcomplex>(AO);
        }();

        // 4C paths may provide operators in a 2-block form while MO coefficients
        // are in a 4-block component basis. Promote to component layout if needed.
        cqmatrix::Matrix<dcomplex> tmp(nAO, coefficients.nColumns());
        cqmatrix::Matrix<dcomplex> MO(coefficients.nColumns());

        // MO = C^dagger * AO * C
        blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans,
                  (int)nAO, (int)coefficients.nColumns(), (int)nAO,
                  dcomplex(1.0, 0.0), AO_work.pointer(), (int)nAO,
                  coefficients.pointer(), (int)coefficients.nRows(),
                  dcomplex(0.0, 0.0), tmp.pointer(), (int)nAO);
        blas::gemm(blas::Layout::ColMajor, blas::Op::ConjTrans, blas::Op::NoTrans,
                  (int)coefficients.nColumns(), (int)coefficients.nColumns(), (int)nAO,
                  dcomplex(1.0, 0.0), coefficients.pointer(), (int)coefficients.nRows(),
                  tmp.pointer(), (int)nAO, dcomplex(0.0, 0.0), MO.pointer(),
                  (int)coefficients.nColumns());
        return MO;
      }
    };

    auto mo_coefficients = [&]() {
      if (this->nC == 1 and not this->iCS and not isROHF) {
        const size_t NB = this->moCoefficients[0].get().nRows();
        cqmatrix::Matrix<dcomplex> coeffs(2 * NB);
        coeffs.clear();

        const auto alpha = cqmatrix::Matrix<dcomplex>(this->moCoefficients[0].get());
        const auto beta = cqmatrix::Matrix<dcomplex>(this->moCoefficients[1].get());

        SetMat('N', NB, NB, dcomplex(1.0, 0.0), alpha.pointer(), NB,
               coeffs.pointer(), 2 * NB);
        SetMat('N', NB, NB, dcomplex(1.0, 0.0), beta.pointer(), NB,
               coeffs.pointer() + NB + NB * (2 * NB), 2 * NB);
        return coeffs;
      }

      if (this->nC == 1) {
        return cqmatrix::Matrix<dcomplex>(this->moCoefficients[0].get())
            .template spatialToSpinBlock<dcomplex>();
      }

      if (this->nC == 4) {
        // Keep coefficients in the native 4C AO metric here.
        // The AO->MO density transformation below applies the proper 4C metric.
        return cqmatrix::Matrix<dcomplex>(this->moCoefficients[0].get());
      }

      return cqmatrix::Matrix<dcomplex>(this->moCoefficients[0].get());
    }();

    if (not has2RDM )
       std::cout << "Warning: 2RDM not provided, using 1RDM to approximate two-body contributions to spin and angular momentum properties. This is exact for a single determinant but approximate for correlated wavefunctions." << std::endl;
    
    // Assemble a single spin-orbital MO density matrix.
    const size_t nMO = mo_coefficients.nRows();
    cqmatrix::Matrix<dcomplex> mo_density(nMO);
    mo_density.clear();

    if (this->nC == 1) {
      // Use the existing overlap-aware AO->MO density utility for 1C cases.
      auto ao_blocks = this->onePDM->template spinGatherToBlocks<MatsT>(true, true);
      // ao_blocks: [AA, AB, BA, BB]
      const auto moDenA = this->generateMODensity(
          ao_blocks[0], this->moCoefficients[0].get());

        const auto moDenB = this->generateMODensity(
          ao_blocks.back(),
          ((this->iCS or isROHF) ? this->moCoefficients[0].get()
                     : this->moCoefficients[1].get()));

      const size_t N = moDenA.nRows();
      for (size_t i = 0; i < N; ++i)
        for (size_t j = 0; j < N; ++j) {
          mo_density(i, j) = dcomplex(moDenA(i, j));
          mo_density(N + i, N + j) = dcomplex(moDenB(i, j));
        }
    } else if (this->nC == 2) {
      // For 2C, reuse the existing overlap-aware density utility.
      const auto moDen = this->generateMODensity(this->onePDM->template spinGather<MatsT>(), this->moCoefficients[0].get());
      for (size_t i = 0; i < moDen.nRows(); ++i)
        for (size_t j = 0; j < moDen.nColumns(); ++j)
          mo_density(i, j) = dcomplex(moDen(i, j));
    } else {
      // For 4C, use the same overlap-aware AO->MO density utility as 2C.
      // This preserves spinor-frame mixing present in the converged 1RDM.
      const auto moDen = this->generateMODensity(this->onePDM->template spinGather<MatsT>(), this->moCoefficients[0].get());
      for (size_t i = 0; i < moDen.nRows(); ++i)
        for (size_t j = 0; j < moDen.nColumns(); ++j)
          mo_density(i, j) = dcomplex(moDen(i, j));

      #ifdef _DEBUG_4C_ANGULAR
            std::cout << std::scientific << std::setprecision(16);
            std::cout << "  4C angular diagnostics" << std::endl;
            std::cout << "    nMO=" << nMO << " nO=" << static_cast<size_t>(this->nO) << std::endl;
            std::cout << "    Tr(mo_density) = " << trace_matrix(mo_density) << std::endl;
            std::cout << "    mo_density diagonal window:" << std::endl;
            const size_t windowEnd = std::min<size_t>(nMO, 8);
            for (size_t i = 0; i < windowEnd; ++i)
              std::cout << "      D(" << i << "," << i << ") = " << mo_density(i, i) << std::endl;
      #endif
    } 

    // Default decomposition values (non-4C path stays zeroed).
    this->SExpectLL = {0., 0., 0.};
    this->SExpectSS = {0., 0., 0.};
    this->SSqLL = 0.0;
    this->SSqSS = 0.0;
    this->SSqCross = 0.0;

    // Handle 4-component relativistic case early
    if (this->nC == 4 and this->aoints_->S_4C != nullptr) {
      if (this->aoints_->S_4C->size() < 3 || this->aoints_->S_4C_LL->size() < 3 ||
          this->aoints_->S_4C_SS->size() < 3)
        CErr("4C Spin operators incomplete - expected 3 Cartesian components");

      auto transformOneBodyLikeOtherComponents = [&](const cqmatrix::Matrix<dcomplex> &aoOp) {
        if constexpr (std::is_same_v<MatsT, dcomplex>) {
          MOIntsTransformer<MatsT, IntsT> moInts(*this);
          const std::vector<std::pair<size_t, size_t>> off_sizes = {
              {0, mo_coefficients.nRows()}, {0, mo_coefficients.nRows()}};
          cqmatrix::Matrix<MatsT> moOpNative(mo_coefficients.nRows());
          cqmatrix::Matrix<MatsT> aoOpNativeMatrix(aoOp);
          OnePInts<MatsT> aoOpNative(aoOpNativeMatrix);
          moInts.subsetTransformOPI(off_sizes, aoOpNative, moOpNative.pointer(), false);
          return cqmatrix::Matrix<dcomplex>(moOpNative);
        } else {
          return AO2MO(aoOp, mo_coefficients);
        }
      };

      const size_t nAO4C = mo_coefficients.nRows();
      if (nAO4C % 4 != 0)
        CErr("4C spin operator build expects AO dimension divisible by 4");

      const size_t nLarge = nAO4C / 4;
      const auto overlapAO = cqmatrix::Matrix<dcomplex>(this->aoints_->overlap->matrix());

      // Build LL with the same overlap-based spin construction as <=2C,
      // embedded in the 4C large-large (blocks 0 and 2) subspace.
      cqmatrix::Matrix<dcomplex> Sx_LL_AO(nAO4C);
      cqmatrix::Matrix<dcomplex> Sy_LL_AO(nAO4C);
      cqmatrix::Matrix<dcomplex> Sz_LL_AO(nAO4C);
      Sx_LL_AO.clear();
      Sy_LL_AO.clear();
      Sz_LL_AO.clear();

      for (size_t mu = 0; mu < nLarge; ++mu)
        for (size_t nu = 0; nu < nLarge; ++nu) {
          const auto s = overlapAO(mu, nu);
          Sx_LL_AO(mu, 2 * nLarge + nu) = dcomplex(0.5, 0.0) * s;
          Sx_LL_AO(2 * nLarge + mu, nu) = dcomplex(0.5, 0.0) * s;

          Sy_LL_AO(mu, 2 * nLarge + nu) = dcomplex(0.0, -0.5) * s;
          Sy_LL_AO(2 * nLarge + mu, nu) = dcomplex(0.0, 0.5) * s;

          Sz_LL_AO(mu, nu) = dcomplex(0.5, 0.0) * s;
          Sz_LL_AO(2 * nLarge + mu, 2 * nLarge + nu) = dcomplex(-0.5, 0.0) * s;
        }

      auto Sx_LL_MO = transformOneBodyLikeOtherComponents(Sx_LL_AO);
      auto Sy_LL_MO = transformOneBodyLikeOtherComponents(Sy_LL_AO);
      auto Sz_LL_MO = transformOneBodyLikeOtherComponents(Sz_LL_AO);

      auto Sx_SS_MO = transformOneBodyLikeOtherComponents((*this->aoints_->S_4C_SS)[0]);
      auto Sy_SS_MO = transformOneBodyLikeOtherComponents((*this->aoints_->S_4C_SS)[1]);
      auto Sz_SS_MO = transformOneBodyLikeOtherComponents((*this->aoints_->S_4C_SS)[2]);

      auto Sx_MO = Sx_LL_MO;
      auto Sy_MO = Sy_LL_MO;
      auto Sz_MO = Sz_LL_MO;
      Sx_MO += Sx_SS_MO;
      Sy_MO += Sy_SS_MO;
      Sz_MO += Sz_SS_MO;

      auto Sx_AO = Sx_LL_AO;
      auto Sy_AO = Sy_LL_AO;
      auto Sz_AO = Sz_LL_AO;
      Sx_AO += (*this->aoints_->S_4C_SS)[0];
      Sy_AO += (*this->aoints_->S_4C_SS)[1];
      Sz_AO += (*this->aoints_->S_4C_SS)[2];

      #ifdef _DEBUG_4C_ANGULAR
            cqmatrix::Matrix<dcomplex> ao_density_check(nMO);
            ao_density_check.clear();
            cqmatrix::Matrix<dcomplex> ao_density_tmp(nMO);
            blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans,
                      nMO, nMO, nMO, dcomplex(1.0, 0.0), mo_coefficients.pointer(), nMO,
                      mo_density.pointer(), nMO, dcomplex(0.0, 0.0), ao_density_tmp.pointer(), nMO);
            blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::ConjTrans,
                      nMO, nMO, nMO, dcomplex(1.0, 0.0), ao_density_tmp.pointer(), nMO,
                      mo_coefficients.pointer(), nMO, dcomplex(0.0, 0.0), ao_density_check.pointer(), nMO);

            auto printTracePair = [&](const char *label, const auto &aoOp, const auto &moOp) {
              const auto aoTrace = std::real(trace_matrix(matrix_product(ao_density_check, aoOp)));
              const auto moTrace = std::real(trace_matrix(matrix_product(mo_density, moOp)));
              std::cout << "    " << label << " AO trace = " << aoTrace
                        << "   MO trace = " << moTrace
                        << "   delta = " << (aoTrace - moTrace) << std::endl;
            };

            std::cout << "    4C angular operator trace checks" << std::endl;
            printTracePair("Sx", Sx_AO, Sx_MO);
            printTracePair("Sy", Sy_AO, Sy_MO);
            printTracePair("Sz", Sz_AO, Sz_MO);
            printTracePair("Sx_LL", Sx_LL_AO, Sx_LL_MO);
            printTracePair("Sy_LL", Sy_LL_AO, Sy_LL_MO);
            printTracePair("Sz_LL", Sz_LL_AO, Sz_LL_MO);
            printTracePair("Sx_SS", (*this->aoints_->S_4C_SS)[0], Sx_SS_MO);
            printTracePair("Sy_SS", (*this->aoints_->S_4C_SS)[1], Sy_SS_MO);
            printTracePair("Sz_SS", (*this->aoints_->S_4C_SS)[2], Sz_SS_MO);

                  // Check the effective 4C metric and LL overlap in the MO occupied block.
                  const size_t n4 = nMO / 4;
                  cqmatrix::Matrix<dcomplex> metric4C(nMO);
                  metric4C.clear();
                  SetMat('N', n4, n4, dcomplex(1.0, 0.0), this->aoints_->overlap->matrix().pointer(), n4,
                    metric4C.pointer(), nMO);
                  SetMat('N', n4, n4, dcomplex(1.0 / (2.0 * SpeedOfLight * SpeedOfLight), 0.0), this->aoints_->kinetic->matrix().pointer(), n4,
                    metric4C.pointer() + nMO * nMO / 4 + nMO / 4, nMO);
                  SetMat('N', n4, n4, dcomplex(1.0, 0.0), this->aoints_->overlap->matrix().pointer(), n4,
                    metric4C.pointer() + nMO * nMO / 2 + nMO / 2, nMO);
                  SetMat('N', n4, n4, dcomplex(1.0 / (2.0 * SpeedOfLight * SpeedOfLight), 0.0), this->aoints_->kinetic->matrix().pointer(), n4,
                    metric4C.pointer() + 3 * nMO * nMO / 4 + 3 * nMO / 4, nMO);

                  cqmatrix::Matrix<dcomplex> llMetric4C(nMO);
                  llMetric4C.clear();
                  SetMat('N', n4, n4, dcomplex(1.0, 0.0), this->aoints_->overlap->matrix().pointer(), n4,
                    llMetric4C.pointer(), nMO);
                  SetMat('N', n4, n4, dcomplex(1.0, 0.0), this->aoints_->overlap->matrix().pointer(), n4,
                    llMetric4C.pointer() + nMO * nMO / 2 + nMO / 2, nMO);

                  auto metric4C_MO = transformOneBodyLikeOtherComponents(metric4C);
                  auto llMetric4C_MO = transformOneBodyLikeOtherComponents(llMetric4C);
                  auto SxSq_MO = matrix_product(Sx_MO, Sx_MO);
                  auto SxLLSq_MO = matrix_product(Sx_LL_MO, Sx_LL_MO);

                  std::cout << "    4C metric/operator consistency checks" << std::endl;
                  std::cout << "      Tr(D * M4C_MO)      = " << trace_matrix(matrix_product(mo_density, metric4C_MO)) << std::endl;
                  std::cout << "      Tr(D * MLL_MO)      = " << trace_matrix(matrix_product(mo_density, llMetric4C_MO)) << std::endl;
                  std::cout << "      Tr(D * Sx^2)        = " << trace_matrix(matrix_product(mo_density, SxSq_MO)) << std::endl;
                  std::cout << "      Tr(D * Sx_LL^2)     = " << trace_matrix(matrix_product(mo_density, SxLLSq_MO)) << std::endl;
      #endif

      // Evaluate total 4C spin expectations from explicit spin operators.
            this->SExpect[0] = std::real(trace_matrix(matrix_product(mo_density, Sx_MO)));
            this->SExpect[1] = std::real(trace_matrix(matrix_product(mo_density, Sy_MO)));
            this->SExpect[2] = std::real(trace_matrix(matrix_product(mo_density, Sz_MO)));
            this->SSq = pair_expectation(mo_density, Sx_MO, Sx_MO)
              + pair_expectation(mo_density, Sy_MO, Sy_MO)
              + pair_expectation(mo_density, Sz_MO, Sz_MO);

            this->SExpectLL[0] = std::real(trace_matrix(matrix_product(mo_density, Sx_LL_MO)));
            this->SExpectLL[1] = std::real(trace_matrix(matrix_product(mo_density, Sy_LL_MO)));
            this->SExpectLL[2] = std::real(trace_matrix(matrix_product(mo_density, Sz_LL_MO)));

            this->SExpectSS[0] = std::real(trace_matrix(matrix_product(mo_density, Sx_SS_MO)));
            this->SExpectSS[1] = std::real(trace_matrix(matrix_product(mo_density, Sy_SS_MO)));
            this->SExpectSS[2] = std::real(trace_matrix(matrix_product(mo_density, Sz_SS_MO)));

            this->SSqLL = pair_expectation(mo_density, Sx_LL_MO, Sx_LL_MO)
                   + pair_expectation(mo_density, Sy_LL_MO, Sy_LL_MO)
                   + pair_expectation(mo_density, Sz_LL_MO, Sz_LL_MO);
            this->SSqSS = pair_expectation(mo_density, Sx_SS_MO, Sx_SS_MO)
                   + pair_expectation(mo_density, Sy_SS_MO, Sy_SS_MO)
                   + pair_expectation(mo_density, Sz_SS_MO, Sz_SS_MO);
      this->SSqCross = this->SSq - this->SSqLL - this->SSqSS;

      // Handle orbital angular momentum
      if (this->aoints_->L_pauli != nullptr) {
        if (this->aoints_->L_pauli->size() < 3)
          CErr("4C Orbital operators incomplete - expected 3 components (x,y,z)");

        // Use L_pauli which has the LL (spatial angular momentum) in the S() component.
        // This properly handles the 4C block structure via spinGather.
        auto gather_spinor = [](const cqmatrix::PauliSpinorMatrices<IntsT> &P) {
          return P.template spinGather<dcomplex>();
        };

        auto Lx_MO = transformOneBodyLikeOtherComponents(gather_spinor((*this->aoints_->L_pauli)[0]));
        auto Ly_MO = transformOneBodyLikeOtherComponents(gather_spinor((*this->aoints_->L_pauli)[1]));
        auto Lz_MO = transformOneBodyLikeOtherComponents(gather_spinor((*this->aoints_->L_pauli)[2]));

        // Compute orbital angular momentum expectations
        this->LExpect[0] = std::real(trace_matrix(matrix_product(mo_density, Lx_MO)));
        this->LExpect[1] = std::real(trace_matrix(matrix_product(mo_density, Ly_MO)));
        this->LExpect[2] = std::real(trace_matrix(matrix_product(mo_density, Lz_MO)));

        // <L^2> = <Lx^2> + <Ly^2> + <Lz^2>
        this->LSq = pair_expectation(mo_density, Lx_MO, Lx_MO)
             + pair_expectation(mo_density, Ly_MO, Ly_MO)
             + pair_expectation(mo_density, Lz_MO, Lz_MO);

          // Use the symmetrized spin-orbit cross term to reduce numerical
          // non-commutativity artifacts after AO->MO transformations.
          const auto SL = pair_expectation(mo_density, Sx_MO, Lx_MO)
            + pair_expectation(mo_density, Sy_MO, Ly_MO)
            + pair_expectation(mo_density, Sz_MO, Lz_MO);
          const auto LS = pair_expectation(mo_density, Lx_MO, Sx_MO)
            + pair_expectation(mo_density, Ly_MO, Sy_MO)
            + pair_expectation(mo_density, Lz_MO, Sz_MO);
          this->SL = SL;
          this->LS = LS;
        // <J> = <L> + <S>, computed via operators
        auto Jx_MO = Lx_MO;
        auto Jy_MO = Ly_MO;
        auto Jz_MO = Lz_MO;
        Jx_MO += Sx_MO;
        Jy_MO += Sy_MO;
        Jz_MO += Sz_MO;

        this->JExpect[0] = std::real(trace_matrix(matrix_product(mo_density, Jx_MO)));
        this->JExpect[1] = std::real(trace_matrix(matrix_product(mo_density, Jy_MO)));
        this->JExpect[2] = std::real(trace_matrix(matrix_product(mo_density, Jz_MO)));

        // <J^2> = <L^2> + <S^2> + <S.L> + <L.S>
        this->JSq = this->LSq + this->SSq + this->SL + this->LS;
      } else {
        // Orbital operators not yet implemented for 4C
        std::cout << "  WARNING: 4C orbital angular momentum operators not yet implemented\n";
        this->LExpect = {0., 0., 0.};
        this->LSq = 0.;
        this->SL = 0.;
        this->LS = 0.;
        this->JExpect = {this->SExpect[0], this->SExpect[1], this->SExpect[2]};
        this->JSq = this->SSq;
      }

      // Clamp tiny values to zero
      for (auto &component : this->LExpect) clampTiny(component);
      for (auto &component : this->SExpect) clampTiny(component);
      for (auto &component : this->JExpect) clampTiny(component);
      for (auto &component : this->SExpectLL) clampTiny(component);
      for (auto &component : this->SExpectSS) clampTiny(component);
      clampTiny(this->LSq);
      clampTiny(this->SSq);
      clampTiny(this->SSqLL);
      clampTiny(this->SSqSS);
      clampTiny(this->SSqCross);
      clampTiny(this->SL);
      clampTiny(this->LS);
      clampTiny(this->JSq);
      
      return;
    }


    cqmatrix::Matrix<dcomplex> LxAO((*this->aoints_->angmom)["X"]->matrix());
    cqmatrix::Matrix<dcomplex> LyAO((*this->aoints_->angmom)["Y"]->matrix());
    cqmatrix::Matrix<dcomplex> LzAO((*this->aoints_->angmom)["Z"]->matrix());

    cqmatrix::Matrix<dcomplex> Lx(mo_coefficients.nRows());
    cqmatrix::Matrix<dcomplex> Ly(mo_coefficients.nRows());
    cqmatrix::Matrix<dcomplex> Lz(mo_coefficients.nRows());

    if (this->nC == 1) {
      // MOIntsTransformer uses only alpha MOs for 1C, so keep explicit spin-block path here.
      Lx = AO2MO(LxAO.template spatialToSpinBlock<dcomplex>(), mo_coefficients);
      Ly = AO2MO(LyAO.template spatialToSpinBlock<dcomplex>(), mo_coefficients);
      Lz = AO2MO(LzAO.template spatialToSpinBlock<dcomplex>(), mo_coefficients);
    } else {
      MOIntsTransformer<MatsT, IntsT> moInts(*this);
      const std::vector<std::pair<size_t, size_t>> off_sizes = {
          {0, mo_coefficients.nRows()}, {0, mo_coefficients.nRows()}};

      cqmatrix::Matrix<MatsT> LxMO_native(mo_coefficients.nRows());
      cqmatrix::Matrix<MatsT> LyMO_native(mo_coefficients.nRows());
      cqmatrix::Matrix<MatsT> LzMO_native(mo_coefficients.nRows());

      OnePInts<MatsT> LxAO_native(cqmatrix::Matrix<MatsT>((*this->aoints_->angmom)["X"]->matrix())
                                      .template spatialToSpinBlock<MatsT>());
      OnePInts<MatsT> LyAO_native(cqmatrix::Matrix<MatsT>((*this->aoints_->angmom)["Y"]->matrix())
                                      .template spatialToSpinBlock<MatsT>());
      OnePInts<MatsT> LzAO_native(cqmatrix::Matrix<MatsT>((*this->aoints_->angmom)["Z"]->matrix())
                                      .template spatialToSpinBlock<MatsT>());

      moInts.subsetTransformOPI(off_sizes, LxAO_native, LxMO_native.pointer(), false);
      moInts.subsetTransformOPI(off_sizes, LyAO_native, LyMO_native.pointer(), false);
      moInts.subsetTransformOPI(off_sizes, LzAO_native, LzMO_native.pointer(), false);

      Lx = cqmatrix::Matrix<dcomplex>(LxMO_native);
      Ly = cqmatrix::Matrix<dcomplex>(LyMO_native);
      Lz = cqmatrix::Matrix<dcomplex>(LzMO_native);
    }


    Lx *= dcomplex(0.0, 1.0);
    Ly *= dcomplex(0.0, 1.0);
    Lz *= dcomplex(0.0, 1.0);

    // Build spin operators from AO overlap blocks, then transform to MO.
    const size_t NB = mo_density.nRows() / 2;
    const auto overlapAO = cqmatrix::Matrix<dcomplex>(this->aoints_->overlap->matrix());
    cqmatrix::Matrix<dcomplex> SxAO(2 * NB);
    cqmatrix::Matrix<dcomplex> SyAO(2 * NB);
    cqmatrix::Matrix<dcomplex> SzAO(2 * NB);
    SxAO.clear();
    SyAO.clear();
    SzAO.clear();
    for (size_t mu = 0; mu < NB; ++mu)
      for (size_t nu = 0; nu < NB; ++nu) {
        const auto s = overlapAO(mu, nu);
        //Setting up the spin operators in the AO basis. The spin operators are defined as follows:
        // Sx is Sx = 0.5 * (Sab + Sba)
        SxAO(mu, NB + nu) = dcomplex(0.5, 0.0) * s; // alpha-beta
        SxAO(NB + mu, nu) = dcomplex(0.5, 0.0) * s; // beta-alpha
        // Sy is Sy = 0.5 * (S ab - S ba) * i
        SyAO(mu, NB + nu) = dcomplex(0.0, -0.5) * s; // alpha-beta
        SyAO(NB + mu, nu) = dcomplex(0.0, 0.5) * s; // beta-alpha
        // Sz is Sz = 0.5 * (Saa - Sbb)
        SzAO(mu, nu) = dcomplex(0.5, 0.0) * s; // alpha-alpha
        SzAO(NB + mu, NB + nu) = dcomplex(-0.5, 0.0) * s; // beta-beta
      }

    // Transform/make spin operators to MO basis
    auto Sx = AO2MO(SxAO, mo_coefficients); 
    auto Sy = AO2MO(SyAO, mo_coefficients);
    auto Sz = AO2MO(SzAO, mo_coefficients);

    auto Jx = Lx + Sx;
    auto Jy = Ly + Sy;
    auto Jz = Lz + Sz;


    //
    // Expectation of the above operators
    //

    // <Sx>, <Sy>, <Sz>
    // std::cout << "  Spin expectation values:" << std::endl;
    auto SxExp = std::real(trace_matrix(matrix_product(mo_density, Sx)));
    auto SyExp = std::real(trace_matrix(matrix_product(mo_density, Sy)));
    auto SzExp = std::real(trace_matrix(matrix_product(mo_density, Sz)));
    
    // <S^2> = <Sx^2> + <Sy^2> + <Sz^2>
    auto SSqExp = pair_expectation(mo_density, Sx, Sx) + pair_expectation(mo_density, Sy, Sy) + pair_expectation(mo_density, Sz, Sz);

    // Compute <S^2> Eigenvalue is S(S+1) or S = (-1 + sqrt(1 + 4<S^2>))/2
    auto SQuantNum = (-1.0 + std::sqrt(1.0 + 4.0 * SSqExp)) / 2.0;

    this->SExpect[0] = SxExp;
    this->SExpect[1] = SyExp;
    this->SExpect[2] = SzExp;
    this->SSq = SSqExp;
    this->SQuantNum = SQuantNum;


    // <Lx>, <Ly>, <Lz>
    // std::cout << "  Orbital angular momentum expectation values:" << std::endl;
    this->LExpect[0] = std::real(trace_matrix(matrix_product(mo_density, Lx)));
    this->LExpect[1] = std::real(trace_matrix(matrix_product(mo_density, Ly)));
    this->LExpect[2] = std::real(trace_matrix(matrix_product(mo_density, Lz)));

    // <L^2> = <Lx^2> + <Ly^2> + <Lz^2>
    this->LSq = pair_expectation(mo_density, Lx, Lx)
           + pair_expectation(mo_density, Ly, Ly)
           + pair_expectation(mo_density, Lz, Lz);

    // Compute <L^2> Eigenvalue is L(L+1) or L = (-1 + sqrt(1 + 4<L^2>))/2
    auto LQuantNum = (-1.0 + std::sqrt(1.0 + 4.0 * this->LSq)) / 2.0;
    this->LQuantNum = LQuantNum;
    
    // Use the SL and LS cross terms just in case non-commutativity
    // std::cout << "  Spin-orbit cross term expectation values:" << std::endl; 
    const auto SL = pair_expectation(mo_density, Sx, Lx)
      + pair_expectation(mo_density, Sy, Ly)
      + pair_expectation(mo_density, Sz, Lz);
    const auto LS = pair_expectation(mo_density, Lx, Sx)
      + pair_expectation(mo_density, Ly, Sy)
      + pair_expectation(mo_density, Lz, Sz);
    this->SL = SL;
    this->LS = LS;

    // <Jx>, <Jy>, <Jz>
    // std::cout << "  Total angular momentum expectation values:" << std::endl;
    this->JExpect[0] = std::real(trace_matrix(matrix_product(mo_density, Jx)));
    this->JExpect[1] = std::real(trace_matrix(matrix_product(mo_density, Jy)));
    this->JExpect[2] = std::real(trace_matrix(matrix_product(mo_density, Jz)));
    
    // <J^2> = <Jx^2> + <Jy^2> + <Jz^2> = <L^2> + <S^2> + <S.L> + <L.S>
    auto temp_JSq  = this->LSq + this->SSq + SL + LS;

    // Also valid approach to J^2 using the J operators directly
    this->JSq = pair_expectation(mo_density, Jx, Jx)
           + pair_expectation(mo_density, Jy, Jy)
           + pair_expectation(mo_density, Jz, Jz);

    // Compute <J^2> Eigenvalue is J(J+1) or J = (-1 + sqrt(1 + 4<J^2>))/2
    auto JQuantNum = (-1.0 + std::sqrt(1.0 + 4.0 * this->JSq)) / 2.0;
    this->JQuantNum = JQuantNum;

    // std::cout << "  J^2 consistency check:  J^2 = " << this->JSq
    //       << "  vs.  L^2 + S^2 + <S.L> + <L.S> = " << temp_JSq
    //       << std::endl;
    
    // Clamp tiny values to zero to clean up numerical noise in the outputs.
    for (auto &component : this->LExpect) clampTiny(component);
    for (auto &component : this->JExpect) clampTiny(component);
    clampTiny(this->LSq);
    clampTiny(this->JSq);
    clampTiny(this->SL);
    clampTiny(this->LS);
    clampTiny(this->SSq);

  };


  template <typename MatsT, typename IntsT>
  typename SingleSlater<MatsT, IntsT>::AngularResults
  SingleSlater<MatsT, IntsT>::computeAngularOrbitalRows(bool includeBeta) {

    AngularResults res;
    if (this->aoints_ == nullptr || this->aoints_->angmom == nullptr) return res;

    auto totalStart = tick();

    // Detect ROHF to avoid duplicating beta orbitals
    auto* fb = this->fockBuilder.get();
    if (auto* neofb = dynamic_cast<NEOFockBuilder<MatsT, IntsT>*>(fb)) {
      fb = neofb->getNonNEOUpstream();
    }
    bool isROHF = (dynamic_cast<ROFock<MatsT, IntsT>*>(fb) != nullptr);

    bool isSpatial = (this->nC == 1);
    auto buildSpinCoefficients = [](const cqmatrix::Matrix<dcomplex> &coeffs,
                                    bool betaBlock) {
      const size_t nBasis = coeffs.nRows();
      const size_t nOrb = coeffs.nColumns();
      cqmatrix::Matrix<dcomplex> spinCoeffs(2 * nBasis, nOrb);
      spinCoeffs.clear();

      const size_t rowOffset = betaBlock ? nBasis : 0;
      for (size_t orb = 0; orb < nOrb; ++orb) {
        const dcomplex *src = coeffs.pointer() + orb * nBasis;
        dcomplex *dst = spinCoeffs.pointer() + orb * (2 * nBasis) + rowOffset;
        for (size_t row = 0; row < nBasis; ++row)
          dst[row] = src[row];
      }

      return spinCoeffs;
    };

    auto matrixProduct = [](const cqmatrix::Matrix<dcomplex> &left,
                            const cqmatrix::Matrix<dcomplex> &right) {
      if (left.nColumns() != right.nRows())
        CErr("computeAngularOrbitalRows: matrixProduct dimension mismatch", std::cout);
      cqmatrix::Matrix<dcomplex> out(left.nRows(), right.nColumns());
      blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans,
          (int)left.nRows(), (int)right.nColumns(), (int)left.nColumns(),
          dcomplex(1.0, 0.0), left.pointer(), (int)left.nRows(),
          right.pointer(), (int)right.nRows(), dcomplex(0.0, 0.0),
          out.pointer(), (int)left.nRows());
      return out;
    };

    auto orbitalExpectations = [&](const cqmatrix::Matrix<dcomplex> &op,
                                   const cqmatrix::Matrix<dcomplex> &coeffs) {
      const size_t nBasis = coeffs.nRows();
      const size_t nOrb = coeffs.nColumns();

      cqmatrix::Matrix<dcomplex> tmp(nBasis, nOrb);
      cqmatrix::Matrix<dcomplex> transformed(nOrb, nOrb);

      blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans,
          (int)nBasis, (int)nOrb, (int)nBasis, dcomplex(1.0, 0.0),
          op.pointer(), (int)nBasis, coeffs.pointer(), (int)nBasis,
          dcomplex(0.0, 0.0), tmp.pointer(), (int)nBasis);

      blas::gemm(blas::Layout::ColMajor, blas::Op::ConjTrans, blas::Op::NoTrans,
          (int)nOrb, (int)nOrb, (int)nBasis, dcomplex(1.0, 0.0),
          coeffs.pointer(), (int)nBasis, tmp.pointer(), (int)nBasis,
          dcomplex(0.0, 0.0), transformed.pointer(), (int)nOrb);

      std::vector<double> values(nOrb);
      for (size_t orb = 0; orb < nOrb; ++orb)
        values[orb] = std::real(transformed(orb, orb));
      return values;
    };

    auto operatorStart = tick();

    cqmatrix::Matrix<MatsT> lxAO((*this->aoints_->angmom)["X"]->matrix());
    cqmatrix::Matrix<MatsT> lyAO((*this->aoints_->angmom)["Y"]->matrix());
    cqmatrix::Matrix<MatsT> lzAO((*this->aoints_->angmom)["Z"]->matrix());

    auto Lx = cqmatrix::Matrix<dcomplex>(this->orthoSpinor->nonortho2ortho(lxAO).template spatialToSpinBlock<MatsT>());
    auto Ly = cqmatrix::Matrix<dcomplex>(this->orthoSpinor->nonortho2ortho(lyAO).template spatialToSpinBlock<MatsT>());
    auto Lz = cqmatrix::Matrix<dcomplex>(this->orthoSpinor->nonortho2ortho(lzAO).template spatialToSpinBlock<MatsT>());

    Lx *= dcomplex(0.0, 1.0);
    Ly *= dcomplex(0.0, 1.0);
    Lz *= dcomplex(0.0, 1.0);

    const size_t NB = Lx.nRows() / 2;
    cqmatrix::Matrix<dcomplex> Sx(2 * NB);
    cqmatrix::Matrix<dcomplex> Sy(2 * NB);
    cqmatrix::Matrix<dcomplex> Sz(2 * NB);
    Sx.clear();
    Sy.clear();
    Sz.clear();

    for (size_t i = 0; i < NB; ++i) {
      Sx(i, NB + i) = dcomplex(0.5, 0.0);
      Sx(NB + i, i) = dcomplex(0.5, 0.0);

      Sy(i, NB + i) = dcomplex(0.0, -0.5);
      Sy(NB + i, i) = dcomplex(0.0, 0.5);

      Sz(i, i) = dcomplex(0.5, 0.0);
      Sz(NB + i, NB + i) = dcomplex(-0.5, 0.0);
    }

    auto Jx = Lx + Sx;
    auto Jy = Ly + Sy;
    auto Jz = Lz + Sz;

    auto L2Op = matrixProduct(Lx, Lx);
    L2Op += matrixProduct(Ly, Ly);
    L2Op += matrixProduct(Lz, Lz);

    auto S2Op = matrixProduct(Sx, Sx);
    S2Op += matrixProduct(Sy, Sy);
    S2Op += matrixProduct(Sz, Sz);

    auto J2Op = matrixProduct(Jx, Jx);
    J2Op += matrixProduct(Jy, Jy);
    J2Op += matrixProduct(Jz, Jz);

    const double operatorBuildTime = tock(operatorStart);

    std::vector<AngularOrbitalRow> rows;
    rows.reserve(this->nOA + this->nVA + (this->iCS ? 0 : this->nOB + this->nVB));

    const double occScale = this->iCS ? 2.0 : 1.0;

    auto collectWithScale = [&](cqmatrix::Matrix<MatsT> coeffs,
                                const double *energies,
                                size_t nOcc,
                                double occValue,
                                bool betaBlock,
                                std::vector<AngularOrbitalRow> &destRows) {
      if (isSpatial) this->orthoSpinor->nonortho2orthoCoeffs(coeffs);
      else this->orthoAB->nonortho2orthoCoeffs(coeffs);

      cqmatrix::Matrix<dcomplex> coeffsC = isSpatial
          ? buildSpinCoefficients(cqmatrix::Matrix<dcomplex>(coeffs), betaBlock)
          : cqmatrix::Matrix<dcomplex>(coeffs);

      const auto lzValues = orbitalExpectations(Lz, coeffsC);
      const auto szValues = orbitalExpectations(Sz, coeffsC);
      const auto jzValues = orbitalExpectations(Jz, coeffsC);
      const auto l2Values = orbitalExpectations(L2Op, coeffsC);
      const auto s2Values = orbitalExpectations(S2Op, coeffsC);
      const auto j2Values = orbitalExpectations(J2Op, coeffsC);

      for (size_t iOrb = 0; iOrb < coeffsC.nColumns(); ++iOrb) {
        const bool occupied = iOrb < nOcc;
        destRows.push_back({energies[iOrb], occupied ? occValue : 0.0,
                            szValues[iOrb], lzValues[iOrb], jzValues[iOrb],
                            s2Values[iOrb], l2Values[iOrb], j2Values[iOrb]});
      }
    };

    auto expectationStart = tick();

    if (isSpatial) {
      collectWithScale(cqmatrix::Matrix<MatsT>(this->mo[0]), this->eps1,
                       this->nOA, occScale, false, rows);

      // For UHF, collect beta orbitals; for ROHF, beta is a copy of alpha so skip it
      if (not this->iCS && not isROHF && includeBeta) {
        collectWithScale(cqmatrix::Matrix<MatsT>(this->mo[1]), this->eps2,
                         this->nOB, 1.0, true, rows);
      }
    } else {
      auto collectSpinor = [&](cqmatrix::Matrix<MatsT> coeffs,
                               const double *energies,
                               size_t nOcc,
                               std::vector<AngularOrbitalRow> &destRows) {
        this->orthoAB->nonortho2orthoCoeffs(coeffs);
        cqmatrix::Matrix<dcomplex> coeffsC(coeffs);

        const auto lzValues = orbitalExpectations(Lz, coeffsC);
        const auto szValues = orbitalExpectations(Sz, coeffsC);
        const auto jzValues = orbitalExpectations(Jz, coeffsC);
        const auto l2Values = orbitalExpectations(L2Op, coeffsC);
        const auto s2Values = orbitalExpectations(S2Op, coeffsC);
        const auto j2Values = orbitalExpectations(J2Op, coeffsC);

        for (size_t iOrb = 0; iOrb < coeffsC.nColumns(); ++iOrb) {
          const bool occupied = iOrb < nOcc;
          destRows.push_back({energies[iOrb], occupied ? 1.0 : 0.0,
                              szValues[iOrb], lzValues[iOrb], jzValues[iOrb],
                              s2Values[iOrb], l2Values[iOrb], j2Values[iOrb]});
        }
      };

      collectSpinor(cqmatrix::Matrix<MatsT>(this->mo[0]), this->eps1, this->nO, rows);
    }

    const double expectationTime = tock(expectationStart);

    std::stable_sort(rows.begin(), rows.end(),
        [](const AngularOrbitalRow &left, const AngularOrbitalRow &right) {
          return left.energy < right.energy;
        });

    res.rows = std::move(rows);
    res.operatorBuildTime = operatorBuildTime;
    res.expectationTime = expectationTime;
    res.totalTime = tock(totalStart);
    return res;
  }
  
template <typename MatsT, typename IntsT>
  void SingleSlater<MatsT,IntsT>::compute4CDipole(EMPerturbation &pert) {
    ROOT_ONLY(comm);

    std::vector<cqmatrix::PauliSpinorMatrices<dcomplex>>
       lenElectric4C = *(this->aoints_->lenElectric->gather4CDipole());

    // Compute elecric contribution to the dipoles
    for(auto iXYZ = 0; iXYZ < 3; iXYZ++) {
      // Scalar
      double dipole_s = -this->template computeOBProperty<DENSITY_TYPE::SCALAR>(
          lenElectric4C[iXYZ].S().pointer());
      // MZ
      double dipole_z = -this->template computeOBProperty<DENSITY_TYPE::MZ>(
          lenElectric4C[iXYZ].Z().pointer());
      // MY
      double dipole_y= -this->template computeOBProperty<DENSITY_TYPE::MY>(
          lenElectric4C[iXYZ].Y().pointer());
      // MX
      double dipole_x= -this->template computeOBProperty<DENSITY_TYPE::MX>(
          lenElectric4C[iXYZ].X().pointer());
    
      this->elecDipole[iXYZ] = dipole_s + dipole_z + dipole_y + dipole_x;
    }

    
    for(auto &atom : this->molecule().atoms){
      if (atom.quantum) continue;  
      MatAdd('N','N',3,1,1.,&this->elecDipole[0],3,atom.nucCharge,
          &atom.coord[0],3,&this->elecDipole[0],3);
    }

  };

template <typename MatsT, typename IntsT>
  void SingleSlater<MatsT,IntsT>::computeFockX2CDipole(EMPerturbation &pert) {
    ROOT_ONLY(comm);

    // Compute elecric contribution to the dipoles
    for(auto iXYZ = 0; iXYZ < 3; iXYZ++) {
      // Scalar
      double dipole_s = -this->template computeOBProperty<DENSITY_TYPE::SCALAR>(
          this->pchgDipole_[iXYZ]->S().pointer());
      // MZ
      double dipole_z = -this->template computeOBProperty<DENSITY_TYPE::MZ>(
          this->pchgDipole_[iXYZ]->Z().pointer());
      // MY
      double dipole_y= -this->template computeOBProperty<DENSITY_TYPE::MY>(
          this->pchgDipole_[iXYZ]->Y().pointer());
      // MX
      double dipole_x= -this->template computeOBProperty<DENSITY_TYPE::MX>(
          this->pchgDipole_[iXYZ]->X().pointer());
    
      this->elecDipole[iXYZ] = dipole_s + dipole_z + dipole_y + dipole_x;
    }

    for(auto &atom : this->molecule().atoms){
      if (atom.quantum) continue;  
      MatAdd('N','N',3,1,1.,&this->elecDipole[0],3,atom.nucCharge,
          &atom.coord[0],3,&this->elecDipole[0],3);
    }

  };

/**
 *  \brief A function to calculate Densitry matrix in MO basis, 
 *  given the density matrix and the coefficient matrix in AO basis
 *
 *  Equation: D^{MO} = C^{T} S D^{AO} S C
 */
template <typename MatsT, typename IntsT>
  cqmatrix::Matrix<MatsT> SingleSlater<MatsT,IntsT>::generateMODensity(const cqmatrix::Matrix<MatsT>& denAO, const cqmatrix::Matrix<MatsT>& coeffAO) {
    
    //ROOT_ONLY(comm);

    size_t NB = coeffAO.nRows();
    cqmatrix::Matrix<MatsT> S(NB);
    
    // Obtaining overlap matrix S
    if(this->nC == 1 ) {
      S = this->aoints_->overlap->matrix();
    } else if(this->nC == 2) {
      std::fill_n(S.pointer(),NB*NB,MatsT(0.0));
      SetMat('N',NB/2,NB/2,MatsT(1.),this->aoints_->overlap->matrix().pointer(), NB/2, S.pointer(),NB);
      SetMat('N',NB/2,NB/2,MatsT(1.),this->aoints_->overlap->matrix().pointer(), NB/2, S.pointer()+NB*NB/2+NB/2,NB);
    } else if(this->nC == 4) {
      std::fill_n(S.pointer(),NB*NB,MatsT(0.0));
      SetMat('N',NB/4,NB/4,MatsT(1.),this->aoints_->overlap->matrix().pointer(), NB/4, S.pointer(),NB);
      SetMat('N',NB/4,NB/4,MatsT(1./(2*SpeedOfLight()*SpeedOfLight())),this->aoints_->kinetic->matrix().pointer(), NB/4, S.pointer()+NB*NB/4+NB/4,NB);
      SetMat('N',NB/4,NB/4,MatsT(1.),this->aoints_->overlap->matrix().pointer(), NB/4, S.pointer()+NB*NB/2+NB/2,NB);
      SetMat('N',NB/4,NB/4,MatsT(1./(2*SpeedOfLight()*SpeedOfLight())),this->aoints_->kinetic->matrix().pointer(), NB/4, S.pointer()+NB*NB*3/4+NB*3/4,NB);
    } else{
      CErr("nC invalid in OrbitalModifierNew<singleSlaterT,MatsT,IntsT>::computeMODensity!");
    }
    
    cqmatrix::Matrix<MatsT> SCR(NB);
    cqmatrix::Matrix<MatsT> SCR1(NB);

    // SCR  = D^{AO} S
    blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans, 
        NB, NB, NB, MatsT(1.), denAO.pointer(), NB,
        S.pointer(), NB, MatsT(0.), SCR.pointer(), NB);
    //  SCR1 = S * SCR = S * D^{AO} S
    blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans, 
        NB, NB, NB, MatsT(1.), S.pointer(), NB,
        SCR.pointer(), NB, MatsT(0.), SCR1.pointer(), NB);
    //  SCR  = C^T * SCR1 = C^T * S * D^{AO} S
    blas::gemm(blas::Layout::ColMajor, blas::Op::ConjTrans, blas::Op::NoTrans, 
        NB, NB, NB, MatsT(1.), coeffAO.pointer(), NB,
        SCR1.pointer(), NB, MatsT(0.), SCR.pointer(), NB);
    //  SCR1 = SCR * C = C^T * S * D^{AO} S C
    blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans, 
        NB, NB, NB, MatsT(1.), SCR.pointer(), NB,
        coeffAO.pointer(), NB, MatsT(0.), SCR1.pointer(), NB);
    
    return SCR1;

  };

template <typename MatsT, typename IntsT>
  void SingleSlater<MatsT,IntsT>::printOrbitalPopulation(std::ostream &out) {

    //ROOT_ONLY(comm);

    // Do Orbital Analysis in AO basis

    std::vector<cqmatrix::Matrix<MatsT>> aoDen;
    //this->ortho2aoDen(this->onePDMSquareOrtho);
    //this->setOnePDMAO(this->onePDMSquareAO.data());

    // Spin-Gather AO Density
    if(this->nC == 1 ){
      aoDen = this->onePDM->template spinGatherToBlocks<MatsT>(false);
    } else {
      aoDen.push_back(this->onePDM->template spinGather<MatsT>());
    }

    // Transform a copy of the ground-state MOs in AO basis
    std::vector<cqmatrix::Matrix<MatsT>> aoMO = this->mo;

    // Transform alpha Density and compute populations
    size_t NB = aoMO[0].nRows();
    std::vector<double> population;
    std::vector<cqmatrix::Matrix<MatsT>> moDen;

    moDen.push_back( this->generateMODensity(aoDen[0], aoMO[0]) );
    for( size_t i=0; i<NB; ++i)
      population.push_back( std::real(moDen[0](i,i)) );

    // UHF Beta populations
    if(this->nC == 1 and not this->iCS ){
      moDen.push_back( this->generateMODensity(aoDen[1], aoMO[1]) );
      for( size_t i=0; i<NB; ++i)
        population.push_back( std::real(moDen[1](i,i)) );
    }

    // Printing
    //if( this->printLevel > 1 ) {
    {
      size_t orbPerRow = 5;
      auto printBlock = [&](std::string header, size_t& start, size_t n){
        std::cout << header << std::endl;
        std::cout << std::fixed << std::setprecision(11);

        for(auto idx = 0; idx < n; idx += orbPerRow) {

          size_t end = idx + orbPerRow < n ? orbPerRow : n - idx;
          for(auto idummy = idx; idummy < idx+end; idummy++) {
            std::cout << std::setw(15) << population[start+idummy];
          }
          std::cout << '\n';
        }
        start += n;
      };


      #if 0
        moDen[0].output(std::cout, "MO Density Matrix", true);
        if(this->nC == 1 and not this->iCS )
          moDen[1].output(std::cout, "MO Beta Density Matrix", true);
      #endif

      size_t start = 0;
      if(this->nC == 1 ) {
        printBlock("Alpha occupied orbitals", start, this->nOA);
        printBlock("Alpha virtual orbitals", start, this->nVA);
        if( not this->iCS ){
          printBlock("Beta occupied orbitals", start, this->nOB);
          printBlock("Beta virtual orbitals", start, this->nVB);
        }
      }
      else if(this->nC == 2 ){
        printBlock("Occupied orbitals", start, this->nO);
        printBlock("Virtual orbitals", start, this->nV);
      } else if(this->nC == 4 ){
        start += NB/2;
        printBlock("Positive Energy Occupied orbitals", start, this->nO);
        printBlock("Positive Energy Virtual orbitals",  start, this->nV);
        start = 0;
        printBlock("Negative Energy Orbitals", start, this->nO+this->nV);

        // Print # of particles for 4C
        MatsT negTotal(MatsT(0.)), posTotal(MatsT(0.));
        for (size_t i = 0; i < NB; ++i) 
          if (i < NB/2) negTotal += moDen[0](i, i);
          else          posTotal += moDen[0](i, i);
        MatsT total = negTotal + posTotal;
        std::cout << std::scientific << std::setprecision(16);
        std::cout << "Negative Particles: " << negTotal << std::endl;
        std::cout << "Positive Particles: " << posTotal << std::endl;
        std::cout << "Total    Particles: " << total << std::endl;

      }
      std::cout << std::flush;
    }

};

template <typename MatsT, typename IntsT>
  void SingleSlater<MatsT,IntsT>::formEWDM(bool equil) {

    // ROOT_ONLY(comm);
    size_t NB = this->basisSet().nBasis;

    if( W != nullptr ) {
      W->clear();
    } else {
      if(not iCS and nC == 1 )//and basisSet().basisType == COMPLEX_GIAO)
        W = std::make_shared<cqmatrix::PauliSpinorMatrices<MatsT>>(NB, false);
      else if(nC == 2 or nC == 4)
        W = std::make_shared<cqmatrix::PauliSpinorMatrices<MatsT>>(NB, true);
      else
        W = std::make_shared<cqmatrix::PauliSpinorMatrices<MatsT>>(NB, false, false);
    }

    formEWDM_impl(equil, this->fockMatrix->S(), this->onePDM->S(), W->S());
    if (this->onePDM->hasZ())
      formEWDM_impl(equil, this->fockMatrix->Z(), this->onePDM->Z(), W->Z());
    if (this->onePDM->hasXY()) {
      formEWDM_impl(equil, this->fockMatrix->Y(), this->onePDM->Y(), W->Y());
      formEWDM_impl(equil, this->fockMatrix->X(), this->onePDM->X(), W->X());
    }

    // Form W strictly by definition for debugging
    // W = C * diag(ε) * Cdagger
#ifdef _DEBUG_EWDM
    //ao2orthoFock();
    //diagOrthoFock();
    //ortho2aoMOs();

    NB = NB * this->nC;

    if (this->nC == 1) {

      cqmatrix::Matrix<MatsT> WA(NB);
      cqmatrix::Matrix<MatsT> scaledMOA = this->mo[0];
      
      // WA = C * diag(ε) * Cdagger
      for(int i = 0; i < this->nOA; i++) {
        blas::scal(NB, this->eps1[i], scaledMOA.pointer() + i * NB, 1);
      }
      blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans, blas::Op::ConjTrans, NB, NB, this->nOA, MatsT(1.), scaledMOA.pointer(), NB,
          this->mo[0].pointer(), NB, MatsT(0.), WA.pointer(), NB);
      
      if (this->iCS) {

        *W = cqmatrix::PauliSpinorMatrices<MatsT>::spinBlockScatterBuild(WA);
      } else {

        cqmatrix::Matrix<MatsT> WB(NB);
        cqmatrix::Matrix<MatsT> scaledMOB = this->mo[1];

        // WB = C * diag(ε) * Cdagger
        for(int i = 0; i < this->nOB; i++) {
          blas::scal(NB, this->eps2[i], scaledMOB.pointer() + i * NB, 1);
        }
        blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans, blas::Op::ConjTrans, NB, NB, this->nOB, MatsT(1.), scaledMOB.pointer(), NB,
            this->mo[1].pointer(), NB, MatsT(1.), WB.pointer(), NB);

        *W = cqmatrix::PauliSpinorMatrices<MatsT>::spinBlockScatterBuild(WA, WB);
      }
    } else {

      if (this->nC == 2) {
        cqmatrix::Matrix<MatsT> W_spinBlockForm(NB);
        cqmatrix::Matrix<MatsT> scaledMO = this->mo[0];
        
        // W = C * diag(ε) * Cdagger
        for(int i = 0; i < this->nO; i++) {
          blas::scal(NB, this->eps1[i], scaledMO.pointer() + i * NB, 1);
        }
        blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans, blas::Op::ConjTrans, NB, NB, this->nO, MatsT(1.), scaledMO.pointer(), NB,
            scaledMO.pointer(), NB, MatsT(1.), W_spinBlockForm.pointer(), NB);

        *W = W_spinBlockForm.template spinScatter<MatsT>();
            
      } else if (this->nC == 4) {
        CErr("4C not implemented in SingleSlater<MatsT,IntsT>::formEWDM!");
      }
    }
#endif

    //W->output(std::cout, "W matrix", true);
};

template <typename MatsT, typename IntsT>
  void SingleSlater<MatsT,IntsT>::formEWDM_impl(bool equil, const cqmatrix::Matrix<MatsT>& F, 
      const cqmatrix::Matrix<MatsT>& P, cqmatrix::Matrix<MatsT>& W) {

    const size_t NB = this->basisSet().nBasis;
    cqmatrix::Matrix<MatsT> SCR(NB, NB);

    if (equil) {
      // W = (1/2) * (P * F * P)
      blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,NB,NB,NB,MatsT(1.),P.pointer(),NB,
          F.pointer(),NB,MatsT(0.),SCR.pointer(),NB);
      blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,NB,NB,NB,MatsT(1.),SCR.pointer(),NB,
          P.pointer(),NB,MatsT(0.0),W.pointer(),NB);
      W *= 0.5;
    } else {
      // For general case, SW + WS = FP + PF
      // RHS = FP + PF
      blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans,
                 NB, NB, NB, MatsT(1), F.pointer(), NB, P.pointer(), NB,
                 MatsT(0), W.pointer(), NB);
      blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans,
                 NB, NB, NB, MatsT(1), P.pointer(), NB, F.pointer(), NB,
                 MatsT(0), SCR.pointer(), NB);
      W += SCR;
    
      // Sylvester equation solver
      // Diagonalize S
      cqmatrix::Matrix<MatsT> X = this->aoints_->overlap->matrix();
      double* eigS = CQMemManager::get().malloc<double>(NB);
      int INFO = HermitianEigen('V','U', NB, X.pointer(), NB, eigS);
      if (INFO) CErr("HermitianEigen failed in formEWDM_impl", std::cout);
    
      // Transform into the eigenspace of the overlap matrix
      blas::gemm(blas::Layout::ColMajor, blas::Op::ConjTrans, blas::Op::NoTrans,
        NB, NB, NB, MatsT(1), X.pointer(), NB, W.pointer(), NB,
        MatsT(0), SCR.pointer(), NB);
      blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans,
              NB, NB, NB, MatsT(1), SCR.pointer(), NB, X.pointer(), NB,
              MatsT(0), W.pointer(), NB);
    
      // Normalize
      double smax = 0.0; for (int i=0;i<NB;++i) smax = std::max(smax, eigS[i]);
      const double eps = 1e-14 * std::max(1.0, smax);
      for (int j=0; j<NB; ++j)
        for (int i=0; i<NB; ++i) {
          double denom = eigS[i] + eigS[j];
          if (std::abs(denom) < eps) denom = (denom < 0 ? -eps : eps);
          W(i,j) = W(i,j) / denom;
        }
    
      // Transform back to AO basis
      blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans,
                 NB, NB, NB, MatsT(1), X.pointer(), NB, W.pointer(), NB,
                 MatsT(0), SCR.pointer(), NB);
      blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::ConjTrans,
                 NB, NB, NB, MatsT(1), SCR.pointer(), NB, X.pointer(), NB,
                 MatsT(0), W.pointer(), NB);
    
      // Hermitize
      for (int j = 0; j < NB; ++j) {
        for (int i = 0; i <= j; ++i) {
          if constexpr (std::is_same<MatsT,double>::value) {
            const MatsT hij = (W(i,j) + W(j,i)) * 0.5;
            W(i,j) = hij; W(j,i) = hij;
          } else {
            const MatsT hij = (W(i,j) + std::conj(W(j,i))) * 0.5;
            W(i,j) = hij; W(j,i) = std::conj(hij);
          }
        }
      }
      CQMemManager::get().free(eigS);
    }
    
    W *= 0.5;
}; // SingleSlater<MatsT,IntsT>::formEWDM_impl

}; // namespace ChronusQ

