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

#include <singleslater.hpp>
#include <singleslater/neoss.hpp>
#include <cqlinalg/blasext.hpp>
#include <cqlinalg/blasutil.hpp>
#include <cqlinalg/blas3.hpp>
#include <quantum/properties.hpp>

namespace ChronusQ {

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

    size_t NB  = this->nAlphaOrbital() * nC;
    auto* fb = this->fockBuilder.get();
    if (auto* neofb = dynamic_cast<NEOFockBuilder<MatsT, IntsT>*>(fb)) {
      fb = neofb->getNonNEOUpstream();              // still a raw pointer
    }
    bool iRO = (dynamic_cast<ROFock<MatsT, IntsT>*>(fb) != nullptr);

    // ROHF copy modified orbitals to redundant set
    if( iRO ){
      std::copy_n(this->mo[0].pointer(),NB*NB,this->mo[1].pointer());
      std::copy_n(this->eps1,NB,this->eps2);
    }

    // Form the 1PDM on the root MPI process as slave processes
    // do not posses the up-to-date MO coefficients
    if( MPIRank(comm) == 0 ) {

      if(nC == 1) {

        cqmatrix::Matrix<MatsT> DA(NB);

        //this->mo[0].output(std::cout, "mo1", true);

        // DA = CA * CA**H
        blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans, blas::Op::ConjTrans, NB, NB, this->nOA, MatsT(1.), this->mo[0].pointer(), NB,
            this->mo[0].pointer(), NB, MatsT(0.), DA.pointer(), NB);

        if(iCS) {

          // DS = 2 * DA
          *this->onePDM = cqmatrix::PauliSpinorMatrices<MatsT>::spinBlockScatterBuild(DA);

        } else {

          cqmatrix::Matrix<MatsT> DB(NB);

          // DB = CB * CB**H
          blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans, blas::Op::ConjTrans, NB, NB, this->nOB, MatsT(1.), this->mo[1].pointer(), NB,
              this->mo[1].pointer(), NB, MatsT(0.), DB.pointer(), NB);

          // DS = DA + DB
          // DZ = DA - DB
          *this->onePDM = cqmatrix::PauliSpinorMatrices<MatsT>::spinBlockScatterBuild(DA,DB);

        }
      } else {

        // 2C or 4C cases
        cqmatrix::Matrix<MatsT> spinBlockForm(NB);

        if( nC == 2 ) {
          blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans, blas::Op::ConjTrans, NB, NB, this->nO, MatsT(1.), this->mo[0].pointer(), NB,
              this->mo[0].pointer(), NB, MatsT(0.), spinBlockForm.pointer(), NB);
        } else if( nC == 4 ) {
          blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans, blas::Op::ConjTrans, NB, NB, this->nO, MatsT(1.), this->mo[0].pointer()+(2*(NB/nC))*NB, NB,
              this->mo[0].pointer()+(2*(NB/nC))*NB, NB, MatsT(0.), spinBlockForm.pointer(), NB);
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

    this->SExpect[0] = 0.5 * this->template computeOBProperty<MX>(
      this->aoints_->overlap->pointer());
    this->SExpect[1] = 0.5 * this->template computeOBProperty<MY>(
      this->aoints_->overlap->pointer());
    this->SExpect[2] = 0.5 * this->template computeOBProperty<MZ>(
      this->aoints_->overlap->pointer());

    if( not this->onePDM->hasZ() ) this->SSq = 0;
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
  

      // SCR2 = D(Z) * S * D(Z)
      blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,NB,NB,NB,MatsT(1.),this->aoints_->overlap->pointer(),NB,
           this->onePDM->Z().pointer(),NB,MatsT(0.),SCR,NB);
      blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::ConjTrans,NB,NB,NB,MatsT(1.),this->aoints_->overlap->pointer(),NB,
           SCR,NB,MatsT(0.),SCR2,NB);
      
      this->SSq += 0.5 * this->template computeOBProperty<MZ>(SCR2);

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

      CQMemManager::get().free(SCR,SCR2);
    }

  };
  
template <typename MatsT, typename IntsT>
  void SingleSlater<MatsT,IntsT>::compute4CDipole(EMPerturbation &pert) {
    ROOT_ONLY(comm);

    // Compute elecric contribution to the dipoles
    for(auto iXYZ = 0; iXYZ < 3; iXYZ++) {
      // Scalar
      double dipole_s = -this->template computeOBProperty<DENSITY_TYPE::SCALAR>(
          (*(this->aoints_->lenElectric4C))[iXYZ].S().pointer());
      // MZ
      double dipole_z = -this->template computeOBProperty<DENSITY_TYPE::MZ>(
          (*(this->aoints_->lenElectric4C))[iXYZ].Z().pointer());
      // MY
      double dipole_y= -this->template computeOBProperty<DENSITY_TYPE::MY>(
          (*(this->aoints_->lenElectric4C))[iXYZ].Y().pointer());
      // MX
      double dipole_x= -this->template computeOBProperty<DENSITY_TYPE::MX>(
          (*(this->aoints_->lenElectric4C))[iXYZ].X().pointer());
    
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
      SetMat('N',NB/4,NB/4,MatsT(1./(2*SpeedOfLight*SpeedOfLight)),this->aoints_->kinetic->matrix().pointer(), NB/4, S.pointer()+NB*NB/4+NB/4,NB);
      SetMat('N',NB/4,NB/4,MatsT(1.),this->aoints_->overlap->matrix().pointer(), NB/4, S.pointer()+NB*NB/2+NB/2,NB);
      SetMat('N',NB/4,NB/4,MatsT(1./(2*SpeedOfLight*SpeedOfLight)),this->aoints_->kinetic->matrix().pointer(), NB/4, S.pointer()+NB*NB*3/4+NB*3/4,NB);
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
    // W = C * diag(ε) * C†
#ifdef _DEBUG_EWDM
    //ao2orthoFock();
    //diagOrthoFock();
    //ortho2aoMOs();

    NB = NB * this->nC;

    if (this->nC == 1) {

      cqmatrix::Matrix<MatsT> WA(NB);
      cqmatrix::Matrix<MatsT> scaledMOA = this->mo[0];
      
      // WA = C * diag(ε) * C†
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

        // WB = C * diag(ε) * C†
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
        
        // W = C * diag(ε) * C†
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

