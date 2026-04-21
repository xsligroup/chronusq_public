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

#include <corehbuilder/nonrel.hpp>
#include <matrix.hpp>

namespace ChronusQ {

  template <>
  void NRCoreH<dcomplex, dcomplex>::addMagPert(EMPerturbation &pert,
      std::shared_ptr<cqmatrix::PauliSpinorMatrices<dcomplex>> coreH, double mass, double charge) {

    //Compute the GIAO non-relativistic core Hamiltonian in the CGTO basis
    //H(S) = 2(T + V) + B * L + sigma * B + 1/4 *(B\timesr)^2

    dcomplex onei = dcomplex(0,1);
    auto magAmp = pert.getDipoleAmp(Magnetic);


    // this part add the angular momentum term
    for ( auto index = 0 ; index < 3 ; index++ ) {
      *coreH += magAmp[index] * onei * charge * (1.0/mass) * (*aoints_.magnetic)[index]->matrix();
    } // for ( auto inde = 0 ; inde < 3 ; inde++ )

    // this part add the length gauge electric quadrupole term
    const std::array<std::string,3> diagindex =
      { "XX","YY","ZZ" };

    double diagcoeff[3];
    diagcoeff[0] = 1.0/8.0*(magAmp[1]*magAmp[1]+magAmp[2]*magAmp[2]);
    diagcoeff[1] = 1.0/8.0*(magAmp[0]*magAmp[0]+magAmp[2]*magAmp[2]);
    diagcoeff[2] = 1.0/8.0*(magAmp[0]*magAmp[0]+magAmp[1]*magAmp[1]);

    // add diagonal part
    for ( size_t index = 0 ; index < 3 ; index++ ) {
      *coreH -= 2.0*diagcoeff[index] * charge * (1.0/mass) * (*aoints_.lenElectric)[diagindex[index]]->matrix();
    }

    const std::array<std::string,3> offindex =
      { "XY","XZ","YZ" };

    double offcoeff[3];
    offcoeff[0] = -1.0/4.0*magAmp[0]*magAmp[1];
    offcoeff[1] = -1.0/4.0*magAmp[0]*magAmp[2];
    offcoeff[2] = -1.0/4.0*magAmp[1]*magAmp[2];

    // add off diagonal part
    for ( size_t index = 0 ; index < 3 ; index++ ) {
      *coreH -= 2.0*offcoeff[index] * charge * (1.0/mass) * (*aoints_.lenElectric)[offindex[index]]->matrix();
    }

    // finally spin Zeeman term
    // z component
    if(coreH->hasZ())
      coreH->Z() = -magAmp[2] * charge * (1.0/mass) * aoints_.overlap->matrix();

    if(coreH->hasXY()) {
      // y component
      coreH->Y() = -magAmp[1] * charge * (1.0/mass) * aoints_.overlap->matrix();

      // x coponent
      coreH->X() = -magAmp[0] * charge * (1.0/mass) * aoints_.overlap->matrix();
    }


  }

  template <>
  void NRCoreH<dcomplex, double>::addMagPert(EMPerturbation &pert,
    std::shared_ptr<cqmatrix::PauliSpinorMatrices<dcomplex>> coreH, double mass, double charge) {

    std::cout<<"GTO magnetic field is experimental and not recommended!"<<std::endl;

    dcomplex onei = dcomplex(0,1);
    auto magAmp = pert.getDipoleAmp(Magnetic);


    // this part add the angular momentum term
    for ( auto index = 0 ; index < 3 ; index++ ) {
      *coreH += magAmp[index] * onei * charge * (1.0/mass) * (*aoints_.magnetic)[index]->matrix();
    } // for ( auto inde = 0 ; inde < 3 ; inde++ )

    // this part add the length gauge electric quadrupole term
    const std::array<std::string,3> diagindex =
      { "XX","YY","ZZ" };

    double diagcoeff[3];
    diagcoeff[0] = 1.0/8.0*(magAmp[1]*magAmp[1]+magAmp[2]*magAmp[2]);
    diagcoeff[1] = 1.0/8.0*(magAmp[0]*magAmp[0]+magAmp[2]*magAmp[2]);
    diagcoeff[2] = 1.0/8.0*(magAmp[0]*magAmp[0]+magAmp[1]*magAmp[1]);

    // add diagonal part
    for ( size_t index = 0 ; index < 3 ; index++ ) {
      *coreH -= 2.0*diagcoeff[index] * charge * (1.0/mass) * (*aoints_.lenElectric)[diagindex[index]]->matrix();
    }

    const std::array<std::string,3> offindex =
      { "XY","XZ","YZ" };

    double offcoeff[3];
    offcoeff[0] = -1.0/4.0*magAmp[0]*magAmp[1];
    offcoeff[1] = -1.0/4.0*magAmp[0]*magAmp[2];
    offcoeff[2] = -1.0/4.0*magAmp[1]*magAmp[2];

    // add off diagonal part
    for ( size_t index = 0 ; index < 3 ; index++ ) {
      *coreH -= 2.0*offcoeff[index] * charge * (1.0/mass) * (*aoints_.lenElectric)[offindex[index]]->matrix();
    }

    // finally spin Zeeman term
    // z component
    if(coreH->hasZ())
      coreH->Z() = -magAmp[2] * charge * (1.0/mass) * aoints_.overlap->matrix();

    if(coreH->hasXY()) {
      // y component
      coreH->Y() = -magAmp[1] * charge * (1.0/mass) * aoints_.overlap->matrix();

      // x coponent
      coreH->X() = -magAmp[0] * charge * (1.0/mass) * aoints_.overlap->matrix();
    }

  }

  template <>
  void NRCoreH<double, double>::addMagPert(EMPerturbation&,
    std::shared_ptr<cqmatrix::PauliSpinorMatrices<double>>, double mass, double charge) {

    CErr("Orbital Zeeman + Real integrals is not a valid option");

  }

  /**
   *  \brief Compute the non-relativistic Core Hamiltonian in the CGTO basis
   *
   *  \f[ H(S) = 2(T + V) \f]
   */
  template <typename MatsT, typename IntsT>
  void NRCoreH<MatsT,IntsT>::computeCoreH(EMPerturbation& emPert,
      std::shared_ptr<cqmatrix::PauliSpinorMatrices<MatsT>> coreH) {

    computeNRCH(emPert, coreH);

  };  // void NRCoreH::computeCoreH(std::vector<MatsT*> &CH)

  /**
   *  \brief Compute the non-relativistic Core Hamiltonian in the CGTO basis
   *
   *  \f[ H(S) = 2(T + V) \f]
   */
  template <typename MatsT, typename IntsT>
  void NRCoreH<MatsT,IntsT>::computeNRCH(EMPerturbation& emPert,
      std::shared_ptr<cqmatrix::PauliSpinorMatrices<MatsT>> coreH) {

    *coreH = 2. * (this->aoints_.kinetic->matrix() + this->aoints_.potential->matrix());

    if( this->hamiltonianOptions_.basisType == COMPLEX_GIAO and pert_has_type(emPert,Magnetic) )
      addMagPert(emPert,coreH,this->hamiltonianOptions_.particle.mass,this->hamiltonianOptions_.particle.charge);


#ifdef _DEBUGGIAOONEE
      this->aoints_.kinetic->matrix().output(std::cout, "Kinetic", true);
      this->aoints_.potential->matrix().output(std::cout, "Potential", true);
      this->aoints_.overlap->matrix().output(std::cout, "Overlap", true);
      coreH->output(std::cout, "Core Ham", true);
#endif

  };  // void NRCoreH::computeCoreH(std::vector<MatsT*> &CH)

  template <typename MatsT, typename IntsT>
  std::vector<double> NRCoreH<MatsT,IntsT>::getGrad(EMPerturbation& pert,
    SingleSlater<MatsT,IntsT>& ss) {

    if (not this->aoints_.gradKinetic)
      CErr("Kinetic energy gradient integrals missing in NRCoreH::getGrad!");
    if (not this->aoints_.gradPotential)
      CErr("Potential energy gradient integrals missing in NRCoreH::getGrad!");

    GradInts<OnePInts,IntsT>& gradK = *this->aoints_.gradKinetic;
    GradInts<OnePInts,IntsT>& gradV = *this->aoints_.gradPotential;

    //gradK.output(std::cout, "gradKinetic", true); 
    //gradV.output(std::cout, "gradPotential", true);

    size_t NB = ss.basisSet().nBasis;
    size_t nGrad = 3*ss.molecule().nAtoms;

    std::vector<double> gradient;
    std::shared_ptr<cqmatrix::PauliSpinorMatrices<MatsT>> coreHGrad;
    if (not ss.iCS and ss.nC == 1)// and ss.basisSet().basisType == COMPLEX_GIAO)
      coreHGrad = std::make_shared<cqmatrix::PauliSpinorMatrices<MatsT>>(NB, false);
    else if (ss.nC == 2 or ss.nC == 4)
      coreHGrad = std::make_shared<cqmatrix::PauliSpinorMatrices<MatsT>>(NB, true);
    else
      coreHGrad = std::make_shared<cqmatrix::PauliSpinorMatrices<MatsT>>(NB, false, false);
    
    // Loop over gradient components
    for ( auto iGrad = 0; iGrad < nGrad; iGrad++ ) {
    
      // Assemble core gradient
      *coreHGrad = 2. * (gradK[iGrad]->matrix() + gradV[iGrad]->matrix());

      // Contract
      double element = ss.template computeOBProperty<DENSITY_TYPE::SCALAR>(
        coreHGrad->S().pointer());
      if (coreHGrad->hasZ()) {
        element += ss.template computeOBProperty<DENSITY_TYPE::MZ>(coreHGrad->Z().pointer());
      }
      if (coreHGrad->hasXY()) { 
        element += ss.template computeOBProperty<DENSITY_TYPE::MY>(coreHGrad->Y().pointer());
        element += ss.template computeOBProperty<DENSITY_TYPE::MX>(coreHGrad->X().pointer());
      }
      gradient.emplace_back(0.5*element);

    }

    return gradient;

  };  // std::vector<double> NRCoreH::getGrad

  template <>
  std::vector<double> NRCoreH<dcomplex,dcomplex>::getGrad(EMPerturbation& pert,
    SingleSlater<dcomplex,dcomplex>& ss) {

    //std::cout << "TangDD: mag core grad" << std::endl;

    if (not this->aoints_.gradKinetic)
      CErr("Kinetic energy gradient integrals missing in NRCoreH::getGrad!");
    if (not this->aoints_.gradPotential)
      CErr("Potential energy gradient integrals missing in NRCoreH::getGrad!");

    GradInts<OnePInts,dcomplex>& gradK = *this->aoints_.gradKinetic;
    GradInts<OnePInts,dcomplex>& gradV = *this->aoints_.gradPotential;

    //gradK.output(std::cout, "gradKinetic", true); 
    //gradV.output(std::cout, "gradPotential", true);

    size_t NB = ss.basisSet().nBasis;
    size_t nGrad = 3*ss.molecule().nAtoms;

    double mass_ = this->hamiltonianOptions_.particle.mass;
    double charge_ = this->hamiltonianOptions_.particle.charge;

    std::vector<double> gradient;
    std::shared_ptr<cqmatrix::PauliSpinorMatrices<dcomplex>> coreHGrad;
    if (not ss.iCS and ss.nC == 1)// and ss.basisSet().basisType == COMPLEX_GIAO)
      coreHGrad = std::make_shared<cqmatrix::PauliSpinorMatrices<dcomplex>>(NB, false);
    else if (ss.nC == 2 or ss.nC == 4)
      coreHGrad = std::make_shared<cqmatrix::PauliSpinorMatrices<dcomplex>>(NB, true);
    else
      coreHGrad = std::make_shared<cqmatrix::PauliSpinorMatrices<dcomplex>>(NB, false, false);

      dcomplex onei = dcomplex(0,1);
      auto magAmp = pert.getDipoleAmp(Magnetic);
 
      if (not this->aoints_.gradOverlap)
        CErr("Overlap gradient integrals missing in NRCoreH::getGrad!");
      if (not this->aoints_.gradLen2)
        CErr("Quadropole gradient integrals missing in NRCoreH::getGrad!");
      if (not this->aoints_.gradMagnetic)
        CErr("Angular gradient integrals missing in NRCoreH::getGrad!");
  
      GradInts<OnePInts,dcomplex>& gradS = *this->aoints_.gradOverlap;
      GradInts<VectorInts,dcomplex>& gradQ = *this->aoints_.gradLen2;
      GradInts<VectorInts,dcomplex>& gradL = *this->aoints_.gradMagnetic;
    
      // Loop over gradient components
      for ( auto iGrad = 0; iGrad < nGrad; iGrad++ ) {
 
        *coreHGrad = 2. * (gradK[iGrad]->matrix() + gradV[iGrad]->matrix());

        // this part add the angular momentum term
        for ( auto index = 0 ; index < 3 ; index++ ) {
          *coreHGrad += magAmp[index] * onei * charge_ * (1.0/mass_) * (*gradL[iGrad])[index]->matrix();
        } // for ( auto inde = 0 ; inde < 3 ; inde++ )
    
        // this part add the length gauge electric quadrupole term
        const std::array<std::string,3> diagindex =
          { "XX","YY","ZZ" };
    
        double diagcoeff[3];
        diagcoeff[0] = 1.0/8.0*(magAmp[1]*magAmp[1]+magAmp[2]*magAmp[2]);
        diagcoeff[1] = 1.0/8.0*(magAmp[0]*magAmp[0]+magAmp[2]*magAmp[2]);
        diagcoeff[2] = 1.0/8.0*(magAmp[0]*magAmp[0]+magAmp[1]*magAmp[1]);
    
        // add diagonal part
        for ( size_t index = 0 ; index < 3 ; index++ ) {
          *coreHGrad -= 2.0*diagcoeff[index] * charge_ * (1.0/mass_) * (*gradQ[iGrad])[diagindex[index]]->matrix();
        }
    
        const std::array<std::string,3> offindex =
          { "XY","XZ","YZ" };
    
        double offcoeff[3];
        offcoeff[0] = -1.0/4.0*magAmp[0]*magAmp[1];
        offcoeff[1] = -1.0/4.0*magAmp[0]*magAmp[2];
        offcoeff[2] = -1.0/4.0*magAmp[1]*magAmp[2];
    
        // add off diagonal part
        for ( size_t index = 0 ; index < 3 ; index++ ) {
          *coreHGrad -= 2.0*offcoeff[index] * charge_ * (1.0/mass_) * (*gradQ[iGrad])[offindex[index]]->matrix();
        }
        
        // finally spin Zeeman term
        // z component
        if(coreHGrad->hasZ())
        coreHGrad->Z() = -magAmp[2] * charge_ * (1.0/mass_) * gradS[iGrad]->matrix();
        
        if(coreHGrad->hasXY()) {
          // y component
          coreHGrad->Y() = -magAmp[1] * charge_ * (1.0/mass_) * gradS[iGrad]->matrix();
          
          // x coponent
          coreHGrad->X() = -magAmp[0] * charge_ * (1.0/mass_) * gradS[iGrad]->matrix();
        }

        //coreHGrad->output(std::cout, "dCore", true); 
        
        // Contract
        double element = ss.template computeOBProperty<DENSITY_TYPE::SCALAR>(
          coreHGrad->S().pointer());
        if (ss.onePDM->hasZ()) {
          element += ss.template computeOBProperty<DENSITY_TYPE::MZ>(coreHGrad->Z().pointer());
        }
        if (ss.onePDM->hasXY()) { 
          element += ss.template computeOBProperty<DENSITY_TYPE::MY>(coreHGrad->Y().pointer());
          element += ss.template computeOBProperty<DENSITY_TYPE::MX>(coreHGrad->X().pointer());
        }
        gradient.emplace_back(0.5*element);
  
      }

    return gradient;

  };  // std::vector<double> NRCoreH::getGrad

}
