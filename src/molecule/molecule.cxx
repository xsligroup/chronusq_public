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
#include <molecule.hpp>
#include <physcon.hpp>
#include <cxxapi/output.hpp>

namespace ChronusQ {


  /**
   *  \brief Compute the internuclear distances for the Atoms contained
   *  int the atoms array
   *
   *  \f[
   *    R_{AB} = \vert R_A - R_B \vert
   *  \f]
   */ 
  void Molecule::computeRIJ() {

    RIJ = dynmat_t(nAtoms, dynvec_t(nAtoms,0.));
    double tmp;

    for(size_t iAtm = 0; iAtm < nAtoms; iAtm++) 
    for(size_t jAtm = 0; jAtm < iAtm  ; jAtm++) {
      tmp = atoms[iAtm].coord[0] - atoms[jAtm].coord[0];
      RIJ[iAtm][jAtm] = tmp*tmp;
      tmp = atoms[iAtm].coord[1] - atoms[jAtm].coord[1];
      RIJ[iAtm][jAtm] += tmp*tmp;
      tmp = atoms[iAtm].coord[2] - atoms[jAtm].coord[2];
      RIJ[iAtm][jAtm] += tmp*tmp;

      RIJ[iAtm][jAtm] = std::sqrt(RIJ[iAtm][jAtm]);
      RIJ[jAtm][iAtm] = RIJ[iAtm][jAtm];
    }
      
  }; // Molecule::computeRIJ

  /**
   *  \brief Compute the nuclear-nuclear repulsion energy for classical
   *  point nuclei using the Atoms contained in the atoms array
   *
   *  \f[
   *    V_{NN} = \sum_{A < B} \frac{Z_A Z_B}{R_{AB}}
   *  \f]
   */ 
  void Molecule::computeNNRep() {
    nucRepEnergy = 0.;
    for(size_t iAtm = 0; iAtm < nAtoms; iAtm++)
    for(size_t jAtm = 0; jAtm < iAtm  ; jAtm++)
      if (not atoms[iAtm].quantum and not atoms[jAtm].quantum )
        nucRepEnergy +=
          atoms[iAtm].nucCharge * atoms[jAtm].nucCharge /
          RIJ[iAtm][jAtm];
  }

  /**
   *  \brief Compute the force contribution from the other
   *  nuclei contained in the atoms array
   *
   *  \f[
   *    \frac{\partial V_{NN}}{\partial X_A} = 
   *       Z_A \sum_B \frac{Z_B (X_B - X_A)}{R_{AB}^3}
   *  \f]
   */ 
  void Molecule::computeNNX() {
    nucRepForce = dynmat_t(nAtoms,dynvec_t(3,0.));

    for(size_t iAtm = 0; iAtm < nAtoms; iAtm++) 
    for(size_t iXYZ = 0; iXYZ < 3     ; iXYZ++) { 
      if( atoms[iAtm].quantum )
        continue;
      for(size_t jAtm = 0; jAtm < nAtoms; jAtm++){ 
        if( atoms[jAtm].quantum )
          continue;
        if(iAtm == jAtm) continue;
        nucRepForce[iAtm][iXYZ] += 
          atoms[jAtm].nucCharge *
          ( atoms[jAtm].coord[iXYZ] - atoms[iAtm].coord[iXYZ] ) /
          RIJ[iAtm][jAtm] / RIJ[iAtm][jAtm] / RIJ[iAtm][jAtm]; 
      }

      nucRepForce[iAtm][iXYZ] *= atoms[iAtm].nucCharge;
    }
  } // Molecule::computeNNX


  /**
   *  \brief Compute the center-of-mass of the Atoms contained in the
   *  atoms array
   *
   *  \f[
   *    \vec{C}_M = \frac{1}{M} \sum_A M_A \vec{R}_A \qquad M = \sum_A M_A
   *  \f]
   */
  void Molecule::computeCOM() {

    std::fill(COM.begin(),COM.end(),0.);

    for(size_t iAtm = 0; iAtm < nAtoms; iAtm++) 
    for(size_t iXYZ = 0; iXYZ < 3     ; iXYZ++) {
      COM[iXYZ] += atoms[iAtm].atomicMass * atoms[iAtm].coord[iXYZ];
    }

    double totalMass(0.);
    for(size_t iAtm = 0; iAtm < nAtoms; iAtm++) 
      totalMass += atoms[iAtm].atomicMass;

    std::transform(COM.begin(),COM.end(),COM.begin(),
      [&](double x){ return x / totalMass; }
    );

  } // Molecule::computeCOM

  /**
   *  \brief Compute the center-of-charges of the Atoms contained in the
   *  atoms array
   *
   *  \f[
   *    \vec{C}_C = \frac{1}{C} \sum_A C_A \vec{R}_A \qquad C = \sum_A C_A
   *  \f]
   */
  void Molecule::computeCOC() {

    std::fill(COC.begin(),COC.end(),0.);

    for(size_t iAtm = 0; iAtm < nAtoms; iAtm++) 
    for(size_t iXYZ = 0; iXYZ < 3     ; iXYZ++) {
      COC[iXYZ] += atoms[iAtm].nucCharge * atoms[iAtm].coord[iXYZ];
    }

    std::transform(COC.begin(),COC.end(),COC.begin(),
      [&](double x){ return x / ( nTotalE - charge ); }
    );

  } // Molecule::computeCOC

  /**
   *  \brief Compute the moment of inertia tensor of the Atoms 
   *  contained in the atoms array
   *
   *  \f[
   *    \mathbf{I} = \sum_A M_A \left(
   *      \left(\vec{R}_A \cdot \vec{R}_A \right)\mathbf{E} - 
   *      \vec{R}_A \otimes \vec{R}_A
   *    \right) \qquad E_{ij} = \delta_{ij}
   *  \f]
   */
  void Molecule::computeMOI() {

    std::fill_n(MOI.begin()->begin(),9,0.);

    double tmp;
    for(size_t iAtm = 0; iAtm < nAtoms; iAtm++) {
      tmp  = atoms[iAtm].coord[0] * atoms[iAtm].coord[0];
      tmp += atoms[iAtm].coord[1] * atoms[iAtm].coord[1];
      tmp += atoms[iAtm].coord[2] * atoms[iAtm].coord[2];

      for(size_t iXYZ = 0; iXYZ < 3; iXYZ++){
        MOI[iXYZ][iXYZ] += atoms[iAtm].atomicMass * tmp;
        for(size_t jXYZ = 0; jXYZ < 3; jXYZ++)
          MOI[iXYZ][jXYZ] -= atoms[iAtm].atomicMass *
            atoms[iAtm].coord[iXYZ] * atoms[iAtm].coord[jXYZ];
      }
    }

  } // Molecule::computeMOI


  /**
   *  \breif Generate gaussian shell definitions for the nuclear
   *  charge distribution.
   *
   *  L. Visscher and K. G. Dyall; Atomic Data and Nuclear Data Tables; 67,
   *    207-224 (1997)
   */ 
  void Molecule::computeCDist() {


    for(auto &atom : atoms) {
      double varience = 
        0.836 * std::pow(atom.massNumber,1.0/3.0) + 0.570; // fm

      varience *= 1e-5; // Ang
      varience /= AngPerBohr; // Bohr

      varience *= varience;

      double zeta = 3. / 2. / varience;

      chargeDist.push_back(
        libint2::Shell {
          { zeta }, 
          {{0,false,{atom.nucCharge}}},
          atom.coord
        }
      );

      // Handle the fact that libint likes to make things square normalized
      // not normalized
      chargeDist.back().contr[0].coeff[0] = 
        atom.nucCharge * std::pow(zeta / M_PI, 1.5);

    } // loop over atoms

  }; // Molecule::computeCDist



  /**
   *  Outputs relevant information for the Molecule object
   *  to a specified output.
   *
   *  \param [in/out] out Ouput device
   *  \param [in]     mol Molecule object to output.
   */ 
  std::ostream& operator<<(std::ostream &out, const Molecule &mol) {
    
/*
    int iProc, nProc ;

    MPI_Comm_size(comm,&nProc);
    MPI_Comm_rank(comm,&iProc);
*/

    out << std::endl << "Molecular Information";
 // if(nProc > 1) out << " on MPI Process " << iProc;
    out << ":" << std::endl << BannerTop << std::endl << std::endl;


    const int fieldNameWidth(28);
    const int fieldValueWidth(15);
    const int shortFieldWidth(10);

    out << std::left << std::setprecision(8) << std::scientific;

    out << "  " << std::setw(fieldNameWidth) << "NAtoms" << mol.nAtoms 
        << std::endl;
    out << "  " << std::setw(fieldNameWidth) << "Total Electrons" 
        << mol.nTotalE << std::endl;
    out << "  " << std::setw(fieldNameWidth) << "Particle Number Charge"
        << std::setw(fieldValueWidth) <<  mol.charge << "e"
        << std::endl;
    double molecCharge=-double(mol.nTotalE);
    for(auto i=0; i<mol.atoms.size(); i++) molecCharge += mol.atoms[i].nucCharge;
    out << "  " << std::setw(fieldNameWidth) << "Molecular Charge"
        << std::setw(fieldValueWidth) << std::fixed << std::setprecision(3) <<  molecCharge << "e"
        << std::endl;
    out << std::left << std::setprecision(8) << std::scientific;
    out << "  " << std::setw(fieldNameWidth) <<  "Nuclear Repulsion Energy"
        << std::setw(fieldValueWidth) << mol.nucRepEnergy << "Eh" << std::endl;

    out << std::endl << std::setprecision(5) << std::scientific;
    out << "  " << std::setw(fieldNameWidth) << std::left << "Center of Mass"
        << "{ " 
        << std::setw(fieldValueWidth) << mol.COM[0] << "," 
        << std::setw(fieldValueWidth) << mol.COM[1] << "," 
        << std::setw(fieldValueWidth) << mol.COM[2] 
        << "}" << std::endl; 
    out << "  " << std::setw(fieldNameWidth) << std::left << "Center of Charges"
        << "{ " 
        << std::setw(fieldValueWidth) << mol.COC[0] << "," 
        << std::setw(fieldValueWidth) << mol.COC[1] << "," 
        << std::setw(fieldValueWidth) << mol.COC[2] 
        << "}" << std::endl; 
    out << std::endl;
    out << "  " << std::setw(fieldNameWidth) << std::left << "Moment of Inertia"
        << "{ " 
        << std::setw(fieldValueWidth) << mol.MOI[0][0] << "," 
        << std::setw(fieldValueWidth) << mol.MOI[0][1] << "," 
        << std::setw(fieldValueWidth) << mol.MOI[0][2] << "," 
        << std::endl << "  " << std::setw(fieldNameWidth) << " " << "  " 
        << std::setw(fieldValueWidth) << mol.MOI[1][0] << "," 
        << std::setw(fieldValueWidth) << mol.MOI[1][1] << "," 
        << std::setw(fieldValueWidth) << mol.MOI[1][2] << ","
        << std::endl << "  " << std::setw(fieldNameWidth) << " " << "  " 
        << std::setw(fieldValueWidth) << mol.MOI[2][0] << "," 
        << std::setw(fieldValueWidth) << mol.MOI[2][1] << "," 
        << std::setw(fieldValueWidth) << mol.MOI[2][2] 
        << "}" << std::endl; 


    out << std::endl;
    out << "  " << std::setw(fieldNameWidth) << "Geometry:" << std::endl;
    out << bannerMid << std::endl;
    out << "    " << std::setw(shortFieldWidth) << std::left << "Element";
    out << std::left << std::setw(shortFieldWidth-1) << "ZNuc";
    out << std::right << std::setw(shortFieldWidth) << "Mass (AMU)";
    out << std::right << std::setw(fieldValueWidth) << "X (Bohr)";
    out << std::right << std::setw(fieldValueWidth) << "Y (Bohr)";
    out << std::right << std::setw(fieldValueWidth) << "Z (Bohr)";
    out << std::endl << std::endl;
    for(auto iAtm = 0; iAtm < mol.atoms.size(); iAtm++){

      out << std::left << std::setprecision(5);
      std::map<std::string,Atom>::const_iterator it = 
      std::find_if(atomicReference.begin(),atomicReference.end(),
        [&](const std::pair<std::string,Atom> &st){ 
          return (st.second.atomicNumber == mol.atoms[iAtm].atomicNumber) and
                 (st.second.massNumber == mol.atoms[iAtm].massNumber);}
         );

      out << "    " << std::setw(shortFieldWidth) << std::left << 
        (it == atomicReference.end() ? "X" : it->first);
    //out << "    " << std::setw(shortFieldWidth) << std::left << "X";

      out << std::left <<  std::setw(shortFieldWidth-1) <<
             std::fixed << std::setprecision(1) << mol.atoms[iAtm].nucCharge;
      out << std::right;
      out << std::setw(shortFieldWidth) <<
             std::setprecision(5) << std::scientific << mol.atoms[iAtm].atomicMass ;

      out << std::setw(fieldValueWidth-1) << mol.atoms[iAtm].coord[0];
      out << std::setw(fieldValueWidth) << mol.atoms[iAtm].coord[1];
      out << std::setw(fieldValueWidth) << mol.atoms[iAtm].coord[2];

      out << std::endl;
    }
    out << std::endl << BannerEnd << std::endl;

    return out; // Return std::ostream reference

  }; // Molecule::operator<<

}; // namespace ChronusQ

