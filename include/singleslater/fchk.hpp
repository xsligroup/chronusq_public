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

#include <singleslater.hpp>
#include <cqlinalg.hpp>
#include <util/matout.hpp>

namespace ChronusQ {

  /**
   *  \brief Map used for Gaussian to CQ ordering of angular
   *  momentum functions.
   *
   *  Takes shell value stored on fchk file and returns where that
   *  function will be in CQ ordering.
   *
   *  spherical d: 0,+1,-1,+2,-2 to -2,-1,0,+1,+2
   *  Cartesian d: xx,yy,zz,xy,xz,yz to xx,xy,xz,yy,yz,zz
   *  spherical f: 0,+1,-1,+2,-2,+3,-3 to -3,-2,-1,0,+1,+2,+3
   *  spherical g: 0,+1,-1,+2,-2,+3,-3,+4,-4 to -4,-3,-2,-1,0,+1,+2,+3,+4
   *  spherical h: 0,+1,-1,+2,-2,+3,-3,+4,-4,+5,-5 to
   *              -5,-4,-3,-2,-1,0,+1,+2,+3,+4,+5
   */
  template <typename MatsT, typename IntsT>
  std::unordered_map<int,std::vector<int>> SingleSlater<MatsT,IntsT>::returnAngReorder(){

    std::unordered_map<int,std::vector<int>> angReorder(
      {
        { 0 , {0} },                          // s function
        { 1 , {0,0,0} },                      // p function
        {-1 , {0,0,0,0} },                    // s and p function
        {-2 , {2,2,-1,1,-4} },                // sph. d function
        { 2 , {0,2,3,-2,-2,-1} },             // Cart. d function
        {-3 , {3,3,0,2,-3,1,-6} },            // sph. f function
        {-4 , {4,4,1,3,-2,2,-5,1,-8} },       // sph. g function
        {-5 , {5,5,2,4,-1,3,-4,2,-7,1,-10} }  // sph. h function
      }
    );

    return angReorder;

  }

  /**
   *  \brief Parses the fchk file and overwrites mo1(mo2) 
   *
   *  If using internally-stored correlation-consistent basis
   *  sets in Gaussian calculation, include IOp(3/60=-1).
   *
   **/
  template <typename MatsT, typename IntsT>
  std::vector<int> SingleSlater<MatsT,IntsT>::fchkToCQMO() {

    std::ifstream fchkFile;
    fchkFile.open(fchkFileName);
    std::vector<int> sl;

    if ( not fchkFile.good() ) CErr("Could not find fchkFile. Use -s flag.");

    // dimension of mo1 and mo2
    auto NB =  this->basisSet().nBasis;
    auto NBC = this->nC * NB;
    auto NB2 = NB*NB;
    auto NBC2 = NBC*NBC;

    // Boolean for if fchk entries are found
    bool isBeta=false, readAlpha=false, readBeta=false;
    bool readShell=false;
    // Various integers
    int mo1Counter = 0, mo2Counter = 0, compCounter = 0, maxLfchk;
    // Various double
    double prevValue = 0.0;
    // Double pointer for handling complex mo1
    double* dptr = reinterpret_cast<double*>(this->mo[0].pointer());

    // Parse fchk file
    while( not fchkFile.eof() ) {

      std::string line;
      std::getline(fchkFile,line);

      // Determine position of first and last non-space character
      size_t firstNonSpace = line.find_first_not_of(" ");

      // Skip blank lines
      if( firstNonSpace == std::string::npos ) continue;

      // Strip trailing spaces
      trim_right(line);
      line =
         line.substr(firstNonSpace,line.length()-firstNonSpace);

      // Split the line into tokens, trim spaces
      std::vector<std::string> tokens;
      split(tokens,line,"   ");
      for(auto &X : tokens) { trim(X); }

      // Some sanity checks on fchk data
      // Checking NB
      if ( tokens.size() > 5 ){
        if ( tokens[2] == "basis" and tokens[3] == "functions" ){
          if ( std::stoi(tokens[5]) != NB ){
            std::cout << "      NB from fchk is " << tokens[5] << " and NB in Chronus is " << NB << "\n";
            CErr("Basis functions do not agree between fchk and Chronus!");
          }
        }
      }

      // Checking highest angular momentum
      if ( tokens.size() > 4 ){
        if ( tokens[0] == "Highest" and tokens[1] == "angular" ){
          maxLfchk = std::stoi(tokens[4]);
          if( maxLfchk > 5 ){
            std::cout << "Max L in fchk is " << maxLfchk << ", but greater than " << 5 << " is NYI" << "\n";
            CErr("Angular momentum too high for fchk parser!");
          }
        }
      }

      // Checking size of mo1
      if ( tokens.size() > 5 ){
        if ( tokens[0] == "Alpha" and tokens[1] == "MO" ){
          if ( this->nC == 1 ){
            if ( std::stoi(tokens[5]) != NBC2 ){
              std::cout << "      MO size from fchk is " << tokens[5] << " and MO size in Chronus is " << NBC2 << "\n";
              CErr("MO coefficient size do not agree between fchk and ChronusQ!");
            }
          } else if ( this->nC == 2 ){
            // Factor of two here since real and imaginary are separate in fchk
            if ( std::stoi(tokens[5]) != 2*NBC2 ){
              std::cout << "      MO size from fchk is " << tokens[5] << " but expected size is " << 2*NBC2 << "\n";
              CErr("MO coefficient size on fchk is not compatible!");
            }
          } else{
            CErr("fchk functionality only implemented for 1c and 2c!");
          }
          // Found alpha MO block
          readAlpha=true;
          continue;
        }
      }

      // Check if unrestricted
      if ( tokens.size() > 5 ){
        if ( tokens[0] == "Beta" and tokens[1] == "MO" ){
          if ( this->nC == 2 or this->iCS ) CErr("Beta MOs present in fchk file but Chronus is 2c or restricted");
          isBeta = true;
          readBeta = true;
          readAlpha = false;
          continue;
        }
      }

      // Check if shell list started
      if ( tokens.size() > 4 ){
        if ( tokens[0] == "Shell" and tokens[1] == "types" ){
          readShell = true;
          continue;
        }
      }

      // Check if shell list ended
      if ( tokens.size() > 5 ){
        if ( tokens[2] == "primitives" and tokens[3] == "per" ){
          readShell = false;
          continue;
        }
      }

      // Check if alpha or beta block ended
      if ( tokens.size() > 4 ){
        if ( (tokens[0] == "Orthonormal" and tokens[1] == "basis") or
             (tokens[0] == "Total" and tokens[1] == "SCF")){
          if( this-> nC == 1){
            if( isBeta ) readBeta = false;
            else readAlpha = false;
            continue;
          }
          if( this-> nC == 2){
            readAlpha = false;
            continue;
          }
        }
      }

      // Read in Alpha MO coeffients
      if ( readAlpha ){
        for(int i=0; i<tokens.size(); i++){
          dptr[mo1Counter]=std::stod(tokens[i]);
          mo1Counter=mo1Counter+1;
        }
      }

      // Read in Beta MO coeffients
      if ( readBeta ){
        for(int i=0; i<tokens.size(); i++){
          this->mo[1].pointer()[mo2Counter]=std::stod(tokens[i]);
          mo2Counter=mo2Counter+1;
        }
      }

      // Read in shell types
      if ( readShell ){
        for(int i=0; i<tokens.size(); i++){
          sl.push_back(std::stoi(tokens[i]));
        }
      }

    }// end of eof loop

    if ( not isBeta and this->nC == 1 and not this->iCS ) CErr("Could not find beta MOs on fchk file!");

    fchkFile.close();

    return sl;

  } // SingleSlater<T>::fchkToCQMO()

  /**
   *  \brief Reorders basis functions for a given l to CQ storage  
   *  Example (d orbitals): 0,+1,-1,+2,-2 to -2,-1,0,+1,+2 
   *
   **/
  template <typename MatsT, typename IntsT>
  void SingleSlater<MatsT,IntsT>::reorderAngMO(std::vector<int> sl, MatsT* tmo, int sp) {

    // Dimension of mo
    auto NB = this->basisSet().nBasis;
    auto NBC = this->nC * NB;
    auto NBC2 = NBC*NBC;

    if( this->nC == 4 ) CErr("reorderAngMO NYI for 4c",std::cout);

    // Skip past beta part if 2c
    int skipl = this->nC;

    // Obtain map for angular momentum ordering between Gaussian and CQ
    std::unordered_map<int,std::vector<int>> angReorder = returnAngReorder();

    // Loop over each MO
    for( auto iMO=0; iMO<NBC; iMO++){

      // Loop over each shell entry in sl
      for( auto iSh=0, iAO=0; iSh<sl.size(); iSh++){

        // Check that shell entry is implemented
        bool contains = angReorder.find(sl[iSh]) != angReorder.end();
        if( not contains ){
          std::cout << "Cannot find shell entry: " << sl[iSh] << std::endl;
          CErr("Shell value for FCHKMO NYI!",std::cout);
        }

        const std::vector<int>& reorder = angReorder[sl[iSh]];

        // Loop over each function per shell
        for( auto iML=0; iML<reorder.size(); iML++){

          this->mo[sp].pointer()[iMO*NBC+iAO+skipl*(iML+reorder[iML])] = tmo[iMO*NBC+iAO+iML*skipl];
          // Beta part
          if( this->nC==2 ) 
            this->mo[sp].pointer()[iMO*NBC+1+iAO+skipl*(iML+reorder[iML])] = tmo[iMO*NBC+1+iAO+iML*skipl];

        }

        iAO=iAO+reorder.size()*skipl;

      }

    }

  } // SingleSlater<T>::reorderAngMO()


  /**
   *  \brief Reorders spin components of basis functions (needed for 2c only) 
   *  Example (1 real 2c MO with 2 basis functions): 
   *  alpha1,beta1,alpha2,beta2 to alpha1,alpha2,beta1,beta2
   *
   **/
  template <typename MatsT, typename IntsT>
  void SingleSlater<MatsT,IntsT>::reorderSpinMO() {

    // Dimension of mo1
    auto NB = this->basisSet().nBasis;
    auto NBC = this->nC * NB;
    auto NBC2 = NBC*NBC;

    // Counters for spin reordering
    int acount = 0, bcount = NB - 1;

    // mo scratch
    MatsT* mo1tmp = CQMemManager::get().malloc<MatsT>(NBC2);
    SetMat('N',NBC,NBC,MatsT(1.),this->mo[0].pointer(),NBC,mo1tmp,NBC);

    // Loop over mo1
    for( int i=0; i<NBC2; i++){

      // Reset acount and bcount for each MO
      if( i % NBC == 0 ){
        acount = 0;
        bcount = NB - 1;
      }

      // If alpha
      if( i % 2 == 0 ){
        this->mo[0].pointer()[i-acount] = mo1tmp[i];
        acount = acount + 1;
      }else{ // If beta
        this->mo[0].pointer()[i+bcount] = mo1tmp[i];
        bcount = bcount - 1;
      }
    }

    CQMemManager::get().free(mo1tmp);

  } // SingleSlater<T>::reorderSpinMO()

}
