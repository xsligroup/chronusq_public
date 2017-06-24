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
#include <basisset/reference.hpp>
#include <atom.hpp>
#include <cxxapi/options.hpp>

namespace ChronusQ {

  // Angular momentum map. Only used in this file, so currently
  // no need to place it in a header. This can be remidied if
  // the need arises

  std::unordered_map<std::string,int> LMap = {
    { "S" , 0  },
    { "P" , 1  },
    { "D" , 2  },
    { "F" , 3  },
    { "G" , 4  },
    { "H" , 5  },
    { "I" , 6  },
    { "J" , 7  },
    { "K" , 8  },
    { "L" , 9  },
    { "M" , 10 },
    { "N" , 11 },
    { "O" , 12 },
    { "Q" , 13 },
    { "R" , 14 },
    { "T" , 15 },
    { "U" , 16 },
    { "V" , 17 },
    { "W" , 18 },
    { "X" , 19 },
    { "Y" , 20 },
    { "Z" , 21 }
  };



  /**
   *  \brief Handles G94 format engineering notation
   *
   *  Custom string to double converter that handles
   *  both D and E exponential number notation.
   */
  double ReferenceBasisSet::str2doub(std::string inp) {

    // Convert to upper case then check for 'D'
    // character and replace with 'E'

    unsigned char newchar = 'E';
    std::transform(inp.begin(),inp.end(),inp.begin(),
                   [ &newchar ]( unsigned char c ){
                       if( std::toupper(c) == 'D')
                         c = newchar;
                       return std::toupper(c); });

    // Return string converted to double
    return std::stod(inp);
  }



  /**
   *  \brief Attempts to find a basis file
   *
   *  Populates internal member data for ReferenceBasisSet for the
   *  basis set file object. Terminates program if it cannot find the
   *  file.
   */
  void ReferenceBasisSet::findBasisFile(bool doPrint){

    // Look for basis sets in BASIS_PATH
    std::string fullBasisPath = "/" + basisPath_;
    fullBasisPath.insert(0,BASIS_PATH);

    // Create a file object for the basis set file
    basisFile_ = std::ifstream(fullBasisPath);

    // If the basis set doesn't exist in the default location, check if they
    //   specified a path
    if ( basisFile_.fail() ) {
      basisFile_ = std::ifstream(basisPath_);
      fullBasisPath = basisPath_;
    }

    // Check if file exists
    if(basisFile_.fail()){
      if( doPrint )
      std::cout << "Cannot open basis set " + basisPath_ << std::endl;
      exit(EXIT_FAILURE);
    } else if (doPrint)
      std::cout << "  *** Reading Basis Set from " + fullBasisPath << " ***" << std::endl;

  }; // ReferenceBasisSet::findBasisFile




  /**
   *  \brief Parses the entire basis set file specified in basisPath
   *
   *  Reads a basis set file in the G94 basis format specified by
   *  the basisPath
   */
  void ReferenceBasisSet::parseBasisFile(std::istream& input) {
    std::string readString;
    std::string nameOfAtom;
    std::string shSymb;
    int         contDepth;
    int atomicNumber;
    std::vector<libint2::Shell> tmpShell;
    std::vector<std::vector<double>> tmpCons;

    bool readRec = false;
    bool newRec  = false;
    bool firstRec = true;
    int nEmpty = 0;
    int nComm  = 0;
    int nRec   = 0;

    // Loop over lines in the file
    while(!input.eof()){

      // Obtain file line
      std::getline(input,readString);

      // Strip trailing spaces
      trim(readString);

      if(readString.size() == 0)    nEmpty++; // line empty
      else if(readString[0] == '!') nComm++;  // comment line
      else if(!readString.compare("****")){   // start or end of a shellset

        std::getline(input,readString);

        // end of the record
        if(readString.size() == 0) { nEmpty++; readRec = false; continue;}

        // Set variables to start working on developing a new shell set
        // record
        nRec++;
        readRec = true;
        newRec  = true;

      }

      // We have decided to parse a new shell set record
      if(readRec){

        // Split the line on white space
        std::istringstream iss(readString);
        std::vector<std::string> tokens(std::istream_iterator<std::string>{iss},
          std::istream_iterator<std::string>{});

        if(newRec){

          // If this is not the first record, append the previous record
          // to the shell sets
          if(!firstRec) refShells[atomicNumber] = { tmpShell, tmpCons };


          // Logic for starting a new record

          atomicNumber = atomicNumMap[tokens[0]]; // Grab the atomic number
          newRec = false;
          firstRec = false;
          tmpShell.clear();
          tmpCons.clear();

        } else {
          // Logic for continuing to parse a shell set record

          // Temporary storage for data set
          std::vector<double> exp;
          std::vector<double> contPrimary;
          std::vector<double> contSecondary;

          libint2::svector< double > exp_libint,
                                     contPrimary_libint,
                                     contSecondary_libint;

          // Determine the contraction depth of the shell set
          contDepth = std::stoi(tokens[1]);

          // Angular momentum symbol
          std::string shSymb = tokens[0];

          // Loop over primitives
          for(auto i = 0; i < contDepth; i++) {
            std::getline(input,readString);

            // Split line on white space
            std::istringstream iss2(readString);
            std::vector<std::string> tokens2(
              std::istream_iterator<std::string>{iss2},
              std::istream_iterator<std::string>{}
            );

            // Temporarily  store shell set data
            exp.push_back(str2doub(tokens2[0]));
            contPrimary.push_back(str2doub(tokens2[1]));

            // Append an extra record for SP shells for the "P" orbital
            // set.
            if(!shSymb.compare("SP"))
              contSecondary.push_back(str2doub(tokens2[2]));
          }


          exp_libint.resize( exp.size() );
          contPrimary_libint.resize( contPrimary.size() );
          contSecondary_libint.resize( contSecondary.size() );

          std::copy(exp.begin(), exp.end(), exp_libint.begin() );
          std::copy(contPrimary.begin(), contPrimary.end(), contPrimary_libint.begin() );
          std::copy(contSecondary.begin(), contSecondary.end(), contSecondary_libint.begin() );

          // Create the libint2::Shell object
          if(!shSymb.compare("SP")) {

            // Create 2 shell sets for "SP" shells

            libint2::svector< libint2::Shell::Contraction > s, p;

            s.emplace_back(libint2::Shell::Contraction{0,false,contPrimary_libint});
            p.emplace_back(libint2::Shell::Contraction{1,false,contSecondary_libint});

            // "S" function
            tmpShell.push_back(
              libint2::Shell{ exp_libint, s, {{0,0,0}} }
            );
            tmpCons.push_back(contPrimary);

            // "P" function
            tmpShell.push_back(
              libint2::Shell{ exp_libint, p, {{0,0,0}} }
            );

            tmpCons.push_back(contSecondary);
          } else {

            // Determine L value for the shell
            auto L = LMap[shSymb];

            // Bug in libint2, need to specify S and P functions as cartesian
            // or nonsense ensues. Default spherical for L > 1
            bool doSph = (L > 1);

            // If we're forcing cartesian, force cartesian...
            if(forceCart_) doSph = false;

            libint2::svector< libint2::Shell::Contraction > cnt;
            cnt.emplace_back(libint2::Shell::Contraction{L,doSph,contPrimary_libint});

            // Append temporary shell set
            tmpShell.push_back(
              libint2::Shell{ exp_libint, cnt, {{0,0,0}} }
            );
            tmpCons.push_back(contPrimary);
          } // end append temp record



        } // end shell record parse
      } // end atomic record parse
    } // end file loop

    // Append the last Rec
    //refShells_.push_back(ReferenceShell{atomicNumber,tmpShell,tmpCons});
    refShells[atomicNumber] = { tmpShell, tmpCons };

  }; // ReferenceBasisSet::parseBasisFile

  /**
   *  \brief Generates the proper shell set for a given molecule
   *
   *  \param [in] mol Molecule object for which to generate the shell set
   *  \return         Shell set and coefficients for the given molecule
   */
  std::pair<std::vector<libint2::Shell>,std::vector<std::vector<double>>>
    ReferenceBasisSet::generateShellSet(const Molecule& mol) {

    std::vector<libint2::Shell> shells;
    std::vector<std::vector<double>> cont;

    // loop over atoms in the molecule
    if (not this->nucBasis_) {
      for(auto &atom : mol.atoms){

        if( refShells.find(atom.atomicNumber) == refShells.end() )
          CErr("Cannot find Z=" + std::to_string(atom.atomicNumber) + " in Basis Definition");

        auto &newSh  = refShells[atom.atomicNumber];

        std::vector<libint2::Shell>::iterator shFront =
          shells.insert(shells.end(),newSh.shells.begin(),
            newSh.shells.end());

        cont.insert(cont.end(),newSh.unNormCont.begin(),
          newSh.unNormCont.end());


        std::for_each(shFront,shells.end(),
          [&](libint2::Shell &sh) { sh.move(atom.coord); }
        );

      }
    } else { // nuclear basis
      for (auto &atom : mol.atoms) {

        if ( atom.quantum ) {
          if (atom.atomicNumber != 1)
            CErr("Find non-hydrogen atoms");

          auto &newSh = refShells[atom.atomicNumber];

          std::vector<libint2::Shell>::iterator shFront =
            shells.insert(shells.end(),newSh.shells.begin(),
              newSh.shells.end());

          cont.insert(cont.end(),newSh.unNormCont.begin(),
            newSh.unNormCont.end());

          std::for_each(shFront,shells.end(),
            [&](libint2::Shell &sh) { sh.move(atom.coord); }
          );
        }
      }
    }

    return { shells, cont };

  }; // ReferenceShellSer::generateShellSet


}; // namespace ChronusQ

