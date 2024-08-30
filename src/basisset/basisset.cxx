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

#include <basisset/basisset_def.hpp>
#include <basisset/reference.hpp>
#include <cxxapi/output.hpp>

#include <libcint.hpp>

#include <util/matout.hpp>
#include <cqlinalg/blas1.hpp>
#include <cqlinalg/blas3.hpp>
#include <algorithm>

namespace ChronusQ {



  /**
   *  Path / Molecule BasisSet constructor 
   *
   *  Constructs a BasisSet object and populates member data given
   *  a basis set name and a Molecule object
   *
   *  \param [in] basisName Either a basis file in G94 format
   *                        or a known keyword which maps to a basis
   *                        file
   *  \param [in] mol       Molecule for which to construct the BasisSet
   */
  BasisSet::BasisSet(std::string _basisName, std::string _basisDef,
    bool doDef, const Molecule &mol, BASIS_FUNCTION_TYPE _basisType, bool
    _forceCart, bool doPrint, bool _nucBasis) {

    basisType = _basisType;
    forceCart = _forceCart;
    basisName = _basisName;
    basisDef  = _basisDef;
    inputDef  = doDef;
    nucBasis  = _nucBasis;


    std::string uppercase(basisName);

    // Convert the passed string to uppercase
    std::transform(basisName.begin(),basisName.end(),uppercase.begin(),
      [](unsigned char c){ return std::toupper(c); });  

    // Possibly find appropriate basis file for keyword
    if( basisKeyword.find(uppercase) != basisKeyword.end() ) {
      basisName = basisKeyword[uppercase];
    }

    // Generate the reference basis set of that keyword
    ReferenceBasisSet ref(basisName, basisDef, inputDef, _forceCart, doPrint, _nucBasis);

    // Update appropriate shell set and coefficients for the Molecule
    // object
    std::tie(shells,unNormCont) = std::move(ref.generateShellSet(mol));

    // Obtain a copy of the basis centers
    std::for_each(mol.atoms.begin(),mol.atoms.end(),
      [&]( const Atom &at ){ centers.emplace_back(at.coord); }
    );

    // Update the BasisSt member data
    update();
    
  }; // BasisSet::BasisSet 


  /**
   *  Updates the member data of a BasisSet object.
   *
   *  Computes nShell, nBasis and nPrimitive. Determines  maxL and maxPrim.
   *  Constructs the shell mappings to basis function number and center number
   */ 
  void BasisSet::update(bool computeShellPairs) {

    // Compute the number of shells
    nShell = shells.size();

    // Compute the number of Basis functions and primitives
    nBasis = std::accumulate(shells.begin(),shells.end(),0,
      [](size_t init, libint2::Shell &sh) -> size_t {
        return init + sh.size();
      }
    );

    // Uncontract Basis functions and get number of primitives
    uncontractShells();
    nPrimitive = std::accumulate(primitives.begin(),primitives.end(),0,
      [](size_t init, auto &sh) -> size_t {
        return init + sh.first.size();
      }
    );

    // Determine the maximum angular momentum
    maxL = std::max_element(shells.begin(), shells.end(),
      [](libint2::Shell &sh1, libint2::Shell &sh2){
        return sh1.contr[0].l < sh2.contr[0].l;
      }
    )->contr[0].l;

    // Determine the maximum contraction depth
    maxPrim = std::max_element(shells.begin(), shells.end(),
      [](libint2::Shell &sh1, libint2::Shell &sh2){
        return sh1.alpha.size() < sh2.alpha.size();
      }
    )->alpha.size();


    // Reset maps
    mapSh2Cen.clear();
    mapSh2Bf.clear();
    mapCen2BfSt.clear();

    // Create basis maps
    // Maps Sh # -> BF #
    // Maps Sh # -> Center index
    size_t run(0);
    std::for_each(shells.begin(),shells.end(),
      [&](libint2::Shell &sh) {
        auto cenIt = std::find_if(centers.begin(),centers.end(),
          [&](cart_t &cen){ 
            //return sh.O == cen; 
            return std::abs(sh.O[0] - cen[0]) < 1e-5 and  
                   std::abs(sh.O[1] - cen[1]) < 1e-5 and  
                   std::abs(sh.O[2] - cen[2]) < 1e-5; 
          }
        );
        mapSh2Cen.emplace_back(std::distance(centers.begin(),cenIt));

        mapSh2Bf.emplace_back(run);
        run += sh.size(); 
        
      }
    );


    for(auto iAtm = 0; iAtm < centers.size(); iAtm++) {

      auto it = std::find_if(mapSh2Cen.begin(),mapSh2Cen.end(),
                  [&](size_t x){ return x == iAtm; });

      if (nucBasis and it == mapSh2Cen.end())
        mapCen2BfSt.emplace_back(0);
      else {
        size_t firstShell = std::distance(mapSh2Cen.begin(),it);
        mapCen2BfSt.emplace_back(mapSh2Bf[firstShell]);
      }

    }

    //std::cout << "mapCen2BfSt is " << std::endl;
    //for(size_t ii = 0; ii < mapCen2BfSt.size(); ii++)
    //  std::cout << mapCen2BfSt[ii] << "  ";
    //std::cout << std::endl;




    // Compute shell pair data
    if (computeShellPairs) {
      const auto max_engine_precision = std::numeric_limits<double>::epsilon();
      const auto ln_prec = std::log(max_engine_precision);
      shellData.computeShellPairs(shells,mapSh2Cen,maxPrim,maxL,1e-12,
          ln_prec);
    }

  }; // BasisSet::update

  /**
   *  Updates the member data of a BasisSet object assuming one or more nuclei
   *  have been moved
   *
   */
  void BasisSet::updateNuclearCoordinates(const Molecule &mol) {

    // Early return if this basis set doesn't actually contain anything
    if( basisName == "" and basisDef == "" and shells.size() == 0 )
      return;

    // Generate the reference basis set of that keyword
    ReferenceBasisSet ref(basisName, basisDef, inputDef, forceCart, false, nucBasis);

    // update appropriate shell set and coefficients for the Molecule
    // object
    std::tie(shells,unNormCont) = std::move(ref.generateShellSet(mol));

    // clear the old list of centers
    centers.clear();

    // Update the copy of the basis centers
    std::for_each(mol.atoms.begin(),mol.atoms.end(),
                  [&]( const Atom &at ){
                    centers.emplace_back(at.coord);
                  }
    );

    // Update the BasisSet member data
    update();

  };


  /**
   *  Outputs relevant information for the BasisSet object
   *  to a specified output.
   *
   *  \param [in/out] out Ouput device
   *  \paral [in]     mol BasisSet object to output.
   */ 
  std::ostream& operator<<(std::ostream &out, const BasisSet &basis) {

    out << std::endl;
    out << "Basis Set Information:" << std::endl << BannerTop
        << std::endl << std::endl;
    
    out << std::left;

    out << "  " << std::setw(25) << "NBasis" << basis.nBasis 
        << std::endl;
    out << "  " << std::setw(25) << "NPrimitive" << basis.nPrimitive 
        << std::endl;
    out << "  " << std::setw(25) << "NShell" << basis.nShell 
        << std::endl;
    out << "  " << std::setw(25) << "Max Primitive" << basis.maxPrim 
        << std::endl;
    out << "  " << std::setw(25) << "Max L" << basis.maxL << std::endl;

    out << std::endl << std::endl;

    out << "  " << "Shell Information:" << std::endl << bannerMid << std::endl;

    out << "  " << "  " << std::left;
    out << std::setw(5) << "#" ;
    out << std::setw(5) << "L" ;
    out << std::setw(15) << std::right << "Exponents";
    out << std::setw(30) << std::right << "Normalized Contraction";
  //out << std::setw(30) << std::right << "Unnormalized Contraction";
    out << std::endl << std::endl;
    
    for(auto iShell = 0; iShell < basis.nShell; iShell++){
      out << "  " << "  " << std::left << std::setprecision(4) 
          << std::scientific;
      out << std::setw(5) << iShell ;
      out << std::setw(5) << basis.shells[iShell].contr[0].l << std::right;
      out << std::setw(15) << basis.shells[iShell].alpha[0];
      out << std::setw(30) << basis.shells[iShell].contr[0].coeff[0];
    //out << std::setw(30) << basis.unNormCont[iShell][0];

      out << std::endl;
      for(auto iPrim = 1; iPrim < basis.shells[iShell].alpha.size(); iPrim++){
        out << "  " << "  " << std::setw(5+5) << " ";
        out << std::setw(15) << basis.shells[iShell].alpha[iPrim];
        out << std::setw(30) << basis.shells[iShell].contr[0].coeff[iPrim];
      //out << std::setw(30) << basis.unNormCont[iShell][iPrim];

        out << std::endl;
      }
      out << std::endl;
    }
    out << std::endl << BannerEnd << std::endl;

    return out; // return std::ostream reference

  }; // BasisSet::operator<<

  /**
   *  \brief Return the uncontracted shell set of the current
   *  contracted shell set
   */ 
  void BasisSet::uncontractShells() {

    primitives.clear();

    size_t count = 0;

    for(auto &shell : shells) {// Loop over shells
      for(auto &con : shell.contr) {// Loop over contractions
        for(auto &a : shell.alpha) {// Loop over primitives
          libint2::Shell newShell{
            { a },
            { { con.l, con.pure, { 1.0 } } },
            { { shell.O[0], shell.O[1], shell.O[2] } }
          };
          if (primitives.find(newShell) == primitives.end()) {
            primitives[newShell] = count;
            count += newShell.size();
          }
        }
      }
    }

  }; // BasisSet::uncontractShells


  /**
   * @brief Return a new BasisSet object with only the primitives
   *
   */
  BasisSet BasisSet::uncontractBasis() const {


    // Copy basis
    BasisSet newBasis(*this);

    newBasis.shells.clear();
    for (const auto &prim : primitives) {
      newBasis.shells.push_back(prim.first);
    }

    std::sort(newBasis.shells.begin(), newBasis.shells.end(),
      [this](const libint2::Shell &a, const libint2::Shell &b)->bool {
        return primitives.at(a) < primitives.at(b);
    });

    newBasis.update();

    return newBasis;

  };


  /**
   * @brief Generates a new BasisSet object,
   *        where basis functions sharing the same set of primitives
   *        are grouped together in one libint2::Shell
   *
   */
  BasisSet BasisSet::groupGeneralContractionBasis() const {


    // Copy basis
    BasisSet GCBasis(*this);

    std::vector<libint2::Shell> shells;

    shells.push_back(*GCBasis.shells.begin());

    for (auto it = ++GCBasis.shells.begin(); it != GCBasis.shells.end(); it++) {

      if (shells.back().O == it->O and
          shells.back().alpha == it->alpha and
          shells.back().contr[0].l == it->contr[0].l) {
        // This shell has the same origin, set of exponents, and angular momentum
        // as the previous shell, merge to previous contraction coefficients.

        shells.back().contr.push_back(it->contr[0]);

      } else {
        // A different shell

        shells.push_back(*it);
      }

    }

    GCBasis.shells = shells;

    GCBasis.update(false); // Update, but do not compute shell pair data.
                           // Computing shell pair data does not work with GC

    return GCBasis;

  };


  /**
   * @brief computes the necessary length of the env vector for libcint computation
   *
   */
  size_t BasisSet::getLibcintEnvLength(const Molecule &mol) const {
    return PTR_ENV_START +
        mol.nAtoms * 4 +
        std::accumulate(shells.begin(),
                        shells.end(),
                        0,
                        [](const size_t &count, const libint2::Shell &sh) {
                          return count + sh.alpha.size() * (1 + sh.contr.size());
                        });
  }


  /**
   * @brief sets up the atm, bas, and env vectors for libcint computation
   */
  void BasisSet::setLibcintEnv(const Molecule &mol, int *atm, int *bas, double *env,
                               bool finiteWidthNuc) const {

    double sNorm;

    int nAtoms = mol.nAtoms;
    int nShells = nShell;
    int off = PTR_ENV_START; // = 20

    for(int iAtom = 0; iAtom < nAtoms; iAtom++) {

      atm(CHARGE_OF, iAtom) = mol.atoms[iAtom].integerNucCharge();

      atm(PTR_COORD, iAtom) = off;
      env[off++] = mol.atoms[iAtom].coord[0]; // x (Bohr)
      env[off++] = mol.atoms[iAtom].coord[1]; // y (Bohr)
      env[off++] = mol.atoms[iAtom].coord[2]; // z (Bohr)

      if (finiteWidthNuc) {
        if (mol.atoms[iAtom].fractionalNucCharge()) {
          CErr("Libcint cannot handle finite nucleus + fractional nuclear charge.");
        }
        atm(NUC_MOD_OF, iAtom) = GAUSSIAN_NUC;
        atm(PTR_ZETA  , iAtom) = off;
        env[off++] = mol.chargeDist[iAtom].alpha[0];
      } else {
        if (mol.atoms[iAtom].fractionalNucCharge()) {
          atm(NUC_MOD_OF     , iAtom) = FRAC_CHARGE_NUC;
          atm(PTR_FRAC_CHARGE, iAtom) = off;
          env[off++] = mol.atoms[iAtom].nucCharge;
        } else {
          atm(NUC_MOD_OF, iAtom) = POINT_NUC;
        }
      }

    }

    for (size_t iAtom : mol.atomsQ) {
      atm(CHARGE_OF , iAtom) = 0;
      atm(NUC_MOD_OF, iAtom) = POINT_NUC;
    }

    for(int iShell = 0; iShell < nShells; iShell++) {

      int nContr = shells[iShell].contr.size();

      bas(ATOM_OF , iShell)  = std::distance(mol.atoms.begin(),
          std::find_if( mol.atoms.begin(),
                        mol.atoms.end(),
                        [this, iShell](const Atom &a){
                          return a.coord == shells[iShell].O;
                        }));
      if (bas(ATOM_OF , iShell) >= nAtoms) {
        CErr("setLibcintEnv: Cannot find corresponding atom for shell");
      }
      bas(ANG_OF  , iShell)  = shells[iShell].contr[0].l;
      bas(NPRIM_OF, iShell)  = shells[iShell].alpha.size();
      bas(NCTR_OF , iShell)  = nContr;
      bas(KAPPA_OF, iShell)  = 0;
      bas(PTR_EXP , iShell)  = off;

      for(int iPrim=0; iPrim < shells[iShell].alpha.size(); iPrim++)
        env[off++] = shells[iShell].alpha[iPrim];

      bas(PTR_COEFF, iShell) = off;

      // Spherical GTO normalization constant missing in Libcint
      sNorm = 2.0*std::sqrt(M_PI)/std::sqrt(2.0*shells[iShell].contr[0].l+1.0);
      for (size_t i = 0; i < nContr; i++) {
        for(int iCoeff=0; iCoeff<shells[iShell].alpha.size(); iCoeff++){
          env[off++] = shells[iShell].contr[i].coeff[iCoeff]*sNorm;
        }
      }

    }
  }


  template <> 
  void BasisSet::makeMapPrim2Cont(const double *SUn, double *MAP) const {

    memset(MAP,0,nPrimitive * nBasis * sizeof(double));

    double *rA = MAP;

    size_t cumeNBf = 0;
    // Compute the unnormalized mapping
    for(size_t iSh = 0; iSh < nShell; iSh++) {

      size_t nContr= shells[iSh].contr.size();
      size_t nPrim = shells[iSh].alpha.size();

      for(size_t iC = 0; iC < nContr; iC++) {// Loop over contractions
        size_t nBf   = shells[iSh].contr[iC].size();
        for(size_t iP = 0; iP < nPrim; iP++) {// Loop over primitives
          libint2::Shell prim{
            { shells[iSh].alpha[iP] },
            { {shells[iSh].contr[iC].l, shells[iSh].contr[iC].pure, { 1.0 } } },
            { { shells[iSh].O[0], shells[iSh].O[1], shells[iSh].O[2] } }
          };
          size_t primIdx = primitives.at(prim);
          for(size_t iB = 0; iB < nBf; iB++) {
            MAP[cumeNBf+iB + (primIdx+iB)*nBasis] = unNormCont[iSh][iP];
          }
        }
        cumeNBf += nBf;
      }

    } // loop over shells

 
    // Compute SUn * MAP
    double *SCR = CQMemManager::get().malloc<double>(nBasis*nPrimitive);
    blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::Trans,nPrimitive,nBasis,nPrimitive,static_cast<double>(1.),SUn,nPrimitive,
      MAP,nBasis,static_cast<double>(0.),SCR,nPrimitive);


    libint2::Engine engine(libint2::Operator::overlap,maxPrim, maxL, 0);
    engine.set_precision(0.);

    const auto& buf_vec = engine.results();

    // Loop over the shells and the basis functions,
    // calculate the diagonals of contracted overlap integrals
    size_t n1;
    for(size_t s1 = 0,  Itot = 0 ; s1<nShell; s1++, Itot += n1 ) {
    // Itot is the beginning basis function index of current shell

      // Get the size of the shell
      n1 = shells[s1].size();  

      // This computes the diagonal shell block for s1
      engine.compute(shells[s1],shells[s1]);

      if( buf_vec[0] == nullptr ) continue;
      const double* buff = buf_vec[0];

      for(size_t i = 0; i < n1; i++) {
      // i loop over the function in the current shell             

        double fact = blas::dot(nPrimitive,MAP + (i+Itot) ,nBasis,
                                                   SCR + (i+Itot)*nPrimitive,1);
        // i+I is the basis function index 
        double alpha = buff[i + i*n1];
        // alpha is the diagonal element of overlap of basis function i+Itot
         
        blas::scal(nPrimitive,std::sqrt(alpha)/std::sqrt(fact),
          MAP + (i+Itot),nBasis);

      } // for size_t i = 0
    } // for size_t s1 = 0
    CQMemManager::get().free(SCR);

  };  // BasisSet::makeMapPrim2Cont


// make tempelate for makeMapPrim2Cont for real and complex
  template<>  
  void BasisSet::makeMapPrim2Cont(const dcomplex *SUn, dcomplex *MAP) const {

    memset(MAP,0,nPrimitive * nBasis * sizeof(dcomplex));

    dcomplex *rA = MAP;

    size_t cumeNBf = 0;
    // Compute the unnormalized mapping
    for(size_t iSh = 0; iSh < nShell; iSh++) {

      size_t nContr= shells[iSh].contr.size();
      size_t nPrim = shells[iSh].alpha.size();

      for(size_t iC = 0; iC < nContr; iC++) {// Loop over contractions
        size_t nBf   = shells[iSh].contr[iC].size();
        for(size_t iP = 0; iP < nPrim; iP++) {// Loop over primitives
          libint2::Shell prim{
            { shells[iSh].alpha[iP] },
            { {shells[iSh].contr[iC].l, shells[iSh].contr[iC].pure, { 1.0 } } },
            { { shells[iSh].O[0], shells[iSh].O[1], shells[iSh].O[2] } }
          };
          size_t primIdx = primitives.at(prim);
          for(size_t iB = 0; iB < nBf; iB++) {
            MAP[cumeNBf+iB + (primIdx+iB)*nBasis] = unNormCont[iSh][iP];
          }
        }
        cumeNBf += nBf;
      }

    } // loop over shells

 
    // Compute SUn * MAP
    dcomplex *SCR = CQMemManager::get().malloc<dcomplex>(nBasis*nPrimitive);
    blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::Trans, 
      nPrimitive,nBasis,nPrimitive,static_cast<dcomplex>(1.),SUn,nPrimitive,
      MAP,nBasis,static_cast<dcomplex>(0.),SCR,nPrimitive);


    libint2::Engine engine(libint2::Operator::overlap,maxPrim, maxL, 0);
    engine.set_precision(0.);

    const auto& buf_vec = engine.results();

    // Loop over the shells and the basis functions,
    // calculate the diagonals of contracted overlap integrals
    size_t n1;
    for(size_t s1 = 0,  Itot = 0 ; s1<nShell; s1++, Itot += n1 ) {
    // Itot is the beginning basis function index of current shell

      // Get the size of the shell
      n1 = shells[s1].size();  

      // This computes the diagonal shell block for s1
      engine.compute(shells[s1],shells[s1]);

      if( buf_vec[0] == nullptr ) continue;
      const double* buff = buf_vec[0];

      for(size_t i = 0; i < n1; i++) {
      // i loop over the function in the current shell             

        dcomplex fact = blas::dot(nPrimitive,MAP + (i+Itot) ,nBasis,
                                  SCR + (i+Itot)*nPrimitive,1);
        // i+I is the basis function index 
        dcomplex alpha = buff[i + i*n1];
        // alpha is the diagonal element of overlap of basis function i+Itot
         
        blas::scal(nPrimitive,std::sqrt(alpha)/std::sqrt(fact),
          MAP + (i+Itot),nBasis);

      } // for size_t i = 0
    } // for size_t s1 = 0
    CQMemManager::get().free(SCR);

  };  // BasisSet::makeMapPrim2Cont<dcomplex>




  void ShellPairData::computeShellPairs(std::vector<libint2::Shell> &shs, 
      std::vector<size_t> &mapSh2Cen, size_t maxNPrim, size_t maxL, 
      double shell_thresh, double ln_prec) {


    // Reset the vars
    sigShellPair.clear();
    shData.clear();

    const size_t nShell = shs.size();
    
    // Setup overlap engine
    libint2::Engine engine(libint2::Operator::overlap, maxNPrim, maxL,0);
    const auto& buf = engine.results();

    for(size_t s1 = 0ul, s12 = 0ul; s1 < nShell; s1++) {

      if( sigShellPair.find(s1) == sigShellPair.end() )
        sigShellPair.insert( std::make_pair(s1, std::vector<size_t>()) );

      size_t n1 = shs[s1].size();

      // Compute significant list
      for(size_t s2 = 0ul; s2 <= s1; s2++, s12++) {

        bool sameCenter  = (mapSh2Cen[s1] == mapSh2Cen[s2]);
        bool significant = sameCenter;

        if( not sameCenter ) {

          size_t n2 = shs[s2].size();
          engine.compute(shs[s1],shs[s2]);
          double norm = blas::nrm2(n1*n2,const_cast<double*>(buf[0]),1);
          significant = (norm >= shell_thresh);

        }

        if( significant ) sigShellPair[s1].emplace_back(s2);

      }

      // Sort list in increasing order
      std::sort(sigShellPair[s1].begin(), sigShellPair[s1].end());

      // Compute the shell pair data
      shData.emplace_back();
      for(const auto& s2 : sigShellPair[s1])
        shData.back().emplace_back(
          std::make_shared<libint2::ShellPair>(shs[s1],shs[s2],ln_prec)
        );


    }

  }
















  std::vector<std::vector<std::array<int,3>>> cart_ang_list;
  void pop_cart_ang_list() {
    // generate angular momentum list, can only be called once. 
    int k,xx,yy,x,y,z;
    // TangDD: expand cart_ang_list for Bottomup ERIs
    for (  k = 0 ; k <= LIBINT2_MAX_AM + LIBINT2_MAX_AM ; k++ ){
      cart_ang_list.emplace_back();  //loop over possible angular momentum
      for (  xx=0 ; xx<k+1 ; xx++ ){ 
        x = k -xx ; 
        for ( yy = 0 ; yy<xx+1 ; yy++ ){ 
          y = xx-yy; 
          z = k - x - y;
          cart_ang_list[k].push_back({x,y,z});
        }
      }
    }
  };


// sperical to cartesian transform matrix elements

//math functions, like factorial, double factorial and bionomials
//--------------------------------------//
// factorial function:  t! = 1*2*3...*t //
//--------------------------------------//
double factorial(int t){
  int i;
  double tmp = 1.0;
  if (t<0 ) std::cout<<"Factorial (t!) only defined on domain t in [0,inf)"<<std::endl;
  if (t==0) return 1.0;
  else {
    for(i=1;i<=t;i++) tmp *= i;
    return tmp;
  }
} //factorial

//-----------------------------------------------//
// double factorial:  (2t-1)!! = 1*3*5*..*(2t-1) //
//-----------------------------------------------//
double doubleFact(int t){
  int i;
  double tmp = 1.0;
  if (t<0 ) std::cout<<"Double factorial (t!!) only defined on domain t in [0,inf)" 
                     <<std::endl;
  if (t==0) return 1.0;
  else  {
    for(i=1;i<=t;i++) tmp *= (2*i-1);
    return tmp;
  }
} // doubleFact

//---------------------------------------------------------//
// polynomial coeff:   (x+a)^l= sum(coeff * x^i * a^(l-i)) //
//---------------------------------------------------------//
double polyCoeff(int l, int i) {

  double dummy;
  if (l >= i) {
    return factorial(l)/( factorial(i)*factorial(l-i) );
  } else { 
    CErr("polyCoeff error");
    return dummy;
  }
}


std::complex <double> cart2sphCoeff(int L, int m, int lx, int ly, int lz){

//calculate the cartesian to spherical transformation coefficient.

  int Ltotal;
  dcomplex coeff(0.0);
  Ltotal = lx+ly+lz;
  double tmp = 0.0;
  if (L!=Ltotal) {
    return  coeff;
  }
  double j;
  j = (double(lx+ly)-std::abs(double(m)))/2;
  if (fmod(j,1)>0) {
    return coeff;
  }
  dcomplex sumval(0.0);
  dcomplex ttmmpp,sumsumval;
  dcomplex pref,absmchooselxm2k,ichoosej;
  int i,k;
  if (Ltotal == L) {
  pref = sqrt(factorial(lx*2)*factorial(2*ly)
          *factorial(2*lz)*factorial(L)
          *factorial(L-std::abs(m))
     /(factorial(2*L)*factorial(lx)*factorial(ly)
      *factorial(lz)*factorial(L+std::abs(m))))
     /(factorial(L)*pow(2,L));
  
  i = 0;
  
  while (i<=double((L-std::abs(m))/2) ) {
    sumsumval = 0.0;
    for ( k = 0 ; k <= j ; k++ ) {
      if (m>=0) {
        ttmmpp = double(std::abs(m)-lx+2*k)/2;
      }
      else {
        ttmmpp = -double(std::abs(m)-lx+2*k)/2;
      }
      
      if ((std::abs(m)>=(lx-2*k))&&((lx-2*k)>=0)) {
        absmchooselxm2k =polyCoeff(std::abs(m),lx-2*k);
      }
      else {
        absmchooselxm2k = 0.0;
      }
      sumsumval = sumsumval + polyCoeff(j,k)
                  *absmchooselxm2k*pow(-1.0,ttmmpp);
    }
    if (i<j||(j<0)) {
       ichoosej = 0.0;
    }
    else {
      ichoosej = polyCoeff(i,j);
    }
    sumval = sumval + polyCoeff(L,i)*ichoosej*pow(-1,i)
                *factorial(2*L-2*i)/
                (factorial(L-std::abs(m)-2*i))*sumsumval;
    i = i + 1;
  }
  coeff = pref * sumval;
  return coeff;
  }

  return coeff;
} // cart2sphCoeff   

std::vector<std::vector<double>> car2sph_matrix;
void pop_car2sph_matrix() {
  int l[3];
  int m;
  int carsize; //cartesian size
  double scalecoeff;

  for ( int L = 0 ; L <= LIBINT2_MAX_AM ; L++ ) {

    carsize = (L+1)*(L+2)/2;
    car2sph_matrix.emplace_back((2*L+1)*carsize, 0.0);
    if ( L==0 ){
      car2sph_matrix[L][0]=1.0;
    } 
    else if (L==1) {
      car2sph_matrix[L][0*3+0] = 1.0;
      car2sph_matrix[L][1*3+1] = 1.0;
      car2sph_matrix[L][2*3+2] = 1.0;
    }
    else {
      for ( int p=0 ; p < L+1 ; p++ )
      for ( int q=0 ; q < carsize ; q++ ){
        for ( int k=0 ; k<3 ; k++ ) l[k]=cart_ang_list[L][q][k];
        m = -L+p;    
        if ( m<0 ) {
          auto cplxcoeff = cart2sphCoeff(L,m,l[0],l[1],l[2]);   // complex coefficient
          
          car2sph_matrix[L][p*carsize+q] = sqrt(2.0)*(-cplxcoeff.imag());

          car2sph_matrix[L][(2*L-p)*carsize+q] = sqrt(2.0)*cplxcoeff.real();

          double scalecoeff = sqrt(doubleFact(L)/
          (doubleFact(l[0])*doubleFact(l[1])*doubleFact(l[2])));
        car2sph_matrix[L][p*carsize+q] = car2sph_matrix[L][p*carsize+q] *scalecoeff;
        car2sph_matrix[L][(2*L-p)*carsize+q] = car2sph_matrix[L][(2*L-p)
                                               *carsize+q]*scalecoeff; 
      } else if ( m ==0 ) {
        car2sph_matrix[L][p*carsize+q] = cart2sphCoeff(L,m,l[0],l[1],l[2]).real();
        scalecoeff = sqrt(doubleFact(L)/
          (doubleFact(l[0])*doubleFact(l[1])*doubleFact(l[2])));
        car2sph_matrix[L][p*carsize+q] = car2sph_matrix[L][p*carsize+q]*scalecoeff;   
      } 
    }
  }

 } //loop over angular momentum

} //pop_car2sph_matrix()

void cart2sph_transform( int l_i, int l_j, std::vector<double> &shell_element_sph, 
  std::vector<double> &shell_element_cart) {

  int cart_i = (l_i+1)*(l_i+2)/2;
  int cart_j = (l_j+1)*(l_j+2)/2;
  int cartsize = cart_i*cart_j;
  int sphsize  = (2*l_i+1)*(2*l_j+1);
  double tempVal;

  if (sphsize != shell_element_sph.size() )
    std::cout<<"spherical dimension doesn't match"<<std::endl;
  if (cartsize != shell_element_cart.size() )
    std::cout<<"cartesian dimension doesn't match"<<std::endl;
 
  for( int i = 0 ; i<2*l_i+1 ; i++  ) {
    for( int j = 0 ; j<2*l_j+1 ; j++ ) {
      tempVal = 0.0;
      for( int p = 0 ; p<cart_i ; p++ ) {
        for( int q = 0 ; q<cart_j ; q++ ) {
          tempVal += car2sph_matrix[l_i][i*cart_i+p]*car2sph_matrix[l_j][j*cart_j+q]
                     *shell_element_cart[p*cart_j+q];
        }
      }

    if (std::abs(tempVal)<1.0e-15) tempVal=0.0;
    shell_element_sph[i*(2*l_j+1)+j] = tempVal;
    }
  }   
} //cart2sph_transform   

 /** 
  *  \brief Pre calculate the shell pair data
  *
  *  \param [in] shells shell set for integral evaluation
  *  \param [in] ln_prec something related with sceening
  *
  *  returns     A vector of libint2 shellpairs
  */

  std::vector<libint2::ShellPair> genShellPairs( std::vector<libint2::Shell> &shells, double ln_prec ) {

    std::vector<libint2::ShellPair> shellpairs; 
    for ( int s1 = 0, s12 = 0 ; s1 < shells.size() ; s1++ )
    for ( int s2 = 0 ; s2 <= s1 ; s2++, s12++ ) {

      libint2::ShellPair tmp_shellpair;
      tmp_shellpair.init( shells[s1], shells[s2], ln_prec );   
      shellpairs.push_back(tmp_shellpair);
       
    } // for ( int s2 = 0 )
    return shellpairs;

  } // genShellPairs 
 
 /** 
  *  \brief Pre calculate the shell pair data with angular momentum reordering
  *
  *  \param [in] shells shell set for integral evaluation
  *  \param [in] ln_prec something related with sceening
  *
  *  returns     A vector of libint2 shellpairs
  */
  std::vector<libint2::ShellPair> genOrderedShellPairs( std::vector<libint2::Shell> &shells, double ln_prec ) {

    std::vector<libint2::ShellPair> shellpairs; 
    for ( int s1 = 0, s12 = 0 ; s1 < shells.size() ; s1++ )
    for ( int s2 = 0 ; s2 <= s1 ; s2++, s12++) {
      libint2::ShellPair tmp_shellpair;

      if ( shells[s1].contr[0].l >= shells[s2].contr[0].l ) 
        tmp_shellpair.init( shells[s1], shells[s2], ln_prec );   
      else 
        tmp_shellpair.init( shells[s2], shells[s1], ln_prec );
      
      shellpairs.push_back(tmp_shellpair);
       
    } // for ( int s2 = 0 )
    return shellpairs;

  } // genorderedShellPairs 











void cart2sph_complex_transform( int l_i, int l_j, 
  std::vector<dcomplex> &shell_element_sph, std::vector<dcomplex> &shell_element_cart) {

  int cart_i = (l_i+1)*(l_i+2)/2;
  int cart_j = (l_j+1)*(l_j+2)/2;
  int cartsize = cart_i*cart_j;
  int sphsize  = (2*l_i+1)*(2*l_j+1);
  dcomplex tempVal;

  if (sphsize != shell_element_sph.size() )
    std::cout<<"spherical dimension doesn't match"<<std::endl;
  if (cartsize != shell_element_cart.size() )
    std::cout<<"cartesian dimension doesn't match"<<std::endl;
 
  for( int i = 0 ; i<2*l_i+1 ; i++  ) {
    for( int j = 0 ; j<2*l_j+1 ; j++ ) {
      tempVal = 0.0;
      for( int p = 0 ; p<cart_i ; p++ ) {
        for( int q = 0 ; q<cart_j ; q++ ) {
          tempVal += car2sph_matrix[l_i][i*cart_i+p]*car2sph_matrix[l_j][j*cart_j+q]
                     *shell_element_cart[p*cart_j+q];
        }
      }

    if (std::abs(tempVal)<1.0e-15) tempVal=0.0;
    shell_element_sph[i*(2*l_j+1)+j] = tempVal;
    }
  }   
} //cart2sph_complex_transform   

void cart2sph_2e_transform( int l_i, int l_j, int l_k, int l_l, 
  std::vector<double> &shell_element_sph, std::vector<double> &shell_element_cart) {

  int cart_i = (l_i+1)*(l_i+2)/2;
  int cart_j = (l_j+1)*(l_j+2)/2;
  int cart_k = (l_k+1)*(l_k+2)/2;
  int cart_l = (l_l+1)*(l_l+2)/2;

  int cartsize = cart_i*cart_j*cart_k*cart_l;
  int sphsize  = (2*l_i+1)*(2*l_j+1)*(2*l_k+1)*(2*l_l+1);
  double tempVal;


  if (sphsize != shell_element_sph.size() )
    std::cout<<"spherical dimension doesn't match"<<std::endl;
  if (cartsize != shell_element_cart.size() )
    std::cout<<"cartesian dimension doesn't match"<<std::endl;

 
  for( int i = 0 ; i<2*l_i+1 ; i++  ) {
    for( int j = 0 ; j<2*l_j+1 ; j++ ) {
      for( int k = 0 ; k<2*l_k+1 ; k++  ) {
        for( int l = 0 ; l<2*l_l+1 ; l++  ) {
          tempVal = 0.0;
          for( int p = 0 ; p<cart_i ; p++ ) {
            for( int q = 0 ; q<cart_j ; q++ ) {
              for( int r = 0 ; r<cart_k ; r++ ) {
                for( int s = 0 ; s<cart_l ; s++ ) {
                  tempVal += car2sph_matrix[l_i][i*cart_i+p]
                            *car2sph_matrix[l_j][j*cart_j+q]
                            *car2sph_matrix[l_k][k*cart_k+r]
                            *car2sph_matrix[l_l][l*cart_l+s]
                            // *shell_element_cart[p*cart_j+q]
                            *shell_element_cart[p*cart_j*cart_k*cart_l
                                          +q*cart_k*cart_l+r*cart_l+s];
                } // for s
              }  // for r
            }   // for q
          }    // for p

          if (std::abs(tempVal)<1.0e-15) tempVal=0.0;
          shell_element_sph[i*(2*l_j+1)*(2*l_k+1)*(2*l_l+1) + j*(2*l_k+1)*(2*l_l+1)
                            + k*(2*l_l+1) + l ] = tempVal;

        }  // for l
      }   // for k
    }    // for j
  }     // for i
} // cart2sph_2e_transform   


void cart2sph_complex_2e_transform( int l_i, int l_j, int l_k, int l_l, 
  std::vector<dcomplex> &shell_element_sph, std::vector<dcomplex> &shell_element_cart) {

  int cart_i = (l_i+1)*(l_i+2)/2;
  int cart_j = (l_j+1)*(l_j+2)/2;
  int cart_k = (l_k+1)*(l_k+2)/2;
  int cart_l = (l_l+1)*(l_l+2)/2;

  int cartsize = cart_i*cart_j*cart_k*cart_l;
  int sphsize  = (2*l_i+1)*(2*l_j+1)*(2*l_k+1)*(2*l_l+1);
  dcomplex tempVal;


  if (sphsize != shell_element_sph.size() )
    std::cout<<"spherical dimension doesn't match"<<std::endl;
  if (cartsize != shell_element_cart.size() )
    std::cout<<"cartesian dimension doesn't match"<<std::endl;

 
  for( int i = 0 ; i<2*l_i+1 ; i++  ) {
    for( int j = 0 ; j<2*l_j+1 ; j++ ) {
      for( int k = 0 ; k<2*l_k+1 ; k++  ) {
        for( int l = 0 ; l<2*l_l+1 ; l++  ) {
          tempVal = 0.0;
          for( int p = 0 ; p<cart_i ; p++ ) {
            for( int q = 0 ; q<cart_j ; q++ ) {
              for( int r = 0 ; r<cart_k ; r++ ) {
                for( int s = 0 ; s<cart_l ; s++ ) {
                  tempVal +=
                            shell_element_cart[p*cart_j*cart_k*cart_l
                                          +q*cart_k*cart_l+r*cart_l+s]
                            *car2sph_matrix[l_i][i*cart_i+p]
                            *car2sph_matrix[l_j][j*cart_j+q]
                            *car2sph_matrix[l_k][k*cart_k+r]
                            *car2sph_matrix[l_l][l*cart_l+s];
                            // *shell_element_cart[p*cart_j+q]
                } // for s
              }  // for r
            }   // for q
          }    // for p

          if (std::abs(tempVal)<1.0e-15) tempVal=0.0;
          shell_element_sph[i*(2*l_j+1)*(2*l_k+1)*(2*l_l+1) + j*(2*l_k+1)*(2*l_l+1)
                            + k*(2*l_l+1) + l ] = tempVal;

        }  // for l
      }   // for k
    }    // for j
  }     // for i
} // cart2sph_complex_2e_transform   

}; // namespace ChronusQ

