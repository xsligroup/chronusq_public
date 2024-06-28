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
#include <util/typedefs.hpp>

#include <libint2/shell.h>

#include <unordered_map>

namespace ChronusQ {

  struct Molecule;

  enum BASIS_FUNCTION_TYPE {
    REAL_GTO,
    COMPLEX_GIAO,
    COMPLEX_GTO
  };

  struct ShellPairData {

    // List of significant shell pairs for a given shell pair
    using libint2_shellpair_list_t = 
      std::unordered_map<size_t,std::vector<size_t>>;

    // Shell pair data from libint in the same order as 
    // libint2_shellpair_list_t
    using libint2_shellpair_data_t = 
      std::vector<std::vector<std::shared_ptr<libint2::ShellPair>>>;


    libint2_shellpair_list_t sigShellPair; ///< Significant shell pairs
    libint2_shellpair_data_t shData;       ///< Shell pair data

    /**
     *  \brief Populate shell pair list and corresponding pair data
     */
    void computeShellPairs(std::vector<libint2::Shell> &shs, 
        std::vector<size_t> &mapSh2Cen, size_t maxNPrim, size_t maxL, 
        double shell_thresh, double ln_prec);


  };

  /**
   *  \brief The PrimShellHash struct. Compute a Hash value for
   *  a primitive shell.
   *
   *  Warning: does not distinguish contracted shells.
   */
  struct PrimShellHash
  {
    std::size_t operator()(libint2::Shell const& s) const noexcept
    {
      std::vector<size_t> hs{
        std::hash<double>{}(s.alpha[0]),
        std::hash<int>{}(s.contr[0].l),
        std::hash<bool>{}(s.contr[0].pure),
        std::hash<double>{}(s.O[0]),
        std::hash<double>{}(s.O[1]),
        std::hash<double>{}(s.O[2])
      };
      size_t h = 0;
      for (size_t hi : hs)
        h = hi ^ (h << 1);
      return h;
    }
  };

  /**
   * \brief The PrimShellEqual struct. Compare two primitive shells for equality.
   *
   * Warning: does not distinguish contracted shells.
   */
  struct PrimShellEqual {
    bool operator()(const libint2::Shell &s1, const libint2::Shell &s2) const {
      return s1.alpha[0] == s2.alpha[0] and
             s1.contr[0].l == s2.contr[0].l and
             s1.contr[0].pure == s2.contr[0].pure and
             s1.O[0] == s2.O[0] and
             s1.O[1] == s2.O[1] and
             s1.O[2] == s2.O[2];
    }
  };


  /**
   *  \brief The BasisSet struct. Contains information pertinant
   *  for a gaussian type orbital (GTO) basis set. 
   *
   *  Precomputes useful intermediates relating to the basis set
   *  as well as stores useful constants (nBasis, nPrimitive) and 
   *  maps (maps shell # -> basis function #, etc). Utilizes
   *  a storage format consistent with the Libint2  integral
   *  package (libint2::Shell)
   */ 
  struct BasisSet {

    std::string basisName;
    std::string basisDef;

    bool inputDef;
    bool forceCart;

    size_t nBasis = 0;  ///< Number of CGTO basis functions
    size_t nPrimitive;  ///< Number of primitive GTO functions
    size_t nShell;      ///< Number of CGTO basis shells
    size_t maxPrim;     ///< Max primitive dimension of basis
    size_t maxL;        ///< Max angular momentum of basis

    bool nucBasis;      ///< Whether this basis is for nuclear

    BASIS_FUNCTION_TYPE basisType = REAL_GTO; ///< Type of basis function

    cartvec_t centers; ///< A list of centers that comprise the BasisSet

    std::vector<libint2::Shell> shells;    ///< Basis shells
    ShellPairData               shellData; ///< Shell pair data

    std::unordered_map<libint2::Shell, size_t, PrimShellHash, PrimShellEqual> primitives;
      ///< Primitive shells

    std::vector<std::vector<double>> unNormCont;
      ///< Unnormalized basis coefficients
        
    std::vector<size_t> mapSh2Bf;  ///< Map Shell # -> BF #
    std::vector<size_t> mapSh2Cen; ///< Map Shell # -> Cen #
    std::vector<size_t> mapCen2BfSt; ///< Map Cen # -> Starting BF #

    // Disable default constructor
    BasisSet() { };


    // Path / Molecule constructor.
    // See src/basisset/basisset.cxx for documentation
    BasisSet(std::string _basisName, std::string _basisDef,
      bool stringDef, const Molecule &mol, 
      BASIS_FUNCTION_TYPE _basisType, 
      bool _forceCart = false, bool doPrint = true, 
      bool _nucBasis = false);
       
    /**
     *  Copy constructor.
     *
     *  Copies the contents of one BasisSet object to another
     */ 
    BasisSet(const BasisSet &) = default;


    /**
     *  Move constructor.
     *
     *  Moves the contents of one BasisSet object to another
     */ 
    BasisSet(BasisSet &&) = default;


    /**
     *  Copy assignment operator
     *
     *  Assigns one BasisSet object to another through a copy
     */ 
    BasisSet& operator=(const BasisSet &) = default;


    /**
     *  Move assignment operator
     *
     *  Assigns one BasisSet object to another through a move (rvalue reference)
     */ 
    BasisSet& operator=(BasisSet &&) = default;


    // Ouputs BasisSet info, see src/basisset/basisset.cxx for 
    // documentation
    friend std::ostream& operator<<(std::ostream &, const BasisSet&);

    // Misc functions
      
    void uncontractShells();
    //void makeMapPrim2Cont(const double*, double*) const;
    // TangDD Make makeMapPrim2Cont available for Real/Complex 
    template <typename IntsT>
    void makeMapPrim2Cont(const IntsT*, IntsT*) const;
    // TangDD End

    BasisSet uncontractBasis() const;
    BasisSet groupGeneralContractionBasis() const;
    size_t getLibcintEnvLength(const Molecule &mol) const;
    void setLibcintEnv(const Molecule &mol, int *atm, int *bas, double *env,
                       bool finiteWidthNuc = false) const;


    // Update BasisSet object member data (nBasis, etc).
    // See src/basisset/basisset.cxx for documentation
    void update(bool computeShellPairs = true);

    void updateNuclearCoordinates(const Molecule &mol);

  }; // BasisSet struct








  // Basis set keyword map
  static std::unordered_map<std::string,std::string> 
    basisKeyword = {
    {  "STO-3G"         , "sto3g.gbs"                                },
    {  "STO-6G"         , "sto6g.gbs"                                },
    {  "3-21G"          , "3-21g.gbs"                                },
    {  "4-31G"          , "4-31g.gbs"                                },
    {  "6-31G"          , "6-31g.gbs"                                },
    {  "6-31G(D)"       , "6-31g*.gbs"                               },
    {  "6-31++G(D)"     , "6-31++g*.gbs"                             },
    {  "6-311G"         , "6-311g.gbs"                               },
    {  "6-311+G(D)"     , "6-311+g*.gbs"                             },
    {  "6-311+G(D,P)"   , "6-311+g**.gbs"                            },
    {  "6-311+G(2D,P)"  , "6-311+g_2d_p.gbs"                         },
    {  "CC-PVDZ"        , "cc-pvdz.gbs"                              },
    {  "CC-PVTZ"        , "cc-pvtz.gbs"                              },
    {  "CC-PVQZ"        , "cc-pvqz.gbs"                              },
    {  "CC-PV5Z"        , "cc-pv5z.gbs"                              },
    {  "CC-PV6Z"        , "cc-pv6z.gbs"                              },
    {  "CC-PVDZ-DK3"    , "cc-pvdz-dk3.gbs"                          },
    {  "CC-PVTZ-DK3"    , "cc-pvtz-dk3.gbs"                          },
    {  "CC-PVQZ-DK3"    , "cc-pvqz-dk3.gbs"                          },
    {  "CC-PVDZ-X2C"    , "cc-pvdz-x2c.gbs"                          },
    {  "CC-PVTZ-X2C"    , "cc-pvtz-x2c.gbs"                          },
    {  "CC-PVQZ-X2C"    , "cc-pvqz-x2c.gbs"                          },
    {  "DEF2-SVP"                  , "def2-svp.gbs"                  },
    {  "DEF2-SVPD"                 , "def2-svpd.gbs"                 },
    {  "DEF2-TZVP"                 , "def2-tzvp.gbs"                 },
    {  "DEF2-QZVP"                 , "def2-qzvp.gbs"                 },
    {  "SAPPORO-DKH3-DZP-2012-ALL" , "sapporo-dkh3-dzp-2012_all.gbs" },
    {  "SAPPORO-DKH3-DZP-2012-NO"  , "sapporo-dkh3-dzp-2012_no.gbs"  },
    {  "SAPPORO-DKH3-DZP-2012-SP"  , "sapporo-dkh3-dzp-2012_sp.gbs"  },
    {  "SAPPORO-DKH3-DZP-ALL"      , "sapporo-dkh3-dzp_all.gbs"      },
    {  "SAPPORO-DKH3-DZP-NO"       , "sapporo-dkh3-dzp_no.gbs"       },
    {  "SAPPORO-DKH3-DZP-SP"       , "sapporo-dkh3-dzp_sp.gbs"       },
    {  "SAPPORO-DKH3-QZP-2012-ALL" , "sapporo-dkh3-qzp-2012_all.gbs" },
    {  "SAPPORO-DKH3-QZP-2012-NO"  , "sapporo-dkh3-qzp-2012_no.gbs"  },
    {  "SAPPORO-DKH3-QZP-2012-SP"  , "sapporo-dkh3-qzp-2012_sp.gbs"  },
    {  "SAPPORO-DKH3-QZP-ALL"      , "sapporo-dkh3-qzp_all.gbs"      },
    {  "SAPPORO-DKH3-QZP-NO"       , "sapporo-dkh3-qzp_no.gbs"       },
    {  "SAPPORO-DKH3-QZP-SP"       , "sapporo-dkh3-qzp_sp.gbs"       },
    {  "SAPPORO-DKH3-TZP-2012-ALL" , "sapporo-dkh3-tzp-2012_all.gbs" },
    {  "SAPPORO-DKH3-TZP-2012-NO"  , "sapporo-dkh3-tzp-2012_no.gbs"  },
    {  "SAPPORO-DKH3-TZP-2012-SP"  , "sapporo-dkh3-tzp-2012_sp.gbs"  },
    {  "SAPPORO-DKH3-TZP-ALL"      , "sapporo-dkh3-tzp_all.gbs"      },
    {  "SAPPORO-DKH3-TZP-NO"       , "sapporo-dkh3-tzp_no.gbs"       },
    {  "SAPPORO-DKH3-TZP-SP"       , "sapporo-dkh3-tzp_sp.gbs"       },
    {  "SAPPORO-DZP-2012-ALL"      , "sapporo-dzp-2012_all.gbs"      },
    {  "SAPPORO-DZP-2012-NO"       , "sapporo-dzp-2012_no.gbs"       },
    {  "SAPPORO-DZP-2012-SP"       , "sapporo-dzp-2012_sp.gbs"       },
    {  "SAPPORO-DZP-ALL"           , "sapporo-dzp_all.gbs"           },
    {  "SAPPORO-DZP-NO"            , "sapporo-dzp_no.gbs"            },
    {  "SAPPORO-DZP-SP"            , "sapporo-dzp_sp.gbs"            },
    {  "SAPPORO-QZP-2012-ALL"      , "sapporo-qzp-2012_all.gbs"      },
    {  "SAPPORO-QZP-2012-NO"       , "sapporo-qzp-2012_no.gbs"       },
    {  "SAPPORO-QZP-2012-SP"       , "sapporo-qzp-2012_sp.gbs"       },
    {  "SAPPORO-QZP-ALL"           , "sapporo-qzp_all.gbs"           },
    {  "SAPPORO-QZP-NO"            , "sapporo-qzp_no.gbs"            },
    {  "SAPPORO-QZP-SP"            , "sapporo-qzp_sp.gbs"            },
    {  "SAPPORO-TZP-2012-ALL"      , "sapporo-tzp-2012_all.gbs"      },
    {  "SAPPORO-TZP-2012-NO"       , "sapporo-tzp-2012_no.gbs"       },
    {  "SAPPORO-TZP-2012-SP"       , "sapporo-tzp-2012_sp.gbs"       },
    {  "SAPPORO-TZP-ALL"           , "sapporo-tzp_all.gbs"           },
    {  "SAPPORO-TZP-NO"            , "sapporo-tzp_no.gbs"            },
    {  "SAPPORO-TZP-SP"            , "sapporo-tzp_sp.gbs"            },
    {  "CC-PV5Z-JKFIT"             , "cc-pv5z-jkfit.gbs"             },
    {  "CC-PV5Z-RIFIT"             , "cc-pv5z-rifit.gbs"             },
    {  "CC-PV6Z-RIFIT"             , "cc-pv6z-rifit.gbs"             },
    {  "CC-PVDZ-RIFIT"             , "cc-pvdz-rifit.gbs"             },
    {  "CC-PVQZ-JKFIT"             , "cc-pvqz-jkfit.gbs"             },
    {  "CC-PVQZ-RIFIT"             , "cc-pvqz-rifit.gbs"             },
    {  "CC-PVTZ-JKFIT"             , "cc-pvtz-jkfit.gbs"             },
    {  "CC-PVTZ-RIFIT"             , "cc-pvtz-rifit.gbs"             },
    {  "CC-PWCV5Z-RIFIT"           , "cc-pwcv5z-rifit.gbs"           },
    {  "CC-PWCVDZ-RIFIT"           , "cc-pwcvdz-rifit.gbs"           },
    {  "CC-PWCVQZ-RIFIT"           , "cc-pwcvqz-rifit.gbs"           },
    {  "CC-PWCVTZ-RIFIT"           , "cc-pwcvtz-rifit.gbs"           },
    {  "AUG-CC-PV5Z-RIFIT"         , "aug-cc-pv5z-rifit.gbs"         },
    {  "AUG-CC-PV6Z-RIFIT"         , "aug-cc-pv6z-rifit.gbs"         },
    {  "AUG-CC-PVDZ-RIFIT"         , "aug-cc-pvdz-rifit.gbs"         },
    {  "AUG-CC-PVQZ-RIFIT"         , "aug-cc-pvqz-rifit.gbs"         },
    {  "AUG-CC-PVTZ-RIFIT"         , "aug-cc-pvtz-rifit.gbs"         },
    {  "AUG-CC-PWCV5Z-RIFIT"       , "aug-cc-pwcv5z-rifit.gbs"       },
    {  "AUG-CC-PWCVDZ-RIFIT"       , "aug-cc-pwcvdz-rifit.gbs"       },
    {  "AUG-CC-PWCVQZ-RIFIT"       , "aug-cc-pwcvqz-rifit.gbs"       },
    {  "AUG-CC-PWCVTZ-RIFIT"       , "aug-cc-pwcvtz-rifit.gbs"       },
    {  "DEF2-QZVP-RIFIT"           , "def2-qzvp-rifit.gbs"           },
    {  "DEF2-QZVPP-RIFIT"          , "def2-qzvpp-rifit.gbs"          },
    {  "DEF2-QZVPPD-RIFIT"         , "def2-qzvppd-rifit.gbs"         },
    {  "DEF2-SVP-RIFIT"            , "def2-svp-rifit.gbs"            },
    {  "DEF2-SVPD-RIFIT"           , "def2-svpd-rifit.gbs"           },
    {  "DEF2-TZVP-RIFIT"           , "def2-tzvp-rifit.gbs"           },
    {  "DEF2-TZVPD-RIFIT"          , "def2-tzvpd-rifit.gbs"          },
    {  "DEF2-TZVPP-RIFIT"          , "def2-tzvpp-rifit.gbs"          },
    {  "DEF2-TZVPPD-RIFIT"         , "def2-tzvppd-rifit.gbs"         },
    {  "DEF2-UNIVERSAL-JFIT"       , "def2-universal-jfit.gbs"       },
    {  "DEF2-UNIVERSAL-JKFIT"      , "def2-universal-jkfit.gbs"      },
    {  "SARC2-QZV-DKH2-JKFIT"      , "sarc2-qzv-dkh2-jkfit.gbs"      },
    {  "SARC2-QZV-ZORA-JKFIT"      , "sarc2-qzv-zora-jkfit.gbs"      },
    {  "SARC2-QZVP-DKH2-JKFIT"     , "sarc2-qzvp-dkh2-jkfit.gbs"     },
    {  "SARC2-QZVP-ZORA-JKFIT"     , "sarc2-qzvp-zora-jkfit.gbs"     },
    {  "X2C-JFIT"                  , "x2c-jfit.gbs"                  },
    {  "PROT-SP"                   , "prot-sp.gbs"                   },
    {  "PROT-PB4-D"                , "prot-pb4-d.gbs"                },
    {  "PROT-PB4-F1"               , "prot-pb4-f1.gbs"               },
    {  "PROT-PB4-F2"               , "prot-pb4-f2.gbs"               },
    {  "PROT-PB5-D"                , "prot-pb5-d.gbs"                },
    {  "PROT-PB5-F"                , "prot-pb5-f.gbs"                },
    {  "PROT-PB5-G"                , "prot-pb5-g.gbs"                },
    {  "PROT-PB6-D"                , "prot-pb6-d.gbs"                },
    {  "PROT-PB6-F"                , "prot-pb6-f.gbs"                },
    {  "PROT-PB6-G"                , "prot-pb6-g.gbs"                },
    {  "PROT-PB6-H"                , "prot-pb6-h.gbs"                }
  };



















// SS start  
  extern std::vector<std::vector<std::array<int,3>>> cart_ang_list; //list of cartesian angular momentum.
  extern std::vector<std::vector<double>> car2sph_matrix; //SS: external variable,
          //vector of transform matrix from cartesian to spherical
  
  void pop_cart_ang_list(); //SS:populate angular momentum list, can only be called once.   

  void pop_car2sph_matrix(); //populate transform matrix up to L=6

  std::complex<double> cart2sphCoeff(int, int, int, int, int); 
          //SS: calculate the matrix elements of cartesian to spherical transform

  double factorial(int);
  
  double doubleFact(int);

  double polyCoeff(int, int);

  void cart2sph_transform( int, int, std::vector<double>&, std::vector<double>& ); 
  void cart2sph_complex_transform( int, int, std::vector<dcomplex>&, std::vector<dcomplex>& ); 
  
  void cart2sph_2e_transform( int, int, int, int, std::vector<double>&, 
    std::vector<double>& ); 

  void cart2sph_complex_2e_transform( int, int, int, int, std::vector<dcomplex>&, 
    std::vector<dcomplex>& ); 

  // pre calculate shellpair data
  std::vector<libint2::ShellPair> genShellPairs( std::vector<libint2::Shell>&, double );
  std::vector<libint2::ShellPair> genOrderedShellPairs( std::vector<libint2::Shell>&, double );

  // pre calculate shellpair data
  std::vector<libint2::ShellPair> genShellPairs( std::vector<libint2::Shell>&, double );
  std::vector<libint2::ShellPair> genOrderedShellPairs( std::vector<libint2::Shell>&, double );

// SS end
 
}; // namespace ChronusQ

