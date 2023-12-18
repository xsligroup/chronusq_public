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
#include <molecule.hpp>
#include <basisset/basisset_def.hpp>
#include <memmanager.hpp>
#include <libint2/engine.h>

#include <util/files.hpp>

namespace ChronusQ {

    // Indexmap
    int indexmap(int, int, int, int);

    double * cross(double*, double*);
    int * cross(int*, int*);
    double dot(double*, double*);
    int dot(int*, int*);

  struct RealGTOIntEngine {
    // overlap integral of a shell pair  
    static std::vector<std::vector<double>> computeOverlapS(libint2::ShellPair&,
                                         libint2::Shell&,libint2::Shell&);

    // horizontal recursion of contracted overlap integral 
    static double hRRSab(libint2::ShellPair&, libint2::Shell&,libint2::Shell&,
                  int,int*,int,int*);

    
    static double Overlapformula(int*, int*, libint2::ShellPair&,libint2::Shell&,libint2::Shell&);

    static double f_k(int, int, int*, int*, libint2::ShellPair::PrimPairData&,libint2::Shell&,libint2::Shell&);

    static double I_w(int, int*, int*, libint2::ShellPair::PrimPairData&,libint2::Shell&,libint2::Shell&);

    static double unconoverlap(int*, int*, libint2::ShellPair::PrimPairData&,libint2::Shell&,libint2::Shell&);

    // horizontal recursion of uncontracted overlap integral
    static double hRRiPPSab(libint2::ShellPair::PrimPairData&,libint2::Shell&,libint2::Shell&,
                  int,int*,int,int*);

    // vertical recursion of uncontracted overlap integral
    static double vRRSa0(libint2::ShellPair::PrimPairData&,libint2::Shell&,int,int*);

    /* Kinetic Integrals */
    
    // kinetic integrals of a shell pair
    static std::vector<std::vector<double>> computeKineticT(libint2::ShellPair&,
                                         libint2::Shell&,libint2::Shell&);

    // kinetic integral by just calling overlap integral
    static double RRTab(libint2::ShellPair&, libint2::Shell&,libint2::Shell&,
                 int,int*,int,int*);

    // vertical recursion of contracted kinetic integral 
    static double vRRTab(libint2::ShellPair&, libint2::Shell&,libint2::Shell&,
                  int,int*,int,int*);

    // vertical recursion of uncontracted kinetic integral
    static double vRRiPPTab(libint2::ShellPair::PrimPairData&,libint2::Shell&,libint2::Shell&,
                  int,int*,int,int*);

    /* Angular Momentum Integrals */

    // angular momentum integrals of a shell pair
    static std::vector<std::vector<double>> computeAngularL(libint2::ShellPair&,
                                        libint2::Shell&,libint2::Shell&);

    // vertical recursion of uncontracted angular momentum integral
    static double Labmu(libint2::ShellPair::PrimPairData&,libint2::Shell&,libint2::Shell&,
                 double*,double*,int,int*,int,int*,int);

    /* Momentum Integrals */

    // electric dipole (velocity gauge) integrals of a shell pair
    static std::vector<std::vector<double>> computeEDipoleE1_vel(libint2::ShellPair&,
                      libint2::Shell&, libint2::Shell&);

    // contracted momentum integral
    static double Momentummu(libint2::ShellPair&,libint2::Shell&,libint2::Shell&,
                      int,int*,int,int*,int);

    /* Electric Dipole Integrals */

    // electric dipole (length gauge) integrals of a shell pair
    static std::vector<std::vector<double>> computeDipoleE1(libint2::ShellPair&,
                                         libint2::Shell&,libint2::Shell&);

    // contracted electric dipole integrals
    static double DipoleE1(libint2::ShellPair&,libint2::Shell&,libint2::Shell&,
                    int,int*,int,int*,int);

    /* Electric Quadrupole Integrals */

    // electric quadrupole integrals of a shell pair
    static std::vector<std::vector<double>> computeEQuadrupoleE2_vel(libint2::ShellPair&,
                                                  libint2::Shell&,libint2::Shell&); 

    // contracted electric quadrupole integrals of a shell pair
    static double QuadrupoleE2_vel( libint2::ShellPair&,libint2::Shell&,libint2::Shell&,
                             int,int*,int,int*,int,int );

    /* Magnetic Dipole Integrals */

    // contracted magnetic dipole integral
    static double MDipoleM1( libint2::ShellPair&,libint2::Shell&,libint2::Shell&,
                      int,int*,int,int*,int );
  
    /* Magnetic Quadrupole Integrals */

    // contracted magnetic quadrupole integrals 
    static double QuadrupoleM2_vel( libint2::ShellPair&,libint2::Shell&,libint2::Shell&,
                             int,int*,int,int*,int,int );

    // magnetic quadrupole integrals of a shell pair
    static std::vector<std::vector<double>> computeMQuadrupoleM2_vel( libint2::ShellPair&,
                                                  libint2::Shell&,libint2::Shell&);

    /* Electric Octupole Integrals */

    // electric octupole integrals of a shell pair
    static std::vector<std::vector<double>> computeEOctupoleE3_vel( libint2::ShellPair&,
                                                 libint2::Shell&,libint2::Shell&);

    // contracted electric octupole integral
    static double OctupoleE3_vel( libint2::ShellPair&,libint2::Shell&,libint2::Shell&,
                           int,int*,int,int*,int,int,int );

    /* Nuclear Potential Integrals */

    // contracted nuclear potential integrals of a shell pair
    static std::vector<std::vector<double>> computePotentialV(
      const std::vector<libint2::Shell> &, libint2::ShellPair&, 
      libint2::Shell&,libint2::Shell&,const Molecule&); 

    static inline std::vector<std::vector<double>> computePotentialV(
      libint2::ShellPair& pair, libint2::Shell &s1, libint2::Shell &s2, 
      const Molecule& molecule) {
    
      std::vector<libint2::Shell> dummy;
      return computePotentialV(dummy,pair,s1,s2,molecule);

    }

    // horizontal recursion of contracted nuclear potential integrals
    static double hRRVab(const std::vector<libint2::Shell>&,libint2::ShellPair&,
                  libint2::Shell&,libint2::Shell&,int,int*,int,int*,const Molecule&);

    // Bra vertical recursion of uncontracted nuclear potential integrals
    static double vRRVa0(const std::vector<libint2::Shell>&,
                  libint2::ShellPair::PrimPairData&,libint2::Shell&,
                  double*,double*,int,int,int*,int);

    // horizontal recursion of uncontracted nuclear potential integrals
    static double hRRiPPVab(const std::vector<libint2::Shell>&,
      libint2::ShellPair::PrimPairData&, libint2::Shell&,
      libint2::Shell&, int,int*,int,int*,double*,int,int,const Molecule&);
    
    // Ket vertical recursion of uncontracted nuclear potential integrals
    static double vRRV0b(const std::vector<libint2::Shell>&,
                  libint2::ShellPair::PrimPairData&,libint2::Shell&,
                  double*,double*,int,int,int*,int,const Molecule&);

    /* Spin-Orbit Integrals */

    // spin orbit integrals of a shell pair
    static std::vector<std::vector<double>> computeSL(
      const std::vector<libint2::Shell>&, libint2::ShellPair&,
      libint2::Shell&,libint2::Shell&,const Molecule&);

    static inline std::vector<std::vector<double>> computeSL(libint2::ShellPair &pair,
      libint2::Shell &s1, libint2::Shell &s2,const Molecule& mol) {

      std::vector<libint2::Shell> dummy;
      return computeSL(dummy,pair,s1,s2,mol);

    }

    // vertical recursion of uncontracted spin orbit integral
    static double Slabmu(const std::vector<libint2::Shell>&, 
      libint2::ShellPair::PrimPairData&,libint2::Shell&,
      libint2::Shell&, double*,double*,int,int*,int,int*,int,int,int,const Molecule&);

    /* pV dot p Integrals */

    // pV dot p integrals of a shell pair
    static std::vector<std::vector<double>> computepVdotp(
      const std::vector<libint2::Shell>&, libint2::ShellPair&,
      libint2::Shell&,libint2::Shell&,const Molecule&);

    static inline std::vector<std::vector<double>> computepVdotp(
      libint2::ShellPair &pair, libint2::Shell &s1, libint2::Shell &s2,
      const Molecule& mol) {

      std::vector<libint2::Shell> dummy;
      return computepVdotp(dummy,pair,s1,s2,mol);

    }

    // vertical recursion of uncontracted pV dot p integrals
    static double pVpab(const std::vector<libint2::Shell>&,
      libint2::ShellPair::PrimPairData&,libint2::Shell&,
      libint2::Shell&, int,int*,int,int*,int,int,const Molecule&); 



    /* ERI Integrals */

    // Bottumup HGP
    static std::vector<double> BottomupHGP( 
      libint2::ShellPair &, libint2::ShellPair &, 
      libint2::Shell &, libint2::Shell &,
      libint2::Shell &, libint2::Shell &,
      int, int, int, int);

    // HBL 4C:  Bottomup HGP for 2e SOC ERI
    static std::vector<std::vector<double>> BottomupHGP_TwoESP( 
      libint2::ShellPair &, libint2::ShellPair &, 
      libint2::Shell &, libint2::Shell &,
      libint2::Shell &, libint2::Shell &);

    // SS shit : gradient of AC atoms for ERI
    static std::vector<std::vector<double>> ACderiv( 
      libint2::ShellPair &, libint2::ShellPair &, 
      libint2::Shell &, libint2::Shell &,
      libint2::Shell &, libint2::Shell &);

    // SS: bottom up integral code 
    static std::vector<double> bottomupERI(libint2::ShellPair&,
      libint2::ShellPair&,libint2::Shell&,libint2::Shell&,libint2::Shell&,libint2::Shell&,
      int, int, int, int );

    // compute ERI of shell pair 1 and 2
    static std::vector<double> computeERIabcd(libint2::ShellPair&,libint2::ShellPair&,
      libint2::Shell&,libint2::Shell&,libint2::Shell&,libint2::Shell&);

    // horizontal recursion (ab||cd)
    static double twoehRRabcd( libint2::ShellPair&, libint2::ShellPair&, libint2::Shell& , 
      libint2::Shell&, libint2::Shell&, libint2::Shell&, 
      std::vector<std::vector<std::vector<double>>>&, int, int*, int, int*, int, 
      int*, int, int* );

    // horizontal recursion (a0||cd)
    static double twoehRRa0cd( libint2::ShellPair&, libint2::ShellPair&, libint2::Shell& , 
      libint2::Shell&, libint2::Shell&, libint2::Shell&, 
      std::vector<std::vector<std::vector<double>>>&,int, int*, int, int*, int, int* );

    // vertical recursion (a0||c0)
    static double twoevRRa0c0( libint2::ShellPair::PrimPairData&, 
      libint2::ShellPair::PrimPairData&, std::vector<double>&, 
      libint2::Shell& , libint2::Shell&, int, int, int*, int, int* );

    // vertical recursion (a0||00)
    static double twoevRRa000( libint2::ShellPair::PrimPairData&, 
      libint2::ShellPair::PrimPairData&, std::vector<double>&, 
      libint2::Shell& , int, int, int* );

    // (ss||ss) type integral with m=0 
    static double twoeSSSS0( libint2::ShellPair&, libint2::ShellPair&, libint2::Shell&,
      libint2::Shell&, libint2::Shell&, libint2::Shell& ); 



    // compute gauge integral of shell pair 1 and 2
    static std::vector<std::vector<double>> computegaugeabcd(libint2::ShellPair&,libint2::ShellPair&,
      libint2::Shell&,libint2::Shell&,libint2::Shell&,libint2::Shell&, int, int, int, int );

    // horizontal recursion (ab||cd)
    static double twoegaugehRRabcd( libint2::ShellPair&, libint2::ShellPair&, int, int, 
      libint2::Shell& , libint2::Shell&, libint2::Shell&, libint2::Shell&, 
      std::vector<std::vector<std::vector<double>>>&, int, int*, int, int*, int, 
      int*, int, int* );


    // horizontal recursion (a0||cd)
    static double twoegaugehRRa0cd( libint2::ShellPair&, libint2::ShellPair&, int, int, 
      libint2::Shell& ,  libint2::Shell&, libint2::Shell&, libint2::Shell&, 
      std::vector<std::vector<std::vector<double>>>&,int, int*, int, int*, int, int* );

    // vertical recursion (a0||c0)
    static double twoegaugevRRa0c0( libint2::ShellPair::PrimPairData&,  
      libint2::ShellPair::PrimPairData&, int, int, std::vector<double>&, 
      libint2::Shell& , libint2::Shell&, int, int, int*, int, int* );

    // vertical recursion (a0||00)
    static double twoegaugevRRa000( libint2::ShellPair::PrimPairData&,  
      libint2::ShellPair::PrimPairData&, int, int, std::vector<double>&, 
      libint2::Shell& , libint2::Shell& ,int, int, int* );

    // (ss||ss) type integral with m=0 
    static double twoegaugeSSSS0( libint2::ShellPair&, libint2::ShellPair&, int, int, 
      libint2::Shell&, libint2::Shell&, libint2::Shell&, libint2::Shell& ); 

    // SS shit : gradient of AC atoms for gauge integral
    static std::vector<std::vector<std::vector<double>>> ACgaugederiv( 
      libint2::ShellPair &, libint2::ShellPair &, 
      libint2::Shell &, libint2::Shell &,
      libint2::Shell &, libint2::Shell &);

    // SS shit : gradient of BC atoms for gauge integral
    static std::vector<std::vector<std::vector<double>>> BCgaugederiv( 
      libint2::ShellPair &, libint2::ShellPair &, 
      libint2::Shell &, libint2::Shell &,
      libint2::Shell &, libint2::Shell &);

    // Taylor intrapolation of Boys function
    static void computeFmTTaylor(double*,double,int,int);

  };





  struct ComplexGIAOIntEngine {
    /* 1e GIAO Integrals */
 
    // local one body GIAO integral start
    static std::vector<std::vector<dcomplex>> computeGIAOOverlapS(libint2::ShellPair&,
      libint2::Shell&, libint2::Shell&, double* );

    // calculate the uncontracted overlap of (s||s) type for a shellpair
    static std::vector<dcomplex> computecompOverlapss( libint2::ShellPair&, 
      libint2::Shell&, double*, libint2::Shell&, double* ); 

    // complex overlap horizontal recursion for contracted case
    static dcomplex comphRRSab( libint2::ShellPair&, libint2::Shell&, 
      libint2::Shell&, double*, std::vector<dcomplex>&, int, int*, int, int*); 

    // complex overlap horizontal recursion iPP specific for uncontracted case
    static dcomplex comphRRiPPSab( libint2::ShellPair::PrimPairData&, libint2::Shell&, 
      libint2::Shell&, double*, dcomplex, int, int*, int, int* );

    // complex overlap vertical recursion for uncontracted case
    static dcomplex compvRRSa0(libint2::ShellPair::PrimPairData&, libint2::Shell&, double*, 
      dcomplex, int, int* );

    // complex overlap vertical recursion for uncontracted case
    static dcomplex compvRRS0b(libint2::ShellPair::PrimPairData&, libint2::Shell&, double*, 
      dcomplex, int, int* );

    // compute a shell pair of GIAO Kinetic integral
    static std::vector<std::vector<dcomplex>> computeGIAOKineticT( libint2::ShellPair&,
      libint2::Shell&, libint2::Shell&, double* );

    // calculate the uncontracted kinetic integrals of (s||s) type for a shellpair
    static std::vector<dcomplex> computecompKineticss( libint2::ShellPair&,
      libint2::Shell&, double*, libint2::Shell&, double*, std::vector<dcomplex>& );

    // complex kinetic integral for contracted case
    static dcomplex compRRTab( libint2::ShellPair&, libint2::Shell&, 
      libint2::Shell&, double*, double*, std::vector<dcomplex>&, 
      std::vector<dcomplex>&, int, int*, int, int* ); 

    // compute a shell pair of GIAO angular momentum integral 
    static std::vector<std::vector<dcomplex>> computeGIAOAngularL( libint2::ShellPair&,
      libint2::Shell&, libint2::Shell&, double* );

    // complex angular momentum integral for uncontracted case
    static dcomplex compLabmu( libint2::ShellPair::PrimPairData&, 
      libint2::Shell&,libint2::Shell&, double*, double*, dcomplex, int, int*, 
      int, int*, int );

    // compute a shell pair of GIAO electric dipole length gauge integral
    static std::vector<std::vector<dcomplex>> computeGIAOEDipoleE1_len( libint2::ShellPair&,
      libint2::Shell&, libint2::Shell&, double* );

    // complex GIAO electric dipole length gauge
    static dcomplex compDipoleD1_len( libint2::ShellPair&, 
      libint2::Shell&,libint2::Shell&, double*, std::vector<dcomplex>&, 
      int, int*, int, int*, int );

    /* Momentum Integrals */
    // electric dipole (velocity gauge) integrals of a shell pair
    static std::vector<std::vector<dcomplex>> computeGIAOEDipoleE1_vel(libint2::ShellPair&,
                      libint2::Shell&, libint2::Shell&, double* );

    // contracted momentum integral
    static dcomplex GIAOMomentummu(libint2::ShellPair&,libint2::Shell&,libint2::Shell&, 
                      double*, double*, std::vector<dcomplex>&, int,int*,int,int*,int);

    // compute a shell pair of GIAO electric quadrupole momentum length gauge integral 
    static std::vector<std::vector<dcomplex>> computeGIAOEQuadrupoleE2_len( libint2::ShellPair&,
      libint2::Shell&, libint2::Shell&, double* );

    // complex GIAO electric quadrupole momentum length gauge 
    static dcomplex compQuadrupoleE2_len( libint2::ShellPair&, 
      libint2::Shell&,libint2::Shell&, double*, std::vector<dcomplex>&, 
      int, int*, int, int*, int,int );

    // compute a shell pair of GIAO electric octupole length gauge integral
    static std::vector<std::vector<dcomplex>> computeGIAOEOctupoleE3_len( libint2::ShellPair&,
      libint2::Shell&, libint2::Shell&, double* );

    // complex GIAO electric octupole length gauge
    static dcomplex compOctupoleE3_len( libint2::ShellPair&, 
      libint2::Shell&,libint2::Shell&, double*, std::vector<dcomplex>&, 
      int, int*, int, int*, int,int,int );

    // compute a shell pair of GIAO electric prp length gauge integral
    static std::vector<std::vector<dcomplex>> computeGIAOEprp_len( libint2::ShellPair&,
      libint2::Shell&, libint2::Shell&, double* );

    // complex GIAO electric prp length gauge
    static dcomplex compGIAOprp_len( libint2::ShellPair&, 
      libint2::Shell&,libint2::Shell&, double*, double*, std::vector<dcomplex>&, 
      int, int*, int, int*, int, int, int );



    // compute a shell pair of GIAO Potential integral
    static std::vector<std::vector<dcomplex>> computeGIAOPotentialV( const std::vector<libint2::Shell>&,
      libint2::ShellPair&, libint2::Shell&, libint2::Shell&, double*, const Molecule&);

    static inline std::vector<std::vector<dcomplex>> computeGIAOPotentialV( libint2::ShellPair& pair, 
      libint2::Shell &s1, libint2::Shell &s2, double *H, const Molecule& mol) {

      std::vector<libint2::Shell> dummy;
      return computeGIAOPotentialV(dummy,pair,s1,s2,H,mol);
    }   

    // complex potential horizontal recursion for contracted case
    static dcomplex comphRRVab( const std::vector<libint2::Shell>&,
      libint2::ShellPair&, libint2::Shell&, 
      libint2::Shell&, double*, std::vector<dcomplex>&, int, int*, int, int*, const Molecule&); 

    // complex potential vertical recursion for uncontracted case
    static dcomplex compvRRVa0( const std::vector<libint2::Shell>&,
      libint2::ShellPair::PrimPairData&, libint2::Shell&, double*, dcomplex*, double*,
      int, int, int*, int);

    // complex potential horizontal recursion for uncontracted case
    static dcomplex comphRRiPPVab( const std::vector<libint2::Shell>&,
      libint2::ShellPair::PrimPairData&, libint2::Shell&, 
      libint2::Shell&, double*, dcomplex, int, int*, int, int*, int, int, const Molecule&); 

    // complex potential vertical recursion for uncontracted case
    static dcomplex compvRRV0b( const std::vector<libint2::Shell>&,
      libint2::ShellPair::PrimPairData&, libint2::Shell&, double*, dcomplex*, double*,
      int, int, int*, int);

    /* Complex Spin-Orbit Integrals */

    // spin orbit integrals of a shell pair  
    static std::vector<std::vector<dcomplex>> computeGIAOSL( 
      const std::vector<libint2::Shell>&, libint2::ShellPair&, libint2::Shell&, 
      libint2::Shell&, double*, const Molecule& ); 
  
    static inline std::vector<std::vector<dcomplex>> computeGIAOSL( 
      libint2::ShellPair &pair, libint2::Shell &s1, libint2::Shell &s2, double* H, 
      const Molecule& mol ) {
      
      std::vector<libint2::Shell> dummy;   
      return computeGIAOSL(dummy,pair,s1,s2,H,mol);

    }  

    // calculate spin orbit integral for an element
    static dcomplex compSLabmu(const std::vector<libint2::Shell>&, 
      libint2::ShellPair::PrimPairData&, libint2::Shell&, libint2::Shell&, double*, 
      double*, dcomplex, int, int*, int, int*, int, int, const Molecule& );   


    /* GIAO pV dot p Integrals */ 

    // spin free integrals of a shell pair  
    static std::vector<std::vector<dcomplex>> computeGIAOpVdotp( 
      const std::vector<libint2::Shell>&, libint2::ShellPair&, libint2::Shell&, 
      libint2::Shell&, double*, const Molecule& ); 
  
    static inline std::vector<std::vector<dcomplex>> computeGIAOpVdotp( 
      libint2::ShellPair &pair, libint2::Shell &s1, libint2::Shell &s2, double* H, 
      const Molecule& mol ) {
      
      std::vector<libint2::Shell> dummy;   
      return computeGIAOpVdotp(dummy,pair,s1,s2,H,mol);

    }  

    // calculate spin free integral for an element
    static dcomplex comppVpab(const std::vector<libint2::Shell>&, 
      libint2::ShellPair::PrimPairData&, libint2::Shell&, libint2::Shell&, double*, 
      double*, dcomplex, int, int*, int, int*, int, const Molecule& );   

    /* GIAO X2C integrals */
    static std::vector<std::vector<dcomplex>> computeGIAOrVr(
      const std::vector<libint2::Shell>&, libint2::ShellPair&, libint2::Shell&, 
      libint2::Shell&, double*, const Molecule& ); 
      
    static inline std::vector<std::vector<dcomplex>> computeGIAOrVrp( 
      libint2::ShellPair &pair, libint2::Shell &s1, libint2::Shell &s2, double* H, 
      const Molecule& mol ) {
      
      std::vector<libint2::Shell> dummy;   
      return computeGIAOrVr(dummy,pair,s1,s2,H,mol);

    }  

    // complex rVr integral for contracted case
    static dcomplex comprVr( const std::vector<libint2::Shell>&,
      libint2::ShellPair&, libint2::Shell&, 
      libint2::Shell&, double*, std::vector<dcomplex>&, int, int*, int, int*, 
      int, int, const Molecule&); 



    // pVr + rVp integrals of a shell pair  
    static std::vector<std::vector<dcomplex>> computeGIAOpVrprVp( 
      const std::vector<libint2::Shell>&, libint2::ShellPair&, libint2::Shell&, 
      libint2::Shell&, double*, const Molecule& ); 
  
    static inline std::vector<std::vector<dcomplex>> computeGIAOpVrprVp( 
      libint2::ShellPair &pair, libint2::Shell &s1, libint2::Shell &s2, double* H, 
      const Molecule& mol ) {
      
      std::vector<libint2::Shell> dummy;   
      return computeGIAOpVrprVp(dummy,pair,s1,s2,H,mol);

    }  

    // calculate pVr + rVp integral for an element in uncontracted case 
    static dcomplex comppVrprVp(const std::vector<libint2::Shell>&, 
      libint2::ShellPair::PrimPairData&, libint2::Shell&, libint2::Shell&, double*, 
      double*, dcomplex, int, int*, int, int*, int, int, int, const Molecule& );   


    // pVr - rVp integrals of a shell pair  
    static std::vector<std::vector<dcomplex>> computeGIAOpVrmrVp( 
      const std::vector<libint2::Shell>&, libint2::ShellPair&, libint2::Shell&, 
      libint2::Shell&, double*, const Molecule& ); 
  
    static inline std::vector<std::vector<dcomplex>> computeGIAOpVrmrVp( 
      libint2::ShellPair &pair, libint2::Shell &s1, libint2::Shell &s2, double* H, 
      const Molecule& mol ) {
      
      std::vector<libint2::Shell> dummy;   
      return computeGIAOpVrmrVp(dummy,pair,s1,s2,H,mol);

    }  

    // calculate pVr - rVp integral for an element in uncontracted case 
    static dcomplex comppVrmrVp(const std::vector<libint2::Shell>&, 
      libint2::ShellPair::PrimPairData&, libint2::Shell&, libint2::Shell&, double*, 
      double*, dcomplex, int, int*, int, int*, int, int, int, const Molecule& );   






    /* GIAO ERI */
    
    // bottom up GIAO ERI of shell pair 1 and 2 
    static std::vector<dcomplex> bottomupcomplexERI(libint2::ShellPair&,libint2::ShellPair&,
      libint2::Shell&,libint2::Shell&,libint2::Shell&,libint2::Shell&,double* );

    // compute GIAO ERI of shell pair 1 and 2 
    static std::vector<dcomplex> computeGIAOERIabcd(libint2::ShellPair&,libint2::ShellPair&,
      libint2::Shell&,libint2::Shell&,libint2::Shell&,libint2::Shell&,double* );

    // complex horizontal recursion (ab||cd)
    static dcomplex twoecomphRRabcd( libint2::ShellPair&, libint2::ShellPair&, libint2::Shell& ,
      libint2::Shell&, libint2::Shell&, libint2::Shell&, 
      std::vector<std::vector<std::vector<dcomplex>>>&, double*, double*, 
      std::vector<dcomplex>&, std::vector<dcomplex>&, int, int*, int, int*, int, int*, int, int* );

    // complex horizontal recursion (a0||cd)
    static dcomplex twoecomphRRa0cd( libint2::ShellPair&, libint2::ShellPair&, libint2::Shell& ,
      libint2::Shell&, libint2::Shell&, libint2::Shell&, 
      std::vector<std::vector<std::vector<dcomplex>>>&, double*, double*, 
      std::vector<dcomplex>&, std::vector<dcomplex>&, int, int*, int, int*, int, int* );

    // complex vertical recursion (a0||c0)
    static dcomplex twoecompvRRa0c0(libint2::ShellPair::PrimPairData&,libint2::ShellPair::PrimPairData&, 
      std::vector<dcomplex>&, double*, double*, dcomplex, dcomplex, 
      libint2::Shell& , libint2::Shell&, int, int, int*, int, int* );

    // complex vertical recursion (a0||00)
    static dcomplex twoecompvRRa000( libint2::ShellPair::PrimPairData&, libint2::ShellPair::PrimPairData&, 
      std::vector<dcomplex>&, double*, double*, dcomplex, dcomplex, 
      libint2::Shell& , int, int, int* );

    // complex SSSS integral (contracted)
    static dcomplex twoecompSSSS0( libint2::ShellPair&, libint2::ShellPair&, libint2::Shell& ,
      libint2::Shell&, libint2::Shell&, libint2::Shell&, 
      std::vector<std::vector<std::vector<dcomplex>>>&, 
      std::vector<dcomplex>&, std::vector<dcomplex>& );

    // compute complex boys function 
    static void computecompFmT(dcomplex*,dcomplex,int,int);

  };
}

