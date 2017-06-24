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

#include <singleslater/kohnsham.hpp>

#include <grid/integrator.hpp>
#include <basisset/basisset_util.hpp>
#include <cqlinalg/blasutil.hpp>
#include <cqlinalg/blasext.hpp>

#include <util/threads.hpp>

namespace ChronusQ {

  template <typename MatsT, typename IntsT>
  void KohnSham<MatsT,IntsT>::loadFXCder(size_t NPTS, double *Den, double *Gamma, double *EpsEval, double *VRhoEval, 
    double *V2RhoEval, double *VgammaEval, double *V2gammaEval, double *V2RhogammaEval, 
    double *EpsSCR, double *VRhoSCR, double *VgammaSCR, double *V2RhoEvalSCR, double *V2gammaEvalSCR,
    double *V2RhogammaEvalSCR) { 

    for(auto iF = 0; iF < functionals.size(); iF++) {

      double *ES = nullptr, *VR = nullptr, *VS = nullptr, *V2R = nullptr, 
        *V2S = nullptr, *V2RS = nullptr;

      if( iF == 0 ) {
        ES = EpsEval;
        VR = VRhoEval;
        VS = VgammaEval;

        V2R  = V2RhoEval;
        V2S  = V2gammaEval;
        V2RS = V2RhogammaEval;
      } else {
        ES = EpsSCR;
        VR = VRhoSCR;
        VS = VgammaSCR;

        V2R  = V2RhoEvalSCR;
        V2S  = V2gammaEvalSCR;
        V2RS = V2RhogammaEvalSCR;
      }

      if( functionals[iF]->isGGA() )
        functionals[iF]->evalEXC_VXC_FXC(NPTS,Den,Gamma,ES,VR,VS,V2R,V2RS,V2S);
      else
        functionals[iF]->evalEXC_VXC_FXC(NPTS,Den,ES,VR,V2R);

      if( std::any_of(V2R,V2R + 3*NPTS,[&](double x){ return std::isnan(x); }) ) std::cerr << "V2R NANS\n";
      if( functionals[iF]->isGGA() ) {
        if( std::any_of(V2S,V2S + 6*NPTS,[&](double x){ return std::isnan(x); }) ) std::cerr << "V2S NANS\n";
        if( std::any_of(V2RS,V2RS + 6*NPTS,[&](double x){ return std::isnan(x); }) ) std::cerr << "V2RS NANS\n";
      }
      
      if( iF != 0 ) {

        blas::axpy(NPTS  ,1.,ES ,1,EpsEval  ,1);
        blas::axpy(2*NPTS,1.,VR ,1,VRhoEval ,1);
        blas::axpy(3*NPTS,1.,V2R,1,V2RhoEval,1);

        if( functionals[iF]->isGGA() ) {
          blas::axpy(3*NPTS,1.,VS  ,1,VgammaEval    ,1);
          blas::axpy(6*NPTS,1.,V2S ,1,V2gammaEval   ,1);
          blas::axpy(6*NPTS,1.,V2RS,1,V2RhogammaEval,1);
        }

      }

    }

  }; // KohnSham<T>::loadFXCder







  template <typename MatsT, typename IntsT>
  template <typename U>
  void KohnSham<MatsT, IntsT>::constructZVarsFXC(DENSITY_TYPE denTyp, bool isGGA, size_t NPts, 
    double* GDenS, double* GDenZ, double* GDenY, double* GDenX,
    U* TS, U* TZ, U* TY, U* TX,
    U* GTS, U* GTZ, U* GTY, U* GTX,
    double *VR, double *VG, 
    double *V2R, double *V2G, double *V2RG, 
    U *ZrhoVar1, U *ZgammaVar1, U *ZgammaVar2, U *ZgammaVar3, U *ZgammaVar4){

    memset(ZrhoVar1,0,NPts);


    if( denTyp == SCALAR ) {

      for(auto iPt = 0ul; iPt < NPts; iPt++)
        ZrhoVar1[iPt] = TS[iPt] * ( V2R[3*iPt] + 2 * V2R[3*iPt + 1] + V2R[3*iPt + 2] ) +
                        TZ[iPt] * ( V2R[3*iPt] -     V2R[3*iPt + 2] );

      if( isGGA )
      for(auto iPt = 0ul; iPt < NPts; iPt++) {

        U gPTss = (GDenS[iPt]          * GTS[iPt]         ) +
                  (GDenS[iPt + NPts]   * GTS[iPt + NPts]  ) +
                  (GDenS[iPt + 2*NPts] * GTS[iPt + 2*NPts]);

        U gPTsz = (GDenS[iPt]          * GTZ[iPt]         ) +
                  (GDenS[iPt + NPts]   * GTZ[iPt + NPts]  ) +
                  (GDenS[iPt + 2*NPts] * GTZ[iPt + 2*NPts]);

        U gPTzz(0.);
        if( this->onePDM->hasZ() ) {

          // Sum of both <SZ> and <ZS>
          gPTsz += (GDenZ[iPt]          * GTS[iPt]         ) +
                   (GDenZ[iPt + NPts]   * GTS[iPt + NPts]  ) +
                   (GDenZ[iPt + 2*NPts] * GTS[iPt + 2*NPts]);

          gPTzz = (GDenZ[iPt]          * GTZ[iPt]         ) +
                  (GDenZ[iPt + NPts]   * GTZ[iPt + NPts]  ) +
                  (GDenZ[iPt + 2*NPts] * GTZ[iPt + 2*NPts]);

        }


        ZrhoVar1[iPt] += ( V2RG[6*iPt    ] + V2RG[6*iPt + 1] + V2RG[6*iPt + 2] + 
                           V2RG[6*iPt + 3] + V2RG[6*iPt + 4] + V2RG[6*iPt + 5]   ) * gPTss +
                         ( V2RG[6*iPt    ] - V2RG[6*iPt + 2] + V2RG[6*iPt + 3] -
                           V2RG[6*iPt + 5]                                       ) * gPTsz +
                         ( V2RG[6*iPt    ] - V2RG[6*iPt + 1] + V2RG[6*iPt + 2] + 
                           V2RG[6*iPt + 3] - V2RG[6*iPt + 4] + V2RG[6*iPt + 5]   ) * gPTzz;



        ZgammaVar1[iPt] = 2 * ( VG[3*iPt] + VG[3*iPt + 1] + VG[3*iPt + 2] );
        ZgammaVar2[iPt] = 2 * ( VG[3*iPt] - VG[3*iPt + 2]                 );

        ZgammaVar3[iPt] = 
          ( V2G[6*iPt    ] + 2*V2G[6*iPt + 1] + 2*V2G[6*iPt + 2] + 
            V2G[6*iPt + 3] + 2*V2G[6*iPt + 4] +   V2G[6*iPt + 5]   ) * gPTss +
          ( V2G[6*iPt    ] +   V2G[6*iPt + 1] -   V2G[6*iPt + 4] -
            V2G[6*iPt + 5]                                         ) * gPTsz +
          ( V2G[6*iPt    ] + 2*V2G[6*iPt + 2] -   V2G[6*iPt + 3] +
            V2G[6*iPt + 5]                                         ) * gPTzz +

         ( V2RG[6*iPt    ] + V2RG[6*iPt + 1] + V2RG[6*iPt + 2] + 
           V2RG[6*iPt + 3] + V2RG[6*iPt + 4] + V2RG[6*iPt + 5]   ) * TS[iPt] +
         ( V2RG[6*iPt    ] + V2RG[6*iPt + 1] + V2RG[6*iPt + 2] - 
           V2RG[6*iPt + 3] - V2RG[6*iPt + 4] - V2RG[6*iPt + 5]   ) * TZ[iPt];

        if( this->onePDM->hasZ() )
          ZgammaVar4[iPt] = 
            ( V2G[6*iPt    ] + 2*V2G[6*iPt + 1] -   V2G[6*iPt + 4] -
              V2G[6*iPt + 5]                                         ) * gPTss +
            ( V2G[6*iPt    ] - 2*V2G[6*iPt + 2] +   V2G[6*iPt + 5]   ) * gPTsz +
            ( V2G[6*iPt    ] -   V2G[6*iPt + 1] +   V2G[6*iPt + 4] -
              V2G[6*iPt + 5]                                         ) * gPTzz +
  
           ( V2RG[6*iPt    ] - V2RG[6*iPt + 2] + V2RG[6*iPt + 3] - V2RG[6*iPt + 5]   ) * TS[iPt] +
           ( V2RG[6*iPt    ] - V2RG[6*iPt + 2] - V2RG[6*iPt + 3] + V2RG[6*iPt + 5]   ) * TZ[iPt];

      }

    } else {

      for(auto iPt = 0ul; iPt < NPts; iPt++)
        ZrhoVar1[iPt] = TZ[iPt] * ( V2R[3*iPt] - 2 * V2R[3*iPt + 1] + V2R[3*iPt + 2] ) +
                        TS[iPt] * ( V2R[3*iPt] -     V2R[3*iPt + 2] );



      if( isGGA )
      for(auto iPt = 0ul; iPt < NPts; iPt++) {

        U gPTss = (GDenS[iPt]          * GTS[iPt]         ) +
                  (GDenS[iPt + NPts]   * GTS[iPt + NPts]  ) +
                  (GDenS[iPt + 2*NPts] * GTS[iPt + 2*NPts]);

        U gPTsz = (GDenS[iPt]          * GTZ[iPt]         ) +
                  (GDenS[iPt + NPts]   * GTZ[iPt + NPts]  ) +
                  (GDenS[iPt + 2*NPts] * GTZ[iPt + 2*NPts]);

        U gPTzz(0.);
        if( this->onePDM->hasZ() ) {

          // Sum of both <SZ> and <ZS>
          gPTsz += (GDenZ[iPt]          * GTS[iPt]         ) +
                   (GDenZ[iPt + NPts]   * GTS[iPt + NPts]  ) +
                   (GDenZ[iPt + 2*NPts] * GTS[iPt + 2*NPts]);

          gPTzz = (GDenZ[iPt]          * GTZ[iPt]         ) +
                  (GDenZ[iPt + NPts]   * GTZ[iPt + NPts]  ) +
                  (GDenZ[iPt + 2*NPts] * GTZ[iPt + 2*NPts]);

        }



        ZrhoVar1[iPt] += ( V2RG[6*iPt    ] + V2RG[6*iPt + 1] + V2RG[6*iPt + 2] - 
                           V2RG[6*iPt + 3] - V2RG[6*iPt + 4] - V2RG[6*iPt + 5]   ) * gPTss +
                         ( V2RG[6*iPt    ] - V2RG[6*iPt + 2] - V2RG[6*iPt + 3] +
                           V2RG[6*iPt + 5]                                       ) * gPTsz + 
                         ( V2RG[6*iPt    ] - V2RG[6*iPt + 1] + V2RG[6*iPt + 2] - 
                           V2RG[6*iPt + 3] + V2RG[6*iPt + 4] - V2RG[6*iPt + 5]   ) * gPTzz;

        ZgammaVar1[iPt] = 2 * ( VG[3*iPt] - VG[3*iPt + 2]                 );
        ZgammaVar2[iPt] = 2 * ( VG[3*iPt] - VG[3*iPt + 1] + VG[3*iPt + 2] );

        ZgammaVar3[iPt] = 
          ( V2G[6*iPt    ] +   V2G[6*iPt + 1] - V2G[6*iPt + 4] - V2G[6*iPt + 5]   ) * gPTss +
          ( V2G[6*iPt    ] - 2*V2G[6*iPt + 2] + V2G[6*iPt + 5]                    ) * gPTsz +
          ( V2G[6*iPt    ] -   V2G[6*iPt + 1] + V2G[6*iPt + 4] - V2G[6*iPt + 5]   ) * gPTzz +

          ( V2RG[6*iPt    ] - V2RG[6*iPt + 2] + V2RG[6*iPt + 3] - V2RG[6*iPt + 5] ) * TS[iPt] +
          ( V2RG[6*iPt    ] - V2RG[6*iPt + 2] - V2RG[6*iPt + 3] + V2RG[6*iPt + 5] ) * TZ[iPt];

        if( this->onePDM->hasZ() )
          ZgammaVar4[iPt] = 
            ( V2G[6*iPt    ] + 2*V2G[6*iPt + 2] -   V2G[6*iPt + 3] +
              V2G[6*iPt + 5]                                         ) * gPTss +
            ( V2G[6*iPt    ] -   V2G[6*iPt + 1] +   V2G[6*iPt + 4] -
              V2G[6*iPt + 5]                                         ) * gPTsz +
            ( V2G[6*iPt    ] - 2*V2G[6*iPt + 1] + 2*V2G[6*iPt + 2] + 
              V2G[6*iPt + 3] - 2*V2G[6*iPt + 4] +   V2G[6*iPt + 5]   ) * gPTzz +
  
           ( V2RG[6*iPt    ] - V2RG[6*iPt + 1] + V2RG[6*iPt + 2] + 
             V2RG[6*iPt + 3] - V2RG[6*iPt + 4] + V2RG[6*iPt + 5]   ) * TS[iPt] +
           ( V2RG[6*iPt    ] - V2RG[6*iPt + 1] + V2RG[6*iPt + 2] - 
             V2RG[6*iPt + 3] + V2RG[6*iPt + 4] - V2RG[6*iPt + 5]   ) * TZ[iPt];

      }

    }

  }

  template <typename MatsT, typename IntsT>
  template <typename U>
  void KohnSham<MatsT, IntsT>::constructZVarsFXC(DENSITY_TYPE denTyp, bool isGGA, size_t NPts, 
    double* GDenS, double* GDenZ, double* GDenY, double* GDenX,
    bool * Msmall, double *Mnorm, 
    double *Kx, double *Ky, double *Kz, 
    double *Hx, double *Hy, double *Hz,
    U* TS, U* TZ, U* TY, U* TX,
    U* GTS, U* GTZ, U* GTY, U* GTX,
    U* gPTssv, U* gPTszv, U* gPTsyv, U* gPTsxv, U* gPTzzv, 
    U* gPTyyv, U* gPTxxv,  
    double *VR, double *VG, 
    double *V2R, double *V2G, double *V2RG, 
    U *ZrhoVar1, U *ZgammaVar1, U *ZgammaVar2, U *ZgammaVar3, U *ZgammaVar4){

    memset(ZrhoVar1,0,NPts);


    if( not this->onePDM->hasXY() ) {
      if( denTyp == SCALAR ) {

        for(auto iPt = 0ul; iPt < NPts; iPt++){
            
          ZrhoVar1[iPt] =  TS[iPt] * ( V2R[3*iPt] + 2 * V2R[3*iPt + 1] + V2R[3*iPt + 2] ) +
                           TZ[iPt] * ( V2R[3*iPt] -     V2R[3*iPt + 2] );
                           
        }
        if( isGGA )
        for(auto iPt = 0ul; iPt < NPts; iPt++) {
  
          U gPTss = (GDenS[iPt]          * GTS[iPt]         ) +
                    (GDenS[iPt + NPts]   * GTS[iPt + NPts]  ) +
                    (GDenS[iPt + 2*NPts] * GTS[iPt + 2*NPts]);
  
          U gPTsz = (GDenS[iPt]          * GTZ[iPt]         ) +
                    (GDenS[iPt + NPts]   * GTZ[iPt + NPts]  ) +
                    (GDenS[iPt + 2*NPts] * GTZ[iPt + 2*NPts]);
  
          U gPTzz(0.);
          if( this->onePDM->hasZ() ) {
  
            // Sum of both <SZ> and <ZS>
            gPTsz += (GDenZ[iPt]          * GTS[iPt]         ) +
                     (GDenZ[iPt + NPts]   * GTS[iPt + NPts]  ) +
                     (GDenZ[iPt + 2*NPts] * GTS[iPt + 2*NPts]);
  
            gPTzz = (GDenZ[iPt]          * GTZ[iPt]         ) +
                    (GDenZ[iPt + NPts]   * GTZ[iPt + NPts]  ) +
                    (GDenZ[iPt + 2*NPts] * GTZ[iPt + 2*NPts]);
  
          }
  
  
          ZrhoVar1[iPt] += ( V2RG[6*iPt    ] + V2RG[6*iPt + 1] + V2RG[6*iPt + 2] + 
                             V2RG[6*iPt + 3] + V2RG[6*iPt + 4] + V2RG[6*iPt + 5]   ) * gPTss +
                           ( V2RG[6*iPt    ] - V2RG[6*iPt + 2] + V2RG[6*iPt + 3] -
                             V2RG[6*iPt + 5]                                       ) * gPTsz +
                           ( V2RG[6*iPt    ] - V2RG[6*iPt + 1] + V2RG[6*iPt + 2] + 
                             V2RG[6*iPt + 3] - V2RG[6*iPt + 4] + V2RG[6*iPt + 5]   ) * gPTzz;
  
  
  
          ZgammaVar1[iPt] = 2 * ( VG[3*iPt] + VG[3*iPt + 1] + VG[3*iPt + 2] );
          ZgammaVar2[iPt] = 2 * ( VG[3*iPt] - VG[3*iPt + 2]                 );
  
          ZgammaVar3[iPt] = 
            ( V2G[6*iPt    ] + 2*V2G[6*iPt + 1] + 2*V2G[6*iPt + 2] + 
              V2G[6*iPt + 3] + 2*V2G[6*iPt + 4] +   V2G[6*iPt + 5]   ) * gPTss +
            ( V2G[6*iPt    ] +   V2G[6*iPt + 1] -   V2G[6*iPt + 4] -
              V2G[6*iPt + 5]                                         ) * gPTsz +
            ( V2G[6*iPt    ] + 2*V2G[6*iPt + 2] -   V2G[6*iPt + 3] +
              V2G[6*iPt + 5]                                         ) * gPTzz +
  
           ( V2RG[6*iPt    ] + V2RG[6*iPt + 1] + V2RG[6*iPt + 2] + 
             V2RG[6*iPt + 3] + V2RG[6*iPt + 4] + V2RG[6*iPt + 5]   ) * TS[iPt] +
           ( V2RG[6*iPt    ] + V2RG[6*iPt + 1] + V2RG[6*iPt + 2] - 
             V2RG[6*iPt + 3] - V2RG[6*iPt + 4] - V2RG[6*iPt + 5]   ) * TZ[iPt];
  
          if( this->onePDM->hasZ() )
            ZgammaVar4[iPt] = 
              ( V2G[6*iPt    ] + V2G[6*iPt + 1] -   V2G[6*iPt + 4] -
                V2G[6*iPt + 5]                                         ) * gPTss +
              ( V2G[6*iPt    ] - 2*V2G[6*iPt + 2] +   V2G[6*iPt + 5]   ) * gPTsz +
              ( V2G[6*iPt    ] -   V2G[6*iPt + 1] +   V2G[6*iPt + 4] -
                V2G[6*iPt + 5]                                         ) * gPTzz +
    
             ( V2RG[6*iPt    ] - V2RG[6*iPt + 2] + V2RG[6*iPt + 3] - V2RG[6*iPt + 5]   ) * TS[iPt] +
             ( V2RG[6*iPt    ] - V2RG[6*iPt + 2] - V2RG[6*iPt + 3] + V2RG[6*iPt + 5]   ) * TZ[iPt];
  
        } //iPts
  
      } else { // MZ 1C
  
        for(auto iPt = 0ul; iPt < NPts; iPt++)
          ZrhoVar1[iPt] = TZ[iPt] * ( V2R[3*iPt] - 2 * V2R[3*iPt + 1] + V2R[3*iPt + 2] ) +
                          TS[iPt] * ( V2R[3*iPt] -     V2R[3*iPt + 2] );
  
  
        if( isGGA )
        for(auto iPt = 0ul; iPt < NPts; iPt++) {
  
          U gPTss = (GDenS[iPt]          * GTS[iPt]         ) +
                    (GDenS[iPt + NPts]   * GTS[iPt + NPts]  ) +
                    (GDenS[iPt + 2*NPts] * GTS[iPt + 2*NPts]);
  
          U gPTsz = (GDenS[iPt]          * GTZ[iPt]         ) +
                    (GDenS[iPt + NPts]   * GTZ[iPt + NPts]  ) +
                    (GDenS[iPt + 2*NPts] * GTZ[iPt + 2*NPts]);
  
          U gPTzz(0.);
          if( this->onePDM->hasZ() ) {
  
            // Sum of both <SZ> and <ZS>
            gPTsz += (GDenZ[iPt]          * GTS[iPt]         ) +
                     (GDenZ[iPt + NPts]   * GTS[iPt + NPts]  ) +
                     (GDenZ[iPt + 2*NPts] * GTS[iPt + 2*NPts]);
  
            gPTzz = (GDenZ[iPt]          * GTZ[iPt]         ) +
                    (GDenZ[iPt + NPts]   * GTZ[iPt + NPts]  ) +
                    (GDenZ[iPt + 2*NPts] * GTZ[iPt + 2*NPts]);
  
          }

  
          ZrhoVar1[iPt] += ( V2RG[6*iPt    ] + V2RG[6*iPt + 1] + V2RG[6*iPt + 2] - 
                             V2RG[6*iPt + 3] - V2RG[6*iPt + 4] - V2RG[6*iPt + 5]   ) * gPTss +
                           ( V2RG[6*iPt    ] - V2RG[6*iPt + 2] - V2RG[6*iPt + 3] +
                             V2RG[6*iPt + 5]                                       ) * gPTsz + 
                           ( V2RG[6*iPt    ] - V2RG[6*iPt + 1] + V2RG[6*iPt + 2] - 
                             V2RG[6*iPt + 3] + V2RG[6*iPt + 4] - V2RG[6*iPt + 5]   ) * gPTzz;




          ZgammaVar1[iPt] = 2 * ( VG[3*iPt] - VG[3*iPt + 2]                 );
          ZgammaVar2[iPt] = 2 * ( VG[3*iPt] - VG[3*iPt + 1] + VG[3*iPt + 2] );
  
          ZgammaVar3[iPt] = 
            ( V2G[6*iPt    ] +   V2G[6*iPt + 1] - V2G[6*iPt + 4] - V2G[6*iPt + 5]   ) * gPTss +
            ( V2G[6*iPt    ] - 2*V2G[6*iPt + 2] + V2G[6*iPt + 5]                    ) * gPTsz +
            ( V2G[6*iPt    ] -   V2G[6*iPt + 1] + V2G[6*iPt + 4] - V2G[6*iPt + 5]   ) * gPTzz +
  
            ( V2RG[6*iPt    ] - V2RG[6*iPt + 2] + V2RG[6*iPt + 3] - V2RG[6*iPt + 5] ) * TS[iPt] +
            ( V2RG[6*iPt    ] - V2RG[6*iPt + 2] - V2RG[6*iPt + 3] + V2RG[6*iPt + 5] ) * TZ[iPt];
  
          if( this->onePDM->hasZ() )
            ZgammaVar4[iPt] = 
              ( V2G[6*iPt    ] + 2*V2G[6*iPt + 2] -   V2G[6*iPt + 3] +
                V2G[6*iPt + 5]                                         ) * gPTss +
              ( V2G[6*iPt    ] -   V2G[6*iPt + 1] +   V2G[6*iPt + 4] -
                V2G[6*iPt + 5]                                         ) * gPTsz +
              ( V2G[6*iPt    ] - 2*V2G[6*iPt + 1] + 2*V2G[6*iPt + 2] + 
                V2G[6*iPt + 3] - 2*V2G[6*iPt + 4] +   V2G[6*iPt + 5]   ) * gPTzz +
    
             ( V2RG[6*iPt    ] - V2RG[6*iPt + 1] + V2RG[6*iPt + 2] + 
               V2RG[6*iPt + 3] - V2RG[6*iPt + 4] + V2RG[6*iPt + 5]   ) * TS[iPt] +
             ( V2RG[6*iPt    ] - V2RG[6*iPt + 1] + V2RG[6*iPt + 2] - 
               V2RG[6*iPt + 3] + V2RG[6*iPt + 4] - V2RG[6*iPt + 5]   ) * TZ[iPt];


        } // iPt
  
      }

    } else { // Scalar 2c

      //  SCALAR 
      if( denTyp == SCALAR ) {
        for(auto iPt = 0ul; iPt < NPts; iPt++){

            ZrhoVar1[iPt] =  TS[iPt] * ( V2R[3*iPt] + 2 * V2R[3*iPt + 1] + V2R[3*iPt + 2] ) +
                              (   Kz[iPt] * TZ[iPt]   
                                + Ky[iPt] * TY[iPt]    
                                + Kx[iPt] * TX[iPt]      ) * ( V2R[3*iPt] -     V2R[3*iPt + 2] );

        }
        if( isGGA ) {
          for(auto iPt = 0ul; iPt < NPts; iPt++) {
  
            ZrhoVar1[iPt] += ( V2RG[6*iPt    ] + V2RG[6*iPt + 1] + V2RG[6*iPt + 2] + 
                               V2RG[6*iPt + 3] + V2RG[6*iPt + 4] + V2RG[6*iPt + 5]   ) 
                              * gPTssv[iPt] 
                     
                           + ( V2RG[6*iPt    ] - V2RG[6*iPt + 2] + V2RG[6*iPt + 3] -
                               V2RG[6*iPt + 5]                                       ) 
                              * (gPTsxv[iPt] * Hx[iPt] + gPTsyv[iPt] * Hy[iPt] + gPTszv[iPt] * Hz[iPt] )

                           + ( V2RG[6*iPt    ] - V2RG[6*iPt + 1] + V2RG[6*iPt + 2] + 
                               V2RG[6*iPt + 3] - V2RG[6*iPt + 4] + V2RG[6*iPt + 5]   )  
                               * (gPTxxv[iPt] + gPTyyv[iPt] + gPTzzv[iPt]);

            ZgammaVar1[iPt] = 2 * ( VG[3*iPt] + VG[3*iPt + 1] + VG[3*iPt + 2] );

            ZgammaVar2[iPt] = 2 * ( VG[3*iPt] - VG[3*iPt + 2]                 );

            ZgammaVar3[iPt] = 
              ( V2G[6*iPt    ] + 2*V2G[6*iPt + 1] + 2*V2G[6*iPt + 2] + 
                V2G[6*iPt + 3] + 2*V2G[6*iPt + 4] +   V2G[6*iPt + 5]   ) 
              * gPTssv[iPt] +

              ( V2G[6*iPt    ] + 2*V2G[6*iPt + 2] -   V2G[6*iPt + 3] +
                V2G[6*iPt + 5]                                         ) 
              * (gPTxxv[iPt] + gPTyyv[iPt] + gPTzzv[iPt]) +
    
             ( V2RG[6*iPt    ] + V2RG[6*iPt + 1] + V2RG[6*iPt + 2] + 
               V2RG[6*iPt + 3] + V2RG[6*iPt + 4] + V2RG[6*iPt + 5]   ) 
              * TS[iPt] 
;

            ZgammaVar4[iPt] = 
              ( V2G[6*iPt    ] +  V2G[6*iPt + 1] -   V2G[6*iPt + 4] -
                V2G[6*iPt + 5]                                         ) 
              * gPTssv[iPt] +

              ( V2G[6*iPt    ] -   V2G[6*iPt + 1] +   V2G[6*iPt + 4] -
                V2G[6*iPt + 5]                                         )
              * (gPTxxv[iPt] + gPTyyv[iPt] + gPTzzv[iPt]) +
            
             ( V2RG[6*iPt    ] - V2RG[6*iPt + 2] + V2RG[6*iPt + 3] - V2RG[6*iPt + 5]   ) 
              * TS[iPt] 
;

               
              ZgammaVar3[iPt]+= ( V2G[6*iPt    ] +   V2G[6*iPt + 1] -   V2G[6*iPt + 4] -
                  V2G[6*iPt + 5]                                         ) 
                * (gPTsxv[iPt] * Hx[iPt] + gPTsyv[iPt] * Hy[iPt] + gPTszv[iPt] * Hz[iPt] ) +

             ( V2RG[6*iPt    ] + V2RG[6*iPt + 1] + V2RG[6*iPt + 2] - 
               V2RG[6*iPt + 3] - V2RG[6*iPt + 4] - V2RG[6*iPt + 5]   ) 
              * (   Kz[iPt] * TZ[iPt]   
                  + Ky[iPt] * TY[iPt]    
                  + Kx[iPt] * TX[iPt]  );
 
             ZgammaVar4[iPt] += ( V2G[6*iPt    ] - 2*V2G[6*iPt + 2] +   V2G[6*iPt + 5]   ) 
               * (gPTsxv[iPt] * Hx[iPt] + gPTsyv[iPt] * Hy[iPt] + gPTszv[iPt] * Hz[iPt] ) + 
               ( V2RG[6*iPt    ] - V2RG[6*iPt + 2] - V2RG[6*iPt + 3] + V2RG[6*iPt + 5]   )
               * (   Kz[iPt] * TZ[iPt] 
                   + Ky[iPt] * TY[iPt] 
                   + Kx[iPt] * TX[iPt]  );  

          } //NPTs
        } // GGA

       //  MX 2c 
       } else if (denTyp == MX ) {
         for(auto iPt = 0ul; iPt < NPts; iPt++) {
           if (!Msmall[iPt]) {
             ZrhoVar1[iPt] = ( Kx[iPt] * TX[iPt] + Ky[iPt] * TY[iPt] + Kz[iPt] * TZ[iPt]) * 
                             ( V2R[3*iPt] - 2 * V2R[3*iPt + 1] + V2R[3*iPt + 2] ) +
                               TS[iPt] * ( V2R[3*iPt] -     V2R[3*iPt + 2] );
             ZrhoVar1[iPt] *=  Kx[iPt] ; 
             ZrhoVar1[iPt] +=  2.*(VR[2*iPt] - VR[2*iPt + 1]) 
                                * (TX[iPt]

                                 - (  Kx[iPt] *Kz[iPt] * TZ[iPt]   
                                    + Kx[iPt] *Ky[iPt] * TY[iPt]   
                                    + Kx[iPt] *Kx[iPt] * TX[iPt]   ) )/(Mnorm[iPt]);

             } else {
               ZrhoVar1[iPt] = ( TX[iPt] ) * 
                               ( V2R[3*iPt] - 2 * V2R[3*iPt + 1] + V2R[3*iPt + 2] ) ;
                               
               ZrhoVar1[iPt]+= TS[iPt] * ( V2R[3*iPt] -     V2R[3*iPt + 2] ); 
      
             }
           }
       } else if (denTyp == MY ) {
        for(auto iPt = 0ul; iPt < NPts; iPt++) {
          if (!Msmall[iPt]) {
            ZrhoVar1[iPt] = ( Kx[iPt] * TX[iPt] + Ky[iPt] * TY[iPt] + Kz[iPt] * TZ[iPt]) * 
                            ( V2R[3*iPt] - 2 * V2R[3*iPt + 1] + V2R[3*iPt + 2] ) +
                              TS[iPt] * ( V2R[3*iPt] -     V2R[3*iPt + 2] );
            ZrhoVar1[iPt] *=  Ky[iPt] ; 
            ZrhoVar1[iPt] +=  2.*(VR[2*iPt] - VR[2*iPt + 1]) 
                               * (TY[iPt]

                                - ( Ky[iPt] * Kz[iPt] * TZ[iPt]   
                                   +Ky[iPt] * Ky[iPt] * TY[iPt]   
                                   +Ky[iPt] * Kx[iPt] * TX[iPt]   ) )/(Mnorm[iPt]);

          } else {
            ZrhoVar1[iPt] = ( TY[iPt] ) * 
                            ( V2R[3*iPt] - 2 * V2R[3*iPt + 1] + V2R[3*iPt + 2] ) ;

            ZrhoVar1[iPt]+= TS[iPt] * ( V2R[3*iPt] -     V2R[3*iPt + 2] ); 

            }
          }
       } else if (denTyp == MZ ) {
        for(auto iPt = 0ul; iPt < NPts; iPt++) {
          if (!Msmall[iPt]) {
             ZrhoVar1[iPt] = ( Kx[iPt] * TX[iPt] + Ky[iPt] * TY[iPt] + Kz[iPt] * TZ[iPt]) * 
                             ( V2R[3*iPt] - 2 * V2R[3*iPt + 1] + V2R[3*iPt + 2] ) +
                               TS[iPt] * ( V2R[3*iPt] -     V2R[3*iPt + 2] );
             ZrhoVar1[iPt]  *=  Kz[iPt] ; 
             ZrhoVar1[iPt] +=  2.*(VR[2*iPt] - VR[2*iPt + 1]) 
                                * (TZ[iPt]
                                
                                 - ( Kz[iPt] * Kz[iPt] * TZ[iPt]   
                                    +Kz[iPt] * Ky[iPt] * TY[iPt]   
                                    +Kz[iPt] * Kx[iPt] * TX[iPt]   ) )/(Mnorm[iPt]);

            } else {
             ZrhoVar1[iPt] = ( TZ[iPt] ) * 
                             ( V2R[3*iPt] - 2 * V2R[3*iPt + 1] + V2R[3*iPt + 2] ) ;

            ZrhoVar1[iPt]+= TS[iPt] * ( V2R[3*iPt] -     V2R[3*iPt + 2] ); 

            }
          }
       }
      
       // GGA 2C for non SCALAR
       if( isGGA && denTyp != SCALAR ) {
         for(auto iPt = 0ul; iPt < NPts; iPt++) {

       //ZrhoVar1 for non scalar
           U tmp   = ( V2RG[6*iPt    ] + V2RG[6*iPt + 1] + V2RG[6*iPt + 2] - 
                       V2RG[6*iPt + 3] - V2RG[6*iPt + 4] - V2RG[6*iPt + 5]   ) 
                       * gPTssv[iPt]; // +


           if (!Msmall[iPt]) {
             tmp+=        ( V2RG[6*iPt    ] - V2RG[6*iPt + 2] - V2RG[6*iPt + 3] +
                       V2RG[6*iPt + 5]                                       ) 

                       * (gPTsxv[iPt] * Hx[iPt] + gPTsyv[iPt] * Hy[iPt] + gPTszv[iPt] * Hz[iPt] );// +
           } else {
             if (denTyp == MX ) { 
               tmp+= ( V2RG[6*iPt    ] - V2RG[6*iPt + 2] - V2RG[6*iPt + 3] + V2RG[6*iPt + 5] )
*  gPTsxv[iPt]  ;
             } else if (denTyp == MY) { 
               tmp+= ( V2RG[6*iPt    ] - V2RG[6*iPt + 2] - V2RG[6*iPt + 3] + V2RG[6*iPt + 5] )
*  gPTsyv[iPt]  ;
             } else if (denTyp == MZ) { 
               tmp+= ( V2RG[6*iPt    ] - V2RG[6*iPt + 2] - V2RG[6*iPt + 3] + V2RG[6*iPt + 5] )
*  gPTszv[iPt]  ;
             }
           }  


           tmp +=
                     ( V2RG[6*iPt    ] - V2RG[6*iPt + 1] + V2RG[6*iPt + 2] - 
                       V2RG[6*iPt + 3] + V2RG[6*iPt + 4] - V2RG[6*iPt + 5]   ) 
                      * (gPTxxv[iPt] + gPTyyv[iPt] + gPTzzv[iPt]);



           if (!Msmall[iPt]) {
             if (denTyp == MX ) {
               ZrhoVar1[iPt] += Kx[iPt]*tmp;
             } else if (denTyp == MY) {
               ZrhoVar1[iPt] += Ky[iPt]*tmp;
             } else if (denTyp == MZ) {
               ZrhoVar1[iPt] += Kz[iPt]*tmp;
             }
           }else {
             ZrhoVar1[iPt] += tmp;
           }

           ZgammaVar1[iPt] = 2 * ( VG[3*iPt] - VG[3*iPt + 2]                 );

           ZgammaVar2[iPt] = 2 * ( VG[3*iPt] - VG[3*iPt + 1] + VG[3*iPt + 2] );

           ZgammaVar3[iPt] = 
             ( V2G[6*iPt    ] +   V2G[6*iPt + 1] - V2G[6*iPt + 4] - V2G[6*iPt + 5]   ) 
              * gPTssv[iPt] +
             ( V2G[6*iPt    ] -   V2G[6*iPt + 1] + V2G[6*iPt + 4] - V2G[6*iPt + 5]   ) 
              * (gPTxxv[iPt] + gPTyyv[iPt] + gPTzzv[iPt]) +
             ( V2RG[6*iPt    ] - V2RG[6*iPt + 2] + V2RG[6*iPt + 3] - V2RG[6*iPt + 5] ) 
              * TS[iPt] ;


           ZgammaVar4[iPt] = 
             ( V2G[6*iPt    ] + 2*V2G[6*iPt + 2] -   V2G[6*iPt + 3] +
               V2G[6*iPt + 5]                                         ) 
              * gPTssv[iPt] +
             ( V2G[6*iPt    ] - 2*V2G[6*iPt + 1] + 2*V2G[6*iPt + 2] + 
               V2G[6*iPt + 3] - 2*V2G[6*iPt + 4] +   V2G[6*iPt + 5]   ) 
              * (gPTxxv[iPt] + gPTyyv[iPt] + gPTzzv[iPt]) +
            ( V2RG[6*iPt    ] - V2RG[6*iPt + 1] + V2RG[6*iPt + 2] + 
              V2RG[6*iPt + 3] - V2RG[6*iPt + 4] + V2RG[6*iPt + 5]   ) 
              * TS[iPt] ;



           if (!Msmall[iPt]) {

             ZgammaVar3[iPt] += ( V2G[6*iPt    ] - 2*V2G[6*iPt + 2] + V2G[6*iPt + 5] )
               * (gPTsxv[iPt] * Hx[iPt] + gPTsyv[iPt] * Hy[iPt] + gPTszv[iPt] * Hz[iPt] ) ; 

             ZgammaVar4[iPt] += ( V2G[6*iPt    ] -   V2G[6*iPt + 1] +   V2G[6*iPt + 4] - V2G[6*iPt + 5] )
               * (gPTsxv[iPt] * Hx[iPt] + gPTsyv[iPt] * Hy[iPt] + gPTszv[iPt] * Hz[iPt] );


               ZgammaVar3[iPt] += ( V2RG[6*iPt    ] - V2RG[6*iPt + 2] - V2RG[6*iPt + 3] + V2RG[6*iPt + 5] ) 
                                   * (   Kz[iPt] * TZ[iPt]   
                                       + Ky[iPt] * TY[iPt]    
                                       + Kx[iPt] * TX[iPt]  );
               ZgammaVar4[iPt] += ( V2RG[6*iPt    ] - V2RG[6*iPt + 1] + V2RG[6*iPt + 2] - 
                                    V2RG[6*iPt + 3] + V2RG[6*iPt + 4] - V2RG[6*iPt + 5]   ) 
                                    *  (   Kz[iPt] * TZ[iPt]   
                                        + Ky[iPt] * TY[iPt]    
                                        + Kx[iPt] * TX[iPt]  );

           } else {
             if (denTyp == MX ) { 
               ZgammaVar3[iPt] += ( V2G[6*iPt    ] - 2*V2G[6*iPt + 2] + V2G[6*iPt + 5] ) 
                 * gPTsxv[iPt] ;

               ZgammaVar4[iPt] += ( V2G[6*iPt    ] -   V2G[6*iPt + 1] +   V2G[6*iPt + 4]
                 - V2G[6*iPt + 5] )* gPTsxv[iPt]; 


               ZgammaVar3[iPt] += ( V2RG[6*iPt    ] - V2RG[6*iPt + 2] - V2RG[6*iPt + 3] + V2RG[6*iPt + 5] )
                                      * TX[iPt];

               ZgammaVar4[iPt] += ( V2RG[6*iPt    ] - V2RG[6*iPt + 1] + V2RG[6*iPt + 2] - 
                                    V2RG[6*iPt + 3] + V2RG[6*iPt + 4] - V2RG[6*iPt + 5]   ) 
                                      * TX[iPt];

             } else if (denTyp == MY) {
               ZgammaVar3[iPt] += ( V2G[6*iPt    ] - 2*V2G[6*iPt + 2] + V2G[6*iPt + 5] ) 
                 * gPTsyv[iPt] ; 


               ZgammaVar4[iPt] += ( V2G[6*iPt    ] -   V2G[6*iPt + 1] +   V2G[6*iPt + 4]
                 - V2G[6*iPt + 5] )* gPTsyv[iPt]; 


               ZgammaVar3[iPt] += ( V2RG[6*iPt    ] - V2RG[6*iPt + 2] - V2RG[6*iPt + 3] + V2RG[6*iPt + 5] )
                                      * TY[iPt];

               ZgammaVar4[iPt] += ( V2RG[6*iPt    ] - V2RG[6*iPt + 1] + V2RG[6*iPt + 2] - 
                                    V2RG[6*iPt + 3] + V2RG[6*iPt + 4] - V2RG[6*iPt + 5]   ) 
                                      * TY[iPt];
             } else if (denTyp == MZ) {
               ZgammaVar3[iPt] += ( V2G[6*iPt    ] - 2*V2G[6*iPt + 2] + V2G[6*iPt + 5] ) 
                 * gPTszv[iPt] ; 
               ZgammaVar4[iPt] += ( V2G[6*iPt    ] -   V2G[6*iPt + 1] +   V2G[6*iPt + 4]
                 - V2G[6*iPt + 5] )* gPTszv[iPt]; 


               ZgammaVar3[iPt] += ( V2RG[6*iPt    ] - V2RG[6*iPt + 2] - V2RG[6*iPt + 3] + V2RG[6*iPt + 5] )
                                      * TZ[iPt];


               ZgammaVar4[iPt] += ( V2RG[6*iPt    ] - V2RG[6*iPt + 1] + V2RG[6*iPt + 2] - 
                                    V2RG[6*iPt + 3] + V2RG[6*iPt + 4] - V2RG[6*iPt + 5]   ) 
                                      * TZ[iPt];
             }
           }   //!Msmall  
         } //NPts 
       } //GGA not SCALAR 
    }  // 2C 
  } // end

  template <typename MatsT, typename IntsT>
  template <typename U>
  void KohnSham<MatsT,IntsT>::evalTransDen(SHELL_EVAL_TYPE typ, size_t NPts,size_t NBE, size_t NB, 
    std::vector<std::pair<size_t,size_t>> &subMatCut, U *SCR1,
    U *SCR2, U *DENMAT, U *Den, U *GDenX, U *GDenY, U *GDenZ,
    U *BasisScr){

    size_t IOff = NPts*NBE;

    SubMatSet(NB,NB,NBE,NBE,DENMAT,NB,SCR1,NBE,subMatCut);           

    // Obtain Sum_nu P_mu_nu Phi_nu
    blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,NBE,NPts,NBE,U(1.),SCR1,NBE,BasisScr,NBE,U(0.),SCR2,NBE);

    if( typ != GRADIENT ){
      for(auto iPt = 0; iPt < NPts; iPt++) {
        Den[iPt] = 0.;
        U *SCR_cur = SCR2 + iPt*NBE;
        U *B_cur   = BasisScr + iPt*NBE;
  
        for (size_t j = 0; j < NBE; j++) 
          Den[iPt] += SCR_cur[j] * B_cur[j];

      }
    } else {

      for(auto iPt = 0; iPt < NPts; iPt++) {
        Den[iPt] = 0.;
        GDenX[iPt] = 0.;
        GDenY[iPt] = 0.;
        GDenZ[iPt] = 0.;
        const size_t NBEiPt = iPt*NBE;
        const U *SCR_cur  = SCR2 + NBEiPt;
        const U *B_cur    = BasisScr + NBEiPt;
        const U *B_curX   = B_cur  + IOff;
        const U *B_curY   = B_curX + IOff;
        const U *B_curZ   = B_curY + IOff;
      
        for(size_t j = 0; j < NBE; j++) { 
          Den[iPt]   += SCR_cur[j] * B_cur[j];
          GDenX[iPt] += SCR_cur[j] * B_curX[j];
          GDenY[iPt] += SCR_cur[j] * B_curY[j];
          GDenZ[iPt] += SCR_cur[j] * B_curZ[j];
        }
      }
      // Transition density not symmetric
      // Obtain Sum_nu P^T_mu_nu Phi_nu
      blas::gemm(blas::Layout::ColMajor,blas::Op::Trans,blas::Op::NoTrans,NBE,NPts,NBE,U(1.),SCR1,NBE,BasisScr,NBE,U(0.),SCR2,NBE);
 
      for(auto iPt = 0; iPt < NPts; iPt++) {
        const size_t NBEiPt = iPt*NBE;
        const U *SCR_cur  = SCR2 + NBEiPt;
        const U *B_cur    = BasisScr + NBEiPt;
        const U *B_curX   = B_cur  + IOff;
        const U *B_curY   = B_curX + IOff;
        const U *B_curZ   = B_curY + IOff;

        U denXtmp=0.0;
        U denYtmp=0.0;
        U denZtmp=0.0;
            
        for(size_t j = 0; j < NBE; j++) { 
          denXtmp += SCR_cur[j] * B_curX[j];
          denYtmp += SCR_cur[j] * B_curY[j];
          denZtmp += SCR_cur[j] * B_curZ[j];
        }

        GDenX[iPt] += denXtmp;
        GDenY[iPt] += denYtmp;
        GDenZ[iPt] += denZtmp;

      } //for(auto iPt = 0; iPt < NPts; iPt++) 
    }

  }; // //KohnSham<MatsT,IntsT>::evalTransDen


  template <typename MatsT, typename IntsT>
  template <typename U>
  void KohnSham<MatsT, IntsT>::formZ_fxc(DENSITY_TYPE denType, bool isGGA, size_t NPts, size_t NBE, size_t IOff,
    double epsScreen, std::vector<double> &weights,
    U *ZrhoVar1, U *ZgammaVar1, U *ZgammaVar2, U *ZgammaVar3, U *ZgammaVar4,
    double* GDenS, double* GDenZ, double* GDenY, double* GDenX, U* GTS, U* GTZ, U* GTY, U* GTX,
    double *BasisScr, U* ZMAT) {

    memset(ZMAT,0,IOff*sizeof(U));

    double* ZMAT_RE = reinterpret_cast<double*>(ZMAT);
    double* ZMAT_IM = ZMAT_RE + 1;
    int INCZMAT = sizeof(U)/sizeof(double);

    for(auto iPt = 0ul; iPt < NPts; iPt++) {

      U Fg = 0.5 * weights[iPt] * ZrhoVar1[iPt];

      double Fg_r = std::real(Fg);
      blas::axpy(NBE,Fg_r,BasisScr + iPt*NBE,1, ZMAT_RE + iPt*NBE*INCZMAT ,INCZMAT);

      if( std::is_same<U,dcomplex>::value ) {
        double Fg_i = std::imag(Fg);
        blas::axpy(NBE,Fg_i,BasisScr + iPt*NBE,1, ZMAT_IM + iPt*NBE*INCZMAT ,INCZMAT);
      }



      if( isGGA ) {

        U FgX = ZgammaVar1[iPt] * GTS[iPt         ] + ZgammaVar2[iPt] * GTZ[iPt         ] + ZgammaVar3[iPt] * GDenS[iPt         ];
        U FgY = ZgammaVar1[iPt] * GTS[iPt + NPts  ] + ZgammaVar2[iPt] * GTZ[iPt + NPts  ] + ZgammaVar3[iPt] * GDenS[iPt + NPts  ];
        U FgZ = ZgammaVar1[iPt] * GTS[iPt + 2*NPts] + ZgammaVar2[iPt] * GTZ[iPt + 2*NPts] + ZgammaVar3[iPt] * GDenS[iPt + 2*NPts];

        if( this->onePDM->hasZ() ) {

          FgX += ZgammaVar4[iPt] * GDenZ[iPt         ];
          FgY += ZgammaVar4[iPt] * GDenZ[iPt + NPts  ];
          FgZ += ZgammaVar4[iPt] * GDenZ[iPt + 2*NPts];

        }


        FgX *= weights[iPt];
        FgY *= weights[iPt];
        FgZ *= weights[iPt];

        double FgX_r = std::real(FgX);
        double FgY_r = std::real(FgY);
        double FgZ_r = std::real(FgZ);

        blas::axpy(NBE,FgX_r,BasisScr + iPt*NBE +   IOff,1,ZMAT_RE + iPt*NBE*INCZMAT,INCZMAT);
        blas::axpy(NBE,FgY_r,BasisScr + iPt*NBE + 2*IOff,1,ZMAT_RE + iPt*NBE*INCZMAT,INCZMAT);
        blas::axpy(NBE,FgZ_r,BasisScr + iPt*NBE + 3*IOff,1,ZMAT_RE + iPt*NBE*INCZMAT,INCZMAT);

      }



    }

  }


//SS:start GIAO
  template<> void KohnSham<dcomplex,dcomplex>::formZ_fxc(DENSITY_TYPE denType, bool isGGA, size_t NPts, size_t NBE, size_t IOff,
    double epsScreen, std::vector<double> &weights,
    dcomplex *ZrhoVar1, dcomplex *ZgammaVar1, dcomplex *ZgammaVar2, dcomplex *ZgammaVar3, dcomplex *ZgammaVar4,
    double* GDenS, double* GDenZ, double* GDenY, double* GDenX, dcomplex* GTS, dcomplex* GTZ, dcomplex* GTY, dcomplex* GTX,
    dcomplex *BasisScr, dcomplex* ZMAT) {

    memset(ZMAT,0,IOff*sizeof(dcomplex));

    for(auto iPt = 0ul; iPt < NPts; iPt++) {

      dcomplex Fg = 0.5 * weights[iPt] * ZrhoVar1[iPt];

      blas::axpy(NBE,Fg,BasisScr + iPt*NBE,1, ZMAT + iPt*NBE ,1);

      if( isGGA ) {

        dcomplex FgX = ZgammaVar1[iPt] * GTS[iPt         ] + ZgammaVar2[iPt] * GTZ[iPt         ] + ZgammaVar3[iPt] * GDenS[iPt         ];
        dcomplex FgY = ZgammaVar1[iPt] * GTS[iPt + NPts  ] + ZgammaVar2[iPt] * GTZ[iPt + NPts  ] + ZgammaVar3[iPt] * GDenS[iPt + NPts  ];
        dcomplex FgZ = ZgammaVar1[iPt] * GTS[iPt + 2*NPts] + ZgammaVar2[iPt] * GTZ[iPt + 2*NPts] + ZgammaVar3[iPt] * GDenS[iPt + 2*NPts];

        if( this->onePDM->hasZ() ) {

          FgX += ZgammaVar4[iPt] * GDenZ[iPt         ];
          FgY += ZgammaVar4[iPt] * GDenZ[iPt + NPts  ];
          FgZ += ZgammaVar4[iPt] * GDenZ[iPt + 2*NPts];

        }


        FgX *= weights[iPt];
        FgY *= weights[iPt];
        FgZ *= weights[iPt];

        blas::axpy(NBE,FgX,BasisScr + iPt*NBE +   IOff,1,ZMAT + iPt*NBE,1);
        blas::axpy(NBE,FgY,BasisScr + iPt*NBE + 2*IOff,1,ZMAT + iPt*NBE,1);
        blas::axpy(NBE,FgZ,BasisScr + iPt*NBE + 3*IOff,1,ZMAT + iPt*NBE,1);

      } // if( isGGA )



    }

  }

//SS:end

  template <typename MatsT, typename IntsT>
  template <typename U>
  void KohnSham<MatsT, IntsT>::formZ_fxc(DENSITY_TYPE denType, bool isGGA, size_t NPts, size_t NBE, size_t IOff,
    double epsScreen, std::vector<double> &weights,
    U *ZrhoVar1, U *ZgammaVar1, U *ZgammaVar2, U *ZgammaVar3, U *ZgammaVar4,
    bool * Msmall, double *Mnorm, 
    double* DSDMnorm, double* signMD,
    double* GDenS, double* GDenZ, double* GDenY, double* GDenX, 
    double *Kx, double *Ky, double *Kz, double *Hx, double *Hy, double *Hz,
    U* GTS, U* GTZ, U* GTY, U* GTX, 
    U* gPTss, U* gPTsz, U* gPTsy, U* gPTsx, U* gPTzz,
    U* gPTyy, U* gPTxx,  
    double *BasisScr, U* ZMAT) {

    memset(ZMAT,0,IOff*sizeof(U));

    //The same for 1 or 2C
    for(auto iPt = 0ul; iPt < NPts; iPt++) {
        U Fg = 0.5 * weights[iPt] * ZrhoVar1[iPt];
        blas::axpy(NBE,Fg,BasisScr + iPt*NBE,1, ZMAT + iPt*NBE ,1);
  
      }

      if( isGGA ) {

        if( not this->onePDM->hasXY() ) {
          for(auto iPt = 0ul; iPt < NPts; iPt++) {
          
            U FgX = ZgammaVar1[iPt] * GTS[iPt         ] + ZgammaVar2[iPt] * GTZ[iPt         ] + ZgammaVar3[iPt] * GDenS[iPt         ];
            U FgY = ZgammaVar1[iPt] * GTS[iPt + NPts  ] + ZgammaVar2[iPt] * GTZ[iPt + NPts  ] + ZgammaVar3[iPt] * GDenS[iPt + NPts  ];
            U FgZ = ZgammaVar1[iPt] * GTS[iPt + 2*NPts] + ZgammaVar2[iPt] * GTZ[iPt + 2*NPts] + ZgammaVar3[iPt] * GDenS[iPt + 2*NPts];
  
            if( this->onePDM->hasZ() ) {
  
              FgX += ZgammaVar4[iPt] * GDenZ[iPt         ];
              FgY += ZgammaVar4[iPt] * GDenZ[iPt + NPts  ];
              FgZ += ZgammaVar4[iPt] * GDenZ[iPt + 2*NPts];
  
            }
  
  
            FgX *= weights[iPt];
            FgY *= weights[iPt];
            FgZ *= weights[iPt];
  
            blas::axpy(NBE,FgX,BasisScr + iPt*NBE +   IOff,1,ZMAT + iPt*NBE,1);
            blas::axpy(NBE,FgY,BasisScr + iPt*NBE + 2*IOff,1,ZMAT + iPt*NBE,1);
            blas::axpy(NBE,FgZ,BasisScr + iPt*NBE + 3*IOff,1,ZMAT + iPt*NBE,1);
  
          } // NPTS

        } else {
          for(auto iPt = 0ul; iPt < NPts; iPt++) {

            U FgX = U(0.);
            U FgY = U(0.);
            U FgZ = U(0.);

            U tmp1 = U(0.0);
            U tmp2 = U(0.0);

            if (denType == SCALAR ) {
             // Fx
               if (!Msmall[iPt]) {
                  U tmp1 =                (   GDenX [iPt         ] * gPTsx[iPt] 
                                            + GDenY [iPt         ] * gPTsy[iPt]  
                                            + GDenZ [iPt         ] * gPTsz[iPt] 
                                          ) / DSDMnorm[iPt]; 
                  U tmp2 =                (  gPTsx[iPt] * Hx[iPt] 
                                           + gPTsy[iPt] * Hy[iPt] 
                                           + gPTsz[iPt] * Hz[iPt] ) * (   GDenX [iPt         ] * Hx[iPt]
                                                                   + GDenY [iPt         ] * Hy[iPt]  
                                                                   + GDenZ [iPt         ] * Hz[iPt] 
                                                                                                    ) / DSDMnorm[iPt]; 
                }
              FgX = 
                    ZgammaVar1[iPt] *    GTS   [iPt         ] 
                  + ZgammaVar2[iPt] * (  GTX   [iPt         ] * Hx[iPt] 
                                       + GTY   [iPt         ] * Hy[iPt]  
                                       + GTZ   [iPt         ] * Hz[iPt] ) 
                  + ZgammaVar3[iPt] *    GDenS [iPt         ]
                  + ZgammaVar4[iPt] * (  GDenX [iPt         ] * Hx[iPt] 
                                       + GDenY [iPt         ] * Hy[iPt]  
                                       + GDenZ [iPt         ] * Hz[iPt] ) 
                  + ZgammaVar2[iPt] * (tmp1 - tmp2) *signMD[iPt];

             // Fy
               tmp1 = U(0.0);
               tmp2 = U(0.0);
               if (!Msmall[iPt]) {
                tmp1 =                (   GDenX [iPt + NPts         ] * gPTsx[iPt] 
                                        + GDenY [iPt + NPts         ] * gPTsy[iPt]  
                                        + GDenZ [iPt + NPts         ] * gPTsz[iPt] 
                                      ) / DSDMnorm[iPt]; 
                tmp2 =                (  gPTsx[iPt] * Hx[iPt] 
                                       + gPTsy[iPt] * Hy[iPt] 
                                       + gPTsz[iPt] * Hz[iPt] ) * (   GDenX [iPt + NPts        ] * Hx[iPt]
                                                               + GDenY [iPt + NPts         ] * Hy[iPt]  
                                                               + GDenZ [iPt + NPts         ] * Hz[iPt] 
                                                                                                ) / DSDMnorm[iPt]; 
                }
              FgY = 
                    ZgammaVar1[iPt] *    GTS   [iPt + NPts     ] 
                  + ZgammaVar2[iPt] * (  GTX   [iPt + NPts     ] * Hx[iPt] 
                                       + GTY   [iPt + NPts     ] * Hy[iPt]  
                                       + GTZ   [iPt + NPts     ] * Hz[iPt] ) 
                  + ZgammaVar3[iPt] *    GDenS [iPt + NPts     ]
                  + ZgammaVar4[iPt] * (  GDenX [iPt + NPts     ] * Hx[iPt] 
                                       + GDenY [iPt + NPts     ] * Hy[iPt]  
                                       + GDenZ [iPt + NPts     ] * Hz[iPt] ) 
                  + ZgammaVar2[iPt] * (tmp1 - tmp2) *signMD[iPt];
             // Fz
               tmp1 = U(0.0);
               tmp2 = U(0.0);
               if (!Msmall[iPt]) {
                tmp1 =                (   GDenX [iPt + 2*NPts         ] * gPTsx[iPt] 
                                        + GDenY [iPt + 2*NPts         ] * gPTsy[iPt]  
                                        + GDenZ [iPt + 2*NPts         ] * gPTsz[iPt] 
                                      ) / DSDMnorm[iPt]; 
                tmp2 =                (  gPTsx[iPt] * Hx[iPt] 
                                       + gPTsy[iPt] * Hy[iPt] 
                                       + gPTsz[iPt] * Hz[iPt] ) * (   GDenX [iPt + 2*NPts        ] * Hx[iPt]
                                                               + GDenY [iPt + 2*NPts         ] * Hy[iPt]  
                                                               + GDenZ [iPt + 2*NPts         ] * Hz[iPt] 
                                                                                                ) / DSDMnorm[iPt]; 
                }
              FgZ = 
                    ZgammaVar1[iPt] *    GTS   [iPt + 2*NPts     ] 
                  + ZgammaVar2[iPt] * (  GTX   [iPt + 2*NPts     ] * Hx[iPt] 
                                       + GTY   [iPt + 2*NPts     ] * Hy[iPt]  
                                       + GTZ   [iPt + 2*NPts     ] * Hz[iPt] ) 
                  + ZgammaVar3[iPt] *    GDenS [iPt + 2*NPts     ]
                  + ZgammaVar4[iPt] * (  GDenX [iPt + 2*NPts     ] * Hx[iPt] 
                                       + GDenY [iPt + 2*NPts     ] * Hy[iPt]  
                                       + GDenZ [iPt + 2*NPts     ] * Hz[iPt] ) 
                  + ZgammaVar2[iPt] * (tmp1 - tmp2) *signMD[iPt];


            } else if (denType == MX ) {
              //Fx
               if (!Msmall[iPt]) {
               tmp1 =                (   GDenS [iPt         ] * gPTsx[iPt] 
                                      ) / DSDMnorm[iPt]; 
               tmp2 =                (  gPTsx[iPt] * Hx[iPt] 
                                       + gPTsy[iPt] * Hy[iPt] 
                                       + gPTsz[iPt] * Hz[iPt] ) * (   GDenS [iPt         ] * Hx[iPt]
                                                                                                ) / DSDMnorm[iPt]; 
                
                FgX = 
                      ZgammaVar1[iPt] *    GTS   [iPt         ] * Hx[iPt]
                    + ZgammaVar2[iPt] * (  GTX   [iPt         ] ) 
                    + ZgammaVar3[iPt] *    GDenS [iPt         ] * Hx[iPt]
                    + ZgammaVar4[iPt] * (  GDenX [iPt         ] )
                  + ZgammaVar1[iPt] * (tmp1 - tmp2) *signMD[iPt];
              } else {
                FgX = ZgammaVar1[iPt] *    GTS   [iPt         ]
                      + ZgammaVar2[iPt] * (  GTX   [iPt         ] )
                      + ZgammaVar3[iPt] *    GDenS [iPt         ] 
                    + ZgammaVar4[iPt] * (  GDenX [iPt         ] ); 
              }
              //Fy
               tmp1 = U(0.0);
               tmp2 = U(0.0);
               if (!Msmall[iPt]) {
                tmp1 =                (   GDenS [iPt + NPts      ] * gPTsx[iPt] 
                                      ) / DSDMnorm[iPt]; 
                tmp2 =                (  gPTsx[iPt] * Hx[iPt] 
                                       + gPTsy[iPt] * Hy[iPt] 
                                       + gPTsz[iPt] * Hz[iPt] ) * (   GDenS [iPt + NPts   ] * Hx[iPt]
                                                                                                ) / DSDMnorm[iPt]; 
                
                 FgY = 
                       ZgammaVar1[iPt] *    GTS   [iPt + NPts      ] * Hx[iPt]
                     + ZgammaVar2[iPt] * (  GTX   [iPt + NPts      ] ) 
                     + ZgammaVar3[iPt] *    GDenS [iPt + NPts      ] * Hx[iPt]
                     + ZgammaVar4[iPt] * (  GDenX [iPt + NPts      ] )
                  + ZgammaVar1[iPt] * (tmp1 - tmp2) *signMD[iPt];
               } else {
                 FgY = ZgammaVar1[iPt] *    GTS   [iPt + NPts      ]
                       + ZgammaVar2[iPt] * (  GTX   [iPt + NPts      ] )  
                       + ZgammaVar3[iPt] *    GDenS [iPt + NPts        ] 
                      + ZgammaVar4[iPt] * (  GDenX [iPt + NPts      ] ) ;
               }

              //Fz
               tmp1 = U(0.0);
               tmp2 = U(0.0);
               if (!Msmall[iPt]) {
                tmp1 =                (   GDenS [iPt + 2*NPts      ] * gPTsx[iPt] 
                                      ) / DSDMnorm[iPt]; 
                tmp2 =                (  gPTsx[iPt] * Hx[iPt] 
                                       + gPTsy[iPt] * Hy[iPt] 
                                       + gPTsz[iPt] * Hz[iPt] ) * (   GDenS [iPt + 2*NPts       ] * Hx[iPt]
                                                                                                ) / DSDMnorm[iPt]; 
                
              FgZ = 
                    ZgammaVar1[iPt] *    GTS   [iPt + 2*NPts      ] * Hx[iPt] 
                  + ZgammaVar2[iPt] * (  GTX   [iPt + 2*NPts      ] )
                  + ZgammaVar3[iPt] *    GDenS [iPt + 2*NPts      ] * Hx[iPt]
                  + ZgammaVar4[iPt] * (  GDenX [iPt + 2*NPts      ] )
                  + ZgammaVar1[iPt] * (tmp1 - tmp2) *signMD[iPt];
               } else {
                 FgZ = ZgammaVar1[iPt] *    GTS   [iPt + 2*NPts      ]
                     + ZgammaVar2[iPt] * (  GTX   [iPt + 2*NPts      ] )
                     + ZgammaVar3[iPt] *    GDenS [iPt + 2*NPts        ] 
                     + ZgammaVar4[iPt] * (  GDenX [iPt + 2*NPts      ] ) ;
               }


            } else if (denType == MY ) {
              //Fx
               if (!Msmall[iPt]) {
               tmp1 =                (   GDenS [iPt         ] * gPTsy[iPt] 
                                          ) / DSDMnorm[iPt]; 
               tmp2 =                (  gPTsx[iPt] * Hx[iPt] 
                                       + gPTsy[iPt] * Hy[iPt] 
                                       + gPTsz[iPt] * Hz[iPt] ) * (   GDenS [iPt         ] * Hy[iPt]
                                                                                                ) / DSDMnorm[iPt]; 
               
              FgX = 
                    ZgammaVar1[iPt] *    GTS   [iPt         ] * Hy[iPt]
                  + ZgammaVar2[iPt] * (  GTY   [iPt         ] )
                  + ZgammaVar3[iPt] *    GDenS [iPt         ] * Hy[iPt]
                  + ZgammaVar4[iPt] * (  GDenY [iPt         ] )
                  + ZgammaVar1[iPt] * (tmp1 - tmp2) *signMD[iPt];

               } else {
                 FgX = ZgammaVar1[iPt] *    GTS   [iPt         ] 
                   + ZgammaVar2[iPt] * (  GTY   [iPt         ] ) 
                   + ZgammaVar3[iPt] *    GDenS [iPt         ]
                   + ZgammaVar4[iPt] * (  GDenY [iPt         ] );
               }
              //Fy
               tmp1 = U(0.0);
               tmp2 = U(0.0);
               if (!Msmall[iPt]) {
                tmp1 =                (   GDenS [iPt + NPts      ] * gPTsy[iPt] 
                                      ) / DSDMnorm[iPt]; 
                tmp2 =                (  gPTsx[iPt] * Hx[iPt] 
                                       + gPTsy[iPt] * Hy[iPt] 
                                       + gPTsz[iPt] * Hz[iPt] ) * (   GDenS [iPt + NPts       ] * Hy[iPt]
                                                                                                ) / DSDMnorm[iPt]; 
                
              FgY = 
                    ZgammaVar1[iPt] *    GTS   [iPt + NPts      ] * Hy[iPt]
                  + ZgammaVar2[iPt] * (  GTY   [iPt + NPts      ] )
                  + ZgammaVar3[iPt] *    GDenS [iPt + NPts      ] * Hy[iPt]
                  + ZgammaVar4[iPt] * (  GDenY [iPt + NPts      ] )
                  + ZgammaVar1[iPt] * (tmp1 - tmp2) *signMD[iPt];
               } else {
                 FgY = ZgammaVar1[iPt] *    GTS   [iPt + NPts      ]
                   + ZgammaVar2[iPt] * (  GTY   [iPt + NPts      ] )  
                   + ZgammaVar3[iPt] *    GDenS [iPt + NPts      ]
                   + ZgammaVar4[iPt] * (  GDenY [iPt + NPts      ] ) ;
               } 
              //Fz
               tmp1 = U(0.0);
               tmp2 = U(0.0);
               if (!Msmall[iPt]) {
                tmp1 =                (   GDenS [iPt + 2*NPts      ] * gPTsy[iPt] 
                                      ) / DSDMnorm[iPt]; 
                tmp2 =                (  gPTsx[iPt] * Hx[iPt] 
                                       + gPTsy[iPt] * Hy[iPt] 
                                       + gPTsz[iPt] * Hz[iPt] ) * (   GDenS [iPt + 2*NPts      ] * Hy[iPt]
                                                                                                ) / DSDMnorm[iPt]; 
               
              FgZ = 
                    ZgammaVar1[iPt] *    GTS   [iPt + 2*NPts      ] * Hy[iPt] 
                  + ZgammaVar2[iPt] * (  GTY   [iPt + 2*NPts      ] )
                  + ZgammaVar3[iPt] *    GDenS [iPt + 2*NPts      ] * Hy[iPt]
                  + ZgammaVar4[iPt] * (  GDenY [iPt + 2*NPts      ] )
                  + ZgammaVar1[iPt] * (tmp1 - tmp2) *signMD[iPt];
                } else {
                  FgZ = ZgammaVar1[iPt] *    GTS   [iPt + 2*NPts      ]
                    + ZgammaVar2[iPt] * (  GTY   [iPt + 2*NPts      ] )
                    + ZgammaVar3[iPt] *    GDenS [iPt + 2*NPts      ]
                    + ZgammaVar4[iPt] * (  GDenY [iPt + 2*NPts      ] ) ;
                }


            } else if (denType == MZ ) {
              //Fx
               if (!Msmall[iPt]) {
               tmp1 =                (   GDenS [iPt         ] * gPTsz[iPt] 
                                      ) / DSDMnorm[iPt]; 
               tmp2 =                (  gPTsx[iPt] * Hx[iPt] 
                                       + gPTsy[iPt] * Hy[iPt] 
                                       + gPTsz[iPt] * Hz[iPt] ) * (   GDenS [iPt         ] * Hz[iPt]
                                                                                                ) / DSDMnorm[iPt]; 
               
              FgX = 
                    ZgammaVar1[iPt] *    GTS   [iPt         ] * Hz[iPt]
                  + ZgammaVar2[iPt] * (  GTZ   [iPt         ] )
                  + ZgammaVar3[iPt] *    GDenS [iPt         ] * Hz[iPt]
                  + ZgammaVar4[iPt] * (  GDenZ [iPt         ] )
                  + ZgammaVar1[iPt] * (tmp1 - tmp2) *signMD[iPt];
               } else {
                FgX = ZgammaVar1[iPt] *    GTS   [iPt         ]
                  + ZgammaVar2[iPt] * (  GTZ   [iPt         ] )
                  + ZgammaVar3[iPt] *    GDenS [iPt      ]
                  + ZgammaVar4[iPt] * (  GDenZ [iPt         ] );
               }

              //Fy
               tmp1 = U(0.0);
               tmp2 = U(0.0);
               if (!Msmall[iPt]) {
                tmp1 =                (   GDenS [iPt + NPts      ] * gPTsz[iPt] 
                                      ) / DSDMnorm[iPt]; 
                tmp2 =                (  gPTsx[iPt] * Hx[iPt] 
                                       + gPTsy[iPt] * Hy[iPt] 
                                       + gPTsz[iPt] * Hz[iPt] ) * (   GDenS [iPt + NPts      ] * Hz[iPt]
                                                                                                ) / DSDMnorm[iPt]; 
               
              FgY = 
                    ZgammaVar1[iPt] *    GTS   [iPt + NPts      ] * Hz[iPt]
                  + ZgammaVar2[iPt] * (  GTZ   [iPt + NPts      ] )
                  + ZgammaVar3[iPt] *    GDenS [iPt + NPts      ] * Hz[iPt]
                  + ZgammaVar4[iPt] * (  GDenZ [iPt + NPts      ] )
                  + ZgammaVar1[iPt] * (tmp1 - tmp2) *signMD[iPt];

               } else {
                 FgY = ZgammaVar1[iPt] *    GTS   [iPt + NPts      ]
                    + ZgammaVar2[iPt] * (  GTZ   [iPt + NPts      ] )
                    + ZgammaVar3[iPt] *    GDenS [iPt + NPts      ]
                    + ZgammaVar4[iPt] * (  GDenZ [iPt + NPts      ] ); 
               } 
              //Fz
               tmp1 = U(0.0);
               tmp2 = U(0.0);
               if (!Msmall[iPt]) {
                tmp1 =                (   GDenS [iPt + 2*NPts      ] * gPTsz[iPt] 
                                      ) / DSDMnorm[iPt]; 
                tmp2 =                (  gPTsx[iPt] * Hx[iPt] 
                                       + gPTsy[iPt] * Hy[iPt] 
                                       + gPTsz[iPt] * Hz[iPt] ) * (   GDenS [iPt + 2*NPts    ] * Hz[iPt]
                                                                                                ) / DSDMnorm[iPt]; 
                
              FgZ = 
                    ZgammaVar1[iPt] *    GTS   [iPt + 2*NPts      ] * Hz[iPt]
                  + ZgammaVar2[iPt] * (  GTZ   [iPt + 2*NPts      ] )
                  + ZgammaVar3[iPt] *    GDenS [iPt + 2*NPts      ] * Hz[iPt]
                  + ZgammaVar4[iPt] * (  GDenZ [iPt + 2*NPts      ] )
                  + ZgammaVar1[iPt] * (tmp1 - tmp2) *signMD[iPt];
                } else {
                  FgZ = ZgammaVar1[iPt] *    GTS   [iPt + 2*NPts      ]
                     + ZgammaVar2[iPt] * (  GTZ   [iPt + 2*NPts      ] )
                     + ZgammaVar3[iPt] *    GDenS [iPt + 2*NPts      ]
                     + ZgammaVar4[iPt] * (  GDenZ [iPt + 2*NPts      ] )  ; 
                }  

            }

            FgX *= weights[iPt];
            FgY *= weights[iPt];
            FgZ *= weights[iPt];
    
            blas::axpy(NBE,FgX,BasisScr + iPt*NBE +   IOff,1,ZMAT + iPt*NBE,1);
            blas::axpy(NBE,FgY,BasisScr + iPt*NBE + 2*IOff,1,ZMAT + iPt*NBE,1);
            blas::axpy(NBE,FgZ,BasisScr + iPt*NBE + 3*IOff,1,ZMAT + iPt*NBE,1);

          } //NPTS 

        } //2C GGA

    }  // isGGA

  } //End


//2c s

  template<> void KohnSham<dcomplex, dcomplex>::formZ_fxc(DENSITY_TYPE denType, bool isGGA, size_t NPts, size_t NBE, size_t IOff,
    double epsScreen, std::vector<double> &weights,
    dcomplex *ZrhoVar1, dcomplex *ZgammaVar1, dcomplex *ZgammaVar2, dcomplex *ZgammaVar3, dcomplex *ZgammaVar4,
    bool * Msmall, double *Mnorm, 
    // double* DenS, double* DenZ, double* DenY, double* DenX, 
    double* DSDMnorm, double* signMD,
    double* GDenS, double* GDenZ, double* GDenY, double* GDenX, 
    double *Kx, double *Ky, double *Kz, double *Hx, double *Hy, double *Hz,
    dcomplex* GTS, dcomplex* GTZ, dcomplex* GTY, dcomplex* GTX, 
    dcomplex* gPTss, dcomplex* gPTsz, dcomplex* gPTsy, dcomplex* gPTsx, dcomplex* gPTzz,
    dcomplex* gPTyy, dcomplex* gPTxx,  
    dcomplex *BasisScr, dcomplex* ZMAT) {

    memset(ZMAT,0,IOff*sizeof(dcomplex));

//    double* ZMAT_RE = reinterpret_cast<double*>(ZMAT);
//    double* ZMAT_IM = ZMAT_RE + 1;
//    int INCZMAT = sizeof(U)/sizeof(double);
    //The same for 1 or 2C
    for(auto iPt = 0ul; iPt < NPts; iPt++) {
        dcomplex Fg = 0.5 * weights[iPt] * ZrhoVar1[iPt];
        //std::cerr << "Fg " << Fg << std::endl;
        //double Fg_r = std::real(Fg);
        blas::axpy(NBE,Fg,BasisScr + iPt*NBE,1, ZMAT + iPt*NBE ,1);
  /*
        if( std::is_same<U,dcomplex>::value ) {
          double Fg_i = std::imag(Fg);
          AXPY(NBE,Fg_i,BasisScr + iPt*NBE,1, ZMAT_IM + iPt*NBE*INCZMAT ,INCZMAT);
        }
  */
      }

//SS debug s

      if( isGGA ) {

        if( not this->onePDM->hasXY() ) {
          for(auto iPt = 0ul; iPt < NPts; iPt++) {
          

            dcomplex FgX = ZgammaVar1[iPt] * GTS[iPt         ] + ZgammaVar2[iPt] * GTZ[iPt         ] + ZgammaVar3[iPt] * GDenS[iPt         ];
            dcomplex FgY = ZgammaVar1[iPt] * GTS[iPt + NPts  ] + ZgammaVar2[iPt] * GTZ[iPt + NPts  ] + ZgammaVar3[iPt] * GDenS[iPt + NPts  ];
            dcomplex FgZ = ZgammaVar1[iPt] * GTS[iPt + 2*NPts] + ZgammaVar2[iPt] * GTZ[iPt + 2*NPts] + ZgammaVar3[iPt] * GDenS[iPt + 2*NPts];
  
            if( this->onePDM->hasZ() ) {
  
              FgX += ZgammaVar4[iPt] * GDenZ[iPt         ];
              FgY += ZgammaVar4[iPt] * GDenZ[iPt + NPts  ];
              FgZ += ZgammaVar4[iPt] * GDenZ[iPt + 2*NPts];
  
            }
  
  
            FgX *= weights[iPt];
            FgY *= weights[iPt];
            FgZ *= weights[iPt];
  
            // double FgX_r = std::real(FgX);
            // double FgY_r = std::real(FgY);
            // double FgZ_r = std::real(FgZ);
  
            blas::axpy(NBE,FgX,BasisScr + iPt*NBE +   IOff,1,ZMAT + iPt*NBE,1);
            blas::axpy(NBE,FgY,BasisScr + iPt*NBE + 2*IOff,1,ZMAT + iPt*NBE,1);
            blas::axpy(NBE,FgZ,BasisScr + iPt*NBE + 3*IOff,1,ZMAT + iPt*NBE,1);
          } // NPTS

        } else {
          for(auto iPt = 0ul; iPt < NPts; iPt++) {

/*
            dcomplex gPTss = (GDenS[iPt]          * GTS[iPt]         ) +
                      (GDenS[iPt + NPts]   * GTS[iPt + NPts]  ) +
                      (GDenS[iPt + 2*NPts] * GTS[iPt + 2*NPts]);
        
            dcomplex gPTsz = (GDenS[iPt]          * GTZ[iPt]         ) +
                      (GDenS[iPt + NPts]   * GTZ[iPt + NPts]  ) +
                      (GDenS[iPt + 2*NPts] * GTZ[iPt + 2*NPts]);
             //since gPTsz + gPTzs are always added 
              gPTsz += (GDenZ[iPt]          * GTS[iPt]         ) +
                       (GDenZ[iPt + NPts]   * GTS[iPt + NPts]  ) +
                       (GDenZ[iPt + 2*NPts] * GTS[iPt + 2*NPts]);
        
            dcomplex gPTzz = (GDenZ[iPt]          * GTZ[iPt]         ) +
                      (GDenZ[iPt + NPts]   * GTZ[iPt + NPts]  ) +
                      (GDenZ[iPt + 2*NPts] * GTZ[iPt + 2*NPts]);
         
            dcomplex gPTsy = (GDenS[iPt]          * GTY[iPt]         ) +
                      (GDenS[iPt + NPts]   * GTY[iPt + NPts]  ) +
                      (GDenS[iPt + 2*NPts] * GTY[iPt + 2*NPts]);
             //since gPTsy + gPTys are always added 
              gPTsy += (GDenY[iPt]          * GTS[iPt]         ) +
                       (GDenY[iPt + NPts]   * GTS[iPt + NPts]  ) +
                       (GDenY[iPt + 2*NPts] * GTS[iPt + 2*NPts]);
        
            dcomplex gPTyy = (GDenY[iPt]          * GTY[iPt]         ) +
                      (GDenY[iPt + NPts]   * GTY[iPt + NPts]  ) +
                      (GDenY[iPt + 2*NPts] * GTY[iPt + 2*NPts]);
         
            dcomplex gPTsx = (GDenS[iPt]          * GTX[iPt]         ) +
                      (GDenS[iPt + NPts]   * GTX[iPt + NPts]  ) +
                      (GDenS[iPt + 2*NPts] * GTX[iPt + 2*NPts]);
             //since gPTsx + gPTxs are always added 
              gPTsx += (GDenX[iPt]          * GTS[iPt]         ) +
                       (GDenX[iPt + NPts]   * GTS[iPt + NPts]  ) +
                       (GDenX[iPt + 2*NPts] * GTS[iPt + 2*NPts]);
        
            dcomplex gPTxx = (GDenX[iPt]          * GTX[iPt]         ) +
                      (GDenX[iPt + NPts]   * GTX[iPt + NPts]  ) +
                      (GDenX[iPt + 2*NPts] * GTX[iPt + 2*NPts]);
*/

            dcomplex FgX = dcomplex(0.);
            dcomplex FgY = dcomplex(0.);
            dcomplex FgZ = dcomplex(0.);
/*
              double tmpnMx       = GDenX[iPt]*GDenS[iPt];
                     tmpnMx      += GDenX[iPt + NPts]*GDenS[iPt + NPts];
                     tmpnMx      += GDenX[iPt + 2*NPts]*GDenS[iPt + 2*NPts];
              
              double tmpnMy       = GDenY[iPt]*GDenS[iPt];
                     tmpnMy      += GDenY[iPt + NPts]*GDenS[iPt + NPts];
                     tmpnMy      += GDenY[iPt + 2*NPts]*GDenS[iPt + 2*NPts];
              
              double tmpnMz       = GDenZ[iPt]*GDenS[iPt];
                     tmpnMz      += GDenZ[iPt + NPts]*GDenS[iPt + NPts];
                     tmpnMz      += GDenZ[iPt + 2*NPts]*GDenS[iPt + 2*NPts];
              double DSDMnorm  = std::sqrt(tmpnMx * tmpnMx + tmpnMy * tmpnMy + tmpnMz * tmpnMz);

          double signMD = 1. ;

          double tmpSign   = tmpnMx * DenX[iPt]; // is M defined here?
          tmpSign  += tmpnMy * DenY[iPt];
          tmpSign  += tmpnMz * DenZ[iPt];
          if ( std::signbit(tmpSign) ) signMD = -1.;
*/

              dcomplex tmp1 = dcomplex(0.0);
              dcomplex tmp2 = dcomplex(0.0);


            if (denType == SCALAR ) {
            //if (false) {
             // Fx
               if (!Msmall[iPt]) {
                  dcomplex tmp1 =         (   GDenX [iPt         ] * gPTsx[iPt] 
                                            + GDenY [iPt         ] * gPTsy[iPt]  
                                            + GDenZ [iPt         ] * gPTsz[iPt] 
                                          ) / DSDMnorm[iPt]; 
                  dcomplex tmp2 =         (  gPTsx[iPt] * Hx[iPt] 
                                           + gPTsy[iPt] * Hy[iPt] 
                                           + gPTsz[iPt] * Hz[iPt] ) * (   GDenX [iPt         ] * Hx[iPt]
                                                                   + GDenY [iPt         ] * Hy[iPt]  
                                                                   + GDenZ [iPt         ] * Hz[iPt] 
                                                                                                    ) / DSDMnorm[iPt]; 
//tmp1 = 0.0;
//tmp2 = 0.0;
                
              }

              FgX = 
                    ZgammaVar1[iPt] *    GTS   [iPt         ] 
                  + ZgammaVar2[iPt] * (  GTX   [iPt         ] * Hx[iPt] 
                                       + GTY   [iPt         ] * Hy[iPt]  
                                       + GTZ   [iPt         ] * Hz[iPt] ) 
                  + ZgammaVar3[iPt] *    GDenS [iPt         ]
                  + ZgammaVar4[iPt] * (  GDenX [iPt         ] * Hx[iPt] 
                                       + GDenY [iPt         ] * Hy[iPt]  
                                       + GDenZ [iPt         ] * Hz[iPt] ) 
                  + ZgammaVar2[iPt] * (tmp1 - tmp2) *signMD[iPt];
/*
               } else {

              FgX = 
                    ZgammaVar1[iPt] *    GTS   [iPt         ] 
                  + ZgammaVar2[iPt] * (  GTX   [iPt         ] * Hx[iPt] 
                                       + GTY   [iPt         ] * Hy[iPt]  
                                       + GTZ   [iPt         ] * Hz[iPt] ) 
                  + ZgammaVar3[iPt] *    GDenS [iPt         ]
                  + ZgammaVar4[iPt] * (  GDenX [iPt         ] * Hx[iPt] 
                                       + GDenY [iPt         ] * Hy[iPt]  
                                       + GDenZ [iPt         ] * Hz[iPt] ); 


               }
*/
             // Fy
               tmp1 = dcomplex(0.0);
               tmp2 = dcomplex(0.0);
               if (!Msmall[iPt]) {
                tmp1 =                (   GDenX [iPt + NPts         ] * gPTsx[iPt] 
                                        + GDenY [iPt + NPts         ] * gPTsy[iPt]  
                                        + GDenZ [iPt + NPts         ] * gPTsz[iPt] 
                                      ) / DSDMnorm[iPt]; 
                tmp2 =                (  gPTsx[iPt] * Hx[iPt] 
                                       + gPTsy[iPt] * Hy[iPt] 
                                       + gPTsz[iPt] * Hz[iPt] ) * (   GDenX [iPt + NPts        ] * Hx[iPt]
                                                               + GDenY [iPt + NPts         ] * Hy[iPt]  
                                                               + GDenZ [iPt + NPts         ] * Hz[iPt] 
                                                                                                ) / DSDMnorm[iPt]; 
//tmp1 = 0.0;
//tmp2 = 0.0;

                }
              FgY = 
                    ZgammaVar1[iPt] *    GTS   [iPt + NPts     ] 
                  + ZgammaVar2[iPt] * (  GTX   [iPt + NPts     ] * Hx[iPt] 
                                       + GTY   [iPt + NPts     ] * Hy[iPt]  
                                       + GTZ   [iPt + NPts     ] * Hz[iPt] ) 
                  + ZgammaVar3[iPt] *    GDenS [iPt + NPts     ]
                  + ZgammaVar4[iPt] * (  GDenX [iPt + NPts     ] * Hx[iPt] 
                                       + GDenY [iPt + NPts     ] * Hy[iPt]  
                                       + GDenZ [iPt + NPts     ] * Hz[iPt] ) 
                  + ZgammaVar2[iPt] * (tmp1 - tmp2) *signMD[iPt];
             // Fz
               tmp1 = dcomplex(0.0);
               tmp2 = dcomplex(0.0);
               if (!Msmall[iPt]) {
                tmp1 =                (   GDenX [iPt + 2*NPts         ] * gPTsx[iPt] 
                                        + GDenY [iPt + 2*NPts         ] * gPTsy[iPt]  
                                        + GDenZ [iPt + 2*NPts         ] * gPTsz[iPt] 
                                      ) / DSDMnorm[iPt]; 
                tmp2 =                (  gPTsx[iPt] * Hx[iPt] 
                                       + gPTsy[iPt] * Hy[iPt] 
                                       + gPTsz[iPt] * Hz[iPt] ) * (   GDenX [iPt + 2*NPts        ] * Hx[iPt]
                                                               + GDenY [iPt + 2*NPts         ] * Hy[iPt]  
                                                               + GDenZ [iPt + 2*NPts         ] * Hz[iPt] 
                                                                                                ) / DSDMnorm[iPt]; 
//tmp1 =0.0;
//tmp2 = 0.0;
                }
              FgZ = 
                    ZgammaVar1[iPt] *    GTS   [iPt + 2*NPts     ] 
                  + ZgammaVar2[iPt] * (  GTX   [iPt + 2*NPts     ] * Hx[iPt] 
                                       + GTY   [iPt + 2*NPts     ] * Hy[iPt]  
                                       + GTZ   [iPt + 2*NPts     ] * Hz[iPt] ) 
                  + ZgammaVar3[iPt] *    GDenS [iPt + 2*NPts     ]
                  + ZgammaVar4[iPt] * (  GDenX [iPt + 2*NPts     ] * Hx[iPt] 
                                       + GDenY [iPt + 2*NPts     ] * Hy[iPt]  
                                       + GDenZ [iPt + 2*NPts     ] * Hz[iPt] ) 
                  + ZgammaVar2[iPt] * (tmp1 - tmp2) *signMD[iPt];


            } else if (denType == MX ) {
            //} else if (!Msmall[iPt] ) {
              //Fx
               if (!Msmall[iPt]) {
               tmp1 =                (   GDenS [iPt         ] * gPTsx[iPt] 
                                      ) / DSDMnorm[iPt]; 
               tmp2 =                (   gPTsx[iPt] * Hx[iPt] 
                                       + gPTsy[iPt] * Hy[iPt] 
                                       + gPTsz[iPt] * Hz[iPt]
                                       ) * (   GDenS [iPt         ] * Hx[iPt]
                                                                                                ) / DSDMnorm[iPt]; 
               //tmp1 = 0.0;
               //tmp2 = 0.0;
 
                FgX = 
                      ZgammaVar1[iPt] *    GTS   [iPt         ] * Hx[iPt]
                    + ZgammaVar2[iPt] * (  GTX   [iPt         ] )
                    + ZgammaVar3[iPt] *    GDenS [iPt         ] * Hx[iPt]
                    + ZgammaVar4[iPt] * (  GDenX [iPt         ] )
                    + ZgammaVar1[iPt] * (tmp1 - tmp2) *signMD[iPt] ; 
              } else {
                FgX = ZgammaVar1[iPt] *    GTS   [iPt         ] 
                      + ZgammaVar2[iPt] * (  GTX   [iPt         ] )
                      + ZgammaVar3[iPt] *    GDenS [iPt         ] 
                    + ZgammaVar4[iPt] * (  GDenX [iPt         ] ); 
              }

              //Fy
               tmp1 = dcomplex(0.0);
               tmp2 = dcomplex(0.0);
               if (!Msmall[iPt]) {
                tmp1 =                (   GDenS [iPt + NPts      ] * gPTsx[iPt] 
                                      ) / DSDMnorm[iPt]; 
                tmp2 =                (  gPTsx[iPt] * Hx[iPt] 
                                       + gPTsy[iPt] * Hy[iPt] 
                                       + gPTsz[iPt] * Hz[iPt] 
                                       ) * (   GDenS [iPt + NPts   ] * Hx[iPt]
                                                                                                ) / DSDMnorm[iPt]; 
                //tmp1 = 0.0;
                //tmp2 = 0.0;

                 FgY = 
                       ZgammaVar1[iPt] *    GTS   [iPt + NPts      ] * Hx[iPt]
                     + ZgammaVar2[iPt] * (  GTX   [iPt + NPts      ] ) 
                     + ZgammaVar3[iPt] *    GDenS [iPt + NPts      ] * Hx[iPt]
                     + ZgammaVar4[iPt] * (  GDenX [iPt + NPts      ] )
                     + ZgammaVar1[iPt] * (tmp1 - tmp2) *signMD[iPt]; 
               } else {
                 FgY = ZgammaVar1[iPt] *    GTS   [iPt + NPts      ]
                       + ZgammaVar2[iPt] * (  GTX   [iPt + NPts      ] )  
                       + ZgammaVar3[iPt] *    GDenS [iPt + NPts        ] 
                      + ZgammaVar4[iPt] * (  GDenX [iPt + NPts      ] ) ;
               }

              //Fz
               tmp1 = dcomplex(0.0);
               tmp2 = dcomplex(0.0);
               if (!Msmall[iPt]) {
                tmp1 =                (   GDenS [iPt + 2*NPts      ] * gPTsx[iPt] 
                                      ) / DSDMnorm[iPt]; 
                tmp2 =                (  gPTsx[iPt] * Hx[iPt] 
                                       + gPTsy[iPt] * Hy[iPt] 
                                       + gPTsz[iPt] * Hz[iPt]
                                        ) * (   GDenS [iPt + 2*NPts       ] * Hx[iPt]
                                                                                                ) / DSDMnorm[iPt]; 
              //tmp1 = 0.0;
              //tmp2 = 0.0;
  
              FgZ = 
                    ZgammaVar1[iPt] *    GTS   [iPt + 2*NPts      ] * Hx[iPt] 
                  + ZgammaVar2[iPt] * (  GTX   [iPt + 2*NPts      ] )
                  + ZgammaVar3[iPt] *    GDenS [iPt + 2*NPts      ] * Hx[iPt]
                  + ZgammaVar4[iPt] * (  GDenX [iPt + 2*NPts      ] )
                  + ZgammaVar1[iPt] * (tmp1 - tmp2) *signMD[iPt];
               } else {
                 FgZ = ZgammaVar1[iPt] *    GTS   [iPt + 2*NPts      ]
                     + ZgammaVar2[iPt] * (  GTX   [iPt + 2*NPts      ] )
                     + ZgammaVar3[iPt] *    GDenS [iPt + 2*NPts        ]
                     + ZgammaVar4[iPt] * (  GDenX [iPt + 2*NPts      ] ) ;
               }

            } else if (denType == MY ) {
            //} else if (!Msmall[iPt] ) {
              //Fx
               if (!Msmall[iPt]) {
               tmp1 =                (   GDenS [iPt         ] * gPTsy[iPt] 
                                          ) / DSDMnorm[iPt]; 
               tmp2 =                (  gPTsx[iPt] * Hx[iPt] 
                                       + gPTsy[iPt] * Hy[iPt] 
                                       + gPTsz[iPt] * Hz[iPt] 
                                       ) * (   GDenS [iPt         ] * Hy[iPt]
                                                                                                ) / DSDMnorm[iPt]; 
               //tmp1 = 0.0 ; 
               //tmp2 = 0.0;

              FgX = 
                    ZgammaVar1[iPt] *    GTS   [iPt         ] * Hy[iPt]
                  + ZgammaVar2[iPt] * (  GTY   [iPt         ] )
                  + ZgammaVar3[iPt] *    GDenS [iPt         ] * Hy[iPt]
                  + ZgammaVar4[iPt] * (  GDenY [iPt         ] )
                  + ZgammaVar1[iPt] * (tmp1 - tmp2) *signMD[iPt]; 
               } else {
                 FgX = ZgammaVar1[iPt] *    GTS   [iPt         ]  
                   + ZgammaVar2[iPt] * (  GTY   [iPt         ] ) 
                   + ZgammaVar3[iPt] *    GDenS [iPt         ]    
                   + ZgammaVar4[iPt] * (  GDenY [iPt         ] );
               }
              //Fy
               tmp1 = dcomplex(0.0);
               tmp2 = dcomplex(0.0);
               if (!Msmall[iPt]) {
                tmp1 =                (   GDenS [iPt + NPts      ] * gPTsy[iPt] 
                                      ) / DSDMnorm[iPt]; 
                tmp2 =                (  gPTsx[iPt] * Hx[iPt] 
                                       + gPTsy[iPt] * Hy[iPt] 
                                       + gPTsz[iPt] * Hz[iPt] 
                                       ) * (   GDenS [iPt + NPts       ] * Hy[iPt]
                                                                                                ) / DSDMnorm[iPt]; 
                //tmp1 = 0.0;
               // tmp2 = 0.0;

              FgY = 
                    ZgammaVar1[iPt] *    GTS   [iPt + NPts      ] * Hy[iPt]
                  + ZgammaVar2[iPt] * (  GTY   [iPt + NPts      ] )
                  + ZgammaVar3[iPt] *    GDenS [iPt + NPts      ] * Hy[iPt]
                  + ZgammaVar4[iPt] * (  GDenY [iPt + NPts      ] )
                  + ZgammaVar1[iPt] * (tmp1 - tmp2) *signMD[iPt];
               } else {
                 FgY = ZgammaVar1[iPt] *    GTS   [iPt + NPts      ]   
                   + ZgammaVar2[iPt] * (  GTY   [iPt + NPts      ] )  
                   + ZgammaVar3[iPt] *    GDenS [iPt + NPts      ]   
                   + ZgammaVar4[iPt] * (  GDenY [iPt + NPts      ] ) ;
               } 
              //Fz
               tmp1 = dcomplex(0.0);
               tmp2 = dcomplex(0.0);
               if (!Msmall[iPt]) {
                tmp1 =                (   GDenS [iPt + 2*NPts      ] * gPTsy[iPt] 
                                      ) / DSDMnorm[iPt]; 
                tmp2 =                (  gPTsx[iPt] * Hx[iPt] 
                                       + gPTsy[iPt] * Hy[iPt] 
                                       + gPTsz[iPt] * Hz[iPt] 
                                       ) * (   GDenS [iPt + 2*NPts      ] * Hy[iPt]
                                                                                                ) / DSDMnorm[iPt]; 
               //tmp1 = 0.0;
               //tmp2 = 0.0;

              FgZ = 
                    ZgammaVar1[iPt] *    GTS   [iPt + 2*NPts      ] * Hy[iPt] 
                  + ZgammaVar2[iPt] * (  GTY   [iPt + 2*NPts      ] )
                  + ZgammaVar3[iPt] *    GDenS [iPt + 2*NPts      ] * Hy[iPt]
                  + ZgammaVar4[iPt] * (  GDenY [iPt + 2*NPts      ] )
                  + ZgammaVar1[iPt] * (tmp1 - tmp2) *signMD[iPt];
                } else {
                  FgZ = ZgammaVar1[iPt] *    GTS   [iPt + 2*NPts      ]   
                    + ZgammaVar2[iPt] * (  GTY   [iPt + 2*NPts      ] )
                    + ZgammaVar3[iPt] *    GDenS [iPt + 2*NPts      ]   
                    + ZgammaVar4[iPt] * (  GDenY [iPt + 2*NPts      ] ) ;
                }

            } else if (denType == MZ ) {
            //} else if (!Msmall[iPt] ) {
              //Fx
               if (!Msmall[iPt]) {
               tmp1 =                (   GDenS [iPt         ] * gPTsz[iPt] 
                                      ) / DSDMnorm[iPt]; 
               tmp2 =                (  gPTsx[iPt] * Hx[iPt] 
                                       + gPTsy[iPt] * Hy[iPt] 
                                       + gPTsz[iPt] * Hz[iPt] ) * (   GDenS [iPt         ] * Hz[iPt]
                                                                                                ) / DSDMnorm[iPt]; 
               //tmp1 = 0.0;
               //tmp2 = 0.0;

              FgX = 
                    ZgammaVar1[iPt] *    GTS   [iPt         ] * Hz[iPt]
                  + ZgammaVar2[iPt] * (  GTZ   [iPt         ] )
                  + ZgammaVar3[iPt] *    GDenS [iPt         ] * Hz[iPt]
                  + ZgammaVar4[iPt] * (  GDenZ [iPt         ] )
                  + ZgammaVar1[iPt] * (tmp1 - tmp2) *signMD[iPt];
               } else {
                FgX = ZgammaVar1[iPt] *    GTS   [iPt         ]   
                  + ZgammaVar2[iPt] * (  GTZ   [iPt         ] )
                  + ZgammaVar3[iPt] *    GDenS [iPt      ]
                  + ZgammaVar4[iPt] * (  GDenZ [iPt         ] );
               }

              //Fy
               tmp1 = dcomplex(0.0);
               tmp2 = dcomplex(0.0);
               if (!Msmall[iPt]) {
                tmp1 =                (   GDenS [iPt + NPts      ] * gPTsz[iPt] 
                                      ) / DSDMnorm[iPt]; 
                tmp2 =                (  gPTsx[iPt] * Hx[iPt] 
                                       + gPTsy[iPt] * Hy[iPt] 
                                       + gPTsz[iPt] * Hz[iPt] ) * (   GDenS [iPt + NPts      ] * Hz[iPt]
                                                                                                ) / DSDMnorm[iPt]; 
               // tmp1 = 0.0;
               // tmp2 = 0.0;

              FgY = 
                    ZgammaVar1[iPt] *    GTS   [iPt + NPts      ] * Hz[iPt]
                  + ZgammaVar2[iPt] * (  GTZ   [iPt + NPts      ] )
                  + ZgammaVar3[iPt] *    GDenS [iPt + NPts      ] * Hz[iPt]
                  + ZgammaVar4[iPt] * (  GDenZ [iPt + NPts      ] )
                  + ZgammaVar1[iPt] * (tmp1 - tmp2) *signMD[iPt];

               } else {
                 FgY = ZgammaVar1[iPt] *    GTS   [iPt + NPts      ]  
                    + ZgammaVar2[iPt] * (  GTZ   [iPt + NPts      ] )
                    + ZgammaVar3[iPt] *    GDenS [iPt + NPts      ]
                    + ZgammaVar4[iPt] * (  GDenZ [iPt + NPts      ] ); 
               } 
              //Fz
               tmp1 = dcomplex(0.0);
               tmp2 = dcomplex(0.0);
               if (!Msmall[iPt]) {
                tmp1 =                (   GDenS [iPt + 2*NPts      ] * gPTsz[iPt] 
                                      ) / DSDMnorm[iPt]; 
                tmp2 =                (  gPTsx[iPt] * Hx[iPt] 
                                       + gPTsy[iPt] * Hy[iPt] 
                                       + gPTsz[iPt] * Hz[iPt] ) * (   GDenS [iPt + 2*NPts    ] * Hz[iPt]
                                                                                                ) / DSDMnorm[iPt]; 
//                tmp1 = 0.0;
//                tmp2 =0.0;

              FgZ = 
                    ZgammaVar1[iPt] *    GTS   [iPt + 2*NPts      ] * Hz[iPt]
                  + ZgammaVar2[iPt] * (  GTZ   [iPt + 2*NPts      ] )
                  + ZgammaVar3[iPt] *    GDenS [iPt + 2*NPts      ] * Hz[iPt]
                  + ZgammaVar4[iPt] * (  GDenZ [iPt + 2*NPts      ] )
                  + ZgammaVar1[iPt] * (tmp1 - tmp2) *signMD[iPt];
                } else {
                  FgZ = ZgammaVar1[iPt] *    GTS   [iPt + 2*NPts      ]
                     + ZgammaVar2[iPt] * (  GTZ   [iPt + 2*NPts      ] )
                     + ZgammaVar3[iPt] *    GDenS [iPt + 2*NPts      ]
                     + ZgammaVar4[iPt] * (  GDenZ [iPt + 2*NPts      ] )  ; 
                }  
            }

            FgX *= weights[iPt];
            FgY *= weights[iPt];
            FgZ *= weights[iPt];
    
           // double FgX_r = std::real(FgX);
           // double FgY_r = std::real(FgY);
           // double FgZ_r = std::real(FgZ);
    
            blas::axpy(NBE,FgX,BasisScr + iPt*NBE +   IOff,1,ZMAT + iPt*NBE,1);
            blas::axpy(NBE,FgY,BasisScr + iPt*NBE + 2*IOff,1,ZMAT + iPt*NBE,1);
            blas::axpy(NBE,FgZ,BasisScr + iPt*NBE + 3*IOff,1,ZMAT + iPt*NBE,1);
    
            // double FgX_i = std::imag(FgX);
            // double FgY_i = std::imag(FgY);
            // double FgZ_i = std::imag(FgZ);

          } //NPTS 

        } //2C GGA

    }  // isGGA


// SS debug end
  } //End
//2c e 
//SS:end

  // Calculate gPTss,sx,sy,sz 
  template <typename MatsT, typename IntsT>
  template <typename U>
  void KohnSham<MatsT,IntsT>::mkgPTVar( 
    size_t NPts, 
    double* GDenS, double* GDenZ, double* GDenY, double* GDenX, 
    U* GTS, U* GTZ, U* GTY, U* GTX,
    U* gPTss, U* gPTsz, U* gPTsy, U* gPTsx, U* gPTzz, 
    U* gPTyy, U* gPTxx  
    ) {
    
       for(auto iPt = 0ul; iPt < NPts; iPt++) {
         gPTss[iPt]  = (GDenS[iPt]          * GTS[iPt]         ) +
                   (GDenS[iPt + NPts]   * GTS[iPt + NPts]  ) +
                   (GDenS[iPt + 2*NPts] * GTS[iPt + 2*NPts]);

         gPTsz[iPt] = (GDenS[iPt]          * GTZ[iPt]         ) +
                 (GDenS[iPt + NPts]   * GTZ[iPt + NPts]  ) +
                 (GDenS[iPt + 2*NPts] * GTZ[iPt + 2*NPts]);
   
         if( this->onePDM->hasZ() ){

           //since gPTsz + gPTzs are always added 
           gPTsz[iPt] += (GDenZ[iPt]          * GTS[iPt]         ) +
                    (GDenZ[iPt + NPts]   * GTS[iPt + NPts]  ) +
                    (GDenZ[iPt + 2*NPts] * GTS[iPt + 2*NPts]);
   
           gPTzz[iPt] = (GDenZ[iPt]          * GTZ[iPt]         ) +
                   (GDenZ[iPt + NPts]   * GTZ[iPt + NPts]  ) +
                   (GDenZ[iPt + 2*NPts] * GTZ[iPt + 2*NPts]);

         }
         if( this->onePDM->hasXY() ) {
    
           gPTsy[iPt] = (GDenS[iPt]          * GTY[iPt]         ) +
                   (GDenS[iPt + NPts]   * GTY[iPt + NPts]  ) +
                   (GDenS[iPt + 2*NPts] * GTY[iPt + 2*NPts]);
           // Since gPTsy + gPTys are always added 
           gPTsy[iPt] += (GDenY[iPt]          * GTS[iPt]         ) +
                    (GDenY[iPt + NPts]   * GTS[iPt + NPts]  ) +
                    (GDenY[iPt + 2*NPts] * GTS[iPt + 2*NPts]);
   
           gPTyy[iPt] = (GDenY[iPt]          * GTY[iPt]         ) +
                   (GDenY[iPt + NPts]   * GTY[iPt + NPts]  ) +
                   (GDenY[iPt + 2*NPts] * GTY[iPt + 2*NPts]);
    
           gPTsx[iPt] = (GDenS[iPt]          * GTX[iPt]         ) +
                   (GDenS[iPt + NPts]   * GTX[iPt + NPts]  ) +
                   (GDenS[iPt + 2*NPts] * GTX[iPt + 2*NPts]);
           // Since gPTsx + gPTxs are always added 
           gPTsx[iPt] += (GDenX[iPt]          * GTS[iPt]         ) +
                    (GDenX[iPt + NPts]   * GTS[iPt + NPts]  ) +
                    (GDenX[iPt + 2*NPts] * GTS[iPt + 2*NPts]);
   
           gPTxx[iPt] = (GDenX[iPt]          * GTX[iPt]         ) +
                   (GDenX[iPt + NPts]   * GTX[iPt + NPts]  ) +
                   (GDenX[iPt + 2*NPts] * GTX[iPt + 2*NPts]);
         }
    
    }; //for(auto iPt = 0ul; iPt < NPts; iPt++)

  } // void mkgPTVar






  template <typename MatsT, typename IntsT>
  template <typename U>
  void KohnSham<MatsT, IntsT>::formFXC( MPI_Comm c,  
      std::vector<TwoBodyContraction<U>> &cList, EMPerturbation &pert) {


    size_t itOff = this->nC == 2 ? 5 : 3;
    size_t nVec = cList.size() / itOff;
    size_t NB     = this->basisSet().nBasis;
    size_t NB2    = NB*NB;
    size_t NPPB   = intParam.nRadPerBatch * intParam.nAng;
    size_t nAtoms = this->molecule().nAtoms;


    // Parallelism
    size_t NT = GetNumThreads();
    size_t LAThreads = GetLAThreads();
    size_t mpiRank = MPIRank(c);
    size_t mpiSize = MPISize(c);

    // Turn off LA threads
    SetLAThreads(1);

    // Split the MPI Comm
    int color = ((mpiSize < nAtoms) or 
                 (mpiRank < nAtoms)) ? 1 : MPI_UNDEFINED;
                                                            
                  
    MPI_Comm intComm = MPICommSplit(c,color,mpiRank);


    bool isGGA = std::any_of(functionals.begin(),functionals.end(),
                   [](std::shared_ptr<DFTFunctional> &x) {
                     return x->isGGA(); 
                   }); 

#ifdef CQ_ENABLE_MPI
    if( intComm != MPI_COMM_NULL ) {
#endif

    U* NBNBSCR = CQMemManager::get().malloc<U>(NB2 * NT);
    U* NBNPSCR = CQMemManager::get().malloc<U>(NB*NPPB * NT);

    double* NBNBSCRD = CQMemManager::get().malloc<double>(NB2 * NT);
    double* NBNPSCRD = CQMemManager::get().malloc<double>(NB*NPPB * NT);


    std::shared_ptr<cqmatrix::PauliSpinorMatrices<double>> Re1PDM
        = std::make_shared<cqmatrix::PauliSpinorMatrices<double>>(
            this->onePDM->real_part());


    std::vector<std::vector<U*>> ReTSymm;
    for(auto iVec = 0; iVec < nVec; iVec++) {
      ReTSymm.emplace_back();

      size_t indx = iVec * itOff;
      for(auto iS = 0; iS < 2*this->nC; iS++) {
        
         ReTSymm.back().emplace_back(
             CQMemManager::get().malloc<U>(NB2));
         SetMat('N',NB,NB,U(1.0),cList[indx + iS + 1].X,NB,
           ReTSymm.back().back(),NB);

      }

    }


    U* GxcT_raw = CQMemManager::get().malloc<U>(2*this->nC*nVec*NT*NB2);
    U* Gxc_first = GxcT_raw;
    std::vector<std::vector<std::vector<U*>>> GxcT;
    for(auto ithread = 0; ithread < NT; ithread++) {
      GxcT.emplace_back();
      for(auto iVec = 0; iVec < nVec; iVec++) { 
        GxcT.back().emplace_back();
        for(auto iS = 0; iS < 2*this->nC; iS++){
          GxcT.back().back().emplace_back(Gxc_first);
          memset(GxcT.back().back().back(),0,NB2*sizeof(U));
          Gxc_first += NB2;
        }
      }
    }






    // Allocation of V Vairables
    double *DenS(nullptr),  *DenZ(nullptr),  *DenY(nullptr),  *DenX(nullptr);
    double *GDenS(nullptr), *GDenZ(nullptr), *GDenY(nullptr), *GDenX(nullptr);

    // Allocation of auxiliary variables used in small m limit
    // for 2c TDKS
    double *Mnorm(nullptr);  // |m| 
    // rho^K/|rho^k| for large m and 1/3 for small m
    // Eq. D.21a Shichao Sun's thesis 
    double *KScratch(nullptr);   
    // Eq. D.21b Shichao Sun's thesis 
    double *HScratch(nullptr);   
    // Checks if m is below a threshold at that grid point 
    bool   *Msmall(nullptr);
    // |del rho * del m| 
    double *DSDMnorm(nullptr);
    // Eq. D.6 in Shichao Sun's thesis
    double *signMD(nullptr);

    // Density
    DenS = CQMemManager::get().malloc<double>(NPPB * NT);
    if( this->nC == 2 or not this->iCS )
      DenZ = CQMemManager::get().malloc<double>(NPPB * NT);
    if( this->nC == 2 ) {
      DenY = CQMemManager::get().malloc<double>(NPPB * NT);
      DenX = CQMemManager::get().malloc<double>(NPPB * NT);
      Mnorm    = CQMemManager::get().malloc<double>(NPPB * NT);
      KScratch = CQMemManager::get().malloc<double>(3 * NPPB * NT); // 3 is for K=x,y,z
      Msmall   = CQMemManager::get().malloc<bool>(NPPB * NT);
    }


    // Density Gradient
    if( isGGA ) {

      GDenS = CQMemManager::get().malloc<double>(3*NPPB * NT);
      if( this->nC == 2 or not this->iCS )
        GDenZ = CQMemManager::get().malloc<double>(3*NPPB * NT);
      if( this->nC == 2 ) {
        GDenY = CQMemManager::get().malloc<double>(3*NPPB * NT);
        GDenX = CQMemManager::get().malloc<double>(3*NPPB * NT);
        HScratch = CQMemManager::get().malloc<double>(3*NPPB * NT);
        DSDMnorm    = CQMemManager::get().malloc<double>(NPPB * NT);
        signMD      = CQMemManager::get().malloc<double>(NPPB * NT);
      }

    }

    // Allocation of T evaluation
    U *TS(nullptr),  *TZ(nullptr),  *TY(nullptr),  *TX(nullptr);
    U *GTS(nullptr), *GTZ(nullptr), *GTY(nullptr), *GTX(nullptr);

    // T
    TS = CQMemManager::get().malloc<U>(NPPB * NT);
    TZ = CQMemManager::get().malloc<U>(NPPB * NT);
    if( this->nC == 2 ) {
      TY = CQMemManager::get().malloc<U>(NPPB * NT);
      TX = CQMemManager::get().malloc<U>(NPPB * NT);
    }

    // Allocation of gPT variables 
    // Products of density gradient with transition-density gradient
    U *gPTss(nullptr), *gPTsz(nullptr), *gPTsy(nullptr), *gPTsx(nullptr), 
     *gPTzz(nullptr), *gPTyy(nullptr), *gPTxx(nullptr); 

    // T Gradient
    if( isGGA ) {

      GTS = CQMemManager::get().malloc<U>(3*NPPB * NT);
      gPTss = CQMemManager::get().malloc<U>(NPPB * NT);
      GTZ = CQMemManager::get().malloc<U>(3*NPPB * NT);
      gPTsz = CQMemManager::get().malloc<U>(NPPB * NT);

      if( this->onePDM->hasZ() ) 
        gPTzz = CQMemManager::get().malloc<U>(NPPB * NT);

      if( this->onePDM->hasXY() ) {
        GTY = CQMemManager::get().malloc<U>(3*NPPB * NT);
        GTX = CQMemManager::get().malloc<U>(3*NPPB * NT);


        gPTsx = CQMemManager::get().malloc<U>(NPPB * NT);
        gPTsy = CQMemManager::get().malloc<U>(NPPB * NT);
        gPTyy = CQMemManager::get().malloc<U>(NPPB * NT);
        gPTxx = CQMemManager::get().malloc<U>(NPPB * NT);
      }
    }





    // U Variables

    double * eps     = CQMemManager::get().malloc<double>(NPPB * NT);
    double * U_n     = CQMemManager::get().malloc<double>(2*NPPB * NT);
    double * U_gamma = isGGA ? 
      CQMemManager::get().malloc<double>(3*NPPB * NT) : nullptr;

    // U First Derivatives

    double * dVU_n     = CQMemManager::get().malloc<double>(2*NPPB * NT);
    double * dVU_gamma = isGGA ? 
      CQMemManager::get().malloc<double>(3*NPPB * NT) : nullptr;

    double * eps_SCR(nullptr), * dVU_n_SCR(nullptr), * dVU_gamma_SCR(nullptr);
    if( functionals.size() > 1 ) {

      eps_SCR       = CQMemManager::get().malloc<double>(NPPB * NT);
      dVU_n_SCR     = CQMemManager::get().malloc<double>(2*NPPB * NT);
      dVU_gamma_SCR = isGGA ? 
        CQMemManager::get().malloc<double>(3*NPPB * NT) : nullptr;

    }

    // U Second Derivatives

    double * d2VU_n = CQMemManager::get().malloc<double>(3*NPPB * NT);
    double * d2VU_gamma   = isGGA ? 
      CQMemManager::get().malloc<double>(6*NPPB * NT) : nullptr;
    double * d2VU_n_gamma = isGGA ? 
      CQMemManager::get().malloc<double>(6*NPPB * NT) : nullptr;

    double * d2VU_n_SCR(nullptr), * d2VU_gamma_SCR(nullptr), 
           * d2VU_n_gamma_SCR(nullptr); 

    if( functionals.size() > 1 ) {
      d2VU_n_SCR = CQMemManager::get().malloc<double>(3*NPPB * NT);
      d2VU_gamma_SCR   = isGGA ? 
        CQMemManager::get().malloc<double>(6*NPPB * NT) : nullptr;
      d2VU_n_gamma_SCR = isGGA ? 
        CQMemManager::get().malloc<double>(6*NPPB * NT) : nullptr;
    }


    // Z Vars
    U * ZrhoVar = CQMemManager::get().malloc<U>(NPPB * NT);

    U * ZgammaVar1 = isGGA ? 
      CQMemManager::get().malloc<U>(NPPB * NT) : nullptr;
    U * ZgammaVar2 = isGGA ? 
      CQMemManager::get().malloc<U>(NPPB * NT) : nullptr;
    U * ZgammaVar3 = isGGA ? 
      CQMemManager::get().malloc<U>(NPPB * NT) : nullptr;
    U * ZgammaVar4 = isGGA ? 
      CQMemManager::get().malloc<U>(NPPB * NT) : nullptr;

    U* ZMAT = CQMemManager::get().malloc<U>(NB*NPPB * NT);

    int NDer = isGGA ? 4 : 1;
    U* Basis_cmplx(nullptr);
    if( std::is_same<U,dcomplex>::value ) 
      Basis_cmplx = CQMemManager::get().malloc<U>(NB*NDer*NPPB*NT);

    double intDen = 0.;


    auto fxcbuild = [&](size_t &res, std::vector<cart_t> &batch, 
      std::vector<double> &weights, std::vector<size_t> NBE_vec, 
      std::vector<double*> BasisEval_vec, 
      std::vector<std::vector<size_t>>& batchEvalShells_vec, 
      std::vector<std::vector<std::pair<size_t,size_t>>>& subMatCut_vec) {


      size_t NBE = NBE_vec[0];
      double * BasisEval = BasisEval_vec[0];
      std::vector<size_t> & batchEvalShells = batchEvalShells_vec[0];
      std::vector<std::pair<size_t,size_t>> & subMatCut = subMatCut_vec[0];

      double epsScreen = intParam.epsilon / nAtoms /
        intParam.nAng / intParam.nRad;

      epsScreen = std::max(epsScreen,
                           std::numeric_limits<double>::epsilon());



      size_t NPts = batch.size();
      size_t IOff = NBE*NPts;
      int tid = GetThreadID();
      size_t nPtOff = tid * NPPB;


      // Get thread Local storage
      

      // Density and T
      double *DenS_loc = DenS + nPtOff;
      double *DenZ_loc = DenZ + nPtOff;
      double *DenY_loc = DenY + nPtOff;
      double *DenX_loc = DenX + nPtOff;

      double *GDenS_loc = GDenS + 3 * nPtOff;
      double *GDenZ_loc = GDenZ + 3 * nPtOff;
      double *GDenY_loc = GDenY + 3 * nPtOff;
      double *GDenX_loc = GDenX + 3 * nPtOff;

      // Aux Var for 2C 
      double * Mnorm_loc    = Mnorm        +   nPtOff;
      double * KScratch_loc = KScratch     + 3*nPtOff;
      bool   * Msmall_loc   = Msmall       +   nPtOff;
      double * HScratch_loc = HScratch     + 3*nPtOff;
      double * DSDMnorm_loc = DSDMnorm     +   nPtOff;
      double * signMD_loc   = signMD       +   nPtOff;

      U *TS_loc = TS + nPtOff;
      U *TZ_loc = TZ + nPtOff;
      U *TY_loc = TY + nPtOff;
      U *TX_loc = TX + nPtOff;

      U *GTS_loc = GTS + 3 * nPtOff;
      U *GTZ_loc = GTZ + 3 * nPtOff;
      U *GTY_loc = GTY + 3 * nPtOff;
      U *GTX_loc = GTX + 3 * nPtOff;

      U *gPTss_loc = gPTss + nPtOff;
      U *gPTsz_loc = gPTsz + nPtOff;
      U *gPTsy_loc = gPTsy + nPtOff;
      U *gPTsx_loc = gPTsx + nPtOff;
      U *gPTxx_loc = gPTxx + nPtOff;
      U *gPTyy_loc = gPTyy + nPtOff;
      U *gPTzz_loc = gPTzz + nPtOff;

      // U Vars

      double *U_n_loc        = U_n        + 2 * nPtOff;     
      double *U_gamma_loc    = U_gamma    + 3 * nPtOff;

      double *eps_loc        = eps        +     nPtOff;
      double *dVU_n_loc      = dVU_n      + 2 * nPtOff;     
      double *dVU_gamma_loc  = dVU_gamma  + 3 * nPtOff;

      double *d2VU_n_loc       = d2VU_n       + 3 * nPtOff;     
      double *d2VU_gamma_loc   = d2VU_gamma   + 6 * nPtOff;
      double *d2VU_n_gamma_loc = d2VU_n_gamma + 6 * nPtOff;

      double *eps_SCR_loc       = eps_SCR       +     nPtOff;
      double *dVU_n_SCR_loc     = dVU_n_SCR     + 2 * nPtOff;     
      double *dVU_gamma_SCR_loc = dVU_gamma_SCR + 3 * nPtOff;

      double *d2VU_n_SCR_loc       = d2VU_n_SCR       + 3 * nPtOff;     
      double *d2VU_gamma_SCR_loc   = d2VU_gamma_SCR   + 6 * nPtOff;
      double *d2VU_n_gamma_SCR_loc = d2VU_n_gamma_SCR + 6 * nPtOff;


      // Z Vars

      U* ZrhoVar_loc    = ZrhoVar    + nPtOff;
      U* ZgammaVar1_loc = ZgammaVar1 + nPtOff;
      U* ZgammaVar2_loc = ZgammaVar2 + nPtOff;
      U* ZgammaVar3_loc = ZgammaVar3 + nPtOff;
      U* ZgammaVar4_loc = ZgammaVar4 + nPtOff;

      U* ZMAT_loc = ZMAT + nPtOff*NB;


      U* NBNBSCR_loc = NBNBSCR + NB * NB   * tid;
      U* NBNPSCR_loc = NBNPSCR + NB * NPPB * tid;

      double* NBNBSCRD_loc = NBNBSCRD + NB * NB   * tid;
      double* NBNPSCRD_loc = NBNPSCRD + NB * NPPB * tid;

      U* NBNBSCR_r = NBNBSCR_loc;
      U* NBNPSCR_r = NBNPSCR_loc;

      double* NBNBSCRD_r = NBNBSCRD_loc;
      double* NBNPSCRD_r = NBNPSCRD_loc;

      // Basis
      
      U* Basis_use = reinterpret_cast<U*>(BasisEval);
      if( std::is_same<U,dcomplex>::value ) {
        Basis_use = Basis_cmplx + nPtOff*NB*NDer;
        std::copy_n(BasisEval, NDer * NPts * NBE, Basis_use);
      } 

      // This evaluates the V variables for all components
      // (Scalar, MZ (UKS) and Mx, MY (2 Comp))
      evalDen((isGGA ? GRADIENT : NOGRAD), NPts, NBE, NB, subMatCut, 
        NBNBSCRD_r, NBNPSCRD_r, Re1PDM->S().pointer(), DenS_loc , 
        GDenS_loc, GDenS_loc + NPts, GDenS_loc + 2*NPts, BasisEval);

      if( this->onePDM->hasZ() )
        evalDen((isGGA ? GRADIENT : NOGRAD), NPts, NBE, NB, subMatCut, 
          NBNBSCRD_r ,NBNPSCRD_r, Re1PDM->Z().pointer(), DenZ_loc, 
          GDenZ_loc, GDenZ_loc + NPts, GDenZ_loc + 2*NPts, BasisEval);

      if( this->onePDM->hasXY() ) {
        evalDen((isGGA ? GRADIENT : NOGRAD), NPts, NBE, NB, subMatCut, 
          NBNBSCRD_r ,NBNPSCRD_r, Re1PDM->Y().pointer(), DenY_loc, 
          GDenY_loc, GDenY_loc + NPts, GDenY_loc + 2*NPts, BasisEval);
        evalDen((isGGA ? GRADIENT : NOGRAD), NPts, NBE, NB, subMatCut, 
          NBNBSCRD_r ,NBNPSCRD_r, Re1PDM->X().pointer(), DenX_loc, 
          GDenX_loc, GDenX_loc + NPts, GDenX_loc + 2*NPts, BasisEval);
      
      }

      // V -> U variables for evaluating the kernel derivatives.
      mkAuxVar(this->onePDM,isGGA,epsScreen,NPts,
        DenS_loc,DenZ_loc,DenY_loc,DenX_loc,
        GDenS_loc,GDenS_loc + NPts,GDenS_loc + 2*NPts,
        GDenZ_loc,GDenZ_loc + NPts,GDenZ_loc + 2*NPts,
        GDenY_loc,GDenY_loc + NPts,GDenY_loc + 2*NPts,
        GDenX_loc,GDenX_loc + NPts,GDenX_loc + 2*NPts,
        Mnorm_loc, 
        KScratch_loc, KScratch_loc + NPts, KScratch_loc + 2* NPts,
        HScratch_loc, HScratch_loc + NPts, HScratch_loc + 2* NPts,
        DSDMnorm_loc, signMD_loc, 
        Msmall_loc,U_n_loc,U_gamma_loc
      );



      loadFXCder(NPts,U_n_loc,U_gamma_loc,eps_loc,dVU_n_loc,d2VU_n_loc,
        dVU_gamma_loc, d2VU_gamma_loc,d2VU_n_gamma_loc,eps_SCR_loc,
        dVU_n_SCR_loc,dVU_gamma_SCR_loc, d2VU_n_SCR_loc,d2VU_gamma_SCR_loc,
        d2VU_n_gamma_SCR_loc);


      for(auto iT = 0; iT < nVec; iT++) {

        size_t indx       = itOff * iT;

        // Transition density build 
        // (Scalar, MZ (UKS) and Mx, MY (2 Comp))
        evalTransDen((isGGA ? GRADIENT : NOGRAD), NPts, NBE, NB, subMatCut, 
          NBNBSCR_r, NBNPSCR_r, ReTSymm[iT][SCALAR] , TS_loc,
          GTS_loc, GTS_loc + NPts, GTS_loc + 2*NPts, Basis_use);

        evalTransDen((isGGA ? GRADIENT : NOGRAD), NPts, NBE, NB, subMatCut, 
          NBNBSCR_r ,NBNPSCR_r, ReTSymm[iT][MZ] , TZ_loc,
          GTZ_loc, GTZ_loc + NPts, GTZ_loc + 2*NPts, Basis_use);

        if( this->onePDM->hasXY() ) {

          evalTransDen((isGGA ? GRADIENT : NOGRAD), NPts, NBE, NB, subMatCut, 
            NBNBSCR_r ,NBNPSCR_r, ReTSymm[iT][MY], TY_loc,
            GTY_loc, GTY_loc + NPts, GTY_loc + 2*NPts, Basis_use);

          evalTransDen((isGGA ? GRADIENT : NOGRAD), NPts, NBE, NB, subMatCut, 
            NBNBSCR_r ,NBNPSCR_r, ReTSymm[iT][MX], TX_loc,
            GTX_loc, GTX_loc + NPts, GTX_loc + 2*NPts, Basis_use);

        }

        //TODO: Remove nC constraint and delete gPT build in constructZVars
        if (isGGA and this->nC == 2) {

          mkgPTVar( NPts, 
            GDenS_loc, GDenZ_loc, GDenY_loc, GDenX_loc, 
            GTS_loc, GTZ_loc, GTY_loc, GTX_loc,
            gPTss_loc, gPTsz_loc, gPTsy_loc, gPTsx_loc, gPTzz_loc, gPTyy_loc, gPTxx_loc
            );

        }

        constructZVarsFXC(SCALAR,isGGA,NPts, GDenS_loc, GDenZ_loc, GDenY_loc, GDenX_loc, 
          Msmall_loc,Mnorm_loc,
          KScratch_loc, KScratch_loc + NPts, KScratch_loc + 2* NPts,
          HScratch_loc, HScratch_loc + NPts, HScratch_loc + 2* NPts,
          TS_loc, TZ_loc, TY_loc, TX_loc, GTS_loc, GTZ_loc, GTY_loc, 
          GTX_loc, 
          gPTss_loc, gPTsz_loc, gPTsy_loc, gPTsx_loc, gPTzz_loc, gPTyy_loc, gPTxx_loc, 
          dVU_n_loc, dVU_gamma_loc, d2VU_n_loc, d2VU_gamma_loc, 
          d2VU_n_gamma_loc, ZrhoVar_loc, ZgammaVar1_loc, ZgammaVar2_loc,
          ZgammaVar3_loc, ZgammaVar4_loc);

        formZ_fxc(SCALAR,isGGA,NPts,NBE,IOff,epsScreen,weights,ZrhoVar_loc,
          ZgammaVar1_loc, ZgammaVar2_loc, ZgammaVar3_loc, ZgammaVar4_loc, 
          Msmall_loc,Mnorm_loc,
          DSDMnorm_loc, signMD_loc, 
          GDenS_loc, GDenZ_loc, GDenY_loc, GDenX_loc, 
          KScratch_loc, KScratch_loc + NPts, KScratch_loc + 2* NPts,
          HScratch_loc, HScratch_loc + NPts, HScratch_loc + 2* NPts,
          GTS_loc, GTZ_loc, GTY_loc, GTX_loc, 
          gPTss_loc, gPTsz_loc, gPTsy_loc, gPTsx_loc, gPTzz_loc, gPTyy_loc, gPTxx_loc, 
          BasisEval, ZMAT_loc);

        blas::syr2k(blas::Layout::ColMajor,blas::Uplo::Lower,blas::Op::NoTrans,NBE,NPts,U(0.5),Basis_use,NBE,ZMAT_loc,NBE,U(0.),NBNBSCR_loc,NBE);

        IncBySubMat(NB,NB,NBE,NBE,GxcT[tid][iT][SCALAR],NB,NBNBSCR_loc,NBE,
            subMatCut);


        //std::cerr << "MZ bit "<< std::endl;
        constructZVarsFXC(MZ,isGGA,NPts, GDenS_loc, GDenZ_loc, GDenY_loc, GDenX_loc, 
          Msmall_loc,Mnorm_loc,
          KScratch_loc, KScratch_loc + NPts, KScratch_loc + 2* NPts,
          HScratch_loc, HScratch_loc + NPts, HScratch_loc + 2* NPts,
          TS_loc, TZ_loc, TY_loc, TX_loc, GTS_loc, GTZ_loc, GTY_loc, 
          GTX_loc, 
          gPTss_loc, gPTsz_loc, gPTsy_loc, gPTsx_loc, gPTzz_loc, gPTyy_loc, gPTxx_loc, 
          dVU_n_loc, dVU_gamma_loc, d2VU_n_loc, d2VU_gamma_loc, 
          d2VU_n_gamma_loc, ZrhoVar_loc, ZgammaVar1_loc, ZgammaVar2_loc,
          ZgammaVar3_loc, ZgammaVar4_loc);

        formZ_fxc(MZ,isGGA,NPts,NBE,IOff,epsScreen,weights,ZrhoVar_loc,
          ZgammaVar1_loc, ZgammaVar2_loc, ZgammaVar3_loc, ZgammaVar4_loc,  
          Msmall_loc,Mnorm_loc,
          DSDMnorm_loc, signMD_loc, 
          GDenS_loc, GDenZ_loc, GDenY_loc, GDenX_loc, 
          KScratch_loc, KScratch_loc + NPts, KScratch_loc + 2* NPts,
          HScratch_loc, HScratch_loc + NPts, HScratch_loc + 2* NPts,
          GTS_loc, GTZ_loc, GTY_loc, GTX_loc, 
          gPTss_loc, gPTsz_loc, gPTsy_loc, gPTsx_loc, gPTzz_loc, gPTyy_loc, gPTxx_loc, 
          BasisEval, ZMAT_loc);

        blas::syr2k(blas::Layout::ColMajor,blas::Uplo::Lower,blas::Op::NoTrans,NBE,NPts,U(0.5),Basis_use,NBE,ZMAT_loc,NBE,U(0.),NBNBSCR_loc,NBE);

        IncBySubMat(NB,NB,NBE,NBE,GxcT[tid][iT][MZ],NB,NBNBSCR_loc,NBE,
            subMatCut);


        if( this->onePDM->hasXY() ) {
          // std::cerr << "MX bit "<< std::endl;
          constructZVarsFXC(MX,isGGA,NPts, GDenS_loc, GDenZ_loc, GDenY_loc, GDenX_loc, 
            Msmall_loc,Mnorm_loc,
            KScratch_loc, KScratch_loc + NPts, KScratch_loc + 2* NPts,
            HScratch_loc, HScratch_loc + NPts, HScratch_loc + 2* NPts,
            TS_loc, TZ_loc, TY_loc, TX_loc, GTS_loc, GTZ_loc, GTY_loc, 
            GTX_loc, 
            gPTss_loc, gPTsz_loc, gPTsy_loc, gPTsx_loc, gPTzz_loc, gPTyy_loc, gPTxx_loc, 
            dVU_n_loc, dVU_gamma_loc, d2VU_n_loc, d2VU_gamma_loc, 
            d2VU_n_gamma_loc, ZrhoVar_loc, ZgammaVar1_loc, ZgammaVar2_loc,
            ZgammaVar3_loc, ZgammaVar4_loc);

          formZ_fxc(MX,isGGA,NPts,NBE,IOff,epsScreen,weights,ZrhoVar_loc,
            ZgammaVar1_loc, ZgammaVar2_loc, ZgammaVar3_loc, ZgammaVar4_loc, 
            Msmall_loc,Mnorm_loc,
            DSDMnorm_loc, signMD_loc, 
            GDenS_loc, GDenZ_loc, GDenY_loc, GDenX_loc, 
            KScratch_loc, KScratch_loc + NPts, KScratch_loc + 2* NPts,
            HScratch_loc, HScratch_loc + NPts, HScratch_loc + 2* NPts,
            GTS_loc, GTZ_loc, GTY_loc, GTX_loc, 
            gPTss_loc, gPTsz_loc, gPTsy_loc, gPTsx_loc, gPTzz_loc, gPTyy_loc, gPTxx_loc, 
            BasisEval, ZMAT_loc);

          blas::syr2k(blas::Layout::ColMajor,blas::Uplo::Lower,blas::Op::NoTrans,NBE,NPts,U(0.5),Basis_use,NBE,ZMAT_loc,NBE,U(0.),NBNBSCR_loc,NBE);

          IncBySubMat(NB,NB,NBE,NBE,GxcT[tid][iT][MX],NB,NBNBSCR_loc,NBE,
              subMatCut);

          //std::cerr << "MY bit "<< std::endl;
          constructZVarsFXC(MY,isGGA,NPts, GDenS_loc, GDenZ_loc, GDenY_loc, GDenX_loc, 
            Msmall_loc,Mnorm_loc,
            KScratch_loc, KScratch_loc + NPts, KScratch_loc + 2* NPts,
            HScratch_loc, HScratch_loc + NPts, HScratch_loc + 2* NPts,
            TS_loc, TZ_loc, TY_loc, TX_loc, GTS_loc, GTZ_loc, GTY_loc, 
            GTX_loc, 
            gPTss_loc, gPTsz_loc, gPTsy_loc, gPTsx_loc, gPTzz_loc, gPTyy_loc, gPTxx_loc, 
            dVU_n_loc, dVU_gamma_loc, d2VU_n_loc, d2VU_gamma_loc, 
            d2VU_n_gamma_loc, ZrhoVar_loc, ZgammaVar1_loc, ZgammaVar2_loc,
            ZgammaVar3_loc, ZgammaVar4_loc);

          formZ_fxc(MY,isGGA,NPts,NBE,IOff,epsScreen,weights,ZrhoVar_loc,
            ZgammaVar1_loc, ZgammaVar2_loc, ZgammaVar3_loc, ZgammaVar4_loc,  
            Msmall_loc,Mnorm_loc,
            DSDMnorm_loc, signMD_loc, 
            GDenS_loc, GDenZ_loc, GDenY_loc, GDenX_loc, 
            KScratch_loc, KScratch_loc + NPts, KScratch_loc + 2* NPts,
            HScratch_loc, HScratch_loc + NPts, HScratch_loc +  NPts,
            GTS_loc, GTZ_loc, GTY_loc, GTX_loc, 
            gPTss_loc, gPTsz_loc, gPTsy_loc, gPTsx_loc, gPTzz_loc, gPTyy_loc, gPTxx_loc, 
            BasisEval, ZMAT_loc);

          blas::syr2k(blas::Layout::ColMajor,blas::Uplo::Lower,blas::Op::NoTrans,NBE,NPts,U(0.5),Basis_use,NBE,ZMAT_loc,NBE,U(0.),NBNBSCR_loc,NBE);

          IncBySubMat(NB,NB,NBE,NBE,GxcT[tid][iT][MY],NB,NBNBSCR_loc,NBE,
              subMatCut);

        } // 2C

      } // iT loop

    }; // VXC

    // Create the BeckeIntegrator object
    BeckeIntegrator<EulerMac> 
      integrator(intComm,this->molecule(),
        this->basisSet(), EulerMac(intParam.nRad), intParam.nAng,
        intParam.nRadPerBatch, (isGGA ? GRADIENT : NOGRAD), intParam.epsilon);

    // Integrate the FXC
    integrator.integrate<size_t>(fxcbuild);

    U* mpiScr = nullptr;
#ifdef CQ_ENABLE_MPI
    if( MPIRank(intComm) == 0 and MPISize(intComm) > 1 )
      mpiScr = CQMemManager::get().malloc<U>(NB*NB);
#endif

    for(auto iT = 0; iT < nVec; iT++) {

      for(auto iS = 0; iS < 2*this->nC; iS++) {

        // Add -Gxc[T] thread contributions into single storage
        for(auto ithread = 0; ithread < NT; ithread++) {

          for(auto j = 0; j < NB; j++)
          for(auto i = j; i < NB; i++)
            GxcT[ithread][iT][iS][ j + i*NB ] =
              GxcT[ithread][iT][iS][ i + j*NB ];

          U fact = (ithread == 0)  ? 0 : 1.;
          MatAdd('N','N',NB,NB,
            fact         , GxcT[0][iT][iS]      ,NB,
            U(-4. * M_PI), GxcT[ithread][iT][iS],NB,
            GxcT[0][iT][iS],NB);
        }

#ifdef CQ_ENABLE_MPI
        // Add MPI Contributions together
        if( MPISize(intComm) > 1 )
          MPIReduce(GxcT[0][iT][iS],NB2,mpiScr,0,intComm);

        if( MPIRank(intComm) == 0 and MPISize(intComm) > 1 )
          std::copy_n(mpiScr,NB2,GxcT[0][iT][iS]);
#endif

        // Add GxcT to K
        if( MPIRank(intComm) == 0 )
          MatAdd('N','N',NB,NB,
            U(functionals.back()->xHFX), cList[iT*itOff + iS + 1].AX, NB,
            U(1.),                       GxcT[0][iT][iS],             NB,
            cList[iT*itOff + iS + 1].AX, NB);

      }

    }


    // Free up the memory
    if( mpiScr ) CQMemManager::get().free(mpiScr);
    CQMemManager::get().free( GxcT_raw, NBNBSCR, NBNPSCR, NBNBSCRD, NBNPSCRD );

    for(auto &Y : ReTSymm) for(auto &X : Y) CQMemManager::get().free(X);
    Re1PDM = nullptr;

    if( DenS ) CQMemManager::get().free( DenS );
    if( DenZ ) CQMemManager::get().free( DenZ );
    if( DenY ) CQMemManager::get().free( DenY );
    if( DenX ) CQMemManager::get().free( DenX );

    if( Mnorm ) CQMemManager::get().free( Mnorm ); 
    if( KScratch ) CQMemManager::get().free( KScratch ); 
    if( Msmall ) CQMemManager::get().free( Msmall ); 
    if( HScratch ) CQMemManager::get().free( HScratch ); 
    if( DSDMnorm ) CQMemManager::get().free( DSDMnorm );
    if( signMD ) CQMemManager::get().free( signMD ); 
    if( Basis_cmplx ) CQMemManager::get().free( Basis_cmplx ); 

    if( GDenS ) CQMemManager::get().free( GDenS );
    if( GDenZ ) CQMemManager::get().free( GDenZ );
    if( GDenY ) CQMemManager::get().free( GDenY );
    if( GDenX ) CQMemManager::get().free( GDenX );

    if( TS ) CQMemManager::get().free( TS );
    if( TZ ) CQMemManager::get().free( TZ );
    if( TY ) CQMemManager::get().free( TY );
    if( TX ) CQMemManager::get().free( TX );

    if( GTS ) CQMemManager::get().free( GTS );
    if( GTZ ) CQMemManager::get().free( GTZ );
    if( GTY ) CQMemManager::get().free( GTY );
    if( GTX ) CQMemManager::get().free( GTX );

    if( gPTss ) CQMemManager::get().free( gPTss );
    if( gPTsx ) CQMemManager::get().free( gPTsx );
    if( gPTsy ) CQMemManager::get().free( gPTsy );
    if( gPTsz ) CQMemManager::get().free( gPTsz );
    if( gPTzz ) CQMemManager::get().free( gPTzz );
    if( gPTyy ) CQMemManager::get().free( gPTyy );
    if( gPTxx ) CQMemManager::get().free( gPTxx );

    CQMemManager::get().free( eps, U_n, dVU_n, d2VU_n );
    if( U_gamma )      CQMemManager::get().free( U_gamma );
    if( dVU_gamma )    CQMemManager::get().free( dVU_gamma );
    if( d2VU_gamma )   CQMemManager::get().free( d2VU_gamma );
    if( d2VU_n_gamma ) CQMemManager::get().free( d2VU_n_gamma );


    if( eps_SCR )       CQMemManager::get().free( eps_SCR );
    if( dVU_n_SCR )     CQMemManager::get().free( dVU_n_SCR );
    if( dVU_gamma_SCR ) CQMemManager::get().free( dVU_gamma_SCR );
    if( d2VU_n_SCR )     CQMemManager::get().free( d2VU_n_SCR );
    if( d2VU_gamma_SCR ) CQMemManager::get().free( d2VU_gamma_SCR );
    if( d2VU_n_gamma_SCR ) CQMemManager::get().free( d2VU_n_gamma_SCR );

    CQMemManager::get().free( ZrhoVar, ZMAT );
    if( isGGA ) {
      CQMemManager::get().free( ZgammaVar1, ZgammaVar2, ZgammaVar3, ZgammaVar4 );
    }

#ifdef CQ_ENABLE_MPI
    MPICommFree(intComm); // Free communicator

    } // End of the MPI

    MPI_Barrier(c);
#endif

    // Turn back on LA threads
    SetLAThreads(LAThreads);


  }

// SS: formFXC for GIAO
// TODO: Merge this into general template
  template <> 
  template <> 
  void KohnSham<dcomplex, dcomplex>::formFXC( MPI_Comm c,  
    std::vector<TwoBodyContraction<dcomplex>> &cList, EMPerturbation &pert ) {

    // CErr("TDDFT + Complex Ints NYI",std::cout);


    size_t itOff = this->nC == 2 ? 5 : 3;
    size_t nVec = cList.size() / itOff;
    size_t NB     = this->basisSet().nBasis;
    size_t NB2    = NB*NB;
    size_t NPPB   = intParam.nRadPerBatch * intParam.nAng;
    size_t nAtoms = this->molecule().nAtoms;


    // Parallelism
    size_t NT = GetNumThreads();
    size_t LAThreads = GetLAThreads();
    size_t mpiRank = MPIRank(c);
    size_t mpiSize = MPISize(c);

    // Turn off LA threads
    SetLAThreads(1);

    // Split the MPI Comm
    int color = ((mpiSize < nAtoms) or 
                 (mpiRank < nAtoms)) ? 1 : MPI_UNDEFINED;
                                                            
                  
    MPI_Comm intComm = MPICommSplit(c,color,mpiRank);


    bool isGGA = std::any_of(functionals.begin(),functionals.end(),
                   [](std::shared_ptr<DFTFunctional> &x) {
                     return x->isGGA(); 
                   }); 

#ifdef CQ_ENABLE_MPI
    if( intComm != MPI_COMM_NULL ) {
#endif

    dcomplex* NBNBSCR = CQMemManager::get().malloc<dcomplex>(NB2 * NT);
    dcomplex* NBNPSCR = CQMemManager::get().malloc<dcomplex>(NB*NPPB * NT);

    std::shared_ptr<cqmatrix::PauliSpinorMatrices<dcomplex>> Re1PDM;
    Re1PDM = std::dynamic_pointer_cast<cqmatrix::PauliSpinorMatrices<dcomplex>>(
      this->onePDM);

//FIXME
    std::vector<std::vector<dcomplex *>> ReTSymm; 
    for(auto iVec = 0; iVec < nVec; iVec++) {
      ReTSymm.emplace_back();

      size_t indx = iVec * itOff;
      for(auto iS = 0; iS < 2*this->nC; iS++) {
        
         ReTSymm.back().emplace_back(
             CQMemManager::get().malloc<dcomplex>(NB2));
         SetMat('N',NB,NB,dcomplex(1.0),cList[indx + iS + 1].X,NB,
           ReTSymm.back().back(),NB);

      }

    }



    dcomplex* GxcT_raw = CQMemManager::get().malloc<dcomplex>(2*this->nC*nVec*NT*NB2);
    dcomplex* Gxc_first = GxcT_raw;
    std::vector<std::vector<std::vector<dcomplex*>>> GxcT;
    for(auto ithread = 0; ithread < NT; ithread++) {
      GxcT.emplace_back();
      for(auto iVec = 0; iVec < nVec; iVec++) { 
        GxcT.back().emplace_back();
        for(auto iS = 0; iS < 2*this->nC; iS++){
          GxcT.back().back().emplace_back(Gxc_first);
          memset(GxcT.back().back().back(),0,NB2*sizeof(dcomplex));
          Gxc_first += NB2;
        }
      }
    }






    // Allocation of V Vairables
    double *DenS(nullptr),  *DenZ(nullptr),  *DenY(nullptr),  *DenX(nullptr);
    double *GDenS(nullptr), *GDenZ(nullptr), *GDenY(nullptr), *GDenX(nullptr);

//2c s
    // Allocation of V small aux var (RKS/UKS, GKS Msmall(false), GKS (Msmall(true)

    double *Mnorm(nullptr);  // |m| 
    // rho^K/|rho^k| for large m and 1/3 for small m
    // Eq. D.21a Shichao Sun's thesis 
    double *KScratch(nullptr);   
    // Eq. D.21b Shichao Sun's thesis 
    double *HScratch(nullptr);
    //double *IScratch(nullptr);   
    // Checks if m is below a threshold at that grid point 
    bool   *Msmall(nullptr);
    // |del rho * del m| 
    double *DSDMnorm(nullptr);
    // Eq. D.6 in Shichao Sun's thesis
    double *signMD(nullptr);

//2c e

    // Density
    DenS = CQMemManager::get().malloc<double>(NPPB * NT);
    if( this->nC == 2 or not this->iCS )
      DenZ = CQMemManager::get().malloc<double>(NPPB * NT);
    if( this->nC == 2 ) {
      DenY = CQMemManager::get().malloc<double>(NPPB * NT);
      DenX = CQMemManager::get().malloc<double>(NPPB * NT);

//2c s

      Mnorm    = CQMemManager::get().malloc<double>(NPPB * NT);
      KScratch = CQMemManager::get().malloc<double>(3 * NPPB * NT); // 3 is for K=x,y,z
      //IScratch = CQMemManager::get().malloc<double>(3 * NPPB * NT); // 3 is for K=x,y,z
      Msmall   = CQMemManager::get().malloc<bool>(NPPB * NT);
//2c e

    }


    // Density Gradient
    if( isGGA ) {

      GDenS = CQMemManager::get().malloc<double>(3*NPPB * NT);
      if( this->nC == 2 or not this->iCS )
        GDenZ = CQMemManager::get().malloc<double>(3*NPPB * NT);
      if( this->nC == 2 ) {
        GDenY = CQMemManager::get().malloc<double>(3*NPPB * NT);
        GDenX = CQMemManager::get().malloc<double>(3*NPPB * NT);

//2c s
        HScratch = CQMemManager::get().malloc<double>(3*NPPB * NT);
        DSDMnorm    = CQMemManager::get().malloc<double>(NPPB * NT);
        signMD      = CQMemManager::get().malloc<double>(NPPB * NT);

//2c e
      }

    }

    // Allocation of T evaluation
    dcomplex *TS(nullptr),  *TZ(nullptr),  *TY(nullptr),  *TX(nullptr);
    dcomplex *GTS(nullptr), *GTZ(nullptr), *GTY(nullptr), *GTX(nullptr);

    // Allocation of gPT variables 
    dcomplex *gPTss(nullptr), *gPTsz(nullptr), *gPTsy(nullptr), *gPTsx(nullptr), 
     *gPTzz(nullptr), *gPTyy(nullptr), *gPTxx(nullptr); 

    // T
    TS = CQMemManager::get().malloc<dcomplex>(NPPB * NT);
    TZ = CQMemManager::get().malloc<dcomplex>(NPPB * NT);
    if( this->nC == 2 ) {
      TY = CQMemManager::get().malloc<dcomplex>(NPPB * NT);
      TX = CQMemManager::get().malloc<dcomplex>(NPPB * NT);
    }


    // T Gradient
    if( isGGA ) {

      GTS = CQMemManager::get().malloc<dcomplex>(3*NPPB * NT);
      GTZ = CQMemManager::get().malloc<dcomplex>(3*NPPB * NT); 
      
      gPTss = CQMemManager::get().malloc<dcomplex>(NPPB * NT);
      gPTsz = CQMemManager::get().malloc<dcomplex>(NPPB * NT);
      
      if( this->onePDM->hasZ() ) 
        gPTzz = CQMemManager::get().malloc<dcomplex>(NPPB * NT);

      if( this->onePDM->hasXY() ) {
        GTY = CQMemManager::get().malloc<dcomplex>(3*NPPB * NT);
        GTX = CQMemManager::get().malloc<dcomplex>(3*NPPB * NT);

        gPTsx = CQMemManager::get().malloc<dcomplex>(NPPB * NT);
        gPTsy = CQMemManager::get().malloc<dcomplex>(NPPB * NT);
        gPTyy = CQMemManager::get().malloc<dcomplex>(NPPB * NT);
        gPTxx = CQMemManager::get().malloc<dcomplex>(NPPB * NT);
      }

    }





    // U Variables

    double * eps     = CQMemManager::get().malloc<double>(NPPB * NT);
    double * U_n     = CQMemManager::get().malloc<double>(2*NPPB * NT);
    double * U_gamma = isGGA ? 
      CQMemManager::get().malloc<double>(3*NPPB * NT) : nullptr;

    // U First Derivatives

    double * dVU_n     = CQMemManager::get().malloc<double>(2*NPPB * NT);
    double * dVU_gamma = isGGA ? 
      CQMemManager::get().malloc<double>(3*NPPB * NT) : nullptr;

    double * eps_SCR(nullptr), * dVU_n_SCR(nullptr), * dVU_gamma_SCR(nullptr);
    if( functionals.size() > 1 ) {

      eps_SCR       = CQMemManager::get().malloc<double>(NPPB * NT);
      dVU_n_SCR     = CQMemManager::get().malloc<double>(2*NPPB * NT);
      dVU_gamma_SCR = isGGA ? 
        CQMemManager::get().malloc<double>(3*NPPB * NT) : nullptr;

    }

    // U Second Derivatives

    double * d2VU_n = CQMemManager::get().malloc<double>(3*NPPB * NT);
    double * d2VU_gamma   = isGGA ? 
      CQMemManager::get().malloc<double>(6*NPPB * NT) : nullptr;
    double * d2VU_n_gamma = isGGA ? 
      CQMemManager::get().malloc<double>(6*NPPB * NT) : nullptr;

    double * d2VU_n_SCR(nullptr), * d2VU_gamma_SCR(nullptr), 
           * d2VU_n_gamma_SCR(nullptr); 

    if( functionals.size() > 1 ) {
      d2VU_n_SCR = CQMemManager::get().malloc<double>(3*NPPB * NT);
      d2VU_gamma_SCR   = isGGA ? 
        CQMemManager::get().malloc<double>(6*NPPB * NT) : nullptr;
      d2VU_n_gamma_SCR = isGGA ? 
        CQMemManager::get().malloc<double>(6*NPPB * NT) : nullptr;
    }


    // Z Vars
// not sure Z var should be complex or real : have to be complex in this case! 
    dcomplex * ZrhoVar = CQMemManager::get().malloc<dcomplex>(NPPB * NT);

    dcomplex * ZgammaVar1 = isGGA ? 
      CQMemManager::get().malloc<dcomplex>(NPPB * NT) : nullptr;
    dcomplex * ZgammaVar2 = isGGA ? 
      CQMemManager::get().malloc<dcomplex>(NPPB * NT) : nullptr;
    dcomplex * ZgammaVar3 = isGGA ? 
      CQMemManager::get().malloc<dcomplex>(NPPB * NT) : nullptr;
    dcomplex * ZgammaVar4 = isGGA ? 
      CQMemManager::get().malloc<dcomplex>(NPPB * NT) : nullptr;

    dcomplex* ZMAT = CQMemManager::get().malloc<dcomplex>(NB*NPPB * NT);

//2c s

    dcomplex*  Basis_cmplx = CQMemManager::get().malloc<dcomplex>(NB*NPPB*NT);

//2c e


    double intDen = 0.;

    auto fxcbuild = [&](size_t &res, std::vector<cart_t> &batch, 
      std::vector<double> &weights, std::vector<size_t> NBE_vec, 
      std::vector<dcomplex*> BasisEval_vec, 
      std::vector<std::vector<size_t>>& batchEvalShells_vec, 
      std::vector<std::vector<std::pair<size_t,size_t>>>& subMatCut_vec) {

      size_t NBE = NBE_vec[0];
      dcomplex * BasisEval = BasisEval_vec[0];
      std::vector<size_t> & batchEvalShells = batchEvalShells_vec[0];
      std::vector<std::pair<size_t,size_t>> & subMatCut = subMatCut_vec[0];

      double epsScreen = intParam.epsilon / nAtoms /
        intParam.nAng / intParam.nRad;

      epsScreen = std::max(epsScreen,
                           std::numeric_limits<double>::epsilon());



      size_t NPts = batch.size();
      size_t IOff = NBE*NPts;
      int tid = GetThreadID();
      size_t nPtOff = tid * NPPB;


      // Get thread Local storage
      

      // Density and T
      double *DenS_loc = DenS + nPtOff;
      double *DenZ_loc = DenZ + nPtOff;
      double *DenY_loc = DenY + nPtOff;
      double *DenX_loc = DenX + nPtOff;

      double *GDenS_loc = GDenS + 3 * nPtOff;
      double *GDenZ_loc = GDenZ + 3 * nPtOff;
      double *GDenY_loc = GDenY + 3 * nPtOff;
      double *GDenX_loc = GDenX + 3 * nPtOff;

//2c s

      // Aux Var for 2C 
      double * Mnorm_loc    = Mnorm        +   nPtOff;
      double * KScratch_loc = KScratch     + 3*nPtOff;
      bool   * Msmall_loc   = Msmall       +   nPtOff;
      double * HScratch_loc = HScratch     + 3*nPtOff;

      double * DSDMnorm_loc = DSDMnorm     +   nPtOff;
      double * signMD_loc   = signMD       +   nPtOff;

//2c e


      dcomplex *TS_loc = TS + nPtOff;
      dcomplex *TZ_loc = TZ + nPtOff;
      dcomplex *TY_loc = TY + nPtOff;
      dcomplex *TX_loc = TX + nPtOff;

// SS s
      dcomplex *gPTss_loc = gPTss + nPtOff;
      dcomplex *gPTsz_loc = gPTsz + nPtOff;
      dcomplex *gPTsy_loc = gPTsy + nPtOff;
      dcomplex *gPTsx_loc = gPTsx + nPtOff;
      dcomplex *gPTxx_loc = gPTxx + nPtOff;
      dcomplex *gPTyy_loc = gPTyy + nPtOff;
      dcomplex *gPTzz_loc = gPTzz + nPtOff;
// SS e

      dcomplex *GTS_loc = GTS + 3 * nPtOff;
      dcomplex *GTZ_loc = GTZ + 3 * nPtOff;
      dcomplex *GTY_loc = GTY + 3 * nPtOff;
      dcomplex *GTX_loc = GTX + 3 * nPtOff;


      // U Vars

      double *U_n_loc        = U_n        + 2 * nPtOff;     
      double *U_gamma_loc    = U_gamma    + 3 * nPtOff;

      double *eps_loc        = eps        +     nPtOff;
      double *dVU_n_loc      = dVU_n      + 2 * nPtOff;     
      double *dVU_gamma_loc  = dVU_gamma  + 3 * nPtOff;

      double *d2VU_n_loc       = d2VU_n       + 3 * nPtOff;     
      double *d2VU_gamma_loc   = d2VU_gamma   + 6 * nPtOff;
      double *d2VU_n_gamma_loc = d2VU_n_gamma + 6 * nPtOff;

      double *eps_SCR_loc       = eps_SCR       +     nPtOff;
      double *dVU_n_SCR_loc     = dVU_n_SCR     + 2 * nPtOff;     
      double *dVU_gamma_SCR_loc = dVU_gamma_SCR + 3 * nPtOff;

      double *d2VU_n_SCR_loc       = d2VU_n_SCR       + 3 * nPtOff;     
      double *d2VU_gamma_SCR_loc   = d2VU_gamma_SCR   + 6 * nPtOff;
      double *d2VU_n_gamma_SCR_loc = d2VU_n_gamma_SCR + 6 * nPtOff;


      // Z Vars

      dcomplex* ZrhoVar_loc    = ZrhoVar    + nPtOff;
      dcomplex* ZgammaVar1_loc = ZgammaVar1 + nPtOff;
      dcomplex* ZgammaVar2_loc = ZgammaVar2 + nPtOff;
      dcomplex* ZgammaVar3_loc = ZgammaVar3 + nPtOff;
      dcomplex* ZgammaVar4_loc = ZgammaVar4 + nPtOff;

      dcomplex* ZMAT_loc = ZMAT + nPtOff*NB;


      dcomplex* NBNBSCR_loc = NBNBSCR + NB * NB   * tid;
      dcomplex* NBNPSCR_loc = NBNPSCR + NB * NPPB * tid;

// they are scratch spaces (don't know what NBNBSCR_r and NBNPSCR_r is) 
      dcomplex * NBNBSCR_r = NBNBSCR_loc;
      dcomplex * NBNPSCR_r = NBNPSCR_loc;


//2c s

      // Basis
      
      dcomplex* Basis_use = reinterpret_cast<dcomplex*>(BasisEval);
      if( Basis_cmplx ) {
        Basis_use = Basis_cmplx + nPtOff*NB;
        std::copy_n(BasisEval, NPts * NBE, Basis_use);
      } 

      dcomplex* Basis2 = BasisEval;
      IMatCopy('N', NBE, NPPB*NT, 1.0, Basis2, NBE, NBE);

//2c e


      // This evaluates the V variables for all components 
      // (Scalar, MZ (UKS) and Mx, MY (2 Comp))
      evalDen((isGGA ? GRADIENT : NOGRAD), NPts, NBE, NB, subMatCut, 
        NBNBSCR_r, NBNPSCR_r, Re1PDM->S().pointer(), DenS_loc , 
        GDenS_loc, GDenS_loc + NPts, GDenS_loc + 2*NPts, BasisEval);

      if( this->onePDM->hasZ() )
        evalDen((isGGA ? GRADIENT : NOGRAD), NPts, NBE, NB, subMatCut, 
          NBNBSCR_r ,NBNPSCR_r, Re1PDM->Z().pointer(), DenZ_loc, 
          GDenZ_loc, GDenZ_loc + NPts, GDenZ_loc + 2*NPts, BasisEval);

      if( this->onePDM->hasXY() ) {
        evalDen((isGGA ? GRADIENT : NOGRAD), NPts, NBE, NB, subMatCut, 
          NBNBSCR_r ,NBNPSCR_r, Re1PDM->Y().pointer(), DenY_loc, 
          GDenY_loc, GDenY_loc + NPts, GDenY_loc + 2*NPts, BasisEval);
        evalDen((isGGA ? GRADIENT : NOGRAD), NPts, NBE, NB, subMatCut, 
          NBNBSCR_r ,NBNPSCR_r, Re1PDM->X().pointer(), DenX_loc, 
          GDenX_loc, GDenX_loc + NPts, GDenX_loc + 2*NPts, BasisEval);
      
      }

      // V -> U variables for evaluating the kernel derivatives.
      mkAuxVar(this->onePDM,isGGA,epsScreen,NPts,
        DenS_loc,DenZ_loc,DenY_loc,DenX_loc,
        GDenS_loc,GDenS_loc + NPts,GDenS_loc + 2*NPts,
        GDenZ_loc,GDenZ_loc + NPts,GDenZ_loc + 2*NPts,
        GDenY_loc,GDenY_loc + NPts,GDenY_loc + 2*NPts,
        GDenX_loc,GDenX_loc + NPts,GDenX_loc + 2*NPts,
//        nullptr, 
//        nullptr, nullptr, nullptr,
//        nullptr, nullptr, nullptr,
//        nullptr,U_n_loc,U_gamma_loc
//2c s
        Mnorm_loc, 
        KScratch_loc, KScratch_loc + NPts, KScratch_loc + 2* NPts,
        HScratch_loc, HScratch_loc + NPts, HScratch_loc + 2* NPts,
        DSDMnorm_loc, signMD_loc, 
        Msmall_loc,U_n_loc,U_gamma_loc
//2c e
      );



      loadFXCder(NPts,U_n_loc,U_gamma_loc,eps_loc,dVU_n_loc,d2VU_n_loc,
        dVU_gamma_loc, d2VU_gamma_loc,d2VU_n_gamma_loc,eps_SCR_loc,
        dVU_n_SCR_loc,dVU_gamma_SCR_loc, d2VU_n_SCR_loc,d2VU_gamma_SCR_loc,
        d2VU_n_gamma_SCR_loc);


      for(auto iT = 0; iT < nVec; iT++) {

        // if( std::is_same<U,dcomplex>::value ) CErr("NO COMPLEX YET!");

        size_t indx       = itOff * iT;
        dcomplex *TS_d      = TS_loc;
        dcomplex *GTS_d     = GTS_loc;
        dcomplex *TZ_d      = TZ_loc;
        dcomplex *GTZ_d     = GTZ_loc;

        // This evaluates the V variables for all components 
        // (Scalar, MZ (UKS) and Mx, MY (2 Comp))
        evalDen((isGGA ? GRADIENT : NOGRAD), NPts, NBE, NB, subMatCut, 
          NBNBSCR_r, NBNPSCR_r, ReTSymm[iT][SCALAR] , TS_d, 
          GTS_d, GTS_d + NPts, GTS_d + 2*NPts, BasisEval);

        evalDen((isGGA ? GRADIENT : NOGRAD), NPts, NBE, NB, subMatCut, 
          NBNBSCR_r ,NBNPSCR_r, ReTSymm[iT][MZ] , TZ_d, 
          GTZ_d, GTZ_d + NPts, GTZ_d + 2*NPts, BasisEval);

#ifdef debugGIAOTDDFT  // SS start
  std::cout<<"print GTS x "<<std::endl;

  for (int ii = 0 ; ii < NPts ; ii++ )
    std::cout<<std::real(GTS_d[ii])<<std::endl;

  std::cout<<"GTS x finished "<<std::endl;

  std::cout<<"print GTS y "<<std::endl;

  for (int ii = 0 ; ii < NPts ; ii++ )
    std::cout<<std::real(GTS_d[NPts+ii])<<std::endl;

  std::cout<<"GTS y finished "<<std::endl;

  std::cout<<"print GTS z "<<std::endl;

  for (int ii = 0 ; ii < NPts ; ii++ )
    std::cout<<std::real(GTS_d[2*NPts+ii])<<std::endl;

  std::cout<<"GTS z finished "<<std::endl;

#endif   // SS end



        if( this->onePDM->hasXY()) {

          dcomplex *TY_d      = reinterpret_cast<dcomplex*>(TY_loc);
          dcomplex *GTY_d     = reinterpret_cast<dcomplex*>(GTY_loc);
          dcomplex *TX_d      = reinterpret_cast<dcomplex*>(TX_loc);
          dcomplex *GTX_d     = reinterpret_cast<dcomplex*>(GTX_loc);

          evalDen((isGGA ? GRADIENT : NOGRAD), NPts, NBE, NB, subMatCut, 
            NBNBSCR_r ,NBNPSCR_r, ReTSymm[iT][MY], TY_d, 
            GTY_d, GTY_d + NPts, GTY_d + 2*NPts, BasisEval);

          evalDen((isGGA ? GRADIENT : NOGRAD), NPts, NBE, NB, subMatCut, 
            NBNBSCR_r ,NBNPSCR_r, ReTSymm[iT][MX], TX_d, 
            GTX_d, GTX_d + NPts, GTX_d + 2*NPts, BasisEval);



        }

   
        // need to template constructZVarsFXC  
// this part is for collinear
/*
        constructZVarsFXC(SCALAR,isGGA,NPts, GDenS_loc, GDenZ_loc, GDenY_loc, 
          GDenX_loc, TS_loc, TZ_loc, TY_loc, TX_loc, GTS_loc, GTZ_loc, GTY_loc, 
          GTX_loc, dVU_n_loc, dVU_gamma_loc, d2VU_n_loc, d2VU_gamma_loc, 
          d2VU_n_gamma_loc, ZrhoVar_loc, ZgammaVar1_loc, ZgammaVar2_loc,
          ZgammaVar3_loc, ZgammaVar4_loc);

        formZ_fxc(SCALAR,isGGA,NPts,NBE,IOff,epsScreen,weights,ZrhoVar_loc,
          ZgammaVar1_loc, ZgammaVar2_loc, ZgammaVar3_loc, ZgammaVar4_loc, 
          GDenS_loc, GDenZ_loc, GDenY_loc, GDenX_loc, GTS_loc, GTZ_loc, 
          GTY_loc, GTX_loc, BasisEval, ZMAT_loc);

*/

//std::cout<<"iT="<<iT<<" SCALAR"<<std::endl;
        if ( isGGA and  this->onePDM->hasXY() ) {  
          mkgPTVar( NPts, 
            GDenS_loc, GDenZ_loc, GDenY_loc, GDenX_loc, 
            GTS_loc, GTZ_loc, GTY_loc, GTX_loc,
            gPTss_loc, gPTsz_loc, gPTsy_loc, gPTsx_loc, gPTzz_loc, gPTyy_loc, gPTxx_loc
            );
        } 

        constructZVarsFXC(SCALAR,isGGA,NPts, GDenS_loc, GDenZ_loc, GDenY_loc, GDenX_loc, 
          Msmall_loc,Mnorm_loc,
          KScratch_loc, KScratch_loc + NPts, KScratch_loc + 2* NPts,
          HScratch_loc, HScratch_loc + NPts, HScratch_loc + 2* NPts,
          TS_loc, TZ_loc, TY_loc, TX_loc, GTS_loc, GTZ_loc, GTY_loc, 
          GTX_loc, 
          gPTss_loc, gPTsz_loc, gPTsy_loc, gPTsx_loc, gPTzz_loc, gPTyy_loc, gPTxx_loc, 
          dVU_n_loc, dVU_gamma_loc, d2VU_n_loc, d2VU_gamma_loc, 
          d2VU_n_gamma_loc, ZrhoVar_loc, ZgammaVar1_loc, ZgammaVar2_loc,
          ZgammaVar3_loc, ZgammaVar4_loc);

        formZ_fxc(SCALAR,isGGA,NPts,NBE,IOff,epsScreen,weights,ZrhoVar_loc,
          ZgammaVar1_loc, ZgammaVar2_loc, ZgammaVar3_loc, ZgammaVar4_loc, 
          Msmall_loc,Mnorm_loc,
//          DenS_loc, DenZ_loc, DenY_loc, DenX_loc,  
          DSDMnorm_loc, signMD_loc, 
          GDenS_loc, GDenZ_loc, GDenY_loc, GDenX_loc, 
          KScratch_loc, KScratch_loc + NPts, KScratch_loc + 2* NPts,
          HScratch_loc, HScratch_loc + NPts, HScratch_loc + 2* NPts,
          GTS_loc, GTZ_loc, GTY_loc, GTX_loc,
          gPTss_loc, gPTsz_loc, gPTsy_loc, gPTsx_loc, gPTzz_loc, gPTyy_loc, gPTxx_loc, 
          BasisEval, ZMAT_loc);



        // double *ZMAT_r = reinterpret_cast<double*>(ZMAT_loc);

        for ( int icount = 0 ; icount < NBE*NPts ; icount++ )
          BasisEval[icount] = std::conj(BasisEval[icount]);

        blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::Trans,NBE,NBE,NPts,dcomplex(0.5),ZMAT_loc,NBE, 
          BasisEval,NBE, dcomplex(0.0),NBNBSCR_r,NBE);

        formZ_fxc(SCALAR,isGGA,NPts,NBE,IOff,epsScreen,weights,ZrhoVar_loc,
          ZgammaVar1_loc, ZgammaVar2_loc, ZgammaVar3_loc, ZgammaVar4_loc, 
          Msmall_loc,Mnorm_loc,
          //DenS_loc, DenZ_loc, DenY_loc, DenX_loc,  
          DSDMnorm_loc, signMD_loc, 
          GDenS_loc, GDenZ_loc, GDenY_loc, GDenX_loc, 
          KScratch_loc, KScratch_loc + NPts, KScratch_loc + 2* NPts,
          HScratch_loc, HScratch_loc + NPts, HScratch_loc + 2* NPts,
          GTS_loc, GTZ_loc, GTY_loc, GTX_loc, 
          gPTss_loc, gPTsz_loc, gPTsy_loc, gPTsx_loc, gPTzz_loc, gPTyy_loc, gPTxx_loc, 
          BasisEval, ZMAT_loc);

        for ( int icount = 0 ; icount < NBE*NPts ; icount++ )
          BasisEval[icount] = std::conj(BasisEval[icount]);

        blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::Trans,NBE,NBE,NPts,dcomplex(0.5),BasisEval,NBE, 
          ZMAT_loc,NBE, dcomplex(1.0),NBNBSCR_r,NBE);
         

        // DSYR2K('L','N',NBE,NPts,0.5,BasisEval,NBE,ZMAT_r,NBE,0.,NBNBSCR_r,NBE);

        IncBySubMat(NB,NB,NBE,NBE,GxcT[tid][iT][SCALAR],NB,NBNBSCR_loc,NBE,
            subMatCut);




// this part is collinear
/*
        constructZVarsFXC(MZ,isGGA,NPts, GDenS_loc, GDenZ_loc, GDenY_loc, 
          GDenX_loc, TS_loc, TZ_loc, TY_loc, TX_loc, GTS_loc, GTZ_loc, GTY_loc, 
          GTX_loc, dVU_n_loc, dVU_gamma_loc, d2VU_n_loc, d2VU_gamma_loc, 
          d2VU_n_gamma_loc, ZrhoVar_loc, ZgammaVar1_loc, ZgammaVar2_loc,
          ZgammaVar3_loc, ZgammaVar4_loc);

        formZ_fxc(MZ,isGGA,NPts,NBE,IOff,epsScreen,weights,ZrhoVar_loc,
          ZgammaVar1_loc, ZgammaVar2_loc, ZgammaVar3_loc, ZgammaVar4_loc, 
          GDenS_loc, GDenZ_loc, GDenY_loc, GDenX_loc, GTS_loc, GTZ_loc, 
          GTY_loc, GTX_loc, BasisEval, ZMAT_loc);
*/

//std::cout<<"iT="<<iT<<" MZ, tid"<<tid<<std::endl;
        constructZVarsFXC(MZ,isGGA,NPts, GDenS_loc, GDenZ_loc, GDenY_loc, GDenX_loc, 
          Msmall_loc,Mnorm_loc,
          KScratch_loc, KScratch_loc + NPts, KScratch_loc + 2* NPts,
          HScratch_loc, HScratch_loc + NPts, HScratch_loc + 2* NPts,
          TS_loc, TZ_loc, TY_loc, TX_loc, GTS_loc, GTZ_loc, GTY_loc, 
          GTX_loc, 
          gPTss_loc, gPTsz_loc, gPTsy_loc, gPTsx_loc, gPTzz_loc, gPTyy_loc, gPTxx_loc, 
          dVU_n_loc, dVU_gamma_loc, d2VU_n_loc, d2VU_gamma_loc, 
          d2VU_n_gamma_loc, ZrhoVar_loc, ZgammaVar1_loc, ZgammaVar2_loc,
          ZgammaVar3_loc, ZgammaVar4_loc);

        formZ_fxc(MZ,isGGA,NPts,NBE,IOff,epsScreen,weights,ZrhoVar_loc,
          ZgammaVar1_loc, ZgammaVar2_loc, ZgammaVar3_loc, ZgammaVar4_loc,  
          Msmall_loc,Mnorm_loc,
          //DenS_loc, DenZ_loc, DenY_loc, DenX_loc,  
          DSDMnorm_loc, signMD_loc,
          GDenS_loc, GDenZ_loc, GDenY_loc, GDenX_loc, 
          KScratch_loc, KScratch_loc + NPts, KScratch_loc + 2* NPts,
          HScratch_loc, HScratch_loc + NPts, HScratch_loc + 2* NPts,
          GTS_loc, GTZ_loc, GTY_loc, GTX_loc, 
          gPTss_loc, gPTsz_loc, gPTsy_loc, gPTsx_loc, gPTzz_loc, gPTyy_loc, gPTxx_loc, 
          BasisEval, ZMAT_loc);


//        DSYR2K('L','N',NBE,NPts,0.5,BasisEval,NBE,ZMAT_r,NBE,0.,NBNBSCR_r,NBE);

        for ( int icount = 0 ; icount < NBE*NPts ; icount++ )
          BasisEval[icount] = std::conj(BasisEval[icount]);

        blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::Trans,NBE,NBE,NPts,dcomplex(0.5),ZMAT_loc,NBE, 
          BasisEval,NBE, dcomplex(0.0),NBNBSCR_r,NBE);

        formZ_fxc(MZ,isGGA,NPts,NBE,IOff,epsScreen,weights,ZrhoVar_loc,
          ZgammaVar1_loc, ZgammaVar2_loc, ZgammaVar3_loc, ZgammaVar4_loc,  
          Msmall_loc,Mnorm_loc,
          //DenS_loc, DenZ_loc, DenY_loc, DenX_loc,  
          DSDMnorm_loc, signMD_loc, 
          GDenS_loc, GDenZ_loc, GDenY_loc, GDenX_loc, 
          KScratch_loc, KScratch_loc + NPts, KScratch_loc + 2* NPts,
          HScratch_loc, HScratch_loc + NPts, HScratch_loc + 2* NPts,
          GTS_loc, GTZ_loc, GTY_loc, GTX_loc, 
          gPTss_loc, gPTsz_loc, gPTsy_loc, gPTsx_loc, gPTzz_loc, gPTyy_loc, gPTxx_loc, 
          BasisEval, ZMAT_loc);

        for ( int icount = 0 ; icount < NBE*NPts ; icount++ )
          BasisEval[icount] = std::conj(BasisEval[icount]);

        blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::Trans,NBE,NBE,NPts,dcomplex(0.5),BasisEval,NBE, 
          ZMAT_loc,NBE, dcomplex(1.0),NBNBSCR_r,NBE);
         




        IncBySubMat(NB,NB,NBE,NBE,GxcT[tid][iT][MZ],NB,NBNBSCR_loc,NBE,
            subMatCut);


//2c s
        if( this->onePDM->hasXY() ) {
          //std::cerr << "MX bit "<< std::endl;
//std::cout<<"iT="<<iT<<" MX, tid"<<tid<<std::endl;
          constructZVarsFXC(MX,isGGA,NPts, GDenS_loc, GDenZ_loc, GDenY_loc, GDenX_loc, 
            Msmall_loc,Mnorm_loc,
            KScratch_loc, KScratch_loc + NPts, KScratch_loc + 2* NPts,
            HScratch_loc, HScratch_loc + NPts, HScratch_loc + 2* NPts,
            TS_loc, TZ_loc, TY_loc, TX_loc, GTS_loc, GTZ_loc, GTY_loc, 
            GTX_loc, 
            gPTss_loc, gPTsz_loc, gPTsy_loc, gPTsx_loc, gPTzz_loc, gPTyy_loc, gPTxx_loc, 
            dVU_n_loc, dVU_gamma_loc, d2VU_n_loc, d2VU_gamma_loc, 
            d2VU_n_gamma_loc, ZrhoVar_loc, ZgammaVar1_loc, ZgammaVar2_loc,
            ZgammaVar3_loc, ZgammaVar4_loc);

          formZ_fxc(MX,isGGA,NPts,NBE,IOff,epsScreen,weights,ZrhoVar_loc,
            ZgammaVar1_loc, ZgammaVar2_loc, ZgammaVar3_loc, ZgammaVar4_loc, 
            Msmall_loc,Mnorm_loc,
            //DenS_loc, DenZ_loc, DenY_loc, DenX_loc,  
            DSDMnorm_loc, signMD_loc,
            GDenS_loc, GDenZ_loc, GDenY_loc, GDenX_loc, 
            KScratch_loc, KScratch_loc + NPts, KScratch_loc + 2* NPts,
            HScratch_loc, HScratch_loc + NPts, HScratch_loc + 2* NPts,
            GTS_loc, GTZ_loc, GTY_loc, GTX_loc, 
            gPTss_loc, gPTsz_loc, gPTsy_loc, gPTsx_loc, gPTzz_loc, gPTyy_loc, gPTxx_loc, 
            BasisEval, ZMAT_loc);


          for ( int icount = 0 ; icount < NBE*NPts ; icount++ )
            BasisEval[icount] = std::conj(BasisEval[icount]);
        
          blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::Trans,NBE,NBE,NPts,dcomplex(0.5),ZMAT_loc,NBE, 
            BasisEval,NBE, dcomplex(0.0),NBNBSCR_r,NBE);
        
          formZ_fxc(MX,isGGA,NPts,NBE,IOff,epsScreen,weights,ZrhoVar_loc,
            ZgammaVar1_loc, ZgammaVar2_loc, ZgammaVar3_loc, ZgammaVar4_loc, 
            Msmall_loc,Mnorm_loc,
            //DenS_loc, DenZ_loc, DenY_loc, DenX_loc,  
            DSDMnorm_loc, signMD_loc,
            GDenS_loc, GDenZ_loc, GDenY_loc, GDenX_loc, 
            KScratch_loc, KScratch_loc + NPts, KScratch_loc + 2* NPts,
            HScratch_loc, HScratch_loc + NPts, HScratch_loc + 2* NPts,
            GTS_loc, GTZ_loc, GTY_loc, GTX_loc, 
            gPTss_loc, gPTsz_loc, gPTsy_loc, gPTsx_loc, gPTzz_loc, gPTyy_loc, gPTxx_loc, 
            BasisEval, ZMAT_loc);
        
          for ( int icount = 0 ; icount < NBE*NPts ; icount++ )
            BasisEval[icount] = std::conj(BasisEval[icount]);
        
          blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::Trans,NBE,NBE,NPts,dcomplex(0.5),BasisEval,NBE, 
            ZMAT_loc,NBE, dcomplex(1.0),NBNBSCR_r,NBE);
         




          IncBySubMat(NB,NB,NBE,NBE,GxcT[tid][iT][MX],NB,NBNBSCR_loc,NBE,
              subMatCut);

          //std::cerr << "MY bit "<< std::endl;
//std::cout<<"iT="<<iT<<" MY, tid"<<tid<<std::endl;
          constructZVarsFXC(MY,isGGA,NPts, GDenS_loc, GDenZ_loc, GDenY_loc, GDenX_loc, 
            Msmall_loc,Mnorm_loc,
            KScratch_loc, KScratch_loc + NPts, KScratch_loc + 2* NPts,
            HScratch_loc, HScratch_loc + NPts, HScratch_loc + 2* NPts,
            TS_loc, TZ_loc, TY_loc, TX_loc, GTS_loc, GTZ_loc, GTY_loc, 
            GTX_loc, 
            gPTss_loc, gPTsz_loc, gPTsy_loc, gPTsx_loc, gPTzz_loc, gPTyy_loc, gPTxx_loc, 
            dVU_n_loc, dVU_gamma_loc, d2VU_n_loc, d2VU_gamma_loc, 
            d2VU_n_gamma_loc, ZrhoVar_loc, ZgammaVar1_loc, ZgammaVar2_loc,
            ZgammaVar3_loc, ZgammaVar4_loc);

          formZ_fxc(MY,isGGA,NPts,NBE,IOff,epsScreen,weights,ZrhoVar_loc,
            ZgammaVar1_loc, ZgammaVar2_loc, ZgammaVar3_loc, ZgammaVar4_loc,  
            Msmall_loc,Mnorm_loc,
            //DenS_loc, DenZ_loc, DenY_loc, DenX_loc,  
            DSDMnorm_loc, signMD_loc,
            GDenS_loc, GDenZ_loc, GDenY_loc, GDenX_loc, 
            KScratch_loc, KScratch_loc + NPts, KScratch_loc + 2* NPts,
            HScratch_loc, HScratch_loc + NPts, HScratch_loc +  NPts,
            GTS_loc, GTZ_loc, GTY_loc, GTX_loc, 
            gPTss_loc, gPTsz_loc, gPTsy_loc, gPTsx_loc, gPTzz_loc, gPTyy_loc, gPTxx_loc, 
            BasisEval, ZMAT_loc);




          for ( int icount = 0 ; icount < NBE*NPts ; icount++ )
            BasisEval[icount] = std::conj(BasisEval[icount]);
        
          blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::Trans,NBE,NBE,NPts,dcomplex(0.5),ZMAT_loc,NBE, 
            BasisEval,NBE, dcomplex(0.0),NBNBSCR_r,NBE);
        
          formZ_fxc(MY,isGGA,NPts,NBE,IOff,epsScreen,weights,ZrhoVar_loc,
            ZgammaVar1_loc, ZgammaVar2_loc, ZgammaVar3_loc, ZgammaVar4_loc,  
            Msmall_loc,Mnorm_loc,
            //DenS_loc, DenZ_loc, DenY_loc, DenX_loc,  
            DSDMnorm_loc, signMD_loc,
            GDenS_loc, GDenZ_loc, GDenY_loc, GDenX_loc, 
            KScratch_loc, KScratch_loc + NPts, KScratch_loc + 2* NPts,
            HScratch_loc, HScratch_loc + NPts, HScratch_loc +  NPts,
            GTS_loc, GTZ_loc, GTY_loc, GTX_loc, 
            gPTss_loc, gPTsz_loc, gPTsy_loc, gPTsx_loc, gPTzz_loc, gPTyy_loc, gPTxx_loc, 
            BasisEval, ZMAT_loc);
        
          for ( int icount = 0 ; icount < NBE*NPts ; icount++ )
            BasisEval[icount] = std::conj(BasisEval[icount]);
        
          blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::Trans,NBE,NBE,NPts,dcomplex(0.5),BasisEval,NBE, 
            ZMAT_loc,NBE, dcomplex(1.0),NBNBSCR_r,NBE);
         

          IncBySubMat(NB,NB,NBE,NBE,GxcT[tid][iT][MY],NB,NBNBSCR_loc,NBE,
              subMatCut);

        } // 2C
//2c e

      } // iT loop

    }; // VXC


    // Create the BeckeIntegrator object
    BeckeIntegrator<EulerMac> 
      integrator(intComm,this->molecule(),
       this->basisSet(), EulerMac(intParam.nRad), intParam.nAng, 
        intParam.nRadPerBatch, (isGGA ? GRADIENT : NOGRAD), intParam.epsilon);

    // Integrate the FXC
    integrator.integrate<size_t>(fxcbuild,pert);

    //std::cerr << "NELEC " << 4.*M_PI*intDen << "\n";

//from here the commented part?   
    //for(auto iT = 0; iT < nVec; iT++) {

    //  // Add -Gxc[T] into K[T]
    //  for(auto iS = 0; iS < 2*this->nC; iS++) {

    //    for(auto ithread = 0; ithread < NT; ithread++) {
    //      HerMat('L',NB,GxcT[ithread][iT][iS],NB);

    //      U fact = (ithread == 0)  ? functionals.back()->xHFX : 1.;
    //      MatAdd('N','N',NB,NB,fact,cList[iT*itOff + iS + 1].AX,NB,
    //        U(-4. * M_PI), GxcT[ithread][iT][iS],NB,
    //        cList[iT*itOff + iS + 1].AX,NB);
    //    }


    //  }

    //}
// commented part end

    dcomplex* mpiScr = nullptr;
#ifdef CQ_ENABLE_MPI
    if( MPIRank(intComm) == 0 and MPISize(intComm) > 1 )
      mpiScr = CQMemManager::get().malloc<dcomplex>(NB*NB);
#endif

    for(auto iT = 0; iT < nVec; iT++) {

      for(auto iS = 0; iS < 2*this->nC; iS++) {

        // Add -Gxc[T] thread contributions into single storage
        for(auto ithread = 0; ithread < NT; ithread++) {
          //HerMat('L',NB,GxcT[ithread][iT][iS],NB);

          dcomplex fact = (ithread == 0)  ? 0 : 1.;
          MatAdd('N','N',NB,NB,
            fact         , GxcT[0][iT][iS]      ,NB,
            dcomplex(-4. * M_PI), GxcT[ithread][iT][iS],NB,
            GxcT[0][iT][iS],NB);
        }

#ifdef CQ_ENABLE_MPI
        // Add MPI Contributions together
        if( MPISize(intComm) > 1 )
          MPIReduce(GxcT[0][iT][iS],NB2,mpiScr,0,intComm);

        if( MPIRank(intComm) == 0 and MPISize(intComm) > 1 )
          std::copy_n(mpiScr,NB2,GxcT[0][iT][iS]);
#endif

#ifdef debugGIAOTDDFT //SS start
prettyPrintSmart(std::cout,"GxcT",GxcT[0][iT][iS],NB,NB,NB);
#endif //SS end

        // Add GxcT to K
        if( MPIRank(intComm) == 0 )
          MatAdd('N','N',NB,NB,
            dcomplex(functionals.back()->xHFX), cList[iT*itOff + iS + 1].AX, NB,
            dcomplex(1.),                       GxcT[0][iT][iS],             NB,
            cList[iT*itOff + iS + 1].AX, NB);

      }

    }

// need to check what is not freed 

    // Free up the memory
    if( mpiScr ) CQMemManager::get().free(mpiScr);
    CQMemManager::get().free( GxcT_raw, NBNBSCR, NBNPSCR );

    for(auto &Y : ReTSymm) for(auto &X : Y) CQMemManager::get().free(X);
    Re1PDM = nullptr;

    if( DenS ) CQMemManager::get().free( DenS );
    if( DenZ ) CQMemManager::get().free( DenZ );
    if( DenY ) CQMemManager::get().free( DenY );
    if( DenX ) CQMemManager::get().free( DenX );

// 2c s 
    if( Mnorm ) CQMemManager::get().free( Mnorm ); 
    if( KScratch ) CQMemManager::get().free( KScratch ); 
    if( Msmall ) CQMemManager::get().free( Msmall ); 
    if( HScratch ) CQMemManager::get().free( HScratch ); 
    if( DSDMnorm ) CQMemManager::get().free( DSDMnorm );
    if( signMD ) CQMemManager::get().free( signMD ); 
    if( Basis_cmplx ) CQMemManager::get().free( Basis_cmplx ); 
// 2c e 

    if( GDenS ) CQMemManager::get().free( GDenS );
    if( GDenZ ) CQMemManager::get().free( GDenZ );
    if( GDenY ) CQMemManager::get().free( GDenY );
    if( GDenX ) CQMemManager::get().free( GDenX );

    if( TS ) CQMemManager::get().free( TS );
    if( TZ ) CQMemManager::get().free( TZ );
    if( TY ) CQMemManager::get().free( TY );
    if( TX ) CQMemManager::get().free( TX );

    if( GTS ) CQMemManager::get().free( GTS );
    if( GTZ ) CQMemManager::get().free( GTZ );
    if( GTY ) CQMemManager::get().free( GTY );
    if( GTX ) CQMemManager::get().free( GTX );

    if( gPTss ) CQMemManager::get().free( gPTss );
    if( gPTsx ) CQMemManager::get().free( gPTsx );
    if( gPTsy ) CQMemManager::get().free( gPTsy );
    if( gPTsz ) CQMemManager::get().free( gPTsz );
    if( gPTzz ) CQMemManager::get().free( gPTzz );
    if( gPTyy ) CQMemManager::get().free( gPTyy );
    if( gPTxx ) CQMemManager::get().free( gPTxx );

    CQMemManager::get().free( eps, U_n, dVU_n, d2VU_n );
    if( U_gamma )      CQMemManager::get().free( U_gamma );
    if( dVU_gamma )    CQMemManager::get().free( dVU_gamma );
    if( d2VU_gamma )   CQMemManager::get().free( d2VU_gamma );
    if( d2VU_n_gamma ) CQMemManager::get().free( d2VU_n_gamma );


    if( eps_SCR )       CQMemManager::get().free( eps_SCR );
    if( dVU_n_SCR )     CQMemManager::get().free( dVU_n_SCR );
    if( dVU_gamma_SCR ) CQMemManager::get().free( dVU_gamma_SCR );
    if( d2VU_n_SCR )     CQMemManager::get().free( d2VU_n_SCR );
    if( d2VU_gamma_SCR ) CQMemManager::get().free( d2VU_gamma_SCR );
    if( d2VU_n_gamma_SCR ) CQMemManager::get().free( d2VU_n_gamma_SCR );

    CQMemManager::get().free( ZrhoVar, ZMAT );
    if( isGGA ) {
      CQMemManager::get().free( ZgammaVar1, ZgammaVar2, ZgammaVar3, ZgammaVar4 );
    }

#ifdef CQ_ENABLE_MPI
    MPICommFree(intComm); // Free communicator

    } // End of the MPI

    MPI_Barrier(c);
#endif

    // Turn back on LA threads
    SetLAThreads(LAThreads);

  } // KohnSham<dcomplex, dcomplex>::formFXC

}; // namespace ChronusQ

