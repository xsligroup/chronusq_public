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

#include <particleintegrals/inhouseaointegral.hpp>

namespace ChronusQ {

  /**
   *  \brief Compute the complex ERI of two shell pairs
   *
   *
   *  \param [in] pair1  bra shell pair data for shell1,shell2
   *  \param [in] pair2  ket shell pair data for shell3,shell4
   *  \param [in] shell1
   *  \param [in] shell2
   *  \param [in] shell3
   *  \param [in] shell4
   *
   *  \return complex ERI of two shell pairs ( shell1 shell2 | shell3 shell4 )
   */

   
  std::vector<dcomplex> ComplexGIAOIntEngine::computeGIAOERIabcd(
    libint2::ShellPair &pair1 , libint2::ShellPair &pair2, 
    libint2::Shell &shell1, libint2::Shell &shell2,
    libint2::Shell &shell3, libint2::Shell &shell4, double *H )  {
    
    dcomplex tmpVal=0.0,sqrPQ,PQ;
    std::vector<dcomplex> ERI_cart;
  
/* 
    if ( std::abs(H[0])<1.0e-15 )
      H[0] = 0.0;
  
    if ( std::abs(H[1])<1.0e-15 )
      H[1] = 0.0;
    if ( std::abs(H[2])<1.0e-15 )
      H[2] = 0.0;
*/

//    std::cout<<"H= "<<H[0]<<" "<<H[1]<<" "<<H[2]<<std::endl;


    int lA[3],lB[3],lC[3],lD[3];

    // compute total angular momentum
    auto lTotal = shell1.contr[0].l + shell2.contr[0].l
                + shell3.contr[0].l + shell4.contr[0].l;

/*
    std::cerr<<"LA "<<shell1.contr[0].l
    <<" LB "<<shell2.contr[0].l
    <<" LC "<<shell3.contr[0].l 
    <<" LD "<<shell4.contr[0].l<<std::endl;
*/

    double ka[3],kb[3],kc[3],kd[3],k1[3],k2[3];
    // here calculate the phase factor in LONDON orbital

    ka[0] = - 0.5*( shell1.O[1]*H[2] - shell1.O[2]*H[1] );
    ka[1] = - 0.5*( shell1.O[2]*H[0] - shell1.O[0]*H[2] );
    ka[2] = - 0.5*( shell1.O[0]*H[1] - shell1.O[1]*H[0] );

    kb[0] = 0.5*( shell2.O[1]*H[2] - shell2.O[2]*H[1] );
    kb[1] = 0.5*( shell2.O[2]*H[0] - shell2.O[0]*H[2] );
    kb[2] = 0.5*( shell2.O[0]*H[1] - shell2.O[1]*H[0] );

    kc[0] = - 0.5*( shell3.O[1]*H[2] - shell3.O[2]*H[1] );
    kc[1] = - 0.5*( shell3.O[2]*H[0] - shell3.O[0]*H[2] );
    kc[2] = - 0.5*( shell3.O[0]*H[1] - shell3.O[1]*H[0] );

    kd[0] = 0.5*( shell4.O[1]*H[2] - shell4.O[2]*H[1] );
    kd[1] = 0.5*( shell4.O[2]*H[0] - shell4.O[0]*H[2] );
    kd[2] = 0.5*( shell4.O[0]*H[1] - shell4.O[1]*H[0] );

    for ( int mu = 0 ; mu < 3 ; mu++ ) {
      k1[mu] = ka[mu] + kb[mu];
      k2[mu] = kc[mu] + kd[mu];
    }

    // pre calculate ss type overlap integrals
    auto ss_shellpair1 = computecompOverlapss( pair1, shell1, ka, shell2, kb );
    auto ss_shellpair2 = computecompOverlapss( pair2, shell3, kc, shell4, kd );

    // pre calculate all the Boys functions 
    // dimension is FmT_2e[shellpair1.prim][shellpair2.prim][lTotal+1]
    std::vector<std::vector<std::vector<dcomplex>>> FmT_2e;
    FmT_2e.resize(pair1.primpairs.size());

    dcomplex *FmT = new dcomplex[lTotal+1];
    // double *nFmT = new double[lTotal+1];
    // double nPQ,nsqrPQ;
    dcomplex onei;
    onei.real(0);
    onei.imag(1);
    int shellpair1_i=0, shellpair2_j ;
    for ( auto &pripair1 : pair1.primpairs ) {
      FmT_2e[shellpair1_i].resize( pair2.primpairs.size() );

      shellpair2_j = 0 ; 
      for ( auto &pripair2 : pair2.primpairs ) {
        sqrPQ = 0.0;
        // nsqrPQ = 0.0;
        for ( int mu=0 ; mu<3 ; mu++ ) {
          PQ = ( pripair1.P[mu] + 0.5*k1[mu]*pripair1.one_over_gamma*onei 
              -( pripair2.P[mu] + 0.5*k2[mu]*pripair2.one_over_gamma*onei ) ); 
          sqrPQ += PQ*PQ;
          
          // nPQ = pripair1.P[mu]-pripair2.P[mu];
          // nsqrPQ += nPQ * nPQ;
        }
        auto Zeta = 1.0/pripair1.one_over_gamma;
        auto Eta  = 1.0/pripair2.one_over_gamma;
        
        auto rho = Zeta*Eta/(Zeta+Eta);
        auto T = rho*sqrPQ;
        // auto nT = rho * nsqrPQ;
        // calculate Fm(T) list
        computecompFmT( FmT, T, lTotal, 0 );
        // computeFmTTaylor(nFmT,nT,lTotal,0);

//std::cout<<"start shellpair FmT "<<std::endl;
        for ( int lcurr = 0 ; lcurr < lTotal+1 ; lcurr++ ) {

           // if ( std::abs(nFmT[lcurr] - FmT[lcurr]) > 1.0e-10 )
           //  std::cout<<"FmT doesn't match"<<std::endl;


           if ( std::abs(FmT[lcurr]) < 1.0e-13 ) 
             FmT_2e[shellpair1_i][shellpair2_j].push_back(0.0);
           else
             FmT_2e[shellpair1_i][shellpair2_j].push_back(FmT[lcurr]);

//std::cout<<" FmT_2e[][] [] = "<<FmT_2e[shellpair1_i][shellpair2_j][lcurr]<<std::endl;
//std::cout<<"FmT complex ";
//std::cout<<std::setprecision(12)<<FmT[lcurr]<<std::endl;


        } // for lcurr
        shellpair2_j ++;
      } // for pripair2
    shellpair1_i++;
    } // for pripair1
    delete[] FmT;
    // delete[] nFmT;



    for(int i = 0; i < cart_ang_list[shell1.contr[0].l].size(); i++) 
    for(int j = 0; j < cart_ang_list[shell2.contr[0].l].size(); j++)
    for(int k = 0; k < cart_ang_list[shell3.contr[0].l].size(); k++)
    for(int l = 0; l < cart_ang_list[shell4.contr[0].l].size(); l++) {
      for (int mu=0 ; mu<3 ; mu++) {
        lA[mu] = cart_ang_list[shell1.contr[0].l][i][mu];
        lB[mu] = cart_ang_list[shell2.contr[0].l][j][mu];
        lC[mu] = cart_ang_list[shell3.contr[0].l][k][mu];
        lD[mu] = cart_ang_list[shell4.contr[0].l][l][mu];
      }  // for mu

      
      tmpVal = twoecomphRRabcd(pair1,pair2,shell1,shell2,shell3,shell4,
                 FmT_2e,k1,k2,ss_shellpair1,ss_shellpair2,shell1.contr[0].l, lA, 
                 shell2.contr[0].l, lB, shell3.contr[0].l, lC, shell4.contr[0].l, lD ); 
      
      ERI_cart.push_back(tmpVal);

    }   // for l


    if ( ( not shell1.contr[0].pure ) and ( not shell2.contr[0].pure ) and 
         ( not shell3.contr[0].pure ) and ( not shell4.contr[0].pure ) ) {  
      // if both sides are cartesian, return cartesian gaussian integrals
      return ERI_cart;
    }

    std::vector<dcomplex> ERI_sph;

    ERI_sph.assign(((2*shell1.contr[0].l+1)*(2*shell2.contr[0].l+1)
                   *(2*shell3.contr[0].l+1)*(2*shell4.contr[0].l+1)),0.0); 

    cart2sph_complex_2e_transform( shell1.contr[0].l,shell2.contr[0].l,
      shell3.contr[0].l,shell4.contr[0].l,ERI_sph,ERI_cart );

    return ERI_sph; 

  }   // computeGIAOERIabcd

  std::vector<dcomplex> ComplexGIAOIntEngine::bottomupcomplexERI(libint2::ShellPair &pair1 ,
    libint2::ShellPair &pair2, libint2::Shell &shell1, libint2::Shell &shell2,
    libint2::Shell &shell3, libint2::Shell &shell4, double *H) {


    // here calculate the phase factor in LONDON orbital
    double ka[3],kb[3],kc[3],kd[3],k1[3],k2[3];
    // here calculate the phase factor in LONDON orbital

    ka[0] = - 0.5*( shell1.O[1]*H[2] - shell1.O[2]*H[1] );
    ka[1] = - 0.5*( shell1.O[2]*H[0] - shell1.O[0]*H[2] );
    ka[2] = - 0.5*( shell1.O[0]*H[1] - shell1.O[1]*H[0] );

    kb[0] = 0.5*( shell2.O[1]*H[2] - shell2.O[2]*H[1] );
    kb[1] = 0.5*( shell2.O[2]*H[0] - shell2.O[0]*H[2] );
    kb[2] = 0.5*( shell2.O[0]*H[1] - shell2.O[1]*H[0] );

    kc[0] = - 0.5*( shell3.O[1]*H[2] - shell3.O[2]*H[1] );
    kc[1] = - 0.5*( shell3.O[2]*H[0] - shell3.O[0]*H[2] );
    kc[2] = - 0.5*( shell3.O[0]*H[1] - shell3.O[1]*H[0] );

    kd[0] = 0.5*( shell4.O[1]*H[2] - shell4.O[2]*H[1] );
    kd[1] = 0.5*( shell4.O[2]*H[0] - shell4.O[0]*H[2] );
    kd[2] = 0.5*( shell4.O[0]*H[1] - shell4.O[1]*H[0] );

    for ( int mu = 0 ; mu < 3 ; mu++ ) {
      k1[mu] = ka[mu] + kb[mu];
      k2[mu] = kc[mu] + kd[mu];
    }

    dcomplex onei;
    onei.real(0.0);
    onei.imag(1.0);
    

    int LA, LB, LC, LD, Lbra, Lket, Ltot;
    
    LA = shell1.contr[0].l; 
    LB = shell2.contr[0].l; 
    LC = shell3.contr[0].l; 
    LD = shell4.contr[0].l; 
    
    Lbra = LA + LB;
    Lket = LC + LD;
    Ltot = Lbra + Lket;

    double A[3], B[3], C[3], D[3], AB[3], CD[3];
    
    for (int ii = 0 ; ii < 3 ; ii++ ) { 
      A[ii] = shell1.O[ii];
      B[ii] = shell2.O[ii];
      C[ii] = shell3.O[ii];
      D[ii] = shell4.O[ii];
      AB[ii] = A[ii] - B[ii];
      CD[ii] = C[ii] - D[ii];
    }  
    
    dcomplex *FT = new dcomplex[Ltot+1];   

   
    // allocate memory for horizontal recursion
    std::vector<std::vector<std::vector<dcomplex>>> Vbraketee;
    int counter = 0; 

    Vbraketee.resize((Lbra-LA+1)*(Lket-LC+1));
    for ( int LAidx = LA ; LAidx <= Lbra ; LAidx++ ) {   
      
      for ( int LCidx = LC ; LCidx <= Lket ; LCidx++ ) {   
     // should loop over C index here instead of B for ( int LBidx = 0 ; LBidx <= Lbra - LAidx ; LBidx++ ) {  
        Vbraketee[(LAidx-LA)*(Lket-LC+1)+LCidx-LC].resize((Lbra-LAidx+1)*(Lket-LCidx+1));
     //   Vbraketee[(LAidx-LA)*(LBidx+1)+LBidx].resize((Lket-LC+1)*(LD+1));
        for ( int LBidx = 0 ; LBidx <= Lbra-LAidx ; LBidx++ ) {
          for ( int LDidx = 0 ; LDidx <= Lket-LCidx ; LDidx++ ) {
            Vbraketee[(LAidx-LA)*(Lket-LC+1)+LCidx-LC][LBidx*(Lket-LCidx+1)+LDidx].assign(
              cart_ang_list[LAidx].size()*cart_ang_list[LBidx].size()
             *cart_ang_list[LCidx].size()*cart_ang_list[LDidx].size(),0.0);
          } // for ( int LDidx = 0 ; LDidx <= LD ; LDidx++ )
        } // for ( int LBidx = 0 ; LBidx <= Lbra-LAidx ; LBidx++ )
      } // for ( int LCidx = LC ; LCidx <= Lket ; LCidx++ ) 
    } // for ( int LAidx = LA ; LAidx <= Lbra ; LAidx++ )      



    for ( auto &pripair1 : pair1.primpairs ) {
      double P[3];
      for (int ii = 0 ; ii < 3 ; ii++ ) {
        P[ii] = pripair1.P[ii]; 
      }

      // GIAO specific 
      dcomplex Ptilde[3];
      for ( int ii = 0 ; ii < 3 ; ii++ ) {
        Ptilde[ii] = P[ii]+0.5*k1[ii]*pripair1.one_over_gamma*onei;
      } 

      double zeta = 1.0/pripair1.one_over_gamma; 

      // here calculate primitive ss1 
      double realpart1 = 0.0; 
      for ( int mu = 0 ; mu < 3 ; mu++ ) { 
        realpart1 -= pow( (ka[mu]+kb[mu]), 2 );
      }
      
      realpart1 *= 0.25*pripair1.one_over_gamma;
      
      double imagpart1 = 0.0 ; 
      for ( int mu = 0 ; mu < 3 ; mu++ ) {
        imagpart1 += ka[mu]*(P[mu] - A[mu]) + kb[mu]*(P[mu]-B[mu]);
      }

      dcomplex z1 = realpart1 + imagpart1*onei; 

      dcomplex ss1 = shell1.contr[0].coeff[pripair1.p1]* shell2.contr[0].coeff[pripair1.p2]
        //* pow(sqrt(M_PI),3) * sqrt(pripair1.one_over_gamma)*pripair1.K ; 
        * exp(z1) * pow(sqrt(M_PI),3) * sqrt(pripair1.one_over_gamma)*pripair1.K ; 
      // calculation of primitive ss1 end 

      for ( auto &pripair2 : pair2.primpairs ) {

        double Q[3]; 
        for (int ii = 0 ; ii < 3 ; ii++ ) {
          Q[ii] = pripair2.P[ii]; 
        }

        dcomplex Qtilde[3];
        for ( int ii = 0 ; ii < 3 ; ii++ ) {
          Qtilde[ii] = Q[ii]+0.5*k2[ii]*pripair2.one_over_gamma*onei;
        } 
        
        // here calculate primitive ss2 
        double realpart2 = 0.0; 
        for ( int mu = 0 ; mu < 3 ; mu++ ) { 
          realpart2 -= pow( (kc[mu]+kd[mu]), 2 );
        }
        
        realpart2 *= 0.25*pripair2.one_over_gamma;
        
        double imagpart2 = 0.0 ; 
        for ( int mu = 0 ; mu < 3 ; mu++ ) {
          imagpart2 += kc[mu]*(Q[mu] - C[mu]) + kd[mu]*(Q[mu]-D[mu]);
        }
       
        dcomplex z2 = realpart2 + imagpart2*onei; 
       
        dcomplex ss2 = shell3.contr[0].coeff[pripair2.p1]* shell4.contr[0].coeff[pripair2.p2]
          * exp(z2) * pow(sqrt(M_PI),3) * sqrt(pripair2.one_over_gamma)*pripair2.K ; 
          //* pow(sqrt(M_PI),3) * sqrt(pripair2.one_over_gamma)*pripair2.K ; 
        // calculation of primitive ss1 end 

        double eta = 1.0/pripair2.one_over_gamma; 
        double zetaG = zeta + eta;
        double rho = zeta*eta/zetaG;

        dcomplex T = 0.0;
        //dcomplex sqrPQ = 0.0; 
        for (int ii = 0 ; ii < 3 ; ii++) {
          T += pow( P[ii]-Q[ii] +0.5*k1[ii]*pripair1.one_over_gamma*onei
                    - 0.5*k2[ii]*pripair2.one_over_gamma*onei,2); 
        }
        T*= rho ; 
        
        //double ss2 = shell3.contr[0].coeff[pripair2.p1]* shell4.contr[0].coeff[pripair2.p2]
        //  * pow(sqrt(M_PI),3) * sqrt(pripair2.one_over_gamma)*pripair2.K ; 
        
        computecompFmT(FT,T,Ltot,0);        

        double W[3];
        for (int ii = 0 ; ii < 3 ; ii++ )
          W[ii] = (zeta*P[ii] + eta*Q[ii])/zetaG;         

        dcomplex Wtilde[3];
        for ( int ii = 0 ; ii < 3 ; ii++ ) { 
          Wtilde[ii] = ( Ptilde[ii] / pripair1.one_over_gamma + 
               Qtilde[ii] / pripair2.one_over_gamma )/zetaG;  
        }

        dcomplex pref = 2*sqrt(rho/M_PI)*ss1*ss2;
        
        // now allocate the memory to store vertical recursion elements
        std::vector<std::vector<std::vector<dcomplex>>> Vtempbraket((Lbra+1)*(Lket+1)); 
        for ( int k = 0 ; k <= Lbra ; k++ ){ 
          for ( int l = 0 ; l <= Lket ; l++ ) { 
            Vtempbraket[k*(Lket+1)+l].resize(cart_ang_list[k].size()*cart_ang_list[l].size());
            for ( int cart_i = 0; cart_i < cart_ang_list[k].size() ; cart_i++ ) {
              for ( int cart_j = 0 ; cart_j < cart_ang_list[l].size(); cart_j++ ) {
                Vtempbraket[k*(Lket+1)+l][cart_i*cart_ang_list[l].size()+cart_j].resize(Ltot-k-l+1);
              }
            }  // for ( int cart_i = 0; cart_i < cart_ang_list[k].size()
          } // for ( int l = 0 ; l <= Lket ; l++ ) 
        } //for ( int k = 0 ; k <= Lbra ; k++ )

        // copy the (ss|ss)^m integrals to the Vtempbraket 
        for ( int ii = 0 ; ii < Ltot+1 ; ii++ ) {
          Vtempbraket[0][0][ii] = FT[ii]*pref ; 
//std::cout<<"Vtempbraket= "<<std::setprecision(12)<<Vtempbraket[0][0][ii]<<std::endl;
        }  
        
        // vertical recursion
        for ( int k = 0 ; k <= Lbra ; k++ ){ 
          int l = 0;
          for ( int cart_i = 0; cart_i < cart_ang_list[k].size() ; cart_i++ ) {
            int lA_xyz[3];
            for ( int ii = 0 ; ii < 3 ; ii++ )
              lA_xyz[ii] = cart_ang_list[k][cart_i][ii];
            
            int mbra = Ltot - k ; 
            int iWork; 
            if ( k > 0 ) {
              if ( lA_xyz[0]>0 )  {
                iWork = 0;
              } else if ( lA_xyz[1]>0 ) {
                iWork = 1;
              } else if ( lA_xyz[2]>0 ) {
                iWork = 2;
              }
              
              // calculate the index of the element with angular momentum lower by 1 
              int lAtemp[3];
              for ( int ii = 0 ; ii < 3 ; ii++ )
                lAtemp[ii] = lA_xyz[ii];

              lAtemp[iWork] = lA_xyz[iWork]-1;
              
              int indexlm1;
              indexlm1 = indexmap(k-1,lAtemp[0],lAtemp[1],lAtemp[2]);
              
              for ( int m = 0 ; m<= mbra ; m++ ) {
                dcomplex ERIscratch = 0.0;
                ERIscratch = (Ptilde[iWork]-A[iWork])* Vtempbraket[(k-1)*(Lket+1)][indexlm1*cart_ang_list[l].size()][m]
                  +(Wtilde[iWork]-Ptilde[iWork])*Vtempbraket[(k-1)*(Lket+1)][indexlm1*cart_ang_list[l].size()][m+1];

                if ( lA_xyz[iWork]>1 ) {
                  lAtemp[iWork] = lA_xyz[iWork]-2;
                  int indexlm2 = indexmap(k-2,lAtemp[0],lAtemp[1],lAtemp[2]);
                  ERIscratch += 1/(2*zeta) * (lA_xyz[iWork]-1) * (
                    Vtempbraket[(k-2)*(Lket+1)][indexlm2*cart_ang_list[l].size()][m]
                    -rho/zeta *Vtempbraket[(k-2)*(Lket+1)][indexlm2*cart_ang_list[l].size()][m+1]);
                } // if ( lA_xyz[iWork]>1 )

                Vtempbraket[k*(Lket+1)][cart_i*cart_ang_list[l].size()][m] = ERIscratch;
              } // for ( int m = 0 ; m<= mbra ; m++ )
            } // if ( k > 0 )  
            
            for ( int l = 0 ; l < Lket+1 ; l++ ) {
              for ( int cart_j = 0 ; cart_j < cart_ang_list[l].size() ; cart_j++ ) {
                int lC_xyz[3];
                for ( int ii = 0 ; ii < 3 ; ii++ ) { 
                  lC_xyz[ii] = cart_ang_list[l][cart_j][ii];  
                }  
                int mbraket = Ltot-k-l;
                if (l>0) { 
                  if ( lC_xyz[0]>0 )  {
                    iWork = 0;
                  } else if ( lC_xyz[1]>0 ) {
                    iWork = 1;
                  } else if ( lC_xyz[2]>0 ) {
                    iWork = 2;
                  }
                  // calculate the index of the element with angular momentum lower by 1 
                  int lCtemp[3];
                  for ( int ii = 0 ; ii < 3 ; ii++ )
                    lCtemp[ii] = lC_xyz[ii];
               
                  lCtemp[iWork] = lC_xyz[iWork]-1;
                  
                  int indexlm1;
                  indexlm1 = indexmap(l-1,lCtemp[0],lCtemp[1],lCtemp[2]);
                 
                  for ( int m = 0 ; m <= mbraket ; m++ ) {
                    dcomplex ERIscratch = 0.0 ; 
                    ERIscratch += (Qtilde[iWork]-C[iWork])*Vtempbraket[(k*(Lket+1)+l-1)][cart_i
                      *cart_ang_list[l-1].size()+indexlm1][m] 
                      +(Wtilde[iWork]-Qtilde[iWork])* Vtempbraket[(k*(Lket+1)+l-1)][cart_i*
                      cart_ang_list[l-1].size()+indexlm1][m+1]; 

                    if (lC_xyz[iWork]>1) {
                      lCtemp[iWork] = lC_xyz[iWork]-2;
                      int indexlm2 = indexmap(l-2,lCtemp[0],lCtemp[1],lCtemp[2]);
                      
                      ERIscratch += 1/(2*eta)*(lC_xyz[iWork]-1)*(Vtempbraket[k*(Lket+1)+l-2]
                        [cart_i*cart_ang_list[l-2].size()+indexlm2][m]-rho/eta*Vtempbraket
                        [k*(Lket+1)+l-2][cart_i*cart_ang_list[l-2].size()+indexlm2][m+1]);
                    } // if (lC_xyz[iWork]>1) 
                    
                    if ( lA_xyz[iWork] > 0 ) {
                      int lAtemp[3];
                      for ( int ii = 0 ; ii < 3 ; ii++ )
                        lAtemp[ii] = lA_xyz[ii];
                     
                      lAtemp[iWork] = lA_xyz[iWork]-1;
                      int indexlAm1 = indexmap(k-1,lAtemp[0],lAtemp[1],lAtemp[2]);
                      ERIscratch += 1/(2*zetaG)*lA_xyz[iWork]*Vtempbraket[(k-1)*(Lket+1)+l-1]
                        [indexlAm1*cart_ang_list[l-1].size()+indexlm1][m+1];
                    } // if ( lA_xyz[iWork] > 0 ) 
                    Vtempbraket[k*(Lket+1)+l][cart_i*cart_ang_list[l].size()+cart_j][m] = ERIscratch; 

                  } // for ( int m = 0 ; m <= mbraket ; m++ )  
                } // if (l>0)
                                
              } // for ( int cart_j = 0 ; cart_j < cart_ang_list[l].size()
            } // for ( int l = 0 ; l < Lket+1 ; l++ )   


          } // for ( int cart_i 
        } // for ( int k = 0 ; k <= Lbra ; k++ )
        // should be some contractions 


           
        // copy the elements from Vtempbraket to Vbraketee  
           
        for ( int LAidx = LA ; LAidx <= Lbra ; LAidx++ ) {
          for ( int LCidx = LC ; LCidx <= Lket ; LCidx++ ) {
            // loop over dimension 
            for ( int cart_i = 0 ; cart_i < cart_ang_list[LAidx].size() ; cart_i++ ) {
              for ( int cart_j = 0 ; cart_j < cart_ang_list[LCidx].size() ; cart_j++ ) {
                Vbraketee[(LAidx-LA)*(Lket-LC+1)+LCidx-LC][0][cart_i*1
                  *cart_ang_list[LCidx].size()*1+cart_j*1] += 
                  Vtempbraket[LAidx*(Lket+1)+LCidx][cart_i*cart_ang_list[LCidx].size()+cart_j][0];
              } // for cart_j
            } // for cart_i
          } // for LCidx 
        } // for LAidx        

        

      } // for ( auto &pripair2 : pair1.primpairs )
    } // for ( auto &pripair1 : pair1.primpairs )  
    delete[] FT;

    int iWork;  
    // horizontal recursion 
    if (LB >0) {
      for ( int lB = 1 ; lB <= LB ; lB++ ) { // which implies that we skip lB=0
        for ( int lA = LA ; lA <= Lbra-lB ; lA++ ) {
          for  ( int Aidx = 0 ; Aidx < cart_ang_list[lA].size(); Aidx++ ) {
            
            int lA_xyz[3];
            for ( int ii = 0 ; ii < 3 ; ii++ ) { 
              lA_xyz[ii] = cart_ang_list[lA][Aidx][ii];
            }
            
            // here loop over elements in lB 
            for ( int Bidx = 0 ; Bidx<cart_ang_list[lB].size();Bidx++ ) {
            
              int lB_xyz[3];
              for ( int ii = 0 ; ii < 3 ; ii++ ) { 
                lB_xyz[ii] = cart_ang_list[lB][Bidx][ii];
              }
              
              int iWork;
              if  (lB_xyz[0]>0) {
                iWork = 0;
              } else if (lB_xyz[1]>0) {
                iWork = 1;
              } else if (lB_xyz[2]>0) {
                iWork = 2 ; 
              }
              
              int lBm1[3] ; 
              for ( int ii = 0 ; ii < 3 ; ii++ ) {
                lBm1[ii] = lB_xyz[ii];
              } 
              lBm1[iWork] = lB_xyz[iWork]-1;
                 
              int lAp1[3] ; 
              for ( int ii = 0 ; ii < 3 ; ii++ ) {
                lAp1[ii] = lA_xyz[ii];
              } 
              lAp1[iWork] = lA_xyz[iWork]+1;

              int idxBtemp = indexmap(lB-1,lBm1[0],lBm1[1],lBm1[2]); 
              int idxAp1temp = indexmap(lA+1,lAp1[0],lAp1[1],lAp1[2]);
              
              for ( int lC = LC ; lC<= Lket ; lC++ ) {
                for ( int Cidx = 0 ; Cidx < cart_ang_list[lC].size() ; Cidx++) {
                  int lC_xyz[3];
                  for ( int ii = 0 ; ii < 3 ; ii++ ) {
                    lC_xyz[ii] = cart_ang_list[lC][Cidx][ii];
                  }
                  
                  Vbraketee[(lA-LA)*(Lket-LC+1)+lC-LC][lB*(Lket-lC+1)][Aidx*cart_ang_list[lB].size()
                    *cart_ang_list[lC].size()*1+Bidx*cart_ang_list[lC].size()+Cidx*1]=
                    Vbraketee[(lA+1-LA)*(Lket-LC+1)+lC-LC][(lB-1)*(Lket-lC+1)][idxAp1temp*cart_ang_list[lB-1].size() 
                    *cart_ang_list[lC].size()*1 +idxBtemp *cart_ang_list[lC].size() *1 +Cidx]
                    + AB[iWork]*Vbraketee[(lA-LA)*(Lket-LC+1)+lC-LC][(lB-1)*(Lket-lC+1)][
                    Aidx*cart_ang_list[lB-1].size()*cart_ang_list[lC].size()*1
                    +idxBtemp*cart_ang_list[lC].size()*1+Cidx];
                    // (ab|c0)= (a+1b-1|c0)+AB(ab-1|c0)  
                } // for ( int Cidx = 0 ; Cidx < cart_ang_list[lC].size() 
              } // for ( int lC = LC ; lC<= Lket ; lC++ )
                 
            } //for ( int Bidx = 0 ; Bidx<cart_ang_list[lB].size();Bidx++ )  

          } // for  ( int Aidx = 0 ; Aidx < cart_ang_list[lA].size(); Aidx++ )   
        }  // for ( int lA = LA ; lA <= Lbra-lB ; lA++ ) 
      } //  for ( int lB = 1 ; lB <= LB ; lB++ )

    } // if  (LB >0)  

    if (LD > 0 ) {
      for (int lD = 1 ; lD<=LD ; lD++ ) {
        for ( int lC = LC ; lC<= Lket - lD ;lC++ ) { 
          for ( int Cidx = 0 ; Cidx < cart_ang_list[lC].size() ; Cidx++ ) {
           
            int lC_xyz[3];
            for ( int ii = 0 ; ii < 3 ; ii++ ) {
              lC_xyz[ii] = cart_ang_list[lC][Cidx][ii];
            }
 
            for ( int Didx = 0 ; Didx<cart_ang_list[lD].size() ; Didx++ ) { 

              int lCp1[3];
              for ( int ii = 0 ; ii < 3 ; ii++ ) {
              //  lC_xyz[ii] = cart_ang_list[lC][Cidx][ii];
                lCp1[ii] = lC_xyz[ii];
              } 


              int lD_xyz[3],lDm1[3];
              for ( int ii = 0 ; ii < 3 ; ii++ ) {
                lD_xyz[ii] = cart_ang_list[lD][Didx][ii];
                lDm1[ii] = lD_xyz[ii]; 
              } 
              if (lD_xyz[0]>0) {
                iWork = 0 ; 
              } else if (lD_xyz[1]>0) {
                iWork = 1 ;
              } else if (lD_xyz[2]>0) {
                iWork = 2;
              } 
              
              lDm1[iWork] = lD_xyz[iWork]-1; 
              lCp1[iWork] = lC_xyz[iWork]+1;
               
              int idxCtemp = indexmap(lC+1,lCp1[0],lCp1[1],lCp1[2]);
              int idxDtemp = indexmap(lD-1,lDm1[0],lDm1[1],lDm1[2]);
              
              for ( int Aidx = 0 ; Aidx < cart_ang_list[LA].size() ; Aidx++ ) {

                int lA_xyz[3] ; 
                for ( int ii = 0 ; ii < 3 ; ii++ ) {
                  lA_xyz[ii] =cart_ang_list[LA][Aidx][ii]; 
                }
                
                for ( int Bidx = 0 ; Bidx < cart_ang_list[LB].size() ; Bidx++ ) {
                  int lB_xyz[3];
                  for (int ii = 0 ; ii < 3 ; ii++) {
                    lB_xyz[ii] = cart_ang_list[LB][Bidx][ii];
                  }
                  
                  Vbraketee[(LA-LA)*(Lket-LC+1)+lC-LC][LB*(Lket-lC+1)+lD][Aidx *cart_ang_list[LB].size()
                    *cart_ang_list[lC].size()*cart_ang_list[lD].size()+ Bidx *cart_ang_list[lC].size()
                    *cart_ang_list[lD].size() + Cidx*cart_ang_list[lD].size() + Didx]  
                  = Vbraketee[(LA-LA)*(Lket-LC+1)+lC+1-LC][LB*(Lket-(lC+1)+1)+lD-1][Aidx*cart_ang_list[LB].size()
                    *cart_ang_list[lC+1].size()*cart_ang_list[lD-1].size()+Bidx*cart_ang_list[lC+1].size()*cart_ang_list[lD-1].size()
                    +idxCtemp*cart_ang_list[lD-1].size()+idxDtemp]
                    +CD[iWork]*Vbraketee[(LA-LA)*(Lket-LC+1)+lC-LC][LB*(Lket-lC+1)+lD-1][
                    Aidx*cart_ang_list[LB].size()*cart_ang_list[lC].size()*cart_ang_list[lD-1].size()
                    + Bidx*cart_ang_list[lC].size()*cart_ang_list[lD-1].size() +
                    Cidx*cart_ang_list[lD-1].size()+idxDtemp];
                 

                  // (ab|cd) = (ab|c+1d-1)+(CD)(ab|cd-1)
                } // for ( int Bidx = 0 ; Bidx < cart_ang_list[LB].size() ; Bidx++ )
              } // for ( int Aidx = 0 ; Aidx < cart_ang_list(LA) ; Aidx++ ) 
              

  
            } // for ( int Didx = 0 ; Didx<cart_ang_list[lD].size() ; Dix++ )
            //lAp1[iWork] = lA_xyz[iWork]+1;
             
          } // for ( int Cidx = 0 ; Cidx < cart_ang_list[lC].size() ; Cidx++ )
        } // for ( int lC = LC ; lC<= Lket - lD ;lC++ )
      }// for (int lD = 1 ; lD<=LD ; lD++ ) 
    } // if (LD > 0 ) 
    
    // here Vbraketee[0][LB*(LD+1)+LD] is the desired ints 

    if ( ( not shell1.contr[0].pure ) and ( not shell2.contr[0].pure ) and 
         ( not shell3.contr[0].pure ) and ( not shell4.contr[0].pure ) ) {  
      // if both sides are cartesian, return cartesian gaussian integrals
      return Vbraketee[0][LB*(LD+1)+LD];
    } 

    std::vector<dcomplex> ERI_sph;

    ERI_sph.assign(((2*shell1.contr[0].l+1)*(2*shell2.contr[0].l+1)
                   *(2*shell3.contr[0].l+1)*(2*shell4.contr[0].l+1)),0.0); 

    cart2sph_complex_2e_transform( shell1.contr[0].l,shell2.contr[0].l,
      shell3.contr[0].l,shell4.contr[0].l,ERI_sph,Vbraketee[0][LB*(LD+1)+LD] );

    return ERI_sph; 


  } //ComplexGIAOIntEngine::bottomupcomplexERI

  /**
   *  \brief horizontal recursion of complex ERI when all angular momentum are nonzero
   *
   *
   *  \param [in] pair1  bra shell pair data for shell1,shell2
   *  \param [in] pair2  ket shell pair data for shell3,shell4
   *  \param [in] shell1
   *  \param [in] shell2
   *  \param [in] shell3
   *  \param [in] shell4
   *  \param [in] FmT_2e Boys function between two shell pairs
   *  \param [in] LA
   *  \param [in] lA
   *  \param [in] LB
   *  \param [in] lB
   *  \param [in] LC
   *  \param [in] lC
   *  \param [in] LD
   *  \param [in] lD
   *
   *  \return complex ERI of two shell pairs ( shell1 shell2 | shell3 shell4 )
   */

  //----------------------------------------------------//
  // two-e horizontal recursion from (ab|cd) to (a0|cd) //
  // (ab|cd)=(a+1,b-1|cd)+(A-B)*(a,b-1|cd)              //
  //----------------------------------------------------//
   
  dcomplex ComplexGIAOIntEngine::twoecomphRRabcd(
    libint2::ShellPair &pair1 ,libint2::ShellPair &pair2 ,
    libint2::Shell &shell1, libint2::Shell &shell2,
    libint2::Shell &shell3, libint2::Shell &shell4,
    std::vector<std::vector<std::vector<dcomplex>>> &FmT_2e, double *k1, double *k2, 
    std::vector<dcomplex> &ss_shellpair1, std::vector<dcomplex> & ss_shellpair2, 
    int LA,int *lA,int LB,int *lB,int LC,int *lC,int LD,int *lD) {

     dcomplex tmpVal = 0.0, tmpVal1=0.0;

  // iWork is used to indicate which Cartesian angular momentum we are reducing (x,y,z)

    int iWork;
    int totalL = LA + LB + LC + LD;


    if(totalL==0) {   // (SS||SS) type 

      return twoecompSSSS0(pair1,pair2,shell1,shell2,shell3,shell4, FmT_2e,ss_shellpair1,ss_shellpair2 );

      } else { // when totalL! = 0

        if( LB>=1 ) {

          int lAp1[3],lBm1[3];
          for( iWork=0 ; iWork<3 ; iWork++ ){
            lAp1[iWork]=lA[iWork];     
            lBm1[iWork]=lB[iWork];     
          } // for iWork

          if (lB[0]>0)      iWork=0;   
          else if (lB[1]>0) iWork=1;
          else if (lB[2]>0) iWork=2;
          lAp1[iWork]+=1;
          lBm1[iWork]-=1;

          tmpVal += twoecomphRRabcd(pair1, pair2, shell1, shell2, shell3, shell4,
                      FmT_2e, k1, k2, ss_shellpair1, ss_shellpair2,
                      LA+1, lAp1, LB-1, lBm1, LC, lC, LD, lD);

          if ( std::abs(pair1.AB[iWork]) > 1.0e-15 ) {
            tmpVal += pair1.AB[iWork] * twoecomphRRabcd(pair1,pair2,shell1,shell2,shell3,
                        shell4,FmT_2e, k1, k2, ss_shellpair1, ss_shellpair2,
                        LA,lA,LB-1,lBm1,LC,lC,LD,lD);
          } // if ( std::abs(pair1.AB[iWork]) > 1.0e-15 )

        } else if ( LB == 0 ) {
          tmpVal = twoecomphRRa0cd(pair1,pair2,shell1,shell2,shell3,shell4,FmT_2e,k1,k2,
                                 ss_shellpair1, ss_shellpair2, LA,lA,LC,lC,LD,lD);

        } // LB == 0
      }  // else ( LTOTAL != 0 )

      return tmpVal;

  }  // twoecomphRRabcd



  /**
   *  \brief horizontal recursion of complex ERI when all LB=0
   *
   *
   *  \param [in] pair1  bra shell pair data for shell1,shell2
   *  \param [in] pair2  ket shell pair data for shell3,shell4
   *  \param [in] shell1
   *  \param [in] shell2
   *  \param [in] shell3
   *  \param [in] shell4
   *  \param [in] FmT_2e Boys function between two shell pairs
   *  \param [in] LA
   *  \param [in] lA
   *  \param [in] LC
   *  \param [in] lC
   *  \param [in] LD
   *  \param [in] lD
   *
   *  \return complex ERI of two shell pairs ( shell1 shell2 | shell3 shell4 )
   */
  //----------------------------------------------------//
  // two-e horizontal recursion from (a0|cd) to (a0|c0) //
  // (a0|cd)=(a,0|c+1,d-1)+(C-D)*(a,0|c,d-1)            //
  //----------------------------------------------------//
   
  dcomplex ComplexGIAOIntEngine::twoecomphRRa0cd(
    libint2::ShellPair &pair1, libint2::ShellPair &pair2,
    libint2::Shell &shell1, libint2::Shell &shell2,
    libint2::Shell &shell3, libint2::Shell &shell4, 
    std::vector<std::vector<std::vector<dcomplex>>> &FmT_2e, double *k1, double *k2, 
    std::vector<dcomplex> &ss_shellpair1, std::vector<dcomplex> &ss_shellpair2,
    int LA,int *lA,int LC,int *lC,int LD,int *lD)  {

    dcomplex tmpVal=0.0;

    if(LD==0) {
      int pair1index=0, pair2index=0;
      // go into the vertical recursion
      for ( auto &pripair1 : pair1.primpairs ) {
        pair2index = 0; 
        for ( auto &pripair2 : pair2.primpairs ) {
/*
          auto norm = 
                 shell1.contr[0].coeff[pripair1.p1]* 
                 shell2.contr[0].coeff[pripair1.p2]* 
                 shell3.contr[0].coeff[pripair2.p1]* 
                 shell4.contr[0].coeff[pripair2.p2];  
*/
          // norm is given in ss overlap
          tmpVal +=  twoecompvRRa0c0( pripair1, pripair2,  
             FmT_2e[pair1index][pair2index],k1,k2,ss_shellpair1[pair1index],
               ss_shellpair2[pair2index], shell1,shell3, 0, LA,lA,LC,lC );

          pair2index++;
        } // for pripair2
        pair1index++;
      } // for pripair1

    } else { // if LD>0, go into horizontal recursion 

      int iWork;
      int lCp1[3],lDm1[3];  
     
      for( int iWork=0 ; iWork<3 ; iWork++ ){
        lCp1[iWork]=lC[iWork];
        lDm1[iWork]=lD[iWork];
      }
     
      if (lD[0]>0) iWork=0;
      else if (lD[1]>0) iWork=1;
      else if (lD[2]>0) iWork=2;
     
      lCp1[iWork]+=1;


    // when LD > 0

      lDm1[iWork] -=1 ;
      tmpVal = twoecomphRRa0cd(pair1, pair2, shell1, shell2, shell3, shell4, FmT_2e, k1,k2,
                 ss_shellpair1,ss_shellpair2, LA,lA, LC+1,lCp1, LD-1,lDm1 );
 
      if ( std::abs(pair2.AB[iWork]) > 1.0e-15 ){
        tmpVal += pair2.AB[iWork] * twoecomphRRa0cd( pair1, pair2, shell1, shell2, 
               shell3, shell4, FmT_2e, k1,k2, ss_shellpair1,ss_shellpair2, 
               LA,lA, LC,lC, LD-1,lDm1 );
      }  

    } // else ( that means LD > 0 )
    return tmpVal;
  }  // twoehRRa0cd



  /**
   *  \brief vertical recursion of complex ERI when all LA, LC > 0
   *
   *
   *  \param [in] pripair1  primitive bra shell pair data for shell1,shell2
   *  \param [in] pripair2  primitive ket shell pair data for shell3,shell4
   *  \param [in] pair1index index of primitive shell pair among the contracted pair
   *  \param [in] pair2index index of primitive shell pair among the contracted pair
   *  \param [in] FmT_2e Boys function between two primitive shell pairs
   *  \param [in] shell1    
   *  \param [in] shell3    
   *  \param [in] m         order of auxiliary function
   *  \param [in] LA
   *  \param [in] lA
   *  \param [in] LC
   *  \param [in] lC
   *
   *  \return complex ERI of two primitive shell pairs ( shell1 shell2 | shell3 shell4 )
   */

//---------------------------------------------------------------//
// two-e vertical recursion from [a0|c0] to [a0|00]              //
// [a0|c0]^m = (Q-C)*[a0|c-1,0]^m                                //
//           + (W-Q)*[a0|c-1,0]^(m+1)                            //
//           + N(a)/(2*(zeta+eta))*[a-1,0|c-1,0]^(m+1)           //
//           + (N(c)-1)/(2*eta)*[a0|c-2,0]^m                     //
//           - (N(c)-1)/(2*eta)*zeta/(zeta+eta)*[a0|c-2,0]^(m+1) //
//---------------------------------------------------------------//
   
  dcomplex ComplexGIAOIntEngine::twoecompvRRa0c0(
    libint2::ShellPair::PrimPairData &pripair1,
    libint2::ShellPair::PrimPairData &pripair2, 
    std::vector<dcomplex> &FmT_2epri, double *k1, double *k2, 
    dcomplex ss_pair1, dcomplex ss_pair2, 
    libint2::Shell &shell1, libint2::Shell &shell3,
    int m, int LA, int *lA, int LC, int *lC ) {
 
 
    if(LC==0) return twoecompvRRa000( pripair1, pripair2, FmT_2epri, k1, k2,
                       ss_pair1, ss_pair2, shell1, m, LA, lA );
 
    int lAm1[3],lCm1[3];  
 
    for ( int iWork=0 ; iWork<3 ; iWork++ ){
      lAm1[iWork]=lA[iWork];     
      lCm1[iWork]=lC[iWork];
    }
 
    dcomplex tmpVal=0.0;
    dcomplex Ptilde_iwork, Qtilde_iwork, Wtilde_iwork;
    double W_iWork,Zeta;
    int iWork;
    dcomplex onei;
    onei.real(0);
    onei.imag(1);
 
    if (lC[0]>0) iWork=0;
    else if (lC[1]>0) iWork=1;
    else if (lC[2]>0) iWork=2;

    lCm1[iWork]-=1;

    Qtilde_iwork = pripair2.P[iWork] + 0.5*k2[iWork]*pripair2.one_over_gamma*onei;
    
    if ( std::abs( Qtilde_iwork - shell3.O[iWork] ) > 1.0e-15 ) {
      tmpVal += ( Qtilde_iwork - shell3.O[iWork] ) * 
                twoecompvRRa0c0( pripair1, pripair2, FmT_2epri, k1,k2,ss_pair1,ss_pair2, 
                  shell1, shell3, m, LA,lA, LC-1,lCm1 );
    }  // if (Q-C) > 1.0e-15 


    Zeta = 1.0/pripair1.one_over_gamma + 1.0/pripair2.one_over_gamma; 
/*
      for ( int mu = 0 ; mu<3 ; mu++ )
      W[mu] = (pripair1.P[mu]/pripair1.one_over_gamma 
              + pripair2.P[mu]/pripair2.one_over_gamma )/Zeta ; 
*/
    Ptilde_iwork = pripair1.P[iWork] + 0.5*k1[iWork]*pripair1.one_over_gamma*onei;

    Wtilde_iwork = ( Ptilde_iwork / pripair1.one_over_gamma + 
               Qtilde_iwork / pripair2.one_over_gamma )/Zeta;  

// here should be W - Qtilde
 
    if( std::abs( Wtilde_iwork- Qtilde_iwork ) > 1.0e-15 ) {
      tmpVal += ( Wtilde_iwork- Qtilde_iwork ) * twoecompvRRa0c0( pripair1, pripair2, 
                 FmT_2epri, k1,k2, ss_pair1,ss_pair2,shell1,shell3, m+1, LA,lA, LC-1,lCm1 );
    } // if( abs( W_iWork- Q[iWork] ) > 1.0e-15 )

    if (lA[iWork]>0) {

      lAm1[iWork] -= 1;
      tmpVal += (lAm1[iWork]+1) / (2.0*Zeta) * twoecompvRRa0c0( pripair1, pripair2, 
                 FmT_2epri, k1,k2,ss_pair1,ss_pair2,shell1, shell3, m+1, 
                 LA-1,lAm1, LC-1,lCm1 );
    } // if (lA[iWork]>0) 

    if ( lC[iWork]>=2 ){

      lCm1[iWork] -=1; // right now lCm1(iWork) = lC[iWork]-2 
      tmpVal += 0.5 * (lCm1[iWork]+1) * pripair2.one_over_gamma * 
        ( twoecompvRRa0c0( pripair1, pripair2, FmT_2epri, k1,k2,  ss_pair1, ss_pair2, 
                       shell1, shell3, m, LA,lA, LC-2,lCm1 )

               - ( 1.0/pripair1.one_over_gamma )/Zeta 
          * twoecompvRRa0c0( pripair1, pripair2, FmT_2epri,k1,k2, ss_pair1,ss_pair2, 
                         shell1, shell3, m+1, LA,lA, LC-2,lCm1 ) );
    } // if ( lC[iWork]>=2 )

    return tmpVal;
  }  // twoecompvRRa0c0



  /**
   *  \brief vertical recursion of complex ERI when all LA > 0, all the others are 0
   *
   *
   *  \param [in] pripair1  primitive bra shell pair data for shell1,shell2
   *  \param [in] pripair2  primitive ket shell pair data for shell3,shell4
   *  \param [in] pair1index index of primitive shell pair among the contracted pair
   *  \param [in] pair2index index of primitive shell pair among the contracted pair
   *  \param [in] FmT_2e Boys function between two primitive shell pairs
   *  \param [in] shell1    
   *  \param [in] m         order of auxiliary function
   *  \param [in] LA
   *  \param [in] lA
   *
   *  \return complex ERI of two primitive shell pairs ( shell1 shell2 | shell3 shell4 )
   */
//---------------------------------------------------------------//
// two-e vertical recursion from [a0|00] to [00|00]              //
// [a0|00]^m = (P-A)*[a-1,0|00]^m                                //
//           + (W-P)*[a-1,0|00]^(m+1)                            //
//           + (N(a)-1)/(2*zeta)*[a-2,0|00]^m                    //
//           - (N(a)-1)/(2*zeta)*eta/(zeta+eta)*[a-2,0|00]^(m+1) //
//---------------------------------------------------------------//
   
  dcomplex ComplexGIAOIntEngine::twoecompvRRa000(
    libint2::ShellPair::PrimPairData &pripair1,
    libint2::ShellPair::PrimPairData &pripair2, std::vector<dcomplex> &FmT_2epri,
    double *k1, double *k2, dcomplex ss_pair1, dcomplex ss_pair2, 
    libint2::Shell &shell1, int m,int LA,int *lA ) {


    if(LA==0) {
      // calculate the (SS||SS) integral  	
      //double zetaG; 
      dcomplex SSSS=0.0 ; 
 
      // zeta+eta is zetaG
      //zetaG = 1.0/pripair1.one_over_gamma + 1.0/pripair2.one_over_gamma;  
      auto rho = 1.0 / ( pripair1.one_over_gamma + pripair2.one_over_gamma ); 

      SSSS += 2.0 * ss_pair1 * ss_pair2 * FmT_2epri[m] * sqrt(rho/M_PI); 
 
      return SSSS;
 
    } // if LA==0
 
    // here LA != 0
    double Zeta;
    dcomplex tmpVal;
    dcomplex Ptilde_iwork, Qtilde_iwork, Wtilde_iwork;
    int iWork;
    int lAm1[3];
    dcomplex onei;
    onei.real(0);
    onei.imag(1);
 
    for( iWork=0 ; iWork<3 ; iWork++ ) lAm1[iWork]=lA[iWork];
 
    if (lA[0]>0) iWork=0;
    else if (lA[1]>0) iWork=1;
    else if (lA[2]>0) iWork=2;
 
    if( LA>=1 ) {
 
      lAm1[iWork]-=1;

      Zeta = 1.0/pripair1.one_over_gamma + 1.0/pripair2.one_over_gamma;

      Ptilde_iwork = pripair1.P[iWork] + 0.5*k1[iWork]*pripair1.one_over_gamma*onei; 
      Qtilde_iwork = pripair2.P[iWork] + 0.5*k2[iWork]*pripair2.one_over_gamma*onei;  
      Wtilde_iwork = ( Ptilde_iwork / pripair1.one_over_gamma + 
                 Qtilde_iwork / pripair2.one_over_gamma )/Zeta;
      
/*
      for ( int mu = 0 ; mu<3 ; mu++ )
        W[mu] = (pripair1.P[mu]/pripair1.one_over_gamma 
              + pripair2.P[mu]/pripair2.one_over_gamma )/Zeta ; 
*/
 
      if ( std::abs( Wtilde_iwork - Ptilde_iwork ) > 1.0e-15 ) {
 
        tmpVal += ( Wtilde_iwork - Ptilde_iwork ) * twoecompvRRa000( pripair1, pripair2,
          FmT_2epri, k1, k2, ss_pair1, ss_pair2, shell1, m+1, LA-1,lAm1 ); 
 
      }  // if ( abs( W_iWork-P[iWork] )>1.0e-15 )

      if ( std::abs( Ptilde_iwork - shell1.O[iWork] )>1.0e-15 ) {  
        tmpVal+= ( Ptilde_iwork - shell1.O[iWork] ) * twoecompvRRa000( pripair1, 
          pripair2, FmT_2epri, k1,k2, ss_pair1, ss_pair2, shell1, m, LA-1,lAm1 );
      } // if ( abs( P[iWork]-A[iWork] )>1.0e-15 ) 

      if ( lA[iWork]>=2 ) {

        lAm1[iWork] -=1; // now lAm1[iWork] == lA[iWork]-2
        tmpVal += 0.5 * ( lAm1[iWork]+1 ) * pripair1.one_over_gamma  
                  * ( twoecompvRRa000( pripair1, pripair2, FmT_2epri, k1, k2, ss_pair1,
                    ss_pair2, shell1, m, LA-2,lAm1 )
                - 1.0/(pripair2.one_over_gamma*Zeta) * twoecompvRRa000( pripair1, pripair2,
                  FmT_2epri, k1,k2,ss_pair1,ss_pair2, shell1, m+1, LA-2,lAm1 ) );
      } // if lA[iWork]>=2

    } // if( LA>=1 ) 

    return tmpVal;

  };  // twoevRRa000


  /**
   *  \brief compute complex ERI when all the angular momentum are 0
   *
   *
   *  \param [in] pair1   bra shell pair data for shell1,shell2
   *  \param [in] pair2   ket shell pair data for shell3,shell4
   *  \param [in] shell1    
   *  \param [in] shell2    
   *  \param [in] shell3    
   *  \param [in] shell4    
   *
   *  \return ERI of two contracted shell pairs ( shell1 shell2 | shell3 shell4 )
   */
   
  dcomplex ComplexGIAOIntEngine::twoecompSSSS0(
    libint2::ShellPair &pair1, libint2::ShellPair &pair2,
    libint2::Shell &shell1, libint2::Shell &shell2,
    libint2::Shell &shell3, libint2::Shell &shell4, 
    std::vector<std::vector<std::vector<dcomplex>>> &FmT_2e,
    std::vector<dcomplex> &ss_shellpair1, std::vector<dcomplex> & ss_shellpair2 ) {

    // in this auxiliary integral, m=0 
    // double norm;
    dcomplex SSSS0=0.0 ;
    int pair1index = 0, pair2index = 0;

    //loop over shellpairs
    for ( auto &pripair1 : pair1.primpairs ) {
      pair2index = 0;
      for ( auto &pripair2 : pair2.primpairs ) {

//std::cout<<"pair2 = "<<pair1index<<" pair2 = "<<pair2index<<std::endl;
       
      // zeta+eta is zetaG
      // zetaG = 1.0/pripair1.one_over_gamma + 1.0/pripair2.one_over_gamma;  
        auto rho = 1.0 / ( pripair1.one_over_gamma + pripair2.one_over_gamma );
/*
        norm = shell1.contr[0].coeff[pripair1.p1]* 
             shell2.contr[0].coeff[pripair1.p2]* 
             shell3.contr[0].coeff[pripair2.p1]* 
             shell4.contr[0].coeff[pripair2.p2];  
*/
        // norm is included in the overlap ss integral

        SSSS0 += 2.0 *  ss_shellpair1[pair1index] * ss_shellpair2[pair2index] 
            * FmT_2e[pair1index][pair2index][0] * sqrt(rho/M_PI); 

        pair2index++;
      } // for pripair2
      pair1index++;
    }// for pripair1

// std::cout<<"calculate ssss0 "<<SSSS0<<std::endl;
    return SSSS0;
    
  } // twoeSSSS0

}  // namespace ChronusQ 

