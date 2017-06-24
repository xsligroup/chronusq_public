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

//SS start

  std::vector<double> RealGTOIntEngine::bottomupERI(libint2::ShellPair &pair1 ,
    libint2::ShellPair &pair2, libint2::Shell &shell1, libint2::Shell &shell2,
    libint2::Shell &shell3, libint2::Shell &shell4, int deria, int derib, int deric, int derid) {

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
    
    double *FT = new double[Ltot+1];   

   
    // allocate memory for horizontal recursion
    std::vector<std::vector<std::vector<double>>> Vbraketee;
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
      double zeta = 1.0/pripair1.one_over_gamma; 
      double ss1 = shell1.contr[0].coeff[pripair1.p1]* shell2.contr[0].coeff[pripair1.p2]
        * pow(sqrt(M_PI),3) * sqrt(pripair1.one_over_gamma)*pripair1.K ; 

      for ( auto &pripair2 : pair2.primpairs ) {

        double Q[3]; 
        for (int ii = 0 ; ii < 3 ; ii++ ) {
          Q[ii] = pripair2.P[ii]; 
        }
        double eta = 1.0/pripair2.one_over_gamma; 
        double zetaG = zeta + eta;
        double rho = zeta*eta/zetaG;

        double T = 0.0;
        for (int ii = 0 ; ii < 3 ; ii++) {
          T += pow( P[ii]-Q[ii] ,2); 
        }
        T*= rho ; 
        
        double ss2 = shell3.contr[0].coeff[pripair2.p1]* shell4.contr[0].coeff[pripair2.p2]
          * pow(sqrt(M_PI),3) * sqrt(pripair2.one_over_gamma)*pripair2.K ; 
        
        RealGTOIntEngine::computeFmTTaylor(FT,T,Ltot,0);        
//FT // FIXME calculate boys function
        double W[3];
        for (int ii = 0 ; ii < 3 ; ii++ )
          W[ii] = (zeta*P[ii] + eta*Q[ii])/zetaG;         

        double pref = 2*sqrt(rho/M_PI)*ss1*ss2;
        // here multiply exponents for derivatives
        pref *= pow( shell1.alpha[pripair1.p1], deria);
        pref *= pow( shell2.alpha[pripair1.p2], derib);
        pref *= pow( shell3.alpha[pripair2.p1], deric);
        pref *= pow( shell4.alpha[pripair2.p2], derid);

        
        // now allocate the memory to store vertical recursion elements
        std::vector<std::vector<std::vector<double>>> Vtempbraket((Lbra+1)*(Lket+1)); 
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
                double ERIscratch = 0.0;
                ERIscratch = (P[iWork]-A[iWork])* Vtempbraket[(k-1)*(Lket+1)][indexlm1*cart_ang_list[l].size()][m]
                  +(W[iWork]-P[iWork])*Vtempbraket[(k-1)*(Lket+1)][indexlm1*cart_ang_list[l].size()][m+1];

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
                    double ERIscratch = 0.0 ; 
                    ERIscratch += (Q[iWork]-C[iWork])*Vtempbraket[(k*(Lket+1)+l-1)][cart_i
                      *cart_ang_list[l-1].size()+indexlm1][m] 
                      +(W[iWork]-Q[iWork])* Vtempbraket[(k*(Lket+1)+l-1)][cart_i*
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

    std::vector<double> ERI_sph;

    ERI_sph.assign(((2*shell1.contr[0].l+1)*(2*shell2.contr[0].l+1)
                   *(2*shell3.contr[0].l+1)*(2*shell4.contr[0].l+1)),0.0); 

    cart2sph_2e_transform( shell1.contr[0].l,shell2.contr[0].l,
      shell3.contr[0].l,shell4.contr[0].l,ERI_sph,Vbraketee[0][LB*(LD+1)+LD] );

    return ERI_sph; 


  } //RealGTOIntEngine::bottomupERI


//SS end




  std::vector<double> RealGTOIntEngine::BottomupHGP( 
    libint2::ShellPair &pair1, libint2::ShellPair &pair2, 
    libint2::Shell &shell1, libint2::Shell &shell2,
    libint2::Shell &shell3, libint2::Shell &shell4,
    int deria, int derib, int deric, int derid ) {

    int LA, LB, LC, LD;
    LA = shell1.contr[0].l;
    LB = shell2.contr[0].l;
    LC = shell3.contr[0].l;
    LD = shell4.contr[0].l;

    
    int L = shell1.contr[0].l + shell2.contr[0].l 
            + shell3.contr[0].l + shell4.contr[0].l;

    int Lbra = shell1.contr[0].l + shell2.contr[0].l;
    int Lket = shell3.contr[0].l + shell4.contr[0].l;
    
    double A[3],B[3],C[3],D[3];

    for (int iWork = 0; iWork < 3; iWork++) {

      A[iWork] = shell1.O[iWork];
      B[iWork] = shell2.O[iWork];
      C[iWork] = shell3.O[iWork];
      D[iWork] = shell4.O[iWork];

    }

    double AB[3], CD[3];

    for (int iWork = 0; iWork < 3; iWork++) {

      AB[iWork] = A[iWork] - B[iWork];
      CD[iWork] = C[iWork] - D[iWork];

    }

    // here allocate the boys function  

    double *FT = new double[L + 1];  
    double P[3],Q[3];
 
    // Allocate the vectors for Vbraket 
    std::vector<std::vector<std::vector<double>>> Vbraket((Lbra - LA + 1) 
        * (Lket - LC + 1));

    for (int lA = LA; lA <= Lbra; lA++) {

      for (int lC = LC; lC <= Lket; lC++) {

        Vbraket[(lA - LA) * (Lket - LC + 1) + (lC - LC)].resize((Lbra - lA + 1)
            * (Lket - lC + 1));

        for (int lB = 0; lB <= Lbra - lA; lB++) {

          for (int lD = 0; lD <= Lket - lC; lD++) {

            Vbraket[(lA - LA) * (Lket - LC + 1) + (lC - LC)][
                    lB * (Lket - lC + 1) + lD
                   ].assign(cart_ang_list[lA].size() 
                            * cart_ang_list[lB].size() 
                            * cart_ang_list[lC].size() 
                            * cart_ang_list[lD].size(), 0.0);

          } // for (int lD = 0; lD <= Lket - lC; lD++)

        } // for (int lB = 0; lB <= Lbra - lA; lB++)

      } // for (int lC = LC; lC <= Lket; lC++)

    } // for (int lA = LA; lA <= Lbra; lA++)

    // loop over primitive shellpairs in bra side 
    for (auto &pripair1 : pair1.primpairs) { 

      for (int iWork = 0; iWork < 3; iWork++) {

        P[iWork] = pripair1.P[iWork];

      } // for (int iWork = 0; iWork < 3; iWork++)

      for (auto &pripair2 : pair2.primpairs) {

        for (int iWork = 0; iWork < 3; iWork++) {

          Q[iWork] = pripair2.P[iWork];

        } //  for (int iWork = 0; iWork < 3; iWork++)

        double sqrPQ = 0.0;

        for (int mu = 0; mu < 3; mu++) {

          double PQ = (pripair1.P[mu] - pripair2.P[mu]); 
          sqrPQ += PQ * PQ;

        } // for (int mu=0; mu < 3; mu++)

        auto zeta = 1.0 / pripair1.one_over_gamma;
        auto eta  = 1.0 / pripair2.one_over_gamma;
        auto zetaG = zeta + eta;
        auto rho = zeta * eta / zetaG;
        auto T = rho * sqrPQ;

        // calculate Fm(T) list
        computeFmTTaylor(FT, T, L, 0);

        // zeta + eta is expoT
        double expoT = 1.0 / pripair1.one_over_gamma 
                       + 1.0 / pripair2.one_over_gamma;
      
        // T= rho * PQ^2
        // T = sqrPQ / (pripair1.one_over_gamma + pripair2.one_over_gamma);
      
        // calculate F0(T)
        // computeFmTTaylor(FT, T, 0, 0);
      
        double Kab = sqrt(2.0) * pow(M_PI, 1.25) * pripair1.K; 
        double Kcd = sqrt(2.0) * pow(M_PI, 1.25) * pripair2.K;
      
        double norm = shell1.contr[0].coeff[pripair1.p1]
                      * shell2.contr[0].coeff[pripair1.p2]
                      * shell3.contr[0].coeff[pripair2.p1]
                      * shell4.contr[0].coeff[pripair2.p2];  
     
        //SS start
        // here multiply exponents for derivatives
        norm *= pow( shell1.alpha[pripair1.p1], deria);
        norm *= pow( shell2.alpha[pripair1.p2], derib);
        norm *= pow( shell3.alpha[pripair2.p1], deric);
        norm *= pow( shell4.alpha[pripair2.p2], derid);

        //SS end
         
        double pref = norm * Kab * Kcd / sqrt(expoT); 

        std::vector<std::vector<std::vector<double>>> Vtempbraket((Lbra + 1) 
            * (Lket + 1));

        for (int k = 0; k <= Lbra; k++) {
        
          for (int l = 0; l <= Lket; l++) {

            Vtempbraket[k * (Lket + 1) + l].resize(cart_ang_list[k].size()
                                                   * cart_ang_list[l].size());

            int mbraket = L - k - l;

            for (int cart_i = 0; cart_i < cart_ang_list[k].size(); cart_i++) {

              for (int cart_j = 0; cart_j < cart_ang_list[l].size(); cart_j++) 
              {

                Vtempbraket[k * (Lket + 1) + l][
                  cart_i * cart_ang_list[l].size() + cart_j
                  ].resize(mbraket + 1); 

              } // for (int cart_j = 0; cart_j < cart_ang_list[l].size(); 
                // cart_j++)

            } // for (int cart_i = 0; cart_i < cart_ang_list[k].size(); 
              // cart_i++)
          
          } // for (int l = 0; l <= Lket; l++)

        } // for (int k = 0; k <= Lbra; k++)

        for (int ii = 0; ii <= L; ii++) {

          Vtempbraket[0][0][ii] = FT[ii] * pref;

          // std::cout << "FT" << ii << "= " << FT[ii] << " pref " << pref << 
          //   std::endl;

        } // for (int ii = 0; ii <= L; ii++)

        double W[3];

        for (int ii = 0; ii < 3; ii++) {

          W[ii] = (zeta * P[ii] + eta * Q[ii]) / zetaG;

        } // for (int ii = 0; ii < 3; ii++)

        for (int k = 0; k <= Lbra; k++) { // Loop over the bra angular momentum
          
          for (int cart_i = 0; cart_i < cart_ang_list[k].size(); cart_i++) {

            int lA_xyz[3];     

            for (int ii = 0; ii < 3; ii++) { 

              lA_xyz[ii] = cart_ang_list[k][cart_i][ii];

            } // for (int ii = 0; ii < 3; ii++)

            int mbra = L - k; // mbra is the highest auxiliary number for a 
            // given k

            if (k > 0) {

              int iwork;

              if (lA_xyz[0] > 0) {

                iwork = 0;

              } // if (lA_xyz[0] > 0)

              else if (lA_xyz[1] > 0) {

                iwork = 1;

              } // else if (lA_xyz[1] > 0)

              else if (lA_xyz[2] > 0) { 

                iwork = 2;

              } // else if (lA_xyz[2] > 0)

              // Calculate the index of a lower angular momentum integral
              int lAtemp[3];

              for (int ii = 0; ii < 3; ii++) { 

                lAtemp[ii] = lA_xyz[ii];

              } // for (int ii = 0; ii < 3; ii++)

              lAtemp[iwork] = lA_xyz[iwork] - 1;

              int indexlm1xyz = indexmap(k - 1, lAtemp[0], lAtemp[1], 
                                         lAtemp[2]);

              // Calculate index of lAm1
              
              for (int m = 0; m <= mbra; m++) {
              
                double ERIscratch = 0.0;

                ERIscratch = (P[iwork] - A[iwork]) 
                             * Vtempbraket[
                             (k - 1) * (Lket + 1)][
                             indexlm1xyz][
                             m]
                             + (W[iwork] - P[iwork])
                             * Vtempbraket[
                             (k - 1) * (Lket + 1)][
                             indexlm1xyz][
                             m + 1];

                // if (std::abs(ERIscratch) > 1.0e-10) {
                //  std::cout << "Line 178. ERIscratch = " << ERIscratch << 
                //    std::endl;
                // }
                // Equation 6. First line HGP

                if (lA_xyz[iwork] > 1) {

                  for (int ii = 0; ii < 3; ii++) {

                    lAtemp[ii] = lA_xyz[ii];

                  } // for (int ii = 0; ii < 3; ii++)

                  lAtemp[iwork] = lA_xyz[iwork] - 2;

                  int indexlm2xyz = indexmap(k - 2, lAtemp[0], lAtemp[1], 
                                             lAtemp[2]);
                  ERIscratch = ERIscratch + 1 / (2 * zeta) * (lA_xyz[iwork] - 1) 
                               * (Vtempbraket[
                                   (k - 2) * (Lket + 1)][
                                   indexlm2xyz][
                                   m]
                               - rho / zeta 
                               * Vtempbraket[
                               (k - 2) * (Lket + 1)][
                               indexlm2xyz][
                               m + 1]);

                  // if (std::abs(ERIscratch) > 1.0e-10) { 
                  //   std::cout << "Line 197. ERIscratch = " << ERIscratch << 
                  //     std::endl;
                  // }
                  // Equation 6. Second line HGP

                } // if (lA_xyz[iwork] > 1)

                Vtempbraket[k * (Lket + 1)][cart_i][m] = ERIscratch;

              } // for (int m = 0; m <= mbra; m++) 

            } // if (k > 0)

            for (int l = 0; l <= Lket; l++) {
            
              for (int cart_j = 0; cart_j < cart_ang_list[l].size(); cart_j++) 
              {
              
                int lC_xyz[3];

                for (int ii = 0; ii < 3; ii++) {
                
                  lC_xyz[ii] = cart_ang_list[l][cart_j][ii];

                } // for (int ii = 0; ii < 3; ii++)

                int mbraket = L - k - l;

                if (l > 0) {

                  int iwork;

                  if (lC_xyz[0] > 0) {

                    iwork = 0;

                  } // if (lC_xyz[0] > 0)

                  else if (lC_xyz[1] > 0) {

                    iwork = 1;

                  } // else if (lC_xyz[1] > 0)

                  else if (lC_xyz[2] > 0) {

                    iwork = 2;

                  } // else if (lC_xyz[2] > 0)

                  // Calculate the index of a lower angular momentum integral
                  int lCtemp[3];

                  for (int ii = 0; ii < 3; ii++) {

                    lCtemp[ii] = lC_xyz[ii];

                  } // for (int ii = 0; ii < 3; ii++)

                  // For ket side. Equation 6. HGP
                  lCtemp[iwork] = lC_xyz[iwork] - 1;

                  int indexlm1xyz = indexmap(l - 1, lCtemp[0], lCtemp[1], 
                                             lCtemp[2]);

                  // Calculate indeix of lCm1
                  for (int m = 0; m <= mbraket; m++) {

                    double ERIscratch = 0.0;

                    ERIscratch = (Q[iwork] - C[iwork]) 
                                 * Vtempbraket[
                                 k * (Lket + 1) + (l - 1)][
                                 (cart_i) * cart_ang_list[l - 1].size() 
                                 + indexlm1xyz][
                                 m] 
                                 + (W[iwork] - Q[iwork]) 
                                 * Vtempbraket[
                                 k * (Lket + 1) + (l - 1)][
                                 cart_i * cart_ang_list[l - 1].size() 
                                 + indexlm1xyz][
                                 m + 1];

                    // if (std::abs(ERIscratch) > 1.0e-10) { 
                    //   std::cout << "Line 264. ERIscratch = " << ERIscratch << 
                    //     std::endl;
                    // }

                    if (lC_xyz[iwork] > 1) {

                      for (int ii = 0; ii < 3; ii++) {

                        lCtemp[ii] = lC_xyz[ii];

                      } // for (int ii = 0; ii < 3; ii++)

                      lCtemp[iwork] = lC_xyz[iwork] - 2;

                      int indexlm2xyz = indexmap(l - 2, lCtemp[0], lCtemp[1], 
                                                 lCtemp[2]);
                      ERIscratch = ERIscratch 
                                   + 1 / (2 * eta) * (lC_xyz[iwork] - 1) 
                                   * (Vtempbraket[
                                   k * (Lket + 1) + (l - 2)][
                                   cart_i * cart_ang_list[l - 2].size() 
                                   + indexlm2xyz][m] 
                                   - rho / eta 
                                   * Vtempbraket[
                                   k * (Lket + 1) + (l - 2)][
                                   cart_i * cart_ang_list[l - 2].size() 
                                   + indexlm2xyz][
                                   m + 1]);

                      // if (std::abs(ERIscratch) > 1.0e-10) { 
                      //   std::cout << "Line 277. ERIscratch = " << ERIscratch << 
                      //     std::endl;
                      // }

                    } // if (lC_xyz[iwork] > 1)

                    if (lA_xyz[iwork] > 0) {

                      int lAtemp[3];

                      for (int ii = 0; ii < 3; ii++) {

                        lAtemp[ii] = lA_xyz[ii];

                      } // for (int ii = 0; ii < 3; ii++)

                      lAtemp[iwork] -= 1; 

                      int indexlAm1 = indexmap(k - 1, lAtemp[0], lAtemp[1], 
                                               lAtemp[2]);
                      ERIscratch += 1.0 / (2.0 * zetaG) * lA_xyz[iwork]
                                    * Vtempbraket[
                                    (k - 1) * (Lket + 1) + (l - 1)][
                                    (indexlAm1) * cart_ang_list[l - 1].size() 
                                    + indexlm1xyz][
                                    m + 1];

                      // if (std::abs(ERIscratch) > 1.0e-10) {
                      //   std::cout << "Line 299. ERIscratch = " << ERIscratch << 
                      //     std::endl;
                      // }

                    } // if (lA_xyz[iwork] > 0)

                    Vtempbraket[k * (Lket + 1) + l][
                      cart_i * cart_ang_list[l].size() + cart_j][
                      m] 
                      = ERIscratch;

                  } // for (int m = 0; m <= mbraket; m++)

                } // if (l > 0)

              } // for (int cart_j = 0; cart_j < cart_ang_list[l].size(); 
                // cart_j++)

            } // for (int l = 0; l <= Lket; l++)

          } // for (int cart_i = 0; cart_i < cart_ang_list[k].size(); cart_i++)

        } // for (int k = 0; k <= Lbra; k++)
        
        // Start Contraction Here
        // Vbracket[]+= Vtempbraket[]
        for (int lA = LA; lA <= Lbra; lA++) {

          for (int lC = LC; lC <= Lket; lC++) {

            for (int i = 0; i < cart_ang_list[lA].size(); i++) {

              for (int k = 0; k < cart_ang_list[lC].size(); k++) {

                Vbraket[(lA - LA) * (Lket - LC + 1) + (lC - LC)][0][
                  i * cart_ang_list[lC].size() + k] 
                  += Vtempbraket[
                  lA * (Lket + 1) + lC][i * cart_ang_list[lC].size() + k][0];

                // if (std::abs(Vbraket[(lA - LA) * (Lket - LC + 1) + lC - LC][0][i * cart_ang_list[lC].size() + k]) > 1.0e-10) { 
                //   std::cout << "Line 322. Vbraket = " << Vbraket[(lA - LA) * (Lket - LC + 1) + lC - LC][0][i * cart_ang_list[lC].size() + k] << std::endl;
                // }

              } // for (int k = 0; k <= cart_ang_list[lC].size(); k++)

            } // for (int i = 0; i <= cart_ang_list[lA].size(); i++)

          } // for (int lC = LC; lC <= Lket; lC++)

        } // for (int lA = LA; lA <= Lbra; lA++)

      } // for ( auto &pripair2 : pair2.primpairs )   

    }  // for ( auto &pripair1 : pair1.primpairs )

    for (int lB = 1; lB <= LB; lB++) {

      // Implies LB > 0
      for (int lA = LA; lA <= Lbra - lB; lA++) {

        // Loop over (lA||lB)
        // Loop over elements in lA
        for (int Aidx = 0; Aidx < cart_ang_list[lA].size(); Aidx++) {

          int lA_xyz[3];

          for (int ii = 0; ii < 3; ii++) {

            lA_xyz[ii] = cart_ang_list[lA][Aidx][ii];

          } // for (int ii = 0; ii < 3; ii++)

          for (int Bidx = 0; Bidx < cart_ang_list[lB].size(); Bidx++) {

            int lB_xyz[3];

            for (int ii = 0; ii < 3; ii++) {

              lB_xyz[ii] = cart_ang_list[lB][Bidx][ii];

            } // for (int ii = 0; ii < 3; ii++)

            int iwork;

            if (lB_xyz[0] > 0) { // Means lA_x > 0

              iwork = 0;

            } // if (lB_xyz[0] > 0)

            else if (lB_xyz[1] > 0) { // Means lA_y > 0

              iwork = 1;

            } // else if (lB_xyz[1] > 0)

            else if (lB_xyz[2] > 0) { // Means lA_z > 0

              iwork = 2;

            } // else if (lB_xyz[2] > 0)

            int lBm1[3];

            for (int ii = 0; ii < 3; ii++) {

              lBm1[ii] = lB_xyz[ii];

            } // for (int ii = 0; ii < 3; ii++)

            lBm1[iwork] = lB_xyz[iwork] - 1;

            int lAp1[3];

            for (int ii = 0; ii < 3; ii++) {

              lAp1[ii] = lA_xyz[ii];

            } // for (int ii = 0; ii < 3; ii++)

            lAp1[iwork] = lA_xyz[iwork] + 1;

            int idxBtemp = indexmap(lB - 1, lBm1[0], lBm1[1], lBm1[2]);
            int lD = 0;

            for (int lC = LC; lC <= Lket; lC++) {

              for (int Cidx = 0; Cidx < cart_ang_list[lC].size(); Cidx++) {

                int lC_xyz[3];

                for (int ii = 0; ii < 3; ii++) {

                  lC_xyz[ii] = cart_ang_list[lC][Cidx][ii];

                } // for (int ii = 0; ii < 3; ii++)

                Vbraket[(lA - LA) * (Lket - LC + 1) + (lC - LC)][
                  lB * (Lket - lC + 1)][
                  Aidx * cart_ang_list[lB].size() * cart_ang_list[lC].size() 
                  + Bidx * cart_ang_list[lC].size()
                  + Cidx] 
                  = Vbraket[(lA - LA + 1) * (Lket - LC + 1) + (lC - LC)][
                  (lB - 1) * (Lket - lC + 1)][
                  indexmap(lA + 1, lAp1[0], lAp1[1], lAp1[2]) 
                  * cart_ang_list[lB - 1].size() * cart_ang_list[lC].size() 
                  + idxBtemp * cart_ang_list[lC].size() + Cidx] 
                  + (AB[iwork]) 
                  * Vbraket[(lA - LA) * (Lket - LC + 1) + lC - LC][
                  (lB - 1) * (Lket - lC + 1)][
                  Aidx * cart_ang_list[lB - 1].size() 
                       * cart_ang_list[lC].size() 
                  + idxBtemp * cart_ang_list[lC].size() 
                  + Cidx];

                // if (std::abs(Vbraket[(lA - LA) * (Lket - LC + 1) + lC - LC][lB * (Lket - lC + 1)][Aidx * cart_ang_list[lB].size() * cart_ang_list[lC].size() + Bidx * cart_ang_list[lC].size() + Cidx]) > 1.0e-10) { 
                //   std::cout << "Line 377. Vbraket = " << Vbraket[(lA - LA) * (Lket - LC + 1) + lC - LC][lB * (Lket - lC + 1)][Aidx * cart_ang_list[lB].size() * cart_ang_list[lC].size() + Bidx * cart_ang_list[lC].size() + Cidx] << std::endl;
                // }

              } // for (int Cidx = 0; Cidx < cart_ang_list[lC]; Cidx++)

            } // for (int lC = LC; lC <= Lket; lC++)

          } // for (int Bidx = 0; Bidx < cart_ang_list[lB].size(); Bidx++)

        } // for (int Aidx = 0; Aidx < cart_ang_list[lA].size(); Aidx++)

      } // for (int lA = LA; lA <= Lbra - lB; lB++)

    } // for (int lB = 1; lB <= LB; lB++)

    // Horizontal of bra side finished
    // Horizontal of ket side start
    for (int lD = 1; lD <= LD; lD++) {

      // Implies LD > 0
      for (int lC = LC; lC <= Lket - lD; lC++) {

        // Loop over (lC||lD)
        // Loop over elements in lC
        for (int Cidx = 0; Cidx < cart_ang_list[lC].size(); Cidx++) {

          // Calculate the lC_x, lC_y, lC_z
          int lC_xyz[3];

          for (int ii = 0; ii < 3; ii++) {

            lC_xyz[ii] = cart_ang_list[lC][Cidx][ii];

          } // for (int ii = 0; ii < 3; ii++)
          // Loop over elements in lD

          for (int Didx = 0; Didx < cart_ang_list[lD].size(); Didx++) {

            int lD_xyz[3];

            for (int ii = 0; ii < 3; ii++) {

              lD_xyz[ii] = cart_ang_list[lD][Didx][ii];

            } // for (int ii = 0; ii < 3; ii++)

            int iwork;

            if (lD_xyz[0] > 0) { // Means lC_x > 0

              iwork = 0;

            } // if (lD_xyz[0] > 0)

            else if (lD_xyz[1] > 0) { // Means lC_y > 0

              iwork = 1;

            } // else if (lD_xyz[1] > 0)

            else if (lD_xyz[2] > 0) { // Means lC_z > 0

              iwork = 2;

            } // else if (lD_xyz[2] > 0)

            int lDm1[3];

            for (int ii = 0; ii < 3; ii++) {

              lDm1[ii] = lD_xyz[ii];

            } // for (int ii = 0; ii < 3; ii++)

            lDm1[iwork] = lD_xyz[iwork] - 1;

            int lCp1[3];

            for (int ii = 0; ii < 3; ii++) {

              lCp1[ii] = lC_xyz[ii];

            } // for (int ii = 0; ii < 3; ii++)

            lCp1[iwork] = lC_xyz[iwork] + 1;

            int idxDtemp = indexmap(lD - 1, lDm1[0], lDm1[1], lDm1[2]);

            for (int Aidx = 0; Aidx < cart_ang_list[LA].size(); Aidx++) {

              // Calculate the lA_x, lA_y, lA_z
              int lA_xyz[3];

              for (int ii = 0; ii < 3; ii++) {

                lA_xyz[ii] = cart_ang_list[LA][Aidx][ii];

              } // for (int ii = 0; ii < 3; ii++)
              for (int Bidx = 0; Bidx < cart_ang_list[LB].size(); Bidx++) {

                int lB_xyz[3];

                for (int ii = 0; ii < 3; ii++) {

                  lB_xyz[ii] = cart_ang_list[LB][Bidx][ii];

                } // for (int ii = 0; ii < 3; ii++)

                Vbraket[lC - LC][LB * (Lket - lC + 1) + lD][
                  Aidx * cart_ang_list[LB].size() 
                  * cart_ang_list[lC].size() 
                  * cart_ang_list[lD].size() 
                  + Bidx * cart_ang_list[lC].size()* cart_ang_list[lD].size() 
                  + Cidx * cart_ang_list[lD].size()
                  + Didx] 
                  = Vbraket[lC - LC + 1][
                  LB * (Lket - (lC + 1) + 1) + (lD - 1)][
                  Aidx * cart_ang_list[LB].size() 
                    * cart_ang_list[lC + 1].size() 
                    * cart_ang_list[lD - 1].size() 
                  + Bidx * cart_ang_list[lC + 1].size() 
                    * cart_ang_list[lD - 1].size() 
                  + indexmap(lC + 1, lCp1[0], lCp1[1], lCp1[2]) 
                    * cart_ang_list[lD - 1].size() 
                  + idxDtemp] 
                  + (CD[iwork]) * Vbraket[lC - LC][
                  LB * (Lket - lC + 1) + (lD - 1)][
                  Aidx * cart_ang_list[LB].size() 
                  * cart_ang_list[lC].size() 
                  * cart_ang_list[lD - 1].size() 
                  + Bidx * cart_ang_list[lC].size() 
                  * cart_ang_list[lD - 1].size() 
                  + Cidx * cart_ang_list[lD - 1].size() 
                  + idxDtemp];

                // if (std::abs(Vbraket[lC - LC][LB * (Lket - lC + 1) + lD][Aidx * cart_ang_list[LB].size() * cart_ang_list[lC].size() * cart_ang_list[lD].size() + Bidx * cart_ang_list[lC].size() * cart_ang_list[lD].size() + Cidx * cart_ang_list[lD].size() + Didx]) > 1.0e-10) { 
                //   std::cout << "Line 449. Vbraket = " << Vbraket[lC - LC][LB * (Lket - lC + 1) + lD][Aidx * cart_ang_list[LB].size() * cart_ang_list[lC].size() * cart_ang_list[lD].size() + Bidx * cart_ang_list[lC].size() * cart_ang_list[lD].size() + Cidx * cart_ang_list[lD].size() + Didx] << std::endl;
                // }

              } // for (int Bidx = 0; Bidx < cart_ang_list[LB].size(); Bidx++)

            } // for (int Aidx = 0; Aidx < cart_ang_list[LA].size(); Aidx++)

          } // for (int Didx = 0; Didx < cart_ang_list[lD].size(); Didx++)

        } // for (int Cidx = 0; Cidx < cart_ang_list[lC].size(); Cidx++)

      } // for (int lC = LC; LC <= Lket - lD; lC++)

    } // for (int lD = 1; lD <= LD; lD++)

    if ( ( not shell1.contr[0].pure ) and ( not shell2.contr[0].pure ) and 
         ( not shell3.contr[0].pure ) and ( not shell4.contr[0].pure ) ) {  
      // if both sides are cartesian, return cartesian gaussian integrals
      return Vbraket[0][LB * (LD + 1) + LD];
    }

    std::vector<double> ERI_sph;

    ERI_sph.assign(((2*shell1.contr[0].l+1)*(2*shell2.contr[0].l+1)
                   *(2*shell3.contr[0].l+1)*(2*shell4.contr[0].l+1)),0.0); 

    cart2sph_2e_transform(shell1.contr[0].l,shell2.contr[0].l,
                          shell3.contr[0].l,shell4.contr[0].l,
                          ERI_sph,Vbraket[0][LB * (LD + 1) + LD]);

    return ERI_sph; 



  }  // BottomupHGP

  int indexmap(int L, int x, int y, int z) {
   int a = L - x;
   int b = a - y;
   int indexinshell = a*(a+1)/2+b;

   return indexinshell;

  }

  /**
   *  \brief Compute the ERI of two shell pairs
   *
   *
   *  \param [in] pair1  bra shell pair data for shell1,shell2
   *  \param [in] pair2  ket shell pair data for shell3,shell4
   *  \param [in] shell1
   *  \param [in] shell2
   *  \param [in] shell3
   *  \param [in] shell4
   *
   *  \return ERI of two shell pairs ( shell1 shell2 | shell3 shell4 )
   */
   
  std::vector<double> RealGTOIntEngine::computeERIabcd(libint2::ShellPair &pair1 ,
    libint2::ShellPair &pair2, libint2::Shell &shell1, libint2::Shell &shell2,
    libint2::Shell &shell3, libint2::Shell &shell4)  {
    
    double tmpVal=0.0,sqrPQ,PQ;
    std::vector<double> ERI_cart;
    
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

    // pre calculate all the Boys functions 
    // dimension is FmT_2e[shellpair1.prim][shellpair2.prim][lTotal+1]
    std::vector<std::vector<std::vector<double>>> FmT_2e;
    FmT_2e.resize(pair1.primpairs.size());

    double *FmT = new double[lTotal+1];
    int shellpair1_i=0, shellpair2_j ;
    for ( auto &pripair1 : pair1.primpairs ) {
      FmT_2e[shellpair1_i].resize( pair2.primpairs.size() );

      shellpair2_j = 0 ; 
      for ( auto &pripair2 : pair2.primpairs ) {
        sqrPQ = 0.0;
        for ( int mu=0 ; mu<3 ; mu++ ) {
          PQ = ( pripair1.P[mu]-pripair2.P[mu] ); 
          sqrPQ += PQ*PQ;
        }
        auto Zeta = 1.0/pripair1.one_over_gamma;
        auto Eta  = 1.0/pripair2.one_over_gamma;
        
        auto rho = Zeta*Eta/(Zeta+Eta);
        auto T = rho*sqrPQ;
        // calculate Fm(T) list
        computeFmTTaylor( FmT, T, lTotal, 0 );

        for ( int lcurr = 0 ; lcurr < lTotal+1 ; lcurr++ ) {
           if ( std::abs(FmT[lcurr]) < 1.0e-15 ) 
             FmT_2e[shellpair1_i][shellpair2_j].push_back(0.0);
           else
             FmT_2e[shellpair1_i][shellpair2_j].push_back(FmT[lcurr]);
        } // for lcurr
        shellpair2_j ++;
      } // for pripair2
    shellpair1_i++;
    } // for pripair1
    delete[] FmT;

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

      
      tmpVal = twoehRRabcd(pair1,pair2,shell1,shell2,shell3,shell4,
                 FmT_2e,shell1.contr[0].l, lA, shell2.contr[0].l, lB, 
                 shell3.contr[0].l, lC, shell4.contr[0].l, lD ); 
      
      ERI_cart.push_back(tmpVal);
    }   // for l


    if ( ( not shell1.contr[0].pure ) and ( not shell2.contr[0].pure ) and 
         ( not shell3.contr[0].pure ) and ( not shell4.contr[0].pure ) ) {  
      // if both sides are cartesian, return cartesian gaussian integrals
      return ERI_cart;
    }

    std::vector<double> ERI_sph;

    ERI_sph.assign(((2*shell1.contr[0].l+1)*(2*shell2.contr[0].l+1)
                   *(2*shell3.contr[0].l+1)*(2*shell4.contr[0].l+1)),0.0); 

    cart2sph_2e_transform( shell1.contr[0].l,shell2.contr[0].l,
      shell3.contr[0].l,shell4.contr[0].l,ERI_sph,ERI_cart );
      
      
    

    return ERI_sph; 

  }   // computeERIabcd


  /**
   *  \brief horizontal recursion of ERI when all angular momentum are nonzero
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
   *  \return ERI of two shell pairs ( shell1 shell2 | shell3 shell4 )
   */
  //----------------------------------------------------//
  // two-e horizontal recursion from (ab|cd) to (a0|cd) //
  // (ab|cd)=(a+1,b-1|cd)+(A-B)*(a,b-1|cd)              //
  //----------------------------------------------------//
   
  double RealGTOIntEngine::twoehRRabcd( 
     libint2::ShellPair &pair1 ,libint2::ShellPair &pair2 ,
     libint2::Shell &shell1, libint2::Shell &shell2,
     libint2::Shell &shell3, libint2::Shell &shell4,
     std::vector<std::vector<std::vector<double>>> &FmT_2e,
     int LA,int *lA,int LB,int *lB,int LC,int *lC,int LD,int *lD) {

     double tmpVal = 0.0, tmpVal1=0.0;

  // iWork is used to indicate which Cartesian angular momentum we are reducing (x,y,z)

    int iWork;
    int totalL = LA + LB + LC + LD;


    if(totalL==0) {   // (SS||SS) type 

      return twoeSSSS0(pair1,pair2,shell1,shell2,shell3,shell4);

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

          tmpVal += twoehRRabcd(pair1,pair2,shell1,shell2,shell3,shell4,
                      FmT_2e, LA+1,lAp1,LB-1,lBm1,LC,lC,LD,lD);

          if ( std::abs(pair1.AB[iWork]) > 1.0e-15 ) {
            tmpVal += pair1.AB[iWork]*twoehRRabcd(pair1,pair2,shell1,shell2,shell3,
                        shell4,FmT_2e, LA,lA,LB-1,lBm1,LC,lC,LD,lD);
          } 

        } else if ( LB == 0 ) {
          tmpVal = twoehRRa0cd(pair1,pair2,shell1,shell2,shell3,shell4,FmT_2e,
                                 LA,lA,LC,lC,LD,lD);

        } // LB == 0
      }  // else ( LTOTAL != 0 )

      return tmpVal;

  }  // twoehRRabcd

  /**
   *  \brief horizontal recursion of ERI when all LB=0
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
   *  \return ERI of two shell pairs ( shell1 shell2 | shell3 shell4 )
   */
  //----------------------------------------------------//
  // two-e horizontal recursion from (a0|cd) to (a0|c0) //
  // (a0|cd)=(a,0|c+1,d-1)+(C-D)*(a,0|c,d-1)            //
  //----------------------------------------------------//
   
  double RealGTOIntEngine::twoehRRa0cd(
    libint2::ShellPair &pair1, libint2::ShellPair &pair2,
    libint2::Shell &shell1, libint2::Shell &shell2,
    libint2::Shell &shell3, libint2::Shell &shell4, 
    std::vector<std::vector<std::vector<double>>> &FmT_2e,
    int LA,int *lA,int LC,int *lC,int LD,int *lD)  {

    double tmpVal=0.0;
    if(LD==0) {
      int pair1index=0, pair2index=0;
      // go into the vertical recursion
      for ( auto &pripair1 : pair1.primpairs ) {
        pair2index = 0; 
        for ( auto &pripair2 : pair2.primpairs ) {

          auto norm = 
                 shell1.contr[0].coeff[pripair1.p1]* 
                 shell2.contr[0].coeff[pripair1.p2]* 
                 shell3.contr[0].coeff[pripair2.p1]* 
                 shell4.contr[0].coeff[pripair2.p2];  

          tmpVal +=  norm * twoevRRa0c0( pripair1, pripair2,  
             FmT_2e[pair1index][pair2index], shell1,shell3, 0, LA,lA,LC,lC);

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
      tmpVal = twoehRRa0cd(pair1, pair2, shell1, shell2, shell3, shell4, 
                             FmT_2e,LA,lA, LC+1,lCp1, LD-1,lDm1 ); 
      if ( std::abs(pair2.AB[iWork]) > 1.0e-15 ){
        tmpVal += pair2.AB[iWork] * twoehRRa0cd( pair1, pair2, shell1, shell2, 
               shell3, shell4, FmT_2e, LA,lA, LC,lC, LD-1,lDm1 );
      }

    } // else ( that means LD > 0 )
    return tmpVal;
  }  // twoehRRa0cd


  /**
   *  \brief vertical recursion of ERI when all LA, LC > 0
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
   *  \return ERI of two primitive shell pairs ( shell1 shell2 | shell3 shell4 )
   */

//---------------------------------------------------------------//
// two-e vertical recursion from [a0|c0] to [a0|00]              //
// [a0|c0]^m = (Q-C)*[a0|c-1,0]^m                                //
//           + (W-Q)*[a0|c-1,0]^(m+1)                            //
//           + N(a)/(2*(zeta+eta))*[a-1,0|c-1,0]^(m+1)           //
//           + (N(c)-1)/(2*eta)*[a0|c-2,0]^m                     //
//           - (N(c)-1)/(2*eta)*zeta/(zeta+eta)*[a0|c-2,0]^(m+1) //
//---------------------------------------------------------------//
   
  double RealGTOIntEngine::twoevRRa0c0(
    libint2::ShellPair::PrimPairData &pripair1,
    libint2::ShellPair::PrimPairData &pripair2, 
    std::vector<double> &FmT_2epri, 
    libint2::Shell &shell1, libint2::Shell &shell3,
    int m, int LA, int *lA, int LC, int *lC ) {
 
 
    if(LC==0) return twoevRRa000( pripair1, pripair2, FmT_2epri,
                                  shell1, m, LA, lA );
 
    int lAm1[3],lCm1[3];  
 
    for ( int iWork=0 ; iWork<3 ; iWork++ ){
      lAm1[iWork]=lA[iWork];     
      lCm1[iWork]=lC[iWork];
    }
 
    double tmpVal=0.0;
    double W_iWork,Zeta;
    int iWork;
 
    if (lC[0]>0) iWork=0;
    else if (lC[1]>0) iWork=1;
    else if (lC[2]>0) iWork=2;

    lCm1[iWork]-=1;

    if ( std::abs(pripair2.P[iWork] - shell3.O[iWork]) > 1.0e-15 ) {
      tmpVal += ( pripair2.P[iWork] - shell3.O[iWork] ) * 
                twoevRRa0c0( pripair1, pripair2, FmT_2epri, 
                  shell1, shell3, m, LA,lA, LC-1,lCm1 );
    }  // if (Q-C) > 1.0e-15 


    Zeta = 1.0/pripair1.one_over_gamma + 1.0/pripair2.one_over_gamma; 
/*
      for ( int mu = 0 ; mu<3 ; mu++ )
      W[mu] = (pripair1.P[mu]/pripair1.one_over_gamma 
              + pripair2.P[mu]/pripair2.one_over_gamma )/Zeta ; 
*/

    W_iWork = ( pripair1.P[iWork]/pripair1.one_over_gamma + 
               pripair2.P[iWork]/pripair2.one_over_gamma )/Zeta;

    if( std::abs( W_iWork- pripair2.P[iWork] ) > 1.0e-15 ) {
      tmpVal += ( W_iWork- pripair2.P[iWork] ) * twoevRRa0c0( pripair1, pripair2, 
                 FmT_2epri, shell1,shell3, m+1, LA,lA, LC-1,lCm1 );
    } // if( abs( W_iWork- Q[iWork] ) > 1.0e-15 )

    if (lA[iWork]>0) {

      lAm1[iWork] -= 1;
      tmpVal += (lAm1[iWork]+1) / (2.0*Zeta) * twoevRRa0c0( pripair1, pripair2, 
                 FmT_2epri, shell1, shell3, m+1, LA-1,lAm1, LC-1,lCm1 );
    } // if (lA[iWork]>0) 

    if ( lC[iWork]>=2 ){

      lCm1[iWork] -=1; // right now lCm1(iWork) = lC[iWork]-2 
      tmpVal += 0.5 * (lCm1[iWork]+1) * pripair2.one_over_gamma * 
        ( twoevRRa0c0( pripair1, pripair2, FmT_2epri, 
                       shell1, shell3, m, LA,lA, LC-2,lCm1 )

               - ( 1.0/pripair1.one_over_gamma )/Zeta 
          * twoevRRa0c0( pripair1, pripair2, FmT_2epri, 
                         shell1, shell3, m+1, LA,lA, LC-2,lCm1 ) );
    } // if ( lC[iWork]>=2 )

    return tmpVal;
  }  // twoevRRa0c0



  /**
   *  \brief vertical recursion of ERI when all LA > 0, all the others are 0
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
   *  \return ERI of two primitive shell pairs ( shell1 shell2 | shell3 shell4 )
   */
//---------------------------------------------------------------//
// two-e vertical recursion from [a0|00] to [00|00]              //
// [a0|00]^m = (P-A)*[a-1,0|00]^m                                //
//           + (W-P)*[a-1,0|00]^(m+1)                            //
//           + (N(a)-1)/(2*zeta)*[a-2,0|00]^m                    //
//           - (N(a)-1)/(2*zeta)*eta/(zeta+eta)*[a-2,0|00]^(m+1) //
//---------------------------------------------------------------//
   
  double RealGTOIntEngine::twoevRRa000(
    libint2::ShellPair::PrimPairData &pripair1,
    libint2::ShellPair::PrimPairData &pripair2, std::vector<double> &FmT_2epri,
    libint2::Shell &shell1, int m,int LA,int *lA ) {


    if(LA==0) {
      // calculate the (SS||SS) integral  	
      double expoT,Kab,Kcd,SSSS=0.0 ;
 
      // zeta+eta is expoT
      expoT = 1.0/pripair1.one_over_gamma + 1.0/pripair2.one_over_gamma;  
 
      Kab = sqrt(2.0)* pow(M_PI,1.25) * pripair1.K; 
      Kcd = sqrt(2.0)* pow(M_PI,1.25) * pripair2.K;
 
      SSSS += Kab * Kcd * FmT_2epri[m] / sqrt(expoT); 
 
      return SSSS;
 
    } // if LA==0
 
    // here LA != 0
    double tmpVal=0.0,W[3],Zeta;
    int iWork;
    int lAm1[3];
 
    for( iWork=0 ; iWork<3 ; iWork++ ) lAm1[iWork]=lA[iWork];
 
    if (lA[0]>0) iWork=0;
    else if (lA[1]>0) iWork=1;
    else if (lA[2]>0) iWork=2;
 
    if( LA>=1 ) {
 
      lAm1[iWork]-=1;

      Zeta = 1.0/pripair1.one_over_gamma + 1.0/pripair2.one_over_gamma;
      for ( int mu = 0 ; mu<3 ; mu++ )
        W[mu] = (pripair1.P[mu]/pripair1.one_over_gamma 
              + pripair2.P[mu]/pripair2.one_over_gamma )/Zeta ; 
 
      if ( std::abs( W[iWork]- pripair1.P[iWork] ) > 1.0e-15 ) {
 
        tmpVal += ( W[iWork]- pripair1.P[iWork] ) * twoevRRa000( pripair1, pripair2,
          FmT_2epri, shell1, m+1, LA-1,lAm1 );
 
      }  // if ( abs( W_iWork-P[iWork] )>1.0e-15 )

      if ( std::abs( pripair1.P[iWork] - shell1.O[iWork] )>1.0e-15 ) {  
        tmpVal+= ( pripair1.P[iWork] - shell1.O[iWork] ) * twoevRRa000( pripair1, 
          pripair2, FmT_2epri, shell1, m, LA-1,lAm1 );
      } // if ( abs( P[iWork]-A[iWork] )>1.0e-15 ) 

      if ( lA[iWork]>=2 ) {

        lAm1[iWork] -=1; // now lAm1[iWork] == lA[iWork]-2
        tmpVal += 0.5 * ( lAm1[iWork]+1 ) * pripair1.one_over_gamma  
                  *(twoevRRa000( pripair1, pripair2, FmT_2epri,
                    shell1, m, LA-2,lAm1 )
                - 1.0/(pripair2.one_over_gamma*Zeta) * twoevRRa000( pripair1, pripair2,
                  FmT_2epri, shell1, m+1, LA-2,lAm1 ) );
      } // if lA[iWork]>=2

    } // if( LA>=1 ) 

    return tmpVal;

  };  // twoevRRa000



  /**
   *  \brief vertical recursion of ERI when all the angular momentum are 0
   *
   *
   *  \param [in] pair1   bra shell pair data for shell1,shell2
   *  \param [in] pair2   ket shell pair data for shell3,shell4
   *  \param [in] shell1    
   *  \param [in] shell2    
   *  \param [in] shell3    
   *  \param [in] shell4    
   *
   *  \return ERI of two primitive shell pairs ( shell1 shell2 | shell3 shell4 )
   */
   
  double RealGTOIntEngine::twoeSSSS0(
    libint2::ShellPair &pair1, libint2::ShellPair &pair2,
    libint2::Shell &shell1, libint2::Shell &shell2,
    libint2::Shell &shell3, libint2::Shell &shell4 ) {

    // in this auxiliary integral, m=0 
    double norm,sqrPQ,PQ,expoT,T,FmT[1],Kab,Kcd,SSSS0=0.0 ;
    for ( auto &pripair1 : pair1.primpairs )
    for ( auto &pripair2 : pair2.primpairs ) {
       
      sqrPQ = 0.0;
      for(int m=0 ; m<3 ; m++ ) {
        PQ = ( pripair1.P[m]-pripair2.P[m] );
        sqrPQ += PQ*PQ;
      } 

      // zeta+eta is expoT
      expoT = 1.0/pripair1.one_over_gamma + 1.0/pripair2.one_over_gamma;  

      // T= \rho * PQ^2
      T = sqrPQ/( pripair1.one_over_gamma + pripair2.one_over_gamma );

      // calculate F0(T)
      computeFmTTaylor( FmT, T, 0, 0 );

      Kab = sqrt(2.0)* pow(M_PI,1.25) * pripair1.K; 
      Kcd = sqrt(2.0)* pow(M_PI,1.25) * pripair2.K;

      norm = shell1.contr[0].coeff[pripair1.p1]* 
             shell2.contr[0].coeff[pripair1.p2]* 
             shell3.contr[0].coeff[pripair2.p1]* 
             shell4.contr[0].coeff[pripair2.p2];  

      SSSS0 += norm * Kab * Kcd * FmT[0] / sqrt(expoT); 

    } // for pripair2

    return SSSS0;
    
  } // twoeSSSS0


  std::vector<std::vector<double>> RealGTOIntEngine::ACderiv( 
    libint2::ShellPair &pair1, libint2::ShellPair &pair2, 
    libint2::Shell &shell1, libint2::Shell &shell2,
    libint2::Shell &shell3, libint2::Shell &shell4
    ) {


    int cart_i_size,cart_j_size,cart_k_size,cart_l_size; 
    int cart_ip[3], cart_im[3], cart_kp[3], cart_km[3];
    int cart_ip_size, cart_kp_size;
    int cart_im_size, cart_km_size;
    int LA = shell1.contr[0].l;
    int LB = shell2.contr[0].l;
    int LC = shell3.contr[0].l;
    int LD = shell4.contr[0].l;
    

    std::vector<std::vector<double>> ERIgrad_cart(9);
    cart_i_size = cart_ang_list[shell1.contr[0].l].size();
    cart_j_size = cart_ang_list[shell2.contr[0].l].size();
    cart_k_size = cart_ang_list[shell3.contr[0].l].size();
    cart_l_size = cart_ang_list[shell4.contr[0].l].size(); 

    int totalcartsize = cart_i_size*cart_j_size*cart_k_size*cart_l_size;

    // assign memory for outcome
    for (int component = 0 ; component < 9 ; component++ ) {
      ERIgrad_cart[component].assign(cart_i_size*cart_j_size*cart_k_size*cart_l_size,0.0);
    }
    std::vector<std::vector<double>> buffgrad;
//    for ( int k = 0 ; k < 9 ; k++ ) {
//      ERIgrad_sph[k].assign(  (2*LA+1)*(2*LB+1) *(2*LC+1)*(2*LD+1), 0.0  );
//    } 



    // define shells with lA and lC plus minus 1

    std::vector<libint2::Shell> shellsApm(2);
    shellsApm[0] = shell1;
    shellsApm[0].contr[0].l +=1; 
    shellsApm[0].contr[0].pure = 0; // cartesian

    cart_ip_size = (shellsApm[0].contr[0].l+1)*(shellsApm[0].contr[0].l+2)/2;  //NBasis for lA+1

    if  ( shell1.contr[0].l>0 ) {
      shellsApm[1] = shell1;
      shellsApm[1].contr[0].l -=1;  
      shellsApm[1].contr[0].pure = 0; // cartesian
      cart_im_size = (shellsApm[1].contr[0].l+1)*(shellsApm[1].contr[0].l+2)/2; //NBasis for la-1
    }


    std::vector<libint2::Shell> shellsCpm(2);
    shellsCpm[0] = shell3;
    shellsCpm[0].contr[0].l +=1; 
    shellsCpm[0].contr[0].pure = 0; // cartesian
    cart_kp_size = (shellsCpm[0].contr[0].l+1)*(shellsCpm[0].contr[0].l+2)/2;

    if  ( shell3.contr[0].l>0 ) {
      shellsCpm[1] = shell3;
      shellsCpm[1].contr[0].l -=1;  
      shellsCpm[1].contr[0].pure = 0; // cartesian
      cart_km_size = (shellsCpm[1].contr[0].l+1)*(shellsCpm[1].contr[0].l+2)/2;
    }
    
    // set pure to cartesian
    libint2::Shell shellB = shell2;
    shellB.contr[0].pure = 0; // cartesian

    libint2::Shell shellD = shell4;
    shellD.contr[0].pure = 0; // cartesian

    // compute integrals with different angular momentum

    // lA+1 lB+1
    //std::vector<double> buffpp  = RealGTOIntEngine::BottomupHGP(pair1,pair2,
    std::vector<double> buffpp  = RealGTOIntEngine::bottomupERI(pair1,pair2,
      shellsApm[0],
      shellB,
      shellsCpm[0],
      shellD,
      1, 0, 1, 0
    );

  //for (int ii = 0 ; ii < cart_ip_size*cart_j_size*cart_kp_size*cart_l_size ; ii++ ) {
  //  std::cout<<"ii = "<<ii<<" pp = "<<buffpp[ii]<<std::endl;
  //}

    buffgrad.push_back(buffpp);

    if  ( shell1.contr[0].l>0 ) {
        //auto buffmp  = RealGTOIntEngine::BottomupHGP(pair1,pair2,
        auto buffmp  = RealGTOIntEngine::bottomupERI(pair1,pair2,
          shellsApm[1],
          shellB,
          shellsCpm[0],
          shellD,
          0, 0, 1, 0
        );
  //for (int ii = 0 ; ii < cart_im_size*cart_j_size*cart_kp_size*cart_l_size ; ii++ ) {
  //  std::cout<<"ii = "<<ii<<" mp = "<<buffmp[ii]<<std::endl;
  //}
  /*
        auto testbuffpm  = RealGTOIntEngine::bottomupERI(pair1,pair2,
          shellsCpm[0],
          shellB,
          shellsApm[1],
          shellD,
          1, 0, 0, 0
        );
  std::cout<<"print testpm"<<std::endl;
  for (int ii = 0 ; ii < cart_im_size*cart_j_size*cart_kp_size*cart_l_size ; ii++ ) {
    std::cout<<"ii = "<<ii<<" testpm = "<<testbuffpm[ii]<<std::endl;
  }
*/
        buffgrad.push_back(buffmp);

    }  else {

        std::vector<double> place;
        buffgrad.push_back(place);
    }

    if  ( shell3.contr[0].l>0 ) {

        //auto buffpm  = RealGTOIntEngine::BottomupHGP(pair1,pair2,
        auto buffpm  = RealGTOIntEngine::bottomupERI(pair1,pair2,
          shellsApm[0],
          shellB,
          shellsCpm[1],
          shellD,
          1, 0, 0, 0
        );
/*
  for (int ii = 0 ; ii < cart_ip_size*cart_j_size*cart_km_size*cart_l_size ; ii++ ) {
    std::cout<<"ii = "<<ii<<" pm = "<<buffpm[ii]<<std::endl;
  }
*/
        buffgrad.push_back(buffpm);

    }  else {
        std::vector<double> place;
        buffgrad.push_back(place);
    }
    

    if ( shell1.contr[0].l>0  and  shell3.contr[0].l>0 )  { 
        //auto buffmm  = RealGTOIntEngine::BottomupHGP(pair1,pair2,
        auto buffmm  = RealGTOIntEngine::bottomupERI(pair1,pair2,
          shellsApm[1],
          shellB,
          shellsCpm[1],
          shellD,
          0, 0, 0, 0
        );
  //for (int ii = 0 ; ii < cart_im_size*cart_j_size*cart_km_size*cart_l_size ; ii++ ) {
  //  std::cout<<"ii = "<<ii<<" mm = "<<buffmm[ii]<<std::endl;
  //}
        buffgrad.push_back(buffmm);

/*
for (int cart_i = 0; cart_i < cart_im_size ; cart_i++) {
  for (int cart_j = 0; cart_j < cart_j_size ; cart_j++) {
    for (int cart_k = 0; cart_k < cart_km_size ; cart_k++) {
      for (int cart_l = 0; cart_l < cart_l_size ; cart_l++) {
        std::cout<<"minus minus: i="<<cart_i<<" j= "<<cart_j<<" k= "<<cart_k<<" l= "<<cart_l<<" "<<buffgrad[3][ cart_i *cart_j_size*cart_km_size*cart_l_size+
              cart_j * cart_km_size*cart_l_size +cart_k* cart_l_size
              + cart_l]<<std::endl;

         
      }
    }
  }
}
*/ 

    } else {
        std::vector<double> place;
        buffgrad.push_back(place);
    }

    // combine into derivatives
    //

    for (int cart_i = 0; cart_i < cart_i_size ; cart_i++) { 
      int lA_xyz[3];     

      for (int ii = 0; ii < 3; ii++) { 

        lA_xyz[ii] = cart_ang_list[shell1.contr[0].l][cart_i][ii];

      } // for (int ii = 0; ii < 3; ii++)
      // here find out the cart_ip ( cart i of lA_x + 1 ) 
      cart_ip[0] = indexmap(shellsApm[0].contr[0].l, lA_xyz[0]+1, lA_xyz[1], lA_xyz[2] );
      cart_ip[1] = indexmap(shellsApm[0].contr[0].l, lA_xyz[0], lA_xyz[1]+1, lA_xyz[2] );
      cart_ip[2] = indexmap(shellsApm[0].contr[0].l, lA_xyz[0], lA_xyz[1], lA_xyz[2]+1 );


      if  ( lA_xyz[0]>0 ) {  
        cart_im[0] = indexmap(shellsApm[1].contr[0].l, lA_xyz[0]-1, lA_xyz[1], lA_xyz[2] ); 
      }

      if  ( lA_xyz[1]>0 ) {  
        cart_im[1] = indexmap(shellsApm[1].contr[0].l, lA_xyz[0], lA_xyz[1]-1, lA_xyz[2] ); 
      }

      if  ( lA_xyz[2]>0 ) {  
        cart_im[2] = indexmap(shellsApm[1].contr[0].l, lA_xyz[0], lA_xyz[1], lA_xyz[2]-1 ); 
      }  
       
      for (int cart_j = 0; cart_j < cart_j_size; cart_j++) { 
        int lB_xyz[3];     

        for (int ii = 0; ii < 3; ii++) { 

          lB_xyz[ii] = cart_ang_list[shell2.contr[0].l][cart_j][ii];

        } // for (int ii = 0; ii < 3; ii++)


        for (int cart_k = 0; cart_k < cart_k_size; cart_k++) { 
          int lC_xyz[3];     

          for (int ii = 0; ii < 3; ii++) { 

            lC_xyz[ii] = cart_ang_list[shell3.contr[0].l][cart_k][ii];
          } // for (int ii = 0; ii < 3; ii++)
          // here find out the cart_kp ( cart i of lC_x + 1 ) 
          cart_kp[0] = indexmap(shellsCpm[0].contr[0].l, lC_xyz[0]+1, lC_xyz[1], lC_xyz[2] );
          cart_kp[1] = indexmap(shellsCpm[0].contr[0].l, lC_xyz[0], lC_xyz[1]+1, lC_xyz[2] );
          cart_kp[2] = indexmap(shellsCpm[0].contr[0].l, lC_xyz[0], lC_xyz[1], lC_xyz[2]+1 );
          
          if  ( lC_xyz[0]>0 ) {  
            cart_km[0] = indexmap(shellsCpm[1].contr[0].l, lC_xyz[0]-1, lC_xyz[1], lC_xyz[2] ); 
          }
          if  ( lC_xyz[1]>0 ) {  
            cart_km[1] = indexmap(shellsCpm[1].contr[0].l, lC_xyz[0], lC_xyz[1]-1, lC_xyz[2] ); 
          }
          if  ( lC_xyz[2]>0 ) {  
            cart_km[2] = indexmap(shellsCpm[1].contr[0].l, lC_xyz[0], lC_xyz[1], lC_xyz[2]-1 ); 
          }  

          for (int cart_l = 0; cart_l < cart_l_size; cart_l++) { 
            int lD_xyz[3];     

            for (int ii = 0; ii < 3; ii++) { 

              lD_xyz[ii] = cart_ang_list[shell4.contr[0].l][cart_l][ii];
            } // for ii
/*
            // (1)
            // place the lAx+1 and lCx+1 into integral
            ERIgrad_cart[cart_i*cart_j_size*cart_k_size*cart_l_size+
              cart_j * cart_k_size*cart_l_size + cart_k* cart_l_size
              + cart_l]
            += 4*buffgrad[0][ cart_ipx *cart_j_size*cart_kp_size*cart_l_size+
              cart_j * cart_kp_size*cart_l_size +cart_kpx* cart_l_size
              + cart_l];

            // lAx-1, lCx+1
            // if  ( basisSet.shells[s1].contr[0].l>0 ) {
            if  ( lA_xyz[0]>0 ) {
              ERIgrad_cart[cart_i*cart_j_size*cart_k_size*cart_l_size+
                cart_j * cart_k_size*cart_l_size + cart_k* cart_l_size
                + cart_l]
                -= 2*lA_xyz[0]*buffgrad[1][ cart_imx *cart_j_size*cart_kp_size*cart_l_size+
                cart_j * cart_kp_size*cart_l_size +cart_kpx* cart_l_size
                + cart_l];
            }

            // lAx+1, lCx-1
            if  ( lC_xyz[0]>0 ) {
              ERIgrad_cart[cart_i*cart_j_size*cart_k_size*cart_l_size+
                cart_j * cart_k_size*cart_l_size + cart_k* cart_l_size
                + cart_l]
                -= 2*lC_xyz[0]*buffgrad[2][ cart_ipx *cart_j_size*cart_km_size*cart_l_size+
                cart_j * cart_km_size*cart_l_size +cart_kmx* cart_l_size
                + cart_l];
            }

            // lAx-1, lCx-1
            if  (lA_xyz[0]>0  and  lC_xyz[0]>0 ) {
              ERIgrad_cart[cart_i*cart_j_size*cart_k_size*cart_l_size+
                cart_j * cart_k_size*cart_l_size + cart_k* cart_l_size
                + cart_l]
                += lA_xyz[0]*lC_xyz[0]*buffgrad[3][ cart_imx *cart_j_size*cart_km_size*cart_l_size+
                cart_j * cart_km_size*cart_l_size +cart_kmx* cart_l_size
                + cart_l];
            }
   */ 
            // (2)
            for ( int aidx = 0 ; aidx < 3 ; aidx++ ) {
            for ( int cidx = 0 ; cidx < 3 ; cidx++ ) {
              int acidx = aidx*3+cidx;
            // place the lA+1 and lC+1 into integral
            ERIgrad_cart[acidx][ cart_i*cart_j_size*cart_k_size*cart_l_size+
              cart_j * cart_k_size*cart_l_size + cart_k* cart_l_size
              + cart_l]
            += 4*buffgrad[0][ cart_ip[aidx] *cart_j_size*cart_kp_size*cart_l_size+
              cart_j * cart_kp_size*cart_l_size +cart_kp[cidx]* cart_l_size
              + cart_l];

              // lA-1, lC+1
              // if  ( basisSet.shells[s1].contr[0].l>0 ) {
              if  ( lA_xyz[aidx]>0 ) {
                ERIgrad_cart[acidx][ cart_i*cart_j_size*cart_k_size*cart_l_size+
                  cart_j * cart_k_size*cart_l_size + cart_k* cart_l_size
                  + cart_l]
                  -= 2*lA_xyz[aidx]*buffgrad[1][ cart_im[aidx] *cart_j_size*cart_kp_size*cart_l_size+
                  cart_j * cart_kp_size*cart_l_size +cart_kp[cidx]* cart_l_size
                  + cart_l];
              }
             
              // lA+1, lC-1
              if  ( lC_xyz[cidx]>0 ) {
                ERIgrad_cart[acidx][ cart_i*cart_j_size*cart_k_size*cart_l_size+
                  cart_j * cart_k_size*cart_l_size + cart_k* cart_l_size
                  + cart_l]
                  -= 2*lC_xyz[cidx]*buffgrad[2][ cart_ip[aidx] *cart_j_size*cart_km_size*cart_l_size+
                  cart_j * cart_km_size*cart_l_size +cart_km[cidx]* cart_l_size
                  + cart_l];
              }
             
              // lA-1, lC-1
              if  (lA_xyz[aidx]>0  and  lC_xyz[cidx]>0 ) {
                ERIgrad_cart[acidx][ cart_i*cart_j_size*cart_k_size*cart_l_size+
                  cart_j * cart_k_size*cart_l_size + cart_k* cart_l_size
                  + cart_l]
                  += lA_xyz[aidx]*lC_xyz[cidx]*buffgrad[3][ cart_im[aidx] *cart_j_size*cart_km_size*cart_l_size+
                  cart_j * cart_km_size*cart_l_size +cart_km[cidx]* cart_l_size
                  + cart_l];
              }
            } // cidx
            } // aidx 

          }  // for cart_l

        } // for cart_k
      } // for cart_j
    } // for cart_i 

/*
  for (int ii = 0 ; ii < totalcartsize ; ii++ ) {
    std::cout<<"ii = "<<ii<<" ERIgradsph = "<<ERIgrad_cart[0][ii]<<std::endl;
  }
*/
     
    if ( ( not shell1.contr[0].pure ) and ( not shell2.contr[0].pure ) and 
         ( not shell3.contr[0].pure ) and ( not shell4.contr[0].pure ) ) {  
      // if both sides are cartesian, return cartesian gaussian integrals
      return ERIgrad_cart;
    }

    std::vector<std::vector<double>> ERIgrad_sph(9);
    int k = 0;

    for ( auto cartmatrix : ERIgrad_cart ) {
      ERIgrad_sph[k].assign(  (2*LA+1)*(2*LB+1)*(2*LC+1)*(2*LD+1), 0.0  );


      cart2sph_2e_transform(LA,LB,LC,LD, ERIgrad_sph[k],cartmatrix );
      k++;

    }   
  

    return ERIgrad_sph; 




  } // ACderiv


 // this part is gauge interaction

 
  /**
   *  \brief Compute the ERI of two shell pairs
   *
   *
   *  \param [in] pair1  bra shell pair data for shell1,shell2
   *  \param [in] pair2  ket shell pair data for shell3,shell4
   *  \param [in] shell1
   *  \param [in] shell2
   *  \param [in] shell3
   *  \param [in] shell4
   *
   *  \return ERI of two shell pairs ( shell1 shell2 | shell3 shell4 )
   */
   
  std::vector<std::vector<double>> RealGTOIntEngine::computegaugeabcd(libint2::ShellPair &pair1 ,
    libint2::ShellPair &pair2, libint2::Shell &shell1, libint2::Shell &shell2,
    libint2::Shell &shell3, libint2::Shell &shell4, 
    int deria, int derib, int deric, int derid )  {
   
//std::cout<<"computegaugeabcd"<<std::endl;
 
    // six components 
    
    double tmpVal=0.0,sqrPQ,PQ;
    std::vector<std::vector<double>> ERI_cart(6);
    
    int lA[3],lB[3],lC[3],lD[3],p,q;

    // compute total angular momentum
    auto lTotal = shell1.contr[0].l + shell2.contr[0].l
                + shell3.contr[0].l + shell4.contr[0].l;

/*
    std::cerr<<"LA "<<shell1.contr[0].l
    <<" LB "<<shell2.contr[0].l
    <<" LC "<<shell3.contr[0].l 
    <<" LD "<<shell4.contr[0].l<<std::endl;
*/

    // pre calculate all the Boys functions 
    // dimension is FmT_2e[shellpair1.prim][shellpair2.prim][lTotal+1]
    std::vector<std::vector<std::vector<double>>> FmT_2e;
    FmT_2e.resize(pair1.primpairs.size());

    double *FmT = new double[lTotal+3];  // 0 to l+2
    int shellpair1_i=0, shellpair2_j ;
    for ( auto &pripair1 : pair1.primpairs ) {
      FmT_2e[shellpair1_i].resize( pair2.primpairs.size() );

      shellpair2_j = 0 ; 
      for ( auto &pripair2 : pair2.primpairs ) {
        sqrPQ = 0.0;
        for ( int mu=0 ; mu<3 ; mu++ ) {
          PQ = ( pripair1.P[mu]-pripair2.P[mu] ); 
          sqrPQ += PQ*PQ;
        }
        auto Zeta = 1.0/pripair1.one_over_gamma;
        auto Eta  = 1.0/pripair2.one_over_gamma;
        
        auto rho = Zeta*Eta/(Zeta+Eta);
        auto T = rho*sqrPQ;
        // calculate Fm(T) list
//std::cout<<"before computetaylog"<<std::endl;
        computeFmTTaylor( FmT, T, lTotal+2, 0 );
//std::cout<<"after computetaylor"<<std::endl;
//  multiply the derivative coeffs

        for ( int lcurr = 0 ; lcurr < lTotal+3 ; lcurr++ ) {
           if ( std::abs(FmT[lcurr]) < 1.0e-15 ) 
             FmT_2e[shellpair1_i][shellpair2_j].push_back(0.0);
           else {
             FmT[lcurr] *= pow( shell1.alpha[pripair1.p1], deria); 
/*
             std::cout<<"lcurr"<<lcurr<<" ltotal "<<lTotal<<std::endl;
             std::cout<<"FmT[lcurr]"<<FmT[lcurr]<<std::endl;
             std::cout<<"shell2.alpha "<<shell2.alpha.size()<<" pripair1.p2 "<<pripair1.p2<<
               "shell1.alpha "<<shell1.alpha.size()<<" pripair1.p1 "<<pripair1.p1<<std::endl;
             std::cout<<" shell2.alpha[pripair1.p2]"<<shell2.alpha[pripair1.p2]<<std::endl;
             std::cout<<"FmT[lcurr]"<<FmT[lcurr]<<" shell2.alpha[pripair1.p2]"<<shell2.alpha[pripair1.p2]<<" derib "<<derib<<std::endl;
*/
             FmT[lcurr] *= pow( shell2.alpha[pripair1.p2], derib); 
//             std::cout<<"FmT[lcurr]"<<FmT[lcurr]<<" shell2.alpha[pripair2.p1]"<<shell3.alpha[pripair2.p1]<<" deric "<<deric<<std::endl;
             FmT[lcurr] *= pow( shell3.alpha[pripair2.p1], deric); 
             FmT[lcurr] *= pow( shell4.alpha[pripair2.p2], derid); 
             FmT_2e[shellpair1_i][shellpair2_j].push_back(FmT[lcurr]);
//             std::cout<<"push back FmT"<<std::endl;
           }  
        } // for lcurr
        shellpair2_j ++;
      } // for pripair2
    shellpair1_i++;
    } // for pripair1
    delete[] FmT;

//std::cout<<"shell size= "<<cart_ang_list[shell1.contr[0].l].size()<<std::endl;

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

      for (int mu = 0 ; mu < 6 ; mu++ ) { 
       
/* 
        if ( mu == 0 )  {
          p = 0;
          q = 0;
        }
        if (mu == 1) {
          p = 0;
          q = 1;
        }
        if (mu==2) {
          p = 0;
          q = 2;
        } 
        if (mu = 3 ) {
          p = 1;
          q = 1;
        }
        if (mu==4) {
          p=1;
          q = 2;
        }
        if (mu ==5) {
          p = 2;
          q = 2;
        }   
*/
        if ( mu < 3) {
          p = 0;
          q = mu;
        } 
        if ( mu == 5 ) {
          p = 2;
          q = 2;
        }
        if ( mu > 2 and mu < 5 ) {
          p = 1;
          q = mu-2;
        }
        
// std::cout<<"mu= "<<mu<<" p= "<<p<<" q= "<<q<<std::endl;

        tmpVal = twoegaugehRRabcd(pair1,pair2,p,q,shell1,shell2,shell3,shell4,
                 FmT_2e,shell1.contr[0].l, lA, shell2.contr[0].l, lB, 
                 shell3.contr[0].l, lC, shell4.contr[0].l, lD ); 
     
        // if (std::abs(tmpVal)>1.0e-6) {
          ERI_cart[mu].push_back(tmpVal);
        // } else {
        //   ERI_cart[mu].push_back(0.0);
        // }

      }  // for mu
    }   // for l


    if ( ( not shell1.contr[0].pure ) and ( not shell2.contr[0].pure ) and 
         ( not shell3.contr[0].pure ) and ( not shell4.contr[0].pure ) ) {  
      // if both sides are cartesian, return cartesian gaussian integrals
      return ERI_cart;
    }

    std::vector<std::vector<double>> ERI_sph(6);

    for ( int mu = 0 ; mu < 6 ; mu++ ) { 
      ERI_sph[mu].assign(((2*shell1.contr[0].l+1)*(2*shell2.contr[0].l+1)
                   *(2*shell3.contr[0].l+1)*(2*shell4.contr[0].l+1)),0.0); 

      cart2sph_2e_transform( shell1.contr[0].l,shell2.contr[0].l,
        shell3.contr[0].l,shell4.contr[0].l,ERI_sph[mu],ERI_cart[mu] );
    }

    return ERI_sph; 

  }   // computeERIabcd
//


// here is the hRR of gauge integral
//

  /**
   *  \brief horizontal recursion of gauge ERI when all angular momentum are nonzero
   *
   *
   *  \param [in] pair1  bra shell pair data for shell1,shell2
   *  \param [in] pair2  ket shell pair data for shell3,shell4
   *  \param [in] p
   *  \param [in] q
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
   *  \return ERI of two shell pairs ( shell1 shell2 | shell3 shell4 )
   */
  //----------------------------------------------------//
  // two-e horizontal recursion from (ab|cd) to (a0|cd) //
  // (ab|cd)=(a+1,b-1|cd)+(A-B)*(a,b-1|cd)              //
  //----------------------------------------------------//
   
  double RealGTOIntEngine::twoegaugehRRabcd( 
     libint2::ShellPair &pair1 ,libint2::ShellPair &pair2 ,
     int p, int q, 
     libint2::Shell &shell1, libint2::Shell &shell2,
     libint2::Shell &shell3, libint2::Shell &shell4,
     std::vector<std::vector<std::vector<double>>> &FmT_2e,
     int LA,int *lA,int LB,int *lB,int LC,int *lC,int LD,int *lD) {

     double tmpVal = 0.0, tmpVal1=0.0;

  // iWork is used to indicate which Cartesian angular momentum we are reducing (x,y,z)

    int iWork;
    int totalL = LA + LB + LC + LD;

//std::cout<<"totalL= "<<totalL<<std::endl;

    if(totalL==0) {   // (SS|1p1q|SS) type 

      return twoegaugeSSSS0(pair1,pair2,p,q,shell1,shell2,shell3,shell4);

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

          tmpVal += twoegaugehRRabcd(pair1,pair2,p,q,shell1,shell2,shell3,shell4,
                      FmT_2e, LA+1,lAp1,LB-1,lBm1,LC,lC,LD,lD);

          if ( std::abs(pair1.AB[iWork]) > 1.0e-15 ) {
            tmpVal += pair1.AB[iWork]*twoegaugehRRabcd(pair1,pair2,p,q,shell1,shell2,shell3,
                        shell4,FmT_2e, LA,lA,LB-1,lBm1,LC,lC,LD,lD);
          } 

        } else if ( LB == 0 ) {
          tmpVal = twoegaugehRRa0cd(pair1,pair2,p,q,shell1,shell2,shell3,shell4,FmT_2e,
                                 LA,lA,LC,lC,LD,lD);

        } // LB == 0
      }  // else ( LTOTAL != 0 )

      return tmpVal;

  }  // twoegaugehRRabcd





  /**
   *  \brief horizontal recursion of ERI when all LB=0
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
   *  \return ERI of two shell pairs ( shell1 shell2 | shell3 shell4 )
   */
  //----------------------------------------------------//
  // two-e horizontal recursion from (a0|cd) to (a0|c0) //
  // (a0|cd)=(a,0|c+1,d-1)+(C-D)*(a,0|c,d-1)            //
  //----------------------------------------------------//
   
  double RealGTOIntEngine::twoegaugehRRa0cd(
    libint2::ShellPair &pair1, libint2::ShellPair &pair2,
    int p, int q,
    libint2::Shell &shell1, libint2::Shell &shell2,
    libint2::Shell &shell3, libint2::Shell &shell4, 
    std::vector<std::vector<std::vector<double>>> &FmT_2e,
    int LA,int *lA,int LC,int *lC,int LD,int *lD)  {

    double tmpVal=0.0;
    if(LD==0) {
      int pair1index=0, pair2index=0;
      // go into the vertical recursion
      for ( auto &pripair1 : pair1.primpairs ) {
        pair2index = 0; 
        for ( auto &pripair2 : pair2.primpairs ) {

          auto norm = 
                 shell1.contr[0].coeff[pripair1.p1]* 
                 shell2.contr[0].coeff[pripair1.p2]* 
                 shell3.contr[0].coeff[pripair2.p1]* 
                 shell4.contr[0].coeff[pripair2.p2];  

          tmpVal +=  norm * twoegaugevRRa0c0( pripair1, pripair2, p, q,  
             FmT_2e[pair1index][pair2index], shell1,shell3, 0, LA,lA,LC,lC);

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
      tmpVal = twoegaugehRRa0cd(pair1, pair2, p, q, shell1, shell2, shell3, shell4, 
                             FmT_2e,LA,lA, LC+1,lCp1, LD-1,lDm1 ); 
      if ( std::abs(pair2.AB[iWork]) > 1.0e-15 ){
        tmpVal += pair2.AB[iWork] * twoegaugehRRa0cd( pair1, pair2, p, q, shell1, shell2, 
               shell3, shell4, FmT_2e, LA,lA, LC,lC, LD-1,lDm1 );
      }

    } // else ( that means LD > 0 )
    return tmpVal;


  }  // twoegaugehRRa0cd



  /**
   *  \brief vertical recursion of ERI when all LA, LC > 0
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
   *  \return ERI of two primitive shell pairs ( shell1 shell2 | shell3 shell4 )
   */


   
  double RealGTOIntEngine::twoegaugevRRa0c0(
    libint2::ShellPair::PrimPairData &pripair1,
    libint2::ShellPair::PrimPairData &pripair2, 
    int p, int q,
    std::vector<double> &FmT_2epri, 
    libint2::Shell &shell1, libint2::Shell &shell3,
    int m, int LA, int *lA, int LC, int *lC ) {
 
 
    if(LC==0) return twoegaugevRRa000( pripair1, pripair2, p, q, FmT_2epri,
                                  shell1, shell3, m, LA, lA );

    // define pq vector
    int pq[3],pqm1[3];
    pq[0]=0;
    pq[1]=0;
    pq[2]=0;
    pq[p]+=1;
    pq[q]+=1;
 
    int lAm1[3],lCm1[3],lAp1[3],lCpm1[3];  
 
    for ( int iWork=0 ; iWork<3 ; iWork++ ){
      lAm1[iWork]=lA[iWork];     
      lAp1[iWork]=lA[iWork];     
      lCm1[iWork]=lC[iWork];
      lCpm1[iWork]=lC[iWork];
    }
 
    double tmpVal=0.0;
    double W_iWork,Zeta;
    int iWork,jWork;
 
    if (lC[0]>0) iWork=0;
    else if (lC[1]>0) iWork=1;
    else if (lC[2]>0) iWork=2;

    lCm1[iWork]-=1;

    if ( std::abs(pripair2.P[iWork] - shell3.O[iWork]) > 1.0e-15 ) {
      tmpVal += ( pripair2.P[iWork] - shell3.O[iWork] ) * 
                twoegaugevRRa0c0( pripair1, pripair2, p, q, FmT_2epri, 
                  shell1, shell3, m, LA,lA, LC-1,lCm1 );
    }  // if (Q-C) > 1.0e-15 


    Zeta = 1.0/pripair1.one_over_gamma + 1.0/pripair2.one_over_gamma; 
/*
      for ( int mu = 0 ; mu<3 ; mu++ )
      W[mu] = (pripair1.P[mu]/pripair1.one_over_gamma 
              + pripair2.P[mu]/pripair2.one_over_gamma )/Zeta ; 
*/

    W_iWork = ( pripair1.P[iWork]/pripair1.one_over_gamma + 
               pripair2.P[iWork]/pripair2.one_over_gamma )/Zeta;

    if( std::abs( W_iWork- pripair2.P[iWork] ) > 1.0e-15 ) {
      tmpVal += ( W_iWork- pripair2.P[iWork] ) * twoegaugevRRa0c0( pripair1, pripair2, p, q,
                 FmT_2epri, shell1,shell3, m+1, LA,lA, LC-1,lCm1 );
    } // if( abs( W_iWork- Q[iWork] ) > 1.0e-15 )



    // here needs some normal ERI 
    if ( pq[iWork] > 0 ) {  
      for (int jj = 0 ; jj < 3 ; jj++)
        pqm1[jj]=pq[jj];
       
      pqm1[iWork] -=1;
      for (int jj = 0 ; jj < 3 ; jj++) {
        if ( pqm1[jj]>0 ) 
          jWork = jj;
        
      }
      lAp1[jWork] += 1;
      lCpm1[iWork]-= 1;
      lCpm1[jWork]+= 1;
      // rho/eta = zeta/Zeta
      tmpVal -= ( 1.0/pripair1.one_over_gamma )/Zeta * pq[iWork] *  ( 
        twoevRRa0c0( pripair1, pripair2, FmT_2epri, 
                       shell1, shell3, m+1, LA+1,lAp1, LC-1,lCm1 ) -
        twoevRRa0c0( pripair1, pripair2, FmT_2epri,shell1, shell3, m+1, LA,lA, LC,lCpm1 ) 
      );
      if( std::abs( shell1.O[jWork]- shell3.O[jWork] ) > 1.0e-15 ) {
        tmpVal -= ( 1.0/pripair1.one_over_gamma )/Zeta * pq[iWork] *( shell1.O[jWork]- shell3.O[jWork] )
        * twoevRRa0c0( pripair1, pripair2, FmT_2epri,shell1, shell3, m+1, LA,lA, LC-1,lCm1 );
      }  //(A-C)jWork not 0
    }   //  ( pq[iWork] > 0 )

    // 0.5*2 = 1


    if (lA[iWork]>0) { 

      lAm1[iWork] -= 1;
      tmpVal += (lAm1[iWork]+1) / (2.0*Zeta) * twoegaugevRRa0c0( pripair1, pripair2, p, q, 
                 FmT_2epri, shell1, shell3, m+1, LA-1,lAm1, LC-1,lCm1 );
    } // if (lA[iWork]>0) 

    if ( lC[iWork]>1 ){

      lCm1[iWork] -=1; // right now lCm1(iWork) = lC[iWork]-2 
      tmpVal += 0.5 * (lCm1[iWork]+1) * pripair2.one_over_gamma * 
        ( twoegaugevRRa0c0( pripair1, pripair2, p, q, FmT_2epri, 
                       shell1, shell3, m, LA,lA, LC-2,lCm1 )

               - ( 1.0/pripair1.one_over_gamma )/Zeta 
          * twoegaugevRRa0c0( pripair1, pripair2, p, q, FmT_2epri, 
                         shell1, shell3, m+1, LA,lA, LC-2,lCm1 ) );
    } // if ( lC[iWork]>1 )
    return tmpVal;
  }  // twoegaugevRRa0c0



   
  double RealGTOIntEngine::twoegaugevRRa000(
    libint2::ShellPair::PrimPairData &pripair1,
    libint2::ShellPair::PrimPairData &pripair2, int p, int q, std::vector<double> &FmT_2epri,
    libint2::Shell &shell1, libint2::Shell &shell3, int m,int LA,int *lA ) {

    
    double PQ[3]; 
    for(int m=0 ; m<3 ; m++ ) {
      PQ[m] = ( pripair1.P[m]-pripair2.P[m] );
    } 
    
    if(LA==0) {
      // calculate the (SS||SS) integral  	
      double expoT,Kab,Kcd,SSSS=0.0 ;
 
      // zeta+eta is expoT
      expoT = 1.0/pripair1.one_over_gamma + 1.0/pripair2.one_over_gamma;  
 
      Kab = sqrt(2.0)* pow(M_PI,1.25) * pripair1.K; 
      Kcd = sqrt(2.0)* pow(M_PI,1.25) * pripair2.K;

      if (p==q) {
        SSSS += 0.5* // (pripair1.one_over_gamma + pripair2.one_over_gamma) * 
           Kab * Kcd * 
          FmT_2epri[m+1] / sqrt(expoT)    *2 ;  // rho cancel with 1/zeta+1/eta 
      }

      SSSS += PQ[p]*PQ[q]*2/( pripair1.one_over_gamma + pripair2.one_over_gamma ) // \rho 
                 *Kab * Kcd * (  FmT_2epri[m+1]    -  FmT_2epri[m+2]) / sqrt(expoT);

 
      return SSSS;
 
    } // if LA==0
 

    // define pq vector
    int pq[3],pqm1[3];
    pq[0]=0;
    pq[1]=0;
    pq[2]=0;
    pq[p]+=1;
    pq[q]+=1;
    for (int jj = 0 ; jj < 3 ; jj++) {
      pqm1[jj]=pq[jj];
    } 
 
    // here LA != 0
    double tmpVal=0.0,W[3],Zeta;
    int iWork,jWork;
    int lAm1[3],lApm1[3],lCp1[3],lC[3];
 
    for( iWork=0 ; iWork<3 ; iWork++ ) {
      lAm1[iWork]=lA[iWork];
      lApm1[iWork]=lA[iWork];
      lCp1[iWork]=0;
      lC[iWork]=0;
    }  
    int LC = 0; 
 
    if (lA[0]>0) iWork=0;
    else if (lA[1]>0) iWork=1;
    else if (lA[2]>0) iWork=2;
 
    if( LA>=1 ) {
 
      lAm1[iWork]-=1;

      // P-A
      if ( std::abs( pripair1.P[iWork] - shell1.O[iWork] )>1.0e-15 ) {  
        tmpVal+= ( pripair1.P[iWork] - shell1.O[iWork] ) * twoegaugevRRa000( pripair1, 
          pripair2, p, q, FmT_2epri, shell1, shell3, m, LA-1,lAm1 );
      } // if ( abs( P[iWork]-A[iWork] )>1.0e-15 ) 

      Zeta = 1.0/pripair1.one_over_gamma + 1.0/pripair2.one_over_gamma;
      for ( int mu = 0 ; mu<3 ; mu++ )
        W[mu] = (pripair1.P[mu]/pripair1.one_over_gamma 
              + pripair2.P[mu]/pripair2.one_over_gamma )/Zeta ; 
 
      if ( std::abs( W[iWork]- pripair1.P[iWork] ) > 1.0e-15 ) {
 
        tmpVal += ( W[iWork]- pripair1.P[iWork] ) * twoegaugevRRa000( pripair1, pripair2,
          p, q, FmT_2epri, shell1, shell3, m+1, LA-1, lAm1 );
 
      }  // if ( abs( W_iWork-P[iWork] )>1.0e-15 )


      // here needs some normal ERI 
      if ( pq[iWork] > 0 ) {  
         
        pqm1[iWork] -=1;
        for (int jj = 0 ; jj < 3 ; jj++) {
          if ( pqm1[jj]>0 ) {
            jWork = jj;
          }
        } // for jj 

        lCp1[jWork] += 1;
        lApm1[iWork]-= 1;
        lApm1[jWork]+= 1;
//std::cout<<"iwork= "<<iWork<<" jWork= "<<jWork<<" lapm1[0,1,2]= "<<lApm1[0]<<" "<<lApm1[1]<<" "<<lApm1[2]<<std::endl;
        // rho/zeta = eta/Zeta

        tmpVal += ( 1.0/pripair2.one_over_gamma )/Zeta * pq[iWork] *  ( 
          twoevRRa0c0( pripair1, pripair2, FmT_2epri,shell1, shell3, m+1, LA,lApm1, LC,lC ) //this part might change to vRRa000
          - twoevRRa0c0( pripair1, pripair2, FmT_2epri, shell1, shell3, m+1, LA-1,lAm1, LC+1,lCp1 ) 
        );
        if( std::abs( shell1.O[jWork]- shell3.O[jWork] ) > 1.0e-15 ) {
          tmpVal += ( 1.0/pripair2.one_over_gamma )/Zeta * pq[iWork] *( shell1.O[jWork]- shell3.O[jWork] )
          * twoevRRa0c0( pripair1, pripair2, FmT_2epri,shell1, shell3, m+1, LA-1,lAm1, LC,lC ); // this part might change to vRRa000
        }  //(A-C)jWork not 0

/*
        double firstnormal = twoevRRa0c0( pripair1, pripair2, FmT_2epri,shell1, shell3, m+1, LA,lApm1, LC,lC );
        double firstnormalprime = twoevRRa0c0( pripair1, pripair2, FmT_2epri,shell1, shell3, m+1, LA,lA, LC,lC );
        double secondnormal = twoevRRa0c0( pripair1, pripair2, FmT_2epri, shell1, shell3, m+1, LA-1,lAm1, LC+1,lCp1 )  ;
        double thirdnormal = twoevRRa0c0( pripair1, pripair2, FmT_2epri,shell1, shell3, m+1, LA-1,lAm1, LC,lC );

        double pqiwork = ( 1.0/pripair2.one_over_gamma )/Zeta * pq[iWork] *  ( 
          twoevRRa0c0( pripair1, pripair2, FmT_2epri,shell1, shell3, m+1, LA,lApm1, LC,lC ) //this part might change to vRRa000
          - twoevRRa0c0( pripair1, pripair2, FmT_2epri, shell1, shell3, m+1, LA-1,lAm1, LC+1,lCp1 ) 
        );
        if( std::abs( shell1.O[jWork]- shell3.O[jWork] ) > 1.0e-15 ) {
          pqiwork += ( 1.0/pripair2.one_over_gamma )/Zeta * pq[iWork] *( shell1.O[jWork]- shell3.O[jWork] )
          * twoevRRa0c0( pripair1, pripair2, FmT_2epri,shell1, shell3, m+1, LA-1,lAm1, LC,lC ); // this part might change to vRRa000
        }  //(A-C)jWork not 0
        std::cout<<"W[iwork]= "<<W[iWork]<<"pqiwork = "<<pqiwork<<"first"<<firstnormal<<"firstprime"<<firstnormalprime<<"second"<<secondnormal<<"third"<<thirdnormal <<std::endl;
        tmpVal += pqiwork;
*/
      }   //  ( pq[iWork] > 0 )

    // 0.5*2 = 1



      if ( lA[iWork]>=2 ) {

        lAm1[iWork] -=1; // now lAm1[iWork] == lA[iWork]-2
        tmpVal += 0.5 * ( lAm1[iWork]+1 ) * pripair1.one_over_gamma  
                  *(twoegaugevRRa000( pripair1, pripair2, p, q, FmT_2epri,
                    shell1, shell3, m, LA-2,lAm1 )
                - 1.0/(pripair2.one_over_gamma*Zeta) * twoegaugevRRa000( pripair1, pripair2,
                  p, q, FmT_2epri, shell1, shell3, m+1, LA-2,lAm1 ) );
      } // if lA[iWork]>=2

    } // if( LA>=1 ) 

    return tmpVal;

  };  // twoegaugevRRa000




  /**
   *  \brief (00|00) integral of gauge ERI when all the angular momentum are 0
   *
   *
   *  \param [in] pair1   bra shell pair data for shell1,shell2
   *  \param [in] pair2   ket shell pair data for shell3,shell4
   *  \param [in] p    
   *  \param [in] q    
   *  \param [in] shell1    
   *  \param [in] shell2    
   *  \param [in] shell3    
   *  \param [in] shell4    
   *
   *  \return ERI of two primitive shell pairs ( shell1 shell2 | shell3 shell4 )
   */
   
  double RealGTOIntEngine::twoegaugeSSSS0(
    libint2::ShellPair &pair1, libint2::ShellPair &pair2,
    int p , int q ,
    libint2::Shell &shell1, libint2::Shell &shell2,
    libint2::Shell &shell3, libint2::Shell &shell4 ) {

    // constants
    // gamma(1/2) = sqrt(pi)
    double gammahalf = sqrt(M_PI);
    // gamma(3/2) = 1/2*sqrt(M_PI);
    double gamma3half = 0.5*gammahalf;
    // gamma(1/2)/gamma(3/2) = 2 
    double gammahalfover3half = 2; 

    double tmpval;

    // in this auxiliary integral, m=0 
    double norm,sqrPQ,PQ[3],expoT,T,FmT[3],Kab,Kcd,SSSS0=0.0 ;
    for ( auto &pripair1 : pair1.primpairs )
    for ( auto &pripair2 : pair2.primpairs ) {
       
      tmpval = 0.0;
      sqrPQ = 0.0;
      for(int m=0 ; m<3 ; m++ ) {
        PQ[m] = ( pripair1.P[m]-pripair2.P[m] );
        sqrPQ += PQ[m]*PQ[m];
      } 

//      double zeta = 1.0/pripair1.one_over_gamma;
//      double eta  = 1.0/pripair2.one_over_gamma;

      // zeta+eta is expoT
      expoT = 1.0/pripair1.one_over_gamma + 1.0/pripair2.one_over_gamma;  

      // T= \rho * PQ^2
      T = sqrPQ/( pripair1.one_over_gamma + pripair2.one_over_gamma );

      // calculate F0(T)
      computeFmTTaylor( FmT, T, 2, 0 );

      Kab = sqrt(2.0)* pow(M_PI,1.25) * pripair1.K; 
      Kcd = sqrt(2.0)* pow(M_PI,1.25) * pripair2.K;

      norm = shell1.contr[0].coeff[pripair1.p1]* 
             shell2.contr[0].coeff[pripair1.p2]* 
             shell3.contr[0].coeff[pripair2.p1]* 
             shell4.contr[0].coeff[pripair2.p2];  

      // here we don't consider norm 
      if (p==q) {
        tmpval += 0.5* // (pripair1.one_over_gamma + pripair2.one_over_gamma) * 
           Kab * Kcd * 
          FmT[1] / sqrt(expoT)    *2 ;  // rho cancel with 1/zeta+1/eta 
// std::cout<<"tmpval first part = "<<0.5*Kab * Kcd *FmT[1] / sqrt(expoT)    *2<<std::endl;
      }
// \rho is 1/(1/zeta+1/eta)
//
//std::cout<<"tmpval second part = "<< PQ[p]*PQ[q]*2/( pripair1.one_over_gamma + pripair2.one_over_gamma ) *Kab * Kcd * (  FmT[1]    -  FmT[2]) / sqrt(expoT)<<std::endl;

      tmpval += PQ[p]*PQ[q]*2/( pripair1.one_over_gamma + pripair2.one_over_gamma ) 
                 *Kab * Kcd * (  FmT[1]    -  FmT[2]) / sqrt(expoT);

//std::cout<<" PQ[p] = "<< PQ[p] <<std::endl;
//std::cout<<" T = "<< T <<std::endl;

// std::cout<<"tmpval = "<<tmpval<<" FmT1= "<<FmT[1]<<" FmT2= "<<FmT[2]<<std::endl;
//std::cout<<"tmpval = "<<tmpval<<"FmT0 = "<<FmT[0]<<" FmT1= "<<FmT[1]<<" FmT2= "<<FmT[2]<<std::endl;

      SSSS0 += norm * tmpval; 

    } // for pripair2

    return SSSS0;
    
  } // twoegaugeSSSS0

  // AC derivative


  std::vector<std::vector<std::vector<double>>> RealGTOIntEngine::ACgaugederiv( 
    libint2::ShellPair &pair1, libint2::ShellPair &pair2, 
    libint2::Shell &shell1, libint2::Shell &shell2,
    libint2::Shell &shell3, libint2::Shell &shell4
    ) {


    int cart_i_size,cart_j_size,cart_k_size,cart_l_size; 
    int cart_ip[3], cart_im[3], cart_kp[3], cart_km[3];
    int cart_ip_size, cart_kp_size;
    int cart_im_size, cart_km_size;
    int LA = shell1.contr[0].l;
    int LB = shell2.contr[0].l;
    int LC = shell3.contr[0].l;
    int LD = shell4.contr[0].l;
    

    std::vector<std::vector<std::vector<double>>> ERIgrad_cart;
    cart_i_size = cart_ang_list[shell1.contr[0].l].size();
    cart_j_size = cart_ang_list[shell2.contr[0].l].size();
    cart_k_size = cart_ang_list[shell3.contr[0].l].size();
    cart_l_size = cart_ang_list[shell4.contr[0].l].size(); 

    int totalcartsize = cart_i_size*cart_j_size*cart_k_size*cart_l_size;

    // assign memory for outcome
    std::vector<std::vector<double>> tempplace(6);
    for (int component = 0 ; component < 6 ; component++ ) {
      tempplace[component].assign(cart_i_size*cart_j_size*cart_k_size*cart_l_size,0.0);
    }
    for ( int component = 0 ; component < 9 ; component++ ) {
      ERIgrad_cart.push_back(tempplace);
    }
    std::vector<std::vector<std::vector<double>>> buffgrad;



    // define shells with lA and lC plus minus 1
    //libint2::Shell test;
    //std::vector<libint2::Shell> shellrandom;

    //std::vector<libint2::Shell> shellsApm(2);
//new start
    libint2::Shell shellsAp;
    libint2::Shell shellsAm;
    shellsAp = shell1;
    shellsAp.contr[0].l +=1; 
    shellsAp.contr[0].pure = 0; // cartesian

    cart_ip_size = (shellsAp.contr[0].l+1)*(shellsAp.contr[0].l+2)/2;  //NBasis for lA+1

    if  ( shell1.contr[0].l>0 ) {
      shellsAm = shell1;
      shellsAm.contr[0].l -=1;  
      shellsAm.contr[0].pure = 0; // cartesian
      cart_im_size = (shellsAm.contr[0].l+1)*(shellsAm.contr[0].l+2)/2; //NBasis for la-1
    }
   // new end 
    
// old
/*
    shellsApm[0] = shell1;
    shellsApm[0].contr[0].l +=1; 
    shellsApm[0].contr[0].pure = 0; // cartesian

    cart_ip_size = (shellsApm[0].contr[0].l+1)*(shellsApm[0].contr[0].l+2)/2;  //NBasis for lA+1

    if  ( shell1.contr[0].l>0 ) {
      shellsApm[1] = shell1;
      shellsApm[1].contr[0].l -=1;  
      shellsApm[1].contr[0].pure = 0; // cartesian
      cart_im_size = (shellsApm[1].contr[0].l+1)*(shellsApm[1].contr[0].l+2)/2; //NBasis for la-1
    }
*/
// old
/*
    std::vector<libint2::Shell> shellsCpm(2);
    shellsCpm[0] = shell3;
    shellsCpm[0].contr[0].l +=1; 
    shellsCpm[0].contr[0].pure = 0; // cartesian
    cart_kp_size = (shellsCpm[0].contr[0].l+1)*(shellsCpm[0].contr[0].l+2)/2;

    if  ( shell3.contr[0].l>0 ) {
      shellsCpm[1] = shell3;
      shellsCpm[1].contr[0].l -=1;  
      shellsCpm[1].contr[0].pure = 0; // cartesian
      cart_km_size = (shellsCpm[1].contr[0].l+1)*(shellsCpm[1].contr[0].l+2)/2;
    }
*/

    libint2::Shell shellsCp;
    libint2::Shell shellsCm;
    shellsCp = shell3;
    shellsCp.contr[0].l +=1; 
    shellsCp.contr[0].pure = 0; // cartesian
    cart_kp_size = (shellsCp.contr[0].l+1)*(shellsCp.contr[0].l+2)/2;

    if  ( shell3.contr[0].l>0 ) {
      shellsCm = shell3;
      shellsCm.contr[0].l -=1;  
      shellsCm.contr[0].pure = 0; // cartesian
      cart_km_size = (shellsCm.contr[0].l+1)*(shellsCm.contr[0].l+2)/2;
    }


    
    // set pure to cartesian
    libint2::Shell shellB = shell2;
    shellB.contr[0].pure = 0; // cartesian

    libint2::Shell shellD = shell4;
    shellD.contr[0].pure = 0; // cartesian

    // compute integrals with different angular momentum

    // lA+1 lC+1
    //std::vector<double> buffpp  = RealGTOIntEngine::BottomupHGP(pair1,pair2,
    std::vector<std::vector<double>> buffpp  = RealGTOIntEngine::computegaugeabcd(pair1,pair2,
//      shellsApm[0],
      shellsAp,
      shellB,
//      shellsCpm[0],
      shellsCp,
      shellD,
      1, 0, 1, 0
    );

  //for (int ii = 0 ; ii < cart_ip_size*cart_j_size*cart_kp_size*cart_l_size ; ii++ ) {
  //  std::cout<<"ii = "<<ii<<" pp = "<<buffpp[ii]<<std::endl;
  //}

    buffgrad.push_back(buffpp);

    if  ( shell1.contr[0].l>0 ) {
        auto buffmp  = RealGTOIntEngine::computegaugeabcd(pair1,pair2,
//          shellsApm[1],
          shellsAm,
          shellB,
          //shellsCpm[0],
          shellsCp,
          shellD,
          0, 0, 1, 0
        );
  //for (int ii = 0 ; ii < cart_im_size*cart_j_size*cart_kp_size*cart_l_size ; ii++ ) {
  //  std::cout<<"ii = "<<ii<<" mp = "<<buffmp[ii]<<std::endl;
  //}
        buffgrad.push_back(buffmp);

    }  else {

        std::vector<std::vector<double>> place;
        buffgrad.push_back(place);
    }

    if  ( shell3.contr[0].l>0 ) {

        //auto buffpm  = RealGTOIntEngine::BottomupHGP(pair1,pair2,
        auto buffpm  = RealGTOIntEngine::computegaugeabcd(pair1,pair2,
//          shellsApm[0],
          shellsAp,
          shellB,
          //shellsCpm[1],
          shellsCm,
          shellD,
          1, 0, 0, 0
        );
/*
  for (int ii = 0 ; ii < cart_ip_size*cart_j_size*cart_km_size*cart_l_size ; ii++ ) {
    std::cout<<"ii = "<<ii<<" pm = "<<buffpm[ii]<<std::endl;
  }
*/
        buffgrad.push_back(buffpm);

    }  else {
        std::vector<std::vector<double>> place;
        buffgrad.push_back(place);
    }
    

    if ( shell1.contr[0].l>0  and  shell3.contr[0].l>0 )  { 
        //auto buffmm  = RealGTOIntEngine::BottomupHGP(pair1,pair2,
        auto buffmm  = RealGTOIntEngine::computegaugeabcd(pair1,pair2,
//          shellsApm[1],
          shellsAm,
          shellB,
          //shellsCpm[1],
          shellsCm,
          shellD,
          0, 0, 0, 0
        );
  //for (int ii = 0 ; ii < cart_im_size*cart_j_size*cart_km_size*cart_l_size ; ii++ ) {
  //  std::cout<<"ii = "<<ii<<" mm = "<<buffmm[ii]<<std::endl;
  //}
        buffgrad.push_back(buffmm);


    } else {
        std::vector<std::vector<double>> place;
        buffgrad.push_back(place);
    }

    // combine into derivatives
    //

    for (int cart_i = 0; cart_i < cart_i_size ; cart_i++) { 
      int lA_xyz[3];     

      for (int ii = 0; ii < 3; ii++) { 

        lA_xyz[ii] = cart_ang_list[shell1.contr[0].l][cart_i][ii];

      } // for (int ii = 0; ii < 3; ii++)
      // here find out the cart_ip ( cart i of lA_x + 1 ) 
//old
//      cart_ip[0] = indexmap(shellsApm[0].contr[0].l, lA_xyz[0]+1, lA_xyz[1], lA_xyz[2] );
//     cart_ip[1] = indexmap(shellsApm[0].contr[0].l, lA_xyz[0], lA_xyz[1]+1, lA_xyz[2] );
//      cart_ip[2] = indexmap(shellsApm[0].contr[0].l, lA_xyz[0], lA_xyz[1], lA_xyz[2]+1 );

      cart_ip[0] = indexmap(shellsAp.contr[0].l, lA_xyz[0]+1, lA_xyz[1], lA_xyz[2] );
      cart_ip[1] = indexmap(shellsAp.contr[0].l, lA_xyz[0], lA_xyz[1]+1, lA_xyz[2] );
      cart_ip[2] = indexmap(shellsAp.contr[0].l, lA_xyz[0], lA_xyz[1], lA_xyz[2]+1 );

      if  ( lA_xyz[0]>0 ) {  
        //cart_im[0] = indexmap(shellsApm[1].contr[0].l, lA_xyz[0]-1, lA_xyz[1], lA_xyz[2] ); 
        cart_im[0] = indexmap(shellsAm.contr[0].l, lA_xyz[0]-1, lA_xyz[1], lA_xyz[2] ); 
      }

      if  ( lA_xyz[1]>0 ) {  
        //cart_im[1] = indexmap(shellsApm[1].contr[0].l, lA_xyz[0], lA_xyz[1]-1, lA_xyz[2] ); 
        cart_im[1] = indexmap(shellsAm.contr[0].l, lA_xyz[0], lA_xyz[1]-1, lA_xyz[2] ); 
      }

      if  ( lA_xyz[2]>0 ) {  
        //cart_im[2] = indexmap(shellsApm[1].contr[0].l, lA_xyz[0], lA_xyz[1], lA_xyz[2]-1 ); 
        cart_im[2] = indexmap(shellsAm.contr[0].l, lA_xyz[0], lA_xyz[1], lA_xyz[2]-1 ); 
      }  
       
      for (int cart_j = 0; cart_j < cart_j_size; cart_j++) { 
        int lB_xyz[3];     

        for (int ii = 0; ii < 3; ii++) { 

          lB_xyz[ii] = cart_ang_list[shell2.contr[0].l][cart_j][ii];

        } // for (int ii = 0; ii < 3; ii++)


        for (int cart_k = 0; cart_k < cart_k_size; cart_k++) { 
          int lC_xyz[3];     

          for (int ii = 0; ii < 3; ii++) { 

            lC_xyz[ii] = cart_ang_list[shell3.contr[0].l][cart_k][ii];
          } // for (int ii = 0; ii < 3; ii++)
          // here find out the cart_kp ( cart i of lC_x + 1 ) 
//          cart_kp[0] = indexmap(shellsCpm[0].contr[0].l, lC_xyz[0]+1, lC_xyz[1], lC_xyz[2] );
//          cart_kp[1] = indexmap(shellsCpm[0].contr[0].l, lC_xyz[0], lC_xyz[1]+1, lC_xyz[2] );
//          cart_kp[2] = indexmap(shellsCpm[0].contr[0].l, lC_xyz[0], lC_xyz[1], lC_xyz[2]+1 );
          cart_kp[0] = indexmap(shellsCp.contr[0].l, lC_xyz[0]+1, lC_xyz[1], lC_xyz[2] );
          cart_kp[1] = indexmap(shellsCp.contr[0].l, lC_xyz[0], lC_xyz[1]+1, lC_xyz[2] );
          cart_kp[2] = indexmap(shellsCp.contr[0].l, lC_xyz[0], lC_xyz[1], lC_xyz[2]+1 );
          
          if  ( lC_xyz[0]>0 ) {  
            //cart_km[0] = indexmap(shellsCpm[1].contr[0].l, lC_xyz[0]-1, lC_xyz[1], lC_xyz[2] ); 
            cart_km[0] = indexmap(shellsCm.contr[0].l, lC_xyz[0]-1, lC_xyz[1], lC_xyz[2] ); 
          }
          if  ( lC_xyz[1]>0 ) {  
            //cart_km[1] = indexmap(shellsCpm[1].contr[0].l, lC_xyz[0], lC_xyz[1]-1, lC_xyz[2] ); 
            cart_km[1] = indexmap(shellsCm.contr[0].l, lC_xyz[0], lC_xyz[1]-1, lC_xyz[2] ); 
          }
          if  ( lC_xyz[2]>0 ) {  
            //cart_km[2] = indexmap(shellsCpm[1].contr[0].l, lC_xyz[0], lC_xyz[1], lC_xyz[2]-1 ); 
            cart_km[2] = indexmap(shellsCm.contr[0].l, lC_xyz[0], lC_xyz[1], lC_xyz[2]-1 ); 
          }  

          for (int cart_l = 0; cart_l < cart_l_size; cart_l++) { 
            int lD_xyz[3];     

            for (int ii = 0; ii < 3; ii++) { 

              lD_xyz[ii] = cart_ang_list[shell4.contr[0].l][cart_l][ii];
            } // for ii
            // (2)
            for ( int aidx = 0 ; aidx < 3 ; aidx++ ) {
            for ( int cidx = 0 ; cidx < 3 ; cidx++ ) {
              int acidx = aidx*3+cidx;
              for ( int mu = 0 ; mu < 6 ; mu++ ) {
                // place the lA+1 and lC+1 into integral
                ERIgrad_cart[acidx][mu][ cart_i*cart_j_size*cart_k_size*cart_l_size+
                  cart_j * cart_k_size*cart_l_size + cart_k* cart_l_size
                  + cart_l]
                += 4*buffgrad[0][mu][ cart_ip[aidx] *cart_j_size*cart_kp_size*cart_l_size+
                  cart_j * cart_kp_size*cart_l_size +cart_kp[cidx]* cart_l_size
                  + cart_l];
 
                  // lA-1, lC+1
                  // if  ( basisSet.shells[s1].contr[0].l>0 ) {
                  if  ( lA_xyz[aidx]>0 ) {
                    ERIgrad_cart[acidx][mu][ cart_i*cart_j_size*cart_k_size*cart_l_size+
                      cart_j * cart_k_size*cart_l_size + cart_k* cart_l_size
                      + cart_l]
                      -= 2*lA_xyz[aidx]*buffgrad[1][mu][ cart_im[aidx] *cart_j_size*cart_kp_size*cart_l_size+
                      cart_j * cart_kp_size*cart_l_size +cart_kp[cidx]* cart_l_size
                      + cart_l];
                  }
                 
                  // lA+1, lC-1
                  if  ( lC_xyz[cidx]>0 ) {
                    ERIgrad_cart[acidx][mu][ cart_i*cart_j_size*cart_k_size*cart_l_size+
                      cart_j * cart_k_size*cart_l_size + cart_k* cart_l_size
                      + cart_l]
                      -= 2*lC_xyz[cidx]*buffgrad[2][mu][ cart_ip[aidx] *cart_j_size*cart_km_size*cart_l_size+
                      cart_j * cart_km_size*cart_l_size +cart_km[cidx]* cart_l_size
                      + cart_l];
                  }
                 
                  // lA-1, lC-1
                  if  (lA_xyz[aidx]>0  and  lC_xyz[cidx]>0 ) {
                    ERIgrad_cart[acidx][mu][ cart_i*cart_j_size*cart_k_size*cart_l_size+
                      cart_j * cart_k_size*cart_l_size + cart_k* cart_l_size
                      + cart_l]
                      += lA_xyz[aidx]*lC_xyz[cidx]*buffgrad[3][mu][ cart_im[aidx] *cart_j_size*cart_km_size*cart_l_size+
                      cart_j * cart_km_size*cart_l_size +cart_km[cidx]* cart_l_size
                      + cart_l];
                  }
              } // mu from 0 to 5 
            } // cidx
            } // aidx 

          }  // for cart_l

        } // for cart_k
      } // for cart_j
    } // for cart_i 

     
    if ( ( not shell1.contr[0].pure ) and ( not shell2.contr[0].pure ) and 
         ( not shell3.contr[0].pure ) and ( not shell4.contr[0].pure ) ) {  
      // if both sides are cartesian, return cartesian gaussian integrals
      return ERIgrad_cart;
    }

    // assign memory for outcome
    std::vector<std::vector<std::vector<double>>> ERIgrad_sph;
    std::vector<std::vector<double>> sphtempplace(6);
    for (int component = 0 ; component < 6 ; component++ ) {
      sphtempplace[component].assign( (2*LA+1)*(2*LB+1)*(2*LC+1)*(2*LD+1) , 0.0 );
    }
    for ( int component = 0 ; component < 9 ; component++ ) {
      ERIgrad_sph.push_back(sphtempplace);
    }

    int k = 0;

    for ( int derivcomp = 0 ; derivcomp < 9 ; derivcomp++ ) {
      for ( int gaugecomp = 0 ; gaugecomp < 6 ; gaugecomp++ ) {
        cart2sph_2e_transform(LA,LB,LC,LD,ERIgrad_sph[derivcomp][gaugecomp],  
          ERIgrad_cart[derivcomp][gaugecomp] );
      }
    }   
  

    return ERIgrad_sph; 
    // first layer : xyz of cartesian, second layer xx, xy, xz component of gauge integral, third layer: pqrs index 



  } // ACgaugederiv





  std::vector<std::vector<std::vector<double>>> RealGTOIntEngine::BCgaugederiv( 
    libint2::ShellPair &pair1, libint2::ShellPair &pair2, 
    libint2::Shell &shell1, libint2::Shell &shell2,
    libint2::Shell &shell3, libint2::Shell &shell4
    ) {


    int cart_i_size,cart_j_size,cart_k_size,cart_l_size; 
    int cart_jp[3], cart_jm[3], cart_kp[3], cart_km[3];
    int cart_jp_size, cart_kp_size;
    int cart_jm_size, cart_km_size;
    int LA = shell1.contr[0].l;
    int LB = shell2.contr[0].l;
    int LC = shell3.contr[0].l;
    int LD = shell4.contr[0].l;
    

    std::vector<std::vector<std::vector<double>>> ERIgrad_cart;
    cart_i_size = cart_ang_list[shell1.contr[0].l].size();
    cart_j_size = cart_ang_list[shell2.contr[0].l].size();
    cart_k_size = cart_ang_list[shell3.contr[0].l].size();
    cart_l_size = cart_ang_list[shell4.contr[0].l].size(); 

    int totalcartsize = cart_i_size*cart_j_size*cart_k_size*cart_l_size;

    // assign memory for outcome
    std::vector<std::vector<double>> tempplace(6);
    for (int component = 0 ; component < 6 ; component++ ) {
      tempplace[component].assign(cart_i_size*cart_j_size*cart_k_size*cart_l_size,0.0);
    }
    for ( int component = 0 ; component < 9 ; component++ ) {
      ERIgrad_cart.push_back(tempplace);
    }
    std::vector<std::vector<std::vector<double>>> buffgrad;



    // define shells with lB and lC plus minus 1

    std::vector<libint2::Shell> shellsBpm(2);
    shellsBpm[0] = shell1;
    shellsBpm[0].contr[0].l +=1; 
    shellsBpm[0].contr[0].pure = 0; // cartesian

    cart_jp_size = (shellsBpm[0].contr[0].l+1)*(shellsBpm[0].contr[0].l+2)/2;  //NBasis for lA+1

    if  ( shell2.contr[0].l>0 ) {
      shellsBpm[1] = shell2;
      shellsBpm[1].contr[0].l -=1;  
      shellsBpm[1].contr[0].pure = 0; // cartesian
      cart_jm_size = (shellsBpm[1].contr[0].l+1)*(shellsBpm[1].contr[0].l+2)/2; //NBasis for lb-1
    }


    std::vector<libint2::Shell> shellsCpm(2);
    shellsCpm[0] = shell3;
    shellsCpm[0].contr[0].l +=1; 
    shellsCpm[0].contr[0].pure = 0; // cartesian
    cart_kp_size = (shellsCpm[0].contr[0].l+1)*(shellsCpm[0].contr[0].l+2)/2;

    if  ( shell3.contr[0].l>0 ) {
      shellsCpm[1] = shell3;
      shellsCpm[1].contr[0].l -=1;  
      shellsCpm[1].contr[0].pure = 0; // cartesian
      cart_km_size = (shellsCpm[1].contr[0].l+1)*(shellsCpm[1].contr[0].l+2)/2;
    }
    
    // set pure to cartesian
    libint2::Shell shellA = shell1;
    shellA.contr[0].pure = 0; // cartesian

    libint2::Shell shellD = shell4;
    shellD.contr[0].pure = 0; // cartesian

    // compute integrals with different angular momentum

    // lB+1 lC+1
    //std::vector<double> buffpp  = RealGTOIntEngine::BottomupHGP(pair1,pair2,
    std::vector<std::vector<double>> buffpp  = RealGTOIntEngine::computegaugeabcd(pair1,pair2,
      shellA,
      shellsBpm[0],
      shellsCpm[0],
      shellD,
      0, 1, 1, 0
    );

  //for (int ii = 0 ; ii < cart_ip_size*cart_j_size*cart_kp_size*cart_l_size ; ii++ ) {
  //  std::cout<<"ii = "<<ii<<" pp = "<<buffpp[ii]<<std::endl;
  //}

    buffgrad.push_back(buffpp);

    if  ( shell2.contr[0].l>0 ) {
        auto buffmp  = RealGTOIntEngine::computegaugeabcd(pair1,pair2,
          shellA,
          shellsBpm[1],
          shellsCpm[0],
          shellD,
          0, 0, 1, 0
        );
  //for (int ii = 0 ; ii < cart_im_size*cart_j_size*cart_kp_size*cart_l_size ; ii++ ) {
  //  std::cout<<"ii = "<<ii<<" mp = "<<buffmp[ii]<<std::endl;
  //}
        buffgrad.push_back(buffmp);

    }  else {

        std::vector<std::vector<double>> place;
        buffgrad.push_back(place);
    }

    if  ( shell3.contr[0].l>0 ) {

        //auto buffpm  = RealGTOIntEngine::BottomupHGP(pair1,pair2,
        auto buffpm  = RealGTOIntEngine::computegaugeabcd(pair1,pair2,
          shellA,
          shellsBpm[0],
          shellsCpm[1],
          shellD,
          0, 1, 0, 0
        );
/*
  for (int ii = 0 ; ii < cart_ip_size*cart_j_size*cart_km_size*cart_l_size ; ii++ ) {
    std::cout<<"ii = "<<ii<<" pm = "<<buffpm[ii]<<std::endl;
  }
*/
        buffgrad.push_back(buffpm);

    }  else {
        std::vector<std::vector<double>> place;
        buffgrad.push_back(place);
    }
    

    if ( shell2.contr[0].l>0  and  shell3.contr[0].l>0 )  { 
        //auto buffmm  = RealGTOIntEngine::BottomupHGP(pair1,pair2,
        auto buffmm  = RealGTOIntEngine::computegaugeabcd(pair1,pair2,
          shellA,
          shellsBpm[1],
          shellsCpm[1],
          shellD,
          0, 0, 0, 0
        );
  //for (int ii = 0 ; ii < cart_im_size*cart_j_size*cart_km_size*cart_l_size ; ii++ ) {
  //  std::cout<<"ii = "<<ii<<" mm = "<<buffmm[ii]<<std::endl;
  //}
        buffgrad.push_back(buffmm);


    } else {
        std::vector<std::vector<double>> place;
        buffgrad.push_back(place);
    }

    // combine into derivatives
    //

    for (int cart_i = 0; cart_i < cart_i_size ; cart_i++) { 
      int lA_xyz[3];     

      for (int ii = 0; ii < 3; ii++) { 

        lA_xyz[ii] = cart_ang_list[shell1.contr[0].l][cart_i][ii];

      } // for (int ii = 0; ii < 3; ii++)


      for (int cart_j = 0; cart_j < cart_j_size; cart_j++) { 
        int lB_xyz[3];     

        for (int ii = 0; ii < 3; ii++) { 

          lB_xyz[ii] = cart_ang_list[shell2.contr[0].l][cart_j][ii];

        } // for (int ii = 0; ii < 3; ii++)
        // here find out the cart_jp ( cart j of lB_x + 1 ) 
        cart_jp[0] = indexmap(shellsBpm[0].contr[0].l, lB_xyz[0]+1, lB_xyz[1], lB_xyz[2] );
        cart_jp[1] = indexmap(shellsBpm[0].contr[0].l, lB_xyz[0], lB_xyz[1]+1, lB_xyz[2] );
        cart_jp[2] = indexmap(shellsBpm[0].contr[0].l, lB_xyz[0], lB_xyz[1], lB_xyz[2]+1 );


        if  ( lB_xyz[0]>0 ) {  
          cart_jm[0] = indexmap(shellsBpm[1].contr[0].l, lB_xyz[0]-1, lB_xyz[1], lB_xyz[2] ); 
        }

        if  ( lB_xyz[1]>0 ) {  
          cart_jm[1] = indexmap(shellsBpm[1].contr[0].l, lB_xyz[0], lB_xyz[1]-1, lB_xyz[2] ); 
        }
 
        if  ( lB_xyz[2]>0 ) {  
          cart_jm[2] = indexmap(shellsBpm[1].contr[0].l, lB_xyz[0], lB_xyz[1], lB_xyz[2]-1 ); 
        }  
       


        for (int cart_k = 0; cart_k < cart_k_size; cart_k++) { 
          int lC_xyz[3];     

          for (int ii = 0; ii < 3; ii++) { 

            lC_xyz[ii] = cart_ang_list[shell3.contr[0].l][cart_k][ii];
          } // for (int ii = 0; ii < 3; ii++)
          // here find out the cart_kp ( cart i of lC_x + 1 ) 
          cart_kp[0] = indexmap(shellsCpm[0].contr[0].l, lC_xyz[0]+1, lC_xyz[1], lC_xyz[2] );
          cart_kp[1] = indexmap(shellsCpm[0].contr[0].l, lC_xyz[0], lC_xyz[1]+1, lC_xyz[2] );
          cart_kp[2] = indexmap(shellsCpm[0].contr[0].l, lC_xyz[0], lC_xyz[1], lC_xyz[2]+1 );
          
          if  ( lC_xyz[0]>0 ) {  
            cart_km[0] = indexmap(shellsCpm[1].contr[0].l, lC_xyz[0]-1, lC_xyz[1], lC_xyz[2] ); 
          }
          if  ( lC_xyz[1]>0 ) {  
            cart_km[1] = indexmap(shellsCpm[1].contr[0].l, lC_xyz[0], lC_xyz[1]-1, lC_xyz[2] ); 
          }
          if  ( lC_xyz[2]>0 ) {  
            cart_km[2] = indexmap(shellsCpm[1].contr[0].l, lC_xyz[0], lC_xyz[1], lC_xyz[2]-1 ); 
          }  

          for (int cart_l = 0; cart_l < cart_l_size; cart_l++) { 
            int lD_xyz[3];     

            for (int ii = 0; ii < 3; ii++) { 

              lD_xyz[ii] = cart_ang_list[shell4.contr[0].l][cart_l][ii];
            } // for ii
            // (2)
            for ( int bidx = 0 ; bidx < 3 ; bidx++ ) {
            for ( int cidx = 0 ; cidx < 3 ; cidx++ ) {
              int bcidx = bidx*3+cidx;
              for ( int mu = 0 ; mu < 6 ; mu++ ) {
                // place the lB+1 and lC+1 into integral
                ERIgrad_cart[bcidx][mu][ cart_i*cart_j_size*cart_k_size*cart_l_size+
                  cart_j * cart_k_size*cart_l_size + cart_k* cart_l_size
                  + cart_l]
                += 4*buffgrad[0][mu][ cart_i *cart_jp_size*cart_kp_size*cart_l_size+
                  cart_jp[bidx] * cart_kp_size*cart_l_size +cart_kp[cidx]* cart_l_size
                  + cart_l];
 
                  // lB-1, lC+1
                  // if  ( basisSet.shells[s2].contr[0].l>0 ) {
                  if  ( lB_xyz[bidx]>0 ) {
                    ERIgrad_cart[bcidx][mu][ cart_i*cart_j_size*cart_k_size*cart_l_size+
                      cart_j * cart_k_size*cart_l_size + cart_k* cart_l_size
                      + cart_l]
                      -= 2*lB_xyz[bidx]*buffgrad[1][mu][ cart_i *cart_jm_size*cart_kp_size*cart_l_size+
                      cart_jm[bidx] * cart_kp_size*cart_l_size +cart_kp[cidx]* cart_l_size
                      + cart_l];
                  }
                 
                  // lB+1, lC-1
                  if  ( lC_xyz[cidx]>0 ) {
                    ERIgrad_cart[bcidx][mu][ cart_i*cart_j_size*cart_k_size*cart_l_size+
                      cart_j * cart_k_size*cart_l_size + cart_k* cart_l_size
                      + cart_l]
                      -= 2*lC_xyz[cidx]*buffgrad[2][mu][ cart_i *cart_jp_size*cart_km_size*cart_l_size+
                      cart_jp[bidx] * cart_km_size*cart_l_size +cart_km[cidx]* cart_l_size
                      + cart_l];
                  }
                 
                  // lB-1, lC-1
                  if  (lB_xyz[bidx]>0  and  lC_xyz[cidx]>0 ) {
                    ERIgrad_cart[bcidx][mu][ cart_i*cart_j_size*cart_k_size*cart_l_size+
                      cart_j * cart_k_size*cart_l_size + cart_k* cart_l_size
                      + cart_l]
                      += lB_xyz[bidx]*lC_xyz[cidx]*buffgrad[3][mu][ cart_i *cart_jm_size*cart_km_size*cart_l_size+
                      cart_jm[bidx] * cart_km_size*cart_l_size +cart_km[cidx]* cart_l_size
                      + cart_l];
                  }
              } // mu from 0 to 5 
            } // cidx
            } // bidx 

          }  // for cart_l

        } // for cart_k
      } // for cart_j
    } // for cart_i 

     
    if ( ( not shell1.contr[0].pure ) and ( not shell2.contr[0].pure ) and 
         ( not shell3.contr[0].pure ) and ( not shell4.contr[0].pure ) ) {  
      // if both sides are cartesian, return cartesian gaussian integrals
      return ERIgrad_cart;
    }

    // assign memory for outcome
    std::vector<std::vector<std::vector<double>>> ERIgrad_sph;
    std::vector<std::vector<double>> sphtempplace(6);
    for (int component = 0 ; component < 6 ; component++ ) {
      sphtempplace[component].assign( (2*LA+1)*(2*LB+1)*(2*LC+1)*(2*LD+1) , 0.0 );
    }
    for ( int component = 0 ; component < 9 ; component++ ) {
      ERIgrad_sph.push_back(sphtempplace);
    }

    int k = 0;

    for ( int derivcomp = 0 ; derivcomp < 9 ; derivcomp++ ) {
      for ( int gaugecomp = 0 ; gaugecomp < 6 ; gaugecomp++ ) {
        cart2sph_2e_transform(LA,LB,LC,LD,ERIgrad_sph[derivcomp][gaugecomp],  
          ERIgrad_cart[derivcomp][gaugecomp] );
      }
    }   
  

    return ERIgrad_sph; 
    // first layer : xyz of cartesian (cartesian components of derivatives), second layer xx, xy, xz component of gauge integral, third layer: pqrs index 



  } // BCgaugederiv



}  // namespace ChronusQ 

