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
#include <molecule.hpp>

namespace ChronusQ {

  static double smallT[21]={
  1.00000000000000,
  0.33333333333333,
  0.20000000000000,
  0.14285714285714,
  0.11111111111111,
  0.09090909090909,
  0.07692307692308,
  0.06666666666667,
  0.05882352941176,
  0.05263157894737,
  0.04761904761905,
  0.04347826086957,
  0.04000000000000,
  0.03703703703704,
  0.03448275862069,
  0.03225806451613,
  0.03030303030303,
  0.02857142857143,
  0.02702702702703,
  0.02564102564103,
  0.02439024390244
  };
  
  static double factTLarge[21] = {
  0.88622692545275800,
  0.44311346272637900,
  0.66467019408956900,
  1.66167548522392000,
  5.81586419828372000,
  26.1713888922768000,
  143.942638907522000,
  935.627152898894000,
  7017.20364674171000,
  59646.2309973045000,
  566639.194474393000,
  5949711.54198112000,
  68421682.7327829000,
  855271034.159787000,
  11546158961.1571000,
  167419304936.778000,
  2594999226520.06000,
  42817487237581.0000,
  749306026657668.000,
  13862161493166900.0,
  270312149116754000.0
  };
  
  std::array<std::array<double,25>,3201> FmTTable;

  void generateFmTTable(){
  
    double intervalFmT = 0.025;
    double T = 0.0;
    int MaxTotalL=25;
    int MaxFmTPt = 3201;
    double critT = 33.0;  // critical value for T. for T>critT, use limit formula
    double expT, factor, term, sum, twoT, Tn;
    for( int i = 0; i < MaxFmTPt ; i++){
      if( std::abs(T) <= 1.0e-10 ) {
        for( int m = 0; m <= MaxTotalL ; m++) 
          FmTTable[i][m] =  1.0/ (2.0 * m + 1);
      } else if ( T > critT ) {
        FmTTable[i][0] = 0.5 * sqrt( M_PI/ T);
        twoT = 2.0 * T;
        Tn = 1.0;
        for(int m = 1; m < MaxTotalL; m++) {
          Tn *= twoT;
          FmTTable[i][m] = FmTTable[i][m-1] * (2*m-1) / twoT ;
        }
      } else {
        expT   = exp(-T);
        factor = MaxTotalL + 0.5;
        term   = 0.5 / factor;
        sum    = term;
        while(term > 1.0e-10) {
          factor += 1.0;
          term   *= T / factor;
          sum    += term;
        };
        FmTTable[i][MaxTotalL] = expT * sum;
        twoT = 2.0 * T;
        for(int m = MaxTotalL - 1; m  >=  0; m--) 
          FmTTable[i][m] = (twoT * FmTTable[i][m+1] + expT)/(2*m + 1);
      } // else
      T += intervalFmT;
    } // for i
  } // generateFmTTable()
  
    
  void RealGTOIntEngine::computeFmTTaylor(
    double *FmT, double T, int maxM, int minM){
    int m,i;
    double intervalFmT = 0.025;
    double critT = 33.0;

//std::cout<<"T value "<<T<<std::endl;
    if(T > critT) {

    // when T> T critical,
    // F_m(T) = (2m-1)!!  (pi)^1/2
    //          --------- ---------
    //          2(2T)^n   (T )^1/2
    //        = factTLarge 
    //          -------------
    //          T^n * T ^1/2
    //
    //  where factTLarge = (2n-1)!! pi^1/2
    //                     --------
    //                      2^(n+1) 
    //
    // S.Obara, A. Saika, J. Chem. Phys. 84,3963
    //
      double Tn = sqrt(T);
      FmT[0] = factTLarge[0] / Tn;
      for(m = minM + 1; m <= maxM; m++) {
        Tn    *= T;
        if ( m<=20 ) {
          FmT[m] = factTLarge[m]/Tn;
        } else {
          // when m exceed the pretabulate limit (m>20) 
          auto tmpfactTLarge = doubleFact(m)*sqrt(M_PI)/pow(2.0,m+1);
          FmT[m] = tmpfactTLarge/Tn ; 
        }          
      }
    } else {
      int    j      = T / intervalFmT;
      double deltaT = j * intervalFmT - T;
      double sum = FmTTable[j][maxM];
      double tmp = 1.0;
      double change;
      for(i = 1; i < 7 ; i++){
      // here seems it use up to i = 4. may have a try with i = 6. 
        tmp   *= deltaT / i;
        change = tmp * FmTTable[j][maxM + i];
        sum   += change;
      }
  //  down-recursion to obtain m<maxM values
      FmT[maxM] = sum;
      if(minM < maxM) {
        double twoT = 2.0 * T;
        double expT = exp(-T);
        for(int m = maxM - 1; m  >=  minM; m--) FmT[m] = (twoT * FmT[m+1] + expT) * smallT[m];
      }
    }
  }
  /**
   *  \brief Computes a shell block of the overlap matrix.
   *
   *
   *  \param [in] pair    Shell pair data for shell1, shell2
   *  \param [in] shell1  Bra shell
   *  \param [in] shell2  Ket shell
   *
   *  \returns Shell block of the overlap matrix for (shell1 | shell2)
   */ 
   
  std::vector<std::vector<double>> RealGTOIntEngine::computeOverlapS(
    libint2::ShellPair &pair, libint2::Shell &shell1, libint2::Shell &shell2 ){
  
    int nPGTOPair = pair.primpairs.size();
    int nElement = cart_ang_list[shell1.contr[0].l].size() 
                   * cart_ang_list[shell2.contr[0].l].size();
                    
    // evaluate the number of integral in the shell pair with cartesian gaussian
    // XXX: Allocates memory
    std::vector<double> S_shellpair(nElement);
    
    int lA[3],lB[3];  
    double S;
  
    if( shell1  ==  shell2 ) {

      for(int i = 0; i < cart_ang_list[shell1.contr[0].l].size() ; i++)  
      for(int j = i; j < cart_ang_list[shell2.contr[0].l].size() ; j++) {

        for(int k = 0; k < 3; k++) {
          lA[k] = cart_ang_list[shell1.contr[0].l][i][k];
          lB[k] = cart_ang_list[shell2.contr[0].l][j][k];  
        };

  // std::cout<<"lA, lB    lA "<<lA[0]<<"\t"<<lA[1]<<"\t"<<lA[2]<<"\t lB"<<lB[0]<<"\t"<<lB[1]<<"\t"<<lB[2]<<std::endl;
  
         S = hRRSab( pair , shell1, shell2, shell1.contr[0].l , lA , shell2.contr[0].l ,lB);
  
         if( std::abs(S) < 1.0e-15 ) S = 0.0;

         S_shellpair[i*cart_ang_list[shell2.contr[0].l].size() + j] = S;
         S_shellpair[j*cart_ang_list[shell2.contr[0].l].size() + i] = S;

      } // loop ij

    } else {

      S_shellpair.resize(0);

      for(int i = 0; i < cart_ang_list[shell1.contr[0].l].size(); i++) 
      for(int j = 0; j < cart_ang_list[shell2.contr[0].l].size(); j++){
        for(int k = 0 ; k < 3 ; k++ ){
          lA[k] = cart_ang_list[shell1.contr[0].l][i][k];
          lB[k] = cart_ang_list[shell2.contr[0].l][j][k];
        };
  
  // std::cout<<"lA, lB    lA "<<lA[0]<<"\t"<<lA[1]<<"\t"<<lA[2]<<"\t lB"<<lB[0]<<"\t"<<lB[1]<<"\t"<<lB[2]<<std::endl;
  
        S = hRRSab( pair , shell1, shell2, shell1.contr[0].l , lA , shell2.contr[0].l ,lB);

        if( std::abs(S) < 1.0e-15 ) S = 0.0;
        S_shellpair.push_back(S);

      } // loop ij

    } // else when the two shells are different 
  
  // p1: index of primitive gaussian in shell1, p2: index of primitive gaussian in shell2.
  // a1: zeta_a, a2:zeta_b, one_over_gamma: 1/(zeta_a+zeta_b) 


    if ( ( not shell1.contr[0].pure ) and ( not shell2.contr[0].pure ) ) {  
      // if both sides are cartesian, return cartesian gaussian integrals
      return {S_shellpair};
    }

    std::vector<std::vector<double>> S_shellpair_sph(1);
    S_shellpair_sph[0].assign(((2*shell1.contr[0].l+1)*(2*shell2.contr[0].l+1)),0.0);
    
    // Convert from cartesian functions to spherical functions
    cart2sph_transform(shell1.contr[0].l,shell2.contr[0].l, S_shellpair_sph[0],S_shellpair);

    return S_shellpair_sph;
  }

  /**
   *  \brief Computes a shell block of the Kinetic matrix.
   *
   *
   *  \param [in] pair    Shell pair data for shell1, shell2
   *  \param [in] shell1  Bra shell
   *  \param [in] shell2  Ket shell
   *
   *  \returns Shell block of the kinetic integral matrix for (shell1 | shell2)
   */ 
   
  std::vector<std::vector<double>> RealGTOIntEngine::computeKineticT(
    libint2::ShellPair &pair, libint2::Shell &shell1, libint2::Shell &shell2 ){
  
//    int nPGTOPair = pair.primpairs.size();
    int nElement = cart_ang_list[shell1.contr[0].l].size() 
                   * cart_ang_list[shell2.contr[0].l].size();
                    
    // evaluate the number of integral in the shell pair with cartesian gaussian
    // XXX: Allocates memory
    std::vector<double> T_shellpair(nElement);
    
    int lA[3],lB[3];  
    double T;
  
/*
    if( shell1  ==  shell2 ) {

      for(int i = 0; i < cart_ang_list[shell1.contr[0].l].size() ; i++)  
      for(int j = i; j < cart_ang_list[shell2.contr[0].l].size() ; j++) {

        for(int k = 0; k < 3; k++) {
          lA[k] = cart_ang_list[shell1.contr[0].l][i][k];
          lB[k] = cart_ang_list[shell2.contr[0].l][j][k];  
        };

  // std::cout<<"lA, lB    lA "<<lA[0]<<"\t"<<lA[1]<<"\t"<<lA[2]<<"\t lB"<<lB[0]<<"\t"<<lB[1]<<"\t"<<lB[2]<<std::endl;
  
         S = hRRSab( pair , shell1, shell2, shell1.contr[0].l , lA , shell2.contr[0].l ,lB);
  
         if( std::abs(S) < 1.0e-15 ) S = 0.0;

         S_shellpair[i*cart_ang_list[shell2.contr[0].l].size() + j] = S;
         S_shellpair[j*cart_ang_list[shell2.contr[0].l].size() + i] = S;

      } // loop ij

    } else {
*/

      T_shellpair.resize(0);

      for(int i = 0; i < cart_ang_list[shell1.contr[0].l].size(); i++) 
      for(int j = 0; j < cart_ang_list[shell2.contr[0].l].size(); j++){
        for(int k = 0 ; k < 3 ; k++ ){
          lA[k] = cart_ang_list[shell1.contr[0].l][i][k];
          lB[k] = cart_ang_list[shell2.contr[0].l][j][k];
        };
  
  // std::cout<<"lA, lB    lA "<<lA[0]<<"\t"<<lA[1]<<"\t"<<lA[2]<<"\t lB"<<lB[0]<<"\t"<<lB[1]<<"\t"<<lB[2]<<std::endl;
  
      //  T = vRRTab( pair , shell1, shell2, shell1.contr[0].l , lA , 
        T = RRTab( pair , shell1, shell2, shell1.contr[0].l , lA , 
              shell2.contr[0].l ,lB);

        if( std::abs(T) < 1.0e-16 ) T = 0.0;
        T_shellpair.push_back(T);

      } // loop ij

//    } // else when the two shells are different 
  
  // p1: index of primitive gaussian in shell1, p2: index of primitive gaussian in shell2.
  // a1: zeta_a, a2:zeta_b, one_over_gamma: 1/(zeta_a+zeta_b) 


    if ( ( not shell1.contr[0].pure ) and ( not shell2.contr[0].pure ) ) {  
      // if both sides are cartesian, return cartesian gaussian integrals
      return {T_shellpair};
    }

    std::vector<std::vector<double>> T_shellpair_sph(1);
    T_shellpair_sph[0].assign(((2*shell1.contr[0].l+1)*(2*shell2.contr[0].l+1)),0.0);
    
    // Convert from cartesian functions to spherical functions
    cart2sph_transform(shell1.contr[0].l,shell2.contr[0].l, T_shellpair_sph[0],T_shellpair);

    return T_shellpair_sph;
  }




  /**
   *  \brief Computes a shell block of the angular momentum (magnetic dipole) matrix which is r\cross \nabla.
   *
   *
   *  \param [in] pair    Shell pair data for shell1, shell2
   *  \param [in] shell1  Bra shell
   *  \param [in] shell2  Ket shell
   *
   *  \returns Shell block of the angular momentum matrix for (shell1 | shell2)
   */  
  //compute angular momentum integrals
    
  std::vector<std::vector<double>> RealGTOIntEngine::computeAngularL(
    libint2::ShellPair &pair, libint2::Shell &shell1 , libint2::Shell &shell2 ) {
  
    int nElement = cart_ang_list[shell1.contr[0].l].size()
                    * cart_ang_list[shell2.contr[0].l].size(); 
                    //number of elements in each dimension
                      
    std::vector<std::vector<double>> L_shellpair(3);
    
    int lA[3],lB[3];
    double L[3];
    int nPGTOPair = pair.primpairs.size(); 
  
  //  RealMatrix OneixBC(3,3);//OneixBC(i,mu)
    double OneixBC[9];
    OneixBC[0*3+0] = 0.0;
    OneixBC[0*3+1] = -shell2.O[2] ; // -Bz
    OneixBC[0*3+2] = shell2.O[1];   //  By
    OneixBC[1*3+0] = shell2.O[2] ;  //  Bz
    OneixBC[1*3+1] = 0.0;
    OneixBC[1*3+2] = -shell2.O[0];  // -Bx
    OneixBC[2*3+0] = -shell2.O[1];  // -By
    OneixBC[2*3+1] = shell2.O[0];   //  Bx
    OneixBC[2*3+2] = 0.0;
   
                                      
  //  RealMatrix OneixAC(3,3);//OneixAC(i,mu)
    double OneixAC[9];
    OneixAC[0*3+0] = 0.0;
    OneixAC[0*3+1] = -shell1.O[2];  // -Az
    OneixAC[0*3+2] = shell1.O[1];   //  Ay
    OneixAC[1*3+0] = shell1.O[2];   //  Az   
    OneixAC[1*3+1] = 0.0;         
    OneixAC[1*3+2] = -shell1.O[0];  // -Ax       
    OneixAC[2*3+0] = -shell1.O[1];  // -Ay    
    OneixAC[2*3+1] = shell1.O[0];   //  Ax   
    OneixAC[2*3+2] = 0.0;         
  
    for(int i = 0; i < cart_ang_list[shell1.contr[0].l].size() ; i++) 
    for(int j = 0; j < cart_ang_list[shell2.contr[0].l].size() ; j++){
      for(int k = 0; k < 3; k++){
        lA[k] = cart_ang_list[shell1.contr[0].l][i][k];
        lB[k] = cart_ang_list[shell2.contr[0].l][j][k];  
      };
  
  //  std::cout<<"lA, lB    lA "<<lA[0]<<"\t"<<lA[1]<<"\t"<<lA[2]<<"\t lB"
  //     <<lB[0]<<"\t"<<lB[1]<<"\t"<<lB[2]<<std::endl;
  
      for( int mu = 0 ; mu < 3 ; mu++ ) {
        L[mu] = 0.0;
        for ( auto pripair : pair.primpairs ){
          L[mu] += shell1.contr[0].coeff[pripair.p1]* shell2.contr[0].coeff[pripair.p2]* 
                   Labmu(pripair,shell1,shell2,OneixAC,OneixBC,
                   shell1.contr[0].l,lA,shell2.contr[0].l,lB,mu) ;
        }
        L_shellpair[mu].push_back(L[mu]);
      }
  
    } // for j


    if ( ( not shell1.contr[0].pure ) and ( not shell2.contr[0].pure ) ) {  
      // if both sides are cartesian, return cartesian gaussian integrals
      return {L_shellpair};
    }
  
  // do cartesian to spherical transform
  
    std::vector<std::vector<double>> Angular_shellpair_sph(3);
    int k = 0;
    for ( auto cartmatrix : L_shellpair ) {
      Angular_shellpair_sph[k].assign(((2*shell1.contr[0].l+1)*(2*shell2.contr[0].l+1)),0.0);
  //    std::cout<<"xyz= "<<k<<std::endl;
      cart2sph_transform(shell1.contr[0].l,shell2.contr[0].l,Angular_shellpair_sph[k],
                           cartmatrix );
  /*
       for (auto elementofshell : Angular_shellpair_sph[k] ) 
        std::cout<<elementofshell<<std::endl;
  */
      k++;
    }
  //  return L_shellpair; 
    return Angular_shellpair_sph;
  }
  
  /**
   *  \brief Computes a shell block of the electric dipole (length gauge) matrix.
   *
   *
   *  \param [in] pair    Shell pair data for shell1, shell2
   *  \param [in] shell1  Bra shell
   *  \param [in] shell2  Ket shell
   *
   *  \returns Shell block of the electric dipole matrix for (shell1 | shell2)
   */ 
    
  std::vector<std::vector<double>> RealGTOIntEngine::computeDipoleE1( 
      libint2::ShellPair &pair, libint2::Shell &shell1 , libint2::Shell &shell2){
  
    double E1[3];
    int lA[3],lB[3];
    int nElement = cart_ang_list[shell1.contr[0].l].size()
                    * cart_ang_list[shell2.contr[0].l].size(); 
                    //number of elements in each dimension
    std::vector<std::vector<double>> tmpED1(3);
  //  std::vector<double> tmpED1x;
  //  std::vector<double> tmpED1y;
  //  std::vector<double> tmpED1z;
  
    int nPGTOPair = pair.primpairs.size(); 
  
    for(int i = 0; i < cart_ang_list[shell1.contr[0].l].size() ; i++)
    for(int j = 0; j < cart_ang_list[shell2.contr[0].l].size() ; j++){
      for(int k = 0 ; k < 3 ; k++ ){ 
        lA[k] = cart_ang_list[shell1.contr[0].l][i][k];
        lB[k] = cart_ang_list[shell2.contr[0].l][j][k];  
      }
  
      for(int mu = 0 ; mu < 3 ; mu++ ) { 
        E1[mu] = 0;
        E1[mu] = DipoleE1( pair,shell1,shell2, shell1.contr[0].l , lA , 
                                          shell2.contr[0].l , lB , mu ); 
        tmpED1[mu].push_back(E1[mu]);
  /*
        if ( RealGTOIntEngine::basisSet_->getforceCart()  ==  true ) {
           edipoleE1[bf1 + bf2*N + mu*N*N] = E1[mu];
           edipoleE1[bf2 + bf1*N + mu*N*N] = E1[mu];
        }
  */
      }
  
    }  // for j


    if ( ( not shell1.contr[0].pure ) and ( not shell2.contr[0].pure ) ) {  
      // if both sides are cartesian, return cartesian gaussian integrals
      return tmpED1;
    }

    std::vector<std::vector<double>> ED1sph(3);
    int k = 0;
    for ( auto cartmatrix : tmpED1 ) {
      ED1sph[k].assign(((2*shell1.contr[0].l+1)*(2*shell2.contr[0].l+1)),0.0);
      cart2sph_transform(shell1.contr[0].l,shell2.contr[0].l,ED1sph[k],cartmatrix);
      k++;
    }
    return ED1sph;   
  }

  /**
   *  \brief Computes a shell block of the electric dipole (velocity gauge) matrix.
   *
   *
   *  \param [in] pair    Shell pair data for shell1, shell2
   *  \param [in] shell1  Bra shell
   *  \param [in] shell2  Ket shell
   *
   *  \returns Shell block of the electric dipole matrix for (shell1 | shell2)
   */ 
    
  std::vector<std::vector<double>> RealGTOIntEngine::computeEDipoleE1_vel( 
      libint2::ShellPair &pair, libint2::Shell &shell1 , libint2::Shell &shell2){
  
    double E1[3];
    int lA[3],lB[3];
    int nElement = cart_ang_list[shell1.contr[0].l].size()
                    * cart_ang_list[shell2.contr[0].l].size(); 
                    //number of elements in each dimension
    std::vector<std::vector<double>> tmpED1(3);
  //  std::vector<double> tmpED1x;
  //  std::vector<double> tmpED1y;
  //  std::vector<double> tmpED1z;
  
    int nPGTOPair = pair.primpairs.size(); 
  
    for(int i = 0; i < cart_ang_list[shell1.contr[0].l].size() ; i++)
    for(int j = 0; j < cart_ang_list[shell2.contr[0].l].size() ; j++){
      for(int k = 0 ; k < 3 ; k++ ){ 
        lA[k] = cart_ang_list[shell1.contr[0].l][i][k];
        lB[k] = cart_ang_list[shell2.contr[0].l][j][k];  
      }
  
      for(int mu = 0 ; mu < 3 ; mu++ ) { 
        E1[mu] = 0;
        E1[mu] = Momentummu( pair,shell1,shell2, shell1.contr[0].l , lA , 
                                          shell2.contr[0].l , lB , mu ); 
        tmpED1[mu].push_back(E1[mu]);
  /*
        if ( RealGTOIntEngine::basisSet_->getforceCart()  ==  true ) {
           edipoleE1[bf1 + bf2*N + mu*N*N] = E1[mu];
           edipoleE1[bf2 + bf1*N + mu*N*N] = E1[mu];
        }
  */
      }
  
    }  // for j


    if ( ( not shell1.contr[0].pure ) and ( not shell2.contr[0].pure ) ) {  
      // if both sides are cartesian, return cartesian gaussian integrals
      return tmpED1;
    }

    std::vector<std::vector<double>> ED1sph(3);
    int k = 0;
    for ( auto cartmatrix : tmpED1 ) {
      ED1sph[k].assign(((2*shell1.contr[0].l+1)*(2*shell2.contr[0].l+1)),0.0);
      cart2sph_transform(shell1.contr[0].l,shell2.contr[0].l,ED1sph[k],cartmatrix);
      k++;
    }
    return ED1sph;   
  }

  
  /**
   *  \brief Computes a shell block of the electric quadrupole (velocity gauge) matrix.
   *
   *
   *  \param [in] pair    Shell pair data for shell1, shell2
   *  \param [in] shell1  Bra shell
   *  \param [in] shell2  Ket shell
   *
   *  \returns Shell block of the electric quadrupole matrix for (shell1 | shell2)
   */ 
   
  std::vector<std::vector<double>> RealGTOIntEngine::computeEQuadrupoleE2_vel( 
         libint2::ShellPair &pair, libint2::Shell &shell1 , libint2::Shell &shell2
         ){
  
    double E2[6]; 
                                  
    int munu[2],lA[3],lB[3];
    
    std::vector<std::vector<double>> tmpEQ2(6); 
  
      for(int i = 0; i < cart_ang_list[shell1.contr[0].l].size(); i++)
      for(int j = 0; j < cart_ang_list[shell2.contr[0].l].size(); j++){
        for(int k = 0 ; k < 3 ; k++ ){ 
          lA[k] = cart_ang_list[shell1.contr[0].l][i][k];
          lB[k] = cart_ang_list[shell2.contr[0].l][j][k];  
        }
  
/*  
        for(int mu = 0 ; mu < 3 ; mu++ ) 
        for(int nu = 0 ; nu < 3 ; nu++ ) { 
          E2[mu*3 + nu] = 0.0;
        
          E2[mu*3 + nu] = QuadrupoleE2_vel( pair,shell1,shell2,shell1.contr[0].l, lA , 
                                                         shell2.contr[0].l, lB , mu, nu ); 

        } // for nu
*/

      /* the ordering of electric quadrupole is 
                 alpha     beta
                  0          0
                  0          1
                  0          2
                  1          1
                  1          2
                  2          2
       */

        for ( int q = 0 ; q < cart_ang_list[2].size() ; q++ ) {
          int totalL = 0;
          for ( int qelement = 0 ; qelement < 3 ; qelement++ ){
            if ( cart_ang_list[2][q][qelement] == 2 ) {
              munu[0] = qelement;
              munu[1] = qelement;
              totalL += 2;
            } else if ( cart_ang_list[2][q][qelement] == 1 ) {
              munu[totalL] = qelement;
              totalL++;
            }  
          } // for qelement

          if ( totalL!= 2 ) std::cerr<<"quadrupole wrong!!"<<std::endl;
          E2[q] = QuadrupoleE2_vel( pair,shell1,shell2,shell1.contr[0].l, lA ,
                                    shell2.contr[0].l, lB ,munu[0],munu[1]);
        }
        for (auto p = 0 ; p<6 ; p++ ){
          tmpEQ2[p].push_back(E2[p]);
        } 
      } //for j


    if ( ( not shell1.contr[0].pure ) and ( not shell2.contr[0].pure ) ) {  
      // if both sides are cartesian, return cartesian gaussian integrals
      return tmpEQ2;
    }

      std::vector<std::vector<double>> EQ2sph(6);
      int kk = 0;
        for (auto cartmatrix : tmpEQ2 ){
          EQ2sph[kk].assign(((2*shell1.contr[0].l+1)*(2*shell2.contr[0].l+1)),0.0);
          cart2sph_transform(shell1.contr[0].l,shell2.contr[0].l,EQ2sph[kk],cartmatrix);
  
          kk++;
        } // for cartmatrix
    return EQ2sph; 
  
  }
  
  /**
   *  \brief Computes a shell block of the magnetic quadrupole matrix.
   *
   *
   *  \param [in] pair    Shell pair data for shell1, shell2
   *  \param [in] shell1  Bra shell
   *  \param [in] shell2  Ket shell
   *
   *  \returns Shell block of the magnetic quadrupole matrix for (shell1 | shell2)
   */ 
   
  std::vector<std::vector<double>> RealGTOIntEngine::computeMQuadrupoleM2_vel(
    libint2::ShellPair &pair, libint2::Shell &shell1 , libint2::Shell &shell2 ){
  
    double M2[9];
    int lA[3],lB[3];
    
    std::vector<std::vector<double>> tmpMQ2(9); 
    
      for(int i = 0; i < cart_ang_list[shell1.contr[0].l].size(); i++)
      for(int j = 0; j < cart_ang_list[shell2.contr[0].l].size(); j++){
        for(int k = 0 ; k < 3 ; k++ ){ 
          lA[k] = cart_ang_list[shell1.contr[0].l][i][k];
          lB[k] = cart_ang_list[shell2.contr[0].l][j][k];  
        }
  
  
        for(int mu = 0 ; mu < 3 ; mu++ ) 
        for(int nu = 0 ; nu < 3 ; nu++ ) { 
          M2[mu*3 + nu] = 0.0;
          M2[mu*3 + nu] = QuadrupoleM2_vel( pair,shell1,shell2,shell1.contr[0].l, lA ,
                                                      shell2.contr[0].l, lB , mu, nu );
  
        } //for nu
        for (auto p = 0 ; p<9 ; p++ ){
          tmpMQ2[p].push_back(M2[p]);
        }
      } //for j


    if ( ( not shell1.contr[0].pure ) and ( not shell2.contr[0].pure ) ) {  
      // if both sides are cartesian, return cartesian gaussian integrals
      return tmpMQ2;
    }

        
      std::vector<std::vector<double>> MQ2sph(9);
      int kk = 0;
      for (auto cartmatrix : tmpMQ2 ){
        MQ2sph[kk].assign(((2*shell1.contr[0].l+1)*(2*shell2.contr[0].l+1)),0.0);
        cart2sph_transform(shell1.contr[0].l,shell2.contr[0].l,MQ2sph[kk],cartmatrix);
        kk++;
      }
      return MQ2sph;
  }
  
  /**
   *  \brief Computes a shell block of the electric octupole (velocity gauge) matrix.
   *
   *
   *  \param [in] pair    Shell pair data for shell1, shell2
   *  \param [in] shell1  Bra shell
   *  \param [in] shell2  Ket shell
   *
   *  \returns Shell block of the electric octupole matrix for (shell1 | shell2)
   */ 
   
  std::vector<std::vector<double>> RealGTOIntEngine::computeEOctupoleE3_vel(
          libint2::ShellPair &pair, libint2::Shell &shell1 , libint2::Shell &shell2 ){
    double E3[10];
    int lA[3],lB[3],alphabetagamma[3];
  
    std::vector<std::vector<double>> tmpEO3(10);
  
      for(int i = 0; i < cart_ang_list[shell1.contr[0].l].size(); i++)
      for(int j = 0; j < cart_ang_list[shell2.contr[0].l].size(); j++){
        for(int k = 0 ; k < 3 ; k++ ){ 
          lA[k] = cart_ang_list[shell1.contr[0].l][i][k];
          lB[k] = cart_ang_list[shell2.contr[0].l][j][k];  
        }
/*
        for( int alpha = 0 ; alpha < 3 ; alpha++ ) 
        for( int beta  = 0 ; beta  < 3 ; beta ++ ) 
        for( int gamma = 0 ; gamma < 3 ; gamma++ ) {
   
          E3[alpha*9 + beta*3 +gamma] = 0.0;
  
          E3[alpha*9 + beta*3 +gamma] = OctupoleE3_vel( pair,shell1,shell2, 
          shell1.contr[0].l , lA , shell2.contr[0].l , lB , alpha, beta, gamma ); 
   
        } // for gamma
*/

/*
 * orderring of octupole
   alpha  beta  gamma
   0      0     0
   0      0     1
   0      0     2
   0      1     1 
   0      1     2
   0      2     2
   1      1     1
   1      1     2 
   1      2     2
   2      2     2
 */

// std::cerr<<"start a cycle"<<std::endl;

        for ( int q = 0 ; q < cart_ang_list[3].size() ; q++ ) {
          int totalL = 0;
          for ( int qelement = 0 ; qelement < 3 ; qelement++ ){
            if ( cart_ang_list[3][q][qelement] == 3 ) {
              alphabetagamma[0] = qelement;
              alphabetagamma[1] = qelement;
              alphabetagamma[2] = qelement;
              totalL += 3;
            } else if ( cart_ang_list[3][q][qelement] == 2 ) {
              alphabetagamma[totalL] = qelement;
              totalL++;
              alphabetagamma[totalL] = qelement;
              totalL++;
            } else if ( cart_ang_list[3][q][qelement] == 1 ) {
              alphabetagamma[totalL] = qelement;
              totalL++;
            }  
          } // for qelement

//std::cerr<<" alpha = "<<alphabetagamma[0]<<" beta = "<<alphabetagamma[1]<<" gamma = "<<alphabetagamma[2]<<std::endl;

          if ( totalL!= 3 ) std::cerr<<"octupole wrong!!"<<std::endl;
          E3[q] = OctupoleE3_vel( pair,shell1,shell2,
                  shell1.contr[0].l , lA , shell2.contr[0].l , lB , 
                  alphabetagamma[0],alphabetagamma[1],alphabetagamma[2] );
        } // for q  

        for ( int p = 0 ; p < 10 ; p++ ){
          tmpEO3[p].push_back(E3[p]);
        }
      } // for j
  

    if ( ( not shell1.contr[0].pure ) and ( not shell2.contr[0].pure ) ) {  
      // if both sides are cartesian, return cartesian gaussian integrals
      return tmpEO3;
    }

      std::vector<std::vector<double>> EO3sph(10);
      int kk = 0;
      for (auto cartmatrix : tmpEO3) {
        EO3sph[kk].assign(((2*shell1.contr[0].l+1)*(2*shell2.contr[0].l+1)),0.0);
        cart2sph_transform(shell1.contr[0].l,shell2.contr[0].l,EO3sph[kk],cartmatrix); 
        kk++;
      }
      return EO3sph;
  }
  
  /**
   *  \brief Computes a shell block of the nuclear potential matrix.
   *
   *
   *  \param [in] nucShell nuclear shell, give the exponents of gaussian function of nuclei
   *  \param [in] pair    Shell pair data for shell1, shell2
   *  \param [in] shell1  Bra shell
   *  \param [in] shell2  Ket shell
   *
   *  \returns Shell block of the potential integral matrix for (shell1 | shell2)
   */ 
   
  std::vector<std::vector<double>> RealGTOIntEngine::computePotentialV(
    const std::vector<libint2::Shell> &nucShell, libint2::ShellPair &pair, 
    libint2::Shell &shell1 , libint2::Shell &shell2 , const Molecule& molecule){
  

    bool useFiniteWidthNuclei = nucShell.size() > 0;


    std::vector<double> potential_shellpair;
    double V,tmpV;
    int lA[3],lB[3];
  
    for(int i = 0; i < cart_ang_list[shell1.contr[0].l].size() ; i++) 
    for(int j = 0; j < cart_ang_list[shell2.contr[0].l].size() ; j++){
      for(int k = 0; k < 3; k++){
        lA[k] = cart_ang_list[shell1.contr[0].l][i][k];
        lB[k] = cart_ang_list[shell2.contr[0].l][j][k];  
      };
      V = 0.0;
  
      V = hRRVab(nucShell,pair,shell1,shell2,shell1.contr[0].l,lA,shell2.contr[0].l,lB,molecule);
        // calculate from hRRVab
      
      int iAtom;
  /*
      tmpV = 0.0;
      for ( auto pripair : pair.primpairs ){
        iAtom = 0;
        for ( auto atom : molecule_.atoms ){
          double C[3];
          for ( int mu = 0 ; mu < 3 ; mu++ ) C[mu] = atom.coord[mu];
          auto norm = shell1.contr[0].coeff[pripair.p1]* 
                      shell2.contr[0].coeff[pripair.p2];
  
          tmpV+= atom.atomicNumber * norm *
                 hRRiPPVab(pripair,shell1,shell2,shell1.contr[0].l,lA,shell2.contr[0].l,lB,C,0,iAtom);
          iAtom++;
        }  
      }
      V = tmpV;
  */  // use hRRiPPVab to calculate 
      potential_shellpair.push_back(-V);
  
    } // for j
   
    if ( ( not shell1.contr[0].pure ) and ( not shell2.contr[0].pure ) ) {  
      // if both sides are cartesian, return cartesian gaussian integrals
      return {potential_shellpair};
    }
 
    std::vector<double> V_shellpair_sph;
  
    V_shellpair_sph.assign(((2*shell1.contr[0].l+1)*(2*shell2.contr[0].l+1)),0.0);
  
    cart2sph_transform(shell1.contr[0].l,shell2.contr[0].l,
                       V_shellpair_sph,potential_shellpair );
  
    return { V_shellpair_sph }; 
  
  }
  
  /**
   *  \brief Computes a shell block of the spin orbit integral matrix.
   *
   *
   *  \param [in] nucShell nuclear shell, give the exponents of gaussian function of nuclei
   *  \param [in] pair    Shell pair data for shell1, shell2
   *  \param [in] shell1  Bra shell
   *  \param [in] shell2  Ket shell
   *
   *  \returns Shell block of the spin orbit integral matrix for (shell1 | shell2)
   */  
   
  std::vector<std::vector<double>> RealGTOIntEngine::computeSL(
    const std::vector<libint2::Shell> &nucShell, libint2::ShellPair &pair, 
    libint2::Shell &shell1 , libint2::Shell &shell2 , const Molecule &molecule){


    bool useFiniteWidthNuclei = nucShell.size() > 0;
  
    std::vector<std::vector<double>> SL_shellpair(3);
    
    int lA[3],lB[3],iAtom;
    double Sl[3],C[3],SlC;
  
    double OneixBC[9];
    double OneixAC[9];
  
    for(int i = 0; 
      i < cart_ang_list[shell1.contr[0].l].size() ; i++) 
    for(int j = 0;
      j < cart_ang_list[shell2.contr[0].l].size() ; j++){
      for(int k = 0; k < 3; k++){
        lA[k] = cart_ang_list[shell1.contr[0].l][i][k];
        lB[k] = cart_ang_list[shell2.contr[0].l][j][k];  
      };
  
  
      for( int mu = 0 ; mu < 3 ; mu++ ) {
        Sl[mu] = 0.0;
        for ( auto pripair : pair.primpairs ){
          SlC = 0.0;
          iAtom = 0;
          for( auto atom: molecule.atoms ) {
            for( int m=0 ; m<3 ; m++ ) C[m] = atom.coord[m];
  
          //  RealMatrix OneixBC(3,3);//OneixBC(i,mu)
            OneixBC[0*3+0] = 0.0;
            OneixBC[0*3+1] = -(shell2.O[2] - C[2]) ; // -Bz
            OneixBC[0*3+2] = (shell2.O[1] - C[1]);   //  By
            OneixBC[1*3+0] = (shell2.O[2] - C[2]);  //  Bz
            OneixBC[1*3+1] = 0.0;
            OneixBC[1*3+2] = -(shell2.O[0] - C[0]);  // -Bx
            OneixBC[2*3+0] = -(shell2.O[1] - C[1]);  // -By
            OneixBC[2*3+1] = (shell2.O[0] - C[0]);   //  Bx
            OneixBC[2*3+2] = 0.0;
       
                                              
                //  RealMatrix OneixAC(3,3);//OneixAC(i,mu)
            OneixAC[0*3+0] = 0.0;
            OneixAC[0*3+1] = -(shell1.O[2] - C[2]);  // -Az
            OneixAC[0*3+2] = (shell1.O[1] - C[1]);   //  Ay
            OneixAC[1*3+0] = (shell1.O[2] - C[2]);   //  Az   
            OneixAC[1*3+1] = 0.0;         
            OneixAC[1*3+2] = -(shell1.O[0] - C[0]);  // -Ax       
            OneixAC[2*3+0] = -(shell1.O[1] - C[1]);  // -Ay    
            OneixAC[2*3+1] = (shell1.O[0] - C[0]);   //  Ax   
            OneixAC[2*3+2] = 0.0;         
          
            SlC += shell1.contr[0].coeff[pripair.p1]* 
              shell2.contr[0].coeff[pripair.p2]* 
              atom.nucCharge
              * Slabmu(nucShell,pripair,shell1,shell2,OneixAC,OneixBC,
              shell1.contr[0].l,lA,shell2.contr[0].l,lB,mu,0,iAtom,molecule) ;
   
            iAtom++; 
  
          } // for atoms
  
          Sl[mu] += SlC; 
        } // for primitive shell pairs

        // Negation takes delV x del -> pV x p
        SL_shellpair[mu].push_back(-Sl[mu]);

      } // for mu
  
    } // for j


    if ( ( not shell1.contr[0].pure ) and ( not shell2.contr[0].pure ) ) {  
      // if both sides are cartesian, return cartesian gaussian integrals
      return SL_shellpair;
    }
  
  // do cartesian to spherical transform
  
    std::vector<std::vector<double>> SOCoupling_shellpair_sph(3);
    int k = 0;
    for ( auto cartmatrix : SL_shellpair ) {
      SOCoupling_shellpair_sph[k].assign(((2*shell1.contr[0].l+1)*(2*shell2.contr[0].l+1)),0.0);
      cart2sph_transform(shell1.contr[0].l,shell2.contr[0].l,
                         SOCoupling_shellpair_sph[k],cartmatrix );
      k++;
    } // loop x,y,z
  
    return SOCoupling_shellpair_sph;
  }
  
  /**
   *  \brief Computes a shell block of the pV dot p matrix.
   *
   *
   *  \param [in] nucShell nuclear shell, give the exponents of gaussian function of nuclei
   *  \param [in] pair    Shell pair data for shell1, shell2
   *  \param [in] shell1  Bra shell
   *  \param [in] shell2  Ket shell
   *
   *  \returns Shell block of the pV dot p matrix for (shell1 | shell2)
   */ 
   
  std::vector<std::vector<double>> RealGTOIntEngine::computepVdotp( 
    const std::vector<libint2::Shell> &nucShell, libint2::ShellPair &pair, 
    libint2::Shell &shell1, libint2::Shell &shell2, const Molecule& molecule){

    bool useFiniteWidthNuclei = nucShell.size() > 0;
  
    std::vector<double> pVdotp_shellpair;
    double pVp,pVpC,C[3];
    int lA[3],lB[3];
  
    for(int i = 0; i < cart_ang_list[shell1.contr[0].l].size() ; i++) 
    for(int j = 0; j < cart_ang_list[shell2.contr[0].l].size() ; j++){
      for(int k = 0; k < 3; k++){
        lA[k] = cart_ang_list[shell1.contr[0].l][i][k];
        lB[k] = cart_ang_list[shell2.contr[0].l][j][k];  
      };
     
      int iAtom;
  
      pVp = 0.0;
      for ( auto pripair : pair.primpairs ){
        iAtom = 0;
        pVpC = 0.0;
        for ( auto atom : molecule.atoms ){
          for ( int mu = 0 ; mu < 3 ; mu++ ) C[mu] = atom.coord[mu];
          auto norm = shell1.contr[0].coeff[pripair.p1]* 
                      shell2.contr[0].coeff[pripair.p2];
  
          pVpC += atom.nucCharge * norm *
                 pVpab(nucShell,pripair,shell1,shell2,shell1.contr[0].l,lA,
                       shell2.contr[0].l,lB,0,iAtom,molecule);
          iAtom++;
        } // atoms
        pVp += pVpC;  
      } // primpairs

      // Negation takes delV.del -> pV.p
      pVdotp_shellpair.push_back(-pVp);
  
    } // for j


    if ( ( not shell1.contr[0].pure ) and ( not shell2.contr[0].pure ) ) {  
      // if both sides are cartesian, return cartesian gaussian integrals
      return {pVdotp_shellpair};
    }
    
    std::vector<std::vector<double>> pVdotp_shellpair_sph(1);
  
    pVdotp_shellpair_sph[0].assign(((2*shell1.contr[0].l+1)*(2*shell2.contr[0].l+1)),0.0);
  
    cart2sph_transform(shell1.contr[0].l,shell2.contr[0].l,
                       pVdotp_shellpair_sph[0],pVdotp_shellpair );
  
    return pVdotp_shellpair_sph; 
  
  }
  
  //------------------------------------//
  // overlap horizontal recursion       //
  // (a|b) = (A-B)(a|b-1) + (a+1|b-1)   //
  //------------------------------------//

  /**
   *  \brief Perform the horizontal recurrence relation for the contracted overlap integral
   *
   *  (a|b) = (A-B)(a|b-1) + (a+1|b-1)
   *
   *  where a,b are the angular momentum, A,B are the nuclear coordinates.
   *
   *  \param [in] pair    Shell pair data for shell1, shell2
   *  \param [in] shell1  Bra shell
   *  \param [in] shell2  Ket shell
   *  \param [in] LA      total Bra angular momentum
   *  \param [in] lA      Bra angular momentum vector (lAx,lAy,lAz)
   *  \param [in] LB      total Ket angular momentum
   *  \param [in] lB      Ket angular momentum vector (lBx,lBy,lBz)
   *
   *  \returns a contracted overlap integral
   *
   */ 
   
  double RealGTOIntEngine::hRRSab(
    libint2::ShellPair &pair, libint2::Shell &shell1, libint2::Shell &shell2, int LA, int *lA ,int LB, int *lB) {
    int iWork,lAp1[3],lBm1[3];
    double tmpVal = 0.0;
  
    if ( LB > LA ) {     // if LB>LA, use horizontal recursion to make LA>LB. 
      for(iWork = 0 ; iWork < 3 ; iWork++ ){
        lAp1[iWork] = lA[iWork];     
        lBm1[iWork] = lB[iWork];     
      };
    
      if ( lB[0] > 0 )      iWork = 0;
      else if ( lB[1] > 0 ) iWork = 1;
      else if ( lB[2] > 0 ) iWork = 2;
      lAp1[iWork]++;
      lBm1[iWork]--;
   
      tmpVal = hRRSab( pair , shell1, shell2, LA+1 , lAp1 , LB-1 , lBm1 );
      tmpVal+= pair.AB[iWork]*hRRSab( pair , shell1, shell2, LA , lA , LB-1 , lBm1 );
      return tmpVal;
    
    }
   
    if(LA == 0) {
      // (s|s)
      for( auto &pripair : pair.primpairs ) { 
        tmpVal += shell1.contr[0].coeff[pripair.p1]* shell2.contr[0].coeff[pripair.p2]*
                  pow(sqrt(M_PI),3) * sqrt(pripair.one_over_gamma)*pripair.K ;
  //tmpVal+= pow(M_PI,3/2) * sqrt(pripair.one_over_gamma)*pripair.K ;
  //std::cout<<"p1 "<<pripair.p1<<" p2 "<<pripair.p2<<" primitive ss"<<shell1.contr[0].coeff[pripair.p1]* shell2.contr[0].coeff[pripair.p2]* 
  // pow(M_PI,3/2) * sqrt(pripair.one_over_gamma)*pripair.K<<"pi^3/2  "<<pow(sqrt(M_PI),3)<<std::endl;
      }
      return tmpVal;
    }  else if(LB == 0) {
      // (|s)
      for( auto pripair : pair.primpairs ) {
      
        tmpVal+= shell1.contr[0].coeff[pripair.p1]* shell2.contr[0].coeff[pripair.p2]*
                  vRRSa0(pripair,shell1,LA,lA) ;
  //      std::cout<<"LA2LB0 "<<shell1.contr[0].coeff[pripair.p1]* shell2.contr[0].coeff[pripair.p2]
  //       <<"\t"<<LA<<"\t"<<vRRSa0(pripair,shell1,LA,lA)<<std::endl;
      }  
      return tmpVal;
    };
  
    for(iWork = 0;iWork < 3;iWork++){
      lAp1[iWork]=lA[iWork];     
      lBm1[iWork]=lB[iWork];     
    };
    if (lB[0] > 0)      iWork = 0;
    else if (lB[1] > 0) iWork=1;
    else if (lB[2] > 0) iWork=2;
    lAp1[iWork]++;
    lBm1[iWork]--;
  
  
    tmpVal = hRRSab(pair,shell1,shell2,LA+1,lAp1,LB-1,lBm1);
    tmpVal+= (shell1.O[iWork]-shell2.O[iWork])*hRRSab(pair,shell1,shell2,LA,lA,LB-1,lBm1);
  //  tmpVal+= pair.AB[iWork]*hRRSab(pair,shell1,shell2,LA,lA,LB-1,lBm1);
    return tmpVal;
   
  }
  
  //----------------------------------------------------------//
  // overlap vertical recursion                               //
  // (a|0) = (P-A)(a-1|0) + halfInvZeta*N_(a-1)*(a-2|0)       //
  //----------------------------------------------------------//
  /**
   *  \brief Perform the vertical recurrence relation for the uncontracted overlap integral 
   *
   *  (a|0) = (P-A)(a-1|0) + 1/2 *1/Zeta * N_(a-1)*(a-2|0)
   *
   *  where a is angular momentum, Zeta=zeta_a+zeta_b, A is bra nuclear coordinate. 
   *  P = (zeta_a*A+zeta_b*B)/Zeta
   *
   *  \param [in] pripair Primitive Shell pair data for shell1, shell2
   *  \param [in] shell1  Bra shell
   *  \param [in] LA      total Bra angular momentum
   *  \param [in] lA      Bra angular momentum vector (lAx,lAy,lAz)
   *
   *  \returns an uncontracted overlap integral
   *
   */ 
   
  double RealGTOIntEngine::vRRSa0(
    libint2::ShellPair::PrimPairData &pripair, libint2::Shell &shell1, int LA, int *lA){

  //notice: vRRSa0 doesn't include contraction coeffs. it is givin in hRRSab.
    double tmpVal=0.0;
    if (LA == 0) {   //[s||s]
      tmpVal = pow(sqrt(M_PI),3) * sqrt(pripair.one_over_gamma)*pripair.K ;
      return tmpVal;
    }
  
    int iWork;
    int lAm1[3];
    for(iWork = 0;iWork < 3;iWork++) lAm1[iWork]=lA[iWork];
    if (lA[0] > 0) iWork = 0;
    else if (lA[1] > 0) iWork=1;
    else if (lA[2] > 0) iWork=2;
  
  /*
    if(LA == 1) {
     tmpVal = (pripair.P[iWork]-shell1.O[iWork])*pow(sqrt(M_PI),3) * sqrt(pripair.one_over_gamma)*pripair.K ; 
     return tmpVal;
    }
  */
  
    lAm1[iWork]--;
    tmpVal = (pripair.P[iWork]-shell1.O[iWork])*vRRSa0(pripair,shell1,LA-1,lAm1); 
  
  //  if(LA == 2&&lA[iWork] == 2) 
  //    tmpVal += 1/2*pripair.one_over_gamma
  //            *pow(sqrt(M_PI),3) * sqrt(pripair.one_over_gamma)*pripair.K;
  //  else if (lA[iWork] >=2) {
      if ( lA[iWork] >=2 ) {
        lAm1[iWork]--; 
        tmpVal += (lA[iWork]-1)*0.5*pripair.one_over_gamma*vRRSa0(pripair,shell1,LA-2,lAm1);
  //  }
      }
    return tmpVal;
   
  }
  
  //---------------------------------------------//
  // overlap horizontal recursion (iPP specific) //
  // (a|b) = (A-B)(a|b-1) + (a+1|b-1)            //
  //   LA  >=  LB                                //
  //---------------------------------------------//
  /**
   *  \brief Perform the horizontal recurrence relation for the uncontracted overlap integral
   *
   *  (a|b) = (A-B)(a|b-1) + (a+1|b-1)
   *
   *  where a,b are the angular momentum, A,B are the nuclear coordinates.
   *
   *  \param [in] pripair Primitive Shell pair data for shell1, shell2
   *  \param [in] shell1  Bra shell
   *  \param [in] shell2  Ket shell
   *  \param [in] LA      total Bra angular momentum
   *  \param [in] lA      Bra angular momentum vector (lAx,lAy,lAz)
   *  \param [in] LB      total Ket angular momentum
   *  \param [in] lB      Ket angular momentum vector (lBx,lBy,lBz)
   *
   *  \returns a uncontracted overlap integral. It doesn't include contraction coeffs.
   *
   */ 
   
  double RealGTOIntEngine::hRRiPPSab(
    libint2::ShellPair::PrimPairData &pripair, libint2::Shell &shell1, 
    libint2::Shell &shell2, int LA,int *lA,int LB,int *lB) {

    double tmpVal = 0.0;
    int iWork,lAp1[3],lBm1[3];
  
  // SS start this part may be redundant
    if (LB>LA) {
      for(iWork = 0;iWork < 3;iWork++){
        lAp1[iWork]=lA[iWork];     
        lBm1[iWork]=lB[iWork];     
      };
      if (lB[0] > 0)      iWork = 0;
      else if (lB[1] > 0) iWork=1;
      else if (lB[2] > 0) iWork=2;
      lAp1[iWork]++;
      lBm1[iWork]--;
  
      tmpVal = RealGTOIntEngine::hRRiPPSab(pripair,shell1,shell2,LA+1,lAp1,LB-1,lBm1);
      tmpVal+= (shell1.O[iWork]-shell2.O[iWork])
                *RealGTOIntEngine::hRRiPPSab(pripair,shell1,shell2,LA,lA,LB-1,lBm1);
      return tmpVal;
    } 
  // SS end
  
    if((LA+LB) == 0) {
      tmpVal = pow(sqrt(M_PI),3) * sqrt(pripair.one_over_gamma)*pripair.K; 
             //doesn't include contraction coefficients
      return tmpVal;
    } else if(LB == 0) {
      tmpVal = RealGTOIntEngine::vRRSa0(pripair,shell1,LA,lA);
      return tmpVal;
    };
  
    for(iWork = 0;iWork < 3;iWork++){
      lAp1[iWork]=lA[iWork];     
      lBm1[iWork]=lB[iWork];     
    };
    if (lB[0] > 0)      iWork = 0;
    else if (lB[1] > 0) iWork=1;
    else if (lB[2] > 0) iWork=2;
    lAp1[iWork]++;
    lBm1[iWork]--;
  
    tmpVal = RealGTOIntEngine::hRRiPPSab(pripair,shell1,shell2,LA+1,lAp1,LB-1,lBm1);
    tmpVal+= (shell1.O[iWork]-shell2.O[iWork])
       *RealGTOIntEngine::hRRiPPSab(pripair,shell1,shell2,LA,lA,LB-1,lBm1);
    return tmpVal;
  };
  


//----------------------------------------------//
//kinetic integrals                             //
//
//[a|T|b]=[a|Tx|b]+[a|Ty|b]+[a|Tz|b]            //
//[a|Ti|b]=beta*(2bi+1)[a|b]-2beta^2 *[a|b+2i]-1/2*bi*(bi-1)[a|b-2i]//
// LA>=LB?                                      //
//S.S.                                          //
//----------------------------------------------//
   
  double RealGTOIntEngine::RRTab(
    libint2::ShellPair &pair, libint2::Shell &shell1, libint2::Shell &shell2, int LA, int *lA ,int LB, int *lB ) {

  double tmpVal = 0.0;
  if ( (LA==0) and (LB==0) ) {
    // (s|T|s)
    for( auto &pripair : pair.primpairs ) {
      
      auto ssS = shell1.contr[0].coeff[pripair.p1]* shell2.contr[0].coeff[pripair.p2]
         * pow(sqrt(M_PI),3) * sqrt(pripair.one_over_gamma)*pripair.K ;
      
      auto Xi = shell1.alpha[pripair.p1] * shell2.alpha[pripair.p2]
                                    *pripair.one_over_gamma; 

      auto ABsquare = 0.0;
      for ( int mu = 0 ; mu < 3 ; mu++ ) {
        ABsquare += pow(pair.AB[mu],2);  
      } // for mu in AB square
      
      auto ssT = Xi*(3.0-2.0*Xi*ABsquare)*ssS;    
      tmpVal += ssT; 
    } // for( auto &pripair : pair.primpairs )
    return tmpVal;
  } // if ( (LA==0) and (LB==0) ) 

  int iWork,lBm[3];
 
  for ( iWork=0 ; iWork<3 ; iWork++ ) { 
    for( int i=0 ; i<3 ; i++ ) lBm[i]=lB[i];     
    lBm[iWork] += 2; 

    for( auto &pripair : pair.primpairs ) {
      auto tmploop = 0.0 ;
      tmploop -= 2.0 * pow( shell2.alpha[pripair.p2], 2 ) 
        * hRRiPPSab( pripair, shell1, shell2, LA, lA, LB+2, lBm );

      tmploop += shell2.alpha[pripair.p2] * ( 2 * lB[iWork] + 1 ) 
        * hRRiPPSab( pripair, shell1, shell2, LA, lA, LB, lB );

      tmpVal += shell1.contr[0].coeff[pripair.p1]* shell2.contr[0].coeff[pripair.p2]
        * tmploop ; 
    }
    if(lB[iWork]>=2) {
      lBm[iWork] -= 4; 
      for( auto &pripair : pair.primpairs ) {

        tmpVal -= shell1.contr[0].coeff[pripair.p1]* shell2.contr[0].coeff[pripair.p2]
          * 0.5 * lB[iWork] * (lB[iWork]-1) 
          * hRRiPPSab( pripair, shell1, shell2, LA, lA, LB-2, lBm );
      } // for( auto &pripair : pair.primpairs )
    } // if(lB[iWork]>=2) 
  } // for(iWork=0;iWork<3;iWork++) 
  return tmpVal;
} // RealGTOIntEngine::RRTab
  

//---------------------------------------------------------//
// kinetic vertical recursion                              //
// (a|T|b) = (P-B)(a|T|b-1) + halfInvZeta*N_(b-1)(a|T|b-2) //
//         + halfInvZeta*N_(a)(a-1|T|b-1)                  //
//         + 2*Xi*[(a|b) - halfInvzeta_b*N_(b-1)(a|b-2)]     //
//                                                         //
// (a|T|0) = (P-A)(a-1|T|0) + halfInvZeta*N_(a-1)(a-2|T|0) //
//         + 2*Xi*[(a|0) - halfInvZeta_a*N_(a-1)(a-2|0)]     //
//---------------------------------------------------------//
   
  double RealGTOIntEngine::vRRTab(
    libint2::ShellPair &pair, libint2::Shell &shell1, libint2::Shell &shell2, int LA, int *lA ,int LB, int *lB ) {

  double tmpVal = 0.0;
  // int i,j,iPP;

  if ( (LA==0) and (LB==0) ) {
    for ( auto &pripair : pair.primpairs ) {

//  tmpVal += ijSP->ssT[iPP];
      auto ssS = shell1.contr[0].coeff[pripair.p1]* shell2.contr[0].coeff[pripair.p2]
         * pow(sqrt(M_PI),3) * sqrt(pripair.one_over_gamma)*pripair.K ;
      
      auto Xi = shell1.alpha[pripair.p1] * shell2.alpha[pripair.p2]
                                    *pripair.one_over_gamma; 

      auto ABsquare = 0.0;
      for ( int mu = 0 ; mu < 3 ; mu++ ) {
        ABsquare += pow(pair.AB[mu],2);  
      } // for mu in AB square
      
      auto ssT = Xi*(3.0-2.0*Xi*ABsquare)*ssS;    
      tmpVal += ssT; 
    } // for pripair
    return tmpVal;
  } else if(LB==0) {
    for(auto &pripair : pair.primpairs ) {
      tmpVal += shell1.contr[0].coeff[pripair.p1]* shell2.contr[0].coeff[pripair.p2]
        * vRRiPPTab( pripair, shell1, shell2, LA, lA, LB, lB );
    } // for pripair
    return tmpVal;
  } // else if LB == 0

  // here LB>0 , which include LA>0 and LA==0

  int iWork,lAm1[3],lBm1[3];
  for( iWork=0 ; iWork<3 ; iWork++ ){
    lAm1[iWork]=lA[iWork];
    lBm1[iWork]=lB[iWork];
  } // for iWork

  if (lB[0]>0)      iWork=0;
  else if (lB[1]>0) iWork=1;
  else if (lB[2]>0) iWork=2;
  lBm1[iWork]--;

  for( auto &pripair : pair.primpairs )
    tmpVal += ( pripair.P[iWork] - shell2.O[iWork] )
      * shell1.contr[0].coeff[pripair.p1]* shell2.contr[0].coeff[pripair.p2]
      * vRRiPPTab( pripair, shell1, shell2, LA, lA, LB-1, lBm1 ); 

  if(lA[iWork]>0) {
    lAm1[iWork]--;
    for( auto &pripair : pair.primpairs )
      tmpVal += (lAm1[iWork]+1) * 0.5 * pripair.one_over_gamma 
        * shell1.contr[0].coeff[pripair.p1]* shell2.contr[0].coeff[pripair.p2]
        * vRRiPPTab( pripair, shell1, shell2, LA-1, lAm1, LB-1, lBm1 );
  } // if lA[iWork]>0

  if(lB[iWork]>=2) {
    lBm1[iWork]--;
    auto tmploop=0.0;
    for( auto &pripair : pair.primpairs ) {
      tmploop = 0.0;
      tmploop += (lBm1[iWork]+1) * 0.5 * pripair.one_over_gamma 
        * vRRiPPTab( pripair, shell1, shell2, LA, lA, LB-2, lBm1 );

      auto Xi = shell1.alpha[pripair.p1] * shell2.alpha[pripair.p2]
                                    *pripair.one_over_gamma; 

      tmploop -= (lBm1[iWork]+1) * 2.0 * Xi * 0.5 /shell2.alpha[pripair.p2] 
      // pripair.one_over_gamma 
        * hRRiPPSab( pripair, shell1, shell2, LA, lA, LB-2, lBm1 );
      tmpVal += tmploop * shell1.contr[0].coeff[pripair.p1]
        * shell2.contr[0].coeff[pripair.p2];
    } // for pripair
    
  } // if(lB[iWork]>=2) 

  for( auto &pripair : pair.primpairs ) {

    auto Xi = shell1.alpha[pripair.p1] * shell2.alpha[pripair.p2]
                                  *pripair.one_over_gamma; 

    tmpVal += 2.0 * Xi * shell1.contr[0].coeff[pripair.p1] 
      * shell2.contr[0].coeff[pripair.p2]
      * hRRiPPSab( pripair, shell1, shell2, LA, lA, LB, lB );
  }
  return tmpVal;
};  // RealGTOIntEngine::vRRTab


//---------------------------------------------------------//
// kinetic vertical recursion (iPP specific)               //
// (a|T|b) = (P-B)(a|T|b-1) + halfInvZeta*N_(b-1)(a|T|b-2) //
//         + halfInvZeta*N_(a)(a-1|T|b-1)                  //
//         + 2*Xi*[(a|b) - halfInvZetai_b*N_(b-1)(a|b-2)]     //
//                                                         //
// (a|T|0) = (P-A)(a-1|T|0) + halfInvZeta*N_(a-1)(a-2|T|0) //
//         + 2*Xi*[(a|0) - halfInvZetai_a*N_(a-1)(a-2|0)]     //
//---------------------------------------------------------//
   
  double RealGTOIntEngine::vRRiPPTab(libint2::ShellPair::PrimPairData &pripair, 
    libint2::Shell &shell1, libint2::Shell &shell2,int LA, int *lA, int LB, int *lB ) {

  double tmpVal = 0.0;
 
  auto Xi = shell1.alpha[pripair.p1] * shell2.alpha[pripair.p2]
                              *pripair.one_over_gamma; 


  if ( (LA==0) and (LB==0) ) {

    // seems the coefficients are given in vRRTab, so we don't include coefficient
    //auto ssS = shell1.contr[0].coeff[pripair.p1]* shell2.contr[0].coeff[pripair.p2]

    auto ssS = pow(sqrt(M_PI),3) * sqrt(pripair.one_over_gamma)*pripair.K ;
         auto ABsquare = 0.0;

    for ( int mu = 0 ; mu < 3 ; mu++ ) {
      ABsquare += pow( shell1.O[mu] - shell2.O[mu], 2 );   
    } // for mu in AB square

    auto ssT = Xi*(3.0-2.0*Xi*ABsquare)*ssS;    

    tmpVal += ssT; 
    return tmpVal;

  } else if(LB==0) {

    int iWork,lAm1[3];
    for(iWork=0;iWork<3;iWork++) lAm1[iWork]=lA[iWork];
    if (lA[0]>0)      iWork=0;
    else if (lA[1]>0) iWork=1;
    else if (lA[2]>0) iWork=2;
    lAm1[iWork]--;
//    cout<<"xsli test 1 :"<<iWork<<" "<<ijSP->PA[iPP][iWork]<<endl;

    tmpVal += ( pripair.P[iWork] - shell1.O[iWork] ) 
      * vRRiPPTab( pripair, shell1, shell2, LA-1, lAm1, LB, lB );

//    cout<<"xsli test 2 :"<<tmpVal<<endl;

//    if(lA[iWork]==2) {
// here we don't distinguish lA == 2 
//      tmpVal += 0.5 * pripair.one_over_gamma * ssT[iPP];
//    cout<<"xsli test 3 :"<<tmpVal<<endl;
//      tmpVal -= math.two*ijSP->Xi[iPP]*ijSP->halfInvZeta[iPP]*ijSP->ss[iPP];
//    cout<<"xsli test 4 :"<<tmpVal<<endl;

//    } else if(lA[iWork]>2) {
    
    if ( lA[iWork] >= 2 ) {
      lAm1[iWork]--;
      tmpVal += (lAm1[iWork]+1) * 0.5 * pripair.one_over_gamma 
        * vRRiPPTab( pripair, shell1, shell2, LA-2, lAm1, LB, lB );
      tmpVal -= (lAm1[iWork]+1) * 2.0 * Xi * 0.5 / shell1.alpha[pripair.p1] 
        * vRRSa0( pripair, shell1, LA-2, lAm1 );

    } // if ( lA[iWork] >= 2 )
    tmpVal += 2.0 * Xi * vRRSa0( pripair, shell1, LA, lA );
//    cout<<"xsli test 5 :"<<tmpVal<<endl;
    return tmpVal;
  } // else if(LB==0) 

  // here LB>0 , which include LA>0 and LA == 0 
  int iWork,lAm1[3],lBm1[3];
  for( iWork=0 ; iWork<3 ; iWork++ ){
    lAm1[iWork]=lA[iWork];
    lBm1[iWork]=lB[iWork];
  };
  if (lB[0]>0)      iWork=0;
  else if (lB[1]>0) iWork=1;
  else if (lB[2]>0) iWork=2;
  lBm1[iWork]--;

  tmpVal += ( pripair.P[iWork] - shell2.O[iWork] ) 
    * vRRiPPTab( pripair, shell1, shell2, LA, lA, LB-1, lBm1 );

  if(lA[iWork]>0) {
    lAm1[iWork]--;
    tmpVal += ( lAm1[iWork] + 1 ) * 0.5 * pripair.one_over_gamma 
      * vRRiPPTab( pripair, shell1, shell2, LA-1, lAm1, LB-1, lBm1 );
  } // if(lA[iWork]>0)

  if(lB[iWork]>=2) {
    lBm1[iWork]--;
    tmpVal += (lBm1[iWork]+1) * 0.5 * pripair.one_over_gamma 
      * vRRiPPTab( pripair, shell1, shell2, LA, lA, LB-2, lBm1 );
    tmpVal -= (lBm1[iWork]+1) * 2.0 * Xi * 0.5 / shell2.alpha[pripair.p2] 
      * hRRiPPSab( pripair, shell1, shell2, LA, lA, LB-2, lBm1 );
  } // if(lB[iWork]>=2)

  tmpVal += 2.0 * Xi * hRRiPPSab( pripair, shell1, shell2, LA, lA, LB, lB );

  return tmpVal;
}  // RealGTOIntEngine::vRRiPPTab

  
  //---------------------------------------------------------------------------------//
  //angular momentum vertical recursion                                              //
  // if LB == 0,then [a|L|0]=(Pi-Ai)[a-1i|L|0]]+halfInvZeta*Ni(a-1i)[a-2i|L|0]         //
  //                  +zeta_b/zeta*{1i cross(B-C)}_mu*[a-1i|b]                       //
  // if LB>0,then  [a|L|b]=(Pi-Bi)[a|L|b-1i]                                         //
  //                  +halfInvZeta*Ni(b-1i)[a|L|b-2i]+halfInvZeta*Ni(a)[a-1i|L|b-1i] //
  //                  -zeta_a/zeta*{1i cross(A-C)}_mu*[a|b-1i]                       //
  //                  -halfInvZeta*Sum_k=x,y,z N_k(a){1i cross 1k}_mu*[a-1k|b-1i]    //
  //---------------------------------------------------------------------------------//
  /**
   *  \brief Perform the vertical recurrence relation for the uncontracted 
   *  angular momentum integral
   *
   *   if LB == 0,then [a|L|0]=(Pi-Ai)[a-1i|L|0]]+halfInvZeta*Ni(a-1i)[a-2i|L|0]      
   *                    +zeta_b/zeta*{1i cross(B-C)}_mu*[a-1i|b]                      
   *   if LB>0,then  [a|L|b]=(Pi-Bi)[a|L|b-1i]                                        
   *                    +halfInvZeta*Ni(b-1i)[a|L|b-2i]+halfInvZeta*Ni(a)[a-1i|L|b-1i]
   *                    -zeta_a/zeta*{1i cross(A-C)}_mu*[a|b-1i]                      
   *                    -halfInvZeta*Sum_k=x,y,z N_k(a){1i cross 1k}_mu*[a-1k|b-1i]   
   *
   *  where a,b are the angular momentum, A,B are the nuclear coordinates.
   *  halInvZeta = 1/2 * 1/Zeta
   *
   *  \param [in] pripair primitive Shell pair data for shell1, shell2
   *  \param [in] shell1  Bra shell
   *  \param [in] shell2  Ket shell
   *  \param [in] OneixAC 1_i cross AC vector is a 3 by 3 tensor
   *  \param [in] OneixBC 1_i cross BC vector is a 3 by 3 tensor
   *  \param [in] LA      total Bra angular momentum
   *  \param [in] lA      Bra angular momentum vector (lAx,lAy,lAz)
   *  \param [in] LB      total Ket angular momentum
   *  \param [in] lB      Ket angular momentum vector (lBx,lBy,lBz)
   *  \param [in] mu      index of the angular momentum integral component being calculated
   *
   *  \returns the mu component of an uncontracted angular momentum integral
   *
   */ 
   
  double RealGTOIntEngine::Labmu(
    libint2::ShellPair::PrimPairData &pripair, libint2::Shell &shell1, libint2::Shell &shell2, 
    double *OneixAC, double *OneixBC, int LA, int *lA, int LB, int *lB, int mu){

    double tmpVal=0.0;
  
    if ((LA+LB) == 0){  // [s|L_mu|s]
      double ACxBCmu;
      if (mu == 0)
        ACxBCmu     = shell1.O[1]*shell2.O[2]-shell1.O[2]*shell2.O[1];  
      else if (mu == 1)
        ACxBCmu     = shell1.O[2]*shell2.O[0]-shell1.O[0]*shell2.O[2];
      else if (mu == 2)
        ACxBCmu     = shell1.O[0]*shell2.O[1]-shell1.O[1]*shell2.O[0];   
   
      tmpVal = pow(sqrt(M_PI),3) * sqrt(pripair.one_over_gamma)*pripair.K ;
                                     //overlap ss integral 
      double Xi = shell1.alpha[pripair.p1] * shell2.alpha[pripair.p2]
                                    *pripair.one_over_gamma; 
      tmpVal = 2.0*Xi*ACxBCmu *tmpVal;
      return tmpVal;
    };
  
    int lAm1[3],lBm1[3];
    int iWork;
   
    if( ( LB  ==  0 ) and ( LA > 0 )) {
      for( iWork = 0; iWork < 3; iWork++) lAm1[iWork] = lA[iWork];
      if( lA[0] > 0 ) iWork = 0;
      else if( lA[1] > 0 ) iWork=1;
      else if( lA[2] > 0 ) iWork=2;
      lAm1[iWork]--;
  
      tmpVal  = (pripair.P[iWork]-shell1.O[iWork]) * 
        Labmu(pripair,shell1,shell2,OneixAC,OneixBC,LA-1,lAm1,LB,lB,mu);
      tmpVal += shell2.alpha[pripair.p2] * pripair.one_over_gamma * OneixBC[iWork*3+mu] *
        hRRiPPSab(pripair,shell1,shell2,LA-1,lAm1,LB,lB); 
  
  
      if(lAm1[iWork] >=1){
        lAm1[iWork]--;
        tmpVal += (lA[iWork] - 1) *0.5* pripair.one_over_gamma * 
          Labmu(pripair,shell1,shell2,OneixAC,OneixBC,LA-2,lAm1,LB,lB,mu);
      };
      return tmpVal;
  
    } else if( ( LB > 0 ) and ( LA  ==  0 ) ) {
  //    std::cout<<"here we are fucked"<<std::endl; 
      for( iWork = 0; iWork < 3; iWork++) lBm1[iWork] = lB[iWork];
      if( lB[0] > 0 ) iWork = 0;
      else if( lB[1] > 0 ) iWork=1;
      else if( lB[2] > 0 ) iWork=2;
      lBm1[iWork]--;
      
      tmpVal  = (pripair.P[iWork]-shell2.O[iWork]) *
        Labmu(pripair,shell1,shell2,OneixAC,OneixBC,LA,lA,LB-1,lBm1,mu);
      tmpVal -= shell1.alpha[pripair.p1] * pripair.one_over_gamma * OneixAC[3*iWork+mu] *
        hRRiPPSab(pripair,shell1,shell2,LA,lA,LB-1,lBm1);
  
      if(lBm1[iWork] >=1){
        lBm1[iWork]--;
        tmpVal += (lB[iWork] - 1) * 0.5*pripair.one_over_gamma *
          Labmu(pripair,shell1,shell2,OneixAC,OneixBC,LA,lA,LB-2,lBm1,mu);
      };
      return tmpVal;
    }
      
  //  else if( LB>LA ) cout <<"fucked again"<<endl;
    else if( ( LB > 0 ) and ( LA > 0 ) ){
  //   cout<<"here we are"<<endl;
      int lBm1[3],lAm1k[3];
      for(iWork = 0;iWork < 3;iWork++) {
        lAm1[iWork]=lA[iWork];
        lBm1[iWork]=lB[iWork];
        lAm1k[iWork]=lA[iWork];
      };
      if( lB[0] > 0 ) iWork = 0;
      else if( lB[1] > 0 ) iWork=1;
      else if( lB[2] > 0 ) iWork=2;
  
      lAm1[iWork]--;
      lBm1[iWork]--;
      tmpVal  = (pripair.P[iWork]-shell2.O[iWork]) * 
        Labmu(pripair,shell1,shell2,OneixAC,OneixBC,LA,lA,LB-1,lBm1,mu);
  
      if( lAm1[iWork]  >=  0 )
       tmpVal += 0.5*pripair.one_over_gamma * lA[iWork] * 
         Labmu(pripair,shell1,shell2,OneixAC,OneixBC,LA-1,lAm1,LB-1,lBm1,mu);
  
      tmpVal -= shell1.alpha[pripair.p1] * pripair.one_over_gamma * OneixAC[3*iWork+mu] *
        hRRiPPSab(pripair,shell1,shell2,LA,lA,LB-1,lBm1);
  
  /*  DBWY
      for(k = 0; k < 3; k++) {
        if(lA[k] > 0){
         lAm1k[k]--;
         tmpVal -= ijSP->halfInvZeta[iPP]*lA[k]*(*OneiOnek)(iWork,k,mu)*RealGTOIntEngine::hRRiPPSab(ijSP,LA-1,lAm1k,LB-1,lBm1,iPP);
         lAm1k[k]++;
        };
      };
  */ 
  
       if( iWork  ==  0 ) {
         if( mu  ==  1 and lA[2] > 0 ) {
           lAm1k[2]--;
           tmpVal += 0.5*pripair.one_over_gamma * lA[2] * 
             hRRiPPSab(pripair,shell1,shell2,LA-1,lAm1k,LB-1,lBm1);
           lAm1k[2]++;
         } else if( mu  ==  2 and lA[1] > 0 ) {
           lAm1k[1]--;
           tmpVal -= 0.5*pripair.one_over_gamma * lA[1] * 
             hRRiPPSab(pripair,shell1,shell2,LA-1,lAm1k,LB-1,lBm1);
           lAm1k[1]++;
         }
       } else if( iWork  ==  1 ) {
         if( mu  ==  2 and lA[0] > 0 ) {
           lAm1k[0]--;
           tmpVal += 0.5*pripair.one_over_gamma * lA[0] * 
             hRRiPPSab(pripair,shell1,shell2,LA-1,lAm1k,LB-1,lBm1);
           lAm1k[0]++;
         } else if( mu  ==  0 and lA[2] > 0 ) {
           lAm1k[2]--;
           tmpVal -= 0.5*pripair.one_over_gamma * lA[2] * 
             hRRiPPSab(pripair,shell1,shell2,LA-1,lAm1k,LB-1,lBm1);
           lAm1k[2]++;
         }
       } else {
         if( mu  ==  0 and lA[1] > 0 ) {
           lAm1k[1]--;
           tmpVal += 0.5*pripair.one_over_gamma * lA[1] * 
             hRRiPPSab(pripair,shell1,shell2,LA-1,lAm1k,LB-1,lBm1);
           lAm1k[1]++;
         } else if( mu  ==  1 and lA[0] > 0 ) {
           lAm1k[0]--;
           tmpVal -= 0.5*pripair.one_over_gamma * lA[0] * 
             hRRiPPSab(pripair,shell1,shell2,LA-1,lAm1k,LB-1,lBm1);
           lAm1k[0]++;
         }
       }
      
      if (lBm1[iWork] > 0){
       lBm1[iWork]--;
       tmpVal += (lBm1[iWork] + 1) * 0.5*pripair.one_over_gamma * 
         Labmu(pripair,shell1,shell2,OneixAC,OneixBC,LA,lA,LB-2,lBm1,mu);
      };
    };
  
    return tmpVal;
  };
  
   
  //-----------------------------------------------------------------//
  //momentum integral                                                //
  //[a|\nabla_\beta|b]=-2*zeta_b)([a+1_\beta||b]+(A-B)_\beta[a||b])   //
  //                    +N_\beta(b)[a||b-1\beta]                     //
  //-----------------------------------------------------------------//
  /**
   *  \brief Decompose the momentum integral into several contracted overlap integral
   *         and uncontracted integrals
   *
   *  [a|\nabla_\mu|b]=-2*zeta_b)([a+1_\mu||b]+(A-B)_\mu[a||b])
   *                     +N_\mu(b)[a||b-1\mu]
   *
   *  where a,b are the angular momentum, A,B are the nuclear coordinates.
   *
   *  \param [in] pair    Shell pair data for shell1, shell2
   *  \param [in] shell1  Bra shell
   *  \param [in] shell2  Ket shell
   *  \param [in] LA      total Bra angular momentum
   *  \param [in] lA      Bra angular momentum vector (lAx,lAy,lAz)
   *  \param [in] LB      total Ket angular momentum
   *  \param [in] lB      Ket angular momentum vector (lBx,lBy,lBz)
   *  \param [in] mu      the component of \nabla operator (x,y,z)
   *
   *  \returns the mu component of contracted momentum integral
   *
   */ 
   
  double RealGTOIntEngine::Momentummu( 
    libint2::ShellPair &pair, libint2::Shell &shell1, libint2::Shell &shell2, int LA, int *lA, int LB, int *lB, int mu) {
  
    int lAp1[3],lBm1[3];
    double tmpVal = 0.0;
  
  //  if(LB == 0) {
      for (auto k = 0 ; k < 3 ; k++ ){
        lAp1[k] = lA[k];
        lBm1[k] = lB[k];
      }
      lAp1[mu] = lAp1[mu]+1;
      
      for ( auto pripair : pair.primpairs ){
  //      tmpVal2 += RealGTOIntEngine::hRRiPPSab(ijSP,LA+1,lAp1,LB,lB,iPP);
  //      tmpVal2 += RealGTOIntEngine::hRRiPPSab(ijSP,LA,lA,LB,lB,iPP)*ijSP->AB[mu]; 
  //      tmpVal += ijSP->zetab[iPP] * tmpVal2 ;
        auto coeffs = shell1.contr[0].coeff[pripair.p1]* shell2.contr[0].coeff[pripair.p2];
        // -2*zeta_b*(a+1||b)
        tmpVal -= coeffs*2.0*shell2.alpha[pripair.p2]
                  *hRRiPPSab(pripair,shell1,shell2,LA+1,lAp1,LB,lB);
        // -2*(A-B)*(a||b)   
        tmpVal -= coeffs*2.0*shell2.alpha[pripair.p2] 
                  *pair.AB[mu]*hRRiPPSab(pripair,shell1,shell2,LA,lA,LB,lB);
  //      tmpVal2 = 0;
      }
      if (lB[mu] == 0) {
        return tmpVal;
      } else if ( lB[mu] > 0 ){
        lBm1[mu] = lBm1[mu]-1;
        tmpVal +=  lB[mu]*hRRSab( pair,shell1,shell2,LA,lA,LB-1,lBm1 );
        return tmpVal;  
      }

      return 0.0;
  }
  
  //--------------------------------------------------------------------//
  //Electric Dipole (E1) integral                                       //
  //[a|r_\alpha|b] = [a+1_\alpha||b] + A_\alpha [a||b]                  //
  //--------------------------------------------------------------------//
  /**
   *  \brief Decompose electric dipole integral into several contracted overlap integral
   *
   *  [a|r_\mu|b] = [a+1_\mu||b] + A_\mu [a||b]
   *
   *  where a,b are the angular momentum, A,B are the nuclear coordinates.
   *
   *  \param [in] pair    Shell pair data for shell1, shell2
   *  \param [in] shell1  Bra shell
   *  \param [in] shell2  Ket shell
   *  \param [in] LA      total Bra angular momentum
   *  \param [in] lA      Bra angular momentum vector (lAx,lAy,lAz)
   *  \param [in] LB      total Ket angular momentum
   *  \param [in] lB      Ket angular momentum vector (lBx,lBy,lBz)
   *  \param [in] mu      the component of r operator (x,y,z)
   *
   *  \returns the mu component of a contracted electric dipole integral
   *
   */ 
   
  double RealGTOIntEngine::DipoleE1( 
    libint2::ShellPair &pair, libint2::Shell &shell1, libint2::Shell &shell2, int LA, int *lA, int LB, int *lB, int mu) {
     
    int lAp1[3];
    double tmpVal = 0.0;
    
    for (auto k = 0 ; k < 3 ; k++ ){ 
      lAp1[k] = lA[k];  
    }
    lAp1[mu] = lAp1[mu] + 1;
    
    tmpVal = hRRSab( pair,shell1,shell2, LA+1, lAp1, LB, lB ) + 
             +shell1.O[mu]*hRRSab( pair,shell1,shell2, LA, lA, LB, lB );
    
    return tmpVal;
  }
  
  //----------------------------------------------------------------------//
  //Electric Quadrupole (E2) integral                                     //
  //[a|\nabla_\alpha r_\beta+r_\alpha\nabla_\beta|b]                      //
  // = \delta_\alpha\eta [a||b]+[a+1_\beta|\nabla_\alpha|b]               //
  //   +A_\beta[a|\nabla_\alpha|b]+[a+1_\alpha|\nabla_\beta|b]            //
  //   +A_\beta[a|\nabla_\beta|b]                                         //
  //----------------------------------------------------------------------//
  /**
   *  \brief Decompose electric quadrupole integral into several contracted overlap integral
   *
   *   [a|\nabla_\alpha r_\beta+r_\alpha\nabla_\beta|b]           
   *    = \delta_\alpha\eta [a||b]+[a+1_\beta|\nabla_\alpha|b]   
   *      +A_\beta[a|\nabla_\alpha|b]+[a+1_\alpha|\nabla_\beta|b]
   *      +A_\beta[a|\nabla_\beta|b]                             
   *
   *  where a,b are the angular momentum, A,B are the nuclear coordinates.
   *
   *  \param [in] pair    Shell pair data for shell1, shell2
   *  \param [in] shell1  Bra shell
   *  \param [in] shell2  Ket shell
   *  \param [in] LA      total Bra angular momentum
   *  \param [in] lA      Bra angular momentum vector (lAx,lAy,lAz)
   *  \param [in] LB      total Ket angular momentum
   *  \param [in] lB      Ket angular momentum vector (lBx,lBy,lBz)
   *  \param [in] alpha   the component of nabla or r operator (x,y,z)
   *  \param [in] beta    the component of nabla or r operator (x,y,z)
   *
   *  \returns the alpha beta component of a contracted electric quadrupole integral
   *
   */ 
   
  double RealGTOIntEngine::QuadrupoleE2_vel( libint2::ShellPair &pair, libint2::Shell &shell1, 
    libint2::Shell &shell2, int LA, int *lA, int LB, int *lB, int alpha, int beta ) {
  
    int lAp1[3],lBp1[3];
    double tmpVal = 0.0;  
  
    for ( auto k = 0 ; k < 3 ; k++ ) {
      lAp1[k] = lA[k];
      lBp1[k] = lB[k];
    }
    lAp1[alpha] = lA[alpha]+1;
    lBp1[beta]  = lB[beta] +1;
    
    tmpVal += Momentummu( pair,shell1,shell2, LA, lA, LB+1, lBp1, alpha );
    tmpVal += shell2.O[beta]* Momentummu( pair,shell1,shell2, LA, lA, LB, lB, alpha );
    tmpVal += Momentummu( pair,shell1,shell2, LA+1, lAp1, LB, lB, beta );
    tmpVal += shell1.O[alpha] * Momentummu( pair,shell1,shell2, LA, lA, LB, lB, beta );
    
    return tmpVal;
  
  
  }  
  
  //Angular momentum integral contracted (magnetic dipole)     
  /**
   *  \brief Calculate contracted angular momenum integral(magnetic dipole integral) 
   *
   *   [a|\nabla_\alpha r_\beta+r_\alpha\nabla_\beta|b]           
   *    = \delta_\alpha\eta [a||b]+[a+1_\beta|\nabla_\alpha|b]   
   *      +A_\beta[a|\nabla_\alpha|b]+[a+1_\alpha|\nabla_\beta|b]
   *      +A_\beta[a|\nabla_\beta|b]                             
   *
   *  where a,b are the angular momentum, A,B are the nuclear coordinates.
   *
   *  \param [in] pair    Shell pair data for shell1, shell2
   *  \param [in] shell1  Bra shell
   *  \param [in] shell2  Ket shell
   *  \param [in] LA      total Bra angular momentum
   *  \param [in] lA      Bra angular momentum vector (lAx,lAy,lAz)
   *  \param [in] LB      total Ket angular momentum
   *  \param [in] lB      Ket angular momentum vector (lBx,lBy,lBz)
   *  \param [in] mu   the component of angular momentum operator (x,y,z)
   *
   *  \returns the mu component of a contracted magnetic dipole integral
   *
   */ 

    
  double RealGTOIntEngine::MDipoleM1( libint2::ShellPair &pair, libint2::Shell &shell1, 
    libint2::Shell &shell2, int LA, int *lA, int LB, int *lB, int mu){
   
    double OneixBC[9];
    OneixBC[0*3+0] = 0.0;
    OneixBC[0*3+1] = -shell2.O[2] ; // -Bz
    OneixBC[0*3+2] = shell2.O[1];   //  By
    OneixBC[1*3+0] = shell2.O[2] ;  //  Bz
    OneixBC[1*3+1] = 0.0;
    OneixBC[1*3+2] = -shell2.O[0];  // -Bx
    OneixBC[2*3+0] = -shell2.O[1];  // -By
    OneixBC[2*3+1] = shell2.O[0];   //  Bx
    OneixBC[2*3+2] = 0.0;
   
    double OneixAC[9];
    OneixAC[0*3+0] = 0.0;
    OneixAC[0*3+1] = -shell1.O[2];  // -Az
    OneixAC[0*3+2] = shell1.O[1];   //  Ay
    OneixAC[1*3+0] = shell1.O[2];   //  Az   
    OneixAC[1*3+1] = 0.0;         
    OneixAC[1*3+2] = -shell1.O[0];  // -Ax       
    OneixAC[2*3+0] = -shell1.O[1];  // -Ay    
    OneixAC[2*3+1] = shell1.O[0];   //  Ax   
    OneixAC[2*3+2] = 0.0;         
  
    double L_mu = 0.0;
  
    int kk = 0;
    for ( auto pripair : pair.primpairs ){
  
      L_mu += shell1.contr[0].coeff[pripair.p1]* shell2.contr[0].coeff[pripair.p2]* 
                   Labmu(pripair,shell1,shell2,OneixAC,OneixBC,
                   LA,lA,LB,lB,mu) ;
    }
    return L_mu;
    
  }
  
  
  //------------------------------------------------------------------//
  //Magnetic Quadrupole                                               //
  //[a|MQ_\alpha\beta|b]=[a+1_\alpha|L_\beta|b]+A_\alpha[a|L_\beta|b] //
  //                    +[a|L_\alpha|b+1_\beta]+B_\beta[a|L_\alpha|b] //
  //------------------------------------------------------------------//
  /**
   *  \brief Calculate contracted angular momenum integral(magnetic dipole integral) 
   *
   *  [a|MQ_\alpha\beta|b]=[a+1_\alpha|L_\beta|b]+A_\alpha[a|L_\beta|b]
   *                      +[a|L_\alpha|b+1_\beta]+B_\beta[a|L_\alpha|b]
   *
   *  where a,b are the angular momentum, A,B are the nuclear coordinates. L is angular
   *  momentum operator
   *
   *  \param [in] pair    Shell pair data for shell1, shell2
   *  \param [in] shell1  Bra shell
   *  \param [in] shell2  Ket shell
   *  \param [in] LA      total Bra angular momentum
   *  \param [in] lA      Bra angular momentum vector (lAx,lAy,lAz)
   *  \param [in] LB      total Ket angular momentum
   *  \param [in] lB      Ket angular momentum vector (lBx,lBy,lBz)
   *  \param [in] alpha   the component of L or r operator (x,y,z)
   *  \param [in] beta    the component of L or r operator (x,y,z)
   *
   *  \returns the alpha beta component of a contracted magnetic quadrupole integral
   *
   */ 

   
  double RealGTOIntEngine::QuadrupoleM2_vel( libint2::ShellPair &pair, libint2::Shell &shell1, 
    libint2::Shell &shell2, int LA, int *lA, int LB, int *lB, int alpha, int beta ) {
  
    int lAp1[3],lBp1[3];
    double tmpVal = 0.0;
  
    for ( auto k = 0 ; k < 3 ; k++ ) {
      lAp1[k] = lA[k];
      lBp1[k] = lB[k];
    }
    lAp1[beta] = lA[beta] +1;
    lBp1[beta] = lB[beta] +1;
  
    tmpVal += MDipoleM1( pair,shell1,shell2, LA+1, lAp1, LB, lB, alpha );
  
  
    tmpVal += MDipoleM1( pair,shell1,shell2, LA, lA, LB+1, lBp1, alpha );
  
  
    tmpVal += (shell1.O[beta]+shell2.O[beta]) * MDipoleM1( pair,shell1,shell2, LA, lA, LB, lB, alpha );
    
    return tmpVal;
  
  }
  
  
  //--------------------------------------------------------------//
  //Electric Octupole integral                                    //
  //[a|Oct(alpha,beta,gamma)|b]=[a+1_\alpha|Qua(\beta,\gamma)|b]  //
  //  +A_\alpha [a|Qua(\beta,\gamma)|b]                           //
  //  +[a|\nabla_\alpha|b+1_\beta+1_\gamma]                       //
  //  +B_\beta [a|\nabla_\alpha|b+1_\gamma]                       //
  //  +B_\gamma [a|\nabla_\alpha|b+1_\beta]                       //
  //  +B_\beta B_\gamma [a|\nabla_\alpha|b]                       //
  //--------------------------------------------------------------//
  /**
   *  \brief Decompose electric octupole integral into several contracted electric
   *   quadrupole, angular momemtum and overlap integral
   *
   *  [a|Oct(alpha,beta,gamma)|b]=[a+1_\alpha|Qua(\beta,\gamma)|b] 
   *    +A_\alpha [a|Qua(\beta,\gamma)|b]                         
   *    +[a|\nabla_\alpha|b+1_\beta+1_\gamma]                     
   *    +B_\beta [a|\nabla_\alpha|b+1_\gamma]                     
   *    +B_\gamma [a|\nabla_\alpha|b+1_\beta]                     
   *    +B_\beta B_\gamma [a|\nabla_\alpha|b]                     
   *     
   *  where a,b are the angular momentum, A,B are the nuclear coordinates.
   *
   *  \param [in] pair    Shell pair data for shell1, shell2
   *  \param [in] shell1  Bra shell
   *  \param [in] shell2  Ket shell
   *  \param [in] LA      total Bra angular momentum
   *  \param [in] lA      Bra angular momentum vector (lAx,lAy,lAz)
   *  \param [in] LB      total Ket angular momentum
   *  \param [in] lB      Ket angular momentum vector (lBx,lBy,lBz)
   *  \param [in] alpha   the component of electric quadrupole or nabla or r operator (x,y,z)
   *  \param [in] beta    the component of electric quadrupole or nabla or r operator (x,y,z)
   *  \param [in] gamma      the component of electric quadrupole or nabla or r operator (x,y,z)
   *
   *  \returns the alpha beta gamma component of a contracted electric octupole integral
   *
   */ 

   
  double RealGTOIntEngine::OctupoleE3_vel( libint2::ShellPair &pair,libint2::Shell &shell1,
    libint2::Shell &shell2, int LA, int *lA, int LB, int *lB, int alpha, int beta, int gamma ) {
   
    int lAp1[3],lBp1_beta[3],lBp1_gamma[3],lBp2[3];
    double tmpVal = 0.0;  
  
    for ( auto k = 0 ; k < 3 ; k++ ) {
      lAp1[k] = lA[k];
      lBp1_beta[k] = lB[k];
      lBp1_gamma[k]= lB[k];
      lBp2[k] = lB[k];
    }
    lAp1[alpha] = lA[alpha]+1;
    lBp1_beta[beta]  = lB[beta] +1;
    lBp1_gamma[gamma] = lB[gamma]+1;
    lBp2[gamma] = lBp2[gamma]+1;
    lBp2[beta] = lBp2[beta]+1;
     
    tmpVal += QuadrupoleE2_vel( pair,shell1,shell2, LA+1, lAp1, LB, lB, beta, gamma );
    tmpVal += shell1.O[alpha] 
              * QuadrupoleE2_vel( pair,shell1,shell2, LA, lA, LB, lB, beta, gamma );
    tmpVal += Momentummu( pair,shell1,shell2, LA, lA, LB+2, lBp2, alpha );
    tmpVal += shell2.O[beta] * Momentummu( pair,shell1,shell2, LA, lA, LB+1, lBp1_gamma, alpha );
    tmpVal += shell2.O[gamma] * Momentummu( pair,shell1,shell2, LA, lA, LB+1, lBp1_beta, alpha );
    tmpVal += shell2.O[beta] *shell2.O[gamma]*Momentummu( pair,shell1,shell2, LA, lA, LB, lB, alpha ) ; 
  
    return tmpVal;
  
    }
  
 
  //---------------------------------------------------------//
  // potential integral horizontal recursion                 //
  //  (a|0_c|b) = (a+1|0_c|b-1) + (A-B)(a|0_c|b-1)           //
  //   LA  >=  LB                                              //
  //   horizontal recursion doesn't increase (m). and it's only used once, so m=0 here. //
  //---------------------------------------------------------//

  /**
   *  \brief Perform the horizontal recurrence relation for the contracted 
   *  nuclear potential integral
   *
   *  (a|0_c|b) = (a+1|0_c|b-1) + (A-B)(a|0_c|b-1)
   *
   *  where a,b are the angular momentum, A,B are the nuclear coordinates.
   *
   *  \param [in] nucShell nuclear shell, give the exponents of gaussian function of nuclei
   *  \param [in] pair    Shell pair data for shell1, shell2
   *  \param [in] shell1  Bra shell
   *  \param [in] shell2  Ket shell
   *  \param [in] LA      total Bra angular momentum
   *  \param [in] lA      Bra angular momentum vector (lAx,lAy,lAz)
   *  \param [in] LB      total Ket angular momentum
   *  \param [in] lB      Ket angular momentum vector (lBx,lBy,lBz)
   *
   *  \returns a contracted nuclear potential integral
   *
   */ 
   
  double RealGTOIntEngine::hRRVab(const std::vector<libint2::Shell> &nucShell, 
    libint2::ShellPair &pair, libint2::Shell &shell1,
    libint2::Shell &shell2, int LA, int *lA, int LB, int *lB, const Molecule& molecule){

    int iWork,iAtom;
    double tmpVal=0.0;
    double PC[3];
    double rho;
    double squarePC=0.0;
    bool useFiniteWidthNuclei = nucShell.size() > 0;
  
    if(LB == 0) {
      // (LA|s)
      for( auto pripair : pair.primpairs ) {
      iAtom = 0;
      for( auto atom : molecule.atoms ) {
  //      std::cerr<<atom.atomicNumber<<std::endl;
  //      std::cerr<<"iAtom = "<<iAtom<<std::endl;
        squarePC=0.0;     
  
        for( int m=0 ; m<3 ; m++ ) {
          PC[m] = pripair.P[m] - atom.coord[m];
          squarePC += PC[m]*PC[m];
        }
        auto lTotal = shell1.contr[0].l + shell2.contr[0].l; 
        double *tmpFmT = new double[lTotal+1];
        if ( useFiniteWidthNuclei ) {
          rho = (1/pripair.one_over_gamma)* nucShell[iAtom].alpha[0]
                       /(1/pripair.one_over_gamma + nucShell[iAtom].alpha[0]);

          // or use double rho = nucShell[iAtom].alpha[0]/
          // (1.0+nucShell[iAtom].alpha[0]*pripair.one_over_gamma)

          RealGTOIntEngine::computeFmTTaylor(tmpFmT,rho*squarePC,lTotal,0);
        }
        else if ( !useFiniteWidthNuclei ) {
          RealGTOIntEngine::computeFmTTaylor(tmpFmT,(shell1.alpha[pripair.p1]+shell2.alpha[pripair.p2])
                     *squarePC,lTotal,0);
        }
        if (LA == 0) {
          auto norm = shell1.contr[0].coeff[pripair.p1]* 
                      shell2.contr[0].coeff[pripair.p2];
          auto ssS = pow(sqrt(M_PI),3) * sqrt(pripair.one_over_gamma)*pripair.K ;

          if ( !useFiniteWidthNuclei ) {
            auto ssV = 2.0*sqrt(1.0/(pripair.one_over_gamma*M_PI))*norm*ssS;
  
            tmpVal += atom.nucCharge * ssV * tmpFmT[0];
//std::cout<<"tmpFmT "<<tmpFmT[0]<<std::endl;
  // std::cerr<<"actual"<<std::endl;
          }
          else if ( useFiniteWidthNuclei ) {
  //          tmpVal += (static_cast<double>(mc->atomZ[iAtom]))*math.two*sqrt(rho/math.pi)*ijSP->ss[iPP]*tmpFmT[0];
            auto ssV = 2.0*sqrt(rho/M_PI)*norm*ssS;
 
            tmpVal += atom.nucCharge * ssV * tmpFmT[0];

//std::cout<<"tmpFmT "<<tmpFmT[0]<<std::endl;
//  std::cerr<<"no finite nuclei"<<std::endl;
          }
        } // LA == 0
        else {   // LA > 0, go into vertical recursion
  
          auto norm = shell1.contr[0].coeff[pripair.p1]* 
                      shell2.contr[0].coeff[pripair.p2];
          auto ssS = pow(sqrt(M_PI),3) * sqrt(pripair.one_over_gamma)*pripair.K ;
 
          if ( !useFiniteWidthNuclei ) {
            auto ssV = 2.0*sqrt(1.0/(pripair.one_over_gamma*M_PI))*norm*ssS;  
   
            tmpVal += atom.nucCharge * ssV * vRRVa0(nucShell,pripair,shell1,
                                                      tmpFmT,PC,0,LA,lA,iAtom);
  //std::cerr<<"actual"<<std::endl;
  
          } else if ( useFiniteWidthNuclei ) {
  //          tmpVal += mc->atomZ[iAtom]*math.two*sqrt(rho/math.pi)*ijSP->ss[iPP]*RealGTOIntEngine::vRRVa0(ijSP,tmpFmT,PC,0,LA,lA,iPP,iAtom);

            auto ssV = 2.0*sqrt(rho/M_PI)*norm*ssS;
            tmpVal += atom.nucCharge * ssV * vRRVa0(nucShell,pripair,shell1,
                                                       tmpFmT,PC,0,LA,lA,iAtom);

//  std::cerr<<"no finite nuclei"<<std::endl;
          }
        } // else
        delete[] tmpFmT;
        iAtom++;
      } // atom
      }  // pripair
    } // if LB  ==  0
  
     else {   // LB>0
      int lAp1[3],lBm1[3];
      for( int m=0 ; m<3 ; m++ ) {
        lAp1[m]=lA[m];
        lBm1[m]=lB[m];
      };
      if (lB[0] > 0) iWork = 0;
      else if (lB[1] > 0) iWork=1;
      else if (lB[2] > 0) iWork=2;
  
      lAp1[iWork]++;
      lBm1[iWork]--;
      tmpVal = hRRVab(nucShell,pair ,shell1,shell2, LA+1,lAp1, LB-1,lBm1,molecule);
      tmpVal+= pair.AB[iWork]*hRRVab(nucShell,pair,shell1, shell2, LA,lA,LB-1,lBm1,molecule);
    }
  
    return tmpVal;
  }
  
  //----------------------------------------------------------------------------//
  // potential integral vertical recursion                                      //
  //  (a|0_c|0)^(m) = (P-A)(a-1|0_c|0)^(m) - (P-C)(a-1|0_c|0)^(m+1)             //
  //                + halfInvZeta*N_(a-1)*[(a-2|0_c|0)^(m)-(a-2|0_c|0)^(m+1)]   //
  //  since LB == 0, we only decrease a                                         //
  //----------------------------------------------------------------------------//

  /**
   *  \brief Perform the vertical recurrence relation for the uncontracted 
   *    nuclear potential integral 
   *
   *   (a|0_c|0)^(m) = (P-A)(a-1|0_c|0)^(m) - (P-C)(a-1|0_c|0)^(m+1) 
   *                 + halfInvZeta*N_(a-1)*[(a-2|0_c|0)^(m)-(a-2|0_c|0)^(m+1)]
   *
   *  where a is angular momentum, Zeta=zeta_a+zeta_b, A is bra nuclear coordinate. 
   *  P = (zeta_a*A+zeta_b*B)/Zeta
   *
   *  \param [in] nucShell nuclear shell, give the exponents of gaussian function of nuclei
   *  \param [in] pripair Primitive Shell pair data for shell1, shell2
   *  \param [in] shell1  Bra shell
   *  \param [in] FmT     table of Boys function with different m
   *  \param [in] PC      vector P-C
   *  \param [in] m       the order of auxiliary function
   *  \param [in] LA      total Bra angular momentum
   *  \param [in] lA      Bra angular momentum vector (lAx,lAy,lAz)
   *  \param [in] iAtom   the index of the atom of the nuclei
   *
   *  \returns an uncontracted nuclear potential integral
   *
   */ 

   
  double RealGTOIntEngine::vRRVa0( const std::vector<libint2::Shell> &nucShell,
    libint2::ShellPair::PrimPairData &pripair, libint2::Shell &shell1,
    double *FmT, double *PC, int m, int LA, int *lA, int iAtom){
  
    double rhoovzeta;
    bool useFiniteWidthNuclei = nucShell.size() > 0;

    if(LA == 0) {
      if ( useFiniteWidthNuclei ) {
//  std::cerr<<"out here"<<std::endl;
        rhoovzeta =  nucShell[iAtom].alpha[0] / ( nucShell[iAtom].alpha[0]
                     +(1/pripair.one_over_gamma) );

//std::cout<<"tmpFmT "<<FmT[m]<<std::endl;
        return (pow(rhoovzeta,m)*FmT[m]);

      } else if ( !useFiniteWidthNuclei ) {       
//std::cout<<"tmpFmT "<<FmT[m]<<std::endl;
        return FmT[m];  //Z*2*sqrt[zeta/pi]*[s|s]is given in hRRVab, in ssV.
      }
    }
  
    double tmpVal=0.0;
    int lAm1[3]; //means lA minus 1_i
    int iWork;
    for( iWork = 0 ; iWork < 3 ; iWork++ ) lAm1[iWork]=lA[iWork];
    if (lA[0] > 0) iWork = 0;
    else if (lA[1] > 0) iWork=1;
    else if (lA[2] > 0) iWork=2;
    lAm1[iWork]--;
  
  if (LA >= -1) {
  }
    tmpVal  = (pripair.P[iWork]-shell1.O[iWork])*
                vRRVa0(nucShell,pripair,shell1,FmT,PC,m,LA-1,lAm1,iAtom);
  
    tmpVal -= PC[iWork] * vRRVa0(nucShell,pripair,shell1,FmT,PC,m+1,LA-1,lAm1,iAtom);
  
    if( lAm1[iWork] >=1 ){
      lAm1[iWork]--;
      tmpVal += (lAm1[iWork]+1)*0.5*pripair.one_over_gamma *
                (vRRVa0(nucShell,pripair,shell1,FmT,PC,m,LA-2,lAm1,iAtom)
                -vRRVa0(nucShell,pripair,shell1,FmT,PC,m+1,LA-2,lAm1,iAtom));
    }
    return tmpVal; 
  }
  
  //----------------------------------------------------------------------------//
  // potential integral horizontal recursion iPP specific                       //
  //  [a|A(0)|b]^(m) = [a+1i|A(0)|b-1i]^(m) + (Ai-Bi)*[a|A(0)|b-1i]^(m)         //
  //  LA  >= LB                                                                   //
  //  m can be nonzero                                                          //
  //----------------------------------------------------------------------------//

  /**
   *  \brief Perform the horizontal recurrence relation for the uncontracted 
   *   nuclear potential integral
   *
   *  [a|A(0)|b]^(m) = [a+1i|A(0)|b-1i]^(m) + (Ai-Bi)*[a|A(0)|b-1i]^(m)
   *
   *  where a,b are the angular momentum, A,B are the nuclear coordinates.
   *
   *  \param [in] nucShell nuclear shell, give the exponents of gaussian function of nuclei
   *  \param [in] pripair Primitive Shell pair data for shell1, shell2
   *  \param [in] shell1  Bra shell
   *  \param [in] shell2  Ket shell
   *  \param [in] LA      total Bra angular momentum
   *  \param [in] lA      Bra angular momentum vector (lAx,lAy,lAz)
   *  \param [in] LB      total Ket angular momentum
   *  \param [in] lB      Ket angular momentum vector (lBx,lBy,lBz)
   *  \param [in] C       coordinate of nuclei
   *  \param [in] m       order of auxiliary function
   *  \param [in] iAtom   index of nuclei
   *
   *  \returns a uncontracted nuclear potential integral. 
   *     It doesn't include contraction coeffs.
   *
   */ 

   
  double RealGTOIntEngine::hRRiPPVab( const std::vector<libint2::Shell> &nucShell,
    libint2::ShellPair::PrimPairData &pripair, libint2::Shell &shell1, 
    libint2::Shell &shell2, int LA, int *lA, int LB, int *lB, double *C, int m, int iAtom, 
    const Molecule& molecule){
  
    double tmpVal=0.0;
    double PC[3];
    double squarePC=0.0;
    double rho,rhoovzeta;
    bool useFiniteWidthNuclei = nucShell.size() > 0 ;
    for(int k = 0 ; k < 3 ; k++ ) {
      PC[k] = pripair.P[k] - C[k];
      squarePC += PC[k]*PC[k];
    }
    int lTotal = shell1.contr[0].l + shell2.contr[0].l; 
  
  
    if(LB == 0) {
      // (LA|s)
        double *tmpFmT = new double[lTotal+m+1];
        if ( useFiniteWidthNuclei ) {
  //        rho = ijSP->Zeta[iPP]*RealGTOIntEngine::molecule_->nucShell(iAtom).alpha[0]/(ijSP->Zeta[iPP]+RealGTOIntEngine::molecule_->nucShell(iAtom).alpha[0]);
  //        RealGTOIntEngine::computeFmTTaylor(tmpFmT,rho*squarePC,ijSP->lTotal+m,0);
          rho = (1/pripair.one_over_gamma)* nucShell[iAtom].alpha[0]
                /(1/pripair.one_over_gamma + nucShell[iAtom].alpha[0]);
          computeFmTTaylor(tmpFmT,rho*squarePC,lTotal+m,0);
  
//  std::cerr<<"no finite nuclei"<<std::endl;
        }
        else if ( !useFiniteWidthNuclei ) {
          computeFmTTaylor(tmpFmT,(shell1.alpha[pripair.p1]+shell2.alpha[pripair.p2])
                     *squarePC,lTotal+m,0);
        }
        if(LA == 0) {
  //          auto norm = shell1.contr[0].coeff[pripair.p1]* 
  //                      shell2.contr[0].coeff[pripair.p2];

          auto ssS = pow(sqrt(M_PI),3) * sqrt(pripair.one_over_gamma)*pripair.K ;

          if ( !useFiniteWidthNuclei ) {
            auto ssV = 2.0*sqrt(1.0/(pripair.one_over_gamma*M_PI))*ssS;
  
            tmpVal = ssV * tmpFmT[m];
//std::cout<<"tmpFmT "<<tmpFmT[m]<<std::endl;
  //          tmpVal = ijSP->ssV[iPP]*tmpFmT[m];
          }
          else if ( useFiniteWidthNuclei ) {
            auto ssV = 2.0*sqrt(rho/M_PI)*ssS;            
            rhoovzeta = nucShell[iAtom].alpha[0] / ( nucShell[iAtom].alpha[0]+
                        1.0/pripair.one_over_gamma );
            
  //          double  rhoovzeta =  RealGTOIntEngine::molecule_->nucShell(iAtom).alpha[0]/(RealGTOIntEngine::molecule_->nucShell(iAtom).alpha[0]+ijSP->Zeta[iPP]);
  //          tmpVal = math.two*sqrt(rho/math.pi)*pow(rhoovzeta,m)*ijSP->ss[iPP]*tmpFmT[m];
            tmpVal = pow(rhoovzeta,m) * ssV * tmpFmT[m]; 
//std::cout<<"tmpFmT "<<tmpFmT[m]<<std::endl;
//  std::cerr<<"no finite nuclei"<<std::endl;
          }
        } // LA == 0
  
        else if (LA>0) {

          auto ssS = pow(sqrt(M_PI),3) * sqrt(pripair.one_over_gamma)*pripair.K ;

          if ( !useFiniteWidthNuclei ) {
  //          auto norm = shell1.contr[0].coeff[pripair.p1]* 
  //                      shell2.contr[0].coeff[pripair.p2];
            auto ssV = 2.0*sqrt(1.0/(pripair.one_over_gamma*M_PI))*ssS;
   
            tmpVal = ssV * vRRVa0(nucShell,pripair,shell1,tmpFmT,PC,m,LA,lA,iAtom);
  
  //          tmpVal = ijSP->ssV[iPP]*RealGTOIntEngine::vRRVa0(ijSP,tmpFmT,PC,m,LA,lA,iPP,iAtom);
          }
          else if ( useFiniteWidthNuclei ) {
            tmpVal = 2.0 * sqrt(rho/M_PI) * ssS * 
                     vRRVa0(nucShell,pripair,shell1,tmpFmT,PC,m,LA,lA,iAtom);

  //          tmpVal = math.two*sqrt(rho/math.pi)*ijSP->ss[iPP]*RealGTOIntEngine::vRRVa0(ijSP,tmpFmT,PC,m,LA,lA,iPP,iAtom);
//  std::cerr<<"no finite nuclei"<<std::endl; 
          } 
        } // else if (LA>0)
        delete[] tmpFmT;
    } // LB == 0
  
    else if ((LA == 0)&&(LB>0)){
  
        double *tmpFmT = new double[lTotal+m+1];
        auto ssS = pow(sqrt(M_PI),3) * sqrt(pripair.one_over_gamma)*pripair.K ;

        if ( useFiniteWidthNuclei ) {
  //        rho = ijSP->Zeta[iPP]*RealGTOIntEngine::molecule_->nucShell(iAtom).alpha[0]/(ijSP->Zeta[iPP]+RealGTOIntEngine::molecule_->nucShell(iAtom).alpha[0]);
  //        RealGTOIntEngine::computeFmTTaylor(tmpFmT,rho*squarePC,ijSP->lTotal+m,0);
  //        tmpVal =  math.two*sqrt(rho/math.pi)*ijSP->ss[iPP]*RealGTOIntEngine::vRRV0b(ijSP,tmpFmT,PC,m,LB,lB,iPP,iAtom);

          rho = (1/pripair.one_over_gamma)* nucShell[iAtom].alpha[0]
                /(1/pripair.one_over_gamma + nucShell[iAtom].alpha[0]);
 
          computeFmTTaylor(tmpFmT,rho*squarePC,lTotal+m,0);
          tmpVal =  2.0*sqrt(rho/M_PI)*ssS*
                    vRRV0b(nucShell,pripair,shell2,tmpFmT,PC,m,LB,lB,iAtom,molecule);

//  std::cerr<<"no finite nuclei"<<std::endl;  
        }
        else if ( !useFiniteWidthNuclei ) {
          computeFmTTaylor(tmpFmT,(shell1.alpha[pripair.p1]+shell2.alpha[pripair.p2])
                     *squarePC,lTotal+m,0);
  
          auto ssV = 2.0*sqrt(1.0/(pripair.one_over_gamma*M_PI))*ssS;
    
          tmpVal = ssV * vRRV0b(nucShell,pripair,shell2,tmpFmT,PC,m,LB,lB,iAtom,molecule); 
        }
        delete[] tmpFmT;
  
      }
    
    else if ((LA>0)&&(LB>0)) {
      int iWork,lAp1[3],lBm1[3];
      for( iWork = 0 ; iWork < 3 ; iWork++ ) {
        lAp1[iWork]=lA[iWork];
        lBm1[iWork]=lB[iWork];
      };
      if (lB[0] > 0) iWork = 0;
      else if (lB[1] > 0) iWork=1;
      else if (lB[2] > 0) iWork=2;
  
      lAp1[iWork]++;
      lBm1[iWork]--;
      tmpVal = hRRiPPVab(nucShell,pripair,shell1,shell2,LA+1,lAp1,LB-1,lBm1,C,m,iAtom,molecule);
      tmpVal+= (shell1.O[iWork]-shell2.O[iWork])
               *hRRiPPVab( nucShell,pripair,shell1,shell2,LA,lA,LB-1,lBm1,C,m,iAtom,molecule);
    }
    return tmpVal;
  }
  
  
  //----------------------------------------------------------------------------//
  // potential integral vertical recursion                                      //
  //  (0|0_c|b)^(m) = (P-B)(0|0_c|b-1i)^(m) - (P-C)(0|0_c|b-1i)^(m+1)           //
  //                + halfInvZeta*N_(b-1)*[(0|0_c|b-2i)^(m)-(0|0_c|b-2i)^(m+1)] //
  //  since LA == 0, we only decrease b                                           //
  //----------------------------------------------------------------------------//
  /**
   *  \brief Perform the vertical recurrence relation for the uncontracted 
   *    nuclear potential integral 
   *
   *   (0|0_c|b)^(m) = (P-B)(0|0_c|b-1i)^(m) - (P-C)(0|0_c|b-1i)^(m+1) 
   *                 + halfInvZeta*N_(b-1)*[(0|0_c|b-2i)^(m)-(0|0_c|b-2i)^(m+1)]
   *
   *  where a is angular momentum, Zeta=zeta_a+zeta_b, A is bra nuclear coordinate. 
   *  P = (zeta_a*A+zeta_b*B)/Zeta
   *
   *  \param [in] nucShell nuclear shell, give the exponents of gaussian function of nuclei
   *  \param [in] pripair Primitive Shell pair data for shell1, shell2
   *  \param [in] shell1  Bra shell
   *  \param [in] FmT     table of Boys function with different m
   *  \param [in] PC      vector P-C
   *  \param [in] m       the order of auxiliary function
   *  \param [in] LB      total Bra angular momentum
   *  \param [in] lB      Bra angular momentum vector (lAx,lAy,lAz)
   *  \param [in] iAtom   the index of the atom of the nuclei
   *
   *  \returns an uncontracted nuclear potential integral
   *
   */ 

   double RealGTOIntEngine::vRRV0b( const std::vector<libint2::Shell> &nucShell, 
    libint2::ShellPair::PrimPairData &pripair, libint2::Shell &shell2,
    double *FmT, double *PC, int m, int LB, int *lB, int iAtom, const Molecule& molecule){
  
    double rhoovzeta;
    bool useFiniteWidthNuclei = nucShell.size() > 0;

    if( LB == 0 ) {
       if ( useFiniteWidthNuclei ) {
//  std::cerr<<"out here"<<std::endl; 
  //      rhoovzeta =  RealGTOIntEngine::molecule_->nucShell(iAtom).alpha[0]/(RealGTOIntEngine::molecule_->nucShell(iAtom).alpha[0]+ijSP->Zeta[iPP]);
         rhoovzeta =  nucShell[iAtom].alpha[0]/( nucShell[iAtom].alpha[0]
                      +(1/pripair.one_over_gamma) );

//std::cout<<"tmpFmT "<<FmT[m]<<std::endl;
         return (pow(rhoovzeta,m)*FmT[m]);

  //     return (pow(rhoovzeta,m)*FmT[m]);
      } else if ( !useFiniteWidthNuclei ) {           
//std::cout<<"tmpFmT "<<FmT[m]<<std::endl;
        return FmT[m];  //Z*2*sqrt[zeta/pi]*[s|s]is given in hRRVab, in ssV.
      }
    }
  
    double tmpVal=0.0;
    int lBm1[3]; //means lA minus 1_i
    int iWork;
    for( iWork = 0 ; iWork < 3 ; iWork++ ) lBm1[iWork]=lB[iWork];
    if (lB[0] > 0) iWork = 0;
    else if (lB[1] > 0) iWork=1;
    else if (lB[2] > 0) iWork=2;
    lBm1[iWork]--;
  
    tmpVal  = (pripair.P[iWork]-shell2.O[iWork])*
                vRRV0b(nucShell,pripair,shell2,FmT,PC,m,LB-1,lBm1,iAtom,molecule);
  
    tmpVal -= PC[iWork]*vRRV0b(nucShell,pripair,shell2,FmT,PC,m+1,LB-1,lBm1,iAtom,molecule);
  
    if(lBm1[iWork] >=1){
      lBm1[iWork]--;
      tmpVal += (lBm1[iWork]+1)*0.5*pripair.one_over_gamma * 
                (vRRV0b(nucShell,pripair,shell2,FmT,PC,m,LB-2,lBm1,iAtom,molecule) 
                -vRRV0b(nucShell,pripair,shell2,FmT,PC,m+1,LB-2,lBm1,iAtom,molecule));
    }
    return tmpVal; 
  }
  
   
  //---------------------------------------------------------------------------------//
  //spin orbital vertical recursion                                                  //
  // if LA == 0,LB == 0, then [0|S|0]^(m)=4*zetaa*zetab*((A-C)X(B-C))_mu*[s|A(0)|s]^(m+1)//
  // if LB == 0,then [a|S|0]^(m)=(Pi-Ai)[a-1i|S|0]^(m)-(Pi-Ci)[a-1i|S|0]^(m+1)         //
  //                           +halfInvZeta*Ni(a-1i)([a-2i|S|0]^(m)-[a-2i|S|0]^(m+1))//
  //                           +2*zeta_b*{1i cross(B-C)}_mu*[a-1i|A(0)|0]^(m+1)      //
  // if LB>0,then  [a|S|b]^(m)=(Pi-Bi)[a|S|b-1i]^(m)-(Pi-Ci)[a|S|b-1i]^(m+1)         //
  //                  +halfInvZeta*Ni(b-1i)([a|S|b-2i]^(m)-[a|S|b-2i]^(m+1))         //
  //                  +halfInvZeta*Ni(a)([a-1i|S|b-1i]^(m)-[a-1i|S|b-1i]^(m+1))      //
  //                  -2*zeta_a*{1i cross(A-C)}_mu*[a|A(0)|b-1i]^(m+1)               //
  //                  -Sum_k=x,y,z N_k(a)*{1i cross 1k}_mu*[a-1k|A(0)|b-1i]^(m+1)    //
  //---------------------------------------------------------------------------------//
  /**
   *  \brief Perform the vertical recurrence relation for the uncontracted 
   *  spin orbit integral
   *
   *  if LA == 0,LB == 0, then [0|S|0]^(m)=4*zetaa*zetab*((A-C)X(B-C))_mu*[s|A(0)|s]^(m+1)
   *  if LB == 0,then [a|S|0]^(m)=(Pi-Ai)[a-1i|S|0]^(m)-(Pi-Ci)[a-1i|S|0]^(m+1)       
   *                            +halfInvZeta*Ni(a-1i)([a-2i|S|0]^(m)-[a-2i|S|0]^(m+1))
   *                            +2*zeta_b*{1i cross(B-C)}_mu*[a-1i|A(0)|0]^(m+1)      
   *  if LB>0,then  [a|S|b]^(m)=(Pi-Bi)[a|S|b-1i]^(m)-(Pi-Ci)[a|S|b-1i]^(m+1)         
   *                   +halfInvZeta*Ni(b-1i)([a|S|b-2i]^(m)-[a|S|b-2i]^(m+1))         
   *                   +halfInvZeta*Ni(a)([a-1i|S|b-1i]^(m)-[a-1i|S|b-1i]^(m+1))      
   *                   -2*zeta_a*{1i cross(A-C)}_mu*[a|A(0)|b-1i]^(m+1)               
   *                   -Sum_k=x,y,z N_k(a)*{1i cross 1k}_mu*[a-1k|A(0)|b-1i]^(m+1)    
   *
   *  where a,b are the angular momentum, A,B are the nuclear coordinates.
   *  halInvZeta = 1/2 * 1/Zeta
   *
   *  \param [in] nucShell nuclear shell, give the exponents of gaussian function of nuclei
   *  \param [in] pripair primitive Shell pair data for shell1, shell2
   *  \param [in] shell1  Bra shell
   *  \param [in] shell2  Ket shell
   *  \param [in] OneixAC 1_i cross AC vector is a 3 by 3 tensor
   *  \param [in] OneixBC 1_i cross BC vector is a 3 by 3 tensor
   *  \param [in] LA      total Bra angular momentum
   *  \param [in] lA      Bra angular momentum vector (lAx,lAy,lAz)
   *  \param [in] LB      total Ket angular momentum
   *  \param [in] lB      Ket angular momentum vector (lBx,lBy,lBz)
   *  \param [in] mu      index of the angular momentum integral component being calculated
   *  \param [in] m       order of auxiliary function
   *  \param [in] iAtom   index of atom
   *
   *  \returns the mu component of an uncontracted spin orbit integral
   *
   */ 

   double RealGTOIntEngine::Slabmu(const std::vector<libint2::Shell> &nucShell,
    libint2::ShellPair::PrimPairData &pripair, libint2::Shell &shell1,
    libint2::Shell &shell2, double *OneixAC, double *OneixBC, int LA, int *lA, 
    int LB, int *lB, int mu, int m, int iAtom, const Molecule& molecule){
  
    double tmpVal=0.0;
    double C[3],AC[3],BC[3],ACxBC[3];
    
    bool useFiniteWidthNuclei = nucShell.size() > 0;
  
    for (int k = 0 ; k < 3 ; k++) C[k] = molecule.atoms[iAtom].coord[k]; 
  
    if( (LA+LB)  == 0 ){ 
      for( int k = 0; k < 3; k++ ) {
        AC[k] = shell1.O[k] - C[k];
        BC[k] = shell2.O[k] - C[k];
      }
    
      ACxBC[0] = AC[1]*BC[2] - AC[2]*BC[1]; 
      ACxBC[1] = AC[2]*BC[0] - AC[0]*BC[2];
      ACxBC[2] = AC[0]*BC[1] - AC[1]*BC[0];
  
     // double u;
      auto u = hRRiPPVab(nucShell,pripair,shell1,shell2,LA,lA,LB,lB,C,m+1,iAtom,molecule);
      tmpVal = 4*shell1.alpha[pripair.p1]*shell2.alpha[pripair.p2]*ACxBC[mu]*u;        
  
      return tmpVal;
    };
  
    int lAm1[3];
    int iWork;
   
    double PC[3];
    for( int k = 0; k < 3; k++ ) PC[k] = pripair.P[k] - C[k];
  
    if( (LB  ==  0) and (LA > 0) ){
      for( iWork = 0; iWork < 3; iWork++) lAm1[iWork] = lA[iWork];
      if( lA[0] > 0 ) iWork = 0;
      else if( lA[1] > 0 ) iWork=1;
      else if( lA[2] > 0 ) iWork=2;
  
      lAm1[iWork]--;
      tmpVal  = (pripair.P[iWork]-shell1.O[iWork]) *
        Slabmu(nucShell,pripair,shell1,shell2,OneixAC,OneixBC,
               LA-1,lAm1,LB,lB,mu,m,iAtom,molecule);
      tmpVal -= PC[iWork]*
        Slabmu(nucShell,pripair,shell1,shell2,OneixAC,OneixBC,
               LA-1,lAm1,LB,lB,mu,m+1,iAtom,molecule);
      tmpVal += 2 * shell2.alpha[pripair.p2] * OneixBC[iWork*3+mu] *
        hRRiPPVab(nucShell,pripair,shell1,shell2,LA-1,lAm1,LB,lB,C,m+1,iAtom,molecule); 
  
      if( lAm1[iWork]  >=  1 ){
        lAm1[iWork]--;
        tmpVal += (lA[iWork]-1) * 0.5* pripair.one_over_gamma * 
          (  
             Slabmu(nucShell,pripair,shell1,shell2,OneixAC,OneixBC,LA-2,lAm1,
               LB,lB,mu,m,iAtom,molecule) -       
             Slabmu(nucShell,pripair,shell1,shell2,OneixAC,OneixBC,LA-2,lAm1,
               LB,lB,mu,m+1,iAtom,molecule)
          );
      }
  
      return tmpVal;
  
    } // else if( (LA > 0) and (LB  ==  0)) cout<<"here"; 
  //  else if( (LB > 0) and (LA > 0) ){
    else if( (LB > 0) and ( LA >= 0)){
      int lBm1[3],lAm1k[3];
      for( iWork = 0; iWork < 3; iWork++ ) {
        lAm1[iWork]=lA[iWork];
        lBm1[iWork]=lB[iWork];
        lAm1k[iWork]=lA[iWork];
      };
      if( lB[0] > 0 ) iWork = 0;
      else if( lB[1] > 0 ) iWork=1;
      else if( lB[2] > 0 ) iWork=2;
  
      lAm1[iWork]--;
      lBm1[iWork]--;
  
      tmpVal  = ( pripair.P[iWork]-shell2.O[iWork]) * 
        Slabmu(nucShell,pripair,shell1,shell2,OneixAC,OneixBC,
               LA,lA,LB-1,lBm1,mu,m,iAtom,molecule);
  
      tmpVal -= PC[iWork] * 
        Slabmu(nucShell,pripair,shell1,shell2,OneixAC,OneixBC,
               LA,lA,LB-1,lBm1,mu,m+1,iAtom,molecule);
  
      if( lAm1[iWork]  >=  0 )
        tmpVal += 0.5* pripair.one_over_gamma * lA[iWork] * 
          (
            Slabmu(nucShell,pripair,shell1,shell2,OneixAC,OneixBC,
                   LA-1,lAm1,LB-1,lBm1,mu,m,
              iAtom,molecule) - 
            Slabmu(nucShell,pripair,shell1,shell2,OneixAC,OneixBC,
                   LA-1,lAm1,LB-1,lBm1,mu,m+1,
              iAtom,molecule)
          );
  
      tmpVal -= 2 * shell1.alpha[pripair.p1] * OneixAC[3*iWork+mu] * 
        hRRiPPVab(nucShell,pripair,shell1,shell2,LA,lA,LB-1,lBm1,C,m+1,iAtom,molecule);
  
  
  /* DBWY
      for( k = 0; k < 3; k++) {
        if(lA[k] > 0){
          lAm1k[k]--;
          tmpVal -= lA[k] * (*OneiOnek)(iWork,k,mu) * 
            hRRiPPVab(ijSP,LA-1,lAm1k,LB-1,lBm1,C,m+1,iPP,iAtom);
          
          lAm1k[k]++;
        };
      };
  */
      if( iWork  ==  0 ) {
        if( mu  ==  1 and lA[2] > 0 ) {
          lAm1k[2]--;
          tmpVal += lA[2] * 
            hRRiPPVab(nucShell,pripair,shell1,shell2,LA-1,lAm1k,
                      LB-1,lBm1,C,m+1,iAtom,molecule);
          lAm1k[2]++;
        } else if( mu  ==  2 and lA[1] > 0 ) {
          lAm1k[1]--;
          tmpVal -= lA[1] * 
            hRRiPPVab(nucShell,pripair,shell1,shell2,LA-1,lAm1k,
                      LB-1,lBm1,C,m+1,iAtom,molecule);
          lAm1k[1]++;
        }
      } else if( iWork  ==  1 ) {
        if( mu  ==  2 and lA[0] > 0 ) {
          lAm1k[0]--;
          tmpVal += lA[0] * 
            hRRiPPVab(nucShell,pripair,shell1,shell2,LA-1,lAm1k,
                      LB-1,lBm1,C,m+1,iAtom,molecule);
          lAm1k[0]++;
        } else if( mu  ==  0 and lA[2] > 0 ) {
          lAm1k[2]--;
          tmpVal -= lA[2] * 
            hRRiPPVab(nucShell,pripair,shell1,shell2,LA-1,lAm1k,
                      LB-1,lBm1,C,m+1,iAtom,molecule);
          lAm1k[2]++;
        }
      } else {
        if( mu  ==  0 and lA[1] > 0 ) {
          lAm1k[1]--;
          tmpVal += lA[1] * 
            hRRiPPVab(nucShell,pripair,shell1,shell2,LA-1,lAm1k,
                      LB-1,lBm1,C,m+1,iAtom,molecule);
          lAm1k[1]++;
        } else if( mu  ==  1 and lA[0] > 0 ) {
          lAm1k[0]--;
          tmpVal -= lA[0] * 
            hRRiPPVab(nucShell,pripair,shell1,shell2,LA-1,lAm1k,
                      LB-1,lBm1,C,m+1,iAtom,molecule);
          lAm1k[0]++;
        }
      }
      
      if( lBm1[iWork] > 0 ) {
        lBm1[iWork]--;
        tmpVal += (lBm1[iWork] + 1) * 0.5* pripair.one_over_gamma * 
          (
             Slabmu(nucShell,pripair,shell1,shell2,OneixAC,OneixBC,
                    LA,lA,LB-2,lBm1,mu,m,iAtom,molecule) - 
             Slabmu(nucShell,pripair,shell1,shell2,OneixAC,OneixBC,
                    LA,lA,LB-2,lBm1,mu,m+1,iAtom,molecule)
          );
      }
    }
  
    return tmpVal;
  };
  
  
  //----------------------------------------------------------------------------------//
  //pVp vertical recursion                                                            //
  //                                                                                  //
  //if LA+LB == 0,then                                                                  //
  //[s|pVp|s]=(6*zeta_a*zeta_b/zeta+4*zeta_a*zeta_b*(P-A)dot(P-B)[0_A|A(0)|0_B]^(m)   //
  //     -[4*zeta_a*zeta_b*(2P-A-B)dot(P-C)+6*zeta_a*zeta_b/zeta]*[0_A|A(0)|0_B]^(m+1)//
  //     +4*zeta_a*zeta_b*(P-C)^2*[0_A|A(0)|0_B]^(m+2)                                //
  //                                                                                  //
  //if LB == 0, LA<>0, then                                                             //
  //[a|pVp|0]^(m)=(P-A)[a-1i|pVp|0]^(m)-(P-C)[a-1i|pVp|0]^(m+1)                       //
  //     +invZeta*Ni(a)*([a-2i|pVp|0]^(m)-[a-2i|pVp|0]^(m+1))                         //
  //     -zeta_b/zeta*2*zeta_b*[a|A|1i]^(m)                                           //
  //     +zeta_b/zeta*{2*zeta_a*[a+1i|A|0]^(m)-Ni(a)[a-1i|A|0]^(m)}                   //
  //     -zeta_a/zeta*2*zeta_b*[a|A|1i]^(m+1)                                         //
  //     -zeta_b/zeta*{2*zeta_a*[a+1i|A|0]^(m+1)-Ni(a)*[a-1i|A|0]^(m+1)}              //
  //                                                                                  //    
  //if LA == 0, LB<>0, then                                                             //
  //[0|pVp|b]^(m)=(P-B)[0|pVp|b-1i]^(m)-(P-C)[0|pVp|b-1i]^(m+1)                       //
  //     +halfinvzeta*Ni(b)*{[0|pVp|b-2i]^(m)-[0|pVp|b-2i]^(m+1)}                     //
  //     -zeta_a/zeta*2*zeta_a*[1i|A|b-1i]^(m)                                        //
  //     +zeta_a/zeta*{2*zeta_b*[0|A|b]^(m)-Ni(b)[0|A|b-2i]^(m)                       //
  //     -zeta_b/zeta*2*zeta_a*[1i|A|b-1i]^(m+1)                                      //
  //     -zeta_a/zeta*{2*zeta_b*[a|A|b]^(m+1)-Ni(b)[a|A|b-2i]^(m+1)}                  //
  //                                                                                  //
  //if LA>0, LB>0, then                                                               //
  //[a|pVp|b]^(m)=(P-B)[a|pVp|b-1i]^(m)-(P-C)[a|pVp|b-1i]^(m+1)                       //
  //     +halfinvzeta*Ni(b)*{[a|pVp|b-2i]^(m)-[a|pVp|b-2i]^(m+1)}                     //
  //     +halfinvzeta*Ni(a)*{[a-1i|pVp|b-1i]^(m)-[a-1i|pVp|b-1i]^(m+1)}               //
  //     -zeta_a/zeta*{2*zeta_a*[a+1i|A|b-1i]^(m)-Ni(a)*[a-1i|A|b-1i]^(m)}            //
  //     +zeta_a/zeta*{2*zeta_b*[a|A|b]^(m)-Ni(b)[a|A|b-2i]^(m)                       //
  //     -zeta_b/zeta*{2*zeta_a*[a+1i|A|b-1i]^(m+1)-Ni(a)[a-1i|A|b-1i]^(m+1)}         //
  //     -zeta_a/zeta*{2*zeta_b*[a|A|b]^(m+1)-Ni(b)[a|A|b-2i]^(m+1)}                  //
  //----------------------------------------------------------------------------------//

  /**
   *  \brief Perform the vertical recurrence relation for the uncontracted 
   *   p V dot p integral
   *   
   *  if LA+LB == 0,then                                                                
   *  [s|pVp|s]=(6*zeta_a*zeta_b/zeta+4*zeta_a*zeta_b*(P-A)dot(P-B)[0_A|A(0)|0_B]^(m)   
   *       -[4*zeta_a*zeta_b*(2P-A-B)dot(P-C)+6*zeta_a*zeta_b/zeta]*[0_A|A(0)|0_B]^(m+1)
   *       +4*zeta_a*zeta_b*(P-C)^2*[0_A|A(0)|0_B]^(m+2)                                
   *                                                                                    
   *  if LB == 0, LA<>0, then                                                           
   *  [a|pVp|0]^(m)=(P-A)[a-1i|pVp|0]^(m)-(P-C)[a-1i|pVp|0]^(m+1)                       
   *       +invZeta*Ni(a)*([a-2i|pVp|0]^(m)-[a-2i|pVp|0]^(m+1))                         
   *       -zeta_b/zeta*2*zeta_b*[a|A|1i]^(m)                                           
   *       +zeta_b/zeta*{2*zeta_a*[a+1i|A|0]^(m)-Ni(a)[a-1i|A|0]^(m)}                   
   *       -zeta_a/zeta*2*zeta_b*[a|A|1i]^(m+1)                                         
   *       -zeta_b/zeta*{2*zeta_a*[a+1i|A|0]^(m+1)-Ni(a)*[a-1i|A|0]^(m+1)}              
   *                                                                                    
   *  if LA == 0, LB<>0, then                                                           
   *  [0|pVp|b]^(m)=(P-B)[0|pVp|b-1i]^(m)-(P-C)[0|pVp|b-1i]^(m+1)                       
   *       +halfinvzeta*Ni(b)*{[0|pVp|b-2i]^(m)-[0|pVp|b-2i]^(m+1)}                     
   *       -zeta_a/zeta*2*zeta_a*[1i|A|b-1i]^(m)                                        
   *       +zeta_a/zeta*{2*zeta_b*[0|A|b]^(m)-Ni(b)[0|A|b-2i]^(m)                       
   *       -zeta_b/zeta*2*zeta_a*[1i|A|b-1i]^(m+1)                                      
   *       -zeta_a/zeta*{2*zeta_b*[a|A|b]^(m+1)-Ni(b)[a|A|b-2i]^(m+1)}                  
   *                                                                                    
   *  if LA>0, LB>0, then                                                               
   *  [a|pVp|b]^(m)=(P-B)[a|pVp|b-1i]^(m)-(P-C)[a|pVp|b-1i]^(m+1)                       
   *       +halfinvzeta*Ni(b)*{[a|pVp|b-2i]^(m)-[a|pVp|b-2i]^(m+1)}                     
   *       +halfinvzeta*Ni(a)*{[a-1i|pVp|b-1i]^(m)-[a-1i|pVp|b-1i]^(m+1)}               
   *       -zeta_a/zeta*{2*zeta_a*[a+1i|A|b-1i]^(m)-Ni(a)*[a-1i|A|b-1i]^(m)}            
   *       +zeta_a/zeta*{2*zeta_b*[a|A|b]^(m)-Ni(b)[a|A|b-2i]^(m)                       
   *       -zeta_b/zeta*{2*zeta_a*[a+1i|A|b-1i]^(m+1)-Ni(a)[a-1i|A|b-1i]^(m+1)}         
   *       -zeta_a/zeta*{2*zeta_b*[a|A|b]^(m+1)-Ni(b)[a|A|b-2i]^(m+1)}                  
   *
   *  where a,b are the angular momentum, A,B are the nuclear coordinates.
   *
   *  \param [in] nucShell nuclear shell, give the exponents of gaussian function of nuclei
   *  \param [in] pripair Primitive Shell pair data for shell1, shell2
   *  \param [in] shell1  Bra shell
   *  \param [in] shell2  Ket shell
   *  \param [in] LA      total Bra angular momentum
   *  \param [in] lA      Bra angular momentum vector (lAx,lAy,lAz)
   *  \param [in] LB      total Ket angular momentum
   *  \param [in] lB      Ket angular momentum vector (lBx,lBy,lBz)
   *  \param [in] m       order of auxiliary function
   *  \param [in] iAtom   index of nuclei
   *
   *  \returns a uncontracted pV dot p integral. It doesn't include contraction coeffs.
   *
   */ 

   double RealGTOIntEngine::pVpab(const std::vector<libint2::Shell> &nucShell,
    libint2::ShellPair::PrimPairData &pripair, libint2::Shell &shell1 , 
    libint2::Shell &shell2, int LA, int *lA, int LB, int *lB, int m, int iAtom, 
    const Molecule& molecule){
  
    double tmpVal=0.0;
    double C[3],PAPB,PABPC,PCsquare;

    bool useFiniteWidthNuclei = nucShell.size() > 0;

    for (int k = 0 ; k < 3 ; k++) C[k] = molecule.atoms[iAtom].coord[k];
  
    if ((LA+LB) == 0){
      PAPB = 0.0;
      PABPC = 0.0;
      PCsquare = 0.0;
      for ( int k = 0; k < 3 ; k++ ) {
        PAPB += ( pripair.P[k]-shell1.O[k] )*( pripair.P[k]- shell2.O[k] );
        PABPC += (2.0*pripair.P[k]- shell1.O[k] - shell2.O[k] )*( pripair.P[k] - C[k] );
        PCsquare += ( pripair.P[k] - C[k] )*( pripair.P[k] - C[k] );
      }
  
    double u,uone,utwo;
    u = hRRiPPVab(nucShell,pripair,shell1,shell2,LA,lA,LB,lB,C,m,iAtom,molecule);
    uone = hRRiPPVab(nucShell,pripair,shell1,shell2,LA,lA,LB,lB,C,m+1,iAtom,molecule);
    utwo = hRRiPPVab(nucShell,pripair,shell1,shell2,LA,lA,LB,lB,C,m+2,iAtom,molecule);
  
    tmpVal = (6 * shell1.alpha[pripair.p1] * shell2.alpha[pripair.p2] * pripair.one_over_gamma 
             +4 * shell1.alpha[pripair.p1] * shell2.alpha[pripair.p2] * PAPB ) * u;
  
    tmpVal -= (PABPC * 4 * shell1.alpha[pripair.p1] * shell2.alpha[pripair.p2] 
              + 6 * shell1.alpha[pripair.p1] * shell2.alpha[pripair.p2] 
              * pripair.one_over_gamma ) * uone;
              
    tmpVal += 4 * shell1.alpha[pripair.p1] * shell2.alpha[pripair.p2] * PCsquare * utwo;
  
    return tmpVal;
    } // if ((LA+LB) == 0)
    
    int lAm1[3],lAp1[3],onei[3];
    int iWork;
  
    double PC[3];
    for ( int k = 0 ; k < 3 ; k++ ) {
      PC[k] = pripair.P[k] - C[k];
    }
  
    if ((LB == 0)&(LA>0)){
      for ( iWork = 0 ; iWork < 3 ; iWork++ ) {
        lAm1[iWork] = lA[iWork];
        onei[iWork] = 0;
      }
  
      if (lA[0] > 0) iWork = 0;
      else if (lA[1] > 0)  iWork = 1;
      else if (lA[2] > 0)  iWork = 2;
      lAm1[iWork]--;
      onei[iWork]=1;
  
      tmpVal = (pripair.P[iWork]-shell1.O[iWork]) 
               * pVpab(nucShell,pripair,shell1,shell2,LA-1,lAm1,LB,lB,m,iAtom,molecule);
  
      tmpVal -= PC[iWork] * 
               pVpab(nucShell,pripair,shell1,shell2,LA-1,lAm1,LB,lB,m+1,iAtom,molecule);
  
      tmpVal -= shell2.alpha[pripair.p2] * pripair.one_over_gamma * 2 * 
                shell2.alpha[pripair.p2] * 
                hRRiPPVab(nucShell,pripair,shell1,shell2,LA-1,lAm1,
                LB+1,onei,C,m,iAtom,molecule);
  
      tmpVal += shell2.alpha[pripair.p2] * pripair.one_over_gamma * 
                ( 2 * shell1.alpha[pripair.p1] * 
                hRRiPPVab(nucShell,pripair,shell1,shell2,LA,lA,LB,lB,C,m,iAtom,molecule));
  
      tmpVal -= shell1.alpha[pripair.p1] * pripair.one_over_gamma * 2 * 
                shell2.alpha[pripair.p2] * 
                hRRiPPVab(nucShell,pripair,shell1,shell2,
                LA-1,lAm1,LB+1,onei,C,m+1,iAtom,molecule);
  
      tmpVal -= shell2.alpha[pripair.p2] * pripair.one_over_gamma * 
                (2 * shell1.alpha[pripair.p1] * 
                hRRiPPVab(nucShell,pripair,shell1,shell2,LA,lA,LB,lB,C,m+1,iAtom,molecule));
  
      if (lAm1[iWork] > 0) {
        lAm1[iWork]--;
        tmpVal += 0.5 * pripair.one_over_gamma * (lA[iWork]-1) * 
                 ( pVpab(nucShell,pripair,shell1,shell2,LA-2,lAm1,LB,lB,m,iAtom,molecule)
                  -pVpab(nucShell,pripair,shell1,shell2,LA-2,lAm1,LB,lB,m+1,iAtom,molecule) );
  
        tmpVal -= shell2.alpha[pripair.p2] * pripair.one_over_gamma * (lAm1[iWork]+1) * 
                  hRRiPPVab(nucShell,pripair,shell1,shell2,LA-2,lAm1,LB,lB,C,m,iAtom,molecule);
  
        tmpVal += shell2.alpha[pripair.p2] * pripair.one_over_gamma * (lAm1[iWork]+1) * 
                  hRRiPPVab(nucShell,pripair,shell1,shell2,
                            LA-2,lAm1,LB,lB,C,m+1,iAtom,molecule);
      } //if (lAm1[iWork] > 0)
    } // if ((LB == 0)&(LA>0))
  
  //  else if ((LA == 0)&(LB>0)) {
  //cout<<"here we are wrong"<<endl;
  //  }
  
  //  else if ((LA>0)&(LB)>0) {
    else if ( (LA >= 0) & (LB>0) ) {
      int lBm1[3],lAp1[3];
      for ( iWork = 0 ; iWork < 3 ; iWork++) {
        lAm1[iWork] = lA[iWork];
        lBm1[iWork] = lB[iWork];
        lAp1[iWork] = lA[iWork];
      }
      if (lB[0] > 0) iWork = 0;
      else if (lB[1] > 0) iWork = 1;
      else if (lB[2] > 0) iWork = 2;
      
      lAm1[iWork]--;
      lBm1[iWork]--;
      lAp1[iWork]++;
  
      tmpVal = ( pripair.P[iWork] - shell2.O[iWork] ) * 
               pVpab(nucShell,pripair,shell1,shell2,LA,lA,LB-1,lBm1,m,iAtom,molecule);
  
      tmpVal -= PC[iWork] * pVpab(nucShell,pripair,shell1,shell2,
                LA,lA,LB-1,lBm1,m+1,iAtom,molecule);
  
      tmpVal -= shell1.alpha[pripair.p1] * pripair.one_over_gamma * 2 * 
                shell1.alpha[pripair.p1] * 
                hRRiPPVab(nucShell,pripair,shell1,shell2,
                LA+1,lAp1,LB-1,lBm1,C,m,iAtom,molecule);
  
      tmpVal += shell1.alpha[pripair.p1] * pripair.one_over_gamma * 2 * 
                shell2.alpha[pripair.p2] * 
                hRRiPPVab(nucShell,pripair,shell1,shell2,LA,lA,LB,lB,C,m,iAtom,molecule);
  
      tmpVal -= shell2.alpha[pripair.p2] * pripair.one_over_gamma * 2 * 
                shell1.alpha[pripair.p1] * 
                hRRiPPVab(nucShell,pripair,shell1,shell2,
                LA+1,lAp1,LB-1,lBm1,C,m+1,iAtom,molecule);
  
      tmpVal -= shell1.alpha[pripair.p1] * pripair.one_over_gamma * 2 * 
                shell2.alpha[pripair.p2] * 
                hRRiPPVab(nucShell,pripair,shell1,shell2,LA,lA,LB,lB,C,m+1,iAtom,molecule);
  
      if (lA[iWork] > 0) {
        tmpVal += 0.5 * pripair.one_over_gamma * lA[iWork] * 
                ( pVpab(nucShell,pripair,shell1,shell2,LA-1,lAm1,LB-1,lBm1,m,iAtom,molecule)
                 -pVpab(nucShell,pripair,shell1,shell2,
                        LA-1,lAm1,LB-1,lBm1,m+1,iAtom,molecule));
  
        tmpVal += shell1.alpha[pripair.p1] * pripair.one_over_gamma * lA[iWork] * 
                  hRRiPPVab(nucShell,pripair,shell1,shell2,
                  LA-1,lAm1,LB-1,lBm1,C,m,iAtom,molecule);
  
        tmpVal += shell2.alpha[pripair.p2] * pripair.one_over_gamma * lA[iWork] * 
                  hRRiPPVab(nucShell,pripair,shell1,shell2,
                  LA-1,lAm1,LB-1,lBm1,C,m+1,iAtom,molecule);
  
      } // if (lA[iWork] > 0)
  
      if (lBm1[iWork] > 0) {
        lBm1[iWork]--;
        tmpVal += 0.5 * pripair.one_over_gamma * (lB[iWork]-1) * 
                  (pVpab(nucShell,pripair,shell1,shell2,LA,lA,LB-2,lBm1,m,iAtom,molecule)
                  -pVpab(nucShell,pripair,shell1,shell2,LA,lA,LB-2,lBm1,m+1,iAtom,molecule));
  
        tmpVal -= shell1.alpha[pripair.p1] * pripair.one_over_gamma * (lB[iWork]-1) * 
                  hRRiPPVab(nucShell,pripair,shell1,shell2,LA,lA,LB-2,lBm1,C,m,iAtom,molecule);
  
        tmpVal += shell1.alpha[pripair.p1] * pripair.one_over_gamma * (lB[iWork]-1) * 
                  hRRiPPVab(nucShell,pripair,shell1,shell2,
                  LA,lA,LB-2,lBm1,C,m+1,iAtom,molecule);
      }  // if (lBm1[iWork] > 0)
  
    } // else if ( (LA >= 0) & (LB)>0 )
    return tmpVal;
  }; // pVpab

}; //namespace ChronusQ 

