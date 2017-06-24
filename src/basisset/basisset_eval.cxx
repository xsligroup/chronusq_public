#include <memmanager.hpp>
#include <basisset/basisset_util.hpp>

// Basis_DEBUG_LEVEL >= 3 - Print EVERYTHING 
#ifndef Basis_DEBUG_LEVEL
#  define Basis_DEBUG_LEVEL 0
#endif

namespace ChronusQ {

  /**
   *  Level 1 Basis Set Evaluation Function (only for debug)
   *  Evaluates a shell set over a specified number of cartesian points.
   *  All distances of the points from each shell origin are computed and used
   *  in the Level 2 Basis Set Evaluation Function
   *  \param [in] typ        Type of evaluation to perform (gradient, etc)
   *  \param [in] shells     Shell set for evaluation(vector of libint2::Shell).
   *  \param [in] pts        Raw storage of the cartesian points (dimension 3 * npts)
   *  \param [in] npts       Number of cartesian points to evaluate.
   *  \param [in/out] fEval  eval Storage for the shell set evaluation, f(ixyz,iSh,ipt). 
   *                         Variable dimensions(npts*NBasisEff*1 or 4) allocated outside.
   *                         This storage will have all values of the functions in the shell, 
   *                         for each shell, for each point. If requested there will be appended 
   *                         f/dx(ixyz,iSh,ipt), f/dy(ixyz,iSh,ipt), f/dz(ixyz,iSh,ipt) values of 
   *                         the functions in the shell, for each shell(nShSize), for each point(npts)
   */ 
  void evalShellSet(SHELL_EVAL_TYPE typ, std::vector<libint2::Shell> &shells, 
    double *pts, size_t npts, double *fEval, bool forceCart, double * SCR) {
    size_t NBasisEff = 0;
    size_t nShSize = shells.size();
    double * r;
    double * rSq;
    if(!SCR)
    {
      r   = CQMemManager::get().malloc<double>(3*npts*nShSize);
      rSq = CQMemManager::get().malloc<double>(npts*nShSize);
    }
    else
    {
      r = SCR;
      rSq = SCR + 3*npts*nShSize;
    }
    // figure the size (number of basis) of the all shells inputed that need to be evaluated
    // and store in NBasisEff. It will be used for defining the pointers later on.
    int LMax = 0;
    for (auto iSh = 0; iSh < nShSize; iSh++){
      LMax = std::max(shells[iSh].contr[0].l,LMax);
    }
    size_t shSizeCar = ((LMax+1)*(LMax+2))/2; 
    size_t NDer = 1;
    if (typ ==GRADIENT) NDer = 4;
    double * SCR_Car;
    if(!SCR)
    {
      SCR_Car = CQMemManager::get().malloc<double>(NDer*shSizeCar);
    }
    else
    {
      SCR_Car = SCR + 3*npts*nShSize + npts*nShSize;
    }
        
    for (auto ipts = 0; ipts < npts; ipts++){
      double xp = *(pts + ipts*3);
      double yp = *(pts + 1 + ipts*3);
      double zp = *(pts + 2 + ipts*3);
      // we need this counter to keep track of the nbasis evaluated for all shells (at a given point)
      for (auto iSh = 0; iSh < nShSize; iSh++){
        r[   iSh*3 + ipts*3*nShSize] =  xp - shells[iSh].O[0];
        r[1+ iSh*3 + ipts*3*nShSize] =  yp - shells[iSh].O[1];
        r[2+ iSh*3 + ipts*3*nShSize] =  zp - shells[iSh].O[2];
        rSq[ iSh + ipts*nShSize] = r[    iSh*3 + ipts*3*nShSize]*r[    iSh*3 + ipts*3*nShSize] + 
                                       r[1 + iSh*3 + ipts*3*nShSize]*r[1 + iSh*3 + ipts*3*nShSize] + 
                                       r[2 + iSh*3 + ipts*3*nShSize]*r[2 + iSh*3 + ipts*3*nShSize]; 
      } // loop over shells
    } // loop over points

    std::vector<size_t> mapSh2Cen(nShSize); 
    std::iota(mapSh2Cen.begin(),mapSh2Cen.end(),0);

    std::vector<bool> evalShell(nShSize,true);
    // Call to Level 2 Basis Set Evaluation
    evalShellSet(typ,shells,evalShell,rSq,r,npts,nShSize,mapSh2Cen,
      NBasisEff,fEval,SCR_Car,shSizeCar,forceCart); 
    if(!SCR)
    {
      CQMemManager::get().free(r,rSq,SCR_Car);
    }

  }; // evalShellSet Level 1

  /**
   *  \brief Level 2 Basis Set Evaluation Function - Used in the KS - DFT
   *  Evaluates a shell set over a specified number of cartesian points. This function requires a precomputed
   *  set of the distances and their x,y,z component for each point from each shell origin in the shells vector.
   *  \param [in] typ        Type of evaluation to perform (gradient, etc)
   *  \param [in] shells     Shell set for evaluation(vector of libint2::Shell).
   *  \param [in] evalshells Vector of bool to know if that Shell is relevant  for evaluation.
   *  \param [in] rSq        Raw storage of overall squared distances between each point and the shell origin (precomputed outside)
   *                         Dimension (nshells * npts)
   *  \param [in] r          Raw storage of overall x,y,z of the vector between each point and the shell origin 
   *                         (precomputed outside). Dimension (3*nshells*npts) 
   *  \param [in] npts       Number of cartesian points to evaluate.
   *  \param [in] nCenter    Number of total distint shell center (is NAtoms).
   *  \param [in] NBasisEff  Number of basis function to be evaluated given all the shells in input.
   *  \param [in] mapSh2Cen  Shell Mapping to atom centers.
   *  \param [in/out] fEval  eval Storage for the shell set evaluation, f(ixyz,iSh,ipt). 
   *                         Variable dimensions(npts*NBasisEff*1 or 4) allocated outside.
   *                         This storage will have all values of the functions in the shell, 
   *                         for each shell, for each point. If requested there will be appended 
   *                         f/dx(ixyz,iSh,ipt), f/dy(ixyz,iSh,ipt), f/dz(ixyz,iSh,ipt) values of 
   *                         the functions in the shell, for each shell(nShSize), for each point(npts)
   *  \param [in] SCR        eval Storage for the shell set evaluation, Cartesian f(MaxShellSize). 
   *  \param [in] IOffSCR    offset for eval Storage for the shell set evaluation, Cartesian f(MaxShellSize). 
   *  \param [in] forceCart  True if force cartesian, otherwise sperical evaluation
   */ 
  void evalShellSet(SHELL_EVAL_TYPE typ, std::vector<libint2::Shell> &shells, 
    std::vector<bool> &evalShell, double* rSq, double *r, size_t npts, size_t nCenter, 
    std::vector<size_t> &mapSh2Cen, size_t NBasisEff, double *fEval, double *SCR, size_t IOffSCR, bool forceCart) {

    assert(shells.size() == evalShell.size());

    size_t nShSize = shells.size();
    size_t IOff =  npts*NBasisEff;
    std::array<double,3> rVal;

    for (auto ipts = 0ul; ipts < npts; ipts++){
      size_t Ic = 0;
    for (auto iSh = 0ul; iSh < nShSize; iSh++){
      if(evalShell[iSh]) {
        double * fStart    = fEval + Ic + ipts*NBasisEff;

        rVal [0]= r[0 + mapSh2Cen[iSh]*3 + ipts*3*nCenter];
        rVal [1]= r[1 + mapSh2Cen[iSh]*3 + ipts*3*nCenter];
        rVal [2]= r[2 + mapSh2Cen[iSh]*3 + ipts*3*nCenter];

        evalShellSet(typ,shells[iSh],rSq[mapSh2Cen[iSh] + ipts*nCenter],
          rVal,SCR,IOffSCR); 

        CarToSpDEval(typ, shells[iSh].contr[0].l, SCR, fStart, IOff, IOffSCR, forceCart);

        Ic += shells[iSh].size(); // Increment offset in basis
      }

    } // loop over shells
    } // loop over points

  }; // evalShellSet Level 2

//SS start 
  // level two evel shellset for GIAO 
  void evalShellSet(SHELL_EVAL_TYPE typ, std::vector<libint2::Shell> &shells, 
    std::vector<bool> &evalShell, double* rSq, double *r, size_t npts, size_t nCenter, 
    std::vector<size_t> &mapSh2Cen, size_t NBasisEff, dcomplex *fEval, dcomplex *SCR, 
    size_t IOffSCR, bool forceCart, EMPerturbation &pert) {

    assert(shells.size() == evalShell.size());

    size_t nShSize = shells.size();
    size_t IOff =  npts*NBasisEff;
    std::array<double,3> rVal;

    for (auto ipts = 0ul; ipts < npts; ipts++){
      size_t Ic = 0;
    for (auto iSh = 0ul; iSh < nShSize; iSh++){
      if(evalShell[iSh]) {
        dcomplex * fStart    = fEval + Ic + ipts*NBasisEff;

        rVal [0]= r[0 + mapSh2Cen[iSh]*3 + ipts*3*nCenter];
        rVal [1]= r[1 + mapSh2Cen[iSh]*3 + ipts*3*nCenter];
        rVal [2]= r[2 + mapSh2Cen[iSh]*3 + ipts*3*nCenter];

        // calculate wave vector

        double k[3];

        auto magAmp = pert.getDipoleAmp(Magnetic);

        k[0] = 0.5*( shells[iSh].O[1]*magAmp[2] - shells[iSh].O[2]*magAmp[1] );
        k[1] = 0.5*( shells[iSh].O[2]*magAmp[0] - shells[iSh].O[0]*magAmp[2] );
        k[2] = 0.5*( shells[iSh].O[0]*magAmp[1] - shells[iSh].O[1]*magAmp[0] );

        evalShellSet(typ,shells[iSh],rSq[mapSh2Cen[iSh] + ipts*nCenter],
          rVal,SCR,IOffSCR,k); 

        CarToSpDEval(typ, shells[iSh].contr[0].l, SCR, fStart, IOff, IOffSCR, forceCart);

        Ic += shells[iSh].size(); // Increment offset in basis
      }

    } // loop over shells
    } // loop over points

  }; // evalShellSet Level 2
//SS end

  /**
   *   \brief Level 3 Basis Set Evaluation Function
   *   Evaluates a single shell over a single cartesian point. This function requires a precomputed
   *   the distance and its x,y,z components for the point from the shell origin. An offset
   *   to properly store the results can be used.
   *   \param [in] typ        Type of evaluation to perform (gradient, etc)
   *   \param [in] shell      Shell for evaluation(libint2::Shell).
   *   \param [in] rSq        Raw storage of square distance between the point and the shell origin (precomputed outside)
   *   \param [in] r          Raw storage of x,y,z of the vector between each point and the shell origin 
   *                          (precomputed outside).  
   *   \param [in/out] SCR    eval Storage for the shell set evaluation, f(ixyz,iSh,ipt). 
   *                          Variable dimensions(npts*NBasisEff*1 or 4) allocated outside.
   *                          This storage will have all values of the functions in the shell, 
   *                          for each shell, for each point. If requested there will be appended 
   *                          f/dx(ixyz,iSh,ipt), f/dy(ixyz,iSh,ipt), f/dz(ixyz,iSh,ipt) values of 
   *                          the functions in the shell, for each shell(nShSize), for each point(npts)
   *   \param [in] IOffSCR       OffSet to properly store the basis set. 
   */ 
  void evalShellSet(SHELL_EVAL_TYPE typ, const libint2::Shell &shell,double rSq, const std::array<double,3> &xyz, 
    double *SCR, size_t IOffSCR) {
    auto L         = shell.contr[0].l;
    auto shSize    = ((L+1)*(L+2))/2; 
    auto shSize_car   = ((L+1)*(L+2))/2; 
    double * f_car     = SCR ;
    double * dx_car = f_car   + IOffSCR;
    double * dy_car = dx_car  + IOffSCR;
    double * dz_car = dy_car  + IOffSCR;
    auto contDepth = shell.alpha.size(); 
    double alpha(0.0);
    double expFactor(0.0);
    double expArg(0);
    double tmpcoef,tmpalpha;
    int lx,ly,lz, ixyz;
    double tmpxyz;
    double tmpdx;
    double tmpdy;
    double tmpdz;
    // Generating the expArgument, expFactotr and the
    // alpha (for derivatives later on) and store them
    // in temp variables
    for(auto k = 0; k < contDepth; k++){
      tmpcoef = shell.contr[0].coeff[k];
      tmpalpha = shell.alpha[k];
      expArg = std::exp(-tmpalpha*rSq);
      expFactor += tmpcoef * expArg;
      if (typ == GRADIENT) { 
        // quantities for derivatives
        tmpcoef *= tmpalpha;
        alpha += tmpcoef * expArg;
      }
    } 

    if (typ == GRADIENT) alpha *= 2;

    for(auto i = 0u, I = 0u; i <= L; i++) {
      lx = L - i;
      for( auto j = 0u; j <= i; j++, I++) {
        ly = i - j;
        lz = L - lx - ly;
        tmpxyz= 1.0;
        tmpdx = 0.0;
        tmpdy = 0.0;
        tmpdz = 0.0;
        for(ixyz = 0; ixyz < lx-1; ixyz++) tmpxyz *= xyz[0];
        for(ixyz = 0; ixyz < ly-1; ixyz++) tmpxyz *= xyz[1];
        for(ixyz = 0; ixyz < lz-1; ixyz++) tmpxyz *= xyz[2];
        f_car[I]  =  tmpxyz;

        if (typ == GRADIENT) {
        // Derivatives
          if(lx> 0) {tmpdx = expFactor * lx;}
          if(ly> 0) {tmpdy = expFactor * ly;}
          if(lz> 0) {tmpdz = expFactor * lz;}
           
          dx_car[I] = tmpxyz*tmpdx;
          dy_car[I] = tmpxyz*tmpdy;
          dz_car[I] = tmpxyz*tmpdz;
    
          // finishing up        
          if(lx> 0) {f_car[I]  *= xyz[0]; dy_car[I] *=xyz[0];dz_car[I] *=xyz[0];}
          if(ly> 0) {f_car[I]  *= xyz[1]; dx_car[I] *=xyz[1];dz_car[I] *=xyz[1];}
          if(lz> 0) {f_car[I]  *= xyz[2]; dx_car[I] *=xyz[2];dy_car[I] *=xyz[2];}
    
          dx_car[I] -= f_car[I] * xyz[0] * alpha;
          dy_car[I] -= f_car[I] * xyz[1] * alpha;
          dz_car[I] -= f_car[I] * xyz[2] * alpha;
          f_car[I]  *= expFactor;
    
        } else{
        // Only basis (not GGA)
          if(lx> 0) {f_car[I]  *= xyz[0];}
          if(ly> 0) {f_car[I]  *= xyz[1];}
          if(lz> 0) {f_car[I]  *= xyz[2];}
          f_car[I]  *= expFactor;
        }

#if Basis_DEBUG_LEVEL >= 3
          // Debug Printing
          std::cerr << I <<" "<< lx << " " 
            << ly << " "<<lz <<"  f(pt) "<< f_car[I] <<" " <<std::endl;
          std::cerr << I<<" "<< lx << " "  
            << ly << " "<<lz <<" dx(pt) "<< dx_car[I] << std::endl;
          std::cerr << I<<" "<< lx << " " 
            << ly << " "<<lz <<" dy(pt) "<< dy_car[I] << std::endl;
          std::cerr << I<<" "<< lx << " " 
            << ly << " "<<lz <<" dz(pt) "<< dz_car[I] << std::endl;
#endif

      } //loop overj, j[0,i]
    } //loop over i, i[0,L] this to loop required to build the lx,ly,lz combination given L

  }; // evalShellSet Level3

// SS start 
  /**
   *   \brief Level 3 Basis Set Evaluation Function for GIAO 
   *   Evaluates a single shell over a single cartesian point. This function requires a precomputed
   *   the distance and its x,y,z components for the point from the shell origin. An offset
   *   to properly store the results can be used.
   *   \param [in] typ        Type of evaluation to perform (gradient, etc)
   *   \param [in] shell      Shell for evaluation(libint2::Shell).
   *   \param [in] rSq        Raw storage of square distance between the point and the shell origin (precomputed outside)
   *   \param [in] r          Raw storage of x,y,z of the vector between each point and the shell origin 
   *                          (precomputed outside).  
   *   \param [in/out] SCR    eval Storage for the shell set evaluation, f(ixyz,iSh,ipt). 
   *                          Variable dimensions(npts*NBasisEff*1 or 4) allocated outside.
   *                          This storage will have all values of the functions in the shell, 
   *                          for each shell, for each point. If requested there will be appended 
   *                          f/dx(ixyz,iSh,ipt), f/dy(ixyz,iSh,ipt), f/dz(ixyz,iSh,ipt) values of 
   *                          the functions in the shell, for each shell(nShSize), for each point(npts)
   *   \param [in] IOffSCR       OffSet to properly store the basis set. 
   */ 
  void evalShellSet(SHELL_EVAL_TYPE typ, const libint2::Shell &shell,double rSq, const std::array<double,3> &xyz, 
    dcomplex *SCR, size_t IOffSCR, double *k) {
    auto L         = shell.contr[0].l;
    auto shSize    = ((L+1)*(L+2))/2; 
    auto shSize_car   = ((L+1)*(L+2))/2; 
    dcomplex * f_car     = SCR ;
    dcomplex * dx_car = f_car   + IOffSCR;
    dcomplex * dy_car = dx_car  + IOffSCR;
    dcomplex * dz_car = dy_car  + IOffSCR;
    auto contDepth = shell.alpha.size(); 
    dcomplex alpha=0.0;
    dcomplex expFactor=0.0;
    dcomplex expArg=0.0;
    double tmpcoef,tmpalpha;
    int lx,ly,lz, ixyz;
    double tmpxyz;
    dcomplex tmpdx;
    dcomplex tmpdy;
    dcomplex tmpdz;
    // Generating the expArgument, expFactotr and the
    // alpha (for derivatives later on) and store them
    // in temp variables
    
    // calculate the phase factor 
    dcomplex phase=0.0;
    dcomplex onei;
    onei.real(0.0);
    onei.imag(1.0);
    double kdotr = 0.0;
    for ( int ii = 0 ; ii < 3 ; ii++ ) kdotr += k[ii]*xyz[ii];
    phase = std::exp(onei*kdotr); 

    for(auto kk = 0; kk < contDepth; kk++){
      tmpcoef = shell.contr[0].coeff[kk];
      tmpalpha = shell.alpha[kk];
      expArg = phase * std::exp(-tmpalpha*rSq);
      expFactor += tmpcoef * expArg;
      if (typ == GRADIENT) { 
        // quantities for derivatives
        tmpcoef *= tmpalpha;
        alpha += tmpcoef * expArg;
      }
    } 

    if (typ == GRADIENT) alpha *= 2;

    for(auto i = 0u, I = 0u; i <= L; i++) {
      lx = L - i;
      for( auto j = 0u; j <= i; j++, I++) {   // I count the number of elements in a shell
        ly = i - j;
        lz = L - lx - ly;
        tmpxyz= 1.0;
        tmpdx = 0.0;
        tmpdy = 0.0;
        tmpdz = 0.0;
        for(ixyz = 0; ixyz < lx-1; ixyz++) tmpxyz *= xyz[0];
        for(ixyz = 0; ixyz < ly-1; ixyz++) tmpxyz *= xyz[1];
        for(ixyz = 0; ixyz < lz-1; ixyz++) tmpxyz *= xyz[2];
        f_car[I]  =  tmpxyz;

        if (typ == GRADIENT) {
        // Derivatives
/*
          if(lx> 0) {tmpdx = -expFactor * static_cast<dcomplex>( lx );}
          if(ly> 0) {tmpdy = -expFactor * static_cast<dcomplex>( ly );}
          if(lz> 0) {tmpdz = -expFactor * static_cast<dcomplex>( lz );}
*/

          if(lx> 0) {tmpdx = expFactor * static_cast<dcomplex>( lx );}
          if(ly> 0) {tmpdy = expFactor * static_cast<dcomplex>( ly );}
          if(lz> 0) {tmpdz = expFactor * static_cast<dcomplex>( lz );}

           
          dx_car[I] = tmpxyz*tmpdx;
          dy_car[I] = tmpxyz*tmpdy;
          dz_car[I] = tmpxyz*tmpdz;
    
          // finishing up        
          if(lx> 0) {f_car[I]  *= xyz[0]; dy_car[I] *=xyz[0];dz_car[I] *=xyz[0];}
          if(ly> 0) {f_car[I]  *= xyz[1]; dx_car[I] *=xyz[1];dz_car[I] *=xyz[1];}
          if(lz> 0) {f_car[I]  *= xyz[2]; dx_car[I] *=xyz[2];dy_car[I] *=xyz[2];}

/*    
          dx_car[I] += f_car[I] * xyz[0] * alpha;
          dy_car[I] += f_car[I] * xyz[1] * alpha;
          dz_car[I] += f_car[I] * xyz[2] * alpha;
*/

          dx_car[I] -= f_car[I] * xyz[0] * alpha;
          dy_car[I] -= f_car[I] * xyz[1] * alpha;
          dz_car[I] -= f_car[I] * xyz[2] * alpha;



          f_car[I]  *= expFactor;
          // extra term for GIAO 
/*
          dx_car[I] += -onei*k[0]*f_car[I];
          dy_car[I] += -onei*k[1]*f_car[I];
          dz_car[I] += -onei*k[2]*f_car[I];
*/

          dx_car[I] += onei*k[0]*f_car[I];
          dy_car[I] += onei*k[1]*f_car[I];
          dz_car[I] += onei*k[2]*f_car[I];
    
        } else{
        // Only basis (not GGA)
          if(lx> 0) {f_car[I]  *= xyz[0];}
          if(ly> 0) {f_car[I]  *= xyz[1];}
          if(lz> 0) {f_car[I]  *= xyz[2];}
          f_car[I]  *= expFactor;
        }

#if Basis_DEBUG_LEVEL >= 3
          // Debug Printing
          std::cout << I <<" "<< lx << " " 
            << ly << " "<<lz <<"  f(pt) "<< std::real(f_car[I]) <<" " <<std::endl;
          std::cout << I<<" "<< lx << " "  
            << ly << " "<<lz <<" dx(pt) "<< std::real(dx_car[I]) << std::endl;
          std::cout << I<<" "<< lx << " " 
            << ly << " "<<lz <<" dy(pt) "<< std::real(dy_car[I]) << std::endl;
          std::cout << I<<" "<< lx << " " 
            << ly << " "<<lz <<" dz(pt) "<< std::real(dz_car[I]) << std::endl;
#endif

      } //loop overj, j[0,i]
    } //loop over i, i[0,L] this to loop required to build the lx,ly,lz combination given L

  }; // evalShellSet Level3

//SS end 


  /**
   *   \brief Basis Set Cartesian to Sperical conversion over a single shell.
   *   It needs to be called also for cartisian evaluation, since this function
   *   handles the final population of the eval Storage for the shell and points.
   *   \param [in] typ       Type of evaluation to perform (gradient, etc)
   *   \param [in] L         Angular quantum momentum
   *   \param [in] fCarEVal  eval Storage for the shell set evaluation, f(ixyz,iSh,ipt). 
   *                         Variable dimensions(npts*NBasisEff*1 or 4) allocated outside.
   *                         This storage will have all values of the functions in the shell, 
   *                         for each shell, for each point. If requested there will be appended 
   *                         f/dx(ixyz,iSh,ipt), f/dy(ixyz,iSh,ipt), f/dz(ixyz,iSh,ipt) values of 
   *                         the functions in the shell, for each shell(nShSize), for each point(npts)
   *  \param [in/out] fEval  eval Storage for the shell set evaluation, f(ixyz,iSh,ipt). 
   *                         Variable dimensions(npts*NBasisEff*1 or 4) allocated outside.
   *                         This storage will have all values of the functions in the shell, 
   *                         for each shell, for each point. If requested there will be appended 
   *                         f/dx(ixyz,iSh,ipt), f/dy(ixyz,iSh,ipt), f/dz(ixyz,iSh,ipt) values of 
   *                         the functions in the shell, for each shell(nShSize), for each point(npts)
   *  \param [in] IOff      OffSet to properly store the basis set and Gradient components 
   *  \param [in] IOffSCR   OffSet to properly read the basis set in cartesian and Gradient components.. 
   *  \param [in] forceCart  True if force cartesian, otherwise sperical evaluation
   */ 
  void CarToSpDEval(SHELL_EVAL_TYPE typ, size_t L, double *fCarEVal, double *FSpEVAl, size_t IOff, size_t IOffSCR, 
    bool forceCart){

    auto shSize_sp  = (2*L+1);
    auto shSize_car   = ((L+1)*(L+2))/2; 
    double * f_sp     = FSpEVAl ;
    double * dx_sp = f_sp   + IOff;
    double * dy_sp = dx_sp  + IOff;
    double * dz_sp = dy_sp  + IOff;
    double * f_car     = fCarEVal ;
    double * dx_car = f_car   + IOffSCR;
    double * dy_car = dx_car  + IOffSCR;
    double * dz_car = dy_car  + IOffSCR;
    double tmp, tmpx, tmpy, tmpz ;
    // No trasformation needed
    //FIXME if (L < 2 or force cart) 
    //if (L < 2 ){
    //bool forceCart = true ;
    if (L < 2 or forceCart){

      if (typ != GRADIENT) {
        for( auto I = 0u ; I<shSize_car ; I++){ 
          f_sp[I] = f_car[I];
          } //loop over car/sp (equal in this case)

      } else {
        for( auto I = 0u ; I<shSize_car ; I++) {
          f_sp[I]  =  f_car[I];
          dx_sp[I] = dx_car[I];
          dy_sp[I] = dy_car[I];
          dz_sp[I] = dz_car[I];
          } //loop over car/sp (equal in this case)

      } //GGA or not

    //We do transform here
    } else {

      if (typ != GRADIENT) {
        for( auto I = 0u ; I<shSize_sp ;I++  ) {
          tmp = 0.0;
          for( auto p = 0u; p<shSize_car ; p++) { 
            tmp += car2sph_matrix[L][I*shSize_car+p] * f_car[p];
            } //loop over cart
          f_sp[I] = tmp;
        } //loop over sp

      } else {
        for( auto I = 0u ; I<shSize_sp ;I++  ) {
          tmp  = 0.0;
          tmpx = 0.0;
          tmpy = 0.0;
          tmpz = 0.0;
          for( auto p = 0u; p<shSize_car ; p++) { 
            tmp  += car2sph_matrix[L][I*shSize_car+p] * f_car[p];
            tmpx += car2sph_matrix[L][I*shSize_car+p] * dx_car[p];
            tmpy += car2sph_matrix[L][I*shSize_car+p] * dy_car[p];
            tmpz += car2sph_matrix[L][I*shSize_car+p] * dz_car[p];
            } //loop over cart
          f_sp[I]  =  tmp;
          dx_sp[I] = tmpx;
          dy_sp[I] = tmpy;
          dz_sp[I] = tmpz;

#if Basis_DEBUG_LEVEL >= 3
          // Debug Printing
          std::cerr << I << 
            "  f(pt) "<< f_sp[I] <<" " <<std::endl;
          std::cerr << I<<" " <<  
            " dx(pt) "<< dx_sp[I] << std::endl;
          std::cerr << I<<" " << 
            " dy(pt) "<< dy_sp[I] << std::endl;
          std::cerr << I<<" " <<  
            " dz(pt) "<< dz_sp[I] << std::endl;
#endif

        } //loop over sp

      } //GGA or not

    } //copy vs transform
  }; // CarToSpDEval

// SS start

  void CarToSpDEval(SHELL_EVAL_TYPE typ, size_t L, dcomplex *fCarEVal, dcomplex *FSpEVAl, size_t IOff, size_t IOffSCR, 
    bool forceCart){

    auto shSize_sp  = (2*L+1);
    auto shSize_car   = ((L+1)*(L+2))/2; 
    dcomplex * f_sp     = FSpEVAl ;
    dcomplex * dx_sp = f_sp   + IOff;
    dcomplex * dy_sp = dx_sp  + IOff;
    dcomplex * dz_sp = dy_sp  + IOff;
    dcomplex * f_car     = fCarEVal ;
    dcomplex * dx_car = f_car   + IOffSCR;
    dcomplex * dy_car = dx_car  + IOffSCR;
    dcomplex * dz_car = dy_car  + IOffSCR;
    dcomplex tmp, tmpx, tmpy, tmpz ;
    // No trasformation needed
    //FIXME if (L < 2 or force cart) 
    //if (L < 2 ){
    //bool forceCart = true ;
    if (L < 2 or forceCart){

      if (typ != GRADIENT) {
        for( auto I = 0u ; I<shSize_car ; I++){ 
          f_sp[I] = f_car[I];
          } //loop over car/sp (equal in this case)

      } else {
        for( auto I = 0u ; I<shSize_car ; I++) {
          f_sp[I]  =  f_car[I];
          dx_sp[I] = dx_car[I];
          dy_sp[I] = dy_car[I];
          dz_sp[I] = dz_car[I];
          } //loop over car/sp (equal in this case)

      } //GGA or not

    //We do transform here
    } else {

      if (typ != GRADIENT) {
        for( auto I = 0u ; I<shSize_sp ;I++  ) {
          tmp = 0.0;
          for( auto p = 0u; p<shSize_car ; p++) { 
            tmp += car2sph_matrix[L][I*shSize_car+p] * f_car[p];
            } //loop over cart
          f_sp[I] = tmp;
        } //loop over sp

      } else {
        for( auto I = 0u ; I<shSize_sp ;I++  ) {
          tmp  = 0.0;
          tmpx = 0.0;
          tmpy = 0.0;
          tmpz = 0.0;
          for( auto p = 0u; p<shSize_car ; p++) { 
            tmp  += car2sph_matrix[L][I*shSize_car+p] * f_car[p];
            tmpx += car2sph_matrix[L][I*shSize_car+p] * dx_car[p];
            tmpy += car2sph_matrix[L][I*shSize_car+p] * dy_car[p];
            tmpz += car2sph_matrix[L][I*shSize_car+p] * dz_car[p];
            } //loop over cart
          f_sp[I]  =  tmp;
          dx_sp[I] = tmpx;
          dy_sp[I] = tmpy;
          dz_sp[I] = tmpz;

#if Basis_DEBUG_LEVEL >= 3
/*
          // Debug Printing
          std::cerr << I << 
            "  f(pt) "<< f_sp[I] <<" " <<std::endl;
          std::cerr << I<<" " <<  
            " dx(pt) "<< dx_sp[I] << std::endl;
          std::cerr << I<<" " << 
            " dy(pt) "<< dy_sp[I] << std::endl;
          std::cerr << I<<" " <<  
            " dz(pt) "<< dz_sp[I] << std::endl;
*/
#endif

        } //loop over sp

      } //GGA or not

    } //copy vs transform
  }; // CarToSpDEval
// SS end

  /**
   *  \brief only for debug 
   *  it tests 5 pts and the evaluation of the f,dx,dy,dz at those points for a vector of shells
   *
   */ 
  void testEval(double *SCR, std::vector<libint2::Shell> &vshells, bool forceCart){
    std::vector<std::array<double,3>>  testpts;
    testpts.push_back({0,0,0});
    testpts.push_back({0.1,0,0});
    testpts.push_back({0,0.1,0});
    testpts.push_back({0,0.,0.1});
    testpts.push_back({1.,0.5,0.1});
    size_t npts = testpts.size();
    std::cout <<"inside testEval" <<std::endl;
    evalShellSet(GRADIENT,vshells,&testpts[0][0],npts,SCR,forceCart);
  }; // testEval

  /**
   *  \brief Level 2 Basis Set Gradient Evaluation Function - Used in the KS - DFT
   *  Evaluates a shell set over a specified number of cartesian points. This function requires a precomputed
   *  set of the distances and their x,y,z component for each point from each shell origin in the shells vector.
   *  \param [in] typ        Type of evaluation to perform (gradient, etc)
   *  \param [in] shells     Shell set for evaluation(vector of libint2::Shell).
   *  \param [in] evalshells Vector of bool to know if that Shell is relevant  for evaluation.
   *  \param [in] rSq        Raw storage of overall squared distances between each point and the shell origin (precomputed outside)
   *                         Dimension (nshells * npts)
   *  \param [in] r          Raw storage of overall x,y,z of the vector between each point and the shell origin 
   *                         (precomputed outside). Dimension (3*nshells*npts) 
   *  \param [in] npts       Number of cartesian points to evaluate.
   *  \param [in] nCenter    Number of total distint shell center (is NAtoms).
   *  \param [in] NBasisEff  Number of basis function to be evaluated given all the shells in input.
   *  \param [in] mapSh2Cen  Shell Mapping to atom centers.
   *  \param [in/out] fEval  eval Storage for the shell set evaluation, f(ixyz,iSh,ipt). 
   *                         Variable dimensions(npts*NBasisEff*1 or 4) allocated outside.
   *                         This storage will have all values of the functions in the shell, 
   *                         for each shell, for each point. If requested there will be appended 
   *                         f/dx(ixyz,iSh,ipt), f/dy(ixyz,iSh,ipt), f/dz(ixyz,iSh,ipt) values of 
   *                         the functions in the shell, for each shell(nShSize), for each point(npts)
   *  \param [in] SCR        eval Storage for the shell set evaluation, Cartesian f(MaxShellSize). 
   *  \param [in] IOffSCR    offset for eval Storage for the shell set evaluation, Cartesian f(MaxShellSize). 
   *  \param [in] forceCart  True if force cartesian, otherwise sperical evaluation
   */ 
  void evalShellSetGrad(SHELL_EVAL_TYPE typ, std::vector<libint2::Shell> &shells, 
    std::vector<bool> &evalShell, double* rSq, double *r, size_t npts, size_t nCenter, 
    std::vector<size_t> &mapSh2Cen, size_t NBasisEff, double *fEval, double *SCR, size_t IOffSCR, bool forceCart) {

    assert(shells.size() == evalShell.size());

    size_t nShSize = shells.size();
    size_t IOff =  npts*NBasisEff;
    std::array<double,3> rVal;

    for (auto ipts = 0ul; ipts < npts; ipts++){
      size_t Ic = 0;
    for (auto iSh = 0ul; iSh < nShSize; iSh++){
      if(evalShell[iSh]) {
        double * fStart    = fEval + Ic + ipts*NBasisEff;

        rVal [0]= r[0 + mapSh2Cen[iSh]*3 + ipts*3*nCenter];
        rVal [1]= r[1 + mapSh2Cen[iSh]*3 + ipts*3*nCenter];
        rVal [2]= r[2 + mapSh2Cen[iSh]*3 + ipts*3*nCenter];

        evalShellSetGrad(typ,shells[iSh],rSq[mapSh2Cen[iSh] + ipts*nCenter],
          rVal,SCR,IOffSCR); 

        CarToSpDGradEval(typ, shells[iSh].contr[0].l, SCR, fStart, IOff, IOffSCR, forceCart);

        Ic += shells[iSh].size(); // Increment offset in basis
      }

    } // loop over shells
    } // loop over points

    //std::cout << "Basis Gradient in evalShellSetGrad" << std::endl;
    //for (size_t i = 0; i < IOff; i++)
    //  std::cout << (fEval+IOff)[i] << std::endl;

  }; // evalShellSetGrad Level 2


  /**
   *   \brief Level 3 Basis Set Nuclear Gradient Evaluation Function
   *   Evaluates a single shell over a single cartesian point. This function requires a precomputed
   *   the distance and its x,y,z components for the point from the shell origin. An offset
   *   to properly store the results can be used.
   *   \param [in] typ        Type of evaluation to perform (gradient, etc)
   *   \param [in] shell      Shell for evaluation(libint2::Shell).
   *   \param [in] rSq        Raw storage of square distance between the point and the shell origin (precomputed outside)
   *   \param [in] r          Raw storage of x,y,z of the vector between each point and the shell origin 
   *                          (precomputed outside).  
   *   \param [in/out] SCR    eval Storage for the shell set evaluation, f(ixyz,iSh,ipt). 
   *                          Variable dimensions(npts*NBasisEff*1 or 4) allocated outside.
   *                          This storage will have all values of the functions in the shell, 
   *                          for each shell, for each point. If requested there will be appended 
   *                          f/dx(ixyz,iSh,ipt), f/dy(ixyz,iSh,ipt), f/dz(ixyz,iSh,ipt) values of 
   *                          the functions in the shell, for each shell(nShSize), for each point(npts)
   *   \param [in] IOffSCR       OffSet to properly store the basis set. 
   */ 
  void evalShellSetGrad(SHELL_EVAL_TYPE typ, const libint2::Shell &shell,double rSq, const std::array<double,3> &xyz, 
    double* SCR, size_t IOffSCR) {
    auto L         = shell.contr[0].l;
    auto shSize    = ((L+1)*(L+2))/2; 
    auto shSize_car   = ((L+1)*(L+2))/2; 

    //std::cout << "x " << xyz[0] << std::endl;
    //std::cout << "y " << xyz[1] << std::endl;
    //std::cout << "z " << xyz[2] << std::endl;

    // pointer positions for nuclear gradients of density 
    double * dX_car     = SCR ;
    double * dY_car = dX_car  + IOffSCR;
    double * dZ_car = dY_car  + IOffSCR;

    // pointer positions for nuclear gradients of density gradient
    double * dxX_car = dZ_car + IOffSCR;
    double * dxY_car = dxX_car + IOffSCR;
    double * dxZ_car = dxY_car + IOffSCR;
    double * dyX_car = dxZ_car + IOffSCR;
    double * dyY_car = dyX_car + IOffSCR;
    double * dyZ_car = dyY_car + IOffSCR;
    double * dzX_car = dyZ_car + IOffSCR;
    double * dzY_car = dzX_car + IOffSCR;
    double * dzZ_car = dzY_car + IOffSCR;

    // contraction length
    auto contDepth = shell.alpha.size(); 
    double alpha(0.0);
    double alpha_s(0.0);
    double expFactor(0.0);
    double expArg(0);
    double tmpcoef,tmpalpha;
    int lx,ly,lz, ixyz;
    double tmpxyz;
    double tmpdx;
    double tmpdy;
    double tmpdz;
    double tmpdxx;
    double tmpdxy;
    double tmpdxz;
    double tmpdyx;
    double tmpdyy;
    double tmpdyz;
    double tmpdzx;
    double tmpdzy;
    double tmpdzz;

    double shift = 1e-8;
    double rSq_x_shift = xyz[0]*xyz[0] + xyz[1]*xyz[1] + (xyz[2]-shift)*(xyz[2]-shift); 
    double expFactor_x_shift(0.0);
    double alpha_x_shift(0.0);
    // Generating the expArgument, expFactotr and the
    // alpha (for derivatives later on) and store them
    // in temp variables
    for(auto k = 0; k < contDepth; k++){
      tmpcoef = shell.contr[0].coeff[k];
      tmpalpha = shell.alpha[k];
      expArg = std::exp(-tmpalpha*rSq);
      expFactor += tmpcoef * expArg;
      expFactor_x_shift += tmpcoef * std::exp(-tmpalpha*rSq_x_shift);
      tmpcoef *= tmpalpha;
      alpha += tmpcoef * expArg;
      alpha_x_shift += tmpcoef * std::exp(-tmpalpha*rSq_x_shift);
      if (typ == GRADIENT) { 
        // quantities for derivatives
        alpha_s += tmpalpha * tmpcoef * expArg;
      }
    } 

    alpha *= 2;
    alpha_x_shift *= 2;
    if (typ == GRADIENT) alpha_s *= 4;

    // sign
    double sign = -1;

    // lambda function that computes x^lx * y^ly * z^lz
    auto compute_power = [&](int x, int y, int z) -> double {
      
      // return zero if any power is less than zero 
      if ( x < 0 or y < 0 or z < 0 ) return 0;

      double result = 1.0;
      for (size_t xi = 0; xi < x; xi++) result *= double(xyz[0]);
      for (size_t yi = 0; yi < y; yi++) result *= double(xyz[1]);
      for (size_t zi = 0; zi < z; zi++) result *= double(xyz[2]);

      return result;

    }; 

    for(auto i = 0u, I = 0u; i <= L; i++) {
      lx = L - i;
      for( auto j = 0u; j <= i; j++, I++) {
        ly = i - j;
        lz = L - lx - ly;

        //std::cout << "lx " << lx << " ly " << ly << " lz " << lz << std::endl;
        //std::cout << "expFactor  " << expFactor << std::endl;
        //std::cout << "alpha " << alpha << std::endl;
        //std::cout << "alpha_check " << alpha_check << std::endl;
        //std::cout << "alpha_s " << alpha_s << std::endl;
        //std::cout << "x " << xyz[0] << std::endl;
        //std::cout << "y " << xyz[1] << std::endl;
        //std::cout << "z " << xyz[2] << std::endl;
        tmpxyz= 1.0;
        tmpdx = 0.0;
        tmpdy = 0.0;
        tmpdz = 0.0;
        tmpdxx = 0.0;
        tmpdxy = 0.0;
        tmpdxz = 0.0;
        tmpdyx = 0.0;
        tmpdyy = 0.0;
        tmpdyz = 0.0;
        tmpdzx = 0.0;
        tmpdzy = 0.0;
        tmpdzz = 0.0;
        for(ixyz = 0; ixyz < lx-2; ixyz++) tmpxyz *= xyz[0];
        for(ixyz = 0; ixyz < ly-2; ixyz++) tmpxyz *= xyz[1];
        for(ixyz = 0; ixyz < lz-2; ixyz++) tmpxyz *= xyz[2];

        double check_xyz = 1.0;
        for(ixyz = 0; ixyz < lx-2; ixyz++) check_xyz *=  xyz[0];
        for(ixyz = 0; ixyz < ly-2; ixyz++) check_xyz *=  xyz[1];
        for(ixyz = 0; ixyz < lz-2; ixyz++) check_xyz *=  (xyz[2]-shift);
        
        double f_car = tmpxyz;
        double f_car_check = check_xyz;
        //f_car[I]  =  tmpxyz;

        if(lx> 0) {tmpdx = sign * expFactor * lx;}
        if(ly> 0) {tmpdy = sign * expFactor * ly;}
        if(lz> 0) {tmpdz = sign * expFactor * lz;}
         
        dX_car[I] = tmpxyz*tmpdx;
        dY_car[I] = tmpxyz*tmpdy;
        dZ_car[I] = tmpxyz*tmpdz;
    
        // finishing up        
        double tmpdx_check = 0.0;
        if(ly> 0) {tmpdx_check = sign * expFactor_x_shift * ly;}
        double dX_car_check = check_xyz*tmpdx_check;

        if(lx> 0) {f_car *= xyz[0]; f_car_check *= xyz[0]; dY_car[I] *=xyz[0];dZ_car[I] *=xyz[0];}
        if(ly> 0) {f_car *= xyz[1]; f_car_check *= xyz[1]; dX_car[I] *=xyz[1];dZ_car[I] *=xyz[1];}
        if(lz> 0) {f_car *= xyz[2]; f_car_check *=(xyz[2]-shift); dX_car[I] *=xyz[2];dY_car[I] *=xyz[2];}

        // debug code
        if(lx> 0) {dX_car_check *=xyz[0];}
        if(lz> 0) {dX_car_check *=xyz[2];}


        if(lx> 1) {f_car *= xyz[0]; f_car_check *= xyz[0]; dX_car[I] *=xyz[0];dY_car[I] *=xyz[0];dZ_car[I] *=xyz[0];}
        if(ly> 1) {f_car *= xyz[1]; f_car_check *= xyz[1]; dY_car[I] *=xyz[1];dX_car[I] *=xyz[1];dZ_car[I] *=xyz[1];}
        if(lz> 1) {f_car *= xyz[2]; f_car_check *=(xyz[2]-shift); dZ_car[I] *=xyz[2];dX_car[I] *=xyz[2];dY_car[I] *=xyz[2];}

        // debug code
        if(lx> 1) {dX_car_check *= xyz[0];}
        if(ly> 1) {dX_car_check *=xyz[1];}
        if(lz> 1) {dX_car_check *=(xyz[2]-shift);}

        dX_car[I] -= sign * f_car * xyz[0] * alpha;
        dY_car[I] -= sign * f_car * xyz[1] * alpha;
        dZ_car[I] -= sign * f_car * xyz[2] * alpha;

        dX_car_check -= sign * f_car_check * xyz[1] * alpha_x_shift;

        // numerical gradient check
        //double check_dx = (f_car_check * expFactor_x_shift - f_car * expFactor) / 1e-6;
        //dX_car[I] = check_dx;
        //std::cout << "unshifted dx " << dX_car[I] << std::endl;
        //std::cout << "shifted dx " << dX_car_check << std::endl;
        double check_dx = (-dX_car_check + dY_car[I]) / shift;
        //std::cout << "Numerical nuclear " << check_dx << std::endl;
        //std::cout << "Analytial nuclear" << dX_car[I] << std::endl;
        //std::cout << std::endl;

        if ( typ == GRADIENT ) {

          // diagonal terms
          dxX_car[I] = f_car * (lx + lx + 1) * alpha;
          dyY_car[I] = f_car * (ly + ly + 1) * alpha;
          dzZ_car[I] = f_car * (lz + lz + 1) * alpha;

          dxX_car[I] -= f_car * xyz[0] * xyz[0] * alpha_s;
          dyY_car[I] -= f_car * xyz[1] * xyz[1] * alpha_s;
          dzZ_car[I] -= f_car * xyz[2] * xyz[2] * alpha_s;

          if(lx> 1) {tmpdxx = -expFactor * lx * (lx-1);}
          if(ly> 1) {tmpdyy = -expFactor * ly * (ly-1);}
          if(lz> 1) {tmpdzz = -expFactor * lz * (lz-1);}

          tmpdxx *= tmpxyz;
          tmpdyy *= tmpxyz;
          tmpdzz *= tmpxyz;

          if(lx> 0) {tmpdyy *=xyz[0];tmpdzz *=xyz[0];}
          if(ly> 0) {tmpdxx *=xyz[1];tmpdzz *=xyz[1];}
          if(lz> 0) {tmpdxx *=xyz[2];tmpdyy *=xyz[2];}

          if(lx> 1) {tmpdyy *=xyz[0];tmpdzz *=xyz[0];}
          if(ly> 1) {tmpdxx *=xyz[1];tmpdzz *=xyz[1];}
          if(lz> 1) {tmpdxx *=xyz[2];tmpdyy *=xyz[2];}

          dxX_car[I] += tmpdxx;
          dyY_car[I] += tmpdyy;
          dzZ_car[I] += tmpdzz;


          // off-diagonal terms
          dxY_car[I] = -f_car * xyz[0] * xyz[1] * alpha_s;
          dxZ_car[I] = -f_car * xyz[0] * xyz[2] * alpha_s;
          dyX_car[I] = -f_car * xyz[1] * xyz[0] * alpha_s;
          dyZ_car[I] = -f_car * xyz[1] * xyz[2] * alpha_s;
          dzX_car[I] = -f_car * xyz[2] * xyz[0] * alpha_s;
          dzY_car[I] = -f_car * xyz[2] * xyz[1] * alpha_s;

          dxY_car[I] -= expFactor * lx * ly * compute_power(lx-1,ly-1,lz);
          dxZ_car[I] -= expFactor * lx * lz * compute_power(lx-1,ly,lz-1);
          dyX_car[I] -= expFactor * ly * lx * compute_power(lx-1,ly-1,lz);
          dyZ_car[I] -= expFactor * ly * lz * compute_power(lx,ly-1,lz-1);
          dzX_car[I] -= expFactor * lz * lx * compute_power(lx-1,ly,lz-1);
          dzY_car[I] -= expFactor * lz * ly * compute_power(lx,ly-1,lz-1);

          dxY_car[I] += alpha * lx * compute_power(lx-1,ly+1,lz);
          dxY_car[I] += alpha * ly * compute_power(lx+1,ly-1,lz);
          //std::cout << "test " << compute_power(lx+1,ly-1,lz) << std::endl;
          dxZ_car[I] += alpha * lx * compute_power(lx-1,ly,lz+1);
          dxZ_car[I] += alpha * lz * compute_power(lx+1,ly,lz-1);

          dyX_car[I] += alpha * ly * compute_power(lx+1,ly-1,lz);
          dyX_car[I] += alpha * lx * compute_power(lx-1,ly+1,lz);
          dyZ_car[I] += alpha * ly * compute_power(lx,ly-1,lz+1);
          dyZ_car[I] += alpha * lz * compute_power(lx,ly+1,lz-1);

          dzX_car[I] += alpha * lz * compute_power(lx+1,ly,lz-1);
          dzX_car[I] += alpha * lx * compute_power(lx-1,ly,lz+1);
          dzY_car[I] += alpha * ly * compute_power(lx,ly-1,lz+1);
          dzY_car[I] += alpha * lz * compute_power(lx,ly+1,lz-1);

          //std::cout << "Analytical nuclear " << dzY_car[I] << std::endl;
          //if (std::abs(dzY_car[I] - check_dx) > 1e-4) 
          //  std::cout << "Incorrect gradient" << std::endl;
          //std::cout << std::endl;


        }

        //std::cout << "dX " << dX_car[I] << "  " << "dY " << dY_car[I] << "  " << "dZ " << dZ_car[I] << std::endl;
    
#if Basis_DEBUG_LEVEL >= 3
          // Debug Printing
          //std::cerr << I <<" "<< lx << " " 
          //  << ly << " "<<lz <<"  f(pt) "<< f_car[I] <<" " <<std::endl;
          std::cerr << I<<" "<< lx << " "  
            << ly << " "<<lz <<" dX(pt) "<< dX_car[I] << std::endl;
          std::cerr << I<<" "<< lx << " " 
            << ly << " "<<lz <<" dY(pt) "<< dY_car[I] << std::endl;
          std::cerr << I<<" "<< lx << " " 
            << ly << " "<<lz <<" dZ(pt) "<< dZ_car[I] << std::endl;
#endif

      } //loop overj, j[0,i]
    } //loop over i, i[0,L] this to loop required to build the lx,ly,lz combination given L

  }; // evalShellSetGrad Level3


  /**
   *   \brief Basis Set Cartesian to Sperical conversion over a single shell.
   *   It needs to be called also for cartisian evaluation, since this function
   *   handles the final population of the eval Storage for the shell and points.
   *   \param [in] typ       Type of evaluation to perform (gradient, etc)
   *   \param [in] L         Angular quantum momentum
   *   \param [in] fCarEVal  eval Storage for the shell set evaluation, f(ixyz,iSh,ipt). 
   *                         Variable dimensions(npts*NBasisEff*1 or 4) allocated outside.
   *                         This storage will have all values of the functions in the shell, 
   *                         for each shell, for each point. If requested there will be appended 
   *                         f/dx(ixyz,iSh,ipt), f/dy(ixyz,iSh,ipt), f/dz(ixyz,iSh,ipt) values of 
   *                         the functions in the shell, for each shell(nShSize), for each point(npts)
   *  \param [in/out] fEval  eval Storage for the shell set evaluation, f(ixyz,iSh,ipt). 
   *                         Variable dimensions(npts*NBasisEff*1 or 4) allocated outside.
   *                         This storage will have all values of the functions in the shell, 
   *                         for each shell, for each point. If requested there will be appended 
   *                         f/dx(ixyz,iSh,ipt), f/dy(ixyz,iSh,ipt), f/dz(ixyz,iSh,ipt) values of 
   *                         the functions in the shell, for each shell(nShSize), for each point(npts)
   *  \param [in] IOff      OffSet to properly store the basis set and Gradient components 
   *  \param [in] IOffSCR   OffSet to properly read the basis set in cartesian and Gradient components.. 
   *  \param [in] forceCart  True if force cartesian, otherwise sperical evaluation
   */ 
  void CarToSpDGradEval(SHELL_EVAL_TYPE typ, size_t L, double *fCarEVal, double *FSpEVAl, size_t IOff, size_t IOffSCR, 
    bool forceCart){

    auto shSize_sp  = (2*L+1);
    auto shSize_car   = ((L+1)*(L+2))/2; 

    double * dX_sp = FSpEVAl;
    double * dY_sp = dX_sp + IOff;
    double * dZ_sp = dY_sp + IOff;
    double * dxX_sp = dZ_sp + IOff;
    double * dxY_sp = dxX_sp + IOff;
    double * dxZ_sp = dxY_sp + IOff;
    double * dyX_sp = dxZ_sp + IOff;
    double * dyY_sp = dyX_sp + IOff;
    double * dyZ_sp = dyY_sp + IOff;
    double * dzX_sp = dyZ_sp + IOff;
    double * dzY_sp = dzX_sp + IOff;
    double * dzZ_sp = dzY_sp + IOff;

    double * dX_car = fCarEVal;
    double * dY_car = dX_car + IOffSCR;
    double * dZ_car = dY_car + IOffSCR;
    double * dxX_car = dZ_car + IOffSCR;
    double * dxY_car = dxX_car + IOffSCR;
    double * dxZ_car = dxY_car + IOffSCR;
    double * dyX_car = dxZ_car + IOffSCR;
    double * dyY_car = dyX_car + IOffSCR;
    double * dyZ_car = dyY_car + IOffSCR;
    double * dzX_car = dyZ_car + IOffSCR;
    double * dzY_car = dzX_car + IOffSCR;
    double * dzZ_car = dzY_car + IOffSCR;

    double tmp, tmpx, tmpy, tmpz ;
    double tmp_dX, tmp_dY, tmp_dZ;
    double tmp_dxX, tmp_dxY, tmp_dxZ;
    double tmp_dyX, tmp_dyY, tmp_dyZ;
    double tmp_dzX, tmp_dzY, tmp_dzZ;
    // No trasformation needed
    //FIXME if (L < 2 or force cart) 
    //bool forceCart = true ;
    if (L < 2 or forceCart){

      if (typ != GRADIENT) {
        for( auto I = 0u ; I<shSize_car ; I++){ 
          dX_sp[I] = dX_car[I];
          dY_sp[I] = dY_car[I];
          dZ_sp[I] = dZ_car[I];
          } //loop over car/sp (equal in this case)

      } else {
        for( auto I = 0u ; I<shSize_car ; I++) {
          dX_sp[I] = dX_car[I];
          dY_sp[I] = dY_car[I];
          dZ_sp[I] = dZ_car[I];
          dxX_sp[I] = dxX_car[I];
          dxY_sp[I] = dxY_car[I];
          dxZ_sp[I] = dxZ_car[I];
          dyX_sp[I] = dyX_car[I];
          dyY_sp[I] = dyY_car[I];
          dyZ_sp[I] = dyZ_car[I];
          dzX_sp[I] = dzX_car[I];
          dzY_sp[I] = dzY_car[I];
          dzZ_sp[I] = dzZ_car[I];
          } //loop over car/sp (equal in this case)

      } //GGA or not

    //We do transform here
    } else {

      if (typ != GRADIENT) {
        for( auto I = 0u ; I<shSize_sp ;I++  ) {
          tmp_dX = 0.0;
          tmp_dY = 0.0;
          tmp_dZ = 0.0;
          for( auto p = 0u; p<shSize_car ; p++) { 
            tmp_dX += car2sph_matrix[L][I*shSize_car+p] * dX_car[p];
            tmp_dY += car2sph_matrix[L][I*shSize_car+p] * dY_car[p];
            tmp_dZ += car2sph_matrix[L][I*shSize_car+p] * dZ_car[p];
            } //loop over cart
          dX_sp[I] = tmp_dX;
          dY_sp[I] = tmp_dY;
          dZ_sp[I] = tmp_dZ;
        } //loop over sp

      } else {
        for( auto I = 0u ; I<shSize_sp ;I++  ) {
          tmp_dX  = 0.0;
          tmp_dY  = 0.0;
          tmp_dZ  = 0.0;
          tmp_dxX = 0.0;
          tmp_dxY = 0.0;
          tmp_dxZ = 0.0;
          tmp_dyX = 0.0;
          tmp_dyY = 0.0;
          tmp_dyZ = 0.0;
          tmp_dzX = 0.0;
          tmp_dzY = 0.0;
          tmp_dzZ = 0.0;
          for( auto p = 0u; p<shSize_car ; p++) { 
            tmp_dX  += car2sph_matrix[L][I*shSize_car+p] * dX_car[p];
            tmp_dY  += car2sph_matrix[L][I*shSize_car+p] * dY_car[p];
            tmp_dZ  += car2sph_matrix[L][I*shSize_car+p] * dZ_car[p];
            tmp_dxX += car2sph_matrix[L][I*shSize_car+p] * dxX_car[p];
            tmp_dxY += car2sph_matrix[L][I*shSize_car+p] * dxY_car[p];
            tmp_dxZ += car2sph_matrix[L][I*shSize_car+p] * dxZ_car[p];
            tmp_dyX += car2sph_matrix[L][I*shSize_car+p] * dyX_car[p];
            tmp_dyY += car2sph_matrix[L][I*shSize_car+p] * dyY_car[p];
            tmp_dyZ += car2sph_matrix[L][I*shSize_car+p] * dyZ_car[p];
            tmp_dzX += car2sph_matrix[L][I*shSize_car+p] * dzX_car[p];
            tmp_dzY += car2sph_matrix[L][I*shSize_car+p] * dzY_car[p];
            tmp_dzZ += car2sph_matrix[L][I*shSize_car+p] * dzZ_car[p];
            } //loop over cart
          dX_sp[I] = tmp_dX;
          dY_sp[I] = tmp_dY;
          dZ_sp[I] = tmp_dZ;
          dxX_sp[I] = tmp_dxX;
          dxY_sp[I] = tmp_dxY;
          dxZ_sp[I] = tmp_dxZ;
          dyX_sp[I] = tmp_dyX;
          dyY_sp[I] = tmp_dyY;
          dyZ_sp[I] = tmp_dyZ;
          dzX_sp[I] = tmp_dzX;
          dzY_sp[I] = tmp_dzY;
          dzZ_sp[I] = tmp_dzZ;

#if Basis_DEBUG_LEVEL >= 3
          // Debug Printing
          std::cerr << I<<" " <<  
            " dx(pt) "<< dX_sp[I] << std::endl;
          std::cerr << I<<" " << 
            " dy(pt) "<< dY_sp[I] << std::endl;
          std::cerr << I<<" " <<  
            " dz(pt) "<< dZ_sp[I] << std::endl;
#endif

        } //loop over sp

      } //GGA or not

    } //copy vs transform
  }; // CarToSpDEval

}; // namespace ChronusQ
