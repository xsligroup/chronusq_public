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

#include <libint2/shell.h>
#include <basisset/basisset_def.hpp>
#include <fields.hpp> 

namespace ChronusQ {


  /**
   *  \brief Specification for shellset evaluation.
   *
   *  Currently supports GRADIENT and NOGRAD (simple eval) evaluation
   *  types.
   */ 
  enum SHELL_EVAL_TYPE {
    GRADIENT,
    NOGRAD
  };
 
  /**
   *  \brief Level 1 Basis Set Evaluation Function
   *  \brief Evaluates a shell set over a specified number of cartesian points.
   */ 
  void evalShellSet(SHELL_EVAL_TYPE, std::vector<libint2::Shell> &, double *, size_t, double *, bool, double * SCR = nullptr);

  inline void evalShellSet(SHELL_EVAL_TYPE typ, std::vector<libint2::Shell> &shells,
    std::vector<std::array<double,3>> &pts, double *eval, bool &forceCart){

    evalShellSet(typ,shells,&(pts[0][0]),pts.size(),eval,forceCart);

  }; // evalShellSet (over vector of arrays)

  /**
   *  \brief Level 2 Basis Set Evaluation Function.
   *  \brief Evaluates a shell set over a specified number of cartesian points. This 
   *  \brief function requires a precomputed set of distances and their x,y,z components 
   *  \brief for each point from each shell origin in the shells vector..
   */ 
  void evalShellSet(SHELL_EVAL_TYPE, std::vector<libint2::Shell> &, std::vector<bool> &,double *, double *, size_t, 
    size_t, std::vector<size_t>&, size_t, double*, double*, size_t, bool );
// SS start
  // define level 2 basis set eval for GIAO 
  void evalShellSet(SHELL_EVAL_TYPE, std::vector<libint2::Shell> &, std::vector<bool> &,double *, double *, size_t, 
    size_t, std::vector<size_t>&, size_t, dcomplex*, dcomplex*, size_t, bool, EMPerturbation& );
// SS end 

  /**
   *  \brief Level 3 Basis Set Evaluation Function
   *  \brief Evaluates a single shell over a single cartesian point. This function requires a precomputed
   *  \brief distance and its x,y,z components for the point from the shell origin. An offset
   *  \brief to properly store the results can be used..
   */ 
  void evalShellSet(SHELL_EVAL_TYPE,const libint2::Shell&,double,const std::array<double,3>&, double *, size_t);
//SS start
  // define level 3 basis set eval for GIAO  
  void evalShellSet(SHELL_EVAL_TYPE,const libint2::Shell&,double,const std::array<double,3>&, dcomplex *, size_t, double *);
//SS end

  /**
   *  \brief Level 2 Basis Set Gradient Evaluation Function
   *  \brief Evaluates a shell set over a specified number of cartesian points. This
   *  \brief function requires a precomputed set of distances and their x,y,z components
   *  \brief for each point from each shell origin in the shells vector..
   */
  void evalShellSetGrad(SHELL_EVAL_TYPE, std::vector<libint2::Shell> &, std::vector<bool> &, double *, double *, size_t, 
    size_t, std::vector<size_t> &, size_t, double*, double*, size_t, bool );

  /**
   *  \brief Level 3 Basis Set Gradient Evaluation Function
   *  \brief Evaluates a single shell over a single cartesian point. This function requires a precomputed
   *  \brief distance and its x,y,z components for the point from the shell origin. An offset
   *  \brief to properly store the results can be used..
   */
  void evalShellSetGrad(SHELL_EVAL_TYPE,const libint2::Shell&,double,const std::array<double,3>&, double*, size_t);

  /**
   *  \brief Basis Set transformation from Cartesian to Spherical
   *  it also provide to copy to the final storage (pointer) and requires that the transformation matrix Sp <-> Cart
   *  is already populated.
   */ 
  void CarToSpDEval(SHELL_EVAL_TYPE, size_t , double *, double*, size_t, size_t, bool);
  void CarToSpDEval(SHELL_EVAL_TYPE, size_t , dcomplex *, dcomplex*, size_t, size_t, bool); //GIAO

  /**
   *  \brief Basis Set Gradient transformation from Cartesian to Spherical
   *  it also provide to copy to the final storage (pointer) and requires that the transformation matrix Sp <-> Cart
   *  is already populated.
   */
  void CarToSpDGradEval(SHELL_EVAL_TYPE, size_t, double *, double *, size_t, size_t, bool);

  void testEval(double *, std::vector<libint2::Shell> &, bool);





















}; // namespace ChronusQ

