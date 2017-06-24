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

#include <particleintegrals/twopints/gtodirectreleri.hpp>
#include <particleintegrals/contract/direct.hpp>
#include <libcint.hpp>

namespace ChronusQ {


  /************************************/
  /* Libcint 4C-direct Implementation */
  /************************************/
 
  template <typename MatsT, typename IntsT>
  void GTODirectRelERIContraction<MatsT,IntsT>::directScaffoldLibcint(
    MPI_Comm comm, const bool screen,
    std::vector<TwoBodyContraction<MatsT>> &matList,
    const bool computeExchange,
    const APPROXIMATION_TYPE_4C approximate4C) const {

    CErr("Please install CQDirect4C library to run direct four-component implementation.");

  }

  template <>
  void GTODirectRelERIContraction<double,double>::directScaffoldLibcint(
    MPI_Comm comm, const bool screen,
    std::vector<TwoBodyContraction<double>> &matList,
    const bool computeExchange,
    const APPROXIMATION_TYPE_4C approximate4C) const {
    CErr("Dirac-Coulomb + Real is an invalid option",std::cout);  
  }

  template <>
  void GTODirectRelERIContraction<dcomplex,dcomplex>::directScaffoldLibcint(
    MPI_Comm comm, const bool screen,
    std::vector<TwoBodyContraction<dcomplex>> &matList,
    const bool computeExchange,
    const APPROXIMATION_TYPE_4C approximate4C) const {
    CErr("Complex integral is is an invalid option",std::cout);  
  }


  template <typename MatsT, typename IntsT>
  void GTODirectRelERIContraction<MatsT,IntsT>::directRelScaffoldLibcintCoulombOnly(
      MPI_Comm comm, const bool screen,
      std::vector<TwoBodyRelContraction<MatsT>> &matList,
      const APPROXIMATION_TYPE_4C approximate4C) const {

    CErr("Please install CQDirect4C library to run direct four-component implementation.");

  }; // GTODirectRelERIContraction::directRelScaffoldLibcintCoulombOnly

  template <>
  void GTODirectRelERIContraction<double,double>::directRelScaffoldLibcintCoulombOnly(
      MPI_Comm comm, const bool screen,
      std::vector<TwoBodyRelContraction<double>> &matList,
      const APPROXIMATION_TYPE_4C approximate4C) const {
    CErr("Dirac-Coulomb + Real is an invalid option",std::cout);
  }

  template <>
  void GTODirectRelERIContraction<dcomplex,dcomplex>::directRelScaffoldLibcintCoulombOnly(
      MPI_Comm comm, const bool screen,
      std::vector<TwoBodyRelContraction<dcomplex>> &matList,
      const APPROXIMATION_TYPE_4C approximate4C) const {
    CErr("Complex integral is is an invalid option",std::cout);
  }

  template <typename MatsT, typename IntsT>
  size_t GTODirectRelERIContraction<MatsT,IntsT>::libcintCacheSize(
      const TWOBODY_CONTRACTION_TYPE & contType, int * atm, const int nAtoms,
      int * bas, const int nShells, double * env) const {

    CErr("Please install CQDirect4C library to run direct four-component implementation.");
    return 0;
  }

  /*******************************************************************************/
  /*                                                                             */
  /* Compute memory requirement for Libcint Batch 4C-direct                      */
  /* Returns:                                                                    */
  /*   size_t SCR size needed for one batch                                      */
  /*   IMPORTANT HERE: size are all in MatsT(dcomplex)                           */
  /*******************************************************************************/

  template <typename MatsT, typename IntsT>
  size_t GTODirectRelERIContraction<MatsT,IntsT>::directRelScaffoldLibcintSCRSize(
      const TWOBODY_CONTRACTION_TYPE & contType, const bool computeExchange) const {

    CErr("Please install CQDirect4C library to run direct four-component implementation.");
    return 0;
  }

  template <>
  size_t GTODirectRelERIContraction<double,double>::directRelScaffoldLibcintSCRSize(
      const TWOBODY_CONTRACTION_TYPE & contType, const bool computeExchange) const {
    CErr("Dirac-Coulomb + Real is an invalid option",std::cout);
    abort();
  }

  template <>
  size_t GTODirectRelERIContraction<dcomplex,dcomplex>::directRelScaffoldLibcintSCRSize(
      const TWOBODY_CONTRACTION_TYPE & contType, const bool computeExchange) const {
    CErr("Complex  is an invalid option",std::cout);
    abort();
  }


  template <typename MatsT, typename IntsT>
  void GTODirectRelERIContraction<MatsT,IntsT>::directRelScaffoldLibcintCoulombOnlySpinFree(
      MPI_Comm comm, const bool screen,
      std::vector<TwoBodyRelContraction<MatsT>> &matList,
      const APPROXIMATION_TYPE_4C approximate4C) const {

    CErr("Please install CQDirect4C library to run direct four-component implementation.");

  }; // GTODirectRelERIContraction::directRelScaffoldLibcintCoulombOnly

  template <>
  void GTODirectRelERIContraction<double,double>::directRelScaffoldLibcintCoulombOnlySpinFree(
      MPI_Comm comm, const bool screen,
      std::vector<TwoBodyRelContraction<double>> &matList,
      const APPROXIMATION_TYPE_4C approximate4C) const {
    CErr("Dirac-Coulomb + Real is an invalid option",std::cout);
  }

  template <>
  void GTODirectRelERIContraction<dcomplex,dcomplex>::directRelScaffoldLibcintCoulombOnlySpinFree(
      MPI_Comm comm, const bool screen,
      std::vector<TwoBodyRelContraction<dcomplex>> &matList,
      const APPROXIMATION_TYPE_4C approximate4C) const {
    CErr("Complex integral is is an invalid option",std::cout);
  }


  template <typename MatsT, typename IntsT>
  void GTODirectRelERIContraction<MatsT,IntsT>::directScaffoldLibcintSpinDependent(
      MPI_Comm comm, const bool screen,
      std::vector<TwoBodyContraction<MatsT>> &matList,
      const bool computeExchange,
      const APPROXIMATION_TYPE_4C approximate4C) const {

    CErr("Please install CQDirect4C library to run direct four-component implementation.");

  }

  template <>
  void GTODirectRelERIContraction<double,double>::directScaffoldLibcintSpinDependent(
      MPI_Comm comm, const bool screen,
      std::vector<TwoBodyContraction<double>> &matList,
      const bool computeExchange,
      const APPROXIMATION_TYPE_4C approximate4C) const {
    CErr("Dirac-Coulomb + Real is an invalid option",std::cout);
  }

  template <>
  void GTODirectRelERIContraction<dcomplex,dcomplex>::directScaffoldLibcintSpinDependent(
      MPI_Comm comm, const bool screen,
      std::vector<TwoBodyContraction<dcomplex>> &matList,
      const bool computeExchange,
      const APPROXIMATION_TYPE_4C approximate4C) const {
    CErr("Complex integral is is an invalid option",std::cout);
  }


  template <typename MatsT, typename IntsT>
  void GTODirectRelERIContraction<MatsT,IntsT>::directScaffoldLibcintSpinFree(
      MPI_Comm comm, const bool screen,
      std::vector<TwoBodyContraction<MatsT>> &matList,
      const bool computeExchange,
      const APPROXIMATION_TYPE_4C approximate4C) const {

    CErr("Please install CQDirect4C library to run direct four-component implementation.");

  }

  template <>
  void GTODirectRelERIContraction<double,double>::directScaffoldLibcintSpinFree(
      MPI_Comm comm, const bool screen,
      std::vector<TwoBodyContraction<double>> &matList,
      const bool computeExchange,
      const APPROXIMATION_TYPE_4C approximate4C) const {
    CErr("Dirac-Coulomb + Real is an invalid option",std::cout);
  }

  template <>
  void GTODirectRelERIContraction<dcomplex,dcomplex>::directScaffoldLibcintSpinFree(
      MPI_Comm comm, const bool screen,
      std::vector<TwoBodyContraction<dcomplex>> &matList,
      const bool computeExchange,
      const APPROXIMATION_TYPE_4C approximate4C) const {
    CErr("Complex integral is is an invalid option",std::cout);
  }


  template <typename MatsT, typename IntsT>
  void GTODirectRelERIContraction<MatsT,IntsT>::directScaffold(
      MPI_Comm comm, const bool screen,
      std::vector<TwoBodyContraction<MatsT>> &matList) const {

    CErr("Please install CQDirect4C library to run direct four-component implementation.");

  }

  template <>
  void GTODirectRelERIContraction<double,double>::directScaffold(
      MPI_Comm comm, const bool screen,
      std::vector<TwoBodyContraction<double>> &matList) const {
    CErr("Dirac-Coulomb + Real is an invalid option",std::cout);
  }

  template <>
  void GTODirectRelERIContraction<dcomplex,dcomplex>::directScaffold(
      MPI_Comm comm, const bool screen,
      std::vector<TwoBodyContraction<dcomplex>> &matList) const {
    CErr("Complex integral is is an invalid option",std::cout);
  }


  /**
   *  \brief Perform various tensor contractions of the full ERI
   *  tensor in core. Wraps other helper functions and provides
   *  loop structure
   *
   *  Currently supports
   *    - Coulomb-type (34,12) contractions
   *    - Exchange-type (23,12) contractions
   *
   *  Works with both real and complex matricies
   *
   *  \param [in/out] list Contains information pertinent to the
   *    matricies to be contracted with. See TwoBodyContraction
   *    for details
   */
  template <typename MatsT, typename IntsT>
  void GTODirectRelERIContraction<MatsT,IntsT>::twoBodyContract3Index(
      MPI_Comm comm, std::vector<TwoBodyContraction<MatsT>> &list) const {

    ROOT_ONLY(comm);

    auto topIncore = tick();

    // Loop over matricies to contract with
    for(auto &C : list) {

      // Coulomb-type (34,12) ERI contraction
      // AX(mn) = (mn | kl) X(kl)
      if( C.contType == COULOMB ) {
        JContract3Index(comm,C);
        // Exchange-type (23,12) ERI contraction
        // AX(mn) = (mk |ln) X(kl)
      } else if( C.contType == EXCHANGE ) {
        KContract3Index(comm,C);
      }

    } // loop over matricies

    auto durIncore = tock(topIncore);


  }; // GTODirectRelERIContraction::twoBodyContractIncore


  /**
   *  \brief Perform a Coulomb-type (34,12) ERI contraction with
   *  a one-body operator.
   */
  template <typename MatsT, typename IntsT>
  void GTODirectRelERIContraction<MatsT,IntsT>::JContract3Index(
      MPI_Comm comm, TwoBodyContraction<MatsT> &C) const {

    CErr("Please install CQDirect4C library to run direct four-component implementation.");

  }; // GTODirectRelERIContraction::JContract3Index


  template <typename MatsT, typename IntsT>
  void GTODirectRelERIContraction<MatsT,IntsT>::KContract3Index(
      MPI_Comm comm, TwoBodyContraction<MatsT> &C) const {

    CErr("Please install CQDirect4C library to run direct four-component implementation.");

  }; // GTODirectRelERIContraction::KContract3Index

}; // namespace ChronusQ











