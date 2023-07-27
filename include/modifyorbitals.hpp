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

/*
 *     Brief: This header defines the inteface between ModifyOrbitals Objects and
 *            the objects that call them. The objects inherit the SCFInterface
 *            base class
 *
 */

namespace ChronusQ {

  // Assign input types to alias for ease of use
template<typename MatsT>
using VecMORef = std::vector<std::reference_wrapper<SquareMatrix<MatsT>>>;
using VecEPtr  = std::vector<double*>;
template<typename MatsT>
using VecShrdPtrMat = std::vector<std::shared_ptr<SquareMatrix<MatsT>>>;
template<typename MatsT>
using VecShrdPtrOrtho = std::vector<std::shared_ptr<Orthogonalization<MatsT>>>;

/*
 *   These are the functions that should be defined in order to use the
 *   modifyOrbitals objects.
 *
 *
 */
template<typename MatsT>
struct ModifyOrbitalsOptions {
    //--------------------------------------------------------------------------------------------
    // I/O  FUNCTIONS
    std::function<void()> printProperties;    ///< Prints the properties after completion of SCF
    std::function<void()> saveCurrentState;   ///< Saves the current State to Disk

    //--------------------------------------------------------------------------------------------

    //--------------------------------------------------------------------------------------------
    // NECESSARY COMPUTATION FUNCTIONS
    std::function<void(EMPerturbation&)> formFock;            ///< Form the Fock Matrix
    std::function<void(EMPerturbation&)> computeProperties;   ///< Compute Implemented Properties
    std::function<void()> formDensity;                        ///< Form the Density

    //--------------------------------------------------------------------------------------------

    //--------------------------------------------------------------------------------------------
    // GETTER FUNCTIONS
    std::function<double()> getTotalEnergy;   ///< Returns the energy for Convergence test
    std::function<VecShrdPtrMat<MatsT>()>
        getFock;   ///< Returns a vector of shared_ptr to the fock Matrix(same basis as MOs)
    std::function<VecShrdPtrMat<MatsT>()>
        getOnePDM;   ///< Returns a vector of shared_ptr to density matrix(same basis as MOs)
    std::function<VecShrdPtrOrtho<MatsT>()>
        getOrtho;   ///< Returns a vector to the orthogonalization object(same basis as MOs)
    //--------------------------------------------------------------------------------------------

    //--------------------------------------------------------------------------------------------
    // SETTER FUNCTIONS
    std::function<void(bool val)> setDenEqCoeff;   ///< Setting the bool flag in SingleSlater whether Coeffs and Density present the same wavefunction
    //--------------------------------------------------------------------------------------------

    //--------------------------------------------------------------------------------------------
    // OPTIONAL FUNCTIONS
    // User defined Functions to overload default gradient functions in the algorithms
    // if you want to use your own functions then bind a function to these, otherwise the
    // code will use the default.

    std::function<void(MatsT*)> computeFullNRStep;   ///< Function to compute the full Newton-Raphson step if Hessian isn't implemented call CErr

    std::function<void(MatsT*)> computeNROrbGrad;   ///< Computes and returns the Orbital Gradient in the Nonredundant Basis
                                                    ///<   this replaces the gradient function in NewtonRaphsonSCF if an
                                                    ///<   approximation is being used

    std::function<void(std::vector<SquareMatrix<MatsT>>&)> computeErrorVector;   ///< Computes the Orbital Gradient in the AO Basis
                                                                                 ///<   this replaces the F,D Commutator function in ConventionalSCF
    //--------------------------------------------------------------------------------------------
};

/*
 *   Brief: Abstract Base class for the modifyOrbitals object. This allows
 *          us to abstract the algorithms that modify the orbitals into one
 *          interface. This is the interface that the object that owns ModifyOrbitals
 *          runs the ModifyOrbitals algorithms.
 */
template<typename MatsT>
class ModifyOrbitals {

  protected:
    MPI_Comm comm;   ///< MPI Communication
    CQMemManager& memManager;
    ModifyOrbitalsOptions<MatsT> modOrbOpt;   ///< Options struct for Modify Orbitals

  public:
    ModifyOrbitals() = delete;
    ModifyOrbitals(MPI_Comm c, ModifyOrbitalsOptions<MatsT> m, CQMemManager& mem): comm(c), modOrbOpt(m), memManager(mem) {
        // Check that Essential functions are bound
        if( not modOrbOpt.printProperties ) CErr("printProperties was not bound in ModifyOrbitalOptions");
        if( not modOrbOpt.saveCurrentState ) CErr("saveCurrentState was not bound in ModifyOrbitalOptions");
        if( not modOrbOpt.formFock ) CErr("formFock was not bound in ModifyOrbitalOptions");
        if( not modOrbOpt.formDensity ) CErr("formDensity was not bound in ModifyOrbitalOptions");
        if( not modOrbOpt.computeProperties ) CErr("computeProperties was not bound in ModifyOrbitalOptions");
        if( not modOrbOpt.getTotalEnergy ) CErr("getTotalEnergy was not bound in ModifyOrbitalOptions");
        if( not modOrbOpt.getFock ) CErr("getFock was not bound in ModifyOrbitalOptions");
        if( not modOrbOpt.getOnePDM ) CErr("getOnePDM was not bound in ModifyOrbitalOptions");
        if( not modOrbOpt.getOrtho ) CErr("getOrtho was not bound in ModifyOrbitalOptions");
        if( not modOrbOpt.setDenEqCoeff ) CErr("setDenEqCoeff was not bound in ModifyOrbitalOptions");
    };

    // Getter/Setter functions for options struct
    ModifyOrbitalsOptions<MatsT> getModOpt() { return modOrbOpt; };
    void setModOpt(const ModifyOrbitalsOptions<MatsT>& mod) { modOrbOpt = mod; };

    // runModifyOrbitals is the driver function that performs
    // the whole optimization/simulation
    virtual void runModifyOrbitals(EMPerturbation& pert, VecMORef<MatsT>&, VecEPtr&) = 0;

    // getNewOrbitals performs only a single step in the
    // optimization/simulation
    virtual void getNewOrbitals(EMPerturbation& pert, VecMORef<MatsT>&, VecEPtr&) = 0;

    // Printing functions
    virtual void printRunHeader(std::ostream& out, EMPerturbation&) const                   = 0;
    virtual void printIteration(std::ostream& out = std::cout, bool printDiff = true) const = 0;
};

};   // Namespace ChronusQ

#include <modifyorbitals/optOrbitals.hpp>
#include <modifyorbitals/conventionalSCF.hpp>
#include <modifyorbitals/newtonRaphsonSCF.hpp>
#include <modifyorbitals/skipSCF.hpp>
