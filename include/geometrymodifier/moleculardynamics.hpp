/*
 *  This file is part of the Chronus Quantum (ChronusQ) software package
 *
 *  Copyright (C) 2014-2026 Li Research Group (University of Washington)
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
#include <cmath>
#include <physcon.hpp>
#include <molecule.hpp>
#include <geometrymodifier.hpp>
namespace ChronusQ {

  struct NEODynamicsOptions {
    bool tpb = false;            // Default NEO-Ehrenfest alg is to fix quantum proton basis center
    bool includeQProtKE = false; // Whether to include the translational KE associated with changes in protonic basis centers         
  }; // struct NEODynamicsOptions

  struct MDOptions {

    size_t nNuclearSteps;                // Number of nuclear steps for molecular dynamics
    size_t nMidpointFockSteps = 10;      // Number of nuclear steps between each gradient calculations (for Ehfenfest calculations)
    size_t nElectronicSteps = 5;         // Number of RT elctronic steps for electron dynamics
    long int restoreFromNuclearStep = 0; // < Restore MD from this nuclear step 
    
    // Save data for all geometry steps, including all the mid-point Fock steps with a nuclear step
    bool saveAllGeometry = false;

    // Print properties at each step
    bool printProperty = false;

    // Molecular Dynamics Options
    double timeStepAU; // Nuclear timestep for molecular dynamics in a.u.
    double timeStepFS; // Nuclear timestep for molecular dynamics in fs

    bool pertFirstAtom = false; ///< Perturb first atom's geometry at t=0
    double pert_val_x = 1e-5;   ///< perturbation value
    double pert_val_y = 1e-5;   ///< perturbation value
    double pert_val_z = 1e-5;   ///< perturbation value

    bool projectOrthoDen = false; // Project the orthonormal density to initial geometry
    bool onlyMoveH = false; // Only move the hydrogen atoms

    MDOptions(double tmax, double deltat)
    {
      timeStepAU = deltat;
      timeStepFS = deltat*FSPerAUTime();

      nNuclearSteps = size_t(ceil(tmax/deltat));
    }

  }; // struct MDOptions
  std::ostream& operator<<(std::ostream&, const MDOptions&);

  struct MDProgress {

    size_t  iStep = 0;         ///< Step index of current time point
    double  time = 0.;         ///< keeps track of the geometry time during all the mid-point Fock steps
    double  ptime = 0.;        ///< keeps track of the velocity time during all the mid-point Fock steps
    double  ptimeHalf = 0.;    ///< keeps track of the half-step velocity time during all the mid-point Fock steps
    double  currentStepSize;   ///< Current step size at which we move the geometry
    size_t  lastSavePoint = 0; ///< last nuclear step that was saved 
    size_t  lastGradientSavePoint = 0; ///< last gradient nuclear step that was saved (restart from here)

  };

  /**
   * \brief The MolecularDynamics class
   */
  class MolecularDynamics : public GeometryModifier {

    MDProgress curState;
    double tMax;
    SafeFile savFile; ///< Data File

  public:

    std::function<std::vector<double>()> gradientGetter;
    std::function<double()> finalMidpointFock;
    std::function<void()> updateBasisIntsHamiltonian;
    std::function<void()> pertFirstAtom;
    std::vector<double> gradient;
    std::vector<double> velocity;          ///< nuclear velocity at the current time (t)
    std::vector<double> velocityHalfTN;    ///< nuclear velocity at the half Δt_N step (t+0.5*Δt_N)
                                           ///< We use this to get full tN step velocity
    std::vector<double> acceleration;      ///< acceleration at the current time (t)

    double   nuclearKineticEnergy;      ///< nuclear kinetic energy
    double   totalEnergy0;              ///< total energy at step 0
    double   currentTotalEnergy;        ///< total energy at the current step
    double   previousTotalEnergy;       ///< total energy at the previous step
    
    MDOptions mdOptions;
    NEODynamicsOptions NEODynamicsOpts;


    // Constructors
    MolecularDynamics() = delete;
    MolecularDynamics(MDOptions mdOptions, Molecule& molecule, SafeFile& rstFile, MPI_Comm mpiComm):
      GeometryModifier(mpiComm),
      mdOptions(mdOptions),
      gradient(3*molecule.nAtoms, 0.),
      velocity(3*molecule.nAtoms, 0.),
      velocityHalfTN(3*molecule.nAtoms, 0.),
      acceleration(3*molecule.nAtoms, 0.)
    {
      curState.currentStepSize = mdOptions.timeStepAU;
      tMax = mdOptions.timeStepAU * mdOptions.nNuclearSteps;
      savFile = rstFile;
    };

    // Different type
    MolecularDynamics(const MolecularDynamics &other):
        GeometryModifier(other),
        mdOptions(other.mdOptions){}
    MolecularDynamics(MolecularDynamics &&other):
        GeometryModifier(other),
        mdOptions(other.mdOptions){}

    const MDOptions& getMDOptions() const {
      return mdOptions;
    }

    // Virtual destructor
    virtual ~MolecularDynamics() {}

    virtual bool hasNext() override;

    virtual void update(bool print, Molecule &molecule, bool firstStep, 
        TDSCFOptions& tdSCFOptions, std::shared_ptr<SingleSlaterBase> ss,
        EMPerturbation& emPert, std::vector<std::shared_ptr<CubeGen>> cubes = {} ) override;

    void velocityVV(Molecule &molecule, std::vector<double>& vIn, std::vector<double>& vOut, std::vector<double> gradient, double timeStep);
    void velocityVV_EXPK1(Molecule &molecule, std::vector<double>& vIn, std::vector<double>& vOut, std::vector<double> gradient, double timeStep, EMPerturbation& emPert);

    void geometryVV(Molecule &molecule, double timeStep);
    
    void computeKineticEnergy(Molecule &molecule);

    void printCurrentGeometry(Molecule &molecule);
    
    void printCurrentVelocity(Molecule &molecule);

    void printMDInfo(Molecule &molecule, double totalEnergy);
    
    void initializeMD(Molecule& molecule, std::shared_ptr<SingleSlaterBase>);

    void createMDDataSets(Molecule& mol, std::shared_ptr<SingleSlaterBase>);

    void saveState(Molecule& mol, std::shared_ptr<SingleSlaterBase>);
    
    void restoreState(Molecule& mol, std::shared_ptr<SingleSlaterBase>);

    template <typename MatsT, typename IntsT>
    void createOnePDM(const std::shared_ptr<SingleSlater<MatsT, IntsT>> ss);

    template <typename MatsT, typename IntsT>
    void writeOnePDM(const std::shared_ptr<SingleSlater<MatsT, IntsT>> ss);

    template <typename MatsT, typename IntsT>
    void readOnePDM(const std::shared_ptr<SingleSlater<MatsT, IntsT>> ss);

    void parseVelocityFromInput(Molecule &mol, std::string &velocityStr, std::ostream &out);
  };
}
