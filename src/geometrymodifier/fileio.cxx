/*
 *  This file is part of the Chronus Quantum (ChronusQ) software package
 *
 *  Copyright (C) 2014-2020 Li Research Group (University of Washington)
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
#include <realtime.hpp>
#include <geometrymodifier/moleculardynamics.hpp>

namespace ChronusQ {

  void MolecularDynamics::createMDDataSets(Molecule& mol, std::shared_ptr<SingleSlaterBase> ss){

    size_t maxPoints = (mdOptions.saveAllGeometry) ? mdOptions.nNuclearSteps*mdOptions.nMidpointFockSteps : mdOptions.nNuclearSteps;

    savFile.createGroup("MD");

    savFile.createDataSet<size_t>("MD/LASTSAVEPOINT", {1});
    savFile.createDataSet<size_t>("MD/LASTGRADIENTSAVEPOINT", {1});
    savFile.createDataSet<size_t>("MD/MAXSAVEPOINTS", {1});

    savFile.createDataSet<size_t>("MD/STEP", {maxPoints});
    savFile.createDataSet<double>("MD/TIME", {maxPoints});

    savFile.createDataSet<double>("MD/ETOT0",   {1});
    savFile.createDataSet<double>("MD/ETOTPREV",{1});
    savFile.createDataSet<double>("MD/ETOT",    {maxPoints});
    savFile.createDataSet<double>("MD/EKIN",    {maxPoints});
    savFile.createDataSet<double>("MD/EPOT",    {maxPoints});

    savFile.createDataSet<double>("MD/TRAJECTORY",        {maxPoints*mol.nAtoms*3});
    savFile.createDataSet<double>("MD/FORCES",            {maxPoints*mol.nAtoms*3});
    savFile.createDataSet<double>("MD/VELOCITY_FULLSTEP", {maxPoints*mol.nAtoms*3});
    savFile.createDataSet<double>("MD/VELOCITY_HALFSTEP", {maxPoints*mol.nAtoms*3});

    if( auto ss_t = std::dynamic_pointer_cast<SingleSlater<double,double>>(ss) )
      createOnePDM<double,double>(ss_t);
    else if( auto ss_t = std::dynamic_pointer_cast<SingleSlater<dcomplex,double>>(ss) )
      createOnePDM<dcomplex,double>(ss_t);
    else if( auto ss_t = std::dynamic_pointer_cast<SingleSlater<dcomplex,dcomplex>>(ss) )
      createOnePDM<dcomplex,dcomplex>(ss_t);
    else
      CErr("Unsuccessful Cast in MolecularDynamics::createMDDataSets!");

  }

  void MolecularDynamics::saveState(Molecule& mol, std::shared_ptr<SingleSlaterBase> ss){

    std::cout << "  *** MD Saving step #"<< curState.iStep <<"( t = "<< curState.time <<" au) to binary file ***" << std::endl;
    size_t maxPoints = (mdOptions.saveAllGeometry) ? mdOptions.nNuclearSteps*mdOptions.nMidpointFockSteps : mdOptions.nNuclearSteps;
    savFile.safeWriteData("MD/MAXSAVEPOINTS", &(maxPoints), {1});
    savFile.safeWriteData("MD/LASTSAVEPOINT", &(curState.lastSavePoint),  {1});
    
    // If this is a gradient step, update "LASTGRADIENTSAVEPOINT" to allow full gradient-step restarting
    if (mdOptions.nMidpointFockSteps == 0 || curState.iStep % mdOptions.nMidpointFockSteps == 0 )
      curState.lastGradientSavePoint = curState.lastSavePoint;
    savFile.safeWriteData("MD/LASTGRADIENTSAVEPOINT", &(curState.lastGradientSavePoint),  {1});
    
    savFile.partialWriteData("MD/STEP", &curState.iStep,            {curState.lastSavePoint},{1},{0},{1});
    savFile.partialWriteData("MD/TIME", &curState.time,             {curState.lastSavePoint},{1},{0},{1});

    savFile.safeWriteData("MD/ETOT0",    &totalEnergy0,              {1});
    savFile.safeWriteData("MD/ETOTPREV", &previousTotalEnergy,       {1});
    savFile.partialWriteData("MD/ETOT",  &currentTotalEnergy,        {curState.lastSavePoint},{1},{0},{1});
    savFile.partialWriteData("MD/EKIN",  &nuclearKineticEnergy,      {curState.lastSavePoint},{1},{0},{1});
    savFile.partialWriteData("MD/EPOT",  &electronicPotentialEnergy, {curState.lastSavePoint},{1},{0},{1});

    size_t len3D    = mol.nAtoms*3;
    size_t offset3D = curState.lastSavePoint*len3D;
    std::vector<double> totalCoordinates = mol.getTotalCoordinates();
    savFile.partialWriteData("MD/TRAJECTORY", &totalCoordinates[0], {offset3D},{len3D},{0},{len3D});

    std::vector<double> forces(gradientCurrent.size());
    std::transform(gradientCurrent.begin(), gradientCurrent.end(), forces.begin(), [](double coord) { return -coord; });
    savFile.partialWriteData("MD/FORCES", &forces[0], {offset3D},{len3D},{0},{len3D});

    savFile.partialWriteData("MD/VELOCITY_FULLSTEP", &velocityCurrent[0], {offset3D},{len3D},{0},{len3D});
    savFile.partialWriteData("MD/VELOCITY_HALFSTEP", &velocityHalfTime[0], {offset3D},{len3D},{0},{len3D});

    if( auto ss_t = std::dynamic_pointer_cast<SingleSlater<double,double>>(ss) )
      writeOnePDM<double,double>(ss_t);
    else if( auto ss_t = std::dynamic_pointer_cast<SingleSlater<dcomplex,double>>(ss) )
      writeOnePDM<dcomplex,double>(ss_t);
    else if( auto ss_t = std::dynamic_pointer_cast<SingleSlater<dcomplex,dcomplex>>(ss) )
      writeOnePDM<dcomplex,dcomplex>(ss_t);
    else
      CErr("Unsuccessful Cast in MolecularDynamics::saveState!");

    curState.lastSavePoint++;
    
  }

  void MolecularDynamics::restoreState(Molecule& mol, std::shared_ptr<SingleSlaterBase> ss) {
   

    hsize_t maxSavePoints, lastSavePoint;
    savFile.readData("MD/MAXSAVEPOINTS", &maxSavePoints);
    size_t expectedMaxPoints = (mdOptions.saveAllGeometry) ? mdOptions.nNuclearSteps*mdOptions.nMidpointFockSteps : mdOptions.nNuclearSteps;
    if ( maxSavePoints != expectedMaxPoints ) CErr("Mismatched requested and saved propagation length!");
    
    // Restart from full gradient steps
    savFile.readData("MD/LASTGRADIENTSAVEPOINT", &lastSavePoint);

    if (mdOptions.restoreFromNuclearStep < 0) {
      curState.lastSavePoint = lastSavePoint;
    } else {
      curState.lastSavePoint = mdOptions.restoreFromNuclearStep;
      if(mdOptions.restoreFromNuclearStep > lastSavePoint) CErr("Cannot restart from a not-yet calculated time-step!");
    } 

    curState.lastGradientSavePoint = curState.lastSavePoint;
    

    savFile.partialReadData("MD/STEP", &curState.iStep, {curState.lastSavePoint}, {1}, {0}, {1});
    savFile.partialReadData("MD/TIME", &curState.time,  {curState.lastSavePoint}, {1}, {0}, {1});

    std::cout << "  *** MD Reading step #"<< curState.iStep <<"( t = "<< curState.time <<" au) to binary file ***" << std::endl;

    savFile.readData("MD/ETOT0",       &totalEnergy0);
    savFile.readData("MD/ETOTPREV",    &previousTotalEnergy);
    savFile.partialReadData("MD/ETOT", &currentTotalEnergy,        {curState.lastSavePoint},{1},{0},{1});
    savFile.partialReadData("MD/EKIN", &nuclearKineticEnergy,      {curState.lastSavePoint},{1},{0},{1});
    savFile.partialReadData("MD/EPOT", &electronicPotentialEnergy, {curState.lastSavePoint},{1},{0},{1});

    size_t len3D    = mol.nAtoms*3;
    size_t offset3D = curState.lastSavePoint*len3D;
    std::vector<double> totalCoordinates(len3D);
    savFile.partialReadData("MD/TRAJECTORY", &totalCoordinates[0], {offset3D},{len3D},{0},{len3D});
    mol.setCoordinates(totalCoordinates);

    std::vector<double> forces(gradientCurrent.size());
    savFile.partialReadData("MD/FORCES", &forces[0], {offset3D},{len3D},{0},{len3D});
    std::transform(forces.begin(), forces.begin(), gradientCurrent.end(), [](double coord) { return -coord; });

    savFile.partialReadData("MD/VELOCITY_FULLSTEP", &velocityCurrent[0],  {offset3D},{len3D},{0},{len3D});
    savFile.partialReadData("MD/VELOCITY_HALFSTEP", &velocityHalfTime[0], {offset3D},{len3D},{0},{len3D});

    // Restore time dependent density
    if( auto ss_t = std::dynamic_pointer_cast<SingleSlater<double,double>>(ss) )
      readOnePDM<double,double>(ss_t);
    else if( auto ss_t = std::dynamic_pointer_cast<SingleSlater<dcomplex,double>>(ss) )
      readOnePDM<dcomplex,double>(ss_t);
    else if( auto ss_t = std::dynamic_pointer_cast<SingleSlater<dcomplex,dcomplex>>(ss) )
      readOnePDM<dcomplex,dcomplex>(ss_t);
    else
      CErr("Unsuccessful Cast in MolecularDynamics::restoreState!");

    mol.update();
    updateBasisIntsHamiltonian();
  
    //curState.lastSavePoint++;

  }; // RealTime::restoreState

  template <typename MatsT, typename IntsT>
  void MolecularDynamics::createOnePDM(const std::shared_ptr<SingleSlater<MatsT, IntsT>> ss) {
    size_t maxPoints = (mdOptions.saveAllGeometry) ? mdOptions.nNuclearSteps*mdOptions.nMidpointFockSteps : mdOptions.nNuclearSteps;
    std::vector<std::shared_ptr<cqmatrix::Matrix<MatsT>>> onePDMs = ss->getOnePDM();
    for( size_t i = 0; i < onePDMs.size(); i++ ) {
      size_t nBasis = onePDMs[i]->dimension();
      savFile.createDataSet<MatsT>("MD/1PDM"+std::to_string(i), {maxPoints*nBasis*nBasis});
    }
  }

  template <typename MatsT, typename IntsT>
  void MolecularDynamics::writeOnePDM(const std::shared_ptr<SingleSlater<MatsT, IntsT>> ss) {
    std::vector<std::shared_ptr<cqmatrix::Matrix<MatsT>>> onePDMs = ss->getOnePDM();
    for( size_t i = 0; i < onePDMs.size(); i++ ) {
      size_t nBasis = onePDMs[i]->dimension();
      savFile.partialWriteData("MD/1PDM"+std::to_string(i), onePDMs[i]->pointer(), 
          {curState.lastSavePoint*nBasis*nBasis},{nBasis * nBasis}, {0}, {nBasis * nBasis});
      //onePDMs[i]->output(std::cout, "Write 1PDM"+std::to_string(i), true);
    }
  }

  template <typename MatsT, typename IntsT>
  void MolecularDynamics::readOnePDM(const std::shared_ptr<SingleSlater<MatsT, IntsT>> ss) {
    std::vector<std::shared_ptr<cqmatrix::Matrix<MatsT>>> onePDMs = ss->getOnePDM();
    for( size_t i = 0; i < onePDMs.size(); i++ ) {
      size_t nBasis = onePDMs[i]->dimension();
      savFile.partialReadData("MD/1PDM"+std::to_string(i), onePDMs[i]->pointer(),
          {curState.lastSavePoint*nBasis*nBasis},{nBasis * nBasis}, {0}, {nBasis * nBasis});
    }
    // Transform shared pointers into raw objects
    std::vector<cqmatrix::Matrix<MatsT>> tempOnePDMs;
    tempOnePDMs.reserve(onePDMs.size());
    std::transform(onePDMs.begin(), onePDMs.end(), std::back_inserter(tempOnePDMs),
                   [](const std::shared_ptr<cqmatrix::Matrix<MatsT>>& ptr) { return *ptr; });
    for( size_t i = 0; i < tempOnePDMs.size(); i++ ) {
      size_t nBasis = tempOnePDMs[i].dimension();
      //tempOnePDMs[i].output(std::cout, "Read 1PDM"+std::to_string(i), true);
    }
    ss->setOnePDMAO(tempOnePDMs.data());
  }


void MolecularDynamics::parseVelocityFromInput( Molecule &mol, std::string &velocityStr, std::ostream &out) {

    std::istringstream velocityStream; velocityStream.str(velocityStr);
    std::vector<std::string> tokens;
    std::vector<Atom> atoms;
    std::vector<double> velocity;
    std::locale loc;

    // Loop over lines of velocity specification
    size_t iAtom = 0;
    velocity.reserve(mol.nAtoms*3);
    for(std::string line; std::getline(velocityStream, line); ){
      split(tokens,line," \t");

      if( tokens.size() == 0 ) continue;

      if( tokens.size() != 4 ) CErr("Error in velocity reader. A line should have 4 entries: Atom Symbol, v_x, v_y, v_z, ");

      for( auto i=0; i<tokens.size(); i++ )
        if( tokens[i].find("NAN") != std::string::npos or tokens[i].find("INF") != std::string::npos ) CErr("Invalid entry for GEOM!");

      std::string atmSymb = tokens[0];

      // Checking if atom specified by symbol or number
      bool hasDig = std::any_of(atmSymb.begin(),atmSymb.end(),
        [&](char a) { return std::isdigit(a,loc); });
      bool hasAlpha = std::any_of(atmSymb.begin(),atmSymb.end(),
        [&](char a) { return std::isalpha(a,loc); });

      // After isotope there should only be a atomic number or symbol left
      if( hasDig ){

        auto it =
        std::find_if(atomicReference.begin(),atomicReference.end(),
          [&](std::pair<std::string,Atom> st){
            return st.second.atomicNumber == std::stoi(atmSymb);}
           );

        std::string parseAtmSymb = it->first.substr(0,it->first.find("-",0));
        atoms.emplace_back((it == atomicReference.end() ? "X" : defaultIsotope[parseAtmSymb]));

      } else if( hasAlpha ){

        atoms.emplace_back(defaultIsotope[atmSymb]);

      }

      if (atoms.back().atomicNumber != mol.atoms[iAtom].atomicNumber)
        CErr("Order of atoms in velocity specification must match the order in geometry specification!");

      // TODO: Handle unit conversion 
      // For now assuming it's already in atomic unit
      velocity.emplace_back(std::stod(tokens[1]));
      velocity.emplace_back(std::stod(tokens[2]));
      velocity.emplace_back(std::stod(tokens[3]));

      iAtom++;
    }

    //if ( atoms.size() == 0 )
    //  CErr("MOLECULE.GEOM must not be empty and must be indented");

    if (velocity.size() != mol.nAtoms * 3) 
        CErr("The size of the velocity vector must be 3 times the number of atoms.");
    
    std::copy(velocity.begin(), velocity.end(), velocityCurrent.begin());

    std::cout << "Read in Velocity:"<<std::endl;
    size_t i = 0;
    for( Atom& atom : mol.atoms ) {

      std::cout << std::right <<"AtomicNumber = " << std::setw(4) << atom.atomicNumber 
                << std::right <<"  X= "<< std::setw(24) <<  velocityCurrent[i  ]
                << std::right <<"  Y= "<< std::setw(24) <<  velocityCurrent[i+1]
                << std::right <<"  Z= "<< std::setw(24) <<  velocityCurrent[i+2]<<std::endl;
      i += 3;
 
    }


  } // parseVelocity


  template void MolecularDynamics::createOnePDM<double, double>(const std::shared_ptr<SingleSlater<double, double>>);
  template void MolecularDynamics::createOnePDM<dcomplex, double>(const std::shared_ptr<SingleSlater<dcomplex, double>>);
  template void MolecularDynamics::createOnePDM<dcomplex, dcomplex>(const std::shared_ptr<SingleSlater<dcomplex, dcomplex>>);

  template void MolecularDynamics::writeOnePDM<double, double>(const std::shared_ptr<SingleSlater<double, double>>);
  template void MolecularDynamics::writeOnePDM<dcomplex, double>(const std::shared_ptr<SingleSlater<dcomplex, double>>);
  template void MolecularDynamics::writeOnePDM<dcomplex, dcomplex>(const std::shared_ptr<SingleSlater<dcomplex, dcomplex>>);

  template void MolecularDynamics::readOnePDM<double, double>(const std::shared_ptr<SingleSlater<double, double>>);
  template void MolecularDynamics::readOnePDM<dcomplex, double>(const std::shared_ptr<SingleSlater<dcomplex, double>>);
  template void MolecularDynamics::readOnePDM<dcomplex, dcomplex>(const std::shared_ptr<SingleSlater<dcomplex, dcomplex>>);
}