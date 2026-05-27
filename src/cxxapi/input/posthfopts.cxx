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
#include <cxxapi/options.hpp>
#include <cerr.hpp>
#include <posthartreefock/base.hpp>
#include <regex>

namespace ChronusQ {
  
  void HandlePostHFProperties(std::ostream &out, CQInputFile &input,
    std::shared_ptr<PostHartreeFockBase> & postHF, 
    std::string postHFSection) {
  
    // Mulliken charge analysis
    OPTOPT( postHF->PopulationAnalysis = input.getData<bool>(postHFSection + ".POPULATION"); )
 
    // Spin analysis
    OPTOPT( postHF->SpinAnalysis = input.getData<bool>(postHFSection + ".PRINTSPIN"); )

    // Oscillator strength
    OPTOPT( postHF->osc_str = input.getData<bool>(postHFSection + ".OSCISTREN"); )
    if ( postHF->osc_str ){
      OPTOPT( postHF->osc_str_order = input.getData<size_t>(postHFSection + ".OSCISTREN_ORDER"); )
      OPTOPT( postHF->NosS1 = input.getData<size_t>(postHFSection + ".OSCISTREN_INITSTATES"); )
    }
    else{
      if (input.containsData(postHFSection + ".OSCISTREN_ORDER") || 
          input.containsData(postHFSection + ".OSCISTREN_INTISTATES")){
        CErr("Cannot Set OSCISTREN_ORDER or OSCISTREN_INITSTATES Without Setting OSCISTREN == True/1");}} 

    // Printing Options
    // MOs
    if ( input.containsData(postHFSection + ".PRINTMOS") ) {
      try { postHF->printMOCoeffs = input.getData<size_t>(postHFSection + ".PRINTMOS"); }
      catch(...) {
        CErr("Invalid PRINTMOS input. Please use number 0 ~ 9.");
      }
    }

    // Print determinant occupations?
    OPTOPT( postHF->printDetailedCICoeffs = input.getData<bool>(postHFSection + ".PRINTDETOCC"); )

    if ( postHF->printMOCoeffs >= 10 ) CErr(postHFSection + " print level is not valid!");
  
  }; // 

  void HandlePostHFRDMPrinting(std::ostream &out, CQInputFile &input,
    std::shared_ptr<PostHartreeFockBase> & postHF, 
    std::string postHFSection) {

    // Parse RDM printing
    std::string printRDMString;
    OPTOPT( printRDMString = input.getData<std::string>(postHFSection + ".PRINTRDMS"));
    if ( not printRDMString.empty() ) {
      std::cout << "  * Printing RDM detected: " << std::endl;

      std::vector<std::string> rdmTokens;
      split(rdmTokens, printRDMString, " \t,");

      if( rdmTokens.size() != 1 and rdmTokens.size() != 2 ) CErr("Need 1 or 2 entries in single line for RDM printing");

      // Parse rdmCut if present
      if( rdmTokens.size() == 2 ) postHF->rdmCut=std::stod(trim(rdmTokens[1]));

      try { postHF->printRDMs = std::stoi(trim(rdmTokens[0])); }
      catch(...) {
        CErr("Invalid PRINTRDMS input. Please use number 0 ~ 2.");
      }
      if (postHF->printRDMs >= 3 ) CErr(postHFSection + " print RDM level is not valid!");

    }

  } // HandleRDMPrinting

  std::unordered_map<std::string,int> OrbSpinMap = {
    { "A" , 0  },
    { "B" , 1  }
  };

  void HandlePostHFOrbitalSwaps(std::ostream &out, CQInputFile &input,
    std::shared_ptr<SingleSlaterBase> &ss, 
    std::shared_ptr<PostHartreeFockBase> &postHF,
    std::string postHFSection) {

    // MO swapping
    std::string swapMOStrings;
    OPTOPT( swapMOStrings = input.getData<std::string>(postHFSection + ".SWAPMO"));
    if ( not swapMOStrings.empty() ) {
      std::cout << "  * Manually MO Swapping Detected: " << std::endl;

      // Pair function for [PostHF] MO swap
      std::vector<std::vector<std::pair<size_t, size_t>>> moPairs;
      moPairs.resize(2, {});

      std::vector<std::string> moTokens;
      //Loop over lines of mo swapping
      std::istringstream moStream(swapMOStrings);

      for( std::string line; std::getline(moStream, line); ) {
        split(moTokens, line, " \t,");

        if( moTokens.size() == 0 ) continue;
        else if( moTokens.size() != 2 and moTokens.size() != 3 ) CErr("Need 2 or 3 entries in single line for swapping");

        // Parse spin if present
        std::string spinDir("A");
        if( moTokens.size() == 3 ) spinDir=moTokens[2];
        trim(spinDir);

        // mo[1] NYI for MCSCF
        if( spinDir == "B" ) CErr("Swapping of beta MOs NYI for " + postHFSection);

        moPairs[OrbSpinMap[spinDir]].emplace_back(std::stoul(moTokens[0]), std::stoul(moTokens[1]));
      }

      postHF->swapMOs(moPairs,isAlpha);

    }

  }; // HandlePostHFOrbitalSwaps
 
  void ConstructActiveSpaces(std::ostream & out, CQInputFile & input,
                             const std::vector<size_t> & nActOs, size_t nActE,
                             size_t MOOffset,
                             int maxInterspaceEX,
                             std::vector<ActiveSpaceParameters> & actS,
                             std::vector<std::vector<size_t>> & refOcc,
                             std::string postHFSection) {

    // Parse DAS input section
    // {[Ne, No, "Label"],[],[], 5+, 6-}
    size_t nDAS;
    std::vector<size_t> iDASOrb;
    std::vector<size_t> iDASEle;
    std::vector<size_t> iDASMinOcc;
    std::vector<size_t> iDASMaxOcc;
    std::vector<size_t> iDASGroup;
    std::vector<int> iDASeLimit;
    std::vector<int> iDAShLimit;
    std::vector<std::string> iDASLabel;
    std::vector<bool> iIntContraction;

    try {
      std::string DASPartition;
      OPTOPT( DASPartition = input.getData<std::string>(postHFSection +".DAS"));
      if (not DASPartition.empty()) {
        std::string iDASDefinition;
        std::istringstream inputDASStream(DASPartition);

        std::smatch nMatchDAS;
        std::regex integers("[[:digit:]]+");

        std::getline(inputDASStream, iDASDefinition);
        std::regex_search(iDASDefinition,nMatchDAS, integers);
        if (nMatchDAS.size()==0) CErr("Missing the number of distributed active space!", std::cout);
        nDAS = std::stoul(nMatchDAS.str());
        if (nDAS==0) CErr("The number of distributed active space cannot be zero.", std::cout);
        std::cout<<"Number of DAS: "<< nDAS <<std::endl;

        // If there are more than 1 DAS partitioning
        if( nDAS > 1) {
          std::regex nEle("\\b([0-9]+)(e|E)\\b"); // an integer number followed by "e" - number of electrons
          std::regex nOrb("\\b([0-9]+)(o|O)\\b"); // an integer number followed by "o" - number of orbitals
          std::regex nLabel("\"(.*?)\"");         // "Label"
          std::regex nMultiply("\\b([0-9]+)(x|X)\\b"); // an integer number followed by "x" - multiply the space

          std::regex iGroup("\\{(.*?)\\}");
          std::regex iSpace("\\[(.*?)\\]");
          std::regex orbRange("\\((.*?)\\)");
          std::regex orbRangeSelection("\\b([0-9]+)(o|O)-([0-9]+)(o|O)\\b");
          std::regex eRestriction("([0-9]+)-$");
          std::regex hRestriction("([0-9]+)\\+");
          std::regex intContraction("\\b(ic)\\b",std::regex_constants::icase);

          std::smatch nMatchGroup;
          std::smatch nMatchSpace;
          std::smatch nMatchOrbs;
          std::smatch eMax;
          std::smatch hMax;

          size_t iGroupNo = 0;
          int eExcitation;
          int hExcitation;
          size_t eGroupMax, hGroupMax;
          int typeOrb = -1;
          size_t nMultiDAS = 1;
          for(size_t iDAS = 0; iDAS < nDAS; ) {
            if(inputDASStream.eof()) CErr("Insufficient number of DAS definitions!", std::cout);
            std::getline(inputDASStream, iDASDefinition);

            // a group of DASs is defined within {} with excitation restrictions defined by
            // N+ and N- for number of additional holes and electrons allowed in the group.
            for (std::sregex_iterator itGroup = std::sregex_iterator(iDASDefinition.begin(),iDASDefinition.end(), iGroup);
            itGroup != std::sregex_iterator(); itGroup++) {

              nMatchGroup = *itGroup;

              if(nMatchGroup.str(1).size()>0) {

                std::string groupValues = nMatchGroup.str(1);
                std::regex_search(groupValues, eMax, eRestriction);
                std::regex_search(groupValues, hMax, hRestriction);

                eExcitation = 0;
                hExcitation = 0;

                std::cout<<std::endl<< "Group No. "<< iGroupNo+1 <<":"<<std::endl;
                std::cout<<"  Excitation Restrictions: ";
                if( eMax.str(1).size()>0 or hMax.str(1).size()>0) {
                  if(eMax.str(1).size()>0) {
                    eExcitation = std::stoul(eMax.str(1));
                    std::cout<<eExcitation<<"(e) ";
                  }
                  if(hMax.str(1).size()>0) {
                    hExcitation = std::stoul(hMax.str(1));
                    std::cout<<hExcitation<<"(h)";
                  }
                  std::cout<<std::endl;
                } else if ( maxInterspaceEX > 0) {
                  eExcitation = maxInterspaceEX;
                  hExcitation = maxInterspaceEX;
                  std::cout<<eExcitation<<"(e) ";
                  std::cout<<hExcitation<<"(h)";
                  std::cout<<std::endl;
                } else std::cout<<" None"<<std::endl;

                // each DAS is defined within []
                for (std::sregex_iterator itSpace = std::sregex_iterator(groupValues.begin(),groupValues.end(), iSpace);
                itSpace != std::sregex_iterator() and iDAS < nDAS; itSpace++) {

                  nMatchSpace = *itSpace;
                  if(nMatchSpace.str(1).size()>0) {

                    std::cout << "    DAS No. "<<iDAS+1<<std::endl;
                    std::string spaceValues = nMatchSpace.str(1);

                    nMultiDAS = 1;
                    // locate the multiply for the current DAS definition
                    std::regex_search(spaceValues,nMatchDAS, nMultiply);
                    if (nMatchDAS.str(1).size() > 0) {
                      std::cout<<"      nMultiply = "+nMatchDAS.str(1)<< std::endl;
                      nMultiDAS = std::stoul(nMatchDAS.str(1));
                    }

                    // locate the definition of space occupation
                    std::regex_search(spaceValues,nMatchDAS, nEle);
                    std::cout<<"      nElectrons = "+nMatchDAS.str(1)<< std::endl;
                    for (auto i = 0; i < nMultiDAS; i++) iDASEle.push_back(std::stoul(nMatchDAS.str(1)));

                    // locate the definition of orbital partitioning
                    std::regex_search(spaceValues, nMatchOrbs, orbRange);
                    if(nMatchOrbs.str(1).size()>0) {
                      if(typeOrb < 0) typeOrb = 0;
                      else if (typeOrb!=0) CErr("Orbital selection type must be the same!", std::cout);
                      std::string orbRangeValues = nMatchOrbs.str(1);
                      std::regex_search(orbRangeValues,nMatchDAS, orbRangeSelection);
                      std::cout << "       nOrbitals = #"+nMatchDAS.str(1)+"-#"+nMatchDAS.str(3)<<std::endl;
                    } else {
                      if(typeOrb < 0) typeOrb = 1;
                      else if (typeOrb!=1) CErr("Orbital selection type must be the same!", std::cout);
                      std::regex_search(spaceValues,nMatchDAS, nOrb);
                      std::cout << "       nOrbitals = "+nMatchDAS.str(1)<< std::endl;
                      for (auto i = 0; i < nMultiDAS; i++) iDASOrb.push_back(std::stoul(nMatchDAS.str(1)));
                    }
                    // locate the name of the space
                    std::regex_search(spaceValues,nMatchDAS, nLabel);
                    if(nMatchDAS.str(1).size() > 0) {
                      std::cout << "           Label = "+nMatchDAS.str(1)<< std::endl;
                      for (auto i = 0; i < nMultiDAS; i++) iDASLabel.push_back(nMatchDAS.str(1));
                    } else {
                      std::cout << "           Label = DAS "<<iDAS+1<< std::endl;
                      for (auto i = 0; i < nMultiDAS; i++) iDASLabel.push_back("DAS "+std::to_string(iDAS+1));
                    }

                    if(hExcitation > 0) {
                      int minOcc = (int) iDASEle[iDAS] - hExcitation;
                      if( minOcc > 0 ) 
                        for (auto i = 0; i < nMultiDAS; i++) 
                          iDASMinOcc.push_back(minOcc);
                      else for (auto i = 0; i < nMultiDAS; i++) 
                          iDASMinOcc.push_back(0ul);
                    } else for (auto i = 0; i < nMultiDAS; i++) 
                        iDASMinOcc.push_back(0ul);
                    //std::cout << "          minOcc = "<<iDASMinOcc[iDAS]<< std::endl;

                    if(eExcitation > 0) {
                      int maxOcc = (int) iDASEle[iDAS] + eExcitation;
                      if( maxOcc < (int) iDASOrb[iDAS] ) 
                        for (auto i = 0; i < nMultiDAS; i++) 
                          iDASMaxOcc.push_back(maxOcc);
                      else for (auto i = 0; i < nMultiDAS; i++) 
                        iDASMaxOcc.push_back(iDASOrb[iDAS]);
                    } else for (auto i = 0; i < nMultiDAS; i++) 
                      iDASMaxOcc.push_back(iDASOrb[iDAS]);
                    //std::cout << "          maxOcc = "<<iDASMaxOcc[iDAS]<< std::endl;

                    for (auto i = 0; i < nMultiDAS; i++) {
                      iIntContraction.push_back(std::regex_search(spaceValues, intContraction));
                      //std::cout <<   "Int. Contraction = "<<std::boolalpha<<iIntContraction[iDAS]<< std::endl;

                      iDASeLimit.push_back(eExcitation);
                      iDAShLimit.push_back(hExcitation);
                      iDASGroup.push_back(iGroupNo);
                      iDAS++;
                    }

                  }
                }
                iGroupNo++;
              }
            }
          }
        }
      } else nDAS = 1;
    } catch(...) {
      CErr("Something is wrong with the DAS definition!", std::cout);
    }

    if(nDAS==1) {
      iDASOrb.push_back((size_t)std::accumulate(nActOs.begin(), nActOs.end(), 0));
      iDASEle.push_back(nActE);
      iDASMinOcc.push_back(0);
      iDASMaxOcc.push_back(nActE);
      iDASGroup.push_back(0);
      iDASeLimit.push_back(0);
      iDAShLimit.push_back(0);
      iIntContraction.push_back(false);
      iDASLabel.push_back("Complete Active Space");
    }

    if(std::accumulate(nActOs.begin(), nActOs.end(), 0) != std::accumulate(iDASOrb.begin(), iDASOrb.end(), 0))
      CErr("Incorrect number of active orbitals defined in DAS!", std::cout);
    if(std::accumulate(iDASEle.begin(), iDASEle.end(), 0) != nActE)
      CErr("Incorrect number of active electrons defined in DAS!", std::cout);

    // generate active space
    actS.resize(nDAS);
    actS.clear();
    size_t offsetDAS = MOOffset;
    for (auto iDAS = 0; iDAS < nDAS; iDAS ++) {
      actS.push_back(ActiveSpaceParameters({offsetDAS, iDASOrb[iDAS], iDASEle[iDAS],
                                            iDASMinOcc[iDAS], iDASMaxOcc[iDAS],
                                            iDASGroup[iDAS], iDASeLimit[iDAS],
                                            iDAShLimit[iDAS], iDASLabel[iDAS],
                                            iIntContraction[iDAS]}));
      offsetDAS += iDASOrb[iDAS];
    }

    refOcc.resize(nDAS);
    refOcc.clear();
    refOcc.push_back(iDASEle);

#if false
    // generate active space
    // TODO: set it for one component  
    size_t off = MOOffset;
    for (const auto & nO: nActOs) {
      actS.push_back(ActiveSpaceParameters({off, nO, 0, nO}));
      off += nO;
    }
    
    // parse restrictions
    // as  each line 
    //   min_occ, max_occ, if_accumulate
    //
    std::string occResStrings;
    std::vector<std::string> occResTokens;
    OPTOPT( occResStrings = input.getData<std::string>(postHFSection + ".OCCRESTRICTIONS"));
    if (not occResStrings.empty()) {
      std::istringstream occResStream(occResStrings);
      int i = -1;
      for(std::string occRes; std::getline(occResStream, occRes); ++i) {
         split(occResTokens, occRes, " \t,");
         if (occResTokens.size() == 0) continue;
         if (occResTokens.size() < 2) CErr("Need to specify min and max Occ for each space");
         actS[i].minOcc = std::stoul(occResTokens[0]);
         actS[i].maxOcc = std::stoul(occResTokens[1]);
         if (occResTokens.size() > 2)
           actS[i].accumulateOcc = occResTokens[2] == "TRUE" 
             or occResTokens[2] == "1" or occResTokens[2] == "T";
      }
    }
#endif

  }; // ConstructActiveSpaces  
  
  void ReadReferenceOcc(std::ostream & out, CQInputFile & input,
    std::vector<std::vector<size_t>> & refOcc, std::string postHFSection) {
    
    refOcc.clear();
    std::string refOccStrings;
    std::vector<std::string> refOccTokens;
    OPTOPT( refOccStrings = input.getData<std::string>(postHFSection + ".REFERENCEOCC"));
    if (not refOccStrings.empty()) {
      std::istringstream refOccStream(refOccStrings);    
      for( std::string refOccStr; std::getline(refOccStream, refOccStr);) {
         split(refOccTokens, refOccStr, " \t,");
         if (refOccTokens.size() == 0) continue;
         
         std::vector<size_t> refOccTmp;
         
         for (const auto & occ: refOccTokens) 
           refOccTmp.push_back(std::stoul(occ));
         
         refOcc.push_back(refOccTmp);
      } 
    } 

  }; // ReadReferenceOcc 

  void HandleSavePDMSPostHF(std::ostream &out, CQInputFile &input,
    std::shared_ptr<PostHartreeFockBase> &posthf, std::string postHFSection) {

    // Parse natural orbital option
    std::string savePDMsString;
    OPTOPT( savePDMsString = input.getData<std::string>(postHFSection + ".SAVEONEPDMS"); )
    if ( not savePDMsString.empty() ) {
      std::cout << "  * Requesting to Save PDMs keyword detected: " << std::endl;

      auto const regexstatestosavepdm = std::regex("\\(.*\\)",std::regex_constants::icase);
      posthf->saveOnePDMS = true;

      if( std::regex_search(savePDMsString, regexstatestosavepdm) ) {
        std::smatch SOI;
        std::regex_search(savePDMsString, SOI, regexstatestosavepdm);
        auto states_of_interest = SOI.str();
        states_of_interest.erase(std::remove(states_of_interest.begin(), states_of_interest.end(), '('), states_of_interest.end());
        states_of_interest.erase(std::remove(states_of_interest.begin(), states_of_interest.end(), ')'), states_of_interest.end());
        std::stringstream states_stream(states_of_interest);
        std::vector<size_t> result;
        posthf->saveOnePDM_states.resize(0);
        while( states_stream.good() )
        {
            std::string substr;
            std::getline( states_stream, substr, ',' );
            if (substr.find('-') != std::string::npos){
              std::stringstream range_stream(substr);
              std::string start_string, end_string;
              std::getline( range_stream, start_string, '-' );
              std::getline( range_stream, end_string, '-' );
              std::istringstream start_ss(start_string);
              std::istringstream end_ss(end_string);
              size_t start, end;
              start_ss >> start;
              end_ss >> end;
              start -= 1;
              end -= 1;
              for (auto i = start; i < end+1; i++){
                posthf->saveOnePDM_states.push_back( i );
              }
            } else {
              std::istringstream orb_ss(substr);
              size_t orb;
              orb_ss >> orb;
              posthf->saveOnePDM_states.push_back(orb-1);
            }
        }
        std::sort(posthf->saveOnePDM_states.begin(), posthf->saveOnePDM_states.end());
        //posthf->noSOI = std::stoi(SOI.str());
      }
    }
  }



}; // namespace ChronusQ

