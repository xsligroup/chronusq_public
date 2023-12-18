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
#include <morspec.hpp>
#include <cerr.hpp>

namespace ChronusQ {

  /**
   *
   *  Check valid keywords in the section.
   *
  */
  void CQRESPONSE_VALID( std::ostream &out, CQInputFile &input ) {


    if( not input.containsSection("RESPONSE") ) return;

    // Allowed keywords
    std::vector<std::string> allowedKeywords = {
      "TYPE",
      "PROPAGATOR",
      "CONV",
      "MAXITER",
      "FULLMAT",
      "DOFULL",
      "TDA",
      "DISTMATFROMROOT",
      "FORMMATDIST",
      "DAMP",
      "FORCEDAMP",
      "NROOTS",
      "DEMIN",
      "GPLHR_M",
      "GPLHR_SIGMA",
      "DOAPBAMB",
      "DOREDUCED",
      "DOSTAB",
      "PPSPINMAT",
      "PPTDAMAT",
      "PPSTAR",
      "AOPS",
      "BOPS",
      "BFREQ",
      "NEO",
    };

    // Specified keywords
    std::vector<std::string> responseKeywords = input.getDataInSection("RESPONSE");

    // Make sure all of basisKeywords in allowedKeywords
    for( auto &keyword : responseKeywords ) {
      auto ipos = std::find(allowedKeywords.begin(),allowedKeywords.end(),keyword);
      if( ipos == allowedKeywords.end() ) 
        CErr("Keyword RESPONSE." + keyword + " is not recognized",std::cout);// Error
    }



    // Check for disallowed combinations (if any)

/*
    // No GIAO + RESPONSE
    if( input.containsData("BASIS.BASISTYPE") ) {

      std::string btype = 
        input.getData<std::string>("BASIS.BASISTYPE");

      if( not btype.compare("GIAO") )
        CErr("GIAO + RESPONSE not allowed");

    } 
*/

#if 0
    // No MPI + Full Residue
    if( MPISize(MPI_COMM_WORLD) > 1 ) {

      bool isRes = false;

      // WARNING: relies on default to RESIDUE if
      // TYPE not specified
      if( input.containsData("RESPONSE.TYPE") ) {

        std::string type = 
          input.getData<std::string>("RESPONSE.TYPE");

        isRes = not type.compare("RESIDUE");

      } else isRes = true;

      bool isDist = false;

      if( input.containsData("RESPONSE.FORMMATDIST") )
        isDist = isDist or 
                   input.getData<bool>("RESPONSE.FORMMATDIST");

      if( input.containsData("RESPONSE.DISTMATFROMROOT") )
        isDist = isDist or 
                   input.getData<bool>("RESPONSE.DISTMATFROMROOT");

      if( isRes and isDist ) 
        CErr("FULL RESIDUE + MPI not allowed");

    }
#endif


    // No GPLHR + A+B/A-B/REDUCED
    if( input.containsData("RESPONSE.DOAPBAMB") or
        input.containsData("RESPONSE.DOREDUCED") ) {


      bool doAPB = false;
      bool doRed = false;

      OPTOPT( doAPB = input.getData<bool>("RESPONSE.DOAPBAMB"); )
      OPTOPT( doRed = input.getData<bool>("RESPONSE.DOREDUCED"); )


      bool isRes = false;

      // WARNING: relies on default to RESIDUE if
      // TYPE not specified
      if( input.containsData("RESPONSE.TYPE") ) {

        std::string type = 
          input.getData<std::string>("RESPONSE.TYPE");

        isRes = not type.compare("RESIDUE");

      } else isRes = true;


      bool doFull = true; // WARNING: default doFULL
      if( input.containsData("RESPONSE.DOFULL") ) {

        doFull = input.getData<bool>("RESPONSE.DOFULL");

      }


      if( (isRes and not doFull) and (doAPB or doRed) )
        CErr("GPLHR + non full dimensional not allowed");

    }

  }

  /**
   *
   *  Check valid keywords in the section.
   *
  */
  void CQMOR_VALID( std::ostream &out, CQInputFile &input ) {

    // Allowed keywords
    std::vector<std::string> allowedKeywords = {
      "NMODEL",
      "REFINE",
      "NMODELMAX",
      "GETEIG",
      "ERRMETH",
    };

    // Specified keywords
    std::vector<std::string> morKeywords = input.getDataInSection("MOR");

    // Make sure all of basisKeywords in allowedKeywords
    for( auto &keyword : morKeywords ) {
      auto ipos = std::find(allowedKeywords.begin(),allowedKeywords.end(),keyword);
      if( ipos == allowedKeywords.end() ) 
        CErr("Keyword MOR." + keyword + " is not recognized",std::cout);// Error
    }
    // Check for disallowed combinations (if any)
  }

  /**
   *  \brief Construct a Response object using the input 
   *  file.
   *
   *  \param [in] out    Output device for data / error output.
   *  \param [in] input  Input file datastructure
   *  \param [in] ss     SingleSlater reference
   *                     
   *
   *  \returns shared_ptr to a ResponseBase object
   *    constructed from the input options.
   *
   */ 
  std::shared_ptr<ResponseBase> CQResponseOptions(std::ostream &out, 
    CQInputFile &input, std::shared_ptr<SingleSlaterBase> &ss,
    EMPerturbation &scfPert ) {


    if( not input.containsSection("RESPONSE") )
      CErr("RESPONSE Section must be specified for RESPONSE job",out);


    out << "  *** Parsing RESPONSE options ***\n";

    std::shared_ptr<ResponseBase> resp     = nullptr;
    SingleSlaterPolarBase* factory  = nullptr;
    SingleSlaterParticleParticleBase* factoryPP  = nullptr;
    std::shared_ptr<MORSpecBase>  mor      = nullptr;


    // Determine  reference and construct RESPONSE object
      

    ResponseType jobTyp = RESIDUE;
    bool doMOR = false;
    bool doPP  = false;
    bool doNEO = false;

    std::string jt = "RESIDUE";

    OPTOPT( jt = input.getData<std::string>("RESPONSE.TYPE"); );

    if( not jt.compare("RESIDUE") )   jobTyp = RESIDUE;
    else if( not jt.compare("FDR") )  jobTyp = FDR;
    else if( not jt.compare("MOR") )  doMOR = true;
    else CErr(jt + " NOT RECOGNIZED RESPONSE.TYPE");

		if ( input.containsData("RESPONSE.NEO") ) doNEO = input.getData<bool>("RESPONSE.NEO");

    // Determine propagator
    try {

      std::string prop = input.getData<std::string>("RESPONSE.PROPAGATOR");


      std::vector< std::string > availProp = 
      {
        "PARTICLEHOLE",
        "PARTICLEPARTICLE"
      };

      bool goodProp = std::any_of(availProp.begin(),availProp.end(),
          [&](std::string &x){ return not x.compare(prop); });

      if( not goodProp )
        CErr(prop + " is not a valid RESPONSE.PROPAGATOR",out);

      doPP = not prop.compare("PARTICLEPARTICLE");

    } catch(...){ }



    if( doMOR and doPP )
      CErr("MOR + ParticleParticlePropagator not valid!",out);


		if (doNEO and doPP )
			CErr("NEO + ParticleParticlePropagator NYI!",out);

		if (doNEO and doMOR )
			CErr("NEO + MOR NYI!",out);


    bool found = false;

    #define CONSTRUCT_PH_RESP(_CLASS) \
      if(not found) try { \
        auto &tmp = dynamic_cast<_CLASS&>(*ss);\
        resp = std::dynamic_pointer_cast<ResponseBase>(\
          std::make_shared<PolarizationPropagator<_CLASS>>(MPI_COMM_WORLD,jobTyp,\
            std::dynamic_pointer_cast<_CLASS>(ss)\
          )\
        );\
        if( std::is_base_of<SingleSlaterPolarBase,PolarizationPropagator<_CLASS>>::value )\
          factory = &dynamic_cast<PolarizationPropagator<_CLASS>&>(*resp);\
        found = true; \
      } catch(...) { }

    #define CONSTRUCT_PP_RESP(_CLASS) \
      if(not found) try { \
        auto &tmp = dynamic_cast<_CLASS&>(*ss);\
        resp = std::dynamic_pointer_cast<ResponseBase>(\
          std::make_shared<ParticleParticlePropagator<_CLASS>>(MPI_COMM_WORLD,jobTyp,\
            std::dynamic_pointer_cast<_CLASS>(ss)\
          )\
        );\
        if( std::is_base_of<SingleSlaterParticleParticleBase,ParticleParticlePropagator<_CLASS>>::value )\
          factoryPP = &dynamic_cast<ParticleParticlePropagator<_CLASS>&>(*resp);\
        found = true; \
      } catch(...) { }

    #define CONSTRUCT_MOR(_CLASS) \
      if(not found) try { \
        auto &tmp = dynamic_cast<_CLASS&>(*ss);\
        mor = std::static_pointer_cast<MORSpecBase>(\
          std::make_shared<MORSpec<_CLASS>>(\
            MPI_COMM_WORLD,\
            std::dynamic_pointer_cast<_CLASS>(ss)\
          )\
        );\
        resp = std::dynamic_pointer_cast<ResponseBase>(\
          std::static_pointer_cast<MORSpec<_CLASS>>(mor)\
        );\
        if( std::is_base_of<SingleSlaterPolarBase,PolarizationPropagator<_CLASS>>::value )\
          factory = & (std::static_pointer_cast<MORSpec<_CLASS>>(mor)->respFactory());\
        found = true; \
      } catch(...) { }



    using HF_dd = HartreeFock<double  ,double>;
    using HF_cd = HartreeFock<dcomplex,double>;
    using HF_cc = HartreeFock<dcomplex,dcomplex>;
    using KS_dd = KohnSham<double  ,double>;
    using KS_cd = KohnSham<dcomplex,double>;
    using KS_cc = KohnSham<dcomplex,dcomplex>;
    using NEO_dd = NEOSS<double  ,double>;
    using NEO_cd = NEOSS<dcomplex,double>;
    using NEO_cc = NEOSS<dcomplex,dcomplex>;


    if( doMOR ) {

      CONSTRUCT_MOR( HF_dd );
      CONSTRUCT_MOR( HF_cd );
      CONSTRUCT_MOR( HF_cc );
      CONSTRUCT_MOR( KS_dd );
      CONSTRUCT_MOR( KS_cd );
      CONSTRUCT_MOR( KS_cc );

    } else if( doPP ) {

      CONSTRUCT_PP_RESP( HF_dd );
      CONSTRUCT_PP_RESP( HF_cd );
      CONSTRUCT_PP_RESP( HF_cc );
      CONSTRUCT_PP_RESP( KS_dd );
      CONSTRUCT_PP_RESP( KS_cd );
      CONSTRUCT_PP_RESP( KS_cc );

    } else {

      CONSTRUCT_PH_RESP( HF_dd );
      CONSTRUCT_PH_RESP( HF_cd );
      CONSTRUCT_PH_RESP( HF_cc );
      CONSTRUCT_PH_RESP( KS_dd );
      CONSTRUCT_PH_RESP( KS_cd );
      CONSTRUCT_PH_RESP( KS_cc );
    if(doNEO){
				CONSTRUCT_PH_RESP( NEO_dd );
      	CONSTRUCT_PH_RESP( NEO_cd );
      	CONSTRUCT_PH_RESP( NEO_cc );
			}
   
    }

    // Copy over SCF Perturbation
    resp->scfPert = scfPert;

    // General settings

    OPTOPT( resp->genSettings.convCrit = 
              input.getData<double>("RESPONSE.CONV") );
    OPTOPT( resp->genSettings.maxIter = 
              input.getData<double>("RESPONSE.MAXITER") );
    OPTOPT( resp->genSettings.formFullMat = 
              input.getData<bool>("RESPONSE.FULLMAT") );
    OPTOPT( resp->genSettings.doFull = 
              input.getData<bool>("RESPONSE.DOFULL") );
    OPTOPT( resp->genSettings.doTDA = 
              input.getData<bool>("RESPONSE.TDA") );




    OPTOPT( resp->genSettings.distMatFromRoot = 
              input.getData<bool>("RESPONSE.DISTMATFROMROOT") );
    OPTOPT( resp->genSettings.formMatDist = 
              input.getData<bool>("RESPONSE.FORMMATDIST") );

    // FDR settings
      
    OPTOPT( resp->fdrSettings.dampFactor = 
              input.getData<double>("RESPONSE.DAMP") );
    OPTOPT( resp->fdrSettings.forceDamp = 
              input.getData<bool>("RESPONSE.FORCEDAMP") );


    // RESIDUE settings

    OPTOPT( resp->resSettings.nRoots =
              input.getData<size_t>("RESPONSE.NROOTS") );
    OPTOPT( resp->resSettings.deMin =
              input.getData<double>("RESPONSE.DEMIN") );
    OPTOPT( resp->resSettings.gplhr_m =
              input.getData<size_t>("RESPONSE.GPLHR_M") );
    OPTOPT( resp->resSettings.gplhr_sigma =
              input.getData<double>("RESPONSE.GPLHR_SIGMA") );
    OPTOPT( resp->resSettings.useGDiag =
              input.getData<bool>("RESPONSE.USEGDIAG") );



    // Handling consistency when applicable
    if( input.containsData("RESPONSE.DEMIN") and 
        not input.containsData("RESPONSE.GPLHR_SIGMA") )
      resp->resSettings.gplhr_sigma = resp->resSettings.deMin;

    if( input.containsData("RESPONSE.GPLHR_SIGMA") and 
        not input.containsData("RESPONSE.DEMIN") )
      resp->resSettings.deMin = resp->resSettings.gplhr_sigma;

    // Forces full diag if the number of roots requested gives a
    // subspace greater than or equal to half of the full problem dimension.
    if( input.containsData("RESPONSE.NROOTS") and 
        input.containsData("RESPONSE.DOFULL") ) {
      if( ((3 + resp->resSettings.gplhr_m) * resp->resSettings.nRoots) >= resp->getNSingleDim() / 2 ) {

        resp->genSettings.doFull = true;

        std::cout << " " << std::endl;
        std::cout << "  ** REQUESTED ITERATIVE SUBSPACE IS >= 1/2 OF THE FULL PROBLEM DIMENSION: " << std::endl;
        std::cout << "    * DEFAULTING TO FULL DIAGONALIZATION " << std::endl;
        std::cout << " " << std::endl;
      }
    }











    // Method Specific
    if( not doPP )

      // Particle-Hole propagator
      try {

        SingleSlaterPolarBase *hf = factory;
      
        OPTOPT( hf->doAPB_AMB = input.getData<bool>("RESPONSE.DOAPBAMB") );
        OPTOPT( hf->doReduced = input.getData<bool>("RESPONSE.DOREDUCED") );
        
        OPTOPT( hf->doStab = input.getData<bool>("RESPONSE.DOSTAB") );

      } catch(...){ }

    else 

      // Particle-Particle propagator
      try {

        SingleSlaterParticleParticleBase *hf = factoryPP;
        OPTOPT( 

          std::string mat = input.getData<std::string>("RESPONSE.PPSPINMAT");
          if( not mat.compare("AA") )
            hf->spinSepProp = PP_AA;
          else if( not mat.compare("AB") )
            hf->spinSepProp = PP_AB;
          else if( not mat.compare("BB") )
            hf->spinSepProp = PP_BB;
          else
            CErr("RESPONSE.PPSPINMAT not recognized");

        );

        OPTOPT( 

          std::string mat = input.getData<std::string>("RESPONSE.PPTDAMAT");
          if( not mat.compare("A") )
            hf->tdaOp = PP_A;
          else if( not mat.compare("C") )
            hf->tdaOp = PP_C;
          else
            CErr("RESPONSE.PPTDAMAT not recognized");

        );

        OPTOPT( hf->doStarRef = input.getData<bool>("RESPONSE.PPSTAR") );

      } catch(...){ }


    // MOR options
    if( doMOR ) {

      OPTOPT( mor->morSettings.nModel = input.getData<size_t>("MOR.NMODEL") );
      OPTOPT( mor->morSettings.doRefine = input.getData<bool>("MOR.REFINE") );
      OPTOPT( mor->morSettings.nModelMax 
          = input.getData<size_t>("MOR.NMODELMAX") );

      OPTOPT( mor->morSettings.getEig = input.getData<bool>("MOR.GETEIG") );

      OPTOPT( 
        std::string errMeth = input.getData<std::string>("MOR.ERRMETH");
        mor->morSettings.doRelErr = not errMeth.compare("RELATIVE");
      );
    }

    // Handle Operator specification
    try {

      std::string aOpsStr, bOpsStr, bFreqStr;
      std::vector<std::string> aOps, bOps, bFreq;

      OPTOPT( aOpsStr  = input.getData<std::string>("RESPONSE.AOPS")  );
      OPTOPT( bOpsStr  = input.getData<std::string>("RESPONSE.BOPS")  );
      OPTOPT( bFreqStr = input.getData<std::string>("RESPONSE.BFREQ") );

      if(not aOpsStr.empty())  split(aOps ,aOpsStr ," ,;\t");
      if(not bOpsStr.empty())  split(bOps ,bOpsStr ," ,;\t");

      if( aOps.empty() and bOps.empty() and bFreqStr.empty()) 
        throw std::runtime_error("ISSUES");

      auto genList = [&](std::vector<std::string> ops) {

        std::vector<ResponseOperator> opsList;

        for(auto &op : ops) {
          trim(op);

          if(not op.compare("EDL")) opsList.push_back(LenElectricDipole);
          if(not op.compare("EQL")) opsList.push_back(LenElectricQuadrupole);
          if(not op.compare("EOL")) opsList.push_back(LenElectricOctupole);

          if(not op.compare("EDV")) opsList.push_back(VelElectricDipole);
          if(not op.compare("EQV")) opsList.push_back(VelElectricQuadrupole);
          if(not op.compare("EOV")) opsList.push_back(VelElectricOctupole);

          if(not op.compare("MD")) opsList.push_back(MagneticDipole);
          if(not op.compare("MQ")) opsList.push_back(MagneticQuadrupole);

          if(not op.compare("ALL")) 
            std::copy(AllOps.begin(),AllOps.end(),std::back_inserter(opsList));

          if(not op.compare("HEROPS"))
            std::copy(HerOps.begin(),HerOps.end(),std::back_inserter(opsList));

          if(not op.compare("ANTIHEROPS"))
            std::copy(AntiHerOps.begin(),AntiHerOps.end(),
              std::back_inserter(opsList));
        }

        return opsList;
      };

      if(not aOps.empty() ) resp->genSettings.aOps = genList(aOps);
      if(not bOps.empty() ) resp->genSettings.bOps = genList(bOps);

      if(not bFreqStr.empty() ) {

        resp->fdrSettings.bFreq.clear();


        if(bFreqStr.find("RANGE") != std::string::npos) {

          auto parStPos = bFreqStr.find("(");
          bFreqStr = bFreqStr.substr(parStPos+1);

          auto parEnPos = bFreqStr.find(")");
          bFreqStr = bFreqStr.substr(0,parEnPos);

          split(bFreq,bFreqStr,",");

          if( bFreq.size() != 3 ) CErr();

          double start = std::stod(bFreq[0]);
          int    count = std::stoi(bFreq[1]);
          double step  = std::stod(bFreq[2]);

          for(auto i = 0; i < count; i++)
            resp->fdrSettings.bFreq.push_back(start + i*step);

        } else {
          split(bFreq ,bFreqStr ," ,;\t");

          for(auto &freq : bFreq) { 
            trim(freq); resp->fdrSettings.bFreq.push_back(std::stod(freq));
          }
        }
      }
      

    } catch(...) { }

    return resp;
  }

};
