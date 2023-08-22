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

#include <unordered_set>

#include <cxxapi/input.hpp>
#include <cerr.hpp>
#include <regex>

namespace ChronusQ {


  /**
   *  \brief Parses a ChronusQ input file
   *
   *  Parses the file and populates the dict_ map which holds the
   *  input data fields to control the ChronusQ calculation
   */
  void CQInputFile::parse() {

    // Check if file actually exists
    if(not inFile_->good()) CErr("Input File Couldn't Be Found!",std::cout);

    // Read in all lines of the file
    std::vector<std::string> lines;
    while( not inFile_->eof() ) {
      std::string line;
      std::getline(*inFile_,line);
      lines.push_back(line);
    }

    // Parse the file
    parse(lines.cbegin(), lines.cend(), "");

  }; // CQInputFile::parse



  /**
   *  \brief Parses a section of the input file
   *
   *  Parses the file and populates the dict_ map which holds the
   *  input data fields to control the ChronusQ calculation
   */
  void CQInputFile::parse(std::vector<std::string>::const_iterator lines_begin,
                          std::vector<std::string>::const_iterator lines_end,
                          const std::string &prefix) {

    bool parseSection(false);
    bool prevLineData(false);
    size_t prevIndent(0);
  
    std::string sectionHeader;
    std::string dataHeader;

    // Keywords that are case sensitive (do *not* transform data to UPPER)
    std::set<std::string> caseSens, caseSensReverse;

    // Add case sensitive data keywords here
    caseSens.insert("BASIS.BASIS");

    // Reverse entries in caseSens
    for (auto &sec : caseSens)
      caseSensReverse.insert(reverse_by_dot(sec));
  
    // Loop over all lines of the file
    for( auto line_iter = lines_begin; line_iter != lines_end; ++line_iter ) {

      std::string line = *line_iter;
  
      // Skip blank lines
      if(line.length() < 1) {
        prevLineData = false;
        continue;
      }
        
  
      // Determine position of first and last non-space character
      size_t firstNonSpace = line.find_first_not_of(" ");
      size_t lastNonSpace  = line.find_last_not_of(" ");
  
      size_t comPos = line.find("#");
  
      // Skip lines in which the first non-space character is #
      // (Comment line)
      if(comPos == firstNonSpace) continue;
  
      // Remove comment portion of the line if it exists
      //  - This is general to when # does not appear in the line
      line = line.substr(0,comPos); 
  
      // Strip trailing spaces
      trim_right(line);
  
  
      auto strToUpper = [](std::string& s){
        std::transform(s.begin(), s.end(), s.begin(),
        [](unsigned char c){ return std::toupper(c);});
      };


      // Check if we have a free-style CQ input
      parseFreeCQInput(line);

      size_t lBrckPos = line.find('[');
      size_t rBrckPos = line.find(']');
  
      size_t eqPos  = line.find('=');
      size_t colPos = line.find(':');
  
      // Determine if this is a line with a section header
      bool sectionLine = 
        lBrckPos == firstNonSpace and rBrckPos == lastNonSpace;
  
      // Determine if this is a line that contains a data field
      bool dataLine    = 
        eqPos != std::string::npos or 
        colPos != std::string::npos;
  
      // Determine if this is a line continuation of a previous data field
      bool multiLine   = prevLineData and firstNonSpace > prevIndent;
  
      // Section line
      if(sectionLine) {
  
        // Strip first spaces
        line = line.substr(firstNonSpace,line.length());
  
        // Obtain the section header name
        sectionHeader = line.substr(1,line.length()-2);
        
        // Convert to UPPER
        strToUpper(sectionHeader);
  
//        // Create a dictionary entry for the section header
//        dict_[sectionHeader] =
//          std::unordered_map<std::string,std::string>();
  
        // XXX: Possibly check if the section is already defined?
  
        parseSection = true;
        prevLineData = false;
        continue;
  
      }
  
  
      // Data line
      if(parseSection and dataLine) {
  
        line = 
          line.substr(firstNonSpace,line.length()-firstNonSpace);
  
        // Split the line into tokens, trim spaces
        std::vector<std::string> tokens;
        split(tokens,line,"=:");
        for(auto &X : tokens) { trim(X); }
  
        dataHeader = tokens[0];
        strToUpper(dataHeader);
        dataHeader = sectionHeader + "." + dataHeader;

        // Capitalize data if not case sensitive
        auto it = caseSensReverse.lower_bound(reverse_by_dot(dataHeader));
        if ((it == caseSensReverse.end()
              or it->find(reverse_by_dot(dataHeader)) != 0)
            and tokens.size() > 1)
          strToUpper(tokens[1]);

        // Create a dictionary entry for the data field in the current
        // section header
        if(tokens.size() > 1) 
          dict_[dataHeader] = tokens[1];
        else 
          dict_[dataHeader] = " ";
  
        prevLineData = true;
        prevIndent = firstNonSpace;
      }
  
      // Multiline data
      else if(parseSection and multiLine) {
        // Capitalize data if not case sensitive
        auto it = caseSensReverse.lower_bound(reverse_by_dot(dataHeader));
        if (it == caseSensReverse.end()
             or it->find(reverse_by_dot(dataHeader)) != 0)
          strToUpper(line);
 
        line = 
          line.substr(firstNonSpace,line.length()-firstNonSpace);
        dict_[dataHeader] += "\n" + line;
      }
      
    };
  
  /* Debug code which prints out the contents of the dict_ map
    for(auto &sec : dict_) {
      std::cout << "Section: " << sec.first << std::endl;
      for(auto &data : sec.second) {
        std::cout << "  DATA: " << data.first << " ; " << data.second << std::endl;
      }
    }
  */
  
  }; // CQInputFile::parse(lines)
  
  void CQInputFile::parseFreeCQInput (std::string &line){

    /********************************************************************************/
    /* CQ Free Format Input                                                         */
    /* example, CQ= HF/STO-3G NEO(EPC17/PROT-BP4-D) SCF(accuracy=1.e-6)             */
    /* example, ChronuQ: X2C-HF/CD-6-31G RT(time=10fs, stepsize = 1as)              */
    /* example, ChronuQ= 4C-HF/ano-rcc hamiltonian(DCB, scalar, atomic)             */
    /********************************************************************************/
    auto const freeCQInput = std::regex("CQ[[:blank:]]*=|CQ[[:blank:]]*:|CHRONUSQ[[:blank:]]*=|CHRONUSQ[[:blank:]]*:",std::regex_constants::icase);
    if(!std::regex_search(line, freeCQInput)) return;
    std::cout<<"xsli test CQ Input"<<std::endl;
    line = std::regex_replace(line, freeCQInput, "");

    // Parse NEO Section
    parseFreeCQInputNEO(line);
    parseFreeCQInputElectron(line);
    parseFreeCQInputSCF(line);

  }; // Free Format Input Parser



  void CQInputFile::parseFreeCQInputNEO (std::string &line){

    auto const freeCQInputHF    = std::regex("((2C)|(X2C)|(4C)|(G)-?)?HF/?",std::regex_constants::icase);
    auto const freeCQInputCCSD  = std::regex("((2C)|(X2C)|(4C)|(G)-?)?CCSD((T)|(\\(T\\)))?/?",std::regex_constants::icase);

    /*************************************/
    /* NEO Input                         */
    /* example, NEO(EPC19/prot-pb6-g)    */
    /* example, NEO(EPC19/CD-prot-pb6-g) */
    /* example, NEO(EPC19/ri-prot-pb6-g) */
    /*************************************/
    auto const freeCQInputEPC17 = std::regex("EPC17",std::regex_constants::icase);
    auto const freeCQInputEPC19 = std::regex("EPC19",std::regex_constants::icase);

    auto const freeCQInputPROTSP    = std::regex("((CD)|(RI)-?)?PROT-SP",std::regex_constants::icase);
    auto const freeCQInputPROTPB4D  = std::regex("((CD)|(RI)-?)?PROT-PB4-D",std::regex_constants::icase);
    auto const freeCQInputPROTPB4F1 = std::regex("((CD)|(RI)-?)?PROT-PB4-F1",std::regex_constants::icase);
    auto const freeCQInputPROTPB4F2 = std::regex("((CD)|(RI)-?)?PROT-PB4-F2",std::regex_constants::icase);
    auto const freeCQInputPROTPB5G  = std::regex("((CD)|(RI)-?)?PROT-PB5-G",std::regex_constants::icase);
    auto const freeCQInputPROTPB6G  = std::regex("((CD)|(RI)-?)?PROT-PB6-G",std::regex_constants::icase);

    auto const freeCQInputNEO = std::regex("NEO(\\((.*?)\\))?",std::regex_constants::icase);
    std::smatch NEOmatch;

    if( std::regex_search(line, NEOmatch, freeCQInputNEO) ){
      // Catch what is inside NEO()
      if(NEOmatch.str(2).size()>0) {
        // Parse user-defined input
        std::cout<<"xsli test NEO Section "<<std::endl;
        std::string NEOInputOptions = NEOmatch.str(2);

        // Methods
        if ( std::regex_search(NEOInputOptions, NEOmatch, freeCQInputHF) ) {
          std::cout<<"xsli test NEO HF"<<std::endl;
        }
        else if ( std::regex_search(NEOInputOptions, NEOmatch, freeCQInputEPC17) ) {
          std::cout<<"xsli test NEO EPC17"<<std::endl;
        }
        else if ( std::regex_search(NEOInputOptions, NEOmatch, freeCQInputEPC19) ) {
          std::cout<<"xsli test NEO EPC19"<<std::endl;
        }

        // Basis Sets
        if ( std::regex_search(NEOInputOptions, NEOmatch, freeCQInputPROTSP) ) {
          if( NEOmatch.str(1).size()==0 ) std::cout<<"xsli test NEO PORT-SP"<<std::endl;
          else if( NEOmatch.str(2).size()>0 ) std::cout<<"xsli test NEO CD-PORT-SP"<<std::endl;
          else if( NEOmatch.str(3).size()>0 ) std::cout<<"xsli test NEO RI-PORT-SP"<<std::endl;
        }
        else if ( std::regex_search(NEOInputOptions, NEOmatch, freeCQInputPROTPB4D) ) {
          if( NEOmatch.str(1).size()==0 ) std::cout<<"xsli test NEO PROT-PB4-D"<<std::endl;
          else if( NEOmatch.str(2).size()>0 ) std::cout<<"xsli test NEO CD-PROT-PB4-D"<<std::endl;
          else if( NEOmatch.str(3).size()>0 ) std::cout<<"xsli test NEO RI-PROT-PB4-D"<<std::endl;
        }
        else if ( std::regex_search(NEOInputOptions, NEOmatch, freeCQInputPROTPB4F1) ) {
          if( NEOmatch.str(1).size()==0 ) std::cout<<"xsli test NEO PROT-PB4-F1"<<std::endl;
          else if( NEOmatch.str(2).size()>0 ) std::cout<<"xsli test NEO CD-PROT-PB4-F1"<<std::endl;
          else if( NEOmatch.str(3).size()>0 ) std::cout<<"xsli test NEO RI-PROT-PB4-F1"<<std::endl;
        }
        else if ( std::regex_search(NEOInputOptions, NEOmatch, freeCQInputPROTPB4F2) ) {
          if( NEOmatch.str(1).size()==0 ) std::cout<<"xsli test NEO PROT-PB4-F2"<<std::endl;
          else if( NEOmatch.str(2).size()>0 ) std::cout<<"xsli test NEO CD-PROT-PB4-F2"<<std::endl;
          else if( NEOmatch.str(3).size()>0 ) std::cout<<"xsli test NEO RI-PROT-PB4-F2"<<std::endl;
        }
        else if ( std::regex_search(NEOInputOptions, NEOmatch, freeCQInputPROTPB5G) ) {
          if( NEOmatch.str(1).size()==0 ) std::cout<<"xsli test NEO PROT-PB4-D"<<std::endl;
          else if( NEOmatch.str(2).size()>0 ) std::cout<<"xsli test NEO CD-PROT-PB4-D"<<std::endl;
          else if( NEOmatch.str(3).size()>0 ) std::cout<<"xsli test NEO RI-PROT-PB4-D"<<std::endl;
        }
        else if ( std::regex_search(NEOInputOptions, NEOmatch, freeCQInputPROTPB6G) ) {
          if( NEOmatch.str(1).size()==0 ) std::cout<<"xsli test NEO PROT-PB5-G"<<std::endl;
          else if( NEOmatch.str(2).size()>0 ) std::cout<<"xsli test NEO CD-PROT-PB5-G"<<std::endl;
          else if( NEOmatch.str(3).size()>0 ) std::cout<<"xsli test NEO RI-PROT-PB5-G"<<std::endl;
        }

      } else {
        // Choose default parameters
      }
      // We need to delete the NEO section so that we can parse the electronic section properly
      line = std::regex_replace(line, freeCQInputNEO, "");
    } // NEO Input

  };

  void CQInputFile::parseFreeCQInputElectron (std::string &line){

    auto const freeCQInputHF    = std::regex("((2C)|(X2C)|(4C)|(G)-?)?HF/?",std::regex_constants::icase);
    // Match strings                             1  2      3       4     5         6 7      8
    auto const freeCQInputCCSD  = std::regex("((2C)|(X2C)|(4C)|(G)-?)?CCSD((T)|(\\(T\\)))?/?",std::regex_constants::icase);

    /****************************/
    /* Electronic Input         */
    /* example, B3LYP/6-31G     */
    /* example, 4C-B3LYP/6-31G  */
    /* example, X2C-B3LYP/6-31G */
    /* example, B3LYP/RI-6-31G  */
    /* example, B3LYP/CD-6-31G  */
    /****************************/

    auto const freeCQInputB3LYP = std::regex("((2C)|(X2C)|(4C)|(G)-?)?B3LYP/?",std::regex_constants::icase);
    auto const freeCQInputPBE   = std::regex("((2C)|(X2C)|(4C)|(G)-?)?PBE/?",std::regex_constants::icase);

    std::smatch methodMatch;
    if ( std::regex_search(line, methodMatch, freeCQInputHF) ) {
      if( methodMatch.str(1).size()==0 ) std::cout<< "xsli test HF " <<std::endl;
      if( methodMatch.str(2).size()>0 ) std::cout<< "xsli test HF type: " <<methodMatch.str(2)<<std::endl;
      if( methodMatch.str(3).size()>0 ) std::cout<< "xsli test HF type: " <<methodMatch.str(3)<<std::endl;
      if( methodMatch.str(4).size()>0 ) std::cout<< "xsli test HF type: " <<methodMatch.str(4)<<std::endl;
      if( methodMatch.str(5).size()>0 ) std::cout<< "xsli test HF type: " <<methodMatch.str(5)<<std::endl;
      line = std::regex_replace(line, freeCQInputHF, "");
    }
    else if ( std::regex_search(line, methodMatch, freeCQInputCCSD) ) {
      if( methodMatch.str(1).size()==0 ) std::cout<< "xsli test CCSD " <<std::endl;
      if( methodMatch.str(2).size()>0 ) std::cout<< "xsli test CCSD type: " <<methodMatch.str(2)<<std::endl;
      if( methodMatch.str(3).size()>0 ) std::cout<< "xsli test CCSD type: " <<methodMatch.str(3)<<std::endl;
      if( methodMatch.str(4).size()>0 ) std::cout<< "xsli test CCSD type: " <<methodMatch.str(4)<<std::endl;
      if( methodMatch.str(5).size()>0 ) std::cout<< "xsli test CCSD type: " <<methodMatch.str(5)<<std::endl;
      if( methodMatch.str(7).size()>0 ) std::cout<< "xsli test CCSDT "<<std::endl;
      if( methodMatch.str(8).size()>0 ) std::cout<< "xsli test CCSD(T) "<<std::endl;
      line = std::regex_replace(line, freeCQInputCCSD, "");
    }


    auto const freeCQInputSTO3G   = std::regex("((CD)|(RI)-?)?STO-3G",std::regex_constants::icase);
    auto const freeCQInput321G    = std::regex("((CD)|(RI)-?)?3-21G",std::regex_constants::icase);
    auto const freeCQInput631G    = std::regex("((CD)|(RI)-?)?6-31G",std::regex_constants::icase);
    auto const freeCQInput6311G   = std::regex("((CD)|(RI)-?)?6-311G",std::regex_constants::icase);

    std::smatch basisMatch;
    if ( std::regex_search(line, basisMatch, freeCQInputSTO3G) ) {
      if( basisMatch.str(1).size()==0 ) std::cout<<"xsli test STO-3G"<<std::endl;
      else if( basisMatch.str(2).size()>0 ) std::cout<<"xsli test CD-STO-3G"<<std::endl;
      else if( basisMatch.str(3).size()>0 ) std::cout<<"xsli test RI-STO-3G"<<std::endl;
      line = std::regex_replace(line, freeCQInputSTO3G, "");
    }
    else if ( std::regex_search(line, basisMatch, freeCQInput321G) ) {
      if( basisMatch.str(1).size()==0 ) std::cout<<"xsli test 3-21G"<<std::endl;
      else if( basisMatch.str(2).size()>0 ) std::cout<<"xsli test CD-3-21G"<<std::endl;
      else if( basisMatch.str(3).size()>0 ) std::cout<<"xsli test RI-3-21G"<<std::endl;
      line = std::regex_replace(line, freeCQInput321G, "");
    }
    else if ( std::regex_search(line, basisMatch, freeCQInput631G) ) {
      if( basisMatch.str(1).size()==0 ) std::cout<<"xsli test 6-31G"<<std::endl;
      else if( basisMatch.str(2).size()>0 ) std::cout<<"xsli test 6-31G"<<std::endl;
      else if( basisMatch.str(3).size()>0 ) std::cout<<"xsli test 6-31G"<<std::endl;
      line = std::regex_replace(line, freeCQInput631G, "");
    }
    else if ( std::regex_search(line, basisMatch, freeCQInput6311G) ) {
      if( basisMatch.str(1).size()==0 ) std::cout<<"xsli test 6-311G"<<std::endl;
      else if( basisMatch.str(2).size()>0 ) std::cout<<"xsli test 6-311GG"<<std::endl;
      else if( basisMatch.str(3).size()>0 ) std::cout<<"xsli test 6-311G"<<std::endl;
      line = std::regex_replace(line, freeCQInput6311G, "");
    }

    auto const freeDividers = std::regex("\\s+|,+",std::regex_constants::icase);
    line = std::regex_replace(line, freeDividers, " ");
    // the second call remove multiple space left after replace ','
    line = std::regex_replace(line, freeDividers, " ");
    if(line.size()>0) std::cout<<"CQ input ignored: "<< line <<std::endl;

  };

  void CQInputFile::parseFreeCQInputSCF (std::string &line) {

    /********************************/
    /* SCF Input                    */
    /* example, SCF(accuracy=1.e-8) */
    /* example, SCF(endiis)         */
    /* example, SCF(energyonly)     */
    /********************************/

  };

  /** 
   *  \brief Splits a query string on a period "."
   * 
   *  This is a helpder function for the getData function which takes a 
   *  formatted string and splits it into a section and data field.
   *
   *  i.e.  "QM.REFERENCE" -> { "QM", "REFERENCE" }
   *
   *  \param [in] query Query string to be split
   *  \return     std::pair containing the two fields separated by a "."
   */
  std::pair<std::string,std::string> CQInputFile::splitQuery(
    const std::string &query) {
  
    std::vector<std::string> tokens;
  
    // Make sure that the query contains a period
  //assert( query.find(".") != query.end() );
  
    split(tokens,query,".");
    for(auto &X : tokens) {
      trim(X);
      std::transform(X.begin(),X.end(),X.begin(),
        [](unsigned char c){ return std::toupper(c);} );
    }
  
    return 
      std::pair<std::string,std::string>(tokens[0],tokens[1]);

  }; // CQInputFile::splitQuery

  
  /**
   *  \brief Custom exception type for handeling the case when
   *  a data field is not found for a query
   */
  class data_not_found : public std::exception {
  
    std::string msg; ///< Error message
  
  public:
  
    // Disable default constructor
    data_not_found() = delete;
  
    /**
     *  Exception constructor. Creates a useful error message
     *  which specifies the failed query
     */ 
    data_not_found(std::string x) { 
      msg = "Data ";
      msg += x; 
      msg += " Not Found\n";
    };
  
    /**
     *  Specialization of std::exception::what. Outputs the error message
     */ 
    virtual const char* what() const throw() {
      return msg.c_str();
    }
  
  }; // data_not_found class
  
  
  /**
   *  \brief Specialization of getData to return std::string of query 
   *  data field
   *
   *  \param [in] query Formatted query string to be parsed
   *  \return     Value of query data field as a std::string
   */
  template<>
  std::string CQInputFile::getData(std::string query) {
      auto kv = dict_.find(query);
  
      if(kv != dict_.end())
        return kv->second;

      else throw data_not_found(query);
  
  }; // CQInputFile::getData<std::string>
  
  /**
   *  \brief Specialization of getData to return int of query 
   *  data field
   *
   *  \param [in] query Formatted query string to be parsed
   *  \return     Value of query data field as a int
   */
  template<>
  int CQInputFile::getData(std::string query) {
  
    return std::stoi(getData<std::string>(query));
  
  }; // CQInputFile::getData<int>
  
  /**
   *  \brief Specialization of getData to return bool of query 
   *  data field
   *
   *  \param [in] query Formatted query string to be parsed
   *  \return     Value of query data field as a bool
   */
  template<>
  bool CQInputFile::getData(std::string query) {
  
    query = getData<std::string>(query);
    if (not query.compare("TRUE") or not query.compare("ON")){
      return true;
    }
      
    if (not query.compare("FALSE") or not query.compare("OFF")){
      return false;
    }
  
    CErr("Invalid Input For Boolean-Type Keyword!");

    return false;
  }; // CQInputFile::getData<bool>
  
  /**
   *  \brief Specialization of getData to return size_t of query 
   *  data field
   *
   *  \param [in] query Formatted query string to be parsed
   *  \return     Value of query data field as a size_t
   */
  template<>
  size_t CQInputFile::getData(std::string query) {
  
    return std::stoul(getData<std::string>(query));
  
  }; // CQInputFile::getData<size_t>
  
  /**
   *  \brief Specialization of getData to return double of query 
   *  data field
   *
   *  \param [in] query Formatted query string to be parsed
   *  \return     Value of query data field as a double
   */
  template<>
  double CQInputFile::getData(std::string query) {
  
    return std::stod(getData<std::string>(query));
  
  }; // CQInputFile::getData<double>

}; // namespace ChronusQ

