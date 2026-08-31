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

#include <chronusq_sys.hpp>
#include <regex>

namespace ChronusQ {

  /**
   * \brief A class to compare two strings in lexicographic order with exceptions:
   *        1. Slash and Dot are before any other character
   *        2. Numbers in brackets are compared before any other character
   *        3. Numbers in brackets are compared numerically
   */
  struct InputKeyCompare {
    static size_t extractNumber(const std::string& str, size_t index) {
      size_t num = 0;
      while (index < str.size() and str[index] != ']') {
        if (std::isdigit(str[index])) {
          num = num * 10 + (str[index] - '0');
        } else {
          return std::numeric_limits<size_t>::max();
        }
        ++index;
      }
      return num;
    }

    bool operator()(const std::string& a, const std::string& b) const {
      size_t minSize = std::min(a.size(), b.size());
      for (size_t i = 0; i < minSize; i++) {
        if (a[i] == '[' && b[i] == '[') {
          size_t numA = extractNumber(a, i+1);
          size_t numB = extractNumber(b, i+1);
          if (numA != numB) return numA < numB;
        } else {
          if (a[i] != b[i]) {
            if (a[i] == '/') return true;
            if (b[i] == '/') return false;
            if (a[i] == '.') return true;
            if (b[i] == '.') return false;
            if (a[i] == '[') return true;
            if (b[i] == '[') return false;
            return a[i] < b[i];
          }
        }
      }
      return a.size() < b.size();
    }
  };

  typedef std::map<std::string,std::map<std::string,std::string>,InputKeyCompare> InputMap;


  /**
   *  \brief A class to handle the parsing a data fetching from a 
   *  ChronusQ input file.
   */
  class CQInputFile {

    // friend function to overload the << operator for the CQInputFile class
    friend std::ostream& operator<<(std::ostream& os, const CQInputFile& inputFile);
  
    std::shared_ptr<std::ifstream> inFile_ = nullptr;  ///< Input file

    InputMap dict_;
    ///< Input data fields partitioned by section headings

    //Required Job Input:
    bool parseFreeCQInputJob(std::string&);
    bool parseFreeCQInputRef(std::string&);
    bool parseFreeCQInputBas(std::string&);
    //Additional input options:
    void parseFreeCQInput(std::string&);
    void parseFreeCQInputNEO(std::string&);
    //void parseFreeCQInputElectron(std::string&);
    void parseFreeCQInputSCF(std::string&);
    void parseFreeCQInputRT(std::string&, const std::regex &freeCQInputRT);
    void parseFreeCQInputField(std::string&);
    void parseFreeCQInputSSGuess(std::string&);
    void parseFreeCQInputCI(std::string&, const std::regex &freeCQInputCI);
    void parseFreeCQInputCC(std::string&, const std::regex &freeCQInputCC);



    /**
     * \brief Add an key-value pair to the input storage
     *  \param [in] prefix Prefix string for the section
     *  \param [in] key    Key of the data
     *  \param [in] path   Formatted query string to be parsed
     */
    void addData(const std::string &prefix, const std::string &key, const std::string &value);
    void addData(const std::string &path, const std::string &value);


    /**
     * \brief Merge a subsection into the input storage
     * \param [in] prefix  Prefix of the data field
     * \param [in] section Section to be merged
     */
    void mergeSection(const std::string &prefix, const std::map<std::string,std::string> &section);


    // Splits query string on last "/"
    // (See src/cxxapi/input/parse/cxxapi.cxx for documentation)
    static std::pair<std::string,std::string>
      splitQuery(const std::string&);
  
    /**
     *  std::ofstream constructor.
     *
     *  Sets and parses input file from std::ofstream object
     *  \param [in] inFile  File object to parse
     */  
    CQInputFile(std::shared_ptr<std::ifstream> inFile) :
      inFile_(inFile){ }
  
  
  
  public:
  
    // Disable default, copy and move constructors and assignment operators
    CQInputFile()                               = delete;
    CQInputFile(const CQInputFile &)            = delete;
    CQInputFile(CQInputFile &&)                 = delete;
    CQInputFile& operator=(const CQInputFile &) = delete; 
    CQInputFile& operator=(CQInputFile &&)      = delete;

    // Parses the input file
    // (See src/cxxapi/input/parse/cxxapi.cxx for documentation)
    void parse();

    // Parses a section of the input file
    // (See src/cxxapi/input/parse/cxxapi.cxx for documentation)
    void parse(std::vector<std::string>::const_iterator lines_begin,
               std::vector<std::string>::const_iterator lines_end,
               const std::string &prefix);
    /**
     *  Filename constructor.
     *
     *  Sets and parses CQ input file given a file name
     *  \param [in] inFileName  Name of CQ input file
     */ 
    CQInputFile(std::string inFileName) :
      CQInputFile(std::make_shared<std::ifstream>(inFileName)){ }
  
  
  
    /**
     *  \brief Template function which returns the value of a data field
     *  from the input file in a specified datatype given a formatted 
     *  query string.
     *
     *  i.e.
     *
     *  Input entry:
     *    [SCF]
     *    DENTOL = 1E-6
     *    
     *  Query
     *    double tol = input.getData<double>("SCF.DENTOL");
     *
     *  This example returns the value of the string data field "SCF.DENTOL"
     *  as a double precision number. Various specializations of this function
     *  exist for various datatypes
     *
     *  \param [in] prefix Prefix string for the section
     *  \param [in] key    Key of the data
     *  \param [in] path   Formatted query string to be parsed
     *  \return            Value of query data field as specified datatype
     */
    template <typename T> T getData(const std::string &prefix,
                                    const std::string &key) const;
    template <typename T> T getData(const std::string &path) const;
  
  
  
  
    /**
     *  Checks whether or not the parsed CQ input file contains
     *  a query section.
     *
     *  \paral  [in] str Query string of a section heading
     *  \return      True if input file contains that heading
     */ 
    bool containsSection(const std::string &prefix) const;
  
    /**
     *  Checks whether or not the parsed CQ input file contains
     *  a query data field.
     *
     *  \paral  [in] str Query string of a data field (includes section heading)
     *  \return      True if input file contains that data field
     */
    bool containsData(const std::string &prefix,
                      const std::string &key) const;
    bool containsData(const std::string &path) const;
  




    std::set<std::string> getDataInSection(const std::string &prefix) const;


    /**
     *  \brief Returns a subsection of data fields from the input file
     *
     *  \param [in] prefix  Section heading
     *  \return             Vector of data fields in section
     */
    const std::map<std::string,std::string>& getSection(const std::string &prefix) const;


    /**
     * \brief Add an key-value pair to the input storage
     *  \param [in] prefix Prefix string for the section
     *  \param [in] key    Key of the data
     *  \param [in] path   Formatted query string to be parsed
     */
    void modifyData(const std::string &prefix, const std::string &key, const std::string &value);
    void modifyData(const std::string &path, const std::string &value);

  }; // CQInputFile class

  /**
   *  \brief Overload the << operator for the CQInputFile class
   *
   *  \param [in] os        Output device for data / error output.
   *  \param [in] inputFile CQInputFile object to be printed
   *
   *  \returns std::ostream object
   */
  std::ostream& operator<<(std::ostream& os, const CQInputFile& inputFile);
  
  
  // Misc string functions
  
  /**
   *  Trim a string of left trailing whitespace
   *
   *  \param [in/out] s std::string to be trimmed
   */
  static inline std::string& trim_left(std::string &s) {
      s.erase(s.begin(), std::find_if_not(s.begin(), s.end(),
              [](auto& x){ return std::isspace(x); }));
      return s;
  }; // trim_left
  
  
  /**
   *  Trim a string of right trailing whitespace
   *
   *  \param [in/out] s std::string to be trimmed
   */
  static inline std::string& trim_right(std::string &s) {
      s.erase(std::find_if(s.rbegin(), s.rend(),
              [](auto& x){ return !std::isspace(x); }).base(), s.end());
      return s;
  }; // trim_right
  
  
  /**
   *  Trim a string of trailing whitespace from both ends
   *
   *  \param [in/out] s std::string to be trimmed
   */
  static inline std::string &trim(std::string &s) {
      return trim_left(trim_right(s));
  }; // trim
  
  /**
   *  Splits a string into tokens  based on a demiliter
   *
   *  \param [out] tokens     std::vector of std::string objects which hold
   *                          the split tokens
   *  \param [in]  str        std::string to split
   *  \param [in]  delimiters Delimiters on which to split str
   */
  static inline void split(std::vector<std::string>& tokens, 
    const std::string& str, const std::string& delimiters = " ") {
  
      tokens.clear();
      // Skip delimiters at beginning.
      std::string::size_type lastPos = str.find_first_not_of(delimiters, 0);
      // Find first "non-delimiter".
      std::string::size_type pos     = str.find_first_of(delimiters, lastPos);
  
      while (std::string::npos != pos || std::string::npos != lastPos) {
          // Found a token, add it to the vector.
          tokens.push_back(str.substr(lastPos, pos - lastPos));
          // Skip delimiters.  Note the "not_of"
          lastPos = str.find_first_not_of(delimiters, pos);
          // Find next "non-delimiter"
          pos = str.find_first_of(delimiters, lastPos);
      }
  }; // split

  /**
   *  Reverse the order of tokens in a string separated by slashs
   *
   *  \param [in]  str        std::string to reverse
   */
  static inline std::string reverse_by_dot(const std::string& str) {

    std::vector<std::string> tokens;

    split(tokens,str,"/");

    // Combine tokens in reverse order by dot
    std::string reversed = "";
    for( auto it = tokens.rbegin(); it != tokens.rend(); ++it ) {
      reversed += *it;
      if( it != tokens.rend()-1 ) reversed += "/";
    }

    return reversed;
  }; // split

  template <typename T>
  T CQStringTo(const std::string &s);

}; // namespace ChronusQ

