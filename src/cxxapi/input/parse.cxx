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
#include <stack>

namespace ChronusQ {


  enum class InputLineType {
    SECTION_HEADER,
    DATA_ENTRY,
    CONTINUATION,
    EMPTY
  };

  bool containsUnenclosedEqualSign(const std::string &s) {
    std::stack<char> st;

    for (char c : s) {
      if (c == '(' or c == '[' or c == '{') {
        st.push(c);
      } else if (c == ')' or c == ']' or c == '}') {
        if (st.empty()) {
          // unmatched closing bracket
          CErr("Unmatched closing bracket in input file line:\n" + s, std::cout);
        }
        char top = st.top();
        st.pop();
        if ((c == ')' and top != '(')
            or (c == ']' and top != '[')
            or (c == '}' and top != '{'))
          CErr("Unmatched bracket in input file line:\n" + s, std::cout);
      } else if (c == '=' or c == ':') {
        if (st.empty()) {
          // unenclosed '=' or ':' sign
          return true;
        }
      }
    }

    return false; // If we don't find an unenclosed '=' by the end
  }

  InputLineType get_input_line_type_and_trim(std::string &line) {

    // Determine position of first and last non-space character
    size_t firstNonSpace = line.find_first_not_of(" ");
    size_t lastNonSpace  = line.find_last_not_of(" ");

    size_t comPos = line.find("#");

    // Skip lines in which the first non-space character is #
    // (Comment line)
    if(comPos == firstNonSpace) return InputLineType::EMPTY;

    // Remove comment portion of the line if it exists
    //  - This is general to when # does not appear in the line
    line = line.substr(0,comPos);

    // Strip trailing spaces
    trim_right(line);

    size_t lBrckPos = line.find('[');
    size_t rBrckPos = line.find(']');

    // Check if we have a section header
    if( lBrckPos == firstNonSpace and rBrckPos == lastNonSpace )
      return InputLineType::SECTION_HEADER;

    // Check if we have a data entry
    if( containsUnenclosedEqualSign(line) ) {
      return InputLineType::DATA_ENTRY;
    }

    // If we get here, we have a continuation line
    return InputLineType::CONTINUATION;

  }; // get_input_line_type_and_trim


  /**
   *  \brief Parses a section of the input file
   *
   *  Parses the file and populates the dict_ map which holds the
   *  input data fields to control the ChronusQ calculation
   */
  void CQInputFile::parse(std::vector<std::string>::const_iterator lines_begin,
                          std::vector<std::string>::const_iterator lines_end,
                          const std::string &prefix) {


    auto strToUpper = [](std::string& s){
      std::for_each(s.begin(), s.end(),
                    [](char &c){ c = std::toupper(c);});
    };
  
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

      InputLineType lineType = get_input_line_type_and_trim(line);

      // Skip empty lines
      if (lineType == InputLineType::EMPTY) continue;

      // Section line
      if(lineType == InputLineType::SECTION_HEADER) {
  
        // Obtain the section header name
        sectionHeader = line.substr(1,line.length()-2);
        
        // Convert to UPPER
        strToUpper(sectionHeader);

        continue;
  
      }
  
  
      // Data line
      if(lineType == InputLineType::DATA_ENTRY) {

        // Find first = or : and get substring before and after it
        // TODO: make sure = or : is not enclosed by brackets
        size_t equalIndex = line.find_first_of("=:");
        dataHeader = line.substr(0,equalIndex);
        trim(dataHeader);
        std::string value = line.substr(equalIndex+1);
        trim(value);

        strToUpper(dataHeader);
        if (sectionHeader != "")
          dataHeader = sectionHeader + "." + dataHeader;

        // Check if the data entry has continuation lines below
        while (line_iter + 1 != lines_end) {
          std::string next_line = *(line_iter + 1);
          InputLineType next_line_type = get_input_line_type_and_trim(next_line);

          // End while loop if next line is not a continuation line or empty
          if (next_line_type != InputLineType::CONTINUATION
              and next_line_type != InputLineType::EMPTY)
            break;

          if (next_line_type == InputLineType::CONTINUATION)
            value += "\n" + next_line;
          ++line_iter;
        }

        // Capitalize data if not case sensitive
        auto it = caseSensReverse.lower_bound(reverse_by_dot(dataHeader));
        if ((it == caseSensReverse.end()
              or it->find(reverse_by_dot(dataHeader)) != 0)
            and not value.empty())
          strToUpper(value);

        // Create a dictionary entry for the data field in the current
        // section header
        if(not value.empty())
          addData(dataHeader,value);
        else
          CErr("No data entry for " + dataHeader + " in input file.");

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
   * \brief Add an key-value pair to the input storage
   * \param [in] key   Key of the data field
   * \param [in] value Value of the data field
   */
  void CQInputFile::addData(const std::string &key, const std::string &value) {
    if (containsData(key))
      CErr("Key " + key + " already exists in the parsed input.", std::cout);
    dict_[key] = value;
  }


  /**
   * \brief Merge a subsection into the input storage
   * \param [in] subsection Subsection to be merged
   * \param [in] prefix Prefix of the data field
   */
  void CQInputFile::mergeSection(const InputMap &subsection,
                                 const std::string &prefix) {
    if (prefix == "")
      for (auto &kv : subsection) {
        addData(kv.first, kv.second);
      }
    else
      for (auto &kv : subsection) {
        addData(prefix + "." + kv.first, kv.second);
      }
  }


  /**
   *  Checks whether or not the parsed CQ input file contains
   *  a query section.
   *
   *  \paral  [in] str Query string of a section heading
   *  \return      True if input file contains that heading
   */
  bool CQInputFile::containsSection(const std::string &str) const {
    auto it = dict_.lower_bound(str);
    if (it == dict_.end()) return false;

    // Check if the query string is identical section heading, this case is a data entry instead of section
    if (it->first == str) it++;
    if (it == dict_.end()) return false;
    // Check if the query string is a substring of the section heading
    return it->first.find(str) == 0 and it->first.size() > str.size() and it->first[str.size()] == '.';
  }


  /**
   *  Checks whether or not the parsed CQ input file contains
   *  a query section.
   *
   *  \paral  [in] str Query string of a section heading
   *  \return      True if input file contains that heading
   */
  bool CQInputFile::containsList(const std::string &str) const {
    auto it = dict_.lower_bound(str);
    if (it == dict_.end()) return false;
    // Check if the query string is a substring of the section heading
    while (it->first.find(str) == 0) {
      if (it->first.size() > str.size() and it->first[str.size()] == '[')
        return true;
      it++;
    }
    return false;
  }


  /**
   *  Checks the size of a query list.
   *
   *  \paral  [in] str Query string of a section heading
   *  \return      Size of the query list.
   */
  size_t CQInputFile::getListSize(const std::string &str) const {
    if (not containsList(str))
      return 0;
    size_t max_index = 0;
    auto it = dict_.lower_bound(str);
    // Check if the query string is a substring of the section heading
    while (it->first.find(str) == 0) {
      if (it->first.size() > str.size() and it->first[str.size()] == '[')
        max_index = std::max(max_index, InputKeyCompare::extractNumber(it->first, str.size() + 1));
      it++;
    }
    return max_index + 1;
  }

  /**
   * \brief Add a data to the end of a list
   * @param str  Query string of a section heading
   * @param data Data to be added to the list
   * @param dict InputMap to be appended
   */
  void CQInputFile::appendList(const std::string &str, const std::string &data) {
    addData(str + "[" + std::to_string(getListSize(str)) + "]", data);
  }
  void CQInputFile::appendList(const std::string &str, const InputMap &dict) {
    mergeSection(dict, str + "[" + std::to_string(getListSize(str)) + "]");
  }


  /**
   *  \brief Returns a subsection of data fields from the input file
   *
   *  \param [in] section Section heading
   *  \return             Vector of data fields in section
   */
  InputMap CQInputFile::getSection(const std::string &section) const {

    if (not containsSection(section))
      CErr("Section " + section + " not found in input file!");

    InputMap sectionData;

    auto it = dict_.lower_bound(section);

    // Check if the query string is identical section heading, this case is a data entry instead of section
    if (it->first == section) it++;
    while (it != dict_.end()) {
      const std::string &key = it->first;
      if (key.find(section) != 0 or key[section.size()] != '.') break;

      if (key.size() > section.size())
        sectionData.emplace(key.substr(section.size() + 1), it->second);

      ++it;
    }

    return sectionData;
  }

  /**
   *  Checks whether or not the parsed CQ input file contains
   *  a query data field.
   *
   *  \paral  [in] str Query string of a data field (includes section heading)
   *  \return      True if input file contains that data field
   */
  bool CQInputFile::containsData(std::string str) const {
    return dict_.find(str) != dict_.end();
  }


  std::vector<std::string> CQInputFile::getDataInSection( std::string section ) const  {

    std::set<std::string> datasets;

    std::string::size_type lenSection = section.size();
    auto it = dict_.lower_bound(section);
    while (it != dict_.end()) {
      const std::string &key = it->first;
      if (key.find(section) != 0) break;

      if (key.size() > lenSection) {
        std::string::size_type nextDotPos = key.find('.', lenSection + 1);

        if (nextDotPos == std::string::npos) {
          datasets.emplace(key.substr(lenSection + 1));
        } else {
          datasets.emplace(key.substr(lenSection + 1, nextDotPos - lenSection - 1));
        }
      }
      ++it;
    }

    return std::vector<std::string>(datasets.begin(), datasets.end());

  }
  
  
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

    std::string freeFormatLine;

    // Parse the free format input
    if (containsData("CQ")) {
      freeFormatLine = getData<std::string>("CQ");
    } else if (containsData("CHRONUSQ")) {
      freeFormatLine = getData<std::string>("CHRONUSQ");
    }

    if (freeFormatLine != "")
      parseFreeCQInput(freeFormatLine);

  }; // CQInputFile::parse

  /**
   *  \brief Overload the << operator for the CQInputFile class
   *
   *  \param [in] os        Output device for data / error output.
   *  \param [in] inputFile CQInputFile object to be printed
   *
   *  \returns std::ostream object
   */
  std::ostream& operator<<(std::ostream& os, const CQInputFile& inputFile) {
    std::vector<std::string> lastSegments;

    for (const auto& [key, value] : inputFile.dict_) {
      std::vector<std::string> currentSegments;
      split(currentSegments, key, ".");

      int commonCount = 0;
      while (commonCount < currentSegments.size()
              and commonCount < lastSegments.size()
              and currentSegments[commonCount] == lastSegments[commonCount]) {
        ++commonCount;
      }

      for (int i = commonCount; i < currentSegments.size() - 1; ++i) {
        os << std::string(i * 4, ' ') << currentSegments[i] << ":" << std::endl;
      }

      os << std::string((currentSegments.size() - 1) * 4, ' ')
          << currentSegments.back() << ": " << value << std::endl;

      lastSegments = std::move(currentSegments);
    }

    return os;
  }

}; // namespace ChronusQ

