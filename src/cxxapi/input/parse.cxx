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
#include <cxxapi/input.hpp>
#include <cxxapi/options.hpp>
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
    std::vector<std::regex> caseSensRegexes;

    // Add case sensitive data keywords here
    caseSensRegexes.emplace_back(".*BASIS/BASIS$");
    caseSensRegexes.emplace_back("^FILES/.*");

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
        std::string dataPath;
        if (sectionHeader == "")
          dataPath = dataHeader;
        else if (sectionHeader == "/")
          dataPath = "/" + dataHeader;
        else
          dataPath = sectionHeader + "/" + dataHeader;
        if (not std::any_of(caseSensRegexes.begin(), caseSensRegexes.end(), [&dataPath](const std::regex& regex) {
              return std::regex_match(dataPath, regex);
            }) and not value.empty())
          strToUpper(value);

        // Create a dictionary entry for the data field in the current
        // section header
        if(not value.empty())
          addData(sectionHeader, dataHeader, value);
        else
          CErr("No data entry for " + dataPath + " in input file.");

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
   *  \brief Splits a query string on the last slash "/" character
   * 
   *  This is a helpder function for the getData function which takes a 
   *  formatted string and splits it into a section and data field.
   *
   *  i.e.  "QM/REFERENCE" -> { "QM", "REFERENCE" }
   *
   *  \param [in] query Query string to be split
   *  \return     std::pair containing the two fields separated by a "/"
   */
  std::pair<std::string, std::string> CQInputFile::splitQuery(const std::string& query) {
    size_t lastSlashPos = query.find_last_of('/');

    // If no slash is found, return the whole string as the first part, and an empty second part
    if (lastSlashPos == std::string::npos) {
      return {"", query};
    }

    // If slash is in the beginning, the section header is "/"
    if (lastSlashPos == 0) {
      return {"/", query.substr(1)};
    }

    // Split the string into two parts: before and after the last slash
    std::string beforeSlash = query.substr(0, lastSlashPos);
    std::string afterSlash = query.substr(lastSlashPos + 1);

    return {beforeSlash, afterSlash};
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
  void CQInputFile::addData(const std::string &prefix, const std::string &key, const std::string &value) {
    if (containsData(prefix + key))
      CErr("Key " + key + " already exists in the parsed input.", std::cout);
    dict_[prefix][key] = value;
  }
  void CQInputFile::addData(const std::string &path, const std::string &value) {
    auto [prefix, key] = splitQuery(path);
    addData(prefix, key, value);
  }
  void CQInputFile::modifyData(const std::string &prefix, const std::string &key, const std::string &value) {
    if (not containsData(prefix + key))
      CErr("Key " + key + " does not exist in the parsed input.", std::cout);
    dict_[prefix][key] = value;
  }
  void CQInputFile::modifyData(const std::string &path, const std::string &value) {
    auto [prefix, key] = splitQuery(path);
    addData(prefix, key, value);
  }


  /**
   * \brief Merge a subsection into the input storage
   * \param [in] subsection Subsection to be merged
   * \param [in] prefix Prefix of the data field
   */
  void CQInputFile::mergeSection(const std::string &prefix,
                                 const std::map<std::string,std::string> &section) {
    if (dict_.find(prefix) == dict_.end())
      dict_[prefix] = section;
    else
      for (const auto &[k, v] : section) {
        addData(prefix, k, v);
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
    return dict_.find(str) != dict_.end();
  }


  /**
   *  \brief Returns a subsection of data fields from the input file
   *
   *  \param [in] section Section heading
   *  \return             Vector of data fields in section
   */
  const std::map<std::string,std::string>& CQInputFile::getSection(const std::string &prefix) const {

    if (not containsSection(prefix))
      CErr("Section " + prefix + " not found in input file!");

    return dict_.at(prefix);
  }


  /**
   *  Checks whether or not the parsed CQ input file contains
   *  a query data field.
   *
   *  \paral  [in] str Query string of a data field (includes section heading)
   *  \return      True if input file contains that data field
   */
  bool CQInputFile::containsData(const std::string &prefix,
                                 const std::string &key) const {
    return containsSection(prefix) and dict_.at(prefix).find(key) != dict_.at(prefix).end();
  }
  bool CQInputFile::containsData(const std::string &path) const {
    const auto [prefix, key] = splitQuery(path);
    return containsData(prefix, key);
  }


  std::set<std::string> CQInputFile::getDataInSection(const std::string &prefix) const  {

    std::set<std::string> keys;

    if (containsSection(prefix))
      for (const auto &[key, value] : dict_.at(prefix))
        keys.insert(key);

    return keys;

  }


  template <typename T>
  T CQStringTo(const std::string &s) {
    std::istringstream iss(s);
    T t;
    iss >> t;
    return t;
  }

  template <>
  bool CQStringTo(const std::string &s) {
    if (s == "TRUE" or s == "ON")
      return true;
    if (s == "FALSE" or s == "OFF")
      return false;
    CErr("Invalid Input For Boolean-Type Keyword!");
    return false;
  }

  template <>
  std::string CQStringTo(const std::string &s) {
    return s;
  }

  template double CQStringTo(const std::string &s);
  template int CQStringTo(const std::string &s);
  template size_t CQStringTo(const std::string &s);

  template <typename T>
  T CQInputFile::getData(const std::string &prefix,
                         const std::string &key) const {

    if (not containsData(prefix, key))
      throw data_not_found(prefix + "/" + key);

    return CQStringTo<T>(dict_.at(prefix).at(key));

  }; // CQInputFile::getData

  template <typename T>
  T CQInputFile::getData(const std::string &path) const {

    const auto &[prefix, key] = splitQuery(path);

    return getData<T>(prefix, key);

  }; // CQInputFile::getData
  /**
   *  \brief Specialization of getData to return std::string of query 
   *  data field
   *
   *  \param [in] query Formatted query string to be parsed
   *  \return     Value of query data field as a std::string
   */
  template
  std::string CQInputFile::getData(const std::string &prefix,
                                   const std::string &key) const;
  template std::string CQInputFile::getData(const std::string &path) const;
  
  /**
   *  \brief Specialization of getData to return int of query 
   *  data field
   *
   *  \param [in] query Formatted query string to be parsed
   *  \return     Value of query data field as a int
   */
  template
  int CQInputFile::getData(const std::string &prefix,
                           const std::string &key) const;
  template int CQInputFile::getData(const std::string &path) const;

  /**
   *  \brief Specialization of getData to return bool of query 
   *  data field
   *
   *  \param [in] query Formatted query string to be parsed
   *  \return     Value of query data field as a bool
   */
  template
  bool CQInputFile::getData(const std::string &prefix,
                            const std::string &key) const;
  template bool CQInputFile::getData(const std::string &path) const;
  
  /**
   *  \brief Specialization of getData to return size_t of query 
   *  data field
   *
   *  \param [in] query Formatted query string to be parsed
   *  \return     Value of query data field as a size_t
   */
  template
  size_t CQInputFile::getData(const std::string &prefix,
                              const std::string &key) const;
  template size_t CQInputFile::getData(const std::string &path) const;
  
  /**
   *  \brief Specialization of getData to return double of query 
   *  data field
   *
   *  \param [in] query Formatted query string to be parsed
   *  \return     Value of query data field as a double
   */
  template
  double CQInputFile::getData(const std::string &prefix,
                              const std::string &key) const;
  template double CQInputFile::getData(const std::string &path) const;


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

    for (const auto& [prefix, section] : inputFile.dict_) {
      os << "[" << prefix << "]" << std::endl;
      for (const auto& [key, value] : section) {
        os << "  " << key << " = " << value << std::endl;
      }
    }

    return os;
  }

}; // namespace ChronusQ

