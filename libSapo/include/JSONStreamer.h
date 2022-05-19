#ifndef JSONSTREAMER_H_
#define JSONSTREAMER_H_

#include <ostream>
#include <string>
#include <algorithm>
#include <list>

namespace JSON
{

enum Command { debug, production };

/**
 * @brief A JSON output stream
 */
class ostream : public std::ostream
{
  bool filter_unnecessary; //!< A flag to filter spaces and new lines

public:
  /**
   * @brief A constructor for `JSON::ostream`
   * 
   * @param os is the output stream to which the JSON output will be redirected
   * @param filter_unnecessary is a flag to filter spaces and new lines
   */
  ostream(std::ostream &os = std::cout, const bool filter_unnecessary = false):
      std::ostream(os.rdbuf()), filter_unnecessary(filter_unnecessary)
  {
  }

  /**
   * @brief Print a JSON command in a JSON stream
   * 
   * @param out is the JSON output stream
   * @param cmd is a JSON command
   * @return a reference to the JSON output stream
   */
  friend JSON::ostream &operator<<(JSON::ostream &out,
                                   const JSON::Command &cmd)
  {
    switch (cmd) {
    case debug:
      out.filter_unnecessary = false;
      break;
    case production:
      out.filter_unnecessary = true;
      break;
    }

    return out;
  }

  /**
   * @brief Print a character in a JSON stream
   * 
   * @param out is the JSON output stream
   * @param ch is a character
   * @return a reference to the JSON output stream
   */
  friend JSON::ostream &operator<<(JSON::ostream &out, const char &ch)
  {
    std::list<char> un_chars{'\n', '\t', ' '};

    if (!out.filter_unnecessary
        || std::find(un_chars.begin(), un_chars.end(), ch) == un_chars.end()) {
      dynamic_cast<std::ostream &>(out) << ch;
    }

    return out;
  }

  /**
   * @brief Print a string in a JSON stream
   * 
   * @param out is the JSON output stream
   * @param str is a string
   * @return a reference to the JSON output stream
   */
  friend JSON::ostream &operator<<(JSON::ostream &out, const std::string &str)
  {
    if (out.filter_unnecessary) {
      std::string new_str(str);
      std::list<char> un_chars{'\n', '\t', ' '};

      for (auto it = std::begin(un_chars); it != std::end(un_chars); ++it) {
        new_str.erase(std::remove(new_str.begin(), new_str.end(), *it),
                      new_str.end());
      }

      dynamic_cast<std::ostream &>(out) << new_str;
    } else {
      dynamic_cast<std::ostream &>(out) << str;
    }

    return out;
  }

  /**
   * @brief Print a value in a JSON stream
   * 
   * @param out is the JSON output stream
   * @param value is a value
   * @return a reference to the JSON output stream
   */
  friend JSON::ostream &operator<<(JSON::ostream &out, const double &value)
  {
    if (value == 0) {
      // This is to avoid "-0"
      dynamic_cast<std::ostream &>(out) << 0;
    } else {
      dynamic_cast<std::ostream &>(out) << value;
    }

    return out;
  }

  /**
   * @brief Print a constant string in a JSON stream
   * 
   * @param out is the JSON output stream
   * @param value is a constant string
   * @return a reference to the JSON output stream
   */
  friend JSON::ostream &operator<<(JSON::ostream &out, const char *value)
  {
    dynamic_cast<std::ostream &>(out) << value;

    return out;
  }
};

}

/**
 * @brief Print a linear system in a JSON stream
 * 
 * @param out is the output JSON stream
 * @param ls is the linear system to be print
 * @return a reference to the output JSON stream
 */
template<typename T>
JSON::ostream &operator<<(JSON::ostream &out, const std::vector<T> &v)
{
  out << "[";
  for (auto el_it = std::begin(v); el_it != std::end(v); ++el_it) {
    if (el_it != std::begin(v)) {
      out << ",";
    }
    out << *el_it;
  }
  out << "]";

  return out;
}

#endif // JSONSTREAMER_H_