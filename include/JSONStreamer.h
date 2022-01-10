#ifndef JSONSTREAMER_H_
#define JSONSTREAMER_H_

#include <ostream>
#include <string>
#include <algorithm>

namespace JSON
{

enum Command { debug, production };

class ostream : public std::ostream
{
  bool filter_unnecessary;

public:
  ostream(std::ostream &os = std::cout, const bool filter_unnecessary = false):
      std::ostream(os.rdbuf()), filter_unnecessary(filter_unnecessary)
  {
  }

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

  friend JSON::ostream &operator<<(JSON::ostream &out, const char &ch)
  {
    std::list<char> un_chars{'\n', '\t', ' '};

    if (!out.filter_unnecessary
        || std::find(un_chars.begin(), un_chars.end(), ch) == un_chars.end()) {
      dynamic_cast<std::ostream &>(out) << ch;
    }

    return out;
  }

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

  friend JSON::ostream &operator<<(JSON::ostream &out, const char *value)
  {
    dynamic_cast<std::ostream &>(out) << value;

    return out;
  }
};

}
#endif // JSONSTREAMER_H_