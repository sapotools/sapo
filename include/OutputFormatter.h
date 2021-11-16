#ifndef OUTPUTFORMATTER_H_
#define OUTPUTFORMATTER_H_

#include <algorithm>
#include <string>

template<typename T>
class OutputFormatter
{
public:
  OutputFormatter() {}

  static std::string field_begin(std::string name)
  {
    (void)name;
    return "";
  }

  static std::string field_end()
  {
    return "";
  }

  static std::string field_separator()
  {
    return "";
  }

  static std::string sequence_begin()
  {
    return "";
  }

  static std::string sequence_end()
  {
    return "";
  }

  static std::string sequence_separator()
  {
    return "";
  }

  static std::string list_begin()
  {
    return "";
  }

  static std::string list_end()
  {
    return "";
  }

  static std::string list_separator()
  {
    return "";
  }

  static std::string short_list_begin()
  {
    return T::list_begin();
  }

  static std::string short_list_end()
  {
    return T::list_end();
  }

  static std::string short_list_separator()
  {
    return T::list_separator();
  }

  static std::string empty_list()
  {
    return "";
  }

  static std::string set_begin()
  {
    return "";
  }

  static std::string set_end()
  {
    return "";
  }

  static std::string set_separator()
  {
    return "";
  }

  static std::string empty_set()
  {
    return "";
  }

  static std::string object_header()
  {
    return "";
  }

  static std::string object_footer()
  {
    return "";
  }
};

template<>
class OutputFormatter<std::ostream>
{
public:
  static std::string field_header(std::string name)
  {
    std::transform(name.begin(), name.end(), name.begin(), ::toupper);
    return name + "\n";
  }

  static std::string field_footer()
  {
    return "";
  }

  static std::string field_separator()
  {
    return "\n\n";
  }

  static std::string sequence_begin()
  {
    return "";
  }

  static std::string sequence_end()
  {
    return "";
  }

  static std::string sequence_separator()
  {
    return "\n" + std::string(20, '-') + "\n";
  }

  static std::string list_begin()
  {
    return "";
  }

  static std::string list_end()
  {
    return "\n";
  }

  static std::string list_separator()
  {
    return "\n" + std::string(20, '=') + "\n";
  }

  static std::string short_list_begin()
  {
    return "[";
  }

  static std::string short_list_end()
  {
    return "]";
  }

  static std::string short_list_separator()
  {
    return ",";
  }

  static std::string empty_list()
  {
    std::string pre_post(10, '=');
    return pre_post + "empty list" + pre_post;
  }

  static std::string set_begin()
  {
    return "";
  }

  static std::string set_end()
  {
    return "";
  }

  static std::string set_separator()
  {
    return "\n\n";
  }

  static std::string empty_set()
  {
    return "empty set";
  }

  static std::string object_header()
  {
    return "";
  }

  static std::string object_footer()
  {
    return "";
  }
};

namespace JSON
{
class ostream;
}

template<>
class OutputFormatter<JSON::ostream>
{
public:
  OutputFormatter() {}

  static std::string field_header(std::string name)
  {
    std::transform(name.begin(), name.end(), name.begin(), ::tolower);
    return "\"" + name + "\":";
  }

  static std::string field_footer()
  {
    return "";
  }

  static std::string field_separator()
  {
    return ",";
  }

  static std::string sequence_begin()
  {
    return OutputFormatter<JSON::ostream>::list_begin();
  }

  static std::string sequence_end()
  {
    return OutputFormatter<JSON::ostream>::list_end();
  }

  static std::string sequence_separator()
  {
    return OutputFormatter<JSON::ostream>::list_separator();;
  }

  static std::string list_begin()
  {
    return "[";
  }

  static std::string list_end()
  {
    return "]";
  }

  static std::string list_separator()
  {
    return ",";
  }

  static std::string short_list_begin()
  {
    return OutputFormatter<JSON::ostream>::list_begin();
  }

  static std::string short_list_end()
  {
    return OutputFormatter<JSON::ostream>::list_end();
  }

  static std::string short_list_separator()
  {
    return OutputFormatter<JSON::ostream>::list_separator();
  }

  static std::string empty_list()
  {
    return (OutputFormatter<JSON::ostream>::list_begin()+
            OutputFormatter<JSON::ostream>::list_end());
  }

  static std::string set_begin()
  {
    return OutputFormatter<JSON::ostream>::list_begin();
  }

  static std::string set_end()
  {
    return OutputFormatter<JSON::ostream>::list_end();
  }

  static std::string set_separator()
  {
    return OutputFormatter<JSON::ostream>::list_separator();;
  }

  static std::string empty_set()
  {
    return OutputFormatter<JSON::ostream>::empty_set();
  }

  static std::string object_header()
  {
    return "{";
  }

  static std::string object_footer()
  {
    return "}";
  }
};

#endif // OUTPUTFORMATTER_H_
