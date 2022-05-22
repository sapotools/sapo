#ifndef OUTPUTFORMATTER_H_
#define OUTPUTFORMATTER_H_

#include <algorithm>
#include <string>

/**
 * @brief A formatter class
 *
 * @tparam T the output stream type
 */
template<typename T>
class OutputFormatter
{
public:
  /**
   * @brief The constructor
   */
  OutputFormatter() {}

  /**
   * @brief Get a field header
   *
   * This method formats the header of a field.
   *
   * @param name is the name of the field to be printed
   * @return the header of the field
   */
  static std::string field_begin(std::string name)
  {
    (void)name;
    return "";
  }

  /**
   * @brief Get a field footer
   *
   * This method formats the footer of a field.
   *
   * @return the footer of the field
   */
  static std::string field_end()
  {
    return "";
  }

  /**
   * @brief Get a field separator
   *
   * @return a field separator
   */
  static std::string field_separator()
  {
    return "";
  }

  /**
   * @brief Get the header of a sequence
   *
   * This method formats the header of a sequence.
   *
   * @return the header of a sequence
   */
  static std::string sequence_begin()
  {
    return "";
  }

  /**
   * @brief Get the footer of a sequence
   *
   * This method formats the footer of a sequence.
   *
   * @return the footer of a sequence
   */
  static std::string sequence_end()
  {
    return "";
  }

  /**
   * @brief Get a sequence separator
   *
   * @return a sequence separator
   */
  static std::string sequence_separator()
  {
    return "";
  }

  /**
   * @brief Get the header of a list
   *
   * This method formats the header of a list.
   *
   * @return the header of a list
   */
  static std::string list_begin()
  {
    return "";
  }

  /**
   * @brief Get the footer of a list
   *
   * This method formats the footer of a list.
   *
   * @return the footer of a list
   */
  static std::string list_end()
  {
    return "";
  }

  /**
   * @brief Get a list separator
   *
   * @return a list separator
   */
  static std::string list_separator()
  {
    return "";
  }

  /**
   * @brief Get the header of a short list
   *
   * This method formats the header of a short list.
   *
   * @return the header of a short list
   */
  static std::string short_list_begin()
  {
    return T::list_begin();
  }

  /**
   * @brief Get the footer of a short list
   *
   * This method formats the footer of a short list.
   *
   * @return the footer of a short list
   */
  static std::string short_list_end()
  {
    return T::list_end();
  }

  /**
   * @brief Get the separator of a short list
   *
   * @return the separator of a short list
   */
  static std::string short_list_separator()
  {
    return T::list_separator();
  }

  /**
   * @brief Get the representation of an empty list
   *
   * This method formats an empty list.
   *
   * @return the representation of an empty list
   */
  static std::string empty_list()
  {
    return "";
  }

  /**
   * @brief Get the header of a set
   *
   * This method formats the header of a set.
   *
   * @return the header of a set
   */
  static std::string set_begin()
  {
    return "";
  }

  /**
   * @brief Get the footer of a set
   *
   * This method formats the footer of a set.
   *
   * @return the footer of a set
   */
  static std::string set_end()
  {
    return "";
  }

  /**
   * @brief Get the separator for a set
   *
   * @return the separator for a set
   */
  static std::string set_separator()
  {
    return "";
  }

  /**
   * @brief Get the representation of an empty set
   *
   * This method formats an empty set.
   *
   * @return the representation of an empty set
   */
  static std::string empty_set()
  {
    return "";
  }

  /**
   * @brief Get the header of an object
   *
   * This method formats the header of an object.
   *
   * @return the header of an object
   */
  static std::string object_header()
  {
    return "";
  }

  /**
   * @brief Get the footer of an object
   *
   * This method formats the footer of an object.
   *
   * @return the footer of an object
   */
  static std::string object_footer()
  {
    return "";
  }
};

/**
 * @brief The formatter class specialization for standard output streams
 */
template<>
class OutputFormatter<std::ostream>
{
public:
  /**
   * @brief Get a field header
   *
   * This method formats the header of a field.
   *
   * @param name is the name of the field to be printed
   * @return the header of the field
   */
  static std::string field_begin(std::string name)
  {
    std::transform(name.begin(), name.end(), name.begin(), ::toupper);
    return name + "\n";
  }

  /**
   * @brief Get a field footer
   *
   * This method formats the footer of a field.
   *
   * @return the footer of the field
   */
  static std::string field_end()
  {
    return "";
  }

  /**
   * @brief Get a field separator
   *
   * @return a separator of the field
   */
  static std::string field_separator()
  {
    return "\n\n";
  }

  /**
   * @brief Get a sequence header
   *
   * This method formats the header of a sequence.
   *
   * @return the header of a sequence
   */
  static std::string sequence_begin()
  {
    return "";
  }

  /**
   * @brief Get a sequence footer
   *
   * This method formats the footer of a sequence.
   *
   * @return the footer of a sequence
   */
  static std::string sequence_end()
  {
    return "";
  }

  /**
   * @brief Get a sequence separator
   *
   * @return a separator of a sequence
   */
  static std::string sequence_separator()
  {
    return "\n" + std::string(20, '-') + "\n";
  }

  /**
   * @brief Get a list header
   *
   * This method formats the header of a list.
   *
   * @return the header of a list
   */
  static std::string list_begin()
  {
    return "";
  }

  /**
   * @brief Get a list footer
   *
   * This method formats the footer of a list.
   *
   * @return the footer of a list
   */
  static std::string list_end()
  {
    return "\n";
  }

  /**
   * @brief Get a list separator
   *
   * @return a separator of a list
   */
  static std::string list_separator()
  {
    return "\n" + std::string(20, '=') + "\n";
  }

  /**
   * @brief Get a short list header
   *
   * This method formats the header of a short list.
   *
   * @return the header of a short list
   */
  static std::string short_list_begin()
  {
    return "[";
  }

  /**
   * @brief Get a short list footer
   *
   * This method formats the footer of a short list.
   *
   * @return the footer of a short list
   */
  static std::string short_list_end()
  {
    return "]";
  }

  /**
   * @brief Get a short list separator
   *
   * @return a separator of a short list
   */
  static std::string short_list_separator()
  {
    return ",";
  }

  /**
   * @brief Get the representation of an empty list
   *
   * This method formats an empty list.
   *
   * @return the representation of an empty list
   */
  static std::string empty_list()
  {
    std::string pre_post(10, '=');
    return pre_post + "empty list" + pre_post;
  }

  /**
   * @brief Get a set header
   *
   * This method formats the header of a set.
   *
   * @return the header of a set
   */
  static std::string set_begin()
  {
    return "";
  }

  /**
   * @brief Get a set footer
   *
   * This method formats the footer of a set.
   *
   * @return the footer of a set
   */
  static std::string set_end()
  {
    return "";
  }

  /**
   * @brief Get the separator for a set
   *
   * @return the separator for a set
   */
  static std::string set_separator()
  {
    return "\n\n";
  }

  /**
   * @brief Get the representation of an empty set
   *
   * This method formats an empty set.
   *
   * @return the representation of an empty set
   */
  static std::string empty_set()
  {
    return "empty set";
  }

  /**
   * @brief Get an object header
   *
   * This method formats the header of an object.
   *
   * @return the header of an object
   */
  static std::string object_header()
  {
    return "";
  }

  /**
   * @brief Get an object footer
   *
   * This method formats the footer of an object.
   *
   * @return the footer of an object
   */
  static std::string object_footer()
  {
    return "";
  }
};

namespace JSON
{
class ostream;
}

/**
 * @brief The formatter class specialization for JSON output streams
 */
template<>
class OutputFormatter<JSON::ostream>
{
public:
  /**
   * @brief Constructor
   */
  OutputFormatter() {}

  /**
   * @brief Get a field header
   *
   * This method formats the header of a field.
   *
   * @param name is the name of the field to be printed
   * @return the header of the field
   */
  static std::string field_begin(std::string name)
  {
    std::transform(name.begin(), name.end(), name.begin(), ::tolower);
    return "\"" + name + "\":";
  }

  /**
   * @brief Get a field footer
   *
   * This method formats the footer of a field.
   *
   * @return the footer of the field
   */
  static std::string field_end()
  {
    return "";
  }

  /**
   * @brief Get a field separator
   *
   * @return a field separator
   */
  static std::string field_separator()
  {
    return ",";
  }

  /**
   * @brief Get the header of a sequence
   *
   * This method formats the header of a sequence.
   *
   * @return the header of a sequence
   */
  static std::string sequence_begin()
  {
    return OutputFormatter<JSON::ostream>::list_begin();
  }

  /**
   * @brief Get the footer of a sequence
   *
   * This method formats the footer of a sequence.
   *
   * @return the footer of a sequence
   */
  static std::string sequence_end()
  {
    return OutputFormatter<JSON::ostream>::list_end();
  }

  /**
   * @brief Get a sequence separator
   *
   * @return a sequence separator
   */
  static std::string sequence_separator()
  {
    return OutputFormatter<JSON::ostream>::list_separator();
    ;
  }

  /**
   * @brief Get the header of a list
   *
   * This method formats the header of a list.
   *
   * @return the header of a list
   */
  static std::string list_begin()
  {
    return "[";
  }

  /**
   * @brief Get the separator of a short list
   *
   * @return the separator of a short list
   */
  static std::string list_end()
  {
    return "]";
  }

  /**
   * @brief Get the separator of a short list
   *
   * @return the separator of a short list
   */
  static std::string list_separator()
  {
    return ",";
  }

  /**
   * @brief Get the header of a short list
   *
   * This method formats the header of a short list.
   *
   * @return the header of a short list
   */
  static std::string short_list_begin()
  {
    return OutputFormatter<JSON::ostream>::list_begin();
  }

  /**
   * @brief Get the footer of a short list
   *
   * This method formats the footer of a short list.
   *
   * @return the footer of a short list
   */
  static std::string short_list_end()
  {
    return OutputFormatter<JSON::ostream>::list_end();
  }

  /**
   * @brief Get the separator of a short list
   *
   * @return the separator of a short list
   */
  static std::string short_list_separator()
  {
    return OutputFormatter<JSON::ostream>::list_separator();
  }

  /**
   * @brief Get the representation of an empty list
   *
   * This method formats an empty list.
   *
   * @return the representation of
   */
  static std::string empty_list()
  {
    return (OutputFormatter<JSON::ostream>::list_begin()
            + OutputFormatter<JSON::ostream>::list_end());
  }

  /**
   * @brief Get the header of a set
   *
   * This method formats the header of a set.
   *
   * @return the header of a set
   */
  static std::string set_begin()
  {
    return OutputFormatter<JSON::ostream>::list_begin();
  }

  /**
   * @brief Get the representation of an empty set
   *
   * This method formats an empty set.
   *
   * @return the representation of an empty set
   */
  static std::string set_end()
  {
    return OutputFormatter<JSON::ostream>::list_end();
  }

  /**
   * @brief Get the separator for a set
   *
   * @return the separator for a set
   */
  static std::string set_separator()
  {
    return OutputFormatter<JSON::ostream>::list_separator();
  }

  /**
   * @brief Get the representation of an empty set
   *
   * This method formats an empty set.
   *
   * @return the representation of an empty set
   */
  static std::string empty_set()
  {
    return OutputFormatter<JSON::ostream>::set_begin()
           + OutputFormatter<JSON::ostream>::set_end();
  }

  /**
   * @brief Get an object header
   *
   * This method formats the header of an object.
   *
   * @return the header of an object
   */
  static std::string object_header()
  {
    return "{";
  }

  /**
   * @brief Get an object footer
   *
   * This method formats the footer of an object.
   *
   * @return the footer of an object
   */
  static std::string object_footer()
  {
    return "}";
  }
};

#endif // OUTPUTFORMATTER_H_
