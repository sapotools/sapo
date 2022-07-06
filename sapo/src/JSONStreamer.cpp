/**
 * @file JSONStreamer.cpp
 * @author Alberto Casagrande <acasagrande@units.it>
 * @brief Functions printing on a JSON stream
 * @version 0.1
 * @date 2022-06-08
 *
 * @copyright Copyright (c) 2022
 */

#include "JSONStreamer.h"

/**
 * @brief Print a linear system in a JSON stream
 *
 * @param out is the output JSON stream
 * @param ls is the linear system to be print
 * @return a reference to the output JSON stream
 */
JSON::ostream &operator<<(JSON::ostream &out, const LinearSystem &ls)
{
  out << "{\"A\":" << ls.A() << ","
      << "\"b\":" << ls.b() << "}";

  return out;
}