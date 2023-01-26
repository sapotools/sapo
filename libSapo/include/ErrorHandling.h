/**
 * @file ErrorHandling.h
 * @author Alberto Casagrande <acasagrande@units.it>
 * @brief Defines macros for error handling
 * @version 0.1
 * @date 2023-01-26
 *
 * @copyright Copyright (c) 2023
 */

#ifndef ERROR_HANDLING_H_
#define ERROR_HANDLING_H_

#include <sstream>

#define SAPO_ERROR(MSG, TYPE)                                                 \
  {                                                                           \
    std::ostringstream os;                                                    \
    os << __PRETTY_FUNCTION__ << ": " << MSG;                                 \
    throw TYPE(os.str().c_str());                                             \
  }

#endif // ERROR_HANDLING_H_