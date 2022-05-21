/**
 * @file Always.cpp
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @author Alberto Casagrande <acasagrande@units.it>
 * @brief Always STL formula
 * @version 0.2
 * @date 2022-05-04
 * 
 * @copyright Copyright (c) 2015-2022
 */

#include <iostream>

#include "Always.h"

namespace STL
{

/**
 * @brief A constructor for STL always formulas
 * 
 * This constructor creates an object of the type 
 * Always to represent the STL formula 
 * \f$G_{[\textrm{begin}, \textrm{end}]} \textrm{formula}\f$.
 *
 * @param[in] begin is the begin of temporal interval
 * @param[in] end is the end of temporal interval
 * @param[in] formula is the formula globally holding
 */
Always::Always(const int begin, const int end,
               const std::shared_ptr<STL> formula):
    STL(ALWAYS), _subformula(formula), _t_interval(begin, end)
{
}

/**
 * @brief Print the STL formula in a stream
 * 
 * This method is publicly wrapped by the function
 * `operator<<(std::ostream&, const STL::STL &)`.
 * 
 * @param os is the output stream
 * @return a reference to the output stream
 */
std::ostream &Always::print(std::ostream &os) const
{
  os << "G " << this->_t_interval << " (" 
     << *(this->_subformula) << ")";

  return os;
}

/**
 * @brief Destroy the STL always formula
 */
Always::~Always() {}

}