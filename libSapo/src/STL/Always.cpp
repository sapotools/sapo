/**
 * @file Always.cpp
 * Always STL formula
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @version 0.1
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
 * \f$G_{[\texrm{begin}, \textrm{end}]} \textrm{formula}\f$.
 *
 * @param[in] begin is the begin of temporal interval
 * @param[in] end is the end of temporal interval
 * @param[in] formula is the formula globally holding
 */
Always::Always(const int begin, const int end,
               const std::shared_ptr<STL> formula):
    STL(ALWAYS), _subformula(formula), _t_itvl(begin, end)
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
  os << "G " << this->_t_itvl << " (" 
     << *(this->_subformula) << ")";

  return os;
}

/**
 * @brief Destroy the STL always formula
 */
Always::~Always() {}

}