/**
 * @file Eventually.cpp
 * Eventually STL formula
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @version 0.1
 */

#include "Eventually.h"

#include <iostream>

namespace STL
{

/**
 * @brief A constructor for STL eventually formulas
 * 
 * This constructor creates an object of the type 
 * Eventually to represent the STL formula 
 * \f$F_{[\texrm{begin}, \textrm{end}]} \textrm{formula}\f$.
 *
 * @param[in] begin is the begin of temporal interval
 * @param[in] end is the end of temporal interval
 * @param[in] formula is the formula eventually holding
 */
Eventually::Eventually(const int begin, const int end,
                       const std::shared_ptr<STL> formula):
    STL(EVENTUALLY), _subformula(formula), _t_itvl(begin, end)
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
std::ostream &Eventually::print(std::ostream &os) const
{
  return os << "F " << this->_t_itvl << " (" << *_subformula << ")";
}

/**
 * @brief Destroy the STL eventually formula
 */
Eventually::~Eventually() {}

}
