/**
 * @file Until.cpp
 * Until STL formula
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @version 0.1
 */

#include "Until.h"

namespace STL
{

/**
 * @brief A constructor for STL until formulas
 * 
 * This constructor creates an object of the type 
 * Until to represent the STL formula 
 * \f$\textrm{left} U_{[\textrm{begin},\textrm{end}]} \textrm{right}\f$.
 *
 * @param[in] left is the left subformula
 * @param[in] begin is the begin of temporal interval
 * @param[in] end is the end of temporal interval
 * @param[in] right is the right subformula
 */
Until::Until(const std::shared_ptr<STL> left,
             const int begin, const int end,
             const std::shared_ptr<STL> right):
    STL(UNTIL), _left(left), _right(right), _t_itvl(begin, end)
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
std::ostream &Until::print(std::ostream &os) const
{
  return os << "(" << *_left << ") U " 
            << _t_itvl << " (" << *_right << ")";
}

/**
 * @brief Destroy the STL until formula
 */
Until::~Until() {}
}
