/**
 * @file Until.cpp
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @author Alberto Casagrande <acasagrande@units.it>
 * @brief Until STL formula
 * @version 0.2
 * @date 2022-05-04
 * 
 * @copyright Copyright (c) 2015-2022
 */

#include "Until.h"

#include <algorithm>
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
 * @brief Get the formula variables
 * 
 * @return the set of formula variables
 */
std::set<SymbolicAlgebra::Symbol<>> Until::get_variables() const
{
  auto variables = _left->get_variables();
  auto right_vars = _right->get_variables();

  variables.insert(std::begin(right_vars), std::end(right_vars));

  return variables;
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
