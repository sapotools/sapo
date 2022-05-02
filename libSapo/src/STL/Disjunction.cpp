/**
 * @file Disjunction.cpp
 * Disjunction STL formula
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @version 0.1
 */

#include "Disjunction.h"

#include <iostream>

namespace STL
{

/**
 * @brief A constructor for STL dicjunction formulas
 * 
 * This constructor creates an object of the type 
 * Disjunction to represent the STL formula 
 * \f$\textrm{left} \lor \textrm{right}\f$.
 *
 * @param[in] left is the left-side disjunction subformula
 * @param[in] right is the right-side disjunction subformula
 */
Disjunction::Disjunction(const std::shared_ptr<STL> left,
                         const std::shared_ptr<STL> right):
    STL(DISJUNCTION), _left(left), _right(right)
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
std::ostream &Disjunction::print(std::ostream &os) const
{
  return os << "(" << *_left << ") || (" << *_right << ")";
}

/**
 * @brief Get the formula time bounds 
 * 
 * @return the time interval affecting the formula sematics
 */
TimeInterval Disjunction::time_bounds() const
{
  const TimeInterval ti1 = _left->time_bounds();
  const TimeInterval ti2 = _right->time_bounds();

  return TimeInterval(std::min(ti1.begin(), ti2.begin()),
                      std::max(ti1.end(), ti2.end()));
}

/**
 * @brief Destroy the STL disjunction formula
 */
Disjunction::~Disjunction() {}

}
