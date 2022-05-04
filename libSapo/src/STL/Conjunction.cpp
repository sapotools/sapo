/**
 * @file Conjunction.cpp
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @author Alberto Casagrande <acasagrande@units.it>
 * @brief Conjunction STL formula
 * @version 0.2
 * @date 2022-05-04
 * 
 * @copyright Copyright (c) 2015-2022
 */

#include "Conjunction.h"

namespace STL
{

/**
 * @brief A constructor for STL conjunction
 * 
 * This constructor creates an object of the type 
 * Conjunction to represent the STL formula 
 * \f$\textrm{left} \land \textrm{right}\f$.
 *
 * @param[in] left is the left-side conjunction subformula
 * @param[in] right is the right-side conjunction subformula
 */
Conjunction::Conjunction(const std::shared_ptr<STL> left,
                         const std::shared_ptr<STL> right):
    STL(CONJUNCTION), _left(left), _right(right)
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
std::ostream &Conjunction::print(std::ostream &os) const
{
  return os << "(" << *_left << ") && (" << *_right << ")";
}

/**
 * @brief Get the formula time bounds 
 * 
 * @return the time interval affecting the formula sematics
 */
TimeInterval Conjunction::time_bounds() const
{
  const TimeInterval ti1 = _left->time_bounds();
  const TimeInterval ti2 = _right->time_bounds();

  return TimeInterval(std::min(ti1.begin(), ti2.begin()),
                      std::max(ti1.end(), ti2.end()));
}

/**
 * @brief Destroy a STL conjunction formula
 */
Conjunction::~Conjunction() {}

}