/**
 * @file Atom.cpp
 * Atomic STL formula
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @version 0.1
 */

#include "Atom.h"

namespace STL
{

unsigned int Atom::_num_of_atoms = 0;

/**
 * @brief A constructor for STL atomic formulas
 * 
 * This constructor creates an object of the type 
 * Atom to represent the STL formula 
 * \f$\textrm{expression} \geq 0\f$.
 *
 * @param[in] expression is the atom expression
 */
Atom::Atom(const SymbolicAlgebra::Expression<> &expression):
    STL(ATOM), _expr(expression), _id(this->_num_of_atoms++)
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
std::ostream &Atom::print(std::ostream &os) const
{
  return os << this->_expr << " <= 0";
}

/**
 * @brief Destroy the STL atomic formula
 */
Atom::~Atom()
{}

}