/**
 * @file Atom.h
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @author Alberto Casagrande <acasagrande@units.it>
 * @brief Atomic STL formula
 * @version 0.2
 * @date 2022-05-04
 * 
 * @copyright Copyright (c) 2015-2022
 */

#ifndef ATOM_H_
#define ATOM_H_

#include "SymbolicAlgebra.h"

#include "STL.h"

namespace STL
{

/**
 * @brief The STL atomic formula
 * 
 * This class represents STL atomic formulas.
 * Every atoms have the form \f$e \geq 0\f$ where 
 * \f$e\f$ is a `SymbolicAlgebra::Expression`.
 */
class Atom : public STL
{
private:
  static unsigned int _num_of_atoms;    //!< Number of instanciated atoms

  SymbolicAlgebra::Expression<> _expr;  //!< Atom expression
  unsigned int _id;                     //!< Atom identifier

protected:

  /**
   * @brief Print the STL formula in a stream
   * 
   * This method is publicly wrapped by the function
   * `operator<<(std::ostream&, const STL::STL &)`.
   * 
   * @param os is the output stream
   * @return a reference to the output stream
   */
  std::ostream &print(std::ostream &os) const;

public:
  /**
   * @brief A constructor for STL atomic formulas
   * 
   * This constructor creates an object of the type 
   * Atom to represent the STL formula 
   * \f$\textrm{expression} \geq 0\f$.
   *
   * @param[in] expression is the atom expression
   */
  Atom(const SymbolicAlgebra::Expression<> &expression);

  /**
   * @brief Get the atom formula expression
   * 
   * Every atomic fomula has the form 
   * \f$\textrm{expression} \geq 0\f$. This method 
   * returns the formula expression.
   * 
   * @return the atom formula expression
   */
  inline const SymbolicAlgebra::Expression<> &get_expression() const
  {
    return _expr;
  }

  /**
   * @brief Get the atom id
   * 
   * Every atom is uniquely identified by a natural number.
   * This method returns the atom identifier.
   * 
   * @return the atom identifier
   */
  inline const unsigned int& get_id() const
  {
    return _id;
  }

  /**
   * @brief Get an equivalent formula in Positive Normal Form (PNF)
   * 
   * A STL formula is in Positive Normal Form (PNF) if it does 
   * not use the negation operator. 
   * For any STL formula \f$\phi\f$ there exists a STL 
   * formula \f$\varphi\f$ in PNF such that 
   * \f[\models \phi \Longleftrightarrow \models \varphi.\f]
   * This method computes the STL formula in PNF of the 
   * current formula. In the case of an atom, this consists 
   * is a copy of the atom itself.
   * 
   * @return The equivalent STL formula in PNF
   */
  inline const std::shared_ptr<STL> get_PNF() const
  {
    return std::make_shared<Atom>(_expr);
  }

  /**
   * @brief Get the formula variables
   * 
   * @return the set of formula variables
   */
  inline std::set<SymbolicAlgebra::Symbol<>> get_variables() const
  {
    return _expr.get_symbols();
  }

  /**
   * @brief Destroy the STL atomic formula
   */
  ~Atom();
};

}

#endif /* ATOM_H_ */
