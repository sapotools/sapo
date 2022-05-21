/**
 * @file Negation.h
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @author Alberto Casagrande <acasagrande@units.it>
 * @brief STL Negation
 * @version 0.2
 * @date 2022-05-04
 * 
 * @copyright Copyright (c) 2015-2022
 */

#ifndef __NEGATION_H__
#define __NEGATION_H__

#include <iostream>
#include <memory>

#include "STL.h"

#include "TimeInterval.h"

namespace STL 
{

/**
 * @brief The STL negation formula
 * 
 * This class represents STL negation formulas.
 */
class Negation : public STL
{

private:
  const std::shared_ptr<STL> _subformula; //!< negation sub-formula

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
   * @brief A constructor for STL negation formulas
   * 
   * This constructor creates an object of the type 
   * Negation to represent the STL formula 
   * \f$\lnot \textrm{formula}\f$.
   *
   * @param[in] formula is the formula to be negated
   */
  Negation(const std::shared_ptr<STL> formula);

  /**
   * @brief Get the negation subformula
   * 
   * @return the negation subformula
   */
  inline const std::shared_ptr<STL> get_subformula() const
  {
    return _subformula;
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
   * current formula.
   * 
   * @return The equivalent STL formula in PNF
   */
  const std::shared_ptr<STL> get_PNF() const;

  /**
   * @brief Get the formula variables
   * 
   * @return the set of formula variables
   */
  inline std::set<SymbolicAlgebra::Symbol<>> get_variables() const
  {
    return _subformula->get_variables();
  }

  /**
   * @brief Get the formula time bounds 
   * 
   * @return the time interval affecting the formula sematics
   */
  inline TimeInterval time_bounds() const
  {
    return _subformula->time_bounds();
  }

  /**
   * @brief Destroy the STL negation formula
   */
  ~Negation();
};

}
#endif /* Eventually_H */
