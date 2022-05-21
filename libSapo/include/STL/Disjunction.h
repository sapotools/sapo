/**
 * @file Disjunction.h
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @author Alberto Casagrande <acasagrande@units.it>
 * @brief Disjunction STL formula
 * @version 0.2
 * @date 2022-05-04
 * 
 * @copyright Copyright (c) 2015-2022
 */

#ifndef DISJUNCTION_H_
#define DISJUNCTION_H_

#include <memory>

#include "STL.h"

namespace STL
{

/**
 * @brief The STL disjunction formula
 * 
 * This class represents STL disjunction formulas.
 */
class Disjunction : public STL
{

private:
  const std::shared_ptr<STL> _left;   //!< left-side disjunction sub-formula
  const std::shared_ptr<STL> _right;  //!< right-side disjunction sub-formula

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
   * @brief A constructor for STL disjunction formulas
   * 
   * This constructor creates an object of the type 
   * Disjunction to represent the STL formula 
   * \f$\textrm{left} \lor \textrm{right}\f$.
   *
   * @param[in] left is the left-side disjunction sub-formula
   * @param[in] right is the right-side disjunction sub-formula
   */
  Disjunction(const std::shared_ptr<STL> left,
              const std::shared_ptr<STL> right);

  /**
   * @brief Get the left-side disjunction sub-formula
   * 
   * @return the left-side disjunction sub-formula
   */
  inline const std::shared_ptr<STL> get_left_subformula() const
  {
    return _left;
  }

  /**
   * @brief Get the right-side disjunction sub-formula
   * 
   * @return the right-side disjunction sub-formula
   */
  inline const std::shared_ptr<STL> get_right_subformula() const
  {
    return _right;
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
  inline const std::shared_ptr<STL> get_PNF() const
  {
    return std::make_shared<Disjunction>(_left->get_PNF(),
                                         _right->get_PNF());
  }

  /**
   * @brief Get the formula variables
   * 
   * @return the set of formula variables
   */
  std::set<SymbolicAlgebra::Symbol<>> get_variables() const;

  /**
   * @brief Get the formula time bounds 
   * 
   * @return the time interval affecting the formula sematics
   */
  TimeInterval time_bounds() const;

  /**
   * @brief Destroy the STL disjunction formula
   */
  ~Disjunction();
};
}

#endif /* DISJUNCTION_H_ */
