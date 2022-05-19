/**
 * @file Eventually.h
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @author Alberto Casagrande <acasagrande@units.it>
 * @brief Eventually STL formula
 * @version 0.2
 * @date 2022-05-04
 * 
 * @copyright Copyright (c) 2015-2022
 */

#ifndef EVENTUALLY_H_
#define EVENTUALLY_H_

#include <memory>

#include "STL.h"
#include "TimeInterval.h"

namespace STL
{

/**
 * @brief The STL eventually formula
 * 
 * This class represents STL eventually formulas.
 */
class Eventually : public STL
{
private:
  const std::shared_ptr<STL> _subformula; //!< eventually subformula
  TimeInterval _t_itvl;                   //!< temporal formula bounds

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
  Eventually(const int begin, const int end, 
             const std::shared_ptr<STL> formula);

  /**
   * @brief Get the eventually subformula
   * 
   * @return the eventually subformula
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
  inline const std::shared_ptr<STL> get_PNF() const
  {
    return std::make_shared<Eventually>(_t_itvl.begin(), _t_itvl.end(),
                                        _subformula->get_PNF());
  }

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
    return _t_itvl;
  }

  /**
   * @brief Destroy the STL eventually formula
   */
  ~Eventually();
};

}

#endif /* Eventually_H */
