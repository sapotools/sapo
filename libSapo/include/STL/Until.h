/**
 * @file Until.h
 * Until STL formula
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @author Alberto Casagrande <acasagrande@units.it>
 * @version 0.2
 */

#ifndef UNTIL_H_
#define UNTIL_H_

#include <iostream>
#include <memory>

#include "STL.h"
#include "TimeInterval.h"

namespace STL {

/**
 * @brief The STL until formula
 * 
 * This class represents STL until formulas.
 */
class Until : public STL
{

private:
  const std::shared_ptr<STL> _left;   //!< left-side until subformula
  const std::shared_ptr<STL> _right;  //!< right-side until subformula
  TimeInterval _t_itvl;               //!< temporal formula bounds

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
   * @brief A constructor for STL until formulas
   * 
   * This constructor creates an object of the type 
   * Until to represent the STL formula 
   * \f$\textrm{left} U_{[\textrm{begin},\textrm{end}]} \textrm{right}\f$.
   *
   * @param[in] left  is the left subformula
   * @param[in] begin is the begin of temporal interval
   * @param[in] end is the end of temporal interval
   * @param[in] right  is the right subformula
   */
  Until(const std::shared_ptr<STL> left, 
        const int begin, const int end,
        const std::shared_ptr<STL> right);

  /**
   * @brief Get the left-side until subformula
   * 
   * @return the left-side until subformula
   */
  inline const std::shared_ptr<STL> get_left_subformula() const
  {
    return _left;
  }

  /**
   * @brief Get the right-side until subformula
   * 
   * @return the right-side until subformula
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
    return std::make_shared<Until>(_left->get_PNF(), _t_itvl.begin(),
                                   _t_itvl.end(), _right->get_PNF());
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
   * @brief Destroy the STL until formula
   */
  ~Until();
};

}
#endif /* UNTIL_H */
