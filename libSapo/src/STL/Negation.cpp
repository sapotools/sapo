/**
 * @file Negation.cpp
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @author Alberto Casagrande <acasagrande@units.it>
 * @brief STL Negation
 * @version 0.2
 * @date 2022-05-04
 *
 * @copyright Copyright (c) 2015-2022
 */

#include "STL/Negation.h"

#include "STL/Atom.h"
#include "STL/Conjunction.h"
#include "STL/Disjunction.h"
#include "STL/Negation.h"
#include "STL/Always.h"
#include "STL/Eventually.h"
#include "STL/Until.h"

namespace STL
{

/**
 * @brief A constructor for STL negation formulas
 *
 * This constructor creates an object of the type
 * Negation to represent the STL formula
 * \f$\lnot \textrm{formula}\f$.
 *
 * @param[in] formula is the formula to be negated
 */
Negation::Negation(const std::shared_ptr<STL> formula):
    STL(NEGATION), _subformula(formula)
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
std::ostream &Negation::print(std::ostream &os) const
{
  return os << "! (" << *_subformula << ")";
}

const std::shared_ptr<STL> Negation::get_PNF() const
{
  switch (_subformula->get_type()) {
  case formula_type::ATOM: {
    std::shared_ptr<Atom> atom = std::dynamic_pointer_cast<Atom>(_subformula);
    return std::make_shared<Atom>(-atom->get_expression());
  }
  case formula_type::CONJUNCTION: {
    std::shared_ptr<Conjunction> conj
        = std::dynamic_pointer_cast<Conjunction>(_subformula);
    std::shared_ptr<STL> left
        = std::make_shared<Negation>(conj->get_left_subformula())->get_PNF();
    std::shared_ptr<STL> right
        = std::make_shared<Negation>(conj->get_right_subformula())->get_PNF();
    return std::make_shared<Disjunction>(left, right);
  }
  case formula_type::DISJUNCTION: {
    std::shared_ptr<Disjunction> disj
        = std::dynamic_pointer_cast<Disjunction>(_subformula);
    std::shared_ptr<STL> left
        = std::make_shared<Negation>(disj->get_left_subformula())->get_PNF();
    std::shared_ptr<STL> right
        = std::make_shared<Negation>(disj->get_right_subformula())->get_PNF();
    return std::make_shared<Conjunction>(left, right);
  }
  case formula_type::NEGATION: {
    std::shared_ptr<Negation> neg
        = std::dynamic_pointer_cast<Negation>(_subformula);
    return neg->get_subformula()->get_PNF();
  }
  case formula_type::ALWAYS: {
    std::shared_ptr<Always> alw
        = std::dynamic_pointer_cast<Always>(_subformula);
    std::shared_ptr<STL> sub
        = std::make_shared<Negation>(alw->get_subformula())->get_PNF();
    return std::make_shared<Eventually>(alw->time_bounds().begin(),
                                        alw->time_bounds().end(), sub);
  }
  case formula_type::EVENTUALLY: {
    std::shared_ptr<Eventually> event
        = std::dynamic_pointer_cast<Eventually>(_subformula);
    std::shared_ptr<STL> sub
        = std::make_shared<Negation>(event->get_subformula())->get_PNF();
    return std::make_shared<Always>(event->time_bounds().begin(),
                                    event->time_bounds().end(), sub);
  }
  case formula_type::UNTIL:
    SAPO_ERROR("until does not support negation", std::logic_error);
  default:
    SAPO_ERROR("unsupported formula type", std::logic_error);
  }
}

/**
 * @brief Destroy the STL negation formula
 */
Negation::~Negation() {}

}
