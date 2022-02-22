/**
 * @file Negation.cpp
 * Negation STL formula
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @version 0.1
 */

#include "Negation.h"

/**
 * Constructor that instantiates a negation formula (not f)
 *
 * @param[in] f subformula
 */
Negation::Negation(const std::shared_ptr<STL> f):
    STL(NEGATION),
    f(f)
{
}

/**
 * Print the formula
 */
std::ostream &Negation::print(std::ostream &os) const
{
  return os << "! (" << *f << ")";
}

Negation::~Negation() {}


const std::shared_ptr<STL> Negation::simplify() const
{
	switch (f->getType()) {
		case formula_type::ATOM: {
			std::shared_ptr<Atom> atom = std::dynamic_pointer_cast<Atom>(f);
			return std::make_shared<Atom>(-atom->getPredicate());
		}
		case formula_type::CONJUNCTION: {
			std::shared_ptr<Conjunction> conj = std::dynamic_pointer_cast<Conjunction>(f);
			std::shared_ptr<STL> left = std::make_shared<Negation>(conj->getLeftSubFormula())->simplify();
			std::shared_ptr<STL> right = std::make_shared<Negation>(conj->getRightSubFormula())->simplify();
			return std::make_shared<Disjunction>(left, right);
		}
		case formula_type::DISJUNCTION: {
			std::shared_ptr<Disjunction> disj = std::dynamic_pointer_cast<Disjunction>(f);
			std::shared_ptr<STL> left = std::make_shared<Negation>(disj->getLeftSubFormula())->simplify();
			std::shared_ptr<STL> right = std::make_shared<Negation>(disj->getRightSubFormula())->simplify();
			return std::make_shared<Conjunction>(left, right);
		}
		case formula_type::NEGATION: {
			std::shared_ptr<Negation> neg = std::dynamic_pointer_cast<Negation>(f);
			return neg->getSubFormula()->simplify();
		}
		case formula_type::ALWAYS: {
			std::shared_ptr<Always> alw = std::dynamic_pointer_cast<Always>(f);
			std::shared_ptr<STL> sub = std::make_shared<Negation>(alw->getSubFormula())->simplify();
			return std::make_shared<Eventually>(alw->time_bounds().begin(), alw->time_bounds().end(), sub);
		}
		case formula_type::EVENTUALLY: {
			std::shared_ptr<Eventually> event = std::dynamic_pointer_cast<Eventually>(f);
			std::shared_ptr<STL> sub = std::make_shared<Negation>(event->getSubFormula())->simplify();
			return std::make_shared<Always>(event->time_bounds().begin(), event->time_bounds().end(), sub);
		}
		case formula_type::UNTIL: {
			throw std::logic_error("Negation of until is illegal");
		}
		default: {
			throw std::logic_error("Unsupported formula type");
		}
	}
}
