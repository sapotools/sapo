#ifndef __FORMULA_H__
#define __FORMULA_H__

#include <memory>

#include "AbsSynIO.h"
#include "Expr.h"
#include "STL.h"
#include "Atom.h"
#include "Conjunction.h"
#include "Disjunction.h"
#include "Until.h"
#include "Always.h"
#include "Eventually.h"
#include "SymbolicAlgebra.h"

namespace AbsSyn
{

class Formula
{
  friend std::ostream &operator<<(std::ostream &os, const Formula &f);

public:
  enum formulaType {
    ATOM,  // atomic formula (ex <= 0)
    CONJ,  // conjunction (f1 && f2)
    DISJ,  // disjunction (f1 || f2)
    F_NEG, // negation (! f1)
    ALW,   // always (G[i] f1)
    EVENT, // eventually (F[i] f1)
    UNTIL  // until (f1 U[i] f2)
  };

  Formula(SymbolicAlgebra::Expression<> e)
  {
    type = formulaType::ATOM;
    ex = e;
    f1 = NULL;
    f2 = NULL;
  } // atomic formula
    //	Formula(Formula *l, Formula *r, formulaType t) { type = t; ex = NULL;
    // f1 = l; f2 = r; }		// state formula (CONJ, DISJ, F_NEG
    //	Formula(Formula *l, Formula *r, formulaType t, pair<int,int> in) { type
    //= t; ex = NULL; f1 = l; f2 = r; i = in; }		// path formula (other
    // types)

  Formula *conj(Formula *f);
  Formula *disj(Formula *f);
  Formula *neg();
  Formula *always(std::pair<int, int> in);
  Formula *eventually(std::pair<int, int> in);
  Formula *until(std::pair<int, int> in, Formula *f);

  ~Formula()
  {
    delete (f1);
    delete (f2);
  }

  SymbolicAlgebra::Expression<> getEx()
  {
    return ex;
  }

  Formula *getLeft()
  {
    return f1;
  }

  Formula *getRight()
  {
    return f2;
  }

  std::pair<int, int> getInterval()
  {
    return i;
  }

  formulaType getType() const
  {
    return type;
  }

  bool
  isLinear(const Context &ctx) const; // checks if the formula is a boolean
                                       // combination of linear inequalities

  /*
   * simplification, removes negations
   */
  bool simplify();

  /**
   * @brief Transforms a Formula into a SAPO STL formula.
   *
   * @param m is the input data object.
   * @param vars is the vector of variable symbols.
   * @param params is the vector of parameter symbols.
   * @return a pointer to the STL formula.
   */
  std::shared_ptr<STL>
  toSTL(const Context &ctx, const std::vector<SymbolicAlgebra::Symbol<>> &vars,
        const std::vector<SymbolicAlgebra::Symbol<>> &params) const;

protected:
  Formula()
  {
    ex = 0;
    f1 = NULL;
    f2 = NULL;
  }

  formulaType type;      // type of formula
  SymbolicAlgebra::Expression<> ex;              // expression of atomic formula (ex >= 0)
  Formula *f1;           // left formula
  Formula *f2;           // right formula
  std::pair<int, int> i; // interval of application

  /*
   * utility for semplification
   * returns 0 if no further calls are needed
   * returns 1 if more are required
   * returns 2 if a negation of an until is found
   */
  int simplifyRec();
};

}

#endif
