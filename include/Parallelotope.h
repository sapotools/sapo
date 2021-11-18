/**
 * @file Parallelotope.h
 * Represent and manipulate a parallelotope
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @version 0.1
 */

#ifndef PARALLELOTOPE_H_
#define PARALLELOTOPE_H_

#include <vector>

#include "Polytope.h"

struct poly_values { // numerical values for polytopes
  std::vector<double> base_vertex;
  std::vector<double> lenghts;
};

class Parallelotope
{
  using Vector = std::vector<double>;
  using Matrix = std::vector<Vector>;

private:
  unsigned int dim;              // dimension of the parallelotope
  std::vector<GiNaC::lst> vars;  // variables appearing in generato function
                                 // vars[0] q: base vertex
                                 // vars[1] alpha : free variables \in [0,1]
                                 // vars[2] beta : generator amplitudes
  GiNaC::lst generator_function; // generator function
  Matrix u;                      // versors

  Vector hyperplaneThroughPts(const std::vector<Vector> &pts)
      const; // find hyper plane passing through pts
  Vector
  lst2vec(const GiNaC::ex &list) const; // convert a lst to a vector of doubles
  double euclidNorm(const Vector &v) const; // compute euclidean norm

  std::vector<double> base_vertex;
  std::vector<double> lenghts;

public:
  Parallelotope(const std::vector<GiNaC::lst> &vars, const Matrix &u);
  Parallelotope(const std::vector<GiNaC::lst> &vars,
                const Matrix &template_matrix, const Vector &offset);
  Parallelotope(const std::vector<GiNaC::lst> &vars, const Polytope &P);

  /**
   * Get the generator functions
   *
   * @returns generator functions
   */
  const GiNaC::lst &getGeneratorFunction() const
  {
    return this->generator_function;
  }

  /**
   * Get variables of base vertex
   *
   * @returns base vertex variables
   */
  const GiNaC::lst &getQ() const
  {
    return this->vars[0];
  }

  /**
   * Get free variables
   *
   * @returns free variables
   */
  const GiNaC::lst &getAlpha() const
  {
    return this->vars[1];
  }

  /**
   * Get variables of generator lengths
   *
   * @returns generator lengths variables
   */
  const GiNaC::lst &getBeta() const
  {
    return this->vars[2];
  }

  /**
   * Get the parallelotope dimension
   *
   * @returns parallelotope dimension
   */
  const unsigned int &getDim() const
  {
    return this->dim;
  }

  // Representation conversion
  Polytope
  gen2const(const Vector &q,
            const Vector &beta) const;       // from generator to constraints
  poly_values const2gen(Polytope *LS) const; // from constraints to generators
  Polytope getPolytope() const
  {
    return this->gen2const(this->base_vertex, this->lenghts);
  }

  const std::vector<double> &getBaseVertex() const
  {
    return this->base_vertex;
  }

  const std::vector<double> &getLenghts() const
  {
    return this->lenghts;
  }

  Matrix getVersors() const
  {
    return this->u;
  }

  friend void swap(Parallelotope &A, Parallelotope &B);

  virtual ~Parallelotope();
};

#endif /* PARALLELOTOPE_H_ */
