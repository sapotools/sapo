/**
 * @file Polytope.h
 * Represent and manipulate polytopes (reached states, parameters, etc.)
 *
 * @author Alberto Casagrande <acasagrande@units.it>
 * @version 0.1
 */

#ifndef POLYTOPE_H_
#define POLYTOPE_H_

#include "LinearSystem.h"

class Polytope : public LinearSystem
{

private:
  std::list<Polytope> split(const std::vector<bool> &bvect_base,
                            const unsigned int cidx,
                            std::list<Polytope> &tmp_covering,
                            std::vector<std::vector<double>> &A,
                            std::vector<double> &b) const;

public:
  /**
   * Constructor that instantiates a unbounded polytope
   */
  Polytope(): LinearSystem() {}

  /**
   * Copy constructor
   *
   * @param[in] orig the original polytope
   */
  Polytope(const Polytope &orig): LinearSystem(orig) {}

  /**
   * Swap constructor
   *
   * @param[in] orig the original polytope
   */
  Polytope(Polytope &&orig): LinearSystem(orig) {}

  /**
   * Constructor
   *
   * @param[in] A template matrix
   * @param[in] b offset vector
   */
  Polytope(const std::vector<std::vector<double>> &A,
           const std::vector<double> &b):
      LinearSystem(A, b)
  {
  }

  /**
   * Constructor from a set of symbolic expressions
   *
   * @param[in] vars list of variables appearing in the constraints
   * @param[in] constraints symbolic constraints
   */
  Polytope(const GiNaC::lst &vars, const GiNaC::lst &constraints):
      LinearSystem(vars, constraints)
  {
  }

  /**
   * Establish whether a polytope is empty
   *
   * Due to approximation errors, it may return false for some empty
   * polytopes too. However, when it returns true, the polytope is certainly
   * empty.
   *
   * @param[in] strict_inequality specifies whether the polytope is
   *         defined by a strict inequality (i.e., Ax < b).
   * @return a Boolean value. If the returned value is true, then the
   *       polytope is empty.
   */
  bool is_empty(const bool strict_inequality = false) const
  {
    return !this->has_solutions(strict_inequality);
  }

  /**
   * Check whether one polytope contains another polytope.
   *
   * This method establishes whether the current Polytope fully
   * contains another polytope. Due to the approximation errors,
   * the method may return false even if this is the case.
   * However, whenever it returns true, the current object
   * certaintly contains the polytope.
   *
   * @param[in] P is the polytope that are compared to the current
   *     object.
   * @return a Boolean value. When the current object does not
   *     contain the parameter, the retured value is false. When
   *     the method returns true, the current polytope contains
   *     the parameter. There are cases in which the current
   *     object contains the parameter and, still, this method
   *     returns false.
   */
  bool contains(const Polytope &P) const;

  /**
   *  Split a polytope in a list of polytopes.
   *
   *  This method splits a polytope in a list of polytopes such
   *  that their set union equals the original polytope.
   *
   *  @return A list of polytopes such that their union equals
   *      the original polytope.
   */
  std::list<Polytope> split() const;

  /**
   * Update a polytope by intersecting it with another one.
   *
   * This method works in-place and changes the calling object.
   *
   * @param[in] P is a polytope.
   * @return a reference to the updated object.
   */
  Polytope &intersect_with(const Polytope &P);

  /**
   * Determine the volume of the bounding box of the polytope
   *
   * @return volume of the bounding box
   */
  double bounding_box_volume() const;

  void plotRegion(std::ostream &os = std::cout, const char color = ' ') const;

  void plotRegionT(std::ostream &os, const double t) const;
  void plotRegion(std::ostream &os, const std::vector<int> &rows,
                  const std::vector<int> &cols) const;

  friend void swap(Polytope &P1, Polytope &P2);

  /**
   * Compute the intersection of two polytopes
   *
   * @param[in] P1 is a polytope
   * @param[in] P2 is a polytope
   * @return the intersection of the two parameters.
   */
  friend Polytope intersect(const Polytope &P1, const Polytope &P2);
};

inline void swap(Polytope &P1, Polytope &P2)
{
  swap(*(static_cast<LinearSystem *>(&P1)),
       *(static_cast<LinearSystem *>(&P2)));
}

#endif /* LINEARSYSTEM_H_ */
