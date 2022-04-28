/**
 * @file Polytope.cpp
 * Represent and manipulate polytopes (reached states, parameters, etc.)
 *
 * @author Alberto Casagrande <acasagrande@units.it>
 * @version 0.1
 */

#include "Polytope.h"
#include "PolytopesUnion.h"

#include "LinearAlgebra.h"

#define MAX_APPROX_ERROR 1e-8 // necessary for double comparison

using namespace std;

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
bool Polytope::contains(const Polytope &P) const
{

  for (unsigned int i = 0; i < this->size(); i++) {
    if (!P.satisfies(this->A[i], this->b[i])) {
      return false;
    }
  }

  return true;
}

template<typename T>
bool are_independent(const std::vector<T> &v1, const std::vector<T> &v2)
{
  if (v1.size() != v2.size()) {
    return true;
  }

  if (v1.size() == 0) {
    return false;
  }

  unsigned int fnz_v1(0);
  while (fnz_v1 < v1.size() && v1[fnz_v1] == 0)
    fnz_v1++;

  unsigned int fnz_v2(0);
  while (fnz_v2 < v2.size() && v2[fnz_v2] == 0)
    fnz_v2++;

  if (fnz_v1 != fnz_v2) {
    return true;
  }

  for (unsigned int i = 1; i < v1.size(); ++i) {
    if (v1[fnz_v1] * v2[i]
        != v2[fnz_v2] * v1[i]) { // this is to avoid numerical errors,
                                 // but it can produce overflows
      return true;
    }
  }

  return false;
}

std::vector<unsigned int>
get_a_polytope_base(const std::vector<std::vector<double>> &A)
{
  if (A.size() == 0) {
    return std::vector<unsigned int>();
  }

  std::vector<unsigned int> base{0};

  unsigned int row_idx = 1;
  while (row_idx < A.size()) {
    bool indep_from_base = true;
    auto b_it = std::begin(base);
    while (indep_from_base && b_it != std::end(base)) {
      indep_from_base = are_independent(A[row_idx], A[*b_it]);
      ++b_it;
    }

    if (indep_from_base) {
      base.push_back(row_idx);
    }
    ++row_idx;
  }

  return base;
}

inline std::vector<bool>
get_a_polytope_base_bit_vector(const std::vector<std::vector<double>> &A)
{
  std::vector<unsigned int> base = get_a_polytope_base(A);

  std::vector<bool> bvect_base(A.size(), false);

  for (auto it = std::begin(base); it != std::end(base); ++it) {
    bvect_base[*it] = true;
  }

  return bvect_base;
}

std::list<Polytope> Polytope::split(const std::vector<bool> &bvect_base,
                                    const unsigned int cidx,
                                    std::list<Polytope> &tmp_covering,
                                    std::vector<std::vector<double>> &A,
                                    std::vector<double> &b,
                                    const unsigned int num_of_splits) const
{
  using namespace LinearAlgebra;
  
  if (this->A.size() == cidx) {
    Polytope ls(A, b);
    ls.simplify();

    tmp_covering.push_back(ls);
    return tmp_covering;
  }

  if (bvect_base[cidx]) {
    A.push_back(this->A[cidx]);
    b.push_back(this->b[cidx]);

    try {
      const double min_value = minimize(this->A[cidx]).optimum();

      A.push_back(-this->A[cidx]);
      if (num_of_splits > 0) {
        const double avg_value = (this->b[cidx] + min_value) / 2;
        b.push_back(-avg_value);

        split(bvect_base, cidx + 1, tmp_covering, A, b, num_of_splits - 1);

        b[b.size() - 1] = -min_value;
        b[b.size() - 2] = avg_value;

        split(bvect_base, cidx + 1, tmp_covering, A, b, num_of_splits - 1);
      } else {
        b.push_back(-min_value);
        split(bvect_base, cidx + 1, tmp_covering, A, b, num_of_splits);
      }

      A.pop_back();
      b.pop_back();
    } catch (std::logic_error &e) {
      std::cerr << "This polytope is open." << std::endl;

      split(bvect_base, cidx + 1, tmp_covering, A, b, num_of_splits);
    }

    A.pop_back();
    b.pop_back();
  } else {
    split(bvect_base, cidx + 1, tmp_covering, A, b, num_of_splits);
  }

  return tmp_covering;
}

std::list<Polytope> Polytope::split(const unsigned int num_of_splits) const
{
  std::list<Polytope> result;

  std::vector<bool> bvect_base = get_a_polytope_base_bit_vector(this->A);

  std::vector<std::vector<double>> A;
  std::vector<double> b;

  split(bvect_base, 0, result, A, b, num_of_splits);

  return result;
}

/**
 * Update a polytope by intersecting it with another one.
 *
 * This method works in-place and changes the calling object.
 *
 * @param[in] P is a polytope.
 * @return a reference to the updated object.
 */
Polytope &Polytope::intersect_with(const Polytope &P)
{
  for (unsigned int i = 0; i < P.size(); i++) {
    if (!this->satisfies(P.A[i], P.b[i])) {
      (this->A).push_back(P.A[i]);
      (this->b).push_back(P.b[i]);
    }
  }

  return *this;
}

Polytope intersect(const Polytope &P1, const Polytope &P2)
{
  Polytope result(P1.A, P1.b);

  result.intersect_with(P2);

  return result;
}

/**
 * Determine the volume of the bounding box of the polytope
 *
 * @return volume of the bounding box
 */
double Polytope::bounding_box_volume() const
{

  vector<double> zeros(this->dim(), 0);
  double vol = 1;

  for (unsigned int i = 0; i < this->dim(); i++) {
    vector<double> facet = zeros;
    facet[i] = 1;
    const double b_plus = maximize(facet).optimum();
    facet[i] = -1;
    const double b_minus = minimize(facet).optimum();
    vol = vol * (b_plus + b_minus);
  }

  return vol;
}

/**
 * Print the polytope in Matlab format (for plotregion script)
 *
 * @param[in] os is the output stream
 * @param[in] color color of the polytope to plot
 */
void Polytope::plotRegion(std::ostream &os, const char color) const
{

  if (this->dim() > 3) {
    std::domain_error("Polytope::plotRegion : maximum 3d sets are allowed");
  }

  os << "Ab = [" << std::endl;
  for (unsigned i = 0; i < A.size(); i++) {
    for (auto el = std::begin(A[i]); el != std::end(A[i]); ++el) {
      os << *el << " ";
    }
    os << " " << this->b[i] << ";" << std::endl;
  }
  os << "];" << std::endl;
  os << "plotregion(-Ab(:,1:" << this->A[0].size() << "),-Ab(:,"
     << this->A[0].size() + 1 << "),[],[],";

  if (color == ' ') {
    os << "color";
  } else {
    os << color;
  }
  os << ");" << std::endl;
}

/**
 * Print the 2d polytope in Matlab format (for plotregion script) over
 * time
 *
 * @param[in] os is the output stream
 * @param[in] t thickness of the set to plot
 */
void Polytope::plotRegionT(std::ostream &os, const double t) const
{
  if (this->dim() > 2) {
    std::domain_error("Polytope::plotRegionT : maximum 2d sets are allowed");
  }

  os << "Ab = [" << std::endl;
  os << " 1 ";
  for (unsigned int j = 0; j < this->A[0].size(); j++) {
    os << " 0 ";
  }
  os << t << ";" << std::endl;
  os << " -1 ";
  for (unsigned int j = 0; j < this->A[0].size(); j++) {
    os << " 0 ";
  }
  os << -t << ";" << std::endl;

  for (unsigned i = 0; i < A.size(); i++) {
    os << " 0 ";
    for (auto el = std::begin(A[i]); el != std::end(A[i]); ++el) {
      os << *el << " ";
    }
    os << this->b[i] << ";" << std::endl;
  }

  os << "];" << std::endl;
  os << "plotregion(-Ab(:,1:3),-Ab(:,4),[],[],color);" << std::endl;
}

/**
 * Print the specified projections in Matlab format (for plotregion script)
 * into a file
 *
 * @param[in] os is the output stream
 * @param[in] rows rows to be plot
 * @param[in] cols colors of the plots
 */
void Polytope::plotRegion(std::ostream &os, const vector<int> &rows,
                          const vector<int> &cols) const
{

  if (cols.size() > 3) {
    std::domain_error(
        "Polytope::plotRegion : cols maximum 3d sets are allowed");
  }

  os << "Ab = [" << std::endl;
  for (auto r_it = std::begin(rows); r_it != std::end(rows); ++r_it) {
    for (auto c_it = std::begin(cols); c_it != std::end(cols); ++c_it) {
      os << this->A[*r_it][*c_it] << " ";
    }
    os << " " << this->b[*r_it] << ";" << std::endl;
  }
  os << "];" << std::endl;
  os << "plotregion(-Ab(:,1:" << cols.size() << "),-Ab(:," << cols.size() + 1
     << "),[],[],color);" << std::endl;
}
