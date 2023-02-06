/**
 * @file Bernstein.h
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @author Alberto Casagrande <acasagrande@units.it>
 * @brief Computing Bernstein coefficients of rational polynomials
 * @version 0.1
 * @date 2022-07-13
 *
 * @copyright Copyright (c) 2015-2022
 */

#ifndef BERNSTEIN_H_
#define BERNSTEIN_H_

#include <vector>
#include <cmath>

#include "ErrorHandling.h"
#include "LinearAlgebra.h"
#include "SymbolicAlgebra.h"

template<typename C>
class Bernstein
{

private:
  /**
   * @brief Initialize the degree shift vector
   *
   * @param degrees is the vector of variable degrees in polynomial
   */
  static std::vector<unsigned int>
  get_shifts(const std::vector<unsigned int> &degrees)
  {
    std::vector<unsigned int> shifts(degrees.size(), 0);

    shifts[degrees.size() - 1] = degrees[degrees.size() - 1] + 1;

    for (int i = degrees.size() - 2; i >= 0; i--) {
      shifts[i] = (degrees[i] + 1) * shifts[i + 1];
    }

    return shifts;
  }

  /**
   * @brief Get the coefficients of a polynomial
   *
   * @param[out] coeffs is the destination coefficient vectors
   * @param[in] polynomial is the polynomial whose coefficients must be
   * evaluated
   * @param[in] vars is the vector of the variables in `polynomial`
   * @param[in] degrees is the vector of the maximal variable degrees
   * @param[in] shifts is the shift vector to store the coefficients
   * @param[in] idx is the index of the considered variable
   * @param[in] position is the position of the current coefficient in `coeffs`
   */
  static void
  get_coeffs_aux(std::vector<SymbolicAlgebra::Expression<C>> &coeffs,
                 const SymbolicAlgebra::Expression<C> &polynomial,
                 const std::vector<SymbolicAlgebra::Symbol<C>> &vars,
                 const std::vector<unsigned int> &degrees,
                 const std::vector<unsigned int> &shifts, unsigned int idx,
                 const unsigned int position)
  {
    if (polynomial == 0) {
      return;
    }

    auto p_coeffs = polynomial.get_coeffs(vars[idx]);

    // Base case, there's only one variable
    if (idx == vars.size() - 1) {

      for (unsigned int i = 0; i <= degrees[idx]; ++i) {
        auto found = p_coeffs.find(i);

        if (found != std::end(p_coeffs)) {
          coeffs[position + i] = std::move(found->second);
        }
      }

      return;
    }

    const unsigned int next_idx = idx + 1;
    SymbolicAlgebra::Expression<C> zero(0);
    SymbolicAlgebra::Expression<C> *p_coeff;
    for (unsigned int i = 0; i <= degrees[idx]; ++i) {
      auto found = p_coeffs.find(i);

      const unsigned int new_pos = position + i * shifts[next_idx];

      if (found == std::end(p_coeffs)) {
        p_coeff = &zero;
      } else {
        p_coeff = &(found->second);
      }

      get_coeffs_aux(coeffs, *p_coeff, vars, degrees, shifts, next_idx,
                     new_pos);
    }
  }

  /**
   * @brief Get the coefficients of a polynomial
   *
   * @param[in] polynomial is the polynomial whose coefficients must be
   * evaluated
   * @param[in] vars is the vector of the variables in `polynomial`
   * @param[in] degrees is the vector of the maximal variable degrees
   * @param[in] shifts is the shift vector to store the coefficients
   * @return the polynomial coefficient vectors
   */
  static std::vector<SymbolicAlgebra::Expression<C>>
  get_coeffs(const SymbolicAlgebra::Expression<C> &polynomial,
             const std::vector<SymbolicAlgebra::Symbol<C>> &vars,
             const std::vector<unsigned int> &degrees,
             const std::vector<unsigned int> &shifts)
  {
    using namespace SymbolicAlgebra;

    std::vector<Expression<C>> coeffs(shifts[0], 0);

    get_coeffs_aux(coeffs, polynomial, vars, degrees, shifts, 0, 0);

    return coeffs;
  }

  /**
   * @brief Compute the degrees of the variables in a polynomial
   *
   * @param p is a polynomial expression
   * @param q is a polynomial expression
   * @param vars is the vector of variables whose degrees must be computed
   * @return the vector of the degree
   */
  static std::vector<unsigned int>
  get_degree(const SymbolicAlgebra::Expression<C> &p,
             const SymbolicAlgebra::Expression<C> &q,
             const std::vector<SymbolicAlgebra::Symbol<C>> &vars)
  {
    std::vector<unsigned int> degrees;

    degrees.reserve(vars.size());
    for (auto var_it = std::begin(vars); var_it != std::end(vars); ++var_it) {
      degrees.push_back(std::max(p.degree(*var_it), q.degree(*var_it)));
    }

    return degrees;
  }

  /**
   * @brief Compute the degrees of the variables in a polynomial
   *
   * @param polynomial is a polynomial expression
   * @param vars is the vector of variables whose degrees must be computed
   * @return the vector of the degree
   */
  static std::vector<unsigned int>
  get_degree(const SymbolicAlgebra::Expression<C> &polynomial,
             const std::vector<SymbolicAlgebra::Symbol<C>> &vars)
  {
    std::vector<unsigned int> degrees;

    degrees.reserve(vars.size());
    for (auto var_it = std::begin(vars); var_it != std::end(vars); ++var_it) {
      degrees.push_back(polynomial.degree(*var_it));
    }

    return degrees;
  }

  /**
   * Decode a position into a multi-index
   *
   * @param[in] degrees is the vector of variable degrees
   * @param[in] shifts is the vector of shifts
   * @param[in] position is position to convert
   * @returns converted position
   */
  static std::vector<unsigned int>
  pos2multi_index(const std::vector<unsigned int> &degrees,
                  const std::vector<unsigned int> &shifts,
                  unsigned int position)
  {

    std::vector<unsigned int> multi_index(degrees.size());

    for (unsigned int i = 0; i < multi_index.size() - 1; i++) {
      multi_index[i] = (int)position / shifts[i + 1];
      position = position % shifts[i + 1];
    }

    multi_index[multi_index.size() - 1] = position;

    return multi_index;
  }

  /**
   * Calculate the binomial coefficient with multiplicative formula
   *
   * @param[in] n n
   * @param[in] k k
   * @returns n choose k
   */
  static unsigned int nChoosek(unsigned int n, unsigned int k)
  {

    if (k > n) {
      SAPO_ERROR("n must be larger equal then k", std::domain_error);
    }

    unsigned int res = 1;
    for (unsigned int i = 1; i <= k; ++i) {
      res = res * (n - (k - i)) / i;
    }
    return res;
  }

  /**
   * Calculate the binomial coefficient of two multi-indices
   *
   * @param[in] n upper multi-index
   * @param[in] k lower multi-index
   * @returns n choose k
   */
  static unsigned int multi_index_nChoosek(const std::vector<unsigned int> &n,
                                           const std::vector<unsigned int> &k)
  {

    if (n.size() != k.size()) {
      SAPO_ERROR("n and k must have same dimension", std::domain_error);
    }

    unsigned int res = 1;
    for (unsigned int i = 0; i < n.size(); i++) {
      res = res * nChoosek(n[i], k[i]);
    }

    return res;
  }

  /**
   * Determines whether b dominates
   *
   * @param[in] a multi-index
   * @param[in] b multi-index
   * @returns true if a <= b
   */
  static bool multi_index_leq(const std::vector<unsigned int> &a,
                              const std::vector<unsigned int> &b)
  {

    bool leq = true;
    unsigned int i = 0;

    while (i < a.size() && leq) {
      leq = leq && (a[i] <= b[i]);
      i++;
    }
    return leq;
  }

  /**
   * Productory of the elements of a vector within an interval
   *
   * @param[in] v vector with elements to multiply
   * @param[in] a beginning of the interval
   * @param[in] b end of the interval
   * @returns product v[a]v[a+1]...v[b]
   */
  static unsigned int prod(const std::vector<unsigned int> &v,
                           const unsigned int &a, const unsigned int &b)
  {
    unsigned int prod = 1;
    for (unsigned int i = a; i < b; i++) {
      prod = prod * v[i];
    }
    return prod;
  }

  /**
   * Shift a vector by one (rotate)
   *
   * @param[in] v vector to shift
   * @returns shifted vector
   */
  static std::vector<unsigned int> shift(const std::vector<unsigned int> &v)
  {
    std::vector<unsigned int> sv(v.size(), 0);

    for (unsigned int i = 1; i < v.size(); i++) {
      sv[i - 1] = v[i];
    }
    sv[v.size() - 1] = v[0];
    return sv;
  }

  /**
   * Convert an nd matrix into a 2d one
   *
   * @param[in] a matrix to convert
   * @param[in] dim is the dimension of the matrix
   * @returns 2d converted matrix
   */
  static std::vector<unsigned int> n2t(const std::vector<unsigned int> &a,
                                       const std::vector<unsigned int> &dim)
  {
    using namespace std;

    if (a.size() != dim.size()) {
      SAPO_ERROR("a and dim must have the same sizes", std::domain_error);
    }

    vector<unsigned int> b(2, 0);
    b[0] = a[0];
    if (a.size() > 1) {
      b[1] = a[1];

      for (unsigned int i = 2; i < a.size(); i++) {
        b[1] = b[1] + (a[i] * prod(dim, 1, i));
      }
    } else {
      b[1] = 0;
    }
    return b;
  }

  /**
   * Get the first component of the nd coordinate transpose
   *
   * @param[in] b_1 is the second component of the nd coordinate to be
   * transposed
   * @param[in] deg_1 dimensions of the second_component
   * @returns the first component of the transposed coordinate
   */
  static inline unsigned int first_transp(const unsigned int &b_1,
                                          const unsigned int &deg_1)
  {
    return b_1 % deg_1;
  }

  /**
   * Get the second component of the nd coordinate transpose
   *
   * @param[in] b_0 is the first component of the nd coordinate to be
   * transposed
   * @param[in] b_1 is the second component of the nd coordinate to be
   * transposed
   * @param[in] deg_1 dimensions of the second_component
   * @param[in] dim_prod product of the dimensions (prod(dim))
   * @returns the second component of the transposed coordinate
   */
  static inline unsigned int second_transp(const unsigned int &b_0,
                                           const unsigned int &b_1,
                                           const unsigned int &deg_1,
                                           const unsigned int &dim_prod)
  {
    return ((b_1 - (b_1 % deg_1)) / deg_1) + b_0 * dim_prod;
  }

  /**
   * Transpose an nd coordinate
   *
   * @param[in] b nd coordinate to transpose
   * @param[in] dim dimensions of the coordinate
   * @param[in] dim_prod product of the dimensions (prod(dim))
   * @returns transposed coordinate
   */
  static inline std::vector<unsigned int>
  transp(const std::vector<unsigned int> &b,
         const std::vector<unsigned int> &dim, const unsigned int &dim_prod)
  {
    return std::vector<unsigned int>{
        first_transp(b[1], dim[1]),
        second_transp(b[0], b[1], dim[1], dim_prod)};
  }

  /**
   * Convert a 2d matrix into a nd one
   *
   * @param[in] c matrix to convert
   * @param[in] dim dimensions of the matrix
   * @returns nd converted matrix
   */
  static std::vector<unsigned int> t2n(const std::vector<unsigned int> &c,
                                       const std::vector<unsigned int> &dim)
  {

    std::vector<unsigned int> a(dim.size(), 0);

    a[dim.size() - 1] = floor(c[1] / prod(dim, 1, dim.size() - 1));
    unsigned int c_value = c[1];
    for (unsigned int i = dim.size() - 1; i > 0; i--) {
      unsigned int div = prod(dim, 1, i);
      a[i] = floor(c_value / div);
      c_value = c_value % div;
    }
    a[0] = c[0];

    return a;
  }

  /**
   * Transpose an 2d matrix
   *
   * @param[in] M matrix to transpose
   * @param[in] dim dimensions of the matrix
   * @returns transposed matrix
   */
  static std::vector<std::vector<SymbolicAlgebra::Expression<C>>>
  transp(const std::vector<std::vector<SymbolicAlgebra::Expression<C>>> &M,
         const std::vector<unsigned int> &dim)
  {
    using namespace std;
    using namespace SymbolicAlgebra;

    const unsigned int prod_dim2n = prod(dim, 2, dim.size());

    const unsigned int rows_t = (dim.size() > 1 ? dim[1] : 1);
    const unsigned int cols_t = prod(dim, 2, dim.size()) * dim[0];

    vector<vector<Expression<C>>> M_transp(rows_t,
                                           vector<Expression<C>>(cols_t, 0.0));

    for (unsigned int i = 0; i < M.size(); i++) {
      const std::vector<SymbolicAlgebra::Expression<C>> &M_row = M[i];
      for (unsigned int j = 0; j < M_row.size(); j++) {
        if (M_row[j] != 0) {
          const unsigned int ij_t_0 = first_transp(j, rows_t);
          const unsigned int ij_t_1 = second_transp(i, j, rows_t, prod_dim2n);

          M_transp[ij_t_0][ij_t_1] = M_row[j];
        }
      }
    }

    return M_transp;
  }

  /**
   * Generate the U tilde matrix for improved matrix method
   *
   * @param[in] n dimension of the matrix
   * @returns U tilde matrix
   */
  static std::vector<std::vector<SymbolicAlgebra::Expression<C>>>
  genUtilde(const unsigned int &n)
  {
    using namespace std;
    using namespace SymbolicAlgebra;

    vector<vector<Expression<C>>> U(n + 1, vector<Expression<C>>(n + 1, 0));

    for (unsigned int i = 0; i < n + 1; i++) {
      U[i][0] = 1;
      U[n][i] = 1;
    }

    for (unsigned int i = 1; i < n; i++) {
      for (unsigned int j = 1; j <= i; j++) {
        U[i][j] = (C)nChoosek(i, i - j);
        U[i][j] /= (C)nChoosek(n, j);
      }
    }

    return U;
  }

  //! @private
  static std::vector<unsigned int>
  increase_by_one(const std::vector<unsigned int> &orig)
  {
    std::vector<unsigned int> increased(orig);
    for (auto it = std::begin(increased); it != std::end(increased); ++it) {
      ++(*it);
    }

    return increased;
  }

  /**
   * @brief Initialize the matrix of Bernstein coefficients
   *
   * @param coeffs is the vector of polynomial coefficients
   * @param degrees is the vector of variable degrees
   * @param degrees_p is the vector of variable degrees plus one
   * @param shifts is the vector of degree shifts
   * @return the initial matrix of Bernstein coefficients
   */
  static LinearAlgebra::Dense::Matrix<SymbolicAlgebra::Expression<C>>
  init_matrix_coeffs(const std::vector<SymbolicAlgebra::Expression<C>> &coeffs,
                     const std::vector<unsigned int> &degrees,
                     const std::vector<unsigned int> &degrees_p,
                     const std::vector<unsigned int> &shifts)
  {
    using namespace SymbolicAlgebra;
    using namespace LinearAlgebra::Dense;

    std::vector<Expression<C>> Ai(prod(degrees_p, 1, degrees_p.size()), 0);
    Matrix<Expression<C>> A(degrees_p[0], Ai);

    for (unsigned int i = 0; i < coeffs.size(); ++i) {
      if (coeffs[i] != 0) {
        auto pos2d = n2t(pos2multi_index(degrees, shifts, i), degrees_p);
        A[pos2d[0]][pos2d[1]] = coeffs[i];
      }
    }

    return A;
  }

  /**
   * @brief Compute the list of Bernstein coefficients
   *
   * @param[in] vars is the variable vector
   * @param[in] polynomial is the polynomial whose coefficient are computed
   * @param[in] degrees is the vector of the variable degrees
   * @param[in] shifts is the degree shift vector
   * @return the Bernstein coefficients of `polynomial`
   */
  static std::vector<SymbolicAlgebra::Expression<C>>
  get_coefficients_from_polynomial(
      const std::vector<SymbolicAlgebra::Symbol<C>> &vars,
      const SymbolicAlgebra::Expression<C> &polynomial,
      const std::vector<unsigned int> &degrees,
      const std::vector<unsigned int> &shifts)
  {
    using namespace SymbolicAlgebra;
    using namespace LinearAlgebra;

    // Initialize the coefficients vector
    auto coeffs = get_coeffs(polynomial, vars, degrees, shifts);

    // degrees increased by one
    auto degrees_p = increase_by_one(degrees);

    // initialize the matrix for the coefficients
    auto A = init_matrix_coeffs(coeffs, degrees, degrees_p, shifts);

    auto UAt = transp(genUtilde(degrees[0]) * A, degrees_p);
    for (long unsigned int i = 1; i < degrees.size(); i++) {
      degrees_p = shift(degrees_p);
      UAt = transp(genUtilde(degrees[i]) * UAt, degrees_p);
    }

    std::vector<Expression<C>> bernCoeffs;
    bernCoeffs.reserve(coeffs.size());

    for (long unsigned int i = 0; i < UAt.size(); i++) {
      for (long unsigned int j = 0; j < UAt[i].size(); j++) {

        // TODO: remove expand to speed-up computation
        if (UAt[i][j] != 0) {
          bernCoeffs.push_back(UAt[i][j].expand());
        } else {
          bernCoeffs.push_back(UAt[i][j]);
        }
      }
    }

    return bernCoeffs;
  }

public:
  /**
   * @brief Compute the Bernstein coefficients of a polynomial
   *
   * @param[in] vars is the vector of the variables
   * @param[in] polynomial is the polynomial whose coefficient are computed
   * @return the Bernstein coefficients of `polynomial`
   */
  static std::vector<SymbolicAlgebra::Expression<C>>
  get_coefficients_from_polynomial(
      const std::vector<SymbolicAlgebra::Symbol<C>> &vars,
      const SymbolicAlgebra::Expression<C> &polynomial)
  {
    if (!polynomial.is_a_polynomial()) {
      SAPO_ERROR("the parameter is not a polynomial", std::domain_error);
    }
    // Put the polynomial in extended form and extract variables degrees
    auto degrees = get_degree(polynomial, vars);

    // Initialize the degree shifts
    auto shifts = get_shifts(degrees);

    return get_coefficients_from_polynomial(vars, polynomial, degrees, shifts);
  }

  /**
   * @brief Compute the Bernstein coefficients of a rational polynomial
   *
   * @param[in] vars is the vector of the variables
   * @param[in] f is the rational polynomial whose coefficient are computed
   * @return the vector of Bernstein coefficients of `f`
   */
  static std::vector<SymbolicAlgebra::Expression<C>>
  get_coefficients(const std::vector<SymbolicAlgebra::Symbol<C>> &vars,
                   const SymbolicAlgebra::Expression<C> &f)
  {
    using namespace SymbolicAlgebra;
    Expression<C> p, q;

    if (!f.is_in_rational_form()) {
      auto f_rational = f.get_rational_form();

      p = f_rational.get_numerator();
      q = f_rational.get_denominator();
    } else {
      p = f.get_numerator();
      q = f.get_denominator();
    }

    // Put the polynomial in extended form and extract variables degrees
    auto degrees = get_degree(p, q, vars);

    // Initialize the degree shifts
    auto shifts = get_shifts(degrees);

    auto coeff = get_coefficients_from_polynomial(vars, p, degrees, shifts);

    if (q.is_a_constant() && q.evaluate() == 1) {
      return coeff;
    }

    auto q_coeff = get_coefficients_from_polynomial(vars, q, degrees, shifts);

    for (size_t i = 0; i < q_coeff.size(); ++i) {
      coeff[i] /= q_coeff[i];
    }

    return coeff;
  }
};

/**
 * @brief Compute the Bernstein coefficients of a polynomial
 *
 * @tparam C is the type of constants
 * @param[in] vars is the vector of the variables
 * @param[in] polynomial is the polynomial whose coefficient are computed
 * @return the Bernstein coefficients of `polynomial`
 */
template<typename C>
inline std::vector<SymbolicAlgebra::Expression<C>>
get_Bernstein_coefficients_from_polynomial(
    const std::vector<SymbolicAlgebra::Symbol<C>> &vars,
    const SymbolicAlgebra::Expression<C> &polynomial)
{
  return Bernstein<C>::get_coefficients_from_polynomial(vars, polynomial);
}

/**
 * @brief Compute the Bernstein coefficients of a rational polynomial
 *
 * @tparam C is the type of constants
 * @param[in] vars is the vector of the variables
 * @param[in] f is the rational polynomial whose coefficient are computed
 * @return the vector of Bernstein coefficients of `f`
 */
template<typename C>
inline std::vector<SymbolicAlgebra::Expression<C>>
get_Bernstein_coefficients(const std::vector<SymbolicAlgebra::Symbol<C>> &vars,
                           const SymbolicAlgebra::Expression<C> &f)
{
  return Bernstein<C>::get_coefficients(vars, f);
}
#endif // BERNSTEIN_H_
