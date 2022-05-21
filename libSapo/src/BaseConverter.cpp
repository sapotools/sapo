/**
 * @file BaseConverter.cpp
 * Convert the basis of the given polynomial
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @version 0.1
 */

#include "BaseConverter.h"

#include <cmath>

/**
 * Constructor that instantiates the base converter
 *
 * @param[in] vars list of variables appearing the current polynomial
 * @param[in] polynomial polynomial to be converted
 */
BaseConverter::BaseConverter(
    const std::vector<SymbolicAlgebra::Symbol<>> &vars,
    const SymbolicAlgebra::Expression<> &polynomial):
    vars(vars),
    polynomial(polynomial)
{
  this->polynomial.expand();

  // Put the polynomial in extended form and extract variables degrees
  for (auto var_it = std::begin(vars); var_it != end(vars); ++var_it) {
    this->degrees.push_back(this->polynomial.degree(*var_it));
  }

  // Initialize the degree shifts
  initShifts();

  this->coeffs
      = std::vector<SymbolicAlgebra::Expression<>>(this->shifts[0], 0);

  // Initialize the coefficients vector
  initCoeffs(this->polynomial);
}

/**
 * Constructor that instantiates the base converter knowing the variable
 * degrees
 *
 * @param[in] vars list of variables appearing the current polynomial
 * @param[in] polynomial polynomial to convert
 * @param[in] degrees of variables
 */
BaseConverter::BaseConverter(
    const std::vector<SymbolicAlgebra::Symbol<>> &vars,
    const SymbolicAlgebra::Expression<> &polynomial,
    const std::vector<unsigned int> &degrees):
    vars(vars),
    polynomial(polynomial), degrees(degrees)
{
  // Put the polynomial in extended form and extract variables degrees
  this->polynomial.expand();

  // Initialize the degree shifts
  initShifts();

  for (unsigned int i = 0; i < this->shifts[0]; i++) {
    this->coeffs.push_back(0);
  }

  // Initialize the coefficients vector
  initCoeffs(this->polynomial);
}

/**
 * @brief Constructor for rational polynomial base converters
 *
 * @param[in] vars list of variables appearing the current polynomial
 * @param[in] numerator of the rational polynomial to convert
 * @param[in] denominator of the rational polynomial to convert
 */
BaseConverter::BaseConverter(
    const std::vector<SymbolicAlgebra::Symbol<>> &vars,
    const SymbolicAlgebra::Expression<> &numerator,
    const SymbolicAlgebra::Expression<> &denominator):
    vars(vars),
    num(numerator), denom(denominator)
{
}

/**
 * Initialize the degree shift vector used to extract the multi-indices
 */
void BaseConverter::initShifts()
{
  this->shifts = std::vector<unsigned int>(this->degrees.size(), 0);

  this->shifts[this->degrees.size() - 1]
      = this->degrees[this->degrees.size() - 1] + 1;

  for (int i = this->degrees.size() - 2; i >= 0; i--) {
    this->shifts[i] = (this->degrees[i] + 1) * this->shifts[i + 1];
  }
}

/**
 * Convert a multi-index in a coefficient position index
 *
 * @param[in] multi_index multi-index to convert
 * @returns converted multi-index
 */
unsigned int BaseConverter::multi_index2pos(
    const std::vector<unsigned int> &multi_index) const
{
  using namespace std;

  if (multi_index.size() != this->degrees.size()) {
    throw std::domain_error("BaseConverter::multi_index2pos : multi_index "
                            "and degrees must have same dimension");
  }

  // Compute the position
  int position = multi_index[multi_index.size() - 1];
  for (unsigned int i = 0; i < multi_index.size() - 1; ++i) {
    position = position + multi_index[i] * this->shifts[i + 1];
  }

  return position;
}

/**
 * Decode a position into a multi-index
 *
 * @param[in] position position to convert
 * @returns converted position
 */
std::vector<unsigned int>
BaseConverter::pos2multi_index(unsigned int position) const
{

  std::vector<unsigned int> multi_index(this->degrees.size(), 0);

  for (unsigned int i = 0; i < multi_index.size() - 1; i++) {
    multi_index[i] = (int)position / this->shifts[i + 1];
    position = position % this->shifts[i + 1];
  }

  multi_index[multi_index.size() - 1] = position;

  return multi_index;
}

// TODO: Parallelize the following method.
void BaseConverter::initCoeffs(const SymbolicAlgebra::Expression<> &polynomial,
                               unsigned int var_idx,
                               const unsigned int position)
{
  if (polynomial == 0) {
    return;
  }

  auto p_coeffs = polynomial.get_coeffs(this->vars[var_idx]);

  // Base case, there's only one variable
  if (var_idx == this->vars.size() - 1) {

    for (unsigned int i = 0; i <= this->degrees[var_idx]; i++) {
      auto found = p_coeffs.find(i);

      if (found != std::end(p_coeffs)) {
        this->coeffs[position + i] = found->second;
      }
    }
  } else {
    const unsigned int next_idx = var_idx + 1;
    for (unsigned int i = 0; i <= this->degrees[var_idx]; i++) {
      auto found = p_coeffs.find(i);

      if (found == std::end(p_coeffs)) {
        SymbolicAlgebra::Expression<> zero(0);
        initCoeffs(zero, next_idx, position + i * this->shifts[next_idx]);
      } else {
        initCoeffs(found->second, next_idx,
                   position + i * this->shifts[next_idx]);
      }
    }
  }
}

/**
 * Calculate the binomial coefficient with multiplicative formula
 *
 * @param[in] n n
 * @param[in] k k
 * @returns n choose k
 */
unsigned int nChoosek(unsigned int n, unsigned int k)
{

  if (k > n) {
    throw std::domain_error("BaseConverter::nChoosek: n must be "
                            "larger equal then k");
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
unsigned int multi_index_nChoosek(const std::vector<unsigned int> &n,
                                  const std::vector<unsigned int> &k)
{

  if (n.size() != k.size()) {
    throw std::domain_error("BaseConverter::multi_index_nChoosek: n "
                            "and k must have same dimension");
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
bool multi_index_leq(const std::vector<unsigned int> &a,
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
 * Compute the mi-th Bernstein coefficient of the current polynomial
 *
 * @param[in] mi multi-index
 * @returns mi-th Bernstein coefficient
 */
SymbolicAlgebra::Expression<>
BaseConverter::bernCoeff(const std::vector<unsigned int> &mi) const
{
  using namespace SymbolicAlgebra;

  unsigned int i = multi_index2pos(mi);
  Expression<> coeff = 0;

  for (unsigned int j = 0; j <= i; j++) {

    if (this->coeffs[j] != 0) {

      std::vector<unsigned int> mj = pos2multi_index(j);
      if (multi_index_leq(mj, mi)) {

        int ichoosej = multi_index_nChoosek(mi, mj);
        int dchoosej = multi_index_nChoosek(this->degrees, mj);

        coeff += ((double)ichoosej * this->coeffs[j] / (double)dchoosej);
      }
    }
  }
  return coeff;
}

/**
 * Compute the list of Bernstein coefficients of the current polynomial
 *
 * @returns list of Bernstein coefficients
 */
std::vector<SymbolicAlgebra::Expression<>> BaseConverter::getBernCoeffs() const
{

  // cout<<"\tComputing Bernstein coefficients...\n";

  std::vector<SymbolicAlgebra::Expression<>> bern_coeffs;

  for (unsigned int i = 0; i < this->coeffs.size(); i++) {
    bern_coeffs.push_back(bernCoeff(pos2multi_index(i)));
  }

  return bern_coeffs;
}

/**
 * Compute the list of multi-indices of the current polynomial
 *
 * @returns list of multi-indices
 */
std::vector<std::vector<unsigned int>> BaseConverter::getMultiIdxList() const
{
  using namespace std;

  vector<vector<unsigned int>> mi_list;

  for (unsigned int i = 0; i < this->coeffs.size(); i++) {
    mi_list.push_back(pos2multi_index(i));
  }
  return mi_list;
}

/**
 * Refine the Bernstein coefficients with splitting algorithm
 *
 * @param[in] direction direction in which to split
 * @param[in] split_point splitting point
 */
void BaseConverter::split(unsigned int direction, double split_point) const
{
  using namespace std;

  if (direction >= this->vars.size()) {
    throw std::domain_error("BaseConverter::split: split direction must "
                            "be between 0 and space dimension excluded");
  }

  if ((split_point < 0) || (split_point > 1)) {
    throw std::domain_error("BaseConverter::split: split_point must be "
                            "between 0 and 1");
  }

  // Extract list of multi-indices and max degree of direction
  vector<vector<unsigned int>> multi_index_list = this->getMultiIdxList();

  std::vector<SymbolicAlgebra::Expression<>> B = this->getBernCoeffs();

  for (long unsigned int mi = 0; mi < B.size(); mi++) {

    vector<unsigned int> i(this->pos2multi_index(mi));

    for (unsigned int k = 1; k < this->degrees[direction]; k++) {

      for (unsigned int j = 0; j < this->degrees.size(); j++) {

        cout << "j:" << j << " dir: " << direction << "\n"
             << " dir: " << direction << "\n";

        if (j != direction) {

          for (unsigned int ij = 0; ij < this->degrees[j]; ij++) {

            if (ij < k) {
              // B[mi] = B[mi];
            } else {
              vector<unsigned int> ni = i;
              ni[j] = ij - 1;
              int bi = this->multi_index2pos(ni);
              B[mi] = (1 - split_point) * B[bi] + split_point * B[mi];
            }
          }
        }
      }
    }
  }
}

/**
 * Productory of the elements of a vector within an interval
 *
 * @param[in] v vector with elements to multiply
 * @param[in] a beginning of the interval
 * @param[in] b end of the interval
 * @returns product v[a]v[a+1]...v[b]
 */
unsigned int prod(const std::vector<unsigned int> &v, const unsigned int &a,
                  const unsigned int &b)
{
  unsigned int prod = 1;
  for (unsigned int i = a; i < b; i++) {
    prod = prod * v[i];
  }
  return prod;
}

// TODO: replace this function by using the one provided by linear algebra
// module.
/**
 * Multiplication of two 2d matrices
 *
 * @param[in] A left matrix to multiply
 * @param[in] B right matrix to multiply
 * @returns product A*B
 */
template<typename T>
std::vector<std::vector<T>> matrixProd(const std::vector<std::vector<T>> &A,
                                       const std::vector<std::vector<T>> &B)
{
  using namespace std;

  if (A[0].size() != B.size()) {
    throw std::domain_error("BaseConverter::matrixProd: matrices dimensions "
                            "must be the same");
  }

  vector<T> product_i(B[0].size(), 0);
  vector<vector<T>> product(A.size(), product_i);

  for (unsigned int i = 0; i < A.size(); i++) {
    for (unsigned int j = 0; j < B[0].size(); j++) {
      T inner_prod(0);
      for (long unsigned int k = 0; k < A[i].size(); k++) {
        inner_prod += A[i][k] * B[k][j];
      }
      std::swap(product[i][j], inner_prod);
    }
  }

  return product;
}

/**
 * Shift a vector by one (rotate)
 *
 * @param[in] v vector to shift
 * @returns shifted vector
 */
std::vector<unsigned int> shift(const std::vector<unsigned int> &v)
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
std::vector<unsigned int> n2t(const std::vector<unsigned int> &a,
                              const std::vector<unsigned int> &dim)
{
  using namespace std;

  if (a.size() != dim.size()) {
    throw std::domain_error("BaseConverter::n2t: a and dim must have "
                            "the same sizes");
  }

  vector<unsigned int> b(2, 0);
  b[0] = a[0];
  if (a.size()>1) {
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
 * @param[in] b_1 is the second component of the nd coordinate to be transposed
 * @param[in] deg_1 dimensions of the second_component
 * @returns the first component of the transposed coordinate
 */
inline unsigned int first_transp(const unsigned int &b_1,
                                 const unsigned int &deg_1)
{
  return b_1 % deg_1;
}

/**
 * Get the second component of the nd coordinate transpose
 *
 * @param[in] b_0 is the first component of the nd coordinate to be transposed
 * @param[in] b_1 is the second component of the nd coordinate to be transposed
 * @param[in] deg_1 dimensions of the second_component
 * @param[in] dim_prod product of the dimensions (prod(dim))
 * @returns the second component of the transposed coordinate
 */
inline unsigned int second_transp(const unsigned int &b_0,
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
inline std::vector<unsigned int> transp(const std::vector<unsigned int> &b,
                                        const std::vector<unsigned int> &dim,
                                        const unsigned int &dim_prod)
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
std::vector<unsigned int> t2n(const std::vector<unsigned int> &c,
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
 * Transpose an 2d coordinate
 *
 * @param[in] b 2d coordinate to transpose
 * @param[in] dim dimensions of the coordinate
 * @returns transposed coordinate
 */
std::vector<unsigned int> transp_naive(const std::vector<unsigned int> &b,
                                       const std::vector<unsigned int> &dim)
{

  return n2t(shift(t2n(b, dim)), shift(dim));
}

/**
 * Transpose an 2d matrix
 *
 * @param[in] M matrix to transpose
 * @param[in] dim dimensions of the matrix
 * @returns transposed matrix
 */
std::vector<std::vector<SymbolicAlgebra::Expression<>>>
transp(const std::vector<std::vector<SymbolicAlgebra::Expression<>>> &M,
       const std::vector<unsigned int> &dim)
{
  using namespace std;
  using namespace SymbolicAlgebra;

  const unsigned int prod_dim2n = prod(dim, 2, dim.size());

  const unsigned int rows_t = (dim.size()>1 ? dim[1]: 1);
  const unsigned int cols_t = prod(dim, 2, dim.size()) * dim[0];

  vector<vector<Expression<>>> M_transp(rows_t,
                                        vector<Expression<>>(cols_t, 0.0));

  for (unsigned int i = 0; i < M.size(); i++) {
    const std::vector<SymbolicAlgebra::Expression<>> &M_row = M[i];
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
std::vector<std::vector<SymbolicAlgebra::Expression<>>>
genUtilde(const unsigned int &n)
{
  using namespace std;
  using namespace SymbolicAlgebra;

  vector<vector<Expression<>>> U(n + 1, vector<Expression<>>(n + 1, 0.0));

  for (unsigned int i = 0; i < n + 1; i++) {
    U[i][0] = 1;
    U[n][i] = 1;
  }

  for (unsigned int i = 1; i < n; i++) {
    for (unsigned int j = 1; j <= i; j++) {
      U[i][j] = ((double)nChoosek(i, i - j)) / nChoosek(n, j);
    }
  }

  return U;
}

/**
 * Compute the list of Bernstein coefficients with improved matrix method
 *
 * @returns list of Bernstein coefficients
 */
std::vector<SymbolicAlgebra::Expression<>>
BaseConverter::getBernCoeffsMatrix() const
{
  using namespace std;
  using namespace SymbolicAlgebra;
  // cout<<"\tComputing Bernstein coefficients...\n";

  // degrees increased by one
  vector<unsigned int> degrees_p(this->degrees.size(), 0);
  for (unsigned int i = 0; i < degrees_p.size(); i++) {
    degrees_p[i] = this->degrees[i] + 1;
  }

  // initialize the matrix for the coefficients
  vector<Expression<>> Ai(prod(degrees_p, 1, degrees_p.size()), 0.0);
  vector<vector<Expression<>>> A(degrees_p[0], Ai);

  for (unsigned int i = 0; i < this->coeffs.size(); i++) {
    if (this->coeffs[i] != 0) {
      vector<unsigned int> pos2d = n2t(this->pos2multi_index(i), degrees_p);
      A[pos2d[0]][pos2d[1]] = this->coeffs[i];
    }
  }

  vector<vector<Expression<>>> UAt
      = transp(matrixProd(genUtilde(this->degrees[0]), A), degrees_p);
  for (long unsigned int i = 1; i < this->degrees.size(); i++) {
    degrees_p = shift(degrees_p);
    UAt = transp(matrixProd(genUtilde(this->degrees[i]), UAt), degrees_p);
  }

  std::vector<Expression<>> bernCoeffs;
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

/**
 * Print the list of computed Bernstein coefficients
 */
void BaseConverter::print() const
{
  using namespace std;

  for (unsigned int i = 0; i < this->coeffs.size(); i++) {
    if (this->coeffs[i] != 0) {
      cout << this->coeffs[i] << ": ";
      vector<unsigned int> multi_index = pos2multi_index(i);
      for (long unsigned int j = 0; j < multi_index.size(); j++) {
        cout << multi_index[j] << " ";
      }
      cout << "\n";
    }
  }
}

/**
 * Print the given matrix
 *
 * @param[in] M matrix to print
 */
template<typename T>
void print(const std::vector<std::vector<SymbolicAlgebra::Expression<>>> &M)
{
  using namespace std;

  for (long unsigned int i = 0; i < M.size(); i++) {
    for (long unsigned int j = 0; j < M[i].size(); j++) {
      cout << M[i][j] << " ";
    }
    cout << "\n";
  }
}

/**
 * @TODO Collection of functions for implicit computation of Bernstein
 * coefficients
 */

std::pair<std::vector<SymbolicAlgebra::Expression<>>,
          std::vector<std::vector<unsigned int>>>
BaseConverter::compressZeroCoeffs() const
{
  using namespace std;
  using namespace SymbolicAlgebra;

  vector<Expression<>> comp_coeffs;
  vector<vector<unsigned int>> comp_dim;
  pair<vector<Expression<>>, vector<vector<unsigned int>>> compression;

  for (unsigned int i = 0; i < this->coeffs.size(); i++) {
    if (this->coeffs[i] != 0) {
      comp_coeffs.push_back(this->coeffs[i]);
      comp_dim.push_back(pos2multi_index(i));
    }
  }

  compression.first = comp_coeffs;
  compression.second = comp_dim;
  return compression;
}

void BaseConverter::implicitMaxIndex() const
{
  using namespace std;
  using namespace SymbolicAlgebra;

  pair<vector<Expression<>>, vector<vector<unsigned int>>> compression
      = this->compressZeroCoeffs();
  // Check the three properties (Uniqueness,Monotonicity, and Dominance)
  vector<Expression<>> coeffs = compression.first;
  vector<vector<unsigned int>> multi_index = compression.second;
  vector<int> implicit_max(this->vars.size(), -1);

  // Check uniqueness
  vector<unsigned int> unique(this->vars.size(), 0);
  vector<unsigned int> unique_deg(this->vars.size(), -1);
  for (unsigned int i = 0; i < this->vars.size(); i++) {
    unsigned int j = 0;
    while ((j < coeffs.size()) && (unique[i] < 2)) {
      if (multi_index[j][i] > 0) {
        unique[i] = unique[i] + 1;
        unique_deg[i] = multi_index[j][i];
      }
      j++;
    }
  }
  for (unsigned int i = 0; i < unique.size(); i++) {
    if (unique[i] == 1) {
      implicit_max[i] = unique_deg[i];
    }
  }

  // Check monotonicity
  vector<int> increase(this->vars.size(), 0);
  for (long unsigned int i = 0; i < this->vars.size(); i++) {
    long unsigned int j = 0;
    while ((j < coeffs.size()) && (increase[i])) {
      if (multi_index[j][i] > 0) {
        increase[i] = increase[i] && (coeffs[j].evaluate() > 0);
      }
      j++;
    }
  }
  vector<int> decrease(this->vars.size(), 0);
  for (long unsigned int i = 0; i < this->vars.size(); i++) {
    long unsigned int j = 0;
    while ((j < coeffs.size()) && (decrease[i])) {
      if (multi_index[j][i] > 0) {
        decrease[i] = decrease[i] && (coeffs[j].evaluate() < 0);
      }
      j++;
    }
  }
  for (unsigned int i = 0; i < increase.size(); i++) {
    if (increase[i] == 1) {
      implicit_max[i] = this->degrees[i];
    }
    if (decrease[i] == 1) {
      implicit_max[i] = 0;
    }
  }

  // Dominance
}

BaseConverter::~BaseConverter()
{
  // TODO Auto-generated destructor stub
}
