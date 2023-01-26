/**
 * @file Simplex.h
 * @author Alberto Casagrande <acasagrande@units.it>
 * @brief A template implementation for the simplex algorithm
 * @version 0.1
 * @date 2023-01-26
 *
 * @copyright Copyright (c) 2023
 */

#ifndef SIMPLEX_H_
#define SIMPLEX_H_

#include "LinearAlgebra.h"

#include "ErrorHandling.h"

/**
 * @brief The result of an optimization process
 *
 * @tparam T is the type of the output value
 */
template<typename T>
class OptimizationResult
{
public:
  enum Status {
    OPTIMUM_AVAILABLE, //!< optimum available
    UNBOUNDED,         //!< the feasible set is unbounded
    INFEASIBLE         //!< the feasible set is empty
  };

private:
  LinearAlgebra::Vector<T> _optimum; //!< The optimum
  Status _status;                    //!< The status of the process

public:
  /**
   * @brief Build a new optimization result
   *
   * @param optimum is the optimum point
   */
  OptimizationResult(const LinearAlgebra::Vector<T> &optimum):
      _optimum(optimum), _status(OPTIMUM_AVAILABLE)
  {
  }

  /**
   * @brief Build a new optimization result
   *
   * @param optimum is the optimum point
   */
  OptimizationResult(LinearAlgebra::Vector<T> &&optimum):
      _optimum(std::move(optimum)), _status(OPTIMUM_AVAILABLE)
  {
  }

  /**
   * @brief Build a new optimization result
   *
   * @param status is the status of the optimization process
   */
  OptimizationResult(const Status status): _status(status)
  {
    if (status == OPTIMUM_AVAILABLE) {
      SAPO_ERROR("this construtor requires a status "
                 "different from OPTIMUM_AVAILABLE",
                 std::domain_error);
    }
  }

  /**
   * @brief Copy constructor
   *
   * @param orig is the template object
   */
  OptimizationResult(const OptimizationResult<T> &orig):
      _optimum(orig._optimum), _status(orig._status)
  {
  }

  /**
   * @brief Move constructor
   *
   * @param orig is the template object
   */
  OptimizationResult(OptimizationResult<T> &&orig):
      _optimum(std::move(orig._optimum)), _status(orig._status)
  {
  }

  /**
   * @brief Get the optimum value
   *
   * @return The optimum value
   */
  inline const LinearAlgebra::Vector<T> &optimum() const
  {
    if (_status != OPTIMUM_AVAILABLE) {
      SAPO_ERROR("this method can be called only if the object "
                 "status is OPTIMUM_AVAILABLE",
                 std::domain_error);
    }

    return _optimum;
  }

  /**
   * @brief Get the status of the optimization process
   *
   * @return The status of the optimization process
   */
  inline const Status &status() const
  {
    return _status;
  }

  /**
   * @brief Assignment operator
   *
   * @param orig is the original optimization result
   * @return a reference to the updated object
   */
  OptimizationResult<T> &operator=(const OptimizationResult<T> &orig)
  {
    _optimum = orig._optimum;
    _status = orig._status;

    return *this;
  }

  /**
   * @brief Assignment operator
   *
   * @param orig is the original optimization result
   * @return a reference to the updated object
   */
  OptimizationResult<T> &operator=(OptimizationResult<T> &&orig)
  {
    _optimum = std::move(orig._optimum);
    _status = orig._status;

    return *this;
  }
};

/**
 * @brief Optimization goal
 */
enum OptimizationGoal { MAXIMIZE, MINIMIZE };

/**
 * @brief An linear optimizer that implement the simplex method
 */
class SimplexMethodOptimizer
{
  template<typename T>
  static inline int sign(const T &value)
  {
    if (0 < value) {
      return 1;
    }

    if (value < 0) {
      return -1;
    }

    return 0;
  }

  /**
   * @brief Build the simplex algorithm tableau for a linear system
   *
   * This method build the tableau used by the simplex algorithm to minimize
   * the objective function \f$\texttt{objective} \cdot x\f$ subject to the
   * linear system \f$A \cdot x \leq b\f$.
   *
   * The tableau will contains num_rows(A)+2 rows: the first num_rows(A) rows
   * will represent the linear system constraints and the last two correspond
   * to the objective function and the omega row, respectively. Each tableau
   * row will be organized as follows:
   * - indices in [1, num_cols(A)-1]: positive values for system variables
   * - indices in [num_cols(A), 2*num_cols(A)-1]: negative values for system
   * variable
   * - indices in [2*num_cols(A)-1, 2*num_cols(A)-1+num_rows(A)]: slack
   * variables
   * - indices in [2*num_cols(A)+num_rows(A), 2*num_cols(A)-1+2*num_rows(A)]:
   * as variables
   * - index 2*(num_cols(A)+num_rows(A)) constant value
   *
   * @tparam T is the type of tableau coefficients
   * @param A is the linear system matrix
   * @param b is the linear system vector
   * @param objective is the objective coefficient vector
   * @return a tableau for the linear
   */
  template<typename T>
  static std::vector<LinearAlgebra::Vector<T>>
  build_tableau(const std::vector<LinearAlgebra::Vector<T>> &A,
                const LinearAlgebra::Vector<T> &b,
                const LinearAlgebra::Vector<T> &objective,
                OptimizationGoal optimization_type)
  {
    using namespace LinearAlgebra;

    size_t num_rows = A.size();

    size_t num_cols = (A.size() == 0 ? 0 : A[0].size());
    if (num_cols != objective.size()) {
      SAPO_ERROR("the objective dimension does not equal the number "
                 "of columns in A",
                 std::domain_error);
    }

    std::vector<LinearAlgebra::Vector<T>> tableau{A};
    switch (optimization_type) {
    case MINIMIZE:
      tableau.push_back(objective);
      break;
    case MAXIMIZE:
      tableau.push_back(-objective);
      break;
    default:
      SAPO_ERROR("unknown optimization type", std::runtime_error);
    }
    size_t row_idx = 0;
    for (auto &row: tableau) {

      // double variable to handle negative solutions
      row.reserve(2 * num_cols);
      for (size_t i = 0; i < num_cols; ++i) {
        row.push_back(-row[i]);
      }

      // add slack variables
      for (size_t i = 0; i < num_rows; ++i) {
        row.push_back(row_idx == i ? 1 : 0);
      }

      // add as variables
      for (size_t i = 0; i < num_rows; ++i) {
        row.push_back(row_idx == i ? sign(b[i]) * 1 : 0);
      }

      row.push_back(row_idx < num_rows ? b[row_idx] : 0);

      ++row_idx;
    }

    // add omega row
    LinearAlgebra::Vector<T> omega_row(tableau[0].size(), 0);
    for (size_t i = 0; i < num_rows; ++i) {
      omega_row[2 * num_cols + num_rows + i] = 1;
    }
    tableau.push_back(std::move(omega_row));

    return tableau;
  }

  /**
   * @brief Perform a pivot operation on a tableau
   *
   * @tparam T is the type of tableau coefficients
   * @param tableau is the tableau
   * @param basic_variables is the vector of basic variables
   * @param pivot_index is the index of the pivot in the basic
   *          variable vector
   */
  template<typename T>
  static void pivot_operation(std::vector<LinearAlgebra::Vector<T>> &tableau,
                              const std::vector<size_t> &basic_variables,
                              const size_t pivot_index)
  {
    using namespace LinearAlgebra;

    const size_t pivot_column_index = basic_variables[pivot_index];
    Vector<T> &pivot_row = tableau[pivot_index];
    const T pivot_value{pivot_row[pivot_column_index]};

    if (pivot_value != 1) {
      pivot_row = pivot_row / pivot_value;
    }

    for (size_t row_index = 0; row_index < tableau.size(); ++row_index) {
      Vector<T> &row = tableau[row_index];
      const T &coeff = row[pivot_column_index];
      if (coeff != 0 && row_index != pivot_index) {
        row = row - coeff * pivot_row;
      }
    }
  }

  /**
   * @brief Test the order of two tableau columns
   *
   * This function tests the order of two tableau columns with respect to their
   * omega and objective coefficients. It returns true if and only if the omega
   * of the first column is lesser than the omega of the second column or the
   * two columns have the same omega value and the objective coefficient of the
   * former is lesser of that of the latter.
   *
   * @tparam T is the type of tableau coefficients
   * @param tableau is a simplex method tableau as build by `build_tableau`
   * @param a is the index in the tableau of the first column to compare
   * @param b is the index in the tableau of the second column to compare
   * @return `true` if and only if the omega of the `a`-th column is lesser
   * than the omega of the `b`-th column or the two columns have the same omega
   *         value and the objective coefficient of the former is lesser of
   * that of the latter
   */
  template<typename T>
  static inline bool
  tableau_obj_less(const std::vector<LinearAlgebra::Vector<T>> &tableau,
                   const size_t &a, const size_t &b)
  {
    const LinearAlgebra::Vector<T> &omega_row = tableau[tableau.size() - 1];
    const LinearAlgebra::Vector<T> &obj_row = tableau[tableau.size() - 2];

    return ((omega_row[a] < omega_row[b])
            || (omega_row[a] == omega_row[b] && obj_row[a] < obj_row[b]));
  }

  /**
   * @brief Choose the entering variable
   *
   * @tparam T is the type of tableau coefficients
   * @param tableau is a simplex method tableau as build by `build_tableau`
   * @return the index in the tableau of the entering variable
   */
  template<typename T>
  static size_t choose_entering_variable(
      const std::vector<LinearAlgebra::Vector<T>> &tableau)
  {
    const size_t num_of_columns = tableau[0].size() - 1;
    size_t min_column = 0;
    for (size_t i = 1; i < num_of_columns; ++i) {
      if (tableau_obj_less(tableau, i, min_column)) {
        min_column = i;
      }
    }

    return min_column;
  }

  /**
   * @brief Choose the leaving variable
   *
   * @param tableau is a simplex method tableau as build by `build_tableau`
   * @return the index in the tableau of the entering variable
   */

  /**
   * @brief Choose the leaving variable
   *
   * @tparam T is the type of tableau coefficients
   * @param tableau is a simplex method tableau as build by `build_tableau`
   * @param entering_variable is the index in the tableau of the entering
   * variable
   * @param basic_variables is the vector of basic variables
   * @return the index in the tableau of the leaving variable
   */
  template<typename T>
  static size_t
  choose_leaving_variable(const std::vector<LinearAlgebra::Vector<T>> &tableau,
                          const size_t &entering_variable,
                          const std::vector<size_t> &basic_variables)
  {
    auto coeff_col = tableau[0].size() - 1;
    size_t min_positive_ratio_row = basic_variables.size();
    T min_positive_ratio = -1;
    for (size_t i = 0; i < basic_variables.size(); ++i) {
      const auto &row{tableau[i]};
      if (sign(row[coeff_col]) * sign(row[entering_variable]) > 0) {
        if (min_positive_ratio < 0
            || (std::abs(row[coeff_col])
                < std::abs(min_positive_ratio * row[entering_variable]))) {
          min_positive_ratio = row[coeff_col] / row[entering_variable];
          min_positive_ratio_row = i;
        }
      }
    }

    if (min_positive_ratio_row == basic_variables.size()) {
      return basic_variables.size();
    }

    for (size_t i = 0; i < basic_variables.size(); ++i) {
      if (tableau[min_positive_ratio_row][basic_variables[i]] == 1) {
        return i;
      }
    }

    return basic_variables.size();
  }

  /**
   * @brief Collext the optimum
   *
   * This method collects and returns the optimum value at the end of the
   * execution of the simplex method.
   *
   * @tparam T is the type of tableau coefficients
   * @param tableau is a simplex method tableau as build by `build_tableau`
   * @param basic_variables is the vector of basic variables
   * @param num_of_system_variables is the number of variables in the original
   * linear system
   * @return the optimal solution of the system
   */
  template<typename T>
  static LinearAlgebra::Vector<T>
  collect_optimum(const std::vector<LinearAlgebra::Vector<T>> &tableau,
                  const std::vector<size_t> &basic_variables,
                  const size_t num_of_system_variables)
  {
    std::vector<T> result(num_of_system_variables, 0);

    const size_t result_column = tableau[0].size() - 1;
    for (size_t i = 0; i < basic_variables.size(); ++i) {
      const size_t var_column = basic_variables[i];
      if (var_column < 2 * num_of_system_variables) {
        size_t j = 0;
        while (j < basic_variables.size() && tableau[j][var_column] != 1) {
          ++j;
        }
        const size_t var_idx = var_column % num_of_system_variables;
        result[var_idx] = (var_column >= num_of_system_variables
                               ? -tableau[j][result_column]
                               : tableau[j][result_column]);
      }
    }
    return result;
  }

  /**
   * @brief Build a initial basic variable vector
   *
   * @tparam T is the type of linear system coefficients
   * @param A is the linear system matrix
   * @return a vector containing the indexes of the basic variables in
   *      the correspong tableau
   */
  template<typename T>
  static std::vector<size_t>
  build_basic_variable_vector(const std::vector<LinearAlgebra::Vector<T>> &A)
  {
    const size_t omega_begin = 2 * A[0].size() + A.size() - 1;
    std::vector<size_t> basic_variables;
    std::generate_n(std::back_inserter(basic_variables), A.size(),
                    [n = omega_begin]() mutable { return ++n; });

    return basic_variables;
  }

  template<typename T>
  static bool
  problem_is_infeasable(const std::vector<LinearAlgebra::Vector<T>> &tableau,
                        const std::vector<size_t> &basic_variables)
  {
    const size_t last_artificial_index = tableau.front().size() - 2;
    const size_t first_artificial_index
        = last_artificial_index - tableau.size() + 3;

    for (size_t i = 0; i < basic_variables.size(); ++i) {
      if (first_artificial_index < basic_variables[i]
          && basic_variables[i] < last_artificial_index) {
        return true;
      }
    }

    return false;
  }

public:
  /**
   * @brief An optimizer constructor
   */
  SimplexMethodOptimizer() {}

  /**
   * @brief Solves a linear programming problem
   *
   * This method solves a linear programming problem by using the simplex
   * method.
   *
   * @tparam T is the type of the linear system coefficients
   * @param A is the linear system matrix
   * @param b is the linear system vector
   * @param objective is the objective coefficient vector
   * @param optimization_type is the required optimization goal, i.e.,
   *                          minimize or maximize
   * @return the result of the optimization process
   */
  template<typename T>
  OptimizationResult<T>
  operator()(const std::vector<LinearAlgebra::Vector<T>> &A,
             const LinearAlgebra::Vector<T> &b,
             const LinearAlgebra::Vector<T> &objective,
             OptimizationGoal optimization_type = MINIMIZE)
  {
    if (A.size() != b.size()) {
      SAPO_ERROR("the number of rows in A and that of "
                 "elements in b differ",
                 std::domain_error);
    }

    if ((A.size() == 0 && objective.size())
        || (A[0].size() != objective.size())) {
      SAPO_ERROR("objective size does not equal the number "
                 "of columns in A: they must be the same",
                 std::domain_error);
    }

    if (A.size() == 0) {
      return OptimizationResult<T>(OptimizationResult<T>::UNBOUNDED);
    }

    using namespace LinearAlgebra;

    auto tableau = build_tableau(A, b, objective, optimization_type);
    auto basic_variables = build_basic_variable_vector(A);

    std::cout << "tableau:" << std::endl
              << tableau << std::endl
              << "basic_variables:" << std::endl
              << basic_variables << std::endl;

    for (size_t pivot_row_index = 0; pivot_row_index < basic_variables.size();
         ++pivot_row_index) {
      pivot_operation(tableau, basic_variables, pivot_row_index);
    }

    const size_t obj_row_idx = tableau.size() - 2;
    const size_t omega_row_idx = tableau.size() - 1;
    const Vector<T> &omega_row = tableau[omega_row_idx];
    const Vector<T> &obj_row = tableau[obj_row_idx];

    while (true) {
      std::cout << "tableau:" << std::endl
                << tableau << std::endl
                << "basic_variables:" << std::endl
                << basic_variables << std::endl;
      size_t entering_variable = choose_entering_variable(tableau);

      if (omega_row[entering_variable] > 0
          || (omega_row[entering_variable] == 0
              && obj_row[entering_variable] >= 0)) {

        if (problem_is_infeasable(tableau, basic_variables)) {
          return OptimizationResult<T>(OptimizationResult<T>::INFEASIBLE);
        }

        auto optimum = collect_optimum(tableau, basic_variables, A[0].size());

        return OptimizationResult<T>(std::move(optimum));
      }

      const size_t leaving_variable = choose_leaving_variable(
          tableau, entering_variable, basic_variables);

      if (leaving_variable == basic_variables.size()) {
        if (problem_is_infeasable(tableau, basic_variables)) {
          return OptimizationResult<T>(OptimizationResult<T>::INFEASIBLE);
        }

        return OptimizationResult<T>(OptimizationResult<T>::UNBOUNDED);
      }

      basic_variables[leaving_variable] = entering_variable;

      pivot_operation(tableau, basic_variables, leaving_variable);
    }
  }
};

#endif // SIMPLEX_H_