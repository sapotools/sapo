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

#include "Approximation.h"
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
  T _obj_value;                      //!< Value of the objective on the optimum
  Status _status;                    //!< The status of the process

public:
  /**
   * @brief Build a new optimization result
   *
   * @param optimum is the optimum point
   */
  OptimizationResult(const LinearAlgebra::Vector<T> &optimum,
                     const T objective_value):
      _optimum(optimum),
      _obj_value(objective_value), _status(OPTIMUM_AVAILABLE)
  {
  }

  /**
   * @brief Build a new optimization result
   *
   * @param optimum is the optimum point
   */
  OptimizationResult(LinearAlgebra::Vector<T> &&optimum,
                     const T objective_value):
      _optimum(std::move(optimum)),
      _obj_value(objective_value), _status(OPTIMUM_AVAILABLE)
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
      SAPO_ERROR("this constructor requires a status "
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
      _optimum(orig._optimum), _obj_value(orig._obj_value),
      _status(orig._status)
  {
  }

  /**
   * @brief Move constructor
   *
   * @param orig is the template object
   */
  OptimizationResult(OptimizationResult<T> &&orig):
      _optimum(std::move(orig._optimum)),
      _obj_value(std::move(orig._obj_value)), _status(orig._status)
  {
  }

  /**
   * @brief Get the optimum
   *
   * @return The optimum
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
   * @brief Get the objective value on the optimum
   *
   * @return the objective function value on the optimum
   */
  inline const T &objective_value() const
  {
    if (_status != OPTIMUM_AVAILABLE) {
      SAPO_ERROR("this method can be called only if the object "
                 "status is OPTIMUM_AVAILABLE",
                 std::domain_error);
    }

    return _obj_value;
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
    _obj_value = orig._obj_value;
    _status = orig._status;

    return *this;
  }

  /**
   * @brief Assignment operator with move semantics
   *
   * @param orig is the original optimization result
   * @return a reference to the updated object
   */
  OptimizationResult<T> &operator=(OptimizationResult<T> &&orig)
  {
    _optimum = std::move(orig._optimum);
    _obj_value = std::move(orig._obj_value);
    _status = orig._status;

    return *this;
  }

  inline bool feasible_set_is_empty() const
  {
    return _status == INFEASIBLE;
  }

  inline bool feasible_set_is_unbounded() const
  {
    return _status == UNBOUNDED;
  }

  inline bool optimum_is_available() const
  {
    return _status == OPTIMUM_AVAILABLE;
  }
};

/**
 * @brief Optimization goal
 */
enum OptimizationGoal { MAXIMIZE, MINIMIZE };

/**
 * @brief An linear optimizer that implements the simplex method
 *
 * @tparam T is the type of the linear system coefficients
 * @tparam APPROX_TYPE is the type of approximation used in place of T
 *      during the computation
 */
template<typename T, typename APPROX_TYPE = T>
class SimplexMethodOptimizer
{
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
   * - indices in [0, num_cols(A)-1]: positive values for system variables
   * - indices in [num_cols(A), 2*num_cols(A)-1]: negative values for system
   * variable
   * - indices in [2*num_cols(A), 2*num_cols(A)+num_rows(A)-1]: slack
   * variables
   * - index 2*num_cols(A)+num_rows(A): the artificial variable
   * - index 2*num_cols(A)+num_rows(A)+1: the constant coefficient
   *
   * @param A is the linear system matrix
   * @param b is the linear system vector
   * @param objective is the objective coefficient vector
   * @param optimization_type is the type of optimization
   *      to achieve, i.e., MAXIMIZE or MINIMIZE
   */
  class Tableau
  {
    std::vector<LinearAlgebra::Vector<APPROX_TYPE>>
        _tableau;                         //!< simplex method tableau
    std::vector<size_t> _basic_variables; //!< basic variable vector

    /**
     * @brief Initialize the initial basic variable vector
     */
    void init_basic_variables(const bool &has_objective)
    {
      const size_t slack_begin = _tableau[0].size()
                                 - (_tableau.size() - (has_objective ? 2 : 1))
                                 - 3;
      _basic_variables.resize(0);
      std::generate_n(std::back_inserter(_basic_variables),
                      _tableau.size() - (has_objective ? 2 : 1),
                      [n = slack_begin]() mutable { return ++n; });
    }

    /**
     * @brief Add omega row
     */
    void add_omega_row()
    {
      LinearAlgebra::Vector<APPROX_TYPE> omega_row(_tableau[0].size(),
                                                   APPROX_TYPE(0));
      omega_row[last_artificial_column_index()] = APPROX_TYPE(1);
      _tableau.push_back(std::move(omega_row));
    }

    /**
     * @brief Initialize the tableau for feasibility test
     *
     * @param A is the linear system matrix
     * @param b is the linear system vector
     */
    void init_tableau(const std::vector<LinearAlgebra::Vector<T>> &A,
                      const LinearAlgebra::Vector<T> &b)
    {
      this->init_tableau(A, b, LinearAlgebra::Vector<APPROX_TYPE>());
    }

    /**
     * @brief Initialize the tableau for minimization
     *
     * @param A is the linear system matrix
     * @param b is the linear system vector
     * @param objective is the coefficient vector of the function to be
     *        minimized
     */
    void init_tableau(const std::vector<LinearAlgebra::Vector<T>> &A,
                      const LinearAlgebra::Vector<T> &b,
                      LinearAlgebra::Vector<APPROX_TYPE> &&objective)
    {
      const size_t num_rows = A.size();
      const size_t num_cols = (num_rows == 0 ? 0 : A[0].size());

      _tableau = std::vector<LinearAlgebra::Vector<APPROX_TYPE>>();
      for (const auto &row: A) {
        LinearAlgebra::Vector<APPROX_TYPE> approx_row;

        std::copy(row.begin(), row.end(), std::back_inserter(approx_row));

        _tableau.push_back(std::move(approx_row));
      }

      // add objective row
      if (objective.size() > 0) {
        _tableau.push_back(std::move(objective));
      }

      // complete tableau rows
      size_t row_idx = 0;
      for (auto &row: _tableau) {
        // double variable to handle negative solutions
        row.reserve(2 * num_cols + num_rows + 2);
        for (size_t i = 0; i < num_cols; ++i) {
          row.push_back(-row[i]);
        }

        // add slack/artificial variables
        for (size_t i = 0; i < num_rows; ++i) {
          row.push_back(APPROX_TYPE(row_idx == i ? 1 : 0));
        }

        // add the last artificial variable
        row.push_back(APPROX_TYPE(b[row_idx] < 0 ? -1 : 0));

        // add the constant coefficient
        row.push_back(APPROX_TYPE(row_idx < num_rows ? b[row_idx] : 0));

        ++row_idx;
      }
    }

    /**
     * @brief Let artificial variable in null constraints be non-basic
     *
     * After the bootstrap phase, some artificial variables may be basic
     * even though they can be replaced by non-artificial ones still
     * preserving the constant coefficient of the omega row. This occurs
     * when the row corresponding to one of artificial basic variable
     * has 0 as constant coefficient. For the sake of example, let us
     * consider the linear minimization problem:
     *
     *  objective: x
     *  subject to: x <= 1 && x >= 1
     *
     * which corresponds to the tableau:
     *
     *                    x -x s1 s2 a0 | b
     *                   ------------------
     *     constraint 1:  1 -1  1  0  0 | 1
     *     constraint 2: -1  1  0  1 -1 |-1
     *        objective: -1  1  0  0  0 | 0
     *                   ------------------
     *        omega row:  0  0  0  0  1 | 0
     *
     * with basic variables [s1, s2]
     *
     * After the bootstrap phase, the tableau becomes:
     *
     *                    x -x s1 s2 a0 | b
     *                   ------------------
     *     constraint 1:  1 -1  1  0  0 | 1
     *     constraint 2:  0  0 -1 -1  1 | 0
     *        objective:  0  0 -1  0  0 |-1
     *                   ------------------
     *        omega row:  0  0  1  1  0 | 0
     *
     * with basic variables [x, a0]. However, the artificial
     * variable a0 corresponds to the constraint 2 whose constant
     * coefficient 0. Thus, all the variables whose coefficients
     * in the row are non-null must equal 0 and the coefficient
     * signs can be reversed. It follows that s2 can replace a0
     * as a basic variable still maintaining the same omega
     * coefficient.
     */
    void let_artificial_variables_in_null_constraints_be_non_basic()
    {
      std::set<size_t> artificials_in_basics;
      std::vector<size_t> basic_vector(num_of_tableau_variables(), false);
      for (size_t i = 0; i < _basic_variables.size(); ++i) {
        if (is_artificial(_basic_variables[i])) {
          artificials_in_basics.insert(i);
        }
        basic_vector[_basic_variables[i]] = true;
      }

      for (auto &artificial_in_basics: artificials_in_basics) {
        const size_t artificial_column
            = _basic_variables[artificial_in_basics];

        auto &artificial_row = _tableau[artificial_in_basics];
        if (!is_false(artificial_row[b_column_index()] == 0)) {

          size_t entering = num_of_columns();
          for (size_t i = 0; i < first_artificial_column_index(); ++i) {
            // `artificial_row[i] == 0` may be tri-Boolean. In that case
            // `is_false(artificial_row[i] == 0)` !=
            // `!is_true(artificial_row[i] == 0)`
            if (!basic_vector[i] && is_false(artificial_row[i] == 0)) {
              entering = i;
            }
          }

          if (entering < num_of_columns()) {
            basic_vector[entering] = true;
            artificial_row[artificial_column] = APPROX_TYPE(0);
            _basic_variables[artificial_in_basics] = entering;
            pivot_operation(artificial_in_basics);
          }
        }
      }
    }

    /**
     * @brief Nullify the artificial variable coefficients
     */
    void nullify_artificial_variable_coefficients()
    {
      for (auto &row: _tableau) {
        const auto begin_it
            = std::begin(row) + first_artificial_column_index();
        const auto end_it
            = std::begin(row) + last_artificial_column_index() + 1;
        for (auto elem_it = begin_it; elem_it != end_it; ++elem_it) {
          *elem_it = APPROX_TYPE(0);
        }
      }
    }

  public:
    /**
     * @brief Get the tableau
     *
     * @return the tableau
     */
    inline std::vector<LinearAlgebra::Vector<APPROX_TYPE>> &get_tableau() const
    {
      return _tableau;
    }

    /**
     * @brief Get the basic variable vector
     *
     * @return the basic variable vector
     */
    inline const std::vector<size_t> &get_basic_variables() const
    {
      return _basic_variables;
    }

    /**
     * @brief Get the index of the first artificial variable column
     *
     * @return the index of the first artificial variable column
     */
    inline size_t first_artificial_column_index() const
    {
      return _tableau[0].size() - 2;
    }

    /**
     * @brief Get the index of the last artificial variable column
     *
     * @return the index of the last artificial variable column
     */
    inline size_t last_artificial_column_index() const
    {
      return _tableau[0].size() - 2;
    }

    /**
     * @brief Test whether a variable is artificial
     *
     * @param variable_index is the index of the variable to be
     *      tested
     * @return `true` if and only if the `variable_index`-th
     *      variable of the tableau is artificial
     */
    inline bool is_artificial(const size_t &variable_index) const
    {
      return variable_index >= first_artificial_column_index();
    }

    /**
     * @brief Get the index of the constant coefficient column
     *
     * @return the index of the constant coefficient column
     */
    inline size_t b_column_index() const
    {
      return _tableau[0].size() - 1;
    }

    /**
     * @brief Get the number of constraints in the linear system
     *
     * @return the number of constraints in the linear system
     */
    inline size_t num_of_system_variables() const
    {
      return (num_of_columns() - _basic_variables.size() - 2) / 2;
    }

    /**
     * @brief Get the number of tableau rows
     *
     * @return the number of tableau rows
     */
    inline size_t num_of_rows() const
    {
      return _tableau.size();
    }

    /**
     * @brief Get the number of tableau columns
     *
     * @return the number of tableau columns
     */
    inline size_t num_of_columns() const
    {
      return _tableau[0].size();
    }

    /**
     * @brief Get the number of tableau columns
     *
     * @return the number of tableau columns
     */
    inline size_t num_of_tableau_variables() const
    {
      return num_of_columns() - 1;
    }

    /**
     * @brief Get the objective row
     *
     * @return a reference to the objective row
     */
    inline LinearAlgebra::Vector<APPROX_TYPE> &objective_row()
    {
      return _tableau[num_of_rows() - 1];
    }

    /**
     * @brief Get the objective row
     *
     * @return a reference to the objective row
     */
    inline const LinearAlgebra::Vector<APPROX_TYPE> &objective_row() const
    {
      return _tableau[num_of_rows() - 1];
    }

    /**
     * @brief A constructor
     *
     * This constructor build a new Tableau object for the feasibility
     * test of the linear system \f$A * x <= b\f$.
     *
     * @param A is the linear system matrix
     * @param b is the linear system vector
     */
    Tableau(const std::vector<LinearAlgebra::Vector<T>> &A,
            const LinearAlgebra::Vector<T> &b)
    {
      using namespace LinearAlgebra;

      if (A.size() != b.size()) {
        SAPO_ERROR("A and b differ in the number of rows", std::domain_error);
      }

      init_tableau(A, b);

      add_omega_row();
      init_basic_variables(false);
    }

    /**
     * @brief A constructor
     *
     * This constructor build a new Tableau object for a the linear
     * optimization problem subject to the system \f$A * x <= b\f$.
     *
     * @param A is the linear system matrix
     * @param b is the linear system vector
     * @param objective is the objective coefficient vector
     * @param optimization_type is the required optimization goal, i.e.,
     *                          MINIMIZE or MAXIMIZE
     */
    Tableau(const std::vector<LinearAlgebra::Vector<T>> &A,
            const LinearAlgebra::Vector<T> &b,
            const LinearAlgebra::Vector<APPROX_TYPE> &objective,
            const OptimizationGoal optimization_type)

    {
      using namespace LinearAlgebra;

      const size_t num_rows = A.size();
      const size_t num_cols = (A.size() == 0 ? 0 : A[0].size());

      if (num_cols != objective.size()) {
        SAPO_ERROR("the objective dimension does not equal the number "
                   "of columns in A",
                   std::domain_error);
      }

      if (num_rows != b.size()) {
        SAPO_ERROR("A and b differ in the number of rows", std::domain_error);
      }

      if (optimization_type != MAXIMIZE && optimization_type != MINIMIZE) {
        SAPO_ERROR("unknown optimization type", std::runtime_error);
      }

      if (optimization_type == MAXIMIZE) {
        init_tableau(A, b, -objective);
      } else {
        init_tableau(A, b, Vector<APPROX_TYPE>(objective));
      }

      add_omega_row();
      init_basic_variables(true);
    }

    /**
     * @brief Perform a pivot operation on a tableau
     *
     * @param pivot_index is the index of the pivot in the basic
     *          variable vector
     */
    void pivot_operation(const size_t pivot_index)
    {
      using namespace LinearAlgebra;

      const size_t pivot_column_index = _basic_variables[pivot_index];
      Vector<APPROX_TYPE> &pivot_row = _tableau[pivot_index];
      const APPROX_TYPE pivot_value{pivot_row[pivot_column_index]};

      for (size_t row_index = 0; row_index < num_of_rows(); ++row_index) {
        if (row_index != pivot_index) {
          Vector<APPROX_TYPE> &row = _tableau[row_index];
          const APPROX_TYPE coeff = row[pivot_column_index];

          if (is_true(coeff != 0)) {

            // row = row - coeff * pivot_row / pivot_value;
            for (size_t i = 0; i < row.size(); ++i) {
              const auto term = coeff * pivot_row[i];

              // row[i] -= coeff * pivot_row[i] / pivot_value;
              if (is_true(row[i] * pivot_value == term)) {
                row[i] = APPROX_TYPE(0);
              } else {
                row[i] = row[i] - (term / pivot_value);
              }
            }
            row[pivot_column_index] = APPROX_TYPE(0);
          }
        }
      }

      pivot_row /= pivot_value;
      pivot_row[pivot_column_index] = APPROX_TYPE(1);
    }

    /**
     * @brief Find the most negative objective coefficient
     *
     * This function implements the "most negative coefficient rule"
     *
     * @return if there exists a negative coefficient in the objective
     *    function, the index of the most negative coefficient in the
     *    objective function. The number of columns in the tableau,
     *    otherwise.
     */
    size_t find_the_most_negative_obj_coefficient() const
    {
      const auto &obj_row = objective_row();
      const size_t num_of_vars = num_of_tableau_variables();

      size_t min_column = 0;
      for (size_t i = 1; i < num_of_vars; ++i) {
        if (obj_row[i] < obj_row[min_column]) {
          min_column = i;
        }
      }

      if (obj_row[min_column] < 0) {
        return min_column;
      }

      return num_of_columns();
    }

    /**
     * @brief Find the first negative objective coefficient
     *
     * This function implements the "Bland's rule"
     *
     * @return if there exists a negative coefficient in the objective
     *    function, the index of the first negative coefficient in the
     *    objective function. The number of columns in the tableau,
     *    otherwise.
     */
    size_t find_the_first_negative_obj_coefficient() const
    {
      const auto &obj_row = objective_row();
      const size_t num_of_vars = num_of_tableau_variables();
      for (size_t i = 0; i < num_of_vars; ++i) {
        if (is_true(obj_row[i] < 0)) {
          return i;
        }
      }

      return num_of_columns();
    }

    /**
     * @brief Choose the leaving variable
     *
     * @param entering_variable is the index in the tableau of the entering
     * variable
     * @return the index in the tableau of the leaving variable
     */
    size_t choose_leaving_variable(const size_t &entering_variable) const
    {
      auto coeff_col = b_column_index();
      size_t min_ratio_row = _basic_variables.size();
      APPROX_TYPE min_ratio{0};
      for (size_t i = 0; i < _basic_variables.size(); ++i) {
        const auto &row{_tableau[i]};

        if (is_true(row[entering_variable] > 0)) {
          const auto candidate_min = row[coeff_col] / row[entering_variable];

          if (min_ratio_row == _basic_variables.size()
              || is_true(candidate_min < min_ratio)) {
            min_ratio = candidate_min;
            min_ratio_row = i;
          }
        }
      }
      if (min_ratio_row == _basic_variables.size()) {
        return _basic_variables.size();
      }

      for (size_t i = 0; i < _basic_variables.size(); ++i) {
        if (is_true(_tableau[min_ratio_row][_basic_variables[i]] == 1)) {
          return i;
        }
      }

      return _basic_variables.size();
    }

    /**
     * @brief Extract the candidate linear problem solution
     *
     * @return the candidate linear problem solution
     */
    LinearAlgebra::Vector<APPROX_TYPE> get_candidate_solution() const
    {
      std::vector<APPROX_TYPE> candidate_solution(num_of_system_variables(),
                                                  APPROX_TYPE(0));

      for (size_t i = 0; i < _basic_variables.size(); ++i) {
        const size_t var_column = _basic_variables[i];
        if (var_column < 2 * num_of_system_variables()) {
          const size_t var_idx = var_column % num_of_system_variables();
          candidate_solution[var_idx]
              = (var_column >= num_of_system_variables()
                     ? -_tableau[i][b_column_index()]
                     : _tableau[i][b_column_index()]);
        }
      }
      return candidate_solution;
    }

    /**
     * @brief Minimize the objective
     *
     * This method minimize the tableau objective by repeating the
     * follow four steps:
     * 1. identify an entering variable
     * 2. identify the corresponding leaving variable
     * 3. replace the leaving variable by the entering variable in
     *    the basic variable vector
     * 4. apply a pivot operation on the entering variable
     *
     * The four steps are repeated until either:
     * - any entering variable are available: in this case, the
     *   method returns `OptimizationResult<T>::OPTIMUM_AVAILABLE`
     * - any leaving variable are available: in this case, the
     *   method returns `OptimizationResult<T>::UNBOUNDED`
     *
     * @return either `OptimizationResult<T>::OPTIMUM_AVAILABLE`
     *      (when no new entering variable are available) or
     *      `OptimizationResult<T>::UNBOUNDED` (when no new
     *      entering variable are available)
     */
    typename OptimizationResult<APPROX_TYPE>::Status minimize_objective()
    {
      while (true) {
        // using the Bald's rule to avoid loops
        size_t entering_variable = find_the_first_negative_obj_coefficient();

        if (entering_variable >= num_of_columns()) {
          return OptimizationResult<APPROX_TYPE>::OPTIMUM_AVAILABLE;
        }

        const size_t leaving_variable
            = choose_leaving_variable(entering_variable);

        if (leaving_variable >= _basic_variables.size()) {
          return OptimizationResult<APPROX_TYPE>::UNBOUNDED;
        }

        _basic_variables[leaving_variable] = entering_variable;

        pivot_operation(leaving_variable);
      }
    }

    /**
     * @brief Get the index of the minimum constant coefficient
     *
     * @return the index of the minimum constant coefficient
     */
    size_t get_the_minimum_coefficient_index() const
    {
      size_t min_coefficient_index = 0;
      auto min_coeff = _tableau[0][b_column_index()];
      for (size_t i = 1; i < _basic_variables.size(); ++i) {
        if (is_true(_tableau[i][b_column_index()] < min_coeff)) {
          min_coeff = _tableau[i][b_column_index()];
          min_coefficient_index = i;
        }
      }

      return min_coefficient_index;
    }

    /**
     * @brief Get the vertex currently considered by the tableau
     *
     * @return the vertex currently considered by the tableau
     */
    LinearAlgebra::Vector<APPROX_TYPE> get_current_vertex() const
    {
      LinearAlgebra::Vector<APPROX_TYPE> vertex(num_of_tableau_variables(), 0);

      size_t row = 0;
      for (const auto &variable: _basic_variables) {
        vertex[variable] = _tableau[row][b_column_index()];
        ++row;
      }

      return vertex;
    }

    /**
     * @brief Simplex method bootstrap phase
     *
     * The simplex method minimization phase jumps from one vertex of the
     * polytope defined by the linear system to another vertex reducing the
     * value of the objective function. Since all the tableau variables are
     * assumed to be non-negative, whenever all the system constant
     * coefficient are non-negative, the initial vertex is the null vector.
     * If instead some of the constants are negative, the bootstrap phase
     * is required.
     * This method implements the bootstrap phase and finds a valid set
     * of basic variables in problems dealing with constraints having
     * negative constant coefficient. At the end of the computation, it
     * also remove the omega row and nullify all the coefficients of
     * the artificial variables.
     *
     * @return if the problem solution set is certaintly not empty, this
     *      method returns `TriBool::TRUE`. When the set is certaintly
     *      empty, the method returns `TriBool::FALSE`. Otherwise,
     *      `TriBool::UNCERTAIN` is returned.
     */
    TriBool bootstrap()
    {
      const size_t entering_variable = last_artificial_column_index();
      size_t leaving_variable = get_the_minimum_coefficient_index();

      if (is_true(_tableau[leaving_variable][b_column_index()] < 0)) {
        _basic_variables[leaving_variable] = entering_variable;

        pivot_operation(leaving_variable);

        minimize_objective();

        let_artificial_variables_in_null_constraints_be_non_basic();
      }

      auto feasible = _tableau[num_of_rows() - 1][b_column_index()] == 0;

      nullify_artificial_variable_coefficients();

      // remove omega row
      _tableau.pop_back();

      return feasible;
    }
  };

public:
  /**
   * @brief An optimizer constructor
   */
  SimplexMethodOptimizer() {}

  /**
   * @brief Check whether a linear system is feasible
   *
   * This method verifies whether a linear system has a non-empty
   * set of solutions by using the simplex method.
   *
   * @param A is the linear system matrix
   * @param b is the linear system vector
   * @return if the problem solution set is certaintly not empty, this
   *      method returns `TriBool::TRUE`. When the set is certaintly
   *      empty, the method returns `TriBool::FALSE`. Otherwise,
   *      `TriBool::UNCERTAIN` is returned.
   */
  TriBool is_feasible(const std::vector<LinearAlgebra::Vector<T>> &A,
                      const LinearAlgebra::Vector<T> &b)
  {
    if (A.size() != b.size()) {
      SAPO_ERROR("the number of rows in A and that of "
                 "elements in b differ",
                 std::domain_error);
    }

    if (A.size() == 0) {
      return true;
    }

    Tableau tableau{A, b};

    return tableau.bootstrap();
  }

  /**
   * @brief Solves a linear programming problem
   *
   * This method solves a linear programming problem by using the simplex
   * method.
   *
   * @param A is the linear system matrix
   * @param b is the linear system vector
   * @param objective is the objective coefficient vector
   * @param optimization_type is the required optimization goal, i.e.,
   *                          minimize or maximize
   * @return the result of the optimization process
   */
  OptimizationResult<APPROX_TYPE>
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
      return OptimizationResult<APPROX_TYPE>(
          OptimizationResult<APPROX_TYPE>::UNBOUNDED);
    }

    LinearAlgebra::Vector<APPROX_TYPE> approx_objective;
    std::copy(objective.begin(), objective.end(),
              std::back_inserter(approx_objective));

    Tableau tableau{A, b, approx_objective, optimization_type};

    if (is_false(tableau.bootstrap())) {
      return OptimizationResult<APPROX_TYPE>(
          OptimizationResult<APPROX_TYPE>::INFEASIBLE);
    }

    const auto status = tableau.minimize_objective();
    switch (status) {
    case OptimizationResult<APPROX_TYPE>::OPTIMUM_AVAILABLE: {
      using namespace LinearAlgebra;

      auto candidate = tableau.get_candidate_solution();
      return {candidate, candidate * approx_objective};
    }
    case OptimizationResult<APPROX_TYPE>::UNBOUNDED:
      return OptimizationResult<APPROX_TYPE>::UNBOUNDED;
    default:
      SAPO_ERROR("unknown \"solve\" result", std::runtime_error);
    }
  }
};

#endif // SIMPLEX_H_