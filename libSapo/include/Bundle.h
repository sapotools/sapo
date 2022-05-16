/**
 * @file Bundle.h
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @author Alberto Casagrande <acasagrande@units.it>
 * @brief Represent and manipulate bundles of parallelotopes 
 *        whose intersection represents a polytope
 * @version 0.2
 * @date 2022-05-04
 * 
 * @copyright Copyright (c) 2016-2022
 */

#ifndef BUNDLE_H_
#define BUNDLE_H_

#include "LinearAlgebra.h"
#include "BaseConverter.h"
#include "Polytope.h"
#include "PolytopesUnion.h"
#include "Parallelotope.h"
#include "VarsGenerator.h"

#include "STL/Atom.h"

#define SPLIT_MAGNITUDE_RATIO 0.75 //!< define the default versor magnitude multiplier for bundle splits

/**
 * @brief A class for parallelotope bundles
 * 
 * A paralletope bundle represents the intersection between 
 * different non-singular parallelotopes. 
 * This class stores all of the parallelotope directions/axes 
 * in one single array and constraints each of them between 
 * an upper and a lower bound. 
 * The original parallelotopes can be rebuild by using the 
 * a array of templates. Each template is a Natural-valued 
 * array corresponding to a parallelotope. It stores which 
 * directions are involved in the corresponding paralletope. 
 */ 
class Bundle
{
public:
  /**
   * @brief Approach to evaluate the image of a bundle
   * 
   * The are two different approches to evalutate the
   * image of a bundle through a polynomial function:
   * 1. One-For-One: the boundaries of any parallelotope
   *                 image in the bundle are evaluated by
   *                 exclusively considering the original
   *                 paralletope itself
   * 2. All-For-One: the boundaries of any parallelotope
   *                 image in the bundle are evaluated by
   *                 exploiting all the bundle templates
   */
  typedef enum {
    OFO,  /* One-For-One */
    AFO   /* All-For-One */
  } transfomation_mode;

private:
  std::vector<LinearAlgebra::Vector<double>> directions; //!< the vector of directions
  LinearAlgebra::Vector<double> lower_bounds;            //!< direction upper bounds
  LinearAlgebra::Vector<double> upper_bounds;            //!< direction lower bounds
  std::vector<LinearAlgebra::Vector<int>> templates;     //!< templates matrix

  /**
   * Compute the edge lengths
   *
   * This method compute the edge lengths, i.e., the distances between the
   * half-spaces bounding the bundle.
   *
   * @returns the vector of the edge lengths
   */
  LinearAlgebra::Vector<double> edge_lengths();

  /**
   * @brief A class to find the minimum and the maximum Bernstein coefficients
   */
  class MinMaxCoeffFinder
  {
    /**
     * @brief Evaluate a Bernstein coefficient
     *
     * @param bernCoeff is the symbolical representation of Bernstein
     *    coefficient
     * @return The numerical evaluation of Bernstein coefficient
     */
    double eval_coeff(const SymbolicAlgebra::Expression<> &bernCoeff) const;

  public:
    /**
     * @brief Create MaxCoeffType objects.
     */
    MinMaxCoeffFinder() {}

    /**
     * @brief Find the interval containing the Bernstein coefficients.
     *
     * @param b_coeffs is a list of symbolical Bernstein coefficients.
     * @return The pair minimum-maximum among all the Bernstein
     *          coefficients in `b_coeffs`.
     */
    virtual std::pair<double, double> find_coeffs_itvl(
        const std::vector<SymbolicAlgebra::Expression<>> &b_coeffs) const;
  };

  /**
   * @brief  A class to find the minimum and the maximum parametric Bernstein
   * coefficients
   */
  class ParamMinMaxCoeffFinder : public MinMaxCoeffFinder
  {
    const std::vector<SymbolicAlgebra::Symbol<>> &params;
    const Polytope &paraSet;
    /**
     * @brief Evaluate the parametric Bernstein coefficient upper-bound
     *
     * @param bernCoeff is the symbolical representation of Bernstein
     *                  coefficient
     * @return The numerical evaluation of parametric Bernstein
     *         coefficient upper-bound
     */
    double
    maximize_coeff(const SymbolicAlgebra::Expression<> &bernCoeff) const;

    /**
     * @brief Evaluate the parametric Bernstein coefficient lower-bound
     *
     * @param bernCoeff is the symbolical representation of Bernstein
     *                  coefficient
     * @return The numerical evaluation of parametric Bernstein coefficient
     *         lower-bound
     */
    double
    minimize_coeff(const SymbolicAlgebra::Expression<> &bernCoeff) const;

  public:
    /**
     * @brief Constructor
     *
     * @param params is the list of parameters
     * @param paraSet is the set of admissible values for parameters
     */
    ParamMinMaxCoeffFinder(
        const std::vector<SymbolicAlgebra::Symbol<>> &params,
        const Polytope &paraSet):
        MinMaxCoeffFinder(),
        params(params), paraSet(paraSet)
    {
    }

    /**
     * @brief Find the interval containing the Bernstein coefficients
     *
     * @param b_coeffs is a list of symbolical Bernstein coefficients
     * @return The pair minimum-maximum among all the Bernstein
     *          coefficients in `b_coeffs`
     */
    std::pair<double, double> find_coeffs_itvl(
        const std::vector<SymbolicAlgebra::Expression<>> &b_coeffs) const;
  };

  /**
   * @brief Transform the bundle according to a dynamic law
   *
   * @param[in] variables is the vector of variables
   * @param[in] dynamics is the vector of dynamic law expressions
   * @param[in] max_finder is a pointer to a MinMaxCoeffFinder object
   * @param[in] mode transformation mode, i.e., OFO or AFO
   * @returns the transformed bundle
   */
  Bundle transform(const std::vector<SymbolicAlgebra::Symbol<>> &variables,
                   const std::vector<SymbolicAlgebra::Expression<>> &dynamics,
                   const MinMaxCoeffFinder *max_finder,
                   Bundle::transfomation_mode mode) const;

  /**
   * @brief Add to the bundle templates for some directions
   *
   * @param missing_template_dirs the indices of the directions whose
   *                              template we are missing
   */
  void add_templates_for(std::set<unsigned int> &missing_template_dirs);

public:
  /**
   * @brief A copy constructor
   *
   * @param orig is the model for the new object
   */
  Bundle(const Bundle &orig);

  /**
   * @brief A move constructor
   *
   * @param orig is the model for the new object
   */
  Bundle(Bundle &&orig);

  /**
   * @brief A constructor
   *
   * @param[in] directions is the direction vector
   * @param[in] lower_bounds is the vector of direction lower bounds
   * @param[in] upper_bounds is the vector of direction upper bounds
   * @param[in] templates is the template vector
   */
  Bundle(const std::vector<LinearAlgebra::Vector<double>> &directions,
         const LinearAlgebra::Vector<double> &lower_bounds,
         const LinearAlgebra::Vector<double> &upper_bounds,
         const std::vector<LinearAlgebra::Vector<int>> &templates);

  /**
   * @brief A move constructor
   *
   * @param[in] directions is the direction vector
   * @param[in] lower_bounds is the vector of direction lower bounds
   * @param[in] upper_bounds is the vector of direction upper bounds
   * @param[in] templates is the template vector
   */
  Bundle(std::vector<LinearAlgebra::Vector<double>> &&directions,
         LinearAlgebra::Vector<double> &&lower_bounds, 
         LinearAlgebra::Vector<double> &&upper_bounds,
         std::vector<LinearAlgebra::Vector<int>> &&templates);

  /**
   * @brief Get the intersection between two bundles
   *
   * This method intersects the current object and another bundle.
   * The result is stored in the current object and a reference
   * to it is returned.
   *
   * @param A is the intersecting bundle
   * @return a reference to the updated object
   */
  Bundle &intersect_with(const Bundle &A);

  /**
   * @brief Get the intersection between this bundle and a linear set
   *
   * This method intersects the current instance of the `Bundle` class
   * and a possible unbounded linear set provided as a `LinearSystem`.
   * The result is stored in the current object and a reference to it
   * is returned.
   *
   * @param ls is the intersecting linear system
   * @return a reference to the updated object
   */
  Bundle &intersect_with(const LinearSystem &ls);

  /**
   * @brief Get the dimension of the bundle space
   *
   * @return the number of dimensions of the bundle space
   */
  unsigned int dim() const
  {
    return (directions.size() == 0 ? 0 : directions.front().size());
  }

  /**
   * @brief Get the number of templates in the bundle
   *
   * @return the number of templates in the bundle
   */
  unsigned int num_of_templates() const
  {
    return templates.size();
  }

  /**
   * @brief Get the number of directions in the bundle
   *
   * @return the number of directions in the bundle
   */
  unsigned int size() const
  {
    return this->directions.size();
  }

  /**
   * @brief Get the template vector
   *
   * @return a reference to the template vector
   */
  const std::vector<LinearAlgebra::Vector<int>> &get_templates() const
  {
    return this->templates;
  }

  /**
   * @brief Get the i-th templatein the bundle
   *
   * @param i is the index of the aimed template
   * @return a reference to the i-th template
   */
  const LinearAlgebra::Vector<int> &get_template(unsigned int i) const
  {
    return this->templates[i];
  }

  /**
   * @brief Get the vector of directions
   *
   * @return a reference to the vector of directions
   */
  const std::vector<LinearAlgebra::Vector<double>> &get_directions() const
  {
    return this->directions;
  }

  /**
   * @brief Get the i-th direction in the bundle
   *
   * @param i is the index of the aimed direction
   * @return  a reference to the i-th direction in the bundle
   */
  const LinearAlgebra::Vector<double> &get_direction(unsigned int i) const
  {
    return this->directions[i];
  }

  /**
   * @brief Get the upper bound of the i-th direction
   *
   * @param i is the index of the direction whose
   *        upper bound must be returned
   * @return a reference to the upper bound of the i-th direction
   */
  const double &get_upper_bound(unsigned int i) const
  {
    return this->upper_bounds[i];
  }

  /**
   * @brief Get the lower bound of the i-th direction
   *
   * @param i is the index of the direction whose
   *        lower bound must be returned
   * @return a reference to the lower bound of the i-th direction
   */
  const double &get_lower_bound(unsigned int i) const
  {
    return this->lower_bounds[i];
  }

  /**
   * @brief Get the vector of the direction upper bounds
   *
   * @return a reference to the direction upper bounds
   */
  const std::vector<double> &get_upper_bounds() const
  {
    return this->upper_bounds;
  }

  /**
   * @brief Get the vector of the direction lower bounds
   *
   * @return a reference to the direction lower bounds
   */
  const std::vector<double> &get_lower_bounds() const
  {
    return this->lower_bounds;
  }

  /**
   * @brief Generate the polytope represented by the bundle
   *
   * @returns polytope represented by the bundle
   */
  operator Polytope() const;

  /**
   * @brief Get the i-th parallelotope of the bundle
   *
   * @param[in] i is the index of the aimed parallelotope
   * @returns i-th parallelotope of the bundle
   */
  Parallelotope get_parallelotope(unsigned int i) const;

  /**
   * @brief Set the bundle templates
   *
   * @param[in] templates is the new template vector
   */
  void set_templates(const std::vector<LinearAlgebra::Vector<int>> &templates)
  {
    this->templates = templates;
  }

  /**
   * @brief Set the direction upper bounds
   *
   * @param upper_bounds is the vector of the direction upper bounds
   */
  void set_upper_bounds(const LinearAlgebra::Vector<double> &upper_bounds)
  {
    this->upper_bounds = upper_bounds;
  }

  /**
   * @brief Set the direction lower bounds
   *
   * @param lower_bounds is the vector of the direction lower bounds
   */
  void set_lower_bounds(const LinearAlgebra::Vector<double> &lower_bounds)
  {
    this->lower_bounds = lower_bounds;
  }

  /**
   * @brief Get a canonical bundle equivalent to the current one
   *
   * Get the canonical form for this bundle by minimizing the
   * difference between lower and upper bounds over all the
   * directions.
   *
   * @returns canonized bundle
   */
  Bundle get_canonical() const;

  /**
   * @brief Canonize the current bundle
   *
   * Turn this bundle in canonical form by minimizing the
   * difference between lower and upper bounds over all the
   * directions.
   *
   * @returns a reference to the canonized bundle
   */
  Bundle &canonize();

  /**
   * @brief Split the bundle in smaller sub-bundles
   * 
   * This method splits the bundles whose maximal magnitude, 
   * the maximal lenght of its generators, is greater than 
   * `max_magnitude` into a list of sub-bundles whose 
   * maximal magnitude is 
   * \f$\textrm{max_magnitude}*\textrm{split_ratio}\f$.
   *
   * If \f$m\f$ in the maximal magnitude of the input bundle 
   * and the bundle itself has \f$d\f$ directions,
   * this method may produces upto
   * \f$(\frac{m}{\textrm{max_magnitude}})^d\f$ 
   * sub-bundles.
   * 
   * @param[in] max_magnitude is the maximal magnitude that
   *                          triggers the split
   * @param[in] split_ratio is the ratio between maximal
   *                        magnitude of the output bundles and
   *                        that that triggers the split
   * @return a list of sub-bundles whose union covers the 
   *         input bundle
   */
  std::list<Bundle> split(const double max_magnitude,
                          const double split_ratio
                          = SPLIT_MAGNITUDE_RATIO) const;

  /**
   * Decompose the current symbolic polytope
   *
   * @param[in] alpha weight parameter in [0,1] for decomposition (0 for
   * distance, 1 for orthogonality)
   * @param[in] max_iters maximum number of randomly generated templates
   * @todo method does not work; it should be reviewed and fixed
   * @returns new bundle decomposing current symbolic polytope
   */
  Bundle decompose(double alpha, int max_iters);

  /**
   * @brief Transform the bundle according to a dynamic law
   *
   * @param[in] variables is the vector of variables
   * @param[in] dynamics is the vector of dynamic law expressions
   * @param[in] mode transformation mode, i.e., OFO or AFO
   * @returns the transformed bundle
   */
  inline Bundle
  transform(const std::vector<SymbolicAlgebra::Symbol<>> &variables,
            const std::vector<SymbolicAlgebra::Expression<>> &dynamics,
            transfomation_mode mode) const
  {
    MinMaxCoeffFinder max_finder;

    return transform(variables, dynamics, &max_finder, mode);
  }

  /**
   * @brief Transform the bundle according to a dynamic law
   *
   * @param[in] variables is the vector of variables
   * @param[in] parameters is the vector of parameter
   * @param[in] dynamics is the vector of dynamic law expressions
   * @param[in] parameter_set is the parameter set
   * @param[in] mode transformation mode, i.e., OFO or AFO
   * @returns the transformed bundle
   */
  inline Bundle
  transform(const std::vector<SymbolicAlgebra::Symbol<>> &variables,
            const std::vector<SymbolicAlgebra::Symbol<>> &parameters,
            const std::vector<SymbolicAlgebra::Expression<>> &dynamics,
            const Polytope &parameter_set, transfomation_mode mode) const
  {
    ParamMinMaxCoeffFinder max_finder(parameters, parameter_set);

    return transform(variables, dynamics, &max_finder, mode);
  }

  /**
   * @brief Perform the parametric synthesis for an atom
   *
   * This method computes a set of parameters such that the
   * transformation from the current bundle satisfies the
   * provided atom.
   *
   * @param[in] variables is the vector of variables
   * @param[in] parameters is the vector of parameter
   * @param[in] dynamics is the vector of dynamic law expressions
   * @param[in] parameter_set is the initial parameter set
   * @param[in] atom is the specification to be satisfied
   * @return a subset of `pSet` such that the transformation
   *         from the current bundle through the dynamic laws
   *         when the parameters are in the returned set satisfies
   *         the atomic formula
   */
  PolytopesUnion
  synthesize(const std::vector<SymbolicAlgebra::Symbol<>> &variables,
             const std::vector<SymbolicAlgebra::Symbol<>> &parameters,
             const std::vector<SymbolicAlgebra::Expression<>> &dynamics,
             const PolytopesUnion &parameter_set,
             const std::shared_ptr<STL::Atom> atom) const;

  /**
   * @brief Copy operator
   *
   * @return a reference to the updated object
   */
  Bundle &operator=(Bundle &&);

  virtual ~Bundle();

  friend void swap(Bundle &A, Bundle &B);
};

/**
 * @brief Swap the content of two bundles
 *
 * @param A is the first bundle to be swapped
 * @param B is the second bundle to be swapped
 */
void swap(Bundle &A, Bundle &B);

/**
 * @brief Get the intersection between two bundles
 * 
 * @param b1 is a bundle
 * @param b2 is a bundle
 * @return A bundle representing the intersection between 
 *         the bundles `b1` and `b2`
 */
inline Bundle intersect(const Bundle& b1, const Bundle& b2)
{
  Bundle res(b1);

  res.intersect_with(b2);

  return res;
}

#endif /* BUNDLE_H_ */
