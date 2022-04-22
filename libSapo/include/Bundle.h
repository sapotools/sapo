/**
 * @file Bundle.h
 * Represent and manipulate bundles of parallelotopes whose intersection
 * represents a polytope
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @version 0.1
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

// define the default versor magnitude multiplier for bundle splits
#define SPLIT_MAGNITUDE_RATIO 0.75

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
  typedef enum {
    AFO, /* the image of any parallelotope in the bundle will
          * be over-approximated by using all the templates
          * of the bundle itself */
    OFO  /* the image of any parallelotope in the bundle will
          * be over-approximated by using exclusively the
          * parallelotope own template */
  } transfomation_mode;

private:
  std::vector<Vector<double>> directions; //!< the vector of directions
  Vector<double> lower_bounds;            //!< direction upper bounds
  Vector<double> upper_bounds;            //!< direction lower bounds
  std::vector<Vector<int>> templates;     //!< templates matrix

  // constraints over directions (assertions)
  // constrainedDirection[i] * vars <= constraintOffset
  std::vector<Vector<double>> constraintDirections;
  Vector<double> constraintOffsets;

  /**
   * Compute the edge lengths
   *
   * This method compute the edge lengths, i.e., the distances between the
   * half-spaces bounding the bundle.
   *
   * @returns the vector of the edge lengths
   */
  Vector<double> edge_lengths();

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
   *            template we are missing
   */
  void add_templates_for(std::set<unsigned int> &missing_template_dirs);

public:
  /**
   * @brief A copy constructor
   *
   * @param orig is the model for the new object.
   */
  Bundle(const Bundle &orig);

  /**
   * @brief A move constructor
   *
   * @param orig is the model for the new object.
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
  Bundle(const std::vector<Vector<double>> &directions,
         const Vector<double> &lower_bounds,
         const Vector<double> &upper_bounds,
         const std::vector<Vector<int>> &templates);

  /**
   * @brief A constructor
   *
   * @param[in] directions is the direction vector
   * @param[in] lower_bounds is the vector of direction lower bounds
   * @param[in] upper_bounds is the vector of direction upper bounds
   * @param[in] templates is the template vector
   * @param[in] constrDirs directions that are constrained by assumptions
   * @param[in] constrOffsets offsets of assumptions
   */
  Bundle(const std::vector<Vector<double>> &directions,
         const Vector<double> &lower_bounds,
         const Vector<double> &upper_bounds,
         const std::vector<Vector<int>> &templates,
         const std::vector<Vector<double>> &constrDirs,
         const Vector<double> &constrOffsets);

  /**
   * @brief A move constructor
   *
   * @param[in] directions is the direction vector
   * @param[in] lower_bounds is the vector of direction lower bounds
   * @param[in] upper_bounds is the vector of direction upper bounds
   * @param[in] templates is the template vector
   */
  Bundle(std::vector<Vector<double>> &&directions,
         Vector<double> &&lower_bounds, Vector<double> &&upper_bounds,
         std::vector<Vector<int>> &&templates);

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
  Bundle &intersect(const Bundle &A);

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
  Bundle &intersect(const LinearSystem &ls);

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
  const std::vector<Vector<int>> &get_templates() const
  {
    return this->templates;
  }

  /**
   * @brief Get the i-th templatein the bundle
   *
   * @param i is the index of the aimed template
   * @return a reference to the i-th template
   */
  const Vector<int> &get_template(unsigned int i) const
  {
    return this->templates[i];
  }

  /**
   * @brief Get the vector of directions
   *
   * @return a reference to the vector of directions
   */
  const std::vector<Vector<double>> &get_directions() const
  {
    return this->directions;
  }

  /**
   * @brief Get the i-th direction in the bundle
   *
   * @param i is the index of the aimed direction
   * @return  a reference to the i-th direction in the bundle
   */
  const Vector<double> &get_direction(unsigned int i) const
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
   * @returns i-th parallelotope
   */
  Parallelotope get_parallelotope(unsigned int i) const;

  /**
   * @brief Set the bundle templates
   *
   * @param[in] templates is the new template vector
   */
  void set_templates(const std::vector<Vector<int>> &templates)
  {
    this->templates = templates;
  }

  /**
   * @brief Set the direction upper bounds
   *
   * @param upper_bounds is the vector of the direction upper bounds
   */
  void set_upper_bounds(const Vector<double> &upper_bounds)
  {
    this->upper_bounds = upper_bounds;
  }

  /**
   * @brief Set the direction lower bounds
   *
   * @param upper_bounds is the vector of the direction lower bounds
   */
  void set_lower_bounds(const Vector<double> &lower_bounds)
  {
    this->lower_bounds = lower_bounds;
  }

  /**
   * @brief Canonize the current bundle
   *
   * Get the canonical form for this bundle by minimizing the
   * difference between lower and upper bounds over all the
   * directions.
   *
   * @returns canonized bundle
   */
  Bundle get_canonical() const;

  /**
   * @brief Split a bundle in a list of smaller bundles
   *
   * @param max_bundle_magnitude is the maximal edge length of the
   *        resulting bundles
   * @param split_magnitude_ratio is the ratio of the `max_bundle_magnitude`
   *        that is used a maximal magnitude of the bundles in output
   * @return A list of bundles whose maximal versor magnitude is
   *        `split_magnitude_ratio`*`max_bundle_magnitude` and whose
   *        union is the current bundle
   */
  std::list<Bundle> split(const double max_bundle_magnitude,
                          const double split_magnitude_ratio
                          = SPLIT_MAGNITUDE_RATIO) const;

  /**
   * @brief Decompose the current symbolic polytope
   *
   * @param[in] alpha weight parameter in [0,1] for decomposition (0 for
   * distance, 1 for orthogonality)
   * @param[in] max_iter maximum number of randomly generated templatess
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
   * transformation from the current bundle satisfy the
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
             const std::shared_ptr<Atom> atom) const;

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

#endif /* BUNDLE_H_ */
