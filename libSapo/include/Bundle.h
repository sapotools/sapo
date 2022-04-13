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
#include "Parallelotope.h"
#include "VarsGenerator.h"

// define the default versor magnitude multiplier for bundle splits
#define SPLIT_MAGNITUDE_RATIO 0.75

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
  DenseLinearAlgebra::Matrix<double> dir_matrix; //!< direction matrix
  Vector<double> upper_offset;                   //!< superior offset
  Vector<double> lower_offset;                   //!< inferior offset
  DenseLinearAlgebra::Matrix<int> t_matrix;      //!< templates matrix

  // constraints over directions (assertions)
  // constrainedDirection[i] * vars <= constraintOffset
  DenseLinearAlgebra::Matrix<double> constraintDirections;
  Vector<double> constraintOffsets;

  /**
   * Compute the distances between the half-spaced of the parallelotopes
   *
   * @returns vector of distances
   */
  Vector<double> offsetDistances();

  /**
   * @brief A finder for Bernstein coefficient upper and lower-bounds.
   */
  class MaxCoeffFinder
  {
    /**
     * @brief Evaluate the Bernstein coefficient upper-bound.
     *
     * @param bernCoeff is the symbolical representation of Bernstein
     * coefficient.
     * @return The numerical evaluation of Bernstein coefficient upper-bound.
     */
    virtual double
    coeff_eval_p(const SymbolicAlgebra::Expression<> &bernCoeff) const;

    /**
     * @brief Evaluate the Bernstein coefficient lower-bound complementary.
     *
     * @param bernCoeff is the symbolical representation of Bernstein
     *                  coefficient.
     * @return The numerical evaluation of Bernstein coefficient
     *         lower-bound complementary.
     */
    virtual double
    coeff_eval_m(const SymbolicAlgebra::Expression<> &bernCoeff) const;

  public:
    typedef struct MaxCoeffType {
      double p; //<! the Bernstein coefficient upper-bound
      double m; //<! the Bernstein coefficient lower-bound complementary
    } MaxCoeffType;

    /**
     * @brief Create MaxCoeffType objects.
     */
    MaxCoeffFinder() {}

    /**
     * @brief Find the maximum of lower-bound complementary and
     *        upper-bound for Bernstein coefficients.
     *
     * @param b_coeffs is a list of symbolical Bernstein coefficients.
     * @return The maximum of both lower-bound complementary and
     *         upper-bound for all the Bernstein coefficients in
     *         `b_coeffs`.
     */
    MaxCoeffType find_max_coeffs(
        const std::vector<SymbolicAlgebra::Expression<>> &b_coeffs) const;
  };

  /**
   * @brief A finder for parametric Bernstein coefficient upper and
   * lower-bounds.
   */
  class ParamMaxCoeffFinder : public MaxCoeffFinder
  {
    const std::vector<SymbolicAlgebra::Symbol<>> &params;
    const Polytope &paraSet;
    /**
     * @brief Evaluate the parametric Bernstein coefficient upper-bound.
     *
     * @param bernCoeff is the symbolical representation of Bernstein
     * coefficient.
     * @return The numerical evaluation of parametric Bernstein
     *         coefficient upper-bound.
     */
    double coeff_eval_p(const SymbolicAlgebra::Expression<> &bernCoeff) const;

    /**
     * @brief Evaluate the parametric Bernstein coefficient lower-bound
     *        complementary.
     *
     * @param bernCoeff is the symbolical representation of Bernstein
     *                  coefficient.
     * @return The numerical evaluation of parametric Bernstein coefficient
     *         lower-bound complementary.
     */
    double coeff_eval_m(const SymbolicAlgebra::Expression<> &bernCoeff) const;

  public:
    /**
     * @brief Create ParamMaxCoeffFinder objects.
     *
     * @param params is the list of parameters.
     * @param paraSet is the set of admissible values for parameters.
     */
    ParamMaxCoeffFinder(const std::vector<SymbolicAlgebra::Symbol<>> &params,
                        const Polytope &paraSet):
        MaxCoeffFinder(),
        params(params), paraSet(paraSet)
    {
    }
  };

  /**
   * Transform the bundle
   *
   * @param[in] vars variables appearing in the transforming function
   * @param[in] f transforming function
   * @param[in] max_finder is a pointer to a MaxCoeffFinder object.
   * @param[in] mode transformation mode, i.e., OFO or AFO.
   * @returns transformed bundle
   */
  Bundle transform(const std::vector<SymbolicAlgebra::Symbol<>> &vars,
                   const std::vector<SymbolicAlgebra::Expression<>> &f,
                   const MaxCoeffFinder *max_finder,
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
   * @brief Copy constructor.
   *
   * @param orig is the model for the new object.
   */
  Bundle(const Bundle &orig);

  /**
   * @brief Swap constructor.
   *
   * @param orig is the model for the new object.
   */
  Bundle(Bundle &&orig);

  /**
   * Constructor
   *
   * @param[in] dir_matrix matrix of directions
   * @param[in] upper_offset upper offsets
   * @param[in] lower_offset lower offsets
   * @param[in] t_matrix t_matrixs matrix
   */
  Bundle(const DenseLinearAlgebra::Matrix<double> &dir_matrix,
         const Vector<double> &upper_offset,
         const Vector<double> &lower_offset,
         const DenseLinearAlgebra::Matrix<int> &t_matrix);

  /**
   * Constructor
   *
   * @param[in] dir_matrix matrix of directions
   * @param[in] upper_offset upper offsets
   * @param[in] lower_offset lower offsets
   * @param[in] t_matrix t_matrixs matrix
   * @param[in] constrDirs directions that are constrained by assumptions
   * @param[in] constrOffsets offsets of assumptions
   */
  Bundle(const DenseLinearAlgebra::Matrix<double> &dir_matrix,
         const Vector<double> &upper_offset,
         const Vector<double> &lower_offset,
         const DenseLinearAlgebra::Matrix<int> &t_matrix,
         const DenseLinearAlgebra::Matrix<double> &constrDirs,
         const Vector<double> &constrOffsets);

  /**
   * Move constructor
   *
   * @param[in] dir_matrix matrix of directions
   * @param[in] upper_offset upper offsets
   * @param[in] lower_offset lower offsets
   * @param[in] t_matrix t_matrixs matrix
   */
  Bundle(DenseLinearAlgebra::Matrix<double> &&dir_matrix,
         Vector<double> &&upper_offset, Vector<double> &&lower_offset,
         DenseLinearAlgebra::Matrix<int> &&t_matrix);

  /**
   * @brief Get the intersection between two bundles
   *
   * This method intersects the current object and another bundle.
   * The result is stored in the current object and a reference
   * to it is returned.
   *
   * @param A is the intersecting bundle.
   * @return a reference to the updated object.
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
   * @param ls is the intersecting bundle.
   * @return a reference to the updated object.
   */
  Bundle &intersect(const LinearSystem &ls);

  /**
   * @brief Return the dimension of the bundle space
   *
   * @return the number of dimensions of the bundle space
   */
  unsigned int dim() const
  {
    return (dir_matrix.size() == 0 ? 0 : dir_matrix.front().size());
  }

  unsigned int num_of_templates() const
  {
    return t_matrix.size();
  }

  unsigned int size() const
  {
    return this->dir_matrix.size();
  }

  const DenseLinearAlgebra::Matrix<int> &get_templates() const
  {
    return this->t_matrix;
  }

  const Vector<int> &get_templates(long unsigned int i) const
  {
    return this->t_matrix[i];
  }

  const DenseLinearAlgebra::Matrix<double> &get_directions() const
  {
    return this->dir_matrix;
  }

  const double &get_offsetp(long unsigned int i) const
  {
    return this->upper_offset[i];
  }

  const double &get_offsetm(long unsigned int i) const
  {
    return this->lower_offset[i];
  }

  const std::vector<double> &get_upper_offset() const
  {
    return this->upper_offset;
  }

  const std::vector<double> &get_lower_offset() const
  {
    return this->lower_offset;
  }

  /**
   * Generate the polytope represented by the bundle
   *
   * @returns polytope represented by the bundle
   */
  operator Polytope() const;

  /**
   * Get the i-th parallelotope of the bundle
   *
   * @param[in] i parallelotope index to fetch
   * @returns i-th parallelotope
   */
  Parallelotope getParallelotope(unsigned int i) const;

  /**
   * Set the bundle t_matrix
   *
   * @param[in] t_matrix new t_matrix
   */
  void setTemplate(const DenseLinearAlgebra::Matrix<int> &t_matrix)
  {
    this->t_matrix = t_matrix;
  }

  void setOffsetP(const Vector<double> &upper_offset)
  {
    this->upper_offset = upper_offset;
  }

  void setOffsetM(const Vector<double> &lower_offset)
  {
    this->lower_offset = lower_offset;
  }

  // operations on bundles
  /**
   * Canonize the current bundle pushing the constraints toward the symbolic
   * polytope
   *
   * @returns canonized bundle
   */
  Bundle get_canonical() const;

  /**
   * @brief Split a bundle in a list of smaller bundles.
   *
   * @param max_bundle_magnitude is the maximal versor magnitude of the
   * resulting bundles.
   * @param split_magnitude_ratio is the ratio of the `max_bundle_magnitude`
   * that is used a maximal magnitude of the bundles in output.
   * @return A list of bundles whose maximal versor magnitude is
   * `split_magnitude_ratio`*`max_bundle_magnitude` and whose union is
   * the current bundle.
   */
  std::list<Bundle> split(const double max_bundle_magnitude,
                          const double split_magnitude_ratio
                          = SPLIT_MAGNITUDE_RATIO) const;

  /**
   * Decompose the current symbolic polytope
   *
   * @param[in] alpha weight parameter in [0,1] for decomposition (0 for
   * distance, 1 for orthogonality)
   * @param[in] max_iter maximum number of randomly generated t_matrixs
   * @returns new bundle decomposing current symbolic polytope
   */
  Bundle decompose(double alpha, int max_iters);

  /**
   * Transform the bundle
   *
   * @param[in] vars variables appearing in the transforming function
   * @param[in] f transforming function
   * @param[in] mode transformation mode, i.e., OFO or AFO
   * @returns transformed bundle
   */
  inline Bundle transform(const std::vector<SymbolicAlgebra::Symbol<>> &vars,
                          const std::vector<SymbolicAlgebra::Expression<>> &f,
                          transfomation_mode mode) const
  {
    MaxCoeffFinder max_finder;

    return transform(vars, f, &max_finder, mode);
  }

  /**
   * Parametric transformation of the bundle
   *
   * @param[in] vars variables appearing in the transforming function
   * @param[in] params parameters appearing in the transforming function
   * @param[in] f transforming function
   * @param[in] paraSet set of parameters
   * @param[in] mode transformation mode, i.e., OFO or AFO
   * @returns transformed bundle
   */
  inline Bundle transform(const std::vector<SymbolicAlgebra::Symbol<>> &vars,
                          const std::vector<SymbolicAlgebra::Symbol<>> &params,
                          const std::vector<SymbolicAlgebra::Expression<>> &f,
                          const Polytope &paraSet,
                          transfomation_mode mode) const
  {
    ParamMaxCoeffFinder max_finder(params, paraSet);

    return transform(vars, f, &max_finder, mode);
  }

  Bundle &operator=(Bundle &&);

  virtual ~Bundle();

  friend void swap(Bundle &A, Bundle &B);
};
void swap(Bundle &A, Bundle &B);

std::vector<SymbolicAlgebra::Expression<>> build_instanciated_generator_functs(
    const std::vector<SymbolicAlgebra::Symbol<>> &alpha,
    const Parallelotope &P);

std::vector<SymbolicAlgebra::Expression<>>
sub_vars(const std::vector<SymbolicAlgebra::Expression<>> &ex_list,
         const std::vector<SymbolicAlgebra::Symbol<>> &vars,
         const std::vector<SymbolicAlgebra::Expression<>> &expressions);

#endif /* BUNDLE_H_ */