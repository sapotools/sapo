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

#include "BaseConverter.h"
#include "Polytope.h"
#include "Parallelotope.h"
#include "VarsGenerator.h"
#include "ControlPointStorage.h"

class Bundle
{
  using Vector = std::vector<double>;
  using Matrix = std::vector<Vector>;
  using CtrlPointType
      = std::map<std::vector<int>, std::pair<GiNaC::lst, GiNaC::lst>>;

private:
  Matrix dir_matrix;                      //!< direction matrix
  Vector offp;                            //!< superior offset
  Vector offm;                            //!< inferior offset
  std::vector<std::vector<int>> t_matrix; //!< templates matrix
  Matrix Theta;                           //!< matrix of orthogonal proximity
  GiNaC::lst q;     //!< variables to represent base vertex
  GiNaC::lst alpha; //!< free variables \in [0,1]
  GiNaC::lst beta;  //!< generator amplitude variables
  std::map<std::vector<int>, GiNaC::lst>
      bernCoeffs; //!< Bernstein coefficients map

  /**
   * Compute the distances between the half-spaced of the parallelotopes
   *
   * @returns vector of distances
   */
  Vector offsetDistances();

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
    virtual double coeff_eval_p(const GiNaC::ex &bernCoeff) const;

    /**
     * @brief Evaluate the Bernstein coefficient lower-bound complementary.
     *
     * @param bernCoeff is the symbolical representation of Bernstein
     *                  coefficient.
     * @return The numerical evaluation of Bernstein coefficient
     *         lower-bound complementary.
     */
    virtual double coeff_eval_m(const GiNaC::ex &bernCoeff) const;

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
     * @param subs is a list of variable assignaments
     * @return The maximum of both lower-bound complementary and
     *         upper-bound for all the Bernstein coefficients in
     *         `b_coeffs`.
     */
    MaxCoeffType find_max_coeffs(const GiNaC::lst &b_coeffs,
                                 const GiNaC::lst &subs) const;
  };

  /**
   * @brief A finder for parametric Bernstein coefficient upper and
   * lower-bounds.
   */
  class ParamMaxCoeffFinder : public MaxCoeffFinder
  {
    const GiNaC::lst &params;
    const Polytope &paraSet;
    /**
     * @brief Evaluate the parametric Bernstein coefficient upper-bound.
     *
     * @param bernCoeff is the symbolical representation of Bernstein
     * coefficient.
     * @return The numerical evaluation of parametric Bernstein
     *         coefficient upper-bound.
     */
    double coeff_eval_p(const GiNaC::ex &bernCoeff) const;

    /**
     * @brief Evaluate the parametric Bernstein coefficient lower-bound
     *        complementary.
     *
     * @param bernCoeff is the symbolical representation of Bernstein
     *                  coefficient.
     * @return The numerical evaluation of parametric Bernstein coefficient
     *         lower-bound complementary.
     */
    double coeff_eval_m(const GiNaC::ex &bernCoeff) const;

  public:
    /**
     * @brief Create ParamMaxCoeffFinder objects.
     *
     * @param params is the list of parameters.
     * @param paraSet is the set of admissible values for parameters.
     */
    ParamMaxCoeffFinder(const GiNaC::lst &params, const Polytope &paraSet):
        MaxCoeffFinder(), params(params), paraSet(paraSet)
    {
    }
  };

  /**
   * Transform the bundle
   *
   * @param[in] vars variables appearing in the transforming function
   * @param[in] f transforming function
   * @param[in,out] controlPts control points computed so far that might be
   *                           updated
   * @param[in] max_finder is a pointer to a MaxCoeffFinder object.
   * @param[in] mode transformation mode (0=OFO,1=AFO)
   * @returns transformed bundle
   */
  Bundle transform(const GiNaC::lst &vars, const GiNaC::lst &f,
                   ControlPointStorage &controlPts,
                   const MaxCoeffFinder *max_finder, int mode) const;

public:
  // constructors
  Bundle(const Bundle &orig);
  Bundle(Bundle &&orig);

  /**
   * Constructor that instantiates the bundle with auto-generated variables
   *
   * @param[in] dir_matrix matrix of directions
   * @param[in] offp upper offsets
   * @param[in] offm lower offsets
   * @param[in] t_matrix t_matrixs matrix
   */
  Bundle(const Matrix &dir_matrix, const Vector &offp, const Vector &offm,
         const std::vector<std::vector<int>> &t_matrix);

  /**
   * Constructor that instantiates the bundle
   *
   * @param[in] q is a list of variables to represent base vertex
   * @param[in] alpha is a list of free variables
   * @param[in] beta are generator amplitude variables
   * @param[in] dir_matrix matrix of directions
   * @param[in] offp upper offsets
   * @param[in] offm lower offsets
   * @param[in] t_matrix t_matrixs matrix
   */
  Bundle(const GiNaC::lst &q, const GiNaC::lst &alpha, const GiNaC::lst &beta,
         const Matrix &dir_matrix, const Vector &offp, const Vector &offm,
         const std::vector<std::vector<int>> &t_matrix);

  unsigned int dim() const
  {
    return (dir_matrix.size() == 0 ? 0 : dir_matrix.front().size());
  }

  unsigned int size() const
  {
    return dir_matrix.size();
  }

  unsigned int num_of_templates() const
  {
    return t_matrix.size();
  }

  unsigned int num_of_dirs() const
  {
    return this->dir_matrix.size();
  }

  const std::vector<int> &get_templates(long unsigned int i) const
  {
    return this->t_matrix[i];
  }

  const std::vector<std::vector<double>> &get_directions() const
  {
    return this->dir_matrix;
  }

  /**
   * Get variables of base vertex
   *
   * @returns base vertex variables
   */
  const GiNaC::lst &get_q() const
  {
    return this->q;
  }

  /**
   * Get free variables
   *
   * @returns free variables
   */
  const GiNaC::lst &get_alpha() const
  {
    return this->alpha;
  }

  /**
   * Get variables of generator lengths
   *
   * @returns generator lengths variables
   */
  const GiNaC::lst &get_beta() const
  {
    return this->beta;
  }

  const double &get_offsetp(long unsigned int i) const
  {
    return this->offp[i];
  }

  const double &get_offsetm(long unsigned int i) const
  {
    return this->offm[i];
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
  void setTemplate(const std::vector<std::vector<int>> &t_matrix)
  {
    this->t_matrix = t_matrix;
  }

  void setOffsetP(Vector offp)
  {
    this->offp = offp;
  }

  void setOffsetM(Vector offm)
  {
    this->offm = offm;
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
   * @param[in,out] controlPts control points computed so far that might be
   * updated
   * @param[in] mode transformation mode (0=OFO,1=AFO)
   * @returns transformed bundle
   */
  inline Bundle transform(const GiNaC::lst &vars, const GiNaC::lst &f,
                          ControlPointStorage &controlPts, int mode) const
  {
    MaxCoeffFinder max_finder;

    return transform(vars, f, controlPts, &max_finder, mode);
  }

  /**
   * Parametric transformation of the bundle
   *
   * @param[in] vars variables appearing in the transforming function
   * @param[in] params parameters appearing in the transforming function
   * @param[in] f transforming function
   * @param[in] paraSet set of parameters
   * @param[in,out] controlPts control points computed so far that might be
   * updated
   * @param[in] mode transformation mode (0=OFO,1=AFO)
   * @returns transformed bundle
   */
  inline Bundle transform(const GiNaC::lst &vars, const GiNaC::lst &params,
                          const GiNaC::lst &f, const Polytope &paraSet,
                          ControlPointStorage &controlPts, int mode) const
  {
    ParamMaxCoeffFinder max_finder(params, paraSet);

    return transform(vars, f, controlPts, &max_finder, mode);
  }

  Bundle &operator=(Bundle &&);

  virtual ~Bundle();

  friend void swap(Bundle &A, Bundle &B);
};

void swap(Bundle &A, Bundle &B);

GiNaC::lst build_generator_functs(const GiNaC::lst &q, const GiNaC::lst &alpha,
                                  const GiNaC::lst &beta,
                                  const Parallelotope &P);

#endif /* BUNDLE_H_ */
