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
  unsigned int dim;                       //!< dimension
  Matrix dir_matrix;                      //!< direction matrix
  Vector offp;                            //!< superior offset
  Vector offm;                            //!< inferior offset
  std::vector<std::vector<int>> t_matrix; //!< templates matrix
  Matrix Theta;                           //!< matrix of orthogonal proximity
  std::vector<GiNaC::lst> vars; //!< variables appearing in generato function
                                //!< vars[0] q: base vertex
                                //!< vars[1] alpha : free variables \in [0,1]
                                //!< vars[2] beta : generator amplitudes

  std::map<std::vector<int>, GiNaC::lst>
      bernCoeffs; //!< Bernstein coefficients map

  /**
   * Compute the distances between the half-spaced of the parallelotopes
   *
   * @returns vector of distances
   */
  Vector offsetDistances();

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
   * @param[in] vars list of variables for parallelotope generator functions
   * @param[in] dir_matrix matrix of directions
   * @param[in] offp upper offsets
   * @param[in] offm lower offsets
   * @param[in] t_matrix t_matrixs matrix
   */
  Bundle(const std::vector<GiNaC::lst> &vars, const Matrix &dir_matrix,
         const Vector &offp, const Vector &offm,
         const std::vector<std::vector<int>> &t_matrix);

  unsigned int getDim() const
  {
    return this->dim;
  }

  unsigned int getSize() const
  {
    return dir_matrix.size();
  }

  unsigned int getCard() const
  {
    return t_matrix.size();
  }

  unsigned int getNumDirs() const
  {
    return this->dir_matrix.size();
  }

  const std::vector<int> &getTemplate(long unsigned int i) const
  {
    return this->t_matrix[i];
  }

  const std::vector<std::vector<double>> &getDirectionMatrix() const
  {
    return this->dir_matrix;
  }

  const double &getOffp(long unsigned int i) const
  {
    return this->offp[i];
  }

  const double &getOffm(long unsigned int i) const
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
  Bundle transform(const GiNaC::lst &vars, const GiNaC::lst &f,
                   ControlPointStorage &controlPts, int mode) const;

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
  Bundle transform(const GiNaC::lst &vars, const GiNaC::lst &params,
                   const GiNaC::lst &f, const Polytope &paraSet,
                   ControlPointStorage &controlPts, int mode) const;

  Bundle &operator=(Bundle &&);

  virtual ~Bundle();

  friend void swap(Bundle &A, Bundle &B);
};

void swap(Bundle &A, Bundle &B);

#endif /* BUNDLE_H_ */
