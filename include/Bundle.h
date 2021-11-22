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
  unsigned int dim;                // dimension
  Matrix L;                        // direction matrix
  Vector offp;                     // superior offset
  Vector offm;                     // inferior offset
  std::vector<std::vector<int>> T; // templates matrix
  Matrix Theta;                    // matrix of orthogonal proximity
  std::vector<GiNaC::lst> vars;    // variables appearing in generato function
                                   // vars[0] q: base vertex
                                   // vars[1] alpha : free variables \in [0,1]
                                   // vars[2] beta : generator amplitudes

  // map with Bernstein coefficients
  std::map<std::vector<int>, GiNaC::lst> bernCoeffs;

  double initTheta();
  Vector offsetDistances();

  // operations on vectors
  double norm(std::vector<double> v);
  double prod(std::vector<double> v1, std::vector<double> v2);
  double angle(std::vector<double> v1, std::vector<double> v2);
  double orthProx(std::vector<double> v1, std::vector<double> v2);
  double maxOrthProx(int vIdx, std::vector<int> dirsIdx);
  double maxOrthProx(std::vector<int> dirsIdx);
  double maxOrthProx(std::vector<std::vector<int>> T);
  double maxOffsetDist(int vIdx, std::vector<int> dirsIdx,
                       std::vector<double> dists);
  double maxOffsetDist(std::vector<int> dirsIdx, std::vector<double> dists);
  double maxOffsetDist(std::vector<std::vector<int>> T,
                       std::vector<double> dists);

  bool validTemp(std::vector<std::vector<int>> T, unsigned int card,
                 std::vector<int> dirs); // check if a template is valid
  std::vector<GiNaC::lst> transformContrPts(GiNaC::lst vars, GiNaC::lst f,
                                            int mode);

public:
  // constructors
  Bundle(const Bundle &orig);
  Bundle(Bundle &&orig);
  Bundle(const Matrix &L, const Vector &offp, const Vector &offm,
         const std::vector<std::vector<int>> &T);
  Bundle(const std::vector<GiNaC::lst> &vars, const Matrix &L,
         const Vector &offp, const Vector &offm,
         const std::vector<std::vector<int>> &T);

  unsigned int getDim() const
  {
    return this->dim;
  }

  unsigned int getSize() const
  {
    return L.size();
  }

  unsigned int getCard() const
  {
    return T.size();
  }

  unsigned int getNumDirs() const
  {
    return this->L.size();
  }

  const std::vector<int> &getTemplate(long unsigned int i) const
  {
    return this->T[i];
  }

  const std::vector<std::vector<double>> &getDirectionMatrix() const
  {
    return this->L;
  }

  const double &getOffp(long unsigned int i) const
  {
    return this->offp[i];
  }

  const double &getOffm(long unsigned int i) const
  {
    return this->offm[i];
  }

  operator Polytope() const;
  Parallelotope getParallelotope(unsigned int i) const;

  void setTemplate(std::vector<std::vector<int>> T);
  void setOffsetP(Vector offp)
  {
    this->offp = offp;
  }

  void setOffsetM(Vector offm)
  {
    this->offm = offm;
  }

  // operations on bundles
  Bundle get_canonical() const;
  Bundle decompose(double alpha, int max_iters);
  Bundle transform(const GiNaC::lst &vars, const GiNaC::lst &f,
                   ControlPointStorage &controlPts, int mode) const;
  Bundle transform(const GiNaC::lst &vars, const GiNaC::lst &params,
                   const GiNaC::lst &f, const Polytope &paraSet,
                   ControlPointStorage &controlPts, int mode) const;

  Bundle &operator=(Bundle &&);

  virtual ~Bundle();

  friend void swap(Bundle &A, Bundle &B);
};

void swap(Bundle &A, Bundle &B);

#endif /* BUNDLE_H_ */
