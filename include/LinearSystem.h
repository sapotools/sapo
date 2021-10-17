/**
 * @file LinearSystem.h
 * Represent and manipulate a linear system
 * It can be used to represent polytopes (reached states, parameters, etc.)
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @version 0.1
 */

#ifndef LINEARSYSTEM_H_
#define LINEARSYSTEM_H_

#include <vector>
#include <ginac/ginac.h>

class LinearSystemSet;

class LinearSystem {

private:
	std::vector< std::vector<double> > A; // matrix A
	std::vector< double > b; 		// vector b

	bool isIn(std::vector< double > Ai, const double bi) const;	// check if a constraint is already in
	bool isRedundant(const std::vector< double >& Ai, const double bi) const;

	LinearSystemSet* get_a_finer_covering(LinearSystemSet *tmp_covering,
								          std::vector<std::vector<double> >& A,
										  std::vector<double>& b) const;

public:
	LinearSystem();
	LinearSystem(const LinearSystem& orig);
	LinearSystem(LinearSystem&& orig);
	LinearSystem(const std::vector< std::vector<double> >& A, const std::vector< double >& b);
	LinearSystem(const GiNaC::lst& vars, const GiNaC::lst& constraints);

	/**
	 * Return the template matrix
	 *
	 * @return template matrix
	 */
	inline const std::vector< std::vector<double> >& getA() const
	{
		return this->A;
	}

	/**
	 * Return the offset vector
	 *
	 * @return offset vector
	 */
	inline const std::vector<double>& getb() const
	{
		return this->b;
	}

	const double& getA(const unsigned int i, const unsigned int j) const;
	const double& getb(const unsigned int i) const;

	// optimization functions
	double minLinearSystem(const GiNaC::lst& vars, const GiNaC::ex& obj_fun) const;
	double maxLinearSystem(const GiNaC::lst& vars, const GiNaC::ex& obj_fun) const;

	double minLinearSystem(const std::vector< double >& obj_fun_coeffs) const;
	double maxLinearSystem(const std::vector< double >& obj_fun_coeffs) const;

    // testing methods 
	bool isEmpty(const bool strict_inequality=false) const;	 // determine this LS is empty
	bool operator==(const LinearSystem& LS) const;  // decide whether two systems are equivalent
	bool solutionsAlsoSatisfy(const LinearSystem& LS) const;

	LinearSystemSet* get_a_finer_covering() const;

	// operations on linear system
	LinearSystem intersectWith(const LinearSystem& ls) const;
	std::vector<bool> redundantCons() const;
	LinearSystem& simplify();

	/**
	+ * Remove redundant constraints
	+ */
	inline LinearSystem get_simplified() const
	{
		return LinearSystem(*this).simplify();
	}

	/**
	 * Return the number of variables
	 */
	inline unsigned int dim() const
	{
		if (this->isEmpty()) {
			return 0;
		}
		return this->A[0].size();
	}

	/**
	 * Return the number of inequalities
	 */
	inline unsigned int size() const
	{
		return this->b.size();
	}

	double volBoundingBox();

	void print() const;
	void plotRegion() const;
	void plotRegionToFile(const char *file_name, const char color) const;
	void plotRegionT(const double t) const;
	void plotRegion(const std::vector<int>& rows, const std::vector<int>& cols) const;
};

/**
 * Compute the complementary of a vector of values.
 *
 * @param[in] orig the vector of values to be complementated.
 * @return A vector containing the complementaries of the parameter values
 */
template<typename T>
std::vector<T> get_complementary(const std::vector<T>& orig)
{
   std::vector<T> res{orig};

   for (typename std::vector<T>::iterator it=std::begin(res); it!=std::end(res); ++it) {
	   *it = -*it;
   }

   return res;
}

#endif /* LINEARSYSTEM_H_ */
