/**
 * @file Parallelotope.h
 * Represent and manipulate a parallelotope
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @version 0.1
 */

#ifndef PARALLELOTOPE_H_
#define PARALLELOTOPE_H_

#include "Common.h"
#include "LinearSystem.h"

class Parallelotope{
	using Vector = vector<double>;
	using Matrix = vector< Vector >;
private:
	unsigned int dim;						// dimension of the parallelotope
	vector<lst> vars;						// variables appearing in generato function
											// vars[0] q: base vertex
											// vars[1] alpha : free variables \in [0,1]
											// vars[2] beta : generator amplitudes
	lst generator_function;				// generator function
	Matrix u;			// versors
	Matrix template_matrix;	// Template matrix


	Vector hyperplaneThroughPts(const vector<Vector>& pts) const;	// find hyper plane passing through pts
	Vector lst2vec(const ex& list) const;						// convert a lst to a vector of doubles
	double euclidNorm(const Vector& v) const;		// compute euclidean norm

	vector< double > actual_base_vertex;
	vector< double > actual_lenghts;

public:

	Parallelotope(const vector<lst>& vars, const Matrix& u);
	Parallelotope(const vector<lst>& vars, const Matrix& template_matrix, const Vector& offset);
	Parallelotope(const vector<lst>& vars, const LinearSystem& ls);

	/**
	 * Get the generator functions
	 *
	 * @returns generator functions
	 */
	inline const lst& getGeneratorFunction() const
	{
		return this->generator_function;
	}

	/**
	 * Get variables of base vertex
	 *
	 * @returns base vertex variables
	 */
	inline const lst& getQ() const
	{
		return this->vars[0];
	}

	/**
	 * Get free variables
	 *
	 * @returns free variables
	 */
	inline const lst& getAlpha() const
	{
		return this->vars[1];
	}

	/**
	 * Get variables of generator lengths
	 *
	 * @returns generator lengths variables
	 */
	inline const lst& getBeta() const
	{
		return this->vars[2];
	}

	/**
	 * Get the parallelotope dimension
	 *
	 * @returns parallelotope dimension
	 */
	inline const unsigned int& getDim() const { return this->dim; }

	/**
	 * Get the template of the parallelotope
	 *
	 * @returns parallelotope template
	 */
	inline const vector< vector< double > >& getTemplate() const
	{
		return this->template_matrix;
	}

	// Representation conversion
	LinearSystem gen2const(const Vector& q, const Vector& beta) const;		// from generator to constraints
	poly_values const2gen(LinearSystem *LS) const;							// from constraints to generators
	inline LinearSystem getLinearSystem() const
	{
		return this->gen2const(this->actual_base_vertex, this->actual_lenghts);
	}

	inline const vector< double >& getBaseVertex() const
	{
		return this->actual_base_vertex; 
	}

	inline const vector< double >& getLenghts() const
	{
		return this->actual_lenghts;
	}

	inline Matrix getVersors() const
	{
		return this->u;
	}

	friend void swap(Parallelotope& A, Parallelotope& B);

	virtual ~Parallelotope();
};

#endif /* PARALLELOTOPE_H_ */
