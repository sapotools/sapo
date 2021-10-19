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

private:
	unsigned int dim;						// dimension of the parallelotope
	vector<lst> vars;						// variables appearing in generato function
											// vars[0] q: base vertex
											// vars[1] alpha : free variables \in [0,1]
											// vars[2] beta : generator amplitudes
	lst generator_function;				// generator function
	vector< vector<double> > u;			// versors
	vector< vector<double> > template_matrix;	// Template matrix


	vector<double> hyperplaneThroughPts(vector< vector<double> > pts);	// find hyper plane passing through pts
	vector<double> lst2vec(ex list);									// convert a lst to a vector of doubles
	double euclidNorm(vector<double> v);								// compute euclidean norm

	vector< double > actual_base_vertex;
	vector< double > actual_lenghts;

public:

	Parallelotope(vector<lst> vars, vector< vector<double> > u);
	Parallelotope(vector<lst> vars, LinearSystem *LS);

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
	LinearSystem* gen2const(vector<double> q, vector<double> beta);		// from generator to constraints
	poly_values const2gen(LinearSystem *LS);							// from constraints to generators
	LinearSystem* getLS(){ return this->gen2const(this->actual_base_vertex, this->actual_lenghts); };

	vector< double > getBaseVertex(){ return this->actual_base_vertex; };
	vector< double > getLenghts(){ return this->actual_lenghts; };

	vector< vector<double> > getVersors(){ return this->u; };

	int getMaxSize(){ return 0; };
	vector<ex> getConvCombs(int i){ vector<ex> c; return c; };

	virtual ~Parallelotope();
};

#endif /* PARALLELOTOPE_H_ */
