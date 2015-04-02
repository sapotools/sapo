/*
 * Parallelotope.cpp
 *
 *  Created on: Oct 30, 2014
 *      Author: Tommaso Dreossi
 */

#include "Parallelotope.h"

Parallelotope::Parallelotope(vector<lst> vars, vector< vector<double> > u) {

	if(vars.size() != 3){
		cout<<"Parallelotope::Parallelotope : vars must contain 3 collections of variable names (q,alpha,beta)";
		exit (EXIT_FAILURE);
	}

	this->vars.push_back(vars[0]);
	this->vars.push_back(vars[1]);
	this->vars.push_back(vars[2]);

	// get the dimension of the parallelotope
	this->dim = vars[0].nops();
	// and store its variable names
	for(int i=0; i<3; i++){
		if((signed)vars[i].nops() != this->dim){
			cout<<"Parallelotope::Parallelotope : vars["<<i<<"] must have "<<this->dim<<" variables";
			exit (EXIT_FAILURE);
		}
	}

	// extract variable names
	lst q = this->vars[0];
	lst alpha = this->vars[1];
	lst beta = this->vars[2];

	// initialize generator function
	for(int i=0; i<this->dim; i++){
		this->generator_function.append(q[i]);
	}

	// create the generation function accumulating the versor values
	for(int i=0; i<this->dim; i++){

		// check dimension of versors
		if((signed)u[i].size() != this->dim){
			cout<<"Parallelotope::Parallelotope : dim and ui dimensions must agree";
			exit (EXIT_FAILURE);
		}

		// Normalize the versors and generatre the generator function
		//double norm = euclidNorm(u[i]);
		//vector< double > norm_versor;

		for(int j=0; j<this->dim; j++){
			//norm_versor.push_back(u[i][j]/norm);
			this->generator_function[j] = this->generator_function[j] + alpha[i]*beta[i]*u[i][j];
		}
		// store the i-th versor
		this->u.push_back(u[i]);
	}

	// Initialize the template matrix
	vector< double > base_vertex (this->dim,0);
	vector<double> lenghts (this->dim,1);
	LinearSystem *LS = this->gen2const(base_vertex,lenghts);

	this->template_matrix = LS->getA();

}

// Convert from generator to constraints representation
// q : numeric base vertex
// beta : numeric lenghts
LinearSystem* Parallelotope::gen2const(vector<double> q, vector<double> beta){

	if((signed)q.size() != this->dim){
		cout<<"Parallelotope::gen2const : q must have dimension "<<this->dim;
		exit (EXIT_FAILURE);
	}


	vector< vector<double> > hps;	// hyperplane equations

	// first set of hyperplanes
	for(int i=0; i<this->dim; i++){

		vector< vector<double> > pts;
		pts.push_back(q);							// add base vertex

		for(int j=0; j<this->dim; j++){ 			// for all the generators u
			if( i != j ){
				vector<double> p;
				for(int k =0; k<this->dim; k++){	// coordinate of the point
					p.push_back(q[k] + this->u[j][k]*beta[j]);
				}
				pts.push_back(p);
			}
		}
		hps.push_back(hyperplaneThroughPts(pts));
	}

	//second set of hyperplanes
	for(int i=0; i<this->dim; i++){

			vector< vector<double> > pts;

			vector<double> qt; // traslate q
			for(int j=0; j<this->dim; j++){
				qt.push_back(q[j] + beta[i]*this->u[i][j]);
			}
			pts.push_back(qt);							// add base vertex

			for(int j=0; j<this->dim; j++){ 			// for all the generators u
				if( i != j ){
					vector<double> p;
					for(int k =0; k<this->dim; k++){	// coordinate of the point
						p.push_back( (q[k] + this->u[j][k]*beta[j]) + beta[i]*this->u[i][k]);
					}
					pts.push_back(p);
				}
			}

			hps.push_back(hyperplaneThroughPts(pts));

	}

	vector< vector<double> > Lambda;
	vector< double > d (this->dim*2,0);

	// initialize template Lambda
	vector<double> Lambda_i (this->dim,0);
	for(int i=0; i<(this->dim)*2; i++){
		Lambda.push_back(Lambda_i);
	}

	for(int i=0; i<this->dim; i++){

		// find hyperplane with smaller direction
		double b1 = -hps[i][hps[i].size()-1];
		double b2 = -hps[i+this->dim][hps[i].size()-1];

		if(b1 > b2){
			d[i] = b1;
			d[i+this->dim] = -b2;
			for(int j=0; j<this->dim; j++){
				Lambda[i][j] = hps[i][j];
				Lambda[i+this->dim][j] = -hps[i+this->dim][j];
			}
		}else{
			d[i] = -b1;
			d[i+this->dim] = b2;
			for(int j=0; j<this->dim; j++){
				Lambda[i][j] = -hps[i][j];
				Lambda[i+this->dim][j] = hps[i+this->dim][j];
			}
		}
	}

	LinearSystem *LS = new LinearSystem(Lambda,d);

	return LS;
}


// Compute the equation of the hyperplane passing through the points in pts
// the equation is res[0]*x_0 + .. + res[n]*x_n + res[n+1] = 0
// pts : matrix containing the points
vector<double> Parallelotope::hyperplaneThroughPts(vector< vector<double> > pts){

	if(pts.size() != pts[0].size()){
		cout<<"Parallelotope::hyperplaneThroughPts : pts must contain n n-dimensional points";
		exit (EXIT_FAILURE);
	}

	vector< vector<double> > A;
	// build the linear system Ax = 0
	for(int i=1; i<(signed)pts.size(); i++){
		vector<double> Ai;
		for(int j=0; j < this->dim; j++){
			Ai.push_back(pts[0][j] - pts[i][j]);
		}
		A.push_back(Ai);
	}

	// Build the linear system to find the normal vector
	lst LS;
	lst a = this->vars[1];

	for(int i=0; i<this->dim-1; i++){

		ex eq = 0;
		lst sub;

		for(int j=0; j<this->dim; j++){
			eq = eq + A[i][j]*a[j];
		}

		eq = eq == 0;
		LS.append(eq);

	}

	// Solve the linear system
	ex solLS = lsolve(LS,a);

	// Build the linear inequality
	lst x = this->vars[0];
	ex eqr = 0;
	vector<double> lambda (this->dim+1,0);

	ex sub;
	int sub_idx;
	// search for the tautology
	for(int i=0; i<(signed)solLS.nops(); i++){
		if(solLS[i].is_equal(a[i] == a[i])){
			sub = a[i] == 1;
			sub_idx = i;
		}
	}

	lambda[sub_idx] = 1;
	for(int i=0; i<this->dim; i++){
		if(i != sub_idx){
			lambda[i] = ex_to<numeric>(evalf(a[i].subs(solLS[i].subs(sub)))).to_double();
		}
		eqr = eqr + a[i].subs(a[i] == lambda[i])*pts[0][i];
	}

	lambda[this->dim] = -ex_to<numeric>(evalf(eqr)).to_double();

	return lambda;
}


// Convert the constrain representation to the generator one
// constr : LinearSystem
poly_values Parallelotope::const2gen(LinearSystem *constr){

	vector< vector<double> > Lambda = constr->getA();
	vector<double> d = constr->getb();
	vector< vector<double> > vertices;

	// find base vertex
	//build the linear system
	ex q = this->vars[0];
	lst LS;

	for(int i=0; i<this->dim; i++){
		ex eq = 0;
		for(int j=0; j<(signed)Lambda[i].size(); j++){
			eq = eq + Lambda[i][j]*q[j];
		}
		eq = eq == d[i];
		LS.append(eq);
	}

	ex solLS = lsolve(LS,q);
	vertices.push_back( lst2vec(q.subs(solLS)) ); // store the base_vertex

	// Compute the vertices v
	for(int k=0; k<this->dim; k++){
		ex a = this->vars[1];
		lst LS;

		for(int i=0; i<this->dim; i++){
			ex eq = 0;
			for(int j=0; j<(signed)Lambda[i].size(); j++){
				eq = eq + Lambda[i][j]*a[j];
			}
			if(i != k){
				eq = eq == d[i];
			}else{
				eq = eq == -d[i+this->dim];
			}
			LS.append(eq);
		}

		ex solLS = lsolve(LS,a);
		vertices.push_back( lst2vec(a.subs(solLS)) ); // store the i-th vertex
	}

	// Compute the generators
	vector< vector<double> > g;
	for(int i=0; i<this->dim; i++){
		vector<double> gi;
		for(int j=0; j<this->dim; j++){
			gi.push_back( vertices[i+1][j] - vertices[0][j] );
		}
		g.push_back(gi);
	}

	// Compute the generators lengths
	vector< double > lengths;
	for(int i=0; i<this->dim; i++){
		lengths.push_back(euclidNorm(g[i]));
	}


//	// Find the versors (optional since versors are specified by the user)
//	vector< vector< double > > versors;
//	for(int i=0; i<this->dim; i++){
//		vector<double> versori;
//		for(int j=0; j<this->dim; j++){
//			versori.push_back( g[i][j]/lengths[i] );
//		}
//		versors.push_back(versori);
//	}

	// Return the conversion
	poly_values result;
	result.base_vertex = vertices[0];
	result.lenghts = lengths;
	//result.versors = versors;

	return result;
}

// Convert a symbolic list to a vector
vector<double> Parallelotope::lst2vec(ex list){

	vector< double > res;

	for(int i=0; i<(signed)list.nops(); i++)
		res.push_back( ex_to<numeric>(evalf( list[i] )).to_double() );

	return res;
}

// Compute the Euclidean norm of the vector v
double Parallelotope::euclidNorm(vector<double> v){

	double norm = 0;

	for(int i=0; i<(signed)v.size(); i++){
		norm = norm + (v[i]*v[i]);
	}

	return sqrt(norm);
}


Parallelotope::~Parallelotope() {
	// TODO Auto-generated destructor stub
}

