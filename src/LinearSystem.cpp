/**
 * @file LinearSystem.cpp
 * Represent and manipulate a linear system
 * It can be used to represent polytopes (reached states, parameters, etc.)
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @version 0.1
 */

#include "LinearSystem.h"

#define MAX_APPROX_ERROR 1e-6  // necessary for double comparison


/**
 * Optimize a linear system
 *
 * @param[in] A template matrix of the system to optimize
 * @param[in] b offset vector of the system to optimize
 * @param[in] obj_fun objective function
 * @param[in] min_max minimize of maximize Ax<=b (GLP_MIN=min, GLP_MAX=max)
 * @return optimum
 */
double solveLinearSystem(const vector< vector< double > > &A, const vector< double > &b, const vector< double > &obj_fun, const int min_max){

	int num_rows = A.size();
	int num_cols = obj_fun.size();
	int size_lp = num_rows*num_cols;

	int ia[size_lp+1], ja[size_lp+1];
	double ar[size_lp+1];

	glp_prob *lp;
	lp = glp_create_prob();
	glp_set_obj_dir(lp, min_max);

	// Turn off verbose mode
	glp_smcp lp_param;
	glp_init_smcp(&lp_param);
	lp_param.msg_lev = GLP_MSG_ERR;

	glp_add_rows(lp, num_rows);
	for(int i=0; i<num_rows; i++){
		glp_set_row_bnds(lp, i+1, GLP_UP, 0.0, b[i]);
	}

	glp_add_cols(lp, num_cols);
	for(int i=0; i<num_cols; i++){
		glp_set_col_bnds(lp, i+1, GLP_FR, 0.0, 0.0);
	}

	for(int i=0; i<num_cols; i++){
		glp_set_obj_coef(lp, i+1, obj_fun[i]);
	}

	int k=1;
	for(int i=0; i<num_rows; i++){
		for(int j=0; j<num_cols; j++){
			ia[k] = i+1, ja[k] = j+1, ar[k] = A[i][j]; /* a[i+1,j+1] = A[i][j] */
			k++;
		}
	}

	glp_load_matrix(lp, size_lp, ia, ja, ar);
	glp_exact(lp, &lp_param);
//	glp_simplex(lp, &lp_param);

	double res = glp_get_obj_val(lp);
	glp_delete_prob(lp);
	glp_free_env();
	return res;
}

/**
 * Check if if a vector is null, i.e.,
 * it's a vector of zeros (used to detected useless constraints)
 *
 * @param[in] line vector to test
 * @return true is the vector is nulle
 */
bool zeroLine(const vector<double> &line) {

	bool zeros = true;
	int i=0;
	while(zeros && i<(signed)line.size()){
		zeros = zeros && (abs(line[i]) < MAX_APPROX_ERROR);
		i++;
	}
	return zeros;
}

/**
 * Constructor that instantiates a linear system
 *
 * @param[in] A template matrix
 * @param[in] b offset vector
 */
LinearSystem::LinearSystem(vector< vector<double> > A, vector< double > b){

	bool smart_insert = false;

	if(!smart_insert){
		this->A = A;
		this->b = b;
	}else{
		for(int i=0; i<(signed)A.size(); i++){
			if(!this->isIn(A[i],b[i]) && (!zeroLine(A[i]))){
				this->A.push_back(A[i]);
				this->b.push_back(b[i]);
			}
		}

	}
	this->n_vars = this->A[0].size();
}

/**
 * Constructor that instantiates an empty linear system
 */
LinearSystem::LinearSystem(){
	vector< vector<double> > A;
	vector< double > b;
	this->A = A;
	this->b = b;
	this->n_vars = 0;

}

/**
 * Check if a constraint belongs to the linear system
 *
 * @param[in] Ai direction
 * @param[in] bi offset
 * @returns true is Ai x <= b is in the linear system
 */
bool LinearSystem::isIn(vector< double > Ai, const double bi) const{
	Ai.push_back(bi);

	for( int i=0; i<(signed)this->A.size(); i++ ){
		vector< double > line = this->A[i];
		line.push_back(this->b[i]);
		bool is_in = true;
		for(int j=0; j<(signed)Ai.size(); j++){
			is_in = is_in && (abs(Ai[j] - line[j]) < MAX_APPROX_ERROR);
		}
		if(is_in){ return true; }
	}
	return false;
}

/**
 * Constructor that instantiates a linear system from a set of symbolic expressions
 *
 * @param[in] vars list of variables appearing in the constraints
 * @param[in] constraints symbolic constraints
 */
LinearSystem::LinearSystem(lst vars, lst constraints) {

	this->vars = vars;
	this->constraints = constraints;
	this->constraints.unique();

	this->n_vars = this->vars.nops();

	initLS();	// initialize Linear System
}

/**
 * Initialize a linear system extracting template and offsets from symbolic expressions
 */
void LinearSystem::initLS(){

	for(int i=0; i<(signed)this->constraints.nops(); i++){

		vector<double> Ai;
		ex const_term = this->constraints[i];

		for(int j=0; j<(signed)this->vars.nops(); j++){

			// Extract the coefficient of the i-th variable (grade 1)
			double coeff = ex_to<numeric>(evalf(this->constraints[i].coeff(this->vars[j],1))).to_double();
			Ai.push_back(coeff);

			// Project to obtain the constant term
			const_term = const_term.coeff(this->vars[j],0);
		}

		double bi = ex_to<numeric>(evalf(const_term)).to_double();

		if(!this->isIn(Ai,-bi)){
			this->A.push_back(Ai);
			this->b.push_back(-bi);
		}
	}

}

/**
 * Return the (i,j) element of the template matrix
 *
 * @param[in] i row index
 * @param[in] j column index
 * @return (i,j) element
 */
const double& LinearSystem::getA(int i, int j) const {
	if(( 0<= i ) && (i < (signed)this->A.size())){
		if(( 0<= j ) && (j < (signed)this->A[j].size())){
			return this->A[i][j];
		}
	}
	cout<<"LinearSystem::getA : i and j must be within the LS->A size";
	exit (EXIT_FAILURE);
}

/**
 * Return the i-th element of the offset vector
 *
 * @param[in] i column index
 * @return i-th element
 */
const double& LinearSystem::getb(int i) const {
	if(( 0<= i ) && (i < (signed)this->b.size())){
			return this->b[i];
	}
	cout<<"LinearSystem::getb : i and j must be within the LS->b size";
	exit (EXIT_FAILURE);
}

/**
 * Determine whether this linear system is empty or not, i.e.,
 * the linear system has solutions
 *
 * @param[in] strict_inequality specifies whether the linear system is a 
 * 						strict inequality (i.e., Ax < b)
 * @return true if the linear system is empty
 */
bool LinearSystem::isEmpty(const bool strict_inequality) const {

	vector< vector< double > > extA = this->A;
	vector< double > obj_fun (this->n_vars, 0);
	obj_fun.push_back(1);

	// Add an extra variable to the linear system
	for(int i=0; i<(signed)extA.size(); i++){
		extA[i].push_back(-1);
	}

	const double z = solveLinearSystem(extA,this->b,obj_fun,GLP_MIN);

	return (z>0)||(strict_inequality&&(z>=0));
}

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

/**
 * Determine all the solutions of this linear system are also 
 * solutions for another linear system.
 *
 * @param[in] LS a linear system
 * @return true if all the solutions of this object are also 
 * 				solutions for the parameter.
 */
bool LinearSystem::solutionsAlsoSatisfy(const LinearSystem& LS) const {

	LinearSystem extLS(LS.A, LS.b);
	extLS.A.push_back(vector<double>(1, 0));
	extLS.b.push_back(0);

	for (int i=0; i<this->size(); i++){
		extLS.A[LS.size()] = get_complementary(this->A[i]);
		extLS.b[LS.size()] = -this->b[i];

		if (!extLS.isEmpty(true)) {
			return false;
		}
	}

	return true;
}


/**
 * Determine whether two linear systems are equivalent
 *
 * @param[in] LS a linear system to be compared this object.
 * @return true if this object is equivalent to the parameter
 */
bool LinearSystem::operator==(const LinearSystem& LS) const {

	return this->solutionsAlsoSatisfy(LS) && LS.solutionsAlsoSatisfy(*this);
}

/**
 * Minimize the linear system
 *
 * @param[in] vars list of variables
 * @param[in] obj_fun objective function
 * @return minimum
 */
double LinearSystem::minLinearSystem(const lst& vars, const ex& obj_fun) const {

	vector< double > obj_fun_coeffs;
	ex const_term = obj_fun;

	// Extract the coefficient of the i-th variable (grade 1)
	for(int i=0; i<(signed)vars.nops(); i++){

		double coeff = ex_to<numeric>(evalf(obj_fun.coeff(vars[i],1))).to_double();

		obj_fun_coeffs.push_back(coeff);
		const_term = const_term.coeff(vars[i],0);

	}

	const double c = ex_to<numeric>(evalf(const_term)).to_double();
	const double min = solveLinearSystem(this->A,this->b,obj_fun_coeffs,GLP_MIN);

	return (min+c);

}

/**
 * Maximize the linear system
 *
 * @param[in] obj_fun objective function
 * @return maximum
 */
double LinearSystem::maxLinearSystem(const vector< double >& obj_fun_coeffs) const {
	return solveLinearSystem(this->A,this->b,obj_fun_coeffs,GLP_MAX);
}

/**
 * Maximize the linear system
 *
 * @param[in] vars list of variables
 * @param[in] obj_fun objective function
 * @return maximum
 */
double LinearSystem::maxLinearSystem(const lst& vars, const ex& obj_fun) const {

	vector< double > obj_fun_coeffs;
	ex const_term = obj_fun;

	// Extract the coefficient of the i-th variable (grade 1)
	for(int i=0; i<(signed)vars.nops(); i++){

		double coeff = ex_to<numeric>(evalf(obj_fun.coeff(vars[i],1))).to_double();
		obj_fun_coeffs.push_back(coeff);
		const_term = const_term.coeff(vars[i],0);

	}

	const double c = ex_to<numeric>(evalf(const_term)).to_double();
	const double max = solveLinearSystem(this->A,this->b,obj_fun_coeffs,GLP_MAX);

	return (max+c);
}

/**
 * Create a new linear system by concatenating this LS and the specified one
 *
 * @param[in] LS linear system to be appended
 * @return linear system obtained by merge
 */
LinearSystem* LinearSystem::intersectWith(const LinearSystem *LS) const {
	LinearSystem *result = new LinearSystem(this->A, this->b);

	for(int i=0; i<LS->size(); i++){
		if( !this->isIn(LS->A[i], LS->b[i]) ){		// check for duplicates
			result->A.push_back( LS->A[i] );
			result->b.push_back( LS->b[i] );
		}
	}

	return result;
}

/**
 * Determine the redundant constraint of the linear system
 *
 * @return boolean vector (true is if i-th constrain is redundant)
 */
vector<bool> LinearSystem::redundantCons() const{

	vector<bool> redun (this->size(), false);

	for(int i=0; i<this->size(); i++){
		double max = this->maxLinearSystem(this->A[i]);
		if (abs(max - this->b[i]) > MAX_APPROX_ERROR) {  /* This should be max != this->b[i], 
								    however, due to double approximation 
								    errors, testing whether the distance 
								    between max and this->b[i] is 
								    greater than a fixed positive 
								    approximation constant is more 
								    conservative */
			redun[i] = true;
		} else {
			auto max_coeff = max_element(std::begin(A[i]), std::end(A[i]));
			auto min_coeff = min_element(std::begin(A[i]), std::end(A[i]));
			redun[i] = ((*max_coeff==*min_coeff)&&(*min_coeff==0)&&(b[i]>=0));
		}
	}
	return redun;
}

/**
+ * Remove redundant constraints
+ */
LinearSystem *LinearSystem::simplify() {
	vector<bool> R = this->redundantCons();

	for (int i = R.size()-1; i>=0; i--) {
		if (R[i])
		{
			(this->A).erase((this->A).begin()+i);
			(this->b).erase((this->b).begin()+i);
		}	
	}
	
	return this;
}

/**
 * Determine the volume of the bounding box of the linear system
 *
 * @return volume of the bounding box
 */
double LinearSystem::volBoundingBox(){

	vector<double> zeros (this->dim(),0);
	double vol = 1;

	for(int i=0; i<this->dim(); i++){
		vector<double> facet = zeros;
		facet[i] = 1;
		const double b_plus = solveLinearSystem(this->A,this->b,facet,GLP_MAX);
		facet[i] = -1;
		const double b_minus = solveLinearSystem(this->A,this->b,facet,GLP_MAX);
		vol = vol*(b_plus+b_minus);
	}

	return vol;
}

/**
 * Print the linear system
 */
void LinearSystem::print() const {
	for(int i=0; i<(signed)this->A.size(); i++){
		for(int j=0; j<(signed)this->A[i].size(); j++){
			cout<<this->A[i][j] << (j == (signed)this->A[i].size()-1 ? "" : " ");
		}
		cout<<" <= "<<this->b[i]<<"\n";
	}
	cout<<"\n";
}

/**
 * Print the linear system in Matlab format (for plotregion script)
 */
void LinearSystem::plotRegion() const {

	if(this->dim() > 3){
		cout<<"LinearSystem::plotRegion : maximum 3d sets are allowed";
		exit (EXIT_FAILURE);
	}

	cout<<"Ab = [\n";
	for(int i=0; i<(signed)this->A.size(); i++){
		for(int j=0; j<(signed)this->A[i].size(); j++){
			cout<<this->A[i][j]<<" ";
		}
		cout<<" "<<this->b[i]<<";\n";
	}
	cout<<"];\n";
	cout<<"plotregion(-Ab(:,1:"<< this->A[0].size() <<"),-Ab(:,"<<this->A[0].size()+1<<"),[],[],color);\n";

}

/**
 * Print the linear system in Matlab format (for plotregion script) into a file
 *
 * @param[in] file_name name of the file
 * @param[in] color color of the polytope to plot
 */
void LinearSystem::plotRegionToFile(const char *file_name, const char color) const {

	if(this->dim() > 3){
		cout<<"LinearSystem::plotRegion : maximum 3d sets are allowed";
		exit (EXIT_FAILURE);
	}


	ofstream matlab_script;
	matlab_script.open (file_name, ios_base::app);

	matlab_script<<"Ab = [\n";
	for(int i=0; i<(signed)this->A.size(); i++){
		for(int j=0; j<(signed)this->A[i].size(); j++){
			matlab_script<<this->A[i][j]<<" ";
		}
		matlab_script<<" "<<this->b[i]<<";\n";
	}
	matlab_script<<"];\n";
	matlab_script<<"plotregion(-Ab(:,1:"<< this->A[0].size() <<"),-Ab(:,"<<this->A[0].size()+1<<"),[],[],'"<<color<<"');\n";
	matlab_script.close();

}

/**
 * Print the 2d linear system in Matlab format (for plotregion script) over time
 *
 * @param[in] t thickness of the set to plot
 */
void LinearSystem::plotRegionT(const double t) const {
	if(this->dim() > 2){
		cout<<"LinearSystem::plotRegionT : maximum 2d sets are allowed";
		exit (EXIT_FAILURE);
	}

	cout<<"Ab = [\n";
	cout<<" 1 ";
	for(int j=0; j<(signed)this->A[0].size(); j++){
		cout<<" 0 ";
	}
	cout<<t<<";\n";
	cout<<" -1 ";
	for(int j=0; j<(signed)this->A[0].size(); j++){
		cout<<" 0 ";
	}
	cout<<-t<<";\n";

	for(int i=0; i<(signed)this->A.size(); i++){
		cout<<" 0 ";
		for(int j=0; j<(signed)this->A[i].size(); j++){
			cout<<this->A[i][j]<<" ";
		}
		cout<<this->b[i]<<";\n";
	}

	cout<<"];\n";
	cout<<"plotregion(-Ab(:,1:3),-Ab(:,4),[],[],color);\n";

}

/**
 * Print the specified projections in Matlab format (for plotregion script) into a file
 *
 * @param[in] rows rows to be plot
 * @param[in] cols colors of the plots
 */
void LinearSystem::plotRegion(const vector<int>& rows, const vector<int>& cols) const{

		if(cols.size() > 3){
			cout<<"LinearSystem::plotRegion : cols maximum 3d sets are allowed";
			exit (EXIT_FAILURE);
		}

		cout<<"Ab = [\n";
		for(int i=0; i<(signed)rows.size(); i++){
			for(int j=0; j<(signed)cols.size(); j++){
				cout<<this->A[rows[i]][cols[j]]<<" ";
			}
			cout<<" "<<this->b[rows[i]]<<";\n";
		}
		cout<<"];\n";
		cout<<"plotregion(-Ab(:,1:"<< cols.size() <<"),-Ab(:,"<<cols.size()+1<<"),[],[],color);\n";
}


LinearSystem::~LinearSystem() {
	// TODO Auto-generated destructor stub
}

