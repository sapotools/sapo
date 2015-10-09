/*
 * BaseConverter.cpp
 *
 *  Created on: Oct 23, 2014
 *      Author: Tommaso Dreossi
 */

#include "BaseConverter.h"

BaseConverter::BaseConverter(lst vars, ex polynomial) {

	this->vars = vars;
	this->polynomial = polynomial;

	// Put the polynomial in extended form and extract variables degrees
	this->polynomial = this->polynomial.expand();
	for(int i=0; i<(signed)this->vars.nops();i++){
		this->degrees.push_back(this->polynomial.degree(this->vars[i]));
	}

//	cout<<"\nPoly: "<<polynomial<<"\n";
//		cout<<"\nDegs: ";
//		for(int i=0; i<this->degrees.size(); i++){
//			cout<<this->degrees[i];
//		}
//		cout<<"\n";


	// Initialize the degree shifts
	initShifts();

	for(int i=0;i<this->shifts[0];i++){
		this->coeffs.push_back(0);
	}

	// Initialize the coefficients vector
	vector<int> multi_index;
	extractCoeffs(this->polynomial,0,multi_index);

}

BaseConverter::BaseConverter(lst vars, ex polynomial, vector<int> degrees) {

	this->vars = vars;
	this->polynomial = polynomial;

	// Put the polynomial in extended form and extract variables degrees
	this->polynomial = this->polynomial.expand();
	this->degrees = degrees;

	// Initialize the degree shifts
	initShifts();

	for(int i=0;i<this->shifts[0];i++){
		this->coeffs.push_back(0);
	}

	// Initialize the coefficients vector
	vector<int> multi_index;
	extractCoeffs(this->polynomial,0,multi_index);

//	cout<<"\nDegs: ";
//	for(int i=0; i<this->degrees.size(); i++){
//		cout<<this->degrees[i];
//	}
//	cout<<"\n";

}

// for rational polynomials
BaseConverter::BaseConverter(lst vars, ex num, ex denom){

	this->vars = vars;
	this->num = num;
	this->denom = denom;

}

// Initialize the degree shift vector used to extract the multi-indices
void BaseConverter::initShifts(){

	vector<int> shifts(this->degrees.size(),0);
	this->shifts = shifts;

	this->shifts[this->degrees.size()-1] = this->degrees[this->degrees.size()-1] + 1;

	for(int i=this->degrees.size()-2; i>=0; i--){
		this->shifts[i] = (this->degrees[i]+1)*this->shifts[i+1];
	}

}


// Convert a multi-index in a coefficient position index
int BaseConverter::multi_index2pos(vector<int> multi_index){

	if( multi_index.size() != this->degrees.size() ){
		cout<<"BaseConverter::multi_index2pos : multi_index and degrees must have same dimension";
		exit (EXIT_FAILURE);
	}

	// Compute the position
	int position = multi_index[multi_index.size()-1];
	for(int i=multi_index.size()-2; i>=0; i--){
		position = position + multi_index[i]*this->shifts[i+1];
	}

	return position;

}

// Convert a position to a multi-index
vector<int> BaseConverter::pos2multi_index(int position){

	vector<int> multi_index(this->degrees.size(),0);

	for(int i=0; i<(signed)multi_index.size() - 1; i++){
		multi_index[i] = (int)position / this->shifts[i+1];
		position = position % this->shifts[i+1];
	}

	multi_index[multi_index.size()-1] = position;

	return multi_index;
}


// Extract recursively the coefficients of the polynomial and populate the vector coeffs
void BaseConverter::extractCoeffs(ex polynomial,int var_idx,vector<int> multi_index){

	// Get the variable index
	multi_index.push_back(0);

	// Base case - there's only one variable
	if( var_idx == (signed)this->vars.nops()-1 ){

		for(int i=0; i<=this->degrees[var_idx]; i++){
			ex coeff = polynomial.coeff(this->vars[var_idx],i); 	// Extract the coefficient
			multi_index[multi_index.size()-1] = i;					// and add it to the coeff table
			this->coeffs[multi_index2pos(multi_index)] = coeff;
		}
	}else{
		for(int i=0; i<=this->degrees[var_idx]; i++){
			ex sub_poly = polynomial.coeff(this->vars[var_idx],i); 	// Extract the sub-polynomials
			multi_index[multi_index.size()-1] = i;
			extractCoeffs(sub_poly,var_idx+1,multi_index);		// Recursive call
		}
	}

}

// Calculate the binomial coefficient with multiplicative formula
int BaseConverter::nChoosek(int n, int k){

	if( k > n ){
		cout<<"BaseConverter::nChoosek : n must be larger equal then k";
		exit (EXIT_FAILURE);
	}

	int res = 1;
	for(int i=1; i<=k; i++){
		res = res*(n-(k-i))/i;
	}
	return res;
}

// Calculate the binomial coefficient of two multi-indeces
int BaseConverter::multi_index_nChoosek(vector<int> n, vector<int> k){

	if(n.size() != k.size()){
		cout<<"BaseConverter::multi_index_nChoosek : n and k must have same dimension";
		exit (EXIT_FAILURE);
	}

	int res = 1;
	for(int i=0; i<(signed)n.size(); i++){
		res = res*nChoosek(n[i],k[i]);
	}

	return res;

}

// Check whether b dominates a
bool BaseConverter::multi_index_leq(vector<int> a, vector<int> b){

	bool leq = true;
	int i = 0;

	while(i<a.size() && leq){
		leq = leq && (a[i] <= b[i]);
		i++;
	}


	return leq;

}

// Compute the mi-th Bernstein coefficient
ex BaseConverter::bernCoeff(vector<int> mi){

	int i = multi_index2pos(mi);
	ex coeff = 0;

	for(int j = 0; j <= i; j++){

		vector<int> mj = pos2multi_index(j);
		if( multi_index_leq(mj,mi) ){

			int ichoosej = multi_index_nChoosek(mi,mj);
			int dchoosej = multi_index_nChoosek(this->degrees,mj);

			coeff = coeff + ((ex)ichoosej/(ex)dchoosej)*this->coeffs[j];
		}
	}

	return coeff;

}

// Compute the the list of Bernstein control points
lst BaseConverter::getBernCoeffs(){

	lst bern_coeffs;

	for(int i=0; i<(signed)this->coeffs.size(); i++){
		bern_coeffs.append(bernCoeff(pos2multi_index(i)));
	}

	return bern_coeffs;

}

// Compute the the list of Bernstein control points for rational polynomials
lst BaseConverter::getRationalBernCoeffs(){

	lst bern_coeffs;

	cout<<"Degrees: ";
	vector<int> degs;
	for(int i=0; i<(signed)this->vars.nops();i++){
		degs.push_back(max(this->num.degree(this->vars[i]),this->denom.degree(this->vars[i])));
		cout<<degs[i]<<", ";
	}

	BaseConverter *num_conv = new BaseConverter(this->vars,this->num,degs);
	BaseConverter *denom_conv = new BaseConverter(this->vars,this->denom,degs);

	lst num_bern_coeffs = num_conv->getBernCoeffs();
	lst denom_bern_coeffs = denom_conv->getBernCoeffs();

	for(int i=0; i<num_bern_coeffs.nops(); i++){
		if(denom_bern_coeffs[i] != 0){ // skip negative denominators
			bern_coeffs.append(num_bern_coeffs[i]/denom_bern_coeffs[i]);
		}
	}

	// eliminate duplicates
	bern_coeffs.unique();
	cout<<"(Total points:"<<bern_coeffs.nops()<<")\n";
	return bern_coeffs;

}

// Compute the the list of Bernstein control points
vector< vector< int > > BaseConverter::getMultiIdxList(){

	vector< vector< int > > mi_list;

	for(int i=0; i<(signed)this->coeffs.size(); i++){
		mi_list.push_back(pos2multi_index(i));
	}
	return mi_list;

}

// Split
void BaseConverter::split(int direction, double split_point){

	if(( direction < 0 ) || (direction >= this->vars.nops())){
		cout<<"BaseConverter::split : split direction must be between 0 and "<<this->vars.nops();
		exit (EXIT_FAILURE);
	}

	if(( split_point < 0 ) || (split_point > 1)){
		cout<<"BaseConverter::split : split_point must be between 0 and 1";
		exit (EXIT_FAILURE);
	}

	// Extract list of multi indices and max degree of direction
	vector< vector< int > > multi_index_list = this->getMultiIdxList();
	int nr = this->degrees[direction];

	lst B;
	B = this->getBernCoeffs();

	for(int mi=0; mi<B.nops(); mi++ ){

		vector< int > i = this->pos2multi_index(mi);

		for(int k=1; k<this->degrees[direction]; k++){

			for(int j=0; j<this->degrees.size(); j++){

				cout<<"j:"<<j<<" dir: "<<direction<<"\n"<<" dir: "<<direction<<"\n";

				if( j!= direction ){

					for(int ij=0; ij < this->degrees[j]; ij++){

						if( ij<k ){
							//B[mi] = B[mi];
						}else{
							vector< int > ni = i;
							ni[j] = ij-1;
							int bi = this->multi_index2pos(ni);
							B[mi] = (1-split_point)*B[bi] + split_point*B[mi];
						}
					}
				}
			}
		}
	}

	cout<<"\n"<<B<<"\n";

}


// remove monomes with zero coefficient
pair< vector<ex>, vector< vector<int> > > BaseConverter::compressZeroCoeffs(){

	vector<ex> comp_coeffs;
	vector< vector<int> > comp_degs;
	pair< vector<ex>, vector< vector<int> > > compression;

	for(int i=0; i<(signed)this->coeffs.size(); i++){
		if(this->coeffs[i] != 0){
			comp_coeffs.push_back(this->coeffs[i]);
			comp_degs.push_back(pos2multi_index(i));
		}
	}

	compression.first = comp_coeffs;
	compression.second = comp_degs;
	return compression;

}

void BaseConverter::implicitMaxIndex(){

	pair< vector<ex>, vector< vector<int> > > compression = this->compressZeroCoeffs();
	// Check the three properties (Uniqueness,Monotonicity, and Dominance)
	vector<ex> coeffs = compression.first;
	vector< vector<int> > multi_index = compression.second;
	vector<int> implicit_max (this->vars.nops(),-1);

	// Check uniqueness
	vector< int > unique (this->vars.nops(), 0);
	vector< int > unique_deg (this->vars.nops(), -1);
	for(int i=0; i<(signed)this->vars.nops(); i++){
		int j=0;
		while(( j< (signed)coeffs.size() ) && ( unique[i] < 2 )){
			if( multi_index[j][i] > 0 ){
				unique[i] = unique[i] + 1;
				unique_deg[i] = multi_index[j][i];
			}
			j++;
		}
	}
	for(int i=0; i<(signed)unique.size(); i++){
		if(unique[i] == 1){
			implicit_max[i] = unique_deg[i];
		}
	}

	// Check monotonicity
	vector< int > increase (this->vars.nops(), 0);
	for(int i=0; i<this->vars.nops(); i++){
		int j=0;
		while(( j< (signed)coeffs.size() ) && ( increase[i]) ){
			if( multi_index[j][i] > 0 ){
				increase[i] = increase[i] && (coeffs[j] > 0);
			}
			j++;
		}
	}
	vector< int > decrease (this->vars.nops(), 0);
	for(int i=0; i<this->vars.nops(); i++){
		int j=0;
		while(( j< (signed)coeffs.size() ) && ( decrease[i]) ){
			if( multi_index[j][i] > 0 ){
				decrease[i] = decrease[i] && (coeffs[j] < 0);
			}
			j++;
		}
	}
	for(int i=0; i<(signed)increase.size(); i++){
		if(increase[i] == 1){
			implicit_max[i] = this->degrees[i];
		}
		if(decrease[i] == 1){
			implicit_max[i] = 0;
		}
	}

	// Dominance

}

lst BaseConverter::getBernCoeffsMatrix(){

	// degrees increased by one
	vector<int> degrees_p (this->degrees.size(),0);
	for(int i=0; i<degrees_p.size(); i++){
		degrees_p[i] = this->degrees[i] + 1;
	}


	// initialize the matrix for the coefficients
	vector< ex > Ai (this->prod(degrees_p,1,degrees_p.size()),0);
	vector< vector< ex > > A (degrees_p[0],Ai);

	for(int i=0; i<this->coeffs.size(); i++){
		if(this->coeffs[i] != 0){
			vector< int > pos2d = this->n2t(this->pos2multi_index(i),degrees_p);
			A[pos2d[0]][pos2d[1]] = this->coeffs[i];
		}
	}

	vector< vector< ex > > UAt = this->transp(this->matrixProd(this->genUtilde(this->degrees[0]),A),degrees_p);
	for(int i=1; i<this->degrees.size(); i++){
		degrees_p = this->shift(degrees_p);
		UAt = this->transp(this->matrixProd(this->genUtilde(this->degrees[i]),UAt),degrees_p);
	}


	lst bernCoeffs;
	for(int i=0; i<UAt.size(); i++){
		for(int j=0; j<UAt[i].size(); j++){
			bernCoeffs.append(UAt[i][j]);
		}
	}

	return bernCoeffs;


}

// shift the vectory by one
vector<int> BaseConverter::shift(vector<int> v){

	vector<int> sv (v.size(),0);

	for(int i=1; i<v.size(); i++){
		sv[i-1] = v[i];
	}
	sv[v.size()-1] = v[0];
	return sv;

}

// multiplication of two matrices
vector< vector<ex> > BaseConverter::matrixProd(vector< vector<ex> > A, vector< vector<ex> > B){

	if( A[0].size() != B.size() ) {
		cout<<"BaseConverter::matrixProd : matrices dimensions must agree";
		exit (EXIT_FAILURE);
	}

	vector<ex> product_i (B[0].size(),0);
	vector< vector<ex> > product (A.size(),product_i);

	for(int i=0; i<A.size(); i++){
		for(int j=0; j<B[0].size(); j++){

			ex inner_prod;
			inner_prod = 0;
			for(int k=0; k<A[i].size(); k++){
				inner_prod = inner_prod + A[i][k]*B[k][j];
			}
			product[i][j] = inner_prod;

		}
	}

	return product;
}


// Operation on multi dimensional matrices

// from n-dimensions to 2
vector< int > BaseConverter::n2t( vector<int> a, vector<int> degs ){

	if( a.size() != degs.size() ){
		cout<<"BaseConverter::n2t : a and degs must have the same sizes";
		exit (EXIT_FAILURE);
	}

	vector< int > b (2,0);
	b[0] = a[0];
	b[1] = a[1];

	for (int i=2; i<a.size(); i++){
		 b[1] = b[1] + (a[i]*this->prod(degs,1,i));
	}
	return b;

}

vector< int > BaseConverter::t2n( vector<int> c, vector<int> degs ){

	vector<int> a(degs.size(),0);

	a[degs.size()-1] = floor( c[1]/this->prod(degs,1,degs.size()-1) );
	for(int i=degs.size()-1; i>0; i--){
		int div = this->prod(degs,1,i);
		a[i] = floor( c[1]/div );
		c[1] = c[1]%div;
	}
	a[0] = c[0];

	return a;
}

// transpose a 2d coordinate
vector< int > BaseConverter::transp( vector<int> b, vector<int> degs, int degs_prod ){

	vector< int > b_transp(2,0);

	b_transp[0] = b[1] % degs[1];
	b_transp[1] = ((b[1] - ( b[1] % degs[1] ))/degs[1]) + b[0]*degs_prod;

	return b_transp;

}

// transpose a 2d coordinate
vector< int > BaseConverter::transp_naive( vector<int> b, vector<int> degs ){

	return this->n2t(this->shift(this->t2n(b,degs)),this->shift(degs));

}


// transpose a 2d matrix
vector< vector< ex > > BaseConverter::transp( vector< vector<ex> > M, vector<int> degs ){

	int prod_degs2n = this->prod(degs,2,degs.size());

	int rows_t = degs[1];
	int cols_t = this->prod(degs,2,degs.size())*degs[0];

	vector< ex > M_transp_i(cols_t,0);
	vector< vector< ex > > M_transp(rows_t,M_transp_i);

	for(int i=0; i<M.size(); i++){
		for(int j=0; j<M[i].size(); j++){
			if(M[i][j] != 0){
				vector< int > ij(2,0);
				ij[0] = i; ij[1] = j;
				//vector< int > ij_t = this->transp_naive(ij,degs);
				vector< int > ij_t = this->transp(ij,degs,prod_degs2n);
				M_transp[ij_t[0]][ij_t[1]] = M[i][j];
			}
		}
	}

	return M_transp;

}

int BaseConverter::prod( vector<int> v, int a, int b ){
	int prod = 1;
	for(int i=a; i<b; i++){
		prod = prod * v[i];
	}
	return prod;

}

int BaseConverter::nchoosek(int n, int k){

	if( k == 0 ){
		return 1;
	}else{
	    return (n * this->nchoosek(n - 1, k - 1)) / k;
	}

}

// generate the U tilde matrix
vector< vector< ex > > BaseConverter::genUtilde(int n){

	vector< ex > Ui(n+1,0);
	vector< vector< ex > > U (n+1,Ui);

	for(int i=0; i<n+1; i++){
		U[i][0] = 1;
		U[n][i] = 1;
	}

	for(int i=1; i<n; i++){
		for(int j=1; j<=i; j++){
			U[i][j] = (double)this->nchoosek(i,i-j)/(double)this->nchoosek(n,j);
		}
	}
	return U;
}





void BaseConverter::print(){

	for(int i=0; i<(signed)this->coeffs.size(); i++){
		if(this->coeffs[i] != 0){
		cout<<this->coeffs[i]<<": ";
		vector<int> multi_index = pos2multi_index(i);
		for(int j=0; j<multi_index.size(); j++){
			cout<<multi_index[j]<<" ";
		}
		cout<<"\n";
		}
	}

}

void BaseConverter::print( vector< vector< ex > > M){

	for(int i=0; i<M.size(); i++){
		for(int j=0; j<M[i].size(); j++){
			cout<<M[i][j]<<" ";
		}
		cout<<"\n";
	}

}


BaseConverter::~BaseConverter() {
	// TODO Auto-generated destructor stub
}

