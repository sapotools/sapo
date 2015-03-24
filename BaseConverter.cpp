/*
 * BaseConverter.cpp
 *
 *  Created on: Oct 23, 2014
 *      Author: dreossi
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


	// Initialize the degree shifts
	initShifts();

	for(int i=0;i<this->shifts[0];i++){
		this->coeffs.push_back(0);
	}

	// Initialize the coefficients vector
	vector<int> multi_index;
	extractCoeffs(this->polynomial,0,multi_index);

//	cout<<"\nPoly: "<<polynomial<<"\n";
//	cout<<"\nDegs: ";
//	for(int i=0; i<this->degrees.size(); i++){
//		cout<<this->degrees[i];
//	}
//	cout<<"\n";

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

	cout<<"Degrees: ";
	for(int i=0; i<(signed)this->vars.nops();i++){
		cout<<this->degrees[i]<<", ";
	}
	cout<<"(Total points:"<<bern_coeffs.nops()<<")\n";

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






BaseConverter::~BaseConverter() {
	// TODO Auto-generated destructor stub
}

