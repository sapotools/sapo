/*
 * MultiParallelotope.cpp
 *
 *  Created on: Oct 30, 2014
 *      Author: Tommaso Dreossi
 */

#include "MultiParallelotope.h"

MultiParallelotope::MultiParallelotope(lst qs, lst as, lst bs, vector< lst > us) {

	if(qs.nops() != as.nops() || qs.nops() != bs.nops() || as.nops() != bs.nops()){
		cout<<"MultiParallelotope::MultiParallelotope : qs, as, and bs must have the same dimensions";
		exit (EXIT_FAILURE);
	}

	this->qs = qs;
	this->as = as;
	this->bs = bs;
	this->us = us;

	// compute the generator function for a single parallelotope
	for( int i=0 ; i<(signed)qs.nops(); i++ ){
		ex gen_fun = this->qs[i];
		for( int j=0 ; j<(signed)qs.nops(); j++ ){
			gen_fun = gen_fun + this->us[i][j]*this->bs[j]*this->as[j];
		}
		this->generator_funs.append(gen_fun);
	}

	this->initAngles();

}

MultiParallelotope::MultiParallelotope(vector< vector< double > > directions, vector< double > upper_bounds, vector< double > lower_bounds, vector< vector< int > > templates) {

	if(directions.size() <= 0){
		cout<<"MultiParallelotope::MultiParallelotope : must specify at least one direction";
		exit (EXIT_FAILURE);
	}

	int n = directions[0].size();
	for( int i=0; i<(signed)directions.size(); i++ ){
		if( (signed)directions[i].size() != n ){
			cout<<"MultiParallelotope::MultiParallelotope : directions must have the same dimensions";
			exit (EXIT_FAILURE);
		}
	}

	if( upper_bounds.size() != directions.size() ){
		cout<<"MultiParallelotope::MultiParallelotope : upper_bounds and directions must have size";
		exit (EXIT_FAILURE);
	}
	if( lower_bounds.size() != directions.size() ){
		cout<<"MultiParallelotope::MultiParallelotope : lower_bounds and directions must have size";
		exit (EXIT_FAILURE);
	}

	if(templates.size() <= 0){
		cout<<"MultiParallelotope::MultiParallelotope : must specify at least one template";
		exit (EXIT_FAILURE);
	}

	for( int i=0; i<(signed)templates.size(); i++ ){
		if( (signed)templates[i].size() != n ){
			cout<<"MultiParallelotope::MultiParallelotope : templates must have the same dimensions";
			exit (EXIT_FAILURE);
		}
	}

	this->directions = directions;
	this->upper_bounds = upper_bounds;
	this->lower_bounds = lower_bounds;
	this->templates = templates;

	this->initAngles();

}

// get the p-th parallelotope of the set
LinearSystem* MultiParallelotope::getParallelotope(int p){

	if( p<0 || p>= (signed)this->templates.size() ){
		cout<<"MultiParallelotope::getParallelotope : p must must be between 0 and "<<templates.size();
		exit (EXIT_FAILURE);
	}

	vector< vector< double > > A;	// template of the parallelotope
	vector< double > b;	// offsets of the parallelotope

	for( int i=0; i<(signed)this->templates[p].size(); i++ ){
		A.push_back(this->directions[this->templates[p][i]]);
		b.push_back(this->upper_bounds[this->templates[p][i]]);
	}
	for( int i=0; i<(signed)this->templates[p].size(); i++ ){
		A.push_back(this->negate(this->directions[this->templates[p][i]]));
		b.push_back(this->lower_bounds[this->templates[p][i]]);
	}

	return new LinearSystem(A,b);

}

Parallelotope* MultiParallelotope::getPar(int p){
	if( p<0 || p>= (signed)this->templates.size() ){
			cout<<"MultiParallelotope::getParallelotopeGenFun : p must must be between 0 and "<<templates.size();
			exit (EXIT_FAILURE);
		}
		LinearSystem *LS = this->getParallelotope(p);	// get p-th parallelotope
		vector< vector< double > > A = LS->getA();	// extract directions of parallelotope
		vector<lst> set_vars;
		set_vars.push_back(this->qs); set_vars.push_back(this->as); set_vars.push_back(this->bs);
		Parallelotope *P = new Parallelotope(set_vars,LS);
		return P;
}


lst MultiParallelotope::getParallelotopeGenFun(int p){
	if( p<0 || p>= (signed)this->templates.size() ){
		cout<<"MultiParallelotope::getParallelotopeGenFun : p must must be between 0 and "<<templates.size();
		exit (EXIT_FAILURE);
	}
	LinearSystem *LS = this->getParallelotope(p);	// get p-th parallelotope
	vector< vector< double > > A = LS->getA();	// extract directions of parallelotope
	vector<lst> set_vars;
	set_vars.push_back(this->qs); set_vars.push_back(this->as); set_vars.push_back(this->bs);
	Parallelotope *P = new Parallelotope(set_vars,LS);
	return P->getGeneratorFunction();
}

// change sign to all elements of a vector
vector< double > MultiParallelotope::negate(vector< double > v){

	vector< double > minus_v;
	for(int i=0; i<v.size(); i++){
		minus_v.push_back(-v[i]);
	}
	return minus_v;
}

// get the i-th direction of the decomposition
vector< double > MultiParallelotope::getDirection(int i){
	if(( i < 0 ) || ( i >= this->getNumDirs()) ){
		cout<<"MultiParallelotope::getDirection : i must must be between 0 and "<<this->getNumDirs();
		exit (EXIT_FAILURE);
	}
	return this->directions[i];
}

// get the i-th template of the decomposition
vector< int > MultiParallelotope::getTemplate(int i){
	if(( i < 0 ) || ( i >= this->getCard()) ){
		cout<<"MultiParallelotope::getTemplate : i must must be between 0 and "<<this->getCard();
		exit (EXIT_FAILURE);
	}
	return this->templates[i];
}

// get the cardinality of the decomposition
int MultiParallelotope::getCard(){
	return this->templates.size();
}

// get the number of directions of the decomposition
int MultiParallelotope::getNumDirs(){
	return this->directions.size();
}

// couple the directions of the decomposition
// wrt the proximity of 90 degs
vector< vector< int > > MultiParallelotope::coupleDirectionsAngle(vector< vector< double > > dirs){


	vector< vector< int > > T;						// template to populate
	vector<bool> marks(dirs.size(),false);	// marks on the directions
	int dirs_dim = dirs[0].size();

	// uni-dimensional system
	// each direction is a good template
	if(dirs_dim <= 1){
		for(int i=0; i<(signed)dirs.size(); i++){
			vector<int> Lambda (1,i);
			T.push_back(Lambda);
		}
		return T;
	}


	// more than one dimension
	vector<double> angles;					// angles between directions
	vector< vector< int > > angles_idx;		// index of the directions

	// compute angles of vectors
	vector< int > angles_idx_i (2,0);
	for(int i=0; i<(signed)dirs.size(); i++){
		for(int j=i+1; j<(signed)dirs.size(); j++){
			// calulcate proximity to 90 degs
			angles.push_back(abs(this->angle(dirs[i],dirs[j])-1.5708));
			angles_idx_i[0] = i; angles_idx_i[1] = j;
			angles_idx.push_back(angles_idx_i);
		}
	}

	//for(int i=0; i<angles.size(); i++){
	//	cout<<"["<<angles_idx[i][0]<<","<<angles_idx[i][1]<<"]"<<angles[i]<<"\n";
	//}


	while( !(this->allMarked(marks)) ){

		vector<int> lambda;

		// pick the best first couple
		pair<int,int> best_first_couple = this->minAngleBothNotMarked(angles,angles_idx,marks);
		// store and mark it
		lambda.push_back(best_first_couple.first);
		lambda.push_back(best_first_couple.second);
		marks[best_first_couple.first] = true;
		marks[best_first_couple.second] = true;

//		for(int j=0; j<(signed)marks.size(); j++){
//					cout<<marks[j]<<" ";
//				}
//		cout<<"\n";

		while( (signed)lambda.size() < dirs_dim ){

			//special case where all become marked
			if(this->allMarked(marks)){
				lambda = this->fillLambda( angles, angles_idx, lambda, dirs_dim);
			}else{
				int best_new_dir = this->minAngleBothNotMarkedWith(angles,angles_idx,marks,lambda);
				lambda.push_back(best_new_dir);
				marks[best_new_dir] = true;
			}
		}
		T.push_back(lambda);
	}

	//cout<<"\nTemplate:\n";
	//for(int i=0; i<(signed)T.size(); i++){
	//	for(int j=0; j<(signed)T[i].size(); j++){
	//		cout<<T[i][j]<<" ";
	//	}
	//	cout<<"\n";
	//}
	return T;
}

void MultiParallelotope::initAngles(){

	for(int i=0; i<(signed)this->directions.size(); i++){
		vector<double> angles_i;
		for(int j=0; j<(signed)this->directions.size(); j++){
			// calculate proximity to 90 degs
			angles_i.push_back(abs(this->angle(this->directions[i],this->directions[j])-1.5708));
		}
		this->angles.push_back(angles_i);
	}

}

// shrink current set of parallelotopes around Hb
pair< vector< double >, vector< double > > MultiParallelotope::shrink(LinearSystem *Hb){

	cout<<"Shrinking...\n";

	pair< vector< double >, vector< double > > bounds;
	vector< double > upper_bounds, lower_bounds;

	for(int i=0; i<this->getNumDirs(); i++){

		// push facets towards Hb
		double b_plus = Hb->maxLinearSystem(this->getDirection(i));
		double b_minus = Hb->maxLinearSystem(this->negate(this->getDirection(i)));
		//cout<<b_plus<<" "<<b_minus<<"\n";

		upper_bounds.push_back(b_plus);
		lower_bounds.push_back(b_minus);

	}

	bounds.first = upper_bounds;
	bounds.second = lower_bounds;

	return bounds;
}

// decompose the current set wrt to H and the given directions
void MultiParallelotope::decompose(LinearSystem *H, vector< vector< double > > directions){

	cout<<"Decomposing set...\n";
	// couple directions and generate the templates
	vector< vector<int> > templates = this->coupleDirectionsAngle(directions);
	this->directions = directions;
	this->templates = templates;
	pair< vector< double >, vector< double > > bounds = this->shrink(H);
	this->upper_bounds = bounds.first;
	this->lower_bounds = bounds.second;

}

// decompose the polytope described by hp_idx and offsets
void MultiParallelotope::decompose(vector<int> hp_idx, vector<bool> offsets){

	int card_dec = ceil(hp_idx.size()/this->dim);
	double alpha = 0.5;
	vector<double> delta (card_dec,1);




}

// compute the intersection of all the parallelotopes
LinearSystem* MultiParallelotope::parallelotopesIntersection(){

	LinearSystem *intersection = new LinearSystem();

	for(int i=0; i<this->getCard(); i++){
		intersection = intersection->appendLinearSystem(this->getParallelotope(i));
	}
	return intersection;

}

// select some directions of a template
pair< vector<int>,  vector<bool> > MultiParallelotope::pickNdirections(vector< vector<double> > directions, vector<bool> marks){
	pair< vector<int>,  vector<bool> > result;
	return result;
}


// check if all the elements of v are true
bool MultiParallelotope::allMarked(vector<bool> v){
	for(int i=0; i<(signed)v.size(); i++){
		if( v[i] == 0 ){ return 0; }
	}
	return 1;
}

// return the couple of indexes with minimum angle and not both coordinates marked
pair< int,int > MultiParallelotope::minAngleBothNotMarked( vector<double> angles, vector< vector< int > > angle_idx, vector<bool> marks){

	double min_angle = DBL_MAX;
	pair< int,int > best_couple;

	for(int i=0; i<(signed)angles.size(); i++){
		if( (!marks[angle_idx[i][0]] || !marks[angle_idx[i][1]]) && angles[i] < min_angle){
			min_angle = angles[i];
			best_couple.first = angle_idx[i][0];
			best_couple.second = angle_idx[i][1];
		}
	}
	return best_couple;
}

// return the index with minimum angle with elements appearing in tmp_lambda
int MultiParallelotope::minAngleBothNotMarkedWith( vector<double> angles, vector< vector< int > > angle_idx, vector<bool> marks, vector<int> tmp_lambda){

	double min_angle = DBL_MAX;
	int best_idx;

	for(int i=0; i<(signed)angles.size(); i++){
		if(this->isIn(angle_idx[i][0],tmp_lambda) && !marks[angle_idx[i][1]] && angles[i] < min_angle){
			min_angle = angles[i];
			best_idx = angle_idx[i][1];
		}else if(this->isIn(angle_idx[i][1],tmp_lambda) && !marks[angle_idx[i][0]] && angles[i] < min_angle){
			min_angle = angles[i];
			best_idx = angle_idx[i][0];
		}
	}

	return best_idx;
}

vector<int> MultiParallelotope::fillLambda( vector<double> angles, vector< vector< int > > angle_idx, vector<int> tmp_lambda, int dim_lambda){

	while( tmp_lambda.size() < dim_lambda ){

		double min_angle = DBL_MAX;
		int best_idx;

		for(int i=0; i<angles.size(); i++){
			if( angles[i] < min_angle && !this->isIn(angle_idx[i][0],tmp_lambda)){
				min_angle = angles[i];
				best_idx = angle_idx[i][0];
			}else if( angles[i] < min_angle && !this->isIn(angle_idx[i][1],tmp_lambda)){
					min_angle = angles[i];
					best_idx = angle_idx[i][1];
			}
		}
		tmp_lambda.push_back(best_idx);
	}

	return tmp_lambda;

}

void MultiParallelotope::setUpperBounds(vector< double > upper_bounds){
	if(upper_bounds.size() == this->upper_bounds.size() ){
		this->upper_bounds = upper_bounds;
		return;
	}
	cout<<"MultiParallelotope::setUpperBounds : upper_bounds must have dimension"<<this->upper_bounds.size();
	exit (EXIT_FAILURE);
}
void MultiParallelotope::setLowerBounds(vector< double > lower_bounds){
	if(lower_bounds.size() == this->lower_bounds.size() ){
		this->lower_bounds = lower_bounds;
		return;
	}
	cout<<"MultiParallelotope::setLowerBounds : lower_bounds must have dimension"<<this->lower_bounds.size();
	exit (EXIT_FAILURE);
}
void MultiParallelotope::setBounds(vector< double > upper_bounds, vector< double > lower_bounds){

	this->setUpperBounds(upper_bounds);
	this->setLowerBounds(lower_bounds);

}

double MultiParallelotope::norm(vector<double> v){
	double sum = 0;
	for(int i=0; i<(signed)v.size(); i++){
		sum = sum + v[i]*v[i];
	}
	return sqrt(sum);
}

double MultiParallelotope::prod(vector<double> v1, vector<double> v2){
	double prod = 0;
	for(int i=0; i<(signed)v1.size(); i++){
		prod = prod + v1[i]*v2[i];
	}
	return prod;
}

double MultiParallelotope::angle(vector<double> v1, vector<double> v2){
	return acos(this->prod(v1,v2)/(this->norm(v1)*this->norm(v2)));
}

bool MultiParallelotope::isIn(int e, vector<int> v){
	for(int i=0; i<(signed)v.size(); i++){
		if( e == v[i] ){
			return true;
		}
	}
	return false;
}




void MultiParallelotope::print() {

	cout<<"b+ : ";
	for(int i=0; i<this->upper_bounds.size(); i++){
		cout<< this->upper_bounds[i]<<" ";
	}

	cout<<"\nb- : ";
	for(int i=0; i<this->lower_bounds.size(); i++){
		cout<< this->lower_bounds[i]<<" ";
	}

	cout<<"\n----------------------\n";

	for(int i=0; i<(signed)this->directions[0].size(); i++){
		for(int j=0; j<(signed)this->directions.size(); j++){
			cout<<this->directions[j][i]<<" ";
		}
		cout<<"\n";
	}
	cout<<"----------------------\n";

	for(int i=0; i<(signed)this->templates[0].size(); i++){
		for(int j=0; j<(signed)this->templates.size(); j++){
			cout<<this->templates[j][i]<<" ";
		}
		cout<<"\n";
	}

}

// print all the parallelotopes
void MultiParallelotope::printParallelotopes(){
	for(int i=0; i<this->getCard(); i++){
		this->getParallelotope(i)->plotRegion();
	}
}

// print all the parallelotopes intersection
void MultiParallelotope::printParallelotopesIntersection(){
	this->parallelotopesIntersection()->print();
}


MultiParallelotope::~MultiParallelotope() {
	// TODO Auto-generated destructor stub
}

