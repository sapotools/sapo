/*
 *
 *  Created on: Sep 18, 2015
 *      Author: dreossi
 */

#include "Bundle.h"
#include <string>


// initial bundle to decompose
Bundle::Bundle(vector<lst> vars, vector< vector< double > > L, vector< double > offp, vector< double > offm, 	vector< vector< int > > T) {

	if( L.size() > 0 ){
		this->dim = L[0].size();
	}else{
		cout<<"Bundle::Bundle : L must be non empty";
	}
	if( L.size() != offp.size() ){
		cout<<"Bundle::Bundle : L and offp must have the same size";
		exit (EXIT_FAILURE);
	}
	if( L.size() != offm.size() ){
		cout<<"Bundle::Bundle : L and offm must have the same size";
		exit (EXIT_FAILURE);
	}
	if( T.size() > 0 ){
		for(int i=0; i<T.size(); i++){
			if( T[i].size() != this->getDim() ){
				cout<<"Bundle::Bundle : T must have "<<this->getDim()<<" columns";
				exit (EXIT_FAILURE);
			}
		}
	}else{
		cout<<"Bundle::Bundle : T must be non empty";
		exit (EXIT_FAILURE);
	}

	this->vars = vars;
	this->L = L;
	this->offp = offp;
	this->offm = offm;
	this->T = T;
	this->coeffMap = new CoeffMap();

	// initialize orthogonal proximity
	for(int i=0; i<this->getNumDirs(); i++){
		vector< double > Thetai (this->getNumDirs(),0);
		for(int j=i; j<this->getNumDirs(); j++){
			this->Theta.push_back(Thetai);
		}
	}
	for(int i=0; i<this->getNumDirs(); i++){
		this->Theta[i][i] = 0;
		for(int j=i+1; j<this->getNumDirs(); j++){
			double prox = this->orthProx(this->L[i],this->L[j]);
			this->Theta[i][j] = prox;
			this->Theta[j][i] = prox;
		}
	}
}

// construct the bundle
LinearSystem* Bundle::getBundle(){
	vector< vector< double> > A;
	vector< double> b;
	for(int i=0; i<this->getSize(); i++){
		A.push_back(this->L[i]);
		b.push_back(this->offp[i]);
	}
	for(int i=0; i<this->getSize(); i++){
		A.push_back(this->negate(this->L[i]));
		b.push_back(this->offm[i]);
	}

	LinearSystem* Ab = new LinearSystem(A,b);

	return Ab;
}

Parallelotope* Bundle::getParallelotope(int i){

	if( i<0 || i>this->T.size() ){
		cout<<"Bundle::getParallelotope : i must be between 0 and "<<T.size();
		exit (EXIT_FAILURE);
	}

	vector<double> d;
	vector< vector< double > > Lambda;

	// upper facets
	for(int j=0; j<this->getDim(); j++){
		Lambda.push_back(this->L[this->T[i][j]]);
		d.push_back(this->offp[this->T[i][j]]);
	}
	// lower facets
	for(int j=0; j<this->getDim(); j++){
		Lambda.push_back(this->negate(this->L[this->T[i][j]]));
		d.push_back(this->offm[this->T[i][j]]);
	}

	LinearSystem *Lambdad = new LinearSystem(Lambda,d);
	Parallelotope *P = new Parallelotope(this->vars, Lambdad);


	return P;

}

// Canonize the current Bundle
void Bundle::canonize(){
	// get current polytope
	LinearSystem *bund = this->getBundle();
	for(int i=0; i<this->getSize(); i++){
		this->offp[i] = bund->maxLinearSystem(this->L[i]);
		this->offm[i] = bund->maxLinearSystem(this->negate(this->L[i]));
	}
}

vector< vector<int> > Bundle::decompose(){

	double alpha = 0.5;

	vector< vector< int > > newT;
	LinearSystem *bundle = this->getBundle();
	vector<bool> redunConstr = bundle->redundantCons();

	vector< double > offDists = this->offsetDistances();		// compute current distances

	int temp_card = ceil(this->getNumDirs()*2/this->getDim());

	vector<int> templ_template_i;
	vector< vector<int> > tmp_templates;

	// initialize templates in construction
	for(int i=0; i<temp_card; i++){
		tmp_templates.push_back(templ_template_i);
	}

	vector< vector<int> > complete_templates;

	while( true ){

		for( int i=0; i<this->getNumDirs(); i++ ){

			if( !redunConstr[i] || !redunConstr[i+this->getDim()] ){

				int candidateParatope;
				double candidateParatopeW = DBL_MAX;

				// find parallelotope that minimizes the distances
				for( int j=0; j<tmp_templates.size(); j++ ){
					if( !this->isIn(i,tmp_templates[j]) ){

						double ortProx = this->maxOrthProx(i,tmp_templates[j]);
						double dist = this->maxOffsetDist(i,tmp_templates[j],offDists);
						double w = alpha*ortProx + (1-alpha)*dist;
						if( w < candidateParatopeW ){
							candidateParatope  = j;
						}
					}
				}

				if( !this->isIn(i,tmp_templates[candidateParatope]) ){
					tmp_templates[candidateParatope].push_back(i);
				}

				// check if the parallelotope is ready
				if(tmp_templates[candidateParatope].size() == this->getDim()){
					// store the complete template
					complete_templates.push_back(tmp_templates[candidateParatope]);
					// and remove it from the list of parallelotopes in construction
					tmp_templates.erase(tmp_templates.begin()+candidateParatope);
				}

				if(tmp_templates.empty()){
					return complete_templates;
				}
			}
		}
	}
}


vector< vector<int> > Bundle::decomposeRand(){

	int max_iters = 100;


	vector< double > offDists = this->offsetDistances();				// compute current distances
	vector<bool> redunConstr = this->getBundle()->redundantCons();		// redundant constraints

	vector< int > dirs;
	for(int i=0; i<this->getSize(); i++){
		if( !redunConstr[i] || !redunConstr[i+this->getSize()] ){
			dirs.push_back(i);
		}
	}

	//int temp_card = ceil(dirs.size()/this->getDim());
	int temp_card = this->getSize();

	// initialize a valid template
	vector< int > Ti (this->getDim(),0);
	vector< vector< int > > T (temp_card,Ti);

	int dirs_i = 0;
	for(int i=0; i<T.size(); i++){
		for(int j=0; j<T[i].size(); j++){
			T[i][j] = dirs[dirs_i%dirs.size()];
			dirs_i++;
		}
	}

	bool goodTemp = true;
	// we have to check that the initialized template is non-signlaur
	for(int i=0; i<T.size(); i++){
		// check if it's non singluar
		ex eq = 0;
		lst LS;
		for(int j=0; j<this->getDim(); j++){
			for(int k=0; k<this->getDim(); k++){
				eq = eq + this->vars[0][k]*this->L[T[i][j]][k];
			}
			LS.append( eq == this->offp[j] );
		}
		ex solLS = lsolve(LS,this->vars[0]);
		goodTemp = goodTemp && solLS.nops() != 0 ;
	}

	while(!goodTemp){

		dirs.clear();
		for(int i=0; i<this->getSize(); i++){
			dirs.push_back(i);
		}
		random_shuffle ( dirs.begin(), dirs.end() );

		temp_card = ceil(this->getSize()/this->getDim());

		random_shuffle ( dirs.begin(), dirs.end() );

		temp_card = ceil(this->getSize()/this->getDim());
		vector< int > newTi (this->getDim(),0);
		vector< vector< int > > newT (temp_card,newTi);
		dirs_i = 0;
		for(int i=0; i<newT.size(); i++){
			for(int j=0; j<newT[i].size(); j++){
				newT[i][j] = dirs_i%this->getSize();
				dirs_i++;
			}
		}

		T = newT;

		for(int i=0; i<T.size(); i++){
				// check if it's non singluar
				ex eq = 0;
				lst LS;
				for(int j=0; j<this->getDim(); j++){
					for(int k=0; k<this->getDim(); k++){
						eq = eq + this->vars[0][k]*this->L[T[i][j]][k];
					}
					LS.append( eq == this->offp[j] );
				}
				ex solLS = lsolve(LS,this->vars[0]);
				goodTemp = solLS.nops() != 0 ;
			}

	}

//	cout<<"Tprima\n";
//	for(int i=0; i<T.size(); i++){
//		for(int j=0; j<T[i].size(); j++){
//			cout<<T[i][j]<<" ";
//		}
//		cout<<"\n";
//	}


	// try to improve the acutal template
	int i=0;
	while( i<max_iters ){

		vector< vector< int > > tmpT = T;

		// generate random coordinates to swap
		int i1 = rand() % temp_card;
		int i2 = rand() % temp_card;
		int j1 = rand() % this->getDim();
		int j2 = rand() % this->getDim();

		// swap them
		if( ( i1 != i2 ) && ( j1 != j2 ) ){
			int d1 = tmpT[i1][j1];
			int d2 = tmpT[i2][j2];
			tmpT[i1][j1] = d2;
			tmpT[i2][j2] = d1;
		}

		// check if it's non singluar
		ex eq1 = 0;
		lst LS1;
		for(int j=0; j<this->getDim(); j++){
			for(int k=0; k<this->getDim(); k++){
				eq1 = eq1 + this->vars[0][k]*this->L[tmpT[i1][j]][k];
			}
			LS1.append( eq1 == this->offp[j] );
		}
		ex solLS1 = lsolve(LS1,this->vars[0]);
		ex eq2 = 0;
		lst LS2;
		for(int j=0; j<this->getDim(); j++){
			for(int k=0; k<this->getDim(); k++){
				eq2 = eq2 + this->vars[0][k]*this->L[tmpT[i2][j]][k];
			}
			LS2.append( eq2 == this->offp[j] );
		}
		ex solLS2 = lsolve(LS2,this->vars[0]);

//		if( tmpT[i1][0] == 0 && tmpT[i1][1] == 1 && tmpT[i1][2] == 3){
//			cout<<LS1;
//			cout<<"sol1 "<<solLS1.nops()<<"\n";
//		}
//		if( tmpT[i2][0] == 0 && tmpT[i2][1] == 1 && tmpT[i2][2] == 3){
//			cout<<LS2;
//			cout<<"sol2 "<<solLS2.nops()<<"\n";
//		}

		double alpha = 1;
		if( solLS1.nops() != 0 && solLS2.nops() != 0 ){

			double w1 = alpha*this->maxOffsetDist(tmpT,offDists) + (1-alpha)*this->maxOrthProx(tmpT);
			double w2 = alpha*this->maxOffsetDist(T,offDists) + (1-alpha)*this->maxOrthProx(T);

			if( w1 < w2 ){
				T = tmpT;
			}
		}
		i++;
	}

//	cout<<"Tdopop\n";
//	for(int i=0; i<T.size(); i++){
//		for(int j=0; j<T[i].size(); j++){
//			cout<<T[i][j]<<" ";
//		}
//		cout<<"\n";
//	}

	return T;

}

vector< vector<int> > Bundle::decomposeTotalRand(){

	vector< double > offDists = this->offsetDistances();

//	for(int i=0; i<offDists.size(); i++){
//		cout<<offDists[i]<<" ";
//	}
//	cout<<"\n";

	vector< vector<int> > curT = this->T;		// get actual template and try to improve it
	vector< vector<int> > bestT = this->T;		// get actual template and try to improve it
	int temp_card = this->T.size();
	int max_iters = 500;

	//cout<<"w entrata:"<<this->maxOffsetDist(this->T,offDists)<<"\n";

	int i=0;
	while( i<max_iters ){

		vector< vector<int> > tmpT = curT;

		// generate random coordinates to swap
		int i1 = rand() % temp_card;
		int j1 = rand() % this->getDim();
		int new_element = rand() % this->getSize();

		// swap them
		tmpT[i1][j1] = new_element;

		bool valid = true;
		// check for duplicates
		vector<int> newTemp1 = tmpT[i1];
		for(int j=0; j<tmpT.size(); j++){
			if( j != i1 ){
				valid = valid && !(this->isPermutation(newTemp1,tmpT[j]));
			}
		}

		if(valid){
			ex eq1 = 0;
			lst LS1;
			for(int j=0; j<this->getDim(); j++){
				for(int k=0; k<this->getDim(); k++){
					eq1 = eq1 + this->vars[0][k]*this->L[tmpT[i1][j]][k];
				}
				LS1.append( eq1 == this->offp[j] );
			}
			ex solLS1 = lsolve(LS1,this->vars[0]);

			if( solLS1.nops() != 0 ){

				double alpha = 0;
				double w1 = alpha*this->maxOffsetDist(tmpT,offDists) + (1-alpha)*this->maxOrthProx(tmpT);
				double w2 = alpha*this->maxOffsetDist(bestT,offDists) + (1-alpha)*this->maxOrthProx(bestT);

				if( w1 < w2 ){
					bestT = tmpT;
				}
				curT = tmpT;
			}
		}
		i++;
	}

	return bestT;

}


// compute the transformation of the bundle
pair< vector<double>, vector<double> > Bundle::transform(lst vars, lst f, bool mode){

	vector<double> newDp (this->getSize(),DBL_MAX);
	vector<double> newDm (this->getSize(),DBL_MAX);

	vector<int> dirs_to_bound;
	if(mode){	// dynamic transformation
		for(int i=0; i<(signed)this->L.size(); i++){
			dirs_to_bound.push_back(i);
		}
	}

	for(int i=0; i<this->getCard(); i++){	// for each parallelotope


		Parallelotope *P = this->getParallelotope(i);
		lst genFun = P->getGeneratorFunction();

		vector< double > base_vertex = P->getBaseVertex();
		vector< double > lengths = P->getLenghts();

		lst subParatope;

		for(int k=0; k<(signed)this->vars[0].nops(); k++){
			subParatope.append(this->vars[0][k] == base_vertex[k]);
			subParatope.append(this->vars[2][k] == lengths[k]);
		}


		if(!mode){	// static mode
			dirs_to_bound = this->T[i];
		}

		for(int j=0; j<(signed)dirs_to_bound.size(); j++){	// for each direction

			vector<int> keyp = this->T[i]; // direction indexes of the actual parallelotope
			keyp.push_back(dirs_to_bound[j]);
			vector<int> keym = this->T[i]; // direction indexes of the actual parallelotope
			keym.push_back(-(dirs_to_bound[j]+1));


			lst bernCoeffsp,bernCoeffsm;

			//if( this->bernCoeffs.find(keyp) == this->bernCoeffs.end() ){

				// the combination parallelotope/direction to bound is not present in hash table
				// compute control points
				lst sub, fog;

				for(int k=0; k<(signed)vars.nops(); k++){
					sub.append(vars[k] == genFun[k]);
				}
				for(int k=0; k<(signed)vars.nops(); k++){
					fog.append(f[k].subs(sub));
				}

				ex Lfogp; Lfogp = 0;
				// upper facets
				for(int k=0; k<this->getDim(); k++){
					Lfogp = Lfogp + this->L[dirs_to_bound[j]][k]*fog[k];
				}

				BaseConverter *BCp = new BaseConverter(this->vars[1],Lfogp);
				bernCoeffsp = BCp->getBernCoeffsMatrix();
				this->bernCoeffs[keyp] = bernCoeffsp;
			//}else{
			//	bernCoeffsp = this->bernCoeffs.find(keyp)->second;
			//}

		//if( this->bernCoeffs.find(keym) == this->bernCoeffs.end() ){

//			lst sub, fog;
//
//			for(int k=0; k<(signed)vars.nops(); k++){
//				sub.append(vars[k] == genFun[k]);
//			}
//			for(int k=0; k<(signed)vars.nops(); k++){
//				fog.append(f[k].subs(sub));
//			}


				ex Lfogm; Lfogm = 0;
				// upper facets
				for(int k=0; k<this->getDim(); k++){
					Lfogm = Lfogm - this->L[dirs_to_bound[j]][k]*fog[k];
				}

				BaseConverter *BCm = new BaseConverter(this->vars[1],Lfogm);
				bernCoeffsm = BCm->getBernCoeffsMatrix();

				this->bernCoeffs[keym] = bernCoeffsm;
		//}else{
		//	bernCoeffsm = this->bernCoeffs.find(keym)->second;
		//}

		double maxCoeff = -DBL_MAX;
		for(int k=0; k<bernCoeffsp.nops(); k++){
			double actCoeff = ex_to<numeric>(evalf(bernCoeffsp[k].subs(subParatope))).to_double();
			maxCoeff = max(maxCoeff,actCoeff);
		}
		newDp[dirs_to_bound[j]] = min(newDp[dirs_to_bound[j]],maxCoeff);


		maxCoeff = -DBL_MAX;
		for(int k=0; k<bernCoeffsm.nops(); k++){
			double actCoeff = ex_to<numeric>(evalf(bernCoeffsm[k].subs(subParatope))).to_double();
			maxCoeff = max(maxCoeff,actCoeff);
		}
		newDm[dirs_to_bound[j]]  = min(newDm[dirs_to_bound[j]],maxCoeff);
		}
	}

	if(!mode){
		this->canonize();
	}

	pair< vector<double>, vector<double> > newOffsets (newDp,newDm);
	return newOffsets;
}


// update the templates matrix
void Bundle::setTemplate(vector< vector< int > > T){
	this->T = T;
}


// compute the distances between parallel half spaces
vector< double > Bundle::offsetDistances(){

	vector< double > dist;
	for(int i=0; i<this->getSize(); i++){
		dist.push_back( abs(this->offp[i] - this->offm[i]) / this->norm(this->L[i]) );
	}
	return dist;

}


double Bundle::norm(vector<double> v){
	double sum = 0;
	for(int i=0; i<(signed)v.size(); i++){
		sum = sum + v[i]*v[i];
	}
	return sqrt(sum);
}

double Bundle::prod(vector<double> v1, vector<double> v2){
	double prod = 0;
	for(int i=0; i<(signed)v1.size(); i++){
		prod = prod + v1[i]*v2[i];
	}
	return prod;
}

double Bundle::angle(vector<double> v1, vector<double> v2){
	return acos(this->prod(v1,v2)/(this->norm(v1)*this->norm(v2)));
}

// orthogonal proximity
double Bundle::orthProx(vector<double> v1, vector<double> v2){
	return abs(this->angle(v1,v2) - (3.14159265/2));
}

// max orthogonal proximity wrt p_set
double Bundle::maxOrthProx(int vIdx, vector<int> dirsIdx){

	if(dirsIdx.empty()){
		return 0;
	}

	double maxProx = 0;
	for( int i=0; i<dirsIdx.size(); i++ ){
		maxProx = max(maxProx, this->orthProx(this->L[vIdx],this->L[dirsIdx[i]]));
	}
	return maxProx;
}

double Bundle::maxOrthProx(vector<int> dirsIdx){
	double maxProx = 0;
	for( int i=0; i<dirsIdx.size(); i++ ){
		for(int j=i+1; j<dirsIdx.size(); j++){
			maxProx = max(maxProx, this->orthProx(this->L[dirsIdx[i]],this->L[dirsIdx[j]]));
		}
	}
	return maxProx;
}

double Bundle::maxOrthProx(vector< vector<int> > T){
	double maxorth = -DBL_MAX;
	for(int i=0; i<T.size(); i++){
		maxorth = max(maxorth,this->maxOrthProx(T[i]));
	}
	return maxorth;
}

double Bundle::maxOffsetDist(int vIdx, vector<int> dirsIdx, vector<double> dists){

	if(dirsIdx.empty()){
		return 0;
	}

	double dist = dists[vIdx];
	for(int i=0; i<dirsIdx.size(); i++){
		dist = dist * dists[dirsIdx[i]];
	}
	return dist;

}

double Bundle::maxOffsetDist(vector<int> dirsIdx, vector<double> dists){

	double dist = 1;
	for(int i=0; i<dirsIdx.size(); i++){
		dist = dist * dists[dirsIdx[i]];
	}
	return dist;

}

double Bundle::maxOffsetDist(vector< vector<int> > T, vector<double> dists){
	double maxdist = -DBL_MAX;
	for(int i=0; i<T.size(); i++){
		maxdist = max(maxdist,this->maxOffsetDist(T[i],dists));
	}
	return maxdist;
}

// change sign to all elements of a vector
vector< double > Bundle::negate(vector< double > v){

	vector< double > minus_v;
	for(int i=0; i<v.size(); i++){
		minus_v.push_back(-v[i]);
	}
	return minus_v;
}

bool Bundle::isIn(int n, vector<int> v){

	for(int i=0; i<v.size(); i++){
		if( n == v[i] ){
			return true;
		}
	}
	return false;
}

bool Bundle::isIn(vector<int> v, vector< vector< int > > vlist){
	for(int i=0; i<vlist.size(); i++){
		if( this->isPermutation(v,vlist[i]) ){
			return true;
		}
	}
	return false;
}

// check if v1 is a permutation of v2
bool Bundle::isPermutation(vector<int> v1, vector<int> v2){
	for( int i=0; i<v1.size(); i++ ){
		if( !this->isIn(v1[i],v2) ){
			return false;
		}
	}
	return true;
}

//check if T is valid
bool Bundle::validTemp(vector< vector<int> > T, int card, vector<int> dirs){

	cout<<"dirs: ";
	for(int i=0; i<dirs.size(); i++){
		cout<<dirs[i]<<" ";
	}
	cout<<"\n";

	cout<<"T:\n";
	for(int i=0; i<T.size(); i++){
		for(int j=0; j<T[i].size(); j++){
			cout<<T[i][j]<<" ";
		}
		cout<<"\n";
	}
	cout<<"\n";

	if( T.size() != card ){
		return false;
	}

	// check if all the directions appear in T
	vector< bool > dirIn (dirs.size(), false);
	for(int i=0; i<dirs.size(); i++){
		for(int j=0; j<T.size(); j++){
			dirIn[i] = this->isIn(dirs[i],T[j]);
		}
	}
	for(int i=0; i<dirIn.size(); i++){
		if( !dirIn[i] ){
			return false;
		}
	}
//
//	// check if all the directions are non null
//	vector< bool > nonNullDir (this->getDim(),false);
//	for( int i=0; i<T.size(); i++ ){
//		for( int j=0; j<T[i].size(); j++ ){
//			for(int k=0; k<this->getDim(); k++){
//				nonNullDir[k] = this->L[T[i][j]][k] != 0;
//			}
//		}
//	}
//	for(int i=0; i<nonNullDir.size(); i++){
//		if( !nonNullDir[i] ){
//			return false;
//		}
//	}
	return true;
}

Bundle::~Bundle() {
	// TODO Auto-generated destructor stub
}

