/**
 * @file Bundle.cpp
 * Represent and manipulate bundles of parallelotopes whose intersection
 * represents a polytope
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @version 0.1
 */

#include "Bundle.h"
#include <string>


/**
 * Constructor that instantiates the bundle
 *
 * @param[in] vars list of variables for parallelotope generator functions
 * @param[in] L matrix of directions
 * @param[in] offp upper offsets
 * @param[in] offm lower offsets
 * @param[in] T templates matrix
 */
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

/**
 * Constructor that instantiates the bundle with auto-generated variables
 *
 * @param[in] L matrix of directions
 * @param[in] offp upper offsets
 * @param[in] offm lower offsets
 * @param[in] T templates matrix
 */
Bundle::Bundle(vector< vector< double > > L, vector< double > offp, vector< double > offm, 	vector< vector< int > > T){

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

	//generate the variables
	VarsGenerator *varsGen = new VarsGenerator(T[0].size());
	lst qs, as, bs, ls;
	vector<lst> us;
	qs = varsGen->getBaseVertex();
	as = varsGen->getFreeVars();
	bs = varsGen->getLenghts();

	vector<lst> paraVars;
	paraVars.push_back(qs);
	paraVars.push_back(as);
	paraVars.push_back(bs);

	this->vars = paraVars;
	this->L = L;
	this->offp = offp;
	this->offm = offm;
	this->T = T;

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

/**
 * Generate the polytope represented by the bundle
 *
 * @returns polytope represented by a linear system
 */
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

/**
 * Get the i-th parallelotope of the bundle
 *
 * @param[in] i parallelotope index to fetch
 * @returns i-th parallelotope
 */
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

/**
 * Canonize the current bundle pushing the constraints toward the symbolic polytope
 *
 * @returns canonized bundle
 */
Bundle* Bundle::canonize(){
	// get current polytope
	LinearSystem *bund = this->getBundle();
	vector<double> canoffp,canoffm;
	for(int i=0; i<this->getSize(); i++){
		canoffp.push_back(bund->maxLinearSystem(this->L[i]));
		canoffm.push_back(bund->maxLinearSystem(this->negate(this->L[i])));
	}
	return new Bundle(this->vars,this->L,canoffp,canoffm,this->T);
}

/**
 * Decompose the current symbolic polytope
 *
 * @param[in] alpha weight parameter in [0,1] for decomposition (0 for distance, 1 for orthogonality)
 * @param[in] max_iter maximum number of randomly generated templates
 * @returns new bundle decomposing current symbolic polytope
 */
Bundle* Bundle::decompose(double alpha, int max_iters){

	vector< double > offDists = this->offsetDistances();

	vector< vector<int> > curT = this->T;		// get actual template and try to improve it
	vector< vector<int> > bestT = this->T;		// get actual template and try to improve it
	int temp_card = this->T.size();

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

	return new Bundle(this->vars,this->L,this->offp,this->offp,bestT);

}

/**
 * Transform the bundle
 *
 * @param[in] vars variables appearing in the transforming function
 * @param[in] f transforming function
 * @param[in,out] controlPts control points computed so far that might be updated
 * @param[in] mode transformation mode (0=OFO,1=AFO)
 * @returns transformed bundle
 */
Bundle* Bundle::transform(lst vars, lst f, map< vector<int>,pair<lst,lst> > &controlPts, int mode){

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

		if(mode == 0){	// static mode
			dirs_to_bound = this->T[i];
		}

		for(int j=0; j<(signed)dirs_to_bound.size(); j++){	// for each direction

			// key of the control points
			vector<int> key = this->T[i];
			key.push_back(dirs_to_bound[j]);

			lst actbernCoeffs;


			if( controlPts.count(key) == 0 || (!controlPts[key].first.is_equal(genFun)) ){	// check if the coefficients were already computed

				// the combination parallelotope/direction to bound is not present in hash table
				// compute control points
				lst sub, fog;

				for(int k=0; k<(signed)vars.nops(); k++){
					sub.append(vars[k] == genFun[k]);
				}
				for(int k=0; k<(signed)vars.nops(); k++){
					fog.append(f[k].subs(sub));
				}

				ex Lfog; Lfog = 0;
				// upper facets
				for(int k=0; k<this->getDim(); k++){
					Lfog = Lfog + this->L[dirs_to_bound[j]][k]*fog[k];
				}

				BaseConverter *BC = new BaseConverter(this->vars[1],Lfog);
				actbernCoeffs = BC->getBernCoeffsMatrix();

				pair<lst,lst> element (genFun,actbernCoeffs);
				controlPts[key] = element;	// store the computed coefficients

			}else{
				actbernCoeffs = controlPts[key].second;
			}

			// find the maximum coefficient
			double maxCoeffp = -DBL_MAX;
			double maxCoeffm = -DBL_MAX;
			
			for (lst::const_iterator c = actbernCoeffs.begin(); c != actbernCoeffs.end(); ++c)
			{
				double actCoeffp = ex_to<numeric>((*c).subs(subParatope)).to_double();
				double actCoeffm = ex_to<numeric>((-(*c)).subs(subParatope)).to_double();
				maxCoeffp = max(maxCoeffp,actCoeffp);
				maxCoeffm = max(maxCoeffm,actCoeffm);
			}
			newDp[dirs_to_bound[j]] = min(newDp[dirs_to_bound[j]],maxCoeffp);
			newDm[dirs_to_bound[j]]  = min(newDm[dirs_to_bound[j]],maxCoeffm);
		}
	}

	Bundle *res = new Bundle(this->vars,this->L,newDp,newDm,this->T);
	if(mode == 0){
		res = res->canonize();
	}

	return res;
}

/**
 * Parametric transformation of the bundle
 *
 * @param[in] vars variables appearing in the transforming function
 * @param[in] params parameters appearing in the transforming function
 * @param[in] f transforming function
 * @param[in] paraSet set of parameters
 * @param[in,out] controlPts control points computed so far that might be updated
 * @param[in] mode transformation mode (0=OFO,1=AFO)
 * @returns transformed bundle
 */
Bundle* Bundle::transform(lst vars, lst params, lst f, LinearSystem *paraSet, map< vector<int>,pair<lst,lst> > &controlPts, int mode){

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


		if(mode == 0){	// static mode
			dirs_to_bound = this->T[i];
		}

		for(int j=0; j<(signed)dirs_to_bound.size(); j++){	// for each direction

			// key of the control points
			vector<int> key = this->T[i];
			key.push_back(dirs_to_bound[j]);

			lst actbernCoeffs;


			if( controlPts.count(key) == 0 || (!controlPts[key].first.is_equal(genFun)) ){	// check if the coefficients were already computed

				// the combination parallelotope/direction to bound is not present in hash table
				// compute control points
				lst sub, fog;

				for(int k=0; k<(signed)vars.nops(); k++){
					sub.append(vars[k] == genFun[k]);
				}
				for(int k=0; k<(signed)vars.nops(); k++){
					fog.append(f[k].subs(sub));
				}

				ex Lfog; Lfog = 0;
				// upper facets
				for(int k=0; k<this->getDim(); k++){
					Lfog = Lfog + this->L[dirs_to_bound[j]][k]*fog[k];
				}

				BaseConverter *BC = new BaseConverter(this->vars[1],Lfog);
				actbernCoeffs = BC->getBernCoeffsMatrix();

				pair<lst,lst> element (genFun,actbernCoeffs);
				controlPts[key] = element;	// store the computed coefficients

			}else{
				actbernCoeffs = controlPts[key].second;
			}

			// find the maximum coefficient
			double maxCoeffp = -DBL_MAX;
			double maxCoeffm = -DBL_MAX;
		

			for (lst::const_iterator c = actbernCoeffs.begin(); c != actbernCoeffs.end(); ++c)
			{
				ex paraBernCoeff;
				paraBernCoeff = (*c).subs(subParatope);
				maxCoeffp = max(maxCoeffp,paraSet->maxLinearSystem(params,paraBernCoeff));
				maxCoeffm = max(maxCoeffm,paraSet->maxLinearSystem(params,-paraBernCoeff));
			}
			newDp[dirs_to_bound[j]] = min(newDp[dirs_to_bound[j]],maxCoeffp);
			newDm[dirs_to_bound[j]]  = min(newDm[dirs_to_bound[j]],maxCoeffm);
		}
	}

	Bundle *res = new Bundle(this->vars,this->L,newDp,newDm,this->T);
	if(mode == 0){
		res = res->canonize();
	}

	return res;
}


/**
 * Set the bundle template
 *
 * @param[in] T new template
 */
void Bundle::setTemplate(vector< vector< int > > T){
	this->T = T;
}


/**
 * Compute the distances between the half-spaced of the parallelotopes
 *
 * @returns vector of distances
 */
vector< double > Bundle::offsetDistances(){

	vector< double > dist;
	for(int i=0; i<this->getSize(); i++){
		dist.push_back( abs(this->offp[i] - this->offm[i]) / this->norm(this->L[i]) );
	}
	return dist;

}

/**
 * Compute the norm of a vector
 *
 * @param[in] v vector to normalize
 * @returns norm of the given vector
 */
double Bundle::norm(vector<double> v){
	double sum = 0;
	for(int i=0; i<(signed)v.size(); i++){
		sum = sum + v[i]*v[i];
	}
	return sqrt(sum);
}

/**
 * Compute the product of two vectors
 *
 * @param[in] v1 left vector
 * @param[in] v2 right vector
 * @returns product v1*v2'
 */
double Bundle::prod(vector<double> v1, vector<double> v2){
	double prod = 0;
	for(int i=0; i<(signed)v1.size(); i++){
		prod = prod + v1[i]*v2[i];
	}
	return prod;
}

/**
 * Compute the angle between two vectors
 *
 * @param[in] v1 vector
 * @param[in] v2 vector
 * @returns angle between v1 and v2
 */
double Bundle::angle(vector<double> v1, vector<double> v2){
	return acos(this->prod(v1,v2)/(this->norm(v1)*this->norm(v2)));
}

/**
 * Orthogonal proximity of v1 and v2, i.e.,
 * how close is the angle between v1 and v2 is to pi/2
 *
 * @param[in] v1 vector
 * @param[in] v2 vector
 * @returns orthogonal proximity
 */
double Bundle::orthProx(vector<double> v1, vector<double> v2){
	return abs(this->angle(v1,v2) - (3.14159265/2));
}

/**
 * Maximum orthogonal proximity of a vector w.r.t. a set of vectors
 *
 * @param[in] vIdx index of the reference vector
 * @param[in] dirsIdx indexes of vectors to be considered
 * @returns maximum orthogonal proximity
 */
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

/**
 * Maximum orthogonal proximity within a set of vectors
 *
 * @param[in] dirsIdx indexes of vectors to be considered
 * @returns maximum orthogonal proximity
 */
double Bundle::maxOrthProx(vector<int> dirsIdx){
	double maxProx = 0;
	for( int i=0; i<dirsIdx.size(); i++ ){
		for(int j=i+1; j<dirsIdx.size(); j++){
			maxProx = max(maxProx, this->orthProx(this->L[dirsIdx[i]],this->L[dirsIdx[j]]));
		}
	}
	return maxProx;
}

/**
 * Maximum orthogonal proximity of all the vectors of a matrix
 *
 * @param[in] T collection of vectors
 * @returns maximum orthogonal proximity
 */
double Bundle::maxOrthProx(vector< vector<int> > T){
	double maxorth = -DBL_MAX;
	for(int i=0; i<T.size(); i++){
		maxorth = max(maxorth,this->maxOrthProx(T[i]));
	}
	return maxorth;
}

/**
 * Maximum distance accumulation of a vector w.r.t. a set of vectors
 *
 * @param[in] vIdx index of the reference vector
 * @param[in] dirsIdx indexes of vectors to be considered
 * @param[in] dists pre-computed distances
 * @returns distance accumulation
 */
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

/**
 * Maximum distance accumulation of a set of vectors
 *
 * @param[in] dirsIdx indexes of vectors to be considered
 * @param[in] dists pre-computed distances
 * @returns distance accumulation
 */
double Bundle::maxOffsetDist(vector<int> dirsIdx, vector<double> dists){

	double dist = 1;
	for(int i=0; i<dirsIdx.size(); i++){
		dist = dist * dists[dirsIdx[i]];
	}
	return dist;

}

/**
 * Maximum distance accumulation of matrix
 *
 * @param[in] T matrix from which fetch the vectors
 * @param[in] dists pre-computed distances
 * @returns distance accumulation
 */
double Bundle::maxOffsetDist(vector< vector<int> > T, vector<double> dists){
	double maxdist = -DBL_MAX;
	for(int i=0; i<T.size(); i++){
		maxdist = max(maxdist,this->maxOffsetDist(T[i],dists));
	}
	return maxdist;
}

/**
 * Change sign to all elements of a vector
 *
 * @param[in] v vector to reverse
 * @returns reversed vector
 */
vector< double > Bundle::negate(vector< double > v){

	vector< double > minus_v;
	for(int i=0; i<v.size(); i++){
		minus_v.push_back(-v[i]);
	}
	return minus_v;
}

/**
 * Determine belonging of an element in a vector
 *
 * @param[in] n element to be searched
 * @param[in] v vector in which to look for
 * @returns true is n belongs to v
 */
bool Bundle::isIn(int n, vector<int> v){

	for(int i=0; i<v.size(); i++){
		if( n == v[i] ){
			return true;
		}
	}
	return false;
}

/**
 * Determine belonging of a vector in a set of vectors
 *
 * @param[in] v vector to be searched
 * @param[in] vlist set of vectors in which to look for
 * @returns true is v belongs to vlist
 */
bool Bundle::isIn(vector<int> v, vector< vector< int > > vlist){
	for(int i=0; i<vlist.size(); i++){
		if( this->isPermutation(v,vlist[i]) ){
			return true;
		}
	}
	return false;
}

/**
 * Check if v1 is a permutation of v2
 *
 * @param[in] v1 first vector
 * @param[in] v2 second vector
 * @returns true is v1 is a permutation of v2
 */
bool Bundle::isPermutation(vector<int> v1, vector<int> v2){
	for( int i=0; i<v1.size(); i++ ){
		if( !this->isIn(v1[i],v2) ){
			return false;
		}
	}
	return true;
}

/**
 * Check if a matrix is a valid template for the current bundle
 *
 * @param[in] T template matrix to be tested
 * @param[in] card cardinality of the bundle
 * @param[in] dirs directions
 * @returns true T is a valid template
 */
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

