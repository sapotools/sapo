/*
 * ParameterSynthesizer.cpp
 *
 *  Created on: Apr 5, 2014
 *      Author: Tommaso Dreossi
 */

#include "ParameterSynthesizer.h"


ParameterSynthesizer::ParameterSynthesizer(DiscreteDynamicalSystem *dynamicalSystem, STL *stl_constraint, synthesizer_opt options){

	this->dynamicalSystem = dynamicalSystem;
	this->dynamicalSystemControlPts = this->dynamicalSystem->getTemplatesControlPts();

	// initialize the control points of the predicates of the constraint
	this->stl_constraint = initConstraintControlPts(stl_constraint);
	this->options = options;

}

// Initialize the control points of the current constraint
STL* ParameterSynthesizer::initConstraintControlPts(STL *formula){

	if( formula->getType() == 0 ){ // atomic formula

		vector< lst > fogs = this->dynamicalSystem->getFogs();			// get the composition (f o gamma)
		lst sys_vars = this->dynamicalSystem->getVars();
		lst sub;

		vector< lst > bernCoeffs;

		for(int j=0; j<(signed)fogs.size(); j++){

			for(int i=0; i<(signed)fogs[j].nops(); i++){
				sub.append(sys_vars[i] == fogs[j][i]);
			}

			// Compute the composition current constraint(f(gamma))
			ex cofog = (formula->getPredicate()).subs(sub);

			// and its control points and store the control points in the atom
			if(!this->dynamicalSystem->isRational()){
				BaseConverter *BC = new BaseConverter((this->dynamicalSystem->getReachSet(j))->getAlpha(), cofog);
				bernCoeffs.push_back(BC->getBernCoeffs());
			}else{
				BaseConverter *BC = new BaseConverter((this->dynamicalSystem->getReachSet(j))->getAlpha(), cofog.numer(),cofog.denom());
				bernCoeffs.push_back(BC->getBernCoeffs());
			}
		}
		formula->setPredicateControlPts(bernCoeffs);
	}else{
		if( (formula->getType() == 4) || (formula->getType() == 5)){ // always or eventually
			initConstraintControlPts(formula->getSubFormula());
		}else{
			initConstraintControlPts(formula->getLeftSubFormula());
			initConstraintControlPts(formula->getRightSubFormula());
		}
	}

	return formula;

}

// Call to the synthesis of STL formula
LinearSystemSet* ParameterSynthesizer::synthesizeSTL(poly_values reach_set, LinearSystemSet *parameterSet) {
	vector< poly_values > reach_sets (1,reach_set);
	return this->synthesize(reach_sets, parameterSet, this->stl_constraint);
}

LinearSystemSet* ParameterSynthesizer::synthesizeSTL(vector<poly_values> reach_sets, LinearSystemSet *parameterSet) {
	return this->synthesize(reach_sets, parameterSet, this->stl_constraint);
}

// Synthesis algorithm for STL
LinearSystemSet* ParameterSynthesizer::synthesize(vector<poly_values> reach_sets, LinearSystemSet *parameterSet, STL *formula) {

	switch( formula->getType() ){

		// Atomic predicate
		case 0:
			return refineParameters(reach_sets, parameterSet, formula->getPredicateControlPts());
		break;

		// Conjunction
		case 1:{
			LinearSystemSet *LS1 = this->synthesize(reach_sets, parameterSet, formula->getLeftSubFormula());
			LinearSystemSet *LS2 = this->synthesize(reach_sets, parameterSet, formula->getRightSubFormula());
			return LS1->intersectWith(LS2);
		}
		break;

		// Disjunction
		case 2:{
			LinearSystemSet *LS1 = this->synthesize(reach_sets, parameterSet, formula->getLeftSubFormula());
			LinearSystemSet *LS2 = this->synthesize(reach_sets, parameterSet, formula->getRightSubFormula());
			// Check result type (all set or largest set)
			if( !this->options.largest_para_set ){
				return LS1->unionWith(LS2);
			}else{	// return the set with max volume
				if(LS1->boundingVol() > LS2->boundingVol()){
					return LS1;
				}else{
					return LS2;
				}
			}
		}
		break;

		// Until
		case 3:
			//return this->synthesizeUntil(base_v, lenghts, parameterSet, formula);
		break;

		// Always
		case 4:
			//return this->synthesizeAlways(base_v, lenghts, parameterSet, formula);
		break;

		// Eventually
		case 5:
			//return this->synthesizeEventually(base_v, lenghts, parameterSet, formula);
		break;

	}

	return parameterSet;

}

// Special procedure for Unitl formulas
LinearSystemSet* ParameterSynthesizer::synthesizeUntil(vector<poly_values> reach_sets, LinearSystemSet *parameterSet, STL *formula) {

	LinearSystemSet* result = new LinearSystemSet();
	int a = formula->getA();
	int b = formula->getB();

	// Until interval far
	if((a > 0) && (b > 0)){
		// Synthesize wrt phi1
		LinearSystemSet *P1 = this->synthesize(reach_sets, parameterSet, formula->getLeftSubFormula());
		if( P1->isEmpty() ){
			return P1;			// false until
		}else{
			// shift until interval
			formula->setA(a-1);
			formula->setB(b-1);

			if(this->options.largest_para_set){		// by assumption there's only one parameter set
				vector<poly_values> new_reach_sets = reachStep(reach_sets, P1->at(0));
				result = synthesizeUntil(new_reach_sets, P1, formula);
			}else{
				// Reach step wrt to the i-th linear system of P1
				for(int i=0; i<P1->size(); i++){
					vector<poly_values> new_reach_sets = reachStep(reach_sets, P1->at(i));
					LinearSystemSet* tmpLSset = new LinearSystemSet(P1->at(i));
					tmpLSset = synthesizeUntil(new_reach_sets, tmpLSset, formula);
					result = result->unionWith(tmpLSset);
				}
			}
			return result;
		}
	}

	// Inside until interval
	if((a == 0) && (b > 0)){

		// Refine wrt phi1 and phi2
		LinearSystemSet *P1 = this->synthesize(reach_sets, parameterSet, formula->getLeftSubFormula());
		LinearSystemSet *P2 = this->synthesize(reach_sets, parameterSet, formula->getRightSubFormula());

		if( P1->isEmpty() ){
			return P2;
		}

		// shift until interval
		formula->setB(b-1);

		if( this->options.largest_para_set ){
			if( P2->boundingVol() > P1->boundingVol() ){
				return P2;
			}else{
				vector<poly_values> new_reach_sets = reachStep(reach_sets, P1->at(0)); 		// by assumption there's only one parameter set
				result = synthesizeUntil(new_reach_sets, P1, formula);
				if( result->boundingVol() > P2->boundingVol() ){
					return result;
				}else{
					return P2;
				}
			}
		}else{
			// Reach step wrt to the i-th linear system of P1
			for(int i=0; i<P1->size(); i++){
				vector<poly_values> new_reach_sets = reachStep(reach_sets, P1->at(i));
				LinearSystemSet* tmpLSset = new LinearSystemSet(P1->at(i));
				tmpLSset = synthesizeUntil(new_reach_sets, tmpLSset, formula);
				result = result->unionWith(tmpLSset);
			}
			return P2->unionWith(result);
		}

	}

	// Base case
	if((a == 0) && (b == 0)){
		return this->synthesize(reach_sets, parameterSet, formula->getRightSubFormula());
	}

	return result;

}


// Special procedure for Always formulas
LinearSystemSet* ParameterSynthesizer::synthesizeAlways(vector<poly_values> reach_sets, LinearSystemSet *parameterSet, STL *formula) {

	LinearSystemSet* result = new LinearSystemSet();
	int a = formula->getA();
	int b = formula->getB();

	// Always interval far
	if((a > 0) && (b > 0)){

		// shift always interval
		formula->setA(a-1);
		formula->setB(b-1);

		// Reach step wrt to the i-th linear system of parameterSet
		for(int i=0; i<parameterSet->size(); i++){
			vector<poly_values> new_reach_sets = reachStep(reach_sets, parameterSet->at(i));
			LinearSystemSet* tmpLSset = new LinearSystemSet(parameterSet->at(i));
			tmpLSset = synthesizeAlways(new_reach_sets, tmpLSset, formula);
			if(this->options.largest_para_set){
				if( result->boundingVol() < tmpLSset->boundingVol() ){	// store the largest result
					result = tmpLSset;
				}
			}else{
				result = result->unionWith(tmpLSset);
			}
		}
		return result;
	}

	// Inside Always interval
	if((a == 0) && (b > 0)){

		// Refine wrt phi
		LinearSystemSet *P = this->synthesize(reach_sets, parameterSet, formula->getSubFormula());

		if(!P->isEmpty()){

			// shift until interval
			formula->setB(b-1);

			// Reach step wrt to the i-th linear system of P
			for(int i=0; i<P->size(); i++){
				vector<poly_values> new_reach_sets = reachStep(reach_sets, P->at(i));
				LinearSystemSet* tmpLSset = new LinearSystemSet(P->at(i));
				tmpLSset = synthesizeAlways(new_reach_sets, tmpLSset, formula);
				if(this->options.largest_para_set){
					if( result->boundingVol() < tmpLSset->boundingVol() ){	// store the largest result
						result = tmpLSset;
					}
					}else{
						result = result->unionWith(tmpLSset);
				}
			}

			return result;
		}

		return P;

	}

	// Base case
	if((a == 0) && (b == 0)){
		return this->synthesize(reach_sets, parameterSet, formula->getSubFormula());
	}

	return result;

}


// Special procedure for Eventually formulas
LinearSystemSet* ParameterSynthesizer::synthesizeEventually(vector<poly_values> reach_sets, LinearSystemSet *parameterSet, STL *formula) {

	LinearSystemSet* result = new LinearSystemSet();
	int a = formula->getA();
	int b = formula->getB();

	// Eventually interval far
	if((a > 0) && (b > 0)){

		// shift until interval
		formula->setA(a-1);
		formula->setB(b-1);

		// Reach step wrt to the i-th linear system of parameterSet
		for(int i=0; i<parameterSet->size(); i++){
			vector<poly_values> new_reach_sets = reachStep(reach_sets, parameterSet->at(i));
			LinearSystemSet* tmpLSset = new LinearSystemSet(parameterSet->at(i));
			tmpLSset = synthesizeEventually(new_reach_sets, tmpLSset, formula);
			if(this->options.largest_para_set){
				if( result->boundingVol() < tmpLSset->boundingVol() ){	// store the largest result
					result = tmpLSset;
				}
			}else{
				result = result->unionWith(tmpLSset);
			}
		}
		return result;
	}

	// Inside eventually interval
	if((a == 0) && (b > 0)){

		// Refine wrt phi
		LinearSystemSet *P = this->synthesize(reach_sets, parameterSet, formula->getSubFormula());

		// shift until interval
		formula->setB(b-1);
		result = P;

		// Reach step wrt to the i-th linear system of P1
		for(int i=0; i<parameterSet->size(); i++){
			vector<poly_values> new_reach_sets = reachStep(reach_sets, parameterSet->at(i));
			LinearSystemSet* tmpLSset = new LinearSystemSet(parameterSet->at(i));
			tmpLSset = synthesizeEventually(new_reach_sets, tmpLSset, formula);
			if(this->options.largest_para_set){
				if( result->boundingVol() < tmpLSset->boundingVol() ){	// store the largest result
					result = tmpLSset;
				}
			}else{
				result = result->unionWith(tmpLSset);
			}
		}


		return result;

	}

	// Base case
	if((a == 0) && (b == 0)){
		return this->synthesize(reach_sets, parameterSet, formula->getSubFormula());
	}

	return result;

}

// Refine parameters (atomic case)
LinearSystemSet* ParameterSynthesizer::refineParameters(vector< poly_values > reach_sets, LinearSystemSet *parameterSet, vector< lst > constraintControlPts){

	LinearSystemSet *result = new LinearSystemSet();

	// for each template
	for( int t=0; t<(signed)constraintControlPts.size(); t++ ){

		lst q = (this->dynamicalSystem->getReachSet(t))->getQ();
		lst beta = (this->dynamicalSystem->getReachSet(t))->getBeta();
		lst sub;

		for(int i=0; i<(signed)q.nops(); i++){
				sub.append(q[i] == reach_sets[t].base_vertex[i]);
				sub.append(beta[i] == reach_sets[t].lenghts[i]);
		}

		// Evaluate numerically the actual constraint control points
		lst num_constraintControlPts;
		for(int i=0; i<(signed)constraintControlPts[t].nops(); i++){
			num_constraintControlPts.append(constraintControlPts[t][i].subs(sub));
		}

		// Eliminate redundant constraints
		num_constraintControlPts.unique();

		// Intersect the control points with the parameter set
		if(!this->dynamicalSystem->isRational()){	// check if there are rational control points
			LinearSystem *num_constraintLS = new LinearSystem(this->dynamicalSystem->getParams(), num_constraintControlPts);
			LinearSystemSet *controlPtsLS = new LinearSystemSet(num_constraintLS);
			if(this->options.largest_para_set){
				if( parameterSet->intersectWith(controlPtsLS)->boundingVol() > result->boundingVol() ){
					result = parameterSet->intersectWith(controlPtsLS);
				}
			}else{
				result = result->unionWith(parameterSet->intersectWith(controlPtsLS));
			}
		}else{
			lst num_constraintControlPts_numer; // numerators
			lst num_constraintControlPts_numer_minus; // numerators with swapped sign
			lst num_constraintControlPts_denom; // denominators
			lst num_constraintControlPts_denom_minus; // denominators with swapped sign
			for(int i=0; i<(signed)num_constraintControlPts.nops();i++){
				num_constraintControlPts_numer.append(num_constraintControlPts[i].numer());
				num_constraintControlPts_numer_minus.append(-num_constraintControlPts_numer[i]);
				num_constraintControlPts_denom.append(num_constraintControlPts[i].denom());
				num_constraintControlPts_denom_minus.append(-num_constraintControlPts_denom[i]);
			}
			// Concatenate the numerator and denominator control points
			for(int i=0; i<(signed)num_constraintControlPts_denom_minus.nops();i++){
				num_constraintControlPts_numer.append(num_constraintControlPts_denom_minus[i]);
			}
			for(int i=0; i<(signed)num_constraintControlPts_denom.nops();i++){
				num_constraintControlPts_numer_minus.append(num_constraintControlPts_denom[i]);
			}

			// Eliminate reduntant constraints
			num_constraintControlPts_numer.unique();
			num_constraintControlPts_numer_minus.unique();

			LinearSystem *num_constraintLS_1 = new LinearSystem(this->dynamicalSystem->getParams(), num_constraintControlPts_numer);
			LinearSystem *num_constraintLS_2 = new LinearSystem(this->dynamicalSystem->getParams(), num_constraintControlPts_numer_minus);

			LinearSystemSet *controlPtsLS_1 = new LinearSystemSet(num_constraintLS_1);
			LinearSystemSet *controlPtsLS_2 = new LinearSystemSet(num_constraintLS_2);

			// return (P \cap Ref1) \cup (P \cap Ref2)
			if(this->options.largest_para_set){
				if( result->boundingVol() < parameterSet->intersectWith(controlPtsLS_1)->boundingVol() ){
					result = parameterSet->intersectWith(controlPtsLS_1);
				}else{
					if( result->boundingVol() < parameterSet->intersectWith(controlPtsLS_2)->boundingVol() ){
						result = parameterSet->intersectWith(controlPtsLS_2);
					}
				}
			}else{
				result = result->unionWith(parameterSet->intersectWith(controlPtsLS_1)->unionWith(parameterSet->intersectWith(controlPtsLS_2)));
			}
		}
	}
	return result;
}

// Perform a single reachability step
// reach_sets is a vector of base vertices and lengths, one for each template
vector< poly_values > ParameterSynthesizer::reachStep(vector< poly_values > reach_sets, LinearSystem *parameterSet){

	LinearSystem* capReachSets = new LinearSystem();		// Linear system containing the intersection of all the sets

	for(int t=0; t<(signed)reach_sets.size(); t++){

		lst q = (this->dynamicalSystem->getReachSet(t))->getQ();
		lst beta = (this->dynamicalSystem->getReachSet(t))->getBeta();
		lst params = this->dynamicalSystem->getParams();
		lst sub;

		// Substitutions for symbloic control points
		for(int i=0; i<(signed)q.nops(); i++){
			sub.append(q[i] == reach_sets[t].base_vertex[i]);
			sub.append(beta[i] == reach_sets[t].lenghts[i]);
		}

		// Compute numerical control points of the dynamics
		vector< lst > num_dynamicalSystemControlPts;
		for(int i=0; i<(signed)this->dynamicalSystemControlPts[t].size(); i++){
			lst num_dynamicalSystemControlPts_i;
			for(int j=0; j<(signed)this->dynamicalSystemControlPts[t][i].nops(); j++){
				num_dynamicalSystemControlPts_i.append(this->dynamicalSystemControlPts[t][i][j].subs(sub));
			}
			num_dynamicalSystemControlPts.push_back(num_dynamicalSystemControlPts_i);
		}

		// And extract the maximums
		vector<double> max_control_pts;
		for(int i=0; i<(signed)num_dynamicalSystemControlPts.size(); i++){
				vector< double > max_control_pts_i;
				for(int j=0; j<(signed)num_dynamicalSystemControlPts[i].nops(); j++){
					max_control_pts_i.push_back(parameterSet->maxLinearSystem(params,num_dynamicalSystemControlPts[i][j]));
				}
				max_control_pts.push_back(this->maxVec(max_control_pts_i));
		}

		// Extract the template matrix of the reach set
		vector< vector< double> > temp_mat = this->dynamicalSystem->getReachSet(t)->getTemplate();
		LinearSystem *reachSet = new LinearSystem(temp_mat,max_control_pts);

		// Intersect actual template with others
		capReachSets = capReachSets->appendLinearSystem(reachSet);
	}

	capReachSets->plotRegion();

	// Shrink the various templates around the computed intersection
	vector< poly_values > new_reach_sets;

	for(int t=0; t<(signed)reach_sets.size(); t++){	// for each template

		vector< vector< double> > temp_mat = this->dynamicalSystem->getReachSet(t)->getTemplate();
		vector< double > shrinked_offset;

		for(int i=0; i<(signed)temp_mat.size(); i++){									// for each facet
			shrinked_offset.push_back(capReachSets->maxLinearSystem(temp_mat[i]));		// compute the new offset
		}

		LinearSystem *shrinkedReachSet = new LinearSystem(temp_mat,shrinked_offset);

		// Convert each shrank template to the generator representation
		poly_values pv = this->dynamicalSystem->getReachSet(t)->const2gen(shrinkedReachSet);
		new_reach_sets.push_back(pv);
	}

	return new_reach_sets;
}

// Reachability up to k steps
vector<poly_values> ParameterSynthesizer::reach(vector<poly_values> reach_sets, LinearSystem *parameterSet, int k){

	for(int i=0; i<k; i++){
		reach_sets = this->reachStep(reach_sets,parameterSet);
	}

	return reach_sets;
}

// Find the maximum of vec
double ParameterSynthesizer::maxVec( vector<double> vec){

	double max = vec[0];

	for(int i=1; i<(signed)vec.size(); i++){
		if(vec[i] > max){
			max = vec[i];
		}
	}

	return max;

}

ParameterSynthesizer::~ParameterSynthesizer() {
	// TODO Auto-generated destructor stub
}

