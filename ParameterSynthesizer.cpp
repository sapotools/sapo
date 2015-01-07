/*
 * ParameterSynthesizer.cpp
 *
 *  Created on: Nov 5, 2014
 *      Author: dreossi
 */

#include "ParameterSynthesizer.h"


ParameterSynthesizer::ParameterSynthesizer(DiscreteDynamicalSystem *dynamicalSystem, STL *stl_constraint){

	this->dynamicalSystem = dynamicalSystem;
	this->dynamicalSystemControlPts = this->dynamicalSystem->getTemplateControlPts();
	// initialize the control points of the predicates of the constraint
	this->stl_constraint = initConstraintControlPts(stl_constraint);

}

// Initialize the control points of the current constraint
STL* ParameterSynthesizer::initConstraintControlPts(STL *formula){

	if( formula->getType() == 0 ){ // atomic formula

		lst fog = this->dynamicalSystem->getFog();			// get the composition (f o gamma)
		lst sys_vars = this->dynamicalSystem->getVars();
		lst sub;

		for(int i=0; i<(signed)fog.nops(); i++){
			sub.append(sys_vars[i] == fog[i]);
		}

		// Compute the composition current constraint(f(gamma))
		ex cofog = (formula->getPredicate()).subs(sub);

		// and its control points
		BaseConverter *BC = new BaseConverter((this->dynamicalSystem->getReachSet())->getAlpha(), cofog);
		// store the control points in the atom

		formula->setPredicateControlPts(BC->getBernCoeffs());

	}else{

		initConstraintControlPts(formula->getLeftSubFormula());
		initConstraintControlPts(formula->getRightSubFormula());

	}

	return formula;

}

// Call to the synthesis of STL formula
LinearSystemSet* ParameterSynthesizer::synthesizeSTL(vector< double > base_v, vector< double > lenghts, LinearSystemSet *parameterSet) {
	return this->synthesize(base_v, lenghts, parameterSet, this->stl_constraint);
}

// Synthesis algorithm for STL
LinearSystemSet* ParameterSynthesizer::synthesize(vector< double > base_v, vector< double > lenghts, LinearSystemSet *parameterSet, STL *formula) {

	switch( formula->getType() ){

		// Atomic predicate
		case 0:
			return refineParameters(base_v, lenghts, parameterSet, formula->getPredicateControlPts());
		break;

		// Conjunction
		case 1:{
			LinearSystemSet *LS1 = this->synthesize(base_v, lenghts, parameterSet, formula->getLeftSubFormula());
			LinearSystemSet *LS2 = this->synthesize(base_v, lenghts, parameterSet, formula->getRightSubFormula());
			return LS1->intersectWith(LS2);
		}
		break;

		// Disjunction
		case 2:{
			LinearSystemSet *LS1 = this->synthesize(base_v, lenghts, parameterSet, formula->getLeftSubFormula());
			LinearSystemSet *LS2 = this->synthesize(base_v, lenghts, parameterSet, formula->getRightSubFormula());
			return LS1->unionWith(LS2);
		}
		break;

		// Until
		case 3:
			return this->synthesizeUntil(base_v, lenghts, parameterSet, formula);
		break;

	}

	return parameterSet;

}

// Special procedure for Unitl formulas
LinearSystemSet* ParameterSynthesizer::synthesizeUntil(vector< double > base_v, vector< double > lenghts, LinearSystemSet *parameterSet, STL *formula) {

	LinearSystemSet* result = new LinearSystemSet();
	int a = formula->getA();
	int b = formula->getB();

	// Until interval far
	if((a > 0) && (b > 0)){
		// Synthesize wrt phi1
		LinearSystemSet *P1 = this->synthesize(base_v, lenghts, parameterSet, formula->getLeftSubFormula());
		if( P1->isEmpty() ){
			return P1;
		}else{

			// shift until interval
			formula->setA(a-1);
			formula->setB(b-1);

			// Reach step wrt to the i-th linear system of P1
			for(int i=0; i<P1->size(); i++){
				poly_values pv = reachStep(base_v, lenghts, P1->at(i));
				LinearSystemSet* tmpLSset = new LinearSystemSet(P1->at(i));
				tmpLSset = synthesizeUntil(pv.base_vertex, pv.lenghts, tmpLSset, formula);
				result = result->unionWith(tmpLSset);
			}
			return result;
		}
	}

	// Inside until interval
	if((a == 0) && (b > 0)){

		// Refine wrt phi1 and phi2
		LinearSystemSet *P1 = this->synthesize(base_v, lenghts, parameterSet, formula->getLeftSubFormula());
		LinearSystemSet *P2 = this->synthesize(base_v, lenghts, parameterSet, formula->getRightSubFormula());

		if( P1->isEmpty() ){
			return P2;
		}

		// shift until interval
		formula->setB(b-1);

		// Reach step wrt to the i-th linear system of P1
		for(int i=0; i<P1->size(); i++){
			poly_values pv = reachStep(base_v, lenghts, P1->at(i));
			LinearSystemSet* tmpLSset = new LinearSystemSet(P1->at(i));
			tmpLSset = synthesizeUntil(pv.base_vertex, pv.lenghts, tmpLSset, formula);
			result = result->unionWith(tmpLSset);
		}

		return P2->unionWith(result);

	}

	// Base case
	if((a == 0) && (b == 0)){
		return this->synthesize(base_v, lenghts, parameterSet, formula->getRightSubFormula());
	}

	return result;

}

// Refine parameters (atomic case)
LinearSystemSet* ParameterSynthesizer::refineParameters(vector< double > base_v, vector< double > lenghts, LinearSystemSet *parameterSet, lst constraintControlPts){

	lst q = (this->dynamicalSystem->getReachSet())->getQ();
	lst beta = (this->dynamicalSystem->getReachSet())->getBeta();
	lst sub;

	for(int i=0; i<(signed)q.nops(); i++){
			sub.append(q[i] == base_v[i]);
	}
	for(int i=0; i<(signed)beta.nops(); i++){
		sub.append(beta[i] == lenghts[i]);
	}

	// Evaluate numerically the actual constraint control points
	lst num_constraintControlPts;
	for(int i=0; i<(signed)constraintControlPts.nops(); i++){
		num_constraintControlPts.append(constraintControlPts[i].subs(sub));
	}

	// Intersect the control points with the parameter set
	LinearSystem *num_constraintLS = new LinearSystem(this->dynamicalSystem->getParams(), num_constraintControlPts);
	LinearSystemSet *controlPtsLS = new LinearSystemSet(num_constraintLS);

	return parameterSet->intersectWith(controlPtsLS);

}

LinearSystem* ParameterSynthesizer::refineParameters(vector< double > base_v, vector< double > lenghts, LinearSystem *parameterSet){

	lst q = (this->dynamicalSystem->getReachSet())->getQ();
	lst beta = (this->dynamicalSystem->getReachSet())->getBeta();
	lst sub;


	for(int i=0; i<(signed)q.nops(); i++){
		sub.append(q[i] == base_v[i]);
	}
	for(int i=0; i<(signed)beta.nops(); i++){
		sub.append(beta[i] == lenghts[i]);
	}

	// Evaluate numerically the actual constraint control points
	lst num_constraintControlPts;
	for(int i=0; i<(signed)this->constraintControlPts.nops(); i++){
		num_constraintControlPts.append(this->constraintControlPts[i].subs(sub));
	}

	LinearSystem *num_constraintLS = new LinearSystem(this->dynamicalSystem->getParams(), num_constraintControlPts);
	parameterSet = parameterSet->appendLinearSystem(num_constraintLS);

	return parameterSet;
}

// Reachability step
poly_values ParameterSynthesizer::reachStep(vector< double > base_v, vector< double > lenghts, LinearSystem *parameterSet){

	lst q = (this->dynamicalSystem->getReachSet())->getQ();
	lst beta = (this->dynamicalSystem->getReachSet())->getBeta();
	lst params = this->dynamicalSystem->getParams();
	lst sub;


	for(int i=0; i<(signed)q.nops(); i++){
		sub.append(q[i] == base_v[i]);
	}
	for(int i=0; i<(signed)beta.nops(); i++){
		sub.append(beta[i] == lenghts[i]);
	}

	// Compute numerical control points of the dynamics
	vector< lst > num_dynamicalSystemControlPts;
	for(int i=0; i<(signed)this->dynamicalSystemControlPts.size(); i++){
		lst num_dynamicalSystemControlPts_i;
		for(int j=0; j<(signed)this->dynamicalSystemControlPts[i].nops(); j++){
			num_dynamicalSystemControlPts_i.append(this->dynamicalSystemControlPts[i][j].subs(sub));
		}
		num_dynamicalSystemControlPts.push_back(num_dynamicalSystemControlPts_i);
	}

	vector<double> max_control_pts;
	for(int i=0; i<(signed)num_dynamicalSystemControlPts.size(); i++){

			vector< double > max_control_pts_i;

			for(int j=0; j<(signed)num_dynamicalSystemControlPts[i].nops(); j++){
				max_control_pts_i.push_back(parameterSet->maxLinearSystem(params,num_dynamicalSystemControlPts[i][j]));
			}

			max_control_pts.push_back(this->maxVec(max_control_pts_i));
	}

	// Extract the template matrix of the reach set
	vector< vector< double> > temp_mat = this->dynamicalSystem->getReachSet()->getTemplate();
	LinearSystem *reachSet = new LinearSystem(temp_mat,max_control_pts);
	poly_values pv = this->dynamicalSystem->getReachSet()->const2gen(reachSet);

	return pv;
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

