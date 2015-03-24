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

		// and its control points and store the control points in the atom
		if(!this->dynamicalSystem->isRational()){
			BaseConverter *BC = new BaseConverter((this->dynamicalSystem->getReachSet())->getAlpha(), cofog);
			formula->setPredicateControlPts(BC->getBernCoeffs());
		}else{
			BaseConverter *BC = new BaseConverter((this->dynamicalSystem->getReachSet())->getAlpha(), cofog.numer(),cofog.denom());
			formula->setPredicateControlPts(BC->getRationalBernCoeffs());
		}
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

		// Always
		case 4:
			return this->synthesizeAlways(base_v, lenghts, parameterSet, formula);
		break;

		// Eventually
		case 5:
			return this->synthesizeEventually(base_v, lenghts, parameterSet, formula);
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

// Special procedure for Always formulas
LinearSystemSet* ParameterSynthesizer::synthesizeAlways(vector< double > base_v, vector< double > lenghts, LinearSystemSet *parameterSet, STL *formula) {

	LinearSystemSet* result = new LinearSystemSet();
	int a = formula->getA();
	int b = formula->getB();

	// Always interval far
	if((a > 0) && (b > 0)){

		// shift until interval
		formula->setA(a-1);
		formula->setB(b-1);

		// Reach step wrt to the i-th linear system of parameterSet
		for(int i=0; i<parameterSet->size(); i++){
			poly_values pv = reachStep(base_v, lenghts, parameterSet->at(i));
			LinearSystemSet* tmpLSset = new LinearSystemSet(parameterSet->at(i));
			tmpLSset = synthesizeAlways(pv.base_vertex, pv.lenghts, tmpLSset, formula);
			result = result->unionWith(tmpLSset);
		}
		return result;
	}

	// Inside Always interval
	if((a == 0) && (b > 0)){

		// Refine wrt phi
		LinearSystemSet *P = this->synthesize(base_v, lenghts, parameterSet, formula->getSubFormula());

		if(!P->isEmpty()){

			// shift until interval
			formula->setB(b-1);

			// Reach step wrt to the i-th linear system of P
			for(int i=0; i<P->size(); i++){
				poly_values pv = reachStep(base_v, lenghts, P->at(i));
				LinearSystemSet* tmpLSset = new LinearSystemSet(P->at(i));
				tmpLSset = synthesizeAlways(pv.base_vertex, pv.lenghts, tmpLSset, formula);
				result = result->unionWith(tmpLSset);
			}

			return result;
		}

		return P;

	}

	// Base case
	if((a == 0) && (b == 0)){
		return this->synthesize(base_v, lenghts, parameterSet, formula->getSubFormula());
	}

	return result;

}

// Special procedure for Eventually formulas
LinearSystemSet* ParameterSynthesizer::synthesizeEventually(vector< double > base_v, vector< double > lenghts, LinearSystemSet *parameterSet, STL *formula) {

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
			poly_values pv = reachStep(base_v, lenghts, parameterSet->at(i));
			LinearSystemSet* tmpLSset = new LinearSystemSet(parameterSet->at(i));
			tmpLSset = synthesizeEventually(pv.base_vertex, pv.lenghts, tmpLSset, formula);
			result = result->unionWith(tmpLSset);
		}
		return result;
	}

	// Inside until interval
	if((a == 0) && (b > 0)){

		// Refine wrt phi
		LinearSystemSet *P = this->synthesize(base_v, lenghts, parameterSet, formula->getSubFormula());

		// shift until interval
		formula->setB(b-1);

		// Reach step wrt to the i-th linear system of P1
		for(int i=0; i<parameterSet->size(); i++){
			poly_values pv = reachStep(base_v, lenghts, parameterSet->at(i));
			LinearSystemSet* tmpLSset = new LinearSystemSet(parameterSet->at(i));
			tmpLSset = synthesizeEventually(pv.base_vertex, pv.lenghts, tmpLSset, formula);
			result = result->unionWith(tmpLSset);
		}

		return P->unionWith(result);

	}

	// Base case
	if((a == 0) && (b == 0)){
		return this->synthesize(base_v, lenghts, parameterSet, formula->getSubFormula());
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

	// Eliminate redundant constraints
	num_constraintControlPts.unique();

	// Intersect the control points with the parameter set
	if(!this->dynamicalSystem->isRational()){	// check if there are rational control points
		LinearSystem *num_constraintLS = new LinearSystem(this->dynamicalSystem->getParams(), num_constraintControlPts);
		LinearSystemSet *controlPtsLS = new LinearSystemSet(num_constraintLS);
		return parameterSet->intersectWith(controlPtsLS);
	}else{
		lst num_constraintControlPts_numer; // numerators
		lst num_constraintControlPts_numer_minus; // numerators with swapped sign
		lst num_constraintControlPts_denom; // denominators
		lst num_constraintControlPts_denom_minus; // denominators with swapped sign
		for(int i=0; i<num_constraintControlPts.nops();i++){
			num_constraintControlPts_numer.append(num_constraintControlPts[i].numer());
			num_constraintControlPts_numer_minus.append(-num_constraintControlPts_numer[i]);
			num_constraintControlPts_denom.append(num_constraintControlPts[i].denom());
			num_constraintControlPts_denom_minus.append(-num_constraintControlPts_denom[i]);
		}
		// Concatenate the numerator and denominator control points
		for(int i=0; i<num_constraintControlPts_denom_minus.nops();i++){
			num_constraintControlPts_numer.append(num_constraintControlPts_denom_minus[i]);
		}
		for(int i=0; i<num_constraintControlPts_denom.nops();i++){
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
		return parameterSet->intersectWith(controlPtsLS_1)->unionWith(parameterSet->intersectWith(controlPtsLS_1));

	}
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


	for(int i=0; i<pv.lenghts.size(); i++){
		cout<<"["<<pv.base_vertex[i]<<","<<pv.base_vertex[i]+pv.versors[i][i]*pv.lenghts[i]<<"]";
	}
	cout<<"\n";

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

