/*
 * CoeffMap.h
 *
 *  Created on: Sep 23, 2015
 *      Author: dreossi
 */

#ifndef COEFFMAP_H_
#define COEFFMAP_H_

#include "Common.h"

class CoeffMap {

private:
	vector< vector<int> > keys;
	vector< lst > elements;

public:
	CoeffMap();
	virtual ~CoeffMap();

	bool isIn(vector<int> key);
	void insert(vector<int> key, lst element);
	int find(vector<int> key);
	lst get(int i);

};

#endif /* COEFFMAP_H_ */
