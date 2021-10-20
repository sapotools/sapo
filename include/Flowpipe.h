/**
 * @file Flowpipe.h
 * Represent and manipulate flowpipes of bundle
 * Used to represent the reachable set of a dynamical system
 *
 * @author Tommaso Dreossi <tommasodreossi@berkeley.edu>
 * @version 0.1
 */

#ifndef FLOWPIPE_H_
#define FLOWPIPE_H_

#include "Common.h"
#include "Bundle.h"

class Flowpipe {

private:
	vector< Bundle* > flowpipe;			// flowpipe

public:

	// constructors
	Flowpipe();
	Flowpipe(vector< Bundle* >);

	const Bundle* get(const unsigned int i) const;	// get i-th bundle

	void append( Bundle* bundle );

	inline const std::size_t size() const { return this->flowpipe.size(); }

	void print() const;
	void plotRegion() const;
	void plotRegionToFile(const char *file_name, const char color) const;
	void plotProjToFile(const unsigned int var, const double time_step,
					            const char *file_name, const char color) const;

	virtual ~Flowpipe();
};

#endif /* BUNDLE_H_ */
