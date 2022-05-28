/**
 * @file PolytopesUnion.h
 * @author Alberto Casagrande (acasagrande@units.it)
 * @brief Simplify polytopes union representations
 * @version 0.1
 * @date 2022-05-28
 * 
 * @copyright Copyright (c) 2022
 */

#ifndef POLYTOPESUNION_H_
#define POLYTOPESUNION_H_

#include "SetsUnion.h"
#include "Polytope.h"

/**
 * @brief Simplify the union representation
 *
 * This method simplifies the union representation 
 * by removing redundant constraints from the polytope
 * representations.
 * 
 * @param[in,out] polytope_union is the union whose representation 
 *                must be simplified
 * @return a reference to the simplified union
 */
SetsUnion<Polytope> &simplify(SetsUnion<Polytope>& polytope_union);

#endif /* POLYTOPESUNION_H_ */
