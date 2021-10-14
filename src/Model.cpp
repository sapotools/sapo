/**
 * @file Model.cpp
 * Model implementation
 *
 * @author Alberto Casagrande <acasagrande@units.it>
 * @version 0.1
 */

#include "Model.h"

Model::~Model()
{
    delete this->reachSet;
    delete this->paraSet;
	delete this->spec;
}