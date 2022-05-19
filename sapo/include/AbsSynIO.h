#ifndef __ABSSYNIO_H__
#define __ABSSYNIO_H__

#include <iostream>
#include <utility> // pair
#include <vector>

//#include "AbsSyn.h"

#include "types.h"

namespace std
{

template<typename T>
std::ostream &operator<<(std::ostream &os, const std::pair<T, T> p)
{
  return os << "[" << p.first << "," << p.second << "]";
}

std::ostream &operator<<(std::ostream &os, const AbsSyn::problemType t);
std::ostream &operator<<(std::ostream &os, const AbsSyn::modeType t);

}
#endif
