#ifndef __ABSSYNIO_H__
#define __ABSSYNIO_H__

#include <iostream>
#include <utility> // pair
#include <vector>

#include "AbsSyn.h"

namespace std
{
template<typename T>
std::ostream &operator<<(std::ostream &os, const std::vector<std::vector<T>> &v)
{
  for (unsigned i = 0; i < v.size(); i++) {
    os << "{" << v[i] << "}";
    if (i+1 == v.size())
      os << endl;
  }
  return os;
}

template<typename T>
std::ostream &operator<<(std::ostream &os, const std::vector<T> v)
{
  for (unsigned i = 0; i < v.size(); i++)
    os << v[i] << (i+1 == v.size() ? "" : ",");

  return os;
}

template<typename T>
std::ostream &operator<<(std::ostream &os, const std::pair<T, T> p)
{
  return os << "[" << p.first << "," << p.second << "]";
}

std::ostream &operator<<(std::ostream &os, const AbsSyn::problemType t);
std::ostream &operator<<(std::ostream &os, const AbsSyn::modeType t);

}
#endif
