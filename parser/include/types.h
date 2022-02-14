#ifndef __TYPES_H__
#define __TYPES_H__

#include <iostream>

namespace AbsSyn
{

// defines the types of problems solved by SAPO
enum problemType {
  P_UNDEF, // undefined
  REACH,   // reachability
  SYNTH    // parameter synthesis
};


// defines modalities for representing variables and parameters
enum modeType {
  M_UNDEF, // undefined
  BOX,     // boxes
  PARAL,   // parallelotopes
  POLY     // polytopes
};

// defines types of transformation
enum transType {
  T_UNDEF, // not yet defined
  AFO,     // All-For-One
  OFO      // One-For-One
};

}

namespace std
{

inline ostream &operator<<(ostream &os, const AbsSyn::problemType t)
{
  if (t == AbsSyn::problemType::P_UNDEF)
    return os << "UNDEF";
  else if (t == AbsSyn::problemType::REACH)
    return os << "reachability";
  else
    return os << "synthesis";
}

inline ostream &operator<<(ostream &os, const AbsSyn::modeType t)
{
  if (t == AbsSyn::modeType::M_UNDEF)
    return os << "UNDEF";
  else if (t == AbsSyn::modeType::BOX)
    return os << "boxes";
  else if (t == AbsSyn::modeType::PARAL)
    return os << "parallelotopes";
  else
    return os << "polytopes";
}

}

#endif
