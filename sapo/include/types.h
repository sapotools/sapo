#ifndef __TYPES_H__
#define __TYPES_H__

#include <iostream>

namespace AbsSyn
{

// defines the types of problems solved by SAPO
enum problemType {
  P_UNDEF, // undefined
  REACH,   // reachability
  SYNTH,   // parameter synthesis
  INVARIANT // invariant validation
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
  switch(t) {
    case AbsSyn::problemType::P_UNDEF:
      os << "Undefined";
      break;
    case AbsSyn::problemType::REACH:
      os << "reachability";
      break;
    case AbsSyn::problemType::SYNTH:
      os << "synthesis";
      break;
    case AbsSyn::problemType::INVARIANT:
      os << "invariant validation";
      break;
    default:
      os << "Unknown problem type ("<< (unsigned int)t << ")";
      break;
  }

  return os;
}

inline ostream &operator<<(ostream &os, const AbsSyn::modeType t)
{
  switch(t) {
    case AbsSyn::modeType::M_UNDEF:
      os << "Undefined";
      break;
    case AbsSyn::modeType::BOX:
      os << "boxes";
      break;
    case AbsSyn::modeType::PARAL:
      os << "parallelotopes";
      break;
    case AbsSyn::modeType::POLY:
      os << "polytopes";
      break;
    default:
      os << "Unknown mode type ("<< (unsigned int)t << ")";
      break;
  }

  return os;
}

}

#endif
