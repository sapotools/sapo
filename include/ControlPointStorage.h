/**
 * @file ControlPointStorage.h
 * @brief Defines the ControlPointStorage class.
 *
 * This file contains the definition of the ControlPointStorage class.
 * This class is meant to represent a control point storage for
 * Bundle transformations.
 *
 * @author Alberto Casagrande <alberto.casagrande@units.it>
 * @version 0.1
 */

#ifndef CONTROLPOINTS_H_
#define CONTROLPOINTS_H_

#include <shared_mutex>
#include <ginac/ginac.h>

class ControlPointStorage
{
  std::map<std::vector<int>, std::pair<GiNaC::lst, GiNaC::lst>> _genF_ctrlP;

  mutable std::shared_timed_mutex
      _mutex; //! A mutex to safely access to the data structure
public:
  ControlPointStorage() {}
  ControlPointStorage(const ControlPointStorage &orig);

  std::pair<GiNaC::lst, GiNaC::lst> get(const std::vector<int> index) const;

  GiNaC::lst get_gen_fun(const std::vector<int> index) const;

  GiNaC::lst get_ctrl_pts(const std::vector<int> index) const;

  bool gen_fun_is_equal_to(const std::vector<int> index,
                           const GiNaC::lst &genFun) const;

  bool contains(const std::vector<int> index) const;

  ControlPointStorage &set(const std::vector<int> index,
                           const GiNaC::lst &genFun, const GiNaC::lst &ctrlPts)
  {
    this->set_gen_fun(index, genFun);
    this->set_ctrl_pts(index, ctrlPts);

    return *this;
  }

  ControlPointStorage &set_gen_fun(const std::vector<int> index,
                                   const GiNaC::lst &genFun);
  ControlPointStorage &set_ctrl_pts(const std::vector<int> index,
                                    const GiNaC::lst &ctrlPts);
};

#endif // CONTROLPOINTS_H_