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
public:
    ControlPointStorage() {}
    ControlPointStorage(const ControlPointStorage& orig)
{
    for (auto it = std::cbegin(orig._genF_ctrlP); it!= std::end(orig._genF_ctrlP); ++it) {
        _genF_ctrlP[it->first] = std::pair<GiNaC::lst, GiNaC::lst>(GiNaC::lst(it->second.first), GiNaC::lst(it->second.second));
    }
}

    std::pair<GiNaC::lst, GiNaC::lst> get(const std::vector<int> index) const
{
    return _genF_ctrlP.at(index);
}

    GiNaC::lst get_gen_fun(const std::vector<int> index) const
{
    return _genF_ctrlP.at(index).first;
}

    GiNaC::lst get_ctrl_pts(const std::vector<int> index) const
{
    return _genF_ctrlP.at(index).second;
}

    bool gen_fun_is_equal_to(const std::vector<int> index, const GiNaC::lst& genFun) const
{
    return _genF_ctrlP.at(index).first.is_equal(genFun);
}

    bool contains(const std::vector<int> index) const
{
    return (_genF_ctrlP.count(index)>0);
}

    ControlPointStorage& set(const std::vector<int> index, const GiNaC::lst& genFun, const GiNaC::lst& ctrlPts)
    {
        this->set_gen_fun(index, genFun);
        this->set_ctrl_pts(index, ctrlPts);

        return *this;
    }

    ControlPointStorage& set_gen_fun(const std::vector<int> index, const GiNaC::lst& genFun)
    {
        _genF_ctrlP[index].first=genFun;

        return *this;
    }

    ControlPointStorage& set_ctrl_pts(const std::vector<int> index, const GiNaC::lst& ctrlPts)
    {
        _genF_ctrlP[index].second=ctrlPts;

        return *this;
    }
};

#endif // CONTROLPOINTS_H_