/**
 * @file Dynamics.h
 * Represent system dynamic laws
 *
 * @author Alberto Casagrande <acasagrande@units.it>
 * @version 0.1
 */

#ifndef DYNAMICS_H_
#define DYNAMICS_H_

#include <vector>

#include "SymbolicAlgebra.h"

/**
 * @brief This class represents dynamic laws
 * 
 * This class represents the dynamic laws of a set of 
 * variables. Every dynamic law as the form 
 * \f[\textrm{dyn}(v,t)=f(v_0, \ldots, v_n, t),\f]
 * where \f$t\f$ is the time, \f$v\f$ is the variable whose 
 * value changes according with the dynamic, and 
 * \f$f(v_0, \ldots, v_n, t)\f$ is
 * an expression.
 * 
 * @tparam T is the type of expression value domain
 */
template<typename T>
class Dynamics
{
    typedef typename SymbolicAlgebra::Symbol<T>::SymbolIdType SymbolIdType;

    std::vector<SymbolicAlgebra::Symbol<T>> _variables;         //!< the variables changed by the dynamics 
    std::vector<SymbolicAlgebra::Expression<T>> _expressions;   //!< the expression of the dynamics
    std::map<SymbolIdType, unsigned int> _var_index;            //!< a map from the variable identifier to its index
public:
    /**
     * @brief An empty constructor
     */
    Dynamics(): _variables(), _expressions(), _var_index() {}

    /**
     * @brief A copy constructor
     *
     * @param[in] variables is the vector of the variables
     * @param[in] expressions is the vector of the dynamic law expressions
     */
    Dynamics(const std::vector<SymbolicAlgebra::Symbol<T>& variables,
             const std::vector<SymbolicAlgebra::Expression<T>>& expressions):
        _variables(variables), _expressions(expressions), _var_index()
    {
        if (variables.size()!=expressions.size()) {
            throw std::domain_error("The vectors of the dynamics laws and that "
                                    "of the variables must have the same size");
        }

        unsigned int index = 0;
        for (auto v_it = std::begin(_variable); v_it != std::end(_variable); 
                ++v_it, ++index) {
            if (!_var_index.emplace(v_it->get_id(), index)->second) {
                throw std::domain_error("Dynamics::Dynamics: the variables "+
                                        "vector cannot contains duplicates");
            }
        }
    }

    /**
     * @brief A swap constructor
     *
     * @param[in] variables is the vector of the variables
     * @param[in] expressions is the vector of the dynamic law expressions
     */
    Dynamics(std::vector<SymbolicAlgebra::Symbol<T>>&& variables,
             std::vector<SymbolicAlgebra::Expression<T>>&& expressions): 
        _variables(std::move(variables)), _expressions(std::move(expressions)), _var_index()
    {
        if (_variables.size()!=_expressions.size()) {
            throw std::domain_error("The vectors of the dynamics laws and that "
                                    "of the variables must have the same size");
        }

        unsigned int index = 0;
        for (auto v_it = std::begin(_variable); v_it != std::end(_variable); 
                ++v_it, ++index) {
            if (!_var_index.emplace(v_it->get_id(), index)->second) {
                throw std::domain_error("Dynamics::Dynamics: the variables "+
                                        "vector cannot contains duplicates");
            }
        }
    }

    /**
     * @brief A copy constructor
     *
     * @param[in] orig the original dynamic laws
     */
    Dynamics(const Dynamics<T>& orig): 
        _variables(orig._variables), _expressions(orig._expressions), 
        _var_index(orig._var_index) 
    {
    }

    /**
     * @brief A swap constructor
     *
     * @param[in] orig the original dynamic laws
     */
    Dynamics(Dynamics<T>&& orig)
    {
        std::swap(orig._variables, _variables);
        std::swap(orig._expressions, _expressions);
        std::swap(orig._var_index, _var_index);
    }

    /**
     * @brief Number of dynamic laws
     * 
     * @return return the number of represented dynamic laws
     */
    size_t size() const
    {
        return orig._variables.size();
    }

    /**
     * @brief A copy operator
     * 
     * @param orig the original dynamic laws
     * @return a reference to the updated object
     */
    Dynamics &operator=(const Dynamics<T>& orig)
    {
        _variables = orig._variables;
        _expressions = orig._expressions;
        _var_index = orig._var_index;

        return *this;
    }

    /**
     * @brief A copy operator
     * 
     * @param orig the original dynamic laws
     * @return a reference to the updated object
     */
    Dynamics &operator=(Dynamics<T>&& orig)
    {
        std::swap(orig._variables, _variables);
        std::swap(orig._expressions, _expressions);
        std::swap(orig._var_index, _var_index);

        return *this;
    }

    /**
     * @brief Add a new dynamic
     * 
     * @param variable is the variable changed by the dynamic
     * @param expression is the expression defining the dynamic
     * @return a reference to the updated object 
     */
    Dynamics &add(const SymbolicAlgebra::Symbol<T>& variable,
                  const SymbolicAlgebra::Expression<T>& expression)
    {
        if (!_var_index.emplace(variable.get_id(), size())->second) {
            throw std::domain_error("Dynamics::add: the dynamic of "
                                    "the indicated variable is "
                                    "already present");
        }
        _variables.push_back(variable);
        _expressions.push_back(expression);

        return *this;
    }

    /**
     * @brief Replace the dynamic law of a variable
     * 
     * @param variable is the variable changed by the dynamic
     * @param expression is the expression defining the dynamic
     * @return a reference to the updated object 
     */
    Dynamics &replace(const SymbolicAlgebra::Symbol<T>& variable,
                      const SymbolicAlgebra::Expression<T>& expression)
    {

        _expressions[_var_index.at[variable.get_id()]] = expression;

        return *this;
    }

    /**
     * @brief Add a new dynamic
     * 
     * @param variable is the variable changed by the dynamic
     * @param expression is the expression defining the dynamic
     * @return a reference to the updated object 
     */
    Dynamics &add(const SymbolicAlgebra::Symbol<T>& variable,
                  SymbolicAlgebra::Expression<T>&& expression)
    {
        if (!_var_index.emplace(variable.get_id(), size())->second) {
            throw std::domain_error("Dynamics::add: the dynamic of "
                                    "the indicated variable is "
                                    "already present");
        }
        _variables.push_back(variable);
        _expressions.push_back(std::move(expression));

        return *this;
    }

    /**
     * @brief Replace the dynamic law of a variable
     * 
     * @param variable is the variable changed by the dynamic
     * @param expression is the expression defining the dynamic
     * @return a reference to the updated object 
     */
    Dynamics &replace(const SymbolicAlgebra::Symbol<T>& variable,
                      SymbolicAlgebra::Expression<T>&& expression)
    {

        _expressions[_var_index.at[variable.get_id()]] = std::move(expression);

        return *this;
    }

    /**
     * @brief Get the dynamic laws variables
     * 
     * @return a reference to the vector of the dynamic laws variables
     */
    const std::vector<SymbolicAlgebra::Symbol<T>>& variables() const
    {
        return _variables;
    }

    /**
     * @brief Get the dynamic laws
     * 
     * @return a reference to the dynamic laws vector
     */
    const std::vector<SymbolicAlgebra::Expression<T>>& expressions() const
    {
        return _expressions;
    }

    /**
     * @brief Get the i-th variables
     * 
     * @return a reference to the i-th variables
     */
    const SymbolicAlgebra::Symbol<T>& variable(const unsigned int i) const
    {
        return _variables[i];
    }

    /**
     * @brief Get the i-th dynamic law
     * 
     * @return a reference to the i-th dynamic law
     */
    const SymbolicAlgebra::Expression<T>& expression(const unsigned int i) const
    {
        return _expressions[i];
    }

    /**
     * @brief Get the i-th dynamic law
     * 
     * @return a reference to the i-th dynamic law
     */
    const SymbolicAlgebra::Expression<T>& operator[](const SymbolicAlgebra::Symbol<T> &variable) const
    {
        return _expressions[_var_index.at[variable.get_id()]];
    }
};

/**
 * @brief A class to represent dicrete dynamic laws 
 * 
 * @tparam T is the type of expression value domain
 */
template<typename T>
class DiscreteDynamics : public Dynamics
{
public:
    /**
     * @brief An empty constructor
     */
    DiscreteDynamics(): Dynamics()
    {}

    /**
     * @brief A copy constructor
     *
     * @param[in] variables is the vector of the variables
     * @param[in] expressions is the vector of the dynamic law expressions
     */
    DiscreteDynamics(const std::vector<SymbolicAlgebra::Symbol<T>& variables,
                     const std::vector<SymbolicAlgebra::Expression<T>>& expressions):
        Dynamics(variable, expressions) 
    {}

    /**
     * @brief A swap constructor
     *
     * @param[in] variables is the vector of the variables
     * @param[in] expressions is the vector of the dynamic law expressions
     */
    DiscreteDynamics(std::vector<SymbolicAlgebra::Symbol<T>>&& variables,
                     std::vector<SymbolicAlgebra::Expression<T>>&& expressions):
        Dynamics(std::move(variable), std::move(expressions))
    {}


    /**
     * @brief A copy constructor
     *
     * @param[in] orig the original dynamic laws
     */
    DiscreteDynamics(const DiscreteDynamics<T>& orig): 
        Dynamics(orig) 
    {
    }

    /**
     * @brief A swap constructor
     *
     * @param[in] orig the original dynamic laws
     */
    DiscreteDynamics(DiscreteDynamics<T>&& orig):
        Dynamics(std::move(orig)) 
    {
    }

    /**
     * @brief Transform a parallelotope accorting with the dynamics
     * 
     * @param parallelotope the parallelotope to be transformed
     * @return an over-approximation of the parallelotope transformed
     *         by the dynamic laws
     */
    Parallelotope transform(const Parallelotope &parallelotope) const;

    /**
     * @brief Transform a bundle accorting with the dynamics
     * 
     * @param parallelotope the parallelotope to be transformed
     * @return an over-approximation of the bundle transformed
     *         by the dynamic laws
     */
    Bundle transform(const Bundle &bundle) const;

    
};

template<typename T>
class ContinousDynamics : public Dynamics
{
public:
    /**
     * @brief An empty constructor
     */
    ContinousDynamics(): Dynamics()
    {}

    /**
     * @brief A copy constructor
     *
     * @param[in] variables is the vector of the variables
     * @param[in] expressions is the vector of the dynamic law expressions
     */
    ContinousDynamics(const std::vector<SymbolicAlgebra::Symbol<T>& time_variable,
                      const std::vector<SymbolicAlgebra::Symbol<T>& variables,
                      const std::vector<SymbolicAlgebra::Expression<T>>& expressions):
        Dynamics(variable, expressions) 
    {}

    /**
     * @brief A swap constructor
     *
     * @param[in] variables is the vector of the variables
     * @param[in] expressions is the vector of the dynamic law expressions
     */
    ContinousDynamics(const std::vector<SymbolicAlgebra::Symbol<T>& time_variable,
                      std::vector<SymbolicAlgebra::Symbol<T>>&& variables,
                      std::vector<SymbolicAlgebra::Expression<T>>&& expressions):
        Dynamics(std::move(variable), std::move(expressions))
    {}


    /**
     * @brief A copy constructor
     *
     * @param[in] orig the original dynamic laws
     */
    ContinousDynamics(const ContinousDynamics<T>& orig): 
        Dynamics(orig) 
    {
    }

    /**
     * @brief A swap constructor
     *
     * @param[in] orig the original dynamic laws
     */
    ContinousDynamics(ContinousDynamics<T>&& orig):
        Dynamics(std::move(orig)) 
    {
    }

    Parallelotope elapse(const Parallelotope &parallelotope,
                         const T timestep) const;

    Bundle elapse(const Bundle &bundle,
                  const T timestep) const;

    Parallelotope flow(const Parallelotope &parallelotope,
                       const T timestep) const;

    Bundle flow(const Bundle &bundle,
                const T timestep) const;
};

#endif // DYNAMICS_H_