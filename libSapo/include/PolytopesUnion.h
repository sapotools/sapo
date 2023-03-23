/**
 * @file PolytopesUnion.h
 * @author Alberto Casagrande <acasagrande@units.it>
 * @brief Polytope union functions
 * @version 0.1
 * @date 2022-05-28
 *
 * @copyright Copyright (c) 2022
 */

#ifndef POLYTOPESUNION_H_
#define POLYTOPESUNION_H_

#include "SetsUnion.h"
#include "Polytope.h"
#include "ThreadPool.h"

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
template<typename T, typename APPROX_TYPE>
SetsUnion<Polytope<T, APPROX_TYPE>> &
simplify(SetsUnion<Polytope<T, APPROX_TYPE>> &polytope_union);

/**
 * @brief Subtract two polytopes and close the result
 *
 * @tparam T is the numeric type used to represent polytopes
 * @tparam APPROX_TYPE is the type of approximation used in place of T
 *      during the computation
 * @param[in] P1 is a polytope
 * @param[in] P2 is a polytope
 * @return a union of polytopes obtained by closing the set
 *         \f$P1\setminus P2\f$
 */
template<typename T, typename APPROX_TYPE>
SetsUnion<Polytope<T, APPROX_TYPE>>
subtract_and_close(const Polytope<T, APPROX_TYPE> &P1,
                   const Polytope<T, APPROX_TYPE> &P2);

template<typename T, typename APPROX_TYPE>
SetsUnion<Polytope<T, APPROX_TYPE>> &
simplify(SetsUnion<Polytope<T, APPROX_TYPE>> &polytope_union)
{
#ifdef WITH_THREADS

  if (polytope_union.size() < 2) { // if the union consists in less
                                   // than two polytopes, use the
                                   // non-threaded version to avoid
                                   // the thread overhead
    for (auto it = std::begin(polytope_union); it != std::end(polytope_union);
         ++it) {
      it->simplify();
    }

    return polytope_union;
  }

  auto simplify_polytope = [](Polytope<T, APPROX_TYPE> &P) { P.simplify(); };

  ThreadPool::BatchId batch_id = thread_pool.create_batch();

  for (auto it = std::begin(polytope_union); it != std::end(polytope_union);
       ++it) {
    // submit the task to the thread pool
    thread_pool.submit_to_batch(batch_id, simplify_polytope, std::ref(*it));
  }

  // join to the pool threads
  thread_pool.join_threads(batch_id);

  // close the batch
  thread_pool.close_batch(batch_id);

#else  // WITH_THREADS
  for (auto it = std::begin(polytope_union); it != std::end(polytope_union);
       ++it) {
    it->simplify();
  }
#endif // WITH_THREADS

  return polytope_union;
}

template<typename T, typename APPROX_TYPE>
SetsUnion<Polytope<T, APPROX_TYPE>>
subtract_and_close(const Polytope<T, APPROX_TYPE> &P1,
                   const Polytope<T, APPROX_TYPE> &P2)
{
  using namespace LinearAlgebra;

  SetsUnion<Polytope<T, APPROX_TYPE>> su;

  if (is_true(P2.includes(P1))) {
    return su;
  }

  if (!is_false(are_disjoint(P1, P2))) {
    su.add(P1);

    return su;
  }

  for (unsigned int i = 0; i < P2.size(); ++i) {
    auto max = P1.maximize(P2.A(i)).objective_value();

    if (max > P2.b(i)) {
      Polytope new_P1 = P1;
      new_P1.add_constraint(-P2.A(i), -P2.b(i));

      su.add(std::move(new_P1));
    }
  }

  return su;
}

#endif /* POLYTOPESUNION_H_ */
