/**
 * @file PolytopesUnion.cpp
 * @author Alberto Casagrande (acasagrande@units.it)
 * @brief Simplify polytopes union representations
 * @version 0.1
 * @date 2022-05-28
 * 
 * @copyright Copyright (c) 2022
 */

#include "Polytope.h"
#include "SetsUnion.h"

SetsUnion<Polytope> &simplify(SetsUnion<Polytope> &polytope_union)
{
#ifdef WITH_THREADS

  if (polytope_union.size() < 2) { // if the union consists in less 
                                   // than two polytopes, use the
                                   // non-threaded version to avoid
                                   // the thread overhead
    for (auto it = std::begin(polytope_union); 
              it != std::end(polytope_union); ++it) {
      it->simplify();
    }

    return polytope_union;
  }

  auto simplify_polytope = [](Polytope &P) { P.simplify(); };

  ThreadPool::BatchId batch_id = thread_pool.create_batch();

  for (auto it = std::begin(polytope_union); 
            it != std::end(polytope_union); ++it) {
    // submit the task to the thread pool
    thread_pool.submit_to_batch(batch_id, simplify_polytope,
                                std::ref(*it));
  }

  // join to the pool threads
  thread_pool.join_threads(batch_id);

  // close the batch
  thread_pool.close_batch(batch_id);

#else  // WITH_THREADS
  for (auto it = std::begin(polytope_union); 
            it != std::end(polytope_union); ++it) {
    it->simplify();
  }
#endif // WITH_THREADS

  return polytope_union;
}