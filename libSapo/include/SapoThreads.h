/**
 * @file SapoThreads.h
 *
 * An header for the Sapo thread pool.
 *
 * @author Alberto Casagrande <acasagrande@units.it>
 * @version 0.1
 */

#ifndef SAPOTHREADS_H_
#define SAPOTHREADS_H_

#ifdef WITH_THREADS
#include "ThreadPool.h"

extern ThreadPool thread_pool;

#endif // WITH_THREADS

#endif // SAPOTHREADS_H_
