/**
 * @file SapoThreads.h
 * @author Alberto Casagrande <acasagrande@units.it>
 * @brief An header for importing the Sapo thread pool
 * @version 0.1
 * @date 2022-02-17
 * 
 * @copyright Copyright (c) 2022
 */

#ifndef SAPOTHREADS_H_
#define SAPOTHREADS_H_

#ifdef WITH_THREADS
#include "ThreadPool.h"

extern ThreadPool thread_pool;

#endif // WITH_THREADS

#endif // SAPOTHREADS_H_
