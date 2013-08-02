#pragma once

#include <unistd.h>
#include "blis.h"
#include <omp.h>

typedef unsigned int rank_t;
typedef unsigned int thread_count_t;
typedef omp_lock_t lock_t;


typedef struct
{
  void*   sent_object;
  thread_count_t num_threads;
  thread_count_t num_groups;
  
  bool_t  barrier_sense;
  lock_t  barrier_lock;
  thread_count_t   barrier_threads_arrived;
} thread_comm_t;

void    th_setup_comm( thread_comm_t *comm, 
				thread_count_t groups_this_level, thread_count_t groups_below );
void   th_broadcast( thread_comm_t *comm, Rank root, void *to_sendRecv, unsigned int size );
void    th_barrier( thread_comm_t *comm );
void    th_set_lock( lock_t *lock );
void    th_unset_lock( lock_t *lock );
void    th_init_lock( lock_t *lock );
void    th_destroy_lock( lock_t *lock );

bool_t  th_am_root( thread_comm_t *comm );

rank_t   th_thread_id( thread_comm_t *comm );
rank_t   th_group_id( thread_comm_t *comm );
rank_t   th_global_thread_id();

void th_shift_start_end(dim_t *start, dim_t *end, thread_comm_t *comm);
