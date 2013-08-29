#include "support.h"
#include "math.h"
#include "string.h"

rank_t th_group_id( thread_comm_t *comm )
{
  return ( th_global_thread_id() / comm->num_threads_in_group )
    % comm->multiplicative_factor_above;
}

rank_t th_global_group_id( thread_comm_t *comm )
{
  return ( th_global_thread_id() / comm->num_threads_in_group );
}

rank_t th_thread_id( thread_comm_t *comm )
{
  return th_global_thread_id() % comm->num_threads_in_group;
}

rank_t th_global_thread_id()
{
    return omp_get_thread_num();
}

//I have armies of threads and need officers to organize them
bool_t th_am_root( thread_comm_t *comm )
{
    return th_thread_id( comm ) == 0;
}

void th_cleanup_communicator( thread_comm_t *comm )
{
    if( comm == NULL ) return;
    th_destroy_lock( &comm->barrier_lock );
}

void    th_setup_comm( thread_comm_t *comm, 
		       thread_count_t threads_in_group, 
		       thread_count_t multiplicative_factor_above )
{
    if( comm == NULL ) return;
    comm->sent_object = NULL;
    comm->num_threads_in_group = threads_in_group;
    comm->multiplicative_factor_above = multiplicative_factor_above;
    comm->barrier_sense = 0;
    th_init_lock( &comm->barrier_lock );
    comm->barrier_threads_arrived = 0;
}

void    th_release_comm( thread_comm_t *comm )
{
  if (comm)
    th_destroy_lock( &comm->barrier_lock );
}

void   th_broadcast( thread_comm_t *comm, rank_t root, void *to_sendRecv, unsigned int size )
{
  if (comm == NULL) return;
  th_broadcast_without_second_barrier(comm, root, to_sendRecv, size);
  th_barrier( comm );
}

void th_broadcast_without_second_barrier( thread_comm_t *comm, rank_t root, void *to_sendRecv, unsigned int size )
{   
  if( comm == NULL ) return;
    bool_t isRoot = th_thread_id( comm ) == root;
    if( isRoot )
        comm->sent_object = to_sendRecv;
    th_barrier( comm );
    if (!isRoot)
      memcpy( to_sendRecv, comm->sent_object, size );
}

void th_init_lock( lock_t *lock )
{
    omp_init_lock( lock );
}
void th_destroy_lock( lock_t *lock )
{
    omp_destroy_lock( lock );
}
void th_set_lock( lock_t *lock )
{
    omp_set_lock( lock );
}
void th_unset_lock( lock_t *lock )
{
    omp_unset_lock( lock );
}

//barrier routine taken from art of multicore programming or something
void th_barrier( thread_comm_t *comm )
{
    if(comm == NULL)
        return;
    bool_t my_sense = comm->barrier_sense;
    dim_t my_threads_arrived;
    th_set_lock(&comm->barrier_lock);
    my_threads_arrived = comm->barrier_threads_arrived + 1;
    comm->barrier_threads_arrived = my_threads_arrived;
    th_unset_lock(&comm->barrier_lock);
    if( my_threads_arrived == comm->num_threads_in_group ) {
        comm->barrier_threads_arrived = 0;
        comm->barrier_sense = !comm->barrier_sense;
        bool_t *mouth = &comm->barrier_sense;
        _Pragma("omp flush ( mouth )")
    }
    else {
        volatile bool_t *listener = &comm->barrier_sense;
        while( *listener == my_sense ) {}
    }
}

void th_shift_start_end(dim_t *start, dim_t *end, thread_comm_t *comm, dim_t round_factor)
{
  if (comm == NULL) return;
  if (comm->multiplicative_factor_above > 1) {
    rank_t group = th_group_id(comm);
    dim_t len = *end - *start;
    dim_t n_pt = len / comm->multiplicative_factor_above;
    n_pt = (n_pt * comm->multiplicative_factor_above < len) ? n_pt + 1 : n_pt;
    n_pt = (n_pt % round_factor == 0) ? n_pt : n_pt + round_factor - (n_pt % round_factor);
    *start = *start + group * n_pt;
    *end   = bli_min( *start + n_pt, *end );
  }
}
