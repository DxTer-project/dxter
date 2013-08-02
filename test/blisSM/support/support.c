#include "support.h"
#include "math.h"

rank_t th_group_id( thread_comm_t *comm )
{
  return floor( (double)th_global_thread_id() / comm->num_threads );
}

rank_t th_thread_id( thread_comm_t *comm )
{
    return th_global_thread_id() % comm->num_threads;
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

void th_setup_communicator( thread_comm_t *comm, dim_t num_threads)
{
    if( comm == NULL ) return;
    comm->sent_object = NULL;
    comm->num_threads = num_threads;
    comm->barrier_sense = 0;
    th_init_lock( &comm->barrier_lock );
    comm->barrier_threads_arrived = 0;
}

void   th_broadcast( thread_comm_t *comm, Rank root, void *to_sendRecv, unsigned int size );
{   
    if( comm == NULL ) return to_send;
    bool_t isRoot = th_thread_id( comm ) == root;
    if( isRoot )
        comm->sent_object = to_sendRecv;
    th_barrier( comm );
    if (!isRoot)
      memcpy( to_sendRecv, comm->sent_object, size );
    th_barrier( comm );
    return object;
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
    if( my_threads_arrived == comm->num_threads ) {
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

void th_shift_start_end(dim_t *start, dim_t *end, thread_comm_t *comm)
{
  rank_t rank = th_group_id( comm );
  dim_t len = *end - *start;
  dim_t n_pt = len / comm->num_groups;
  n_pt = (n_pt * comm->num_groups < len) ? n_pt + 1 : n_pt;
  n_pt = (n_pt % 4 == 0) ? n_pt : n_pt + 4 - (n_pt % 4);
  *start = *start + rank * n_pt;
  *end   = bli_min( *start + n_pt, *start + len );
}
