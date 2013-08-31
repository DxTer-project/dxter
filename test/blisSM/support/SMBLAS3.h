#include "support.h"

#include "bli_trmm_ll_ker_var2_par.h"
#include "bli_trmm_lu_ker_var2_par.h"
#include "bli_trmm_rl_ker_var2_par.h"
#include "bli_trmm_ru_ker_var2_par.h"
//#include "bli_trsm_rl_ker_var2_par.h"
//#include "bli_trsm_ru_ker_var2_par.h"
//#include "bli_trsm_ll_ker_var2_par.h"
//#include "bli_trsm_lu_ker_var2_par.h"
#include "bli_packm_blk_var2_par.h"
#include "bli_packm_blk_var3_par.h"
#include "bli_gemm_ker_var2_par.h"
#include "bli_herk_l_ker_var2_par.h"
#include "bli_herk_u_ker_var2_par.h"


extern blksz_t *gemm_mc;
extern blksz_t *gemm_kc;
extern blksz_t *gemm_nc;
extern blksz_t *gemm_mr;
extern blksz_t *gemm_kr;
extern blksz_t *gemm_nr;
extern blksz_t *gemm_extmr;
extern blksz_t *gemm_extkr;
extern blksz_t *gemm_extnr;
extern blksz_t *trmm_mr;

#define NUMTHREADSPERL1 1
#define NUML1PERL2 1
#define NUML2PERPROC 2
#define NUMPROCS 2

#define NUML1 (NUMPROCS*NUML2PERPROC*NUML1PERL2)
#define NUML2 (NUMPROCS*NUML2PERPROC)
#define NUMTHREADS (NUML1*NUMTHREADSPERL1)
#define NUMTHREADSPERL2 (NUML1PERL2*NUMTHREADSPERL1)
#define NUMTHREADSPERPROC (NUMTHREADSPERL2*NUML2PERPROC)

void SetupComms(thread_comm_t *global_comm, thread_comm_t *proc_comms,
		thread_comm_t *l2_comms, thread_comm_t *l1_comms,
		thread_comm_t **GlobalComm, thread_comm_t **ProcComm,
		thread_comm_t **L2Comm, thread_comm_t **L1Comm)
{
  rank_t rank = omp_get_thread_num();
  *GlobalComm = global_comm;
  *ProcComm = proc_comms + (rank / NUMTHREADSPERPROC);
  *L2Comm = l2_comms + (rank / NUMTHREADSPERL2);
  *L1Comm = l1_comms + (rank / NUMTHREADSPERL1);

  if ((rank % (NUMTHREADSPERPROC)) == 0) {
    th_setup_comm(*ProcComm, NUMTHREADSPERPROC, NUMPROCS);
  }

  if ((rank % NUMTHREADSPERL2) == 0) {
    th_setup_comm(*L2Comm, NUMTHREADSPERL2, NUML2PERPROC);
  }

  if ((rank % NUMTHREADSPERL1) == 0) {
    th_setup_comm(*L1Comm, NUMTHREADSPERL1, NUML1PERL2);
  }

  th_barrier(*GlobalComm);
}

#define DEFINECOMMS \



#define FUNCTIONSTART \
  obj_t packed_A_blk;\
  obj_t packed_B_pan;\
  bli_obj_init_pack( &packed_A_blk );\
  bli_obj_init_pack( &packed_B_pan );\
  thread_comm_t *GlobalComm;\
  thread_comm_t *ProcComm;\
  thread_comm_t *L2Comm;\
  thread_comm_t *L1Comm;\
  SetupComms(global_comm, proc_comms,\
	     l2_comms, l1_comms, \
	     &GlobalComm, &ProcComm, \
	     &L2Comm, &L1Comm);

#define FUNCTIONEND \
  th_barrier(GlobalComm);\
  if (th_am_root(L1Comm)) {\
    th_release_comm(L1Comm);\
  }\
  if (th_am_root(L2Comm)) {\
    bli_obj_release_pack( &packed_A_blk );\
    th_release_comm(L2Comm);\
  }\
  if (th_am_root(ProcComm)) {\
    bli_obj_release_pack( &packed_B_pan );\
    th_release_comm(ProcComm);\
  }\
  if (th_am_root(GlobalComm)) {\
    th_release_comm(GlobalComm);\
  }



