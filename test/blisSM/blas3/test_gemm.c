#include <unistd.h>
#include <omp.h>
#include "blis.h"
#include "SMBLAS3.h"

thread_comm_t global_comm[1];
thread_comm_t all_l2_comms[1];
thread_comm_t proc_comms[NUMPROCS];
thread_comm_t l2_comms[NUML2];
thread_comm_t l2_comms_sub[NUML2];
thread_comm_t l1_comms[NUML1];


void DxT_GemmNN( obj_t *alpha,
		 obj_t *A,
		 obj_t *B,
		 obj_t *beta,
		 obj_t *C )
{
  FUNCTIONSTART
  bli_scalm(beta, C);

  dim_t idx1, dimLen1, bs1;
  dimLen1 = bli_obj_width_after_trans( *C );
  idx1 = 0;
  th_shift_start_end(&idx1, &dimLen1, ProcComm, bli_blksz_for_obj( C, gemm_nr));
  for ( ; idx1 < dimLen1; idx1 += bs1 ) {
    bs1 = bli_determine_blocksize_f( idx1, dimLen1, C, gemm_nc );
    dim_t idx2, dimLen2, bs2;
    //****
    obj_t B_1;
    bli_acquire_mpart_l2r( BLIS_SUBPART1, idx1, bs1, B, &B_1 );
    obj_t C_1;
    bli_acquire_mpart_l2r( BLIS_SUBPART1, idx1, bs1, C, &C_1 );
    //------------------------------------//

    dimLen2 = bli_obj_width_after_trans( *A );
    for ( idx2 = 0; idx2 < dimLen2; idx2 += bs2 ) {
      bs2 = bli_determine_blocksize_f( idx2, dimLen2, A, gemm_kc );
      dim_t idx3, dimLen3, bs3;
      //****
      obj_t A_1;
      bli_acquire_mpart_l2r( BLIS_SUBPART1, idx2, bs2, A, &A_1 );
      obj_t B_1_1;
      bli_acquire_mpart_t2b( BLIS_SUBPART1, idx2, bs2, &B_1, &B_1_1 );
      //------------------------------------//

      th_barrier( ProcComm );
      if (th_am_root(ProcComm)) {
	alloced_B = TRUE;
	bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_COL_PANELS, 
			     BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			     BLIS_BUFFER_FOR_B_PANEL,
			     gemm_kr, gemm_nr, 
			     &B_1_1, &packed_B_pan_local_alloc );
	th_broadcast_without_second_barrier(ProcComm, 0,
					    (void*)(&packed_B_pan_local_alloc),
					    (void*)(&packed_B_pan), sizeof(packed_B_pan));
      }
      else {
	th_broadcast_without_second_barrier(ProcComm, 0, (void*)NULL,
					    (void*)(&packed_B_pan), sizeof(packed_B_pan));
      }
      bli_packm_blk_var2_par( &BLIS_ONE, &B_1_1, &packed_B_pan, ProcComm );
      //// ***Parallelized with communicator ProcComm
      dimLen3 = bli_obj_length_after_trans( C_1 );
      idx3 = 0;
      th_shift_start_end(&idx3, &dimLen3, L2Comm, bli_blksz_for_obj( &C_1, gemm_mr));
      for ( ; idx3 < dimLen3; idx3 += bs3 ) {
	bs3 = bli_determine_blocksize_f( idx3, dimLen3, &C_1, gemm_mc );
	dim_t idx4, dimLen4, bs4;
	//****
	obj_t A_1_1;
	bli_acquire_mpart_t2b( BLIS_SUBPART1, idx3, bs3, &A_1, &A_1_1 );
	obj_t C_1_1;
	bli_acquire_mpart_t2b( BLIS_SUBPART1, idx3, bs3, &C_1, &C_1_1 );
	//------------------------------------//

	th_barrier( L2Comm );
	if (th_am_root(L2Comm)) {
	  alloced_A = TRUE;
	  bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_ROW_PANELS, 
			       BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			       BLIS_BUFFER_FOR_A_BLOCK,
			       gemm_mr, gemm_kr, 
			       &A_1_1, &packed_A_blk_local_alloc );
	  th_broadcast_without_second_barrier(L2Comm, 0,
					      (void*)(&packed_A_blk_local_alloc),
					      (void*)(&packed_A_blk), sizeof(packed_A_blk));
	}
	else {
	  th_broadcast_without_second_barrier(L2Comm, 0, (void*)NULL,
					      (void*)(&packed_A_blk), sizeof(packed_A_blk));
	}
	bli_packm_blk_var2_par( &BLIS_ONE, &A_1_1, &packed_A_blk, L2Comm );
	bli_gemm_ker_var2_par( &BLIS_ONE, &packed_A_blk, &packed_B_pan, 
			       &BLIS_ONE, &C_1_1, (gemm_t*)NULL, L1Comm );

	//------------------------------------//

	//****
      }

      //------------------------------------//

      //****
    }

    //------------------------------------//

    //****
  }


  FUNCTIONEND
}

void DxT_GemmTN( obj_t *alpha,
	  obj_t *A,
	  obj_t *B,
	  obj_t *beta,
	  obj_t *C )
{
  FUNCTIONSTART
  bli_scalm(beta, C);

  dim_t idx1, dimLen1, bs1;
  dimLen1 = bli_obj_width_after_trans( *C );
  idx1 = 0;
  th_shift_start_end(&idx1, &dimLen1, ProcComm, bli_blksz_for_obj( C, gemm_nr));
  for ( ; idx1 < dimLen1; idx1 += bs1 ) {
    bs1 = bli_determine_blocksize_f( idx1, dimLen1, C, gemm_nc );
    dim_t idx2, dimLen2, bs2;
    //****
    obj_t B_1;
    bli_acquire_mpart_l2r( BLIS_SUBPART1, idx1, bs1, B, &B_1 );
    obj_t C_1;
    bli_acquire_mpart_l2r( BLIS_SUBPART1, idx1, bs1, C, &C_1 );
    //------------------------------------//

    dimLen2 = bli_obj_length_after_trans( *A );
    for ( idx2 = 0; idx2 < dimLen2; idx2 += bs2 ) {
      bs2 = bli_determine_blocksize_f( idx2, dimLen2, A, gemm_kc );
      dim_t idx3, dimLen3, bs3;
      //****
      obj_t A_1;
      bli_acquire_mpart_t2b( BLIS_SUBPART1, idx2, bs2, A, &A_1 );
      obj_t B_1_1;
      bli_acquire_mpart_t2b( BLIS_SUBPART1, idx2, bs2, &B_1, &B_1_1 );
      //------------------------------------//

      th_barrier( ProcComm );
      if (th_am_root(ProcComm)) {
	alloced_B = TRUE;
	bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_COL_PANELS, 
			     BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			     BLIS_BUFFER_FOR_B_PANEL,
			     gemm_kr, gemm_nr, 
			     &B_1_1, &packed_B_pan_local_alloc );
	th_broadcast_without_second_barrier(ProcComm, 0,
					    (void*)(&packed_B_pan_local_alloc),
					    (void*)(&packed_B_pan), sizeof(packed_B_pan));
      }
      else {
	th_broadcast_without_second_barrier(ProcComm, 0, (void*)NULL,
					    (void*)(&packed_B_pan), sizeof(packed_B_pan));
      }
      bli_packm_blk_var2_par( &BLIS_ONE, &B_1_1, &packed_B_pan, ProcComm );
      //// ***Parallelized with communicator ProcComm
      dimLen3 = bli_obj_length_after_trans( C_1 );
      idx3 = 0;
      th_shift_start_end(&idx3, &dimLen3, L2Comm, bli_blksz_for_obj( &C_1, gemm_mr));
      for ( ; idx3 < dimLen3; idx3 += bs3 ) {
	bs3 = bli_determine_blocksize_f( idx3, dimLen3, &C_1, gemm_mc );
	dim_t idx4, dimLen4, bs4;
	//****
	obj_t A_1_1;
	bli_acquire_mpart_l2r( BLIS_SUBPART1, idx3, bs3, &A_1, &A_1_1 );
	obj_t C_1_1;
	bli_acquire_mpart_t2b( BLIS_SUBPART1, idx3, bs3, &C_1, &C_1_1 );
	//------------------------------------//

	obj_t A_1_1T;
	bli_obj_alias_with_trans( BLIS_TRANSPOSE, A_1_1, A_1_1T);
	th_barrier( L2Comm );
	if (th_am_root(L2Comm)) {
	  alloced_A = TRUE;
	  bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_ROW_PANELS, 
			       BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			       BLIS_BUFFER_FOR_A_BLOCK,
			       gemm_mr, gemm_kr, 
			       &A_1_1T, &packed_A_blk_local_alloc );
	  th_broadcast_without_second_barrier(L2Comm, 0,
					      (void*)(&packed_A_blk_local_alloc),
					      (void*)(&packed_A_blk), sizeof(packed_A_blk));
	}
	else {
	  th_broadcast_without_second_barrier(L2Comm, 0, (void*)NULL,
					      (void*)(&packed_A_blk), sizeof(packed_A_blk));
	}
	bli_packm_blk_var2_par( &BLIS_ONE, &A_1_1T, &packed_A_blk, L2Comm );
	bli_gemm_ker_var2_par( &BLIS_ONE, &packed_A_blk, &packed_B_pan, 
			       &BLIS_ONE, &C_1_1, (gemm_t*)NULL, L1Comm );

	//------------------------------------//

	//****
      }

      //------------------------------------//

      //****
    }

    //------------------------------------//

    //****
  }

  FUNCTIONEND
}

void DxT_GemmNT( obj_t *alpha,
	  obj_t *A,
	  obj_t *B,
	  obj_t *beta,
	  obj_t *C )
{
  FUNCTIONSTART
  bli_scalm(beta, C);

  dim_t idx1, dimLen1, bs1;
  dimLen1 = bli_obj_width_after_trans( *C );
  idx1 = 0;
  th_shift_start_end(&idx1, &dimLen1, ProcComm, bli_blksz_for_obj( C, gemm_nr));
  for ( ; idx1 < dimLen1; idx1 += bs1 ) {
    bs1 = bli_determine_blocksize_f( idx1, dimLen1, C, gemm_nc );
    dim_t idx2, dimLen2, bs2;
    //****
    obj_t B_1;
    bli_acquire_mpart_t2b( BLIS_SUBPART1, idx1, bs1, B, &B_1 );
    obj_t C_1;
    bli_acquire_mpart_l2r( BLIS_SUBPART1, idx1, bs1, C, &C_1 );
    //------------------------------------//

    dimLen2 = bli_obj_width_after_trans( *A );
    for ( idx2 = 0; idx2 < dimLen2; idx2 += bs2 ) {
      bs2 = bli_determine_blocksize_f( idx2, dimLen2, A, gemm_kc );
      dim_t idx3, dimLen3, bs3;
      //****
      obj_t A_1;
      bli_acquire_mpart_l2r( BLIS_SUBPART1, idx2, bs2, A, &A_1 );
      obj_t B_1_1;
      bli_acquire_mpart_l2r( BLIS_SUBPART1, idx2, bs2, &B_1, &B_1_1 );
      //------------------------------------//

      obj_t B_1_1T;
      bli_obj_alias_with_trans( BLIS_TRANSPOSE, B_1_1, B_1_1T);
      th_barrier( ProcComm );
      if (th_am_root(ProcComm)) {
	alloced_B = TRUE;
	bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_COL_PANELS, 
			     BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			     BLIS_BUFFER_FOR_B_PANEL,
			     gemm_kr, gemm_nr, 
			     &B_1_1T, &packed_B_pan_local_alloc );
	th_broadcast_without_second_barrier(ProcComm, 0,
					    (void*)(&packed_B_pan_local_alloc),
					    (void*)(&packed_B_pan), sizeof(packed_B_pan));
      }
      else {
	th_broadcast_without_second_barrier(ProcComm, 0, (void*)NULL,
					    (void*)(&packed_B_pan), sizeof(packed_B_pan));
      }
      bli_packm_blk_var2_par( &BLIS_ONE, &B_1_1T, &packed_B_pan, ProcComm );
      //// ***Parallelized with communicator ProcComm
      dimLen3 = bli_obj_length_after_trans( C_1 );
      idx3 = 0;
      th_shift_start_end(&idx3, &dimLen3, L2Comm, bli_blksz_for_obj( &C_1, gemm_mr));
      for ( ; idx3 < dimLen3; idx3 += bs3 ) {
	bs3 = bli_determine_blocksize_f( idx3, dimLen3, &C_1, gemm_mc );
	dim_t idx4, dimLen4, bs4;
	//****
	obj_t A_1_1;
	bli_acquire_mpart_t2b( BLIS_SUBPART1, idx3, bs3, &A_1, &A_1_1 );
	obj_t C_1_1;
	bli_acquire_mpart_t2b( BLIS_SUBPART1, idx3, bs3, &C_1, &C_1_1 );
	//------------------------------------//

	th_barrier( L2Comm );
	if (th_am_root(L2Comm)) {
	  alloced_A = TRUE;
	  bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_ROW_PANELS, 
			       BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			       BLIS_BUFFER_FOR_A_BLOCK,
			       gemm_mr, gemm_kr, 
			       &A_1_1, &packed_A_blk_local_alloc );
	  th_broadcast_without_second_barrier(L2Comm, 0,
					      (void*)(&packed_A_blk_local_alloc),
					      (void*)(&packed_A_blk), sizeof(packed_A_blk));
	}
	else {
	  th_broadcast_without_second_barrier(L2Comm, 0, (void*)NULL,
					      (void*)(&packed_A_blk), sizeof(packed_A_blk));
	}
	bli_packm_blk_var2_par( &BLIS_ONE, &A_1_1, &packed_A_blk, L2Comm );
	bli_gemm_ker_var2_par( &BLIS_ONE, &packed_A_blk, &packed_B_pan, 
			       &BLIS_ONE, &C_1_1, (gemm_t*)NULL, L1Comm );

	//------------------------------------//

	//****
      }

      //------------------------------------//

      //****
    }

    //------------------------------------//

    //****
  }


  FUNCTIONEND
}

void DxT_GemmTT( obj_t *alpha,
	  obj_t *A,
	  obj_t *B,
	  obj_t *beta,
	  obj_t *C )
{
  FUNCTIONSTART
  bli_scalm(beta, C);

  dim_t idx1, dimLen1, bs1;
  dimLen1 = bli_obj_width_after_trans( *C );
  idx1 = 0;
  th_shift_start_end(&idx1, &dimLen1, ProcComm, bli_blksz_for_obj( C, gemm_nr));
  for ( ; idx1 < dimLen1; idx1 += bs1 ) {
    bs1 = bli_determine_blocksize_f( idx1, dimLen1, C, gemm_nc );
    dim_t idx2, dimLen2, bs2;
    //****
    obj_t B_1;
    bli_acquire_mpart_t2b( BLIS_SUBPART1, idx1, bs1, B, &B_1 );
    obj_t C_1;
    bli_acquire_mpart_l2r( BLIS_SUBPART1, idx1, bs1, C, &C_1 );
    //------------------------------------//

    dimLen2 = bli_obj_length_after_trans( *A );
    for ( idx2 = 0; idx2 < dimLen2; idx2 += bs2 ) {
      bs2 = bli_determine_blocksize_f( idx2, dimLen2, A, gemm_kc );
      dim_t idx3, dimLen3, bs3;
      //****
      obj_t A_1;
      bli_acquire_mpart_t2b( BLIS_SUBPART1, idx2, bs2, A, &A_1 );
      obj_t B_1_1;
      bli_acquire_mpart_l2r( BLIS_SUBPART1, idx2, bs2, &B_1, &B_1_1 );
      //------------------------------------//

      obj_t B_1_1T;
      bli_obj_alias_with_trans( BLIS_TRANSPOSE, B_1_1, B_1_1T);
      th_barrier( ProcComm );
      if (th_am_root(ProcComm)) {
	alloced_B = TRUE;
	bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_COL_PANELS, 
			     BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			     BLIS_BUFFER_FOR_B_PANEL,
			     gemm_kr, gemm_nr, 
			     &B_1_1T, &packed_B_pan_local_alloc );
	th_broadcast_without_second_barrier(ProcComm, 0,
					    (void*)(&packed_B_pan_local_alloc),
					    (void*)(&packed_B_pan), sizeof(packed_B_pan));
      }
      else {
	th_broadcast_without_second_barrier(ProcComm, 0, (void*)NULL,
					    (void*)(&packed_B_pan), sizeof(packed_B_pan));
      }
      bli_packm_blk_var2_par( &BLIS_ONE, &B_1_1T, &packed_B_pan, ProcComm );
      //// ***Parallelized with communicator ProcComm
      dimLen3 = bli_obj_length_after_trans( C_1 );
      idx3 = 0;
      th_shift_start_end(&idx3, &dimLen3, L2Comm, bli_blksz_for_obj( &C_1, gemm_mr));
      for ( ; idx3 < dimLen3; idx3 += bs3 ) {
	bs3 = bli_determine_blocksize_f( idx3, dimLen3, &C_1, gemm_mc );
	dim_t idx4, dimLen4, bs4;
	//****
	obj_t A_1_1;
	bli_acquire_mpart_l2r( BLIS_SUBPART1, idx3, bs3, &A_1, &A_1_1 );
	obj_t C_1_1;
	bli_acquire_mpart_t2b( BLIS_SUBPART1, idx3, bs3, &C_1, &C_1_1 );
	//------------------------------------//

	obj_t A_1_1T;
	bli_obj_alias_with_trans( BLIS_TRANSPOSE, A_1_1, A_1_1T);
	th_barrier( L2Comm );
	if (th_am_root(L2Comm)) {
	  alloced_A = TRUE;
	  bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_ROW_PANELS, 
			       BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			       BLIS_BUFFER_FOR_A_BLOCK,
			       gemm_mr, gemm_kr, 
			       &A_1_1T, &packed_A_blk_local_alloc );
	  th_broadcast_without_second_barrier(L2Comm, 0,
					      (void*)(&packed_A_blk_local_alloc),
					      (void*)(&packed_A_blk), sizeof(packed_A_blk));
	}
	else {
	  th_broadcast_without_second_barrier(L2Comm, 0, (void*)NULL,
					      (void*)(&packed_A_blk), sizeof(packed_A_blk));
	}
	bli_packm_blk_var2_par( &BLIS_ONE, &A_1_1T, &packed_A_blk, L2Comm );
	bli_gemm_ker_var2_par( &BLIS_ONE, &packed_A_blk, &packed_B_pan, 
			       &BLIS_ONE, &C_1_1, (gemm_t*)NULL, L1Comm );

	//------------------------------------//

	//****
      }

      //------------------------------------//

      //****
    }

    //------------------------------------//

    //****
  }

  
  FUNCTIONEND
}

int main( int argc, char** argv )
{
  obj_t a, b, c1, c2;
  obj_t c_save;
  obj_t alpha, beta, negOne, normVal;
  dim_t m, n, k;
  dim_t p;
  dim_t p_begin, p_end, p_inc;
  int   m_input, n_input, k_input;
  num_t dt_a, dt_b, dt_c;
  num_t dt_alpha, dt_beta;
  int   r, n_repeats;

  double dtimeDxT;
  double dtime_saveDxT;
  double gflopsDxT;
  double dtimeBLI;
  double dtime_saveBLI;
  double gflopsBLI;

  int transA=0, transB=0;

  dt_alpha = BLIS_DOUBLE;
  dt_beta = BLIS_DOUBLE;

  bli_obj_create( dt_alpha, 1, 1, 0, 0, &negOne );
  bli_setsc(  -1.0, 0.0, &negOne );		
			
  if (argc != 3) {
    printf("test T/N T/N\n");
    return 0;
  }

  if (*(argv[1]) == 'T')
    transA = 1;
  else if (*(argv[1]) != 'N') {
    printf("transA not correct\n");
    return 0;
  }

  if (*(argv[2]) == 'T')
    transB = 1;
  else if (*(argv[2]) != 'N') {
    printf("transB not correct\n");
    return 0;
  }

  bli_init();

  n_repeats = 1;

  p_begin = 40;
  p_end   = 600;
  p_inc   = 40;

  m_input = -1;
  //m_input = 384;
  n_input = -1;
  k_input = -1;
  //k_input = 200;

  dt_a = BLIS_DOUBLE;
  dt_b = BLIS_DOUBLE;
  dt_c = BLIS_DOUBLE;
  dt_alpha = BLIS_DOUBLE;
  dt_beta = BLIS_DOUBLE;

  printf("%% %d Procs, %d L2 Per Proc, %d L1 Per L2, %d Per L1\n",
	 NUMPROCS, NUML2PERPROC, NUML1PERL2, NUMTHREADSPERL1);

	  bli_error_checking_level_set( BLIS_NO_ERROR_CHECKING );

  for ( p = p_begin; p <= p_end; p += p_inc )
    {

      if ( m_input < 0 ) m = p * ( dim_t )abs(m_input);
      else               m =     ( dim_t )    m_input;
      if ( n_input < 0 ) n = p * ( dim_t )abs(n_input);
      else               n =     ( dim_t )    n_input;
      if ( k_input < 0 ) k = p * ( dim_t )abs(k_input);
      else               k =     ( dim_t )    k_input;


      bli_obj_create( dt_alpha, 1, 1, 0, 0, &alpha );
      bli_obj_create( dt_beta,  1, 1, 0, 0, &beta );

      bli_obj_create( dt_a, m, k, 0, 0, &a );
      bli_obj_create( dt_b, k, n, 0, 0, &b );
      bli_obj_create( dt_c, m, n, 0, 0, &c1 );
      bli_obj_create( dt_c, m, n, 0, 0, &c2 );
      bli_obj_create( dt_c, m, n, 0, 0, &c_save );

      bli_randm( &a );
      bli_randm( &b );
      bli_randm( &c1 );


      bli_setsc( 1.0, 0.0, &alpha );
      bli_setsc( 1.0/1.0, 0.0, &beta );

      bli_copym( &c1, &c_save );
	
      dtime_saveDxT = 1.0e9;
      dtime_saveBLI = 1.0e9;

      for ( r = 0; r < n_repeats; ++r )
	{
	  bli_copym( &c_save, &c1 );
	  bli_copym( &c_save, &c2 );




	  dtimeDxT = bli_clock();

	  th_setup_comm(&global_comm[0], NUMTHREADS, 1);

	  _Pragma( "omp parallel num_threads(NUMTHREADS)" ) 
	    {
	      if (!transA && !transB) {

		{
		  DxT_GemmNN( &alpha,
			      &a,
			      &b,
			      &beta,
			      &c2 );
		}
			
	      }
	      else if (transA && !transB) {
		DxT_GemmTN( &alpha,
			    &a,
			    &b,
			    &beta,
			    &c2 );

	      }
	      else if (!transA) {// &&transB 
		DxT_GemmNT( &alpha,
			    &a,
			    &b,
			    &beta,
			    &c2 );
	      }
	      else {//transA && transB
		DxT_GemmTT( &alpha,
			    &a,
			    &b,
			    &beta,
			    &c2 );
	      }
	    }

	  dtime_saveDxT = bli_clock_min_diff( dtime_saveDxT, dtimeDxT );


	  dtimeBLI = bli_clock();
	  if (transA)
	    bli_obj_set_conjtrans( BLIS_TRANSPOSE, a);
	  if (transB)
	    bli_obj_set_conjtrans( BLIS_TRANSPOSE, b);
	  bli_gemm( &alpha,
		    &a,
		    &b,
		    &beta,
		    &c1);
	  if (transA)
	    bli_obj_toggle_trans( a );
	  if (transB)
	    bli_obj_toggle_trans( b );

	  dtime_saveBLI = bli_clock_min_diff( dtime_saveBLI, dtimeBLI );
	  
	  //	  bli_printm( "c_save", &c_save, "%4.1f", "" );
	  //	  bli_printm( "c1", &c1, "%4.1f", "" );
	  //	  bli_printm( "c2", &c2, "%4.1f", "" );

	  bli_axpym( &negOne, &c1, &c2 );

	  bli_obj_create( dt_alpha, 1, 1, 0, 0, &normVal );

	  bli_fnormm( &c2, &normVal );

	}

      gflopsDxT = ( 2.0 * m * k * n ) / ( dtime_saveDxT * 1.0e9 );
      gflopsBLI = ( 2.0 * m * k * n ) / ( dtime_saveBLI * 1.0e9 );

      printf( "data_gemm_DxT" );
      printf( "( %2ld, 1:5 ) = [ %4lu %4lu  %10.3e  %6.3f %2.4f ];\n",
	      (p - p_begin + 1)/p_inc + 1, m, n, dtime_saveDxT, gflopsDxT, gflopsDxT/gflopsBLI );

      printf( "data_gemm_blis" );
      printf( "( %2ld, 1:5 ) = [ %4lu %4lu %4lu  %10.3e  %6.3f ];\n",
	      (p - p_begin + 1)/p_inc + 1, m, k, n, dtime_saveBLI, gflopsBLI );

      bli_printm( "%NORM", &normVal, "%4.1f", "" );

      bli_obj_free( &alpha );
      bli_obj_free( &beta );

      bli_obj_free( &a );
      bli_obj_free( &b );
      bli_obj_free( &c1 );
      bli_obj_free( &c2 );
      bli_obj_free( &c_save );
    }

  bli_finalize();

  return 0;
}

