#include <unistd.h>
#include "blis.h"
#include "SMBLAS3.h"

#define TESTLIB 

thread_comm_t global_comm[1];			
thread_comm_t all_l2_comms[1];			
thread_comm_t proc_comms[NUMPROCS];		
thread_comm_t l2_comms[NUML2];
thread_comm_t l2_comms_sub[NUML2];
thread_comm_t l1_comms[NUML1];

void DxT_TrsmLLN( obj_t *alpha,
                 obj_t *L,
                 obj_t *X )
{
  FUNCTIONSTART
  bli_scalm(alpha, X);

  //// ***Parallelized with communicator GlobalComm
  dim_t idx1, dimLen1, bs1;
  dimLen1 = bli_obj_width_after_trans( *X );
  idx1 = 0;
  th_shift_start_end(&idx1, &dimLen1, ProcComm, bli_blksz_for_obj( X, gemm_nr));
  for ( ; idx1 < dimLen1; idx1 += bs1 ) {
    bs1 = bli_determine_blocksize_f( idx1, dimLen1, X, gemm_nc );
    dim_t idx2, dimLen2, bs2;
    //****
    obj_t X_1;
    bli_acquire_mpart_l2r( BLIS_SUBPART1, idx1, bs1, X, &X_1 );
    //------------------------------------//

    dimLen2 = bli_obj_length_after_trans( X_1 );
    for ( idx2 = 0; idx2 < dimLen2; idx2 += bs2 ) {
      bs2 = bli_determine_blocksize_f( idx2, dimLen2, &X_1, gemm_kc );
      dim_t idx3, dimLen3, bs3;
      //****
      obj_t L_11;
      bli_acquire_mpart_tl2br( BLIS_SUBPART11, idx2, bs2, L, &L_11 );
      obj_t L_21;
      bli_acquire_mpart_tl2br( BLIS_SUBPART21, idx2, bs2, L, &L_21 );
      obj_t X_1_1;
      bli_acquire_mpart_t2b( BLIS_SUBPART1, idx2, bs2, &X_1, &X_1_1 );
      obj_t X_1_2;
      bli_acquire_mpart_t2b( BLIS_SUBPART2, idx2, bs2, &X_1, &X_1_2 );
      //------------------------------------//

      th_barrier( ProcComm );
      if (th_am_root(ProcComm)) {
	alloced_B = TRUE;
	bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_COL_PANELS, 
			     BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			     BLIS_BUFFER_FOR_B_PANEL,
			     gemm_mr, gemm_nr, 
			     &X_1_1, &packed_B_pan_local_alloc );
	th_broadcast_without_second_barrier(ProcComm, 0,
					    (void*)(&packed_B_pan_local_alloc),
					    (void*)(&packed_B_pan), sizeof(packed_B_pan));
      }
      else {
	th_broadcast_without_second_barrier(ProcComm, 0, (void*)NULL,
					    (void*)(&packed_B_pan), sizeof(packed_B_pan));
      }
      bli_packm_blk_var2_par( &BLIS_ONE, &X_1_1, &packed_B_pan, ProcComm );
      dimLen3 = bli_obj_length_after_trans( L_11 );
      for ( idx3 = 0; idx3 < dimLen3; idx3 += bs3 ) {
	bs3 = bli_determine_blocksize_f( idx3, dimLen3, &L_11, gemm_mc );
	dim_t idx4, dimLen4, bs4;
	//****
	obj_t L_11_1;
	bli_acquire_mpart_t2b( BLIS_SUBPART1, idx3, bs3, &L_11, &L_11_1 );
	obj_t X_1_1_1;
	bli_acquire_mpart_t2b( BLIS_SUBPART1, idx3, bs3, &X_1_1, &X_1_1_1 );
	//------------------------------------//

	bli_obj_set_struc( BLIS_TRIANGULAR, L_11_1 );
	bli_obj_set_uplo( BLIS_LOWER, L_11_1 );
	th_barrier( ProcComm );
	if (th_am_root(ProcComm)) {
	  alloced_A = TRUE;
	  bli_packm_init_pack( TRUE, BLIS_INVERT_DIAG, BLIS_PACKED_ROW_PANELS, 
			       BLIS_PACK_REV_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			       BLIS_BUFFER_FOR_A_BLOCK,
			       gemm_mr, gemm_mr, 
			       &L_11_1, &packed_A_blk_local_alloc );
	  th_broadcast_without_second_barrier(ProcComm, 0,
					      (void*)(&packed_A_blk_local_alloc),
					      (void*)(&packed_A_blk), sizeof(packed_A_blk));
	}
	else {
	  th_broadcast_without_second_barrier(ProcComm, 0, (void*)NULL,
					      (void*)(&packed_A_blk), sizeof(packed_A_blk));
	}
	bli_packm_blk_var3_par( &BLIS_ONE, &L_11_1, &packed_A_blk, ProcComm );
	bli_trsm_ll_ker_var2_par( &BLIS_ONE, &packed_A_blk, &packed_B_pan, 
				  &BLIS_ZERO, &X_1_1_1, (trsm_t*)NULL, ProcComm);

	//------------------------------------//

	//****
      }
      th_barrier( ProcComm ); //barrier for dependency
      //// ***Parallelized with communicator ProcComm
      dimLen3 = bli_obj_length_after_trans( X_1_2 );
      idx3 = 0;
      th_shift_start_end(&idx3, &dimLen3, L2Comm, bli_blksz_for_obj( &X_1_2, gemm_mr));
      for ( ; idx3 < dimLen3; idx3 += bs3 ) {
	bs3 = bli_determine_blocksize_f( idx3, dimLen3, &X_1_2, gemm_mc );
	dim_t idx4, dimLen4, bs4;
	//****
	obj_t L_21_1;
	bli_acquire_mpart_t2b( BLIS_SUBPART1, idx3, bs3, &L_21, &L_21_1 );
	obj_t X_1_2_1;
	bli_acquire_mpart_t2b( BLIS_SUBPART1, idx3, bs3, &X_1_2, &X_1_2_1 );
	//------------------------------------//

	th_barrier( L2Comm );
	if (th_am_root(L2Comm)) {
	  alloced_A = TRUE;
	  bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_ROW_PANELS, 
			       BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			       BLIS_BUFFER_FOR_A_BLOCK,
			       gemm_mr, gemm_kr, 
			       &L_21_1, &packed_A_blk_local_alloc );
	  th_broadcast_without_second_barrier(L2Comm, 0,
					      (void*)(&packed_A_blk_local_alloc),
					      (void*)(&packed_A_blk), sizeof(packed_A_blk));
	}
	else {
	  th_broadcast_without_second_barrier(L2Comm, 0, (void*)NULL,
					      (void*)(&packed_A_blk), sizeof(packed_A_blk));
	}
	bli_packm_blk_var2_par( &BLIS_ONE, &L_21_1, &packed_A_blk, L2Comm );
	bli_gemm_ker_var2_par( &BLIS_MINUS_ONE, &packed_A_blk, &packed_B_pan, 
			       &BLIS_ONE, &X_1_2_1, (gemm_t*)NULL, L1Comm );

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

void DxT_TrsmLLT( obj_t *alpha,
		  obj_t *L,
		  obj_t *X )
{
  FUNCTIONSTART
  bli_scalm(alpha, X);

  dim_t idx1, dimLen1, bs1;
  dimLen1 = bli_obj_width_after_trans( *X );
  idx1 = 0;
  th_shift_start_end(&idx1, &dimLen1, ProcComm, bli_blksz_for_obj( X, gemm_nr));
  for ( ; idx1 < dimLen1; idx1 += bs1 ) {
    bs1 = bli_determine_blocksize_f( idx1, dimLen1, X, gemm_nc );
    dim_t idx2, dimLen2, bs2;
    //****
    obj_t X_1;
    bli_acquire_mpart_l2r( BLIS_SUBPART1, idx1, bs1, X, &X_1 );
    //------------------------------------//

    dimLen2 = bli_obj_length_after_trans( X_1 );
    for ( idx2 = 0; idx2 < dimLen2; idx2 += bs2 ) {
      bs2 = bli_determine_blocksize_b( idx2, dimLen2, &X_1, gemm_kc );
      dim_t idx3, dimLen3, bs3;
      //****
      obj_t L_10;
      bli_acquire_mpart_br2tl( BLIS_SUBPART10, idx2, bs2, L, &L_10 );
      obj_t L_11;
      bli_acquire_mpart_br2tl( BLIS_SUBPART11, idx2, bs2, L, &L_11 );
      obj_t X_1_0;
      bli_acquire_mpart_b2t( BLIS_SUBPART0, idx2, bs2, &X_1, &X_1_0 );
      obj_t X_1_1;
      bli_acquire_mpart_b2t( BLIS_SUBPART1, idx2, bs2, &X_1, &X_1_1 );
      //------------------------------------//

      th_barrier( ProcComm );
      if (th_am_root(ProcComm)) {
	alloced_B = TRUE;
	bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_COL_PANELS, 
			     BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			     BLIS_BUFFER_FOR_B_PANEL,
			     gemm_mr, gemm_nr, 
			     &X_1_1, &packed_B_pan_local_alloc );
	th_broadcast_without_second_barrier(ProcComm, 0,
					    (void*)(&packed_B_pan_local_alloc),
					    (void*)(&packed_B_pan), sizeof(packed_B_pan));
      }
      else {
	th_broadcast_without_second_barrier(ProcComm, 0, (void*)NULL,
					    (void*)(&packed_B_pan), sizeof(packed_B_pan));
      }
      bli_packm_blk_var2_par( &BLIS_ONE, &X_1_1, &packed_B_pan, ProcComm );
      dimLen3 = bli_obj_width_after_trans( L_11 );
      for ( idx3 = 0; idx3 < dimLen3; idx3 += bs3 ) {
	bs3 = bli_determine_blocksize_b( idx3, dimLen3, &L_11, gemm_mc );
	dim_t idx4, dimLen4, bs4;
	//****
	obj_t L_11_1;
	bli_acquire_mpart_r2l( BLIS_SUBPART1, idx3, bs3, &L_11, &L_11_1 );
	obj_t X_1_1_1;
	bli_acquire_mpart_b2t( BLIS_SUBPART1, idx3, bs3, &X_1_1, &X_1_1_1 );
	//------------------------------------//

	obj_t L_11_1T;
	bli_obj_alias_with_trans( BLIS_TRANSPOSE, L_11_1, L_11_1T);
	bli_obj_set_struc( BLIS_TRIANGULAR, L_11_1T );
	bli_obj_set_uplo( BLIS_LOWER, L_11_1T );
	th_barrier( ProcComm );
	if (th_am_root(ProcComm)) {
	  alloced_A = TRUE;
	  bli_packm_init_pack( TRUE, BLIS_INVERT_DIAG, BLIS_PACKED_ROW_PANELS, 
			       BLIS_PACK_REV_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			       BLIS_BUFFER_FOR_A_BLOCK,
			       gemm_mr, gemm_mr, 
			       &L_11_1T, &packed_A_blk_local_alloc );
	  th_broadcast_without_second_barrier(ProcComm, 0,
					      (void*)(&packed_A_blk_local_alloc),
					      (void*)(&packed_A_blk), sizeof(packed_A_blk));
	}
	else {
	  th_broadcast_without_second_barrier(ProcComm, 0, (void*)NULL,
					      (void*)(&packed_A_blk), sizeof(packed_A_blk));
	}
	bli_packm_blk_var3_par( &BLIS_ONE, &L_11_1T, &packed_A_blk, ProcComm );
	bli_trsm_lu_ker_var2_par( &BLIS_ONE, &packed_A_blk, &packed_B_pan, 
				  &BLIS_ZERO, &X_1_1_1, (trsm_t*)NULL, ProcComm);

	//------------------------------------//

	//****
      }
      th_barrier( ProcComm );	 //barrier for dependency
      //// ***Parallelized with communicator ProcComm
      dimLen3 = bli_obj_length_after_trans( X_1_0 );
      idx3 = 0;
      th_shift_start_end(&idx3, &dimLen3, L2Comm, bli_blksz_for_obj( &X_1_0, gemm_mr));
      for ( ; idx3 < dimLen3; idx3 += bs3 ) {
	bs3 = bli_determine_blocksize_f( idx3, dimLen3, &X_1_0, gemm_mc );
	dim_t idx4, dimLen4, bs4;
	//****
	obj_t L_10_1;
	bli_acquire_mpart_l2r( BLIS_SUBPART1, idx3, bs3, &L_10, &L_10_1 );
	obj_t X_1_0_1;
	bli_acquire_mpart_t2b( BLIS_SUBPART1, idx3, bs3, &X_1_0, &X_1_0_1 );
	//------------------------------------//

	obj_t L_10_1T;
	bli_obj_alias_with_trans( BLIS_TRANSPOSE, L_10_1, L_10_1T);
	th_barrier( L2Comm );
	if (th_am_root(L2Comm)) {
	  alloced_A = TRUE;
	  bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_ROW_PANELS, 
			       BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			       BLIS_BUFFER_FOR_A_BLOCK,
			       gemm_mr, gemm_kr, 
			       &L_10_1T, &packed_A_blk_local_alloc );
	  th_broadcast_without_second_barrier(L2Comm, 0,
					      (void*)(&packed_A_blk_local_alloc),
					      (void*)(&packed_A_blk), sizeof(packed_A_blk));
	}
	else {
	  th_broadcast_without_second_barrier(L2Comm, 0, (void*)NULL,
					      (void*)(&packed_A_blk), sizeof(packed_A_blk));
	}
	bli_packm_blk_var2_par( &BLIS_ONE, &L_10_1T, &packed_A_blk, L2Comm );
	bli_gemm_ker_var2_par( &BLIS_MINUS_ONE, &packed_A_blk, &packed_B_pan, 
			       &BLIS_ONE, &X_1_0_1, (gemm_t*)NULL, L1Comm );

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


void DxT_TrsmLUN( obj_t *alpha,
		  obj_t *U,
		  obj_t *X )
{
  FUNCTIONSTART
  bli_scalm(alpha, X);

  dim_t idx1, dimLen1, bs1;
  dimLen1 = bli_obj_width_after_trans( *X );
  idx1 = 0;
  th_shift_start_end(&idx1, &dimLen1, ProcComm, bli_blksz_for_obj( X, gemm_nr));
  for ( ; idx1 < dimLen1; idx1 += bs1 ) {
    bs1 = bli_determine_blocksize_f( idx1, dimLen1, X, gemm_nc );
    dim_t idx2, dimLen2, bs2;
    //****
    obj_t X_1;
    bli_acquire_mpart_l2r( BLIS_SUBPART1, idx1, bs1, X, &X_1 );
    //------------------------------------//

    dimLen2 = bli_obj_length_after_trans( X_1 );
    for ( idx2 = 0; idx2 < dimLen2; idx2 += bs2 ) {
      bs2 = bli_determine_blocksize_b( idx2, dimLen2, &X_1, gemm_kc );
      dim_t idx3, dimLen3, bs3;
      //****
      obj_t U_01;
      bli_acquire_mpart_br2tl( BLIS_SUBPART01, idx2, bs2, U, &U_01 );
      obj_t U_11;
      bli_acquire_mpart_br2tl( BLIS_SUBPART11, idx2, bs2, U, &U_11 );
      obj_t X_1_0;
      bli_acquire_mpart_b2t( BLIS_SUBPART0, idx2, bs2, &X_1, &X_1_0 );
      obj_t X_1_1;
      bli_acquire_mpart_b2t( BLIS_SUBPART1, idx2, bs2, &X_1, &X_1_1 );
      //------------------------------------//

      th_barrier( ProcComm );
      if (th_am_root(ProcComm)) {
	alloced_B = TRUE;
	bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_COL_PANELS, 
			     BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			     BLIS_BUFFER_FOR_B_PANEL,
			     gemm_mr, gemm_nr, 
			     &X_1_1, &packed_B_pan_local_alloc );
	th_broadcast_without_second_barrier(ProcComm, 0,
					    (void*)(&packed_B_pan_local_alloc),
					    (void*)(&packed_B_pan), sizeof(packed_B_pan));
      }
      else {
	th_broadcast_without_second_barrier(ProcComm, 0, (void*)NULL,
					    (void*)(&packed_B_pan), sizeof(packed_B_pan));
      }
      bli_packm_blk_var2_par( &BLIS_ONE, &X_1_1, &packed_B_pan, ProcComm );
      dimLen3 = bli_obj_length_after_trans( U_11 );
      for ( idx3 = 0; idx3 < dimLen3; idx3 += bs3 ) {
	bs3 = bli_determine_blocksize_b( idx3, dimLen3, &U_11, gemm_mc );
	dim_t idx4, dimLen4, bs4;
	//****
	obj_t U_11_1;
	bli_acquire_mpart_b2t( BLIS_SUBPART1, idx3, bs3, &U_11, &U_11_1 );
	obj_t X_1_1_1;
	bli_acquire_mpart_b2t( BLIS_SUBPART1, idx3, bs3, &X_1_1, &X_1_1_1 );
	//------------------------------------//

	bli_obj_set_struc( BLIS_TRIANGULAR, U_11_1 );
	bli_obj_set_uplo( BLIS_UPPER, U_11_1 );
	th_barrier( ProcComm );
	if (th_am_root(ProcComm)) {
	  alloced_A = TRUE;
	  bli_packm_init_pack( TRUE, BLIS_INVERT_DIAG, BLIS_PACKED_ROW_PANELS, 
			       BLIS_PACK_REV_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			       BLIS_BUFFER_FOR_A_BLOCK,
			       gemm_mr, gemm_mr, 
			       &U_11_1, &packed_A_blk_local_alloc );
	  th_broadcast_without_second_barrier(ProcComm, 0,
					      (void*)(&packed_A_blk_local_alloc),
					      (void*)(&packed_A_blk), sizeof(packed_A_blk));
	}
	else {
	  th_broadcast_without_second_barrier(ProcComm, 0, (void*)NULL,
					      (void*)(&packed_A_blk), sizeof(packed_A_blk));
	}
	bli_packm_blk_var3_par( &BLIS_ONE, &U_11_1, &packed_A_blk, ProcComm );
	bli_trsm_lu_ker_var2_par( &BLIS_ONE, &packed_A_blk, &packed_B_pan, 
				  &BLIS_ZERO, &X_1_1_1, (trsm_t*)NULL, ProcComm);

	//------------------------------------//

	//****
      }
      th_barrier( ProcComm );	 //barrier for dependency
      //// ***Parallelized with communicator ProcComm
      dimLen3 = bli_obj_length_after_trans( X_1_0 );
      idx3 = 0;
      th_shift_start_end(&idx3, &dimLen3, L2Comm, bli_blksz_for_obj( &X_1_0, gemm_mr));
      for ( ; idx3 < dimLen3; idx3 += bs3 ) {
	bs3 = bli_determine_blocksize_f( idx3, dimLen3, &X_1_0, gemm_mc );
	dim_t idx4, dimLen4, bs4;
	//****
	obj_t U_01_1;
	bli_acquire_mpart_t2b( BLIS_SUBPART1, idx3, bs3, &U_01, &U_01_1 );
	obj_t X_1_0_1;
	bli_acquire_mpart_t2b( BLIS_SUBPART1, idx3, bs3, &X_1_0, &X_1_0_1 );
	//------------------------------------//

	th_barrier( L2Comm );
	if (th_am_root(L2Comm)) {
	  alloced_A = TRUE;
	  bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_ROW_PANELS, 
			       BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			       BLIS_BUFFER_FOR_A_BLOCK,
			       gemm_mr, gemm_mr, 
			       &U_01_1, &packed_A_blk_local_alloc );
	  th_broadcast_without_second_barrier(L2Comm, 0,
					      (void*)(&packed_A_blk_local_alloc),
					      (void*)(&packed_A_blk), sizeof(packed_A_blk));
	}
	else {
	  th_broadcast_without_second_barrier(L2Comm, 0, (void*)NULL,
					      (void*)(&packed_A_blk), sizeof(packed_A_blk));
	}
	bli_packm_blk_var2_par( &BLIS_ONE, &U_01_1, &packed_A_blk, L2Comm );
	bli_gemm_ker_var2_par( &BLIS_MINUS_ONE, &packed_A_blk, &packed_B_pan, 
			       &BLIS_ONE, &X_1_0_1, (gemm_t*)NULL, L1Comm );

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

void DxT_TrsmLUT( obj_t *alpha,
		  obj_t *U,
		  obj_t *X )
{
  FUNCTIONSTART
  bli_scalm(alpha, X);

  dim_t idx1, dimLen1, bs1;
  dimLen1 = bli_obj_width_after_trans( *X );
  idx1 = 0;
  th_shift_start_end(&idx1, &dimLen1, ProcComm, bli_blksz_for_obj( X, gemm_nr));
  for ( ; idx1 < dimLen1; idx1 += bs1 ) {
    bs1 = bli_determine_blocksize_f( idx1, dimLen1, X, gemm_nc );
    dim_t idx2, dimLen2, bs2;
    //****
    obj_t X_1;
    bli_acquire_mpart_l2r( BLIS_SUBPART1, idx1, bs1, X, &X_1 );
    //------------------------------------//

    dimLen2 = bli_obj_length_after_trans( X_1 );
    for ( idx2 = 0; idx2 < dimLen2; idx2 += bs2 ) {
      bs2 = bli_determine_blocksize_f( idx2, dimLen2, &X_1, gemm_kc );
      dim_t idx3, dimLen3, bs3;
      //****
      obj_t U_11;
      bli_acquire_mpart_tl2br( BLIS_SUBPART11, idx2, bs2, U, &U_11 );
      obj_t U_12;
      bli_acquire_mpart_tl2br( BLIS_SUBPART12, idx2, bs2, U, &U_12 );
      obj_t X_1_1;
      bli_acquire_mpart_t2b( BLIS_SUBPART1, idx2, bs2, &X_1, &X_1_1 );
      obj_t X_1_2;
      bli_acquire_mpart_t2b( BLIS_SUBPART2, idx2, bs2, &X_1, &X_1_2 );
      //------------------------------------//

      th_barrier( ProcComm );
      if (th_am_root(ProcComm)) {
	alloced_B = TRUE;
	bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_COL_PANELS, 
			     BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			     BLIS_BUFFER_FOR_B_PANEL,
			     gemm_mr, gemm_nr, 
			     &X_1_1, &packed_B_pan_local_alloc );
	th_broadcast_without_second_barrier(ProcComm, 0,
					    (void*)(&packed_B_pan_local_alloc),
					    (void*)(&packed_B_pan), sizeof(packed_B_pan));
      }
      else {
	th_broadcast_without_second_barrier(ProcComm, 0, (void*)NULL,
					    (void*)(&packed_B_pan), sizeof(packed_B_pan));
      }
      bli_packm_blk_var2_par( &BLIS_ONE, &X_1_1, &packed_B_pan, ProcComm );
      dimLen3 = bli_obj_width_after_trans( U_11 );
      for ( idx3 = 0; idx3 < dimLen3; idx3 += bs3 ) {
	bs3 = bli_determine_blocksize_f( idx3, dimLen3, &U_11, gemm_mc );
	dim_t idx4, dimLen4, bs4;
	//****
	obj_t U_11_1;
	bli_acquire_mpart_l2r( BLIS_SUBPART1, idx3, bs3, &U_11, &U_11_1 );
	obj_t X_1_1_1;
	bli_acquire_mpart_t2b( BLIS_SUBPART1, idx3, bs3, &X_1_1, &X_1_1_1 );
	//------------------------------------//

	obj_t U_11_1T;
	bli_obj_alias_with_trans( BLIS_TRANSPOSE, U_11_1, U_11_1T);
	bli_obj_set_struc( BLIS_TRIANGULAR, U_11_1T );
	bli_obj_set_uplo( BLIS_UPPER, U_11_1T );
	th_barrier( ProcComm );
	if (th_am_root(ProcComm)) {
	  alloced_A = TRUE;
	  bli_packm_init_pack( TRUE, BLIS_INVERT_DIAG, BLIS_PACKED_ROW_PANELS, 
			       BLIS_PACK_REV_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			       BLIS_BUFFER_FOR_A_BLOCK,
			       gemm_mr, gemm_mr, 
			       &U_11_1T, &packed_A_blk_local_alloc );
	  th_broadcast_without_second_barrier(ProcComm, 0,
					      (void*)(&packed_A_blk_local_alloc),
					      (void*)(&packed_A_blk), sizeof(packed_A_blk));
	}
	else {
	  th_broadcast_without_second_barrier(ProcComm, 0, (void*)NULL,
					      (void*)(&packed_A_blk), sizeof(packed_A_blk));
	}
	bli_packm_blk_var3_par( &BLIS_ONE, &U_11_1T, &packed_A_blk, ProcComm );
	bli_trsm_ll_ker_var2_par( &BLIS_ONE, &packed_A_blk, &packed_B_pan, 
				  &BLIS_ZERO, &X_1_1_1, (trsm_t*)NULL, ProcComm);

	//------------------------------------//

	//****
      }
      th_barrier( ProcComm );	 //barrier for dependency
      //// ***Parallelized with communicator ProcComm
      dimLen3 = bli_obj_length_after_trans( X_1_2 );
      idx3 = 0;
      th_shift_start_end(&idx3, &dimLen3, L2Comm, bli_blksz_for_obj( &X_1_2, gemm_mr));
      for ( ; idx3 < dimLen3; idx3 += bs3 ) {
	bs3 = bli_determine_blocksize_f( idx3, dimLen3, &X_1_2, gemm_mc );
	dim_t idx4, dimLen4, bs4;
	//****
	obj_t U_12_1;
	bli_acquire_mpart_l2r( BLIS_SUBPART1, idx3, bs3, &U_12, &U_12_1 );
	obj_t X_1_2_1;
	bli_acquire_mpart_t2b( BLIS_SUBPART1, idx3, bs3, &X_1_2, &X_1_2_1 );
	//------------------------------------//

	obj_t U_12_1T;
	bli_obj_alias_with_trans( BLIS_TRANSPOSE, U_12_1, U_12_1T);
	th_barrier( L2Comm );
	if (th_am_root(L2Comm)) {
	  alloced_A = TRUE;
	  bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_ROW_PANELS, 
			       BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			       BLIS_BUFFER_FOR_A_BLOCK,
			       gemm_mr, gemm_mr, 
			       &U_12_1T, &packed_A_blk_local_alloc );
	  th_broadcast_without_second_barrier(L2Comm, 0,
					      (void*)(&packed_A_blk_local_alloc),
					      (void*)(&packed_A_blk), sizeof(packed_A_blk));
	}
	else {
	  th_broadcast_without_second_barrier(L2Comm, 0, (void*)NULL,
					      (void*)(&packed_A_blk), sizeof(packed_A_blk));
	}
	bli_packm_blk_var2_par( &BLIS_ONE, &U_12_1T, &packed_A_blk, L2Comm );
	bli_gemm_ker_var2_par( &BLIS_MINUS_ONE, &packed_A_blk, &packed_B_pan, 
			       &BLIS_ONE, &X_1_2_1, (gemm_t*)NULL, L1Comm );

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

void DxT_TrsmRLN( obj_t *alpha,
		  obj_t *L,
		  obj_t *X )
{
  FUNCTIONSTART
  bli_scalm(alpha, X);

  dim_t idx1, dimLen1, bs1;
  dimLen1 = bli_obj_width_after_trans( *X );
  for ( idx1 = 0; idx1 < dimLen1; idx1 += bs1 ) {
    bs1 = bli_determine_blocksize_b( idx1, dimLen1, X, gemm_nc );
    dim_t idx2, dimLen2, bs2;
    //****
    obj_t L_11;
    bli_acquire_mpart_br2tl( BLIS_SUBPART11, idx1, bs1, L, &L_11 );
    obj_t L_21;
    bli_acquire_mpart_br2tl( BLIS_SUBPART21, idx1, bs1, L, &L_21 );
    obj_t X_1;
    bli_acquire_mpart_r2l( BLIS_SUBPART1, idx1, bs1, X, &X_1 );
    obj_t X_2;
    bli_acquire_mpart_r2l( BLIS_SUBPART2, idx1, bs1, X, &X_2 );
    //------------------------------------//

    dimLen2 = bli_obj_width_after_trans( X_2 );
    for ( idx2 = 0; idx2 < dimLen2; idx2 += bs2 ) {
      bs2 = bli_determine_blocksize_b( idx2, dimLen2, &X_2, gemm_kc );
      dim_t idx3, dimLen3, bs3;
      //****
      obj_t X_2_1;
      bli_acquire_mpart_r2l( BLIS_SUBPART1, idx2, bs2, &X_2, &X_2_1 );
      obj_t L_21_1;
      bli_acquire_mpart_b2t( BLIS_SUBPART1, idx2, bs2, &L_21, &L_21_1 );
      //------------------------------------//

      th_barrier( AllL2Comm );
      if (th_am_root(AllL2Comm)) {
	alloced_B = TRUE;
	bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_COL_PANELS, 
			     BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			     BLIS_BUFFER_FOR_B_PANEL,
			     gemm_kr, gemm_nr, 
			     &L_21_1, &packed_B_pan_local_alloc );
	th_broadcast_without_second_barrier(AllL2Comm, 0,
					    (void*)(&packed_B_pan_local_alloc),
					    (void*)(&packed_B_pan), sizeof(packed_B_pan));
      }
      else {
	th_broadcast_without_second_barrier(AllL2Comm, 0, (void*)NULL,
					    (void*)(&packed_B_pan), sizeof(packed_B_pan));
      }
      bli_packm_blk_var2_par( &BLIS_ONE, &L_21_1, &packed_B_pan, AllL2Comm );
      //// ***Parallelized with communicator AllL2Comm
      dimLen3 = bli_obj_length_after_trans( X_1 );
      idx3 = 0;
      th_shift_start_end(&idx3, &dimLen3, L2SubAllL2Comm, bli_blksz_for_obj( &X_1, gemm_mr));
      for ( ; idx3 < dimLen3; idx3 += bs3 ) {
	bs3 = bli_determine_blocksize_f( idx3, dimLen3, &X_1, gemm_mc );
	dim_t idx4, dimLen4, bs4;
	//****
	obj_t X_2_1_1;
	bli_acquire_mpart_t2b( BLIS_SUBPART1, idx3, bs3, &X_2_1, &X_2_1_1 );
	obj_t X_1_1;
	bli_acquire_mpart_t2b( BLIS_SUBPART1, idx3, bs3, &X_1, &X_1_1 );
	//------------------------------------//

	th_barrier( L2SubAllL2Comm );
	if (th_am_root(L2SubAllL2Comm)) {
	  alloced_A = TRUE;
	  bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_ROW_PANELS, 
			       BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			       BLIS_BUFFER_FOR_A_BLOCK,
			       gemm_mr, gemm_kr, 
			       &X_2_1_1, &packed_A_blk_local_alloc );
	  th_broadcast_without_second_barrier(L2SubAllL2Comm, 0,
					      (void*)(&packed_A_blk_local_alloc),
					      (void*)(&packed_A_blk), sizeof(packed_A_blk));
	}
	else {
	  th_broadcast_without_second_barrier(L2SubAllL2Comm, 0, (void*)NULL,
					      (void*)(&packed_A_blk), sizeof(packed_A_blk));
	}
	bli_packm_blk_var2_par( &BLIS_ONE, &X_2_1_1, &packed_A_blk, L2SubAllL2Comm );
	bli_gemm_ker_var2_par( &BLIS_MINUS_ONE, &packed_A_blk, &packed_B_pan, 
			       &BLIS_ONE, &X_1_1, (gemm_t*)NULL, L1Comm );

	//------------------------------------//

	//****
      }

      //------------------------------------//

      //****
    }
    dimLen2 = bli_obj_width_after_trans( X_1 );
    for ( idx2 = 0; idx2 < dimLen2; idx2 += bs2 ) {
      bs2 = bli_determine_blocksize_b( idx2, dimLen2, &X_1, gemm_kc );
      dim_t idx3, dimLen3, bs3;
      //****
      obj_t L_11_10;
      bli_acquire_mpart_br2tl( BLIS_SUBPART10, idx2, bs2, &L_11, &L_11_10 );
      obj_t L_11_11;
      bli_acquire_mpart_br2tl( BLIS_SUBPART11, idx2, bs2, &L_11, &L_11_11 );
      obj_t X_1_0;
      bli_acquire_mpart_r2l( BLIS_SUBPART0, idx2, bs2, &X_1, &X_1_0 );
      obj_t X_1_1;
      bli_acquire_mpart_r2l( BLIS_SUBPART1, idx2, bs2, &X_1, &X_1_1 );
      //------------------------------------//

      //****
      //------------------------------------//

      bli_obj_set_struc( BLIS_TRIANGULAR, L_11_11 );
      bli_obj_set_uplo( BLIS_LOWER, L_11_11 );
      th_barrier( AllL2Comm );
      if (th_am_root(AllL2Comm)) {
	alloced_B = TRUE;
	bli_packm_init_pack( TRUE, BLIS_INVERT_DIAG, BLIS_PACKED_COL_PANELS, 
			     BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_REV_IF_LOWER, 
			     BLIS_BUFFER_FOR_B_PANEL,
			     gemm_mr, gemm_mr, 
			     &L_11_11, &packed_B_pan_local_alloc );
	th_broadcast_without_second_barrier(AllL2Comm, 0,
					    (void*)(&packed_B_pan_local_alloc),
					    (void*)(&packed_B_pan), sizeof(packed_B_pan));
      }
      else {
	th_broadcast_without_second_barrier(AllL2Comm, 0, (void*)NULL,
					    (void*)(&packed_B_pan), sizeof(packed_B_pan));
      }
      bli_packm_blk_var3_par( &BLIS_ONE, &L_11_11, &packed_B_pan, AllL2Comm );
      dim_t idx4, dimLen4, bs4;
      dimLen4 = bli_obj_length_after_trans( X_1_1 );
      for ( idx4 = 0; idx4 < dimLen4; idx4 += bs4 ) {
	bs4 = bli_determine_blocksize_b( idx4, dimLen4, &X_1_1, gemm_mc );
	dim_t idx5, dimLen5, bs5;
	//****
	obj_t X_1_1_1;
	bli_acquire_mpart_b2t( BLIS_SUBPART1, idx4, bs4, &X_1_1, &X_1_1_1 );
	//------------------------------------//

	th_barrier( AllL2Comm );
	if (th_am_root(AllL2Comm)) {
	  alloced_A = TRUE;
	  bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_ROW_PANELS, 
			       BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			       BLIS_BUFFER_FOR_A_BLOCK,
			       gemm_nr, gemm_mr, 
			       &X_1_1_1, &packed_A_blk_local_alloc );
	  th_broadcast_without_second_barrier(AllL2Comm, 0,
					      (void*)(&packed_A_blk_local_alloc),
					      (void*)(&packed_A_blk), sizeof(packed_A_blk));
	}
	else {
	  th_broadcast_without_second_barrier(AllL2Comm, 0, (void*)NULL,
					      (void*)(&packed_A_blk), sizeof(packed_A_blk));
	}
	bli_packm_blk_var2_par( &BLIS_ONE, &X_1_1_1, &packed_A_blk, AllL2Comm );
	bli_trsm_rl_ker_var2_par( &BLIS_ONE, &packed_A_blk, &packed_B_pan, 
				  &BLIS_ZERO, &X_1_1_1, (trsm_t*)NULL, AllL2Comm);

	//------------------------------------//

	//****
      }

      //------------------------------------//

      //****
      //****
      //------------------------------------//

      th_barrier( AllL2Comm );
      if (th_am_root(AllL2Comm)) {
	alloced_B = TRUE;
	bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_COL_PANELS, 
			     BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			     BLIS_BUFFER_FOR_B_PANEL,
			     gemm_kr, gemm_nr, 
			     &L_11_10, &packed_B_pan_local_alloc );
	th_broadcast_without_second_barrier(AllL2Comm, 0,
					    (void*)(&packed_B_pan_local_alloc),
					    (void*)(&packed_B_pan), sizeof(packed_B_pan));
      }
      else {
	th_broadcast_without_second_barrier(AllL2Comm, 0, (void*)NULL,
					    (void*)(&packed_B_pan), sizeof(packed_B_pan));
      }
      bli_packm_blk_var2_par( &BLIS_ONE, &L_11_10, &packed_B_pan, AllL2Comm );
      //// ***Parallelized with communicator AllL2Comm
      dimLen4 = bli_obj_length_after_trans( X_1_0 );
      idx4 = 0;
      th_shift_start_end(&idx4, &dimLen4, L2SubAllL2Comm, bli_blksz_for_obj( &X_1_0, gemm_mr));
      for ( ; idx4 < dimLen4; idx4 += bs4 ) {
	bs4 = bli_determine_blocksize_f( idx4, dimLen4, &X_1_0, gemm_mc );
	dim_t idx5, dimLen5, bs5;
	//****
	obj_t X_1_1_1;
	bli_acquire_mpart_t2b( BLIS_SUBPART1, idx4, bs4, &X_1_1, &X_1_1_1 );
	obj_t X_1_0_1;
	bli_acquire_mpart_t2b( BLIS_SUBPART1, idx4, bs4, &X_1_0, &X_1_0_1 );
	//------------------------------------//

	th_barrier( L2SubAllL2Comm );
	if (th_am_root(L2SubAllL2Comm)) {
	  alloced_A = TRUE;
	  bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_ROW_PANELS, 
			       BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			       BLIS_BUFFER_FOR_A_BLOCK,
			       gemm_mr, gemm_kr, 
			       &X_1_1_1, &packed_A_blk_local_alloc );
	  th_broadcast_without_second_barrier(L2SubAllL2Comm, 0,
					      (void*)(&packed_A_blk_local_alloc),
					      (void*)(&packed_A_blk), sizeof(packed_A_blk));
	}
	else {
	  th_broadcast_without_second_barrier(L2SubAllL2Comm, 0, (void*)NULL,
					      (void*)(&packed_A_blk), sizeof(packed_A_blk));
	}
	bli_packm_blk_var2_par( &BLIS_ONE, &X_1_1_1, &packed_A_blk, L2SubAllL2Comm );
	bli_gemm_ker_var2_par( &BLIS_MINUS_ONE, &packed_A_blk, &packed_B_pan, 
			       &BLIS_ONE, &X_1_0_1, (gemm_t*)NULL, L1Comm );

	//------------------------------------//

	//****
      }

      //------------------------------------//

      //****

      //------------------------------------//

      //****
    }

    //------------------------------------//

    //****
  }




  FUNCTIONEND
}

void DxT_TrsmRLT( obj_t *alpha,
		  obj_t *L,
		  obj_t *X )
{
  FUNCTIONSTART
  bli_scalm(alpha, X);

  dim_t idx1, dimLen1, bs1;
  dimLen1 = bli_obj_width_after_trans( *X );
  for ( idx1 = 0; idx1 < dimLen1; idx1 += bs1 ) {
    bs1 = bli_determine_blocksize_f( idx1, dimLen1, X, gemm_nc );
    dim_t idx2, dimLen2, bs2;
    //****
    obj_t L_10;
    bli_acquire_mpart_tl2br( BLIS_SUBPART10, idx1, bs1, L, &L_10 );
    obj_t L_11;
    bli_acquire_mpart_tl2br( BLIS_SUBPART11, idx1, bs1, L, &L_11 );
    obj_t X_0;
    bli_acquire_mpart_l2r( BLIS_SUBPART0, idx1, bs1, X, &X_0 );
    obj_t X_1;
    bli_acquire_mpart_l2r( BLIS_SUBPART1, idx1, bs1, X, &X_1 );
    //------------------------------------//

    dimLen2 = bli_obj_width_after_trans( X_0 );
    for ( idx2 = 0; idx2 < dimLen2; idx2 += bs2 ) {
      bs2 = bli_determine_blocksize_b( idx2, dimLen2, &X_0, gemm_kc );
      dim_t idx3, dimLen3, bs3;
      //****
      obj_t X_0_1;
      bli_acquire_mpart_r2l( BLIS_SUBPART1, idx2, bs2, &X_0, &X_0_1 );
      obj_t L_10_1;
      bli_acquire_mpart_r2l( BLIS_SUBPART1, idx2, bs2, &L_10, &L_10_1 );
      //------------------------------------//

      obj_t L_10_1T;
      bli_obj_alias_with_trans( BLIS_TRANSPOSE, L_10_1, L_10_1T);
      th_barrier( AllL2Comm );
      if (th_am_root(AllL2Comm)) {
	alloced_B = TRUE;
	bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_COL_PANELS, 
			     BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			     BLIS_BUFFER_FOR_B_PANEL,
			     gemm_kr, gemm_nr, 
			     &L_10_1T, &packed_B_pan_local_alloc );
	th_broadcast_without_second_barrier(AllL2Comm, 0,
					    (void*)(&packed_B_pan_local_alloc),
					    (void*)(&packed_B_pan), sizeof(packed_B_pan));
      }
      else {
	th_broadcast_without_second_barrier(AllL2Comm, 0, (void*)NULL,
					    (void*)(&packed_B_pan), sizeof(packed_B_pan));
      }
      bli_packm_blk_var2_par( &BLIS_ONE, &L_10_1T, &packed_B_pan, AllL2Comm );
      //// ***Parallelized with communicator AllL2Comm
      dimLen3 = bli_obj_length_after_trans( X_1 );
      idx3 = 0;
      th_shift_start_end(&idx3, &dimLen3, L2SubAllL2Comm, bli_blksz_for_obj( &X_1, gemm_mr));
      for ( ; idx3 < dimLen3; idx3 += bs3 ) {
	bs3 = bli_determine_blocksize_f( idx3, dimLen3, &X_1, gemm_mc );
	dim_t idx4, dimLen4, bs4;
	//****
	obj_t X_0_1_1;
	bli_acquire_mpart_t2b( BLIS_SUBPART1, idx3, bs3, &X_0_1, &X_0_1_1 );
	obj_t X_1_1;
	bli_acquire_mpart_t2b( BLIS_SUBPART1, idx3, bs3, &X_1, &X_1_1 );
	//------------------------------------//

	th_barrier( L2SubAllL2Comm );
	if (th_am_root(L2SubAllL2Comm)) {
	  alloced_A = TRUE;
	  bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_ROW_PANELS, 
			       BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			       BLIS_BUFFER_FOR_A_BLOCK,
			       gemm_mr, gemm_kr, 
			       &X_0_1_1, &packed_A_blk_local_alloc );
	  th_broadcast_without_second_barrier(L2SubAllL2Comm, 0,
					      (void*)(&packed_A_blk_local_alloc),
					      (void*)(&packed_A_blk), sizeof(packed_A_blk));
	}
	else {
	  th_broadcast_without_second_barrier(L2SubAllL2Comm, 0, (void*)NULL,
					      (void*)(&packed_A_blk), sizeof(packed_A_blk));
	}
	bli_packm_blk_var2_par( &BLIS_ONE, &X_0_1_1, &packed_A_blk, L2SubAllL2Comm );
	bli_gemm_ker_var2_par( &BLIS_MINUS_ONE, &packed_A_blk, &packed_B_pan, 
			       &BLIS_ONE, &X_1_1, (gemm_t*)NULL, L1Comm );

	//------------------------------------//

	//****
      }

      //------------------------------------//

      //****
    }
    dimLen2 = bli_obj_width_after_trans( X_1 );
    for ( idx2 = 0; idx2 < dimLen2; idx2 += bs2 ) {
      bs2 = bli_determine_blocksize_f( idx2, dimLen2, &X_1, gemm_kc );
      dim_t idx3, dimLen3, bs3;
      //****
      obj_t L_11_11;
      bli_acquire_mpart_tl2br( BLIS_SUBPART11, idx2, bs2, &L_11, &L_11_11 );
      obj_t L_11_21;
      bli_acquire_mpart_tl2br( BLIS_SUBPART21, idx2, bs2, &L_11, &L_11_21 );
      obj_t X_1_1;
      bli_acquire_mpart_l2r( BLIS_SUBPART1, idx2, bs2, &X_1, &X_1_1 );
      obj_t X_1_2;
      bli_acquire_mpart_l2r( BLIS_SUBPART2, idx2, bs2, &X_1, &X_1_2 );
      //------------------------------------//

      obj_t L_11_11T;
      bli_obj_alias_with_trans( BLIS_TRANSPOSE, L_11_11, L_11_11T);
      obj_t L_11_21T;
      bli_obj_alias_with_trans( BLIS_TRANSPOSE, L_11_21, L_11_21T);
      //****
      //------------------------------------//

      bli_obj_set_struc( BLIS_TRIANGULAR, L_11_11T );
      bli_obj_set_uplo( BLIS_LOWER, L_11_11T );
      th_barrier( AllL2Comm );
      if (th_am_root(AllL2Comm)) {
	alloced_B = TRUE;
	bli_packm_init_pack( TRUE, BLIS_INVERT_DIAG, BLIS_PACKED_COL_PANELS, 
			     BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_REV_IF_LOWER, 
			     BLIS_BUFFER_FOR_B_PANEL,
			     gemm_mr, gemm_mr, 
			     &L_11_11T, &packed_B_pan_local_alloc );
	th_broadcast_without_second_barrier(AllL2Comm, 0,
					    (void*)(&packed_B_pan_local_alloc),
					    (void*)(&packed_B_pan), sizeof(packed_B_pan));
      }
      else {
	th_broadcast_without_second_barrier(AllL2Comm, 0, (void*)NULL,
					    (void*)(&packed_B_pan), sizeof(packed_B_pan));
      }
      bli_packm_blk_var3_par( &BLIS_ONE, &L_11_11T, &packed_B_pan, AllL2Comm );
      	dim_t idx4, dimLen4, bs4;
dimLen4 = bli_obj_length_after_trans( X_1_1 );
      for ( idx4 = 0; idx4 < dimLen4; idx4 += bs4 ) {
	bs4 = bli_determine_blocksize_f( idx4, dimLen4, &X_1_1, gemm_mc );
	dim_t idx5, dimLen5, bs5;
	//****
	obj_t X_1_1_1;
	bli_acquire_mpart_t2b( BLIS_SUBPART1, idx4, bs4, &X_1_1, &X_1_1_1 );

	//------------------------------------//

	th_barrier( AllL2Comm );
	if (th_am_root(AllL2Comm)) {
	  alloced_A = TRUE;
	  bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_ROW_PANELS, 
			       BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			       BLIS_BUFFER_FOR_A_BLOCK,
			       gemm_nr, gemm_mr, 
			       &X_1_1_1, &packed_A_blk_local_alloc );
	  th_broadcast_without_second_barrier(AllL2Comm, 0,
					      (void*)(&packed_A_blk_local_alloc),
					      (void*)(&packed_A_blk), sizeof(packed_A_blk));
	}
	else {
	  th_broadcast_without_second_barrier(AllL2Comm, 0, (void*)NULL,
					      (void*)(&packed_A_blk), sizeof(packed_A_blk));
	}
	bli_packm_blk_var2_par( &BLIS_ONE, &X_1_1_1, &packed_A_blk, AllL2Comm );
	bli_trsm_ru_ker_var2_par( &BLIS_ONE, &packed_A_blk, &packed_B_pan, 
				  &BLIS_ZERO, &X_1_1_1, (trsm_t*)NULL, AllL2Comm);

	//------------------------------------//

	//****
      }

      //------------------------------------//

      //****
      //****
      //------------------------------------//

      th_barrier( AllL2Comm );
      if (th_am_root(AllL2Comm)) {
	alloced_B = TRUE;
	bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_COL_PANELS, 
			     BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			     BLIS_BUFFER_FOR_B_PANEL,
			     gemm_kr, gemm_nr, 
			     &L_11_21T, &packed_B_pan_local_alloc );
	th_broadcast_without_second_barrier(AllL2Comm, 0,
					    (void*)(&packed_B_pan_local_alloc),
					    (void*)(&packed_B_pan), sizeof(packed_B_pan));
      }
      else {
	th_broadcast_without_second_barrier(AllL2Comm, 0, (void*)NULL,
					    (void*)(&packed_B_pan), sizeof(packed_B_pan));
      }
      bli_packm_blk_var2_par( &BLIS_ONE, &L_11_21T, &packed_B_pan, AllL2Comm );
      //// ***Parallelized with communicator AllL2Comm
      dimLen4 = bli_obj_length_after_trans( X_1_2 );
      idx4 = 0;
      th_shift_start_end(&idx4, &dimLen4, L2SubAllL2Comm, bli_blksz_for_obj( &X_1_2, gemm_mr));
      for ( ; idx4 < dimLen4; idx4 += bs4 ) {
	bs4 = bli_determine_blocksize_f( idx4, dimLen4, &X_1_2, gemm_mc );
	dim_t idx5, dimLen5, bs5;
	//****
	obj_t X_1_1_1;
	bli_acquire_mpart_t2b( BLIS_SUBPART1, idx4, bs4, &X_1_1, &X_1_1_1 );
	obj_t X_1_2_1;
	bli_acquire_mpart_t2b( BLIS_SUBPART1, idx4, bs4, &X_1_2, &X_1_2_1 );
	//------------------------------------//

	th_barrier( L2SubAllL2Comm );
	if (th_am_root(L2SubAllL2Comm)) {
	  alloced_A = TRUE;
	  bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_ROW_PANELS, 
			       BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			       BLIS_BUFFER_FOR_A_BLOCK,
			       gemm_mr, gemm_kr, 
			       &X_1_1_1, &packed_A_blk_local_alloc );
	  th_broadcast_without_second_barrier(L2SubAllL2Comm, 0,
					      (void*)(&packed_A_blk_local_alloc),
					      (void*)(&packed_A_blk), sizeof(packed_A_blk));
	}
	else {
	  th_broadcast_without_second_barrier(L2SubAllL2Comm, 0, (void*)NULL,
					      (void*)(&packed_A_blk), sizeof(packed_A_blk));
	}
	bli_packm_blk_var2_par( &BLIS_ONE, &X_1_1_1, &packed_A_blk, L2SubAllL2Comm );
	bli_gemm_ker_var2_par( &BLIS_MINUS_ONE, &packed_A_blk, &packed_B_pan, 
			       &BLIS_ONE, &X_1_2_1, (gemm_t*)NULL, L1Comm );

	//------------------------------------//

	//****
      }

      //------------------------------------//

      //****

      //------------------------------------//

      //****
    }

    //------------------------------------//

    //****
  }


  FUNCTIONEND
}

void DxT_TrsmRUN( obj_t *alpha,
		  obj_t *U,
		  obj_t *X )
{
  FUNCTIONSTART
  bli_scalm(alpha, X);

  dim_t idx1, dimLen1, bs1;
  dimLen1 = bli_obj_width_after_trans( *X );
  for ( idx1 = 0; idx1 < dimLen1; idx1 += bs1 ) {
    bs1 = bli_determine_blocksize_f( idx1, dimLen1, X, gemm_nc );
    dim_t idx2, dimLen2, bs2;
    //****
    obj_t U_01;
    bli_acquire_mpart_tl2br( BLIS_SUBPART01, idx1, bs1, U, &U_01 );
    obj_t U_11;
    bli_acquire_mpart_tl2br( BLIS_SUBPART11, idx1, bs1, U, &U_11 );
    obj_t X_0;
    bli_acquire_mpart_l2r( BLIS_SUBPART0, idx1, bs1, X, &X_0 );
    obj_t X_1;
    bli_acquire_mpart_l2r( BLIS_SUBPART1, idx1, bs1, X, &X_1 );
    //------------------------------------//

    dimLen2 = bli_obj_width_after_trans( X_0 );
    for ( idx2 = 0; idx2 < dimLen2; idx2 += bs2 ) {
      bs2 = bli_determine_blocksize_b( idx2, dimLen2, &X_0, gemm_kc );
      dim_t idx3, dimLen3, bs3;
      //****
      obj_t X_0_1;
      bli_acquire_mpart_r2l( BLIS_SUBPART1, idx2, bs2, &X_0, &X_0_1 );
      obj_t U_01_1;
      bli_acquire_mpart_b2t( BLIS_SUBPART1, idx2, bs2, &U_01, &U_01_1 );
      //------------------------------------//

      th_barrier( AllL2Comm );
      if (th_am_root(AllL2Comm)) {
	alloced_B = TRUE;
	bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_COL_PANELS, 
			     BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			     BLIS_BUFFER_FOR_B_PANEL,
			     gemm_kr, gemm_nr, 
			     &U_01_1, &packed_B_pan_local_alloc );
	th_broadcast_without_second_barrier(AllL2Comm, 0,
					    (void*)(&packed_B_pan_local_alloc),
					    (void*)(&packed_B_pan), sizeof(packed_B_pan));
      }
      else {
	th_broadcast_without_second_barrier(AllL2Comm, 0, (void*)NULL,
					    (void*)(&packed_B_pan), sizeof(packed_B_pan));
      }
      bli_packm_blk_var2_par( &BLIS_ONE, &U_01_1, &packed_B_pan, AllL2Comm );
      //// ***Parallelized with communicator AllL2Comm
      dimLen3 = bli_obj_length_after_trans( X_1 );
      idx3 = 0;
      th_shift_start_end(&idx3, &dimLen3, L2SubAllL2Comm, bli_blksz_for_obj( &X_1, gemm_mr));
      for ( ; idx3 < dimLen3; idx3 += bs3 ) {
	bs3 = bli_determine_blocksize_f( idx3, dimLen3, &X_1, gemm_mc );
	dim_t idx4, dimLen4, bs4;
	//****
	obj_t X_0_1_1;
	bli_acquire_mpart_t2b( BLIS_SUBPART1, idx3, bs3, &X_0_1, &X_0_1_1 );
	obj_t X_1_1;
	bli_acquire_mpart_t2b( BLIS_SUBPART1, idx3, bs3, &X_1, &X_1_1 );
	//------------------------------------//

	th_barrier( L2SubAllL2Comm );
	if (th_am_root(L2SubAllL2Comm)) {
	  alloced_A = TRUE;
	  bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_ROW_PANELS, 
			       BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			       BLIS_BUFFER_FOR_A_BLOCK,
			       gemm_mr, gemm_kr, 
			       &X_0_1_1, &packed_A_blk_local_alloc );
	  th_broadcast_without_second_barrier(L2SubAllL2Comm, 0,
					      (void*)(&packed_A_blk_local_alloc),
					      (void*)(&packed_A_blk), sizeof(packed_A_blk));
	}
	else {
	  th_broadcast_without_second_barrier(L2SubAllL2Comm, 0, (void*)NULL,
					      (void*)(&packed_A_blk), sizeof(packed_A_blk));
	}
	bli_packm_blk_var2_par( &BLIS_ONE, &X_0_1_1, &packed_A_blk, L2SubAllL2Comm );
	bli_gemm_ker_var2_par( &BLIS_MINUS_ONE, &packed_A_blk, &packed_B_pan, 
			       &BLIS_ONE, &X_1_1, (gemm_t*)NULL, L1Comm );

	//------------------------------------//

	//****
      }

      //------------------------------------//

      //****
    }
    dimLen2 = bli_obj_width_after_trans( X_1 );
    for ( idx2 = 0; idx2 < dimLen2; idx2 += bs2 ) {
      bs2 = bli_determine_blocksize_f( idx2, dimLen2, &X_1, gemm_kc );
      dim_t idx3, dimLen3, bs3;
      //****
      obj_t U_11_11;
      bli_acquire_mpart_tl2br( BLIS_SUBPART11, idx2, bs2, &U_11, &U_11_11 );
      obj_t U_11_12;
      bli_acquire_mpart_tl2br( BLIS_SUBPART12, idx2, bs2, &U_11, &U_11_12 );
      obj_t X_1_1;
      bli_acquire_mpart_l2r( BLIS_SUBPART1, idx2, bs2, &X_1, &X_1_1 );
      obj_t X_1_2;
      bli_acquire_mpart_l2r( BLIS_SUBPART2, idx2, bs2, &X_1, &X_1_2 );
      //------------------------------------//

      //****
      //------------------------------------//

      bli_obj_set_struc( BLIS_TRIANGULAR, U_11_11 );
      bli_obj_set_uplo( BLIS_UPPER, U_11_11 );
      th_barrier( AllL2Comm );
      if (th_am_root(AllL2Comm)) {
	alloced_B = TRUE;
	bli_packm_init_pack( TRUE, BLIS_INVERT_DIAG, BLIS_PACKED_COL_PANELS, 
			     BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_REV_IF_LOWER, 
			     BLIS_BUFFER_FOR_B_PANEL,
			     gemm_mr, gemm_mr, 
			     &U_11_11, &packed_B_pan_local_alloc );
	th_broadcast_without_second_barrier(AllL2Comm, 0,
					    (void*)(&packed_B_pan_local_alloc),
					    (void*)(&packed_B_pan), sizeof(packed_B_pan));
      }
      else {
	th_broadcast_without_second_barrier(AllL2Comm, 0, (void*)NULL,
					    (void*)(&packed_B_pan), sizeof(packed_B_pan));
      }
      bli_packm_blk_var3_par( &BLIS_ONE, &U_11_11, &packed_B_pan, AllL2Comm );
      	dim_t idx4, dimLen4, bs4;
dimLen4 = bli_obj_length_after_trans( X_1_1 );
      for ( idx4 = 0; idx4 < dimLen4; idx4 += bs4 ) {
	bs4 = bli_determine_blocksize_f( idx4, dimLen4, &X_1_1, gemm_mc );
	dim_t idx5, dimLen5, bs5;
	//****
	obj_t X_1_1_1;
	bli_acquire_mpart_t2b( BLIS_SUBPART1, idx4, bs4, &X_1_1, &X_1_1_1 );
	
	//------------------------------------//
	th_barrier( AllL2Comm );
	if (th_am_root(AllL2Comm)) {
	  alloced_A = TRUE;
	  bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_ROW_PANELS, 
			       BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			       BLIS_BUFFER_FOR_A_BLOCK,
			       gemm_nr, gemm_mr, 
			       &X_1_1_1, &packed_A_blk_local_alloc );
	  th_broadcast_without_second_barrier(AllL2Comm, 0,
					      (void*)(&packed_A_blk_local_alloc),
					      (void*)(&packed_A_blk), sizeof(packed_A_blk));
	}
	else {
	  th_broadcast_without_second_barrier(AllL2Comm, 0, (void*)NULL,
					      (void*)(&packed_A_blk), sizeof(packed_A_blk));
	}
	bli_packm_blk_var2_par( &BLIS_ONE, &X_1_1_1, &packed_A_blk, AllL2Comm );
	bli_trsm_ru_ker_var2_par( &BLIS_ONE, &packed_A_blk, &packed_B_pan, 
				  &BLIS_ZERO, &X_1_1_1, (trsm_t*)NULL, AllL2Comm);

	//------------------------------------//

	//****
      }

      //------------------------------------//

      //****
      //****
      //------------------------------------//

      th_barrier( AllL2Comm );
      if (th_am_root(AllL2Comm)) {
	alloced_B = TRUE;
	bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_COL_PANELS, 
			     BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			     BLIS_BUFFER_FOR_B_PANEL,
			     gemm_kr, gemm_nr, 
			     &U_11_12, &packed_B_pan_local_alloc );
	th_broadcast_without_second_barrier(AllL2Comm, 0,
					    (void*)(&packed_B_pan_local_alloc),
					    (void*)(&packed_B_pan), sizeof(packed_B_pan));
      }
      else {
	th_broadcast_without_second_barrier(AllL2Comm, 0, (void*)NULL,
					    (void*)(&packed_B_pan), sizeof(packed_B_pan));
      }
      bli_packm_blk_var2_par( &BLIS_ONE, &U_11_12, &packed_B_pan, AllL2Comm );
      //// ***Parallelized with communicator AllL2Comm
      dimLen4 = bli_obj_length_after_trans( X_1_2 );
      idx4 = 0;
      th_shift_start_end(&idx4, &dimLen4, L2SubAllL2Comm, bli_blksz_for_obj( &X_1_2, gemm_mr));
      for ( ; idx4 < dimLen4; idx4 += bs4 ) {
	bs4 = bli_determine_blocksize_f( idx4, dimLen4, &X_1_2, gemm_mc );
	dim_t idx5, dimLen5, bs5;
	//****
	obj_t X_1_1_1;
	bli_acquire_mpart_t2b( BLIS_SUBPART1, idx4, bs4, &X_1_1, &X_1_1_1 );
	obj_t X_1_2_1;
	bli_acquire_mpart_t2b( BLIS_SUBPART1, idx4, bs4, &X_1_2, &X_1_2_1 );
	//------------------------------------//

	th_barrier( L2SubAllL2Comm );
	if (th_am_root(L2SubAllL2Comm)) {
	  alloced_A = TRUE;
	  bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_ROW_PANELS, 
			       BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			       BLIS_BUFFER_FOR_A_BLOCK,
			       gemm_mr, gemm_kr, 
			       &X_1_1_1, &packed_A_blk_local_alloc );
	  th_broadcast_without_second_barrier(L2SubAllL2Comm, 0,
					      (void*)(&packed_A_blk_local_alloc),
					      (void*)(&packed_A_blk), sizeof(packed_A_blk));
	}
	else {
	  th_broadcast_without_second_barrier(L2SubAllL2Comm, 0, (void*)NULL,
					      (void*)(&packed_A_blk), sizeof(packed_A_blk));
	}
	bli_packm_blk_var2_par( &BLIS_ONE, &X_1_1_1, &packed_A_blk, L2SubAllL2Comm );
	bli_gemm_ker_var2_par( &BLIS_MINUS_ONE, &packed_A_blk, &packed_B_pan, 
			       &BLIS_ONE, &X_1_2_1, (gemm_t*)NULL, L1Comm );

	//------------------------------------//

	//****
      }

      //------------------------------------//

      //****

      //------------------------------------//

      //****
    }

    //------------------------------------//

    //****
  }


  FUNCTIONEND
}


void DxT_TrsmRUT( obj_t *alpha,
		  obj_t *U,
		  obj_t *X )
{
  FUNCTIONSTART
  bli_scalm(alpha, X);

  dim_t idx1, dimLen1, bs1;
  dimLen1 = bli_obj_width_after_trans( *X );
  for ( idx1 = 0; idx1 < dimLen1; idx1 += bs1 ) {
    bs1 = bli_determine_blocksize_b( idx1, dimLen1, X, gemm_nc );
    dim_t idx2, dimLen2, bs2;
    //****
    obj_t U_11;
    bli_acquire_mpart_br2tl( BLIS_SUBPART11, idx1, bs1, U, &U_11 );
    obj_t U_12;
    bli_acquire_mpart_br2tl( BLIS_SUBPART12, idx1, bs1, U, &U_12 );
    obj_t X_1;
    bli_acquire_mpart_r2l( BLIS_SUBPART1, idx1, bs1, X, &X_1 );
    obj_t X_2;
    bli_acquire_mpart_r2l( BLIS_SUBPART2, idx1, bs1, X, &X_2 );
    //------------------------------------//

    dimLen2 = bli_obj_width_after_trans( X_2 );
    for ( idx2 = 0; idx2 < dimLen2; idx2 += bs2 ) {
      bs2 = bli_determine_blocksize_b( idx2, dimLen2, &X_2, gemm_kc );
      dim_t idx3, dimLen3, bs3;
      //****
      obj_t X_2_1;
      bli_acquire_mpart_r2l( BLIS_SUBPART1, idx2, bs2, &X_2, &X_2_1 );
      obj_t U_12_1;
      bli_acquire_mpart_r2l( BLIS_SUBPART1, idx2, bs2, &U_12, &U_12_1 );
      //------------------------------------//

      obj_t U_12_1T;
      bli_obj_alias_with_trans( BLIS_TRANSPOSE, U_12_1, U_12_1T);
      th_barrier( AllL2Comm );
      if (th_am_root(AllL2Comm)) {
	alloced_B = TRUE;
	bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_COL_PANELS, 
			     BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			     BLIS_BUFFER_FOR_B_PANEL,
			     gemm_kr, gemm_nr, 
			     &U_12_1T, &packed_B_pan_local_alloc );
	th_broadcast_without_second_barrier(AllL2Comm, 0,
					    (void*)(&packed_B_pan_local_alloc),
					    (void*)(&packed_B_pan), sizeof(packed_B_pan));
      }
      else {
	th_broadcast_without_second_barrier(AllL2Comm, 0, (void*)NULL,
					    (void*)(&packed_B_pan), sizeof(packed_B_pan));
      }
      bli_packm_blk_var2_par( &BLIS_ONE, &U_12_1T, &packed_B_pan, AllL2Comm );
      //// ***Parallelized with communicator AllL2Comm
      dimLen3 = bli_obj_length_after_trans( X_1 );
      idx3 = 0;
      th_shift_start_end(&idx3, &dimLen3, L2SubAllL2Comm, bli_blksz_for_obj( &X_1, gemm_mr));
      for ( ; idx3 < dimLen3; idx3 += bs3 ) {
	bs3 = bli_determine_blocksize_f( idx3, dimLen3, &X_1, gemm_mc );
	dim_t idx4, dimLen4, bs4;
	//****
	obj_t X_2_1_1;
	bli_acquire_mpart_t2b( BLIS_SUBPART1, idx3, bs3, &X_2_1, &X_2_1_1 );
	obj_t X_1_1;
	bli_acquire_mpart_t2b( BLIS_SUBPART1, idx3, bs3, &X_1, &X_1_1 );
	//------------------------------------//

	th_barrier( L2SubAllL2Comm );
	if (th_am_root(L2SubAllL2Comm)) {
	  alloced_A = TRUE;
	  bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_ROW_PANELS, 
			       BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			       BLIS_BUFFER_FOR_A_BLOCK,
			       gemm_mr, gemm_kr, 
			       &X_2_1_1, &packed_A_blk_local_alloc );
	  th_broadcast_without_second_barrier(L2SubAllL2Comm, 0,
					      (void*)(&packed_A_blk_local_alloc),
					      (void*)(&packed_A_blk), sizeof(packed_A_blk));
	}
	else {
	  th_broadcast_without_second_barrier(L2SubAllL2Comm, 0, (void*)NULL,
					      (void*)(&packed_A_blk), sizeof(packed_A_blk));
	}
	bli_packm_blk_var2_par( &BLIS_ONE, &X_2_1_1, &packed_A_blk, L2SubAllL2Comm );
	bli_gemm_ker_var2_par( &BLIS_MINUS_ONE, &packed_A_blk, &packed_B_pan, 
			       &BLIS_ONE, &X_1_1, (gemm_t*)NULL, L1Comm );

	//------------------------------------//

	//****
      }

      //------------------------------------//

      //****
    }
    dimLen2 = bli_obj_width_after_trans( X_1 );
    for ( idx2 = 0; idx2 < dimLen2; idx2 += bs2 ) {
      bs2 = bli_determine_blocksize_b( idx2, dimLen2, &X_1, gemm_kc );
      dim_t idx3, dimLen3, bs3;
      //****
      obj_t U_11_01;
      bli_acquire_mpart_br2tl( BLIS_SUBPART01, idx2, bs2, &U_11, &U_11_01 );
      obj_t U_11_11;
      bli_acquire_mpart_br2tl( BLIS_SUBPART11, idx2, bs2, &U_11, &U_11_11 );
      obj_t X_1_0;
      bli_acquire_mpart_r2l( BLIS_SUBPART0, idx2, bs2, &X_1, &X_1_0 );
      obj_t X_1_1;
      bli_acquire_mpart_r2l( BLIS_SUBPART1, idx2, bs2, &X_1, &X_1_1 );
      //------------------------------------//

      obj_t U_11_01T;
      bli_obj_alias_with_trans( BLIS_TRANSPOSE, U_11_01, U_11_01T);
      obj_t U_11_11T;
      bli_obj_alias_with_trans( BLIS_TRANSPOSE, U_11_11, U_11_11T);
      //****
      //------------------------------------//

      bli_obj_set_struc( BLIS_TRIANGULAR, U_11_11T );
      bli_obj_set_uplo( BLIS_UPPER, U_11_11T );
      th_barrier( AllL2Comm );
      if (th_am_root(AllL2Comm)) {
	alloced_B = TRUE;
	bli_packm_init_pack( TRUE, BLIS_INVERT_DIAG, BLIS_PACKED_COL_PANELS, 
			     BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_REV_IF_LOWER, 
			     BLIS_BUFFER_FOR_B_PANEL,
			     gemm_mr, gemm_mr, 
			     &U_11_11T, &packed_B_pan_local_alloc );
	th_broadcast_without_second_barrier(AllL2Comm, 0,
					    (void*)(&packed_B_pan_local_alloc),
					    (void*)(&packed_B_pan), sizeof(packed_B_pan));
      }
      else {
	th_broadcast_without_second_barrier(AllL2Comm, 0, (void*)NULL,
					    (void*)(&packed_B_pan), sizeof(packed_B_pan));
      }
      bli_packm_blk_var3_par( &BLIS_ONE, &U_11_11T, &packed_B_pan, AllL2Comm );
      	dim_t idx4, dimLen4, bs4;
dimLen4 = bli_obj_length_after_trans( X_1_1 );
      for ( idx4 = 0; idx4 < dimLen4; idx4 += bs4 ) {
	bs4 = bli_determine_blocksize_b( idx4, dimLen4, &X_1_1, gemm_mc );
	dim_t idx5, dimLen5, bs5;
	//****
	obj_t X_1_1_1;
	bli_acquire_mpart_b2t( BLIS_SUBPART1, idx4, bs4, &X_1_1, &X_1_1_1 );
	
	//------------------------------------//

	th_barrier( AllL2Comm );
	if (th_am_root(AllL2Comm)) {
	  alloced_A = TRUE;
	  bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_ROW_PANELS, 
			       BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			       BLIS_BUFFER_FOR_A_BLOCK,
			       gemm_nr, gemm_mr, 
			       &X_1_1_1, &packed_A_blk_local_alloc );
	  th_broadcast_without_second_barrier(AllL2Comm, 0,
					      (void*)(&packed_A_blk_local_alloc),
					      (void*)(&packed_A_blk), sizeof(packed_A_blk));
	}
	else {
	  th_broadcast_without_second_barrier(AllL2Comm, 0, (void*)NULL,
					      (void*)(&packed_A_blk), sizeof(packed_A_blk));
	}
	bli_packm_blk_var2_par( &BLIS_ONE, &X_1_1_1, &packed_A_blk, AllL2Comm );
	bli_trsm_rl_ker_var2_par( &BLIS_ONE, &packed_A_blk, &packed_B_pan, 
				  &BLIS_ZERO, &X_1_1_1, (trsm_t*)NULL, AllL2Comm);

	//------------------------------------//

	//****
      }

      //------------------------------------//

      //****
      //****
      //------------------------------------//

      th_barrier( AllL2Comm );
      if (th_am_root(AllL2Comm)) {
	alloced_B = TRUE;
	bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_COL_PANELS, 
			     BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			     BLIS_BUFFER_FOR_B_PANEL,
			     gemm_kr, gemm_nr, 
			     &U_11_01T, &packed_B_pan_local_alloc );
	th_broadcast_without_second_barrier(AllL2Comm, 0,
					    (void*)(&packed_B_pan_local_alloc),
					    (void*)(&packed_B_pan), sizeof(packed_B_pan));
      }
      else {
	th_broadcast_without_second_barrier(AllL2Comm, 0, (void*)NULL,
					    (void*)(&packed_B_pan), sizeof(packed_B_pan));
      }
      bli_packm_blk_var2_par( &BLIS_ONE, &U_11_01T, &packed_B_pan, AllL2Comm );
      //// ***Parallelized with communicator AllL2Comm
      dimLen4 = bli_obj_length_after_trans( X_1_0 );
      idx4 = 0;
      th_shift_start_end(&idx4, &dimLen4, L2SubAllL2Comm, bli_blksz_for_obj( &X_1_0, gemm_mr));
      for ( ; idx4 < dimLen4; idx4 += bs4 ) {
	bs4 = bli_determine_blocksize_f( idx4, dimLen4, &X_1_0, gemm_mc );
	dim_t idx5, dimLen5, bs5;
	//****
	obj_t X_1_1_1;
	bli_acquire_mpart_t2b( BLIS_SUBPART1, idx4, bs4, &X_1_1, &X_1_1_1 );
	obj_t X_1_0_1;
	bli_acquire_mpart_t2b( BLIS_SUBPART1, idx4, bs4, &X_1_0, &X_1_0_1 );
	//------------------------------------//

	th_barrier( L2SubAllL2Comm );
	if (th_am_root(L2SubAllL2Comm)) {
	  alloced_A = TRUE;
	  bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_ROW_PANELS, 
			       BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			       BLIS_BUFFER_FOR_A_BLOCK,
			       gemm_mr, gemm_kr, 
			       &X_1_1_1, &packed_A_blk_local_alloc );
	  th_broadcast_without_second_barrier(L2SubAllL2Comm, 0,
					      (void*)(&packed_A_blk_local_alloc),
					      (void*)(&packed_A_blk), sizeof(packed_A_blk));
	}
	else {
	  th_broadcast_without_second_barrier(L2SubAllL2Comm, 0, (void*)NULL,
					      (void*)(&packed_A_blk), sizeof(packed_A_blk));
	}
	bli_packm_blk_var2_par( &BLIS_ONE, &X_1_1_1, &packed_A_blk, L2SubAllL2Comm );
	bli_gemm_ker_var2_par( &BLIS_MINUS_ONE, &packed_A_blk, &packed_B_pan, 
			       &BLIS_ONE, &X_1_0_1, (gemm_t*)NULL, L1Comm );

	//------------------------------------//

	//****
      }

      //------------------------------------//

      //****

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
  obj_t a, c1, c2;
  obj_t c_save;
  obj_t alpha, beta, negOne;
#ifdef TESTLIB
  obj_t normVal;
#endif
  dim_t m, n;
  dim_t p;
  dim_t p_begin, p_end, p_inc;
  int   m_input, n_input;
  num_t dt_a, dt_b, dt_c;
  num_t dt_alpha, dt_beta;
  int   r, n_repeats;

  double dtimeDxT;
  double dtime_saveDxT;
  double gflopsDxT;

  double dtimeBLI;
  double dtime_saveBLI;
  double gflopsBLI;

  int left=0, lower=0, trans=0;


  dt_alpha = BLIS_DOUBLE;
  dt_beta = BLIS_DOUBLE;


  bli_obj_create( dt_alpha, 1, 1, 0, 0, &negOne );
  bli_setsc(  -1.0, 0.0, &negOne );		

  if (argc != 4) {
    printf("test L/R L/U N/T\n");
    fflush(stdout);
    return 0;
  }

  if (*(argv[1]) == 'L')
    left = 1;
  else if (*(argv[1]) != 'R') {
    printf("left/right not correct\n");
    return 0;
  }
  else
    left = 0;

  if (*(argv[2]) == 'L')
    lower = 1;
  else if (*(argv[2]) != 'U') {
    printf("lower/upper not correct\n");
    return 0;
  }
  else
    lower = 0;

  if (*(argv[3]) == 'T')
    trans = 1;
  else if (*(argv[3]) != 'N') {
    printf("trans/normal not correct\n");
    return 0;
  }
  else
    trans = 0;



  bli_init();

  n_repeats = 3;

#if 1
  p_begin = 1000;
  p_end   = 1000;
  p_inc   = 40;
#else
  p_begin = 512;
  p_end   = 6144;
  p_inc   = 512;
#endif

  m_input = -1;
  n_input = -1;

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


		bli_obj_create( dt_alpha, 1, 1, 0, 0, &alpha );
		bli_obj_create( dt_beta,  1, 1, 0, 0, &beta );
#ifdef TESTLIB
		bli_obj_create( dt_alpha, 1, 1, 0, 0, &normVal );
#endif
	
		if ( left )
			bli_obj_create( dt_a, m, m, 0, 0, &a );
		else
			bli_obj_create( dt_a, n, n, 0, 0, &a );
		bli_obj_create( dt_c, m, n, 0, 0, &c1 );
		bli_obj_create( dt_c, m, n, 0, 0, &c2 );
		bli_obj_create( dt_c, m, n, 0, 0, &c_save );


		bli_randm( &a );
		double *buff = (double*)bli_obj_buffer_at_off( a );
		dim_t rstride = bli_obj_row_stride( a );
		dim_t cstride = bli_obj_col_stride( a );
		for(unsigned int i=0;i<m;++i) {
		  buff[i*rstride+i*cstride] += m;
		}
		bli_randm( &c1 );

		bli_setsc(  (1.0/1.0), 0.0, &alpha );
		bli_setsc(  (1.0/1.0), 0.0, &beta );
#ifdef TESTLIB
		bli_setsc(  (1.0/1.0), 0.0, &normVal );
#endif

		

		bli_copym( &c1, &c_save );
	
		dtime_saveDxT = 1.0e9;
		dtime_saveBLI = 1.0e9;

		for ( r = 0; r < n_repeats; ++r )
		{
			bli_copym( &c_save, &c1 );
			bli_copym( &c_save, &c2 );

			dtimeBLI = bli_clock();

			if (trans)
			  bli_obj_set_conjtrans( BLIS_TRANSPOSE, a);

			bli_obj_set_struc( BLIS_TRIANGULAR, a );
			if (lower) {
			  bli_obj_set_uplo( BLIS_LOWER, a );
			}
			else {
			  bli_obj_set_uplo( BLIS_UPPER, a );
			}

			bli_trsm( left ? BLIS_LEFT : BLIS_RIGHT,
			          &alpha,
			          &a,
			          &c2 );
			
			if (trans)
			  bli_obj_set_conjtrans(BLIS_NO_TRANSPOSE, a);
			  // bli_obj_toggle_trans( a );
			bli_obj_set_struc( BLIS_GENERAL, a );
			bli_obj_set_uplo( BLIS_DENSE, a );
			
			dtime_saveBLI = bli_clock_min_diff( dtime_saveBLI, dtimeBLI );


			dtimeDxT = bli_clock();

			th_setup_comm(&global_comm[0], NUMTHREADS, 1);
			_Pragma( "omp parallel num_threads(NUMTHREADS)" ) 
			  {
			    if (left) {
			      if (lower) {
				if (trans) {
				  DxT_TrsmLLT(&alpha,
					      &a,
					      &c1 );
				}
				else { //normal
				  DxT_TrsmLLN(&alpha,
					      &a,
					      &c1 );
				}
			      }
			      else { //upper
				if (trans) {
				  DxT_TrsmLUT(&alpha,
					      &a,
					      &c1 );
				}
				else {//normal 
				  DxT_TrsmLUN(&alpha,
					      &a,
					      &c1 );
				}
			      }
			    }
			    else { //right
			      if (lower) {
				if (trans) {
				  DxT_TrsmRLT(&alpha,
					      &a,
					      &c1 );
				}
				else { //normal
				  DxT_TrsmRLN(&alpha,
					      &a,
					      &c1 );
				}

			      }
			      else { //upper
				if (trans) {
				  DxT_TrsmRUT(&alpha,
					      &a,
					      &c1 );
				}
				else { //normal
				  DxT_TrsmRUN(&alpha,
					      &a,
					      &c1 );
				}
			      }
			    }
			  }

			dtime_saveDxT = bli_clock_min_diff( dtime_saveDxT, dtimeDxT );
			/*
			bli_fnormm( &c1, &normVal );
			bli_printm( "c1 norm", &normVal, "%4.3e", "" );
			bli_fnormm( &c2, &normVal );
			bli_printm( "c2 norm", &normVal, "%4.3e", "" );
			*/

#ifdef TESTLIB
			bli_axpym( &negOne, &c1, &c2 );
			bli_fnormm( &c2, &normVal );
#endif

		}

		if ( left ) {
		  gflopsDxT = ( 1.0 * m * m * n ) / ( dtime_saveDxT * 1.0e9 );
		  gflopsBLI = ( 1.0 * m * m * n ) / ( dtime_saveBLI * 1.0e9 );
		}
		else {
		  gflopsDxT = ( 1.0 * m * n * n ) / ( dtime_saveDxT * 1.0e9 );
		  gflopsBLI = ( 1.0 * m * n * n ) / ( dtime_saveBLI * 1.0e9 );
		}

		printf( "data_trsm_DxT" );
		printf( "( %2ld, 1:5 ) = [ %4lu %4lu  %10.3e  %6.3f %2.4f ];\n",
			(p - p_begin + 1)/p_inc + 1, m, n, dtime_saveDxT, gflopsDxT, 
			gflopsDxT/gflopsBLI );
		
		printf( "data_trsm_BLIS" );
		printf( "( %2ld, 1:4 ) = [ %4lu %4lu  %10.3e  %6.3f ];\n\n",
		        (p - p_begin + 1)/p_inc + 1, m, n, dtime_saveBLI, gflopsBLI );

#ifdef TESTLIB
		bli_printm( "% NORM", &normVal, "%4.1f", "" );
#endif
      

		bli_obj_free( &alpha );
		bli_obj_free( &beta );

		bli_obj_free( &a );
		bli_obj_free( &c1 );
		bli_obj_free( &c2 );
		bli_obj_free( &c_save );
	}

	bli_finalize();

	return 0;
}

