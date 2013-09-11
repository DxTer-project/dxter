#include <unistd.h>
#include "blis.h"
#include "SMBLAS3.h"

thread_comm_t global_comm[1];
thread_comm_t all_l2_comms[1];
thread_comm_t proc_comms[NUMPROCS];
thread_comm_t l2_comms[NUML2];
thread_comm_t l2_comms_sub[NUML2];
thread_comm_t l1_comms[NUML1];

void DxT_Her2kNL( obj_t *alpha,
		  obj_t *A,
		  obj_t *B,
		  obj_t *beta,
		  obj_t *C )
{
  FUNCTIONSTART

    obj_t AT;
  bli_obj_alias_with_trans( BLIS_TRANSPOSE, *A, AT);
  obj_t BT;
  bli_obj_alias_with_trans( BLIS_TRANSPOSE, *B, BT);
  //// ***Parallelized with communicator GlobalComm
  dim_t idx1, dimLen1, bs1;
  dimLen1 = bli_obj_width_after_trans( *C );
  idx1 = 0;
  th_shift_start_end(&idx1, &dimLen1, ProcComm, bli_blksz_for_obj( C, gemm_nr));
  for ( ; idx1 < dimLen1; idx1 += bs1 ) {
    bs1 = bli_determine_blocksize_f( idx1, dimLen1, C, gemm_nc );
    dim_t idx2, dimLen2, bs2;
    //****
    obj_t BT_1;
    bli_acquire_mpart_l2r( BLIS_SUBPART1, idx1, bs1, &BT, &BT_1 );
    obj_t AT_1;
    bli_acquire_mpart_l2r( BLIS_SUBPART1, idx1, bs1, &AT, &AT_1 );
    obj_t C_0;
    bli_acquire_mpart_l2r( BLIS_SUBPART0, idx1, bs1, C, &C_0 );
    obj_t C_1;
    bli_acquire_mpart_l2r( BLIS_SUBPART1, idx1, bs1, C, &C_1 );
    //------------------------------------//

    obj_t C_1_B, A_B, B_B;
    dim_t offB = bli_max( 0, -bli_obj_diag_offset_after_trans( C_1 ) );
    dim_t mB = bli_obj_length_after_trans( C_1 ) - offB;
    bli_acquire_mpart_t2b( BLIS_SUBPART1,
			   offB, mB, &C_1, &C_1_B );
    bli_acquire_mpart_t2b( BLIS_SUBPART1,
			   offB, mB, A, &A_B );
    bli_acquire_mpart_t2b( BLIS_SUBPART1,
			   offB, mB, B, &B_B );
    bli_obj_set_struc( BLIS_TRIANGULAR, C_1_B );
    bli_obj_set_uplo( BLIS_LOWER, C_1_B );
    bli_scalm( &BLIS_ONE, &C_1_B );
    bli_obj_set_struc( BLIS_GENERAL, C_1_B );
    bli_obj_set_uplo( BLIS_DENSE, C_1_B );
    dimLen2 = bli_obj_width_after_trans( A_B );
    for ( idx2 = 0; idx2 < dimLen2; idx2 += bs2 ) {
      bs2 = bli_determine_blocksize_f( idx2, dimLen2, &A_B, gemm_kc );
      dim_t idx3, dimLen3, bs3;
      //****
      obj_t A_B_1;
      bli_acquire_mpart_l2r( BLIS_SUBPART1, idx2, bs2, &A_B, &A_B_1 );
      obj_t BT_1_1;
      bli_acquire_mpart_t2b( BLIS_SUBPART1, idx2, bs2, &BT_1, &BT_1_1 );
      obj_t B_B_1;
      bli_acquire_mpart_l2r( BLIS_SUBPART1, idx2, bs2, &B_B, &B_B_1 );
      obj_t AT_1_1;
      bli_acquire_mpart_t2b( BLIS_SUBPART1, idx2, bs2, &AT_1, &AT_1_1 );
      //------------------------------------//

      //****
      //------------------------------------//

      th_barrier( ProcComm );
      if (th_am_root(ProcComm)) {
	alloced_B = TRUE;
	bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_COL_PANELS, 
			     BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			     BLIS_BUFFER_FOR_B_PANEL,
			     gemm_kr, gemm_nr, 
			     &BT_1_1, &packed_B_pan );
      }
      th_broadcast_without_second_barrier(ProcComm, 0, (void*)(&packed_B_pan), sizeof(packed_B_pan));
      bli_packm_blk_var2_par( &BLIS_ONE, &BT_1_1, &packed_B_pan, ProcComm );
      //// ***Parallelized with communicator ProcComm
      			dim_t idx4, dimLen4, bs4;
			dimLen4 = bli_obj_length_after_trans( C_1_B );
      idx4 = 0;
      th_shift_start_end(&idx4, &dimLen4, L2Comm, bli_blksz_for_obj( &C_1_B, gemm_mr));
      for ( ; idx4 < dimLen4; idx4 += bs4 ) {
	bs4 = bli_determine_blocksize_b( idx4, dimLen4, &C_1_B, gemm_mc );
	dim_t idx5, dimLen5, bs5;
	//****
	obj_t A_B_1_1;
	bli_acquire_mpart_b2t( BLIS_SUBPART1, idx4, bs4, &A_B_1, &A_B_1_1 );
	obj_t C_1_B_0;
	bli_acquire_mpart_b2t( BLIS_SUBPART0, idx4, bs4, &C_1_B, &C_1_B_0 );
	obj_t C_1_B_1;
	bli_acquire_mpart_b2t( BLIS_SUBPART1, idx4, bs4, &C_1_B, &C_1_B_1 );
	//------------------------------------//

	th_barrier( L2Comm );
	if (th_am_root(L2Comm)) {
	  alloced_A = TRUE;
	  bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_ROW_PANELS, 
			       BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			       BLIS_BUFFER_FOR_A_BLOCK,
			       gemm_mr, gemm_kr, 
			       &A_B_1_1, &packed_A_blk );
	}
	th_broadcast_without_second_barrier(L2Comm, 0, (void*)(&packed_A_blk), sizeof(packed_A_blk));
	bli_obj_set_struc( BLIS_SYMMETRIC, C_1_B_1 );
	bli_obj_set_uplo( BLIS_LOWER, C_1_B_1 );
	obj_t C_1_B_1_L, packed_B_pan_L;
	dim_t offL = 0;
	dim_t nL = bli_min( bli_obj_width_after_trans( C_1_B_1 ), 
			    bli_obj_diag_offset_after_trans( C_1_B_1 ) + bs4 );
	bli_acquire_mpart_l2r( BLIS_SUBPART1,
			       offL, nL, &C_1_B_1, &C_1_B_1_L );
	bli_acquire_mpart_l2r( BLIS_SUBPART1,
			       offL, nL, &packed_B_pan, &packed_B_pan_L );
	bli_packm_blk_var2_par( &BLIS_ONE, &A_B_1_1, &packed_A_blk, L2Comm );
	bli_herk_l_ker_var2_par( &BLIS_ONE, &packed_A_blk, &packed_B_pan_L, 
				 &BLIS_ONE, &C_1_B_1_L, (herk_t*)NULL , L1Comm );

	//------------------------------------//

	//****
      }

      //------------------------------------//

      //****
      //****
      //------------------------------------//

      th_barrier( ProcComm );
      if (th_am_root(ProcComm)) {
	alloced_B = TRUE;
	bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_COL_PANELS, 
			     BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			     BLIS_BUFFER_FOR_B_PANEL,
			     gemm_kr, gemm_nr, 
			     &AT_1_1, &packed_B_pan );
      }
      th_broadcast_without_second_barrier(ProcComm, 0, (void*)(&packed_B_pan), sizeof(packed_B_pan));
      bli_packm_blk_var2_par( &BLIS_ONE, &AT_1_1, &packed_B_pan, ProcComm );
      //// ***Parallelized with communicator ProcComm
      dimLen4 = bli_obj_length_after_trans( C_1_B );
      idx4 = 0;
      th_shift_start_end(&idx4, &dimLen4, L2Comm, bli_blksz_for_obj( &C_1_B, gemm_mr));
      for ( ; idx4 < dimLen4; idx4 += bs4 ) {
	bs4 = bli_determine_blocksize_b( idx4, dimLen4, &C_1_B, gemm_mc );
	dim_t idx5, dimLen5, bs5;
	//****
	obj_t B_B_1_1;
	bli_acquire_mpart_b2t( BLIS_SUBPART1, idx4, bs4, &B_B_1, &B_B_1_1 );
	obj_t C_1_B_0;
	bli_acquire_mpart_b2t( BLIS_SUBPART0, idx4, bs4, &C_1_B, &C_1_B_0 );
	obj_t C_1_B_1;
	bli_acquire_mpart_b2t( BLIS_SUBPART1, idx4, bs4, &C_1_B, &C_1_B_1 );
	//------------------------------------//

	th_barrier( L2Comm );
	if (th_am_root(L2Comm)) {
	  alloced_A = TRUE;
	  bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_ROW_PANELS, 
			       BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			       BLIS_BUFFER_FOR_A_BLOCK,
			       gemm_mr, gemm_kr, 
			       &B_B_1_1, &packed_A_blk );
	}
	th_broadcast_without_second_barrier(L2Comm, 0, (void*)(&packed_A_blk), sizeof(packed_A_blk));
	bli_obj_set_struc( BLIS_SYMMETRIC, C_1_B_1 );
	bli_obj_set_uplo( BLIS_LOWER, C_1_B_1 );
	obj_t C_1_B_1_L, packed_B_pan_L;
	dim_t offL = 0;
	dim_t nL = bli_min( bli_obj_width_after_trans( C_1_B_1 ), 
			    bli_obj_diag_offset_after_trans( C_1_B_1 ) + bs4 );
	bli_acquire_mpart_l2r( BLIS_SUBPART1,
			       offL, nL, &C_1_B_1, &C_1_B_1_L );
	bli_acquire_mpart_l2r( BLIS_SUBPART1,
			       offL, nL, &packed_B_pan, &packed_B_pan_L );
	bli_packm_blk_var2_par( &BLIS_ONE, &B_B_1_1, &packed_A_blk, L2Comm );
	bli_herk_l_ker_var2_par( &BLIS_ONE, &packed_A_blk, &packed_B_pan_L, 
				 &BLIS_ONE, &C_1_B_1_L, (herk_t*)NULL , L1Comm );

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

void DxT_Her2kNU( obj_t *alpha,
		  obj_t *A,
		  obj_t *B,
		  obj_t *beta,
		  obj_t *C )
{
  FUNCTIONSTART

    obj_t AT;
  bli_obj_alias_with_trans( BLIS_TRANSPOSE, *A, AT);
  obj_t BT;
  bli_obj_alias_with_trans( BLIS_TRANSPOSE, *B, BT);
  //// ***Parallelized with communicator GlobalComm
  dim_t idx1, dimLen1, bs1;
  dimLen1 = bli_obj_width_after_trans( *C );
  idx1 = 0;
  th_shift_start_end(&idx1, &dimLen1, ProcComm, bli_blksz_for_obj( C, gemm_nr));
  for ( ; idx1 < dimLen1; idx1 += bs1 ) {
    bs1 = bli_determine_blocksize_b( idx1, dimLen1, C, gemm_nc );
    dim_t idx2, dimLen2, bs2;
    //****
    obj_t BT_1;
    bli_acquire_mpart_r2l( BLIS_SUBPART1, idx1, bs1, &BT, &BT_1 );
    obj_t AT_1;
    bli_acquire_mpart_r2l( BLIS_SUBPART1, idx1, bs1, &AT, &AT_1 );
    obj_t C_1;
    bli_acquire_mpart_r2l( BLIS_SUBPART1, idx1, bs1, C, &C_1 );
    obj_t C_2;
    bli_acquire_mpart_r2l( BLIS_SUBPART2, idx1, bs1, C, &C_2 );
    //------------------------------------//

    obj_t C_1_T, A_T, B_T;
    dim_t offT = 0;
    dim_t mT = bli_min( bli_obj_length_after_trans( C_1 ), 
			-bli_obj_diag_offset_after_trans( C_1 ) + bs1 );
    bli_acquire_mpart_t2b( BLIS_SUBPART1,
			   offT, mT, &C_1, &C_1_T );
    bli_acquire_mpart_t2b( BLIS_SUBPART1,
			   offT, mT, A, &A_T );
    bli_acquire_mpart_t2b( BLIS_SUBPART1,
			   offT, mT, B, &B_T );
    bli_obj_set_struc( BLIS_TRIANGULAR, C_1_T );
    bli_obj_set_uplo( BLIS_LOWER, C_1_T );
    bli_scalm( &BLIS_ONE, &C_1_T );
    bli_obj_set_struc( BLIS_GENERAL, C_1_T );
    bli_obj_set_uplo( BLIS_DENSE, C_1_T );
    dimLen2 = bli_obj_width_after_trans( A_T );
    for ( idx2 = 0; idx2 < dimLen2; idx2 += bs2 ) {
      bs2 = bli_determine_blocksize_f( idx2, dimLen2, &A_T, gemm_kc );
      dim_t idx3, dimLen3, bs3;
      //****
      obj_t A_T_1;
      bli_acquire_mpart_l2r( BLIS_SUBPART1, idx2, bs2, &A_T, &A_T_1 );
      obj_t BT_1_1;
      bli_acquire_mpart_t2b( BLIS_SUBPART1, idx2, bs2, &BT_1, &BT_1_1 );
      obj_t B_T_1;
      bli_acquire_mpart_l2r( BLIS_SUBPART1, idx2, bs2, &B_T, &B_T_1 );
      obj_t AT_1_1;
      bli_acquire_mpart_t2b( BLIS_SUBPART1, idx2, bs2, &AT_1, &AT_1_1 );
      //------------------------------------//

      //****
      //------------------------------------//

      th_barrier( ProcComm );
      if (th_am_root(ProcComm)) {
	alloced_B = TRUE;
	bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_COL_PANELS, 
			     BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			     BLIS_BUFFER_FOR_B_PANEL,
			     gemm_kr, gemm_nr, 
			     &BT_1_1, &packed_B_pan );
      }
      th_broadcast_without_second_barrier(ProcComm, 0, (void*)(&packed_B_pan), sizeof(packed_B_pan));
      bli_packm_blk_var2_par( &BLIS_ONE, &BT_1_1, &packed_B_pan, ProcComm );
      //// ***Parallelized with communicator ProcComm
      			dim_t idx4, dimLen4, bs4;
dimLen4 = bli_obj_length_after_trans( C_1_T );
      idx4 = 0;
      th_shift_start_end(&idx4, &dimLen4, L2Comm, bli_blksz_for_obj( &C_1_T, gemm_mr));
      for ( ; idx4 < dimLen4; idx4 += bs4 ) {
	bs4 = bli_determine_blocksize_f( idx4, dimLen4, &C_1_T, gemm_mc );
	dim_t idx5, dimLen5, bs5;
	//****
	obj_t A_T_1_1;
	bli_acquire_mpart_t2b( BLIS_SUBPART1, idx4, bs4, &A_T_1, &A_T_1_1 );
	obj_t C_1_T_1;
	bli_acquire_mpart_t2b( BLIS_SUBPART1, idx4, bs4, &C_1_T, &C_1_T_1 );
	obj_t C_1_T_2;
	bli_acquire_mpart_t2b( BLIS_SUBPART2, idx4, bs4, &C_1_T, &C_1_T_2 );
	//------------------------------------//

	th_barrier( L2Comm );
	if (th_am_root(L2Comm)) {
	  alloced_A = TRUE;
	  bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_ROW_PANELS, 
			       BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			       BLIS_BUFFER_FOR_A_BLOCK,
			       gemm_mr, gemm_kr, 
			       &A_T_1_1, &packed_A_blk );
	}
	th_broadcast_without_second_barrier(L2Comm, 0, (void*)(&packed_A_blk), sizeof(packed_A_blk));
	bli_obj_set_struc( BLIS_SYMMETRIC, C_1_T_1 );
	bli_obj_set_uplo( BLIS_UPPER, C_1_T_1 );
	obj_t C_1_T_1_R, packed_B_pan_R;
	dim_t offR = bli_max( 0, bli_obj_diag_offset_after_trans( C_1_T_1 ) );
	dim_t nR = bli_obj_width_after_trans( C_1_T_1 ) - offR;
	bli_acquire_mpart_l2r( BLIS_SUBPART1,
			       offR, nR, &C_1_T_1, &C_1_T_1_R );
	bli_acquire_mpart_l2r( BLIS_SUBPART1,
			       offR, nR, &packed_B_pan, &packed_B_pan_R );
	bli_packm_blk_var2_par( &BLIS_ONE, &A_T_1_1, &packed_A_blk, L2Comm );
	bli_herk_u_ker_var2_par( &BLIS_ONE, &packed_A_blk, &packed_B_pan_R, 
				 &BLIS_ONE, &C_1_T_1_R, (herk_t*)NULL , L1Comm );

	//------------------------------------//

	//****
      }

      //------------------------------------//

      //****
      //****
      //------------------------------------//

      th_barrier( ProcComm );
      if (th_am_root(ProcComm)) {
	alloced_B = TRUE;
	bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_COL_PANELS, 
			     BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			     BLIS_BUFFER_FOR_B_PANEL,
			     gemm_kr, gemm_nr, 
			     &AT_1_1, &packed_B_pan );
      }
      th_broadcast_without_second_barrier(ProcComm, 0, (void*)(&packed_B_pan), sizeof(packed_B_pan));
      bli_packm_blk_var2_par( &BLIS_ONE, &AT_1_1, &packed_B_pan, ProcComm );
      //// ***Parallelized with communicator ProcComm
      dimLen4 = bli_obj_length_after_trans( C_1_T );
      idx4 = 0;
      th_shift_start_end(&idx4, &dimLen4, L2Comm, bli_blksz_for_obj( &C_1_T, gemm_mr));
      for ( ; idx4 < dimLen4; idx4 += bs4 ) {
	bs4 = bli_determine_blocksize_f( idx4, dimLen4, &C_1_T, gemm_mc );
	dim_t idx5, dimLen5, bs5;
	//****
	obj_t B_T_1_1;
	bli_acquire_mpart_t2b( BLIS_SUBPART1, idx4, bs4, &B_T_1, &B_T_1_1 );
	obj_t C_1_T_1;
	bli_acquire_mpart_t2b( BLIS_SUBPART1, idx4, bs4, &C_1_T, &C_1_T_1 );
	obj_t C_1_T_2;
	bli_acquire_mpart_t2b( BLIS_SUBPART2, idx4, bs4, &C_1_T, &C_1_T_2 );
	//------------------------------------//

	th_barrier( L2Comm );
	if (th_am_root(L2Comm)) {
	  alloced_A = TRUE;
	  bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_ROW_PANELS, 
			       BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			       BLIS_BUFFER_FOR_A_BLOCK,
			       gemm_mr, gemm_kr, 
			       &B_T_1_1, &packed_A_blk );
	}
	th_broadcast_without_second_barrier(L2Comm, 0, (void*)(&packed_A_blk), sizeof(packed_A_blk));
	bli_obj_set_struc( BLIS_SYMMETRIC, C_1_T_1 );
	bli_obj_set_uplo( BLIS_UPPER, C_1_T_1 );
	obj_t C_1_T_1_R, packed_B_pan_R;
	dim_t offR = bli_max( 0, bli_obj_diag_offset_after_trans( C_1_T_1 ) );
	dim_t nR = bli_obj_width_after_trans( C_1_T_1 ) - offR;
	bli_acquire_mpart_l2r( BLIS_SUBPART1,
			       offR, nR, &C_1_T_1, &C_1_T_1_R );
	bli_acquire_mpart_l2r( BLIS_SUBPART1,
			       offR, nR, &packed_B_pan, &packed_B_pan_R );
	bli_packm_blk_var2_par( &BLIS_ONE, &B_T_1_1, &packed_A_blk, L2Comm );
	bli_herk_u_ker_var2_par( &BLIS_ONE, &packed_A_blk, &packed_B_pan_R, 
				 &BLIS_ONE, &C_1_T_1_R, (herk_t*)NULL , L1Comm );

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

void DxT_Her2kTL( obj_t *alpha,
		  obj_t *A,
		  obj_t *B,
		  obj_t *beta,
		  obj_t *C )
{
  FUNCTIONSTART

    obj_t AT;
  bli_obj_alias_with_trans( BLIS_TRANSPOSE, *A, AT);
  obj_t BT;
  bli_obj_alias_with_trans( BLIS_TRANSPOSE, *B, BT);
  //// ***Parallelized with communicator GlobalComm
  dim_t idx1, dimLen1, bs1;
  dimLen1 = bli_obj_width_after_trans( *C );
  idx1 = 0;
  th_shift_start_end(&idx1, &dimLen1, ProcComm, bli_blksz_for_obj( C, gemm_nr));
  for ( ; idx1 < dimLen1; idx1 += bs1 ) {
    bs1 = bli_determine_blocksize_f( idx1, dimLen1, C, gemm_nc );
    dim_t idx2, dimLen2, bs2;
    //****
    obj_t A_1;
    bli_acquire_mpart_l2r( BLIS_SUBPART1, idx1, bs1, A, &A_1 );
    obj_t B_1;
    bli_acquire_mpart_l2r( BLIS_SUBPART1, idx1, bs1, B, &B_1 );
    obj_t C_0;
    bli_acquire_mpart_l2r( BLIS_SUBPART0, idx1, bs1, C, &C_0 );
    obj_t C_1;
    bli_acquire_mpart_l2r( BLIS_SUBPART1, idx1, bs1, C, &C_1 );
    //------------------------------------//

    obj_t C_1_B, BT_B, AT_B;
    dim_t offB = bli_max( 0, -bli_obj_diag_offset_after_trans( C_1 ) );
    dim_t mB = bli_obj_length_after_trans( C_1 ) - offB;
    bli_acquire_mpart_t2b( BLIS_SUBPART1,
			   offB, mB, &C_1, &C_1_B );
    bli_acquire_mpart_t2b( BLIS_SUBPART1,
			   offB, mB, &BT, &BT_B );
    bli_acquire_mpart_t2b( BLIS_SUBPART1,
			   offB, mB, &AT, &AT_B );
    bli_obj_set_struc( BLIS_TRIANGULAR, C_1_B );
    bli_obj_set_uplo( BLIS_LOWER, C_1_B );
    bli_scalm( &BLIS_ONE, &C_1_B );
    bli_obj_set_struc( BLIS_GENERAL, C_1_B );
    bli_obj_set_uplo( BLIS_DENSE, C_1_B );
    dimLen2 = bli_obj_width_after_trans( BT_B );
    for ( idx2 = 0; idx2 < dimLen2; idx2 += bs2 ) {
      bs2 = bli_determine_blocksize_f( idx2, dimLen2, &BT_B, gemm_kc );
      dim_t idx3, dimLen3, bs3;
      //****
      obj_t BT_B_1;
      bli_acquire_mpart_l2r( BLIS_SUBPART1, idx2, bs2, &BT_B, &BT_B_1 );
      obj_t A_1_1;
      bli_acquire_mpart_t2b( BLIS_SUBPART1, idx2, bs2, &A_1, &A_1_1 );
      obj_t AT_B_1;
      bli_acquire_mpart_l2r( BLIS_SUBPART1, idx2, bs2, &AT_B, &AT_B_1 );
      obj_t B_1_1;
      bli_acquire_mpart_t2b( BLIS_SUBPART1, idx2, bs2, &B_1, &B_1_1 );
      //------------------------------------//

      //****
      //------------------------------------//

      th_barrier( ProcComm );
      if (th_am_root(ProcComm)) {
	alloced_B = TRUE;
	bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_COL_PANELS, 
			     BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			     BLIS_BUFFER_FOR_B_PANEL,
			     gemm_kr, gemm_nr, 
			     &A_1_1, &packed_B_pan );
      }
      th_broadcast_without_second_barrier(ProcComm, 0, (void*)(&packed_B_pan), sizeof(packed_B_pan));
      bli_packm_blk_var2_par( &BLIS_ONE, &A_1_1, &packed_B_pan, ProcComm );
      //// ***Parallelized with communicator ProcComm
      dim_t idx4, dimLen4, bs4;
      dimLen4 = bli_obj_length_after_trans( C_1_B );
      idx4 = 0;
      th_shift_start_end(&idx4, &dimLen4, L2Comm, bli_blksz_for_obj( &C_1_B, gemm_mr));
      for ( ; idx4 < dimLen4; idx4 += bs4 ) {
	bs4 = bli_determine_blocksize_b( idx4, dimLen4, &C_1_B, gemm_mc );
	dim_t idx5, dimLen5, bs5;
	//****
	obj_t BT_B_1_1;
	bli_acquire_mpart_b2t( BLIS_SUBPART1, idx4, bs4, &BT_B_1, &BT_B_1_1 );
	obj_t C_1_B_0;
	bli_acquire_mpart_b2t( BLIS_SUBPART0, idx4, bs4, &C_1_B, &C_1_B_0 );
	obj_t C_1_B_1;
	bli_acquire_mpart_b2t( BLIS_SUBPART1, idx4, bs4, &C_1_B, &C_1_B_1 );
	//------------------------------------//

	th_barrier( L2Comm );
	if (th_am_root(L2Comm)) {
	  alloced_A = TRUE;
	  bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_ROW_PANELS, 
			       BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			       BLIS_BUFFER_FOR_A_BLOCK,
			       gemm_mr, gemm_kr, 
			       &BT_B_1_1, &packed_A_blk );
	}
	th_broadcast_without_second_barrier(L2Comm, 0, (void*)(&packed_A_blk), sizeof(packed_A_blk));
	bli_obj_set_struc( BLIS_SYMMETRIC, C_1_B_1 );
	bli_obj_set_uplo( BLIS_LOWER, C_1_B_1 );
	obj_t C_1_B_1_L, packed_B_pan_L;
	dim_t offL = 0;
	dim_t nL = bli_min( bli_obj_width_after_trans( C_1_B_1 ), 
			    bli_obj_diag_offset_after_trans( C_1_B_1 ) + bs4 );
	bli_acquire_mpart_l2r( BLIS_SUBPART1,
			       offL, nL, &C_1_B_1, &C_1_B_1_L );
	bli_acquire_mpart_l2r( BLIS_SUBPART1,
			       offL, nL, &packed_B_pan, &packed_B_pan_L );
	bli_packm_blk_var2_par( &BLIS_ONE, &BT_B_1_1, &packed_A_blk, L2Comm );
	bli_herk_l_ker_var2_par( &BLIS_ONE, &packed_A_blk, &packed_B_pan_L, 
				 &BLIS_ONE, &C_1_B_1_L, (herk_t*)NULL , L1Comm );

	//------------------------------------//

	//****
      }

      //------------------------------------//

      //****
      //****
      //------------------------------------//

      th_barrier( ProcComm );
      if (th_am_root(ProcComm)) {
	alloced_B = TRUE;
	bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_COL_PANELS, 
			     BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			     BLIS_BUFFER_FOR_B_PANEL,
			     gemm_kr, gemm_nr, 
			     &B_1_1, &packed_B_pan );
      }
      th_broadcast_without_second_barrier(ProcComm, 0, (void*)(&packed_B_pan), sizeof(packed_B_pan));
      bli_packm_blk_var2_par( &BLIS_ONE, &B_1_1, &packed_B_pan, ProcComm );
      //// ***Parallelized with communicator ProcComm
      dimLen4 = bli_obj_length_after_trans( C_1_B );
      idx4 = 0;
      th_shift_start_end(&idx4, &dimLen4, L2Comm, bli_blksz_for_obj( &C_1_B, gemm_mr));
      for ( ; idx4 < dimLen4; idx4 += bs4 ) {
	bs4 = bli_determine_blocksize_b( idx4, dimLen4, &C_1_B, gemm_mc );
	dim_t idx5, dimLen5, bs5;
	//****
	obj_t AT_B_1_1;
	bli_acquire_mpart_b2t( BLIS_SUBPART1, idx4, bs4, &AT_B_1, &AT_B_1_1 );
	obj_t C_1_B_0;
	bli_acquire_mpart_b2t( BLIS_SUBPART0, idx4, bs4, &C_1_B, &C_1_B_0 );
	obj_t C_1_B_1;
	bli_acquire_mpart_b2t( BLIS_SUBPART1, idx4, bs4, &C_1_B, &C_1_B_1 );
	//------------------------------------//

	th_barrier( L2Comm );
	if (th_am_root(L2Comm)) {
	  alloced_A = TRUE;
	  bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_ROW_PANELS, 
			       BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			       BLIS_BUFFER_FOR_A_BLOCK,
			       gemm_mr, gemm_kr, 
			       &AT_B_1_1, &packed_A_blk );
	}
	th_broadcast_without_second_barrier(L2Comm, 0, (void*)(&packed_A_blk), sizeof(packed_A_blk));
	bli_obj_set_struc( BLIS_SYMMETRIC, C_1_B_1 );
	bli_obj_set_uplo( BLIS_LOWER, C_1_B_1 );
	obj_t C_1_B_1_L, packed_B_pan_L;
	dim_t offL = 0;
	dim_t nL = bli_min( bli_obj_width_after_trans( C_1_B_1 ), 
			    bli_obj_diag_offset_after_trans( C_1_B_1 ) + bs4 );
	bli_acquire_mpart_l2r( BLIS_SUBPART1,
			       offL, nL, &C_1_B_1, &C_1_B_1_L );
	bli_acquire_mpart_l2r( BLIS_SUBPART1,
			       offL, nL, &packed_B_pan, &packed_B_pan_L );
	bli_packm_blk_var2_par( &BLIS_ONE, &AT_B_1_1, &packed_A_blk, L2Comm );
	bli_herk_l_ker_var2_par( &BLIS_ONE, &packed_A_blk, &packed_B_pan_L, 
				 &BLIS_ONE, &C_1_B_1_L, (herk_t*)NULL , L1Comm );

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

void DxT_Her2kTU( obj_t *alpha,
		  obj_t *A,
		  obj_t *B,
		  obj_t *beta,
		  obj_t *C )
{
  FUNCTIONSTART

    obj_t AT;
  bli_obj_alias_with_trans( BLIS_TRANSPOSE, *A, AT);
  obj_t BT;
  bli_obj_alias_with_trans( BLIS_TRANSPOSE, *B, BT);
  //// ***Parallelized with communicator GlobalComm
  dim_t idx1, dimLen1, bs1;
  dimLen1 = bli_obj_width_after_trans( *C );
  idx1 = 0;
  th_shift_start_end(&idx1, &dimLen1, ProcComm, bli_blksz_for_obj( C, gemm_nr));
  for ( ; idx1 < dimLen1; idx1 += bs1 ) {
    bs1 = bli_determine_blocksize_b( idx1, dimLen1, C, gemm_nc );
    dim_t idx2, dimLen2, bs2;
    //****
    obj_t A_1;
    bli_acquire_mpart_r2l( BLIS_SUBPART1, idx1, bs1, A, &A_1 );
    obj_t B_1;
    bli_acquire_mpart_r2l( BLIS_SUBPART1, idx1, bs1, B, &B_1 );
    obj_t C_1;
    bli_acquire_mpart_r2l( BLIS_SUBPART1, idx1, bs1, C, &C_1 );
    obj_t C_2;
    bli_acquire_mpart_r2l( BLIS_SUBPART2, idx1, bs1, C, &C_2 );
    //------------------------------------//

    obj_t C_1_T, BT_T, AT_T;
    dim_t offT = 0;
    dim_t mT = bli_min( bli_obj_length_after_trans( C_1 ), 
			-bli_obj_diag_offset_after_trans( C_1 ) + bs1 );
    bli_acquire_mpart_t2b( BLIS_SUBPART1,
			   offT, mT, &C_1, &C_1_T );
    bli_acquire_mpart_t2b( BLIS_SUBPART1,
			   offT, mT, &BT, &BT_T );
    bli_acquire_mpart_t2b( BLIS_SUBPART1,
			   offT, mT, &AT, &AT_T );
    bli_obj_set_struc( BLIS_TRIANGULAR, C_1_T );
    bli_obj_set_uplo( BLIS_LOWER, C_1_T );
    bli_scalm( &BLIS_ONE, &C_1_T );
    bli_obj_set_struc( BLIS_GENERAL, C_1_T );
    bli_obj_set_uplo( BLIS_DENSE, C_1_T );
    dimLen2 = bli_obj_width_after_trans( BT_T );
    for ( idx2 = 0; idx2 < dimLen2; idx2 += bs2 ) {
      bs2 = bli_determine_blocksize_f( idx2, dimLen2, &BT_T, gemm_kc );
      dim_t idx3, dimLen3, bs3;
      //****
      obj_t BT_T_1;
      bli_acquire_mpart_l2r( BLIS_SUBPART1, idx2, bs2, &BT_T, &BT_T_1 );
      obj_t A_1_1;
      bli_acquire_mpart_t2b( BLIS_SUBPART1, idx2, bs2, &A_1, &A_1_1 );
      obj_t AT_T_1;
      bli_acquire_mpart_l2r( BLIS_SUBPART1, idx2, bs2, &AT_T, &AT_T_1 );
      obj_t B_1_1;
      bli_acquire_mpart_t2b( BLIS_SUBPART1, idx2, bs2, &B_1, &B_1_1 );
      //------------------------------------//

      //****
      //------------------------------------//

      th_barrier( ProcComm );
      if (th_am_root(ProcComm)) {
	alloced_B = TRUE;
	bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_COL_PANELS, 
			     BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			     BLIS_BUFFER_FOR_B_PANEL,
			     gemm_kr, gemm_nr, 
			     &A_1_1, &packed_B_pan );
      }
      th_broadcast_without_second_barrier(ProcComm, 0, (void*)(&packed_B_pan), sizeof(packed_B_pan));
      bli_packm_blk_var2_par( &BLIS_ONE, &A_1_1, &packed_B_pan, ProcComm );
      //// ***Parallelized with communicator ProcComm
      			dim_t idx4, dimLen4, bs4;
dimLen4 = bli_obj_length_after_trans( C_1_T );
      idx4 = 0;
      th_shift_start_end(&idx4, &dimLen4, L2Comm, bli_blksz_for_obj( &C_1_T, gemm_mr));
      for ( ; idx4 < dimLen4; idx4 += bs4 ) {
	bs4 = bli_determine_blocksize_f( idx4, dimLen4, &C_1_T, gemm_mc );
	dim_t idx5, dimLen5, bs5;
	//****
	obj_t BT_T_1_1;
	bli_acquire_mpart_t2b( BLIS_SUBPART1, idx4, bs4, &BT_T_1, &BT_T_1_1 );
	obj_t C_1_T_1;
	bli_acquire_mpart_t2b( BLIS_SUBPART1, idx4, bs4, &C_1_T, &C_1_T_1 );
	obj_t C_1_T_2;
	bli_acquire_mpart_t2b( BLIS_SUBPART2, idx4, bs4, &C_1_T, &C_1_T_2 );
	//------------------------------------//

	th_barrier( L2Comm );
	if (th_am_root(L2Comm)) {
	  alloced_A = TRUE;
	  bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_ROW_PANELS, 
			       BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			       BLIS_BUFFER_FOR_A_BLOCK,
			       gemm_mr, gemm_kr, 
			       &BT_T_1_1, &packed_A_blk );
	}
	th_broadcast_without_second_barrier(L2Comm, 0, (void*)(&packed_A_blk), sizeof(packed_A_blk));
	bli_obj_set_struc( BLIS_SYMMETRIC, C_1_T_1 );
	bli_obj_set_uplo( BLIS_UPPER, C_1_T_1 );
	obj_t C_1_T_1_R, packed_B_pan_R;
	dim_t offR = bli_max( 0, bli_obj_diag_offset_after_trans( C_1_T_1 ) );
	dim_t nR = bli_obj_width_after_trans( C_1_T_1 ) - offR;
	bli_acquire_mpart_l2r( BLIS_SUBPART1,
			       offR, nR, &C_1_T_1, &C_1_T_1_R );
	bli_acquire_mpart_l2r( BLIS_SUBPART1,
			       offR, nR, &packed_B_pan, &packed_B_pan_R );
	bli_packm_blk_var2_par( &BLIS_ONE, &BT_T_1_1, &packed_A_blk, L2Comm );
	bli_herk_u_ker_var2_par( &BLIS_ONE, &packed_A_blk, &packed_B_pan_R, 
				 &BLIS_ONE, &C_1_T_1_R, (herk_t*)NULL , L1Comm );

	//------------------------------------//

	//****
      }

      //------------------------------------//

      //****
      //****
      //------------------------------------//

      th_barrier( ProcComm );
      if (th_am_root(ProcComm)) {
	alloced_B = TRUE;
	bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_COL_PANELS, 
			     BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			     BLIS_BUFFER_FOR_B_PANEL,
			     gemm_kr, gemm_nr, 
			     &B_1_1, &packed_B_pan );
      }
      th_broadcast_without_second_barrier(ProcComm, 0, (void*)(&packed_B_pan), sizeof(packed_B_pan));
      bli_packm_blk_var2_par( &BLIS_ONE, &B_1_1, &packed_B_pan, ProcComm );
      //// ***Parallelized with communicator ProcComm
      dimLen4 = bli_obj_length_after_trans( C_1_T );
      idx4 = 0;
      th_shift_start_end(&idx4, &dimLen4, L2Comm, bli_blksz_for_obj( &C_1_T, gemm_mr));
      for ( ; idx4 < dimLen4; idx4 += bs4 ) {
	bs4 = bli_determine_blocksize_f( idx4, dimLen4, &C_1_T, gemm_mc );
	dim_t idx5, dimLen5, bs5;
	//****
	obj_t AT_T_1_1;
	bli_acquire_mpart_t2b( BLIS_SUBPART1, idx4, bs4, &AT_T_1, &AT_T_1_1 );
	obj_t C_1_T_1;
	bli_acquire_mpart_t2b( BLIS_SUBPART1, idx4, bs4, &C_1_T, &C_1_T_1 );
	obj_t C_1_T_2;
	bli_acquire_mpart_t2b( BLIS_SUBPART2, idx4, bs4, &C_1_T, &C_1_T_2 );
	//------------------------------------//

	th_barrier( L2Comm );
	if (th_am_root(L2Comm)) {
	  alloced_A = TRUE;
	  bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_ROW_PANELS, 
			       BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			       BLIS_BUFFER_FOR_A_BLOCK,
			       gemm_mr, gemm_kr, 
			       &AT_T_1_1, &packed_A_blk );
	}
	th_broadcast_without_second_barrier(L2Comm, 0, (void*)(&packed_A_blk), sizeof(packed_A_blk));
	bli_obj_set_struc( BLIS_SYMMETRIC, C_1_T_1 );
	bli_obj_set_uplo( BLIS_UPPER, C_1_T_1 );
	obj_t C_1_T_1_R, packed_B_pan_R;
	dim_t offR = bli_max( 0, bli_obj_diag_offset_after_trans( C_1_T_1 ) );
	dim_t nR = bli_obj_width_after_trans( C_1_T_1 ) - offR;
	bli_acquire_mpart_l2r( BLIS_SUBPART1,
			       offR, nR, &C_1_T_1, &C_1_T_1_R );
	bli_acquire_mpart_l2r( BLIS_SUBPART1,
			       offR, nR, &packed_B_pan, &packed_B_pan_R );
	bli_packm_blk_var2_par( &BLIS_ONE, &AT_T_1_1, &packed_A_blk, L2Comm );
	bli_herk_u_ker_var2_par( &BLIS_ONE, &packed_A_blk, &packed_B_pan_R, 
				 &BLIS_ONE, &C_1_T_1_R, (herk_t*)NULL , L1Comm );

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
  obj_t a, b, c1, c2;
	obj_t c_save;
	obj_t alpha, beta, normVal;
	dim_t m, k;
	dim_t p;
	dim_t p_begin, p_end, p_inc;
	int   m_input, k_input;
	num_t dt_a, dt_b, dt_c;
	num_t dt_alpha, dt_beta;
	int   r, n_repeats;

	double dtime;
	double dtime_save;
	double gflops;

	int lower=0, trans=0;

	if (argc != 3) {
	  printf("test L/U N/T\n");
	  fflush(stdout);
	  return 0;
	}

	if (*(argv[1]) == 'L')
	  lower = 1;
	else if (*(argv[1]) != 'U') {
	  printf("lower/upper not correct\n");
	  return 0;
	}
	else
	  lower = 0;

	if (*(argv[2]) == 'T')
	  trans = 1;
	else if (*(argv[2]) != 'N') {
	  printf("trans/normal not correct\n");
	  return 0;
	}
	else
	  trans = 0;

	bli_init();

	n_repeats = 3;

	p_begin = 40;
	p_end   = 400;
	p_inc   = 40;

	m_input = -1;
	k_input = -1;

	dt_a = BLIS_DOUBLE;
	dt_b = BLIS_DOUBLE;
	dt_c = BLIS_DOUBLE;
	dt_alpha = BLIS_DOUBLE;
	dt_beta = BLIS_DOUBLE;

	for ( p = p_begin; p <= p_end; p += p_inc )
	{

		if ( m_input < 0 ) m = p * ( dim_t )abs(m_input);
		else               m =     ( dim_t )    m_input;
		if ( k_input < 0 ) k = p * ( dim_t )abs(k_input);
		else               k =     ( dim_t )    k_input;


		bli_obj_create( dt_alpha, 1, 1, 0, 0, &alpha );
		bli_obj_create( dt_beta,  1, 1, 0, 0, &beta );
		bli_obj_create( dt_beta,  1, 1, 0, 0, &normVal );

		bli_obj_create( dt_a, m, k, 0, 0, &a );
		bli_obj_create( dt_a, m, k, 0, 0, &b );
		bli_obj_create( dt_c, m, m, 0, 0, &c1 );
		bli_obj_create( dt_c, m, m, 0, 0, &c2 );
		bli_obj_create( dt_c, m, m, 0, 0, &c_save );

		bli_randm( &a );
		bli_randm( &b );
		bli_randm( &c1 );

		bli_setsc(  (1.0/1.0), 0.0, &alpha );
		bli_setsc(  (1.0/1.0), 0.0, &beta );
		bli_setsc(  (1.0/1.0), 0.0, &beta );

		bli_copym( &c1, &c_save );
	
		dtime_save = 1.0e9;

		for ( r = 0; r < n_repeats; ++r )
		{
			bli_copym( &c_save, &c1 );
			bli_copym( &c_save, &c2 );


			dtime = bli_clock();

			th_setup_comm(&global_comm[0], NUMTHREADS, 1);
			_Pragma( "omp parallel num_threads(NUMTHREADS)" ) 
			  {
			    if (trans) {
			      if (lower) {
				DxT_Her2kTL(&alpha,
					    &a,
					    &b,
					    &beta,
					    &c1);				       
			      }
			      else { //upper
				DxT_Her2kTU(&alpha,
					    &a,
					    &b,
					    &beta,
					    &c1);				       
			      }
			    }
			    else { //normal
			      if (lower) {
				DxT_Her2kNL(&alpha,
					    &a,
					    &b,
					    &beta,
					    &c1);				       
			      }
			      else { //upper
				DxT_Her2kNU(&alpha,
					    &a,
					    &b,
					    &beta,
					    &c1);				       
			      }
			    }
			  }



			if (trans) {
			  bli_obj_set_conjtrans( BLIS_TRANSPOSE, a);
			  bli_obj_set_conjtrans( BLIS_TRANSPOSE, b);
			}
			  
			bli_obj_set_struc( BLIS_SYMMETRIC, c2 );
			if (lower) {
			  bli_obj_set_uplo( BLIS_LOWER, c2 );
			}
			else {
			  bli_obj_set_uplo( BLIS_UPPER, c2 );
			}
			
			bli_syr2k( &alpha,
				   &a,
				   &b,
				   &beta,
				   &c2);

			if (trans) {
			  bli_obj_toggle_trans( a );
			  bli_obj_toggle_trans( b );
			}

			bli_obj_set_struc( BLIS_GENERAL, c2 );
			bli_obj_set_uplo( BLIS_DENSE, c2 );
			
			
			bli_axpym( &BLIS_MINUS_ONE, &c1, &c2 );
			
		    bli_fnormm( &c2, &normVal );			
		    
		    dtime_save = bli_clock_min_diff( dtime_save, dtime );
		}

		gflops = ( 2.0 * m * k * m ) / ( dtime_save * 1.0e9 );

		printf( "data_her2k_blis" );
		printf( "( %2ld, 1:4 ) = [ %4lu %4lu  %10.3e  %6.3f ];\n",
		        (p - p_begin + 1)/p_inc + 1, m, k, dtime_save, gflops );

		bli_printm( "NORM", &normVal, "%4.1f", "" );

		bli_obj_free( &alpha );
		bli_obj_free( &beta );
		bli_obj_free( &normVal );

		bli_obj_free( &a );
		bli_obj_free( &b );
		bli_obj_free( &c1 );
		bli_obj_free( &c2 );
		bli_obj_free( &c_save );
	}

	bli_finalize();

	return 0;
}

