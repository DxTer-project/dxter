#include <unistd.h>
#include "blis.h"

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


void DxT_HemmLL( obj_t *alpha,
		 obj_t *A,
		 obj_t *B,
		 obj_t *beta,
		 obj_t *C )
{
  obj_t B_1_packed;
  obj_t A_11_1_packed;
  obj_t A_21_1_packed;
  obj_t A_10_1T_packed;

  bli_obj_init_pack( &B_1_packed );
  bli_obj_init_pack( &A_11_1_packed );
  bli_obj_init_pack( &A_21_1_packed );
  bli_obj_init_pack( &A_10_1T_packed );
  
  bli_scalm(beta, C);

  dim_t idx1, dimLen1, bs1;
///// Blocksize = 256
dimLen1 = bli_obj_length_after_trans( *C );
for ( idx1 = 0; idx1 < dimLen1; idx1 += bs1 ) {
	bs1 = bli_determine_blocksize_f( idx1, dimLen1, C, gemm_kc );
	dim_t idx2, dimLen2, bs2;
//****
	obj_t A_10;
	bli_acquire_mpart_tl2br( BLIS_SUBPART10, idx1, bs1, A, &A_10 );
	obj_t A_11;
	bli_acquire_mpart_tl2br( BLIS_SUBPART11, idx1, bs1, A, &A_11 );
	obj_t A_21;
	bli_acquire_mpart_tl2br( BLIS_SUBPART21, idx1, bs1, A, &A_21 );
	obj_t B_1;
	bli_acquire_mpart_t2b( BLIS_SUBPART1, idx1, bs1, B, &B_1 );
	obj_t C_0;
	bli_acquire_mpart_t2b( BLIS_SUBPART0, idx1, bs1, C, &C_0 );
	obj_t C_1;
	bli_acquire_mpart_t2b( BLIS_SUBPART1, idx1, bs1, C, &C_1 );
	obj_t C_2;
	bli_acquire_mpart_t2b( BLIS_SUBPART2, idx1, bs1, C, &C_2 );
	//------------------------------------//

	bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_COL_PANELS, 
			BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			BLIS_BUFFER_FOR_B_PANEL,
			gemm_kr,  gemm_nr,  
			&B_1, &B_1_packed );
	bli_packm_blk_var2( &BLIS_ONE, &B_1, &B_1_packed );
	///// Blocksize = 128
	dimLen2 = bli_obj_length_after_trans( C_1 );
	for ( idx2 = 0; idx2 < dimLen2; idx2 += bs2 ) {
		bs2 = bli_determine_blocksize_f( idx2, dimLen2, &C_1, gemm_mc );
		dim_t idx3, dimLen3, bs3;
	//****
		obj_t A_11_1;
		bli_acquire_mpart_t2b( BLIS_SUBPART1, idx2, bs2, &A_11, &A_11_1 );
		obj_t C_1_1;
		bli_acquire_mpart_t2b( BLIS_SUBPART1, idx2, bs2, &C_1, &C_1_1 );
		//------------------------------------//

		bli_obj_set_struc( BLIS_SYMMETRIC, A_11_1 );
		bli_obj_set_uplo( BLIS_LOWER, A_11_1 );
		bli_packm_init_pack( TRUE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_ROW_PANELS, 
				BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
				BLIS_BUFFER_FOR_A_BLOCK,
				gemm_mr,  gemm_kr,  
				&A_11_1, &A_11_1_packed );
		bli_packm_blk_var2( &BLIS_ONE, &A_11_1, &A_11_1_packed );
		bli_gemm_ker_var2( &BLIS_ONE, &A_11_1_packed, &B_1_packed, 
				&BLIS_ONE, &C_1_1, (gemm_t*)NULL );

		//------------------------------------//

	//****
	}
	///// Blocksize = 128
	dimLen2 = bli_obj_length_after_trans( C_0 );
	for ( idx2 = 0; idx2 < dimLen2; idx2 += bs2 ) {
		bs2 = bli_determine_blocksize_f( idx2, dimLen2, &C_0, gemm_mc );
		dim_t idx3, dimLen3, bs3;
	//****
		obj_t A_10_1;
		bli_acquire_mpart_l2r( BLIS_SUBPART1, idx2, bs2, &A_10, &A_10_1 );
		obj_t C_0_1;
		bli_acquire_mpart_t2b( BLIS_SUBPART1, idx2, bs2, &C_0, &C_0_1 );
		//------------------------------------//

		obj_t A_10_1T;
		bli_obj_alias_with_trans( BLIS_TRANSPOSE, A_10_1, A_10_1T);
		bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_ROW_PANELS, 
				BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
				BLIS_BUFFER_FOR_A_BLOCK,
				gemm_mr,  gemm_kr,  
				&A_10_1T, &A_10_1T_packed );
		bli_packm_blk_var2( &BLIS_ONE, &A_10_1T, &A_10_1T_packed );
		bli_gemm_ker_var2( &BLIS_ONE, &A_10_1T_packed, &B_1_packed, 
				&BLIS_ONE, &C_0_1, (gemm_t*)NULL );

		//------------------------------------//

	//****
	}
	///// Blocksize = 128
	dimLen2 = bli_obj_length_after_trans( C_2 );
	for ( idx2 = 0; idx2 < dimLen2; idx2 += bs2 ) {
		bs2 = bli_determine_blocksize_f( idx2, dimLen2, &C_2, gemm_mc );
		dim_t idx3, dimLen3, bs3;
	//****
		obj_t A_21_1;
		bli_acquire_mpart_t2b( BLIS_SUBPART1, idx2, bs2, &A_21, &A_21_1 );
		obj_t C_2_1;
		bli_acquire_mpart_t2b( BLIS_SUBPART1, idx2, bs2, &C_2, &C_2_1 );
		//------------------------------------//

		bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_ROW_PANELS, 
				BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
				BLIS_BUFFER_FOR_A_BLOCK,
				gemm_mr,  gemm_kr,  
				&A_21_1, &A_21_1_packed );
		bli_packm_blk_var2( &BLIS_ONE, &A_21_1, &A_21_1_packed );
		bli_gemm_ker_var2( &BLIS_ONE, &A_21_1_packed, &B_1_packed, 
				&BLIS_ONE, &C_2_1, (gemm_t*)NULL );

		//------------------------------------//

	//****
	}

	//------------------------------------//

//****
}
  
 

  bli_obj_release_pack( &B_1_packed );
  bli_obj_release_pack( &A_11_1_packed );
  bli_obj_release_pack( &A_21_1_packed );
  bli_obj_release_pack( &A_10_1T_packed );
}

void DxT_HemmLU( obj_t *alpha,
		 obj_t *A,
		 obj_t *B,
		 obj_t *beta,
		 obj_t *C )
{
  obj_t B_1_packed;
  obj_t A_11_1_packed;
  obj_t A_12_1T_packed;
  obj_t A_01_1_packed;
  
  bli_obj_init_pack( &A_12_1T_packed );
  bli_obj_init_pack( &A_01_1_packed );
  bli_obj_init_pack( &B_1_packed );
  bli_obj_init_pack( &A_11_1_packed );

  bli_scalm(beta, C);

dim_t idx1, dimLen1, bs1;
///// Blocksize = 256
dimLen1 = bli_obj_length_after_trans( *C );
for ( idx1 = 0; idx1 < dimLen1; idx1 += bs1 ) {
	bs1 = bli_determine_blocksize_f( idx1, dimLen1, C, gemm_kc );
	dim_t idx2, dimLen2, bs2;
//****
	obj_t A_01;
	bli_acquire_mpart_tl2br( BLIS_SUBPART01, idx1, bs1, A, &A_01 );
	obj_t A_11;
	bli_acquire_mpart_tl2br( BLIS_SUBPART11, idx1, bs1, A, &A_11 );
	obj_t A_12;
	bli_acquire_mpart_tl2br( BLIS_SUBPART12, idx1, bs1, A, &A_12 );
	obj_t B_1;
	bli_acquire_mpart_t2b( BLIS_SUBPART1, idx1, bs1, B, &B_1 );
	obj_t C_0;
	bli_acquire_mpart_t2b( BLIS_SUBPART0, idx1, bs1, C, &C_0 );
	obj_t C_1;
	bli_acquire_mpart_t2b( BLIS_SUBPART1, idx1, bs1, C, &C_1 );
	obj_t C_2;
	bli_acquire_mpart_t2b( BLIS_SUBPART2, idx1, bs1, C, &C_2 );
	//------------------------------------//

	bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_COL_PANELS, 
			BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			BLIS_BUFFER_FOR_B_PANEL,
			gemm_kr,  gemm_nr,  
			&B_1, &B_1_packed );
	bli_packm_blk_var2( &BLIS_ONE, &B_1, &B_1_packed );
	///// Blocksize = 128
	dimLen2 = bli_obj_length_after_trans( C_2 );
	for ( idx2 = 0; idx2 < dimLen2; idx2 += bs2 ) {
		bs2 = bli_determine_blocksize_f( idx2, dimLen2, &C_2, gemm_mc );
		dim_t idx3, dimLen3, bs3;
	//****
		obj_t A_12_1;
		bli_acquire_mpart_l2r( BLIS_SUBPART1, idx2, bs2, &A_12, &A_12_1 );
		obj_t C_2_1;
		bli_acquire_mpart_t2b( BLIS_SUBPART1, idx2, bs2, &C_2, &C_2_1 );
		//------------------------------------//

		obj_t A_12_1T;
		bli_obj_alias_with_trans( BLIS_TRANSPOSE, A_12_1, A_12_1T);
		bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_ROW_PANELS, 
				BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
				BLIS_BUFFER_FOR_A_BLOCK,
				gemm_mr,  gemm_kr,  
				&A_12_1T, &A_12_1T_packed );
		bli_packm_blk_var2( &BLIS_ONE, &A_12_1T, &A_12_1T_packed );
		bli_gemm_ker_var2( &BLIS_ONE, &A_12_1T_packed, &B_1_packed, 
				&BLIS_ONE, &C_2_1, (gemm_t*)NULL );

		//------------------------------------//

	//****
	}
	///// Blocksize = 128
	dimLen2 = bli_obj_length_after_trans( C_1 );
	for ( idx2 = 0; idx2 < dimLen2; idx2 += bs2 ) {
		bs2 = bli_determine_blocksize_f( idx2, dimLen2, &C_1, gemm_mc );
		dim_t idx3, dimLen3, bs3;
	//****
		obj_t A_11_1;
		bli_acquire_mpart_t2b( BLIS_SUBPART1, idx2, bs2, &A_11, &A_11_1 );
		obj_t C_1_1;
		bli_acquire_mpart_t2b( BLIS_SUBPART1, idx2, bs2, &C_1, &C_1_1 );
		//------------------------------------//

		bli_obj_set_struc( BLIS_SYMMETRIC, A_11_1 );
		bli_obj_set_uplo( BLIS_UPPER, A_11_1 );
		bli_packm_init_pack( TRUE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_ROW_PANELS, 
				BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
				BLIS_BUFFER_FOR_A_BLOCK,
				gemm_mr,  gemm_kr,  
				&A_11_1, &A_11_1_packed );
		bli_packm_blk_var2( &BLIS_ONE, &A_11_1, &A_11_1_packed );
		bli_gemm_ker_var2( &BLIS_ONE, &A_11_1_packed, &B_1_packed, 
				&BLIS_ONE, &C_1_1, (gemm_t*)NULL );

		//------------------------------------//

	//****
	}
	///// Blocksize = 128
	dimLen2 = bli_obj_length_after_trans( C_0 );
	for ( idx2 = 0; idx2 < dimLen2; idx2 += bs2 ) {
		bs2 = bli_determine_blocksize_f( idx2, dimLen2, &C_0, gemm_mc );
		dim_t idx3, dimLen3, bs3;
	//****
		obj_t A_01_1;
		bli_acquire_mpart_t2b( BLIS_SUBPART1, idx2, bs2, &A_01, &A_01_1 );
		obj_t C_0_1;
		bli_acquire_mpart_t2b( BLIS_SUBPART1, idx2, bs2, &C_0, &C_0_1 );
		//------------------------------------//

		bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_ROW_PANELS, 
				BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
				BLIS_BUFFER_FOR_A_BLOCK,
				gemm_mr,  gemm_kr,  
				&A_01_1, &A_01_1_packed );
		bli_packm_blk_var2( &BLIS_ONE, &A_01_1, &A_01_1_packed );
		bli_gemm_ker_var2( &BLIS_ONE, &A_01_1_packed, &B_1_packed, 
				&BLIS_ONE, &C_0_1, (gemm_t*)NULL );

		//------------------------------------//

	//****
	}

	//------------------------------------//

//****
}

  bli_obj_release_pack( &A_12_1T_packed );
  bli_obj_release_pack( &A_11_1_packed );
  bli_obj_release_pack( &B_1_packed );
  bli_obj_release_pack( &A_01_1_packed );
}

void DxT_HemmRL( obj_t *alpha,
		 obj_t *A,
		 obj_t *B,
		 obj_t *beta,
		 obj_t *C )
{
  obj_t BT_1_packed;
  obj_t AC_11_1_packed;
  obj_t AC_21_1_packed;
  obj_t AC_10_1T_packed;

  bli_obj_init_pack( &AC_21_1_packed );
  bli_obj_init_pack( &AC_10_1T_packed );
  bli_obj_init_pack( &BT_1_packed );
  bli_obj_init_pack( &AC_11_1_packed );
  
  bli_scalm(beta, C);

dim_t idx1, dimLen1, bs1;
obj_t AC;
bli_obj_alias_with_conj( BLIS_CONJUGATE, *A, AC);
obj_t BT;
bli_obj_alias_with_trans( BLIS_TRANSPOSE, *B, BT);
bli_obj_induce_trans( *C );
///// Blocksize = 256
dimLen1 = bli_obj_length_after_trans( *C );
for ( idx1 = 0; idx1 < dimLen1; idx1 += bs1 ) {
	bs1 = bli_determine_blocksize_f( idx1, dimLen1, C, gemm_kc );
	dim_t idx2, dimLen2, bs2;
//****
	obj_t AC_10;
	bli_acquire_mpart_tl2br( BLIS_SUBPART10, idx1, bs1, &AC, &AC_10 );
	obj_t AC_11;
	bli_acquire_mpart_tl2br( BLIS_SUBPART11, idx1, bs1, &AC, &AC_11 );
	obj_t AC_21;
	bli_acquire_mpart_tl2br( BLIS_SUBPART21, idx1, bs1, &AC, &AC_21 );
	obj_t BT_1;
	bli_acquire_mpart_t2b( BLIS_SUBPART1, idx1, bs1, &BT, &BT_1 );
	obj_t C_0;
	bli_acquire_mpart_t2b( BLIS_SUBPART0, idx1, bs1, C, &C_0 );
	obj_t C_1;
	bli_acquire_mpart_t2b( BLIS_SUBPART1, idx1, bs1, C, &C_1 );
	obj_t C_2;
	bli_acquire_mpart_t2b( BLIS_SUBPART2, idx1, bs1, C, &C_2 );
	//------------------------------------//

	bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_COL_PANELS, 
			BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			BLIS_BUFFER_FOR_B_PANEL,
			gemm_kr,  gemm_nr,  
			&BT_1, &BT_1_packed );
	bli_packm_blk_var2( &BLIS_ONE, &BT_1, &BT_1_packed );
	///// Blocksize = 128
	dimLen2 = bli_obj_length_after_trans( C_0 );
	for ( idx2 = 0; idx2 < dimLen2; idx2 += bs2 ) {
		bs2 = bli_determine_blocksize_f( idx2, dimLen2, &C_0, gemm_mc );
		dim_t idx3, dimLen3, bs3;
	//****
		obj_t AC_10_1;
		bli_acquire_mpart_l2r( BLIS_SUBPART1, idx2, bs2, &AC_10, &AC_10_1 );
		obj_t C_0_1;
		bli_acquire_mpart_t2b( BLIS_SUBPART1, idx2, bs2, &C_0, &C_0_1 );
		//------------------------------------//

		obj_t AC_10_1T;
		bli_obj_alias_with_trans( BLIS_TRANSPOSE, AC_10_1, AC_10_1T);
		bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_ROW_PANELS, 
				BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
				BLIS_BUFFER_FOR_A_BLOCK,
				gemm_mr,  gemm_kr,  
				&AC_10_1T, &AC_10_1T_packed );
		bli_packm_blk_var2( &BLIS_ONE, &AC_10_1T, &AC_10_1T_packed );
		bli_gemm_ker_var2( &BLIS_ONE, &AC_10_1T_packed, &BT_1_packed, 
				&BLIS_ONE, &C_0_1, (gemm_t*)NULL );

		//------------------------------------//

	//****
	}
	///// Blocksize = 128
	dimLen2 = bli_obj_length_after_trans( C_1 );
	for ( idx2 = 0; idx2 < dimLen2; idx2 += bs2 ) {
		bs2 = bli_determine_blocksize_f( idx2, dimLen2, &C_1, gemm_mc );
		dim_t idx3, dimLen3, bs3;
	//****
		obj_t AC_11_1;
		bli_acquire_mpart_t2b( BLIS_SUBPART1, idx2, bs2, &AC_11, &AC_11_1 );
		obj_t C_1_1;
		bli_acquire_mpart_t2b( BLIS_SUBPART1, idx2, bs2, &C_1, &C_1_1 );
		//------------------------------------//

		bli_obj_set_struc( BLIS_SYMMETRIC, AC_11_1 );
		bli_obj_set_uplo( BLIS_LOWER, AC_11_1 );
		bli_packm_init_pack( TRUE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_ROW_PANELS, 
				BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
				BLIS_BUFFER_FOR_A_BLOCK,
				gemm_mr,  gemm_kr,  
				&AC_11_1, &AC_11_1_packed );
		bli_packm_blk_var2( &BLIS_ONE, &AC_11_1, &AC_11_1_packed );
		bli_gemm_ker_var2( &BLIS_ONE, &AC_11_1_packed, &BT_1_packed, 
				&BLIS_ONE, &C_1_1, (gemm_t*)NULL );

		//------------------------------------//

	//****
	}
	///// Blocksize = 128
	dimLen2 = bli_obj_length_after_trans( C_2 );
	for ( idx2 = 0; idx2 < dimLen2; idx2 += bs2 ) {
		bs2 = bli_determine_blocksize_f( idx2, dimLen2, &C_2, gemm_mc );
		dim_t idx3, dimLen3, bs3;
	//****
		obj_t AC_21_1;
		bli_acquire_mpart_t2b( BLIS_SUBPART1, idx2, bs2, &AC_21, &AC_21_1 );
		obj_t C_2_1;
		bli_acquire_mpart_t2b( BLIS_SUBPART1, idx2, bs2, &C_2, &C_2_1 );
		//------------------------------------//

		bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_ROW_PANELS, 
				BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
				BLIS_BUFFER_FOR_A_BLOCK,
				gemm_mr,  gemm_kr,  
				&AC_21_1, &AC_21_1_packed );
		bli_packm_blk_var2( &BLIS_ONE, &AC_21_1, &AC_21_1_packed );
		bli_gemm_ker_var2( &BLIS_ONE, &AC_21_1_packed, &BT_1_packed, 
				&BLIS_ONE, &C_2_1, (gemm_t*)NULL );

		//------------------------------------//

	//****
	}

	//------------------------------------//

//****
}
bli_obj_induce_trans( *C );


  bli_obj_release_pack( &AC_21_1_packed );
  bli_obj_release_pack( &AC_11_1_packed );
  bli_obj_release_pack( &BT_1_packed );
  bli_obj_release_pack( &AC_10_1T_packed );
}

void DxT_HemmRU( obj_t *alpha,
		 obj_t *A,
		 obj_t *B,
		 obj_t *beta,
		 obj_t *C )
{
  obj_t BT_1_packed;
  obj_t AC_12_1T_packed;
  obj_t AC_01_1_packed;
  obj_t AC_11_1_packed;
  
  bli_obj_init_pack( &BT_1_packed );
  bli_obj_init_pack( &AC_12_1T_packed );
  bli_obj_init_pack( &AC_01_1_packed );
  bli_obj_init_pack( &AC_11_1_packed );
    
  bli_scalm(beta, C);

dim_t idx1, dimLen1, bs1;
obj_t AC;
bli_obj_alias_with_conj( BLIS_CONJUGATE, *A, AC);
obj_t BT;
bli_obj_alias_with_trans( BLIS_TRANSPOSE, *B, BT);
bli_obj_induce_trans( *C );
///// Blocksize = 256
dimLen1 = bli_obj_length_after_trans( *C );
for ( idx1 = 0; idx1 < dimLen1; idx1 += bs1 ) {
	bs1 = bli_determine_blocksize_f( idx1, dimLen1, C, gemm_kc );
	dim_t idx2, dimLen2, bs2;
//****
	obj_t AC_01;
	bli_acquire_mpart_tl2br( BLIS_SUBPART01, idx1, bs1, &AC, &AC_01 );
	obj_t AC_11;
	bli_acquire_mpart_tl2br( BLIS_SUBPART11, idx1, bs1, &AC, &AC_11 );
	obj_t AC_12;
	bli_acquire_mpart_tl2br( BLIS_SUBPART12, idx1, bs1, &AC, &AC_12 );
	obj_t BT_1;
	bli_acquire_mpart_t2b( BLIS_SUBPART1, idx1, bs1, &BT, &BT_1 );
	obj_t C_0;
	bli_acquire_mpart_t2b( BLIS_SUBPART0, idx1, bs1, C, &C_0 );
	obj_t C_1;
	bli_acquire_mpart_t2b( BLIS_SUBPART1, idx1, bs1, C, &C_1 );
	obj_t C_2;
	bli_acquire_mpart_t2b( BLIS_SUBPART2, idx1, bs1, C, &C_2 );
	//------------------------------------//

	bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_COL_PANELS, 
			BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			BLIS_BUFFER_FOR_B_PANEL,
			gemm_kr,  gemm_nr,  
			&BT_1, &BT_1_packed );
	bli_packm_blk_var2( &BLIS_ONE, &BT_1, &BT_1_packed );
	///// Blocksize = 128
	dimLen2 = bli_obj_length_after_trans( C_0 );
	for ( idx2 = 0; idx2 < dimLen2; idx2 += bs2 ) {
		bs2 = bli_determine_blocksize_f( idx2, dimLen2, &C_0, gemm_mc );
		dim_t idx3, dimLen3, bs3;
	//****
		obj_t AC_01_1;
		bli_acquire_mpart_t2b( BLIS_SUBPART1, idx2, bs2, &AC_01, &AC_01_1 );
		obj_t C_0_1;
		bli_acquire_mpart_t2b( BLIS_SUBPART1, idx2, bs2, &C_0, &C_0_1 );
		//------------------------------------//

		bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_ROW_PANELS, 
				BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
				BLIS_BUFFER_FOR_A_BLOCK,
				gemm_mr,  gemm_kr,  
				&AC_01_1, &AC_01_1_packed );
		bli_packm_blk_var2( &BLIS_ONE, &AC_01_1, &AC_01_1_packed );
		bli_gemm_ker_var2( &BLIS_ONE, &AC_01_1_packed, &BT_1_packed, 
				&BLIS_ONE, &C_0_1, (gemm_t*)NULL );

		//------------------------------------//

	//****
	}
	///// Blocksize = 128
	dimLen2 = bli_obj_length_after_trans( C_1 );
	for ( idx2 = 0; idx2 < dimLen2; idx2 += bs2 ) {
		bs2 = bli_determine_blocksize_f( idx2, dimLen2, &C_1, gemm_mc );
		dim_t idx3, dimLen3, bs3;
	//****
		obj_t AC_11_1;
		bli_acquire_mpart_t2b( BLIS_SUBPART1, idx2, bs2, &AC_11, &AC_11_1 );
		obj_t C_1_1;
		bli_acquire_mpart_t2b( BLIS_SUBPART1, idx2, bs2, &C_1, &C_1_1 );
		//------------------------------------//

		bli_obj_set_struc( BLIS_SYMMETRIC, AC_11_1 );
		bli_obj_set_uplo( BLIS_UPPER, AC_11_1 );
		bli_packm_init_pack( TRUE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_ROW_PANELS, 
				BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
				BLIS_BUFFER_FOR_A_BLOCK,
				gemm_mr,  gemm_kr,  
				&AC_11_1, &AC_11_1_packed );
		bli_packm_blk_var2( &BLIS_ONE, &AC_11_1, &AC_11_1_packed );
		bli_gemm_ker_var2( &BLIS_ONE, &AC_11_1_packed, &BT_1_packed, 
				&BLIS_ONE, &C_1_1, (gemm_t*)NULL );

		//------------------------------------//

	//****
	}
	///// Blocksize = 128
	dimLen2 = bli_obj_length_after_trans( C_2 );
	for ( idx2 = 0; idx2 < dimLen2; idx2 += bs2 ) {
		bs2 = bli_determine_blocksize_f( idx2, dimLen2, &C_2, gemm_mc );
		dim_t idx3, dimLen3, bs3;
	//****
		obj_t AC_12_1;
		bli_acquire_mpart_l2r( BLIS_SUBPART1, idx2, bs2, &AC_12, &AC_12_1 );
		obj_t C_2_1;
		bli_acquire_mpart_t2b( BLIS_SUBPART1, idx2, bs2, &C_2, &C_2_1 );
		//------------------------------------//

		obj_t AC_12_1T;
		bli_obj_alias_with_trans( BLIS_TRANSPOSE, AC_12_1, AC_12_1T);
		bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_ROW_PANELS, 
				BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
				BLIS_BUFFER_FOR_A_BLOCK,
				gemm_mr,  gemm_kr,  
				&AC_12_1T, &AC_12_1T_packed );
		bli_packm_blk_var2( &BLIS_ONE, &AC_12_1T, &AC_12_1T_packed );
		bli_gemm_ker_var2( &BLIS_ONE, &AC_12_1T_packed, &BT_1_packed, 
				&BLIS_ONE, &C_2_1, (gemm_t*)NULL );

		//------------------------------------//

	//****
	}

	//------------------------------------//

//****
}
bli_obj_induce_trans( *C );


      bli_obj_release_pack( &BT_1_packed );
      bli_obj_release_pack( &AC_11_1_packed );
      bli_obj_release_pack( &AC_12_1T_packed );
      bli_obj_release_pack( &AC_01_1_packed );

}

int main( int argc, char** argv )
{
  obj_t a, b, c1, c2;
	obj_t c_save;
	obj_t alpha, beta, normVal;
	dim_t m, n;
	dim_t p;
	dim_t p_begin, p_end, p_inc;
	int   m_input, n_input;
	num_t dt_a, dt_b, dt_c;
	num_t dt_alpha, dt_beta;
	int   r, n_repeats;
	side_t side;

	double dtime;
	double dtime_save;
	double gflops;

	int left, lower;

	bli_init();

	n_repeats = 3;

	p_begin = 40;
	p_end   = 600;
	p_inc   = 40;

	m_input = -1;
	n_input = -1;

	dt_a = BLIS_DOUBLE;
	dt_b = BLIS_DOUBLE;
	dt_c = BLIS_DOUBLE;
	dt_alpha = BLIS_DOUBLE;
	dt_beta = BLIS_DOUBLE;

	if (argc != 3) {
	  printf("test L/R L/U\n");
	  fflush(stdout);
	  return;
	}

	if (*(argv[1]) == 'L')
	  left = 1;
	else if (*(argv[1]) != 'R') {
	  printf("left/right not correct\n");
	  return;
	}
	else
	  left = 0;

	if (*(argv[2]) == 'L')
	  lower = 1;
	else if (*(argv[2]) != 'U') {
	  printf("lower/upper not correct\n");
	  return;
	}
	else
	  lower = 0;


	for ( p = p_begin; p <= p_end; p += p_inc )
	{

		if ( m_input < 0 ) m = p * ( dim_t )abs(m_input);
		else               m =     ( dim_t )    m_input;
		if ( n_input < 0 ) n = p * ( dim_t )abs(n_input);
		else               n =     ( dim_t )    n_input;


		bli_obj_create( dt_alpha, 1, 1, 0, 0, &alpha );
		bli_obj_create( dt_beta,  1, 1, 0, 0, &beta );
		bli_obj_create( dt_alpha, 1, 1, 0, 0, &normVal );
	

		if ( left )
		  bli_obj_create( dt_a, m, m, 0, 0, &a );
		else
		  bli_obj_create( dt_a, n, n, 0, 0, &a );
		bli_obj_create( dt_b, m, n, 0, 0, &b );
		bli_obj_create( dt_c, m, n, 0, 0, &c1 );
		bli_obj_create( dt_c, m, n, 0, 0, &c2 );
		bli_obj_create( dt_c, m, n, 0, 0, &c_save );

		bli_randm( &a );
		bli_randm( &b );
		bli_randm( &c1 );


		bli_setsc( 1.0, 0.0, &alpha );
		bli_setsc( 1.0, 0.0, &beta );

		bli_copym( &c1, &c_save );
	
		dtime_save = 1.0e9;

		for ( r = 0; r < n_repeats; ++r )
		{
			bli_copym( &c_save, &c1 );
			bli_copym( &c_save, &c2 );


			dtime = bli_clock();

			bli_error_checking_level_set( BLIS_NO_ERROR_CHECKING );
			
			bli_obj_set_struc( BLIS_HERMITIAN, a );
			if (lower) {
			  bli_obj_set_uplo( BLIS_LOWER, a );
			}
			else {
			  bli_obj_set_uplo( BLIS_UPPER, a );
			}
			
			bli_hemm( left ? BLIS_LEFT : BLIS_RIGHT,
				  &alpha,
				  &a,
				  &b,
				  &beta,
				  &c1);

			bli_obj_set_struc( BLIS_GENERAL, a );
			bli_obj_set_uplo( BLIS_DENSE, a );


			if (left) {
			  if (lower) {
			    DxT_HemmLL( &alpha,
				  &a,
				  &b,
				  &beta,
				  &c2);
			  }
			  else { //upper
			    DxT_HemmLU( &alpha,
					&a,
					&b,
					&beta,
					&c2);
			  }
			}
			else { //right
			  if (lower) {
			    DxT_HemmRL( &alpha,
				  &a,
				  &b,
				  &beta,
				  &c2);
			  }
			  else { //upper
			    DxT_HemmRU( &alpha,
				  &a,
				  &b,
				  &beta,
				  &c2);
			  }
			}

			bli_axpym( &BLIS_MINUS_ONE, &c1, &c2 );
			
			bli_fnormm( &c2, &normVal );

			dtime_save = bli_clock_min_diff( dtime_save, dtime );
		}

		if ( left )
			gflops = ( 2.0 * m * m * n ) / ( dtime_save * 1.0e9 );
		else
			gflops = ( 2.0 * m * n * n ) / ( dtime_save * 1.0e9 );

		printf( "data_hemm_blis" );
		printf( "( %2ld, 1:4 ) = [ %4lu %4lu  %10.3e  %6.3f ];\n",
		        (p - p_begin + 1)/p_inc + 1, m, n, dtime_save, gflops );

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

