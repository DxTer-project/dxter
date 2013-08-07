#include <unistd.h>
#include "blis.h"
#include "support.h"
#include "bli_trmm_l_ker_var2_par.h"
#include "bli_trmm_u_ker_var2_par.h"


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



void DxT_TrmmLLN( obj_t *alpha,
		  obj_t *L,
		  obj_t *X )
{
  obj_t L_11_packed;
  obj_t L_21_1_packed;
  obj_t X_1_packed;
    bli_obj_init_pack( &L_11_packed );
    bli_obj_init_pack( &L_21_1_packed );
    bli_obj_init_pack( &X_1_packed );

  bli_scalm(alpha, X);
dim_t idx1, dimLen1, bs1;
///// Blocksize = 256
dimLen1 = bli_obj_length_after_trans( *X );
for ( idx1 = 0; idx1 < dimLen1; idx1 += bs1 ) {
	bs1 = bli_determine_blocksize_b( idx1, dimLen1, X, gemm_kc );
	dim_t idx2, dimLen2, bs2;
//****
	obj_t L_11;
	bli_acquire_mpart_br2tl( BLIS_SUBPART11, idx1, bs1, L, &L_11 );
	obj_t L_21;
	bli_acquire_mpart_br2tl( BLIS_SUBPART21, idx1, bs1, L, &L_21 );
	obj_t X_1;
	bli_acquire_mpart_b2t( BLIS_SUBPART1, idx1, bs1, X, &X_1 );
	obj_t X_2;
	bli_acquire_mpart_b2t( BLIS_SUBPART2, idx1, bs1, X, &X_2 );
	//------------------------------------//

	bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_COL_PANELS, 
			BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			BLIS_BUFFER_FOR_B_PANEL,
			gemm_mr,  gemm_nr,  
			&X_1, &X_1_packed );
	bli_packm_blk_var2( &BLIS_ONE, &X_1, &X_1_packed );
	///// Blocksize = 128
	dimLen2 = bli_obj_length_after_trans( X_2 );
	for ( idx2 = 0; idx2 < dimLen2; idx2 += bs2 ) {
		bs2 = bli_determine_blocksize_f( idx2, dimLen2, &X_2, gemm_mc );
		dim_t idx3, dimLen3, bs3;
	//****
		obj_t L_21_1;
		bli_acquire_mpart_t2b( BLIS_SUBPART1, idx2, bs2, &L_21, &L_21_1 );
		obj_t X_2_1;
		bli_acquire_mpart_t2b( BLIS_SUBPART1, idx2, bs2, &X_2, &X_2_1 );
		//------------------------------------//

		bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_ROW_PANELS, 
				BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
				BLIS_BUFFER_FOR_A_BLOCK,
				gemm_mr,  gemm_mr,  
				&L_21_1, &L_21_1_packed );
		bli_packm_blk_var2( &BLIS_ONE, &L_21_1, &L_21_1_packed );
		bli_gemm_ker_var2( &BLIS_ONE, &L_21_1_packed, &X_1_packed, 
				&BLIS_ONE, &X_2_1, (gemm_t*)NULL );

		//------------------------------------//

	//****
	}
	///// Blocksize = 128
	dimLen2 = bli_obj_length_after_trans( L_11 );
	for ( idx2 = 0; idx2 < dimLen2; idx2 += bs2 ) {
		bs2 = bli_determine_blocksize_f( idx2, dimLen2, &L_11, gemm_mc );
		dim_t idx3, dimLen3, bs3;
	//****
		obj_t L_11_1;
		bli_acquire_mpart_t2b( BLIS_SUBPART1, idx2, bs2, &L_11, &L_11_1 );
		obj_t X_1_1;
		bli_acquire_mpart_t2b( BLIS_SUBPART1, idx2, bs2, &X_1, &X_1_1 );
		//------------------------------------//

		bli_obj_set_struc( BLIS_TRIANGULAR, L_11_1 );
		bli_obj_set_uplo( BLIS_LOWER, L_11_1 );
		bli_packm_init_pack( TRUE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_ROW_PANELS, 
				BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
				BLIS_BUFFER_FOR_A_BLOCK,
				gemm_mr,  gemm_mr,  
				&L_11_1, &L_11_packed );
		bli_packm_blk_var3( &BLIS_ONE, &L_11_1, &L_11_packed );
		bli_trmm_l_ker_var2( &BLIS_ONE, &L_11_packed, &X_1_packed, 
				&BLIS_ZERO, &X_1_1, (trmm_t*)NULL );

		//------------------------------------//

	//****
	}

	//------------------------------------//

//****
}



    bli_obj_release_pack( &L_11_packed );
    bli_obj_release_pack( &L_21_1_packed );
    bli_obj_release_pack( &X_1_packed );
}

void DxT_TrmmLLT( obj_t *alpha,
		  obj_t *L,
		  obj_t *X )
{
  obj_t L_11_packed;
  obj_t L_10_1T_packed;
  obj_t X_1_packed;
  bli_obj_init_pack( &L_11_packed );
  bli_obj_init_pack( &L_10_1T_packed );
  bli_obj_init_pack( &X_1_packed );
  
  bli_scalm(alpha, X);
dim_t idx1, dimLen1, bs1;
///// Blocksize = 256
dimLen1 = bli_obj_length_after_trans( *X );
for ( idx1 = 0; idx1 < dimLen1; idx1 += bs1 ) {
	bs1 = bli_determine_blocksize_f( idx1, dimLen1, X, gemm_kc );
	dim_t idx2, dimLen2, bs2;
//****
	obj_t L_10;
	bli_acquire_mpart_tl2br( BLIS_SUBPART10, idx1, bs1, L, &L_10 );
	obj_t L_11;
	bli_acquire_mpart_tl2br( BLIS_SUBPART11, idx1, bs1, L, &L_11 );
	obj_t X_0;
	bli_acquire_mpart_t2b( BLIS_SUBPART0, idx1, bs1, X, &X_0 );
	obj_t X_1;
	bli_acquire_mpart_t2b( BLIS_SUBPART1, idx1, bs1, X, &X_1 );
	//------------------------------------//

	bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_COL_PANELS, 
			BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			BLIS_BUFFER_FOR_B_PANEL,
			gemm_mr,  gemm_nr,  
			&X_1, &X_1_packed );
	bli_packm_blk_var2( &BLIS_ONE, &X_1, &X_1_packed );
	///// Blocksize = 128
	dimLen2 = bli_obj_width_after_trans( L_11 );
	for ( idx2 = 0; idx2 < dimLen2; idx2 += bs2 ) {
		bs2 = bli_determine_blocksize_f( idx2, dimLen2, &L_11, gemm_mc );
		dim_t idx3, dimLen3, bs3;
	//****
		obj_t L_11_1;
		bli_acquire_mpart_l2r( BLIS_SUBPART1, idx2, bs2, &L_11, &L_11_1 );
		obj_t X_1_1;
		bli_acquire_mpart_t2b( BLIS_SUBPART1, idx2, bs2, &X_1, &X_1_1 );
		//------------------------------------//

		obj_t L_11_1T;
		bli_obj_alias_with_trans( BLIS_TRANSPOSE, L_11_1, L_11_1T);
		bli_obj_set_struc( BLIS_TRIANGULAR, L_11_1T );
		bli_obj_set_uplo( BLIS_LOWER, L_11_1T );
		bli_packm_init_pack( TRUE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_ROW_PANELS, 
				BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
				BLIS_BUFFER_FOR_A_BLOCK,
				gemm_mr,  gemm_mr,  
				&L_11_1T, &L_11_packed );
		bli_packm_blk_var3( &BLIS_ONE, &L_11_1T, &L_11_packed );
		bli_trmm_u_ker_var2( &BLIS_ONE, &L_11_packed, &X_1_packed, 
				&BLIS_ZERO, &X_1_1, (trmm_t*)NULL );

		//------------------------------------//

	//****
	}
	///// Blocksize = 128
	dimLen2 = bli_obj_length_after_trans( X_0 );
	for ( idx2 = 0; idx2 < dimLen2; idx2 += bs2 ) {
		bs2 = bli_determine_blocksize_f( idx2, dimLen2, &X_0, gemm_mc );
		dim_t idx3, dimLen3, bs3;
	//****
		obj_t L_10_1;
		bli_acquire_mpart_l2r( BLIS_SUBPART1, idx2, bs2, &L_10, &L_10_1 );
		obj_t X_0_1;
		bli_acquire_mpart_t2b( BLIS_SUBPART1, idx2, bs2, &X_0, &X_0_1 );
		//------------------------------------//

		obj_t L_10_1T;
		bli_obj_alias_with_trans( BLIS_TRANSPOSE, L_10_1, L_10_1T);
		bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_ROW_PANELS, 
				BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
				BLIS_BUFFER_FOR_A_BLOCK,
				gemm_mr,  gemm_mr,  
				&L_10_1T, &L_10_1T_packed );
		bli_packm_blk_var2( &BLIS_ONE, &L_10_1T, &L_10_1T_packed );
		bli_gemm_ker_var2( &BLIS_ONE, &L_10_1T_packed, &X_1_packed, 
				&BLIS_ONE, &X_0_1, (gemm_t*)NULL );

		//------------------------------------//

	//****
	}

	//------------------------------------//

//****
}



    bli_obj_release_pack( &L_11_packed );
    bli_obj_release_pack( &L_10_1T_packed );
    bli_obj_release_pack( &X_1_packed );
}


void DxT_TrmmLUN( obj_t *alpha,
		  obj_t *U,
		  obj_t *X )
{
  obj_t X_1_packed;
  obj_t U_11_packed;
  obj_t U_01_1_packed;
    bli_obj_init_pack( &U_11_packed );
    bli_obj_init_pack( &U_01_1_packed );
    bli_obj_init_pack( &X_1_packed );

dim_t idx1, dimLen1, bs1;
///// Blocksize = 256
dimLen1 = bli_obj_length_after_trans( *X );
for ( idx1 = 0; idx1 < dimLen1; idx1 += bs1 ) {
	bs1 = bli_determine_blocksize_f( idx1, dimLen1, X, gemm_kc );
	dim_t idx2, dimLen2, bs2;
//****
	obj_t U_01;
	bli_acquire_mpart_tl2br( BLIS_SUBPART01, idx1, bs1, U, &U_01 );
	obj_t U_11;
	bli_acquire_mpart_tl2br( BLIS_SUBPART11, idx1, bs1, U, &U_11 );
	obj_t X_0;
	bli_acquire_mpart_t2b( BLIS_SUBPART0, idx1, bs1, X, &X_0 );
	obj_t X_1;
	bli_acquire_mpart_t2b( BLIS_SUBPART1, idx1, bs1, X, &X_1 );
	//------------------------------------//

	bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_COL_PANELS, 
			BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			BLIS_BUFFER_FOR_B_PANEL,
			gemm_mr,  gemm_nr,  
			&X_1, &X_1_packed );
	bli_packm_blk_var2( &BLIS_ONE, &X_1, &X_1_packed );
	///// Blocksize = 128
	dimLen2 = bli_obj_length_after_trans( X_0 );
	for ( idx2 = 0; idx2 < dimLen2; idx2 += bs2 ) {
		bs2 = bli_determine_blocksize_f( idx2, dimLen2, &X_0, gemm_mc );
		dim_t idx3, dimLen3, bs3;
	//****
		obj_t U_01_1;
		bli_acquire_mpart_t2b( BLIS_SUBPART1, idx2, bs2, &U_01, &U_01_1 );
		obj_t X_0_1;
		bli_acquire_mpart_t2b( BLIS_SUBPART1, idx2, bs2, &X_0, &X_0_1 );
		//------------------------------------//

		bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_ROW_PANELS, 
				BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
				BLIS_BUFFER_FOR_A_BLOCK,
				gemm_mr,  gemm_mr,  
				&U_01_1, &U_01_1_packed );
		bli_packm_blk_var2( &BLIS_ONE, &U_01_1, &U_01_1_packed );
		bli_gemm_ker_var2( &BLIS_ONE, &U_01_1_packed, &X_1_packed, 
				&BLIS_ONE, &X_0_1, (gemm_t*)NULL );

		//------------------------------------//

	//****
	}
	///// Blocksize = 128
	dimLen2 = bli_obj_length_after_trans( U_11 );
	for ( idx2 = 0; idx2 < dimLen2; idx2 += bs2 ) {
		bs2 = bli_determine_blocksize_f( idx2, dimLen2, &U_11, gemm_mc );
		dim_t idx3, dimLen3, bs3;
	//****
		obj_t U_11_1;
		bli_acquire_mpart_t2b( BLIS_SUBPART1, idx2, bs2, &U_11, &U_11_1 );
		obj_t X_1_1;
		bli_acquire_mpart_t2b( BLIS_SUBPART1, idx2, bs2, &X_1, &X_1_1 );
		//------------------------------------//

		bli_obj_set_struc( BLIS_TRIANGULAR, U_11_1 );
		bli_obj_set_uplo( BLIS_UPPER, U_11_1 );
		bli_packm_init_pack( TRUE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_ROW_PANELS, 
				BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
				BLIS_BUFFER_FOR_A_BLOCK,
				gemm_mr,  gemm_mr,  
				&U_11_1, &U_11_packed );
		bli_packm_blk_var3( &BLIS_ONE, &U_11_1, &U_11_packed );
		bli_trmm_u_ker_var2( &BLIS_ONE, &U_11_packed, &X_1_packed, 
				&BLIS_ZERO, &X_1_1, (trmm_t*)NULL );

		//------------------------------------//

	//****
	}

	//------------------------------------//

//****
}


    bli_obj_release_pack( &U_11_packed );
    bli_obj_release_pack( &U_01_1_packed );
    bli_obj_release_pack( &X_1_packed );
}

void DxT_TrmmLUT( obj_t *alpha,
		  obj_t *U,
		  obj_t *X )
{
  obj_t X_1_packed;
  obj_t U_11_packed;
  obj_t U_12_1T_packed;
    bli_obj_init_pack( &U_11_packed );
    bli_obj_init_pack( &U_12_1T_packed );
    bli_obj_init_pack( &X_1_packed );

dim_t idx1, dimLen1, bs1;
///// Blocksize = 256
dimLen1 = bli_obj_length_after_trans( *X );
for ( idx1 = 0; idx1 < dimLen1; idx1 += bs1 ) {
	bs1 = bli_determine_blocksize_b( idx1, dimLen1, X, gemm_kc );
	dim_t idx2, dimLen2, bs2;
//****
	obj_t U_11;
	bli_acquire_mpart_br2tl( BLIS_SUBPART11, idx1, bs1, U, &U_11 );
	obj_t U_12;
	bli_acquire_mpart_br2tl( BLIS_SUBPART12, idx1, bs1, U, &U_12 );
	obj_t X_1;
	bli_acquire_mpart_b2t( BLIS_SUBPART1, idx1, bs1, X, &X_1 );
	obj_t X_2;
	bli_acquire_mpart_b2t( BLIS_SUBPART2, idx1, bs1, X, &X_2 );
	//------------------------------------//

	bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_COL_PANELS, 
			BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			BLIS_BUFFER_FOR_B_PANEL,
			gemm_mr,  gemm_nr,  
			&X_1, &X_1_packed );
	bli_packm_blk_var2( &BLIS_ONE, &X_1, &X_1_packed );
	///// Blocksize = 128
	dimLen2 = bli_obj_length_after_trans( X_2 );
	for ( idx2 = 0; idx2 < dimLen2; idx2 += bs2 ) {
		bs2 = bli_determine_blocksize_f( idx2, dimLen2, &X_2, gemm_mc );
		dim_t idx3, dimLen3, bs3;
	//****
		obj_t U_12_1;
		bli_acquire_mpart_l2r( BLIS_SUBPART1, idx2, bs2, &U_12, &U_12_1 );
		obj_t X_2_1;
		bli_acquire_mpart_t2b( BLIS_SUBPART1, idx2, bs2, &X_2, &X_2_1 );
		//------------------------------------//

		obj_t U_12_1T;
		bli_obj_alias_with_trans( BLIS_TRANSPOSE, U_12_1, U_12_1T);
		bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_ROW_PANELS, 
				BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
				BLIS_BUFFER_FOR_A_BLOCK,
				gemm_mr,  gemm_mr,  
				&U_12_1T, &U_12_1T_packed );
		bli_packm_blk_var2( &BLIS_ONE, &U_12_1T, &U_12_1T_packed );
		bli_gemm_ker_var2( &BLIS_ONE, &U_12_1T_packed, &X_1_packed, 
				&BLIS_ONE, &X_2_1, (gemm_t*)NULL );

		//------------------------------------//

	//****
	}
	///// Blocksize = 128
	dimLen2 = bli_obj_width_after_trans( U_11 );
	for ( idx2 = 0; idx2 < dimLen2; idx2 += bs2 ) {
		bs2 = bli_determine_blocksize_f( idx2, dimLen2, &U_11, gemm_mc );
		dim_t idx3, dimLen3, bs3;
	//****
		obj_t U_11_1;
		bli_acquire_mpart_l2r( BLIS_SUBPART1, idx2, bs2, &U_11, &U_11_1 );
		obj_t X_1_1;
		bli_acquire_mpart_t2b( BLIS_SUBPART1, idx2, bs2, &X_1, &X_1_1 );
		//------------------------------------//

		obj_t U_11_1T;
		bli_obj_alias_with_trans( BLIS_TRANSPOSE, U_11_1, U_11_1T);
		bli_obj_set_struc( BLIS_TRIANGULAR, U_11_1T );
		bli_obj_set_uplo( BLIS_UPPER, U_11_1T );
		bli_packm_init_pack( TRUE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_ROW_PANELS, 
				BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
				BLIS_BUFFER_FOR_A_BLOCK,
				gemm_mr,  gemm_mr,  
				&U_11_1T, &U_11_packed );
		bli_packm_blk_var3( &BLIS_ONE, &U_11_1T, &U_11_packed );
		bli_trmm_l_ker_var2( &BLIS_ONE, &U_11_packed, &X_1_packed, 
				&BLIS_ZERO, &X_1_1, (trmm_t*)NULL );

		//------------------------------------//

	//****
	}

	//------------------------------------//

//****
}


  

    bli_obj_release_pack( &U_11_packed );
    bli_obj_release_pack( &U_12_1T_packed );
    bli_obj_release_pack( &X_1_packed );
}

void DxT_TrmmRLN( obj_t *alpha,
		  obj_t *L,
		  obj_t *X )
{
  obj_t L_11_packed;
  obj_t L_10_1T_packed;
  obj_t X_1_packed;
    bli_obj_init_pack( &L_11_packed );
    bli_obj_init_pack( &L_10_1T_packed );
    bli_obj_init_pack( &X_1_packed );

  bli_scalm(alpha, X);
dim_t idx1, dimLen1, bs1;
bli_obj_induce_trans( *X );
///// Blocksize = 256
dimLen1 = bli_obj_length_after_trans( *X );
for ( idx1 = 0; idx1 < dimLen1; idx1 += bs1 ) {
	bs1 = bli_determine_blocksize_f( idx1, dimLen1, X, gemm_kc );
	dim_t idx2, dimLen2, bs2;
//****
	obj_t L_10;
	bli_acquire_mpart_tl2br( BLIS_SUBPART10, idx1, bs1, L, &L_10 );
	obj_t L_11;
	bli_acquire_mpart_tl2br( BLIS_SUBPART11, idx1, bs1, L, &L_11 );
	obj_t X_0;
	bli_acquire_mpart_t2b( BLIS_SUBPART0, idx1, bs1, X, &X_0 );
	obj_t X_1;
	bli_acquire_mpart_t2b( BLIS_SUBPART1, idx1, bs1, X, &X_1 );
	//------------------------------------//

	bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_COL_PANELS, 
			BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			BLIS_BUFFER_FOR_B_PANEL,
			gemm_mr,  gemm_nr,  
			&X_1, &X_1_packed );
	bli_packm_blk_var2( &BLIS_ONE, &X_1, &X_1_packed );
	///// Blocksize = 128
	dimLen2 = bli_obj_width_after_trans( L_11 );
	for ( idx2 = 0; idx2 < dimLen2; idx2 += bs2 ) {
		bs2 = bli_determine_blocksize_f( idx2, dimLen2, &L_11, gemm_mc );
		dim_t idx3, dimLen3, bs3;
	//****
		obj_t L_11_1;
		bli_acquire_mpart_l2r( BLIS_SUBPART1, idx2, bs2, &L_11, &L_11_1 );
		obj_t X_1_1;
		bli_acquire_mpart_t2b( BLIS_SUBPART1, idx2, bs2, &X_1, &X_1_1 );
		//------------------------------------//

		obj_t L_11_1T;
		bli_obj_alias_with_trans( BLIS_TRANSPOSE, L_11_1, L_11_1T);
		bli_obj_set_struc( BLIS_TRIANGULAR, L_11_1T );
		bli_obj_set_uplo( BLIS_LOWER, L_11_1T );
		bli_packm_init_pack( TRUE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_ROW_PANELS, 
				BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
				BLIS_BUFFER_FOR_A_BLOCK,
				gemm_mr,  gemm_mr,  
				&L_11_1T, &L_11_packed );
		bli_packm_blk_var3( &BLIS_ONE, &L_11_1T, &L_11_packed );
		bli_trmm_u_ker_var2( &BLIS_ONE, &L_11_packed, &X_1_packed, 
				&BLIS_ZERO, &X_1_1, (trmm_t*)NULL );

		//------------------------------------//

	//****
	}
	///// Blocksize = 128
	dimLen2 = bli_obj_length_after_trans( X_0 );
	for ( idx2 = 0; idx2 < dimLen2; idx2 += bs2 ) {
		bs2 = bli_determine_blocksize_f( idx2, dimLen2, &X_0, gemm_mc );
		dim_t idx3, dimLen3, bs3;
	//****
		obj_t L_10_1;
		bli_acquire_mpart_l2r( BLIS_SUBPART1, idx2, bs2, &L_10, &L_10_1 );
		obj_t X_0_1;
		bli_acquire_mpart_t2b( BLIS_SUBPART1, idx2, bs2, &X_0, &X_0_1 );
		//------------------------------------//

		obj_t L_10_1T;
		bli_obj_alias_with_trans( BLIS_TRANSPOSE, L_10_1, L_10_1T);
		bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_ROW_PANELS, 
				BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
				BLIS_BUFFER_FOR_A_BLOCK,
				gemm_mr,  gemm_mr,  
				&L_10_1T, &L_10_1T_packed );
		bli_packm_blk_var2( &BLIS_ONE, &L_10_1T, &L_10_1T_packed );
		bli_gemm_ker_var2( &BLIS_ONE, &L_10_1T_packed, &X_1_packed, 
				&BLIS_ONE, &X_0_1, (gemm_t*)NULL );

		//------------------------------------//

	//****
	}

	//------------------------------------//

//****
}
bli_obj_induce_trans( *X );


    bli_obj_release_pack( &L_11_packed );
    bli_obj_release_pack( &L_10_1T_packed );
    bli_obj_release_pack( &X_1_packed );
}

void DxT_TrmmRLT( obj_t *alpha,
		  obj_t *L,
		  obj_t *X )
{
  obj_t L_11_packed;
  obj_t L_21_1_packed;
  obj_t X_1_packed;
    bli_obj_init_pack( &L_11_packed );
    bli_obj_init_pack( &L_21_1_packed );
    bli_obj_init_pack( &X_1_packed );  

  bli_scalm(alpha, X);
dim_t idx1, dimLen1, bs1;
bli_obj_induce_trans( *X );
///// Blocksize = 256
dimLen1 = bli_obj_length_after_trans( *X );
for ( idx1 = 0; idx1 < dimLen1; idx1 += bs1 ) {
	bs1 = bli_determine_blocksize_b( idx1, dimLen1, X, gemm_kc );
	dim_t idx2, dimLen2, bs2;
//****
	obj_t L_11;
	bli_acquire_mpart_br2tl( BLIS_SUBPART11, idx1, bs1, L, &L_11 );
	obj_t L_21;
	bli_acquire_mpart_br2tl( BLIS_SUBPART21, idx1, bs1, L, &L_21 );
	obj_t X_1;
	bli_acquire_mpart_b2t( BLIS_SUBPART1, idx1, bs1, X, &X_1 );
	obj_t X_2;
	bli_acquire_mpart_b2t( BLIS_SUBPART2, idx1, bs1, X, &X_2 );
	//------------------------------------//

	bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_COL_PANELS, 
			BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			BLIS_BUFFER_FOR_B_PANEL,
			gemm_mr,  gemm_nr,  
			&X_1, &X_1_packed );
	bli_packm_blk_var2( &BLIS_ONE, &X_1, &X_1_packed );
	///// Blocksize = 128
	dimLen2 = bli_obj_length_after_trans( X_2 );
	for ( idx2 = 0; idx2 < dimLen2; idx2 += bs2 ) {
		bs2 = bli_determine_blocksize_f( idx2, dimLen2, &X_2, gemm_mc );
		dim_t idx3, dimLen3, bs3;
	//****
		obj_t L_21_1;
		bli_acquire_mpart_t2b( BLIS_SUBPART1, idx2, bs2, &L_21, &L_21_1 );
		obj_t X_2_1;
		bli_acquire_mpart_t2b( BLIS_SUBPART1, idx2, bs2, &X_2, &X_2_1 );
		//------------------------------------//

		bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_ROW_PANELS, 
				BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
				BLIS_BUFFER_FOR_A_BLOCK,
				gemm_mr,  gemm_mr,  
				&L_21_1, &L_21_1_packed );
		bli_packm_blk_var2( &BLIS_ONE, &L_21_1, &L_21_1_packed );
		bli_gemm_ker_var2( &BLIS_ONE, &L_21_1_packed, &X_1_packed, 
				&BLIS_ONE, &X_2_1, (gemm_t*)NULL );

		//------------------------------------//

	//****
	}
	///// Blocksize = 128
	dimLen2 = bli_obj_length_after_trans( L_11 );
	for ( idx2 = 0; idx2 < dimLen2; idx2 += bs2 ) {
		bs2 = bli_determine_blocksize_f( idx2, dimLen2, &L_11, gemm_mc );
		dim_t idx3, dimLen3, bs3;
	//****
		obj_t L_11_1;
		bli_acquire_mpart_t2b( BLIS_SUBPART1, idx2, bs2, &L_11, &L_11_1 );
		obj_t X_1_1;
		bli_acquire_mpart_t2b( BLIS_SUBPART1, idx2, bs2, &X_1, &X_1_1 );
		//------------------------------------//

		bli_obj_set_struc( BLIS_TRIANGULAR, L_11_1 );
		bli_obj_set_uplo( BLIS_LOWER, L_11_1 );
		bli_packm_init_pack( TRUE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_ROW_PANELS, 
				BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
				BLIS_BUFFER_FOR_A_BLOCK,
				gemm_mr,  gemm_mr,  
				&L_11_1, &L_11_packed );
		bli_packm_blk_var3( &BLIS_ONE, &L_11_1, &L_11_packed );
		bli_trmm_l_ker_var2( &BLIS_ONE, &L_11_packed, &X_1_packed, 
				&BLIS_ZERO, &X_1_1, (trmm_t*)NULL );

		//------------------------------------//

	//****
	}

	//------------------------------------//

//****
}
bli_obj_induce_trans( *X );



    bli_obj_release_pack( &L_11_packed );
    bli_obj_release_pack( &L_21_1_packed );
    bli_obj_release_pack( &X_1_packed );
}

void DxT_TrmmRUN( obj_t *alpha,
		  obj_t *U,
		  obj_t *X )
{
  obj_t X_1_packed;
  obj_t U_11_packed;
  obj_t U_12_1T_packed;
    bli_obj_init_pack( &U_11_packed );
    bli_obj_init_pack( &U_12_1T_packed );
    bli_obj_init_pack( &X_1_packed );
  
  bli_scalm(alpha, X);
dim_t idx1, dimLen1, bs1;
bli_obj_induce_trans( *X );
///// Blocksize = 256
dimLen1 = bli_obj_length_after_trans( *X );
for ( idx1 = 0; idx1 < dimLen1; idx1 += bs1 ) {
	bs1 = bli_determine_blocksize_b( idx1, dimLen1, X, gemm_kc );
	dim_t idx2, dimLen2, bs2;
//****
	obj_t U_11;
	bli_acquire_mpart_br2tl( BLIS_SUBPART11, idx1, bs1, U, &U_11 );
	obj_t U_12;
	bli_acquire_mpart_br2tl( BLIS_SUBPART12, idx1, bs1, U, &U_12 );
	obj_t X_1;
	bli_acquire_mpart_b2t( BLIS_SUBPART1, idx1, bs1, X, &X_1 );
	obj_t X_2;
	bli_acquire_mpart_b2t( BLIS_SUBPART2, idx1, bs1, X, &X_2 );
	//------------------------------------//

	bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_COL_PANELS, 
			BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			BLIS_BUFFER_FOR_B_PANEL,
			gemm_mr,  gemm_nr,  
			&X_1, &X_1_packed );
	bli_packm_blk_var2( &BLIS_ONE, &X_1, &X_1_packed );
	///// Blocksize = 128
	dimLen2 = bli_obj_width_after_trans( U_11 );
	for ( idx2 = 0; idx2 < dimLen2; idx2 += bs2 ) {
		bs2 = bli_determine_blocksize_f( idx2, dimLen2, &U_11, gemm_mc );
		dim_t idx3, dimLen3, bs3;
	//****
		obj_t U_11_1;
		bli_acquire_mpart_l2r( BLIS_SUBPART1, idx2, bs2, &U_11, &U_11_1 );
		obj_t X_1_1;
		bli_acquire_mpart_t2b( BLIS_SUBPART1, idx2, bs2, &X_1, &X_1_1 );
		//------------------------------------//

		obj_t U_11_1T;
		bli_obj_alias_with_trans( BLIS_TRANSPOSE, U_11_1, U_11_1T);
		bli_obj_set_struc( BLIS_TRIANGULAR, U_11_1T );
		bli_obj_set_uplo( BLIS_UPPER, U_11_1T );
		bli_packm_init_pack( TRUE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_ROW_PANELS, 
				BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
				BLIS_BUFFER_FOR_A_BLOCK,
				gemm_mr,  gemm_mr,  
				&U_11_1T, &U_11_packed );
		bli_packm_blk_var3( &BLIS_ONE, &U_11_1T, &U_11_packed );
		bli_trmm_l_ker_var2( &BLIS_ONE, &U_11_packed, &X_1_packed, 
				&BLIS_ZERO, &X_1_1, (trmm_t*)NULL );

		//------------------------------------//

	//****
	}
	///// Blocksize = 128
	dimLen2 = bli_obj_length_after_trans( X_2 );
	for ( idx2 = 0; idx2 < dimLen2; idx2 += bs2 ) {
		bs2 = bli_determine_blocksize_f( idx2, dimLen2, &X_2, gemm_mc );
		dim_t idx3, dimLen3, bs3;
	//****
		obj_t U_12_1;
		bli_acquire_mpart_l2r( BLIS_SUBPART1, idx2, bs2, &U_12, &U_12_1 );
		obj_t X_2_1;
		bli_acquire_mpart_t2b( BLIS_SUBPART1, idx2, bs2, &X_2, &X_2_1 );
		//------------------------------------//

		obj_t U_12_1T;
		bli_obj_alias_with_trans( BLIS_TRANSPOSE, U_12_1, U_12_1T);
		bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_ROW_PANELS, 
				BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
				BLIS_BUFFER_FOR_A_BLOCK,
				gemm_mr,  gemm_mr,  
				&U_12_1T, &U_12_1T_packed );
		bli_packm_blk_var2( &BLIS_ONE, &U_12_1T, &U_12_1T_packed );
		bli_gemm_ker_var2( &BLIS_ONE, &U_12_1T_packed, &X_1_packed, 
				&BLIS_ONE, &X_2_1, (gemm_t*)NULL );

		//------------------------------------//

	//****
	}

	//------------------------------------//

//****
}
bli_obj_induce_trans( *X );


    bli_obj_release_pack( &U_11_packed );
    bli_obj_release_pack( &U_12_1T_packed );
    bli_obj_release_pack( &X_1_packed );
}


void DxT_TrmmRUT( obj_t *alpha,
		  obj_t *U,
		  obj_t *X )
{
  obj_t X_1_packed;
  obj_t U_11_packed;
  obj_t U_01_1_packed;
    bli_obj_init_pack( &U_11_packed );
    bli_obj_init_pack( &U_01_1_packed );
    bli_obj_init_pack( &X_1_packed );  

  bli_scalm(alpha, X);
dim_t idx1, dimLen1, bs1;
bli_obj_induce_trans( *X );
///// Blocksize = 256
dimLen1 = bli_obj_length_after_trans( *X );
for ( idx1 = 0; idx1 < dimLen1; idx1 += bs1 ) {
	bs1 = bli_determine_blocksize_f( idx1, dimLen1, X, gemm_kc );
	dim_t idx2, dimLen2, bs2;
//****
	obj_t U_01;
	bli_acquire_mpart_tl2br( BLIS_SUBPART01, idx1, bs1, U, &U_01 );
	obj_t U_11;
	bli_acquire_mpart_tl2br( BLIS_SUBPART11, idx1, bs1, U, &U_11 );
	obj_t X_0;
	bli_acquire_mpart_t2b( BLIS_SUBPART0, idx1, bs1, X, &X_0 );
	obj_t X_1;
	bli_acquire_mpart_t2b( BLIS_SUBPART1, idx1, bs1, X, &X_1 );
	//------------------------------------//

	bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_COL_PANELS, 
			BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			BLIS_BUFFER_FOR_B_PANEL,
			gemm_mr,  gemm_nr,  
			&X_1, &X_1_packed );
	bli_packm_blk_var2( &BLIS_ONE, &X_1, &X_1_packed );
	///// Blocksize = 128
	dimLen2 = bli_obj_length_after_trans( X_0 );
	for ( idx2 = 0; idx2 < dimLen2; idx2 += bs2 ) {
		bs2 = bli_determine_blocksize_f( idx2, dimLen2, &X_0, gemm_mc );
		dim_t idx3, dimLen3, bs3;
	//****
		obj_t U_01_1;
		bli_acquire_mpart_t2b( BLIS_SUBPART1, idx2, bs2, &U_01, &U_01_1 );
		obj_t X_0_1;
		bli_acquire_mpart_t2b( BLIS_SUBPART1, idx2, bs2, &X_0, &X_0_1 );
		//------------------------------------//

		bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_ROW_PANELS, 
				BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
				BLIS_BUFFER_FOR_A_BLOCK,
				gemm_mr,  gemm_mr,  
				&U_01_1, &U_01_1_packed );
		bli_packm_blk_var2( &BLIS_ONE, &U_01_1, &U_01_1_packed );
		bli_gemm_ker_var2( &BLIS_ONE, &U_01_1_packed, &X_1_packed, 
				&BLIS_ONE, &X_0_1, (gemm_t*)NULL );

		//------------------------------------//

	//****
	}
	///// Blocksize = 128
	dimLen2 = bli_obj_length_after_trans( U_11 );
	for ( idx2 = 0; idx2 < dimLen2; idx2 += bs2 ) {
		bs2 = bli_determine_blocksize_f( idx2, dimLen2, &U_11, gemm_mc );
		dim_t idx3, dimLen3, bs3;
	//****
		obj_t U_11_1;
		bli_acquire_mpart_t2b( BLIS_SUBPART1, idx2, bs2, &U_11, &U_11_1 );
		obj_t X_1_1;
		bli_acquire_mpart_t2b( BLIS_SUBPART1, idx2, bs2, &X_1, &X_1_1 );
		//------------------------------------//

		bli_obj_set_struc( BLIS_TRIANGULAR, U_11_1 );
		bli_obj_set_uplo( BLIS_UPPER, U_11_1 );
		bli_packm_init_pack( TRUE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_ROW_PANELS, 
				BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
				BLIS_BUFFER_FOR_A_BLOCK,
				gemm_mr,  gemm_mr,  
				&U_11_1, &U_11_packed );
		bli_packm_blk_var3( &BLIS_ONE, &U_11_1, &U_11_packed );
		bli_trmm_u_ker_var2( &BLIS_ONE, &U_11_packed, &X_1_packed, 
				&BLIS_ZERO, &X_1_1, (trmm_t*)NULL );

		//------------------------------------//

	//****
	}

	//------------------------------------//

//****
}
bli_obj_induce_trans( *X );


  bli_obj_release_pack( &U_11_packed );
  bli_obj_release_pack( &U_01_1_packed );
  bli_obj_release_pack( &X_1_packed );
}



int main( int argc, char** argv )
{
  obj_t a, c1, c2;
  obj_t c_save;
  obj_t alpha, beta, negOne, normVal;
  dim_t m, n;
  dim_t p;
  dim_t p_begin, p_end, p_inc;
  int   m_input, n_input;
  num_t dt_a, dt_b, dt_c;
  num_t dt_alpha, dt_beta;
  int   r, n_repeats;

  double dtime;
  double dtime_save;
  double gflops;

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
      bli_obj_create( dt_c, m, n, 0, 0, &c1 );
      bli_obj_create( dt_c, m, n, 0, 0, &c2 );
      bli_obj_create( dt_c, m, n, 0, 0, &c_save );


      bli_randm( &a );
      bli_randm( &c1 );

      bli_setsc(  (1.0/1.0), 0.0, &alpha );
      bli_setsc(  (1.0/1.0), 0.0, &beta );
      bli_setsc(  (1.0/1.0), 0.0, &normVal );
      bli_copym( &c1, &c_save );
	
      dtime_save = 1.0e9;

      for ( r = 0; r < n_repeats; ++r )
	{
	  bli_copym( &c_save, &c1 );
	  bli_copym( &c_save, &c2 );

	  dtime = bli_clock();

	  if (left) {
	    if (lower) {
	      if (trans) {
		DxT_TrmmLLT(&alpha,
			    &a,
			    &c1 );
	      }
	      else { //normal
		DxT_TrmmLLN(&alpha,
			    &a,
			    &c1 );
	      }
	    }
	    else { //upper
	      if (trans) {
		DxT_TrmmLUT(&alpha,
			    &a,
			    &c1 );
	      }
	      else {//normal 
		DxT_TrmmLUN(&alpha,
			    &a,
			    &c1 );
	      }
	    }
	  }
	  else { //right
	    if (lower) {
	      if (trans) {
		DxT_TrmmRLT(&alpha,
			    &a,
			    &c1 );
	      }
	      else { //normal
		DxT_TrmmRLN(&alpha,
			    &a,
			    &c1 );
	      }

	    }
	    else { //upper
	      if (trans) {
		DxT_TrmmRUT(&alpha,
			    &a,
			    &c1 );
	      }
	      else { //normal
		DxT_TrmmRUN(&alpha,
			    &a,
			    &c1 );
	      }
	    }
	  }

	  if (trans)
	    bli_obj_set_conjtrans( BLIS_TRANSPOSE, a);

	  bli_obj_set_struc( BLIS_TRIANGULAR, a );
	  if (lower) {
	    bli_obj_set_uplo( BLIS_LOWER, a );
	  }
	  else {
	    bli_obj_set_uplo( BLIS_UPPER, a );
	  }

	  bli_trmm( left ? BLIS_LEFT : BLIS_RIGHT,
		    &alpha,
		    &a,
		    &c2 );

	  if (trans)
	    bli_obj_toggle_trans( a );
	  bli_obj_set_struc( BLIS_GENERAL, a );
	  bli_obj_set_uplo( BLIS_DENSE, a );
			
	  bli_axpym( &negOne, &c1, &c2 );
			
	  bli_fnormm( &c2, &normVal );

	  dtime_save = bli_clock_min_diff( dtime_save, dtime );
	}

      if ( left )
	gflops = ( 1.0 * m * m * n ) / ( dtime_save * 1.0e9 );
      else
	gflops = ( 1.0 * m * n * n ) / ( dtime_save * 1.0e9 );

      printf( "data_trmm_blis" );
      printf( "( %2ld, 1:4 ) = [ %4lu %4lu  %10.3e  %6.3f ];\n",
	      (p - p_begin + 1)/p_inc + 1, m, n, dtime_save, gflops );

      bli_printm( "NORM", &normVal, "%4.1f", "" );
      

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

