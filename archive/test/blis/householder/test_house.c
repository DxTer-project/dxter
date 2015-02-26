#include <unistd.h>
#include "blis.h"
#include "FLAME.h"

//#define TESTLIB

extern blksz_t *gemm_mc;
extern blksz_t *gemm_kc;
extern blksz_t *gemm_nc;
extern blksz_t *gemm_mr;
extern blksz_t *gemm_kr;
extern blksz_t *gemm_nr;

extern fla_apqut_t* fla_apqut_cntl_leaf;
extern fla_trmm_t*  fla_trmm_cntl_blas;
extern fla_gemm_t*  fla_gemm_cntl_blas;
extern fla_trsm_t*  fla_trsm_cntl_blas;
extern fla_copyt_t* fla_copyt_cntl_blas;
extern fla_axpyt_t* fla_axpyt_cntl_blas;

fla_qrut_t *qrut_cntl;

FLA_Error wrap_blis_obj( obj_t *blis,
			 FLA_Obj *fla )
{
  FLA_Error err;

  err = FLA_Obj_create_without_buffer( FLA_DOUBLE, 
				       bli_obj_length(*blis),
				       bli_obj_width(*blis),
				       fla);

  if (err != FLA_SUCCESS) {
    return err;
  }

  err = FLA_Obj_attach_buffer( bli_obj_buffer_at_off( *blis ),
			       bli_obj_row_stride( *blis ),
			       bli_obj_col_stride( *blis ),
			       fla );
  return err;
}

void wrap_fla_obj( FLA_Obj *fla,
		   obj_t *blis )
{
  bli_obj_create_with_attached_buffer( BLIS_DOUBLE,
				       FLA_Obj_length(*fla),
				       FLA_Obj_width(*fla),
				       FLA_Obj_buffer_at_view(*fla),
				       FLA_Obj_row_stride(*fla),
				       FLA_Obj_col_stride(*fla),
				       blis );
}

void blis_obj_merge_vert( obj_t *top,
			  obj_t *bott,
			  obj_t *merged )
{
  if (bli_obj_width(*top) != bli_obj_width(*bott)) {
    printf("bad objects in merge_vert 1\n");
    return;
  }
  if (bli_obj_row_stride(*top) != bli_obj_row_stride(*bott)) {
    printf("bad objects in merge_vert 2\n");
    return;
  }
  if (bli_obj_col_stride(*top) != bli_obj_col_stride(*bott)){
    printf("bad objects in merge_vert 3\n");
    return;
  }
  if (bli_obj_width(*bott) != 0 && bli_obj_length(*bott) != 0 &&
      (((double*)bli_obj_buffer_at_off(*top) + (bli_obj_length(*top) * bli_obj_col_stride(*top)))
       != (double*)bli_obj_buffer_at_off(*bott)))
    {
      printf("top %p\n",(double*)bli_obj_buffer_at_off(*top));
      printf("bott %p\n",(double*)bli_obj_buffer_at_off(*bott));
      printf("top+offset %p\n",((double*)bli_obj_buffer_at_off(*top) + (bli_obj_length(*top) * bli_obj_col_stride(*top))));
      printf("length %lu\n", bli_obj_length(*top));
      printf("width %lu\n", bli_obj_width(*top));
      printf("bad objects in merge_vert 4\n");
      return;
    }
  /* 
     printf("top length %lu\n", bli_obj_length(*top));
     printf("top width %lu\n", bli_obj_width(*top));
     printf("bott length %lu\n", bli_obj_length(*bott));
     printf("bott width %lu\n", bli_obj_width(*bott));
  */
  bli_obj_create_with_attached_buffer( BLIS_DOUBLE,
				       bli_obj_length(*top) + bli_obj_length(*bott),
				       bli_obj_width(*top),
				       bli_obj_buffer_at_off(*top),
				       bli_obj_row_stride(*top),
				       bli_obj_col_stride(*top),
				       merged );
}


void DxT_AppHouse( FLA_Obj Ain,
		   FLA_Obj Tin,
		   FLA_Obj Win,
		   FLA_Obj Bin,
		   blksz_t *bs_obj )
{
  obj_t U, T, W, B;

  wrap_fla_obj( &Ain, &U);
  wrap_fla_obj( &Tin, &T );
  wrap_fla_obj( &Win, &W );
  wrap_fla_obj( &Bin, &B );

  obj_t packed_A_blk;
  bli_obj_init_pack( &packed_A_blk );
  obj_t packed_B_pan;
  bli_obj_init_pack( &packed_B_pan );

  if (bli_blksz_for_obj(&U,bs_obj) <= bli_blksz_for_obj(&U,gemm_kc)) {
    dim_t idx1, dimLen1, bs1;
    dimLen1 = min( bli_obj_length_after_trans( U ), bli_obj_width_after_trans( U ) );
    for ( idx1 = 0; idx1 < dimLen1; idx1 += bs1 ) {
      bs1 = bli_determine_blocksize_f( idx1, dimLen1, &U, bs_obj );
      dim_t idx2, dimLen2, bs2;
      //****
      obj_t B_1;
      bli_acquire_mpart_t2b( BLIS_SUBPART1, idx1, bs1, &B, &B_1 );
      obj_t B_2;
      bli_acquire_mpart_t2b( BLIS_SUBPART2, idx1, bs1, &B, &B_2 );
      obj_t U_11;
      bli_acquire_mpart_tl2br( BLIS_SUBPART11, idx1, bs1, &U, &U_11 );
      obj_t U_21;
      bli_acquire_mpart_tl2br( BLIS_SUBPART21, idx1, bs1, &U, &U_21 );
      obj_t T_1;
      bli_acquire_mpart_l2r( BLIS_SUBPART1, idx1, bs1, &T, &T_1 );
      //------------------------------------//

      obj_t T_1TL, T_1TLtmp;
      bli_acquire_mpart_l2r( BLIS_SUBPART1, 0, bli_obj_width( T_1 ), &T_1, &T_1TLtmp );
      bli_acquire_mpart_t2b( BLIS_SUBPART1, 0, bli_obj_length( B_1 ), &T_1TLtmp, &T_1TL );
      obj_t WTL, WTLtmp;
      bli_acquire_mpart_l2r( BLIS_SUBPART1, 0, bli_obj_width( B_1 ), &W, &WTLtmp );
      bli_acquire_mpart_t2b( BLIS_SUBPART1, 0, bli_obj_length( B_1 ), &WTLtmp, &WTL );
      dimLen2 = bli_obj_length_after_trans( WTL );
      for ( idx2 = 0; idx2 < dimLen2; idx2 += bs2 ) {
	bs2 = bli_determine_blocksize_f( idx2, dimLen2, &WTL, gemm_kc );
	dim_t idx3, dimLen3, bs3;
	//****
	obj_t U_11_10;
	bli_acquire_mpart_tl2br( BLIS_SUBPART10, idx2, bs2, &U_11, &U_11_10 );
	obj_t U_11_11;
	bli_acquire_mpart_tl2br( BLIS_SUBPART11, idx2, bs2, &U_11, &U_11_11 );
	obj_t B_1_1;
	bli_acquire_mpart_t2b( BLIS_SUBPART1, idx2, bs2, &B_1, &B_1_1 );
	obj_t WTL_0;
	bli_acquire_mpart_t2b( BLIS_SUBPART0, idx2, bs2, &WTL, &WTL_0 );
	obj_t WTL_1;
	bli_acquire_mpart_t2b( BLIS_SUBPART1, idx2, bs2, &WTL, &WTL_1 );
	//------------------------------------//

	bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_COL_PANELS, 
			     BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			     BLIS_BUFFER_FOR_B_PANEL,
			     gemm_mr, gemm_nr, 
			     &B_1_1, &packed_B_pan );
	bli_packm_blk_var2( &BLIS_ONE, &B_1_1, &packed_B_pan );
	dimLen3 = bli_obj_width_after_trans( U_11_11 );
	for ( idx3 = 0; idx3 < dimLen3; idx3 += bs3 ) {
	  bs3 = bli_determine_blocksize_f( idx3, dimLen3, &U_11_11, gemm_mc );
	  dim_t idx4, dimLen4, bs4;
	  //****
	  obj_t U_11_11_1;
	  bli_acquire_mpart_l2r( BLIS_SUBPART1, idx3, bs3, &U_11_11, &U_11_11_1 );
	  obj_t WTL_1_1;
	  bli_acquire_mpart_t2b( BLIS_SUBPART1, idx3, bs3, &WTL_1, &WTL_1_1 );
	  //------------------------------------//

	  obj_t U_11_11_1H;
	  bli_obj_alias_with_trans( BLIS_CONJ_TRANSPOSE, U_11_11_1, U_11_11_1H);
	  bli_obj_set_struc( BLIS_TRIANGULAR, U_11_11_1H );
	  bli_obj_set_uplo( BLIS_LOWER, U_11_11_1H );
	  bli_obj_set_diag( BLIS_UNIT_DIAG, U_11_11_1H );
	  bli_packm_init_pack( TRUE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_ROW_PANELS, 
			       BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			       BLIS_BUFFER_FOR_A_BLOCK,
			       gemm_mr, gemm_mr, 
			       &U_11_11_1H, &packed_A_blk );
	  bli_packm_blk_var3( &BLIS_ONE, &U_11_11_1H, &packed_A_blk );
	  bli_trmm_u_ker_var2( &BLIS_ONE, &packed_A_blk, &packed_B_pan, 
			       &BLIS_ZERO, &WTL_1_1, (trmm_t*)NULL );

	  //------------------------------------//

	  //****
	}
	dimLen3 = bli_obj_length_after_trans( WTL_0 );
	for ( idx3 = 0; idx3 < dimLen3; idx3 += bs3 ) {
	  bs3 = bli_determine_blocksize_f( idx3, dimLen3, &WTL_0, gemm_mc );
	  dim_t idx4, dimLen4, bs4;
	  //****
	  obj_t U_11_10_1;
	  bli_acquire_mpart_l2r( BLIS_SUBPART1, idx3, bs3, &U_11_10, &U_11_10_1 );
	  obj_t WTL_0_1;
	  bli_acquire_mpart_t2b( BLIS_SUBPART1, idx3, bs3, &WTL_0, &WTL_0_1 );
	  //------------------------------------//

	  obj_t U_11_10_1H;
	  bli_obj_alias_with_trans( BLIS_CONJ_TRANSPOSE, U_11_10_1, U_11_10_1H);
	  bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_ROW_PANELS, 
			       BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			       BLIS_BUFFER_FOR_A_BLOCK,
			       gemm_mr, gemm_mr, 
			       &U_11_10_1H, &packed_A_blk );
	  bli_packm_blk_var2( &BLIS_ONE, &U_11_10_1H, &packed_A_blk );
	  bli_gemm_ker_var2( &BLIS_ONE, &packed_A_blk, &packed_B_pan, 
			     &BLIS_ZERO, &WTL_0_1, (gemm_t*)NULL );

	  //------------------------------------//

	  //****
	}

	//------------------------------------//

	//****
      }
      dimLen2 = bli_obj_length_after_trans( U_21 );
      for ( idx2 = 0; idx2 < dimLen2; idx2 += bs2 ) {
	bs2 = bli_determine_blocksize_b( idx2, dimLen2, &U_21, gemm_kc );
	dim_t idx3, dimLen3, bs3;
	//****
	obj_t U_21_1;
	bli_acquire_mpart_b2t( BLIS_SUBPART1, idx2, bs2, &U_21, &U_21_1 );
	obj_t B_2_1;
	bli_acquire_mpart_b2t( BLIS_SUBPART1, idx2, bs2, &B_2, &B_2_1 );
	//------------------------------------//

	bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_COL_PANELS, 
			     BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			     BLIS_BUFFER_FOR_B_PANEL,
			     gemm_kr, gemm_nr, 
			     &B_2_1, &packed_B_pan );
	bli_packm_blk_var2( &BLIS_ONE, &B_2_1, &packed_B_pan );
	dimLen3 = bli_obj_length_after_trans( WTL );
	for ( idx3 = 0; idx3 < dimLen3; idx3 += bs3 ) {
	  bs3 = bli_determine_blocksize_f( idx3, dimLen3, &WTL, gemm_mc );
	  dim_t idx4, dimLen4, bs4;
	  //****
	  obj_t U_21_1_1;
	  bli_acquire_mpart_l2r( BLIS_SUBPART1, idx3, bs3, &U_21_1, &U_21_1_1 );
	  obj_t WTL_1;
	  bli_acquire_mpart_t2b( BLIS_SUBPART1, idx3, bs3, &WTL, &WTL_1 );
	  //------------------------------------//

	  obj_t U_21_1_1H;
	  bli_obj_alias_with_trans( BLIS_CONJ_TRANSPOSE, U_21_1_1, U_21_1_1H);
	  bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_ROW_PANELS, 
			       BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			       BLIS_BUFFER_FOR_A_BLOCK,
			       gemm_mr, gemm_kr, 
			       &U_21_1_1H, &packed_A_blk );
	  bli_packm_blk_var2( &BLIS_ONE, &U_21_1_1H, &packed_A_blk );
	  bli_gemm_ker_var2( &BLIS_ONE, &packed_A_blk, &packed_B_pan, 
			     &BLIS_ONE, &WTL_1, (gemm_t*)NULL );

	  //------------------------------------//

	  //****
	}

	//------------------------------------//

	//****
      }
      dimLen2 = bli_obj_length_after_trans( WTL );
      for ( idx2 = 0; idx2 < dimLen2; idx2 += bs2 ) {
	bs2 = bli_determine_blocksize_f( idx2, dimLen2, &WTL, gemm_kc );
	dim_t idx3, dimLen3, bs3;
	//****
	obj_t T_1TL_11;
	bli_acquire_mpart_tl2br( BLIS_SUBPART11, idx2, bs2, &T_1TL, &T_1TL_11 );
	obj_t T_1TL_12;
	bli_acquire_mpart_tl2br( BLIS_SUBPART12, idx2, bs2, &T_1TL, &T_1TL_12 );
	obj_t WTL_1;
	bli_acquire_mpart_t2b( BLIS_SUBPART1, idx2, bs2, &WTL, &WTL_1 );
	obj_t WTL_2;
	bli_acquire_mpart_t2b( BLIS_SUBPART2, idx2, bs2, &WTL, &WTL_2 );
	//------------------------------------//

	bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_COL_PANELS, 
			     BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			     BLIS_BUFFER_FOR_B_PANEL,
			     gemm_mr, gemm_nr, 
			     &WTL_1, &packed_B_pan );
	bli_packm_blk_var2( &BLIS_ONE, &WTL_1, &packed_B_pan );
	dimLen3 = bli_obj_width_after_trans( T_1TL_11 );
	for ( idx3 = 0; idx3 < dimLen3; idx3 += bs3 ) {
	  bs3 = bli_determine_blocksize_f( idx3, dimLen3, &T_1TL_11, gemm_mc );
	  dim_t idx4, dimLen4, bs4;
	  //****
	  obj_t T_1TL_11_1;
	  bli_acquire_mpart_l2r( BLIS_SUBPART1, idx3, bs3, &T_1TL_11, &T_1TL_11_1 );
	  obj_t WTL_1_1;
	  bli_acquire_mpart_t2b( BLIS_SUBPART1, idx3, bs3, &WTL_1, &WTL_1_1 );
	  //------------------------------------//

	  obj_t T_1TL_11_1H;
	  bli_obj_alias_with_trans( BLIS_CONJ_TRANSPOSE, T_1TL_11_1, T_1TL_11_1H);
	  bli_obj_set_struc( BLIS_TRIANGULAR, T_1TL_11_1H );
	  bli_obj_set_uplo( BLIS_UPPER, T_1TL_11_1H );
	  bli_packm_init_pack( TRUE, BLIS_INVERT_DIAG, BLIS_PACKED_ROW_PANELS, 
			       BLIS_PACK_REV_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			       BLIS_BUFFER_FOR_A_BLOCK,
			       gemm_mr, gemm_mr, 
			       &T_1TL_11_1H, &packed_A_blk );
	  bli_packm_blk_var3( &BLIS_ONE, &T_1TL_11_1H, &packed_A_blk );
	  bli_trsm_l_ker_var2( &BLIS_ONE, &packed_A_blk, &packed_B_pan, 
			       &BLIS_ZERO, &WTL_1_1, (trsm_t*)NULL );

	  //------------------------------------//

	  //****
	}
	dimLen3 = bli_obj_length_after_trans( WTL_2 );
	for ( idx3 = 0; idx3 < dimLen3; idx3 += bs3 ) {
	  bs3 = bli_determine_blocksize_f( idx3, dimLen3, &WTL_2, gemm_mc );
	  dim_t idx4, dimLen4, bs4;
	  //****
	  obj_t T_1TL_12_1;
	  bli_acquire_mpart_l2r( BLIS_SUBPART1, idx3, bs3, &T_1TL_12, &T_1TL_12_1 );
	  obj_t WTL_2_1;
	  bli_acquire_mpart_t2b( BLIS_SUBPART1, idx3, bs3, &WTL_2, &WTL_2_1 );
	  //------------------------------------//

	  obj_t T_1TL_12_1H;
	  bli_obj_alias_with_trans( BLIS_CONJ_TRANSPOSE, T_1TL_12_1, T_1TL_12_1H);
	  bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_ROW_PANELS, 
			       BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			       BLIS_BUFFER_FOR_A_BLOCK,
			       gemm_mr, gemm_mr, 
			       &T_1TL_12_1H, &packed_A_blk );
	  bli_packm_blk_var2( &BLIS_ONE, &T_1TL_12_1H, &packed_A_blk );
	  bli_gemm_ker_var2( &BLIS_MINUS_ONE, &packed_A_blk, &packed_B_pan, 
			     &BLIS_ONE, &WTL_2_1, (gemm_t*)NULL );

	  //------------------------------------//

	  //****
	}

	//------------------------------------//

	//****
      }
      dimLen2 = bli_obj_width_after_trans( U_21 );
      for ( idx2 = 0; idx2 < dimLen2; idx2 += bs2 ) {
	bs2 = bli_determine_blocksize_b( idx2, dimLen2, &U_21, gemm_kc );
	dim_t idx3, dimLen3, bs3;
	//****
	obj_t U_21_1;
	bli_acquire_mpart_r2l( BLIS_SUBPART1, idx2, bs2, &U_21, &U_21_1 );
	obj_t WTL_1;
	bli_acquire_mpart_b2t( BLIS_SUBPART1, idx2, bs2, &WTL, &WTL_1 );
	obj_t U_11_11;
	bli_acquire_mpart_br2tl( BLIS_SUBPART11, idx2, bs2, &U_11, &U_11_11 );
	obj_t U_11_21;
	bli_acquire_mpart_br2tl( BLIS_SUBPART21, idx2, bs2, &U_11, &U_11_21 );
	obj_t B_1_1;
	bli_acquire_mpart_b2t( BLIS_SUBPART1, idx2, bs2, &B_1, &B_1_1 );
	obj_t B_1_2;
	bli_acquire_mpart_b2t( BLIS_SUBPART2, idx2, bs2, &B_1, &B_1_2 );
	//------------------------------------//

	bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_COL_PANELS, 
			     BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			     BLIS_BUFFER_FOR_B_PANEL,
			     gemm_mr, gemm_nr, 
			     &WTL_1, &packed_B_pan );
	bli_packm_blk_var2( &BLIS_ONE, &WTL_1, &packed_B_pan );
	dimLen3 = bli_obj_length_after_trans( B_2 );
	for ( idx3 = 0; idx3 < dimLen3; idx3 += bs3 ) {
	  bs3 = bli_determine_blocksize_f( idx3, dimLen3, &B_2, gemm_mc );
	  dim_t idx4, dimLen4, bs4;
	  //****
	  obj_t U_21_1_1;
	  bli_acquire_mpart_t2b( BLIS_SUBPART1, idx3, bs3, &U_21_1, &U_21_1_1 );
	  obj_t B_2_1;
	  bli_acquire_mpart_t2b( BLIS_SUBPART1, idx3, bs3, &B_2, &B_2_1 );
	  //------------------------------------//

	  bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_ROW_PANELS, 
			       BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			       BLIS_BUFFER_FOR_A_BLOCK,
			       gemm_mr, gemm_mr, 
			       &U_21_1_1, &packed_A_blk );
	  bli_packm_blk_var2( &BLIS_ONE, &U_21_1_1, &packed_A_blk );
	  bli_gemm_ker_var2( &BLIS_MINUS_ONE, &packed_A_blk, &packed_B_pan, 
			     &BLIS_ONE, &B_2_1, (gemm_t*)NULL );

	  //------------------------------------//

	  //****
	}
	dimLen3 = bli_obj_length_after_trans( U_11_11 );
	for ( idx3 = 0; idx3 < dimLen3; idx3 += bs3 ) {
	  bs3 = bli_determine_blocksize_f( idx3, dimLen3, &U_11_11, gemm_mc );
	  dim_t idx4, dimLen4, bs4;
	  //****
	  obj_t U_11_11_1;
	  bli_acquire_mpart_t2b( BLIS_SUBPART1, idx3, bs3, &U_11_11, &U_11_11_1 );
	  obj_t B_1_1_1;
	  bli_acquire_mpart_t2b( BLIS_SUBPART1, idx3, bs3, &B_1_1, &B_1_1_1 );
	  //------------------------------------//

	  bli_obj_set_struc( BLIS_TRIANGULAR, U_11_11_1 );
	  bli_obj_set_uplo( BLIS_LOWER, U_11_11_1 );
	  bli_obj_set_diag( BLIS_UNIT_DIAG, U_11_11_1 );
	  bli_packm_init_pack( TRUE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_ROW_PANELS, 
			       BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			       BLIS_BUFFER_FOR_A_BLOCK,
			       gemm_mr, gemm_mr, 
			       &U_11_11_1, &packed_A_blk );
	  bli_packm_blk_var3( &BLIS_ONE, &U_11_11_1, &packed_A_blk );
	  bli_trmm_l_ker_var2( &BLIS_MINUS_ONE, &packed_A_blk, &packed_B_pan, 
			       &BLIS_ONE, &B_1_1_1, (trmm_t*)NULL );

	  //------------------------------------//

	  //****
	}
	dimLen3 = bli_obj_length_after_trans( B_1_2 );
	for ( idx3 = 0; idx3 < dimLen3; idx3 += bs3 ) {
	  bs3 = bli_determine_blocksize_f( idx3, dimLen3, &B_1_2, gemm_mc );
	  dim_t idx4, dimLen4, bs4;
	  //****
	  obj_t U_11_21_1;
	  bli_acquire_mpart_t2b( BLIS_SUBPART1, idx3, bs3, &U_11_21, &U_11_21_1 );
	  obj_t B_1_2_1;
	  bli_acquire_mpart_t2b( BLIS_SUBPART1, idx3, bs3, &B_1_2, &B_1_2_1 );
	  //------------------------------------//

	  bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_ROW_PANELS, 
			       BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			       BLIS_BUFFER_FOR_A_BLOCK,
			       gemm_mr, gemm_mr, 
			       &U_11_21_1, &packed_A_blk );
	  bli_packm_blk_var2( &BLIS_ONE, &U_11_21_1, &packed_A_blk );
	  bli_gemm_ker_var2( &BLIS_MINUS_ONE, &packed_A_blk, &packed_B_pan, 
			     &BLIS_ONE, &B_1_2_1, (gemm_t*)NULL );

	  //------------------------------------//

	  //****
	}

	//------------------------------------//

	//****
      }

      //------------------------------------//

      //****
    }

  }
  else {
dim_t idx1, dimLen1, bs1;
dimLen1 = min( bli_obj_length_after_trans( U ), bli_obj_width_after_trans( U ) );
for ( idx1 = 0; idx1 < dimLen1; idx1 += bs1 ) {
	bs1 = bli_determine_blocksize_f( idx1, dimLen1, &U, bs_obj );
	dim_t idx2, dimLen2, bs2;
//****
	obj_t B_1;
	bli_acquire_mpart_t2b( BLIS_SUBPART1, idx1, bs1, &B, &B_1 );
	obj_t B_2;
	bli_acquire_mpart_t2b( BLIS_SUBPART2, idx1, bs1, &B, &B_2 );
	obj_t U_11;
	bli_acquire_mpart_tl2br( BLIS_SUBPART11, idx1, bs1, &U, &U_11 );
	obj_t U_21;
	bli_acquire_mpart_tl2br( BLIS_SUBPART21, idx1, bs1, &U, &U_21 );
	obj_t T_1;
	bli_acquire_mpart_l2r( BLIS_SUBPART1, idx1, bs1, &T, &T_1 );
	//------------------------------------//

	obj_t T_1TL, T_1TLtmp;
	bli_acquire_mpart_l2r( BLIS_SUBPART1, 0, bli_obj_width( T_1 ), &T_1, &T_1TLtmp );
	bli_acquire_mpart_t2b( BLIS_SUBPART1, 0, bli_obj_length( B_1 ), &T_1TLtmp, &T_1TL );
	obj_t WTL, WTLtmp;
	bli_acquire_mpart_l2r( BLIS_SUBPART1, 0, bli_obj_width( B_1 ), &W, &WTLtmp );
	bli_acquire_mpart_t2b( BLIS_SUBPART1, 0, bli_obj_length( B_1 ), &WTLtmp, &WTL );
	dimLen2 = bli_obj_length_after_trans( WTL );
	for ( idx2 = 0; idx2 < dimLen2; idx2 += bs2 ) {
		bs2 = bli_determine_blocksize_f( idx2, dimLen2, &WTL, gemm_kc );
		dim_t idx3, dimLen3, bs3;
	//****
		obj_t U_11_10;
		bli_acquire_mpart_tl2br( BLIS_SUBPART10, idx2, bs2, &U_11, &U_11_10 );
		obj_t U_11_11;
		bli_acquire_mpart_tl2br( BLIS_SUBPART11, idx2, bs2, &U_11, &U_11_11 );
		obj_t B_1_1;
		bli_acquire_mpart_t2b( BLIS_SUBPART1, idx2, bs2, &B_1, &B_1_1 );
		obj_t WTL_0;
		bli_acquire_mpart_t2b( BLIS_SUBPART0, idx2, bs2, &WTL, &WTL_0 );
		obj_t WTL_1;
		bli_acquire_mpart_t2b( BLIS_SUBPART1, idx2, bs2, &WTL, &WTL_1 );
		//------------------------------------//

		bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_COL_PANELS, 
				BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
				BLIS_BUFFER_FOR_B_PANEL,
				gemm_mr, gemm_nr, 
				&B_1_1, &packed_B_pan );
		bli_packm_blk_var2( &BLIS_ONE, &B_1_1, &packed_B_pan );
		dimLen3 = bli_obj_length_after_trans( WTL_0 );
		for ( idx3 = 0; idx3 < dimLen3; idx3 += bs3 ) {
			bs3 = bli_determine_blocksize_f( idx3, dimLen3, &WTL_0, gemm_mc );
			dim_t idx4, dimLen4, bs4;
		//****
			obj_t U_11_10_1;
			bli_acquire_mpart_l2r( BLIS_SUBPART1, idx3, bs3, &U_11_10, &U_11_10_1 );
			obj_t WTL_0_1;
			bli_acquire_mpart_t2b( BLIS_SUBPART1, idx3, bs3, &WTL_0, &WTL_0_1 );
			//------------------------------------//

			obj_t U_11_10_1H;
			bli_obj_alias_with_trans( BLIS_CONJ_TRANSPOSE, U_11_10_1, U_11_10_1H);
			bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_ROW_PANELS, 
					BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
					BLIS_BUFFER_FOR_A_BLOCK,
					gemm_mr, gemm_mr, 
					&U_11_10_1H, &packed_A_blk );
			bli_packm_blk_var2( &BLIS_ONE, &U_11_10_1H, &packed_A_blk );
			bli_gemm_ker_var2( &BLIS_ONE, &packed_A_blk, &packed_B_pan, 
					&BLIS_ONE, &WTL_0_1, (gemm_t*)NULL );

			//------------------------------------//

		//****
		}
		dimLen3 = bli_obj_width_after_trans( U_11_11 );
		for ( idx3 = 0; idx3 < dimLen3; idx3 += bs3 ) {
			bs3 = bli_determine_blocksize_f( idx3, dimLen3, &U_11_11, gemm_mc );
			dim_t idx4, dimLen4, bs4;
		//****
			obj_t U_11_11_1;
			bli_acquire_mpart_l2r( BLIS_SUBPART1, idx3, bs3, &U_11_11, &U_11_11_1 );
			obj_t WTL_1_1;
			bli_acquire_mpart_t2b( BLIS_SUBPART1, idx3, bs3, &WTL_1, &WTL_1_1 );
			//------------------------------------//

			obj_t U_11_11_1H;
			bli_obj_alias_with_trans( BLIS_CONJ_TRANSPOSE, U_11_11_1, U_11_11_1H);
			bli_obj_set_struc( BLIS_TRIANGULAR, U_11_11_1H );
			bli_obj_set_uplo( BLIS_LOWER, U_11_11_1H );
			bli_obj_set_diag( BLIS_UNIT_DIAG, U_11_11_1H );
			bli_packm_init_pack( TRUE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_ROW_PANELS, 
					BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
					BLIS_BUFFER_FOR_A_BLOCK,
					gemm_mr, gemm_mr, 
					&U_11_11_1H, &packed_A_blk );
			bli_packm_blk_var3( &BLIS_ONE, &U_11_11_1H, &packed_A_blk );
			bli_trmm_u_ker_var2( &BLIS_ONE, &packed_A_blk, &packed_B_pan, 
					&BLIS_ZERO, &WTL_1_1, (trmm_t*)NULL );

			//------------------------------------//

		//****
		}

		//------------------------------------//

	//****
	}
	dimLen2 = bli_obj_length_after_trans( U_21 );
	for ( idx2 = 0; idx2 < dimLen2; idx2 += bs2 ) {
		bs2 = bli_determine_blocksize_f( idx2, dimLen2, &U_21, gemm_kc );
		dim_t idx3, dimLen3, bs3;
	//****
		obj_t U_21_1;
		bli_acquire_mpart_t2b( BLIS_SUBPART1, idx2, bs2, &U_21, &U_21_1 );
		obj_t B_2_1;
		bli_acquire_mpart_t2b( BLIS_SUBPART1, idx2, bs2, &B_2, &B_2_1 );
		//------------------------------------//

		bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_COL_PANELS, 
				BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
				BLIS_BUFFER_FOR_B_PANEL,
				gemm_kr, gemm_nr, 
				&B_2_1, &packed_B_pan );
		bli_packm_blk_var2( &BLIS_ONE, &B_2_1, &packed_B_pan );
		dimLen3 = bli_obj_length_after_trans( WTL );
		for ( idx3 = 0; idx3 < dimLen3; idx3 += bs3 ) {
			bs3 = bli_determine_blocksize_f( idx3, dimLen3, &WTL, gemm_mc );
			dim_t idx4, dimLen4, bs4;
		//****
			obj_t U_21_1_1;
			bli_acquire_mpart_l2r( BLIS_SUBPART1, idx3, bs3, &U_21_1, &U_21_1_1 );
			obj_t WTL_1;
			bli_acquire_mpart_t2b( BLIS_SUBPART1, idx3, bs3, &WTL, &WTL_1 );
			//------------------------------------//

			obj_t U_21_1_1H;
			bli_obj_alias_with_trans( BLIS_CONJ_TRANSPOSE, U_21_1_1, U_21_1_1H);
			bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_ROW_PANELS, 
					BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
					BLIS_BUFFER_FOR_A_BLOCK,
					gemm_mr, gemm_kr, 
					&U_21_1_1H, &packed_A_blk );
			bli_packm_blk_var2( &BLIS_ONE, &U_21_1_1H, &packed_A_blk );
			bli_gemm_ker_var2( &BLIS_ONE, &packed_A_blk, &packed_B_pan, 
					&BLIS_ONE, &WTL_1, (gemm_t*)NULL );

			//------------------------------------//

		//****
		}

		//------------------------------------//

	//****
	}
	dimLen2 = bli_obj_length_after_trans( WTL );
	for ( idx2 = 0; idx2 < dimLen2; idx2 += bs2 ) {
		bs2 = bli_determine_blocksize_f( idx2, dimLen2, &WTL, gemm_kc );
		dim_t idx3, dimLen3, bs3;
	//****
		obj_t T_1TL_11;
		bli_acquire_mpart_tl2br( BLIS_SUBPART11, idx2, bs2, &T_1TL, &T_1TL_11 );
		obj_t T_1TL_12;
		bli_acquire_mpart_tl2br( BLIS_SUBPART12, idx2, bs2, &T_1TL, &T_1TL_12 );
		obj_t WTL_1;
		bli_acquire_mpart_t2b( BLIS_SUBPART1, idx2, bs2, &WTL, &WTL_1 );
		obj_t WTL_2;
		bli_acquire_mpart_t2b( BLIS_SUBPART2, idx2, bs2, &WTL, &WTL_2 );
		//------------------------------------//

		bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_COL_PANELS, 
				BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
				BLIS_BUFFER_FOR_B_PANEL,
				gemm_mr, gemm_nr, 
				&WTL_1, &packed_B_pan );
		bli_packm_blk_var2( &BLIS_ONE, &WTL_1, &packed_B_pan );
		dimLen3 = bli_obj_width_after_trans( T_1TL_11 );
		for ( idx3 = 0; idx3 < dimLen3; idx3 += bs3 ) {
			bs3 = bli_determine_blocksize_f( idx3, dimLen3, &T_1TL_11, gemm_mc );
			dim_t idx4, dimLen4, bs4;
		//****
			obj_t T_1TL_11_1;
			bli_acquire_mpart_l2r( BLIS_SUBPART1, idx3, bs3, &T_1TL_11, &T_1TL_11_1 );
			obj_t WTL_1_1;
			bli_acquire_mpart_t2b( BLIS_SUBPART1, idx3, bs3, &WTL_1, &WTL_1_1 );
			//------------------------------------//

			obj_t T_1TL_11_1H;
			bli_obj_alias_with_trans( BLIS_CONJ_TRANSPOSE, T_1TL_11_1, T_1TL_11_1H);
			bli_obj_set_struc( BLIS_TRIANGULAR, T_1TL_11_1H );
			bli_obj_set_uplo( BLIS_UPPER, T_1TL_11_1H );
			bli_packm_init_pack( TRUE, BLIS_INVERT_DIAG, BLIS_PACKED_ROW_PANELS, 
					BLIS_PACK_REV_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
					BLIS_BUFFER_FOR_A_BLOCK,
					gemm_mr, gemm_mr, 
					&T_1TL_11_1H, &packed_A_blk );
			bli_packm_blk_var3( &BLIS_ONE, &T_1TL_11_1H, &packed_A_blk );
			bli_trsm_l_ker_var2( &BLIS_ONE, &packed_A_blk, &packed_B_pan, 
					&BLIS_ZERO, &WTL_1_1, (trsm_t*)NULL );

			//------------------------------------//

		//****
		}
		dimLen3 = bli_obj_length_after_trans( WTL_2 );
		for ( idx3 = 0; idx3 < dimLen3; idx3 += bs3 ) {
			bs3 = bli_determine_blocksize_f( idx3, dimLen3, &WTL_2, gemm_mc );
			dim_t idx4, dimLen4, bs4;
		//****
			obj_t T_1TL_12_1;
			bli_acquire_mpart_l2r( BLIS_SUBPART1, idx3, bs3, &T_1TL_12, &T_1TL_12_1 );
			obj_t WTL_2_1;
			bli_acquire_mpart_t2b( BLIS_SUBPART1, idx3, bs3, &WTL_2, &WTL_2_1 );
			//------------------------------------//

			obj_t T_1TL_12_1H;
			bli_obj_alias_with_trans( BLIS_CONJ_TRANSPOSE, T_1TL_12_1, T_1TL_12_1H);
			bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_ROW_PANELS, 
					BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
					BLIS_BUFFER_FOR_A_BLOCK,
					gemm_mr, gemm_mr, 
					&T_1TL_12_1H, &packed_A_blk );
			bli_packm_blk_var2( &BLIS_ONE, &T_1TL_12_1H, &packed_A_blk );
			bli_gemm_ker_var2( &BLIS_MINUS_ONE, &packed_A_blk, &packed_B_pan, 
					&BLIS_ONE, &WTL_2_1, (gemm_t*)NULL );

			//------------------------------------//

		//****
		}

		//------------------------------------//

	//****
	}
	dimLen2 = bli_obj_width_after_trans( U_21 );
	for ( idx2 = 0; idx2 < dimLen2; idx2 += bs2 ) {
		bs2 = bli_determine_blocksize_b( idx2, dimLen2, &U_21, gemm_kc );
		dim_t idx3, dimLen3, bs3;
	//****
		obj_t U_21_1;
		bli_acquire_mpart_r2l( BLIS_SUBPART1, idx2, bs2, &U_21, &U_21_1 );
		obj_t WTL_1;
		bli_acquire_mpart_b2t( BLIS_SUBPART1, idx2, bs2, &WTL, &WTL_1 );
		obj_t U_11_11;
		bli_acquire_mpart_br2tl( BLIS_SUBPART11, idx2, bs2, &U_11, &U_11_11 );
		obj_t U_11_21;
		bli_acquire_mpart_br2tl( BLIS_SUBPART21, idx2, bs2, &U_11, &U_11_21 );
		obj_t B_1_1;
		bli_acquire_mpart_b2t( BLIS_SUBPART1, idx2, bs2, &B_1, &B_1_1 );
		obj_t B_1_2;
		bli_acquire_mpart_b2t( BLIS_SUBPART2, idx2, bs2, &B_1, &B_1_2 );
		//------------------------------------//

		bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_COL_PANELS, 
				BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
				BLIS_BUFFER_FOR_B_PANEL,
				gemm_mr, gemm_nr, 
				&WTL_1, &packed_B_pan );
		bli_packm_blk_var2( &BLIS_ONE, &WTL_1, &packed_B_pan );
		dimLen3 = bli_obj_length_after_trans( B_2 );
		for ( idx3 = 0; idx3 < dimLen3; idx3 += bs3 ) {
			bs3 = bli_determine_blocksize_f( idx3, dimLen3, &B_2, gemm_mc );
			dim_t idx4, dimLen4, bs4;
		//****
			obj_t U_21_1_1;
			bli_acquire_mpart_t2b( BLIS_SUBPART1, idx3, bs3, &U_21_1, &U_21_1_1 );
			obj_t B_2_1;
			bli_acquire_mpart_t2b( BLIS_SUBPART1, idx3, bs3, &B_2, &B_2_1 );
			//------------------------------------//

			bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_ROW_PANELS, 
					BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
					BLIS_BUFFER_FOR_A_BLOCK,
					gemm_mr, gemm_mr, 
					&U_21_1_1, &packed_A_blk );
			bli_packm_blk_var2( &BLIS_ONE, &U_21_1_1, &packed_A_blk );
			bli_gemm_ker_var2( &BLIS_MINUS_ONE, &packed_A_blk, &packed_B_pan, 
					&BLIS_ONE, &B_2_1, (gemm_t*)NULL );

			//------------------------------------//

		//****
		}
		dimLen3 = bli_obj_length_after_trans( U_11_11 );
		for ( idx3 = 0; idx3 < dimLen3; idx3 += bs3 ) {
			bs3 = bli_determine_blocksize_f( idx3, dimLen3, &U_11_11, gemm_mc );
			dim_t idx4, dimLen4, bs4;
		//****
			obj_t U_11_11_1;
			bli_acquire_mpart_t2b( BLIS_SUBPART1, idx3, bs3, &U_11_11, &U_11_11_1 );
			obj_t B_1_1_1;
			bli_acquire_mpart_t2b( BLIS_SUBPART1, idx3, bs3, &B_1_1, &B_1_1_1 );
			//------------------------------------//

			bli_obj_set_struc( BLIS_TRIANGULAR, U_11_11_1 );
			bli_obj_set_uplo( BLIS_LOWER, U_11_11_1 );
			bli_obj_set_diag( BLIS_UNIT_DIAG, U_11_11_1 );
			bli_packm_init_pack( TRUE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_ROW_PANELS, 
					BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
					BLIS_BUFFER_FOR_A_BLOCK,
					gemm_mr, gemm_mr, 
					&U_11_11_1, &packed_A_blk );
			bli_packm_blk_var3( &BLIS_ONE, &U_11_11_1, &packed_A_blk );
			bli_trmm_l_ker_var2( &BLIS_MINUS_ONE, &packed_A_blk, &packed_B_pan, 
					&BLIS_ONE, &B_1_1_1, (trmm_t*)NULL );

			//------------------------------------//

		//****
		}
		dimLen3 = bli_obj_length_after_trans( B_1_2 );
		for ( idx3 = 0; idx3 < dimLen3; idx3 += bs3 ) {
			bs3 = bli_determine_blocksize_f( idx3, dimLen3, &B_1_2, gemm_mc );
			dim_t idx4, dimLen4, bs4;
		//****
			obj_t U_11_21_1;
			bli_acquire_mpart_t2b( BLIS_SUBPART1, idx3, bs3, &U_11_21, &U_11_21_1 );
			obj_t B_1_2_1;
			bli_acquire_mpart_t2b( BLIS_SUBPART1, idx3, bs3, &B_1_2, &B_1_2_1 );
			//------------------------------------//

			bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_ROW_PANELS, 
					BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
					BLIS_BUFFER_FOR_A_BLOCK,
					gemm_mr, gemm_mr, 
					&U_11_21_1, &packed_A_blk );
			bli_packm_blk_var2( &BLIS_ONE, &U_11_21_1, &packed_A_blk );
			bli_gemm_ker_var2( &BLIS_MINUS_ONE, &packed_A_blk, &packed_B_pan, 
					&BLIS_ONE, &B_1_2_1, (gemm_t*)NULL );

			//------------------------------------//

		//****
		}

		//------------------------------------//

	//****
	}

	//------------------------------------//

//****
}

   }

  bli_obj_release_pack( &packed_A_blk );
  bli_obj_release_pack( &packed_B_pan );
}

void DxT_QR( obj_t *a1,
	     obj_t *t1,
	     blksz_t *bs )
{
  FLA_Error err;
  FLA_Obj A, T;

  err = wrap_blis_obj( a1, &A );
  if (err != FLA_SUCCESS)
    return;

  err = wrap_blis_obj( t1, &T );
  if (err != FLA_SUCCESS)
    return;

  FLA_Obj ATL,   ATR,      A00, A01, A02, 
          ABL,   ABR,      A10, A11, A12,
                           A20, A21, A22;

  FLA_Obj TL,    TR,       T0,  T1,  W12;
  FLA_Obj T1T,   T2B;

  FLA_Obj AB1,   AB2;

  dim_t   b_alg, b;

  // Query the algorithmic blocksize by inspecting the length of T.
  b_alg = FLA_Obj_length( T );

  FLA_Part_2x2( A,    &ATL, &ATR,
                      &ABL, &ABR,     0, 0, FLA_TL );

  FLA_Part_1x2( T,    &TL,  &TR,      0, FLA_LEFT );

  while ( FLA_Obj_min_dim( ABR ) > 0 ){

    b = min( b_alg, FLA_Obj_min_dim( ABR ) );

    FLA_Repart_2x2_to_3x3( ATL, /**/ ATR,       &A00, /**/ &A01, &A02,
                        /* ************* */   /* ******************** */
                                                &A10, /**/ &A11, &A12,
                           ABL, /**/ ABR,       &A20, /**/ &A21, &A22,
                           b, b, FLA_BR );

    FLA_Repart_1x2_to_1x3( TL,  /**/ TR,        &T0, /**/ &T1, &W12,
                           b, FLA_RIGHT );

    /*------------------------------------------------------------*/

    FLA_Part_2x1( T1,   &T1T, 
                        &T2B,    b, FLA_TOP );

    FLA_Merge_2x1( A11,
                   A21,   &AB1 );

    // Perform a QR factorization via the UT transform on AB1:
    //
    //   / A11 \ -> QB1 R11
    //   \ A21 /
    //
    // where:
    //  - QB1 is formed from UB1 (which is stored column-wise below the
    //    diagonal of AB1) and T11 (which is stored to the upper triangle
    //    of T11).
    //  - R11 is stored to the upper triangle of AB1.
    
    FLA_QR_UT_internal( AB1, T1T, qrut_cntl );

    if ( FLA_Obj_width( A12 ) > 0 )
    {
      FLA_Merge_2x1( A12,
                     A22,   &AB2 );

      // Apply the Householder transforms associated with UB1 and T11 to 
      // AB2:
      //
      //   / A12 \ := QB1' / A12 \
      //   \ A22 /         \ A22 /
      //
      // where QB1 is formed from UB1 and T11.

      DxT_AppHouse( AB1, T1T, W12, AB2, bs );

      /*
      FLA_Apply_Q_UT_internal( FLA_LEFT, FLA_CONJ_TRANSPOSE, FLA_FORWARD, FLA_COLUMNWISE,
                               AB1, T1T, W12, AB2,
			       fla_apqut_cntl_leaf);
      */

    }

    /*------------------------------------------------------------*/

    FLA_Cont_with_3x3_to_2x2( &ATL, /**/ &ATR,       A00, A01, /**/ A02,
                                                     A10, A11, /**/ A12,
                            /* ************** */  /* ****************** */
                              &ABL, /**/ &ABR,       A20, A21, /**/ A22,
                              FLA_TL );

    FLA_Cont_with_1x3_to_1x2( &TL,  /**/ &TR,        T0, T1, /**/ W12,
                              FLA_LEFT );
  }

  FLA_Obj_free_without_buffer(&A);
  FLA_Obj_free_without_buffer(&T);
}


void FLA_QR( obj_t *a2, 
	     obj_t *t2 )
{
  FLA_Error err;
  FLA_Obj A, T;

  err = wrap_blis_obj( a2, &A );
  if (err != FLA_SUCCESS)
    return;

  err = wrap_blis_obj( t2, &T );
  if (err != FLA_SUCCESS)
    return;

  err = FLA_QR_UT( A, T );

  if (err != FLA_SUCCESS) {
    double *tmp = NULL;
    *tmp = 0;
    printf("FLA_QR_UT error\n");
  }

  FLA_Obj_free_without_buffer( &A );
  FLA_Obj_free_without_buffer( &T );
}


int main( int argc, char** argv )
{
  obj_t a1, a2;
  obj_t a_save;
  obj_t t1, t2;
  obj_t t_save;
#ifdef TESTLIB
  obj_t normValA;
#endif
  dim_t m;
  dim_t p;
  dim_t p_begin, p_end, p_inc;
  dim_t bs_begin, bs_end, bs_inc;
  int   m_input, n_input, k_input;
  int   r, n_repeats;
  dim_t bsDxTBest, bsFLABest;
  double dtimeDxT;
  double dtime_saveDxT;
  double gflopsDxT;
  double gflopsDxTBest;
  double dtimeFLA;
  double dtime_saveFLA;
  double gflopsFLA;
  double gflopsFLABest;

  blksz_t *bs_obj;

  bli_init();
  FLA_Init();

  qrut_cntl = FLA_Cntl_qrut_obj_create( FLA_FLAT,
					FLA_UNB_OPT_VARIANT2,
					NULL,
					NULL,
					NULL ); 

  bli_error_checking_level_set( BLIS_NO_ERROR_CHECKING );

  n_repeats = 3;

  p_begin = 200;
  p_end   = 5000;
  p_inc   = 200;
  bs_begin = 64;
  bs_end   = 256;
  bs_inc   = 16;

  m_input = -1;
#ifdef TESTLIB
  bli_obj_create( BLIS_DOUBLE, 1, 1, 0, 0, &normValA );
#endif

  for ( p = p_begin; p <= p_end; p += p_inc )
    {

      if ( m_input < 0 ) m = p * ( dim_t )abs(m_input);
      else               m =     ( dim_t )    m_input;

      bli_obj_create( BLIS_DOUBLE, m, m, 0, 0, &a1 );
      bli_obj_create( BLIS_DOUBLE, m, m, 0, 0, &a2 );
      bli_obj_create( BLIS_DOUBLE, m, m, 0, 0, &a_save );

      bli_randm( &a1 );

      bli_copym( &a1, &a_save );

      gflopsDxTBest = 0;
      gflopsFLABest = 0;

      for ( dim_t bs = bs_begin; bs <= bs_end; bs += bs_inc ) {
	bs_obj = bli_blksz_obj_create(bs, 0,
				      bs, 0,
				      bs, 0,
				      bs, 0);

	bli_obj_create( BLIS_DOUBLE, bs, m, 0, 0, &t1 );
	bli_obj_create( BLIS_DOUBLE, bs, m, 0, 0, &t2 );
	bli_obj_create( BLIS_DOUBLE, bs, m, 0, 0, &t_save );
	bli_randm( &t1 );
	bli_copym( &t1, &t_save );
	
	
	dtime_saveDxT = 1.0e9;
	dtime_saveFLA = 1.0e9;

	for ( r = 0; r < n_repeats; ++r )
	  {
	    bli_copym( &a_save, &a1 );
	    bli_copym( &a_save, &a2 );
	    bli_copym( &t_save, &t1 );
	    bli_copym( &t_save, &t2 );

	    dtimeFLA = bli_clock();
	    FLA_QR( &a2, &t2 );
	    dtime_saveFLA = bli_clock_min_diff( dtime_saveFLA, dtimeFLA );

	    dtimeDxT = bli_clock();
	    DxT_QR( &a1, &t1, bs_obj );
	    dtime_saveDxT = bli_clock_min_diff( dtime_saveDxT, dtimeDxT );

	    //We don't expect the T's to be the same because of the
	    // differences in the way DxT and FLAME deal with the upper 
	    //triangular portions
	  
#ifdef TESTLIB
	    bli_axpym( &BLIS_MINUS_ONE, &a1, &a2 );
	    bli_fnormm( &a2, &normValA );
#endif
	  }

	gflopsDxT = ( (4.0 / 3.0) * m * m * m ) / ( dtime_saveDxT * 1.0e9 );
	if (gflopsDxT > gflopsDxTBest) {
	  gflopsDxTBest = gflopsDxT;
	  bsDxTBest = bs;
	}
	gflopsFLA = ( (4.0 / 3.0) * m * m * m ) / ( dtime_saveFLA * 1.0e9 );
	if (gflopsFLA > gflopsFLABest) {
	  gflopsFLABest = gflopsFLA;
	  bsFLABest = bs;
	}

	printf( "data_QR_DxT" );
	printf( "( %2lu, %2lu, 1:4 ) = [ %4lu %4lu %10.3e  %6.3f ];\n",
		(p - p_begin + 1)/p_inc + 1, (bs - bs_begin + 1)/bs_inc + 1,
		m, bs, dtime_saveDxT, gflopsDxT );
	printf( "data_QR_FLA" );
	printf( "( %2lu, %2lu, 1:4 ) = [ %4lu %4lu %10.3e  %6.3f ];\n\n",
		(p - p_begin + 1)/p_inc + 1, (bs - bs_begin + 1)/bs_inc + 1,
		m, bs, dtime_saveFLA, gflopsFLA );

#ifdef TESTLIB
      bli_printm( "%%%%NORM A ", &normValA, "%4.1f", "" );
#endif

	bli_obj_free( &t1 );
	bli_obj_free( &t2 );
	bli_obj_free( &t_save );

	bli_blksz_obj_free( bs_obj );
      }

      printf( "data_QR_DxT_Best" );
      printf( "( %2lu, 1:2 ) = [ %4lu  %6.3f ];  %%%% %lu\n",
	      (p - p_begin + 1)/p_inc + 1,
	      m, gflopsDxTBest, bsDxTBest );
      
      printf( "data_QR_FLA_Best" );
      printf( "( %2lu, 1:2 ) = [ %4lu  %6.3f ];  %%%% %lu\n",
	      (p - p_begin + 1)/p_inc + 1,
	      m, gflopsFLABest, bsFLABest);

      printf("%%%%%%%%%%\n");
      fflush(stdout);

      bli_obj_free( &a1 );
      bli_obj_free( &a2 );
      bli_obj_free( &a_save );
    }

  FLA_Finalize();
  bli_finalize();

  return 0;
}

