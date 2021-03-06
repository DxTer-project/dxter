#include <unistd.h>
#include "blis.h"
#include "FLAME.h"

//#define TESTLIB

extern fla_axpy_t*  fla_axpy_cntl_blas;
extern fla_gemm_t*  fla_gemm_cntl_blas;
extern fla_hemm_t*  fla_hemm_cntl_blas;
extern fla_her2k_t* fla_her2k_cntl_blas;
extern fla_trmm_t*  fla_trmm_cntl_blas;
extern fla_trsm_t*  fla_trsm_cntl_blas;
extern blksz_t *gemm_mc;
extern blksz_t *gemm_kc;
extern blksz_t *gemm_nc;
extern blksz_t *gemm_mr;
extern blksz_t *gemm_kr;
extern blksz_t *gemm_nr;
/*
fla_eig_gest_t*     fla_eig_gest_ix_cntl;
fla_eig_gest_t*     fla_eig_gest_nx_cntl;
*/
fla_eig_gest_t*     fla_eig_gest_ix_cntl_leaf;
fla_eig_gest_t*     fla_eig_gest_nx_cntl_leaf;

FLA_Error wrap_blis_obj( obj_t *blis,
		    FLA_Obj *fla )
{
  FLA_Error err;

  FLA_Datatype type;
  
  switch( bli_obj_datatype(*blis) )
    {
    case (BLIS_DOUBLE):
      type = FLA_DOUBLE;
      break;
    case (BLIS_INT):
      type = FLA_INT;
      break;
    default:
      return 0;
      break;
    }


  err = FLA_Obj_create_without_buffer( type, 
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

void blis_obj_merge_vert( obj_t *top,
			  obj_t *bott,
			  obj_t *merged )
{
  if (bli_obj_datatype(*top) != BLIS_DOUBLE) {
    printf("bad datatype in merge_vert\n");
    return;
  }
	
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
      (((double*)bli_obj_buffer_at_off(*top) + (bli_obj_length(*top) * bli_obj_row_stride(*top)))
       != (double*)bli_obj_buffer_at_off(*bott)))
    {
      printf("top %p\n",(double*)bli_obj_buffer_at_off(*top));
      printf("bot %p\n",(double*)bli_obj_buffer_at_off(*bott));
      printf("top col stride %lu\n", bli_obj_col_stride(*top));
      printf("top row stride %lu\n", bli_obj_row_stride(*top));
      printf("top+offset %p\n",((double*)bli_obj_buffer_at_off(*top) + (bli_obj_length(*top) * bli_obj_row_stride(*top))));
      printf("top length %lu\n", bli_obj_length(*top));
      printf("top width %lu\n", bli_obj_width(*top));
      printf("bott length %lu\n", bli_obj_length(*bott));
      printf("bott width %lu\n", bli_obj_width(*bott));
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

void libflame_Hegst(obj_t *A11,
		    obj_t *B11 )
{
  FLA_Obj A, B;
  wrap_blis_obj( A11, &A );
  wrap_blis_obj( B11, &B );
  
  FLA_Obj tmp;
  FLA_Obj_create( FLA_DOUBLE, 	 
		  bli_obj_length(*A11), bli_obj_width(*A11),
		  0, 0, &tmp );

#if 1
  FLA_Eig_gest_il_opt_var3(A, tmp, B);
#elif 1
  FLA_Eig_gest_internal( FLA_INVERSE, FLA_LOWER_TRIANGULAR,
  			 A, tmp, B,
  			 fla_eig_gest_ix_cntl_leaf);
#else
  FLA_Eig_gest(FLA_INVERSE, FLA_LOWER_TRIANGULAR, A, B);
#endif  

  FLA_Obj_free(&tmp);
  FLA_Obj_free_without_buffer(&A);
  FLA_Obj_free_without_buffer(&B);
}

//var 2 in libflame
void DxT_TwoSidedTrsm( obj_t A,
		       obj_t L,
		       blksz_t *bs_obj)
{
  obj_t packed_A_blk;
  obj_t packed_B_pan;

  bli_obj_init_pack( &packed_A_blk );
  bli_obj_init_pack( &packed_B_pan );

/* 	obj_t Y21; */
/* 	bli_obj_create(BLIS_DOUBLE, bli_obj_length(L_21), bli_obj_width(L_21), */
/* 		   0, 0, &Y21); */

  if (bli_blksz_for_obj(&A,bs_obj) <= bli_blksz_for_obj(&A,gemm_kc)) {
    dim_t idx1, dimLen1, bs1;
    dimLen1 = min( bli_obj_length_after_trans( A ), bli_obj_width_after_trans( A ) );
    for ( idx1 = 0; idx1 < dimLen1; idx1 += bs1 ) {
      bs1 = bli_determine_blocksize_f( idx1, dimLen1, &A, bs_obj );
      dim_t idx2, dimLen2, bs2;
      //****
      obj_t A_10;
      bli_acquire_mpart_tl2br( BLIS_SUBPART10, idx1, bs1, &A, &A_10 );
      obj_t A_20;
      bli_acquire_mpart_tl2br( BLIS_SUBPART20, idx1, bs1, &A, &A_20 );
      obj_t A_11;
      bli_acquire_mpart_tl2br( BLIS_SUBPART11, idx1, bs1, &A, &A_11 );
      obj_t A_21;
      bli_acquire_mpart_tl2br( BLIS_SUBPART21, idx1, bs1, &A, &A_21 );
      obj_t A_22;
      bli_acquire_mpart_tl2br( BLIS_SUBPART22, idx1, bs1, &A, &A_22 );
      obj_t L_11;
      bli_acquire_mpart_tl2br( BLIS_SUBPART11, idx1, bs1, &L, &L_11 );
      obj_t L_21;
      bli_acquire_mpart_tl2br( BLIS_SUBPART21, idx1, bs1, &L, &L_21 );
      //------------------------------------//

	obj_t Y21;
	bli_obj_create(BLIS_DOUBLE, bli_obj_length(L_21), bli_obj_width(L_21),
		   0, 0, &Y21);

      libflame_Hegst( &A_11, &L_11 );
      obj_t A_11C;
      bli_obj_alias_with_conj( BLIS_CONJUGATE, A_11, A_11C);
      obj_t L_21T;
      bli_obj_alias_with_trans( BLIS_TRANSPOSE, L_21, L_21T);
      bli_obj_induce_trans( Y21 );
      bli_obj_induce_trans( A_21 );
      bli_obj_set_struc( BLIS_TRIANGULAR, A_22 );
      bli_obj_set_uplo( BLIS_LOWER, A_22 );
      bli_scalm( &BLIS_ONE, &A_22 );
      bli_obj_set_struc( BLIS_GENERAL, A_22 );
      bli_obj_set_uplo( BLIS_DENSE, A_22 );
      dimLen2 = bli_obj_length_after_trans( A_21 );
      for ( idx2 = 0; idx2 < dimLen2; idx2 += bs2 ) {
	bs2 = bli_determine_blocksize_b( idx2, dimLen2, &A_21, gemm_kc );
	dim_t idx3, dimLen3, bs3;
	//****
	obj_t L_11_10;
	bli_acquire_mpart_br2tl( BLIS_SUBPART10, idx2, bs2, &L_11, &L_11_10 );
	obj_t L_11_11;
	bli_acquire_mpart_br2tl( BLIS_SUBPART11, idx2, bs2, &L_11, &L_11_11 );
	obj_t A_21_0;
	bli_acquire_mpart_b2t( BLIS_SUBPART0, idx2, bs2, &A_21, &A_21_0 );
	obj_t A_21_1;
	bli_acquire_mpart_b2t( BLIS_SUBPART1, idx2, bs2, &A_21, &A_21_1 );
	//------------------------------------//

	bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_COL_PANELS, 
			     BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			     BLIS_BUFFER_FOR_B_PANEL,
			     gemm_mr, gemm_nr, 
			     &A_21_1, &packed_B_pan );
	bli_packm_blk_var2( &BLIS_ONE, &A_21_1, &packed_B_pan );
	dimLen3 = bli_obj_length_after_trans( L_11_11 );
	for ( idx3 = 0; idx3 < dimLen3; idx3 += bs3 ) {
	  bs3 = bli_determine_blocksize_f( idx3, dimLen3, &L_11_11, gemm_mc );
	  dim_t idx4, dimLen4, bs4;
	  //****
	  obj_t L_11_11_1;
	  bli_acquire_mpart_t2b( BLIS_SUBPART1, idx3, bs3, &L_11_11, &L_11_11_1 );
	  obj_t A_21_1_1;
	  bli_acquire_mpart_t2b( BLIS_SUBPART1, idx3, bs3, &A_21_1, &A_21_1_1 );
	  //------------------------------------//

	  obj_t L_11_11_1C;
	  bli_obj_alias_with_conj( BLIS_CONJUGATE, L_11_11_1, L_11_11_1C);
	  bli_obj_set_struc( BLIS_TRIANGULAR, L_11_11_1C );
	  bli_obj_set_uplo( BLIS_LOWER, L_11_11_1C );
	  bli_packm_init_pack( TRUE, BLIS_INVERT_DIAG, BLIS_PACKED_ROW_PANELS, 
			       BLIS_PACK_REV_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			       BLIS_BUFFER_FOR_A_BLOCK,
			       gemm_mr, gemm_mr, 
			       &L_11_11_1C, &packed_A_blk );
	  bli_packm_blk_var3( &BLIS_ONE, &L_11_11_1C, &packed_A_blk );
	  bli_trsm_l_ker_var2( &BLIS_ONE, &packed_A_blk, &packed_B_pan, 
			       &BLIS_ZERO, &A_21_1_1, (trsm_t*)NULL );

	  //------------------------------------//

	  //****
	}
	dimLen3 = bli_obj_length_after_trans( A_21_0 );
	for ( idx3 = 0; idx3 < dimLen3; idx3 += bs3 ) {
	  bs3 = bli_determine_blocksize_f( idx3, dimLen3, &A_21_0, gemm_mc );
	  dim_t idx4, dimLen4, bs4;
	  //****
	  obj_t L_11_10_1;
	  bli_acquire_mpart_l2r( BLIS_SUBPART1, idx3, bs3, &L_11_10, &L_11_10_1 );
	  obj_t A_21_0_1;
	  bli_acquire_mpart_t2b( BLIS_SUBPART1, idx3, bs3, &A_21_0, &A_21_0_1 );
	  //------------------------------------//

	  obj_t L_11_10_1C;
	  bli_obj_alias_with_conj( BLIS_CONJUGATE, L_11_10_1, L_11_10_1C);
	  bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_ROW_PANELS, 
			       BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			       BLIS_BUFFER_FOR_A_BLOCK,
			       gemm_mr, gemm_kr, 
			       &L_11_10_1C, &packed_A_blk );
	  bli_packm_blk_var2( &BLIS_ONE, &L_11_10_1C, &packed_A_blk );
	  bli_gemm_ker_var2( &BLIS_MINUS_ONE, &packed_A_blk, &packed_B_pan, 
			     &BLIS_ONE, &A_21_0_1, (gemm_t*)NULL );

	  //------------------------------------//

	  //****
	}

	//------------------------------------//

	//****
      }
      //****
      //------------------------------------//

      bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_COL_PANELS, 
			   BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			   BLIS_BUFFER_FOR_B_PANEL,
			   gemm_kr, gemm_nr, 
			   &L_21T, &packed_B_pan );
      bli_packm_blk_var2( &BLIS_ONE, &L_21T, &packed_B_pan );
	dim_t idx3, dimLen3, bs3;
      dimLen3 = bli_obj_length_after_trans( Y21 );
      for ( idx3 = 0; idx3 < dimLen3; idx3 += bs3 ) {
	bs3 = bli_determine_blocksize_f( idx3, dimLen3, &Y21, gemm_mc );
	dim_t idx4, dimLen4, bs4;
	//****
	obj_t A_11C_1;
	bli_acquire_mpart_t2b( BLIS_SUBPART1, idx3, bs3, &A_11C, &A_11C_1 );
	obj_t Y21_1;
	bli_acquire_mpart_t2b( BLIS_SUBPART1, idx3, bs3, &Y21, &Y21_1 );
	//------------------------------------//

	bli_obj_set_struc( BLIS_HERMITIAN, A_11C_1 );
	bli_obj_set_uplo( BLIS_LOWER, A_11C_1 );
	bli_packm_init_pack( TRUE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_ROW_PANELS, 
			     BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			     BLIS_BUFFER_FOR_A_BLOCK,
			     gemm_mr, gemm_kr, 
			     &A_11C_1, &packed_A_blk );
	bli_packm_blk_var2( &BLIS_ONE, &A_11C_1, &packed_A_blk );
	bli_gemm_ker_var2( &BLIS_ONE, &packed_A_blk, &packed_B_pan, 
			   &BLIS_ZERO, &Y21_1, (gemm_t*)NULL );

	//------------------------------------//

	//****
      }

      //------------------------------------//

      //****
      bli_obj_induce_trans( Y21 );
      bli_obj_induce_trans( A_21 );
      bli_axpym( &BLIS_ONE_HALF, &Y21,
		 &A_21);
      dimLen2 = bli_obj_width_after_trans( L_21 );
      for ( idx2 = 0; idx2 < dimLen2; idx2 += bs2 ) {
	bs2 = bli_determine_blocksize_f( idx2, dimLen2, &L_21, gemm_kc );
	dim_t idx3, dimLen3, bs3;
	//****
	obj_t L_21_1;
	bli_acquire_mpart_l2r( BLIS_SUBPART1, idx2, bs2, &L_21, &L_21_1 );
	obj_t A_21_1;
	bli_acquire_mpart_l2r( BLIS_SUBPART1, idx2, bs2, &A_21, &A_21_1 );
	obj_t L_11_11;
	bli_acquire_mpart_tl2br( BLIS_SUBPART11, idx2, bs2, &L_11, &L_11_11 );
	obj_t L_11_21;
	bli_acquire_mpart_tl2br( BLIS_SUBPART21, idx2, bs2, &L_11, &L_11_21 );
	obj_t A_10_1;
	bli_acquire_mpart_t2b( BLIS_SUBPART1, idx2, bs2, &A_10, &A_10_1 );
	obj_t A_10_2;
	bli_acquire_mpart_t2b( BLIS_SUBPART2, idx2, bs2, &A_10, &A_10_2 );
	//------------------------------------//

	obj_t L_21_1H;
	bli_obj_alias_with_trans( BLIS_CONJ_TRANSPOSE, L_21_1, L_21_1H);
	obj_t A_21_1H;
	bli_obj_alias_with_trans( BLIS_CONJ_TRANSPOSE, A_21_1, A_21_1H);
	//****
	//------------------------------------//

	bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_COL_PANELS, 
			     BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			     BLIS_BUFFER_FOR_B_PANEL,
			     gemm_mr, gemm_nr, 
			     &A_10_1, &packed_B_pan );
	bli_packm_blk_var2( &BLIS_ONE, &A_10_1, &packed_B_pan );
	dim_t idx4, dimLen4, bs4;
	dimLen4 = bli_obj_length_after_trans( L_11_11 );
	for ( idx4 = 0; idx4 < dimLen4; idx4 += bs4 ) {
	  bs4 = bli_determine_blocksize_f( idx4, dimLen4, &L_11_11, gemm_mc );
	  dim_t idx5, dimLen5, bs5;
	  //****
	  obj_t L_11_11_1;
	  bli_acquire_mpart_t2b( BLIS_SUBPART1, idx4, bs4, &L_11_11, &L_11_11_1 );
	  obj_t A_10_1_1;
	  bli_acquire_mpart_t2b( BLIS_SUBPART1, idx4, bs4, &A_10_1, &A_10_1_1 );
	  //------------------------------------//

	  bli_obj_set_struc( BLIS_TRIANGULAR, L_11_11_1 );
	  bli_obj_set_uplo( BLIS_LOWER, L_11_11_1 );
	  bli_packm_init_pack( TRUE, BLIS_INVERT_DIAG, BLIS_PACKED_ROW_PANELS, 
			       BLIS_PACK_REV_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			       BLIS_BUFFER_FOR_A_BLOCK,
			       gemm_mr, gemm_mr, 
			       &L_11_11_1, &packed_A_blk );
	  bli_packm_blk_var3( &BLIS_ONE, &L_11_11_1, &packed_A_blk );
	  bli_trsm_l_ker_var2( &BLIS_ONE, &packed_A_blk, &packed_B_pan, 
			       &BLIS_ZERO, &A_10_1_1, (trsm_t*)NULL );

	  //------------------------------------//

	  //****
	}
	dimLen4 = bli_obj_length_after_trans( A_10_2 );
	for ( idx4 = 0; idx4 < dimLen4; idx4 += bs4 ) {
	  bs4 = bli_determine_blocksize_f( idx4, dimLen4, &A_10_2, gemm_mc );
	  dim_t idx5, dimLen5, bs5;
	  //****
	  obj_t L_11_21_1;
	  bli_acquire_mpart_t2b( BLIS_SUBPART1, idx4, bs4, &L_11_21, &L_11_21_1 );
	  obj_t A_10_2_1;
	  bli_acquire_mpart_t2b( BLIS_SUBPART1, idx4, bs4, &A_10_2, &A_10_2_1 );
	  //------------------------------------//

	  bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_ROW_PANELS, 
			       BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			       BLIS_BUFFER_FOR_A_BLOCK,
			       gemm_mr, gemm_kr, 
			       &L_11_21_1, &packed_A_blk );
	  bli_packm_blk_var2( &BLIS_ONE, &L_11_21_1, &packed_A_blk );
	  bli_gemm_ker_var2( &BLIS_MINUS_ONE, &packed_A_blk, &packed_B_pan, 
			     &BLIS_ONE, &A_10_2_1, (gemm_t*)NULL );

	  //------------------------------------//

	  //****
	}
	dimLen4 = bli_obj_length_after_trans( A_20 );
	for ( idx4 = 0; idx4 < dimLen4; idx4 += bs4 ) {
	  bs4 = bli_determine_blocksize_f( idx4, dimLen4, &A_20, gemm_mc );
	  dim_t idx5, dimLen5, bs5;
	  //****
	  obj_t L_21_1_1;
	  bli_acquire_mpart_t2b( BLIS_SUBPART1, idx4, bs4, &L_21_1, &L_21_1_1 );
	  obj_t A_20_1;
	  bli_acquire_mpart_t2b( BLIS_SUBPART1, idx4, bs4, &A_20, &A_20_1 );
	  //------------------------------------//

	  bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_ROW_PANELS, 
			       BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			       BLIS_BUFFER_FOR_A_BLOCK,
			       gemm_mr, gemm_mr, 
			       &L_21_1_1, &packed_A_blk );
	  bli_packm_blk_var2( &BLIS_ONE, &L_21_1_1, &packed_A_blk );
	  bli_gemm_ker_var2( &BLIS_MINUS_ONE, &packed_A_blk, &packed_B_pan, 
			     &BLIS_ONE, &A_20_1, (gemm_t*)NULL );

	  //------------------------------------//

	  //****
	}

	//------------------------------------//

	//****
	//****
	//------------------------------------//

	bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_COL_PANELS, 
			     BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			     BLIS_BUFFER_FOR_B_PANEL,
			     gemm_kr, gemm_nr, 
			     &A_21_1H, &packed_B_pan );
	bli_packm_blk_var2( &BLIS_ONE, &A_21_1H, &packed_B_pan );
	dimLen4 = bli_obj_length_after_trans( A_22 );
	for ( idx4 = 0; idx4 < dimLen4; idx4 += bs4 ) {
	  bs4 = bli_determine_blocksize_b( idx4, dimLen4, &A_22, gemm_mc );
	  dim_t idx5, dimLen5, bs5;
	  //****
	  obj_t L_21_1_1;
	  bli_acquire_mpart_b2t( BLIS_SUBPART1, idx4, bs4, &L_21_1, &L_21_1_1 );
	  obj_t A_22_0;
	  bli_acquire_mpart_b2t( BLIS_SUBPART0, idx4, bs4, &A_22, &A_22_0 );
	  obj_t A_22_1;
	  bli_acquire_mpart_b2t( BLIS_SUBPART1, idx4, bs4, &A_22, &A_22_1 );
	  //------------------------------------//

	  bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_ROW_PANELS, 
			       BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			       BLIS_BUFFER_FOR_A_BLOCK,
			       gemm_mr, gemm_kr, 
			       &L_21_1_1, &packed_A_blk );
	  bli_obj_set_struc( BLIS_HERMITIAN, A_22_1 );
	  bli_obj_set_uplo( BLIS_LOWER, A_22_1 );
	  obj_t A_22_1L, packed_B_panL;
	  dim_t offL = 0;
	  dim_t nL = bli_min( bli_obj_width_after_trans( A_22_1 ), 
			      bli_obj_diag_offset_after_trans( A_22_1 ) + bs4 );
	  bli_acquire_mpart_l2r( BLIS_SUBPART1,
				 offL, nL, &A_22_1, &A_22_1L );
	  bli_acquire_mpart_l2r( BLIS_SUBPART1,
				 offL, nL, &packed_B_pan, &packed_B_panL );
	  bli_packm_blk_var2( &BLIS_ONE, &L_21_1_1, &packed_A_blk );
	  bli_herk_l_ker_var2( &BLIS_MINUS_ONE, &packed_A_blk, &packed_B_panL, 
			       &BLIS_ONE, &A_22_1L, (herk_t*)NULL );

	  //------------------------------------//

	  //****
	}

	//------------------------------------//

	//****
	//****
	//------------------------------------//

	bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_COL_PANELS, 
			     BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			     BLIS_BUFFER_FOR_B_PANEL,
			     gemm_kr, gemm_nr, 
			     &L_21_1H, &packed_B_pan );
	bli_packm_blk_var2( &BLIS_ONE, &L_21_1H, &packed_B_pan );
	dimLen4 = bli_obj_length_after_trans( A_22 );
	for ( idx4 = 0; idx4 < dimLen4; idx4 += bs4 ) {
	  bs4 = bli_determine_blocksize_b( idx4, dimLen4, &A_22, gemm_mc );
	  dim_t idx5, dimLen5, bs5;
	  //****
	  obj_t A_21_1_1;
	  bli_acquire_mpart_b2t( BLIS_SUBPART1, idx4, bs4, &A_21_1, &A_21_1_1 );
	  obj_t A_22_0;
	  bli_acquire_mpart_b2t( BLIS_SUBPART0, idx4, bs4, &A_22, &A_22_0 );
	  obj_t A_22_1;
	  bli_acquire_mpart_b2t( BLIS_SUBPART1, idx4, bs4, &A_22, &A_22_1 );
	  //------------------------------------//

	  bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_ROW_PANELS, 
			       BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			       BLIS_BUFFER_FOR_A_BLOCK,
			       gemm_mr, gemm_kr, 
			       &A_21_1_1, &packed_A_blk );
	  bli_obj_set_struc( BLIS_HERMITIAN, A_22_1 );
	  bli_obj_set_uplo( BLIS_LOWER, A_22_1 );
	  obj_t A_22_1L, packed_B_panL;
	  dim_t offL = 0;
	  dim_t nL = bli_min( bli_obj_width_after_trans( A_22_1 ), 
			      bli_obj_diag_offset_after_trans( A_22_1 ) + bs4 );
	  bli_acquire_mpart_l2r( BLIS_SUBPART1,
				 offL, nL, &A_22_1, &A_22_1L );
	  bli_acquire_mpart_l2r( BLIS_SUBPART1,
				 offL, nL, &packed_B_pan, &packed_B_panL );
	  bli_packm_blk_var2( &BLIS_ONE, &A_21_1_1, &packed_A_blk );
	  bli_herk_l_ker_var2( &BLIS_MINUS_ONE, &packed_A_blk, &packed_B_panL, 
			       &BLIS_ONE, &A_22_1L, (herk_t*)NULL );

	  //------------------------------------//

	  //****
	}

	//------------------------------------//

	//****

	//------------------------------------//

	//****
      }
      bli_axpym( &BLIS_ONE_HALF, &Y21,
		 &A_21);

      //------------------------------------//
      bli_obj_free(&Y21);
      //****
    }



  }
  else {
    dim_t idx1, dimLen1, bs1;
    dimLen1 = min( bli_obj_length_after_trans( A ), bli_obj_width_after_trans( A ) );
    for ( idx1 = 0; idx1 < dimLen1; idx1 += bs1 ) {
      bs1 = bli_determine_blocksize_f( idx1, dimLen1, &A, bs_obj );
      dim_t idx2, dimLen2, bs2;
      //****
      obj_t A_10;
      bli_acquire_mpart_tl2br( BLIS_SUBPART10, idx1, bs1, &A, &A_10 );
      obj_t A_20;
      bli_acquire_mpart_tl2br( BLIS_SUBPART20, idx1, bs1, &A, &A_20 );
      obj_t A_11;
      bli_acquire_mpart_tl2br( BLIS_SUBPART11, idx1, bs1, &A, &A_11 );
      obj_t A_21;
      bli_acquire_mpart_tl2br( BLIS_SUBPART21, idx1, bs1, &A, &A_21 );
      obj_t A_22;
      bli_acquire_mpart_tl2br( BLIS_SUBPART22, idx1, bs1, &A, &A_22 );
      obj_t L_11;
      bli_acquire_mpart_tl2br( BLIS_SUBPART11, idx1, bs1, &L, &L_11 );
      obj_t L_21;
      bli_acquire_mpart_tl2br( BLIS_SUBPART21, idx1, bs1, &L, &L_21 );
      //------------------------------------//


	obj_t Y21;
	bli_obj_create(BLIS_DOUBLE, bli_obj_length(L_21), bli_obj_width(L_21),
		   0, 0, &Y21);

      libflame_Hegst( &A_11, &L_11 );
      obj_t A_11C;
      bli_obj_alias_with_conj( BLIS_CONJUGATE, A_11, A_11C);
      obj_t L_21T;
      bli_obj_alias_with_trans( BLIS_TRANSPOSE, L_21, L_21T);
      bli_obj_induce_trans( Y21 );
      bli_obj_induce_trans( A_21 );
      bli_scalm( &BLIS_ZERO, &Y21 );
      bli_obj_set_struc( BLIS_TRIANGULAR, A_22 );
      bli_obj_set_uplo( BLIS_LOWER, A_22 );
      bli_scalm( &BLIS_ONE, &A_22 );
      bli_obj_set_struc( BLIS_GENERAL, A_22 );
      bli_obj_set_uplo( BLIS_DENSE, A_22 );
      dimLen2 = bli_obj_length_after_trans( A_21 );
      for ( idx2 = 0; idx2 < dimLen2; idx2 += bs2 ) {
	bs2 = bli_determine_blocksize_b( idx2, dimLen2, &A_21, gemm_kc );
	dim_t idx3, dimLen3, bs3;
	//****
	obj_t L_11_10;
	bli_acquire_mpart_br2tl( BLIS_SUBPART10, idx2, bs2, &L_11, &L_11_10 );
	obj_t L_11_11;
	bli_acquire_mpart_br2tl( BLIS_SUBPART11, idx2, bs2, &L_11, &L_11_11 );
	obj_t A_21_0;
	bli_acquire_mpart_b2t( BLIS_SUBPART0, idx2, bs2, &A_21, &A_21_0 );
	obj_t A_21_1;
	bli_acquire_mpart_b2t( BLIS_SUBPART1, idx2, bs2, &A_21, &A_21_1 );
	//------------------------------------//

	bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_COL_PANELS, 
			     BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			     BLIS_BUFFER_FOR_B_PANEL,
			     gemm_mr, gemm_nr, 
			     &A_21_1, &packed_B_pan );
	bli_packm_blk_var2( &BLIS_ONE, &A_21_1, &packed_B_pan );
	dimLen3 = bli_obj_length_after_trans( L_11_11 );
	for ( idx3 = 0; idx3 < dimLen3; idx3 += bs3 ) {
	  bs3 = bli_determine_blocksize_f( idx3, dimLen3, &L_11_11, gemm_mc );
	  dim_t idx4, dimLen4, bs4;
	  //****
	  obj_t L_11_11_1;
	  bli_acquire_mpart_t2b( BLIS_SUBPART1, idx3, bs3, &L_11_11, &L_11_11_1 );
	  obj_t A_21_1_1;
	  bli_acquire_mpart_t2b( BLIS_SUBPART1, idx3, bs3, &A_21_1, &A_21_1_1 );
	  //------------------------------------//

	  obj_t L_11_11_1C;
	  bli_obj_alias_with_conj( BLIS_CONJUGATE, L_11_11_1, L_11_11_1C);
	  bli_obj_set_struc( BLIS_TRIANGULAR, L_11_11_1C );
	  bli_obj_set_uplo( BLIS_LOWER, L_11_11_1C );
	  bli_packm_init_pack( TRUE, BLIS_INVERT_DIAG, BLIS_PACKED_ROW_PANELS, 
			       BLIS_PACK_REV_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			       BLIS_BUFFER_FOR_A_BLOCK,
			       gemm_mr, gemm_mr, 
			       &L_11_11_1C, &packed_A_blk );
	  bli_packm_blk_var3( &BLIS_ONE, &L_11_11_1C, &packed_A_blk );
	  bli_trsm_l_ker_var2( &BLIS_ONE, &packed_A_blk, &packed_B_pan, 
			       &BLIS_ZERO, &A_21_1_1, (trsm_t*)NULL );

	  //------------------------------------//

	  //****
	}
	dimLen3 = bli_obj_length_after_trans( A_21_0 );
	for ( idx3 = 0; idx3 < dimLen3; idx3 += bs3 ) {
	  bs3 = bli_determine_blocksize_f( idx3, dimLen3, &A_21_0, gemm_mc );
	  dim_t idx4, dimLen4, bs4;
	  //****
	  obj_t L_11_10_1;
	  bli_acquire_mpart_l2r( BLIS_SUBPART1, idx3, bs3, &L_11_10, &L_11_10_1 );
	  obj_t A_21_0_1;
	  bli_acquire_mpart_t2b( BLIS_SUBPART1, idx3, bs3, &A_21_0, &A_21_0_1 );
	  //------------------------------------//

	  obj_t L_11_10_1C;
	  bli_obj_alias_with_conj( BLIS_CONJUGATE, L_11_10_1, L_11_10_1C);
	  bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_ROW_PANELS, 
			       BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			       BLIS_BUFFER_FOR_A_BLOCK,
			       gemm_mr, gemm_kr, 
			       &L_11_10_1C, &packed_A_blk );
	  bli_packm_blk_var2( &BLIS_ONE, &L_11_10_1C, &packed_A_blk );
	  bli_gemm_ker_var2( &BLIS_MINUS_ONE, &packed_A_blk, &packed_B_pan, 
			     &BLIS_ONE, &A_21_0_1, (gemm_t*)NULL );

	  //------------------------------------//

	  //****
	}

	//------------------------------------//

	//****
      }
      bli_obj_induce_trans( A_21 );
      dimLen2 = bli_obj_length_after_trans( Y21 );
      for ( idx2 = 0; idx2 < dimLen2; idx2 += bs2 ) {
	bs2 = bli_determine_blocksize_f( idx2, dimLen2, &Y21, gemm_kc );
	dim_t idx3, dimLen3, bs3;
	//****
	obj_t A_11C_10;
	bli_acquire_mpart_tl2br( BLIS_SUBPART10, idx2, bs2, &A_11C, &A_11C_10 );
	obj_t A_11C_11;
	bli_acquire_mpart_tl2br( BLIS_SUBPART11, idx2, bs2, &A_11C, &A_11C_11 );
	obj_t A_11C_21;
	bli_acquire_mpart_tl2br( BLIS_SUBPART21, idx2, bs2, &A_11C, &A_11C_21 );
	obj_t L_21T_1;
	bli_acquire_mpart_t2b( BLIS_SUBPART1, idx2, bs2, &L_21T, &L_21T_1 );
	obj_t Y21_0;
	bli_acquire_mpart_t2b( BLIS_SUBPART0, idx2, bs2, &Y21, &Y21_0 );
	obj_t Y21_1;
	bli_acquire_mpart_t2b( BLIS_SUBPART1, idx2, bs2, &Y21, &Y21_1 );
	obj_t Y21_2;
	bli_acquire_mpart_t2b( BLIS_SUBPART2, idx2, bs2, &Y21, &Y21_2 );
	//------------------------------------//

	bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_COL_PANELS, 
			     BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			     BLIS_BUFFER_FOR_B_PANEL,
			     gemm_kr, gemm_nr, 
			     &L_21T_1, &packed_B_pan );
	bli_packm_blk_var2( &BLIS_ONE, &L_21T_1, &packed_B_pan );
	dimLen3 = bli_obj_length_after_trans( Y21_0 );
	for ( idx3 = 0; idx3 < dimLen3; idx3 += bs3 ) {
	  bs3 = bli_determine_blocksize_f( idx3, dimLen3, &Y21_0, gemm_mc );
	  dim_t idx4, dimLen4, bs4;
	  //****
	  obj_t A_11C_10_1;
	  bli_acquire_mpart_l2r( BLIS_SUBPART1, idx3, bs3, &A_11C_10, &A_11C_10_1 );
	  obj_t Y21_0_1;
	  bli_acquire_mpart_t2b( BLIS_SUBPART1, idx3, bs3, &Y21_0, &Y21_0_1 );
	  //------------------------------------//

	  obj_t A_11C_10_1H;
	  bli_obj_alias_with_trans( BLIS_CONJ_TRANSPOSE, A_11C_10_1, A_11C_10_1H);
	  bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_ROW_PANELS, 
			       BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			       BLIS_BUFFER_FOR_A_BLOCK,
			       gemm_mr, gemm_kr, 
			       &A_11C_10_1H, &packed_A_blk );
	  bli_packm_blk_var2( &BLIS_ONE, &A_11C_10_1H, &packed_A_blk );
	  bli_gemm_ker_var2( &BLIS_ONE, &packed_A_blk, &packed_B_pan, 
			     &BLIS_ONE, &Y21_0_1, (gemm_t*)NULL );

	  //------------------------------------//

	  //****
	}
	dimLen3 = bli_obj_length_after_trans( Y21_1 );
	for ( idx3 = 0; idx3 < dimLen3; idx3 += bs3 ) {
	  bs3 = bli_determine_blocksize_f( idx3, dimLen3, &Y21_1, gemm_mc );
	  dim_t idx4, dimLen4, bs4;
	  //****
	  obj_t A_11C_11_1;
	  bli_acquire_mpart_t2b( BLIS_SUBPART1, idx3, bs3, &A_11C_11, &A_11C_11_1 );
	  obj_t Y21_1_1;
	  bli_acquire_mpart_t2b( BLIS_SUBPART1, idx3, bs3, &Y21_1, &Y21_1_1 );
	  //------------------------------------//

	  bli_obj_set_struc( BLIS_HERMITIAN, A_11C_11_1 );
	  bli_obj_set_uplo( BLIS_LOWER, A_11C_11_1 );
	  bli_packm_init_pack( TRUE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_ROW_PANELS, 
			       BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			       BLIS_BUFFER_FOR_A_BLOCK,
			       gemm_mr, gemm_kr, 
			       &A_11C_11_1, &packed_A_blk );
	  bli_packm_blk_var2( &BLIS_ONE, &A_11C_11_1, &packed_A_blk );
	  bli_gemm_ker_var2( &BLIS_ONE, &packed_A_blk, &packed_B_pan, 
			     &BLIS_ONE, &Y21_1_1, (gemm_t*)NULL );

	  //------------------------------------//

	  //****
	}
	dimLen3 = bli_obj_length_after_trans( Y21_2 );
	for ( idx3 = 0; idx3 < dimLen3; idx3 += bs3 ) {
	  bs3 = bli_determine_blocksize_f( idx3, dimLen3, &Y21_2, gemm_mc );
	  dim_t idx4, dimLen4, bs4;
	  //****
	  obj_t A_11C_21_1;
	  bli_acquire_mpart_t2b( BLIS_SUBPART1, idx3, bs3, &A_11C_21, &A_11C_21_1 );
	  obj_t Y21_2_1;
	  bli_acquire_mpart_t2b( BLIS_SUBPART1, idx3, bs3, &Y21_2, &Y21_2_1 );
	  //------------------------------------//

	  bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_ROW_PANELS, 
			       BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			       BLIS_BUFFER_FOR_A_BLOCK,
			       gemm_mr, gemm_kr, 
			       &A_11C_21_1, &packed_A_blk );
	  bli_packm_blk_var2( &BLIS_ONE, &A_11C_21_1, &packed_A_blk );
	  bli_gemm_ker_var2( &BLIS_ONE, &packed_A_blk, &packed_B_pan, 
			     &BLIS_ONE, &Y21_2_1, (gemm_t*)NULL );

	  //------------------------------------//

	  //****
	}

	//------------------------------------//

	//****
      }
      bli_obj_induce_trans( Y21 );
      bli_axpym( &BLIS_ONE_HALF, &Y21,
		 &A_21);
      dimLen2 = bli_obj_length_after_trans( A_10 );
      for ( idx2 = 0; idx2 < dimLen2; idx2 += bs2 ) {
	bs2 = bli_determine_blocksize_f( idx2, dimLen2, &A_10, gemm_kc );
	dim_t idx3, dimLen3, bs3;
	//****
	obj_t L_11_11;
	bli_acquire_mpart_tl2br( BLIS_SUBPART11, idx2, bs2, &L_11, &L_11_11 );
	obj_t L_11_21;
	bli_acquire_mpart_tl2br( BLIS_SUBPART21, idx2, bs2, &L_11, &L_11_21 );
	obj_t A_10_1;
	bli_acquire_mpart_t2b( BLIS_SUBPART1, idx2, bs2, &A_10, &A_10_1 );
	obj_t A_10_2;
	bli_acquire_mpart_t2b( BLIS_SUBPART2, idx2, bs2, &A_10, &A_10_2 );
	obj_t L_21_1;
	bli_acquire_mpart_l2r( BLIS_SUBPART1, idx2, bs2, &L_21, &L_21_1 );
	obj_t A_21_1;
	bli_acquire_mpart_l2r( BLIS_SUBPART1, idx2, bs2, &A_21, &A_21_1 );
	//------------------------------------//

	obj_t L_21_1H;
	bli_obj_alias_with_trans( BLIS_CONJ_TRANSPOSE, L_21_1, L_21_1H);
	obj_t A_21_1H;
	bli_obj_alias_with_trans( BLIS_CONJ_TRANSPOSE, A_21_1, A_21_1H);
	//****
	//------------------------------------//

	bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_COL_PANELS, 
			     BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			     BLIS_BUFFER_FOR_B_PANEL,
			     gemm_mr, gemm_nr, 
			     &A_10_1, &packed_B_pan );
	bli_packm_blk_var2( &BLIS_ONE, &A_10_1, &packed_B_pan );
  dim_t idx4, dimLen4, bs4;
	dimLen4 = bli_obj_length_after_trans( L_11_11 );
	for ( idx4 = 0; idx4 < dimLen4; idx4 += bs4 ) {
	  bs4 = bli_determine_blocksize_f( idx4, dimLen4, &L_11_11, gemm_mc );
	  dim_t idx5, dimLen5, bs5;
	  //****
	  obj_t L_11_11_1;
	  bli_acquire_mpart_t2b( BLIS_SUBPART1, idx4, bs4, &L_11_11, &L_11_11_1 );
	  obj_t A_10_1_1;
	  bli_acquire_mpart_t2b( BLIS_SUBPART1, idx4, bs4, &A_10_1, &A_10_1_1 );
	  //------------------------------------//

	  bli_obj_set_struc( BLIS_TRIANGULAR, L_11_11_1 );
	  bli_obj_set_uplo( BLIS_LOWER, L_11_11_1 );
	  bli_packm_init_pack( TRUE, BLIS_INVERT_DIAG, BLIS_PACKED_ROW_PANELS, 
			       BLIS_PACK_REV_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			       BLIS_BUFFER_FOR_A_BLOCK,
			       gemm_mr, gemm_mr, 
			       &L_11_11_1, &packed_A_blk );
	  bli_packm_blk_var3( &BLIS_ONE, &L_11_11_1, &packed_A_blk );
	  bli_trsm_l_ker_var2( &BLIS_ONE, &packed_A_blk, &packed_B_pan, 
			       &BLIS_ZERO, &A_10_1_1, (trsm_t*)NULL );

	  //------------------------------------//

	  //****
	}
	dimLen4 = bli_obj_length_after_trans( A_10_2 );
	for ( idx4 = 0; idx4 < dimLen4; idx4 += bs4 ) {
	  bs4 = bli_determine_blocksize_f( idx4, dimLen4, &A_10_2, gemm_mc );
	  dim_t idx5, dimLen5, bs5;
	  //****
	  obj_t L_11_21_1;
	  bli_acquire_mpart_t2b( BLIS_SUBPART1, idx4, bs4, &L_11_21, &L_11_21_1 );
	  obj_t A_10_2_1;
	  bli_acquire_mpart_t2b( BLIS_SUBPART1, idx4, bs4, &A_10_2, &A_10_2_1 );
	  //------------------------------------//

	  bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_ROW_PANELS, 
			       BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			       BLIS_BUFFER_FOR_A_BLOCK,
			       gemm_mr, gemm_mr, 
			       &L_11_21_1, &packed_A_blk );
	  bli_packm_blk_var2( &BLIS_ONE, &L_11_21_1, &packed_A_blk );
	  bli_gemm_ker_var2( &BLIS_MINUS_ONE, &packed_A_blk, &packed_B_pan, 
			     &BLIS_ONE, &A_10_2_1, (gemm_t*)NULL );

	  //------------------------------------//

	  //****
	}
	dimLen4 = bli_obj_length_after_trans( A_20 );
	for ( idx4 = 0; idx4 < dimLen4; idx4 += bs4 ) {
	  bs4 = bli_determine_blocksize_f( idx4, dimLen4, &A_20, gemm_mc );
	  dim_t idx5, dimLen5, bs5;
	  //****
	  obj_t L_21_1_1;
	  bli_acquire_mpart_t2b( BLIS_SUBPART1, idx4, bs4, &L_21_1, &L_21_1_1 );
	  obj_t A_20_1;
	  bli_acquire_mpart_t2b( BLIS_SUBPART1, idx4, bs4, &A_20, &A_20_1 );
	  //------------------------------------//

	  bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_ROW_PANELS, 
			       BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			       BLIS_BUFFER_FOR_A_BLOCK,
			       gemm_mr, gemm_mr, 
			       &L_21_1_1, &packed_A_blk );
	  bli_packm_blk_var2( &BLIS_ONE, &L_21_1_1, &packed_A_blk );
	  bli_gemm_ker_var2( &BLIS_MINUS_ONE, &packed_A_blk, &packed_B_pan, 
			     &BLIS_ONE, &A_20_1, (gemm_t*)NULL );

	  //------------------------------------//

	  //****
	}

	//------------------------------------//

	//****
	//****
	//------------------------------------//

	bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_COL_PANELS, 
			     BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			     BLIS_BUFFER_FOR_B_PANEL,
			     gemm_kr, gemm_nr, 
			     &A_21_1H, &packed_B_pan );
	bli_packm_blk_var2( &BLIS_ONE, &A_21_1H, &packed_B_pan );
	dimLen4 = bli_obj_length_after_trans( A_22 );
	for ( idx4 = 0; idx4 < dimLen4; idx4 += bs4 ) {
	  bs4 = bli_determine_blocksize_b( idx4, dimLen4, &A_22, gemm_mc );
	  dim_t idx5, dimLen5, bs5;
	  //****
	  obj_t L_21_1_1;
	  bli_acquire_mpart_b2t( BLIS_SUBPART1, idx4, bs4, &L_21_1, &L_21_1_1 );
	  obj_t A_22_0;
	  bli_acquire_mpart_b2t( BLIS_SUBPART0, idx4, bs4, &A_22, &A_22_0 );
	  obj_t A_22_1;
	  bli_acquire_mpart_b2t( BLIS_SUBPART1, idx4, bs4, &A_22, &A_22_1 );
	  //------------------------------------//

	  bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_ROW_PANELS, 
			       BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			       BLIS_BUFFER_FOR_A_BLOCK,
			       gemm_mr, gemm_kr, 
			       &L_21_1_1, &packed_A_blk );
	  bli_obj_set_struc( BLIS_HERMITIAN, A_22_1 );
	  bli_obj_set_uplo( BLIS_LOWER, A_22_1 );
	  obj_t A_22_1L, packed_B_panL;
	  dim_t offL = 0;
	  dim_t nL = bli_min( bli_obj_width_after_trans( A_22_1 ), 
			      bli_obj_diag_offset_after_trans( A_22_1 ) + bs4 );
	  bli_acquire_mpart_l2r( BLIS_SUBPART1,
				 offL, nL, &A_22_1, &A_22_1L );
	  bli_acquire_mpart_l2r( BLIS_SUBPART1,
				 offL, nL, &packed_B_pan, &packed_B_panL );
	  bli_packm_blk_var2( &BLIS_ONE, &L_21_1_1, &packed_A_blk );
	  bli_herk_l_ker_var2( &BLIS_MINUS_ONE, &packed_A_blk, &packed_B_panL, 
			       &BLIS_ONE, &A_22_1L, (herk_t*)NULL );

	  //------------------------------------//

	  //****
	}

	//------------------------------------//

	//****
	//****
	//------------------------------------//

	bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_COL_PANELS, 
			     BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			     BLIS_BUFFER_FOR_B_PANEL,
			     gemm_kr, gemm_nr, 
			     &L_21_1H, &packed_B_pan );
	bli_packm_blk_var2( &BLIS_ONE, &L_21_1H, &packed_B_pan );
	dimLen4 = bli_obj_length_after_trans( A_22 );
	for ( idx4 = 0; idx4 < dimLen4; idx4 += bs4 ) {
	  bs4 = bli_determine_blocksize_b( idx4, dimLen4, &A_22, gemm_mc );
	  dim_t idx5, dimLen5, bs5;
	  //****
	  obj_t A_21_1_1;
	  bli_acquire_mpart_b2t( BLIS_SUBPART1, idx4, bs4, &A_21_1, &A_21_1_1 );
	  obj_t A_22_0;
	  bli_acquire_mpart_b2t( BLIS_SUBPART0, idx4, bs4, &A_22, &A_22_0 );
	  obj_t A_22_1;
	  bli_acquire_mpart_b2t( BLIS_SUBPART1, idx4, bs4, &A_22, &A_22_1 );
	  //------------------------------------//

	  bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_ROW_PANELS, 
			       BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			       BLIS_BUFFER_FOR_A_BLOCK,
			       gemm_mr, gemm_kr, 
			       &A_21_1_1, &packed_A_blk );
	  bli_obj_set_struc( BLIS_HERMITIAN, A_22_1 );
	  bli_obj_set_uplo( BLIS_LOWER, A_22_1 );
	  obj_t A_22_1L, packed_B_panL;
	  dim_t offL = 0;
	  dim_t nL = bli_min( bli_obj_width_after_trans( A_22_1 ), 
			      bli_obj_diag_offset_after_trans( A_22_1 ) + bs4 );
	  bli_acquire_mpart_l2r( BLIS_SUBPART1,
				 offL, nL, &A_22_1, &A_22_1L );
	  bli_acquire_mpart_l2r( BLIS_SUBPART1,
				 offL, nL, &packed_B_pan, &packed_B_panL );
	  bli_packm_blk_var2( &BLIS_ONE, &A_21_1_1, &packed_A_blk );
	  bli_herk_l_ker_var2( &BLIS_MINUS_ONE, &packed_A_blk, &packed_B_panL, 
			       &BLIS_ONE, &A_22_1L, (herk_t*)NULL );

	  //------------------------------------//

	  //****
	}

	//------------------------------------//

	//****

	//------------------------------------//

	//****
      }
      bli_axpym( &BLIS_ONE_HALF, &Y21,
		 &A_21);

      //------------------------------------//
      bli_obj_free(&Y21);
      //****
    }

  }

  bli_obj_release_pack( &packed_A_blk );
  bli_obj_release_pack( &packed_B_pan );
}

void FLA_TwoSidedTrsm( obj_t *a2, 
		       obj_t *b2,
		       fla_eig_gest_t *cntl )
{
  FLA_Error err;
  FLA_Obj A, B, Y;

  err = wrap_blis_obj( a2, &A );
  if (err != FLA_SUCCESS)
    return;

  err = wrap_blis_obj( b2, &B );
  if (err != FLA_SUCCESS)
    return;

  FLA_Obj_create( FLA_DOUBLE, 	 
		  bli_obj_length(*a2), bli_obj_width(*a2),
		  0, 0, &Y );
  
#if 1
  err = FLA_Eig_gest_il( A, Y, B, cntl );
#elif 1
  err = FLA_Eig_gest(FLA_INVERSE, FLA_LOWER_TRIANGULAR, A, B);
#else
  FLA_Eig_gest_il_opt_var2( A, Y, B);
#endif

  FLA_Obj_free(&Y);

  if (err != FLA_SUCCESS) {
    double *tmp = NULL;
    *tmp = 0;
    printf("FLA_LU error\n");
  }
				  
  FLA_Obj_free_without_buffer( &A );
  FLA_Obj_free_without_buffer( &B );
}


int main( int argc, char** argv )
{
  obj_t a1, a2, b1, b2;
  obj_t a_save, b_save;
#ifdef TESTLIB
  obj_t normVal1, normVal2;
  obj_t x, y, z, normVal3;
#endif
  dim_t m, bs;
  dim_t p;
  dim_t p_begin, p_end, p_inc;
  dim_t bs_begin, bs_end, bs_inc;
  int   m_input, n_input, k_input;
  num_t dt_a, dt_b, dt_c;
  num_t dt_alpha, dt_beta;
  int   r, n_repeats;
  dim_t bsDxTBest, bsFLABest;

  double dtimeDxT;
  double dtime_saveDxT;
  double gflopsDxT, gflopsDxTBest;
  double dtimeFLA;
  double dtime_saveFLA;
  double gflopsFLA, gflopsFLABest;

  blksz_t *bs_obj;

  fla_eig_gest_t*     cntl;

  bli_init();
  FLA_Init();


  n_repeats = 1;

#if 0
#define VAL 130
  p_begin = VAL;
  p_end   = VAL;
  p_inc   = VAL;
#else
  p_begin = 200;
  p_end   = 5000;
  p_inc   = 200;
  bs_begin = 64;
  bs_end   = 256;
  bs_inc   = 16;
#endif

  m_input = -1;

  dt_a = BLIS_DOUBLE;
  dt_b = BLIS_DOUBLE;
  dt_c = BLIS_DOUBLE;
  dt_alpha = BLIS_DOUBLE;
  dt_beta = BLIS_DOUBLE;

#ifdef TESTLIB
      bli_obj_create( dt_alpha, 1, 1, 0, 0, &normVal1 );
      bli_obj_create( dt_alpha, 1, 1, 0, 0, &normVal2 );
      bli_obj_create( dt_alpha, 1, 1, 0, 0, &normVal3 );
#endif



      bli_error_checking_level_set( BLIS_NO_ERROR_CHECKING );

      for ( p = p_begin; p <= p_end; p += p_inc )
    {

      if ( m_input < 0 ) m = p * ( dim_t )abs(m_input);
      else               m =     ( dim_t )    m_input;

      bli_obj_create( dt_a, m, m, 0, 0, &a1 );
      bli_obj_create( dt_a, m, m, 0, 0, &a2 );
      bli_obj_create( dt_c, m, m, 0, 0, &a_save );
      bli_obj_create( dt_a, m, m, 0, 0, &b1 );
      bli_obj_create( dt_a, m, m, 0, 0, &b2 );
      bli_obj_create( dt_c, m, m, 0, 0, &b_save );
#ifdef TESTLIB
      bli_obj_create( dt_a, m, m, 0, 0, &x );
      bli_obj_create( dt_a, m, m, 0, 0, &y );
      bli_obj_create( dt_a, m, m, 0, 0, &z );
#endif
      
      bli_randm( &a1 );
      bli_obj_set_struc( BLIS_TRIANGULAR, b1 );
      bli_obj_set_uplo( BLIS_LOWER, b1 );
      bli_randm( &b1 );
      bli_obj_set_struc( BLIS_GENERAL, b1 );
      bli_obj_set_uplo( BLIS_DENSE, b1 );
      double *buf = (double*)(bli_obj_buffer_at_off(b1));
      for (dim_t loc = 0; loc < m; ++loc) {
	buf[loc*bli_obj_col_stride(b1)+loc*bli_obj_row_stride(b1)] += m;
      }

      bli_copym( &a1, &a_save );
      bli_copym( &b1, &b_save );
	
      gflopsDxTBest = 0;
      gflopsFLABest = 0;
      
      for ( bs = bs_begin; bs <= bs_end; bs += bs_inc ) {

	bs_obj = bli_blksz_obj_create(bs, 0,
				      bs, 0,
				      bs, 0,
				      bs, 0);

	fla_blocksize_t *bs_fla = FLA_Blocksize_create(bs, bs, bs, bs);

	cntl = FLA_Cntl_eig_gest_obj_create( FLA_FLAT,
					     FLA_BLOCKED_VARIANT4,
					     bs_fla,
					     fla_eig_gest_ix_cntl_leaf,
					     fla_axpy_cntl_blas,
					     fla_axpy_cntl_blas,
					     fla_gemm_cntl_blas,
					     fla_gemm_cntl_blas,
					     fla_gemm_cntl_blas,
					     fla_hemm_cntl_blas,
					     fla_her2k_cntl_blas,
					     fla_trmm_cntl_blas,
					     fla_trmm_cntl_blas,
					     fla_trsm_cntl_blas,
					     fla_trsm_cntl_blas );

	dtime_saveDxT = 1.0e9;
	dtime_saveFLA = 1.0e9;

	for ( r = 0; r < n_repeats; ++r )
	  {
	    bli_copym( &a_save, &a1 );
	    bli_copym( &a_save, &a2 );
	    bli_copym( &b_save, &b1 );
	    bli_copym( &b_save, &b2 );
#ifdef TESTLIB
	    bli_randm( &x );
	    bli_copym( &x, &y );
#endif

	    dtimeDxT = bli_clock();
	    DxT_TwoSidedTrsm( a1, b1, bs_obj );
	    dtime_saveDxT = bli_clock_min_diff( dtime_saveDxT, dtimeDxT );

	    dtimeFLA = bli_clock();
	    FLA_TwoSidedTrsm( &a2, &b2, cntl);
	    dtime_saveFLA = bli_clock_min_diff( dtime_saveFLA, dtimeFLA );


	  /* Trsm( LEFT, LOWER, ADJOINT, diag, F(1), B, Y ); */
/* 	  Hemm( LEFT, LOWER, F(1), AOrig, Y, F(0), Z ); */
/* 	  Trsm( LEFT, LOWER, NORMAL, diag, F(1), B, Z ); */
/* 	  Hemm( LEFT, LOWER, F(-1), A, X, F(1), Z ); */

#ifdef TESTLIB
	  bli_obj_set_struc( BLIS_TRIANGULAR, b2 );
	  bli_obj_set_uplo( BLIS_LOWER, b2 );
	  bli_obj_set_conjtrans( BLIS_CONJ_TRANSPOSE, b2);
	  bli_trsm(BLIS_LEFT, &BLIS_ONE, &b2, &y);
	  bli_obj_set_conjtrans( BLIS_NO_TRANSPOSE, b2);

	  bli_obj_set_struc( BLIS_HERMITIAN, a_save );
	  bli_obj_set_uplo( BLIS_LOWER, a_save );
	  bli_hemm(BLIS_LEFT, &BLIS_ONE, &a_save, &y, &BLIS_ZERO, &z);

	  bli_trsm(BLIS_LEFT, &BLIS_ONE, &b2, &z);

	  bli_obj_set_struc( BLIS_HERMITIAN, a2 );
	  bli_obj_set_uplo( BLIS_LOWER, a2 );
	  bli_hemm(BLIS_LEFT, &BLIS_MINUS_ONE, &a2, &x, &BLIS_ONE, &z);
	  bli_obj_set_struc( BLIS_GENERAL, a2 );
	  bli_obj_set_uplo( BLIS_DENSE, a2 );

	  bli_obj_set_struc( BLIS_GENERAL, b2 );
	  bli_obj_set_uplo( BLIS_DENSE, b2 );
	  bli_obj_set_struc( BLIS_GENERAL, a_save );
	  bli_obj_set_uplo( BLIS_DENSE, a_save );

	  bli_fnormm( &z, &normVal3 );

	  bli_axpym( &BLIS_MINUS_ONE, &a1, &a2 );
	  bli_axpym( &BLIS_MINUS_ONE, &b1, &b2 );

	  bli_fnormm( &a2, &normVal1 );
	  bli_fnormm( &b2, &normVal2 );

	  //	  bli_printm( "a1-a2 ", &a2, "%4.1e", "" );
#endif


	  }

	gflopsDxT = ( ((1.0 * m) * m) * m ) / ( dtime_saveDxT * 1.0e9 );
	if (gflopsDxT > gflopsDxTBest) {
	  gflopsDxTBest = gflopsDxT;
	  bsDxTBest = bs;
	}
	gflopsFLA = ( ((1.0 * m) * m) * m ) / ( dtime_saveFLA * 1.0e9 );
	if (gflopsFLA > gflopsFLABest) {
	  gflopsFLABest = gflopsFLA;
	  bsFLABest = bs;
	}
	
	printf( "data_TwoSidedTrsm_DxT" );
	printf( "( %2lu, %2lu, 1:4 ) = [ %4lu %4lu %10.3e  %6.3f ];\n",
		(p - p_begin + 1)/p_inc + 1, (bs - bs_begin + 1)/bs_inc + 1,
		m, bs, dtime_saveDxT, gflopsDxT );

	printf( "data_TwoSidedTrsm_FLA" );
	printf( "( %2lu, %2lu, 1:4 ) = [ %4lu %4lu %10.3e  %6.3f ];\n\n",
		(p - p_begin + 1)/p_inc + 1, (bs - bs_begin + 1)/bs_inc + 1,
		m, bs, dtime_saveFLA, gflopsFLA );

	fflush(stdout);
	
#ifdef TESTLIB
	bli_printm( "%%NORM1 ", &normVal1, "%4.1f", "" );
	bli_printm( "%%NORM2 ", &normVal2, "%4.1f", "" );
	bli_printm( "%%NORM3 ", &normVal3, "%4.1f", "" );
#endif

	fflush(stdout);

	FLA_Cntl_obj_free( cntl );
	FLA_Blocksize_free( bs_fla );
	bli_blksz_obj_free( bs_obj );
      }
	printf( "data_TwoSidedTrsm_DxT_Best" );
	printf( "( %2lu, 1:2 ) = [ %4lu  %6.3f ];  %%%% %lu\n",
		(p - p_begin + 1)/p_inc + 1,
		m, gflopsDxTBest, bsDxTBest );

	printf( "data_TwoSidedTrsm_FLA_Best" );
	printf( "( %2lu, 1:2 ) = [ %4lu  %6.3f ];  %%%% %lu\n",
		(p - p_begin + 1)/p_inc + 1,
		m, gflopsFLABest, bsFLABest);
      printf("%%%%%%%%%%\n");
      fflush(stdout);

      bli_obj_free( &a1 );
      bli_obj_free( &a2 );
      bli_obj_free( &b1 );
      bli_obj_free( &b2 );
      bli_obj_free( &a_save );
      bli_obj_free( &b_save );
#ifdef TESTLIB
      bli_obj_free( &x );
      bli_obj_free( &y );
      bli_obj_free( &z );
#endif
    }

  FLA_Finalize();
  bli_finalize();

  return 0;
}

