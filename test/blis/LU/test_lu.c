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

extern fla_gemm_t*  fla_gemm_cntl_blas;
extern fla_trsm_t*  fla_trsm_cntl_blas;
extern fla_appiv_t* fla_appiv_cntl_leaf;

fla_lu_t*           fla_lu_piv_cntl;
fla_lu_t*           fla_lu_piv_cntl_in;
fla_lu_t*           fla_lu_piv_cntl_leaf;
fla_blocksize_t*    fla_lu_piv_var5_bsize;
fla_blocksize_t*    fla_lu_piv_var5_bsize_in;

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

void PanelLU( obj_t *A_10,
	      obj_t *A_20,
	      obj_t *A_11,
	      obj_t *A_21,
	      obj_t *A_12,
	      obj_t *A_22,
	      obj_t *P1_bl )
{
  obj_t AB0_bl, AB1_bl, AB2_bl;
  FLA_Obj AB0, AB1, AB2, p1;
  
  wrap_blis_obj( P1_bl, &p1 );

  //  printf("merge 1\n");
  //  fflush(stdout);
  blis_obj_merge_vert( A_10, 
		       A_20, &AB0_bl );
  wrap_blis_obj( &AB0_bl, &AB0 );

  //  printf("merge 2\n");
  //  fflush(stdout);
  blis_obj_merge_vert( A_11, 
		       A_21, &AB1_bl );
  wrap_blis_obj( &AB1_bl, &AB1 );

  //  printf("merge 3\n");
  //  fflush(stdout);
  blis_obj_merge_vert( A_12, 
		       A_22, &AB2_bl );
  wrap_blis_obj( &AB2_bl, &AB2 );

  fla_appiv_t *pivCntl = FLA_Cntl_appiv_obj_create( FLA_FLAT, FLA_UNB_OPT_VARIANT1, NULL, NULL );

  fla_lu_t *luCntl = FLA_Cntl_lu_obj_create( FLA_FLAT, FLA_UNB_OPT_VARIANT5, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL );

  // AB1, p1 = LU_piv( AB1 )                                                                   
  FLA_LU_piv_internal( AB1, p1, 
		       fla_lu_piv_cntl_in );
		       //		       luCntl );

  // Apply computed pivots to AB0                                                              
  FLA_Apply_pivots_internal( FLA_LEFT, FLA_NO_TRANSPOSE, p1, AB0,
			     fla_appiv_cntl_leaf );
			     //			     pivCntl );

  // Apply computed pivots to AB2                                                              
  FLA_Apply_pivots_internal( FLA_LEFT, FLA_NO_TRANSPOSE, p1, AB2,
			     fla_appiv_cntl_leaf );
			     //			     pivCntl );

  FLA_Obj_free_without_buffer(&p1);
  FLA_Obj_free_without_buffer(&AB0);
  FLA_Obj_free_without_buffer(&AB1);
  FLA_Obj_free_without_buffer(&AB2);
  FLA_Cntl_obj_free( pivCntl );
  FLA_Cntl_obj_free( luCntl );
}

void DxT_LU( obj_t A,
	     obj_t P,
	     blksz_t *bs_obj)
{
  obj_t packed_A_blk;
  obj_t packed_B_pan;

  bli_obj_init_pack( &packed_A_blk );
  bli_obj_init_pack( &packed_B_pan );

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
    obj_t A_12;
    bli_acquire_mpart_tl2br( BLIS_SUBPART12, idx1, bs1, &A, &A_12 );
    obj_t A_22;
    bli_acquire_mpart_tl2br( BLIS_SUBPART22, idx1, bs1, &A, &A_22 );
    obj_t P_1;
    bli_acquire_mpart_t2b( BLIS_SUBPART1, idx1, bs1, &P, &P_1 );
    //------------------------------------//

    PanelLU( &A_10, &A_20,
	     &A_11, &A_21,
	     &A_12, &A_22,
	     &P_1 );
    dimLen2 = bli_obj_length_after_trans( A_12 );
    for ( idx2 = 0; idx2 < dimLen2; idx2 += bs2 ) {
      bs2 = bli_determine_blocksize_f( idx2, dimLen2, &A_12, gemm_kc );
      dim_t idx3, dimLen3, bs3;
      //****
      obj_t A_11_11;
      bli_acquire_mpart_tl2br( BLIS_SUBPART11, idx2, bs2, &A_11, &A_11_11 );
      obj_t A_11_21;
      bli_acquire_mpart_tl2br( BLIS_SUBPART21, idx2, bs2, &A_11, &A_11_21 );
      obj_t A_12_1;
      bli_acquire_mpart_t2b( BLIS_SUBPART1, idx2, bs2, &A_12, &A_12_1 );
      obj_t A_12_2;
      bli_acquire_mpart_t2b( BLIS_SUBPART2, idx2, bs2, &A_12, &A_12_2 );
      obj_t A_21_1;
      bli_acquire_mpart_l2r( BLIS_SUBPART1, idx2, bs2, &A_21, &A_21_1 );
      //------------------------------------//

      bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_COL_PANELS, 
			   BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			   BLIS_BUFFER_FOR_B_PANEL,
			   gemm_mr, gemm_nr, 
			   &A_12_1, &packed_B_pan );
      bli_packm_blk_var2( &BLIS_ONE, &A_12_1, &packed_B_pan );
      dimLen3 = bli_obj_length_after_trans( A_11_11 );
      for ( idx3 = 0; idx3 < dimLen3; idx3 += bs3 ) {
	bs3 = bli_determine_blocksize_f( idx3, dimLen3, &A_11_11, gemm_mc );
	dim_t idx4, dimLen4, bs4;
	//****
	obj_t A_11_11_1;
	bli_acquire_mpart_t2b( BLIS_SUBPART1, idx3, bs3, &A_11_11, &A_11_11_1 );
	obj_t A_12_1_1;
	bli_acquire_mpart_t2b( BLIS_SUBPART1, idx3, bs3, &A_12_1, &A_12_1_1 );
	//------------------------------------//

	bli_obj_set_struc( BLIS_TRIANGULAR, A_11_11_1 );
	bli_obj_set_uplo( BLIS_LOWER, A_11_11_1 );
	bli_obj_set_diag( BLIS_UNIT_DIAG, A_11_11_1 );
	bli_packm_init_pack( TRUE, BLIS_INVERT_DIAG, BLIS_PACKED_ROW_PANELS, 
			     BLIS_PACK_REV_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			     BLIS_BUFFER_FOR_A_BLOCK,
			     gemm_mr, gemm_mr, 
			     &A_11_11_1, &packed_A_blk );
	bli_packm_blk_var3( &BLIS_ONE, &A_11_11_1, &packed_A_blk );
	bli_trsm_l_ker_var2( &BLIS_ONE, &packed_A_blk, &packed_B_pan, 
			     &BLIS_ZERO, &A_12_1_1, (trsm_t*)NULL );

	//------------------------------------//

	//****
      }
      dimLen3 = bli_obj_length_after_trans( A_12_2 );
      for ( idx3 = 0; idx3 < dimLen3; idx3 += bs3 ) {
	bs3 = bli_determine_blocksize_f( idx3, dimLen3, &A_12_2, gemm_mc );
	dim_t idx4, dimLen4, bs4;
	//****
	obj_t A_11_21_1;
	bli_acquire_mpart_t2b( BLIS_SUBPART1, idx3, bs3, &A_11_21, &A_11_21_1 );
	obj_t A_12_2_1;
	bli_acquire_mpart_t2b( BLIS_SUBPART1, idx3, bs3, &A_12_2, &A_12_2_1 );
	//------------------------------------//

	bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_ROW_PANELS, 
			     BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			     BLIS_BUFFER_FOR_A_BLOCK,
			     gemm_mr, gemm_kr, 
			     &A_11_21_1, &packed_A_blk );
	bli_packm_blk_var2( &BLIS_ONE, &A_11_21_1, &packed_A_blk );
	bli_gemm_ker_var2( &BLIS_MINUS_ONE, &packed_A_blk, &packed_B_pan, 
			   &BLIS_ONE, &A_12_2_1, (gemm_t*)NULL );

	//------------------------------------//

	//****
      }
      dimLen3 = bli_obj_length_after_trans( A_22 );
      for ( idx3 = 0; idx3 < dimLen3; idx3 += bs3 ) {
	bs3 = bli_determine_blocksize_f( idx3, dimLen3, &A_22, gemm_mc );
	dim_t idx4, dimLen4, bs4;
	//****
	obj_t A_21_1_1;
	bli_acquire_mpart_t2b( BLIS_SUBPART1, idx3, bs3, &A_21_1, &A_21_1_1 );
	obj_t A_22_1;
	bli_acquire_mpart_t2b( BLIS_SUBPART1, idx3, bs3, &A_22, &A_22_1 );
	//------------------------------------//

	bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_ROW_PANELS, 
			     BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			     BLIS_BUFFER_FOR_A_BLOCK,
			     gemm_mr, gemm_kr, 
			     &A_21_1_1, &packed_A_blk );
	bli_packm_blk_var2( &BLIS_ONE, &A_21_1_1, &packed_A_blk );
	bli_gemm_ker_var2( &BLIS_MINUS_ONE, &packed_A_blk, &packed_B_pan, 
			   &BLIS_ONE, &A_22_1, (gemm_t*)NULL );

	//------------------------------------//

	//****
      }

      //------------------------------------//

      //****
    }

    //------------------------------------//

    //****
  }



  bli_obj_release_pack( &packed_A_blk );
  bli_obj_release_pack( &packed_B_pan );
}

void FLA_LU( obj_t *a2, 
	     obj_t *pVec,
	     fla_lu_t* fla_lu_piv_cntl2)
{
  FLA_Error err;
  FLA_Obj A, p;

  err = wrap_blis_obj( a2, &A );
  if (err != FLA_SUCCESS)
    return;

  err = wrap_blis_obj( pVec, &p );
  if (err != FLA_SUCCESS)
    return;

  FLA_LU_piv_internal( A, p, fla_lu_piv_cntl2 );

  if (err != FLA_SUCCESS) {
    double *tmp = NULL;
    *tmp = 0;
    printf("FLA_LU error\n");
  }

  FLA_Obj_free_without_buffer( &A );
  FLA_Obj_free_without_buffer( &p );
}


int main( int argc, char** argv )
{
  obj_t pVec, a1, a2;
  obj_t a_save;
  obj_t negOne, normVal;
  dim_t m, bs;
  dim_t p;
  dim_t p_begin, p_end, p_inc;
  dim_t bs_begin, bs_end, bs_inc;
  int   m_input, n_input, k_input;
  num_t dt_a, dt_b, dt_c;
  num_t dt_alpha, dt_beta;
  int   r, n_repeats;

  double dtimeDxT;
  double dtime_saveDxT;
  double gflopsDxT;
  double gflopsDxTBest;
  double dtimeFLA;
  double dtime_saveFLA;
  double gflopsFLA;
  double gflopsFLABest;
  dim_t bsDxTBest, bsFLABest;

  fla_lu_t*           fla_lu_piv_cntl2;

  bli_error_checking_level_set( BLIS_NO_ERROR_CHECKING );

  bli_obj_create( BLIS_DOUBLE, 1, 1, 0, 0, &negOne );
  bli_setsc(  -1.0, 0.0, &negOne );		
	
  blksz_t *bs_obj;
		
  bli_init();
  FLA_Init();

  n_repeats = 3;

  p_begin = 200;
  p_end   = 5000;
  p_inc   = 200;
  bs_begin = 64;
  bs_end   = 256;
  bs_inc   = 16;

  m_input = -1;

  dt_a = BLIS_DOUBLE;
  dt_b = BLIS_DOUBLE;
  dt_c = BLIS_DOUBLE;
  dt_alpha = BLIS_DOUBLE;
  dt_beta = BLIS_DOUBLE;

  for ( p = p_begin; p <= p_end; p += p_inc )
    {

      if ( m_input < 0 ) m = p * ( dim_t )abs(m_input);
      else               m =     ( dim_t )    m_input;

      bli_obj_create( dt_a, m, m, 0, 0, &a1 );
      bli_obj_create( dt_a, m, m, 0, 0, &a2 );
      bli_obj_create( BLIS_INT, m, 1, 0, 0, &pVec );
      bli_obj_create( dt_c, m, m, 0, 0, &a_save );

      bli_randm( &a1 );
      //      bli_randm( &pVec );

      bli_copym( &a1, &a_save );
	
      gflopsDxTBest = 0;
      gflopsFLABest = 0;

      for ( bs = bs_begin; bs <= bs_end; bs += bs_inc ) {

	fla_blocksize_t *bs_fla = FLA_Blocksize_create(bs, bs, bs, bs);


	fla_lu_piv_cntl2       = FLA_Cntl_lu_obj_create( FLA_FLAT,
                                                         FLA_BLOCKED_VARIANT5,
							 bs_fla,
                                                         fla_lu_piv_cntl_in,
                                                         fla_gemm_cntl_blas,
                                                         fla_gemm_cntl_blas,
                                                         fla_gemm_cntl_blas,
                                                         fla_trsm_cntl_blas,
                                                         fla_trsm_cntl_blas,
                                                         fla_appiv_cntl_leaf,
                                                         fla_appiv_cntl_leaf );



	bs_obj = bli_blksz_obj_create(bs, 0,
				      bs, 0,
				      bs, 0,
				      bs, 0);

	dtime_saveDxT = 1.0e9;
	dtime_saveFLA = 1.0e9;

	for ( r = 0; r < n_repeats; ++r )
	  {
	    bli_copym( &a_save, &a1 );
	    bli_copym( &a_save, &a2 );

	    dtimeDxT = bli_clock();
	    DxT_LU( a1, pVec, bs_obj );
	    dtime_saveDxT = bli_clock_min_diff( dtime_saveDxT, dtimeDxT );

	    dtimeFLA = bli_clock();
	    FLA_LU( &a2, &pVec, fla_lu_piv_cntl2 );
	    dtime_saveFLA = bli_clock_min_diff( dtime_saveFLA, dtimeFLA );

	    bli_axpym( &negOne, &a1, &a2 );

	    bli_obj_create( dt_alpha, 1, 1, 0, 0, &normVal );

	    bli_fnormm( &a2, &normVal );
	  }

	gflopsDxT = ( (2.0 / 3.0) * m * m * m ) / ( dtime_saveDxT * 1.0e9 );
	if (gflopsDxT > gflopsDxTBest) {
	  gflopsDxTBest = gflopsDxT;
	  bsDxTBest = bs;
	}
	gflopsFLA = ( (2.0 / 3.0) * m * m * m ) / ( dtime_saveFLA * 1.0e9 );
	if (gflopsFLA > gflopsFLABest) {
	  gflopsFLABest = gflopsFLA;
	  bsFLABest = bs;
	}

	printf( "data_LU_DxT" );
	printf( "( %2lu, 1:4 ) = [ %4lu %4lu %10.3e  %6.3f ];\n",
		(p - p_begin + 1)/p_inc + 1, m, bs, dtime_saveDxT, gflopsDxT );
	printf( "data_LU_FLA" );
	printf( "( %2lu, 1:4 ) = [ %4lu %4lu %10.3e  %6.3f ];\n",
		(p - p_begin + 1)/p_inc + 1, m, bs, dtime_saveFLA, gflopsFLA );

#ifdef TESTLIB
	bli_printm( "%%%%NORM", &normVal, "%4.1f", "" );
#endif

	FLA_Cntl_obj_free( fla_lu_piv_cntl2 );
	FLA_Blocksize_free( bs_fla );
	bli_blksz_obj_free( bs_obj );
      }

	printf( "data_LU_DxT_Best" );
	printf( "( %2lu, 1:2 ) = [ %4lu  %6.3f ];  %%%% %lu\n",
		(p - p_begin + 1)/p_inc + 1,
		m, gflopsDxTBest, bsDxTBest );

	printf( "data_LU_FLA_Best" );
	printf( "( %2lu, 1:2 ) = [ %4lu  %6.3f ];  %%%% %lu\n",
		(p - p_begin + 1)/p_inc + 1,
		m, gflopsFLABest, bsFLABest);
      printf("%%%%%%%%%%\n\n\n");



      bli_obj_free( &a1 );
      bli_obj_free( &a2 );
      bli_obj_free( &pVec );
      bli_obj_free( &a_save );
    }

  FLA_Finalize();
  bli_finalize();

  return 0;
}

