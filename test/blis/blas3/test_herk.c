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


void DxT_HerkNL( obj_t *alpha,
		 obj_t *A,
		 obj_t *beta,
		 obj_t *C )
{
  obj_t packed_A_blk;
  obj_t packed_B_pan;
  bli_obj_init_pack( &packed_A_blk );
  bli_obj_init_pack( &packed_B_pan );


  dim_t idx1, dimLen1, bs1;
  obj_t AT;
  bli_obj_alias_with_trans( BLIS_TRANSPOSE, *A, AT);
  dimLen1 = bli_obj_width_after_trans( *C );
  for ( idx1 = 0; idx1 < dimLen1; idx1 += bs1 ) {
    bs1 = bli_determine_blocksize_f( idx1, dimLen1, C, gemm_nc );
    dim_t idx2, dimLen2, bs2;
    //****
    obj_t AT_1;
    bli_acquire_mpart_l2r( BLIS_SUBPART1, idx1, bs1, &AT, &AT_1 );
    obj_t C_0;
    bli_acquire_mpart_l2r( BLIS_SUBPART0, idx1, bs1, C, &C_0 );
    obj_t C_1;
    bli_acquire_mpart_l2r( BLIS_SUBPART1, idx1, bs1, C, &C_1 );
    //------------------------------------//

    obj_t C_1B, AB;
    printf("off %u\n", bli_obj_diag_offset_after_trans(C_1));
    dim_t offB = bli_max( 0, -bli_obj_diag_offset_after_trans( C_1 ) );
    dim_t mB = bli_obj_length_after_trans( C_1 ) - offB;
    printf("offB %u, mB %u\n",offB, mB);
    printf("A %u x %u\n", bli_obj_length(*A),bli_obj_width(*A));
    printf("C_1 %u x %u\n", bli_obj_length(C_1),bli_obj_width(C_1));
    bli_acquire_mpart_t2b( BLIS_SUBPART1,
			   offB, mB, &C_1, &C_1B );
    bli_acquire_mpart_t2b( BLIS_SUBPART1,
			   offB, mB, A, &AB );
    printf("C_1B %u x %u\n", bli_obj_length(C_1B),bli_obj_width(C_1B));
    printf("AB %u x %u\n", bli_obj_length(AB),bli_obj_width(AB));
    bli_obj_set_struc( BLIS_TRIANGULAR, C_1B );
    bli_obj_set_uplo( BLIS_LOWER, C_1B );
    bli_scalm( &BLIS_ONE, &C_1B );
    bli_obj_set_struc( BLIS_GENERAL, C_1B );
    bli_obj_set_uplo( BLIS_DENSE, C_1B );
    dimLen2 = bli_obj_width_after_trans( AB );
    for ( idx2 = 0; idx2 < dimLen2; idx2 += bs2 ) {
      bs2 = bli_determine_blocksize_f( idx2, dimLen2, &AB, gemm_kc );
      dim_t idx3, dimLen3, bs3;
      //****
      obj_t AB_1;
      bli_acquire_mpart_l2r( BLIS_SUBPART1, idx2, bs2, &AB, &AB_1 );
      obj_t AT_1_1;
      bli_acquire_mpart_t2b( BLIS_SUBPART1, idx2, bs2, &AT_1, &AT_1_1 );
      //------------------------------------//

      bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_COL_PANELS, 
			   BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			   BLIS_BUFFER_FOR_B_PANEL,
			   gemm_kr, gemm_nr, 
			   &AT_1_1, &packed_B_pan );
      bli_packm_blk_var2( &BLIS_ONE, &AT_1_1, &packed_B_pan );
      dimLen3 = bli_obj_length_after_trans( C_1B );
      for ( idx3 = 0; idx3 < dimLen3; idx3 += bs3 ) {
	bs3 = bli_determine_blocksize_b( idx3, dimLen3, &C_1B, gemm_mc );
	dim_t idx4, dimLen4, bs4;
	//****
	obj_t AB_1_1;
	bli_acquire_mpart_b2t( BLIS_SUBPART1, idx3, bs3, &AB_1, &AB_1_1 );
	obj_t C_1B_0;
	bli_acquire_mpart_b2t( BLIS_SUBPART0, idx3, bs3, &C_1B, &C_1B_0 );
	obj_t C_1B_1;
	bli_acquire_mpart_b2t( BLIS_SUBPART1, idx3, bs3, &C_1B, &C_1B_1 );
	//------------------------------------//

	bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_ROW_PANELS, 
			     BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			     BLIS_BUFFER_FOR_A_BLOCK,
			     gemm_mr, gemm_kr, 
			     &AB_1_1, &packed_A_blk );
	bli_obj_set_struc( BLIS_SYMMETRIC, C_1B_1 );
	bli_obj_set_uplo( BLIS_LOWER, C_1B_1 );
	obj_t C_1B_1L, packed_B_panL;
	dim_t offL = 0;
	dim_t nL = bli_min( bli_obj_width_after_trans( C_1B_1 ), 
			    bli_obj_diag_offset_after_trans( C_1B_1 ) + bs3 );
	bli_acquire_mpart_l2r( BLIS_SUBPART1,
			       offL, nL, &C_1B_1, &C_1B_1L );
	bli_acquire_mpart_l2r( BLIS_SUBPART1,
			       offL, nL, &packed_B_pan, &packed_B_panL );
	bli_packm_blk_var2( &BLIS_ONE, &AB_1_1, &packed_A_blk );
	bli_herk_l_ker_var2( &BLIS_ONE, &packed_A_blk, &packed_B_panL, 
			     &BLIS_ONE, &C_1B_1L, (herk_t*)NULL );

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

void DxT_HerkNU( obj_t *alpha,
		 obj_t *A,
		 obj_t *beta,
		 obj_t *C )
{
  obj_t packed_B_pan;
  obj_t packed_A_blk;
  bli_obj_init_pack( &packed_B_pan );
  bli_obj_init_pack( &packed_A_blk );

obj_t AT;
bli_obj_alias_with_trans( BLIS_TRANSPOSE, *A, AT);
dim_t dimLen1 = bli_obj_width_after_trans( *C );
 dim_t bs1, idx1;
for ( idx1 = 0; idx1 < dimLen1; idx1 += bs1 ) {
	bs1 = bli_determine_blocksize_b( idx1, dimLen1, C, gemm_nc );
	dim_t idx2, dimLen2, bs2;
//****
	obj_t AT_1;
	bli_acquire_mpart_r2l( BLIS_SUBPART1, idx1, bs1, &AT, &AT_1 );
	obj_t C_1;
	bli_acquire_mpart_r2l( BLIS_SUBPART1, idx1, bs1, C, &C_1 );
	obj_t C_2;
	bli_acquire_mpart_r2l( BLIS_SUBPART2, idx1, bs1, C, &C_2 );
	//------------------------------------//

	obj_t C_1_T, A_T;
	dim_t offT = 0;
	dim_t mT = bli_min( bli_obj_length_after_trans( C_1 ), 
			-bli_obj_diag_offset_after_trans( C_1 ) + bs1 );
	bli_acquire_mpart_t2b( BLIS_SUBPART1,
			offT, mT, &C_1, &C_1_T );
	bli_acquire_mpart_t2b( BLIS_SUBPART1,
			offT, mT, A, &A_T );
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
		obj_t AT_1_1;
		bli_acquire_mpart_t2b( BLIS_SUBPART1, idx2, bs2, &AT_1, &AT_1_1 );
		//------------------------------------//

		bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_COL_PANELS, 
				BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
				BLIS_BUFFER_FOR_B_PANEL,
				gemm_kr, gemm_nr, 
				&AT_1_1, &packed_B_pan );
		bli_packm_blk_var2( &BLIS_ONE, &AT_1_1, &packed_B_pan );
		dimLen3 = bli_obj_length_after_trans( C_1_T );
		for ( idx3 = 0; idx3 < dimLen3; idx3 += bs3 ) {
			bs3 = bli_determine_blocksize_f( idx3, dimLen3, &C_1_T, gemm_mc );
			dim_t idx4, dimLen4, bs4;
		//****
			obj_t A_T_1_1;
			bli_acquire_mpart_t2b( BLIS_SUBPART1, idx3, bs3, &A_T_1, &A_T_1_1 );
			obj_t C_1_T_1;
			bli_acquire_mpart_t2b( BLIS_SUBPART1, idx3, bs3, &C_1_T, &C_1_T_1 );
			obj_t C_1_T_2;
			bli_acquire_mpart_t2b( BLIS_SUBPART2, idx3, bs3, &C_1_T, &C_1_T_2 );
			//------------------------------------//

			bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_ROW_PANELS, 
					BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
					BLIS_BUFFER_FOR_A_BLOCK,
					gemm_mr, gemm_kr, 
					&A_T_1_1, &packed_A_blk );
			bli_obj_set_struc( BLIS_SYMMETRIC, C_1_T_1 );
			bli_obj_set_uplo( BLIS_UPPER, C_1_T_1 );
			obj_t C_1_T_1_R, packed_B_pan_R;
			dim_t offR = bli_max( 0, bli_obj_diag_offset_after_trans( C_1_T_1 ) );
			dim_t nR = bli_obj_width_after_trans( C_1_T_1 ) - offR;
			bli_acquire_mpart_l2r( BLIS_SUBPART1,
					offR, nR, &C_1_T_1, &C_1_T_1_R );
			bli_acquire_mpart_l2r( BLIS_SUBPART1,
					offR, nR, &packed_B_pan, &packed_B_pan_R );
			bli_packm_blk_var2( &BLIS_ONE, &A_T_1_1, &packed_A_blk );
			bli_herk_u_ker_var2( &BLIS_ONE, &packed_A_blk, &packed_B_pan_R, 
					&BLIS_ONE, &C_1_T_1_R, (herk_t*)NULL );

			//------------------------------------//

		//****
		}

		//------------------------------------//

	//****
	}

	//------------------------------------//

//****
}


  bli_obj_release_pack( &packed_B_pan );
  bli_obj_release_pack( &packed_A_blk );
}

void DxT_HerkTL( obj_t *alpha,
		 obj_t *A,
		 obj_t *beta,
		 obj_t *C )
{
  obj_t packed_B_pan;
  obj_t packed_A_blk;
  bli_obj_init_pack( &packed_B_pan );
  bli_obj_init_pack( &packed_A_blk );

  dim_t idx1, dimLen1, bs1;
  obj_t AT;
  bli_obj_alias_with_trans( BLIS_TRANSPOSE, *A, AT);
  dimLen1 = bli_obj_width_after_trans( *C );
  for ( idx1 = 0; idx1 < dimLen1; idx1 += bs1 ) {
    bs1 = bli_determine_blocksize_f( idx1, dimLen1, C, gemm_nc );
    dim_t idx2, dimLen2, bs2;
    //****
    obj_t A_1;
    bli_acquire_mpart_l2r( BLIS_SUBPART1, idx1, bs1, A, &A_1 );
    obj_t C_0;
    bli_acquire_mpart_l2r( BLIS_SUBPART0, idx1, bs1, C, &C_0 );
    obj_t C_1;
    bli_acquire_mpart_l2r( BLIS_SUBPART1, idx1, bs1, C, &C_1 );
    //------------------------------------//

    obj_t C_1B, ATB;
    dim_t offB = bli_max( 0, -bli_obj_diag_offset_after_trans( C_1 ) );
    dim_t mB = bli_obj_length_after_trans( C_1 ) - offB;
    bli_acquire_mpart_t2b( BLIS_SUBPART1,
			   offB, mB, &C_1, &C_1B );
    bli_acquire_mpart_t2b( BLIS_SUBPART1,
			   offB, mB, &AT, &ATB );
    bli_obj_set_struc( BLIS_TRIANGULAR, C_1B );
    bli_obj_set_uplo( BLIS_LOWER, C_1B );
    bli_scalm( &BLIS_ONE, &C_1B );
    bli_obj_set_struc( BLIS_GENERAL, C_1B );
    bli_obj_set_uplo( BLIS_DENSE, C_1B );
    dimLen2 = bli_obj_width_after_trans( ATB );
    for ( idx2 = 0; idx2 < dimLen2; idx2 += bs2 ) {
      bs2 = bli_determine_blocksize_f( idx2, dimLen2, &ATB, gemm_kc );
      dim_t idx3, dimLen3, bs3;
      //****
      obj_t ATB_1;
      bli_acquire_mpart_l2r( BLIS_SUBPART1, idx2, bs2, &ATB, &ATB_1 );
      obj_t A_1_1;
      bli_acquire_mpart_t2b( BLIS_SUBPART1, idx2, bs2, &A_1, &A_1_1 );
      //------------------------------------//

      bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_COL_PANELS, 
			   BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			   BLIS_BUFFER_FOR_B_PANEL,
			   gemm_kr, gemm_nr, 
			   &A_1_1, &packed_B_pan );
      bli_packm_blk_var2( &BLIS_ONE, &A_1_1, &packed_B_pan );
      dimLen3 = bli_obj_length_after_trans( C_1B );
      for ( idx3 = 0; idx3 < dimLen3; idx3 += bs3 ) {
	bs3 = bli_determine_blocksize_b( idx3, dimLen3, &C_1B, gemm_mc );
	dim_t idx4, dimLen4, bs4;
	//****
	obj_t ATB_1_1;
	bli_acquire_mpart_b2t( BLIS_SUBPART1, idx3, bs3, &ATB_1, &ATB_1_1 );
	obj_t C_1B_0;
	bli_acquire_mpart_b2t( BLIS_SUBPART0, idx3, bs3, &C_1B, &C_1B_0 );
	obj_t C_1B_1;
	bli_acquire_mpart_b2t( BLIS_SUBPART1, idx3, bs3, &C_1B, &C_1B_1 );
	//------------------------------------//

	bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_ROW_PANELS, 
			     BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			     BLIS_BUFFER_FOR_A_BLOCK,
			     gemm_mr, gemm_kr, 
			     &ATB_1_1, &packed_A_blk );
	bli_obj_set_struc( BLIS_SYMMETRIC, C_1B_1 );
	bli_obj_set_uplo( BLIS_LOWER, C_1B_1 );
	obj_t C_1B_1L, packed_B_panL;
	dim_t offL = 0;
	dim_t nL = bli_min( bli_obj_width_after_trans( C_1B_1 ), 
			    bli_obj_diag_offset_after_trans( C_1B_1 ) + bs3 );
	bli_acquire_mpart_l2r( BLIS_SUBPART1,
			       offL, nL, &C_1B_1, &C_1B_1L );
	bli_acquire_mpart_l2r( BLIS_SUBPART1,
			       offL, nL, &packed_B_pan, &packed_B_panL );
	bli_packm_blk_var2( &BLIS_ONE, &ATB_1_1, &packed_A_blk );
	bli_herk_l_ker_var2( &BLIS_ONE, &packed_A_blk, &packed_B_panL, 
			     &BLIS_ONE, &C_1B_1L, (herk_t*)NULL );

	//------------------------------------//

	//****
      }

      //------------------------------------//

      //****
    }

    //------------------------------------//

    //****
  }


   bli_obj_release_pack( &packed_B_pan );
  bli_obj_release_pack( &packed_A_blk );
}

void DxT_HerkTU( obj_t *alpha,
		 obj_t *A,
		 obj_t *beta,
		 obj_t *C )
{
  obj_t packed_B_pan;
  obj_t packed_A_blk;
  bli_obj_init_pack( &packed_B_pan );
  bli_obj_init_pack( &packed_A_blk );
  
  dim_t idx1, dimLen1, bs1;
  obj_t AT;
  bli_obj_alias_with_trans( BLIS_TRANSPOSE, *A, AT);
  dimLen1 = bli_obj_width_after_trans( *C );
  for ( idx1 = 0; idx1 < dimLen1; idx1 += bs1 ) {
    bs1 = bli_determine_blocksize_b( idx1, dimLen1, C, gemm_nc );
    dim_t idx2, dimLen2, bs2;
    //****
    obj_t A_1;
    bli_acquire_mpart_r2l( BLIS_SUBPART1, idx1, bs1, A, &A_1 );
    obj_t C_1;
    bli_acquire_mpart_r2l( BLIS_SUBPART1, idx1, bs1, C, &C_1 );
    obj_t C_2;
    bli_acquire_mpart_r2l( BLIS_SUBPART2, idx1, bs1, C, &C_2 );
    //------------------------------------//

    obj_t C_1T, ATT;
    dim_t offT = 0;
    dim_t mT = bli_min( bli_obj_length_after_trans( C_1 ), 
			-bli_obj_diag_offset_after_trans( C_1 ) + bs1 );
    bli_acquire_mpart_t2b( BLIS_SUBPART1,
			   offT, mT, &C_1, &C_1T );
    bli_acquire_mpart_t2b( BLIS_SUBPART1,
			   offT, mT, &AT, &ATT );
    bli_obj_set_struc( BLIS_TRIANGULAR, C_1T );
    bli_obj_set_uplo( BLIS_LOWER, C_1T );
    bli_scalm( &BLIS_ONE, &C_1T );
    bli_obj_set_struc( BLIS_GENERAL, C_1T );
    bli_obj_set_uplo( BLIS_DENSE, C_1T );
    dimLen2 = bli_obj_width_after_trans( ATT );
    for ( idx2 = 0; idx2 < dimLen2; idx2 += bs2 ) {
      bs2 = bli_determine_blocksize_f( idx2, dimLen2, &ATT, gemm_kc );
      dim_t idx3, dimLen3, bs3;
      //****
      obj_t ATT_1;
      bli_acquire_mpart_l2r( BLIS_SUBPART1, idx2, bs2, &ATT, &ATT_1 );
      obj_t A_1_1;
      bli_acquire_mpart_t2b( BLIS_SUBPART1, idx2, bs2, &A_1, &A_1_1 );
      //------------------------------------//

      bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_COL_PANELS, 
			   BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			   BLIS_BUFFER_FOR_B_PANEL,
			   gemm_kr, gemm_nr, 
			   &A_1_1, &packed_B_pan );
      bli_packm_blk_var2( &BLIS_ONE, &A_1_1, &packed_B_pan );
      dimLen3 = bli_obj_length_after_trans( C_1T );
      for ( idx3 = 0; idx3 < dimLen3; idx3 += bs3 ) {
	bs3 = bli_determine_blocksize_f( idx3, dimLen3, &C_1T, gemm_mc );
	dim_t idx4, dimLen4, bs4;
	//****
	obj_t ATT_1_1;
	bli_acquire_mpart_t2b( BLIS_SUBPART1, idx3, bs3, &ATT_1, &ATT_1_1 );
	obj_t C_1T_1;
	bli_acquire_mpart_t2b( BLIS_SUBPART1, idx3, bs3, &C_1T, &C_1T_1 );
	obj_t C_1T_2;
	bli_acquire_mpart_t2b( BLIS_SUBPART2, idx3, bs3, &C_1T, &C_1T_2 );
	//------------------------------------//

	bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_ROW_PANELS, 
			     BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			     BLIS_BUFFER_FOR_A_BLOCK,
			     gemm_mr, gemm_kr, 
			     &ATT_1_1, &packed_A_blk );
	bli_obj_set_struc( BLIS_SYMMETRIC, C_1T_1 );
	bli_obj_set_uplo( BLIS_UPPER, C_1T_1 );
	obj_t C_1T_1R, packed_B_panR;
	dim_t offR = bli_max( 0, bli_obj_diag_offset_after_trans( C_1T_1 ) );
	dim_t nR = bli_obj_width_after_trans( C_1T_1 ) - offR;
	bli_acquire_mpart_l2r( BLIS_SUBPART1,
			       offR, nR, &C_1T_1, &C_1T_1R );
	bli_acquire_mpart_l2r( BLIS_SUBPART1,
			       offR, nR, &packed_B_pan, &packed_B_panR );
	bli_packm_blk_var2( &BLIS_ONE, &ATT_1_1, &packed_A_blk );
	bli_herk_u_ker_var2( &BLIS_ONE, &packed_A_blk, &packed_B_panR, 
			     &BLIS_ONE, &C_1T_1R, (herk_t*)NULL );

	//------------------------------------//

	//****
      }

      //------------------------------------//

      //****
    }

    //------------------------------------//

    //****
  }


   bli_obj_release_pack( &packed_B_pan );
  bli_obj_release_pack( &packed_A_blk );
}

int main( int argc, char** argv )
{
  obj_t a, ah, c1, c2;
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

	//bli_error_checking_level_set( BLIS_NO_ERROR_CHECKING );

	n_repeats = 3;

	p_begin = 20;
	p_end   = 20;
	p_inc   = 20;

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
		bli_obj_create( dt_alpha, 1, 1, 0, 0, &normVal );

		bli_obj_create( dt_a, m, k, 0, 0, &a );
		bli_obj_create( dt_c, m, m, 0, 0, &c1 );
		bli_obj_create( dt_c, m, m, 0, 0, &c2 );
		bli_obj_create( dt_c, m, m, 0, 0, &c_save );

		bli_randm( &a );
		bli_randm( &c1 );

		//		bli_obj_alias_with_trans( BLIS_CONJ_TRANSPOSE, a, ah );

		bli_setsc( (1.0/1.0), 0.0, &alpha );
		bli_setsc( (1.0/1.0), 0.0, &beta );

		bli_copym( &c1, &c_save );
	
		dtime_save = 1.0e9;

		for ( r = 0; r < n_repeats; ++r )
		  {
		    bli_copym( &c_save, &c1 );
		    bli_copym( &c_save, &c2 );


		    dtime = bli_clock();

		    if (trans) {
		      if (lower) {
			DxT_HerkTL(&alpha,
				   &a,
				   &beta,
				   &c1);				       
		      }
		      else { //upper
			DxT_HerkTU(&alpha,
				   &a,
				   &beta,
				   &c1);				       
		      }
		    }
		    else { //normal
		      if (lower) {
			DxT_HerkNL(&alpha,
				   &a,
				   &beta,
				   &c1);				       
		      }
		      else { //upper
			DxT_HerkNU(&alpha,
				   &a,
				   &beta,
				   &c1);				       
		      }
		    }


		    if (trans)
		      bli_obj_set_conjtrans( BLIS_TRANSPOSE, a);
			
		    bli_obj_set_struc( BLIS_SYMMETRIC, c2 );
		    if (lower) {
		      bli_obj_set_uplo( BLIS_LOWER, c2 );
		    }
		    else {
		      bli_obj_set_uplo( BLIS_UPPER, c2 );
		    }

		    bli_syrk( &alpha,
			      &a,
			      &beta,
			      &c2);


		    if (trans)
		      bli_obj_set_conjtrans( BLIS_NO_TRANSPOSE, a );
		    bli_obj_set_struc( BLIS_GENERAL, c2 );
		    bli_obj_set_uplo( BLIS_DENSE, c2 );
			
		    bli_axpym( &BLIS_MINUS_ONE, &c1, &c2 );
			
		    bli_fnormm( &c2, &normVal );

		    dtime_save = bli_clock_min_diff( dtime_save, dtime );
		  }

		gflops = ( 1.0 * m * k * m ) / ( dtime_save * 1.0e9 );

		printf( "data_herk_blis" );
		printf( "( %2ld, 1:4 ) = [ %4lu %4lu  %10.3e  %6.3f ];\n",
		        (p - p_begin + 1)/p_inc + 1, m, k, dtime_save, gflops );

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

