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
  obj_t A_1T_packed;
  obj_t A_1_1_packed;
    bli_obj_init_pack( &A_1T_packed );
    bli_obj_init_pack( &A_1_1_packed );

dim_t idx1, dimLen1, bs1;
bli_obj_set_struc( BLIS_TRIANGULAR, *C );
bli_obj_set_uplo( BLIS_LOWER, *C );
bli_scalm( &BLIS_ONE, C );
bli_obj_set_struc( BLIS_GENERAL, *C );
bli_obj_set_uplo( BLIS_DENSE, *C );
///// Blocksize = 256
dimLen1 = bli_obj_width_after_trans( *A );
for ( idx1 = 0; idx1 < dimLen1; idx1 += bs1 ) {
	bs1 = bli_determine_blocksize_f( idx1, dimLen1, A, gemm_kc );
	dim_t idx2, dimLen2, bs2;
//****
	obj_t A_1;
	bli_acquire_mpart_l2r( BLIS_SUBPART1, idx1, bs1, A, &A_1 );
	//------------------------------------//

	obj_t A_1T;
	bli_obj_alias_with_trans( BLIS_TRANSPOSE, A_1, A_1T);
	bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_COL_PANELS, 
			BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			BLIS_BUFFER_FOR_B_PANEL,
			gemm_kr,  gemm_nr,  
			&A_1T, &A_1T_packed );
	bli_packm_blk_var2( &BLIS_ONE, &A_1T, &A_1T_packed );
	///// Blocksize = 128
	dimLen2 = bli_obj_length_after_trans( *C );
	for ( idx2 = 0; idx2 < dimLen2; idx2 += bs2 ) {
		bs2 = bli_determine_blocksize_b( idx2, dimLen2, C, gemm_mc );
		dim_t idx3, dimLen3, bs3;
	//****
		obj_t A_1_1;
		bli_acquire_mpart_b2t( BLIS_SUBPART1, idx2, bs2, &A_1, &A_1_1 );
		obj_t C_0;
		bli_acquire_mpart_b2t( BLIS_SUBPART0, idx2, bs2, C, &C_0 );
		obj_t C_1;
		bli_acquire_mpart_b2t( BLIS_SUBPART1, idx2, bs2, C, &C_1 );
		//------------------------------------//

		bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_ROW_PANELS, 
				BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
				BLIS_BUFFER_FOR_A_BLOCK,
				gemm_mr,  gemm_kr,  
				&A_1_1, &A_1_1_packed );
		bli_obj_set_struc( BLIS_SYMMETRIC, C_1 );
		bli_obj_set_uplo( BLIS_LOWER, C_1 );
		obj_t C_1L, A_1T_packedL;
		dim_t offL = 0;
		dim_t nL = bli_min( bli_obj_width_after_trans( C_1 ), 
				bli_obj_diag_offset_after_trans( C_1 ) + bs2 );
		bli_acquire_mpart_l2r( BLIS_SUBPART1,
				offL, nL, &C_1, &C_1L );
				bli_acquire_mpart_l2r( BLIS_SUBPART1,
				offL, nL, &A_1T_packed, &A_1T_packedL );
		bli_packm_blk_var2( &BLIS_ONE, &A_1_1, &A_1_1_packed );
		bli_herk_l_ker_var2( &BLIS_ONE, &A_1_1_packed, &A_1T_packedL, 
				&BLIS_ONE, &C_1L, (herk_t*)NULL );

		//------------------------------------//

	//****
	}

	//------------------------------------//

//****
}

    bli_obj_release_pack( &A_1T_packed );
    bli_obj_release_pack( &A_1_1_packed );
}

void DxT_HerkNU( obj_t *alpha,
		 obj_t *A,
		 obj_t *beta,
		 obj_t *C )
{
  obj_t A_1T_packed;
  obj_t A_1_1_packed;
    bli_obj_init_pack( &A_1T_packed );
    bli_obj_init_pack( &A_1_1_packed );



dim_t idx1, dimLen1, bs1;
bli_obj_set_struc( BLIS_TRIANGULAR, *C );
bli_obj_set_uplo( BLIS_LOWER, *C );
bli_scalm( &BLIS_ONE, C );
bli_obj_set_struc( BLIS_GENERAL, *C );
bli_obj_set_uplo( BLIS_DENSE, *C );
///// Blocksize = 256
dimLen1 = bli_obj_width_after_trans( *A );
for ( idx1 = 0; idx1 < dimLen1; idx1 += bs1 ) {
	bs1 = bli_determine_blocksize_f( idx1, dimLen1, A, gemm_kc );
	dim_t idx2, dimLen2, bs2;
//****
	obj_t A_1;
	bli_acquire_mpart_l2r( BLIS_SUBPART1, idx1, bs1, A, &A_1 );
	//------------------------------------//

	obj_t A_1T;
	bli_obj_alias_with_trans( BLIS_TRANSPOSE, A_1, A_1T);
	bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_COL_PANELS, 
			BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			BLIS_BUFFER_FOR_B_PANEL,
			gemm_kr,  gemm_nr,  
			&A_1T, &A_1T_packed );
	bli_packm_blk_var2( &BLIS_ONE, &A_1T, &A_1T_packed );
	///// Blocksize = 128
	dimLen2 = bli_obj_length_after_trans( *C );
	for ( idx2 = 0; idx2 < dimLen2; idx2 += bs2 ) {
		bs2 = bli_determine_blocksize_f( idx2, dimLen2, C, gemm_mc );
		dim_t idx3, dimLen3, bs3;
	//****
		obj_t A_1_1;
		bli_acquire_mpart_t2b( BLIS_SUBPART1, idx2, bs2, &A_1, &A_1_1 );
		obj_t C_1;
		bli_acquire_mpart_t2b( BLIS_SUBPART1, idx2, bs2, C, &C_1 );
		obj_t C_2;
		bli_acquire_mpart_t2b( BLIS_SUBPART2, idx2, bs2, C, &C_2 );
		//------------------------------------//

		bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_ROW_PANELS, 
				BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
				BLIS_BUFFER_FOR_A_BLOCK,
				gemm_mr,  gemm_kr,  
				&A_1_1, &A_1_1_packed );
		bli_obj_set_struc( BLIS_SYMMETRIC, C_1 );
		bli_obj_set_uplo( BLIS_UPPER, C_1 );
		obj_t C_1R, A_1T_packedR;
		dim_t offR = bli_max( 0, bli_obj_diag_offset_after_trans( C_1 ) );
		dim_t nR = bli_obj_width_after_trans( C_1 ) - offR;
		bli_acquire_mpart_l2r( BLIS_SUBPART1,
				offR, nR, &C_1, &C_1R );
				bli_acquire_mpart_l2r( BLIS_SUBPART1,
				offR, nR, &A_1T_packed, &A_1T_packedR );
		bli_packm_blk_var2( &BLIS_ONE, &A_1_1, &A_1_1_packed );
		bli_herk_u_ker_var2( &BLIS_ONE, &A_1_1_packed, &A_1T_packedR, 
				&BLIS_ONE, &C_1R, (herk_t*)NULL );

		//------------------------------------//

	//****
	}

	//------------------------------------//

//****
}


    bli_obj_release_pack( &A_1T_packed );
    bli_obj_release_pack( &A_1_1_packed );
}

void DxT_HerkTL( obj_t *alpha,
		 obj_t *A,
		 obj_t *beta,
		 obj_t *C )
{
  obj_t A_1_packed, A_1T_1_packed;
    bli_obj_init_pack( &A_1_packed );
    bli_obj_init_pack( &A_1T_1_packed );

dim_t idx1, dimLen1, bs1;
bli_obj_set_struc( BLIS_TRIANGULAR, *C );
bli_obj_set_uplo( BLIS_LOWER, *C );
bli_scalm( &BLIS_ONE, C );
bli_obj_set_struc( BLIS_GENERAL, *C );
bli_obj_set_uplo( BLIS_DENSE, *C );
///// Blocksize = 256
dimLen1 = bli_obj_length_after_trans( *A );
for ( idx1 = 0; idx1 < dimLen1; idx1 += bs1 ) {
	bs1 = bli_determine_blocksize_f( idx1, dimLen1, A, gemm_kc );
	dim_t idx2, dimLen2, bs2;
//****
	obj_t A_1;
	bli_acquire_mpart_t2b( BLIS_SUBPART1, idx1, bs1, A, &A_1 );
	//------------------------------------//

	obj_t A_1T;
	bli_obj_alias_with_trans( BLIS_TRANSPOSE, A_1, A_1T);
	bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_COL_PANELS, 
			BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			BLIS_BUFFER_FOR_B_PANEL,
			gemm_kr,  gemm_nr,  
			&A_1, &A_1_packed );
	bli_packm_blk_var2( &BLIS_ONE, &A_1, &A_1_packed );
	///// Blocksize = 128
	dimLen2 = bli_obj_length_after_trans( *C );
	for ( idx2 = 0; idx2 < dimLen2; idx2 += bs2 ) {
		bs2 = bli_determine_blocksize_b( idx2, dimLen2, C, gemm_mc );
		dim_t idx3, dimLen3, bs3;
	//****
		obj_t A_1T_1;
		bli_acquire_mpart_b2t( BLIS_SUBPART1, idx2, bs2, &A_1T, &A_1T_1 );
		obj_t C_0;
		bli_acquire_mpart_b2t( BLIS_SUBPART0, idx2, bs2, C, &C_0 );
		obj_t C_1;
		bli_acquire_mpart_b2t( BLIS_SUBPART1, idx2, bs2, C, &C_1 );
		//------------------------------------//

		bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_ROW_PANELS, 
				BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
				BLIS_BUFFER_FOR_A_BLOCK,
				gemm_mr,  gemm_kr,  
				&A_1T_1, &A_1T_1_packed );
		bli_obj_set_struc( BLIS_SYMMETRIC, C_1 );
		bli_obj_set_uplo( BLIS_LOWER, C_1 );
		obj_t C_1L, A_1_packedL;
		dim_t offL = 0;
		dim_t nL = bli_min( bli_obj_width_after_trans( C_1 ), 
				bli_obj_diag_offset_after_trans( C_1 ) + bs2 );
		bli_acquire_mpart_l2r( BLIS_SUBPART1,
				offL, nL, &C_1, &C_1L );
				bli_acquire_mpart_l2r( BLIS_SUBPART1,
				offL, nL, &A_1_packed, &A_1_packedL );
		bli_packm_blk_var2( &BLIS_ONE, &A_1T_1, &A_1T_1_packed );
		bli_herk_l_ker_var2( &BLIS_ONE, &A_1T_1_packed, &A_1_packedL, 
				&BLIS_ONE, &C_1L, (herk_t*)NULL );

		//------------------------------------//

	//****
	}

	//------------------------------------//

//****
}

  
  bli_obj_release_pack( &A_1_packed );
    bli_obj_release_pack( &A_1T_1_packed );
}

void DxT_HerkTU( obj_t *alpha,
		 obj_t *A,
		 obj_t *beta,
		 obj_t *C )
{
  obj_t A_1_packed, A_1T_1_packed;
  bli_obj_init_pack( &A_1_packed );
  bli_obj_init_pack( &A_1T_1_packed );
  
dim_t idx1, dimLen1, bs1;
bli_obj_set_struc( BLIS_TRIANGULAR, *C );
bli_obj_set_uplo( BLIS_LOWER, *C );
bli_scalm( &BLIS_ONE, C );
bli_obj_set_struc( BLIS_GENERAL, *C );
bli_obj_set_uplo( BLIS_DENSE, *C );
///// Blocksize = 256
dimLen1 = bli_obj_length_after_trans( *A );
for ( idx1 = 0; idx1 < dimLen1; idx1 += bs1 ) {
	bs1 = bli_determine_blocksize_f( idx1, dimLen1, A, gemm_kc );
	dim_t idx2, dimLen2, bs2;
//****
	obj_t A_1;
	bli_acquire_mpart_t2b( BLIS_SUBPART1, idx1, bs1, A, &A_1 );
	//------------------------------------//

	obj_t A_1T;
	bli_obj_alias_with_trans( BLIS_TRANSPOSE, A_1, A_1T);
	bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_COL_PANELS, 
			BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			BLIS_BUFFER_FOR_B_PANEL,
			gemm_kr,  gemm_nr,  
			&A_1, &A_1_packed );
	bli_packm_blk_var2( &BLIS_ONE, &A_1, &A_1_packed );
	///// Blocksize = 128
	dimLen2 = bli_obj_length_after_trans( *C );
	for ( idx2 = 0; idx2 < dimLen2; idx2 += bs2 ) {
		bs2 = bli_determine_blocksize_f( idx2, dimLen2, C, gemm_mc );
		dim_t idx3, dimLen3, bs3;
	//****
		obj_t A_1T_1;
		bli_acquire_mpart_t2b( BLIS_SUBPART1, idx2, bs2, &A_1T, &A_1T_1 );
		obj_t C_1;
		bli_acquire_mpart_t2b( BLIS_SUBPART1, idx2, bs2, C, &C_1 );
		obj_t C_2;
		bli_acquire_mpart_t2b( BLIS_SUBPART2, idx2, bs2, C, &C_2 );
		//------------------------------------//

		bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_ROW_PANELS, 
				BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
				BLIS_BUFFER_FOR_A_BLOCK,
				gemm_mr,  gemm_kr,  
				&A_1T_1, &A_1T_1_packed );
		bli_obj_set_struc( BLIS_SYMMETRIC, C_1 );
		bli_obj_set_uplo( BLIS_UPPER, C_1 );
		obj_t C_1R, A_1_packedR;
		dim_t offR = bli_max( 0, bli_obj_diag_offset_after_trans( C_1 ) );
		dim_t nR = bli_obj_width_after_trans( C_1 ) - offR;
		bli_acquire_mpart_l2r( BLIS_SUBPART1,
				offR, nR, &C_1, &C_1R );
				bli_acquire_mpart_l2r( BLIS_SUBPART1,
				offR, nR, &A_1_packed, &A_1_packedR );
		bli_packm_blk_var2( &BLIS_ONE, &A_1T_1, &A_1T_1_packed );
		bli_herk_u_ker_var2( &BLIS_ONE, &A_1T_1_packed, &A_1_packedR, 
				&BLIS_ONE, &C_1R, (herk_t*)NULL );

		//------------------------------------//

	//****
	}

	//------------------------------------//

//****
}


    bli_obj_release_pack( &A_1_packed );
    bli_obj_release_pack( &A_1T_1_packed );
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
	  return;
	}

	if (*(argv[1]) == 'L')
	  lower = 1;
	else if (*(argv[1]) != 'U') {
	  printf("lower/upper not correct\n");
	  return;
	}
	else
	  lower = 0;

	if (*(argv[2]) == 'T')
	  trans = 1;
	else if (*(argv[2]) != 'N') {
	  printf("trans/normal not correct\n");
	  return;
	}
	else
	  trans = 0;


	bli_init();

	n_repeats = 3;

	p_begin = 40;
	p_end   = 600;
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
		      bli_obj_set_trans( BLIS_TRANSPOSE, a);
			
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
		      bli_obj_toggle_trans( a );
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

