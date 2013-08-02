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

void DxT_GemmNN( obj_t *alpha,
	  obj_t *A,
	  obj_t *B,
	  obj_t *beta,
	  obj_t *C )
{
  obj_t A_1_1_packed;
  obj_t B_1_packed;
  bli_obj_init_pack( &A_1_1_packed );
  bli_obj_init_pack( &B_1_packed );

  bli_scalm(beta, C);

  dimLen1 = bli_obj_width_after_trans( C );
  idx = 0;
  th_shift_start_end(&idx1, &dimLen1, GlobalComm);
  for ( ; idx1 < dimLen1; idx1 += bs1 ) {
    bs1 = bli_determine_blocksize_f( idx1, dimLen1, &C, gemm_nc );
    dim_t idx2, dimLen2, bs2;
    //****
    obj_t B_1;
    bli_acquire_mpart_l2r( BLIS_SUBPART1, idx1, bs1, &B, &B_1 );
    obj_t C_1;
    bli_acquire_mpart_l2r( BLIS_SUBPART1, idx1, bs1, &C, &C_1 );
    //------------------------------------//

    dimLen2 = bli_obj_width_after_trans( A );
    for ( idx2 = 0; idx2 < dimLen2; idx2 += bs2 ) {
      bs2 = bli_determine_blocksize_f( idx2, dimLen2, &A, gemm_kc );
      dim_t idx3, dimLen3, bs3;
      //****
      obj_t A_1;
      bli_acquire_mpart_l2r( BLIS_SUBPART1, idx2, bs2, &A, &A_1 );
      obj_t B_1_1;
      bli_acquire_mpart_t2b( BLIS_SUBPART1, idx2, bs2, &B_1, &B_1_1 );
      //------------------------------------//

      if (th_am_root(ProcComm)) {
	bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_COL_PANELS, 
			     BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			     BLIS_BUFFER_FOR_B_PANEL,
			     gemm_kr, gemm_nr, 
			     &B_1_1, &packed_B_pan );
      }
      Broadcast(ProcComm, 0, (void*)(&packed_B_pan), sizeof(packed_B_pan);
		bli_packm_blk_var2_par( ProcComm, &BLIS_ONE, &B_1_1, &packed_B_pan );
		//// ***Parallelized with communicator ProcComm; need correct output code
		dimLen3 = bli_obj_length_after_trans( C_1 );
		idx = 0;
		th_shift_start_end(&idx3, &dimLen3, ProcComm);
		for ( ; idx3 < dimLen3; idx3 += bs3 ) {
		  bs3 = bli_determine_blocksize_f( idx3, dimLen3, &C_1, gemm_mc );
		  dim_t idx4, dimLen4, bs4;
		  //****
		  obj_t A_1_1;
		  bli_acquire_mpart_t2b( BLIS_SUBPART1, idx3, bs3, &A_1, &A_1_1 );
		  obj_t C_1_1;
		  bli_acquire_mpart_t2b( BLIS_SUBPART1, idx3, bs3, &C_1, &C_1_1 );
		  //------------------------------------//

		  if (th_am_root(L2Comm)) {
		    bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_ROW_PANELS, 
					 BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
					 BLIS_BUFFER_FOR_A_BLOCK,
					 gemm_mr, gemm_kr, 
					 &A_1_1, &packed_A_blk );
		  }
		  Broadcast(L2Comm, 0, (void*)(&packed_A_blk), sizeof(packed_A_blk);
			    bli_packm_blk_var2_par( L2Comm, &BLIS_ONE, &A_1_1, &packed_A_blk );
			    bli_gemm_ker_var2_par( L2Comm, &BLIS_ONE, &packed_A_blk, &packed_B_pan, 
						   &BLIS_ONE, &C_1_1, (gemm_t*)NULL );

			    //------------------------------------//

			    //****
			    }

		    //------------------------------------//

		    //****
		}

		//------------------------------------//

		//****
		}

  bli_obj_release_pack( &A_1_1_packed );
  bli_obj_release_pack( &B_1_packed );
}

void DxT_GemmTN( obj_t *alpha,
	  obj_t *A,
	  obj_t *B,
	  obj_t *beta,
	  obj_t *C )
{
  obj_t A_1_1T_packed;
  obj_t B_1_packed;
  bli_obj_init_pack( &A_1_1T_packed );
  bli_obj_init_pack( &B_1_packed );

  bli_scalm(beta, C);

dim_t idx1, dimLen1, bs1;
///// Blocksize = 256
dimLen1 = bli_obj_length_after_trans( *A );
for ( idx1 = 0; idx1 < dimLen1; idx1 += bs1 ) {
	bs1 = bli_determine_blocksize_f( idx1, dimLen1, A, gemm_kc );
	dim_t idx2, dimLen2, bs2;
//****
	obj_t A_1;
	bli_acquire_mpart_t2b( BLIS_SUBPART1, idx1, bs1, A, &A_1 );
	obj_t B_1;
	bli_acquire_mpart_t2b( BLIS_SUBPART1, idx1, bs1, B, &B_1 );
	//------------------------------------//

	bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_COL_PANELS, 
			BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			BLIS_BUFFER_FOR_B_PANEL,
			gemm_kr,  gemm_nr,  
			&B_1, &B_1_packed );
	bli_packm_blk_var2( &BLIS_ONE, &B_1, &B_1_packed );
	///// Blocksize = 128
	dimLen2 = bli_obj_length_after_trans( *C );
	for ( idx2 = 0; idx2 < dimLen2; idx2 += bs2 ) {
		bs2 = bli_determine_blocksize_f( idx2, dimLen2, C, gemm_mc );
		dim_t idx3, dimLen3, bs3;
	//****
		obj_t A_1_1;
		bli_acquire_mpart_l2r( BLIS_SUBPART1, idx2, bs2, &A_1, &A_1_1 );
		obj_t C_1;
		bli_acquire_mpart_t2b( BLIS_SUBPART1, idx2, bs2, C, &C_1 );
		//------------------------------------//

		obj_t A_1_1T;
		bli_obj_alias_with_trans( BLIS_TRANSPOSE, A_1_1, A_1_1T);
		bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_ROW_PANELS, 
				BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
				BLIS_BUFFER_FOR_A_BLOCK,
				gemm_mr,  gemm_kr,  
				&A_1_1T, &A_1_1T_packed );
		bli_packm_blk_var2( &BLIS_ONE, &A_1_1T, &A_1_1T_packed );
		bli_gemm_ker_var2( &BLIS_ONE, &A_1_1T_packed, &B_1_packed, 
				&BLIS_ONE, &C_1, (gemm_t*)NULL );

		//------------------------------------//

	//****
	}

	//------------------------------------//

//****
}


  bli_obj_release_pack( &A_1_1T_packed );
  bli_obj_release_pack( &B_1_packed );
}

void DxT_GemmNT( obj_t *alpha,
	  obj_t *A,
	  obj_t *B,
	  obj_t *beta,
	  obj_t *C )
{
  obj_t A_1_1_packed;
  obj_t B_1T_packed;
  bli_obj_init_pack( &A_1_1_packed );
  bli_obj_init_pack( &B_1T_packed );

  bli_scalm(beta, C);

dim_t idx1, dimLen1, bs1;
///// Blocksize = 256
dimLen1 = bli_obj_width_after_trans( *A );
for ( idx1 = 0; idx1 < dimLen1; idx1 += bs1 ) {
	bs1 = bli_determine_blocksize_f( idx1, dimLen1, A, gemm_kc );
	dim_t idx2, dimLen2, bs2;
//****
	obj_t A_1;
	bli_acquire_mpart_l2r( BLIS_SUBPART1, idx1, bs1, A, &A_1 );
	obj_t B_1;
	bli_acquire_mpart_l2r( BLIS_SUBPART1, idx1, bs1, B, &B_1 );
	//------------------------------------//

	obj_t B_1T;
	bli_obj_alias_with_trans( BLIS_TRANSPOSE, B_1, B_1T);
	bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_COL_PANELS, 
			BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			BLIS_BUFFER_FOR_B_PANEL,
			gemm_kr,  gemm_nr,  
			&B_1T, &B_1T_packed );
	bli_packm_blk_var2( &BLIS_ONE, &B_1T, &B_1T_packed );
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
		//------------------------------------//

		bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_ROW_PANELS, 
				BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
				BLIS_BUFFER_FOR_A_BLOCK,
				gemm_mr,  gemm_kr,  
				&A_1_1, &A_1_1_packed );
		bli_packm_blk_var2( &BLIS_ONE, &A_1_1, &A_1_1_packed );
		bli_gemm_ker_var2( &BLIS_ONE, &A_1_1_packed, &B_1T_packed, 
				&BLIS_ONE, &C_1, (gemm_t*)NULL );

		//------------------------------------//

	//****
	}

	//------------------------------------//

//****
}




  bli_obj_release_pack( &A_1_1_packed );
  bli_obj_release_pack( &B_1T_packed );
}

void DxT_GemmTT( obj_t *alpha,
	  obj_t *A,
	  obj_t *B,
	  obj_t *beta,
	  obj_t *C )
{
  obj_t A_1_1T_packed;
  obj_t B_1T_packed;
  bli_obj_init_pack( &A_1_1T_packed );
  bli_obj_init_pack( &B_1T_packed );

  bli_scalm(beta, C);

dim_t idx1, dimLen1, bs1;
///// Blocksize = 256
dimLen1 = bli_obj_length_after_trans( *A );
for ( idx1 = 0; idx1 < dimLen1; idx1 += bs1 ) {
	bs1 = bli_determine_blocksize_f( idx1, dimLen1, A, gemm_kc );
	dim_t idx2, dimLen2, bs2;
//****
	obj_t A_1;
	bli_acquire_mpart_t2b( BLIS_SUBPART1, idx1, bs1, A, &A_1 );
	obj_t B_1;
	bli_acquire_mpart_l2r( BLIS_SUBPART1, idx1, bs1, B, &B_1 );
	//------------------------------------//

	obj_t B_1T;
	bli_obj_alias_with_trans( BLIS_TRANSPOSE, B_1, B_1T);
	bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_COL_PANELS, 
			BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
			BLIS_BUFFER_FOR_B_PANEL,
			gemm_kr,  gemm_nr,  
			&B_1T, &B_1T_packed );
	bli_packm_blk_var2( &BLIS_ONE, &B_1T, &B_1T_packed );
	///// Blocksize = 128
	dimLen2 = bli_obj_length_after_trans( *C );
	for ( idx2 = 0; idx2 < dimLen2; idx2 += bs2 ) {
		bs2 = bli_determine_blocksize_f( idx2, dimLen2, C, gemm_mc );
		dim_t idx3, dimLen3, bs3;
	//****
		obj_t A_1_1;
		bli_acquire_mpart_l2r( BLIS_SUBPART1, idx2, bs2, &A_1, &A_1_1 );
		obj_t C_1;
		bli_acquire_mpart_t2b( BLIS_SUBPART1, idx2, bs2, C, &C_1 );
		//------------------------------------//

		obj_t A_1_1T;
		bli_obj_alias_with_trans( BLIS_TRANSPOSE, A_1_1, A_1_1T);
		bli_packm_init_pack( FALSE, BLIS_NO_INVERT_DIAG, BLIS_PACKED_ROW_PANELS, 
				BLIS_PACK_FWD_IF_UPPER, BLIS_PACK_FWD_IF_LOWER, 
				BLIS_BUFFER_FOR_A_BLOCK,
				gemm_mr,  gemm_kr,  
				&A_1_1T, &A_1_1T_packed );
		bli_packm_blk_var2( &BLIS_ONE, &A_1_1T, &A_1_1T_packed );
		bli_gemm_ker_var2( &BLIS_ONE, &A_1_1T_packed, &B_1T_packed, 
				&BLIS_ONE, &C_1, (gemm_t*)NULL );

		//------------------------------------//

	//****
	}

	//------------------------------------//

//****
}


 bli_obj_release_pack( &A_1_1T_packed );
  bli_obj_release_pack( &B_1T_packed );
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

  double dtime;
  double dtime_save;
  double gflops;

  int transA=0, transB=0;

  dt_alpha = BLIS_DOUBLE;
  dt_beta = BLIS_DOUBLE;

  bli_obj_create( dt_alpha, 1, 1, 0, 0, &negOne );
  bli_setsc(  -1.0, 0.0, &negOne );		
			
  if (argc != 3) {
    printf("test T/N T/N\n");
    return;
  }

  if (*(argv[1]) == 'T')
    transA = 1;
  else if (*(argv[1]) != 'N') {
    printf("transA not correct\n");
    return;
  }

  if (*(argv[2]) == 'T')
    transB = 1;
  else if (*(argv[2]) != 'N') {
    printf("transB not correct\n");
    return;
  }

  bli_init();

  n_repeats = 3;

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
	
      dtime_save = 1.0e9;

      for ( r = 0; r < n_repeats; ++r )
	{
	  bli_copym( &c_save, &c1 );
	  bli_copym( &c_save, &c2 );


	  dtime = bli_clock();

	  //bli_error_checking_level_set( BLIS_NO_ERROR_CHECKING );
	  if (!transA && !transB) {
	    DxT_GemmNN( &alpha,
			&a,
			&b,
			&beta,
			&c2 );
			
	    bli_gemm( &alpha,
		      &a,
		      &b,
		      &beta,
		      &c1);
	  }
	  else if (transA && !transB) {
	    DxT_GemmTN( &alpha,
			&a,
			&b,
			&beta,
			&c2 );

	    bli_obj_set_conjtrans( BLIS_TRANSPOSE, a);
	    bli_gemm( &alpha,
		      &a,
		      &b,
		      &beta,
		      &c1);
	    bli_obj_toggle_trans( a );
	  }
	  else if (!transA) {// &&transB 
	    DxT_GemmNT( &alpha,
			&a,
			&b,
			&beta,
			&c2 );

	    bli_obj_set_conjtrans( BLIS_TRANSPOSE, b);
	    bli_gemm( &alpha,
		      &a,
		      &b,
		      &beta,
		      &c1);
	    bli_obj_toggle_trans( b );
	  }
	  else {//transA && transB
	    DxT_GemmTT( &alpha,
			&a,
			&b,
			&beta,
			&c2 );

	    bli_obj_set_conjtrans( BLIS_TRANSPOSE, a);
	    bli_obj_set_conjtrans( BLIS_TRANSPOSE, b);
	    bli_gemm( &alpha,
		      &a,
		      &b,
		      &beta,
		      &c1);
	    bli_obj_toggle_trans( a );
	    bli_obj_toggle_trans( b );
	  }


	  bli_axpym( &negOne, &c1, &c2 );

	  bli_obj_create( dt_alpha, 1, 1, 0, 0, &normVal );

	  bli_fnormm( &c2, &normVal );

	  dtime_save = bli_clock_min_diff( dtime_save, dtime );
	}

      gflops = ( 2.0 * m * k * n ) / ( dtime_save * 1.0e9 );

      printf( "data_gemm_blis" );
      printf( "( %2ld, 1:5 ) = [ %4lu %4lu %4lu  %10.3e  %6.3f ];\n",
	      (p - p_begin + 1)/p_inc + 1, m, k, n, dtime_save, gflops );

      bli_printm( "NORM", &normVal, "%4.1f", "" );

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

