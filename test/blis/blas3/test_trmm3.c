#include <unistd.h>
#include "blis.h"

extern blksz_t *gemm_mc;
extern blksz_t *gemm_kc;
extern blksz_t *gemm_nc;
extern blksz_t *gemm_mr;
extern blksz_t *gemm_kr;
extern blksz_t *gemm_nr;
extern blksz_t *trmm_mr;



void DxT_TrmmLLN( obj_t *alpha,
		  obj_t *L,
		  obj_t *X,
		  obj_t *C )
{
  obj_t packed_A_blk;
  obj_t packed_B_pan;
  bli_obj_init_pack( &packed_A_blk );
  bli_obj_init_pack( &packed_B_pan );

  bli_scalm(alpha, X);


  bli_obj_release_pack( &packed_A_blk );
  bli_obj_release_pack( &packed_B_pan );
}

void DxT_TrmmLLT( obj_t *alpha,
		  obj_t *L,
		  obj_t *X, obj_t *C )
{
  obj_t packed_A_blk;
  obj_t packed_B_pan;
  bli_obj_init_pack( &packed_A_blk );
  bli_obj_init_pack( &packed_B_pan );

  bli_scalm(alpha, X);


  bli_obj_release_pack( &packed_A_blk );
  bli_obj_release_pack( &packed_B_pan );
}


void DxT_TrmmLUN( obj_t *alpha,
		  obj_t *U,
		  obj_t *X, obj_t *C )
{
  obj_t packed_A_blk;
  obj_t packed_B_pan;
  bli_obj_init_pack( &packed_A_blk );
  bli_obj_init_pack( &packed_B_pan );

  bli_scalm(alpha, X);


  bli_obj_release_pack( &packed_A_blk );
  bli_obj_release_pack( &packed_B_pan );
}

void DxT_TrmmLUT( obj_t *alpha,
		  obj_t *U,
		  obj_t *X, obj_t *C )
{
  obj_t packed_A_blk;
  obj_t packed_B_pan;
  bli_obj_init_pack( &packed_A_blk );
  bli_obj_init_pack( &packed_B_pan );

  bli_scalm(alpha, X);


  bli_obj_release_pack( &packed_A_blk );
  bli_obj_release_pack( &packed_B_pan );
}

void DxT_TrmmRLN( obj_t *alpha,
		  obj_t *L,
		  obj_t *X, obj_t *C )
{
  obj_t packed_A_blk;
  obj_t packed_B_pan;
  bli_obj_init_pack( &packed_A_blk );
  bli_obj_init_pack( &packed_B_pan );

  bli_scalm(alpha, X);


  bli_obj_release_pack( &packed_A_blk );
  bli_obj_release_pack( &packed_B_pan );
}

void DxT_TrmmRLT( obj_t *alpha,
		  obj_t *L,
		  obj_t *X, obj_t *C )
{
  obj_t packed_A_blk;
  obj_t packed_B_pan;
  bli_obj_init_pack( &packed_A_blk );
  bli_obj_init_pack( &packed_B_pan );

  bli_scalm(alpha, X);


  bli_obj_release_pack( &packed_A_blk );
  bli_obj_release_pack( &packed_B_pan );
}

void DxT_TrmmRUN( obj_t *alpha,
		  obj_t *U,
		  obj_t *X, obj_t *C )
{
  obj_t packed_A_blk;
  obj_t packed_B_pan;
  bli_obj_init_pack( &packed_A_blk );
  bli_obj_init_pack( &packed_B_pan );

  bli_scalm(alpha, X);


  bli_obj_release_pack( &packed_A_blk );
  bli_obj_release_pack( &packed_B_pan );
}


void DxT_TrmmRUT( obj_t *alpha,
		  obj_t *U,
		  obj_t *X, obj_t *C )
{
  obj_t packed_A_blk;
  obj_t packed_B_pan;
  bli_obj_init_pack( &packed_A_blk );
  bli_obj_init_pack( &packed_B_pan );

  bli_scalm(alpha, X);


  bli_obj_release_pack( &packed_A_blk );
  bli_obj_release_pack( &packed_B_pan );
}



int main( int argc, char** argv )
{
  obj_t a, x, x2, c1, c2, b;
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

  if (*(argv[3]) == 'T')
    trans = 1;
  else if (*(argv[3]) != 'N') {
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
      bli_obj_create( dt_c, m, n, 0, 0, &b );
      bli_obj_create( dt_c, m, n, 0, 0, &x );
      bli_obj_create( dt_c, m, n, 0, 0, &c1 );
      bli_obj_create( dt_c, m, n, 0, 0, &c2 );
      bli_obj_create( dt_c, m, n, 0, 0, &c_save );


      bli_randm( &a );
      bli_randm( &b );
      bli_randm( &x );
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
			    &a, &x,
			    &c1 );
	      }
	      else { //normal
		DxT_TrmmLLN(&alpha,
			    &a, &x,
			    &c1 );
	      }
	    }
	    else { //upper
	      if (trans) {
		DxT_TrmmLUT(&alpha,
			    &a, &x,
			    &c1 );
	      }
	      else {//normal 
		DxT_TrmmLUN(&alpha,
			    &a, &x,
			    &c1 );
	      }
	    }
	  }
	  else { //right
	    if (lower) {
	      if (trans) {
		DxT_TrmmRLT(&alpha,
			    &a, &x,
			    &c1 );
	      }
	      else { //normal
		DxT_TrmmRLN(&alpha,
			    &a, &x,
			    &c1 );
	      }

	    }
	    else { //upper
	      if (trans) {
		DxT_TrmmRUT(&alpha,
			    &a, &x,
			    &c1 );
	      }
	      else { //normal
		DxT_TrmmRUN(&alpha,
			    &a, &x,
			    &c1 );
	      }
	    }
	  }

	  if (trans)
	    bli_obj_set_trans( BLIS_TRANSPOSE, a);

	  bli_obj_set_struc( BLIS_TRIANGULAR, a );
	  if (lower) {
	    bli_obj_set_uplo( BLIS_LOWER, a );
	  }
	  else {
	    bli_obj_set_uplo( BLIS_UPPER, a );
	  }

	  bli_trmm3( left ? BLIS_LEFT : BLIS_RIGHT,
		     &alpha,
		     &a, &x,
		     &BLIS_ONE,
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

