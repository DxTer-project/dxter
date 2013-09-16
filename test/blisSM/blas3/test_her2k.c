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




  FUNCTIONEND
}

void DxT_Her2kNU( obj_t *alpha,
		  obj_t *A,
		  obj_t *B,
		  obj_t *beta,
		  obj_t *C )
{
  FUNCTIONSTART



  FUNCTIONEND
}

void DxT_Her2kTL( obj_t *alpha,
		  obj_t *A,
		  obj_t *B,
		  obj_t *beta,
		  obj_t *C )
{
  FUNCTIONSTART




  FUNCTIONEND
}

void DxT_Her2kTU( obj_t *alpha,
		  obj_t *A,
		  obj_t *B,
		  obj_t *beta,
		  obj_t *C )
{
  FUNCTIONSTART




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

