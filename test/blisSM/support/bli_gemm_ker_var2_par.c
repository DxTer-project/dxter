/*

   BLIS    
   An object-based framework for developing high-performance BLAS-like
   libraries.

   Copyright (C) 2013, The University of Texas

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are
   met:
    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    - Neither the name of The University of Texas nor the names of its
      contributors may be used to endorse or promote products derived
      from this software without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
   HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
   LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
   DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
   THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
   OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#include "blis.h"
#include "bli_gemm_ker_var2_par.h"
#include "support.h"

#define FUNCPTR_T gemm_fp

typedef void (*FUNCPTR_T)(
                           dim_t   m,
                           dim_t   n,
                           dim_t   k,
                           void*   alpha,
                           void*   a, inc_t rs_a, inc_t cs_a, inc_t ps_a,
                           void*   b, inc_t rs_b, inc_t cs_b, inc_t ps_b,
                           void*   beta,
                           void*   c, inc_t rs_c, inc_t cs_c,
                           dim_t   l2_num_threads,
                           dim_t   l2_thread_id,
                           dim_t   l1_num_threads,
                           dim_t   l1_thread_id
                         );
			   //                           dim_t   l0_thread_id

static FUNCPTR_T GENARRAY(ftypes,gemm_ker_var2_par);


void bli_gemm_ker_var2_par( obj_t*  alpha,
			    obj_t*  a,
			    obj_t*  b,
			    obj_t*  beta,
			    obj_t*  c,
			    gemm_t* cntl,
			    thread_comm_t* l2_comm,
			    thread_comm_t* l1_comm)
{
	num_t     dt_exec   = bli_obj_execution_datatype( *c );

	dim_t     m         = bli_obj_length( *c );
	dim_t     n         = bli_obj_width( *c );
	dim_t     k         = bli_obj_width( *a );

	void*     buf_a     = bli_obj_buffer_at_off( *a );
	inc_t     rs_a      = bli_obj_row_stride( *a );
	inc_t     cs_a      = bli_obj_col_stride( *a );
	inc_t     ps_a      = bli_obj_panel_stride( *a );

	void*     buf_b     = bli_obj_buffer_at_off( *b );
	inc_t     rs_b      = bli_obj_row_stride( *b );
	inc_t     cs_b      = bli_obj_col_stride( *b );
	inc_t     ps_b      = bli_obj_panel_stride( *b );

	void*     buf_c     = bli_obj_buffer_at_off( *c );
	inc_t     rs_c      = bli_obj_row_stride( *c );
	inc_t     cs_c      = bli_obj_col_stride( *c );

	num_t     dt_alpha;
	void*     buf_alpha;

	num_t     dt_beta;
	void*     buf_beta;

	dim_t l2_num_threads = l2_comm->num_threads_in_group;
	dim_t l2_thread_id   = th_thread_id(l2_comm);
	dim_t l1_num_threads;
	dim_t l1_thread_id;
	if (l1_comm) {
	  l1_num_threads = l1_comm->num_threads_in_group;
	  l1_thread_id = th_thread_id(l1_comm);
	}
	else {
	  l1_num_threads = 1;
	  l1_thread_id = 0;
	}


	FUNCPTR_T f;

/*
	// Handle the special case where c and a are complex and b is real.
	// Note that this is the ONLY case allowed by the inner kernel whereby
	// the datatypes of a and b differ. In this situation, the execution
	// datatype is real, so we need to inflate (by a factor of two):
	//  - the m dimension,
	//  - the column stride of c,
	//  - the column stride (ie: the panel length) of a, and
	//  - the panel stride of a.
	if ( bli_obj_is_complex( *a ) && bli_obj_is_real( *b ) )
	{
		m    *= 2;
		cs_c *= 2;
		cs_a *= 2;
		ps_a *= 2;
	}
*/

	// If alpha is a scalar constant, use dt_exec to extract the address of the
	// corresponding constant value; otherwise, use the datatype encoded
	// within the alpha object and extract the buffer at the alpha offset.
	bli_set_scalar_dt_buffer( alpha, dt_exec, dt_alpha, buf_alpha );

	// If beta is a scalar constant, use dt_exec to extract the address of the
	// corresponding constant value; otherwise, use the datatype encoded
	// within the beta object and extract the buffer at the beta offset.
	bli_set_scalar_dt_buffer( beta, dt_exec, dt_beta, buf_beta );

	// Index into the type combination array to extract the correct
	// function pointer.
	f = ftypes[dt_exec];
	
    // Invoke the function.
	f( m,
	   n,
	   k,
	   buf_alpha,
	   buf_a, rs_a, cs_a, ps_a,
	   buf_b, rs_b, cs_b, ps_b,
	   buf_beta,
	   buf_c, rs_c, cs_c,
	   l2_num_threads,
	   l2_thread_id,
	   l1_num_threads,
	   l1_thread_id );
}


#undef  GENTFUNC
#define GENTFUNC( ctype, ch, varname, ukrname ) \
\
void PASTEMAC(ch,varname)( \
                           dim_t   m, \
                           dim_t   n, \
                           dim_t   k, \
                           void*   alpha, \
                           void*   a, inc_t rs_a, inc_t cs_a, inc_t ps_a, \
                           void*   b, inc_t rs_b, inc_t cs_b, inc_t ps_b, \
                           void*   beta, \
                           void*   c, inc_t rs_c, inc_t cs_c, \
                           dim_t   l2_num_threads, \
                           dim_t   l2_thread_id, \
                           dim_t   l1_num_threads, \
                           dim_t   l1_thread_id \
                         ) \
{ \
	/* Temporary buffer for duplicating elements of B. */ \
	ctype           bd[ PASTEMAC(ch,maxkc) * \
	                    PASTEMAC(ch,nr) * \
	                    PASTEMAC(ch,ndup) ] \
	                    __attribute__((aligned(BLIS_STACK_BUF_ALIGN_SIZE))); \
	ctype* restrict bp; \
\
	/* Temporary C buffer for edge cases. */ \
	ctype           ct[ PASTEMAC(ch,mr) * \
	                    PASTEMAC(ch,nr) ] \
	                    __attribute__((aligned(BLIS_STACK_BUF_ALIGN_SIZE))); \
	const inc_t     rs_ct      = 1; \
	const inc_t     cs_ct      = PASTEMAC(ch,mr); \
\
	/* Alias some constants to shorter names. */ \
	const dim_t     MR         = PASTEMAC(ch,mr); \
	const dim_t     NR         = PASTEMAC(ch,nr); \
	const dim_t     NDUP       = PASTEMAC(ch,ndup); \
	const bool_t    DUPB       = NDUP != 1; \
\
	ctype* restrict zero       = PASTEMAC(ch,0); \
	ctype* restrict a_cast     = a; \
	ctype* restrict b_cast     = b; \
	ctype* restrict c_cast     = c; \
	ctype* restrict alpha_cast = alpha; \
	ctype* restrict beta_cast  = beta; \
	ctype* restrict a1; \
	ctype* restrict b1; \
	ctype* restrict c1; \
	ctype* restrict c11; \
	ctype* restrict a2; \
	ctype* restrict b2; \
\
	dim_t           k_nr; \
	dim_t           m_iter, m_left; \
	dim_t           n_iter, n_left; \
	dim_t           i, j; \
	inc_t           rstep_a; \
	inc_t           cstep_b; \
	inc_t           rstep_c, cstep_c; \
\
	/*
	   Assumptions/assertions:
	     rs_a == 1
	     cs_a == GEMM_MR
	     ps_a == stride to next row panel of A
	     rs_b == GEMM_NR
	     cs_b == 1
	     ps_b == stride to next column panel of B
	     rs_c == (no assumptions)
	     cs_c == (no assumptions)
	*/ \
\
	/* If any dimension is zero, return immediately. */ \
	if ( bli_zero_dim3( m, n, k ) ) return; \
\
	/* Clear the temporary C buffer in case it has any infs or NaNs. */ \
	PASTEMAC(ch,set0s_mxn)( MR, NR, \
	                        ct, rs_ct, cs_ct ); \
\
	/* Compute number of primary and leftover components of the m and n
	   dimensions. */ \
	n_iter = n / NR; \
	n_left = n % NR; \
\
	m_iter = m / MR; \
	m_left = m % MR; \
\
	/* Compute the number of elements in B to duplicate per iteration. */ \
	k_nr = k * NR; \
\
	/* Determine some increments used to step through A, B, and C. */ \
	rstep_a = ps_a; \
\
	cstep_b = ps_b; \
\
	rstep_c = rs_c * MR; \
	cstep_c = cs_c * NR; \
\
	/* If the micro-kernel needs elements of B duplicated, set bp to
	   point to the duplication buffer. If no duplication is called for,
	   bp will be set to the current column panel of B for each iteration
	   of the outer loop below. */ \
	if ( DUPB ) bp = bd; \
\
\
	/* Loop over the n dimension (NR columns at a time). */ \
	for ( j = l2_thread_id; j < n_iter; j += l2_num_threads  ) \
	{ \
        b1 = b_cast + j * cstep_b; \
        c1 = c_cast + j * cstep_c; \
\
		/* If duplication is needed, copy the current iteration's NR
		   columns of B to a local buffer with each value duplicated. */ \
		if ( DUPB ) PASTEMAC(ch,dupl)( k_nr, b1, bp ); \
		else        bp = b1; \
\
		/* Initialize our next panel of B to be the current panel of B. */ \
		b2 = b1; \
\
		/* Interior loop over the m dimension (MR rows at a time). */ \
		for ( i = l1_thread_id; i < m_iter; i += l1_num_threads ) \
		{ \
			a1  = a_cast + i * rstep_a; \
			c11 = c1 + i * rstep_c; \
\
			/* Compute the addresses of the next panels of A and B. */ \
			a2 = a1 + rstep_a; \
			if ( i == m_iter - 1 && m_left == 0 ) \
			{ \
				a2 = a_cast; \
				b2 = b1 + cstep_b; \
				if ( j == n_iter - 1 && n_left == 0 ) \
					b2 = b_cast; \
			} \
\
			/* Invoke the gemm micro-kernel. */ \
			PASTEMAC(ch,ukrname)( k, \
			                      alpha_cast, \
			                      a1, \
			                      bp, \
			                      beta_cast, \
			                      c11, rs_c, cs_c, \
			                      a2, b2); \
\
            /*if( ++n_ukernels == 4 ){ \
                bli_barrier( comm ); \
                n_ukernels = 0; \
            }\*/\
		} \
\
		/* Bottom edge handling. */ \
		if ( m_left && !l1_thread_id ) \
		{ \
			a1  = a_cast + m_iter * rstep_a; \
			c11 = c1 + m_iter * rstep_c; \
			/* Compute the addresses of the next panels of A and B. */ \
			a2 = a_cast; \
			b2 = b1 + cstep_b; \
			if ( j == n_iter - 1 && n_left == 0 ) \
				b2 = b_cast; \
\
\
			/* Invoke the gemm micro-kernel. */ \
			PASTEMAC(ch,ukrname)( k, \
			                      alpha_cast, \
			                      a1, \
			                      bp, \
			                      zero, \
			                      ct, rs_ct, cs_ct, \
			                      a2, b2); \
\
			/* Scale the bottom edge of C and add the result from above. */ \
			PASTEMAC(ch,xpbys_mxn)( m_left, NR, \
			                        ct,  rs_ct, cs_ct, \
			                        beta_cast, \
			                        c11, rs_c,  cs_c ); \
\
		} \
\
	} \
\
	if ( n_left && !l2_thread_id ) \
	{ \
		b1 = b_cast + n_iter * cstep_b; \
		c1 = c_cast + n_iter * cstep_c; \
\
		/* If duplication is needed, copy the n_left (+ padding) columns
		   of B to a local buffer with each value duplicated. */ \
		if ( DUPB ) PASTEMAC(ch,dupl)( k_nr, b1, bp ); \
		else        bp = b1; \
\
		/* Initialize our next panel of B to be the current panel of B. */ \
		b2 = b1; \
\
		/* Right edge loop over the m dimension (MR rows at a time). */ \
		for ( i = l1_thread_id; i < m_iter; i += l1_num_threads ) \
		{ \
			a1  = a_cast + i * rstep_a; \
			c11 = c1 + i * rstep_c; \
\
			/* Compute the addresses of the next panels of A and B. */ \
			a2 = a1 + rstep_a; \
			if ( i == m_iter - 1 && m_left == 0 ) \
			{ \
				a2 = a_cast; \
				b2 = b_cast; \
			} \
\
			/* Invoke the gemm micro-kernel. */ \
			PASTEMAC(ch,ukrname)( k, \
			                      alpha_cast, \
			                      a1, \
			                      bp, \
			                      zero, \
			                      ct, rs_ct, cs_ct, \
			                      a2, b2 ); \
\
			/* Scale the right edge of C and add the result from above. */ \
			PASTEMAC(ch,xpbys_mxn)( MR, n_left, \
			                        ct,  rs_ct, cs_ct, \
			                        beta_cast, \
			                        c11, rs_c,  cs_c ); \
\
		} \
\
		/* Bottom-right corner handling. */ \
		if ( m_left && !l1_thread_id ) \
		{ \
			a1  = a_cast + m_iter * rstep_a; \
			c11 = c1 + m_iter * rstep_c; \
			/* Compute the address of the next panel of A. */ \
			a2 = a_cast; \
			b2 = b_cast; \
\
			/* Invoke the gemm micro-kernel. */ \
			PASTEMAC(ch,ukrname)( k, \
			                      alpha_cast, \
			                      a1, \
			                      bp, \
			                      zero, \
			                      ct, rs_ct, cs_ct, \
			                      a2, b2 ); \
\
			/* Scale the bottom-right corner of C and add the result from above. */ \
			PASTEMAC(ch,xpbys_mxn)( m_left, n_left, \
			                        ct,  rs_ct, cs_ct, \
			                        beta_cast, \
			                        c11, rs_c,  cs_c ); \
		} \
	} \
\
/*PASTEMAC(ch,fprintm)( stdout, "gemm_ker_var2_par: b1", k, NR, b1, NR, 1, "%4.1f", "" ); \
PASTEMAC(ch,fprintm)( stdout, "gemm_ker_var2_par: bd", k, NR*NDUP, bp, NR*NDUP, 1, "%4.1f", "" );*/ \
/*PASTEMAC(ch,fprintm)( stdout, "gemm_ker_var2_par: a1", MR, k, a1, 1, MR, "%4.1f", "" );*/ \
}

INSERT_GENTFUNC_BASIC( gemm_ker_var2_par, GEMM_UKERNEL )

