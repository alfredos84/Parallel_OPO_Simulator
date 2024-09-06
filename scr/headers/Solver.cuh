/*---------------------------------------------------------------------------*/
// * This file contains functions to solve the Split-Step Fourier method (SSMF)
// * needed to calculate the electric fields evolution through the nonlinear Cr.
// * 
// * In particular, this file should be used when only three equation describes the 
// * problem, i.e., sum or difference frequency generation (SFG or DFG).
// *
// * For any specific process, please check the form of the Couple wave equations
// * in the first function called dAdz(). 
/*---------------------------------------------------------------------------*/



#ifndef _SOLVERCUH
#define _SOLVERCUH


/** Size-scaling vectors after Fourier transform */
__global__ void realScaling ( complex_t *A )
{
	real_t s = static_cast<real_t>(SIZE);
	
	uint idx = threadIdx.x + blockIdx.x * blockDim.x;
	
	if (idx < SIZE){
		A[idx] = A[idx] / s ;
	}
	
	return ;
}



/** Scales all vectors after Fourier transform */
__global__ void CUFFTScale ( complex_t *Awp, complex_t *Aws, complex_t *Awi )
{
	real_t s = static_cast<real_t>(SIZE);
	
	uint idx = threadIdx.x + blockIdx.x * blockDim.x;
	
	if (idx < SIZE){
		Awp[idx] = Awp[idx] / s ; Aws[idx] = Aws[idx] / s ;	Awi[idx] = Awi[idx] / s ;
	}
	
	return ;
}


/** This function compensates the GVD after a single-pass */
__global__ void GDDKernel(complex_t *A, complex_t *aux, real_t *w, real_t GDD)
{
	
	uint idx = threadIdx.x + blockIdx.x * blockDim.x;
	
	if (idx < SIZE){
		aux[idx] = A[idx] * CpxExp( 0.5*w[idx]*w[idx]*GDD );
	}
	if (idx < SIZE){
		A[idx] = aux[idx];
	}
	
	return ;
}

/** Applies an electro optical modulator to an electric field after one round trip. */
__global__ void EOMKernel(complex_t *A, complex_t *aux, real_t beta, real_t fpm, real_t *t)
{
	
	uint idx = threadIdx.x + blockIdx.x * blockDim.x;
	if (idx < SIZE){
		aux[idx] = A[idx] * CpxExp(beta*sinf(2*PI*fpm*t[idx]));
	}
	if (idx < SIZE){
		A[idx] = aux[idx];}
		
	return ;
}

__global__ void satAbsKernel( complex_t *A_ptr, complex_t *aux_ptr, real_t Isat, real_t alpha_0, real_t sat_len )
{
	uint idx = threadIdx.x + blockIdx.x * blockDim.x;
	
	if (idx < SIZE){
		aux_ptr[idx] = A_ptr[idx] * ( (Isat*(1-sat_len*alpha_0) + CpxAbs2(A_ptr[idx])) / (Isat + CpxAbs2(A_ptr[idx])) );
	}
	if (idx < SIZE){
		A_ptr[idx] = aux_ptr[idx] ;
	}

	return ; 
}


/** Computes the nonlinear part: dA/dz=i.κ.Ax.Ay.exp(i.Δk.L) and saves the result in dAx (x represents the different fields) */
__global__ void dAdz( complex_t *dAp, complex_t *dAs,  complex_t *dAi, complex_t *Ap, complex_t *As, complex_t *Ai, 
	real_t kp, real_t ks, real_t ki, real_t dk, real_t z )
{
	complex_t Im; Im.x = 0; Im.y = 1;
	
	uint idx = threadIdx.x + blockIdx.x * blockDim.x;
	
	if (idx < SIZE){
		dAp[idx]  = Im * kp * As[idx] * Ai[idx] * CpxExp(-0.*dk*z) ;
		dAs[idx]  = Im * ks * Ap[idx] * CpxConj(Ai[idx]) * CpxExp(+dk*0.*z);
		dAi[idx]  = Im * ki * Ap[idx] * CpxConj(As[idx]) * CpxExp(+dk*0.*z);
	}
	
	return ;
}


/** Computes a linear combination Ax + s.kx and saves the result in aux_x */
__global__ void LinealCombination( complex_t *auxp, complex_t *auxs, complex_t *auxi,
								   complex_t *Ap, complex_t *As, complex_t *Ai, 
								   complex_t *kp, complex_t *ks, complex_t *ki, real_t s )
{
	
	uint idx = threadIdx.x + blockIdx.x * blockDim.x;

	if (idx < SIZE){
		auxp[idx] = Ap[idx] + kp[idx] * s;
		auxs[idx] = As[idx] + ks[idx] * s;
		auxi[idx] = Ai[idx] + ki[idx] * s;
	}
	
	return ;
}


/** This kernel computes the final sum after appling the Rounge-Kutta algorithm */
__global__ void rk4(complex_t *Ap, complex_t *As,  complex_t *Ai, complex_t *k1p, complex_t *k1s, complex_t *k1i, 
					complex_t *k2p, complex_t *k2s, complex_t *k2i, complex_t *k3p, complex_t *k3s, complex_t *k3i,
					complex_t *k4p, complex_t *k4s, complex_t *k4i, real_t dz )
{
	
	uint idx = threadIdx.x + blockIdx.x * blockDim.x;
	
	if (idx < SIZE){
		Ap[idx] = Ap[idx] + (k1p[idx] + 2*k2p[idx] + 2*k3p[idx] + k4p[idx]) * dz / 6.0;
		As[idx] = As[idx] + (k1s[idx] + 2*k2s[idx] + 2*k3s[idx] + k4s[idx]) * dz / 6.0;
		Ai[idx] = Ai[idx] + (k1i[idx] + 2*k2i[idx] + 2*k3i[idx] + k4i[idx]) * dz / 6.0;
	}
	
	return ;
}


/** Computes the linear part: Ax = Ax.exp(i.f(Ω)*z), where f(Ω) is a frequency dependant functions
 * including the group velocity and the group velocity dispersion parameters. */
__global__ void LinearOperator(complex_t *auxp, complex_t *auxs, complex_t *auxi, complex_t *Awp, complex_t* Aws, complex_t* Awi,
								real_t *w, real_t vp, real_t vs, real_t vi, real_t b2p, real_t b2s, real_t b2i, 
								real_t b3p, real_t b3s, real_t b3i, real_t alpha_crp,
								real_t alpha_crs, real_t alpha_cri, real_t z)
{
	
	uint idx = threadIdx.x + blockIdx.x * blockDim.x;
	
	real_t attenp = expf(-0.5*alpha_crp*z);
	real_t attens = expf(-0.5*alpha_crs*z);
	real_t atteni = expf(-0.5*alpha_cri*z);
	
	if (idx < SIZE){
		auxp[idx] = Awp[idx] * CpxExp(z*w[idx]*((1/vp-1/vi)+0.5*w[idx]*b2p + w[idx]*w[idx]*b3p/6));
		auxs[idx] = Aws[idx] * CpxExp(z*w[idx]*((1/vs-1/vs)+0.5*w[idx]*b2s + w[idx]*w[idx]*b3s/6));
		auxi[idx] = Awi[idx] * CpxExp(z*w[idx]*((1/vi-1/vs)+0.5*w[idx]*b2i + w[idx]*w[idx]*b3i/6));
	}
	if (idx < SIZE){
		Awp[idx] = auxp[idx] * attenp;
		Aws[idx] = auxs[idx] * attens;
		Awi[idx] = auxi[idx] * atteni;
	}
	
	return ;
}


template<typename Crystal>
class Solver
{	// Difine the class Solver for modelling the fields propagation and heating
public:	

	cVecd_t k1p, k2p, k3p, k4p;
	cVecd_t k1s, k2s, k3s, k4s;
	cVecd_t k1i, k2i, k3i, k4i;
	cVecd_t auxp, auxs, auxi;

	Crystal *Cr;
	EFields<Crystal> *A;
	Cavity<Crystal> *Cav;

	Solver(Crystal *_Cr, EFields<Crystal> *_A) : Cr(_Cr), A(_A)
	{	// Constructor

		k1p.resize(SIZE); k2p.resize(SIZE); k3p.resize(SIZE); k4p.resize(SIZE);
		k1s.resize(SIZE); k2s.resize(SIZE); k3s.resize(SIZE); k4s.resize(SIZE);
		k1i.resize(SIZE); k2i.resize(SIZE); k3i.resize(SIZE); k4i.resize(SIZE);
		auxp.resize(SIZE); auxs.resize(SIZE); auxi.resize(SIZE);

		printf("\nInstance of the class Solver.\n");
	}


	Solver(Crystal *_Cr, Cavity<Crystal> *_Cav, EFields<Crystal> *_A) : Cr(_Cr), Cav(_Cav), A(_A)
	{	// Constructor

		k1p.resize(SIZE); k2p.resize(SIZE); k3p.resize(SIZE); k4p.resize(SIZE);
		k1s.resize(SIZE); k2s.resize(SIZE); k3s.resize(SIZE); k4s.resize(SIZE);
		k1i.resize(SIZE); k2i.resize(SIZE); k3i.resize(SIZE); k4i.resize(SIZE);
		auxp.resize(SIZE); auxs.resize(SIZE); auxi.resize(SIZE);

		printf("\nInstance of the class Solver.\n");
	}
	
	~Solver(){ printf("Solver Destructor.\n"); }

	// Methods definition
	void EOM ( cVecd_t &Ax );
	void GDD( cVecd_t &Ax );
	void satAbs( cVecd_t &Ax );
	void dispersion ( );
	void solverRK4( real_t z );
	void SSFM ( );
	void runOPO(const std::vector<int>& save_roundtrips);
	void runOPO();
	void runSinglePass();

};


// Methods declaration

template<typename Crystal>
void Solver<Crystal>::EOM ( cVecd_t &Ax )
{
	if(Cav->eom){
		dim3 block(BLKX); dim3 grid((SIZE+BLKX-1)/BLKX);

		complex_t *aux_ptr;
		CHECK( cudaMalloc((void **)&aux_ptr, Ax.size()*sizeof(complex_t)) );
		real_t* t_ptr = thrust::raw_pointer_cast( A->t.data() );
		complex_t* A_ptr = thrust::raw_pointer_cast( Ax.data() );
				
		EOMKernel<<<grid,block>>>( A_ptr, aux_ptr, Cav->beta, Cav->fpm, t_ptr );
		CHECK(cudaDeviceSynchronize());
		CHECK(cudaFree(aux_ptr));
	}
	return ;
}


template<typename Crystal>
void Solver<Crystal>::GDD( cVecd_t &Ax )
{
	if(Cav->gdd){
		dim3 block(BLKX);	dim3 grid((SIZE+BLKX-1)/BLKX);
		cufftHandle plan;	cufftPlan1d(&plan, SIZE, CUFFT_C2C, 1);

		complex_t *A_ptr = thrust::raw_pointer_cast(Ax.data());
		real_t* w_ptr = thrust::raw_pointer_cast(A->w.data());
		
		complex_t *Aw_ptr , *aux_ptr;
		CHECK(cudaMalloc((void **)&Aw_ptr, Ax.size()*sizeof(complex_t) ));
		CHECK(cudaMalloc((void **)&aux_ptr, Ax.size()*sizeof(complex_t) ));

		cufftExecC2C(plan, (complex_t *)A_ptr, (complex_t *)Aw_ptr, CUFFT_INVERSE);
		CHECK(cudaDeviceSynchronize());

		realScaling<<<grid,block>>>( Aw_ptr );
		CHECK(cudaDeviceSynchronize());
		
		real_t GDD = -(Cav->gamma)*(Cr->b2s)*(Cr->Lcr);
		GDDKernel<<<grid,block>>>(Aw_ptr, aux_ptr, w_ptr, GDD);
		CHECK(cudaDeviceSynchronize());
		
		cufftExecC2C(plan, (complex_t *)Aw_ptr, (complex_t *)A_ptr, CUFFT_FORWARD);
		CHECK(cudaDeviceSynchronize());
		
		cufftDestroy(plan);
		CHECK(cudaFree(Aw_ptr)); CHECK(cudaFree(aux_ptr));
	}

	return ;
}


template<typename Crystal>
void Solver<Crystal>::satAbs ( cVecd_t &Ax )
{
	if(Cav->satabs){
		dim3 block(BLKX); dim3 grid((SIZE+BLKX-1)/BLKX);

		complex_t *aux_ptr;
		CHECK( cudaMalloc((void **)&aux_ptr, Ax.size()*sizeof(complex_t)) );
		complex_t* A_ptr = thrust::raw_pointer_cast( Ax.data() );
				
		satAbsKernel<<<grid,block>>>( A_ptr, aux_ptr, Cav->Isat, Cav->alpha_0, Cav->sat_len );
		CHECK(cudaDeviceSynchronize());
		CHECK(cudaFree(aux_ptr));
	}
	return ;
}

template<typename Crystal>
void Solver<Crystal>::dispersion ( )
{	// A->Applies the dispersion term to the electric fields
	// Parameters for kernels
	dim3 block(BLKX);	dim3 grid((SIZE+BLKX-1)/BLKX);
	// Set plan for cuFFT //
	cufftHandle plan;	cufftPlan1d(&plan, SIZE, CUFFT_C2C, 1);
	// Define pointers to use them in kernels
	complex_t * Ap_ptr   = thrust::raw_pointer_cast(A->Ap.data());
	complex_t * Awp_ptr  = thrust::raw_pointer_cast(A->Awp.data());
	complex_t * As_ptr   = thrust::raw_pointer_cast(A->As.data());
	complex_t * Aws_ptr  = thrust::raw_pointer_cast(A->Aws.data());
	complex_t * Ai_ptr   = thrust::raw_pointer_cast(A->Ai.data());
	complex_t * Awi_ptr  = thrust::raw_pointer_cast(A->Awi.data());
	real_t    * w_ptr    = thrust::raw_pointer_cast(A->w.data());
	complex_t * auxp_ptr = thrust::raw_pointer_cast(this->auxp.data());
	complex_t * auxs_ptr = thrust::raw_pointer_cast(this->auxs.data());
	complex_t * auxi_ptr = thrust::raw_pointer_cast(this->auxi.data());
	
	real_t vp  = Cr->vp,  vs  = Cr->vs,  vi  = Cr->vs;
	real_t b2p = Cr->b2p, b2s = Cr->b2s, b2i = Cr->b2i; 
	real_t b3p = Cr->b3p, b3s = Cr->b3s, b3i = Cr->b3i;
	real_t alpha_crp = Cr->alpha_crp, alpha_crs = Cr->alpha_crs, alpha_cri = Cr->alpha_cri;
	real_t dz = Cr->dz;

	// Linear operator for dz
	cufftExecC2C(plan, (complex_t *)Ap_ptr, (complex_t *)Awp_ptr, CUFFT_INVERSE);
	CHECK(cudaDeviceSynchronize());
	
	cufftExecC2C(plan, (complex_t *)As_ptr, (complex_t *)Aws_ptr, CUFFT_INVERSE);
	CHECK(cudaDeviceSynchronize());
	
	cufftExecC2C(plan, (complex_t *)Ai_ptr, (complex_t *)Awi_ptr, CUFFT_INVERSE);
	CHECK(cudaDeviceSynchronize());
	
	CUFFTScale<<<grid,block>>>(Awp_ptr, Aws_ptr, Awi_ptr);
	CHECK(cudaDeviceSynchronize());

	LinearOperator<<<grid,block>>>( auxp_ptr, auxs_ptr, auxi_ptr, Awp_ptr, Aws_ptr, Awi_ptr,
									w_ptr, vp, vs, vi, b2p, b2s, b2i, b3p, b3s, b3i, alpha_crp,
									alpha_crs, alpha_cri ,dz);
	CHECK(cudaDeviceSynchronize());
	
	cufftExecC2C(plan, (complex_t *)Awp_ptr, (complex_t *)Ap_ptr, CUFFT_FORWARD);
	CHECK(cudaDeviceSynchronize());
	cufftExecC2C(plan, (complex_t *)Aws_ptr, (complex_t *)As_ptr, CUFFT_FORWARD);
	CHECK(cudaDeviceSynchronize());
	cufftExecC2C(plan, (complex_t *)Awi_ptr, (complex_t *)Ai_ptr, CUFFT_FORWARD);
	CHECK(cudaDeviceSynchronize());
	
	cufftDestroy(plan);
	
	return ;
}


template<typename Crystal>
void Solver<Crystal>::solverRK4( real_t z )
{	
	// A->Applies the Fourh-order Runge-Kutta Method with fixed step size dz
	// This function apply the fourth-order Runge-Kutta method	
	
	// Parameters for kernels
	dim3 block(BLKX);	dim3 grid((SIZE+BLKX-1)/BLKX);

	// Define pointers to use them in kernels
	complex_t * k1p_ptr = thrust::raw_pointer_cast(this->k1p.data());
	complex_t * k2p_ptr = thrust::raw_pointer_cast(this->k2p.data());
	complex_t * k3p_ptr = thrust::raw_pointer_cast(this->k3p.data());
	complex_t * k4p_ptr = thrust::raw_pointer_cast(this->k4p.data());
	
	complex_t * k1s_ptr = thrust::raw_pointer_cast(this->k1s.data());
	complex_t * k2s_ptr = thrust::raw_pointer_cast(this->k2s.data());
	complex_t * k3s_ptr = thrust::raw_pointer_cast(this->k3s.data());
	complex_t * k4s_ptr = thrust::raw_pointer_cast(this->k4s.data());

	complex_t * k1i_ptr = thrust::raw_pointer_cast(this->k1i.data());
	complex_t * k2i_ptr = thrust::raw_pointer_cast(this->k2i.data());
	complex_t * k3i_ptr = thrust::raw_pointer_cast(this->k3i.data());
	complex_t * k4i_ptr = thrust::raw_pointer_cast(this->k4i.data());
	
	complex_t * Ap_ptr  = thrust::raw_pointer_cast(A->Ap.data());
	complex_t * As_ptr  = thrust::raw_pointer_cast(A->As.data());
	complex_t * Ai_ptr  = thrust::raw_pointer_cast(A->Ai.data());

	complex_t * auxp_ptr = thrust::raw_pointer_cast(this->auxp.data());
	complex_t * auxs_ptr = thrust::raw_pointer_cast(this->auxs.data());
	complex_t * auxi_ptr = thrust::raw_pointer_cast(this->auxi.data());
	
	real_t dz = Cr->dz;
	real_t dk = A->dk; 
	real_t kp = A->kp, ks = A->ks, ki = A->ki;

	//k1 = dAdz(kappas,dk,z,A)
	dAdz<<<grid,block>>>( k1p_ptr, k1s_ptr, k1i_ptr, Ap_ptr, As_ptr, Ai_ptr, kp, ks, ki, dk, z );
	CHECK(cudaDeviceSynchronize()); 

	//k2 = dAdz(kappas,dk,z+dz/2,A+k1/2) -> aux = A+k1/2
	LinealCombination<<<grid,block>>>( auxp_ptr, auxs_ptr, auxi_ptr, Ap_ptr, As_ptr, Ai_ptr, k1p_ptr, k1s_ptr, k1i_ptr, 0.5 );
	CHECK(cudaDeviceSynchronize());   
	dAdz<<<grid,block>>>( k2p_ptr, k2s_ptr, k2i_ptr, auxp_ptr, auxs_ptr, auxi_ptr, kp, ks, ki, dk, z+dz/4 );
	CHECK(cudaDeviceSynchronize());

	// k3 = dAdz(kappas,dk,z+dz/2,A+k2/2)
	LinealCombination<<<grid,block>>>( auxp_ptr, auxs_ptr, auxi_ptr, Ap_ptr, As_ptr, Ai_ptr, k2p_ptr, k2s_ptr, k2i_ptr, 0.5 );
	CHECK(cudaDeviceSynchronize());   
	dAdz<<<grid,block>>>( k3p_ptr, k3s_ptr, k3i_ptr, auxp_ptr, auxs_ptr, auxi_ptr, kp, ks, ki, dk, z+dz/4 );
	CHECK(cudaDeviceSynchronize());

	// k4 = dAdz(kappas,dk,z+dz,A+k3)
	LinealCombination<<<grid,block>>>( auxp_ptr, auxs_ptr, auxi_ptr, Ap_ptr, As_ptr, Ai_ptr, k3p_ptr, k3s_ptr, k3i_ptr, 1.0 );
	CHECK(cudaDeviceSynchronize());   
	dAdz<<<grid,block>>>( k4p_ptr, k4s_ptr, k4i_ptr, auxp_ptr, auxs_ptr, auxi_ptr, kp, ks, ki, dk, z+dz/2 );
	CHECK(cudaDeviceSynchronize());

	// A = A + (k1+k2+k3+k4)/6
	rk4<<<grid,block>>>( Ap_ptr, As_ptr, Ai_ptr, k1p_ptr, k1s_ptr, k1i_ptr, k2p_ptr, k2s_ptr, k2i_ptr, 
						k3p_ptr, k3s_ptr, k3i_ptr, k4p_ptr, k4s_ptr, k4i_ptr, dz/2 );
	CHECK(cudaDeviceSynchronize());

	return ;
	
}


template<typename Crystal>
void Solver<Crystal>::SSFM ( )
{
	for (real_t z = 0; z <= Cr->Lcr; z+=(Cr->dz))
	{	
		solverRK4(z); // RK4 in dz/2
		dispersion(); // Dispersion in dz
		solverRK4(z); // RK4 in dz/2
	}

	return ;
}


template<typename Crystal>
void Solver<Crystal>::runOPO(const std::vector<int>& save_roundtrips)
{
	
	std::unordered_set<int> save_set(save_roundtrips.begin(), save_roundtrips.end());

	real_t iStart = seconds();	// Timing code
	
	for( uint r = 0; r < NRT; r++ )
	{

		runOPO_status (r, NRT/10); // print simulation status on screen
		SSFM(); // compute couple-wave equations

		
		if(Cav->eom){EOM(A->As); EOM(A->Ai);}			// intracavity elements: EOM
		if(Cav->satabs){satAbs(A->As); satAbs(A->Ai);}	// saturable absorber SESAM
		if(Cav->gdd){GDD(A->As); GDD(A->Ai);}			// GDD
				
		// Apply boundary conditions to resonant fields	
		if(std::get<1>(A->resFields)){Cav->applyReflectivity(A->As, 0.0f);}
		if(std::get<2>(A->resFields)){Cav->applyReflectivity(A->Ai, 0.0f);}
		
		// save specific round trips
		if (save_set.find(r) != save_set.end()) {
			SaveVectorComplexGPU(A->Ap, "output_pump_"+std::to_string(r));
			SaveVectorComplexGPU(A->As, "output_signal_"+std::to_string(r));
			SaveVectorComplexGPU(A->Ai, "output_idler_"+std::to_string(r));
        }
		
		A->cwField( A->Power );

		if( r%NRT/10 == 0 ){ // Check NaN subrutine
			cVech_t t_vec = A->As;	real_t test_nan =  t_vec[0].x;
			if(isnan(test_nan)){
				std::cout << "Simulation failed. Error in NaN check a r = " << r << std::endl;
				r = NRT;
			}
		}
	}

	real_t iElaps = seconds() - iStart;	// finish timing
	TimingCode( iElaps); // print time
	
	return ;
}


// Overloaded runOPO() without saving data
template<typename Crystal>
void Solver<Crystal>::runOPO() {runOPO(std::vector<int>());}


template<typename Crystal>
void Solver<Crystal>::runSinglePass()
{
	
	double iStart = seconds();	// Timing code

	SSFM();

	double iElaps = seconds() - iStart;	// finish timing
	TimingCode( iElaps); // print time
	
	return ;
}



#endif // -> #ifdef _SOLVERCUH