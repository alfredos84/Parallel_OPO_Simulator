/*---------------------------------------------------------------------------*/
// * This file contains functions to -----
/*---------------------------------------------------------------------------*/



#ifndef _CAVITYCUH
#define _CAVITYCUH


/** Add phase (phase) and mirror losses (R) after a single-pass */
__global__ void roundTripChange(complex_t *A, complex_t *aux, real_t R, real_t phase)
{
	
	uint idx = threadIdx.x + blockIdx.x * blockDim.x;
	
	if (idx < SIZE){
		aux[idx] = (sqrtf(R) * CpxExp(phase)) * A[idx];
	}
	if (idx < SIZE){
		A[idx] = aux[idx];
	}
	
	return ;
}


template<typename Crystal>
class Cavity
{	// Difine the class Cavity....

public:	
	real_t R, Lcav, trt, fsr; // cavity 
	real_t gamma;             // GDD comp
	real_t alpha_0, Isat, sat_len;  // saturable absorber
	real_t beta, fpm;		  // EOM	
	bool gdd, satabs, eom;    // presence of intracavity elements
	Crystal *Cr;

	Cavity( Crystal *_Cr, real_t _Lcav, real_t _R ) : Cr(_Cr), Lcav(_Lcav), R(_R)
	{	// Constructor
		trt = ( Lcav + (Cr->Lcr)*((Cr->ns) - 1) )/C;
		fsr = 1./trt;

		printf("\nInstance of the class Cavity.\n");
	}

	~Cavity()
	{	// Destructor
		printf("Cavity Destructor.\n");
	}

	// Methods definition
	void applyReflectivity(cVecd_t &A, real_t phase);
	void setEOM(real_t _beta, real_t _fpm);
	void setGDD (real_t _gamma);
	void setSatAbs (real_t _Isat, real_t _alpha_0, real_t sat_len);
};


// Methods declaration

template<typename Crystal>
void Cavity<Crystal>::applyReflectivity(cVecd_t &A, real_t phase)
{
	dim3 block(BLKX); dim3 grid((SIZE+BLKX-1)/BLKX);
	
	complex_t * aux_ptr;
	CHECK(cudaMalloc((void **)&aux_ptr, A.size()*sizeof(complex_t) ));
	complex_t * A_ptr = thrust::raw_pointer_cast(A.data());
	roundTripChange<<<grid,block>>>(A_ptr, aux_ptr, this->R, phase);
	CHECK(cudaDeviceSynchronize());
	CHECK(cudaFree(aux_ptr));
	
	return ;
}

template<typename Crystal>
void Cavity<Crystal>::setEOM(real_t _beta, real_t _fpm)
{	
	this->eom = true; 
	this->fpm = _fpm; this->beta = _beta; 
	std::cout << "Set intracavity EOM with \u03B2 = " << this->beta << "- fpm = " << this->fpm << std::endl;
	return ;
}

template<typename Crystal>
void Cavity<Crystal>::setGDD (real_t _gamma)
{
	// set gamma between 0 and 1
	this->gdd = true;
	this->gamma = _gamma;
	std::cout << "Set " << static_cast<int>(100*this->gamma)<< "% intracavity GDD compensation.\n" << std::endl;
	return ;
}

template<typename Crystal>
void Cavity<Crystal>::setSatAbs (real_t _Isat, real_t _alpha_0, real_t _sat_len)
{
	this->satabs = true;
	this->Isat = _Isat; this->alpha_0 = _alpha_0; this->sat_len = _sat_len;
	std::cout << "Set saturable absorber.\n" << std::endl;
	return ;
}


#endif // -> #ifdef _CAVITYCUH