/*---------------------------------------------------------------------------*/
// * This file contains the class Twm which models the three-wave mixing procces
// * in a nonlinear Cr
/*---------------------------------------------------------------------------*/


#ifndef _EFIELDSCUH
#define _EFIELDSCUH

#pragma once

///////////////////////////////////////////////////////////////////////////////////////////////////
/** Define an input field (usually the pump field) in the time
 * domain. This function is overloaded and its use depends on
 * the employed regime (cw, nanosecond or user-defined)*/


/** Flips a vector for Fourier transforms */
template<typename T>
void fftshift( T& V_flip, T V )
{
	int i, c = V.size()/2;
	for ( i = 0; i < V.size()/2; i++ ){
		V_flip[i+c] = V[i];
		V_flip[i]   = V[i+c];
	}
	
	return ;
}


template<typename Crystal>
class EFields
{
public:
    cVecd_t Ap, As, Ai, Awp, Aws, Awi;
    rVecd_t t, F, w;

	real_t lp, ls, li, waist, Power;
	real_t np, ns, ni;
	real_t kp, ks, ki;
	real_t dk, dkp;     // mismatch and group-velocity mismatch
	std::tuple <bool, bool, bool> resFields;  // is resonant? <pump,signal,idler>

    // Constructor
    EFields(real_t _lp, real_t _ls, real_t _li, real_t _Power, real_t _waist, Crystal *Cr) :
			lp(_lp), ls(_ls), li(_li), Power(_Power), waist(_waist)
    {
        // Initialization or other constructor logic if needed
		this->Ap.resize(SIZE); this->As.resize(SIZE); this->Ai.resize(SIZE);
		this->Awp.resize(SIZE); this->Aws.resize(SIZE); this->Awi.resize(SIZE);
		this->t.resize(SIZE); this->F.resize(SIZE); this->w.resize(SIZE);
		printf("\nInstance of the class Efields.\n");
		np  = Cr->np; ns = Cr->ns; ni = Cr->ni;
		kp  = 2*PI*Cr->dQ/(np*lp);   // pump   kappa [1/V]
		ks  = 2*PI*Cr->dQ/(ns*ls);   // signal kappa [1/V]
		ki  = 2*PI*Cr->dQ/(ni*li);   // idler  kappa [1/V]
		dk  = 2*PI*( np/lp-ns/ls-ni/li-1/Cr->Lambda ); // mismatch factor
		dkp = 1/(Cr->GV(lp, Cr->T))-1/(Cr->GV(ls, Cr->T));
    }


    // Destructor
	~EFields(){	printf("Efields Destructor.\n"); }


	// Methods definition
	void cwField( real_t Power );
	void gaussianField( real_t peakpower, real_t FWHM );
	void gaussianField( real_t peakpower, real_t FWHM, real_t Chirp );
	void noiseGenerator ( cVecd_t& Vec );
	void setTimeFreqVectors(real_t trt);

};


// Methods declaration

template<typename Crystal>
void EFields<Crystal>::cwField( real_t Power )
{
	real_t Inten = Power/(PI*waist*waist);
	real_t Ap0 = sqrtf(2*Inten/((this->np)*EPS0*C));
	cVech_t A(SIZE);
	for (int i = 0; i < A.size(); i++){
		A[i].x = static_cast<real_t>(Ap0); // cw field
		A[i].y = static_cast<real_t>(0.0f);
	}
	thrust::copy(A.begin(), A.end(), this->Ap.begin());

	return ;
}


template<typename Crystal>
void EFields<Crystal>::gaussianField( real_t peakpower, real_t FWHM )
{
	real_t t0    = FWHM*sqrtf(2)/(2*sqrtf(2*logf(2)));
	real_t Inten = Power/(PI*waist*waist);
	real_t Ap0 = sqrtf(2*Inten/((this->np)*EPS0*C));
	cVech_t A(SIZE);
	for (int i = 0; i < A.size(); i++){
		A[i].x = static_cast<real_t>(Ap0*expf(-powf(this->t[i]/(2*t0),2))); // gaussian field
		A[i].y = static_cast<real_t>(0.0f);
	}
	thrust::copy(A.begin(), A.end(), this->Ap.begin());

	return ;
}


template<typename Crystal>
void EFields<Crystal>::gaussianField( real_t peakpower, real_t FWHM, real_t Chirp )
{
	real_t t0    = FWHM*sqrtf(2)/(2*sqrtf(2*logf(2)));
	real_t Inten = Power/(PI*waist*waist);
	real_t Ap0 = sqrtf(2*Inten/((this->np)*EPS0*C));
	cVech_t A(SIZE);
	for (int i = 0; i < A.size(); i++){
		A[i].x = static_cast<real_t>(Ap0*expf(-0.5*(1+Chirp)*powf(this->t[i]/t0,2))); // gaussian field
		A[i].y = static_cast<real_t>(0.0f);
	}
	thrust::copy(A.begin(), A.end(), this->Ap.begin());

	return ;
}


template<typename Crystal>
void EFields<Crystal>::noiseGenerator ( cVecd_t& Vec )
{	// Noise generator for initial signal/idler vectors 
	uint seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator(seed);
	std::normal_distribution<real_t> distribution(0.0,1.0e-15);
	
	cVech_t A(Vec.size());
	real_t nsx, nsy;    
	for (int i=0; i<Vec.size(); ++i) {
		nsx = distribution(generator); A[i].x = static_cast<real_t>(nsx);
		nsy = distribution(generator); A[i].y = static_cast<real_t>(nsy);
	}
	thrust::copy(A.begin(), A.end(), Vec.begin());

	return ;	
}


template<typename Crystal>
void EFields<Crystal>::setTimeFreqVectors(real_t trt)
{
	std::cout << "\t--->Setting time and frequency vectors...";
	this->t = linspace<decltype(this->t)>( -trt*0.5, trt*0.5, SIZE );
	this->F = linspace<decltype(this->F)>( -0.5*SIZE/trt, +0.5*SIZE/trt, SIZE );
	fftshift<decltype(this->F)>(this->w, this->F) ;
	this->w *= (2.0*PI);
	std::cout << "done.\n" << std::flush;
	return ;
}


#endif // -> #ifdef _EFIELDSCUH