/*---------------------------------------------------------------------------*/
// * This file contains a set of functions based on the 
// * Sellmeier equations for the MgO:PPLN nonlinear crystal and other 
// * properties of the χ⁽²⁾ material. Sellmeier equations from reference 
// * O. Gayer: Temperature and wavelength dependent refractive index 
// * equations for MgO-doped congruent and stoichiometric LiNbO3.
/*---------------------------------------------------------------------------*/

// All the functions have two input arguments:
// *     L: wavelenght in um
// *     T: temperature in degrees


#ifndef _MGOPPLN 
#define _MGOPPLN 

#pragma once


class MgOPPLN
{
public:
	real_t Lcr;			// crystal length [μm]
	real_t T; 			// crystal temperature [ºC]
	real_t Lambda;	   	// grating period for QPM [um]
	real_t lp, ls, li;  // wavelenghts

	// Define constants
	real_t d33;				// Eff. second-order susceptibility (d33) [um/V] [Ref]
	real_t dQ;				// Eff. second-order susceptibility for QPM [um/V]
	real_t alpha_crp; 		// pump linear absorption [1/μm]
	real_t alpha_crs;  		// signal linear absorption [1/μm]
	real_t alpha_cri;  		// idler linear absorption [1/μm]
	real_t beta_crs;		// signal 2-photons absorption [μm/W]
	real_t chi3p;			// χ⁽³⁾ in [um²/V²]
	real_t chi3s;			// χ⁽³⁾ in [um²/V²]
	real_t chi3i;			// χ⁽³⁾ in [um²/V²]  
	real_t np, ns, ni;  	// relevant refractive indexes
	real_t vp, vs, vi;  	// relevant group-velocities
	real_t b2p, b2s, b2i;  	// relevant group-velocities
	real_t b3p, b3s, b3i;  	// relevant group-velocities
	real_t dz; 			   	// step size
	std::string name;		// crystal name


	MgOPPLN(real_t _Lcr, real_t _T, real_t _Lambda, real_t _lp, real_t _ls, real_t _li) :
			Lcr(_Lcr), T(_T), Lambda(_Lambda), lp(_lp), ls(_ls), li(_li)
	{	// Constructor
		printf("\nInstance of the class MgOPPLN.\n\n");

		name = "MgO:PPLN";
		d33 = 25.20e-6;	 dQ = 2.0*d33/PI;
		alpha_crp = 0.025e-4; alpha_crs = 0.002e-4; alpha_cri = 0.002e-4;
		dz = static_cast<real_t> (Lcr/NZ);
		this->Lcr = Lcr; this->T = T; this->Lambda = Lambda;
		this->lp = lp; this->ls = ls; this->li = li;
		
		np = this->n(lp, T); ns = this->n(ls, T); ni = this->n(li, T);
		vp = this->GV(lp, T); vs = this->GV(ls, T); vi = this->GV(li, T);
		b2p = this->GVD(lp, T); b2s = this->GVD(ls, T); b2i = this->GVD(li, T);
		b3p = this->TOD(lp, T); b3s = this->TOD(ls, T); b3i = this->TOD(li, T);
	}
	

	~MgOPPLN(){printf("MgOPPLN Destructor\n");}
	

	/** This function returns the MgO:PPLN extraordinary refractive index */
	__host__ __device__ 
	real_t n(real_t L, real_t T) 
	{
		
		real_t f = (T - 24.5) * (T + 570.82);
		real_t a1 = 5.756;
		real_t a2 = 0.0983;
		real_t a3 = 0.2020;
		real_t a4 = 189.32;
		real_t a5 = 12.52;
		real_t a6 =  1.32e-2;
		real_t b1 =  2.860e-6;
		real_t b2 =  4.700e-8;
		real_t b3 =  6.113e-8;
		real_t b4 =  1.516e-4;
		real_t G1 = a1 + b1*f;
		real_t G2 = a2 + b2*f;
		real_t G3 = a3 + b3*f;
		real_t G4 = a4 + b4*f;
		return sqrtf(G1+G2/(powf(L,2) - powf(G3,2))+G4/(powf(L,2) - powf(a5,2))-a6*L*L);
		
	}


	/** This function is an auxiliary function related with the resonances */
	__host__ __device__ 
	real_t resonances(real_t L, real_t T, int p) 
	{
		
		real_t f = (T - 24.5) * (T + 570.82);
		real_t a1 = 5.756;
		real_t a2 = 0.0983;
		real_t a3 = 0.2020;
		real_t a4 = 189.32;
		real_t a5 = 12.52;
		real_t a6 =  1.32e-2;
		real_t b1 =  2.860e-6;
		real_t b2 =  4.700e-8;
		real_t b3 =  6.113e-8;
		real_t b4 =  1.516e-4;
		real_t G1 = a1 + b1*f;
		real_t G2 = a2 + b2*f;
		real_t G3 = a3 + b3*f;
		real_t G4 = a4 + b4*f;

		return G2/powf((powf(L,2) - powf(G3,2)), p) + G4/powf((powf(L,2) - powf(a5,2)), p);
		
	}


	/** Returns the first-order derivative of the 
	* refractive index respect to the wavelength dn/dλ. */
	__host__ __device__ 
	real_t dndl(real_t L,real_t T) 
	{
		
		real_t f = (T - 24.5) * (T + 570.82);
		real_t a1 = 5.756;
		real_t a2 = 0.0983;
		real_t a3 = 0.2020;
		real_t a4 = 189.32;
		real_t a5 = 12.52;
		real_t a6 =  1.32e-2;
		real_t b1 =  2.860e-6;
		real_t b2 =  4.700e-8;
		real_t b3 =  6.113e-8;
		real_t b4 =  1.516e-4;
		real_t G1 = a1 + b1*f;
		real_t G2 = a2 + b2*f;
		real_t G3 = a3 + b3*f;
		real_t G4 = a4 + b4*f;

		return -L*(resonances(L, T, 2) + a6)/n(L, T);
		
	}


	/** Returns the second-order derivative of the
	* refractive index respect to the wavelength d²n/dλ². */
	__host__ __device__ 
	real_t d2ndl2(real_t L,real_t T) 
	{
		
		real_t f = (T - 24.5) * (T + 570.82);
		real_t a1 = 5.756;
		real_t a2 = 0.0983;
		real_t a3 = 0.2020;
		real_t a4 = 189.32;
		real_t a5 = 12.52;
		real_t a6 =  1.32e-2;
		real_t b1 =  2.860e-6;
		real_t b2 =  4.700e-8;
		real_t b3 =  6.113e-8;
		real_t b4 =  1.516e-4;
		real_t G1 = a1 + b1*f;
		real_t G2 = a2 + b2*f;
		real_t G3 = a3 + b3*f;
		real_t G4 = a4 + b4*f;


		real_t A  = (L*dndl(L,T)/powf(n(L,T),2)-1/n(L,T))*(resonances(L, T, 2) + a6);
		real_t B  = 4*L*L/n(L,T) * resonances(L, T, 3);
		
		return A+B;
		
	}


	/** Returns the third-order derivative of the
	* refractive index respect to the wavelength d³n/dλ³. */
	__host__ __device__ 
	real_t d3ndl3(real_t L,real_t T) 
	{
		
		real_t f = (T - 24.5) * (T + 570.82);
		real_t a1 = 5.756;
		real_t a2 = 0.0983;
		real_t a3 = 0.2020;
		real_t a4 = 189.32;
		real_t a5 = 12.52;
		real_t a6 =  1.32e-2;
		real_t b1 =  2.860e-6;
		real_t b2 =  4.700e-8;
		real_t b3 =  6.113e-8;
		real_t b4 =  1.516e-4;
		real_t G1 = a1 + b1*f;
		real_t G2 = a2 + b2*f;
		real_t G3 = a3 + b3*f;
		real_t G4 = a4 + b4*f;

		real_t A1 = (2*dndl(L,T)+L*d2ndl2(L,T))/powf(n(L,T),2);
		real_t A2 = -2*L*powf(dndl(L,T),2)/powf(n(L,T),3);
		real_t AA = (A1 + A2)*(resonances(L,T,2)+a6);
		real_t B1 = 12*L/n(L,T);
		real_t B2 = - 4*L*L/n(L,T)*d2ndl2(L,T)*(1-1/n(L,T));
		real_t BB = (B1+B2)*resonances(L,T,3);
		real_t CC = -24*L*L*L/n(L,T)*resonances(L,T,4);

		return AA + BB + CC;	
	}


	/** Returns the group-velocity vg(λ) = c/(n(λ)-λdn/dλ). */
	__host__ __device__ 
	real_t GV(real_t L,real_t T) 
	{
		
		return C/(n(L,T)-L*dndl(L,T));
	}


	/** Returns the group-velocity β2(λ)=λ^3/(2πc²)(d²n/dλ²). */
	__host__ __device__ 
	real_t GVD(real_t L,real_t T) 
	{
		return powf(L,3)*d2ndl2(L, T)/(2*PI*C*C);
	}


	/** Returns the TOD β3(λ)=-λ^4/(4π²c³)[3.d²n/dλ² + λ.d³n/dλ³]. */
	__host__ __device__ 
	real_t TOD(real_t L,real_t T) 
	{
		return -powf(L,4)/(4*PI*PI*C*C*C)*(3*d2ndl2(L, T)+L*d3ndl3(L, T));
	}


	void getCrystalProp () 
	{
		std::cout << "Crystal name    = " << name << std::endl;
		std::cout << "Temp            = " << T << " ºC" << std::endl;
		std::cout << "lp              = " << this->lp << " um" << std::endl;
		std::cout << "ls              = " << this->ls << " um" << std::endl;
		std::cout << "li              = " << this->li << " um" << std::endl;
		std::cout << "vgp             = " << this->vp << " um/ps" << std::endl;
		std::cout << "vgs             = " << this->vs << " um/ps" << std::endl;
		std::cout << "vgi             = " << this->vi << " um/ps" << std::endl;
		std::cout << "b2p             = " << this->b2p << " ps²/um" << std::endl;
		std::cout << "b2s             = " << this->b2s << " ps²/um" << std::endl;
		std::cout << "b2i             = " << this->b2i << " ps²/um" << std::endl;
		std::cout << "dQ              = " << dQ*1e6 << " pm/V"  << std::endl;
		std::cout << "\u039B               = " << Lambda << " \u03BCm"  << std::endl;
		std::cout << "\u03B1cp             = " << alpha_crp << " \u03BCm⁻¹"  << std::endl;
		std::cout << "\u03B1cs             = " << alpha_crs << " \u03BCm⁻¹" << std::endl;
		std::cout << "\u03B1ci             = " << alpha_cri << " \u03BCm⁻¹" << std::endl;
		std::cout << "dz              = " << dz << " \u03BCm"  << std::endl;
		std::cout << "Crystal length  = " << Lcr*1e-3 << " mm\n"  << std::endl;
		return ;
		
	}


};

#endif // -> #ifdef _MGOPPLN