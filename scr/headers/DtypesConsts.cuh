#ifndef _DTYPESCONSTSCUH
#define _DTYPESCONSTSCUH

/** 
 * Complex data type: a set of datatypes are
 * defined to make the code more readable.
 *
 * Definitions for numbers
 * real_t : datatype for real numbers
 * complex_t     : datatype for complex numbers
 * 
 * Definitions for vectors:
 * 
 * rVech_t  : real vector host
 * rVecd_t  : real vector device
 * cVech_t  : complex vector host
 * cVecd_t  : complex vector device
 */

using real_t = float;
using complex_t = cufftComplex;
using rVech_t = thrust::host_vector<real_t>;
using rVecd_t = thrust::device_vector<real_t>;
using cVech_t = thrust::host_vector<complex_t>;
using cVecd_t = thrust::device_vector<complex_t>;


// Define global constants
const uint SIZE   = 1 << 14;  // vector size
const uint NZ     = 150;		  // size discretization
const uint NRT    = 10000;		// number of round trips    
const uint BLKX   = 1 << 7;		// block dimensions for kernels

const real_t PI   = 3.14159265358979323846;		  // pi
const real_t C    = 299792458*1E6/1E12;		    	// speed of ligth in vacuum [um/ps]
const real_t EPS0 = 8.8541878128E-12*1E12/1E6;	// vacuum pertivity [W.ps/V²μm] 

#endif // -> #ifdef _DTYPESCONSTS
