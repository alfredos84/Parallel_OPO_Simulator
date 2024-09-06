#include <sys/time.h>

#ifndef _COMMON_H
#define _COMMON_H

#define CHECK(call)                                                            \
{                                                                              \
    const cudaError_t error = call;                                            \
    if (error != cudaSuccess)                                                  \
    {                                                                          \
        fprintf(stderr, "Error: %s:%d, ", __FILE__, __LINE__);                 \
        fprintf(stderr, "code: %d, reason: %s\n", error,                       \
                cudaGetErrorString(error));                                    \
        exit(1);                                                               \
    }                                                                          \
}

#define CHECK_CUBLAS(call)                                                     \
{                                                                              \
    cublasStatus_t err;                                                        \
    if ((err = (call)) != CUBLAS_STATUS_SUCCESS)                               \
    {                                                                          \
        fprintf(stderr, "Got CUBLAS error %d at %s:%d\n", err, __FILE__,       \
                __LINE__);                                                     \
        exit(1);                                                               \
    }                                                                          \
}

#define CHECK_CURAND(call)                                                     \
{                                                                              \
    curandStatus_t err;                                                        \
    if ((err = (call)) != CURAND_STATUS_SUCCESS)                               \
    {                                                                          \
        fprintf(stderr, "Got CURAND error %d at %s:%d\n", err, __FILE__,       \
                __LINE__);                                                     \
        exit(1);                                                               \
    }                                                                          \
}

#define CHECK_CUFFT(call)                                                      \
{                                                                              \
    cufftResult err;                                                           \
    if ( (err = (call)) != CUFFT_SUCCESS)                                      \
    {                                                                          \
        fprintf(stderr, "Got CUFFT error %d at %s:%d\n", err, __FILE__,        \
                __LINE__);                                                     \
        exit(1);                                                               \
    }                                                                          \
}

#define CHECK_CUSPARSE(call)                                                   \
{                                                                              \
    cusparseStatus_t err;                                                      \
    if ((err = (call)) != CUSPARSE_STATUS_SUCCESS)                             \
    {                                                                          \
        fprintf(stderr, "Got error %d at %s:%d\n", err, __FILE__, __LINE__);   \
        cudaError_t cuda_err = cudaGetLastError();                             \
        if (cuda_err != cudaSuccess)                                           \
        {                                                                      \
            fprintf(stderr, "  CUDA error \"%s\" also detected\n",             \
                    cudaGetErrorString(cuda_err));                             \
        }                                                                      \
        exit(1);                                                               \
    }                                                                          \
}


double seconds() {
    static auto start = std::chrono::high_resolution_clock::now();
    auto now = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = now - start;
    return elapsed.count();
}


void TimingCode( double iElaps)
{
	if ( iElaps < 60. ){std::cout << "\n\nTime elapsed " <<  iElaps << " seconds\n\n " << std::endl;}
	else if ( iElaps >= 60. and iElaps < 3600. ){std::cout << "\n\nTime elapsed " <<  iElaps/60 << " minutes\n\n " << std::endl;}
	else{std::cout << "\n\nTime elapsed " <<  iElaps/3600 << " hours\n\n " << std::endl;}
	
	return ;
}


// Linear spacing for time vectors
template<typename T>
void linspace( T& Vec, real_t xmin, real_t xmax)
{
    uint size = Vec.size();
	for (int i = 0; i < Vec.size(); i++)
		Vec[i] = xmin + i * (xmax - xmin)/(size-1);
	
	return ;
}


template<typename T>
T linspace( real_t xmin, real_t xmax, uint size)
{
    T Vec (size);
	for (int i = 0; i < Vec.size(); i++)
		Vec[i] = xmin + i * (xmax - xmin)/(size-1);
	
	return Vec ;
}


struct PrintComplex
{
    __host__ __device__
    void operator()(const cufftComplex& x) const
    {
        printf("(%f, %f)\n", x.x, x.y);
    }
};


struct PrintReal
{
    __host__ __device__
    void operator()(const real_t& x) const
    {
        printf("%f\n", x);
    }
};


void runOPO_status( uint r, uint print_each )
{
    if( (r%print_each == 0) or (r == NRT-1) )
        std::cout << "# Round trip: " << r << " - Completed " << r*100/NRT << "%" << "\t\r" << std::flush;
    return ;
}

template<typename T>
void printVarOnScreen( std::string text,  T var)
{
    std::cout << text << var << std::endl;
    return ;
}


// Function to verify the availability of a GPU
void is_gpu_available() {
    int dev = cudaGetDevice(&dev);	// Set up device (GPU)
    if(dev>0)
    {cudaDeviceProp deviceProp;
    cudaGetDeviceProperties(&deviceProp, dev);
    std::cout << "\n\nUsing Device " << dev << ": GPU " << deviceProp.name << std::endl;
    cudaSetDevice(dev);}
    return ;
}

#endif // _COMMON_H