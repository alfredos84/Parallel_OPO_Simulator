/*---------------------------------------------------------------------------*/
// * This file contains four functions that save files in .dat extension
// * 1 - SaveFileVectorReal()       : save CPU real vectors 
// * 2 - SaveFileVectorRealGPU()    : save GPU real vectors
// * 3 - SaveFileVectorComplex()    : save CPU complex vectors
// * 4 - SaveFileVectorComplexGPU() : save GPU complex vectors

// Inputs:
// - Vector   : vector to save (stored on CPU or GPU)
// - N        : vector size
// - Filename : name of the saved file
/*---------------------------------------------------------------------------*/


#ifndef _FILESCUH
#define _FILESCUH

#pragma once


void SaveVectorReal (rVech_t V, std::string Filename)
{
	std::ofstream myfile;
	std::string extension = ".dat";
	myfile.open(Filename+extension);
	for (int i = 0; i < V.size(); i++)
		myfile << std::setprecision(20) << V[i] << "\n";
	myfile.close();
	
	return;
}


void SaveVectorRealGPU (rVecd_t V, std::string Filename)
{
	rVech_t Vh = V;
	SaveVectorReal ( Vh, Filename );
	
	return;
}


void SaveVectorComplex (cVech_t V, std::string Filename)
{
	std::ofstream myfile;
	std::string extension_r = "_r.dat", extension_i = "_i.dat";
	myfile.open(Filename+extension_r);
	for (int iy = 0; iy < V.size(); iy++)
		myfile << std::setprecision(20) << V[iy].x << "\n";
	myfile.close();
	myfile.open(Filename+extension_i);
	for (int iy = 0; iy < V.size(); iy++)
		myfile << std::setprecision(20) << V[iy].y << "\n";
	myfile.close();
	
	return;
	
}


void SaveVectorComplexGPU (cVecd_t V, std::string Filename)
{

	cVech_t Vh = V;
	SaveVectorComplex ( Vh, Filename );
	
	return;
}


#endif // -> #ifdef _FILESCUH