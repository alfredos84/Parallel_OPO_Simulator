// * Author name: Alfredo Daniel Sanchez
// * email:       alfredo.daniel.sanchez@gmail.com

#include "headers/Libraries.cuh"		// Required libraries
#include "headers/PackageLibraries.cuh"	// Required package libraries

int main(int argc, char *argv[]){
	
	///////////////////////////////////////////////////////////////////////////////////
	// 1. Build the OPO

	// a. set parameters
	real_t lp=0.532, ls=2*lp, li=ls*lp/(ls-lp), waist=31.0f, reflectivity = 0.98f;
	real_t Lcr = 20.e3f, T = 27.0f, Lambda = 6.97f;


	// b. set crystal
	MgOsPPLT * cr1 = new MgOsPPLT(Lcr, T, Lambda, lp, ls, li);
	// cr1->getCrystalProp();


	// c. set cavity
	Cavity<MgOsPPLT> * cav1 = new Cavity<MgOsPPLT>( cr1, 163.0e3f, reflectivity );
	// cav1->setEOM( 0.8f*PI, cav1->fsr ); 
	cav1->setGDD( 1.0f );

	real_t alphas = 0.5*( (1-cav1->R)+cr1->alpha_crs*cr1->Lcr ), alphai = alphas; 
	real_t Ith   = EPS0*C*cr1->np*cr1->ns*cr1->ni*cr1->ls*cr1->li*powf((1/cr1->dQ/cr1->Lcr/PI),2)*alphas*alphai/8;
	real_t N     = atof(argv[1]);
	printVarOnScreen("N = ", N);
	real_t Power = N*Ith*(waist*waist*PI); 
	
	// real_t Power = 50.0e-3f;
	printVarOnScreen("Power = ", Power);
	
	// cav1->setSatAbs( 0.1f*Ith/(0.5f*EPS0*C*cr1->np), 100*(cr1->alpha_crp), 20.0f );
	
	// d. set electric fields
	EFields<MgOsPPLT> * efs = new EFields<MgOsPPLT>(lp, ls, li, Power, waist, cr1);
	efs->resFields = std::make_tuple(false, true, true); //are (pump,signal,idler) resonant?
	efs->cwField( Power ); efs->noiseGenerator( efs->As ); efs->Ai = efs->As;
	efs->setTimeFreqVectors( cav1->trt );
	// real_t FWHM = 0.2; efs->setTimeFreqVectors( 30*FWHM );	// set time and freq vectors
	// efs->gaussianField( Power, FWHM );	efs->noiseGenerator(efs->As); efs->noiseGenerator(efs->Ai); 


	///////////////////////////////////////////////////////////////////////////////////
	// 2. Run the OPO
	
	Solver<MgOsPPLT> * solv1 = new Solver<MgOsPPLT>(cr1, cav1, efs);
	std::vector<int> save_roundtrips(64); // set of roundtrips to save
	for (uint i = 0; i < save_roundtrips.size(); i++)
		save_roundtrips[i] = NRT - save_roundtrips.size() + i;
	solv1->runOPO(save_roundtrips);
	
	SaveVectorReal(solv1->A->t, "time");
	
	// delete objects
	delete solv1; delete efs, delete cr1, delete cav1;

	///////////////////////////////////////////////////////////////////////////////////
	
	return 0;
	
}