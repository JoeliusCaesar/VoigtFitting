
#include <string>
#include <vector>

using std::string;
using std::vector;

class VoigtFit{

public:
	//isoName: Name of isotope w/ weight and charge, 
	//runInfo: Run number, laser frequency, isotope mass
	//qNumbers: Quantum Numbers and atomic coefficients [I,F,J,A1,A2,B1,B2]
	
	//params: Fitting parameters (member variables below)
	//intensities: peak intensities
	
	//VoigtFit(string isoName, double* runInfo, double* params, double* intensities);
	//VoigtFit(string isoName, double* runInfo, double* params, double* intensities, double* qNumbers);

	VoigtFit(string isoName, double* runInfo, double* qNumbers);
	VoigtFit(string filename);

	~VoigtFit(){
		delete [] peakPositions;
		delete [] peakIntensities;
	}


	double Gaussian(double* x, double* par);
	double Lorentzian(double* x, double* par);
	double Voigt(double* x, double* par);
	double Asymmetric_Voigt(double *x, int i);
	double MultiVoigt(double *x, double *par);

	void SetParams(double*);

	double aCoef(const double, const double, const double);
	double bCoef(const double, const double, const double);
	void deltaE(const double, const double, const double);

	double V_to_Freq(double mass, double laser, double beam_energy, int charge);

	string name_;
	int runNumber_;
	double mass_;
	int numpks_;
	
	double* peakPositions;
	double* peakIntensities;

//member variables
protected:
        double centroid;	//Centroid from where splitting occurs
        double a1;		//Ground state 
        double a2;
        double b1;		//Excited State
        double b2;		
        double fwhm;		//Full width at half max (FWHM)
        double lorentz;		//lorentz fraction
        double asym;		//Assymetry
        double baseline;	//Background
	double qI;
	double qF;
	double qJ;

public:
	double getCentroid(){return centroid;}
	double getA1(){return a1;}
	double getA2(){return a2;}
	double getB1(){return b1;}
	double getB2(){return b2;}
	double getFWHM(){return fwhm;}
	double getLorentz(){return lorentz;}
	double getAsym(){return asym;}
	double getBaseline(){return baseline;}
};
