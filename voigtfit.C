
#include "voigtfit.h"
#include <fstream>
/*
This code is a work in progress. Currently for all fittings, there needs to be a custom Fit[Isotope#] function in voigt_fit.h
The purpose of this code is to replace the need to add additional functions for new isotopes.
This should also reduce the run time of the fitting program.

-For example:
-Currently the TF1 in fit45.C is defined as:
specStart = fspec->GetBinCenter(fspec->FindFirstBinAbove(0));
specEnd = fspec->GetBinCenter(fspec->FindLastBinAbove(0));
asymmvoigt = new TF1("Sc45", Sc45, spectStart, specEnd, 19);

-Now it can be defined as:
VoigtFit sc_45 = new VoigtFit(Name&Charge,{Run#,Mass,#ofPeaks},{I,F,J,a1,a2,b1,b2})
i.e.
VoigtFit sc_45 = new VoigtFit("Sc45I",{6538,44.955912,9},{3.5,1,2,368,-480,-61.7,-12.6})

-And run as
asymmvoigt = new TF1("Sc45", sc_45, &VoigtFit::MultiVoight, specStart, specEnd, 19, "VoigtFit", "MultiVoigt")

*/

//VoigtFit::VoigtFit(string isoName, double* runInfo, double* params, double* intensities);

//VoigtFit object consturctor
oigtFit::VoigtFit(string isoName, double* runInfo, double* qNumbers){
	name_ = isoName;		//optional, for later implamentation
	runNumber_ = int(runInfo[0]);	//optional
	mass_ = runInfo[1];
	numpks_ = int(runInfo[2]);
	peakPositions = new double[numpks_];
	peakIntensities = new double[numpks_];	
	
	//Quantum Numbers
	qI = qNumbers[0];
	qF = qNumbers[1];
	qJ = qNumbers[2];
	
	//Initial guesses for atomic coeff
	a1 = qNumbers[3];
	a2 = qNumbers[4];
	b1 = qNumbers[5];
	b2 = qNumbers[6];

	//Sets peak positions
	deltaE(qI,qF,qJ);
}

//Load fit from file
VoigtFit::VoigtFit(string filename){
	//Open and Read Infile
	int l = 0;
	double arr[30];
	ifstream File;
	File.open(filename);
	while(!File.eof()){
		File >> arr[l];
		l++;
	}	
	File.close();
	//Close Infile
	
	name_ = filename;
	runNumber_ = int(arr[0]);
	mass_ = arr[1];
	numpks_ = int(arr[2]);
	peakPositions = new double[numpks_]; //Position of Signal Peaks on X axis
	peakIntensities = new double[numpks_];	//Height each Peak on Y axis
	
	a1 = arr[3];
	a2 = arr[4];
	b1 = arr[5];
	b2 = arr[6];

	qI = arr[7];
	qF = arr[8];
	qJ = arr[9];		
}


double VoigtFit::Gaussian(double* x, double* par){
	double A = par[0];
	double V0 = par[1];
	double gamma = par[2]/2;
	double coeff = A/gamma;
	coeff*=sqrt(2*log(2)/TMath::Pi());
	double ex = -log(2)*((x[0]-V0)/gamma)*((x[0]-V0)/gamma);
	return coeff*exp(ex);
}

double VoigtFit::Lorentzian(double* x, double* par){
	double A = par[0];
	double V0 = par[1];
	double gamma = par[2]/2;
	double num = A/(TMath::Pi()*gamma);
	double den = 1+((x[0]-V0)/gamma)*((x[0]-V0)/gamma);	
	return num/den;
}

double VoigtFit::Voigt(double* x, double* par){
	double L = lorentz*Lorentzian(x,par);
	double G = (1-lorentz)*Gaussian(x,par);
	return L+G;
}

double VoigtFit::Asymmetric_Voigt(double *x, int i){
	double V0 = peakPositions[i] + centroid;
	double gamma = 2*fwhm/(1+exp(asym*(x[0]-V0)));
	double newpar[3];
	newpar[0] = peakIntensities[i];
	newpar[1] = V0;
	newpar[2] = gamma;
	return Voigt(x,newpar)+baseline/double(numpks_);
}


double VoigtFit::MultiVoigt(double* x, double* par){
	SetParams(par);
	double result = 0;
	for(int i=0; i<numpks_; i++){
		result+=Asymmetric_Voigt(x,i);
	}
	return result;
}

void VoigtFit::SetParams(double* par){		
	//index 0 reserved for later use
	fwhm = par[1];		//Full width at half max (FWHM)
	lorentz =  par[2];	//lorentz fraction
	asym = par[3];		//Assymetry
	baseline = par[4];	//Background
	a2 = par[5];		//Excited state A 
	a1 = par[6];		//Ground state A
	b2 = par[7];		//Excited State B
	b1 = par[8];		//Ground state B
	centroid = par[9];	//Centroid from where splitting occurs

	deltaE(qI,qF,qJ);	//Calculate splittings from new a & b values
	for(int i = 0; i<numpks_; i++){
		//cout<<"Setting Peak#"<<i<<" to: "<<par[i+10]<<endl;
		peakIntensities[i] = par[i+10];
	}
}


double VoigtFit::aCoef(const double I, const double F, const double J){
	//uble F = I+J;
	double K = F*(F+1.0)-I*(I+1.0)-J*(J+1.0);
	return (K/2.0);
}
double VoigtFit::bCoef(const double I, const double F, const double J){
	//double F = I+J;
	double K = F*(F+1.0)-I*(I+1.0)-J*(J+1.0);
	return (3.0*K*(K+1.0)-4*I*(I+1.0)*J*(J+1.0))/(8.0*I*(2.0*I-1.0)*J*(2.0*J-1.0));
}


/* This Functions calculates the */
void VoigtFit::deltaE(const double I, const double J1, const double J2){
	//A1 A2 B1 B2
	//5  6  7  8

	//Debug log title
	//cout<<setw(6)<<"Peak#"<<" f1"<<"  f2"<<"  j1"<<" j2  E        atomic coefs-- A1:"<<atomicCoef[5]<<" A2:"<<atomicCoef[6]<<" B1:"<<atomicCoef[7]<<" B2:"<<atomicCoef[8]<<endl;
	int numPos = 0;	
	for(double f1 = I-J1; f1<=(I+J1); f1++){
		for(double f2 = I-J2; f2<=(I+J2); f2++){
			if((f1 >= (f2-1))&&(f1<=(f2+1))){ //if delta j = -1,0,1
				double E1 =aCoef(I,f1,J1)*a1+bCoef(I,f1,J1)*b1;//A1 B1 for grnd state
				double E2 =aCoef(I,f2,J2)*a2+bCoef(I,f2,J2)*b2;//A2 B2 for excit state
				double dE = E2-E1;
				peakPositions[numPos] = dE;
				numPos++;
			//debug	cout<<"Peak#"<<numPos<<" positioned at: "<<dE<<endl;
			//	cout<<setw(6)<<numPos<<" "<<f1<<" "<<f2<<" "<<J1<<"  "<<J2<<"  "<<E1<<" - "<<E2<<" = "<<dE<<endl;
			//	if(numPos > 99){cout<<"Too many peaks!"<<endl;}
			}
		}
	}
	/*
	if(true){
		for(int i=0; i<numPos; i++){
			cout<<"Peak#"<<i<<" positioned at: "<<peakPositions[numPos]<<endl;
		}
	}
	*/
}


double VoigtFit::V_to_Freq(double mass, double laser, double beam_energy, int charge){
//mass in eV, laser in 1/cm, beam_energy in Volts, charge in e
	long double beta = sqrt(1 - 1/(1+charge*charge*TMath::Power((beam_energy)/(mass),2)+2*charge*(beam_energy)/(mass)));
	long double doppler = sqrt((1-beta)/(1+beta));
	long double freqobs = laser*TMath::C()*100*doppler;
	long double beta0 = sqrt(1- 1/(1+pow((REF_VOLTAGE)/(REF_MASS*AMU),2)+2*(REF_VOLTAGE)/(REF_MASS*AMU)));
	long double doppler0 = sqrt((1-beta0)/(1+beta0));
	long double basefreq = LASER_REF*TMath::C()*100*doppler0;
	basefreq = 7.6190545017635E+14;
	//cout<<setprecision(20)<<basefreq<<endl;
	return double((freqobs-basefreq)*1E-6);
	
	/*
	double laser_freq = laser * 100. * TMath::C();
	double gamma = beam_energy/mass + 1;
	double beta = TMath::Sqrt(1-1/TMath::Power(gamma,2));
	double f_doppler = TMath::Sqrt((1-beta) / (1+beta));
	double freq_dop = laser_freq * f_doppler*1E-6;
	
	double gamma0 = 29850.0/mass + 1;
	double beta0 = TMath::Sqrt(1-1/TMath::Power(gamma0,2));
	double f_doppler0 = TMath::Sqrt((1-beta0) / (1+beta0));
	double freq_dop0 = laser_freq * f_doppler0*1E-6;
	return freq_dop-freq_dop0;
	*/
}



