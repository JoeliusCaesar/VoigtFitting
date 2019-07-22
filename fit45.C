//Excerpt from fit45.C class

#include "voigt_fits.h"
#include <TH1F.h>

//Make fit algorithm object
        double runInfo[] = {6538,44.955912,9};
        double qNumbers[] = {3.5,1,2,-480,368,-12.6,-61.7};
        VoigtFit* sc_45 = new VoigtFit("Sc45I",runInfo,qNumbers);
//      VoigtFit* sc_45 = new VoigtFit("sc45fit.txt"); //Import info from file

        //Create Fit
asymmvoigt = new TF1("Sc45", sc_45, &VoigtFit::MultiVoigt, fspec->GetBinCenter(fspec->FindFirstBinAbove(0)), fspec->GetBinCenter(fspec->FindLastBinAbove(0)),19, "VoigtFit", "MultiVoigt");
