# VoigtFitting
Custom Voigt fit for machine learning fit algorithms using Lorentzian and Gaussian fits for spectroscopy experiments.


Previously static values and properties of the fit were recalculated with every iteration of fitting. This decreased preformace as complexe calculations were unnessesarily repeated.
Additionally a new class had to be added for every isotope analyzed, this decreases scalability.


I created a ROOT TF1 fitting object so that values can be passed to the object for calculation and then stored in the object. The stored values can then easily be retrived on subsequent runs.
Additionally, it eliminates to need to create a custom fit[mass#].C class for every isotope analyzed. They can now all be created as VoigtFit objects with all of the same variables and information conserved.
Also allows for VoigtFit objects to be saved and loaded from text files for easy use, reproduction, and sharing.

Unfortunatly, due to export procedure constraints, experiment data can not be provided here for fitting.
