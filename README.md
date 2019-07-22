# VoigtFitting
Custom Voigt fit using Lorentzian and Gaussian fits for spectroscopy experiments.

Previously static values and properties of the fit were recalculated with every iteration of fitting. This decreased preformace as complexe calculations were unnessesarily repeated.

I created a ROOT TF1 fitting object so that values can be passed to the object for calculation and then stored in the object, so the values can then easily be retrived on subsequent runs.

