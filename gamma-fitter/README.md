Load the macros like in "Load.C" in order to use the objects in a root session.

The fit functions are defined in "GammaFitFunctions.h" (the `GammaFitFunction` is a wrapper for TF1 functions that also contains the actual function implementation, like `GammaProtoGauss` with some added functionality, like parameter naming etc.).

The fit can be done with the `GammaLineFit` object: it needs to be created with a name, a TH1D histogram, and the expected resolution as value or TF1 resolution curve. The fit is then performed using `GammaLineFit::Fit` where one needs to define the expected line positions, the range, which type of background (if I remember a step like background is not really working), and the MCMC precision. The fit is then performed in steps, first it fits only the background range with the background model, then it fits the full model with fixed position/resolution, this fits are simple root fits and used to define the parameter ranges for the bat fit, and only then the bayesian bat fit is preformed using the BCHistogramFitter .

A possible example is the following:
``` cpp
GammaLineFit fitter(TString, TH1, TF1 or double); //name, histo, resolution curve or resolution value
fitter.Fit("K42_1525", {1524.7}, {1505.,1545.}, kLinear,   0.5, 0.2, BCEngineMCMC::kMedium); //name, nominal line pos, range, background model (do not use kStep), prior width on pos, prior width on resolution, MCMC precision
```

The code should also create a lot of output (a log file with the result of each fit, a png picture of the fit, a root file with the data and best fit model, and a pdf with all the posterior distributions).
