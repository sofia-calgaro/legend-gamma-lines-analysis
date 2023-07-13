Load the macros like in "Load.C" in order to use the objects in a root session.

The fit functions are defined in "GammaFitFunctions.h" (the `GammaFitFunction` is a wrapper for TF1 functions that also contains the actual function implementation, like `G`ammaProtoGauss`_` with some added functionality, like parameter naming etc.).

The fit can be done with the `GammaLineFit` object: it needs to be created with a name, a TH1D histogram, and the expected resolution as value or TF1 resolution curve. The fit is then performed using `GammaLineFit::Fit` where one needs to define the expected line positions, the range, which type of background (if I remember a step like background is not really working), and the MCMC precision. The fit is then performed in steps, first it fits only the background range with the background model, then it fits the full model with fixed position/resolution, this fits are simple root fits and used to define the parameter ranges for the bat fit, and only then the bayesian bat fit is preformed using the BCHistogramFitter .

A possible example is the following:
``` cpp
GammaLineFit fitter(TString, TH1, TF1 or double); //name, histo, resolution curve or resolution value
fitter.Fit("K42_1525", {1524.7}, {1505.,1545.}, kLinear,   0.5, 0.2, BCEngineMCMC::kMedium); //name, nominal line pos, range, background model (do not use kStep), prior width on pos, prior width on resolution, MCMC precision
```

The code should also create a lot of output (a log file with the result of each fit, a png picture of the fit, a root file with the data and best fit model, and a pdf with all the posterior distributions).


# Code

## "Load.C"
Macro to _automatically_ run a gamma analysis.

## "GammaFitFunctions.h"
Defines a list of classes, each one for a specific fit function. Available classes are:
- `GammaProtoFunction`
- `GammaProtoGauss` (inheriting from `GammaProtoFunction`)
- `GammaProtoPolynom` (inheriting from `GammaProtoFunction`)
- `GammaProtoFunctionSum` (inheriting from `GammaProtoFunction`)
- `GammaProtoInterruptedFunction` (inheriting from `GammaProtoFunction`)
- `GammaFitFunction`

## "GammaLineFit.cxx"
Macro to perform the fit. **`BAT` dependency is needed!** The main function is `int GammaLineFit::Fit( TString name, vector<double> lines, pair<double,double> range, GammaBackground backgroundType, double linePosPrior, double fwhmPrior, BCEngineMCMC::Precision precision )`:
1. it calls the first function `int GammaLineFit::Fit( vector<double> lines, pair<double,double> range, GammaBackground backgroundType, double linePosPrior, double fwhmPrior, BCEngineMCMC::Precision precision )`
2. here, it creates proto fit function
3. it creates proto background function
4. it collects proto functions
5. it creates functions
6. it creates fit histo
7. it evaluates preliminary background parameters
8. it evaluates preliminary line intensities
9. it sets parameter ranges using preliminaries
10. it creates hist fitter and set priors
11. it sets mc options
12. it performs the fit

The macro provide some log output with best fit parameters and errors. The fit is saved to output files. Marginalized posteriors are saved in a pdf too.




## "runGammaAnalysis.cxx"
Macro to _automatically_ run a series of gamma fits. Just run the function `runGammaAnalysis`. This will call the function `GammaAnalysis` that inspects K lines, Co60 lines, Th-chain lines, U-chain lines, overlapped lines (i.e. e+e- annihilation & Kr85 lines, Pb212 & Pb214, Ac228, Pb214)

## "extractSpectra.cxx"
Macro to extract spectra, with cuts already applied too (ie isTP, isBL, isMuVetoed, multiplicity, is LArVetoed). It returns the exposure values too. **`CLEHP` dependency is needed!**