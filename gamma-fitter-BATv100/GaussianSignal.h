#include <BAT/BCModel.h>
#include <BAT/BCMath.h>
#include <BAT/BCMTF.h>
#include <BAT/BCDataSet.h>
#include <BAT/BCDataPoint.h>

#include <TFile.h>
#include <TMath.h>

#include <string>
#include <vector>

#include "Operations.h"



// ---------------------------------------------------------------------------------------------- Gaus + pol1
#ifndef __BAT__GAUSSIANSIGNAL__H
#define __BAT__GAUSSIANSIGNAL__H

class GaussianSignal : public BCModel
{

public:

    // Constructor
    GaussianSignal(const std::string& name, std::vector<double> bin_content, double E0, double xL, double xR, double E1, double E2, int outputK, int *rng);

    // Destructor
    ~GaussianSignal();

    // Overload LogLikelihood to implement model
    double LogLikelihood(const std::vector<double>& pars);

};
#endif
