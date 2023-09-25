#ifndef GammaLineFit_H
#define GammaLineFit_H

//root
#include "TH1D.h"
#include "TF1.h"
#include "TGraph.h"
#include "TString.h"
#include "TNamed.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TSystem.h"

// c++
#include <vector>
#include <fstream>

// BAT
#include <BAT/BCModel.h>
#include <BAT/BCFitter.h>
#include <BAT/BCHistogramFitter.h>
#include <BAT/BCParameter.h>
#include <BAT/BCLog.h>
#include <BAT/BCAux.h>
#include <BAT/BCH1D.h>

#include "GammaFitFunctions.h"

namespace std {

enum GammaBackground{ kFlat = 0, kStep = 1, kLinear = 2, kQuadratic = 3 };

class GammaLineFit : public TNamed
{
   public:
     //! con(de)structors
     GammaLineFit( TString name, TH1D* hist, TF1* res, TString outputDir="output" ) 
     : TNamed(name.Data(),""), fHist(hist), fRes(res), fOutputDir(outputDir), fFitPerformed(false) {
       gSystem->Exec(Form("mkdir -p %s",outputDir.Data()));
     };
     GammaLineFit( TString name, TH1D* hist, double res, TString outputDir="output" )
     : GammaLineFit(name,hist,new TF1("res",Form("%f",res),0,10000),outputDir) {};
     ~GammaLineFit() {
       if(fFitPerformed) Reset();
     };    

     //! getters
     TH1D*              GetHist()       { return fHist;       };
     TF1*               GetRes()        { return fRes;        };
     BCHistogramFitter* GetHistFitter() { return fHistFitter; };

     //! setters
     void SetHist( TH1D* hist ) { fHist = hist; };
     void SetRes( TF1* res )    { fRes = res;   };

     //! methods
     int Fit( vector<double> lines, pair<double,double> range,
       GammaBackground backgroundType = kLinear,
       double linePosPrior = 0.2, //keV 
       double fwhmPrior    = 0.2, //keV 
       BCEngineMCMC::Precision precision=BCEngineMCMC::kLow );
     int Fit( TString name, vector<double> lines, pair<double,double> range,
       GammaBackground backgroundType=kLinear,
       double linePosPrior = 0.2, //keV
       double fwhmPrior    = 0.2, //keV
       BCEngineMCMC::Precision precision=BCEngineMCMC::kMedium );

   protected:
     void Reset() {
       for(auto kv : fProtoGarbage) delete kv;
       delete fFitFunction;
       delete fBkgFunction;
       delete fHistFitter;
       //delete fFitHist;
     };

   private:
     //! input
     TH1D* fHist;
     TF1*  fRes;

     //! internals
     GammaFitFunction*  fFitFunction;
     GammaFitFunction*  fBkgFunction;
     TH1D*              fFitHist;
     BCHistogramFitter* fHistFitter;
     bool               fFitPerformed;
     vector<GammaProtoFunction*> fProtoGarbage;

     //! output
     TString fOutputDir;
};

} // namespace std

#endif
