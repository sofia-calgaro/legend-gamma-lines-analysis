#include "GammaLineFit.h"

namespace std {

int GammaLineFit::Fit( vector<double> lines, pair<double,double> range,
  GammaBackground backgroundType, double linePosPrior, double fwhmPrior,
  BCEngineMCMC::Precision precision )
{
  if(fFitPerformed) Reset();

  //! create proto fit function
  GammaProtoPolynom* protoPolynom;
  switch( backgroundType) {
    case kStep:      protoPolynom = new GammaProtoPolynom(0,0,lines);     break;
    case kLinear:    protoPolynom = new GammaProtoPolynom(1,range.first); break;
    case kQuadratic: protoPolynom = new GammaProtoPolynom(2,range.first); break;
    default:         protoPolynom = new GammaProtoPolynom(0);             break;
  }
  GammaProtoFunctionSum* protoFitFunction = new GammaProtoFunctionSum({protoPolynom});
  for(unsigned int i=0;i<lines.size();i++) protoFitFunction->AddFunction(new GammaProtoGauss(i));

  //! create proto background function
  vector<pair<double,double>> skip;
  for(auto kv : lines) skip.push_back({kv-fRes->Eval(kv)/2.35*3,kv+fRes->Eval(kv)/2.35*3}); //+-3 sigma
  GammaProtoInterruptedFunction* protoBkgFunction = new GammaProtoInterruptedFunction(protoPolynom,skip);

  //! collect proto functions
  fProtoGarbage = protoFitFunction->GetFunctions();
  fProtoGarbage.push_back(protoFitFunction);
  fProtoGarbage.push_back(protoBkgFunction);

  //! create functions
  fFitFunction = new GammaFitFunction("fitfunc",protoFitFunction,range);
  fBkgFunction = new GammaFitFunction("bkgfunc",protoBkgFunction,range);
  int nBkgPars = fBkgFunction->GetNpar();

  //! create fit histo
  pair<int,int> rangeBins = { fHist->GetXaxis()->FindBin(range.first),
                              fHist->GetXaxis()->FindBin(range.second) };
  int nBins = rangeBins.second - rangeBins.first;
  TH1D* fFitHist = new TH1D("fitHist","",nBins,range.first,range.second);
  for(int i=1;i<=nBins;i++) fFitHist->SetBinContent(i,fHist->GetBinContent(rangeBins.first+i-1));
  double binning = fFitHist->GetBinWidth(1); //BAT treats as cts / 1 keV

  //! evaluate preliminary background parameters
  fFitHist->Fit(fBkgFunction,"RL");
 
  //! evaluate preliminary line intensities
  for(int i=0;i<nBkgPars;i++) fFitFunction->FixParameter(i,fBkgFunction->GetParameter(i));
  for(unsigned int i=0;i<lines.size();i++) {
    fFitFunction->FixParameter(nBkgPars+3*i+0,lines.at(i));
    fFitFunction->FixParameter(nBkgPars+3*i+1,fRes->Eval(lines.at(i)));
    fFitFunction->SetParLimits(nBkgPars+3*i+2,0.,max(5.,1.5*( //counts in +-3sigma w/o bkg +50%
      fFitHist->Integral( fFitHist->FindBin(lines.at(i)-fRes->Eval(lines.at(i))/2.35*3) ,
                          fFitHist->FindBin(lines.at(i)+fRes->Eval(lines.at(i))/2.35*3) ) -
      fBkgFunction->Integral( lines.at(i)-fRes->Eval(lines.at(i))/2.35*3 ,
                              lines.at(i)+fRes->Eval(lines.at(i))/2.35*3 ))));
  }
  fFitHist->Fit(fFitFunction,"RL");

  //! set parameter ranges using preliminaries
  for(int i=0;i<nBkgPars;i++) {
    fFitFunction->SetParLimits(i, //+-4*uncertainty
      (fBkgFunction->GetParameter(i)-4*fBkgFunction->GetParError(i))/binning, 
      (fBkgFunction->GetParameter(i)+4*fBkgFunction->GetParError(i))/binning);
  }
  for(unsigned int i=0;i<lines.size();i++) {
    fFitFunction->SetParLimits(nBkgPars+3*i+0, //+-4*prior width
     lines.at(i)-linePosPrior*3 ,
     lines.at(i)+linePosPrior*3 );
    fFitFunction->SetParLimits(nBkgPars+3*i+1, //+-4*prior width
     fRes->Eval(lines.at(i))-fwhmPrior*3, 
     fRes->Eval(lines.at(i))+fwhmPrior*3 );
    fFitFunction->SetParLimits(nBkgPars+3*i+2, //preliminary parameters +-4*uncertainty
     max(0.,(fFitFunction->GetParameter(nBkgPars+3*i+2)-
     4*fFitFunction->GetParError(nBkgPars+3*i+2))/binning),
     (fFitFunction->GetParameter(nBkgPars+3*i+2)+
     4*fFitFunction->GetParError(nBkgPars+3*i+2))/binning); 
  }

  //! create hist fitter and set priors
  fHistFitter = new BCHistogramFitter(fFitHist,fFitFunction);
  for(unsigned int i=0;i<lines.size();i++) {
    fHistFitter->SetPriorGauss(nBkgPars+3*i+0,lines.at(i),linePosPrior);
    fHistFitter->SetPriorGauss(nBkgPars+3*i+1,fRes->Eval(lines.at(i)),fwhmPrior);
  };

  //! set mc options
  fHistFitter->SetOptimizationMethod(BCModel::kOptMinuit);
  fHistFitter->MCMCSetPrecision(precision);
  for(unsigned int i=0;i<lines.size();i++) {  
    fHistFitter->GetParameter(nBkgPars+3*i+0)->SetNbins(100);   
    fHistFitter->GetParameter(nBkgPars+3*i+1)->SetNbins(100);   
    fHistFitter->GetParameter(nBkgPars+3*i+2)->SetNbins(100);
  }
  for(int i=0;i<nBkgPars;i++) fHistFitter->GetParameter(i)->SetNbins(100);

  //! perform fit
  BCLog::SetLogLevel(BCLog::error);
  BCLog::OpenLog(Form("%s/%s.bat.log",fOutputDir.Data(),GetName()));
  BCAux::SetStyle();
  fHistFitter->Fit();
  BCLog::CloseLog();

  fFitPerformed = true;
  return 0;
}

int GammaLineFit::Fit( TString name, vector<double> lines, pair<double,double> range,
  GammaBackground backgroundType, double linePosPrior, double fwhmPrior,
  BCEngineMCMC::Precision precision ) 
{
  Fit(lines,range,backgroundType,linePosPrior,fwhmPrior,precision);

  //! provide log output
  vector<double> bestFitPars      = fHistFitter->GetBestFitParameters();
  vector<double> bestFitParErrors = fHistFitter->GetBestFitParameterErrors();
  int nBkgPars = fBkgFunction->GetNpar();
  ofstream log;
  log.open(Form("%s/%s.gamma.log",fOutputDir.Data(),GetName()), ios::app);
  //-----------------------------------------------------------------------------------------------
   log << "------------------------------------------------------------"                   << "\n";
   log << setw(60)          << setfill('*') << Form(" %s *",name.Data())                   << "\n"; 
   log << "-------------- " << setfill(' ')                                                << "\n";
   log << Form(" range:%27s",Form("[%.1f,%.f] keV",range.first,range.second))              << "\n";
  for(int i=0;i<fBkgFunction->GetNpar();i++)                               // background parameters
   log << setw(15) << left  << Form(" %s:",fBkgFunction->GetParName(i))
       << setw(15) << right << Form("(%.1f+-%.1f)",bestFitPars[i],bestFitParErrors[i])
       << " cts" << (i>0&&backgroundType!=kStep?"/keV":"") << (i>1?Form("^%d",i):"")       << "\n";
  for(unsigned int i=0;i<lines.size();i++) {                                     // line parameters
   log << "--------------"                                                                 << "\n";
   log << setw(15) << left  << Form(" %s:",fFitFunction->GetParName(nBkgPars+3*i+0))
       << setw(15) << right << Form("(%.1f+-%.1f)",bestFitPars[nBkgPars+3*i+0],bestFitParErrors[nBkgPars+3*i+0])
       << Form(" keV (nom. %.1f keV)",lines.at(i))                                         << "\n";
   log << setw(15) << left  << Form(" %s:",fFitFunction->GetParName(nBkgPars+3*i+1))
       << setw(15) << right << Form("(%.2f+-%.2f)",bestFitPars[nBkgPars+3*i+1],bestFitParErrors[nBkgPars+3*i+1])
       << Form(" keV (nom. %.1f keV)",fRes->Eval(lines.at(i)))                             << "\n";
   log << setw(15) << left << Form(" %s:",fFitFunction->GetParName(nBkgPars+3*i+2));
  //************************************************************************
    BCH1D* intensity = fHistFitter->GetMarginalized(nBkgPars+3*i+2);
    double mode = intensity->GetMode();
    double low, high;
    intensity->GetSmallestInterval(low,high,0.682689);
  //************************************************************************
    if(low>0) {                                                                        // intensity
   log << Form("(%*.1f [%*.1f,%*.1f])",5,mode,5,low,5,high) << " cts (mode [68%])"         << "\n";
    } else                                                                           // upper limit
   log << Form(" <%*.1f",5,intensity->GetQuantile(0.90)) << " cts (90%)"                   << "\n";
   log << "-------------- " << setfill(' ')                                                << "\n";
  }                                                                                      // p-value
   log << " p-value:      " << setw(15) << fHistFitter->GetPValue()                        << "\n";    
   log << "------------------------------------------------------------"                   << "\n";
  //-----------------------------------------------------------------------------------------------
  log.close();

  //! save fit to files
  TCanvas canvas;
   fHistFitter->DrawFit("HIST",false);
   canvas.Print(Form("%s/%s.%s.png",fOutputDir.Data(),GetName(),name.Data()));
  TFile* file = new TFile(Form("%s/%s.root",fOutputDir.Data(),GetName()),"UPDATE");
   fHistFitter->GetHistogram()->Write(Form("%s_hist",name.Data()));
   //fHistFitter->GetFitFunction()->Write(Form("%s_bestfit",name.Data()));
   fHistFitter->GetFitFunctionGraph(bestFitPars)->Write(Form("%s_bestfit",name.Data()));
   fHistFitter->GetErrorBand()/*Graph(0.16,0.84)*/->Write(Form("%s_band",name.Data()));
  file->Close();

  //! provide pdf file with posterior distributions
  fHistFitter->PrintAllMarginalized(Form("%s/%s.%s.pdf",fOutputDir.Data(),GetName(),name.Data()));

  return 0;
}

} // namespace
