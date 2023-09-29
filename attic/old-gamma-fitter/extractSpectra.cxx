//CLHEP
#include <CLHEP/Units/SystemOfUnits.h>

//root 
#include "TH1D.h"

using namespace gada;

void extractSpectra(TString keyList, TString rootDir="")
{
  struct dataset {
    TString name;
    set<int> datasetIDs;
    TH1D* hRaw; 
    TH1D* hLAr;
    TH1D* hLArC;
    double exposure;

    dataset(TString thisName,set<int> thisDatasetIDs,double binning=1.0)
    : name(thisName), 
      datasetIDs(thisDatasetIDs), 
      hRaw( new TH1D(Form("h%sRaw", thisName.Data()),"",(int)8000/binning,0,8000)),
      hLAr( new TH1D(Form("h%sLAr", thisName.Data()),"",(int)8000/binning,0,8000)),
      hLArC(new TH1D(Form("h%sLArC",thisName.Data()),"",(int)8000/binning,0,8000)),
      exposure(0.) {};
  };

  vector<dataset> datasets = {dataset("EnrBEGe",{0,4}),dataset("EnrCoax",{1,5})};

  FileMap* fileMap = new FileMap();
   if(rootDir!="") fileMap->SetRootDir(rootDir.Data());
   fileMap->BuildFromListOfKeys(keyList.Data());
  DataLoader* dataLoader = new DataLoader();
   dataLoader->AddFileMap(fileMap);
   dataLoader->BuildTier4();
  TChain* master = dataLoader->GetMasterChain(0);

  vector<int>* datasetID = new vector<int>;
  vector<double>* energy = new vector<double>;
  int isTP, isBL, isMuVetoed, multiplicity, isLArVetoed;
  ULong64_t timestamp;
   master->SetBranchAddress("tier4_all.datasetID",   &datasetID   ); 
   master->SetBranchAddress("tier4_all.energy",      &energy      );
   master->SetBranchAddress("tier4_all.isTP",        &isTP        );
   master->SetBranchAddress("tier4_all.isBL",        &isBL        );
   master->SetBranchAddress("tier4_all.isMuVetoed",  &isMuVetoed  );
   master->SetBranchAddress("tier4_all.multiplicity",&multiplicity);
   master->SetBranchAddress("tier4_all.isLArVetoed", &isLArVetoed );
   master->SetBranchAddress("tier4_all.timestamp",   &timestamp   );

  GERunConfigurationManager runConfManager;
   runConfManager.SetVerbosity(0);
   runConfManager.AllowRunConfigurationSwitch(true);
  GETRunConfiguration* runConf = new GETRunConfiguration();

  int entries = master->GetEntries();
  for(int i=0;i<entries;i++) {
    if(i%100000==0) cout << i << " of " << entries << endl;
    master->GetEntry(i);
    if(isTP) {
      runConf = runConfManager.GetRunConfiguration(timestamp);
      for(int ch=0; ch<runConf->GetNDetectors(); ch++) {
        if(! runConf->IsOn(ch)) continue;
        for(auto& kv : datasets) {
          if(kv.datasetIDs.count(datasetID->at(ch))) {
            kv.exposure += runConf->GetDetector(ch)->GetMass()/CLHEP::kg  
            / runConf->GetPulserRate() / (60.*60.*24.*365.25);
    } } } }
    else if(!isBL&&!isMuVetoed&&multiplicity==1) {
      for(int ch=0; ch<runConf->GetNDetectors(); ch++) {
        if(!energy->at(ch)) continue;
        for(auto kv : datasets) {
          if(kv.datasetIDs.count(datasetID->at(ch))) {
            kv.hRaw->Fill(energy->at(ch));
            if(!isLArVetoed) kv.hLAr->Fill(energy->at(ch));
            else kv.hLArC->Fill(energy->at(ch));
    } } } }
  }

  TFile* file = new TFile("spectra.root","RECREATE");
  for(auto kv : datasets) {
    kv.hRaw->SetTitle(Form("spectrum after AC, exposure = %.3f kg*yr",kv.exposure));
    kv.hLAr->SetTitle(Form("spectrum after LAr, exposure = %.3f kg*yr",kv.exposure));
    kv.hLArC->SetTitle(Form("spectrum in coincidence with LAr, exposure = %.3f kg*yr",kv.exposure));
    kv.hRaw->Write();
    kv.hLAr->Write();
    kv.hLArC->Write();
  }
  file->Close();
}
