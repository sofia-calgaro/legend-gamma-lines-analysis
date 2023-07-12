#ifndef GammaFitFunctions_H
#define GammaFitFunctions_H

//root
#include "TMath.h"
#include "TString.h"
#include "TF1.h"

// c++
#include <vector>
#include <array>
#include <iostream>

namespace std {

//! virtual parameter container, operators
class GammaProtoFunction
{
   public:
     GammaProtoFunction() {};
     GammaProtoFunction(vector<TString> parNames) : fParNames(parNames) {};
     virtual ~GammaProtoFunction() {};

     virtual double operator() ( double *x, double *par ) = 0;
     friend ostream& operator<<(ostream& os, const GammaProtoFunction& function) {
       os << "[";
       for(int i=0; i<(int) function.fParNames.size(); i++)
         os << " " << function.fParNames.at(i).Data();
       os << " ]";
       return os;
     };

     int             GetNPar()         { return (int) fParNames.size(); };
     TString         GetParName(int i) { return fParNames.at(i);        };
     vector<TString> GetParNames()     { return fParNames;              };

     void SetParName(int i,TString name) { fParNames.at(i) = name;    };

   protected:
     void AddPar(TString name)           { fParNames.push_back(name); };
     void AddPars(vector<TString> names) { 
       fParNames.insert(fParNames.end(),names.begin(),names.end());
     };

   private:
     vector<TString> fParNames;
};

//! normalized gauss (with fwhm, not sigma)
class GammaProtoGauss : virtual public GammaProtoFunction
{
   public:
     GammaProtoGauss( array<TString,3> names={"mean","fwhm","intensity"} )
     : GammaProtoFunction(vector<TString>(begin(names),end(names))) {}; 
     GammaProtoGauss( int number )
     : GammaProtoGauss({Form("mean%d",number),Form("fwhm%d",number),Form("intensity%d",number)}) {};
     virtual ~GammaProtoGauss() {};

     double operator() ( double *x, double *par ) {
       return par[2]/sqrt(2*TMath::Pi()*pow(par[1]/2.35,2))
           *exp(-0.5*pow(x[0]-par[0],2)/pow(par[1]/2.35,2));
     };
};

//! positive part of n-th order polynomial (with offset)
class GammaProtoPolynom : virtual public GammaProtoFunction
{
   public:
     static GammaProtoPolynom constant(
       TString name="background" ) {
       return GammaProtoPolynom({name},0.);
     };
     static GammaProtoPolynom linear(
       double offset=0., array<TString,2> names={"background","slope"} ) {
       return GammaProtoPolynom(vector<TString>(begin(names),end(names)),offset);
     };
     static GammaProtoPolynom quadratic(
       double offset=0., array<TString,3> names={"background","slope","curvature"} ) {
       return GammaProtoPolynom(vector<TString>(begin(names),end(names)),offset);
     };

   protected:
     GammaProtoPolynom( vector<TString> names, double offset ) 
     : GammaProtoFunction(names), fOrder((int)names.size()-1), fOffset(offset) {};

   public:
     GammaProtoPolynom( int order, double offset=0., vector<double> steps={} )
     : fOrder(order), fOffset(offset), fSteps(steps) {
       fSteps.insert(fSteps.begin(),0.);
       for(int i=0;i<(int)fSteps.size();i++) 
         GammaProtoFunction::AddPar(Form("par0.%d",i));
       for(int i=1;i<=fOrder;i++) 
          GammaProtoFunction::AddPar(Form("par%d",i));
     };
     virtual ~GammaProtoPolynom() {};

     double operator() ( double *x, double *par ) {
       double y=0.;
       for(unsigned int i=0;i<fSteps.size();i++) 
         if(x[0]>fSteps[i]) y=par[i];
       for(int i=1;i<=fOrder;i++) 
         y+=par[fSteps.size()-1+i]*pow(x[0]-fOffset,i);
       return y < 0. ? 0. : y;
     };

     int            GetOrder()  { return fOrder;  };
     double         GetOffset() { return fOffset; };
     vector<double> GetSteps()  { return fSteps;  };

     void SetOffset( double offset) { fOffset=offset; };

   private:
     int            fOrder;
     double         fOffset;
     vector<double> fSteps;
};

//! sum of proto functions
class GammaProtoFunctionSum : virtual public GammaProtoFunction
{
   public:
     GammaProtoFunctionSum() {};
     GammaProtoFunctionSum(vector<GammaProtoFunction*> functions) {
       for(auto kv : functions) AddFunction(kv);
     };
     virtual ~GammaProtoFunctionSum() {};

     double operator() ( double *x, double *par ) {
       double  y=0;
       double* thisPar=par;
       for(auto kv : fFunctions) {
         y+=(*kv)(x,thisPar);
         thisPar=&thisPar[kv->GetNPar()];
       }
       return y;
     };

     vector<GammaProtoFunction*> GetFunctions() { return  fFunctions; };

     void AddFunction( GammaProtoFunction* function ) {
       fFunctions.push_back(function);
       AddPars(function->GetParNames());
     };

   private:
     vector<GammaProtoFunction*> fFunctions;
};

//! protofunction with holes
class GammaProtoInterruptedFunction : virtual public GammaProtoFunction
{
   public:
     GammaProtoInterruptedFunction( GammaProtoFunction* function,
       vector<pair<double,double>> skip={} ) 
     : GammaProtoFunction(*function), fFunction(function), fSkip(skip) {};
     virtual ~GammaProtoInterruptedFunction() {};

     double operator() ( double *x, double *par ) {
       for(auto kv : fSkip) { 
         if(x[0]>kv.first&&x[0]<kv.second) {
           TF1::RejectPoint();
           return 0;
       }}
       return (*fFunction)(x,par);
     };

     GammaProtoFunction*         GetFunction() { return fFunction; };
     vector<pair<double,double>> GetSkip()     { return fSkip;     };

     void SetFunction(GammaProtoFunction* function) { fFunction=function; };
     void SetSkip(vector<pair<double,double>> skip) { fSkip=skip; };

   private:     
     GammaProtoFunction*         fFunction;
     vector<pair<double,double>> fSkip;
};

//! fit function
class GammaFitFunction : public TF1
{
   public:
     GammaFitFunction( TF1* function ) : TF1(*function) {};
     GammaFitFunction( TString name, GammaProtoFunction* protoFunction, pair<double,double> range )
     : TF1(name.Data(),protoFunction,range.first,range.second,protoFunction->GetNPar()) {
       for(int i=0; i<protoFunction->GetNPar();i++)
         TF1::SetParName(i,protoFunction->GetParName(i));
     };
     virtual ~GammaFitFunction() {};
};

} // namespace std

#endif
