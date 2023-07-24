#include "GaussianSignal.h"
#include "Operations.h"
#include <TMath.h>
#include <BAT/BCTH1Prior.h>
#include <cmath>



// ----------------------------------------------------------------------------------------------------- CONSTRUCTOR
GaussianSignal::GaussianSignal(const std::string& name, std::vector<double> bin_content, double E0, double xL, double xR, double E1, double E2, int outputK, int *rng)
    : BCModel(name)
{
	    int *max_height = FindMaximumSignalHeight( E0, E1, E2, bin_content, xL, xR, outputK);
	    
	    int xE0=0; // E0
	    if ( rng[3]==0 ) { xE0=0; }
	    if ( rng[3]==1 ) { xE0=100; }
	    if ( rng[3]==2 ) { xE0=200; }
	    if ( rng[3]==3 ) { xE0=400; }
	    if ( rng[3]==4 ) { xE0=800; }
	    if ( rng[3]==5 ) { xE0=1000; }
	    if ( rng[3]==6 ) { xE0=1200; }
	    if ( rng[3]==7 ) { xE0=1400; }
	    if ( rng[3]==8 ) { xE0=1600; }
	    if ( rng[3]==9 ) { xE0=1800; }
	    if ( rng[3]==10 ) { xE0=2000; }
	    
	    // 1) Signal yield (index 0)
	    AddParameter("E0_height", 0, max_height[3]+xE0, "", "[events]");
	    GetParameters().Back().SetPriorConstant();

}

// ----------------------------------------------------------------------------------------------------- DESTRUCTOR
GaussianSignal::~GaussianSignal() { }

// ----------------------------------------------------------------------------------------------------- MY MODEL
double GaussianSignal::LogLikelihood(const std::vector<double>& pars)
{

            double LP = 0.;
            
            double E0 = GetDataSet()->GetDataPoint(5201).GetValue(0);
            double xL = GetDataSet()->GetDataPoint(5202).GetValue(0);
            double xR = GetDataSet()->GetDataPoint(5203).GetValue(0);
            int outputK = GetDataSet()->GetDataPoint(5205).GetValue(0);
            
            // start values for other data points
            const int explow  = 5207;
            const int expup   = 6760;
            const int binslow = 8313;
            const int binsup  = 9866;
            const int ctslow = 11419;
            const int ctssup  = 12972;
            
            xL = (floor((xL*2)+0.5)/2);
            xR = (floor((xR*2)+0.5)/2);
            int line_low = 0;
            while ( (GetDataSet()->GetDataPoint(binslow+line_low).GetValue(0))!=(xL+0.5) ) { line_low++; if (line_low==1553) {break;} }
            if ( line_low==1553 ) {
            	line_low = 0;
            	while ( (GetDataSet()->GetDataPoint(binslow+line_low).GetValue(0))!=(xL) ) { line_low++; }
            }
            int line_up = 0;
            while ( (GetDataSet()->GetDataPoint(binsup+line_up).GetValue(0))!=(xL+0.5) ) { line_up++; if (line_up==1553) {break;} }
            if ( line_up==1553 ) {
            	line_up = 0;
            	while ( (GetDataSet()->GetDataPoint(binsup+line_up).GetValue(0))!=(xL) ) { line_up++; }
            }
            //std::cout << xL << "\t" << xL+0.5 << "\n";
            //std::cout << line_low << "\t" << line_up << "\n";
	    
	    int j = 0;
            for ( double i=xL; i<xR; i++ ) {
                    double y_obs = GetDataSet()->GetDataPoint(ctssup+line_up+j).GetValue(0); //GetDataSet()->GetDataPoint(i).GetValue(0); // observed value ( 0 = 1st column )
                    double y_exp = pars[0]*TMath::Gaus(i, E0, FindSigma(E0), true); // modelling of signal 
                    
                    // modelling of bkg (no more gamma peaks to take into account since they're already included in bkg-modelling)
                    if ( i<195 ) y_exp += GetDataSet()->GetDataPoint(explow+line_low+j).GetValue(0); 
                    if ( i>=195 ) y_exp += GetDataSet()->GetDataPoint(expup+line_up+j).GetValue(0); 
                    
                    j += 1;

                    LP += BCMath::LogPoisson(y_obs, y_exp); // log of conditional probability, p(data|pars)         
            }
            return LP;
 }	    
