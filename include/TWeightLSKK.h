/***********************************************************************
*                         TWeightLSKK                                 
*                                                                           
* Returns weight corresponding to matrix element for process pp-->ppK+K- 
* Lebiedowicz & Szczurek Phys.Rev.D85(2012)014026                                           
*		            
*   									       
***********************************************************************/

#ifndef TWeightLSKK_H
#define TWeightLSKK_H


/*!
  @class TWeightLSKK

  @brief * Returns weight corresponding to matrix element for process pp-->ppK+K- Lebiedowicz & Szczurek Phys.Rev.D85(2012)014026 
  
 */


#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <math.h>
#include <complex>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <assert.h>

#include <TSystem.h>
#include <TROOT.h>
#include <TRandom.h>
#include <TApplication.h>
#include <TSystem.h>
#include <TFile.h>
#include <TObjArray.h>
#include <TLorentzVector.h>
#include <TParticlePDG.h>
#include <TDatabasePDG.h>

#include "Global.h"
#include "TConfigReader.h"
#include "TPolicyReader.h"
#include "TEvent.h"
#include "TWeightStrategy.h"



using namespace std;

////////////////////////////////////////////////////////////////////////
//#define DEBUG
////////////////////////////////////////////////////////////////////////

#define POMERON 0
#define F2	1
#define RHO	2
#define A2	3
#define OMEGA	4



////////////////////////////////////////////////////////////////////////
class TWeightLSKK: public TWeightStrategy
{
private:
  
	///Use reggeized propagators or not
    bool reggeize;
  
     /// regge intercept
    double rintc[5];
    
    /// regge  traj. slope
    double rslope[5];
    
    /// traj. phase 
    complex<double> eta[5];
    
    /// pomeron,f,rho couplings [mb]
    double CNN[5];
    double CKN[5];
    double CKK[5];
    
    /// formfactor slope
    double B0[5];
    
    /// regge scale
    double S0;
    
    /// low mass  cutoff slope
    double aa;
    
    /// low mass threshold cut
    double W0;
    
    /// proton proton total x-section GeV-2
    double sigma_tot;
    
    /// slope of eleastic pp x-section
    double B_elastic;
    
    double Beff;
    
    // GeV
    double Lambda2;
    
    // 1-2 GeV
    double LambdaInt;
    
    // GeV-2
    double alphaKPrim;
   
    /// Imaginary unit 0 + i
    const complex<double> Imag;
   
    /// Pi = 3.14...
    static const double PI;    

    double dazim34,dazim12,coll12;
    
    double * azim;
    double * theta;
    double * pt;
    double * prap;
    double dz1,dz2;
    
	/// four-body event invariants
	double t1,t2,ta,tb,s13,s24,s14,s23,cmass,cmass2;
	
	/// Proton mass
    double mProton;
    
    /// Pion mass
    double mPion;
    
    /// Kaon+- mass
	double mKPlus;
    
    /// number of outgoing particles (2 protons + ( nop-2 ) pions)
    int nop;
    
    /// total energy in CM 
    double tecm ;		

    /// final particles four-vectors
    TLorentzVector * pf;
        
    /// beam 1 && 2 fourvectors 
    TLorentzVector pb1,pb2; 

	///configuration file reader
	TConfigReader * ConfigReader;

	/// Particle Database PDG
	TDatabasePDG * PDGDatabese;

	///reads configuration form the file
	void ReadConfigFile( const string & filename );
      
    ///Calculate matrix element
    complex<double> MatrixElement( void );
    
    ///Sets event 
	void SetEvent( TEvent * event, double evenWeight = 1.0 );
   
	/// calculate weight based on Cylindrical Phase Space and scattering amplitude
	double GetWeight( double eventWeight = 1.0 );
 
 	/// Helper function - normalizes angle to the [-180; 180] range.
	double fdelta(double);
 
 	/// Calculate invariants 2->4 body event.
	void Kinematics2to4(TEvent * event);
    
   
	/// Helper function to calculate scattering amplitude
	double FJT(void);		
  
	/// Complex regge amplitute
	complex<double> reggeamp(int isign,double tt, double ss);
  
	/// Off-shell formfactor. 
	double FFJT(double tt); 

	/// Helper function to calculate ingredients for scattering matrix
	complex <double> A( int particle, double s, double t );
  
	/// Calculates rescattering correction.
	double Rescatt(); 
  
	/// Reggeized propagator - see: Phys. Letters B 701 (2011) 434 - 444:
	double alphaK( double t );
	double betaM( double s );
	double betaR( double s );
	complex<double> PKReggeized( double t, double s, double tt );
	complex<double> ReggeizedKaon( double t, double s, double tt );
  
	/// Logging
	TLog * Logger;
  
public:
    
    ///Constructor
    TWeightLSKK(void);                    
    
    ///Destructor
    virtual ~TWeightLSKK(void);           
  
	/// @returns 0 or 1 in accordance to the fact that event is accepted by detector or not.
	virtual double getAcceptance( TEvent * event );
  
	/// @returns true if event was accepted or false if not.
	virtual bool getIsAccepted( TEvent * event );
 
	/// @returns weight of the event
	virtual double GetWeight( TEvent * event, double eventWeight );
  
	///Calculate matrix element
	virtual complex<double> getMatrixElement( TEvent * event );

};
#endif
