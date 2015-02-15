/***********************************************************************
*                        TWeightGaussPt                                
*                                                                                                                           
* Returns weight corresponding to gaussian distribution of transverse momenta
* of final state particles.                                               
*		       
*   									       
************************************************************************/

#ifndef TWeightGaussPt_H
#define TWeightGaussPt_H


/*!
  @class TWeightGaussPt

  @brief  Returns weight corresponding to gaussian distribution of transverse momenta of final state particles from the paper [1].

  Weight exp( -\alpha \sum_{i} (p^{T})^{2} )
  
  [1] S. Jadach ,,Rapidity Generator for Monte-Carlo Calculations of Cylindrical Phase Space'' Comp. Phys. Comm. 9 (1975) 297
  

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


////////////////////////////////////////////////////////////////////////
class TWeightGaussPt: public TWeightStrategy
{
private:
   
	///Gauss alpha parameter
	double alpha;
    
    ///Norm of the Gauss
    double norm;
      
    /// acceptance cut
    double ypicut;
    
    /// acceptance cut
    double dzcut;

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
    
    /// number of outgoing particles (2 protons + ( nop-2 ) pions)
    int nop;
    
    /// final particles four-vectors
    TLorentzVector * pf;
       
	///configuration file reader
	TConfigReader * ConfigReader;

	/// Particle Database PDG
	TDatabasePDG * PDGDatabese;

	///reads configuration form the file
	void ReadConfigFile( const string & filename );
    
    ///Sets event 
	void SetEvent( TEvent * event, double evenWeight = 1.0 );
   
	/// calculate weight based on Cylindrical Phase Space and scattering amplitude
	double GetWeight( double eventWeight = 1.0 );
	
	/// Logging
	TLog * Logger;
    
public:
    
    ///Constructor
    TWeightGaussPt(void);                    
    
    ///Destructor
    virtual ~TWeightGaussPt(void);           
  
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
